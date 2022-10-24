import enum
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, groupby
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, List, Literal, Optional, Tuple

from attrs import define

from .log import get_logger
from .model import AlignInfo, MonomerReadSeq, Walk
from .settings import WALK_TAG
from .utils import AlignCategory, FileCollection, SamFlags

logger = get_logger()

if TYPE_CHECKING:
    pass


@define
class MappingFileCollection(FileCollection):
    bam: Path = Path("{prefix}.bam")
    temp_bam: Path = Path("{prefix}.temp.bam")
    concatemers: Path = Path("{prefix}.concatemers.parquet")
    contacts: Path = Path("{prefix}.contacts.parquet")
    unmapped: Path = Path("{prefix}.unmapped.fastq")


def group_aligns_by_concatemers(
    aligns: Iterable[MonomerReadSeq], sort: bool = True
) -> Iterable[Tuple[str, List[MonomerReadSeq]]]:
    seen = set()
    for concat_id, aligns in groupby(aligns, lambda x: x.concatemer_id):
        if concat_id in seen:
            raise ValueError(
                f"Concatemer '{concat_id}' has already been seen, "
                "these alignments should be sorted by MI tag"
            )
        aligns = list(aligns)
        if sort:
            aligns = sorted(aligns, key=lambda x: x.coords.subread_idx)
        yield (concat_id, aligns)
        seen.add(concat_id)


def pick_walk_aligns(aligns: List[MonomerReadSeq], concat_id: str) -> List[bool]:
    expected_monomers = aligns[0].coords.subread_total
    keep = [False] * len(aligns)

    # pick one alignment per-monomer to use in the walk
    for monomer_id, cats in groupby(
        [(_.monomer_id, x, _.read_seq.flags.category) for x, _ in enumerate(aligns)],
        lambda x: x[0],
    ):
        # sort by category enum (int enum so get the order):
        #  primary,unmapped,secondary,supplementary
        cats = sorted(cats, key=lambda x: x[2])
        for y, (_, idx, cat) in enumerate(cats):
            if y == 0:  # the pass alignment
                aligns[idx].read_seq.flags.qcfail = False
                keep[idx] = True
                if cat not in (AlignCategory.primary, AlignCategory.unmapped):
                    logger.warning(
                        f"Warning: best alignment for monomer: {monomer_id} "
                        f"has category {cat.name}"
                    )
            else:
                aligns[idx].read_seq.flags.qcfail = True
    num_pass = sum(keep)
    if num_pass != expected_monomers:
        logger.warning(
            f"Expected to see {expected_monomers} alignments for "
            f"concatemer {concat_id}, found {num_pass}"
        )
    return keep


def annotate_monomer_alignments(
    mi_sorted_aligns: Iterable[MonomerReadSeq],
    remove_qcfail: bool = True,
) -> Iterable[Tuple[str, List[MonomerReadSeq]]]:
    for concat_id, aligns in group_aligns_by_concatemers(mi_sorted_aligns, sort=True):
        keep = pick_walk_aligns(aligns, concat_id)
        pass_aligns = [a for x, a in enumerate(aligns) if keep[x]]
        walk = Walk.from_aligns(pass_aligns)
        walk_tag = walk.to_tag()
        for a in pass_aligns:
            a.read_seq.tags[WALK_TAG] = walk_tag
        yield (concat_id, pass_aligns if remove_qcfail else aligns)


def sort_aligns_by_concatemer_idx(
    aligns: List[MonomerReadSeq], validate: bool = True
) -> List[MonomerReadSeq]:
    subread_idx, orig_idx = zip(
        *sorted([(a.coords.subread_idx, x) for x, a in enumerate(aligns)])
    )
    new_aligns = [aligns[i] for i in orig_idx]
    if validate:
        expected_subreads = aligns[0].coords.subread_total
        align_counts = Counter()
        for a in new_aligns:
            if a.read_seq.align_info:
                align_counts[
                    SamFlags.from_int(a.read_seq.align_info.flag).category
                ] += 1
            else:
                align_counts["unmapped"] += 1
        if align_counts["primary"] + align_counts["unmapped"] != expected_subreads:
            raise ValueError(
                f"Unexpected number of alignments: {align_counts} {expected_subreads}"
            )
    return new_aligns


def get_pairs(
    sorted_aligns: List[MonomerReadSeq], direct_only: bool = False
) -> Iterable[Tuple[MonomerReadSeq, MonomerReadSeq, "PairData"]]:
    if direct_only:
        pairs = zip(sorted_aligns[:-1], sorted_aligns[1:])
    else:
        pairs = combinations(sorted_aligns, 2)
    for x, (left, right) in enumerate(pairs):
        yield (left, right, PairData.from_monomer_pair(left, right, pair_idx=x))


def set_walk_tag(align: MonomerReadSeq, walk_str: str):
    align.read_seq.tags[WALK_TAG] = f"{WALK_TAG}:Z:{walk_str}"


def set_next_read(left: MonomerReadSeq, right: MonomerReadSeq):
    return
    # if right.ref_name != "*":
    #    if right.ref_name == left.ref_name:
    #        left.next_ref_name = "="
    #    else:
    #        left.next_ref_name = right.ref_name
    #    left.next_ref_pos = right.ref_pos


class PairAlignState(enum.IntEnum):
    both = 0
    left = 1
    right = 2
    neither = 4

    @staticmethod
    @lru_cache(maxsize=256)
    def from_flags(left_mapped: bool, right_mapped: bool) -> "PairAlignState":
        if left_mapped and right_mapped:
            return PairAlignState.both
        elif not (left_mapped or right_mapped):
            return PairAlignState.neither
        elif left_mapped:
            return PairAlignState.left
        elif right_mapped:
            return PairAlignState.right
        else:
            raise ValueError


@dataclass
class PairData:
    concatemer_id: str
    pair_idx: int
    is_direct: bool
    read_distance: int
    align_state: PairAlignState
    is_cis: Optional[bool] = None
    genome_distance: Optional[int] = None
    relative_ori: Optional[bool] = None

    def to_flags(
        self, align_flags: Tuple[SamFlags, SamFlags]
    ) -> Tuple[SamFlags, SamFlags]:
        if self.align_state == PairAlignState.both:
            flags = (
                SamFlags(proper_pair=self.is_cis is True),
                SamFlags(proper_pair=self.is_cis is True),
            )
        elif self.align_state == PairAlignState.left:
            flags = (SamFlags(munmap=True), SamFlags(unmap=True))
        elif self.align_state == PairAlignState.right:
            flags = (SamFlags(unmap=True), SamFlags(munmap=True))
        elif self.align_state == PairAlignState.neither:
            flags = (
                SamFlags(unmap=True, munmap=True),
                SamFlags(unmap=True, munmap=True),
            )
        else:
            raise ValueError(self.align_state)
        for x, f in enumerate(flags):
            f.paired = True
            if x == 0:
                f.mreverse = align_flags[1].mreverse
                f.read1 = True
            else:
                f.mreverse = align_flags[0].mreverse
                f.read2 = True
            flags[x].supplementary = align_flags[x].supplementary
            flags[x].secondary = align_flags[x].secondary
        return flags

        # l_flag, r_flag = left.flags.copy(), right.flags.copy()
        # l_flag.read2, r_flag.read1 = False, False
        # l_flag.read1, r_flag.read2 = True, True
        # l_flag.proper_pair, r_flag.proper_pair = proper_pair, proper_pair

        # if right.align_info:
        #    l_flag.munmap = False
        #    l_flag.mreverse = right.align_info.strand == "-"
        #    if left.align_info is not None:
        #        next_ref = (
        #            "="
        #            if (left.align_info.ref_name == right.align_info.ref_name)
        #            else right.align_info.ref_name
        #        )
        #    else:
        #        next_ref = right.align_info.ref_name
        #    next_pos = right.align_info.ref_pos + 1
        #    self.write_record(
        #        left,
        #        read_name=l_read_name,
        #        flag=l_flag,
        #        next_reference_name=next_ref,
        #        next_reference_start=next_pos,
        #    )
        # else:
        #    l_flag.munmap = True

        # raise NotImplementedError

    @classmethod
    def from_monomer_pair(
        cls, left: MonomerReadSeq, right: MonomerReadSeq, pair_idx: int = -1
    ) -> "PairData":
        concatemer_id = left.concatemer_id
        md_l, md_r = (
            left.coords,
            right.coords,
        )
        is_direct = (md_r.subread_idx - md_l.subread_idx) == 1
        read_distance = md_r.start - md_l.end
        align_state = PairAlignState.from_flags(
            left.read_seq.flags.unmap is False, right.read_seq.flags.unmap is False
        )
        if not (align_state == PairAlignState.both):
            return PairData(
                concatemer_id, pair_idx, is_direct, read_distance, align_state
            )
        else:
            l_align = left.read_seq.align_info
            r_align = right.read_seq.align_info
            is_cis = l_align.ref_name == r_align.ref_name
            if is_cis:
                strands = (
                    l_align.strand,
                    r_align.strand,
                )
                if "." in strands:
                    relative_ori = None
                else:
                    relative_ori = strands[0] == strands[1]
                genome_distance = calculate_genomic_distance(
                    (l_align.strand, l_align.ref_pos, l_align.ref_end),
                    (r_align.strand, r_align.ref_pos, r_align.ref_end),
                )
            else:
                relative_ori = None
                genome_distance = None
            return PairData(
                concatemer_id,
                pair_idx,
                is_direct,
                read_distance,
                align_state,
                is_cis,
                genome_distance,
                relative_ori,
            )


def calculate_genomic_distance(
    left: Tuple[Literal["+", "-"], int, int],
    right: Tuple[Literal["+", "-"], int, int],
):
    l_coord = left[2] if left[0] == "+" else left[1]
    r_coord = right[1] if right[0] == "+" else right[2]
    return r_coord - l_coord


def get_concatemer_align_string(sorted_aligns: List[MonomerReadSeq]) -> str:
    _aligns = [a.read_seq.align_info for a in sorted_aligns]

    def _to_str(a: Optional[AlignInfo]) -> str:
        if a is None or a.flag == 4:
            return "*"
        else:
            return f"{a.ref_name}:{a.ref_pos}_{a.ref_pos + a.length}"

    return ",".join([_to_str(_) for _ in _aligns])


# TODO: delete this once we're sure we don't need it
# def map_concatemer_read(
#    *,
#    aligner: mp.Aligner,
#    read: MonomerReadSeq,
#    cutter: Cutter,
#    thread_buf: Optional[mp.ThreadBuffer] = None,
# ) -> Tuple[MonomerReadSeq, List[ReadFragment], List[mp.Alignment]]:
#
#    read_frags = sequence_to_read_fragments(cutter, read, store_sequence=False)
#    hits = [
#        list(aligner.map(read.seq[_.read_start : _.read_end], buf=thread_buf))
#        for _ in read_frags
#    ]
#    return (read, read_frags, hits)
