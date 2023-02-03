"""Aligns."""
from collections import Counter
from dataclasses import dataclass
import enum
from functools import lru_cache
from itertools import combinations, groupby
from typing import Iterable, List, Literal, Optional, Tuple

from pore_c_py.log import get_logger
from pore_c_py.model import AlignInfo, MonomerReadSeq, Walk
from pore_c_py.sam_utils import AlignCategory, MOLECULE_TAG, SamFlags, WALK_TAG

logger = get_logger()

PairedMonomers = Tuple[
    Optional[MonomerReadSeq], Optional[MonomerReadSeq], Optional["PairData"]
]


def is_colinear(
        align1: Optional[AlignInfo],
        align2: Optional[AlignInfo], tol: int = 1):
    """Check if two alignments are co-linear in the genome."""
    # if either is is unmapped then they can't be co-linear
    if (
        align1 is None
        or align2 is None
        or align1.sam_flags.unmap
        or align2.sam_flags.unmap
    ):
        return False
    if align1.ref_name != align2.ref_name:
        return False
    if align1.strand != align2.strand:
        return False
    delta = min(align1.ref_end, align2.ref_end) - max(align1.ref_pos, align2.ref_pos) # noqa
    # overlaps
    if delta > 0:
        return True
    # doesn't overlap
    if delta < tol:
        return True
    return False


def group_colinear(aligns: List[AlignInfo], tol: int = 1) -> List[List[int]]:
    """Group alignments into co-linear blocks."""
    if len(aligns) < 2:
        return [list(range(len(aligns)))]
    res = []
    block = []
    for x, a in enumerate(aligns):
        if x == 0:
            block = [x]
        elif is_colinear(aligns[block[-1]], a, tol=tol):
            block.append(x)
        else:
            res.append(block)
            block = [x]
    if block:
        res.append(block)
    return res


def group_aligns_by_concatemers(
    aligns: Iterable[MonomerReadSeq], sort: bool = True
) -> Iterable[Tuple[str, List[MonomerReadSeq]]]:
    """Group aligns by concatemers."""
    seen = set()
    for concat_id, aligns in groupby(aligns, lambda x: x.concatemer_id):
        if concat_id in seen:
            raise ValueError(
                f"Concatemer '{concat_id}' has already been seen, "
                "these alignments should be sorted by read name or "
                f"{MOLECULE_TAG} tag"
            )
        aligns = list(aligns)
        if sort:
            aligns = sorted(aligns, key=lambda x: x.coords.start)
        yield (concat_id, aligns)
        seen.add(concat_id)


def pick_walk_aligns(
        aligns: List[MonomerReadSeq], concat_id: str) -> List[bool]:
    """Pick walk aligns."""
    expected_monomers = aligns[0].coords.subread_total
    keep = [False] * len(aligns)

    # pick one alignment per-monomer to use in the walk
    for monomer_id, cats in groupby(
        [
            (_.monomer_id, x, _.read_seq.flags.category) for x, _ in enumerate(aligns)],  # noqa
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
    """Annotate monomer alignments."""
    for concat_id, aligns in group_aligns_by_concatemers(
            mi_sorted_aligns, sort=True):
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
    """Sort aligns by concetemer idx."""
    _, orig_idx = zip(
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
        if align_counts["primary"] + align_counts["unmapped"] != expected_subreads: # noqa
            raise ValueError(
                f"Unexpected number of alignments: {align_counts} {expected_subreads}" # noqa
            )
    return new_aligns


def get_pairs(
    sorted_aligns: List[MonomerReadSeq],
    direct_only: bool = False,
) -> Iterable[PairedMonomers]:
    """Get pairs."""
    if len(sorted_aligns) == 0:
        yield (None, None, None)
    elif len(sorted_aligns) == 1:
        yield (sorted_aligns[0], None, None)
    else:
        if direct_only:
            pairs = zip(sorted_aligns[:-1], sorted_aligns[1:])
        else:
            pairs = combinations(sorted_aligns, 2)
        for x, (left, right) in enumerate(pairs):
            yield (left, right, PairData.from_monomer_pair(
                left, right, pair_idx=x))


def set_walk_tag(align: MonomerReadSeq, walk_str: str):
    """Set walk tag."""
    align.read_seq.tags[WALK_TAG] = f"{WALK_TAG}:Z:{walk_str}"


class PairAlignState(enum.IntEnum):
    """Pair align state."""

    both = 0
    left = 1
    right = 2
    neither = 4

    @staticmethod
    @lru_cache(maxsize=256)
    def from_flags(left_mapped: bool, right_mapped: bool) -> "PairAlignState":
        """From flags."""
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
    """Pair data."""

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
        """To flags."""
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

    @classmethod
    def from_monomer_pair(
        cls, left: MonomerReadSeq, right: MonomerReadSeq, pair_idx: int = -1
    ) -> "PairData":
        """From monomer pair."""
        concatemer_id = left.concatemer_id
        md_l, md_r = (
            left.coords,
            right.coords,
        )
        is_direct = (md_r.subread_idx - md_l.subread_idx) == 1
        read_distance = md_r.start - md_l.end
        align_state = PairAlignState.from_flags(
            left.read_seq.flags.unmap is False,
            right.read_seq.flags.unmap is False
        )
        if not (align_state == PairAlignState.both):
            return PairData(
                concatemer_id, pair_idx, is_direct, read_distance, align_state
            )
        else:
            l_align = left.read_seq.align_info
            r_align = right.read_seq.align_info
            is_cis = l_align.ref_name == r_align.ref_name  # type: ignore
            if is_cis:
                strands = (
                    l_align.strand,  # type: ignore
                    r_align.strand,  # type: ignore
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
    left: Tuple[Literal["+", "-", "."], int, int],
    right: Tuple[Literal["+", "-", "."], int, int],
):
    """Calculate genomics distance."""
    l_coord = left[2] if (left[0] == "+" or left[0] == ".") else left[1]
    r_coord = right[1] if (right[0] == "+" or right[0] == ".") else right[2]
    return r_coord - l_coord


def get_concatemer_align_string(sorted_aligns: List[MonomerReadSeq]) -> str:
    """Get concatemer align string."""
    _aligns = [a.read_seq.align_info for a in sorted_aligns]

    def _to_str(a: Optional[AlignInfo]) -> str:
        if a is None or a.flag == 4:
            return "*"
        else:
            return f"{a.ref_name}:{a.ref_pos}_{a.ref_pos + a.length}"

    return ",".join([_to_str(_) for _ in _aligns])
