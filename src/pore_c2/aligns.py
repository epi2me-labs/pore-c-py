from collections import Counter
from itertools import combinations, groupby
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, List, Optional, Tuple

from attrs import define

from .log import get_logger
from .model import AlignData
from .settings import WALK_TAG
from .utils import FileCollection, SamFlags

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


def get_concatemer_id(align: AlignData) -> str:
    # TODO: use the MI tag for this
    return align.name.split(":", 1)[0]


def group_aligns_by_concatemers(
    aligns: Iterable[AlignData],
) -> Iterable[Tuple[str, List[AlignData]]]:
    seen = set()
    for concat_id, aligns in groupby(aligns, get_concatemer_id):
        if concat_id in seen:
            raise ValueError(
                f"Concatemer {concat_id} has already been seen, "
                "these alignments should be sorted by MI tag"
            )
        yield (concat_id, list(aligns))
        seen.add(concat_id)


def annotate_monomer_alignments(
    mi_sorted_aligns: Iterable[AlignData],
) -> Iterable[Tuple[str, List[AlignData]]]:
    for concat_id, aligns in group_aligns_by_concatemers(mi_sorted_aligns):
        # sort by read position so that we can correctly pair
        sorted_aligns = sort_aligns_by_concatemer_idx(aligns)
        align_str = get_concatemer_align_string(sorted_aligns)
        for a in aligns:
            set_walk_tag(a, align_str)
        for (left, right) in get_pairs(sorted_aligns, direct_only=True):
            set_next_read(left, right)
        yield (concat_id, sorted_aligns)


def sort_aligns_by_concatemer_idx(
    aligns: List[AlignData], validate: bool = True
) -> List[AlignData]:
    md = [a.concatemer_metadata for a in aligns]
    subread_idx, orig_idx = zip(*sorted([(a.subread_idx, x) for x, a in enumerate(md)]))
    new_aligns = [aligns[i] for i in orig_idx]
    if validate:
        expected_subreads = md[0].subread_total
        align_counts = Counter([SamFlags.from_int(a.flag).category for a in new_aligns])
        if align_counts["primary"] + align_counts["unmapped"] != expected_subreads:
            raise ValueError(
                f"Unexpected number of alignments: {align_counts} {expected_subreads}"
            )
    return new_aligns


def get_pairs(
    sorted_aligns: List[AlignData], direct_only: bool = False
) -> Iterable[Tuple[AlignData, AlignData]]:
    if direct_only:
        pairs = zip(sorted_aligns[:-1], sorted_aligns[1:])
    else:
        pairs = combinations(sorted_aligns, 2)
    for left, right in pairs:
        yield (left, right)


def set_walk_tag(align: AlignData, walk_str: str):
    align.tags.append(f"{WALK_TAG}:Z:{walk_str}")


def set_next_read(left: AlignData, right: AlignData):
    if right.ref_name != "*":
        if right.ref_name == left.ref_name:
            left.next_ref_name = "="
        else:
            left.next_ref_name = right.ref_name
        left.next_ref_pos = right.ref_pos


@define
class PairData:
    is_direct: bool
    read_distance: int
    both_aligned: bool
    is_cis: Optional[bool] = None
    genome_distance: int = -1


def get_pair_data(left: AlignData, right: AlignData) -> PairData:
    md_l, md_r = (
        left.concatemer_metadata,
        right.concatemer_metadata,
    )
    is_direct = (md_r.subread_idx - md_l.subread_idx) == 1

    aligned_l, aligned_r = left.flag != 4, right.flag != 4
    both_aligned = aligned_l & aligned_r
    read_distance = md_r.start - md_l.end
    if not both_aligned:
        return PairData(is_direct, read_distance, both_aligned)
    else:
        is_cis = left.ref_name == right.ref_name
        if is_cis:
            if left.ref_pos < right.ref_pos:
                genome_distance = right.ref_pos - (left.ref_pos + left.length)
            else:
                genome_distance = right.ref_pos - (right.ref_pos + right.length)
        else:
            genome_distance = -1
        return PairData(is_direct, read_distance, both_aligned, is_cis, genome_distance)


def get_concatemer_align_string(sorted_aligns: List[AlignData]) -> str:
    res = [
        f"{a.ref_name}:{a.ref_pos}_{a.ref_pos + a.length}" if a.flag != 4 else "*"
        for a in sorted_aligns
    ]
    return ",".join(res)


# TODO: delete this once we're sure we don't need it
# def map_concatemer_read(
#    *,
#    aligner: mp.Aligner,
#    read: AlignData,
#    cutter: Cutter,
#    thread_buf: Optional[mp.ThreadBuffer] = None,
# ) -> Tuple[AlignData, List[ReadFragment], List[mp.Alignment]]:
#
#    read_frags = sequence_to_read_fragments(cutter, read, store_sequence=False)
#    hits = [
#        list(aligner.map(read.seq[_.read_start : _.read_end], buf=thread_buf))
#        for _ in read_frags
#    ]
#    return (read, read_frags, hits)
