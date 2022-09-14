from itertools import groupby
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, List, Tuple

from attrs import define

from .log import get_logger
from .model import AlignData
from .utils import FileCollection

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


def sort_aligns_by_concatmer_idx(
    aligns: List[AlignData], validate: bool = True
) -> List[AlignData]:
    md = [a.concatemer_metadata for a in aligns]
    subread_idx, orig_idx = zip(*sorted([(a.subread_idx, x) for x, a in enumerate(md)]))
    new_aligns = [aligns[i] for i in orig_idx]
    if validate:
        expected_subreads = md[0].subread_total
        if len(subread_idx) != expected_subreads:
            raise ValueError(f"Unexpected number of alignments: {md}")
    return new_aligns


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
