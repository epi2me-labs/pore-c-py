from pathlib import Path
from typing import Any, Dict, List

import mappy as mp
from attr import define
from pysam import FastxFile, FastxRecord

from .digest import ReadFragment, _get_enzyme, sequence_to_read_fragments
from .overlaps import FragmentOverlapper, OverlapStats


@define(kw_only=True)
class ReadFragmentAligns:
    read_fragment: ReadFragment
    aligns: List[mp.Alignment]
    overlaps: List[OverlapStats]


@define(kw_only=True)
class ConcatemerAlignData:
    read_id: str
    aligns: List[ReadFragmentAligns]


def split_and_map_read(
    *, aligner: mp.Aligner, rec: FastxRecord, enzyme, overlapper: FragmentOverlapper
) -> ConcatemerAlignData:
    aligns = []
    for read_frag in sequence_to_read_fragments(enzyme, rec):
        seq, _ = read_frag.slice_fastq(rec)
        hits = list(aligner.map(seq))
        overlaps = [overlapper.overlaps(h.ctg, h.r_st, h.r_en) for h in hits]
        aligns.append(
            ReadFragmentAligns(read_fragment=read_frag, aligns=hits, overlaps=overlaps)
        )
    return ConcatemerAlignData(read_id=rec.name, aligns=aligns)


def map_concatemers(
    *,
    enzyme: str,
    fastq: Path,
    mmi: Path,
    minimap_settings: Dict[Any, Any],
    fragment_pq: Path
):
    overlapper = FragmentOverlapper.from_parquet(fragment_pq)
    aligner = mp.Aligner(fn_idx_in=str(mmi), **minimap_settings)
    enzyme = _get_enzyme(enzyme)
    for rec in FastxFile(str(fastq)):
        _ = split_and_map_read(
            aligner=aligner, rec=rec, enzyme=enzyme, overlapper=overlapper
        )
        # print(c_align_data)
        # print(
        #    read_frag.read_fragment_id,
        #    hit.ctg,
        #    hit.r_st,
        #    hit.r_en,
        #    hit.q_st,
        #    hit.q_en,
        # )
