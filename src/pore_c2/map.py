from pathlib import Path
from typing import Any, Dict, List

import mappy as mp
from attr import define
from pysam import FastxFile, FastxRecord

from .digest import ReadFragment, _get_enzyme, sequence_to_read_fragments


@define
class ReadFragmentAligns:
    pass


@define(kw_only=True)
class ConcatemerAlignData:
    read_id: str
    fragments: List[ReadFragment]
    aligns: List[List[mp.Alignment]]


def split_and_map_read(
    aligner: mp.Aligner, rec: FastxRecord, enzyme
) -> ConcatemerAlignData:
    frags, aligns = [], []
    for read_frag in sequence_to_read_fragments(enzyme, rec):
        seq, _ = read_frag.slice_fastq(rec)
        frags.append(read_frag)
        aligns.append(list(aligner.map(seq)))
    return ConcatemerAlignData(read_id=rec.name, fragments=frags, aligns=aligns)


def map_concatemers(
    *, enzyme: str, fastq: Path, mmi: Path, minimap_settings: Dict[Any, Any]
):
    aligner = mp.Aligner(fn_idx_in=str(mmi), **minimap_settings)
    enzyme = _get_enzyme(enzyme)
    for rec in FastxFile(str(fastq)):
        c_align_data = split_and_map_read(aligner, rec, enzyme)
        print(c_align_data)
        # print(
        #    read_frag.read_fragment_id,
        #    hit.ctg,
        #    hit.r_st,
        #    hit.r_en,
        #    hit.q_st,
        #    hit.q_en,
        # )
