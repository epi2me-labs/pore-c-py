from pathlib import Path
from typing import Any, Dict

import mappy as mp
from pysam import FastxFile

from .digest import _get_enzyme, sequence_to_read_fragments


def map_concatemers(
    *, enzyme: str, fastq: Path, mmi: Path, minimap_settings: Dict[Any, Any]
):
    aligner = mp.Aligner(fn_idx_in=mmi, **minimap_settings)
    enzyme = _get_enzyme(enzyme)
    for rec in FastxFile(str(fastq)):
        for read_frag in sequence_to_read_fragments(enzyme, rec):
            seq, _ = read_frag.slice_fastq(rec)
            for hit in aligner.map(seq):
                if not hit.is_primary:
                    continue
                # print(
                #    read_frag.read_fragment_id,
                #    hit.ctg,
                #    hit.r_st,
                #    hit.r_en,
                #    hit.q_st,
                #    hit.q_en,
                # )
