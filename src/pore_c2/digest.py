from pathlib import Path
from typing import Dict, List

from loguru import logger
from pysam import FastaFile


def find_cut_sites(enzyme: str, seq: str) -> List[int]:
    from Bio import Restriction
    from Bio.Seq import Seq

    enz = getattr(Restriction, enzyme, None)
    if enz is None:
        raise ValueError(f"Enzyme not found: {enzyme}")
    if enz.cut_twice():
        raise NotImplementedError(
            f"Enzyme cuts twice, not currently supported: {enzyme}"
        )
    s = Seq(seq)
    positions = [_ - 1 for _ in enz.search(s)]
    return positions


def virtual_digest(enzyme: str, fasta: Path) -> Dict[str, int]:
    ff = FastaFile(fasta)
    chrom_lengths = [(t[0], t[1]) for t in zip(ff.references, ff.lengths)]
    logger.debug(chrom_lengths)
    res = {}
    for chrom, length in chrom_lengths:
        seq = ff.fetch(chrom)
        res[chrom] = find_cut_sites(enzyme, seq)
    return res
