from typing import List

from Bio import Restriction
from loguru import logger
from numpy.random import default_rng


def simulate_sequence_with_cut_sites(
    enzyme: str = "EcoRI", seq_length: int = 1000, cut_sites: List[int] = [50, 300, 750]
):
    seq = default_rng(42).choice(list("ACGT"), size=seq_length)
    enz = getattr(Restriction, enzyme, None)
    offset = enz.fst5
    length = len(enz.site)
    for p in cut_sites:
        start = p - offset
        logger.debug(enz.site)
        seq[start : start + length] = list(enz.site)
    return cut_sites, "".join(seq)
