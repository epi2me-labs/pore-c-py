from typing import List

import pytest
from Bio import Restriction
from loguru import logger
from numpy.random import default_rng


def find_site_positions_biopython(enzyme: str, seq: str) -> List[int]:
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


@pytest.mark.parametrize("enzyme", ["EcoRI", "HindIII", "AloI"])
def test_enzyme_digest(enzyme):
    true_pos, seq = simulate_sequence_with_cut_sites(enzyme)
    logger.debug(seq)
    if enzyme != "AloI":
        positions = find_site_positions_biopython(enzyme, seq)
        assert true_pos == positions
    else:
        with pytest.raises(NotImplementedError):
            find_site_positions_biopython(enzyme, seq)
