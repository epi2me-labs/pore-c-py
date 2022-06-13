from pathlib import Path
from typing import Dict, List, Tuple

from Bio import Restriction
from numpy.random import default_rng


def simulate_sequence_with_cut_sites(
    enzyme: str = "EcoRI", seq_length: int = 1000, cut_sites: List[int] = [50, 300, 750]
) -> Tuple[List[int], str]:
    seq = default_rng(42).choice(list("ACGT"), size=seq_length)
    enz = getattr(Restriction, enzyme, None)
    offset = enz.fst5
    length = len(enz.site)
    for p in cut_sites:
        start = p - offset
        seq[start : start + length] = list(enz.site)
    return cut_sites, "".join(seq)


def simulate_fasta_with_cut_sites(
    fasta: Path,
    seq_length: Dict[str, int],
    cut_sites: Dict[str, List[int]],
    enzyme: str = "EcoRI",
) -> Path:
    fh = fasta.open("w")
    for chrom, length in seq_length.items():
        positions = cut_sites[chrom]
        _, seq = simulate_sequence_with_cut_sites(
            enzyme, cut_sites=positions, seq_length=length
        )
        fh.write(f">{chrom}\n{seq}\n")
    fh.close()
    return fasta
