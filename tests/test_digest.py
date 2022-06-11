import pytest
from loguru import logger
from pysam import faidx

from pore_c2.digest import find_cut_sites, virtual_digest
from pore_c2.testing import simulate_sequence_with_cut_sites


@pytest.mark.parametrize("enzyme", ["EcoRI", "HindIII", "AloI"])
def test_enzyme_digest(enzyme):
    true_pos, seq = simulate_sequence_with_cut_sites(enzyme)
    logger.debug(seq)
    if enzyme != "AloI":
        positions = find_cut_sites(enzyme, seq)
        assert true_pos == positions
    else:
        with pytest.raises(NotImplementedError):
            find_cut_sites(enzyme, seq)


@pytest.mark.parametrize("enzyme", ["EcoRI", "HindIII"])
def test_fastq_digest(enzyme, tmp_path):
    fasta = tmp_path / "genome.fa"
    fh = fasta.open("w")
    expected = {}
    for chrom, length, positions in [("chr1", 1000, [10, 50]), ("chr2", 500, [100])]:
        expected[chrom] = positions
        _, seq = simulate_sequence_with_cut_sites(
            enzyme, cut_sites=positions, seq_length=length
        )
        fh.write(f">{chrom}\n{seq}\n")
    fh.close()
    faidx(str(fasta))
    logger.info(fasta)
    observed = virtual_digest(enzyme, fasta)
    assert observed == expected
