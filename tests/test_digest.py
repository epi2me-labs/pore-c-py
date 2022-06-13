import pytest
from loguru import logger
from pandera.typing import DataFrame
from pysam import faidx

from pore_c2.cli import app
from pore_c2.digest import (
    DigestFragment,
    VirtualDigestSchema,
    find_cut_sites,
    virtual_digest,
)
from pore_c2.testing import (
    simulate_fasta_with_cut_sites,
    simulate_sequence_with_cut_sites,
)


@pytest.mark.parametrize("enzyme", ["EcoRI", "HindIII", "AloI"])
def test_enzyme_digest(enzyme):
    true_pos, seq = simulate_sequence_with_cut_sites(enzyme)
    if enzyme != "AloI":
        positions = find_cut_sites(enzyme, seq)
        assert true_pos == positions
    else:
        with pytest.raises(NotImplementedError):
            find_cut_sites(enzyme, seq)


@pytest.mark.parametrize("enzyme", ["EcoRI", "HindIII"])
def test_fastq_digest(enzyme, tmp_path):
    fasta = tmp_path / "genome.fa"
    lengths = {"chr1": 1000, "chr2": 500}
    expected = {"chr1": [10, 50], "chr2": [100]}
    simulate_fasta_with_cut_sites(fasta, lengths, expected, enzyme=enzyme)
    faidx(str(fasta))
    observed = virtual_digest(enzyme, fasta)
    expected = [
        DigestFragment(*_)
        for _ in [
            ("chr1", 0, 10, 1),
            ("chr1", 10, 50, 2),
            ("chr1", 50, 1000, 3),
            ("chr2", 0, 100, 4),
            ("chr2", 100, 500, 5),
        ]
    ]

    assert observed == expected


def test_serde():
    expected = [
        DigestFragment(*_)
        for _ in [
            ("chr1", 0, 10, 1),
            ("chr1", 10, 50, 2),
            ("chr1", 50, 1000, 3),
            ("chr2", 0, 100, 4),
            ("chr2", 100, 500, 5),
        ]
    ]
    df = DataFrame[VirtualDigestSchema](expected)
    logger.debug(df)


@pytest.mark.parametrize("enzyme", ["EcoRI", "HindIII"])
def test_digest_cli(runner, enzyme, tmp_path):
    fasta = tmp_path / "genome.fa"
    outfile = tmp_path / "results.pq"
    lengths = {"chr1": 1000, "chr2": 500}
    expected = {"chr1": [10, 50], "chr2": [100]}
    simulate_fasta_with_cut_sites(fasta, lengths, expected, enzyme=enzyme)
    faidx(str(fasta))
    result = runner.invoke(app, ["digest", str(fasta), enzyme, str(outfile)])
    print(result.stdout)
    assert result.exit_code == 0
