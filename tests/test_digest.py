import pytest
from pysam import faidx

from pore_c2.cli import app
from pore_c2.digest import _get_cut_sites, _get_enzyme, digest_fastq, digest_genome
from pore_c2.testing import (
    simulate_fasta_with_cut_sites,
    simulate_sequence_with_cut_sites,
)


@pytest.mark.parametrize("enzyme_id", ["EcoRI", "HindIII", "AloI"])
def test_enzyme_digest(enzyme_id):
    true_pos, seq = simulate_sequence_with_cut_sites(enzyme_id)
    if enzyme_id != "AloI":
        enzyme = _get_enzyme(enzyme_id)
        positions = _get_cut_sites(enzyme, seq)
        assert true_pos == positions
    else:
        with pytest.raises(NotImplementedError):
            enzyme = _get_enzyme(enzyme_id)
            _get_cut_sites(enzyme, seq)


def test_digest_genome(scenario, tmp_path):
    frag_table = tmp_path / "fragments.pq"
    res = digest_genome(scenario.enzyme, scenario.reference_fasta, frag_table)
    assert res.shape == scenario.fragments_df.shape


def test_digest_fastq(scenario, tmp_path):
    output_fastq = tmp_path / "read_fragments.fastq"
    res = digest_fastq(
        scenario.enzyme, scenario.concatemer_fastq, output_fastq, return_dataframe=True
    )
    print(res)
    # TODO: check that the fragments match the expectation


@pytest.mark.parametrize("enzyme", ["EcoRI", "HindIII"])
def test_digest_cli(runner, enzyme, tmp_path):
    fasta = tmp_path / "genome.fa"
    outfile = tmp_path / "results.pq"
    lengths = {"chr1": 1000, "chr2": 500}
    expected = {"chr1": [10, 50], "chr2": [100]}
    simulate_fasta_with_cut_sites(fasta, lengths, expected, enzyme=enzyme)
    faidx(str(fasta))
    result = runner.invoke(app, ["index", str(fasta), enzyme, str(outfile)])
    assert result.exit_code == 0
