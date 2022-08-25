from pathlib import Path

import pytest

from pore_c2.cli import app
from pore_c2.digest import EnzymeCutter, digest_fastq, digest_genome
from pore_c2.reads import get_reads
from pore_c2.testing import Scenario, simulate_sequence_with_cut_sites


@pytest.mark.parametrize("enzyme_id", ["EcoRI", "HindIII", "AloI"])
def test_enzyme_digest(enzyme_id):
    true_pos, seq = simulate_sequence_with_cut_sites(enzyme_id)
    if enzyme_id != "AloI":
        cutter = EnzymeCutter.from_name(enzyme_id)
        positions = cutter.get_cut_sites(seq)
        assert true_pos == positions
    else:
        with pytest.raises(NotImplementedError):
            cutter = EnzymeCutter.from_name(enzyme_id)
            cutter.get_cut_sites(seq)


def test_digest_genome(scenario: Scenario):
    cutter = EnzymeCutter.from_name(scenario.enzyme)
    res = digest_genome(cutter=cutter, fasta=scenario.reference_fasta)
    assert res.shape == scenario.fragments_df.shape


def test_digest_fastq(scenario: Scenario, tmp_path):
    cutter = EnzymeCutter.from_name(scenario.enzyme)
    output_fastq = tmp_path / "read_fragments.fastq"
    _ = digest_fastq(
        cutter=cutter,
        read_iter=get_reads(scenario.concatemer_fastq),
        fastq_out=output_fastq,
        return_dataframe=True,
    )
    # TODO: check that the fragments match the expectation


def test_index_cli(runner, scenario: Scenario, tmp_path):
    prefix = Path(tmp_path / "test_digest")
    result = runner.invoke(
        app,
        [
            "index",
            str(scenario.reference_fasta),
            str(scenario.enzyme),
            "--prefix",
            str(prefix),
        ],
    )
    assert result.exit_code == 0


def test_digest_cli(runner, scenario: Scenario, tmp_path):
    output_fastq = Path(tmp_path / "monomers.fastq")
    result = runner.invoke(
        app,
        [
            "utils",
            "digest-concatemers",
            str(scenario.concatemer_fastq),
            str(scenario.enzyme),
            str(output_fastq),
        ],
    )
    assert result.exit_code == 0
