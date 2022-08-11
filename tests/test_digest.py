from pathlib import Path

import pytest

from pore_c2.cli import app
from pore_c2.digest import _get_cut_sites, _get_enzyme, digest_fastq, digest_genome
from pore_c2.testing import Scenario, simulate_sequence_with_cut_sites


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


def test_digest_genome(scenario):
    res = digest_genome(enzyme_id=scenario.enzyme, fasta=scenario.reference_fasta)
    assert res.shape == scenario.fragments_df.shape


def test_digest_fastq(scenario, tmp_path):
    output_fastq = tmp_path / "read_fragments.fastq"
    _ = digest_fastq(
        enzyme_id=scenario.enzyme,
        fastq=scenario.concatemer_fastq,
        fastq_out=output_fastq,
        return_dataframe=True,
    )
    # TODO: check that the fragments match the expectation


def test_digest_cli(runner, scenario: Scenario, tmp_path):
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
