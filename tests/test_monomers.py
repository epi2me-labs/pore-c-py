import pytest

from pore_c2.model import EnzymeCutter
from pore_c2.monomers import digest_genome
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


# @pytest.mark.monitor_test
# @pytest.mark.parametrize("n_procs", [0, 1, 2, 3, 4, 5])
# def test_digest_nproc(large_scenario, tmp_path, n_procs):
#    suffix = ".fastq"
#    output_file = tmp_path / f"read_fragments{suffix}"
#    digest_concatemers(
#        large_scenario.concatemer_fastq, large_scenario.enzyme, output_file,
#        n_proc=n_procs
#    )
