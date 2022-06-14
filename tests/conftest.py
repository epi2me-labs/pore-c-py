import pytest
from numpy.random import default_rng
from typer.testing import CliRunner

from pore_c2.testing import Scenario


@pytest.fixture
def runner():
    runner = CliRunner()
    return runner


@pytest.fixture(params=["EcoRI", "HindIII"])
def scenario(request, tmp_path):
    s = Scenario(
        chrom_lengths={"chr1": 2000, "chr2": 1000}, random_state=default_rng(421)
    )
    return s
