from pathlib import Path

import pytest
from numpy.random import default_rng
from typer.testing import CliRunner

from pore_c2.testing import Scenario


@pytest.fixture
def runner():
    runner = CliRunner()
    return runner


@pytest.fixture(scope="session")
def scenario_dir(tmp_path_factory):
    base_dir = tmp_path_factory.mktemp("scenarios")
    return base_dir


@pytest.fixture(params=["EcoRI", "HindIII"], scope="session")
def scenario(request, scenario_dir):
    enzyme = request.param
    data_dir = scenario_dir / f"enzyme_{enzyme}"
    data_dir = data_dir.mkdir()
    s = Scenario(
        chrom_lengths={"chr1": 2000, "chr2": 1000},
        enzyme=enzyme,
        random_state=default_rng(421),
    )

    return s


@pytest.fixture(scope="session")
def default_scenario(scenario_dir):
    data_dir = scenario_dir / "default"
    data_dir.mkdir()
    s = Scenario(
        chrom_lengths={"chr1": 2000, "chr2": 1000},
        random_state=default_rng(421),
        temp_path=Path(str(data_dir)),
    )
    s.reference_fasta
    s.concatemer_fastq

    return s


@pytest.fixture(scope="session")
def het_scenario(scenario_dir):
    data_dir = scenario_dir / "default"
    data_dir.mkdir()
    s = Scenario(
        chrom_lengths={"chr1": 2000, "chr2": 1000},
        random_state=default_rng(421),
        temp_path=Path(str(data_dir)),
        num_haplotypes=2,
        variant_density=0.01,
    )
    s.reference_fasta
    s.concatemer_fastq
    s.phased_vcf
    return s


@pytest.fixture(scope="session")
def large_scenario(scenario_dir):
    data_dir = scenario_dir / "default"
    data_dir.mkdir()
    s = Scenario(
        chrom_lengths={"chr1": 20_000, "chr2": 1000},
        random_state=default_rng(421),
        temp_path=Path(str(data_dir)),
        num_concatemers=10_000,
    )
    s.reference_fasta
    s.concatemer_fastq
    return s
