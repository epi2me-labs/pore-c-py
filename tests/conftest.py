import subprocess as sp
from pathlib import Path
from typing import Dict

import pytest
from numpy.random import default_rng
from typer.testing import CliRunner

from pore_c2.model import ReadSeq, TagData
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


@pytest.fixture(scope="session")
def name_sorted_bam(default_scenario: Scenario):
    ns_bam = default_scenario.temp_path / "name_sorted.bam"
    _ = sp.check_output(
        f"minimap2 -y -ax map-ont "
        f"{default_scenario.reference_fasta} {default_scenario.monomer_fastq} "
        f"| samtools sort -t MI > {ns_bam}",
        stderr=sp.STDOUT,
        shell=True,
    )
    return ns_bam


@pytest.fixture(scope="session")
def mock_read_no_qual() -> ReadSeq:
    return ReadSeq(name="read_no_qual", sequence="ACTG", quality=None)


@pytest.fixture(scope="session")
def mock_read_w_qual() -> ReadSeq:
    return ReadSeq(name="read_w_qual", sequence="ACTGACTG", quality="!" * 8)


@pytest.fixture(scope="session")
def mock_read_w_tags() -> ReadSeq:
    return ReadSeq(
        name="read_w_tags",
        sequence="AACGTTCGAAC",
        quality="!!00{}22[]]",
        tags={
            "RG": TagData.from_string("RG:Z:RG01"),
            "Mm": TagData.from_string("Mm:Z:C+m,0,1;"),
            "Ml": TagData.from_string("Ml:B:C,122,128"),
        },
    )


@pytest.fixture(scope="session")
def mock_reads(
    mock_read_no_qual: ReadSeq, mock_read_w_qual: ReadSeq, mock_read_w_tags
) -> Dict[str, ReadSeq]:
    reads = [mock_read_no_qual, mock_read_w_qual, mock_read_w_tags]
    # reads = [
    #    AlignData(name="read_no_qual", seq="ACTG"),
    #    AlignData(name="read_w_qual", seq="ACTGACTG", qual="!" * 8),
    #    AlignData(name="read_w_rg", seq="ACTGACTG", qual="!" * 8, tags=["RG:Z:RG01"]),
    #    AlignData(
    #        name="read_w_mods",
    #        seq="AACGTTCGAAC",
    #        qual="!!00{}22[]]",
    #        tags=["RG:Z:RG01", "Mm:Z:C+m,0,1;", "Ml:B:C,122,128"],
    #    ),
    # ]
    return {r.name: r for r in reads}
