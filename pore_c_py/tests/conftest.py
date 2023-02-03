"""Conf test."""
from pathlib import Path
import subprocess as sp
from typing import Dict

import pytest
from typer.testing import CliRunner

from pore_c_py.model import (
    AlignInfo,
    ConcatemerCoords,
    MonomerReadSeq,
    ReadSeq
)
from pore_c_py.sam_utils import MOLECULE_TAG
from pore_c_py.testing import Scenario


@pytest.fixture
def runner():
    """Runner."""
    runner = CliRunner()
    return runner


@pytest.fixture(scope="session")
def scenario_dir(tmp_path_factory):
    """Scenario dir."""
    base_dir = tmp_path_factory.mktemp("scenarios")
    return base_dir


@pytest.fixture(params=["EcoRI", "HindIII"], scope="session")
def scenario(request, scenario_dir):
    """Scenario."""
    enzyme = request.param
    data_dir = scenario_dir / f"enzyme_{enzyme}"
    data_dir = data_dir.mkdir()
    s = Scenario(
        enzyme=enzyme,
    )
    return s


@pytest.fixture(scope="session")
def default_scenario(scenario_dir):
    """Def scenario."""
    data_dir = scenario_dir / "default"
    data_dir.mkdir()
    s = Scenario(
        temp_path=Path(str(data_dir)),
    )
    s.reference_fasta
    s.concatemer_fastq
    return s


@pytest.fixture(scope="session")
def het_scenario(scenario_dir):
    """Het scenario."""
    data_dir = scenario_dir / "het"
    data_dir.mkdir()
    s = Scenario(
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
    """Large scenario."""
    data_dir = scenario_dir / "large"
    data_dir.mkdir()
    s = Scenario(
        temp_path=Path(str(data_dir)),
        num_concatemers=10_000,
    )
    s.reference_fasta
    s.concatemer_fastq
    return s


@pytest.fixture(scope="session")
def name_sorted_bam(default_scenario: Scenario):
    """Name sorted bam."""
    ns_bam = default_scenario.temp_path / "name_sorted.bam"
    _ = sp.check_output(
        f"minimap2 -y -ax map-ont "
        f"{default_scenario.reference_fasta} {default_scenario.monomer_fastq} "
        f"| samtools sort -t {MOLECULE_TAG} -o {ns_bam}",
        stderr=sp.STDOUT,
        shell=True,
    )
    return ns_bam


@pytest.fixture(scope="session")
def coord_sorted_bam(het_scenario: Scenario):
    """Coord sorted bam."""
    cs_bam = het_scenario.temp_path / "name_sorted.bam"
    _ = sp.check_output(
        f"minimap2 -y -ax map-ont "
        f"{het_scenario.reference_fasta} {het_scenario.monomer_fastq} "
        f"| samtools sort -o {cs_bam}",
        stderr=sp.STDOUT,
        shell=True,
    )
    return cs_bam


@pytest.fixture(scope="session")
def mock_read_no_qual() -> ReadSeq:
    """Mock read no qual."""
    return ReadSeq(name="read_no_qual", sequence="ACTG", quality=None)


@pytest.fixture(scope="session")
def mock_read_w_qual() -> ReadSeq:
    """Mock read with qual."""
    return ReadSeq(name="read_w_qual", sequence="ACTGACTG", quality="!" * 8)


@pytest.fixture(scope="session")
def mock_read_w_tags() -> ReadSeq:
    """Mock read with tags."""
    return ReadSeq(
        name="read_w_tags",
        sequence="AACGTTCGAAC",
        quality="!!00{}22[]]",
        tags={
            "RG": "RG:Z:RG01",
            "Mm": "Mm:Z:C+m,0,1;",
            "Ml": "Ml:B:C,122,128",
        },
    )


@pytest.fixture(scope="session")
def mock_read_w_rg() -> ReadSeq:
    """Mock read w_rg."""
    return ReadSeq(
        name="read_w_rg",
        sequence="AACGTTCGAAC",
        quality="!!00{}22[]]",
        tags={
            "RG": "RG:Z:RG01",
        },
    )


@pytest.fixture(scope="session")
def mock_reads(
    mock_read_no_qual: ReadSeq,
    mock_read_w_qual: ReadSeq,
    mock_read_w_tags,
    mock_read_w_rg,
) -> Dict[str, ReadSeq]:
    """Mock reads."""
    reads = [
        mock_read_no_qual, mock_read_w_qual, mock_read_w_tags, mock_read_w_rg]
    return {r.name: r for r in reads}


@pytest.fixture
def monomer_read_seqs():
    """Monomer read sequences."""
    num_concatemers, aligns_per_concatemer = 2, 10
    res = []
    for concat_idx in range(num_concatemers):
        concatemer_id = f"CONCAT{concat_idx}"
        for monomer_idx in range(aligns_per_concatemer):
            start = monomer_idx * 10
            end = start + 10
            monomer_id = f"{concatemer_id}:{start}:{end}"
            m = MonomerReadSeq(
                concatemer_id=concatemer_id,
                monomer_id=monomer_id,
                coords=ConcatemerCoords(
                    start=start,
                    end=end,
                    read_length=aligns_per_concatemer * 10,
                    subread_idx=monomer_idx,
                    subread_total=aligns_per_concatemer,
                ),
                read_seq=ReadSeq(
                    name=monomer_id,
                    sequence="A" * 10,
                    quality="!" * 10,
                    align_info=AlignInfo(
                        ref_name=f"chr{monomer_idx + 1}",
                        ref_pos=100 * monomer_idx,
                        cigar="10M",
                        map_quality=20,
                        length=10,
                    ),
                    tags={MOLECULE_TAG: f"{MOLECULE_TAG}:Z:{concatemer_id}"},
                ),
            )
            res.append(m)
    return res
