"""Test cli."""
from collections import Counter
import json
from pathlib import Path
import shutil

import polars as pl
import pysam
import pytest
from typer.testing import CliRunner

from pore_c_py.main import (
    app,
    create_chunked_ubam,
    create_test_data,
    digest,
    export_bam,
    parse_bam,
)
from pore_c_py.sam_utils import pysam_verbosity
from pore_c_py.testing import Scenario


@pytest.mark.parametrize("command", ["utils process-monomer-alignments"])
def test_help(runner: CliRunner, command: str):
    """Test help."""
    result = runner.invoke(app, [command, "--help"])
    assert result.exit_code == 2


@pytest.mark.parametrize("suffix", [".bam", ".fastq"])
def test_digest_concatemers(default_scenario: Scenario, tmp_path, suffix):
    """Test digest concatemers."""
    output_file = tmp_path / f"read_fragments{suffix}"
    digest(
        default_scenario.concatemer_fastq,
        default_scenario.enzyme,
        output_file,
        remove_tags=None,
        max_reads=0,
    )
    expected = default_scenario.concatemer_metadata.select(
        ["concatemer_id", "num_monomers"]
    )
    observed = Counter()
    if suffix == ".fastq":
        for rec in pysam.FastxFile(output_file):
            tags = {_.split(":", 1)[0]: _ for _ in rec.comment.split("\t")}
            concatemer_id = tags["MI"].split(":")[-1]
            observed[concatemer_id] += 1
    elif suffix == ".bam":
        with pysam_verbosity(0):
            af = pysam.AlignmentFile(output_file, check_sq=False)
        for rec in af:
            concatemer_id = rec.get_tag("MI")
            observed[concatemer_id] += 1
    observed = pl.DataFrame(
        list(observed.items()), columns=["concatemer_id", "observed_mononmers"]
    )
    assert len(observed) == len(expected)
    oe_df = expected.join(observed, on="concatemer_id", how="outer")
    diff = oe_df.filter(pl.col("num_monomers") != pl.col("observed_mononmers"))
    assert len(diff) == 0


def test_create_test_data(tmp_path):
    """Test create test data."""
    scenario: Scenario = create_test_data(tmp_path, seed=421)
    existing_files = scenario.fc.existing()
    for k in ["reference_fasta", "concatemer_fastq"]:
        assert existing_files[k].exists() is True


@pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 is not in path")
def test_process_monomer_alignments(name_sorted_bam, tmp_path):
    """Test process monomer alignments."""
    prefix = tmp_path / "processed"
    monomer_bam = tmp_path / "processed.ns.bam"
    pe_bam = tmp_path / "processed.pe.bam"
    writer = parse_bam(
        name_sorted_bam,
        prefix,
        paired_end=True,
        chromunity=True,
        summary=True,
        chromunity_merge_distance=10,
        paired_end_minimum_distance=100,
        paired_end_maximum_distance=None,
    )
    writer.close()

    with pysam_verbosity(0):
        input_aligns = len(
            [a for a in pysam.AlignmentFile(str(name_sorted_bam))])

    assert writer.ns_writer is not None
    assert monomer_bam.exists()
    with pysam_verbosity(0):
        num_aligns = len([a for a in pysam.AlignmentFile(str(monomer_bam))])
    assert num_aligns == sum(writer.ns_writer.counter.values())
    assert num_aligns == input_aligns

    assert writer.pe_writer is not None
    with pysam_verbosity(0):
        num_pe_aligns = len([a for a in pysam.AlignmentFile(str(pe_bam))])
    assert num_pe_aligns == sum(writer.pe_writer.counter.values())
    assert writer.pe_writer.counter["primary"] > 0

    assert writer.pq_writer is not None
    _df = pl.read_parquet(str(writer.pq_writer.path))
    expected_fragments = sum(
        [
            writer.ns_writer.counter[k] for k in [
                "primary", "secondary", "supplementary"]]
    )
    assert len(_df) == writer.pq_writer.counter
    assert len(_df) <= expected_fragments
    assert _df["num_fragments"].sum() == expected_fragments

    assert writer.stats_writer is not None
    data = json.loads(writer.stats_writer.path.read_text())
    assert sum(data["cardinality"].values()) == _df["cid"].n_unique()


# TODO: this might be redundant
def test_digest_cli(runner, scenario: Scenario, tmp_path):
    """Test digest cli."""
    output_fastq = Path(tmp_path / "monomers.fastq")
    result = runner.invoke(
        app,
        [
            "digest",
            str(scenario.concatemer_fastq),
            str(scenario.enzyme),
            str(output_fastq),
        ],
    )
    assert result.exit_code == 0


@pytest.fixture
def processed_bam(name_sorted_bam, tmp_path):
    """Process bam."""
    prefix = tmp_path / "processed"
    monomer_bam = tmp_path / "processed.ns.bam"
    writer = parse_bam(
        name_sorted_bam,
        prefix,
        paired_end_minimum_distance=None,
        paired_end_maximum_distance=None,
        chromunity_merge_distance=None,
    )
    writer.close()
    return monomer_bam


@pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 is not in path")
def test_export_paired_end(processed_bam, tmp_path):
    """Test export paired end."""
    counts = {}
    defaults = dict(
        paired_end_maximum_distance=None,
        paired_end_minimum_distance=None,
    )

    for setting, params in [
        ("default", {}),
        (
            "no_trans",
            {"paired_end_minimum_distance": 0},
        ),  # setting the minimum distance to 0 gets rid of trans contacts
        ("no_short", {"paired_end_minimum_distance": 100}),
        (
            "no_short_or_long",
            {"paired_end_minimum_distance": 100,
             "paired_end_maximum_distance": 500},
        ),
        ("direct_only", {"direct_only": True}),
    ]:
        kwds = {**defaults, **params}
        writer = export_bam(
            processed_bam, tmp_path / setting, paired_end=True, **kwds)
        with pysam_verbosity(0):
            counts[setting] = len(
                [a for a in pysam.AlignmentFile(str(writer.pe_writer.path))]
            )

    # sanity check that increasing levels of
    # filters leads to fewer reads in the
    # results
    assert counts["no_trans"] < counts["default"]
    assert counts["no_short"] < counts["no_trans"]
    assert counts["no_short_or_long"] <= counts["no_short"]
    assert counts["direct_only"] < counts["default"]


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 is not in path") # noqa
def test_export_chromunity(processed_bam, tmp_path):
    """Test export chromunity."""
    counts = {}
    fragments = {}
    for setting, params in [
        ("default", {}),
        ("direct_only", {"direct_only": True}),
        ("merge_bookend", {"chromunity_merge_distance": 0}),
        ("merge_long", {"chromunity_merge_distance": 1000}),
    ]:
        writer = export_bam(
            processed_bam, tmp_path / setting, chromunity=True, **params
        )
        assert writer.pq_writer is not None
        _df = pl.read_parquet(str(writer.pq_writer.path))
        counts[setting] = len(_df)
        fragments[setting] = _df["num_fragments"].sum()

    # chromunity based on indirect
    assert counts["direct_only"] == counts["default"]
    assert counts["merge_bookend"] < counts["default"]
    assert counts["merge_long"] < counts["merge_bookend"]
    _frags = list(fragments.values())
    assert all([f == _frags[0] for f in _frags])


def test_create_chunked_bam(default_scenario: Scenario, tmp_path):
    """Test create chunked bam."""
    output_prefix = tmp_path / "chunked"
    res = create_chunked_ubam(
        default_scenario.concatemer_ubam,
        output_prefix,
        chunk_size=10,
        max_reads=0,
    )

    for path, expected in res:
        with pysam_verbosity(0):
            observed = len(
                [a for a in pysam.AlignmentFile(str(path), check_sq=False)])
        assert observed == expected
