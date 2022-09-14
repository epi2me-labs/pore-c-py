from pathlib import Path
from typing import Dict

import pytest
from pysam import AlignmentHeader, fqimport

from pore_c2.io import find_files
from pore_c2.io.reads import FastqReadIter, ReadIter
from pore_c2.monomers import AlignData, EnzymeCutter, digest_genome
from pore_c2.testing import Scenario, simulate_sequence_with_cut_sites


@pytest.fixture(scope="session")
def mock_reads() -> Dict[str, AlignData]:
    reads = [
        AlignData(name="read_no_qual", seq="ACTG"),
        AlignData(name="read_w_qual", seq="ACTGACTG", qual="!" * 8),
        AlignData(name="read_w_rg", seq="ACTGACTG", qual="!" * 8, tags=["RG:Z:RG01"]),
        AlignData(
            name="read_w_mods",
            seq="AACGTTCGAAC",
            qual="!!00{}22[]]",
            tags=["RG:Z:RG01", "Mm:Z:C+m,0,1;", "Ml:B:C,122,128"],
        ),
    ]
    return {r.name: r for r in reads}


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


def test_mods(mock_reads):
    test_read = mock_reads["read_w_mods"]
    header = AlignmentHeader.from_dict(
        {"RG": [{"ID": "RG1", "SM": "sample01", "LB": "lib01"}]}
    )
    seg = test_read.to_pysam(header)
    assert seg.modified_bases[("C", 0, "m")] == [(2, 122), (10, 128)]
    # assert seg.get_tag("RG") == "RG1"


def test_split_read(mock_reads):
    test_read = mock_reads["read_w_mods"]
    reads = test_read.split([5])
    # ps_reads = [r.to_pysam() for r in reads]
    for idx, s in [(0, slice(None, 5)), (1, slice(5, None))]:
        assert reads[idx].seq == test_read.seq[s]
        assert reads[idx].qual == test_read.qual[s]

    assert reads[0].to_pysam().modified_bases == {("C", 0, "m"): [(2, 122)]}
    assert reads[1].to_pysam().modified_bases == {("C", 0, "m"): [(5, 128)]}


@pytest.mark.parametrize(
    "layout,glob,recursive,expected",
    [
        (["A/a/a.fastq"], None, None, None),
        (["a.fastq"], None, None, None),
        (["a.fastq.gz"], None, None, []),
        (["a.fastq", "A/a/a.fastq"], None, False, ["a.fastq"]),
    ],
)
def test_find_files(tmp_path, layout, glob, recursive, expected):
    root = Path(str(tmp_path))
    if expected is None:
        expected = layout
    kwds = {}
    if glob is not None:
        kwds["glob"] = glob
    if recursive is not None:
        kwds["recursive"] = recursive
    for _ in layout:
        p = root / _
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(_)
    res = list(
        [str(_).replace(str(tmp_path), "")[1:] for _ in find_files(root, **kwds)]
    )
    assert res == expected


def test_fastq_reader(tmp_path, mock_reads):
    outfile = tmp_path / "test.fastq"
    written = []
    with outfile.open("w") as fh:
        for r in mock_reads.values():
            if r.name != "read_no_qual":
                fh.write(r.to_fastq())
                written.append(r)
        fh.flush()
    for x, r in enumerate(FastqReadIter(outfile)):
        assert r == written[x]


def test_ubam_reader(tmp_path, mock_reads):
    fastq = tmp_path / "test.fastq"
    ubam = tmp_path / "test.bam"
    written = []
    with fastq.open("w") as fh:
        for r in mock_reads.values():
            if r.name != "read_no_qual":
                fh.write(r.to_fastq())
                written.append(r)
        fh.flush()
    fqimport("-o", str(ubam), "-T", "*", str(fastq))

    for x, r in enumerate(ReadIter.load(ubam)):
        assert r == written[x]
