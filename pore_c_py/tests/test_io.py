"""Test IO."""
from dataclasses import dataclass
from pathlib import Path

from pysam import fqimport
import pytest

from pore_c_py.io import FileCollection, iter_reads


def test_file_collection(tmp_path):
    """Test file collection."""
    @dataclass
    class TestFC(FileCollection):
        a: Path = Path("{prefix}.A.txt")
        b: Path = Path("{prefix}.B.txt")

    test_fc = TestFC.with_prefix(Path(tmp_path))

    assert test_fc.exists_all() is False
    assert test_fc.exists_any() is False

    test_fc.a.write_text("a")

    assert test_fc.exists_all() is False
    assert test_fc.exists_any() is True

    test_fc.b.write_text("b")
    assert test_fc.exists_all() is True


@pytest.fixture
def mock_fastq(tmp_path, mock_reads):
    """Mock fastq."""
    outfile = tmp_path / "test.fastq"
    written = []
    with outfile.open("w") as fh:
        for r in mock_reads.values():
            if r.name != "read_no_qual":
                fh.write(r.to_fastq_str())
                written.append(r)
        fh.flush()
    return outfile


def test_fastq_reader(mock_fastq, mock_reads):
    """Test fastq reader."""
    written = [v for k, v in mock_reads.items() if k != "read_no_qual"]
    reads = list(iter_reads(mock_fastq))
    assert len(reads) == len(written)


def test_ubam_reader(mock_fastq, mock_reads):
    """Test ubam reader."""
    ubam = mock_fastq.with_suffix(".bam")
    written = [v for k, v in mock_reads.items() if k != "read_no_qual"]
    fqimport("-o", str(ubam), "-T", "*", str(mock_fastq))
    reads = list(iter_reads(ubam))
    assert len(reads) == len(written)


def test_remove_tags(mock_fastq, mock_reads):
    """Test remove tags."""
    initial_tags = [
        set(
            v.tags.keys()) for k, v in mock_reads.items() if k != "read_no_qual" # noqa
    ]
    assert initial_tags == [set(), {"RG", "Ml", "Mm"}, {"RG"}]

    fastq_tags = [
        set(
            v.tags.keys()) for v in iter_reads(
                mock_fastq, remove_tags=["Ml", "Mm"])
    ]
    assert fastq_tags == [set(), {"RG"}, {"RG"}]

    ubam = mock_fastq.with_suffix(".bam")
    fqimport("-o", str(ubam), "-T", "*", str(mock_fastq))
    bam_tags = [
        set(v.tags.keys()) for v in iter_reads(ubam, remove_tags=["Ml", "Mm"])]
    assert bam_tags == [set(), {"RG"}, {"RG"}]
