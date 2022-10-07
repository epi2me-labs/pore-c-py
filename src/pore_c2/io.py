from pathlib import Path
from typing import Iterable, List, Literal, Mapping, Optional, Tuple, TypeVar, Union

from pysam import AlignmentFile, AlignmentHeader, FastaFile, FastxFile

from .log import get_logger
from .model import ConcatemerReadSeq, MonomerReadSeq, ReadSeq
from .settings import DEFAULT_ALIGN_HEADER

T = TypeVar("T", ReadSeq, ConcatemerReadSeq, MonomerReadSeq)


FASTQ_EXTS = [".fastq", ".fastq.gz", ".fq", ".fastq.gz"]
FASTA_EXTS = [".fasta", ".fa", ".fasta.gz", ".fasta.gz"]
SAM_EXTS = [".sam", ".bam", ".cram"]
logger = get_logger()


def get_concatemer_seqs(paths: List[Path]) -> Iterable[ConcatemerReadSeq]:
    for p in paths:
        for read_seq in iter_reads(p, primary_only=True, as_unaligned=True):
            yield ConcatemerReadSeq.from_readseq(read_seq)


def get_monomer_aligns(paths: List[Path]) -> Iterable[MonomerReadSeq]:
    for p in paths:
        for read_seq in iter_reads(p, primary_only=True, as_unaligned=True):
            yield MonomerReadSeq.from_readseq(read_seq)


class Writer:
    def __init__(self):
        self.align_count = 0

    @classmethod
    def get_writer(
        cls,
        path: Path,
        as_unaligned: bool = False,
        header: Optional[AlignmentHeader] = None,
    ):
        p = str(path)
        for e in SAM_EXTS:
            if p.endswith(e):
                return SamWriter(path, as_unaligned=as_unaligned, header=header)
        for e in FASTQ_EXTS[:2]:
            if p.endswith(e):
                return FastqWriter(path)
        raise ValueError(f"Error guessing format from path: {path}")

    def write_record(self, rec: ReadSeq):
        ...

    def consume(self, read_seqs: Iterable[ReadSeq]):
        for seq in read_seqs:
            self.write_record(seq)

    def close(self):
        ...


class FastqWriter(Writer):
    def __init__(self, path: Path):
        super().__init__()
        self.path = path
        self.fh = self.path.open("w")

    def write_record(self, rec: ReadSeq):
        self.fh.write(rec.to_fastq_str())

    def close(self):
        self.fh.close()


class SamWriter(Writer):
    def __init__(
        self,
        path: Path,
        header: AlignmentHeader = DEFAULT_ALIGN_HEADER,
        as_unaligned: bool = False,
    ):
        super().__init__()
        self.path = path
        self.header = header
        self.as_unaligned = as_unaligned

        if self.path.suffix == ".bam":
            mode = "wb"
        elif self.path.suffix == ".sam":
            mode = "w"
        elif self.path.suffix == ".cram":
            mode = "wc"
        else:
            raise NotImplementedError(self.path.suffix)
        self.writer = AlignmentFile(str(path), mode=mode, header=header)

    def write_record(self, rec: ReadSeq):
        self.writer.write(
            rec.to_align(header=self.header, as_unaligned=self.as_unaligned)
        )


def get_monomer_writer(
    path: Path, header: Optional[AlignmentHeader] = None, as_unaligned: bool = False
):
    writer = Writer.get_writer(path, header=header, as_unaligned=as_unaligned)
    return writer


def get_raw_reader(
    src: Path, reader_kwds: Optional[Mapping] = None
) -> Tuple[Literal["sam", "fastq", "fasta"], Union[AlignmentFile, FastxFile]]:

    if reader_kwds is None:
        reader_kwds = {}

    p = str(src)
    format = None
    for (fmt, exts) in [
        ("sam", SAM_EXTS),
        ("fastq", FASTQ_EXTS),
        ("fasta", FASTA_EXTS),
    ]:
        for e in exts:
            if p.endswith(e):
                format = fmt
                break
        if format is not None:
            break
    if not format:
        raise ValueError(f"Couldn't guess format of input: {src}")
    if format in ["fastq", "fasta"]:
        reader = FastxFile(p)
    else:
        reader = AlignmentFile(p, check_sq=False, **reader_kwds)

    return (format, reader)


def iter_reads(
    src: Path,
    reader_kwds: Optional[Mapping] = None,
    primary_only: bool = False,
    as_unaligned: bool = False,
) -> Iterable[ReadSeq]:
    format, reader = get_raw_reader(src, reader_kwds=reader_kwds)
    if format == "sam":
        for rec in reader:
            if primary_only and (rec.is_supplementary or rec.is_secondary):
                continue
            yield ReadSeq.from_align(rec, as_unaligned=as_unaligned)
    else:
        for rec in reader:
            yield ReadSeq.from_fastq(rec)


def find_files(
    root: Path, glob: str = "*.fastq", recursive: bool = True
) -> Iterable[Path]:

    if not root.is_dir():
        yield root
    else:
        if recursive and not glob.startswith("**/"):
            glob = f"**/{glob}"
        for f in root.glob(glob):
            yield (f)


def get_alignment_header(
    *, source_files: Optional[List[Path]] = None, reference_fasta: Optional[Path] = None
) -> AlignmentHeader:
    # TODO: add tests for this
    data = {}
    if source_files:
        bams = [f for f in source_files if f.suffix in SAM_EXTS]
        if len(bams) == 1:
            source_header = AlignmentFile(str(bams[0])).header
            data = {**data, **source_header.to_dict()}
        elif len(bams) > 1:
            raise NotImplementedError(f"Too many bams {bams}")
        else:
            pass
    if reference_fasta:
        ff = FastaFile(str(reference_fasta))
        header = AlignmentHeader.from_references(list(ff.references), list(ff.lengths))
        data = {**data, **header.to_dict()}

    if not data:
        header = DEFAULT_ALIGN_HEADER
    else:
        header = AlignmentHeader.from_dict(data)
    return header
