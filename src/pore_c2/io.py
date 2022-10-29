from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import (
    Dict,
    Iterable,
    List,
    Literal,
    Mapping,
    Optional,
    Tuple,
    TypeVar,
    Union,
)

from attrs import define, fields
from pysam import AlignmentFile, AlignmentHeader, FastaFile, FastxFile

from .aligns import get_pairs
from .log import get_logger
from .model import ConcatemerReadSeq, MonomerReadSeq, ReadSeq
from .sam_tags import SamFlags, downgrade_mm_tag, pysam_verbosity
from .settings import DEFAULT_ALIGN_HEADER, DOWNGRADE_MM

T = TypeVar("T", ReadSeq, ConcatemerReadSeq, MonomerReadSeq)


FASTQ_EXTS = [".fastq", ".fastq.gz", ".fq", ".fastq.gz"]
FASTA_EXTS = [".fasta", ".fa", ".fasta.gz", ".fasta.gz"]
SAM_EXTS = [".sam", ".bam", ".cram"]
logger = get_logger()


@define
class FileCollection:
    _path_attrs: List[str]

    @classmethod
    def with_prefix(cls, prefix: Path, drop: Optional[List[str]] = None):
        path_attrs = []
        kwds = {}
        for f in fields(cls):  # pyright: ignore
            if f.name.startswith("_"):
                continue
            if drop and f.name in drop:
                kwds[f.name] = None
            else:
                kwds[f.name] = Path(str(f.default).format(prefix=str(prefix)))
            path_attrs.append(f.name)

        return cls(
            path_attrs=path_attrs, **kwds  # pyright: ignore [reportGeneralTypeIssues]
        )

    def __iter__(self):
        for a in self._path_attrs:
            yield getattr(self, a)

    def items(self) -> List[Tuple[str, Optional[Path]]]:
        return [(a, getattr(self, a)) for a in self._path_attrs]

    def existing(self) -> Dict[str, Path]:
        return {
            key: val for key, val in self.items() if val is not None and val.exists()
        }

    def exists_any(self) -> bool:
        return len(self.existing()) > 0

    def exists_all(self) -> bool:
        for p in self:
            if p is not None and not p.exists():
                return False
        return True


@define
class MappingFileCollection(FileCollection):
    bam: Path = Path("{prefix}.bam")
    temp_bam: Path = Path("{prefix}.temp.bam")
    concatemers: Path = Path("{prefix}.concatemers.parquet")
    contacts: Path = Path("{prefix}.contacts.parquet")
    unmapped: Path = Path("{prefix}.unmapped.fastq")


@define
class AnnotatedMonomerFC(FileCollection):
    namesorted_bam: Path = Path("{prefix}.ns.bam")
    paired_end_bam: Path = Path("{prefix}.pe.bam")
    concatemers: Path = Path("{prefix}.concatemers.parquet")


def get_concatemer_seqs(paths: List[Path]) -> Iterable[ConcatemerReadSeq]:
    for p in paths:
        for read_seq in iter_reads(p, primary_only=True, as_unaligned=True):
            yield ConcatemerReadSeq.from_readseq(read_seq)


def get_monomer_aligns(paths: List[Path]) -> Iterable[MonomerReadSeq]:
    for p in paths:
        for read_seq in iter_reads(p, primary_only=True, as_unaligned=False):
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
                assert header is not None
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
        self.counter = Counter()
        if self.path.suffix == ".bam":
            mode = "wb"
        elif self.path.suffix == ".sam":
            mode = "w"
        elif self.path.suffix == ".cram":
            mode = "wc"
        else:
            raise NotImplementedError(self.path.suffix)
        self.writer = AlignmentFile(str(path), mode=mode, header=header)

    def write_record(self, rec: ReadSeq, **kwds):
        align = rec.to_align(header=self.header, as_unaligned=self.as_unaligned, **kwds)
        self.counter[SamFlags.from_int(align.flag).category.name] += 1
        self.writer.write(align)

    def close(self):
        self.writer.close()


class PairedEndWriter(SamWriter):
    def write_records(self, recs: List[MonomerReadSeq], direct_only: bool = True):
        if len(recs) == 1:
            self.write_record(recs[0].read_seq)
        else:
            walk = [_ for _ in recs if _.read_seq.flags.qcfail is False]
            if len(walk) == 1:
                self.write_record(walk[0].read_seq)
            else:
                for (left, right, pair_data) in get_pairs(
                    walk, direct_only=direct_only
                ):
                    l_flag, r_flag = pair_data.to_flags(
                        align_flags=(left.read_seq.flags, right.read_seq.flags)
                    )
                    l_sam_kwds, r_sam_kwds = {}, {}
                    if not l_flag.munmap:
                        assert right.read_seq.align_info is not None
                        l_sam_kwds["next_reference_name"] = (
                            "="
                            if pair_data.is_cis is True
                            else right.read_seq.align_info.ref_name
                        )
                        l_sam_kwds[
                            "next_reference_start"
                        ] = right.read_seq.align_info.ref_pos
                    if not r_flag.munmap:
                        assert left.read_seq.align_info is not None
                        r_sam_kwds["next_reference_name"] = (
                            "="
                            if pair_data.is_cis is True
                            else left.read_seq.align_info.ref_name
                        )
                        r_sam_kwds[
                            "next_reference_start"
                        ] = left.read_seq.align_info.ref_pos
                    if pair_data.is_cis:
                        template_length = max(
                            left.read_seq.align_info.ref_end,  # type: ignore
                            right.read_seq.align_info.ref_end,  # type: ignore
                        ) - min(
                            left.read_seq.align_info.ref_pos,  # type: ignore
                            right.read_seq.align_info.ref_pos,  # type: ignore
                        )
                    else:
                        template_length = 0

                    self.write_record(
                        left.read_seq,
                        flag=l_flag,
                        template_length=template_length,
                        **l_sam_kwds,
                    )
                    self.write_record(
                        right.read_seq,
                        flag=r_flag,
                        template_length=template_length,
                        **r_sam_kwds,
                    )
            # write any qc fail records
            if len(walk) != len(recs):
                for rec in recs:
                    if rec.read_seq.flags.qcfail:
                        self.write_record(rec.read_seq)


@dataclass
class AnnotatedMonomerWriter(Writer):
    ns_writer: Optional[SamWriter] = None
    pe_writer: Optional[PairedEndWriter] = None

    @classmethod
    def from_file_collection(
        cls, fc: AnnotatedMonomerFC, header: AlignmentHeader = DEFAULT_ALIGN_HEADER
    ):
        ns_writer = (
            SamWriter(fc.namesorted_bam, header=header) if fc.namesorted_bam else None
        )
        pe_writer = (
            PairedEndWriter(fc.paired_end_bam, header=header)
            if fc.paired_end_bam
            else None
        )
        return cls(ns_writer=ns_writer, pe_writer=pe_writer)

    def consume(self, annotated_stream: Iterable[Tuple[str, List[MonomerReadSeq]]]):
        for _, read_seqs in annotated_stream:
            if self.ns_writer:
                [self.ns_writer.write_record(_.read_seq) for _ in read_seqs]
            if self.pe_writer:
                self.pe_writer.write_records(read_seqs)

    def close(self):
        if self.ns_writer:
            self.ns_writer.close()
        if self.pe_writer:
            self.pe_writer.close()


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
        with pysam_verbosity(0):
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

            if primary_only and (
                rec.is_supplementary or rec.is_secondary  # type: ignore
            ):
                continue
            if DOWNGRADE_MM:
                rec = downgrade_mm_tag(rec)  # type: ignore
            yield ReadSeq.from_align(rec, as_unaligned=as_unaligned)  # type: ignore
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
            with pysam_verbosity(0):
                source_header = AlignmentFile(str(bams[0]), check_sq=False).header
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
