"""IO."""
from collections import Counter, defaultdict
from dataclasses import dataclass, fields
from itertools import chain, islice
import json
from pathlib import Path
import sys
import time
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

import pyarrow as pa
import pyarrow.parquet as pq
from pysam import AlignmentFile, AlignmentHeader, FastaFile, FastxFile

from pore_c_py.aligns import get_pairs, group_colinear, PairedMonomers
from pore_c_py.log import get_logger
from pore_c_py.model import ConcatemerReadSeq, MonomerReadSeq, ReadSeq
from pore_c_py.sam_utils import downgrade_mm_tag, pysam_verbosity, SamFlags
from pore_c_py.settings import (
    DEFAULT_ALIGN_HEADER,
    DOWNGRADE_MM,
    EXE_NAME,
    VERSION
)

T = TypeVar("T", ReadSeq, ConcatemerReadSeq, MonomerReadSeq)


FASTQ_EXTS = [".fastq", ".fastq.gz", ".fq", ".fastq.gz"]
FASTA_EXTS = [".fasta", ".fa", ".fasta.gz", ".fasta.gz"]
SAM_EXTS = [".sam", ".bam", ".cram"]
logger = get_logger()


@dataclass
class FileCollection:
    """File collection."""

    # _path_attrs: List[str]

    @classmethod
    def with_prefix(cls, prefix: Path, drop: Optional[List[str]] = None):
        """With prefix."""
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

        return cls(**kwds)  # pyright: ignore [reportGeneralTypeIssues]

    @property
    def _path_attrs(self):
        """Path attrs."""
        return [f.name for f in fields(self) if f.type is Path]

    def __iter__(self):
        """Iterate."""
        for a in self._path_attrs:
            yield getattr(self, a)

    def items(self) -> List[Tuple[str, Optional[Path]]]:
        """Items."""
        return [(a, getattr(self, a)) for a in self._path_attrs]

    def existing(self) -> Dict[str, Path]:
        """EXisting."""
        return {
            key: val for key, val in self.items()
            if val is not None and val.exists()
        }

    def exists_any(self) -> bool:
        """Exist any."""
        return len(self.existing()) > 0

    def exists_all(self) -> bool:
        """Exists all."""
        for p in self:
            if p is not None and not p.exists():
                return False
        return True


@dataclass
class MappingFileCollection(FileCollection):
    """Mapping file collection."""

    bam: Path = Path("{prefix}.bam")
    temp_bam: Path = Path("{prefix}.temp.bam")
    concatemers: Path = Path("{prefix}.concatemers.parquet")
    contacts: Path = Path("{prefix}.contacts.parquet")
    unmapped: Path = Path("{prefix}.unmapped.fastq")


@dataclass
class AnnotatedMonomerFC(FileCollection):
    """Annotated monomer file collection."""

    namesorted_bam: Path = Path("{prefix}.ns.bam")
    paired_end_bam: Path = Path("{prefix}.pe.bam")
    summary_json: Path = Path("{prefix}.summary.json")
    chromunity_parquet: Path = Path("{prefix}.chromunity.parquet")


def get_concatemer_seqs(
    paths: List[Path], remove_tags: Optional[List[str]] = None
) -> Iterable[ConcatemerReadSeq]:
    """Get concatemer seq."""
    for p in paths:
        for read_seq in iter_reads(
            p, primary_only=True, as_unaligned=True, remove_tags=remove_tags
        ):
            yield ConcatemerReadSeq.from_readseq(read_seq)


def get_monomer_aligns(paths: List[Path]) -> Iterable[MonomerReadSeq]:
    """Get monomer aligns."""
    for p in paths:
        for read_seq in iter_reads(p, primary_only=True, as_unaligned=False):
            yield MonomerReadSeq.from_readseq(read_seq)


class Writer:
    """Write."""

    def __init__(self):
        """Init."""
        self.align_count = 0

    @classmethod
    def get_writer(
        cls,
        path: Path,
        as_unaligned: bool = False,
        header: Optional[AlignmentHeader] = None,
    ):
        """Get writer."""
        p = str(path)
        for e in SAM_EXTS:
            if p.endswith(e):
                assert header is not None
                return SamWriter(
                    path, as_unaligned=as_unaligned, header=header)
        for e in FASTQ_EXTS[:2]:
            if p.endswith(e):
                return FastqWriter(path)
        raise ValueError(f"Error guessing format from path: {path}")

    def write_record(self, rec: ReadSeq):
        """Write record."""
        ...

    def consume(self, read_seqs: Iterable[ReadSeq]):
        """Consume."""
        for seq in read_seqs:
            self.write_record(seq)

    def close(self):
        """Close."""
        ...


class FastqWriter(Writer):
    """Fastq writer."""

    def __init__(self, path: Path):
        """Init."""
        super().__init__()
        self.path = path
        self.fh = self.path.open("w")

    def write_record(self, rec: ReadSeq):
        """Write record."""
        self.fh.write(rec.to_fastq_str())

    def close(self):
        """Close."""
        self.fh.close()


class SamWriter(Writer):
    """Sam writer."""

    def __init__(
        self,
        path: Path,
        header: AlignmentHeader = DEFAULT_ALIGN_HEADER,
        as_unaligned: bool = False,
    ):
        """Init."""
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
        """Write record."""
        align = rec.to_align(
            header=self.header, as_unaligned=self.as_unaligned, **kwds)
        self.counter[SamFlags.from_int(align.flag).category.name] += 1
        self.writer.write(align)

    def close(self):
        """Close."""
        self.writer.close()


class PairedEndWriter(SamWriter):
    """Paired end writer."""

    def __init__(
        self,
        *args,
        min_distance: Optional[int] = None,
        max_distance: Optional[int] = None,
        keep_unpaired: bool = False,
        **kwds,
    ):
        """Init."""
        self.min_distance = min_distance if min_distance is not None else -1
        self.max_distance = max_distance if max_distance is not \
            None else float("inf")
        self.keep_unpaired = keep_unpaired
        self._filter_cis = (min_distance is not None) \
            | (max_distance is not None)

        super().__init__(*args, **kwds)

    def write_records(
        self, recs=List[PairedMonomers],
        fail_aligns=Optional[List[MonomerReadSeq]]
    ):
        """Write records."""
        for (left, right, pair_data) in recs:
            if left is None:  # no valid pairs
                continue
            if right is None:  # singleton
                if self.keep_unpaired:
                    self.write_record(left.read_seq)
                continue
            if self._filter_cis:
                if not pair_data.is_cis:
                    continue
                if not (
                    self.min_distance <= pair_data.genome_distance <= self.max_distance # noqa
                ):
                    continue

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
                l_sam_kwds["next_reference_start"] = (
                    right.read_seq.align_info.ref_pos + 1
                )
            if not r_flag.munmap:
                assert left.read_seq.align_info is not None
                r_sam_kwds["next_reference_name"] = (
                    "="
                    if pair_data.is_cis is True
                    else left.read_seq.align_info.ref_name
                )
                r_sam_kwds["next_reference_start"] = (
                    left.read_seq.align_info.ref_pos + 1
                )
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
        if fail_aligns:
            for rec in fail_aligns:
                self.write_record(rec.read_seq)


class ChromunityWriter:
    """Chromunity writer."""

    # TODO: see if chromunity can just parse this from the BAM in the future
    def __init__(self, path: Path, merge_distance: Optional[int] = None):
        """Init."""
        self.path = path
        self.schema = pa.schema(
            [
                ("cid", pa.string()),
                ("chrom", pa.string()),
                ("start", pa.uint32()),
                ("end", pa.uint32()),
                ("num_fragments", pa.uint32()),
            ]
        )
        self.writer = pq.ParquetWriter(str(self.path), self.schema)
        self.counter = 0
        self.merge_distance = merge_distance

    def write_records(self, recs: List[MonomerReadSeq]):
        """Write records."""
        if self.merge_distance is None:
            pylist = [
                {
                    "cid": r.concatemer_id,
                    "chrom": r.read_seq.align_info.ref_name,
                    "start": r.read_seq.align_info.ref_pos,
                    "end": r.read_seq.align_info.ref_end,
                    "num_fragments": 1,
                }
                for r in recs
                if r.read_seq.align_info is not None
            ]
        else:
            pylist = []
            align_blocks = group_colinear(
                [r.read_seq.align_info for r in recs], tol=self.merge_distance
            )
            for a in align_blocks:
                if len(a) == 1:
                    r = recs[a[0]]
                    if r.read_seq.align_info is not None:
                        pylist.append(
                            {
                                "cid": r.concatemer_id,
                                "chrom": r.read_seq.align_info.ref_name,
                                "start": r.read_seq.align_info.ref_pos,
                                "end": r.read_seq.align_info.ref_end,
                                "num_fragments": 1,
                            }
                        )
                else:
                    for x, idx in enumerate(a):
                        r = recs[idx]
                        # multi-alignment blocks should all be alignments
                        assert r.read_seq.align_info is not None
                        if x == 0:
                            cid = r.concatemer_id
                            chrom = r.read_seq.align_info.ref_name
                            start = r.read_seq.align_info.ref_pos
                            end = r.read_seq.align_info.ref_end
                        else:
                            start = min(start, r.read_seq.align_info.ref_pos)
                            end = max(end, r.read_seq.align_info.ref_end)
                    pylist.append(
                        {
                            "cid": cid,
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "num_fragments": len(a),
                        }
                    )
        if pylist:
            batch = pa.RecordBatch.from_pylist(pylist, schema=self.schema)
            self.writer.write_batch(batch)
            self.counter += len(pylist)

    def close(self):
        """Close."""
        self.writer.close()


class StatsWriter:
    """Stats writer."""

    def __init__(self, path: Path):
        """Init."""
        self.path = path
        self.concatemer_count = 0
        self.cardinality_count = defaultdict(int)
        self.pair_count = defaultdict(int)
        self.cis_trans = defaultdict(int)

    def write_records(self, pairs: List[PairedMonomers]):
        """Write records."""
        self.concatemer_count += 1
        for x, (left, right, pair_data) in enumerate(pairs):
            if x == 0:
                if left is None:
                    cardinality = 0
                elif right is None:
                    cardinality = 1
                else:
                    cardinality = left.coords.subread_total
                self.cardinality_count[cardinality] += 1
            if pair_data is not None:
                self.pair_count[pair_data.align_state.name] += 1
                if pair_data.align_state.name == "both":
                    if pair_data.is_cis:
                        self.cis_trans["cis"] += 1
                    else:
                        self.cis_trans["trans"] += 1

    def close(self):
        """Close."""
        d = {
            "cardinality": {k: v for k, v in self.cardinality_count.items()},
            "pair_count": {k: v for k, v in self.pair_count.items()},
            "cis_trans": {k: v for k, v in self.cis_trans.items()},
        }
        with self.path.open("w") as fh:
            fh.write(json.dumps(d))


@dataclass
class AnnotatedMonomerWriter(Writer):
    """Annotated Monomer writer."""

    ns_writer: Optional[SamWriter] = None
    pe_writer: Optional[PairedEndWriter] = None
    pq_writer: Optional[ChromunityWriter] = None
    stats_writer: Optional[StatsWriter] = None

    @classmethod
    def from_file_collection(
        cls,
        fc: AnnotatedMonomerFC,
        header: AlignmentHeader = DEFAULT_ALIGN_HEADER,
        chromunity_merge_distance: Optional[int] = None,
        paired_end_minimum_distance: Optional[int] = None,
        paired_end_maximum_distance: Optional[int] = None,
    ):
        """From file collection."""
        ns_writer = (
            SamWriter(
                fc.namesorted_bam, header=header) if fc.namesorted_bam else None # noqa
        )
        pe_writer = (
            PairedEndWriter(
                fc.paired_end_bam,
                header=header,
                min_distance=paired_end_minimum_distance,
                max_distance=paired_end_maximum_distance,
            )
            if fc.paired_end_bam
            else None
        )
        pq_writer = (
            ChromunityWriter(
                fc.chromunity_parquet, merge_distance=chromunity_merge_distance
            )
            if fc.chromunity_parquet
            else None
        )
        stats_writer = StatsWriter(fc.summary_json) if fc.summary_json else None # noqa
        return cls(
            ns_writer=ns_writer,
            pe_writer=pe_writer,
            pq_writer=pq_writer,
            stats_writer=stats_writer,
        )

    def consume(
        self,
        annotated_stream: Iterable[Tuple[str, List[MonomerReadSeq]]],
        direct_only: bool,
    ):
        """Consume."""
        for _, read_seqs in annotated_stream:
            if self.ns_writer:
                [self.ns_writer.write_record(_.read_seq) for _ in read_seqs]
            if self.pe_writer or self.stats_writer:
                pass_aligns = [
                    _ for _ in read_seqs if _.read_seq.flags.qcfail is False]
                fail_aligns = (
                    []
                    if len(pass_aligns) == len(read_seqs)
                    else [
                        _ for _ in read_seqs if _.read_seq.flags.qcfail is True] # noqa
                )
                pairs = list(get_pairs(pass_aligns, direct_only=direct_only))
                if self.pe_writer:
                    self.pe_writer.write_records(
                        pairs, fail_aligns=fail_aligns)
                if self.stats_writer:
                    self.stats_writer.write_records(pairs)
            if self.pq_writer:
                self.pq_writer.write_records(read_seqs)

    def close(self):
        """Close."""
        if self.ns_writer:
            self.ns_writer.close()
        if self.pe_writer:
            self.pe_writer.close()
        if self.pq_writer:
            self.pq_writer.close()
        if self.stats_writer:
            self.stats_writer.close()


def get_monomer_writer(
    path: Path, header: Optional[AlignmentHeader] = None,
    as_unaligned: bool = False
):
    """Get monomer writer."""
    writer = Writer.get_writer(path, header=header, as_unaligned=as_unaligned)
    return writer


def get_raw_reader(
    src: Path, reader_kwds: Optional[Mapping] = None
) -> Tuple[Literal["sam", "fastq", "fasta"], Union[AlignmentFile, FastxFile]]:
    """Get raw reader."""
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
    remove_tags: Optional[List[str]] = None,
) -> Iterable[ReadSeq]:
    """Iterate reads."""
    format, reader = get_raw_reader(src, reader_kwds=reader_kwds)
    if format == "sam":
        for rec in reader:
            if primary_only and (
                rec.is_supplementary or rec.is_secondary  # type: ignore
            ):
                continue
            if DOWNGRADE_MM:
                rec = downgrade_mm_tag(rec)  # type: ignore
            yield ReadSeq.from_align(
                rec, as_unaligned=as_unaligned, remove_tags=remove_tags
            )
    else:
        for rec in reader:
            yield ReadSeq.from_fastq(rec, remove_tags=remove_tags)


def find_files(
    root: Path, glob: str = "*.fastq", recursive: bool = True
) -> Iterable[Path]:
    """Find files."""
    if not root.is_dir():
        yield root
    else:
        if recursive and not glob.startswith("**/"):
            glob = f"**/{glob}"
        for f in root.glob(glob):
            yield (f)


def get_alignment_header(
    *,
    source_files: Optional[List[Path]] = None,
    reference_fasta: Optional[Path] = None,
    add_pg: bool = True,
) -> AlignmentHeader:
    """Get alignment header."""
    # TODO: add tests for this
    data = {}
    if source_files:
        bams = [f for f in source_files if f.suffix in SAM_EXTS]
        if len(bams) >= 1:
            with pysam_verbosity(0):
                source_header = AlignmentFile(
                    str(bams[0]), check_sq=False).header
            data = {**data, **source_header.to_dict()}
            if len(bams) > 1:
                logger.warning(
                    f"BAM header derived from first BAM only: {bams[0]}")
        else:
            pass
    if reference_fasta:
        ff = FastaFile(str(reference_fasta))
        header = AlignmentHeader.from_references(
            list(ff.references), list(ff.lengths))
        data = {**data, **header.to_dict()}

    if add_pg:
        _pg = data.pop("PG", [])
        pg_data = {
            "ID": f"{EXE_NAME}-{len(_pg) + 1}",
            "PN": EXE_NAME,
            "VN": VERSION,
            "CL": " ".join(sys.argv),
        }
        if len(_pg) > 0:
            if "ID" in _pg[-1]:
                pg_data["PP"] = _pg[-1]["ID"]
        _pg.append(pg_data)
        data["PG"] = _pg

    if not data:
        header = DEFAULT_ALIGN_HEADER
    else:
        header = AlignmentHeader.from_dict(data)

    return header


def fastq_to_ubam(source_fastqs: List[Path], output_ubam: Path) -> Path:
    """Fastq to ubam."""
    writer = SamWriter(output_ubam)
    src = [iter_reads(s) for s in source_fastqs]
    writer.consume(chain(*src))
    return output_ubam


# taken from recipe at https://docs.python.org/3/library/itertools.html
def batched(iterable, n):
    """Batch data into lists of length n. The last batch may be shorter."""
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := list(islice(it, n)):
        yield batch


def create_chunked_bam(
    input_files: List[Path],
    output_prefix: Path,
    chunk_size: int,
    max_reads: Optional[int] = None,
) -> List[Tuple[Path, int]]:
    """Create chunked bam."""
    header = get_alignment_header(source_files=input_files)
    output = []
    with pysam_verbosity(0):
        read_stream = chain(
            *(AlignmentFile(str(_), check_sq=False) for _ in input_files)
        )
        if max_reads:
            read_stream = islice(read_stream, max_reads)

        for batch_idx, reads in enumerate(batched(read_stream, chunk_size)):
            outfile = output_prefix.with_suffix(f".batch_{batch_idx}.bam")
            counter = 0
            start = time.perf_counter()
            with AlignmentFile(
                    str(outfile), mode="wb", header=header) as writer:
                for read in reads:
                    writer.write(read)
                    counter += 1
            elapsed = time.perf_counter() - start
            reads_per_minute = int(counter * 60 / elapsed)
            logger.info(
                f"Wrote batch {batch_idx}: {counter} reads "
                f"[{reads_per_minute:,d} reads/min]"
            )
            output.append((outfile, counter))
    return output
