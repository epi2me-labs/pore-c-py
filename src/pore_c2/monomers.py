from abc import ABCMeta, abstractmethod
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, TextIO, Tuple

import polars as pl
from attrs import Factory, asdict, define, field, frozen
from Bio.Seq import Seq
from pysam import AlignmentHeader, FastaFile

from .io import ReadIter, ReadWriter
from .model import AlignData


def digest_read(cutter: "Cutter", align: AlignData) -> List["AlignData"]:
    positions = cutter.get_cut_sites(align.seq)
    return align.split(positions)


def get_reads(paths: Iterable[Path]) -> Iterable[AlignData]:
    for p in paths:
        reader = ReadIter.load(p)
        for read in reader:
            yield read


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


def find_reads(root: Path, glob: str = "*.fastq", recursive: bool = True):

    if root.is_dir():
        files = find_files(root, glob=glob, recursive=recursive)
    else:
        files = (root,)

    for f in files:

        print(f)


def get_writer(
    path: Path, align_header: Optional[AlignmentHeader] = None
) -> ReadWriter:
    return ReadWriter.load(path, header=align_header)


class Cutter(metaclass=ABCMeta):
    @abstractmethod
    def get_cut_sites(self, seq: str) -> List[int]:
        ...

    def get_cut_intervals(self, seq: str) -> Tuple[List[int], List[int]]:
        sites = self.get_cut_sites(seq)
        seq_len = len(seq)
        if not sites:
            return [0], [seq_len]
        if sites[-1] == seq_len:
            sites = sites[:-1]
        starts = sites
        if starts[0] != 0:
            starts.insert(0, 0)
        ends = starts[1:] + [len(seq)]
        return starts, ends


class EnzymeCutter(Cutter):
    def __init__(self, enzyme: Any):
        self.enzyme = enzyme

    @classmethod
    def from_name(cls, enzyme_id: str):
        from Bio import Restriction

        enz = getattr(Restriction, enzyme_id, None)
        if enz is None:
            raise ValueError(f"Enzyme not found: {enzyme_id}")
        if enz.cut_twice():
            raise NotImplementedError(
                f"Enzyme cuts twice, not currently supported: {enzyme_id}"
            )
        return cls(enz)

    def get_cut_sites(self, seq: str) -> List[int]:
        return [_ - 1 for _ in self.enzyme.search(Seq(seq))]


@frozen
class SeqComposition:
    total_bases: int = -1
    gc_bases: int = 0
    at_bases: int = 0
    other_bases: int = 0
    masked_bases: int = 0

    @classmethod
    def from_seq(cls, seq: str) -> "SeqComposition":
        total = len(seq)
        char_count = Counter(list(seq))
        gc = char_count["G"] + char_count["g"] + char_count["C"] + char_count["c"]
        at = char_count["A"] + char_count["a"] + char_count["T"] + char_count["t"]
        other = total - (gc + at)
        masked = sum([v for k, v in char_count.items() if k.islower()])
        return SeqComposition(total, gc, at, other, masked)


@define(kw_only=True)
class GenomicFragment:
    chrom: str
    start: int
    end: int
    fragment_id: str = field(init=False)
    fragment_idx: int = -1
    composition: SeqComposition = Factory(SeqComposition)

    def __attrs_post_init__(self):
        self.fragment_id = f"{self.chrom}:{self.start}_{self.end}"

    def _flat_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        composition = d.pop("composition")
        return {**d, **{f"composition.{k}": v for k, v in composition.items()}}

    @staticmethod
    def to_dataframe(data: "List[GenomicFragment]") -> pl.DataFrame:
        df = pl.DataFrame([_._flat_dict() for _ in data], orient="row")
        return df

    @staticmethod
    def from_cuts(
        chrom: str, chrom_length: int, positions: List[int], id_iter: Iterable[int]
    ) -> List["GenomicFragment"]:
        if len(positions) == 0:
            p = [0, chrom_length]
        else:
            start_zero = positions[0] == 0
            end_length = positions[-1] == chrom_length
            if start_zero and end_length:
                p = positions
            elif start_zero:
                p = positions + [chrom_length]
            elif end_length:
                p = [0] + positions
            else:
                p = [0] + positions + [chrom_length]

        return [
            GenomicFragment(chrom=chrom, start=start, end=end, fragment_idx=frag_id)
            for start, end, frag_id in zip(p[0:-1], p[1:], id_iter)
        ]


@define(kw_only=True)
class ReadFragment:
    read_fragment_id: str = field(init=False)
    read_start: int
    read_end: int
    read_id: str
    read_fragment_idx: int
    total_read_fragments: int
    sequence: Optional[str] = None

    def __attrs_post_init__(self):
        self.read_fragment_id = (
            f"{self.read_id}:{self.read_start}_{self.read_end}"
            f":{self.read_fragment_idx}_{self.total_read_fragments}"
        )

    def slice_fastq(self, read: AlignData) -> Tuple[str, str]:
        seq = read.seq[self.read_start : self.read_end]
        qual = read.qual[self.read_start : self.read_end]
        return (seq, qual)


def sequence_to_read_fragments(
    cutter: Cutter, read: AlignData, store_sequence: bool = True
) -> List[ReadFragment]:
    starts, ends = cutter.get_cut_intervals(read.seq)
    num_fragments = len(starts)
    read_fragments = [
        ReadFragment(
            read_id=read.name,
            read_start=start,
            read_end=end,
            read_fragment_idx=x,
            total_read_fragments=num_fragments,
            sequence=read.seq[start:end] if store_sequence else None,
        )
        for x, (start, end) in enumerate(zip(starts, ends))
    ]
    return read_fragments


def sequence_to_genomic_fragments(
    cutter: Cutter, seq: str, chrom: str, fasta_fh: Optional[TextIO]
) -> List[GenomicFragment]:
    starts, ends = cutter.get_cut_intervals(seq)
    res = []
    for _, (start, end) in enumerate(zip(starts, ends)):
        frag_seq = seq[start:end]
        frag = GenomicFragment(
            chrom=chrom,
            start=start,
            end=end,
            composition=SeqComposition.from_seq(frag_seq),
        )
        if fasta_fh:
            fasta_fh.write(f">{frag.fragment_id}\n{frag_seq}\n")
        res.append(frag)
    return res


def digest_genome(
    *,
    cutter: Cutter,
    fasta: Path,
    bed_file: Optional[Path] = None,
    fasta_out: Optional[Path] = None,
) -> pl.DataFrame:

    if fasta_out:
        fasta_fh = fasta_out.open("w")
    else:
        fasta_fh = None

    ff = FastaFile(str(fasta))
    fragments = []
    for chrom in ff.references:
        seq = ff.fetch(chrom)  # , 0, 100_000)
        fragments.extend(
            sequence_to_genomic_fragments(cutter, seq, chrom, fasta_fh=fasta_fh)
        )

    df = (
        GenomicFragment.to_dataframe(fragments)
        .drop("fragment_idx")
        .with_row_count(name="fragment_idx", offset=1)
    )
    if bed_file:
        df.select(["chrom", "start", "end", "fragment_id"]).write_csv(
            bed_file, has_header=False, sep="\t"
        )
    return df
