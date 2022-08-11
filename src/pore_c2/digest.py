from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, TextIO, Tuple

import polars as pl
from attrs import Factory, asdict, define, field, frozen
from Bio.Seq import Seq
from pysam import (  # pyright: reportGeneralTypeIssue=false
    FastaFile,
    FastxFile,
    FastxRecord,
)


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
        return [
            GenomicFragment(chrom=chrom, start=start, end=end, fragment_idx=frag_id)
            for start, end, frag_id in zip(
                [0] + positions, positions + [chrom_length], id_iter
            )
        ]


@define(kw_only=True)
class ReadFragment:
    read_fragment_id: str = field(init=False)
    read_start: int
    read_end: int
    read_id: str
    read_fragment_idx: int
    total_read_fragments: int

    def __attrs_post_init__(self):
        self.read_fragment_id = (
            f"{self.read_id}:{self.read_start}_{self.read_end}"
            f":{self.read_fragment_idx}_{self.total_read_fragments}"
        )

    def slice_fastq(self, read: FastxRecord) -> Tuple[str, str]:
        seq = read.sequence[self.read_start : self.read_end]
        qual = read.quality[self.read_start : self.read_end]
        return (seq, qual)


def _get_enzyme(enzyme_id: str):
    from Bio import Restriction

    enz = getattr(Restriction, enzyme_id, None)
    if enz is None:
        raise ValueError(f"Enzyme not found: {enzyme_id}")
    if enz.cut_twice():
        raise NotImplementedError(
            f"Enzyme cuts twice, not currently supported: {enzyme_id}"
        )
    return enz


def _get_cut_sites(enzyme, seq: str) -> List[int]:
    return [_ - 1 for _ in enzyme.search(Seq(seq))]


def _get_cut_intervals(enzyme, seq: str) -> Tuple[List[int], List[int]]:
    sites = _get_cut_sites(enzyme, seq)
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


def sequence_to_read_fragments(enzyme, read: FastxRecord) -> List[ReadFragment]:
    starts, ends = _get_cut_intervals(enzyme, read.sequence)
    num_fragments = len(starts)
    read_fragments = [
        ReadFragment(
            read_id=read.name,
            read_start=start,
            read_end=end,
            read_fragment_idx=x,
            total_read_fragments=num_fragments,
        )
        for x, (start, end) in enumerate(zip(starts, ends))
    ]
    return read_fragments


def sequence_to_genomic_fragments(
    enzyme, seq: str, chrom: str, fasta_fh: Optional[TextIO]
) -> List[GenomicFragment]:
    starts, ends = _get_cut_intervals(enzyme, seq)
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
    enzyme_id: str,
    fasta: Path,
    bed_file: Optional[Path] = None,
    fasta_out: Optional[Path] = None,
):
    enzyme = _get_enzyme(enzyme_id)

    if fasta_out:
        fasta_fh = fasta_out.open("w")
    else:
        fasta_fh = None

    ff = FastaFile(str(fasta))
    fragments = []
    for chrom in ff.references:
        seq = ff.fetch(chrom)  # , 0, 100_000)
        fragments.extend(
            sequence_to_genomic_fragments(enzyme, seq, chrom, fasta_fh=fasta_fh)
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


def digest_fastq(
    *, enzyme_id: str, fastq: Path, fastq_out: Path, return_dataframe: bool = False
):
    ff = FastxFile(str(fastq))
    outfh = fastq_out.open("w")

    enzyme = _get_enzyme(enzyme_id)
    data = []
    for read in ff:
        read_frags = sequence_to_read_fragments(enzyme, read)
        if return_dataframe:
            data.extend(read_frags)
        for frag in read_frags:
            seq, qual = frag.slice_fastq(read)
            outfh.write(f"@{frag.read_fragment_id}\n{seq}\n+\n{qual}\n")
    outfh.close()

    if return_dataframe:
        return pl.DataFrame([asdict(_) for _ in data], orient="row")
    else:
        return None
