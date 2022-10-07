from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, TextIO

import polars as pl
from attrs import Factory, asdict, define, field, frozen
from pysam import FastaFile

from .model import Cutter


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
