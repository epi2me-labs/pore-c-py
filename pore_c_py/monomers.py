"""Monomers."""
from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Optional, TextIO

import polars as pl

from pore_c_py.model import Cutter


@dataclass
class SeqComposition:
    """Seq composition."""

    total_bases: int = -1
    gc_bases: int = 0
    at_bases: int = 0
    other_bases: int = 0
    masked_bases: int = 0

    @classmethod
    def from_seq(cls, seq: str) -> "SeqComposition":
        """From seq."""
        total = len(seq)
        char_count = Counter(list(seq))
        gc = char_count["G"] + char_count["g"] + \
            char_count["C"] + char_count["c"]
        at = char_count["A"] + char_count["a"] + \
            char_count["T"] + char_count["t"]
        other = total - (gc + at)
        masked = sum([v for k, v in char_count.items() if k.islower()])
        return SeqComposition(total, gc, at, other, masked)


@dataclass()
class GenomicFragment:
    """Genomic Fragment."""

    chrom: str
    start: int
    end: int
    fragment_id: str = field(init=False)
    fragment_idx: int = -1
    composition: SeqComposition = field(default_factory=SeqComposition)

    def __post_init__(self):
        """Post init."""
        self.fragment_id = f"{self.chrom}:{self.start}_{self.end}"

    def _flat_dict(self) -> Dict[str, Any]:
        """Flat dic."""
        d = asdict(self)
        composition = d.pop("composition")
        return {**d, **{f"composition.{k}": v for k, v in composition.items()}}

    @staticmethod
    def to_dataframe(data: "List[GenomicFragment]") -> pl.DataFrame:
        """To dataframe."""
        df = pl.DataFrame([_._flat_dict() for _ in data], orient="row")
        return df

    @staticmethod
    def from_cuts(
        chrom: str, chrom_length: int,
        positions: List[int], id_iter: Iterable[int]
    ) -> List["GenomicFragment"]:
        """From cuts."""
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
            GenomicFragment(
                chrom=chrom, start=start, end=end, fragment_idx=frag_id)
            for start, end, frag_id in zip(p[0:-1], p[1:], id_iter)
        ]


def sequence_to_genomic_fragments(
    cutter: Cutter, seq: str, chrom: str, fasta_fh: Optional[TextIO]
) -> List[GenomicFragment]:
    """Sequence to genomics fragments."""
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
