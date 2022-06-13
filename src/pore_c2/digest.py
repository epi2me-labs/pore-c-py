from dataclasses import dataclass, field
from itertools import count
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import pandera as pa
from pandera.typing import DataFrame, Series
from pysam import FastaFile


@dataclass
class DigestFragment:
    chrom: str
    start: int
    end: int
    fragment_id: int
    fragment_length: int = field(init=False)

    def __post_init__(self):
        self.fragment_length = self.end - self.start


class VirtualDigestSchema(pa.SchemaModel):
    chrom: Series[pd.CategoricalDtype] = pa.Field(dtype_kwargs={"ordered": True})
    start: Series[np.uint32] = pa.Field(ge=0, nullable=False)
    end: Series[np.uint32] = pa.Field(gt=0, nullable=False)
    fragment_id: Series[np.uint32] = pa.Field(unique=True, ge=1)
    fragment_length: Series[np.uint32] = pa.Field(ge=1)

    class Config:
        coerce = True


def save_digest(df: DataFrame[VirtualDigestSchema], outfile: Path):
    df.to_parquet(outfile, engine="pyarrow")


def load_digest(infile: Path) -> DataFrame[VirtualDigestSchema]:
    pd.read_parquet(infile, engine="pyarrow")


def find_cut_sites(enzyme: str, seq: str) -> List[int]:
    from Bio import Restriction
    from Bio.Seq import Seq

    enz = getattr(Restriction, enzyme, None)
    if enz is None:
        raise ValueError(f"Enzyme not found: {enzyme}")
    if enz.cut_twice():
        raise NotImplementedError(
            f"Enzyme cuts twice, not currently supported: {enzyme}"
        )
    s = Seq(seq)
    positions = [_ - 1 for _ in enz.search(s)]
    return positions


def virtual_digest(enzyme: str, fasta: Path) -> List[DigestFragment]:
    ff = FastaFile(fasta)
    chrom_lengths = [(t[0], t[1]) for t in zip(ff.references, ff.lengths)]
    res = []
    id_counter = count(1)
    for chrom, length in chrom_lengths:
        seq = ff.fetch(chrom)
        positions = find_cut_sites(enzyme, seq)
        res.extend(
            [
                DigestFragment(chrom, start, end, frag_id)
                for start, end, frag_id in zip(
                    [0] + positions, positions + [length], id_counter
                )
            ]
        )
    return res


def virtual_digest_df(enzyme: str, fasta: Path) -> DataFrame[VirtualDigestSchema]:
    df = DataFrame[VirtualDigestSchema](virtual_digest(enzyme, fasta))
    return df
