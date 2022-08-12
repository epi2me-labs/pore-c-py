from bisect import bisect_right
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import polars as pl
from attrs import define

from .log import get_logger

logger = get_logger()


@define(kw_only=True)
class OverlapStats:
    num_frags: int
    frag_ids: List[str]
    start_frag: str
    start_offset: int
    start_offset_perc: float
    end_frag: str
    end_offset: int
    end_offset_perc: float


@define(kw_only=True)
class FragmentOverlapper:
    left: Dict[str, List[int]]
    ids: Dict[str, List[str]]

    @classmethod
    def from_dataframe(cls, df: pl.DataFrame):
        left = {}
        ids = {}
        logger.info("Creating fragments from dataframe")
        for chrom, _df in df.partition_by(
            "chrom", maintain_order=True, as_dict=True
        ).items():
            logger.info(f"Fragments for {chrom}, {len(_df)}")
            left[chrom] = list(_df["start"]) + [_df["end"][-1] + 1]
            ids[chrom] = list(_df["fragment_id"]) + [_df["fragment_id"][-1]]
        return cls(left=left, ids=ids)

    @classmethod
    def from_parquet(cls, pq: Path):
        df = pl.read_parquet(pq, ["chrom", "start", "end", "fragment_id"])
        return cls.from_dataframe(df)

    def __attrs_post_init__(self):
        assert len(self.left) == len(self.ids)
        for k, v in self.left.items():
            assert len(v) == len(self.ids[k])
            increasing = np.all(np.diff(v) > 0)
            if not increasing:
                logger.debug(np.diff(v) > 0)
                raise ValueError

    def snap(self, chrom: str, coord: int) -> Tuple[int, int, float]:
        """Snap a coordinate to the closest junction"""
        # TODO: more comments
        idx = bisect_right(self.left[chrom], coord)
        try:
            f_start, f_end = self.left[chrom][idx - 1], self.left[chrom][idx]
        except IndexError:
            raise ValueError(coord, idx, len(self.left[chrom]), self.left[chrom][-1])
        d_left = coord - f_start
        d_right = f_end - coord
        frag_idx, delta = -1, -1
        if d_left <= d_right:
            frag_idx, delta = idx - 1, -d_left
        else:
            frag_idx, delta = idx, d_right
        perc = 100.0 * abs(delta) / (f_end - f_start)
        return frag_idx, delta, perc

    def overlaps(self, chrom, start: int, end: int):
        start_idx, start_offset, start_offset_perc = self.snap(chrom, start)
        end_idx, end_offset, end_offset_perc = self.snap(chrom, end)
        # span will be zero if both 'snap' to the same junction
        ids = self.ids[chrom][start_idx:end_idx]
        if isinstance(ids, str):
            ids = [ids]
        span = end_idx - start_idx
        return OverlapStats(
            num_frags=span,
            frag_ids=ids,  # pyright: ignore
            start_frag=self.ids[chrom][start_idx],
            start_offset=start_offset,
            start_offset_perc=start_offset_perc,
            end_frag=self.ids[chrom][end_idx],
            end_offset=end_offset,
            end_offset_perc=end_offset_perc,
        )
