from dataclasses import dataclass
from typing import Any, Dict, Tuple

import numpy as np
from pysam import AlignedSegment


@dataclass
class FragmentLoc:
    """Relative position within fragment"""

    fragment_id: int
    fp_dist: int
    tp_dist: int
    idx: int


@dataclass
class FragmentSpan:
    start: FragmentLoc
    end: FragmentLoc


@dataclass
class FragmentAssignment:
    fragment_id: int
    fragment_start: int
    fragment_end: int
    overlap_length: int = 0
    perc_of_alignment: float = 0.0
    perc_of_fragment: float = 0.0
    is_contained: bool = False


@dataclass
class FragmentsOverlap:
    num_fragments: int
    fragment_id_start: int
    fragment_id_end: int
    relative_start: Tuple[int, int]
    relative_end: Tuple[int, int]
    frag_start: Tuple[int, int]
    frag_end: Tuple[int, int]


@dataclass
class FragmentOverlapper:

    intervals: Dict[str, Any]
    fragment_ids: Dict[str, Any]

    @classmethod
    def from_fragments_df(cls, fragments_df):
        intervals = {}
        fragment_ids = {}
        for chrom, g in fragments_df.groupby("chrom"):
            intervals[chrom] = g["end"].values  # np.insert(g["end"].values, 0, [0])
            fragment_ids[chrom] = g["fragment_id"].values
            # np.insert(
            #    g["fragment_id"].values, 0, g["fragment_id"].values[0:1]
            # )
        return cls(intervals, fragment_ids)

    def find_overlap_pos(self, chrom: str, pos: int):
        """
        0 1 2 3 4
        | | | | |
        |A|C|G|T|
        ---------
        |0|1|2|3|
        -----      [0,2)
            ----   [2,4)

        https://cseducators.stackexchange.com/questions/5023/why-do-we-count-starting-from-zero/5026#5026
        """
        if not (0 <= pos <= self.intervals[chrom][-1]):
            raise ValueError(f"Position {pos} is out of bounds for chromosome {chrom}")
        # to the left of this interval
        idx = np.searchsorted(self.intervals[chrom], pos, "right")
        tp_dist = self.intervals[chrom][idx] - pos
        if idx == 0:
            fp_dist = pos
        else:
            fp_dist = pos - self.intervals[chrom][idx - 1]
        return FragmentLoc(self.fragment_ids[chrom][idx], fp_dist, tp_dist, idx)

    def find_overlap_span(self, chrom: str, start: int, end: int):

        start_loc = self.find_overlap_pos(chrom, start)
        end_loc = self.find_overlap_pos(chrom, end)
        return FragmentSpan(start_loc, end_loc)

    def assign_fragments(self, align: AlignedSegment):

        span = self.find_overlap_span(
            align.reference_name, align.reference_start, align.reference_end
        )
        sites = self.intervals[align.reference_name][span.start.idx : span.end.idx]

        num_frags = span.end.fragment_id - span.start.fragment_id
        return num_frags

        return sites

    def find_overlap(
        self, chrom: str, start: int, end: int, min_overlap_length: int = 5
    ):
        ivals = self.intervals[chrom]
        ids = self.fragment_ids[chrom]
        if end is None:
            raise ZeroDivisionError
        if start < self.intervals[chrom][0]:
            raise ValueError
        elif end > self.intervals[chrom][-1]:
            raise ValueError

        coords = np.array([start, end])
        # index of the interval the start,end is contained in
        indices = np.searchsorted(ivals, coords)
        frag_starts = ivals[indices - 1]
        if frag_starts[0] < 0:
            frag_starts[0] = 0
        frag_ends = ivals[indices]
        # distance to the fragment start site for each coord
        start_offset = coords - ivals[(indices - 1)]
        if indices[0] == 0:
            start_offset[0] = start
        # distance to the fragment end site for each coord
        end_offset = ivals[indices] - coords
        relative_start = (start_offset[0], end_offset[0])
        relative_end = (start_offset[1], end_offset[1])
        return FragmentsOverlap(
            indices[1] - indices[0] + 1,
            ids[indices[0]],
            ids[indices[1]],
            relative_start,
            relative_end,
            (frag_starts[0], frag_ends[0]),
            (frag_starts[1], frag_ends[1]),
        )
