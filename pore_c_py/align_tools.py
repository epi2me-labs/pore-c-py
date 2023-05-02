"""Things for operating on alignments."""
import copy
import enum
from itertools import combinations
import sys
from typing import List, Sequence

import pysam

import pore_c_py


class PairState(enum.Enum):
    """Pair State."""

    both = 0
    left = 1
    right = 2
    neither = 3
    singleton = 4


class PairedSegments:
    """Paired Segment class."""

    __align_state__ = {
        (False, False): PairState.both,
        (False, True): PairState.left,
        (True, False): PairState.right,
        (True, True): PairState.neither}

    def __init__(
            self,
            left: pysam.AlignedSegment,
            right: pysam.AlignedSegment):
        """Init."""
        # don't set these as self.left, self.right, just to save
        # some typing below
        left = copy.copy(left)
        right = copy.copy(right)

        if left is None:
            raise ValueError("Left read in pair cannot be None.")
        if right is None:  # singleton
            self.left = left
            self.right = None
            self.state = PairState.singleton
            return

        self.state = self.__align_state__[
            left.is_unmapped, right.is_unmapped]

        left.is_paired = True
        left.is_read1 = True
        left.mate_is_reverse = right.mate_is_reverse
        left.mate_is_unmapped = right.is_unmapped

        right.is_paired = True
        right.is_read2 = True
        right.mate_is_reverse = left.mate_is_reverse
        right.mate_is_unmapped = left.is_unmapped

        if self.state == PairState.both:
            if left.reference_name == right.reference_name:
                # TODO: check this: the spec says:
                # "each segment is properly aligned according to the aligner"
                # that doesn't strictly necessitate same reference for both
                left.is_proper_pair = True
                right.is_proper_pair = True
                left.next_reference_name = "="
                right.next_reference_name = "="
            else:
                left.next_reference_name = right.reference_name
                right.next_reference_name = left.reference_name
            left.next_reference_start = right.reference_start
            right.next_reference_start = left.reference_start
            mx = max(
                left.reference_end, right.reference_end)
            mn = min(
                left.reference_start, right.reference_start)
            template_length = mx - mn
            left.template_length = template_length
            right.template_length = -1 * template_length
        if self.state == PairState.left:
            left.mate_is_unmapped = True
            right.next_reference_name = left.reference_name
            right.next_reference_start = left.reference_start
        if self.state == PairState.right:
            right.mate_is_unmapped = True
            left.next_reference_name = right.reference_name
            left.next_reference_start = right.reference_start
        if self.state == PairState.neither:
            left.mate_is_unmapped = True
            right.mate_is_unmapped = True

        self.left = left
        self.right = right


def get_pairs(sorted_aligns, direct_only, keep_unpaired=True):
    """Get pairs."""
    if len(sorted_aligns) == 0:
        return
    elif len(sorted_aligns) == 1 and keep_unpaired:
        yield PairedSegments(sorted_aligns[0], None)
    else:
        if direct_only:
            pairs = zip(sorted_aligns[:-1], sorted_aligns[1:])
        else:
            pairs = combinations(sorted_aligns, 2)
        for left, right in pairs:
            yield PairedSegments(left, right)


def filter_pairs(
        pairs: Sequence[PairedSegments],
        min_distance=0, max_distance=float("inf"),
        allow_singletons=False, allow_improper=False, allow_unmapped=False):
    """Filter alignment pairs."""
    # TODO: something is weird here: filtering on distance, but then
    # allowing unmapped and improper is a bit weird. Should we stop
    # the caller from doing this?
    for pair in pairs:
        if pair.state == PairState.singleton:
            if allow_singletons:
                yield pair
            else:
                continue
        elif pair.state is not PairState.both:
            if allow_unmapped:
                yield pair
            else:
                continue
        # state == PairState.both
        distance = genomic_distance(pair.left, pair.right)
        if min_distance <= distance <= max_distance:
            yield pair


def genomic_distance(
        align1: pysam.AlignedSegment, align2: pysam.AlignedSegment):
    """Calculate distance between two alignments."""
    if (
            align1.is_unmapped or align2.is_unmapped
            or (align1.reference_name != align2.reference_name)):
        return float("inf")

    delta = min(align1.reference_end, align2.reference_end) \
        - max(align1.reference_start, align2.reference_start)
    return max(0, -delta)


def is_colinear(
        align1: pysam.AlignedSegment,
        align2: pysam.AlignedSegment, tol=0):
    """Check if two alignments are co-linear (and close) in the genome."""
    if align1.is_reverse != align2.is_reverse:
        return False
    distance = genomic_distance(align1, align2)
    return distance <= tol


def group_colinear(
        aligns: List[pysam.AlignedSegment], tol=0):
    """Group alignments into co-linear blocks."""
    if len(aligns) == 0:
        return []  # an empty list of blocks
    elif len(aligns) == 1:
        return [aligns]  # one block with single alignment
    res = []
    block = []
    last = aligns[0]
    block = [last]
    for aln in aligns[1:]:
        if is_colinear(last, aln, tol=tol):
            block.append(aln)
        else:
            res.append(block)
            block = [aln]
        last = aln
    res.append(block)
    return res


def update_header(header):
    """Add PG tag to existing BAM header.

    :param bamfile: a pysam.AlignmentFile.header
    """
    header = header.to_dict()
    pg = header.pop("PG", [])
    name = __package__.replace("_", "-")
    pg_data = {
        "ID": f"{name}-{len(pg) + 1}",
        "PN": name,
        "VN": pore_c_py.__version__,
        "CL": " ".join(sys.argv)}
    if len(pg) > 0:
        if "ID" in pg[-1]:
            pg_data["PP"] = pg[-1]["ID"]
    pg.append(pg_data)
    header["PG"] = pg
    return header
