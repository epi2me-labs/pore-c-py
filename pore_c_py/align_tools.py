"""Things for operating on alignments."""
import sys
from typing import List

import pysam

import pore_c_py


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


def is_colinear(
        align1: pysam.AlignedSegment,
        align2: pysam.AlignedSegment, tol: int = 1):
    """Check if two alignments are co-linear in the genome."""
    # if either is is unmapped then they can't be co-linear
    if align1.is_unmapped or align2.is_unmapped:
        return False
    if align1.reference_name != align2.reference_name:
        return False
    if align1.is_reverse != align2.is_reverse:
        return False
    delta = min(align1.reference_end, align2.reference_end) \
        - max(align1.reference_start, align2.reference_start)
    # overlaps or within a distance
    return delta > 0 or delta < tol


def group_colinear(
        aligns: List[pysam.AlignedSegment], tol: int = 1
        ) -> List[List[pysam.AlignedSegment]]:
    """Group alignments into co-linear blocks."""
    if len(aligns) < 2:
        return [list(range(len(aligns)))]
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
