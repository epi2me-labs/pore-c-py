"""Tests to align_tools."""
import pysam
import pytest

from pore_c_py import align_tools

HEADER = pysam.AlignmentHeader.from_references(
    ["chr1", "chr2"], [1000, 2000])


def align_from_tuple(t):
    """Align from tuple."""
    rec = pysam.AlignedSegment(header=HEADER)
    rec.reference_id = ["chr1", "chr2"].index(t[0])
    rec.reference_start = t[1]
    rec.query_sequence = "A" * (t[2] - t[1])
    rec.cigartuples = [(pysam.CMATCH, (t[2] - t[1]))]
    rec.flag = 0
    if t[3] == "-":
        rec.flag |= pysam.FREVERSE
    if t[3] == ".":
        rec.flag |= pysam.FUNMAP
    return rec


@pytest.mark.parametrize(
    "left,right,tol,expected",
    [
        # overlapping
        (("chr1", 0, 10, "+"), ("chr1", 0, 20, "+"), 0, True),
        (("chr1", 0, 10, "+"), ("chr1", 0, 20, "-"), 0, False),
        (("chr1", 0, 10, "-"), ("chr1", 0, 20, "-"), 0, True),
        (("chr1", 0, 10, "+"), ("chr2", 0, 20, "+"), 0, False),
        # need tol=1 for bookended features
        (("chr1", 0, 10, "+"), ("chr1", 10, 20, "+"), 1, True),
        # allow some wiggle
        (("chr1", 0, 10, "+"), ("chr1", 5, 20, "+"), 5, True),
        # unaligned not colinear
        (("chr1", 0, 10, "+"), ("chr1", 0, 0, "."), 0, False)
    ],
)
def test_is_colinear(left, right, tol, expected):
    """Test is colinear."""
    res = align_tools.is_colinear(
        align_from_tuple(left),
        align_from_tuple(right),
        tol=tol,)
    assert res == expected


@pytest.mark.parametrize(
    "aligns,expected",
    [
        [
            [("chr1", 0, 10, "+"), ("chr1", 9, 20, "+"), ("chr2", 0, 10, "+")],
            [[0, 1], [2]],
        ],
        [
            [("chr1", 0, 8, "+"), ("chr1", 8, 9, "."), ("chr1", 20, 30, "+")],
            [[0], [1], [2]],
        ],
    ],
)
def test_group_colinear(aligns, expected):
    """Test group colinear."""
    alns = [align_from_tuple(t) for t in aligns]
    exp = list()
    for grp in expected:
        exp.append([alns[x] for x in grp])
    res = align_tools.group_colinear(alns)
    assert res == exp
