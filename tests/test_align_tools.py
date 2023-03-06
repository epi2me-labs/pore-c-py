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


# @pytest.mark.parametrize(
#     "left,right,expected",
#     [
#         (("+", 0, 5), ("+", 10, 15), 5),  # span 5-10
#         (("-", 10, 15), ("-", 0, 5), -5),  # span 5-10
#         (("-", 0, 5), ("+", 10, 15), 10),  # span 0-10
#         (("-", 10, 15), ("+", 0, 5), -10),  # span 0-10
#         (("+", 0, 5), ("-", 10, 15), 10),  # span 5-15
#         (("+", 10, 15), ("-", 0, 5), -10),  # span 5-15
#         (("-", 0, 5), ("-", 10, 15), 15),  # span 0-15
#         (("+", 10, 15), ("+", 0, 5), -15),  # span 0-15
#     ],
# )
# def test_genomic_distance(left, right, expected):
#     """Test genomic distance."""
#     assert calculate_genomic_distance(left, right) == expected
#
#
# def test_short_walk(m1, caplog):
#     """Test short walk."""
#     m1[2].read_seq.flags.secondary = True
#     with caplog.at_level(logging.WARNING):
#         m = list(
#             annotate_monomer_alignments(m1[:-1], remove_qcfail=False))[0][1]
#         qcfail = [_.read_seq.flags.qcfail for _ in m]
#         assert qcfail == [False, False, True]
#         assert (
#             caplog.record_tuples[-1][2]
#             == "Expected to see 3 alignments for concatemer C, found 2"
#         )
#
#
# @pytest.fixture
# def m2():
#     """Make mock monomer aligns m2."""
#     return (
#         _make_mock_monomer_aligns(
#             [
#                 (0, 5, "chr1", 0, None),
#                 (5, 10, "chr2", 100, None),
#                 (20, 25, "chr1", 50, None),
#                 (25, 30, "chr1", 500, None),
#             ]
#         ),
#         [
#             PairData(*_)
#             for _ in [
#                 ("C", 0, True, 0, PairAlignState.both, False, None, None),
#                 ("C", 1, True, 10, PairAlignState.both, False, None, None),
#                 ("C", 2, True, 0, PairAlignState.both, True, 445, True),
#             ]
#         ],
#     )
#
#
# def test_pairs_direct(m2):
#     """Test pairs direct."""
#     aligns, expected = m2
#     _, _, pair_data = zip(*get_pairs(aligns, direct_only=True))
#     assert list(pair_data) == expected
#
#
# @pytest.fixture
# def m3():
#     """Make mock monomer aligns m3."""
#     return (
#         _make_mock_monomer_aligns(
#             [
#                 (0, 5, "chr1", 0, None),
#                 (5, 10, "chr2", 100, None),
#                 (20, 25, "chr1", 50, None),
#             ]
#         ),
#         [
#             PairData(*_)
#             for _ in [
#                 ("C", 0, True, 0, PairAlignState.both, False, None, None),
#                 ("C", 1, False, 15, PairAlignState.both, True, 45, True),
#                 ("C", 2, True, 10, PairAlignState.both, False, None, None),
#             ]
#         ],
#     )
#
#
# def test_pairs_indirect(m3):
#     """Test pairs indirect."""
#     aligns, expected = m3
#     _, _, pair_data = zip(*get_pairs(aligns, direct_only=False))
#     assert list(pair_data) == expected
