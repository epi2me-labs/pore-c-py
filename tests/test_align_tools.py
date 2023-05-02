"""Tests to align_tools."""
import pysam
import pytest

from pore_c_py import align_tools


HEADER = pysam.AlignmentHeader.from_references(
    ["chr1", "chr2", "*"], [1000, 2000, 3000])

def align_from_tuple(t):
    """Align from tuple."""
    rec = pysam.AlignedSegment(header=HEADER)
    rec.reference_id = ["chr1", "chr2", "*"].index(t[0])
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
        [
            [("chr1", 0, 8, "+"), ("chr1", 8, 9, ".")],
            [[0], [1]],
        ],
        [
            [("chr1", 0, 8, "+")],
            [[0]],
        ],
        [
            [],
            []
        ]
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


@pytest.mark.parametrize(
    "left,right,expected",
    [   # overlapping
        (("chr1", 0, 5, "+"), ("chr1", 2, 15, "+"), 0),
        # not overlapping
        (("chr1", 0, 5, "+"), ("chr1", 10, 15, "+"), 5),
        (("chr1", 10, 15, "-"), ("chr1", 0, 5, "+"), 5),
        (("chr1", 0, 5, "-"), ("chr1", 10, 15, "+"), 5),
        (("chr1", 10, 15, "-"), ("chr1", 0, 5, "+"), 5),
        (("chr1", 0, 10, "+"), ("chr1", 11, 15, "-"), 1),
        (("chr1", 0, 5, "+"), ("chr1", 10, 15, "-"), 5),
        (("chr1", 10, 15, "+"), ("chr1", 0, 5, "-"), 5),
        (("chr1", 0, 5, "-"), ("chr1", 10, 15, "-"), 5),
        (("chr1", 10, 15, "+"), ("chr1", 0, 5, "+"), 5),
        # unmapped
        (("chr1", 0, 10, "+"), ("chr1", 11, 15, "."), float("inf")),
        (("chr1", 0, 10, "."), ("chr1", 11, 15, "."), float("inf")),

    ],
)
def test_genomic_distance(left, right, expected):
    """Test genomic distance."""

    res = align_tools.genomic_distance(
        align_from_tuple(left), align_from_tuple(right))
    assert res == expected


@pytest.mark.parametrize(
    "left, right, expected_left_flags, expected_right_flags",
    [
        # Neither
        (("chr1", 0, 5, "."), ("chr1", 0, 5, "."),
         {"mate_is_unmapped": True}, {"mate_is_unmapped": True}),
        # Both same ref
        (("chr1", 0, 5, "+"), ("chr1", 10, 20, "-"),
         {"is_proper_pair": True, "next_reference_name": "chr1",
         "next_reference_start": 10, "template_length": 20},
         {"is_proper_pair": True, "next_reference_name": "chr1",
         'next_reference_start': 0, 'template_length': -20}),
        # Both different ref
        (("chr1", 0, 5, "+"), ("chr2", 10, 20, "+"),
         {"is_proper_pair": False, "next_reference_name": "chr2",
          "next_reference_start": 10},
         {"is_proper_pair": False, "next_reference_name": "chr1",
          "next_reference_start": 0}),
        # right
        (("chr1", 0, 5, "."), ("chr2", 10, 20, "+"),
         {"mate_is_unmapped": False, "next_reference_name": "chr2",
          "next_reference_start": 10},
         {"mate_is_unmapped": True, "next_reference_name": None,
          "next_reference_start": -1}),
        # left
        (("chr1", 0, 5, "+"),("chr2", 10, 20, "."),
         {"mate_is_unmapped": True, "next_reference_name": None,
         "next_reference_start": -1},
         {"mate_is_unmapped": False, "next_reference_name": "chr1",
         "next_reference_start": 0}),
        # singleton
        (("chr1", 0, 5, "+"), None, {}, {})
    ],
)
def test_paired_segments(
        left, right, expected_left_flags, expected_right_flags):
    """Test paired segments class."""
    if right is not None:
        right = align_from_tuple(right)
    left = align_from_tuple(left)
    paired_segments = align_tools.PairedSegments(left, right)
    for k, v in expected_left_flags.items():
        assert getattr(paired_segments.left, k) == v
    for k, v in expected_right_flags.items():
        assert getattr(paired_segments.right, k) == v
    if right is None:
        assert paired_segments.state == align_tools.PairState.singleton
        assert paired_segments.right is None


@pytest.mark.parametrize(
    "aligns,  direct, expected, expected_segments",
    [   
        ([("chr1", 0, 10, "+"), ("chr1", 9, 20, "+"), ("chr2", 0, 10, "+")], False,
         align_tools.PairState.both, 3),
        ([("chr1", 0, 10, "+"), ("chr1", 0, 10, "+"), ("chr1", 9, 20, "+"),
          ("chr1", 9, 20, "+"), ("chr2", 0, 10, "+")], False,
         align_tools.PairState.both, 10),
        ([("chr1", 0, 10, "+"), ("chr1", 0, 10, "+"), ("chr1", 9, 20, "+"),
          ("chr1", 9, 20, "+"), ("chr2", 0, 10, "+")], True,
        align_tools.PairState.both, 4),
        ([], False, None, None),
        ([("chr2", 0, 10, "+")], False, align_tools.PairState.singleton, 1),
        ],
)
def test_get_pairs(aligns,  direct, expected, expected_segments):
    """Test get pairs."""
    aligns = [align_from_tuple(t) for t in aligns]
    pair_segments = align_tools.get_pairs(aligns, direct_only=direct)
    if len(aligns) != 0:
        count = 0
        for i in pair_segments:
            count += 1
            assert i.state == expected
        assert count == expected_segments
    else:
        assert next(pair_segments, None) == expected


@pytest.mark.parametrize(
    "expected",
    [
        ("pore-c-py")
    ],
)
def test_update_header(expected):
    """Test update header."""
    header = align_tools.update_header(HEADER)
    assert header['PG'][-1]['PN'] == expected


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
