"""Test align."""
from itertools import combinations
import logging
from typing import List, Optional, Tuple

import pytest

from pore_c_py.aligns import (
    annotate_monomer_alignments,
    calculate_genomic_distance,
    get_pairs,
    group_colinear,
    is_colinear,
    PairAlignState,
    PairData
)
from pore_c_py.model import (
    AlignInfo,
    ConcatemerCoords,
    MonomerReadSeq,
    ReadSeq
)
from pore_c_py.sam_utils import SamFlags


def _make_mock_monomer_aligns(
    data: List[Tuple[int, int, str, int, Optional[SamFlags]]]
) -> List[MonomerReadSeq]:
    """Make mock monomer aligns."""
    monomers = []
    concatemer_id = "C"
    subread_total = len(set([_[0] for _ in data]))
    read_length = sum([_[1] - _[0] for _ in data])
    for x, (r_start, r_end, chrom, g_start, flag) in enumerate(data):
        if flag is None:
            flag = SamFlags(unmap=chrom is None)
        length = r_end - r_start
        monomer_id = f"{concatemer_id}:{r_start}:{r_end}"
        m = MonomerReadSeq(
            concatemer_id=concatemer_id,
            monomer_id=monomer_id,
            coords=ConcatemerCoords(
                start=r_start,
                end=r_end,
                read_length=read_length,
                subread_idx=x,
                subread_total=subread_total,
            ),
            read_seq=ReadSeq(
                name=monomer_id,
                sequence="",
                quality="",
                flags=flag,
                align_info=AlignInfo(
                    ref_name=chrom,
                    ref_pos=g_start,
                    length=length,
                    flag=flag.to_int(),
                )
                if chrom
                else None,
            ),
        )
        m._update_tags()
        monomers.append(m)

    return monomers


@pytest.fixture
def m1():
    """Make mock monomer aligns m1."""
    return _make_mock_monomer_aligns(
        [
            (0, 5, "chr1", 0, None),
            (5, 10, "chr1", 100, None),
            (5, 10, "chr2", 100, None),
            (10, 15, "chr1", 50, None),
        ]
    )


def test_qcfail_supplementary(m1):
    """Test QC fail supplementary."""
    m1[2].read_seq.flags.supplementary = True
    m = list(annotate_monomer_alignments(m1, remove_qcfail=False))[0][1]
    qcfail = [_.read_seq.flags.qcfail for _ in m]
    assert qcfail == [False, False, True, False]


def test_qcfail_secondary(m1, caplog):
    """Test QC fail secondary."""
    m1[1].read_seq.flags.supplementary = True
    m1[2].read_seq.flags.secondary = True
    with caplog.at_level(logging.WARNING):
        m = list(annotate_monomer_alignments(m1, remove_qcfail=False))[0][1]
        qcfail = [_.read_seq.flags.qcfail for _ in m]
        assert qcfail == [False, False, True, False]
        assert (
            caplog.record_tuples[-1][2]
            == "Warning: best alignment for monomer: C:5:10 has category supplementary" # noqa
        )


def test_qcfail_two_supplementary(m1, caplog):
    """Test QCfail two supplementary."""
    m1[1].read_seq.flags.secondary = True
    m1[2].read_seq.flags.secondary = True
    with caplog.at_level(logging.WARNING):
        m = list(
            annotate_monomer_alignments(m1, remove_qcfail=False))[0][1]
        qcfail = [_.read_seq.flags.qcfail for _ in m]
        assert qcfail == [False, False, True, False]
        assert (
            caplog.record_tuples[-1][2]
            == "Warning: best alignment for monomer: C:5:10 has category secondary" # noqa
        )


def test_short_walk(m1, caplog):
    """Test short walk."""
    m1[2].read_seq.flags.secondary = True
    with caplog.at_level(logging.WARNING):
        m = list(
            annotate_monomer_alignments(m1[:-1], remove_qcfail=False))[0][1]
        qcfail = [_.read_seq.flags.qcfail for _ in m]
        assert qcfail == [False, False, True]
        assert (
            caplog.record_tuples[-1][2]
            == "Expected to see 3 alignments for concatemer C, found 2"
        )


def test_to_sam(m1):
    """Test to sam."""
    m = list(annotate_monomer_alignments(m1))[0][1]
    # baseline
    assert (
        m[0].read_seq.to_sam(strip_tags=True)
        == "C:0:5\t0\tchr1\t1\t0\t*\t*\t0\t0\t\t\t"
    )
    # alter read name at write time
    assert (
        m[0].read_seq.to_sam(read_name="blah", strip_tags=True)
        == "blah\t0\tchr1\t1\t0\t*\t*\t0\t0\t\t\t"
    )
    # strip aligmnent info
    assert (
        m[0].read_seq.to_sam(as_unaligned=True, strip_tags=True)
        == "C:0:5\t4\t*\t0\t0\t*\t*\t0\t0\t\t\t"
    )
    # strip aligmnent info
    assert (
        m[0].read_seq.to_sam(as_unaligned=True, strip_tags=True)
        == "C:0:5\t4\t*\t0\t0\t*\t*\t0\t0\t\t\t"
    )
    # alter flag
    assert (
        m[0].read_seq.to_sam(flag=SamFlags.from_int(8), strip_tags=True)
        == "C:0:5\t8\tchr1\t1\t0\t*\t*\t0\t0\t\t\t"
    )


@pytest.fixture
def m2():
    """Make mock monomer aligns m2."""
    return (
        _make_mock_monomer_aligns(
            [
                (0, 5, "chr1", 0, None),
                (5, 10, "chr2", 100, None),
                (20, 25, "chr1", 50, None),
                (25, 30, "chr1", 500, None),
            ]
        ),
        [
            PairData(*_)
            for _ in [
                ("C", 0, True, 0, PairAlignState.both, False, None, None),
                ("C", 1, True, 10, PairAlignState.both, False, None, None),
                ("C", 2, True, 0, PairAlignState.both, True, 445, True),
            ]
        ],
    )


def test_pairs_direct(m2):
    """Test pairs direct."""
    aligns, expected = m2
    _, _, pair_data = zip(*get_pairs(aligns, direct_only=True))
    assert list(pair_data) == expected


@pytest.fixture
def m3():
    """Make mock monomer aligns m3."""
    return (
        _make_mock_monomer_aligns(
            [
                (0, 5, "chr1", 0, None),
                (5, 10, "chr2", 100, None),
                (20, 25, "chr1", 50, None),
            ]
        ),
        [
            PairData(*_)
            for _ in [
                ("C", 0, True, 0, PairAlignState.both, False, None, None),
                ("C", 1, False, 15, PairAlignState.both, True, 45, True),
                ("C", 2, True, 10, PairAlignState.both, False, None, None),
            ]
        ],
    )


def test_pairs_indirect(m3):
    """Test pairs indirect."""
    aligns, expected = m3
    _, _, pair_data = zip(*get_pairs(aligns, direct_only=False))
    assert list(pair_data) == expected


def test_relative_ori():
    """Test relative ori."""
    flags = [SamFlags(unmap=True), SamFlags(reverse=True), SamFlags()]
    for f1, f2 in combinations(flags, 2):
        m = _make_mock_monomer_aligns([
            (0, 5, "chr1", 0, f1), (5, 10, "chr1", 10, f2)])
        _, _, pair_data = zip(*get_pairs(m, direct_only=False))
        ori = pair_data[0].relative_ori  # type: ignore
        strand_pair = (f1.strand, f2.strand)
        if strand_pair in [("+", "+"), ("-", "-")]:
            assert ori is True
        elif strand_pair in [("+", "-"), ("-", "+")]:
            assert ori is False
        else:
            assert ori is None


@pytest.mark.parametrize(
    "left,right,expected",
    [
        (("+", 0, 5), ("+", 10, 15), 5),  # span 5-10
        (("-", 10, 15), ("-", 0, 5), -5),  # span 5-10
        (("-", 0, 5), ("+", 10, 15), 10),  # span 0-10
        (("-", 10, 15), ("+", 0, 5), -10),  # span 0-10
        (("+", 0, 5), ("-", 10, 15), 10),  # span 5-15
        (("+", 10, 15), ("-", 0, 5), -10),  # span 5-15
        (("-", 0, 5), ("-", 10, 15), 15),  # span 0-15
        (("+", 10, 15), ("+", 0, 5), -15),  # span 0-15
    ],
)
def test_genomic_distance(left, right, expected):
    """Test genomic distance."""
    assert calculate_genomic_distance(left, right) == expected


def test_read_pair_to_sam():
    """Test read pair to sam."""
    for right_chrom in ["chr1", "chr2", None]:
        left, right = _make_mock_monomer_aligns(
            [(0, 5, "chr1", 0, None), (5, 10, right_chrom, 100, None)]
        )
        pair_data = PairData.from_monomer_pair(left, right)
        l_flags, r_flags = pair_data.to_flags(
            (left.read_seq.flags, right.read_seq.flags)
        )
        assert l_flags.paired is True
        assert r_flags.paired is True
        if right_chrom == "chr1":
            assert l_flags == SamFlags(
                proper_pair=True, paired=True, read1=True)
            assert r_flags == SamFlags(
                proper_pair=True, paired=True, read2=True)
            assert pair_data.is_cis is True
            assert pair_data.genome_distance == 95
        elif right_chrom is None:
            assert l_flags == SamFlags(
                proper_pair=False, paired=True, read1=True, munmap=True
            )
            assert r_flags == SamFlags(
                proper_pair=False, paired=True, read2=True, unmap=True
            )
            assert pair_data.is_cis is None
            assert pair_data.genome_distance is None
        else:
            assert l_flags == SamFlags(
                proper_pair=False, paired=True, read1=True)
            assert r_flags == SamFlags(
                proper_pair=False, paired=True, read2=True)
            assert pair_data.is_cis is False
            assert pair_data.genome_distance is None

    # print(pair_data)
    # raise ValueError(pair_data)


def _align_from_tuple(t) -> AlignInfo:
    """Align from tuple."""
    return AlignInfo(
        ref_name=t[0],
        ref_pos=t[1],
        length=t[2] - t[1],
        flag=SamFlags(
            reverse=t[3] == "-", unmap=len(t) > 4 and t[4] is False).to_int(),
    )


@pytest.mark.parametrize(
    "left,right,tol,expected",
    [
        (("chr1", 0, 10, "+"), ("chr1", 0, 20, "+"), 0, True),  # overlapping
        (("chr1", 0, 10, "+"), ("chr1", 0, 20, "-"), 0, False),
        (("chr1", 0, 10, "-"), ("chr1", 0, 20, "-"), 0, True),
        (("chr1", 0, 10, "+"), ("chr2", 0, 20, "+"), 0, False),
        (
            ("chr1", 0, 10, "+"),
            ("chr1", 10, 20, "+"),
            1,  # need tol=1 for bookended features
            True,
        ),
        (
            ("chr1", 0, 10, "+"),
            ("chr1", 5, 20, "+"),
            5,  # allow some wiggle
            True,
        ),
        (
            ("chr1", 0, 10, "+"),
            ("chr1", 0, 0, ".", False),
            0,
            False,
        ),  # unaligned not colinear
    ],
)
def test_is_colinear(left, right, tol, expected):
    """Test is colinear."""
    res = is_colinear(
        _align_from_tuple(left),
        _align_from_tuple(right),
        tol=tol,
    )
    assert res == expected


@pytest.mark.parametrize(
    "aligns,expected",
    [
        [
            [("chr1", 0, 10, "+"), ("chr1", 9, 20, "+"), ("chr2", 0, 10, "+")],
            [[0, 1], [2]],
        ],
        [
            [
                ("chr1", 0, 10, "+"),
                ("chr1", 9, 20, "+", False),
                ("chr1", 20, 30, "+")],
            [[0], [1], [2]],
        ],
    ],
)
def test_group_colinear(aligns, expected):
    """Test group colinear."""
    res = group_colinear([_align_from_tuple(t) for t in aligns])
    assert res == expected
