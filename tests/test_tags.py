from typing import List

import pytest
from attrs import fields_dict
from pysam import AlignedSegment, AlignmentHeader

from pore_c2.model import ConcatemerCoords
from pore_c2.sam_utils import (
    CONCATEMER_TAG,
    SamEnum,
    SamFlags,
    downgrade_mm_tag,
    pysam_verbosity,
    tag_tuple_to_str,
)


def make_mock_aligned_segment(tags: List[str]):
    return AlignedSegment.from_dict(
        dict(
            name="read_w_mods",
            seq="AACGTTCGAAC",
            qual="!!00{}22[]]",
            tags=tags,
            ref_name="*",
            ref_pos="0",
            map_quality="0",
            flag="4",
            next_ref_pos="0",
            next_ref_name="*",
            cigar="*",
            length="0",
        ),
        header=AlignmentHeader.from_dict({}),
    )


def test_concat_coords_tag():
    coords = ConcatemerCoords(
        start=0, end=10, read_length=20, subread_idx=5, subread_total=10
    )
    tag = coords.to_tag()
    assert tag == f"{CONCATEMER_TAG}:B:i,0,10,20,5,10"
    coords1 = coords.from_tag(tag)
    assert coords == coords1


def test_mod_tags():
    align1 = make_mock_aligned_segment(["Mm:Z:C+m,0,1;", "Ml:B:C,122,128"])
    assert align1.modified_bases == {("C", 0, "m"): [(2, 122), (10, 128)]}
    align2 = make_mock_aligned_segment(["MM:Z:C+m,0,1;", "ML:B:C,122,128"])
    assert align1.modified_bases == align2.modified_bases


@pytest.mark.parametrize(
    "tag_str",
    ["Ml:B:C,1,2,3", "mv:B:c,5,1,1,1", "xy:Z:blah", "MM:Z:C+m?,0,1;", "NM:i:0"],
)
def test_tag_roundtrip(tag_str):
    align = make_mock_aligned_segment([tag_str])
    tag = align.get_tags(with_value_type=True)[0]
    res = tag_tuple_to_str(*tag)
    assert res == tag_str


@pytest.mark.parametrize("tag_str", ["Ml:B:C,1,2,3", "mv:B:c,5,1,1,1", "xy:Z:blah"])
def test_tag_roundtrip_w_downgrade(tag_str):
    align = downgrade_mm_tag(make_mock_aligned_segment([tag_str]))
    tag = align.get_tags(with_value_type=True)[0]
    res = tag_tuple_to_str(*tag)
    assert res == tag_str


@pytest.mark.parametrize(
    "mm_tag", ["MM:Z:C+m?,0,1;", "Mm:Z:C+m?,0,1;", "Mm:Z:C+m,0,1;"]
)
def test_fix_mod_tags(mm_tag):
    align = make_mock_aligned_segment([mm_tag, "ML:B:C,122,128"])
    align1 = downgrade_mm_tag(align)
    assert align1.modified_bases == {  # type: ignore
        ("C", 0, "m"): [(2, 122), (10, 128)]
    }


@pytest.mark.parametrize(
    "mm_tag",
    ["MM:Z:C+m?,0,1;", "Mm:Z:C+m?,0,1;"],
)
def test_pysam_still_broken(mm_tag):
    with pysam_verbosity(0):
        align = make_mock_aligned_segment([mm_tag, "ML:B:C,122,128"])
        assert align.modified_bases is None


def test_sam_flags():
    for key, _ in fields_dict(SamFlags).items():  # pyright: ignore
        for tf in (True, False):
            flags = SamFlags(**{key: tf})
            assert getattr(flags, key) == tf
            integer_val = flags.to_int()
            _ = SamFlags.from_int(integer_val)
            assert _ == flags

    flags = SamFlags.from_int(3844)
    assert flags == SamFlags(
        unmap=True, secondary=True, qcfail=True, dup=True, supplementary=True
    )


# TODO: find a home for this
def is_primary(flag: int):
    return bool(flag & ~(SamEnum.supplementary | SamEnum.secondary | SamEnum.unmap))


def test_mask():
    f = SamFlags(unmap=True, secondary=True).to_int()
    assert bool(f & SamEnum.unmap)
    # assert(bool(SamEnum['unmap'].value & f.to_int()) is True)


def test_align_categories():
    flags = [
        SamFlags(unmap=True),
        SamFlags(secondary=True),
        SamFlags(supplementary=True),
        SamFlags(),
    ]
    categories = [_.category for _ in flags]
    assert [c.name for c in categories] == [
        "unmapped",
        "secondary",
        "supplementary",
        "primary",
    ]
    assert [c.name for c in sorted(categories)] == [
        "primary",
        "unmapped",
        "supplementary",
        "secondary",
    ]


def test_strand():
    flags = [
        SamFlags(),
        SamFlags(reverse=True),
        SamFlags(unmap=True),
        SamFlags(secondary=True),
    ]
    strands = [SamFlags.int_to_strand(_.to_int()) for _ in flags]
    assert strands == ["+", "-", ".", "+"]
