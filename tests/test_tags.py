from typing import List

import pytest
from pysam import AlignedSegment, AlignmentHeader

from pore_c2.model import ConcatemerCoords, downgrade_mm_tag, tag_tuple_to_str
from pore_c2.utils import pysam_verbosity


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
    coords = ConcatemerCoords(start=0, end=10, subread_idx=5, subread_total=10)
    tag = coords.to_tag()
    assert tag == "Xc:B:i,0,10,5,10"
    coords1 = coords.from_tag(tag)
    assert coords == coords1
    # align1 = make_mock_aligned_segment(["Mm:Z:C+m,0,1;", "Ml:B:C,122,128"])


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
    assert align1.modified_bases == {("C", 0, "m"): [(2, 122), (10, 128)]}


@pytest.mark.parametrize(
    "mm_tag",
    ["MM:Z:C+m?,0,1;", "Mm:Z:C+m?,0,1;"],
)
def test_pysam_still_broken(mm_tag):
    with pysam_verbosity(0):
        align = make_mock_aligned_segment([mm_tag, "ML:B:C,122,128"])
        align.modified_bases is None
