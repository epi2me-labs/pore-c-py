import pytest
from pysam import AlignedSegment, AlignmentHeader, FastxRecord

from pore_c2.aligns import annotate_monomer_alignments, group_aligns_by_concatemers
from pore_c2.model import (
    AlignInfo,
    ConcatemerCoords,
    ConcatemerReadSeq,
    MonomerReadSeq,
    ReadSeq,
    Walk,
    WalkSegment,
)
from pore_c2.settings import DEFAULT_ALIGN_HEADER, WALK_TAG
from pore_c2.utils import SamFlags


@pytest.mark.parametrize(
    "idx,total,expected",
    [
        (0, 1, "concat001:0"),
        (1, 9, "concat001:1"),
        (1, 10, "concat001:01"),
        (0, 100, "concat001:000"),
        (99, 100, "concat001:099"),
    ],
)
def test_monomer_id(idx, total, expected):
    coords = ConcatemerCoords(start=0, end=1, subread_idx=idx, subread_total=total)
    res = MonomerReadSeq.generate_id("concat001", coords)
    assert res == expected


def test_flow(concatemer_unaligned: AlignedSegment):
    concatemer = ConcatemerReadSeq.from_align(concatemer_unaligned, as_unaligned=True)
    monomers = concatemer.split([5])
    assert len(monomers) == 2
    assert [_.read_seq.tags["MI"].rsplit(":", 1)[1] for _ in monomers] == [
        "read_w_mods",
        "read_w_mods",
    ]
    assert monomers[0].coords == ConcatemerCoords(
        start=0, end=5, subread_idx=0, subread_total=2
    )
    assert monomers[1].coords == ConcatemerCoords(
        start=5, end=11, subread_idx=1, subread_total=2
    )
    monomers[0].read_seq.align_info = AlignInfo(
        ref_name="chr1", ref_pos=0, length=5, flag=0
    )
    monomers[1].read_seq.align_info = AlignInfo(
        ref_name="chr1", ref_pos=10, length=6, flag=SamFlags(reverse=True).to_int()
    )


def test_fastq_round_trip():
    src = FastxRecord(
        name="read01",
        sequence="ACTG",
        quality="????",
        comment="AB:Z:sdfsdfsd\tMI:Z:read1",
    )
    rs = ReadSeq.from_fastq(src)
    assert rs.name == src.name
    assert rs.sequence == src.sequence
    assert rs.quality == src.quality

    dest = rs.to_fastq()
    assert dest.name == src.name
    assert dest.sequence == src.sequence
    assert dest.quality == src.quality
    assert dest.comment == src.comment


@pytest.fixture
def concatemer_unaligned() -> AlignedSegment:
    src = AlignedSegment.from_dict(
        dict(
            name="read_w_mods",
            seq="AACGTTCGAAC",
            qual="!!00{}22[]]",
            tags=["RG:Z:RG01", "Mm:Z:C+m,0,1;", "Ml:B:C,122,128"],
            ref_name="*",
            ref_pos="0",
            map_quality="0",
            flag="4",
            next_ref_pos="0",
            next_ref_name="*",
            cigar="*",
            length="0",
        ),
        header=DEFAULT_ALIGN_HEADER,
    )
    return src


@pytest.fixture
def concatemer_aligned() -> AlignedSegment:
    src = AlignedSegment.from_dict(
        dict(
            name="read_w_mods",
            seq="AACGTTCGAAC",
            qual="!!00{}22[]]",
            tags=["RG:Z:RG01", "Mm:Z:C+m,0,1;", "Ml:B:C,122,128"],
            ref_name="chr1",
            ref_pos="100",
            map_quality="30",
            flag="0",
            next_ref_pos="0",
            next_ref_name="*",
            cigar="11M",
            length="0",
        ),
        header=AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]}),
    )
    return src


def test_unaligned_round_trip(concatemer_unaligned: AlignedSegment):
    rs = ReadSeq.from_align(concatemer_unaligned)
    assert rs.name == concatemer_unaligned.query_name
    assert rs.sequence == concatemer_unaligned.query_sequence
    assert rs.quality == concatemer_unaligned.qual
    dest = rs.to_align()
    assert dest.tostring() == concatemer_unaligned.tostring()


def test_aligned_round_trip(concatemer_aligned: AlignedSegment):
    rs = ReadSeq.from_align(concatemer_aligned)
    assert rs.name == concatemer_aligned.query_name
    assert rs.sequence == concatemer_aligned.query_sequence
    assert rs.quality == concatemer_aligned.qual
    dest = rs.to_align(
        header=AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]})
    )
    assert dest.tostring() == concatemer_aligned.tostring()


def test_group_monomer_aligns(monomer_read_seqs):
    num_concatemers, aligns_per_concatemer = 2, 10
    concatemer_ids = set()
    for concat_id, monomers in group_aligns_by_concatemers(monomer_read_seqs):
        assert len(list(monomers)) == aligns_per_concatemer
        concatemer_ids.add(concat_id)
    assert len(concatemer_ids) == num_concatemers
    shuffled = sorted(monomer_read_seqs, key=lambda x: x.read_seq.align_info.ref_name)
    with pytest.raises(ValueError):
        for concat_id, monomer_read_seqs in group_aligns_by_concatemers(shuffled):
            pass


def test_annotate_monomer_aligns(monomer_read_seqs):
    res = {
        concat_id: aligns
        for (concat_id, aligns) in annotate_monomer_alignments(monomer_read_seqs)
    }
    assert len(res) == 2
    for aligns in res.values():
        assert len(aligns) == 10
        for a in aligns:
            a._update_tags()
            tags = a.read_seq.tags
            assert "MI" in tags
            assert WALK_TAG in tags
            assert "Xc" in tags


def test_mods(mock_reads):
    test_read = mock_reads["read_w_tags"]
    header = AlignmentHeader.from_dict(
        {"RG": [{"ID": "RG1", "SM": "sample01", "LB": "lib01"}]}
    )
    seg = test_read.to_align(header)
    assert seg.modified_bases[("C", 0, "m")] == [(2, 122), (10, 128)]
    # assert seg.get_tag("RG") == "RG1"


def test_split_read(concatemer_unaligned: AlignedSegment):
    concatemer = ConcatemerReadSeq.from_align(
        concatemer_unaligned, as_unaligned=True, init_mod_bases=True
    )
    monomers = concatemer.split([5])
    for idx, s in [(0, slice(None, 5)), (1, slice(5, None))]:
        assert monomers[idx].read_seq.sequence == concatemer_unaligned.seq[s]
        assert monomers[idx].read_seq.quality == concatemer_unaligned.qual[s]
    assert monomers[0].to_align().modified_bases == {("C", 0, "m"): [(2, 122)]}
    assert monomers[1].to_align().modified_bases == {("C", 0, "m"): [(5, 128)]}


@pytest.mark.parametrize(
    "args,expected",
    [((0, 5), "*:0-5"), ((0, 5, "chr1", 0, 10, "+"), "chr1:+:0-10:0-5")],
)
def test_walk_segment(args, expected):
    seg = WalkSegment(*args)
    seg_str = seg.to_string()
    assert seg_str == expected
    seg1 = WalkSegment.from_string(seg_str)
    assert seg == seg1


def test_walk():
    segments = [
        WalkSegment(0, 5),
        WalkSegment(5, 10, "chr1", 10, 15, "+"),
        WalkSegment(10, 30, "chr2", 30, 50, "-"),
        WalkSegment(30, 40),
    ]
    w = Walk(segments)
    tag = "Xc:Z:*:0-5;chr1:+:10-15:5-10;chr2:-:30-50:10-30;*:30-40"
    assert w.to_tag() == tag
    w1 = Walk.from_tag(tag)
    assert w1 == w
