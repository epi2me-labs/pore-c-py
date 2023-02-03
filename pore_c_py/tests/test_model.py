"""Test model."""
from pysam import AlignedSegment, AlignmentHeader, FastxRecord
import pytest

from pore_c_py.aligns import (
    annotate_monomer_alignments,
    group_aligns_by_concatemers
)
from pore_c_py.model import (
    AlignInfo,
    ConcatemerCoords,
    ConcatemerReadSeq,
    MonomerReadSeq,
    ReadSeq,
    splits_to_intervals,
    Walk,
    WalkSegment
)
from pore_c_py.sam_utils import (
    CONCATEMER_TAG,
    MOLECULE_TAG,
    SamFlags,
    WALK_TAG
)
from pore_c_py.settings import DEFAULT_ALIGN_HEADER


@pytest.mark.parametrize(
    "length,positions,expected",
    [
        (10, [5], [(0, 5), (5, 10)]),
        (10, [0, 5], [(0, 5), (5, 10)]),
        (10, [0, 5, 10], [(0, 5), (5, 10)]),
        (10, [], [(0, 10)]),
        (10, [0], [(0, 10)]),
        (10, [0, 10], [(0, 10)]),
    ],
)
def test_splits_to_intervals(length, positions, expected):
    """Test splits to intervals."""
    intervals = splits_to_intervals(positions, length)
    assert intervals == expected


@pytest.mark.parametrize(
    "start,end,read_length,expected",
    [
        (0, 1, 9, "concat001:0:1"),
        (1, 9, 10, "concat001:01:09"),
        (1, 10, 10, "concat001:01:10"),
        # (0, 100, 1000,  "concat001:000"),
        # (99, 100,10, "concat001:099"),
    ],
)
def test_monomer_id(start, end, read_length, expected):
    """Test monomer id."""
    coords = ConcatemerCoords(
        start=start, end=end,
        read_length=read_length, subread_idx=0, subread_total=1
    )
    res = MonomerReadSeq.generate_id("concat001", coords)
    assert res == expected


def test_flow(concatemer_unaligned: AlignedSegment):
    """Test flow."""
    concatemer = ConcatemerReadSeq.from_align(
        concatemer_unaligned, as_unaligned=True)
    monomers = concatemer.split([5])
    assert len(monomers) == 2
    assert [_.read_seq.tags[MOLECULE_TAG].rsplit(":", 1)[1] for _ in monomers] == [ # noqa
        "read_w_mods",
        "read_w_mods",
    ]
    assert monomers[0].coords == ConcatemerCoords(
        start=0, end=5, read_length=11, subread_idx=0, subread_total=2
    )
    assert monomers[1].coords == ConcatemerCoords(
        start=5, end=11, read_length=11, subread_idx=1, subread_total=2
    )
    monomers[0].read_seq.align_info = AlignInfo(
        ref_name="chr1", ref_pos=0, length=5, flag=0
    )
    monomers[1].read_seq.align_info = AlignInfo(
        ref_name="chr1", ref_pos=10,
        length=6, flag=SamFlags(reverse=True).to_int()
    )


def test_fastq_round_trip():
    """Test fastq round trip."""
    src = FastxRecord(
        name="read01",
        sequence="ACTG",
        quality="????",
        comment=f"AB:Z:sdfsdfsd\t{MOLECULE_TAG}:Z:read1",
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
    """Concatemer unaglined."""
    src = AlignedSegment.from_dict(
        dict(
            name="read_w_mods",
            seq="AACGTTCGAAC",  # lenght 11
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
    """Concatemer aligned."""
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


@pytest.fixture
def concatemer_aligned_strand() -> AlignedSegment:
    """Concatemer aligned strand."""
    src = AlignedSegment.from_dict(
        dict(
            name="read_w_mods",
            seq="AACGTTCGAAC",
            qual="!!00{}22[]]",
            tags=["RG:Z:RG01", "Mm:Z:C+m,0,1;", "Ml:B:C,122,128"],
            ref_name="chr1",
            ref_pos="100",
            map_quality="30",
            flag="16",
            next_ref_pos="0",
            next_ref_name="*",
            cigar="11M",
            length="0",
        ),
        header=AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]}),
    )
    return src


def test_unaligned_round_trip(concatemer_unaligned: AlignedSegment):
    """Test unaligned round trip."""
    rs = ReadSeq.from_align(concatemer_unaligned)
    dest = rs.to_align()
    assert dest.tostring() == concatemer_unaligned.tostring()  # type: ignore


def test_aligned_round_trip(concatemer_aligned: AlignedSegment):
    """Test aligned round trip."""
    rs = ReadSeq.from_align(concatemer_aligned)
    dest = rs.to_align(
        header=AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]})
    )
    assert dest.tostring() == concatemer_aligned.tostring()


def test_aligned_strand_round_trip(concatemer_aligned_strand: AlignedSegment):
    """Test aligned strand round trip."""
    rs = ReadSeq.from_align(concatemer_aligned_strand)
    dest = rs.to_align(
        header=AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]})
    )
    assert dest.tostring() == concatemer_aligned_strand.tostring()


def test_group_monomer_aligns(monomer_read_seqs):
    """Test group monomer aligns."""
    num_concatemers, aligns_per_concatemer = 2, 10
    concatemer_ids = set()
    for concat_id, monomers in group_aligns_by_concatemers(monomer_read_seqs):
        assert len(list(monomers)) == aligns_per_concatemer
        concatemer_ids.add(concat_id)
    assert len(concatemer_ids) == num_concatemers
    shuffled = sorted(
        monomer_read_seqs, key=lambda x: x.read_seq.align_info.ref_name)
    with pytest.raises(ValueError):
        for concat_id, monomer_read_seqs in group_aligns_by_concatemers(shuffled): # noqa
            pass


def test_annotate_monomer_aligns(monomer_read_seqs):
    """Test annotate monomer aligns."""
    res = {
        concat_id: aligns
        for (concat_id, aligns) in annotate_monomer_alignments(monomer_read_seqs)   # noqa
    }
    assert len(res) == 2
    for aligns in res.values():
        assert len(aligns) == 10
        for a in aligns:
            a._update_tags()
            tags = a.read_seq.tags
            assert MOLECULE_TAG in tags
            assert WALK_TAG in tags
            assert CONCATEMER_TAG in tags


def test_mods(mock_reads):
    """Test mods."""
    test_read = mock_reads["read_w_tags"]
    header = AlignmentHeader.from_dict(
        {"RG": [{"ID": "RG1", "SM": "sample01", "LB": "lib01"}]}
    )
    seg = test_read.to_align(header)
    assert seg.modified_bases[("C", 0, "m")] == [(2, 122), (10, 128)]
    # assert seg.get_tag("RG") == "RG1"


def test_split_read(concatemer_unaligned: AlignedSegment):
    """Test split read."""
    concatemer = ConcatemerReadSeq.from_align(
        concatemer_unaligned, as_unaligned=True, init_mod_bases=True
    )
    monomers = concatemer.split([5])
    for idx, s in [(0, slice(None, 5)), (1, slice(5, None))]:
        assert (
            monomers[idx].read_seq.sequence
            == concatemer_unaligned.seq[s]  # type: ignore
        )
        assert (
            monomers[idx].read_seq.quality
            == concatemer_unaligned.qual[s]  # type: ignore
        )
    assert monomers[0].to_align().modified_bases == {  # type: ignore
        ("C", 0, "m"): [(2, 122)]
    }
    assert monomers[1].to_align().modified_bases == {  # type: ignore
        ("C", 0, "m"): [(5, 128)]
    }


@pytest.mark.parametrize(
    "args,expected",
    [((0, 5), "*:0-5"), ((0, 5, "chr1", 0, 10, "+"), "chr1:+:0-10:0-5")],
)
def test_walk_segment(args, expected):
    """Test walk segment."""
    seg = WalkSegment(*args)
    seg_str = seg.to_string()
    assert seg_str == expected
    seg1 = WalkSegment.from_string(seg_str)
    assert seg == seg1


def test_walk():
    """Test walk."""
    segments = [
        WalkSegment(0, 5),
        WalkSegment(5, 10, "chr1", 10, 15, "+"),
        WalkSegment(10, 30, "chr2", 30, 50, "-"),
        WalkSegment(30, 40),
    ]
    w = Walk(segments)
    tag = f"{WALK_TAG}:Z:*:0-5;chr1:+:10-15:5-10;chr2:-:30-50:10-30;*:30-40"
    assert w.to_tag() == tag
    w1 = Walk.from_tag(tag)
    assert w1 == w
