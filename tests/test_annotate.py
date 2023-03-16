"""Test functions from annotate."""
import contextlib
import pysam
import pytest
import os
from pore_c_py import annotate, utils


HEADER = pysam.AlignmentHeader.from_references(
    ["chr1", "chr2"], [1000, 2000])


def align_from_tuple(t, name="test", molecule_tag='unknown', walk_tag=None):
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
    rec = utils.MonomerData(molecule_tag, t[4],  t[5], 0,0,0).to_pysam(rec)
    rec.query_name = name
    rec.set_tag(utils.WALK_TAG, walk_tag)
    return rec


def get_data_set(t, name="test", molecule_tag=None, walk_tag=None):
    """Get data set with different categories."""
    primary = align_from_tuple(t, name=name, molecule_tag=molecule_tag)
    unmapped = align_from_tuple(t, name=name, molecule_tag=molecule_tag)
    unmapped.is_unmapped = True
    secondary = align_from_tuple(t, name=name, molecule_tag=molecule_tag)
    secondary.is_secondary = True
    supplementary = align_from_tuple(t, name=name, molecule_tag=molecule_tag)
    supplementary.is_supplementary = True
    alns = [primary, secondary, supplementary, unmapped]  
    return (alns)


@pytest.mark.parametrize(
    "monomer,expected",
    [
        (("chr1", 0, 10, "+", 0, 10), "chr1:+:0-10:0-10"),
        (("chr2", 0, 20, "-", 0, 10), "chr2:-:0-20:0-10"),
        (("chr1", 0, 10, ".", 0, 10), "*:0-10"),
    ],
)
def test_get_walk_component(monomer, expected):
    """Test get monomer tag."""

    res = utils.MonomerData.from_pysam(align_from_tuple(monomer)).name
    assert res == expected


@pytest.mark.parametrize(
    "t",
    [
        (("chr1", 0, 10, "+", 0, 10))
    ],
)
def test_sort_by_category(t):
    """Test sort by category."""
    alns = get_data_set(t)
    primary = alns[0]
    res = annotate.sort_by_category(alns, "test")
    assert next(res) == primary


@pytest.mark.parametrize(
    "alns, expected_monomers, expected_tag",
    [
        ([("chr1", 0, 20, "+", 0, 10), ("chr1", 0, 5, "-", 0, 10)], 2,
         'chr1:+:0-20:0-10;chr1:-:0-5:0-10'),
        ([("chr1", 0, 10, "+", 0, 10), ("chr1", 9, 20, "+", 0, 20),
          ("chr1", 0, 10, ".", 0, 10)], 3, 
          "chr1:+:0-10:0-10;chr1:+:9-20:0-20;*:0-10")
    ]
)
def test_get_walk(alns, expected_monomers, expected_tag):
    """Test get walk."""
    aligns = []
    expected_walk = []
    for i, t in enumerate(alns):
        data_set = get_data_set(t, "test"+str(i), "test"+str(i))
        expected_walk += [data_set[0]]
        aligns += data_set
    walk, n_monomers = annotate.get_walk(aligns)
    assert walk == expected_walk
    assert n_monomers == expected_monomers
    assert walk[0].get_tag(utils.WALK_TAG) == expected_tag


@pytest.mark.parametrize(
    "alns, tag",
    [(
        [("chr1", 0, 20, "+", 0, 10),
         ("chr1", 0, 5, "-", 0, 10)],
        'chr1:+:0-20:0-10;chr1:-:0-5:0-10'),
     (
        [("chr1", 0, 20, "+", 0, 10),
         ("chr1", 0, 0, ".", 0, 10),
         ("chr1", 0, 10, "+", 0, 10)],
        'chr1:+:0-20:0-10;*:0-10;chr1:+:0-10:0-10')]
)
def test_annotate_alignments(alns, tag, tmpdir):
    """Test annotate alignments."""
    aligns = []
    expected_walk = []
    for i, t in enumerate(alns):
        data_set = get_data_set(t, "test"+str(i), "test")
        expected_walk += [align_from_tuple(t, "test"+str(i), "test", tag)]
        aligns += data_set
    bam_path = os.path.join(str(tmpdir), 'ex1.bam')
    with pysam.AlignmentFile(bam_path, mode='wb0', header=HEADER) as f:
        for i in aligns:
            f.write(i)
    with pysam.AlignmentFile(bam_path, "r", check_sq=False) as inbam:
        walks = list(annotate.annotate_alignments(inbam))
        assert len(walks) == 1
        assert walks[0] == expected_walk

#@pytest.mark.parametrize(
#    "start,end,read_length,expected",
#    [
#        (0, 1, 9, "concat001:0:1"),
#        (1, 9, 10, "concat001:01:09"),
#        (1, 10, 10, "concat001:01:10"),
#        # (0, 100, 1000,  "concat001:000"),
#        # (99, 100,10, "concat001:099"),
#    ],
#)
#def test_monomer_id(start, end, read_length, expected):
#    """Test monomer id."""
#    coords = ConcatemerCoords(
#        start=start, end=end,
#        read_length=read_length, subread_idx=0, subread_total=1
#    )
#    res = MonomerReadSeq.generate_id("concat001", coords)
#    assert res == expected
#
#
#def test_group_monomer_aligns(monomer_read_seqs):
#    """Test group monomer aligns."""
#    num_concatemers, aligns_per_concatemer = 2, 10
#    concatemer_ids = set()
#    for concat_id, monomers in group_aligns_by_concatemers(monomer_read_seqs):
#        assert len(list(monomers)) == aligns_per_concatemer
#        concatemer_ids.add(concat_id)
#    assert len(concatemer_ids) == num_concatemers
#    shuffled = sorted(
#        monomer_read_seqs, key=lambda x: x.read_seq.align_info.ref_name)
#    with pytest.raises(ValueError):
#        for concat_id, monomer_read_seqs in group_aligns_by_concatemers(shuffled): # noqa
#            pass
#
#
#def test_annotate_monomer_aligns(monomer_read_seqs):
#    """Test annotate monomer aligns."""
#    res = {
#        concat_id: aligns
#        for (concat_id, aligns) in annotate_monomer_alignments(monomer_read_seqs)   # noqa
#    }
#    assert len(res) == 2
#    for aligns in res.values():
#        assert len(aligns) == 10
#        for a in aligns:
#            a._update_tags()
#            tags = a.read_seq.tags
#            assert MOLECULE_TAG in tags
#            assert WALK_TAG in tags
#            assert CONCATEMER_TAG in tags
#
#
#def test_mods(mock_reads):
#    """Test mods."""
#    test_read = mock_reads["read_w_tags"]
#    header = AlignmentHeader.from_dict(
#        {"RG": [{"ID": "RG1", "SM": "sample01", "LB": "lib01"}]}
#    )
#    seg = test_read.to_align(header)
#    assert seg.modified_bases[("C", 0, "m")] == [(2, 122), (10, 128)]
#    # assert seg.get_tag("RG") == "RG1"
#
#
#def test_split_read(concatemer_unaligned: AlignedSegment):
#    """Test split read."""
#    concatemer = ConcatemerReadSeq.from_align(
#        concatemer_unaligned, as_unaligned=True, init_mod_bases=True
#    )
#    monomers = concatemer.split([5])
#    for idx, s in [(0, slice(None, 5)), (1, slice(5, None))]:
#        assert (
#            monomers[idx].read_seq.sequence
#            == concatemer_unaligned.seq[s]  # type: ignore
#        )
#        assert (
#            monomers[idx].read_seq.quality
#            == concatemer_unaligned.qual[s]  # type: ignore
#        )
#    assert monomers[0].to_align().modified_bases == {  # type: ignore
#        ("C", 0, "m"): [(2, 122)]
#    }
#    assert monomers[1].to_align().modified_bases == {  # type: ignore
#        ("C", 0, "m"): [(5, 128)]
#    }