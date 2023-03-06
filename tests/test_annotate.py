"""Test functions from annotate."""
# flake8: noqa

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