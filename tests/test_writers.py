"""Test functions in writers."""
import pytest
import os
from pathlib import Path
from pore_c_py import align_tools, writers, utils
from tests.test_align_tools import align_from_tuple


def monomer_tagged_align(t, name="*", molecule_tag=None):
    """Align from tuple."""
    rec = align_from_tuple(t)
    rec = utils.MonomerData(
       molecule_tag, t[4],  t[5], 0, 0, 0).to_pysam(rec)
    rec.query_name = name
    rec.set_tag(utils.CONCATEMER_ID_TAG, molecule_tag)
    return rec


@pytest.mark.parametrize(
    "aligns, expected",
    [
        ([("chr1", 0, 10, "+", 10, 20), ("chr1", 9, 20, "+", 0, 11),
          ("chr2", 0, 10, "+", 5, 15 )],
         {'cardinality': 3, 'concatemer_count': 3,
          'pair_count': {'both': 3}, 'cis_trans':{'trans': 2,'cis': 1}}),
        ([("chr2", 0, 10, "+", 10, 20)],
         {'cardinality': 1, 'concatemer_count': 1,
          'pair_count':{'singleton':1}, 'cis_trans': {}}),
    ],
)
def test_stats_writer_append(aligns, expected):
    "Test stats writer append."
    aligns = [monomer_tagged_align(t, molecule_tag="test") for t in aligns]
    pair_segments = align_tools.get_pairs(aligns, direct_only=False)
    st = Path("tmp.json")
    stats = writers.StatsWriter(st)
    for i in pair_segments:
        stats.append(i)
    assert stats.cardinality_count[0] == expected['cardinality']
    assert stats.concatemer_count == expected['concatemer_count']
    assert stats.pair_count == expected['pair_count']
    assert stats.cis_trans == expected['cis_trans']


@pytest.mark.parametrize(
    "aligns, expected",
    [
        ([("chr1", 0, 10, "+", 10,20), ("chr1", 0, 20, "+", 0,20), ("chr1", 9, 20, "+", 0,11)], 3),
        ([("chr1", 0, 10, "+", 10,20), ("chr1", 0, 20, ".", 0,20)], 1),
        ([("chr1", 0, 10, "+", 10,20), ("chr1", 0, 20, "+", 0,20)], 2)
    ],
)
def test_chromunity_get_pylist(aligns, expected, tmpdir):
    "Test Chromunity get pylist."
    aligns = [monomer_tagged_align(t, molecule_tag="test") for t in aligns]
    chr_json = os.path.join(str(tmpdir), 'tmp.json')
    chrom = writers.ChromunityWriter(chr_json, merge_distance=5)
    test = chrom.get_pylist(aligns)
    assert test[0]['num_fragments'] == expected
    

@pytest.mark.parametrize(
    "aligns, expected",
    [
        ([("chr1", 0, 8, "+", 0, 8), ("chr1", 50, 60, "+",50, 60),
          ("chr1", 20, 30, "+",20, 30)], 3),
        ([("chr1", 0, 10, "+", 0, 10), ("chr2", 0, 10, "+",8 , 9)], 2),
    ],
)
def test_chromunity_writer(aligns, expected, tmpdir):
    "Test Chromunity writer."
    aligns = [monomer_tagged_align(t, molecule_tag="test") for t in aligns]
    chr_json = os.path.join(str(tmpdir), 'tmp.json')
    chrom = writers.ChromunityWriter(chr_json, merge_distance=5)
    chrom.write(aligns)
    assert chrom.counter == expected