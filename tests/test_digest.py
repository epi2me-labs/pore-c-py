"""Test functions from digest."""
import array
import copy

import pysam
import pytest

from pore_c_py import digest
from array import array
from tests.test_align_tools import align_from_tuple


def align_with_sequence(t, query_sequence=None, name="test"):
    """Align from tuple with specific sequence."""
    rec = align_from_tuple(t)
    if query_sequence is not None:
        rec.query_sequence = query_sequence
    rec.query_name = name
    rec.qual = "&" * (t[2] - t[1])
    return rec

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
    intervals = digest.splits_to_intervals(positions, length)
    assert intervals == expected


# https://github.com/pysam-developers/pysam/blob/b0f1eb7f71b8f5305b112d26da5f72ce79dc33b3/tests/pysam_data/MM-chebi.sam
example1 = (
    ("*", 0, 36, "+"),
    # . = m
    # , = h (76792)
    # ; = m and h
    # | = N
    #      .        | .|,.          .  ;
    "AGCACTCCAGAGTCGNACGCCATYCGCGCGCCACCA",
    "C+m,2,2,1,4,1;C+76792,6,7;N+n,15,2;",
    # .   .   .   .   .   ,   ,   |   |
    [102,128,153,179,204,161,187,212,169])

@pytest.mark.parametrize(
    "name, sam_data, start, end, mm_exp, ml_exp",
    [
        # the whole read
        ("whole", example1, 0, 36, "C+m,2,2,1,4,1;C+76792,6,7;N+n,15,2;",
         [102, 128, 153, 179, 204, 161, 187, 212, 169]),

        # trim two non-C bases.
        # Doesn't effect C mods, rejigs coords of N
        ("trim-0C", example1, 2, 36, "C+m,2,2,1,4,1;C+76792,6,7;N+n,13,2;",
         [102, 128, 153, 179, 204, 161, 187, 212, 169]),

        # trim first C
        # Need to skip one less C for C mod subtag and for N (since its any base)
        ("trim-1C", example1, 3, 36, "C+m,1,2,1,4,1;C+76792,5,7;N+n,12,2;",
         [102, 128, 153, 179, 204, 161, 187, 212, 169]),

        # trim first two C
        # We now skip 0 Cs before first mod.
        ("trim-2C", example1, 5, 36, "C+m,0,2,1,4,1;C+76792,4,7;N+n,10,2;",
         [102, 128, 153, 179, 204, 161, 187, 212, 169]),

        # trim first three C
        # We've trimmed a mod base so drop the first entry from ML. Also
        ("trim-3C", example1, 7, 36, "C+m,2,1,4,1;C+76792,3,7;N+n,8,2;",
         [128, 153, 179, 204, 161, 187, 212, 169]),

        # trim last C
        # Affects both m and h codes
        ("trim+1C", example1, 0, 34, "C+m,2,2,1,4;C+76792,6;N+n,15,2;",
         [102, 128, 153, 179, 161, 212, 169]),

        # trim past first N
        # Affects both m and h codes also
        ("trim-5C", example1, 17, 36, "C+m,0,1,4,1;C+76792,1,7;N+n,1;",
         [128,153,179,204,161,187,169]),

        # trim past first N and last C
        # Affects both m and h codes also
        ("trim-5C+1C", example1, 17, 34, "C+m,0,1,4;C+76792,1;N+n,1;",
         [128,153,179,161,169]),

    ],
)
def test_get_subread_modified_bases(
        name, sam_data, start, end, mm_exp, ml_exp):
    """Test get subread."""
    # note: the representation isn't unique and the code doesn't produce the same
    #       representation as the test data. So we parse both using pysam and compare

    # create original
    aln_tuple, qseq, mm, ml = sam_data
    aln = align_with_sequence(aln_tuple, query_sequence=qseq, name="*")
    aln.set_tag("Mm", mm)
    aln.set_tag("Ml", ml)

    # create trimmed
    read = copy.copy(aln)
    mm, ml = digest.get_subread_modified_bases(aln, start, end)
    read.set_tag("Mm", mm)
    read.set_tag("Ml", ml)

    # original with expected tags
    expected_aln = align_with_sequence(aln_tuple, query_sequence=qseq, name="*")
    expected_aln.set_tag("Mm", mm_exp)
    expected_aln.set_tag("Ml", ml_exp)

    # compare
    for k, v in read.modified_bases.items():
        assert expected_aln.modified_bases[k] == v

    # check that ML tag is a byte array
    assert isinstance(read.get_tag("Ml"), array)
    found = sum(x.startswith('Ml:B:C') for x in read.to_dict()["tags"])
    assert found == 1


@pytest.mark.parametrize(
    "align, sequence, enzyme, expected",
    [
        (("chr1", 0, 19, "+"), 'AAACATGTTTGGCATGAAA', "NlaIII",
         ['AAACATG', 'TTTGGCATG', 'AAA']),
        (("chr1", 0, 10, "+"), 'GGATCTGATC', "DpnII",
         ['G', 'GATCT', 'GATC'])
    ],
)
def test_digest_sequence(align, sequence, enzyme, expected):
    aln = align_with_sequence(align, query_sequence=sequence)
    enz = digest.get_enzyme(enzyme)
    seqs = digest.digest_sequence(aln, enz)
    for i, (monomer, gen) in enumerate(seqs):
        assert monomer
        assert gen.query_sequence == expected[i]


@pytest.mark.parametrize(
    "align, sequence, MM, ML, expected_MN, expected_ML, expected_MM",
    [(
         ("*", 0, 43, "+"),
         "AGCACTCCAGAGTCGNACGCCATYCGCGCGCCACCACATGTTT",
         "C+m,2,2,1,4,1;C+76792,6,7;N+n,15,2;",
         [102,128,153,179,204,161,187,212,169],
         [40, 0], 
         [[102, 128, 153, 179, 204, 212, 169, 161, 187]],
         ['C+m,2,2,1,4,1;N+n,15,2;C+76792,6,7;']),
     (
         ("*", 0, 43, "+"),
         "AGCACTCCAGAGTCGNACGCCATYCGCGCGCCACCACATGTTC",
         "C+m,2,2,1,4,1,1;C+76792,6,7;N+n,15,2;",
         [102,128,153,179,204,161,187,212,169,100],
         [40, 3],
         [[102, 128, 153, 179, 204, 169, 100, 187, 212],
          [161]],
         ['C+m,2,2,1,4,1;N+n,15,2;C+76792,6,7;', 'C+m,0;'])
    ],
)
def test_digest_sequence_mods(align, sequence, MM, ML, expected_MN, expected_ML, expected_MM):
    # Create aln
    aln = align_with_sequence(align, query_sequence=sequence)
    aln.set_tag("MN", len(sequence))
    aln.set_tag("Mn", len(sequence))
    aln.set_tag("Ml", ML)
    aln.set_tag("Mm", MM)

    # Get digested with modified bases
    enz = digest.get_enzyme("NlaIII")
    seqs = digest.digest_sequence(aln, enz)
 
    # Compare 
    for i, (monomer, gen) in enumerate(seqs):
        assert gen.has_tag("Mn") == False
        if gen.has_tag("MN"):
            assert gen.get_tag("MN") == expected_MN[i]
            assert gen.get_tag("ML") == array('B', expected_ML[i])
            assert gen.get_tag("MM") == expected_MM[i]
        else:
            assert gen.has_tag("MN") == False
        
        


@pytest.mark.parametrize(
    "align, sequence, enzyme",
    [
        (("chr1", 0, 19, "+"), 'AAACATGTTTGGCATGAAA', "NlaIII"),
    ],
)
def test_remove_mv_tag(align, sequence, enzyme):
    aln = align_with_sequence(align, query_sequence=sequence)
    aln.set_tag('mv', [0,1,0,1,0,1])
    enz = digest.get_enzyme(enzyme)
    seqs = digest.digest_sequence(aln, enz)
    for (monomer, gen) in seqs:
        assert not gen.has_tag('mv')

    # check also we don't blow up removing tag that doesn't exist
    seqs = digest.digest_sequence(aln, enz, remove_tags={'thing'})
    for (monomer, gen) in seqs:
        assert not gen.has_tag('mv')

