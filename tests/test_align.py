import mappy as mp
import pytest

from pore_c2.aligns import annotate_monomer_alignments, group_aligns_by_concatemers
from pore_c2.model import AlignData
from pore_c2.overlaps import FragmentOverlapper
from pore_c2.testing import Scenario, simulate_read_sequence


@pytest.mark.parametrize(
    "chrom,start,end,cigar,NM,offsets,alleles",
    [
        ("chr1", 10, 100, "90M", 0, None, None),  # perfect match
        (
            "chr1",
            10,
            100,
            "90M",
            1,
            [10],
            ["G"],
        ),  # SNP  TODO: this test might fail if genome has a G
        (
            "chr1",
            10,
            100,
            "10M1I80M",
            2,
            [10],
            ["GG"],
        ),  # SNP  TODO: this test might fail if genome has a G
    ],
)
def test_mappy_snvs(
    default_scenario: Scenario, chrom, start, end, cigar, NM, offsets, alleles
):
    aligner = mp.Aligner(str(default_scenario.reference_fasta), preset="map-ont")
    # chrom, start, end, cigar, NM = "chr1", 10, 100, "90M", 1
    read, _ = simulate_read_sequence(
        ff=default_scenario.ff,
        chrom=chrom,
        start=start,
        end=end,
        offsets=offsets,
        allele_values=alleles,
    )
    hits = list(aligner.map(read))
    assert len(hits) == 1
    hit = hits[0]
    assert hit.ctg == chrom
    assert hit.r_st == start
    assert hit.cigar_str == cigar
    assert hit.NM == NM


def test_group_aligns():
    num_concatemers, aligns_per_concatemer = 2, 10
    aligns = [
        AlignData(
            name=f"CONCAT{y}:{x}",
            seq="AATGC",
            ref_name=f"chrom{x+1}",
            ref_pos=10,
            tags=[f"MI:Z:CONCAT{y}"],
        )
        for y in range(num_concatemers)
        for x in range(aligns_per_concatemer)
    ]
    concatemer_ids = set()
    for concat_id, monomers in group_aligns_by_concatemers(aligns):
        assert len(list(monomers)) == aligns_per_concatemer
        concatemer_ids.add(concat_id)
    assert len(concatemer_ids) == num_concatemers
    shuffled = sorted(aligns, key=lambda x: x.ref_name)
    with pytest.raises(ValueError):
        for concat_id, aligns in group_aligns_by_concatemers(shuffled):
            pass


def test_aligns_to_junctions():
    monomer_length = 10

    # TODO: annotate fragment overlaps
    overlapper = FragmentOverlapper(
        left={
            "chr1": [0, 10, 100, 110, 120],
            "chr2": [0, 10, 100, 110, 120],
        },
        ids={
            "chr1": ["f11", "f12", "f13", "f14", "f15"],
            "chr2": ["f21", "f22", "f23", "f24", "f25"],
        },
    )  # noqa
    genome_pos = [("chr1", 0), ("chr1", 100), ("chr1", 110), ("chr2", 0), ("chr1", 50)]
    expected_walk = ",".join(
        [f"{chrom}:{pos}_{pos+monomer_length}" for chrom, pos in genome_pos]
    )
    num_monomers = len(genome_pos)
    aligns = [
        AlignData(
            name=f"read1:{x}",
            flag=0,
            seq="ATGC",
            ref_name=chrom,
            ref_pos=start,
            length=monomer_length,
            tags=[
                "MI:Z:read1",
                f"Xc:B:i,{x*monomer_length},{(x+1)*monomer_length},{x},{num_monomers}",
            ],
        )
        for x, (chrom, start) in enumerate(genome_pos)
    ]
    res = list(annotate_monomer_alignments(aligns[::-1]))
    assert len(res) == 1
    concat_id, annotated_aligns = res[0]
    assert concat_id == "read1"
    # check that they're sorted
    assert [(a.ref_name, a.ref_pos) for a in annotated_aligns] == genome_pos
    # they all have the correct walk
    assert [a.get_walk() for a in annotated_aligns] == [expected_walk] * len(aligns)
    # they all point to the next monomer in the walk
    for x, (chrom, pos) in enumerate(genome_pos[1:]):
        a = annotated_aligns[x]
        assert a.next_ref_name in (chrom, "=")
        assert a.next_ref_pos == pos
    # the last monomer in the concatemer should not have a next read
    assert annotated_aligns[-1].next_ref_name == "*"
    assert annotated_aligns[-1].next_ref_pos == 0

    # for concat_id, annotated_aligns in annotate_monomer_alignments(aligns):
    #    print(concat_id, annotated_aligns)
    #    print([a.next_ref_pos for a in annotated_aligns])
