from typing import List

import mappy as mp
import pytest

from pore_c2.aligns import group_aligns_by_concatemers
from pore_c2.model import AlignData
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


def pair_to_junction(align1: AlignData, align2: AlignData):
    raise ValueError(align1, align2)


def aligns_to_pairs(aligns: List[AlignData]):
    from itertools import combinations

    for pair in combinations(aligns):
        print(pair)


def test_aligns_to_junctions():
    monomer_length = 10
    genome_pos = [("chr1", 0), ("chr1", 100), ("chr1", 110), ("chr2", 0), ("chr1", 50)]
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
    # sort by read coords
    aligns = sorted(aligns, key=lambda x: x.concatemer_metadata.subread_idx)
    # look at all combinations of monomers to create junctions
    #
    for a in aligns:
        print(a)
