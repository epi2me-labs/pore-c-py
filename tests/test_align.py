import mappy as mp
import pytest

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
