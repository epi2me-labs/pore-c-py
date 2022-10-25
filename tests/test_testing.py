import numpy as np
import polars as pl
import pytest
from numpy.random import default_rng
from numpy.testing import assert_allclose, assert_array_equal

from pore_c2.model import EnzymeCutter
from pore_c2.testing import (
    Scenario,
    assign_snps_to_fragments,
    random_concatemer_generator,
    simulate_concatemers,
    simulate_contact_prob_matrix,
    simulate_sequence_with_cut_sites,
)


@pytest.mark.parametrize(
    "pos,buffer,snps,offset,fragments",
    [([1, 5, 19], 0, [1, 5, 19], [1, 5, 9], [1, 1, 2]), ([1, 5, 19], 2, [5], [5], [1])],
)
def test_assign_snps_to_fragments(pos, buffer, snps, offset, fragments):
    fragment_df = pl.DataFrame(
        [
            ("chr1", 0, 10, 1),
            ("chr1", 10, 20, 2),
        ],
        columns=["chrom", "start", "end", "fragment_id"],
    )
    _snps, _offset, _fragments = assign_snps_to_fragments(
        "chr1", pos, fragment_df, buffer
    )
    assert _snps == snps
    assert _offset == offset
    assert _fragments == fragments


def test_simualte_sequence_with_cut_sites():
    cutter = EnzymeCutter.from_name("EcoRI")
    enzyme = cutter.enzyme
    site = enzyme.site

    expected = [10]

    observed, seq = simulate_sequence_with_cut_sites(
        seq_length=20,
        random_state=default_rng(421),
        cut_sites=expected,
    )
    assert expected == observed
    assert site in seq


def test_simulate_concatemer(default_scenario: Scenario):
    rng = default_rng(421)
    concat_gen = random_concatemer_generator(
        fragments_df=default_scenario.fragments_df,
        reference_fasta=default_scenario.ff,
        contact_probs=default_scenario.contact_prob_matrix,
        haplotype_df=default_scenario.haplotype_df,
        random_state=rng,
        max_concatemers=20,
    )
    cutter = EnzymeCutter.from_name(default_scenario.enzyme)
    for (concatemer, expected, _) in concat_gen:
        observed = concatemer.cut(cutter)
        assert [_.coords for _ in observed] == [_.coords for _ in expected]


def test_prob_matrix():
    # chrom, start, end, fragment_id, fragment_idx, ...
    frags = [
        ("chr1", 0, 100),
        ("chr1", 100, 200),
        ("chr2", 0, 50),
        ("chr2", 50, 100),
        ("chr2", 100, 150),
    ]
    frag_df = pl.DataFrame(
        [
            {"chrom": _[0], "start": _[1], "end": _[2], "fragment_idx": x + 1}
            for x, _ in enumerate(frags)
        ],
        orient="row",
    )
    p_cis = 1.0
    m = simulate_contact_prob_matrix(frag_df, p_cis=p_cis)
    expected = [
        [
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
        ],  # only 2 frags on chr1, so 100% prob of contact with other frag
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.5, 0.5],  # prob is split evenly among other two frags on chr2
        [0.0, 0.0, 0.5, 0.0, 0.5],
        [0.0, 0.0, 0.5, 0.5, 0.0],
    ]
    assert_array_equal(m, np.array(expected))
    p_cis = 0.0
    m1 = simulate_contact_prob_matrix(frag_df, p_cis=p_cis)
    expected1 = [
        [
            0.0,
            0.0,
            0.5,
            0.5,
            0.5,
        ],  # each of 3 frags on chr2 have equal chance of contacting
        [0.0, 0.0, 0.5, 0.5, 0.5],  # either frag on chrom 1
        [
            1.0 / 3,
            1.0 / 3,
            0.0,
            0.0,
            0.0,
        ],  # prob is split evenly among other 3 frags on chr2
        [1.0 / 3, 1.0 / 3, 0.0, 0.0, 0.0],
        [1.0 / 3, 1.0 / 3, 0.0, 0.0, 0.0],
    ]
    assert_array_equal(m1, np.array(expected1))


def test_random_walk():
    rng = default_rng(seed=42)
    probs = np.array(
        [
            [
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
            ],  # only 2 frags on chr1, so 100% prob of contact with other frag
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [
                0.0,
                0.0,
                0.0,
                0.5,
                0.5,
            ],  # prob is split evenly among other two frags on chr2
            [0.0, 0.0, 0.5, 0.0, 0.5],
            [0.0, 0.0, 0.5, 0.5, 0.0],
        ]
    )
    c = simulate_concatemers(probs, np.ones(100) * 100, rng)
    counts = np.zeros_like(probs)
    for _ in c:
        if len(_) > 1:
            for left, right in zip(_[:-1], _[1:]):
                counts[left, right] += 1
                counts[right, left] += 1
    counts = counts / counts.sum(axis=1)
    assert_allclose(probs, counts, atol=0.1)


def test_scenario(tmp_path):
    seed = 42
    genome_size: int = 5_000
    num_chroms: int = 2
    cut_rate: float = 0.005
    enzyme: str = "NlaIII"
    num_concatemers: int = 100
    num_haplotypes: int = 0
    variant_density: float = 0.0
    rng = default_rng(seed=int(seed))
    chrom_lengths = {
        f"chr{x+1}": v
        for x, v in enumerate(
            sorted(rng.choice(genome_size, size=num_chroms, replace=False))
        )
    }
    scenario = Scenario(
        rng,
        chrom_lengths,
        cut_rate=cut_rate,
        enzyme=enzyme,
        num_concatemers=num_concatemers,
        num_haplotypes=num_haplotypes,
        variant_density=variant_density,
        temp_path=tmp_path,
    )
    # TODO figure out actual test
    print(scenario.monomer_metadata)
    print(scenario.concatemer_metadata)
