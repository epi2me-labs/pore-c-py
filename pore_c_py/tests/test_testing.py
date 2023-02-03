"""Test testing."""
import numpy as np
from numpy.random import default_rng
from numpy.testing import assert_allclose, assert_array_equal
import polars as pl
import pytest

from pore_c_py.model import EnzymeCutter
from pore_c_py.testing import (
    assign_snps_to_fragments,
    random_concatemer_generator,
    Scenario,
    simulate_concatemers,
    simulate_contact_prob_matrix,
    simulate_sequence_with_cut_sites
)


@pytest.mark.parametrize(
    "pos,buffer,snps,offset,fragments",
    [
        ([1, 5, 19], 0, [1, 5, 19], [1, 5, 9], [1, 1, 2]),
        ([1, 5, 19], 2, [5], [5], [1])
    ],
)
def test_assign_snps_to_fragments(pos, buffer, snps, offset, fragments):
    """Test assign snps to fragments."""
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


def test_simulate_sequence_with_cut_sites():
    """Test simulate sequence with cut sites."""
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
    """Test simulate concatemer."""
    rng = default_rng(421)
    concat_gen = random_concatemer_generator(
        fragments_df=default_scenario.fragments_df,
        reference_fasta=default_scenario.ff,
        contact_probs=default_scenario.contact_prob_matrix,
        fragment_to_snps=default_scenario.fragment_to_snps,
        random_state=rng,
        max_concatemers=20,
    )
    cutter = EnzymeCutter.from_name(default_scenario.enzyme)
    for (concatemer, expected, walk, _) in concat_gen:
        walk_length = len(walk.segments)
        observed = concatemer.cut(cutter)
        assert observed[0].coords.subread_total == walk_length
        assert expected[0].coords.subread_total == walk_length
        assert len(observed) == walk_length
        assert len(expected) == len(observed)
        assert [_.coords for _ in observed] == [_.coords for _ in expected]


def test_prob_matrix():
    """Test probability matrix."""
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
        [0.0, 0.0, 0.0, 0.5, 0.5],  # prob split evenly among other two frags on chr2   # noqa
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
    """Test random walk."""
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
    c = simulate_concatemers(probs, np.ones(100, dtype=int) * 100, rng)
    counts = np.zeros_like(probs)
    for _ in c:
        if len(_) > 1:
            for left, right in zip(_[:-1], _[1:]):
                counts[left, right] += 1
                counts[right, left] += 1
    counts = counts / counts.sum(axis=1)
    assert_allclose(probs, counts, atol=0.1)


def test_scenario(tmp_path):
    """Test scenario."""
    s1 = Scenario(
        genome_size=2_000,
        num_concatemers=10,
        temp_path=tmp_path / "s1",
    )
    s2 = Scenario.from_json(s1.fc.params_json, temp_path=tmp_path / "s2")
    # TODO figure out actual test
    assert len(s1.monomer_metadata) == len(s2.monomer_metadata)
