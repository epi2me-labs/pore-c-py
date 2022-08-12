from pore_c2.overlaps import FragmentOverlapper
from pore_c2.testing import Scenario


def test_load(scenario: Scenario):
    _ = FragmentOverlapper.from_dataframe(scenario.fragments_df)


def test_overlapper():
    overlapper = FragmentOverlapper(
        left={"chr1": [0, 10, 20]}, ids={"chr1": ["f01", "f02", "f03"]}
    )
    snaps = [overlapper.snap("chr1", x) for x in range(20)]
    assert (snaps[0]) == (
        0,
        0,
        0.0,
    )  # first fragment, 0 bases from junction, 0% of fragment length
    assert (snaps[4]) == (
        0,
        -4,
        40.0,
    )  # first fragment, 4 bases to the right of junction, 40% of fragment length
    assert (snaps[5]) == (
        0,
        -5,
        50.0,
    )  # half-way across, snaps left
    assert (snaps[6]) == (
        1,
        4,
        40.0,
    )  # second junction, 4 base to the left of it, 40% of the fragment
    # print(overlapper.overlaps("chr1", 0, 19))
