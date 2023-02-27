"""Test functions from digest."""
import pytest

from pore_c_py import digest


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
