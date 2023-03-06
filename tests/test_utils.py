"""Test functions in utils."""
from pathlib import Path

import pytest

from pore_c_py import utils


@pytest.mark.parametrize(
    "layout,glob,recursive,expected",
    [
        (["A/a/a.fastq"], None, None, None),
        (["a.fastq"], None, None, None),
        (["a.fastq.gz"], None, None, []),
        (["a.fastq", "A/a/a.fastq"], None, False, ["a.fastq"]),
    ],
)
def test_find_files(tmp_path, layout, glob, recursive, expected):
    """Test find files."""
    root = Path(str(tmp_path))
    if expected is None:
        expected = layout
    kwds = {}
    if glob is not None:
        kwds["glob"] = glob
    if recursive is not None:
        kwds["recursive"] = recursive
    for _ in layout:
        p = root / _
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(_)
    res = list([
        str(x).replace(str(tmp_path), "")[1:]
        for x in utils.find_files(root, **kwds)]
    )
    assert res == expected
