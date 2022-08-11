from pathlib import Path
from typing import Optional

from attrs import define

from pore_c2.utils import FileCollection


def test_file_collection(tmp_path):
    @define
    class TestFC(FileCollection):
        a: Optional[Path] = Path("{prefix}.A.txt")
        b: Optional[Path] = Path("{prefix}.B.txt")

    test_fc = TestFC.with_prefix(Path(tmp_path))

    assert test_fc.exists_all() is False
    assert test_fc.exists_any() is False

    test_fc.a.write_text("a")

    assert test_fc.exists_all() is False
    assert test_fc.exists_any() is True

    test_fc.b.write_text("b")
    assert test_fc.exists_all() is True
