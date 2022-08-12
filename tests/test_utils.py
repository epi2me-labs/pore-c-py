from pathlib import Path

from attrs import define, fields_dict

from pore_c2.utils import FileCollection, SamFlags


def test_file_collection(tmp_path):
    @define
    class TestFC(FileCollection):
        a: Path = Path("{prefix}.A.txt")
        b: Path = Path("{prefix}.B.txt")

    test_fc = TestFC.with_prefix(Path(tmp_path))

    assert test_fc.exists_all() is False
    assert test_fc.exists_any() is False

    test_fc.a.write_text("a")

    assert test_fc.exists_all() is False
    assert test_fc.exists_any() is True

    test_fc.b.write_text("b")
    assert test_fc.exists_all() is True


def test_sam_flags():

    for key, _ in fields_dict(SamFlags).items():
        for tf in (True, False):
            flags = SamFlags(**{key: tf})
            assert getattr(flags, key) == tf
            integer_val = flags.to_int()
            _ = SamFlags.from_int(integer_val)
            assert _ == flags

    flags = SamFlags.from_int(3844)
    assert flags == SamFlags(
        unmap=True, secondary=True, qcfail=True, dup=True, supplementary=True
    )
