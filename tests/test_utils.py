from pore_c2.utils import PrefixedFileCollection


def test_file_collection(tmp_path):
    test_fc = PrefixedFileCollection(tmp_path, {"a": ".A.txt", "b": ".B.txt"})

    assert test_fc.exists_all() is False
    assert test_fc.exists_any() is False

    test_fc.p.a.write_text("a")

    assert test_fc.exists_all() is False
    assert test_fc.exists_any() is True

    test_fc.p.b.write_text("b")
    assert test_fc.exists_all() is True
