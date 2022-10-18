from enum import IntFlag
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

    for key, _ in fields_dict(SamFlags).items():  # pyright: ignore
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


# TODO: move into main file
class SamBits(IntFlag):
    paired = 1  # template having multiple segments in sequencing
    proper_pair = 2  # each segment properly aligned according to the aligner
    unmap = 4  # segment unmapped
    munmap = 8  # next segment in the template unmapped
    reverse = 16  # SEQ being reverse complemented
    mreverse = 32  # SEQ of the next segment in the template being reverse complemented
    read1 = 64  # the first segment in the template
    read2 = 128  # the last segment in the template
    secondary = 256  # secondary alignment
    qcfail = 512  # not passing filters, such as platform/vendor quality controls
    dup = 1024  # PCR or optical duplicate
    supplementary = 2048  # supplementary alignment


def is_primary(flag: int):
    return bool(flag & ~(SamBits.supplementary | SamBits.secondary | SamBits.unmap))


def test_mask():
    f = SamFlags(unmap=True, secondary=True).to_int()
    assert bool(f & SamBits.unmap)
    # assert(bool(SamEnum['unmap'].value & f.to_int()) is True)


def test_align_categories():
    flags = [
        SamFlags(unmap=True),
        SamFlags(secondary=True),
        SamFlags(supplementary=True),
        SamFlags(),
    ]
    categories = [_.category for _ in flags]
    assert [c.name for c in categories] == [
        "unmapped",
        "secondary",
        "supplementary",
        "primary",
    ]
    assert [c.name for c in sorted(categories)] == [
        "primary",
        "unmapped",
        "supplementary",
        "secondary",
    ]


def test_strand():
    flags = [
        SamFlags(),
        SamFlags(reverse=True),
        SamFlags(unmap=True),
        SamFlags(secondary=True),
    ]
    strands = [SamFlags.int_to_strand(_.to_int()) for _ in flags]
    assert strands == ["+", "-", ".", "+"]
