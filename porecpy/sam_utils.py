"""Sam utils."""
import array
from contextlib import contextmanager
from dataclasses import asdict, dataclass, fields
import enum
from functools import lru_cache
import re
from typing import Any, Literal

import pysam
from pysam import AlignedSegment

TAG_MI_RE = re.compile(r"MI\:Z\:(\S+)")
# XC==concatemer metadata
TAG_XC_RE = re.compile(
    r"Xc:B:i,(?P<start>\d+),(?P<end>\d+),(?P<read_length>\d+),"
    r"(?P<subread_idx>\d+),(?P<subread_total>\d+)"
)
# XW==walk metadata
TAG_XW_RE = re.compile(r"Xw\:Z\:(.+)")  # TODO fill this out


# MM:Z:([ACGTUN][-+]([a-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*
MOD_RE = re.compile(
    r"(?P<canonical>[ACGTUN])"
    "(?P<strand>[-+])"
    "(?P<mod>[a-z]+|[0-9]+)"
    "(?P<skip_scheme>[.?]?)"
    ","
    r"(?P<deltas>[\d,]+)"
)
MOD_TAGS = {"Ml", "ML", "Mm", "MM"}
MM_TAG_SKIP_SCHEME_RE = re.compile(r"[.?]")

MOLECULE_TAG = "MI"
WALK_TAG = "Xw"
CONCATEMER_TAG = "Xc"
#  <alpha-num><alpha-num>:<type>:<data>
FASTQ_TAG_RE = re.compile(r"(?P<tag>\w\w):(?P<type>[ABfHiZ]):(?P<data>\S+)")


WALK_SEGMENT_RE = re.compile(
    r"(?P<chrom>(\S+?)):(?P<orientation>[+-]):"
    r"(?P<genome_start>\d+)-(?P<genome_end>\d+):"
    r"(?P<read_start>\d+)-(?P<read_end>\d+)"
)

# copied from pysam.libcalignedsegment.pyx
SAM_TYPES = "iiiiiif"
HTSLIB_TYPES = "cCsSiIf"
PARRAY_TYPES = "bBhHiIf"
ARRAY_TO_HTSLIB_TRANS = str.maketrans(PARRAY_TYPES, HTSLIB_TYPES)
PYSAM_TO_SAM_TRANS = str.maketrans(HTSLIB_TYPES, SAM_TYPES)


@contextmanager
def pysam_verbosity(level: int = 0):
    """Pysam verbosity."""
    current = pysam.set_verbosity(level)
    yield
    pysam.set_verbosity(current)


def tag_tuple_to_str(key: str, val: Any, value_type: str):
    """Tag tuple to string."""
    # TODO not sure if this works for H
    value_type = value_type.translate(PYSAM_TO_SAM_TRANS)
    if value_type == "B":
        assert isinstance(val, array.array)
        elemtype = val.typecode.translate(ARRAY_TO_HTSLIB_TRANS)
        val = ",".join(map(str, val))
        res = f"{key}:{value_type}:{elemtype},{val}"
    elif value_type in "AfiZ":
        res = f"{key}:{value_type}:{val}"
    else:
        # hope that python string formatting is correct
        # TODO: warn about this?
        res = f"{key}:{value_type}:{val}"
    if ":S:" in res:
        raise ValueError(key, val, value_type)
    return res


def downgrade_mm_tag(align: AlignedSegment) -> AlignedSegment:
    """Downgraed mm tag."""
    tag, orig_tag = None, None
    for tag in ["MM", "Mm"]:
        try:
            orig_tag = align.get_tag(tag)
        except KeyError:
            continue
        if orig_tag:
            break
    # alignment doesn't have tag, return
    if not orig_tag:
        return align
    assert isinstance(orig_tag, str)
    m = MM_TAG_SKIP_SCHEME_RE.search(orig_tag)
    # tag does't contain skip_scheme character
    if not m or not tag:
        return align
    new_tag = f"{tag}:Z:" + MM_TAG_SKIP_SCHEME_RE.sub("", orig_tag)
    d = align.to_dict()
    d["tags"] = [t for t in d["tags"] if not t.startswith(tag)] + [new_tag]
    return AlignedSegment.from_dict(d, header=align.header)


@enum.unique
class SamEnum(enum.IntFlag):
    """Sam enum."""

    paired = 1  # template having multiple segments in sequencing
    proper_pair = 2  # each segment properly aligned according to the aligner
    unmap = 4  # segment unmapped
    munmap = 8  # next segment in the template unmapped
    reverse = 16  # SEQ being reverse complemented
    mreverse = 32  # SEQ of the next segment in the template being reverse
    # complemented
    read1 = 64  # the first segment in the template
    read2 = 128  # the last segment in the template
    secondary = 256  # secondary alignment
    qcfail = 512  # not passing filters,
    # such as platform/vendor quality controls
    dup = 1024  # PCR or optical duplicate
    supplementary = 2048  # supplementary alignment


class AlignCategory(enum.IntEnum):
    """Align category."""

    primary = 0
    unmapped = 1
    supplementary = 2
    secondary = 3


@dataclass()
class SamFlags:
    """Sam flags."""

    paired: bool = False
    proper_pair: bool = False
    unmap: bool = False
    munmap: bool = False
    reverse: bool = False
    mreverse: bool = False
    read1: bool = False
    read2: bool = False
    secondary: bool = False
    qcfail: bool = False
    dup: bool = False
    supplementary: bool = False

    def to_int(self):
        """To int."""
        res = 0
        for key, val in asdict(self).items():
            if val is True:
                res = res | SamEnum[key].value
        return res

    def copy(self):
        """Copy."""
        settings = asdict(self)
        return SamFlags(**settings)

    @classmethod
    def from_int(cls, val: int):
        """From int."""
        kwds = {}
        for f in fields(cls):
            kwds[f.name] = (val & SamEnum[f.name].value) > 0
        return cls(**kwds)

    @property
    def primary(self):
        """Primary."""
        return not (self.secondary | self.supplementary)

    @property
    def category(self) -> AlignCategory:
        """Category."""
        if self.secondary:
            return AlignCategory.secondary
        elif self.supplementary:
            return AlignCategory.supplementary
        elif self.unmap:
            return AlignCategory.unmapped
        else:
            return AlignCategory.primary

    @property
    def strand(self) -> Literal["+", "-", "."]:
        """Strand."""
        if self.unmap:
            return "."
        elif self.reverse:
            return "-"
        else:
            return "+"

    @staticmethod
    @lru_cache
    def int_to_strand(flag: int) -> Literal["+", "-", "."]:
        """Int to strand."""
        return SamFlags.from_int(flag).strand

    @staticmethod
    @lru_cache
    def int_to_category(
        flag: int,
    ) -> AlignCategory:
        """Int to category."""
        return SamFlags.from_int(flag).category
