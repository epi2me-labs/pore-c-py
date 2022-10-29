import array
import re
from typing import Any

from pysam import AlignedSegment

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

WALK_TAG = "Xw"
CONCATEMER_TAG = "Xc"
#  <alpha-num><alpha-num>:<type>:<data>
FASTQ_TAG_RE = re.compile(r"(?P<tag>\w\w):(?P<type>[ABfHiZ]):(?P<data>\S+)")

# copied from pysam.libcalignedsegment.pyx
SAM_TYPES = "iiiiiif"
HTSLIB_TYPES = "cCsSiIf"
PARRAY_TYPES = "bBhHiIf"


ARRAY_TO_HTSLIB_TRANS = str.maketrans(PARRAY_TYPES, HTSLIB_TYPES)
PYSAM_TO_SAM_TRANS = str.maketrans(HTSLIB_TYPES, SAM_TYPES)


TAG_MI_RE = re.compile(r"MI\:Z\:(\S+)")
# XC==concatemer metadata
TAG_XC_RE = re.compile(
    r"Xc:B:i,(?P<start>\d+),(?P<end>\d+),(?P<subread_idx>\d+),(?P<subread_total>\d+)"
)
# XW==walk metadata
TAG_XW_RE = re.compile(r"Xw\:Z\:(.+)")  # TODO fill this out

WALK_SEGMENT_RE = re.compile(
    r"(?P<chrom>(\S+?)):(?P<orientation>[+-]):"
    r"(?P<genome_start>\d+)-(?P<genome_end>\d+):"
    r"(?P<read_start>\d+)-(?P<read_end>\d+)"
)


def tag_tuple_to_str(key: str, val: Any, value_type: str):
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
