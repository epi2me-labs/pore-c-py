import re

from pysam import AlignmentHeader

MINIMAP2_SETTINGS = {"preset": "ont"}
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

WALK_TAG = "Xw"
MOLECULE_TAG = "MI"
TAG_XW_RE = re.compile(WALK_TAG + r"\:Z\:(.+)")  # TODO fill this out

CONCATEMER_TAG = "Xc"

#  <alpha-num><alpha-num>:<type>:<data>
FASTQ_TAG_RE = re.compile(r"(?P<tag>\w\w):(?P<type>[ABfHiZ]):(?P<data>\S+)")
DEFAULT_ALIGN_HEADER = AlignmentHeader.from_dict({"CO": [""]})
DOWNGRADE_MM = True

TAG_MI_RE = re.compile(MOLECULE_TAG + r"\:Z\:(\S+)")
# XC==concatemer metadata
TAG_XC_RE = re.compile(
    CONCATEMER_TAG
    + r":B:i,(?P<start>\d+),(?P<end>\d+),(?P<subread_idx>\d+),(?P<subread_total>\d+)"
)
