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
CONCATEMER_TAG = "Xc"
#  <alpha-num><alpha-num>:<type>:<data>
FASTQ_TAG_RE = re.compile(r"(?P<tag>\w\w):(?P<type>[ABfHiZ]):(?P<data>\S+)")
DEFAULT_ALIGN_HEADER = AlignmentHeader.from_dict({"CO": [""]})
DOWNGRADE_MM = True
