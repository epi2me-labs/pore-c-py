"""Settings."""
from pysam import AlignmentHeader

from porecpy import __version__


DEFAULT_ALIGN_HEADER = AlignmentHeader.from_dict({"CO": [""]})
DOWNGRADE_MM = True

EXE_NAME = "porecpy"
VERSION = __version__
