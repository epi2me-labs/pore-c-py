"""Settings."""
from pysam import AlignmentHeader

from pore_c_py import __version__


DEFAULT_ALIGN_HEADER = AlignmentHeader.from_dict({"CO": [""]})
DOWNGRADE_MM = True

EXE_NAME = "porecpy"
VERSION = __version__
