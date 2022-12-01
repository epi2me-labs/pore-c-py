from pysam import AlignmentHeader

from pore_c2 import __version__

DEFAULT_ALIGN_HEADER = AlignmentHeader.from_dict({"CO": [""]})
DOWNGRADE_MM = True

EXE_NAME = "pore-c2"
VERSION = __version__
