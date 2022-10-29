from pysam import AlignmentHeader

MINIMAP2_SETTINGS = {"preset": "ont"}
DEFAULT_ALIGN_HEADER = AlignmentHeader.from_dict({"CO": [""]})
DOWNGRADE_MM = True
