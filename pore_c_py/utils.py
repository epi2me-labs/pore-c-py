"""Log."""
import argparse
from dataclasses import dataclass
import logging
import os
from pathlib import Path
import stat
import sys

import pysam


def log_level():
    """Parser to set logging level and acquire software version/commit."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    modify_log_level = parser.add_mutually_exclusive_group()
    modify_log_level.add_argument(
        '--debug', action='store_const',
        dest='log_level', const=logging.DEBUG, default=logging.INFO,
        help='Verbose logging of debug information.')
    modify_log_level.add_argument(
        '--quiet', action='store_const',
        dest='log_level', const=logging.WARNING, default=logging.INFO,
        help='Minimal logging; warnings only.')
    parser.add_argument("--logfile", type=Path, help="Specify a log file.")

    return parser


def get_named_logger(name):
    """Create a logger with a name.

    :param name: name of logger.
    """
    name = name.ljust(10)[:10]  # so logging is aligned
    logger = logging.getLogger('{}.{}'.format(__package__, name))
    logger.name = name
    return logger


def find_files(root: Path, glob="*.fastq", recursive=True):
    """Find files."""
    if not root.is_dir():
        yield root
    else:
        if recursive and not glob.startswith("**/"):
            glob = f"**/{glob}"
        for f in root.glob(glob):
            yield (f)


def stdout_is_regular_file() -> bool:
    """
    Detect if standard output is a regular file (or say a pipe).

    :return: True if stdout is a regular file, else False.
    """
    mode = os.fstat(sys.stdout.buffer.fileno()).st_mode
    return stat.S_ISREG(mode)


CONCATEMER_ID_TAG = "MI"
MONOMER_DATA_TAG = "Xc"
WALK_TAG = "Xw"

# NOTE: the class below exists only because we don't store individual bits in
#       separate tags; which would be simpler and more direct to access.


@dataclass
class MonomerData:
    """Concatemer Data."""

    concatemer_id: str
    start: int
    end: int
    read_length: int
    subread_idx: int
    subread_total: int

    @classmethod
    def from_pysam(cls, align: pysam.AlignedSegment):
        """Extract concatemer meta information from alignment."""
        cls.concatemer_id = align.get_tag(CONCATEMER_ID_TAG)
        tags = align.get_tag(MONOMER_DATA_TAG)
        cls.start = tags[0]
        cls.end = tags[1]
        cls.read_length = tags[2]
        cls.subread_idx = tags[3]
        cls.subread_total = tags[4]

        if align.is_unmapped:
            cls.name = (
                f"*:"
                f"{cls.start}-{cls.end}")
        else:
            orientation = "+-"[align.is_reverse]
            cls.name = (
                f"{align.reference_name}:{orientation}:"
                f"{align.reference_start}-{align.reference_end}:"
                f"{cls.start}-{cls.end}")
        return cls

    def to_pysam(self, align: pysam.AlignedSegment):
        """Set concatemer metadata on alignment."""
        align.set_tag(CONCATEMER_ID_TAG, self.concatemer_id)
        align.set_tag(MONOMER_DATA_TAG, [
            self.start, self.end, self.read_length,
            self.subread_idx, self.subread_total])
        return align

    @staticmethod
    def set_monomer_data(
            align: pysam.AlignedSegment,
            start, end, read_length, idx, num_intervals):
        """Set the monomer data on an alignment."""
        align.set_tag(
            MONOMER_DATA_TAG,
            [start, end, read_length, idx, num_intervals])

    @staticmethod
    def get_subread_total(align: pysam.AlignedSegment):
        """Return the number of sibling monomers."""
        return align.get_tag(MONOMER_DATA_TAG)[4]
