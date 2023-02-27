"""Log."""
import argparse
import logging
import os
from pathlib import Path
import stat
import sys


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
