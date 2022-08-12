import logging
from pathlib import Path
from typing import Optional

CONSOLE = logging.StreamHandler()
CONSOLE.setLevel(logging.INFO)
ROOT_LOGGER = logging.getLogger("pore-c")
ROOT_LOGGER.setLevel(logging.DEBUG)
ROOT_LOGGER.addHandler(CONSOLE)


def init_logger(quiet: bool = False, logfile: Optional[Path] = None):

    if logfile:
        fh = logging.FileHandler(logfile, "w")
        fh.setLevel(logging.DEBUG)
        ROOT_LOGGER.addHandler(fh)
    if quiet:
        CONSOLE.setLevel(logging.WARNING)


def get_logger():
    return logging.getLogger("pore-c")
