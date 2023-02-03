"""Log."""
import logging
from pathlib import Path
from typing import Optional


# taken from Remora
class CustomFormatter(logging.Formatter):
    """Custom formatter for python logging module.

    Creates separate logging format for error,
    warning, info and debug logging
    levels.
    """

    err_fmt = "*" * 100 + "\n\tERROR: %(msg)s\n" + "*" * 100
    warn_fmt = "*" * 20 + " WARNING: %(msg)s " + "*" * 20
    info_fmt = "[%(asctime)s %(module)s] %(message)s"
    dbg_fmt = (
        "DBG %(asctime)s : %(msg)s --- %(processName)s-"
        + "%(threadName)s %(module)s.py:%(lineno)d"
    )

    def __init__(self, fmt="[%(asctime)s] %(message)s"):
        """Init."""
        super().__init__(fmt=fmt, datefmt="%H:%M:%S", style="%")

    def format(self, record):
        """Format."""
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = self.dbg_fmt
        elif record.levelno == logging.INFO:
            self._style._fmt = self.info_fmt
        elif record.levelno == logging.WARNING:
            self._style._fmt = self.warn_fmt
        elif record.levelno == logging.ERROR:
            self._style._fmt = self.err_fmt
        result = logging.Formatter.format(self, record)

        self._fmt = format_orig

        return result


CONSOLE = logging.StreamHandler()
CONSOLE.setLevel(logging.INFO)
CONSOLE.setFormatter(CustomFormatter())
ROOT_LOGGER = logging.getLogger("porecpy")
ROOT_LOGGER.setLevel(logging.DEBUG)
ROOT_LOGGER.addHandler(CONSOLE)


def init_logger(quiet: bool = False, logfile: Optional[Path] = None):
    """Initiate logger."""
    if logfile:
        fh = logging.FileHandler(logfile, "w")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(CustomFormatter())
        ROOT_LOGGER.addHandler(fh)
    if quiet:
        CONSOLE.setLevel(logging.WARNING)


def get_logger():
    """Get logger."""
    return logging.getLogger("porecpy")
