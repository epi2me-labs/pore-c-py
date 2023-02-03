"""Init."""

__version__ = "1.0.0"

import sys

if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional)
    # when `python_requires = >= 3.8`
    from importlib.metadata import (
        # pyright: reportMissingImports=false
        PackageNotFoundError,
        version,
    )
else:
    from importlib_metadata import PackageNotFoundError, version
    # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = "pore_c_py"
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
