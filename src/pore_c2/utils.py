import enum
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam
from attrs import asdict, define, fields, fields_dict


@contextmanager
def pysam_verbosity(level: int = 0):
    current = pysam.set_verbosity(level)
    yield
    pysam.set_verbosity(current)


@define
class FileCollection:
    _path_attrs: List[str]

    @classmethod
    def with_prefix(cls, prefix: Path, drop: Optional[List[str]] = None):
        path_attrs = []
        kwds = {}
        for f in fields(cls):  # pyright: ignore
            if f.name.startswith("_"):
                continue
            if drop and f.name in drop:
                kwds[f.name] = None
            else:
                kwds[f.name] = Path(str(f.default).format(prefix=str(prefix)))
            path_attrs.append(f.name)

        return cls(
            path_attrs=path_attrs, **kwds  # pyright: ignore [reportGeneralTypeIssues]
        )

    def __iter__(self):
        for a in self._path_attrs:
            yield getattr(self, a)

    def items(self) -> List[Tuple[str, Optional[Path]]]:
        return [(a, getattr(self, a)) for a in self._path_attrs]

    def existing(self) -> Dict[str, Path]:
        return {
            key: val for key, val in self.items() if val is not None and val.exists()
        }

    def exists_any(self) -> bool:
        return len(self.existing()) > 0

    def exists_all(self) -> bool:
        for p in self:
            if p is not None and not p.exists():
                return False
        return True


@enum.unique
class SamEnum(enum.IntFlag):
    paired = 1
    proper_pair = 2
    unmap = 4
    munmap = 8
    reverse = 16
    mreverse = 32
    read1 = 64
    read2 = 128
    secondary = 256
    qcfail = 512
    dup = 1024
    supplementary = 2048


@define(kw_only=True)
class SamFlags:
    paired: bool = False
    proper_pair: bool = False
    unmap: bool = False
    munmap: bool = False
    reverse: bool = False
    mreverse: bool = False
    read1: bool = False
    read2: bool = False
    secondary: bool = False
    qcfail: bool = False
    dup: bool = False
    supplementary: bool = False

    def to_int(self):
        res = 0
        for key, val in asdict(self).items():
            if val is True:
                res = res | SamEnum[key].value
        return res

    @classmethod
    def from_int(cls, val: int):
        kwds = {}
        for key, _ in fields_dict(cls).items():  # pyright: ignore
            kwds[key] = (val & SamEnum[key].value) > 0
        return cls(**kwds)

    @property
    def primary(self):
        return not (self.secondary | self.supplementary)

    @property
    def category(self):
        if self.secondary:
            return "secondary"
        elif self.supplementary:
            return "supplementary"
        elif self.unmap:
            return "unmapped"
        else:
            return "primary"
