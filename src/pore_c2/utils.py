from pathlib import Path
from typing import Dict, List, Optional, Tuple

from attrs import define, fields


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
