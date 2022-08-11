from collections import namedtuple
from pathlib import Path
from typing import Dict, NamedTuple

from attrs import define, field


@define
class PrefixedFileCollection:
    prefix: Path
    suffixes: Dict[str, str]
    p: NamedTuple = field(init=False)

    def __attrs_post_init__(self):
        nt = namedtuple("PrefixedPaths", list(self.suffixes.keys()))
        p = nt(
            **{key: self.prefix.with_suffix(val) for key, val in self.suffixes.items()}
        )
        self.p = p

    def __iter__(self):
        for key in self.suffixes:
            yield getattr(self.p, key)

    def exists_all(self):
        for path in self.p:
            if path is None:
                continue
            if not path.exists():
                return False
        return True

    def exists_any(self):
        for path in self.p:
            if path is None:
                continue
            if path.exists():
                return True
        return False
