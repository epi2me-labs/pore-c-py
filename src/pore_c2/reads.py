from pathlib import Path
from typing import Any, Dict, Iterator

from attrs import Factory, define
from pysam import FastqFile


@define(frozen=True)
class Read:
    name: str
    sequence: str
    quality: str
    tags: Dict[str, Any] = Factory(dict)


def get_reads(path: Path) -> Iterator[Read]:
    # TODO: handle glob
    # TOOD  handle ubam
    # TODO  model tags
    if path.is_dir():
        raise NotImplementedError
    else:
        fq = FastqFile(str(path))
        for rec in fq:
            yield Read(rec.name, rec.sequence, rec.quality)
