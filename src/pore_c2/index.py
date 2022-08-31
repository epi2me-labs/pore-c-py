from pathlib import Path
from typing import Any, Dict, List, Optional

from attrs import define
from cattr.preconf.json import make_converter

from .log import get_logger
from .utils import FileCollection

logger = get_logger()


@define(kw_only=True)
class IndexMetadata:
    enzyme: str
    chrom_order: List[str]
    chrom_lengths: Dict[str, int]
    reference_path: str
    pore_c_version: str
    mappy_settings: Dict[Any, Any]


@define
class IndexFileCollection(FileCollection):
    metadata: Path = Path("{prefix}.metadata.json")
    fragments: Path = Path("{prefix}.fragments.parquet")
    mmi: Path = Path("{prefix}.mmi")
    fasta: Optional[Path] = Path("{prefix}.fragments.fasta")
    bed: Optional[Path] = Path("{prefix}.bed.txt")

    def save_metadata(self, md: IndexMetadata):
        converter = make_converter()
        with self.metadata.open("w") as fh:
            fh.write(converter.dumps(md))

    def load_metadata(self) -> IndexMetadata:
        converter = make_converter()
        return converter.loads(self.metadata.read_text(), IndexMetadata)
