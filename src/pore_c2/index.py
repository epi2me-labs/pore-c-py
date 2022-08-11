from pathlib import Path
from typing import Any, Dict, List, Optional

import mappy
from attrs import define
from cattr.preconf.json import make_converter
from loguru import logger
from pysam import FastaFile, faidx  # pyright: ignore

from pore_c2 import __version__

from .digest import _get_enzyme, digest_genome
from .settings import MINIMAP2_SETTINGS
from .utils import FileCollection


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


def create_index(
    *, fasta: Path, enzyme: str, prefix: Optional[Path] = None, force: bool = False
) -> IndexFileCollection:
    try:
        _ = _get_enzyme(enzyme)
    except Exception:
        logger.error(f"Error loading enzyme {enzyme}", exc_info=True)
        raise
    if prefix is None:
        prefix = fasta.parent / f"{fasta.stem}.porec.{enzyme}"
    index_files = IndexFileCollection.with_prefix(prefix)
    if index_files.exists_any() and not force:
        logger.error(
            "Some of the outputs already exist, please remove before continuing"
        )
        raise IOError
    idx_path = Path(str(fasta) + ".fai")
    if not idx_path.exists():
        logger.info(f"Creating a .fai for {fasta}")
        faidx(str(fasta))
    df = digest_genome(
        enzyme_id=enzyme,
        fasta=fasta,
        bed_file=index_files.bed,
        fasta_out=index_files.fasta,
    )
    if index_files.fragments:
        df.write_parquet(index_files.fragments)
        logger.info(f"Wrote {len(df)} fragments to {index_files.fragments}")
    logger.debug(
        f"Creating minimap index of {fasta} at {index_files.mmi} "
        f"using preset '{MINIMAP2_SETTINGS}'"
    )
    mappy.Aligner(
        fn_idx_in=str(fasta), fn_idx_out=str(index_files.mmi), **MINIMAP2_SETTINGS
    )
    ff = FastaFile(str(fasta))
    metadata = IndexMetadata(
        enzyme=enzyme,
        reference_path=str(fasta.absolute()),
        chrom_order=list(ff.references),
        chrom_lengths={c: ff.get_reference_length(c) for c in ff.references},
        pore_c_version=__version__,
        mappy_settings=MINIMAP2_SETTINGS,
    )
    index_files.save_metadata(metadata)
    return index_files
