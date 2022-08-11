from pathlib import Path
from typing import Optional

import mappy
from attrs import define
from loguru import logger
from pysam import faidx

from .digest import _get_enzyme, digest_genome
from .settings import MINIMAP2_PRESET
from .utils import FileCollection


@define
class IndexFileCollection(FileCollection):
    fragments: Optional[Path] = Path("{prefix}.fragments.parquet")
    fasta: Optional[Path] = Path("{prefix}.fragments.fasta")
    bed: Optional[Path] = Path("{prefix}.bed.txt")
    mmi: Optional[Path] = Path("{prefix}.mmi")


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
        f"using preset '{MINIMAP2_PRESET}'"
    )
    mappy.Aligner(
        fn_idx_in=str(fasta), fn_idx_out=str(index_files.mmi), preset=MINIMAP2_PRESET
    )

    return index_files
