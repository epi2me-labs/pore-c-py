from pathlib import Path
from typing import Optional

import typer
from loguru import logger
from pysam import faidx
from rich.console import Console

from .digest import _get_enzyme, digest_genome
from .map import map_concatemers
from .utils import PrefixedFileCollection

app = typer.Typer()
console = Console()


@app.command()
def index(fasta: Path, enzyme: str, prefix: Optional[Path] = None, force: bool = False):

    try:
        _ = _get_enzyme(enzyme)
    except Exception:
        logger.error(f"Error loading enzyme {enzyme}", exc_info=True)
        raise
    if prefix is None:
        prefix = fasta.parent / f"{fasta.stem}.porec.{enzyme}"
    else:
        if enzyme not in prefix.name:
            prefix = Path(str(prefix) + f".{enzyme}")

    index_files = PrefixedFileCollection(
        prefix,
        {
            "mmi": ".mmi",
            "fragments": ".digest.parquet",
            "bed": ".bed",
            "fasta": ".fasta",
        },
    )
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
        bed_file=index_files.p.bed,
        fasta_out=index_files.p.fasta,
    )
    df.write_parquet(index_files.p.fragments)
    logger.debug(df)
    logger.info(f"Wrote {len(df)} fragments to {index_files.p.fragments}")


@app.command()
def map(fastq: Path, enzyme: str, reference_fasta: Path, outfile: Path):
    _ = map_concatemers(enzyme=enzyme, fastq=fastq, reference_fasta=reference_fasta)
    raise NotImplementedError


@app.command()
def merge():
    raise NotImplementedError


if __name__ == "__main__":
    app()
