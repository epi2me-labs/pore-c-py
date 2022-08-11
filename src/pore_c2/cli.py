from pathlib import Path
from typing import Optional

import typer
from loguru import logger
from rich.console import Console

from .index import IndexFileCollection, create_index
from .map import map_concatemers

app = typer.Typer()
console = Console()


@app.command()
def index(fasta: Path, enzyme: str, prefix: Optional[Path] = None, force: bool = False):
    index_files = create_index(fasta=fasta, enzyme=enzyme, prefix=prefix, force=force)
    logger.debug(index_files.metadata.read_text())


@app.command()
def map(fastq: Path, index_metadata: Path, prefix: Optional[Path]):

    index_prefix = Path(str(index_metadata).replace(".metadata.json", ""))
    index_files = IndexFileCollection.with_prefix(index_prefix)
    if not index_files.exists_any():
        raise IOError(f"Couldn't find index files {index_files}")

    md = index_files.load_metadata()
    map_concatemers(
        enzyme=md.enzyme,
        fastq=fastq,
        mmi=index_files.mmi,
        minimap_settings=md.mappy_settings,
    )


@app.command()
def merge():
    raise NotImplementedError


if __name__ == "__main__":
    app()
