from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from .index import IndexFileCollection, create_index
from .log import get_logger, init_logger
from .map import map_concatemers

app = typer.Typer()
console = Console()


@app.callback()
def main(quiet: bool = False, logfile: Optional[Path] = None):
    init_logger(quiet=quiet, logfile=logfile)


@app.command()
def index(fasta: Path, enzyme: str, prefix: Optional[Path] = None, force: bool = False):
    logger = get_logger()
    index_files = create_index(fasta=fasta, enzyme=enzyme, prefix=prefix, force=force)
    logger.debug(index_files.metadata.read_text())


@app.command()
def map(fastq: Path, index_metadata: Path, prefix: Optional[Path]):
    index_prefix = Path(str(index_metadata).replace(".metadata.json", ""))
    index_files = IndexFileCollection.with_prefix(index_prefix)
    if not index_files.exists_any():
        raise IOError(f"Couldn't find index files {index_files}")

    index_md = index_files.load_metadata()
    map_concatemers(
        enzyme=index_md.enzyme,
        fastq=fastq,
        mmi=index_files.mmi,
        minimap_settings=index_md.mappy_settings,
        fragment_pq=index_files.fragments,
    )


@app.command()
def merge():
    raise NotImplementedError


if __name__ == "__main__":
    app()
