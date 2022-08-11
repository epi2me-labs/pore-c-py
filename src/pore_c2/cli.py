from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from .index import create_index
from .map import map_concatemers

app = typer.Typer()
console = Console()


@app.command()
def index(fasta: Path, enzyme: str, prefix: Optional[Path] = None, force: bool = False):
    create_index(fasta=fasta, enzyme=enzyme, prefix=prefix, force=force)


@app.command()
def map(fastq: Path, enzyme: str, reference_fasta: Path, outfile: Path):
    _ = map_concatemers(enzyme=enzyme, fastq=fastq, reference_fasta=reference_fasta)
    raise NotImplementedError


@app.command()
def merge():
    raise NotImplementedError


if __name__ == "__main__":
    app()
