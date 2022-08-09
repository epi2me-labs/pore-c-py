from pathlib import Path

import typer
from rich.console import Console

app = typer.Typer()
console = Console()


@app.command()
def index(reference_fasta: Path, enzyme: str):
    raise NotImplementedError


@app.command()
def map():
    raise NotImplementedError


@app.command()
def merge():
    raise NotImplementedError


if __name__ == "__main__":
    app()
