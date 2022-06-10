from pathlib import Path

import typer
from rich.console import Console

app = typer.Typer()
console = Console()


@app.command()
def digest(reference_fasta: Path):
    raise NotImplementedError


@app.command()
def map():
    raise NotImplementedError


@app.command()
def merge():
    raise NotImplementedError


@app.command()
def haplotag():
    raise NotImplementedError


@app.command()
def to_cooler():
    raise NotImplementedError


if __name__ == "__main__":
    app()
