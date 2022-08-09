from pathlib import Path

import typer
from loguru import logger
from pysam import faidx
from rich.console import Console

from pore_c2.digest import digest_genome

app = typer.Typer()
console = Console()


@app.command()
def index(fasta: Path, enzyme: str, outfile: Path):
    idx_path = Path(str(fasta) + ".fai")
    if not idx_path.exists():
        faidx(str(fasta))
    df = digest_genome(enzyme, fasta, outfile)
    df.write_parquet(outfile)
    logger.debug(df)
    logger.info(f"Wrote {len(df)} fragments to {outfile}")


@app.command()
def map():
    raise NotImplementedError


@app.command()
def merge():
    raise NotImplementedError


if __name__ == "__main__":
    app()
