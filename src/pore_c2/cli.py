from functools import partial
from pathlib import Path
from typing import Optional

import mappy
import typer
from pysam import FastaFile, faidx  # pyright: ignore [reportGeneralTypeIssues]
from rich.console import Console

from pore_c2 import __version__

from .index import IndexFileCollection, IndexMetadata
from .io import get_alignment_header
from .log import get_logger, init_logger
from .monomers import (
    EnzymeCutter,
    digest_genome,
    digest_read,
    find_files,
    get_reads,
    get_writer,
)
from .settings import MINIMAP2_SETTINGS

app = typer.Typer()
console = Console()


@app.callback()
def main(quiet: bool = False, logfile: Optional[Path] = None):
    init_logger(quiet=quiet, logfile=logfile)


@app.command()
def index(fasta: Path, enzyme: str, prefix: Optional[Path] = None, force: bool = False):
    logger = get_logger()
    try:
        cutter = EnzymeCutter.from_name(enzyme)
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
        cutter=cutter,
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
    logger.debug(index_files.metadata.read_text())
    return index_files


@app.command()
def align():
    raise NotImplementedError


@app.command()
def merge():
    raise NotImplementedError


utils = typer.Typer()


@utils.command()
def digest_concatemers(
    file_or_root: Path,
    enzyme: str,
    output_path: Path,
    glob: str = "*.fastq",
    recursive: bool = True,
):

    logger = get_logger()
    logger.info("Digesting concatemers")
    input_files = list(find_files(file_or_root, glob=glob, recursive=recursive))
    header = get_alignment_header(source_files=input_files)
    read_stream = get_reads(input_files)
    cutter = EnzymeCutter.from_name(enzyme)
    writer = get_writer(output_path, align_header=header)
    cut_read = partial(digest_read, cutter)

    writer.consume(map(cut_read, read_stream))
    logger.info(
        f"Wrote {writer.base_counter:,} bases in "
        f"{writer.read_counter:,} reads to {output_path}"
    )
    return writer


@utils.command()
def process_monomer_alignments(bam: Path, output_path: Path):
    # logger = get_logger()
    # read_stream = get_reads(bam)

    raise NotImplementedError


app.add_typer(utils, name="utils")

if __name__ == "__main__":
    app()
