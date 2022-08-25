from itertools import islice
from pathlib import Path
from typing import Optional

import mappy as mp
import typer
from pysam import AlignmentHeader
from rich.console import Console

from .digest import EnzymeCutter, digest_fastq
from .index import IndexFileCollection, create_index
from .io import MappingFileCollection
from .log import get_logger, init_logger
from .multiprocessing import MappyThreadPool
from .reads import get_reads

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
def map(
    fastq: Path, index_metadata: Path, prefix: Path, max_reads: Optional[int] = None
):
    # load index
    index_prefix = Path(str(index_metadata).replace(".metadata.json", ""))
    index_files = IndexFileCollection.with_prefix(index_prefix)
    if not index_files.exists_any():
        raise IOError(f"Couldn't find index files {index_files}")
    index_md = index_files.load_metadata()
    cutter = EnzymeCutter.from_name(index_md.enzyme_id)

    # output files
    map_fc = MappingFileCollection.with_prefix(prefix)
    header = AlignmentHeader.from_references(
        index_md.chrom_order,
        [index_md.chrom_lengths[c] for c in index_md.chrom_order],
    )

    # create a read stream
    read_stream = get_reads(fastq)
    if max_reads is not None:
        assert max_reads > 0
        read_stream = islice(read_stream, 0, max_reads)

    # initialise the Mappy thread pool
    # overlapper = FragmentOverlapper.from_parquet(index_files.fragments)
    aligner = mp.Aligner(fn_idx_in=str(index_files.mmi))
    mapping_stream = MappyThreadPool(
        read_iter=read_stream, aligner=aligner, cutter=cutter, n_threads=1
    )

    # writer = MapWriter(
    #    mapping_stream,
    #    fc=map_fc,
    #    sam_header=header,jj
    #    reference_filename=Path(index_md.reference_path),
    # )

    # map_concatemers(
    #    enzyme=index_md.enzyme,
    #    fastq=fastq,
    #    mmi=index_files.mmi,
    #    minimap_settings=index_md.mappy_settings,
    #    fragment_pq=index_files.fragments,
    #    writer=writer,
    # )


@app.command()
def merge():
    raise NotImplementedError


utils = typer.Typer()


@utils.command()
def digest_concatemers(fastq: Path, enzyme: str, output_path: Path):
    logger = get_logger()
    read_stream = get_reads(fastq)
    cutter = EnzymeCutter.from_name(enzyme)
    digest_fastq(
        cutter=cutter,
        read_iter=read_stream,
        fastq_out=output_path,
        return_dataframe=False,
    )
    logger.info(f"Split reads written to {output_path}")


app.add_typer(utils, name="utils")

if __name__ == "__main__":
    app()
