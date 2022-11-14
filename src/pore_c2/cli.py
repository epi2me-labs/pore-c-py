from contextlib import closing
from itertools import islice
from pathlib import Path
from typing import Optional

import mappy
import typer
from pysam import FastaFile, faidx  # pyright: ignore [reportGeneralTypeIssues]

from pore_c2 import __version__

from .aligns import annotate_monomer_alignments
from .index import IndexFileCollection, IndexMetadata
from .io import (
    AnnotatedMonomerFC,
    AnnotatedMonomerWriter,
    find_files,
    get_alignment_header,
    get_concatemer_seqs,
    get_monomer_aligns,
    get_monomer_writer,
)
from .log import get_logger, init_logger
from .model import EnzymeCutter
from .monomers import digest_genome
from .settings import MINIMAP2_SETTINGS
from .testing import Scenario

# from rich.console import Console


app = typer.Typer(pretty_exceptions_enable=False)
# console = Console()


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
    max_reads: int = 0,
):

    logger = get_logger()
    logger.info("Digesting concatemers")
    input_files = list(find_files(file_or_root, glob=glob, recursive=recursive))
    header = get_alignment_header(source_files=input_files)
    read_stream = get_concatemer_seqs(input_files)
    if max_reads:
        read_stream = islice(read_stream, max_reads)
    cutter = EnzymeCutter.from_name(enzyme)
    monomer_stream = (
        monomer.read_seq for read in read_stream for monomer in read.cut(cutter)
    )

    with closing(get_monomer_writer(output_path, header=header)) as writer:
        writer.consume(monomer_stream)
    # logger.info(
    #    f"Wrote {writer.base_counter:,} bases in "
    #    f"{writer.read_counter:,} reads to {output_path}"
    # )
    return writer


@utils.command()
def create_test_data(
    base_dir: Path,
    genome_size: int = 5_000,
    num_chroms: int = 2,
    cut_rate: float = 0.005,
    enzyme: str = "NlaIII",
    num_concatemers: int = 100,
    num_haplotypes: int = 0,
    variant_density: float = 0.05,
    p_cis: float = 0.8,
    mean_frags_per_concatemer: int = 5,
    max_frags_per_concatemer: int = 10,
    seed: int = 42,
    file_prefix: str = "porec_test",
    create_truth_files: bool = False,
):
    logger = get_logger()
    logger.info(f"Creating test data at: {base_dir}")
    temp_path = base_dir / file_prefix
    scenario = Scenario(
        seed=seed,
        genome_size=genome_size,
        num_chroms=num_chroms,
        cut_rate=cut_rate,
        enzyme=enzyme,
        num_concatemers=num_concatemers,
        num_haplotypes=num_haplotypes,
        variant_density=variant_density,
        temp_path=temp_path,
        p_cis=p_cis,
        mean_frags_per_concatemer=mean_frags_per_concatemer,
        max_frags_per_concatemer=max_frags_per_concatemer,
    )
    logger.info(f"Creating scenario: {scenario}")

    logger.info(f"Genome fasta: {scenario.reference_fasta}")
    logger.info(f"Concatemer fastq: {scenario.concatemer_fastq}")
    logger.info(f"Concatemer ubam: {scenario.concatemer_ubam}")
    if num_haplotypes >= 2 and variant_density > 0:
        logger.info(f"Phased VCF: {scenario.phased_vcf}")
    if create_truth_files:
        for file_id, path in scenario.fc.truth_files():
            logger.info(f"Truth data {file_id}: {path}")
    return scenario


@utils.command()
def process_monomer_alignments(
    bam: Path,
    output_prefix: Path,
    force: bool = False,
    monomers: bool = True,
    paired_end: bool = False,
    chromunity: bool = False,
):
    logger = get_logger()
    logger.info(f"Processing reads from {bam}")
    input_files = [bam]
    drop_outputs = []
    for flag, filekey in [
        (monomers, "namesorted_bam"),
        (paired_end, "paired_end_bam"),
        (chromunity, "chromunity_parquet"),
    ]:
        if not flag:
            drop_outputs.append(filekey)
    output_files = AnnotatedMonomerFC.with_prefix(output_prefix, drop=drop_outputs)
    if output_files.exists_any() and not force:
        logger.error(
            "Some of the outputs already exist, please remove before continuing"
        )
        raise IOError
    header = get_alignment_header(source_files=input_files)
    monomer_aligns = get_monomer_aligns(input_files)
    annotated_stream = annotate_monomer_alignments(monomer_aligns)
    with closing(
        AnnotatedMonomerWriter.from_file_collection(output_files, header=header)
    ) as writer:
        writer.consume(annotated_stream)
    # logger.info(
    #    f"Wrote {writer.base_counter:,} bases in "
    #    f"{writer.read_counter:,} reads to {output_path}"
    # )
    return writer


app.add_typer(utils, name="utils")

if __name__ == "__main__":
    app()
