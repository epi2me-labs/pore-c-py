from contextlib import closing
from itertools import islice
from pathlib import Path
from typing import List, Optional

import typer

from .aligns import annotate_monomer_alignments
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
from .testing import Scenario

# from rich.console import Console


app = typer.Typer(pretty_exceptions_enable=False)
# console = Console()


@app.callback()
def main(quiet: bool = False, logfile: Optional[Path] = None):
    init_logger(quiet=quiet, logfile=logfile)


@app.command()
def align():
    raise NotImplementedError


@app.command()
def merge():
    raise NotImplementedError


@app.command()
def digest(
    file_or_root: Path,
    enzyme: str,
    output_path: Path,
    glob: str = "*.fastq",
    recursive: bool = True,
    max_reads: int = 0,
    remove_tags: Optional[List[str]] = typer.Option(
        None, help="Optionally remove SAM tags from input"
    ),
):

    logger = get_logger()
    logger.info("Digesting concatemers")
    input_files = list(find_files(file_or_root, glob=glob, recursive=recursive))
    header = get_alignment_header(source_files=input_files)
    remove_tags = list(set(remove_tags)) if remove_tags else None
    read_stream = get_concatemer_seqs(input_files, remove_tags=remove_tags)
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


utils = typer.Typer()


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
    summary: bool = False,
    direct_only: bool = False,
):
    logger = get_logger()
    logger.info(f"Processing reads from {bam}")
    input_files = [bam]
    drop_outputs = []
    for flag, filekey in [
        (monomers, "namesorted_bam"),
        (paired_end, "paired_end_bam"),
        (chromunity, "chromunity_parquet"),
        (summary, "summary_json"),
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
        writer.consume(annotated_stream, direct_only=direct_only)

    if summary:
        logger.info(f"Summary information at {output_files.summary_json}")
    # logger.info(
    #    f"Wrote {writer.base_counter:,} bases in "
    #    f"{writer.read_counter:,} reads to {output_path}"
    # )
    return writer


app.add_typer(utils, name="utils")

if __name__ == "__main__":
    app()
