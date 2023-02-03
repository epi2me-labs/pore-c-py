"""Cli."""
from contextlib import closing
from itertools import islice
from pathlib import Path
from typing import List, Optional

import typer

from pore_c_py.aligns import (
    annotate_monomer_alignments,
    group_aligns_by_concatemers
)
from pore_c_py.io import (
    AnnotatedMonomerFC,
    AnnotatedMonomerWriter,
    create_chunked_bam,
    find_files,
    get_alignment_header,
    get_concatemer_seqs,
    get_monomer_aligns,
    get_monomer_writer,
)
from pore_c_py.log import get_logger, init_logger
from pore_c_py.model import EnzymeCutter
from pore_c_py.testing import Scenario

# from rich.console import Console


app = typer.Typer(pretty_exceptions_enable=False)
# console = Console()


@app.callback()
def main(quiet: bool = False, logfile: Optional[Path] = None):
    """Entry point."""
    init_logger(quiet=quiet, logfile=logfile)


@app.command()
def digest(
    input: Path = typer.Argument(
        ..., help="An unaligned BAM file or FASTQ file or a directory of same"
    ),
    enzyme: str = typer.Argument(
        ..., help="A restriction enzyme name eg. NlaIII"),
    output_bam: Path = typer.Argument(
        ...,
        help="An unaligned BAM file with a separate record for each monomer"
    ),
    glob: str = typer.Option(
        default="*.fastq",
        help="If INPUT is a directory use this glob to search for files",
    ),
    recursive: bool = typer.Option(
        default=True,
        help="If INPUT is a directory search recusively for files matching 'glob'", # noqa
    ),
    max_reads: int = typer.Option(
        default=0, help="Take the first n reads (useful for testing)"
    ),
    remove_tags: Optional[List[str]] = typer.Option(
        None, help="Optionally remove SAM tags from input"
    ),
):
    """
    Digest concatemer sequences into monomers using a restriction enzyme.

    The resulting BAM file has a record for each of the monomer subreads. The
    BAM file has the following properties:

    1. Read ID format: <original_read_id>:<subread_start>:<subread_end>

        The subread start and end are left-padded with zeros so that \
        when lexographically sorted the monomers are in the same order as \
        they originally appeared in the concatemer.

    2. Molecule Identifier:

        The MI SAM tag is set to the originial read/concatemer id.

    3. Concatemer Coordinates (Xc):

        The Xc tag contains the relative coordinates of the monomer on the
        concatemer. Xc:B:i,start,end,concatemer_length,
        subread_idx,subread_total

    4. Modified bases (MM/ML):

        If modified base tags are found in the input then they will processed
        so that each subread has the correct ML/MM tags. Note that the 'mv' SAM
        tag is not handled in the same way and should be removed using the
        --remove-tags option if present.

    """
    logger = get_logger()
    logger.info("Digesting concatemers")
    input_files = list(find_files(input, glob=glob, recursive=recursive))
    header = get_alignment_header(source_files=input_files)
    remove_tags = list(set(remove_tags)) if remove_tags else None
    read_stream = get_concatemer_seqs(input_files, remove_tags=remove_tags)
    if max_reads:
        read_stream = islice(read_stream, max_reads)
    cutter = EnzymeCutter.from_name(enzyme)
    monomer_stream = (
        monomer.read_seq for read in read_stream for monomer in read.cut(cutter)    # noqa
    )

    with closing(get_monomer_writer(output_bam, header=header)) as writer:
        writer.consume(monomer_stream)
    # logger.info(
    #    f"Wrote {writer.base_counter:,} bases in "
    #    f"{writer.read_counter:,} reads to {output_path}"
    # )
    return writer


@app.command()
def parse_bam(
    bam: Path = typer.Argument(
        ..., help="A namesorted BAM of aligned + unaligned monomers"
    ),
    output_prefix: Path = typer.Argument(
        ..., help="Output files will share this prefix"
    ),
    force: bool = typer.Option(
        default=False, help="Overwite any existing files"),
    monomers: bool = typer.Option(
        default=True,
        help="Create a namesorted BAM file with pore-c annotations"
    ),
    chromunity: bool = typer.Option(
        default=False, help="Create a chromunity-compatible parquet"
    ),
    chromunity_merge_distance: Optional[int] = typer.Option(
        default=None,
        help=(
            "Merge co-linear monomers that are separated by less than this"
            "distance into a single monomer"
        ),
    ),
    summary: bool = False,
    direct_only: bool = typer.Option(
        default=False,
        help=(
            "Only output monomer pairs that are adjacent on the concatemer"
            ", don't do combinatorial expansion"
        ),
    ),
    paired_end: bool = typer.Option(
        default=False, help="Create a mock paired-end BAM"),
    paired_end_minimum_distance: Optional[int] = typer.Option(
        None,
        help=(
            "Filter out any pairs shorter than this distance."
            "Note setting this removes all trans pairs"
        ),
    ),
    paired_end_maximum_distance: Optional[int] = typer.Option(
        None,
        help=(
            "Filter out any pairs longer than this distance. "
            "Note setting this removes all trans pairs"
        ),
    ),
):
    """Parse a BAM file of aligned monomers."""
    # TODO: fill out documentation
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
    output_files = AnnotatedMonomerFC.with_prefix(
        output_prefix, drop=drop_outputs)
    if output_files.exists_any() and not force:
        logger.error(
            "Some of the outputs already exist,"
            "please remove before continuing"
        )
        raise IOError
    header = get_alignment_header(source_files=input_files)
    monomer_aligns = get_monomer_aligns(input_files)
    annotated_stream = annotate_monomer_alignments(monomer_aligns)
    with closing(
        AnnotatedMonomerWriter.from_file_collection(
            output_files,
            header=header,
            chromunity_merge_distance=chromunity_merge_distance,
            paired_end_minimum_distance=paired_end_minimum_distance,
            paired_end_maximum_distance=paired_end_maximum_distance,
        )
    ) as writer:
        writer.consume(annotated_stream, direct_only=direct_only)

    if summary:
        logger.info(f"Summary information at {output_files.summary_json}")
    # logger.info(
    #    f"Wrote {writer.base_counter:,} bases in "
    #    f"{writer.read_counter:,} reads to {output_path}"
    # )
    return writer


utils = typer.Typer()


@utils.command()
def export_bam(
    bam: Path,
    output_prefix: Path,
    force: bool = False,
    monomers: bool = False,
    paired_end: bool = False,
    chromunity: bool = False,
    summary: bool = False,
    direct_only: bool = False,
    chromunity_merge_distance: Optional[int] = None,
    paired_end_minimum_distance: Optional[int] = None,
    paired_end_maximum_distance: Optional[int] = None,
):
    """Export bam."""
    logger = get_logger()
    logger.info(f"Processing reads from {bam}")
    input_files = [bam]
    drop_outputs = []
    at_least_one = False
    for flag, filekey in [
        (monomers, "namesorted_bam"),
        (paired_end, "paired_end_bam"),
        (chromunity, "chromunity_parquet"),
        (summary, "summary_json"),
    ]:
        if not flag:
            drop_outputs.append(filekey)
        else:
            at_least_one = True

    if not at_least_one:
        logger.error(
            "You need to supply at least one of the "
            "--paired_end, --chromunity or --summary flags"
        )
        raise typer.Exit(code=1)
    output_files = AnnotatedMonomerFC.with_prefix(
        output_prefix, drop=drop_outputs)
    if not ("namesorted_bam" in drop_outputs):
        logger.warning(
            f"The output file {output_files.namesorted_bam} "
            "will have the same content as input"
        )

    if output_files.exists_any() and not force:
        logger.error(
            "Some of the outputs already exist,"
            "please remove before continuing"
        )
        raise IOError
    header = get_alignment_header(source_files=input_files)
    annotated_stream = group_aligns_by_concatemers(get_monomer_aligns([bam]))
    with closing(
        AnnotatedMonomerWriter.from_file_collection(
            output_files,
            header=header,
            chromunity_merge_distance=chromunity_merge_distance,
            paired_end_minimum_distance=paired_end_minimum_distance,
            paired_end_maximum_distance=paired_end_maximum_distance,
        )
    ) as writer:
        writer.consume(annotated_stream, direct_only=direct_only)

    if summary:
        logger.info(f"Summary information at {output_files.summary_json}")

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
    """Create test data."""
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
def create_chunked_ubam(
    input: Path = typer.Argument(
        ..., help="An unaligned BAM file or directory of same"
    ),
    output_prefix: Path = typer.Argument(
        ..., help="An unaligned BAM file with a"
        "separate record for each monomer"
    ),
    chunk_size: int = typer.Argument(
        ..., help="The number of reads in each chunk"),
    glob: str = typer.Option(
        default="*.bam",
        help="If INPUT is a directory use this glob to search for files",
    ),
    recursive: bool = typer.Option(
        default=True,
        help="If INPUT is a directory search recusively"
        "for files matching 'glob'",
    ),
    max_reads: int = typer.Option(
        default=0, help="Take the first n reads (useful for testing)"
    ),
):
    """Create chunked ubam."""
    logger = get_logger()
    input_files = list(find_files(input, glob=glob, recursive=recursive))
    if len(input_files) == 0:
        logger.error("No input files found at {input}")
        raise typer.Exit(code=1)

    if output_prefix.is_dir():
        output_prefix = output_prefix / "chunked"
    output_files = create_chunked_bam(
        input_files, output_prefix, chunk_size, max_reads=max_reads
    )
    logger.info(f"Wrote {len(output_files)} chunks")
    return output_files


app.add_typer(utils, name="utils")


def run_main():
    """Entry point."""
    app()


if __name__ == "__main__":
    run_main()
