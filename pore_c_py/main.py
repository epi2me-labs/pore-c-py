"""Command-line interface."""
import argparse
import contextlib
import logging
import os
from pathlib import Path
import stat
import sys

import pysam

from pore_c_py import annotate, writers
from pore_c_py import digests
from pore_c_py.annotate import update_header
from pore_c_py.io import (
    create_chunked_bam,
    find_files
)
from pore_c_py.log import get_named_logger, log_level


def porec_parser():
    """Create CLI parser."""
    from pore_c_py import __version__
    parser = argparse.ArgumentParser(
        "pore-c-py", parents=[log_level()],
        description="Tools for processing Pore-C data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(
        title='subcommands', description='valid commands',
        help='additional help', dest='command')
    subparsers.required = True

    parser.add_argument(
        '--version', action='version',
        version='%(prog)s {}'.format(__version__))

    digest_parse = subparsers.add_parser(
        'digest',
        description=(
            "Digest concatemer sequences into monomers "
            "using a restriction enzyme."),
        help=(
            "Resulting BAM file has query names of the format "
            "<original_read_id>:<subread_start>:<subread_end> with "
            "coordinates left padded with zeroes such that they are "
            "sortable as strings. The MI SAM tag is set to the original "
            "read/concatemer ID. The Xc tag contains the relative "
            "coordinates of the monomer on the concatemer. "
            "Xc:B:i,start,end,concatemer_length,subread_idx,subread_total "
            "If modified base tags are found in the input then they will "
            "processed so that each subread has the correct ML/MM tags. "
            "Note that the 'mv' SAM tag is not handled in the same way "
            "and is removed."),
        parents=[log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    digest_parse.set_defaults(func=digest)
    digest_parse.add_argument(
        "input", type=Path,
        help="An unaligned BAM file or FASTQ file or a directory of same")
    digest_parse.add_argument(
        "enzyme",
        help="A restriction enzyme name eg. NlaIII")
    digest_parse.add_argument(
        "--output", type=Path, default=sys.stdout,
        help="An unaligned BAM file with a separate record for each monomer,"
        "default will write to standard out.")
    digest_parse.add_argument(
        "--recursive", default=False, action="store_true",
        help=(
            "If INPUT is a directory search recursively "
            "for files matching 'glob'")),
    digest_parse.add_argument(
        "--glob", default="*.bam",
        help="If INPUT is a directory use this glob to search for files")
    digest_parse.add_argument(
        "--max_reads", type=int, default=0,  # should be inf?
        help="Take the first n reads (useful for testing)")
    digest_parse.add_argument(
        "--remove_tags", nargs="+", default=list(),
        help="Optionally remove SAM tags from input")
    digest_parse.add_argument(
        "--threads", default=1, type=int,
        help="Compute threads of bam compression.")

    annotate_parse = subparsers.add_parser(
        "annotate",
        description=(
            "Filter a BAM of concatemer-grouped monomer alignments, "
            "and add 'walk' tag."),
        parents=[log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    annotate_parse.set_defaults(func=annotate_bam)
    annotate_parse.add_argument(
        "bam", type=Path,
        help=(
            "A (concatemer) namesorted BAM of aligned "
            "and unaligned monomers."))
    annotate_parse.add_argument(
        "output_prefix", type=Path,
        help="Output files will share this prefix")
    annotate_parse.add_argument(
        "--threads", default=1, type=int,
        help="Compute threads of bam compression.")
    annotate_parse.add_argument(
        "--force", default=False, action="store_true",
        help="Overwite any existing files")
    annotate_parse.add_argument(
        "--monomers", default=False, action="store_true",
        help="Create a namesorted BAM file with Pore-C annotations")
    annotate_parse.add_argument(
        "--chromunity", default=False, action="store_true",
        help="Create a chromunity-compatible parquet")
    annotate_parse.add_argument(
        "--chromunity_merge_distance", type=int,
        help=(
            "Merge co-linear monomers that are separated by less than this"
            "distance into a single monomer"))
    annotate_parse.add_argument(
        "--summary", default=False, action="store_true")
    annotate_parse.add_argument(
        "--direct_only", default=False, action="store_true",
        help=(
            "Only output monomer pairs that are adjacent on the concatemer"
            ", don't do combinatorial expansion."))
    annotate_parse.add_argument(
        "--paired_end", default=False, action="store_true",
        help="Create a mock paired-end BAM")
    annotate_parse.add_argument(
        "--paired_end_minimum_distance", type=int,
        help=(
            "Filter out any pairs shorter than this distance."
            "Note setting this removes all trans pairs."))
    annotate_parse.add_argument(
        "--paired_end_maximum_distance", type=int,
        help=(
            "Filter out any pairs longer than this distance. "
            "Note setting this removes all trans pairs."))

    chunk_bam_parse = subparsers.add_parser(
        "chunk-bam", help="Chunk a BAM into pieces.",
        parents=[log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    chunk_bam_parse.set_defaults(func=chunk_bam)
    chunk_bam_parse.add_argument(
        "input", type=Path,
        help="An unaligned BAM file or directory of same")
    chunk_bam_parse.add_argument(
        "output_prefix", type=Path,
        help="An unaligned BAM file with a separate record for each monomer")
    chunk_bam_parse.add_argument(
        "chunk_size", type=int,
        help="The number of reads in each chunk"),
    chunk_bam_parse.add_argument(
        "--recursive", default=False, action="store_true",
        help=(
            "If INPUT is a directory search recursively for "
            "files matching 'glob'"))
    chunk_bam_parse.add_argument(
        "--glob", default="*.fastq",
        help="If INPUT is a directory use this glob to search for files")
    chunk_bam_parse.add_argument(
        "max_reads", type=int, default=0,
        help="Take the first n reads (useful for testing)")

    return parser


def stdout_is_regular_file() -> bool:
    """
    Detect if standard output is a regular file (or say a pipe).

    :return: True if stdout is a regular file, else False.
    """
    mode = os.fstat(sys.stdout.buffer.fileno()).st_mode
    return stat.S_ISREG(mode)


def digest(args):
    """Digest entry point."""
    logger = get_named_logger("Digest")
    logger.info(f"Digesting concatemers from {args.input}.")
    input_files = list(
        find_files(
            args.input, glob=args.glob, recursive=args.recursive))
    with pysam.AlignmentFile(
            str(input_files[0]), "r", check_sq=False) as inputfile:
        header = update_header(inputfile.header)
    mode = "wb"
    if (args.output is not sys.stdout) and (not stdout_is_regular_file()):
        mode = "wb0"
    with pysam.AlignmentFile(
            args.output, header=header,
            threads=args.threads, mode=mode) as outbam:
        for input_file in input_files:
            with pysam.AlignmentFile(
                    str(input_file), "r", check_sq=False) as inputfile:
                for monomer in digests.get_concatemer_seqs(
                    inputfile, enzyme=args.enzyme,
                        remove_tags=args.remove_tags):
                    outbam.write(monomer)
    logger.info("Finished digestion.")


def annotate_bam(args):
    """Parse BAM entry point."""
    logger = get_named_logger("AnntateBAM")
    logger.info(f"Processing reads from {args.bam}")
    if all(
        x is False for x in (
            args.monomers, args.chromunity, args.summary, args.paired_end)):
        logger.error(
            "You must specify at least one of --monomers, --chromunity, "
            "--summary, or --paired_end")
        sys.exit(1)

    namesorted_bam = Path(f"{args.output_prefix}.ns.bam")
    paired_end_bam = Path(f"{args.output_prefix}.pe.bam")
    summary_json = Path(f"{args.output_prefix}.summary.json")
    chromunity_parquet = Path(f"{args.output_prefix}.chromunity.parquet")
    outputs = (
        namesorted_bam, paired_end_bam, summary_json, chromunity_parquet)

    if any(x.exists() for x in outputs) and not args.force:
        logger.error(
            "Some of the outputs exist, please remove before continuing.")
        sys.exit(1)

    # TODO: other outputs
    with contextlib.ExitStack() as manager:
        inbam = pysam.AlignmentFile(args.bam, "r", check_sq=False)
        manager.enter_context(inbam)

        if args.monomers:
            header = annotate.update_header(inbam.header)
            outbam = pysam.AlignmentFile(
                namesorted_bam, "wb", header=header, threads=args.threads)
            manager.enter_context(outbam)
        if args.chromunity:
            chrom_writer = writers.ChromunityWriter(
                chromunity_parquet,
                merge_distance=args.chromunity_merge_distance)
            manager.enter_context(chrom_writer)

        for concat_walk in annotate.annotate_alignments(inbam):
            if args.monomers:
                for aln in concat_walk:
                    outbam.write(aln)
            if args.chromunity:
                chrom_writer.write(concat_walk)
            if args.summary or args.paired_end:
                # make pairs, send to writers.
                if args.summary:
                    pass
                if args.paired_end:
                    pass
    logger.info("Finished BAM parsing.")


def chunk_bam(args):
    """Create chunked ubam."""
    logger = get_named_logger()
    input_files = list(find_files(
        input, glob=args.glob, recursive=args.recursive))
    if len(input_files) == 0:
        logger.error("No input files found at {args.input}")
    else:
        if args.output_prefix.is_dir():
            args.output_prefix = args.output_prefix / "chunked"
        output_files = create_chunked_bam(
            input_files, args.output_prefix, args.chunk_size,
            max_reads=args.max_reads)
        logger.info(f"Wrote {len(output_files)} chunks")


def run_main():
    """Entry point."""
    parser = porec_parser()
    args = parser.parse_args()

    logging.basicConfig(
        format='[%(asctime)s - %(name)s] %(message)s',
        datefmt='%H:%M:%S', level=logging.INFO)
    logger = logging.getLogger(__package__)
    logger.setLevel(args.log_level)
    if args.logfile:
        fh = logging.FileHandler(args.logfile, "w")
        fh.setLevel(args.log_level)
        logger.addHandler(fh)

    args.func(args)
