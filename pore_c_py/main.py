"""Command-line interface."""
import argparse
from contextlib import closing
from itertools import islice
import logging
from pathlib import Path
import sys

from pore_c_py.aligns import annotate_monomer_alignments
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
from pore_c_py.log import get_named_logger, log_level
from pore_c_py.model import EnzymeCutter


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
        "output", type=Path,
        help="An unaligned BAM file with a separate record for each monomer")
    digest_parse.add_argument(
        "--recursive", default=False, action="store_true",
        help=(
            "If INPUT is a directory search recursively "
            "for files matching 'glob'")),
    digest_parse.add_argument(
        "--glob", default="*.fastq",
        help="If INPUT is a directory use this glob to search for files")
    digest_parse.add_argument(
        "--max_reads", type=int, default=0,  # should be inf?
        help="Take the first n reads (useful for testing)")
    digest_parse.add_argument(
        "--remove_tags", nargs="+", default=list(),
        help="Optionally remove SAM tags from input")

    bam_parse = subparsers.add_parser(
        "parse-bam",
        description=(
            "Filter a BAM of concatemer-grouped monomer alignments, "
            "and add 'walk' tag."),
        parents=[log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    bam_parse.set_defaults(func=parse_bam)
    bam_parse.add_argument(
        "bam", type=Path,
        help="A namesorted BAM of aligned + unaligned monomers.")
    bam_parse.add_argument(
        "output_prefix", type=Path,
        help="Output files will share this prefix")
    bam_parse.add_argument(
        "--force", default=False, action="store_true",
        help="Overwite any existing files")
    bam_parse.add_argument(
        "--monomers", default=False, action="store_true",
        help="Create a namesorted BAM file with pore-c annotations")
    bam_parse.add_argument(
        "--chromunity", default=False, action="store_true",
        help="Create a chromunity-compatible parquet")
    bam_parse.add_argument(
        "--chromunity_merge_distance", type=int,
        help=(
            "Merge co-linear monomers that are separated by less than this"
            "distance into a single monomer"))
    bam_parse.add_argument(
        "--summary", default=False, action="store_true")
    bam_parse.add_argument(
        "--direct_only", default=False, action="store_true",
        help=(
            "Only output monomer pairs that are adjacent on the concatemer"
            ", don't do combinatorial expansion."))
    bam_parse.add_argument(
        "--paired_end", default=False, action="store_true",
        help="Create a mock paired-end BAM")
    bam_parse.add_argument(
        "--paired_end_minimum_distance", type=int,
        help=(
            "Filter out any pairs shorter than this distance."
            "Note setting this removes all trans pairs."))
    bam_parse.add_argument(
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


def digest(args):
    """Digest entry point."""
    logger = get_named_logger("Digest")
    logger.info("Digesting concatemers.")
    input_files = list(find_files(
        args.input, glob=args.glob, recursive=args.recursive))
    header = get_alignment_header(source_files=args.input_files)
    args.remove_tags.append('mv')  # always remove this, as we don't chop it up
    read_stream = get_concatemer_seqs(
        input_files, remove_tags=args.remove_tags)
    if args.max_reads:
        read_stream = islice(read_stream, args.max_reads)
    cutter = EnzymeCutter.from_name(args.enzyme)
    monomer_stream = (
        monomer.read_seq
        for read in read_stream
        for monomer in read.cut(cutter))

    with closing(get_monomer_writer(args.output_bam, header=header)) as writer:
        writer.consume(monomer_stream)
    logger.info("Finished digestion.")


def parse_bam(args):
    """Parse BAM entry point."""
    logger = get_named_logger("ParseBAM")
    logger.info(f"Processing reads from {args.bam}")
    if all(
        x is False for x in (
            args.monomers, args.chromunity, args.summary, args.paired_end)):
        logger.error(
            "You must specify at least one of --monomers, --chromunity, "
            "--summary, or --paired_end")
        sys.exit(1)
    input_files = [args.bam]
    drop_outputs = []
    for flag, filekey in [
        (args.monomers, "namesorted_bam"),
        (args.paired_end, "paired_end_bam"),
        (args.chromunity, "chromunity_parquet"),
        (args.summary, "summary_json"),
    ]:
        if not flag:
            drop_outputs.append(filekey)
    output_files = AnnotatedMonomerFC.with_prefix(
        args.output_prefix, drop=drop_outputs)
    if output_files.exists_any() and not args.force:
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
            chromunity_merge_distance=args.chromunity_merge_distance,
            paired_end_minimum_distance=args.paired_end_minimum_distance,
            paired_end_maximum_distance=args.paired_end_maximum_distance,
        )
    ) as writer:
        writer.consume(annotated_stream, direct_only=args.direct_only)

    if args.summary:
        logger.info(f"Summary information at {args.output_files.summary_json}")
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
