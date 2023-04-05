"""Command-line interface."""
import argparse
import contextlib
import itertools
import logging
from pathlib import Path
import sys
import time

import pysam

from pore_c_py import align_tools, annotate, digest, utils, writers


def porec_parser():
    """Create CLI parser."""
    from pore_c_py import __version__
    parser = argparse.ArgumentParser(
        "pore-c-py", parents=[utils.log_level()],
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
            f"sortable as strings. The f{utils.CONCATEMER_ID_TAG} SAM tag "
            "is set to the original "
            f"read/concatemer ID. The f{utils.MONOMER_DATA_TAG} tag "
            "contains the relative "
            "coordinates of the monomer on the concatemer. "
            f"{utils.MONOMER_DATA_TAG}:B:i,start,end,concatemer_length,subread_idx,subread_total "  # noqa:E501
            "If modified base tags are found in the input then they will "
            "processed so that each subread has the correct ML/MM tags. "
            "Note that the 'mv' SAM tag is not handled in the same way "
            "and is removed."),
        parents=[utils.log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    digest_parse.set_defaults(func=digest_bam)
    digest_parse.add_argument(
        "input", type=Path, default=[sys.stdin], nargs='*',
        help="An unaligned BAM file or FASTQ file or a directory of same")
    digest_parse.add_argument(
        "enzyme",
        help="A restriction enzyme name eg. NlaIII")
    digest_parse.add_argument(
        "--output", type=Path, default=sys.stdout,
        help="An unaligned BAM file with a separate record for each monomer,"
        "default will write to standard out.")
    digest_parse.add_argument(
        "--header", type=Path,
        help="The output bam requires a header."
        "If INPUT parameter is stdin any existing headers will not be"
        "accessible so optionally provide an existing bam that"
        "contains headers to be copied to the output bam."
        "Otherwise the header will be 'pore-c-py'.")
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
        parents=[utils.log_level()],
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
        "--stdout", default=False,
        help="If True output name sorted bam to stdout")
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
        "--filter_pairs", action="store_true",
        help=(
            "Filter out pairs using min and max distance"))
    annotate_parse.add_argument(
        "--paired_end_minimum_distance", default=0, type=int,
        help=(
            "Filter out any pairs shorter than this distance."
            "Note setting this removes all trans pairs."))
    annotate_parse.add_argument(
        "--paired_end_maximum_distance", default=float("inf"), type=int,
        help=(
            "Filter out any pairs longer than this distance. "
            "Note setting this removes all trans pairs."))
    annotate_parse.add_argument(
        "--allow_singletons", default=False, action="store_true",
        help=(
            "If filter pairs is true, allow singletons in paired end bam"))
    annotate_parse.add_argument(
        "--allow_improper", default=False, action="store_true",
        help=(
            "If filter pairs is true, allow improper pairs in paired end bam"))
    annotate_parse.add_argument(
        "--allow_unmapped", default=False, action="store_true",
        help=(
            "If filter pairs is true, allow unmapped in paired end bam"))

    chunk_bam_parse = subparsers.add_parser(
        "chunk-bam", help="Chunk a BAM into pieces.",
        parents=[utils.log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    chunk_bam_parse.set_defaults(func=chunk_bam)
    chunk_bam_parse.add_argument(
        "input", type=Path,
        help="An unaligned BAM file or directory of same")
    chunk_bam_parse.add_argument(
        "output_prefix", type=Path,
        help="An unaligned BAM file with a separate record for each monomer")
    chunk_bam_parse.add_argument(
        "--chunk_size", type=int, default=8000,
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
        "--max_reads", type=int,
        help="Take the first n reads (useful for testing)")
    chunk_bam_parse.add_argument(
        "--threads", default=1, type=int,
        help="Compute threads of bam compression.")

    return parser


def digest_bam(args):
    """Digest entry point.

    Read in concatemer sequences from unaligned bam or fastq.
    Cuts each read with enzyme provided,
    and tags with relevant concatemer info.
    Outputs a bam containing each monomer sequence.
    """
    logger = utils.get_named_logger("Digest")
    logger.info(f"Digesting concatemers from {args.input}.")
    # Default header will be pore-c-py.
    # If input type is a file or directory, Header will be copied from that.
    # If input type is stdin, there is an optional header parameter available
    # to supply a preferred header present in an existing bam.
    get_header = args.header
    input_files = [args.input[0]]
    input_mode = "rb"
    pysam_kwargs = {"text": "pore-c-py"}
    if sys.stdin not in args.input:
        input_files = list(
            utils.find_files(
                args.input[0], glob=args.glob, recursive=args.recursive))
        get_header = str(input_files[0])
        input_mode = "r"
    if args.header or sys.stdin not in args.input:
        with pysam.AlignmentFile(
                get_header, input_mode, check_sq=False) as inputfile:
            header = align_tools.update_header(inputfile.header)
        pysam_kwargs = {"header": header}
    mode = "wb"
    if (args.output is not sys.stdout) and \
            (not utils.stdout_is_regular_file()):
        mode = "wb0"
    with pysam.AlignmentFile(
            args.output,
            threads=args.threads, mode=mode, **pysam_kwargs) as outbam:
        for input_file in input_files:
            with pysam.AlignmentFile(
                    input_file, input_mode, check_sq=False) as inputfile:
                for monomer in digest.get_concatemer_seqs(
                    inputfile, enzyme=args.enzyme,
                        remove_tags=args.remove_tags):
                    outbam.write(monomer)
    logger.info("Finished digestion.")


def annotate_bam(args):
    """Parse BAM entry point."""
    logger = utils.get_named_logger("AnntateBAM")
    logger.info(f"Processing reads from {args.bam}")
    if all(
        x is False for x in (
            args.monomers, args.chromunity, args.summary, args.paired_end)):
        logger.error(
            "You must specify at least one of --monomers, --chromunity, "
            "--summary, or --paired_end")
        sys.exit(1)
    paired_end_bam = Path(f"{args.output_prefix}.pe.bam")
    summary_json = Path(f"{args.output_prefix}.summary.json")
    chromunity_parquet = Path(f"{args.output_prefix}.chromunity.parquet")
    outputs = [paired_end_bam, summary_json, chromunity_parquet]
    mode = "wb"
    if args.stdout:
        namesorted_bam = sys.stdout
        if not utils.stdout_is_regular_file():
            mode = "wb0"
    else:
        namesorted_bam = Path(f"{args.output_prefix}.ns.bam")
        outputs += namesorted_bam
    if any(x.exists() for x in outputs) and not args.force:
        logger.error(
            "Some of the outputs exist, please remove before continuing.")
        sys.exit(1)

    # TODO: other outputs
    with contextlib.ExitStack() as manager:
        inbam = pysam.AlignmentFile(args.bam, "r", check_sq=False)
        manager.enter_context(inbam)

        if args.monomers:
            header = align_tools.update_header(inbam.header)
            outbam = pysam.AlignmentFile(
                namesorted_bam, mode=mode, header=header, threads=args.threads)
            manager.enter_context(outbam)
        if args.chromunity:
            chrom_writer = writers.ChromunityWriter(
                chromunity_parquet,
                merge_distance=args.chromunity_merge_distance)
            manager.enter_context(chrom_writer)
        if args.paired_end:
            pairs_bam = pysam.AlignmentFile(
                paired_end_bam, "wb", header=header, threads=args.threads)
            manager.enter_context(pairs_bam)
        if args.summary:
            summary = writers.StatsWriter(summary_json)
            manager.enter_context(summary)

        for concat_walk in annotate.annotate_alignments(inbam):
            if args.monomers:
                for aln in concat_walk:
                    outbam.write(aln)
            if args.chromunity:
                chrom_writer.write(concat_walk)
            if args.summary or args.paired_end:
                pairs = align_tools.get_pairs(concat_walk, args.direct_only)
                if args.filter_pairs:
                    pairs = align_tools.filter_pairs(
                        pairs,
                        allow_singletons=args.allow_singletons,
                        allow_improper=args.allow_improper,
                        allow_unmapped=args.allow_unmapped)
                for p in pairs:
                    if args.paired_end:
                        pairs_bam.write(p.left)
                        if p.right is not None:
                            pairs_bam.write(p.right)
                    if args.summary:
                        summary.append(p)
    logger.info("Finished BAM parsing.")


def chunk_bam(args):
    """Create chunked ubam."""
    logger = utils.get_named_logger("ChunkedBAM")
    input_files = list(utils.find_files(
        args.input, glob=args.glob, recursive=args.recursive))
    if len(input_files) == 0:
        logger.error("No input files found at {args.input}")
        sys.exit(1)

    if args.output_prefix.is_dir():
        args.output_prefix = args.output_prefix / "chunked"

    with pysam.AlignmentFile(
            str(input_files[0]), "r", check_sq=False) as inbam:
        header = align_tools.update_header(inbam.header)

    read_stream = itertools.chain(
        *(pysam.AlignmentFile(str(x), check_sq=False) for x in input_files))
    if args.max_reads is not None:
        read_stream = itertools.islice(read_stream, args.max_reads)

    def _new_file(index, bam):
        if bam is not None:
            bam.close()
        return index + 1, pysam.AlignmentFile(
            str(args.output_prefix.with_suffix(f".batch_{index}.bam")),
            threads=args.threads, mode="wb", header=header)

    def _log_time(start, index):
        elapsed = time.perf_counter() - start
        reads_per_minute = int(args.chunk_size * 60 / elapsed)
        logger.info(
            f"Wrote batch {index}: {reads_per_minute} reads/min.")
        return time.perf_counter()

    bam = None
    index, written = 0, 0
    start = time.perf_counter()
    for align in read_stream:
        if written == 0:
            index, bam = _new_file(index, bam)
        bam.write(align)
        written += 1
        if written % args.chunk_size == 0:
            written = 0
            start = _log_time(start, index)
    if written != 0:
        _log_time(start, index)
    if bam is not None:
        bam.close()
    logger.info("Finished repacking BAMs.")


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
