"""Annotation of monomer alignments."""

from itertools import groupby
import sys

import pysam

import pore_c_py
from pore_c_py.log import get_named_logger

logger = get_named_logger("AnntateAln")


def get_monomer_tag(monomer):
    """Get monomer tag."""
    read_start, read_end, _, _, _ = monomer.get_tag("Xc")
    if monomer.is_unmapped:
        monomer_tag = (
            f"*:"
            f"{read_start}-{read_end}")
    else:
        orientation = "+-"[monomer.is_reverse]
        monomer_tag = (
            f"{monomer.reference_name}:{orientation}:"
            f"{monomer.reference_start}-{monomer.reference_end}:"
            f"{read_start}-{read_end}")
    return monomer_tag


def sort_by_category(alns, monomer_id):
    """Sort by category."""
    reordered = list()
    for aln in alns:
        # original claims this as "primary", thats not technically accurate
        category = 0
        if aln.is_secondary:
            category = 3
        elif aln.is_supplementary:
            category = 2
        elif aln.is_unmapped:
            category = 1
        reordered.append(list([aln, category]))
    # create a generator to be consistent with input
    reordered.sort(key=lambda a: a[1])
    if reordered[0][1] not in {0, 1}:
        logger.warning(
            f"Warning: best alignment for monomer: {monomer_id} "
            f"is secondary or supplementary")
    yield from (x[0] for x in reordered)


def update_header(header):
    """Add PG tag to existing header.

    :param bamfile: a pysam.AlignmentFile.header

    """
    header = header.to_dict()
    pg = header.pop("PG", [])
    name = __package__.replace("_", "-")
    pg_data = {
        "ID": f"{name}-{len(pg) + 1}",
        "PN": name,
        "VN": pore_c_py.__version__,
        "CL": " ".join(sys.argv)}
    if len(pg) > 0:
        if "ID" in pg[-1]:
            pg_data["PP"] = pg[-1]["ID"]
    pg.append(pg_data)
    header["PG"] = pg
    return header


def annotate_alignments(input_bam):
    """Annotate (and filter) monomer alignments.

    :param input_bam: The input BAM file.
    :yields: the file header following by lists of alignments by concatemer.

    Alignments are group by concatemer (the input should be sorted by
    concatemer stored in the MI tag), and one alignment is selected per
    monomer. Alignments are annotated with a "walk", which simply enumerates
    how the monomers were linked together.
    """
    logger.info(f"Annotating alignments from file {input_bam}.")

    n_monomers = 0
    seen = set()
    with pysam.AlignmentFile(input_bam, "r", check_sq=False) as bamfile:
        # first yield header, this allows the input to be stdin and for us
        # to still be able to create a new BAM from the stream.
        yield bamfile.header
        # group alignments by concatemer and select one alignment per monomer
        for concat_id, aligns in groupby(
                bamfile.fetch(until_eof=True),
                lambda x: x.get_tag("MI")):
            if concat_id in seen:
                logger.warning(
                    f"Input file appears unsorted, found {concat_id} "
                    "more than once.")
            seen.add(concat_id)
            walk = list()
            for monomer_id, alns in groupby(aligns, lambda x: x.query_name):
                n_monomers += 1
                sorted_aligns = sort_by_category(alns, monomer_id)
                aln = next(sorted_aligns)
                walk.append(aln)
            expected_monomers = walk[-1].get_tag("Xc")[4]
            concat_id = walk[-1].get_tag("MI")
            if expected_monomers != len(walk):
                logger.warning(
                    f"Expected to see {expected_monomers} alignments for "
                    f"concatemer {concat_id}, found {len(walk)}"
                )
            walk_tag = ";".join(get_monomer_tag(aln) for aln in walk)
            for aln in walk:
                aln.set_tag("Xw", walk_tag)
            yield walk
    logger.info(f"Found {n_monomers} monomers in {len(seen)} concatemers.")
