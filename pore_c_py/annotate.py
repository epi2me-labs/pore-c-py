"""Annotation of monomer alignments."""
from itertools import groupby

from pore_c_py.utils import get_named_logger

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


def annotate_alignments(bamfile):
    """Annotate (and filter) monomer alignments.

    :param bamfile: pysam.AlignmentFile
    :yields: the file header following by lists of alignments by concatemer.

    Alignments are grouped by concatemer (the input should be sorted by
    concatemer stored in the MI tag), and one alignment is selected per
    monomer. Alignments are annotated with a "walk", which simply enumerates
    how the monomers were linked together.
    """
    n_monomers = 0
    seen = set()
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
