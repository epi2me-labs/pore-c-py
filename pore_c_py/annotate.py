"""Annotation of monomer alignments."""
from itertools import groupby

from pore_c_py import utils

logger = utils.get_named_logger("AnntateAln")


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


def get_walk(aligns):
    """Get walk."""
    n_monomers = 0
    walk = list()
    for monomer_id, alns in groupby(aligns, lambda x: x.query_name):
        n_monomers += 1
        sorted_aligns = sort_by_category(alns, monomer_id)
        aln = next(sorted_aligns)
        walk.append(aln)
    expected_monomers = utils.MonomerData.get_subread_total(walk[-1])
    if expected_monomers != len(walk):
        concat_id = walk[-1].get_tag(utils.CONCATEMER_ID_TAG)
        logger.warning(
            f"Expected to see {expected_monomers} alignments for "
            f"concatemer {concat_id}, found {len(walk)}"
        )
    walk_tag = ";".join(utils.MonomerData.from_pysam(aln).name for aln in walk)
    for aln in walk:
        aln.set_tag(utils.WALK_TAG, walk_tag)
    return (walk, n_monomers)


def annotate_alignments(bamfile):
    """Annotate (and filter) monomer alignments.

    :param bamfile: pysam.AlignmentFile
    :yields: the file header following by lists of alignments by concatemer.

    Alignments are grouped by concatemer (the input should be sorted by
    concatemer stored in the concatemer ID tag), and one alignment is selected
    per monomer. Alignments are annotated with a "walk", which simply
    enumerates the alignment coordinates of the monomers comprising the
    concatemer.
    """
    total_monomers = 0
    seen = set()
    for concat_id, aligns in groupby(
            bamfile.fetch(until_eof=True),
            lambda x: x.get_tag(utils.CONCATEMER_ID_TAG)):
        if concat_id in seen:
            logger.warning(
                f"Input file appears unsorted, found {concat_id} "
                "more than once.")
        seen.add(concat_id)
        walk, n_monomers = get_walk(aligns)
        total_monomers += n_monomers
        yield walk
    logger.info(f"Found {total_monomers} monomers in {len(seen)} concatemers.")
