"""Digestion of unaligned concatemers."""
import copy

from Bio import Restriction
from Bio.Seq import Seq

from pore_c_py import utils


logger = utils.get_named_logger("Digest")


def get_subread_modified_bases(align, start, end):
    """Get modified bases subread.

    :param align: pysam.AlignedSegment.
    :param start: start coordinate to trim to.
    :param end: exclusive end co-ordinate.
    """
    mm_str, ml_str = "", ""
    base_indices = {}
    seq = align.query_sequence[start:end]
    for mod_key, mod_data in align.modified_bases.items():
        # find the modifications that overlap the subread
        idx = [
            x for x in range(len(mod_data)) if
            start <= mod_data[x][0] < end]
        probs_dic = dict(mod_data)
        if not idx:  # no mod bases (of this type) in the subread
            continue
        try:
            canonical_base, strand, skip_scheme, mod_type = mod_key
        except ValueError:
            canonical_base, strand, mod_type = mod_key
            skip_scheme = ""
        if canonical_base == "N":
            base_indices[canonical_base] = list(range(len(seq)))
        elif canonical_base not in base_indices:
            base_indices[canonical_base] = [
                x for x, b in enumerate(seq) if b.upper() == canonical_base
            ]
        base_offsets, probs = zip(*[mod_data[i] for i in idx])
        strand = "+" if strand == 0 else "-"
        deltas = []
        counter = 0
        probs = []
        for seq_idx in base_indices[canonical_base]:
            orig_idx = seq_idx + start
            if orig_idx in base_offsets:  # is modified
                probs += [probs_dic[orig_idx]]
                deltas.append(str(counter))
                counter = 0
            else:
                counter += 1
            deltas_formatted = ','.join(deltas)
            prob_str = ",".join(map(str, probs))
        mm_str += (
                f"{canonical_base}{strand}{mod_type}{skip_scheme}"
                f",{deltas_formatted};"
            )
        ml_str += f"{prob_str},"
    return mm_str, ml_str


def splits_to_intervals(positions, length):
    """Split to intervals."""
    if len(positions) == 0:
        return [(0, length)]
    prefix, suffix = [], []
    if positions[0] != 0:
        prefix = [0]
    if positions[-1] != length:
        suffix = [length]
    breaks = prefix + positions + suffix
    return [(start, end) for start, end in zip(breaks[:-1], breaks[1:])]


def get_enzyme(enzyme):
    """Get enzyme."""
    enz = getattr(Restriction, enzyme, None)
    if enz is None:
        raise ValueError(f"Enzyme not found: {enzyme}")
    if enz.cut_twice():
        raise NotImplementedError(
            f"Enzyme cuts twice, not currently supported: {enzyme}"
        )
    return enz


def digest_sequence(align, enzyme, tags_remove=None):
    """Digest sequence."""
    # the move tag massively bloats files, and we don't care for
    # it or handle it in trimming, so force its removal by default.
    if tags_remove is None:
        tags_remove = {'mv'}

    concatemer_id = align.query_name
    cut_points = [x - 1 for x in enzyme.search(Seq(align.query_sequence))]
    read_length = len(align.query_sequence)
    num_digits = len(str(read_length))
    intervals = splits_to_intervals(cut_points, read_length)
    num_intervals = len(intervals)

    for idx, (start, end) in enumerate(intervals):
        read = copy.copy(align)
        # trim the sequence and quality
        seq = align.query_sequence[start:end]
        qual = None
        if align.query_qualities:
            qual = align.query_qualities[start:end]
        read.query_sequence = seq
        read.query_qualities = qual
        # deal with mods, upgrading tag from interim to approved spec
        if ('Mm' in tags_remove) or ('MM' in tags_remove):
            # Setting to None effectively deletes,
            # or does nothing if not present
            for tag in ('Mm', 'Ml', 'MM', 'ML'):
                read.set_tag(tag, None)
        else:
            mm, ml = get_subread_modified_bases(align, start, end)
            for tag in ('Mm', 'Ml'):
                read.set_tag(tag, None)
            read.set_tag("MM", mm)
            read.set_tag("ML", ml)
        # lexograpically sortable monomer ID
        read.query_name = \
            f"{concatemer_id}:{start:0{num_digits}d}:{end:0{num_digits}d}"
        read.set_tag(
            utils.MONOMER_DATA_TAG,
            [start, end, read_length, idx, num_intervals])
        utils.MonomerData.set_monomer_data(
                read, start, end, read_length, idx, num_intervals)
        read.set_tag(utils.CONCATEMER_ID_TAG, concatemer_id, "Z")
        yield read


def get_concatemer_seqs(input_file, enzyme, remove_tags=None):
    """Digest concatemers in to unaligned monomers.

    :param input_file: pysam.AlignmentFile input
    :param enzyme: Name of the digestion enzyme used in sample preperation
    :param remove_tags: Comma seperated list of additional
                        tags to remove from file

    Concatemers are split with chosen digestion enzyme in to monomers.
    Monomers are tagged with unique monomer id and concatemer info.
    """
    logger.info(f"Digesting unaligned sequences from {input_file}")
    n_concatemers = 0
    n_monomers = 0
    enzyme = get_enzyme(enzyme)
    tags_remove = {"mv"}
    if remove_tags:
        tags_remove.update(set(remove_tags))
    for align in input_file.fetch(until_eof=True):
        n_concatemers += 1
        reads = digest_sequence(
            align, enzyme, tags_remove)
        for read in reads:
            n_monomers += 1
            yield read
    logger.info(
        f"Found {n_monomers} monomers in {n_concatemers} concatemers.")
