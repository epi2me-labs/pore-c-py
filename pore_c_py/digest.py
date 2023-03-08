"""Digestion of unaligned concatemers."""
import copy

from Bio import Restriction
from Bio.Seq import Seq

from pore_c_py import utils

logger = utils.get_named_logger("Digest")


def get_subread(align, start, end):
    """Get sub read."""
    seq = align.query_sequence[start:end]
    if align.query_qualities:
        qual = align.query_qualities[start:end]
    else:
        qual = None
    if align.modified_bases:
        mm_str, ml_str = "", ""
        base_indices = {}
        for mod_key, mod_data in align.modified_bases.items():
            # find the modifications that overlap the subread
            idx = [
                x for x in range(len(mod_data)) if
                start <= mod_data[x][0] < end]
            if not idx:  # no mods in this subread
                continue
            try:
                canonical_base, strand, skip_scheme, mod_type = mod_key
            except ValueError:
                canonical_base, strand, mod_type = mod_key
                skip_scheme = ""

            canonical_base = mod_key[0]
            # find the positions of that canonical base in the subread
            if canonical_base not in base_indices:
                base_indices[canonical_base] = [
                    x for x, b in enumerate(seq) if b.upper() == canonical_base
                ]

            base_offsets, probs = zip(*[mod_data[_] for _ in idx])
            deltas = []
            counter = 0
            for seq_idx in base_indices[canonical_base]:
                orig_idx = seq_idx + start
                if orig_idx in base_offsets:  # is modified
                    deltas.append(str(counter))
                    counter = 0
                else:
                    counter += 1
            assert len(deltas) == len(probs)
            prob_str = ",".join(map(str, probs))
            strand = "+" if strand == 0 else "-"
            mm_str += (
                f"{canonical_base}{strand}{mod_type}{skip_scheme}"
                f",{','.join(deltas)};"
            )
            ml_str += f"{canonical_base},{prob_str};"
    else:
        mm_str, ml_str = None, None

    return seq, qual, mm_str, ml_str


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


def get_concatemer_seqs(input_file, enzyme, remove_tags):
    """Digest concatemers in to unaligned monomers.

    :param input_file: The input BAM or Fastq file.
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
    tags_remove = {"Ml", "ML", "Mm", "MM", "mv"}
    if remove_tags:
        tags_remove.update(set(remove_tags))
    for align in input_file.fetch(until_eof=True):
        n_concatemers += 1
        concatemer_id = align.query_name
        cut_points = [x - 1 for x in enzyme.search(Seq(align.query_sequence))]
        read_length = len(align.query_sequence)
        num_digits = len(str(read_length))
        intervals = splits_to_intervals(cut_points, read_length)
        num_intervals = len(intervals)
        for idx, (start, end) in enumerate(intervals):
            n_monomers += 1
            read = copy.copy(align)
            seq, qual, mm_str, ml_str = get_subread(
                align=align, start=start, end=end)
            for item in tags_remove:
                read.set_tag(item, None)
            read.query_sequence = seq
            read.query_qualities = qual
            # lexograpically sortable monomer ID
            read.query_name = \
                f"{concatemer_id}:{start:0{num_digits}d}:{end:0{num_digits}d}"
            utils.MonomerData.set_monomer_data(
                read, start, end, read_length, idx, num_intervals)
            read.set_tag(utils.CONCATEMER_ID_TAG, concatemer_id, "Z")
            read.set_tag("MM", mm_str)
            read.set_tag("ML", ml_str)
            yield read

    logger.info(f"Found {n_monomers} monomers in {n_concatemers} concatemers.")
