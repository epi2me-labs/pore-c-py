import subprocess as sp
from dataclasses import dataclass
from itertools import groupby
from pathlib import Path
from typing import List

import pytest
from loguru import logger
from pysam import AlignedSegment, AlignmentFile, FastxFile

from pore_c2.digest import _get_enzyme, sequence_to_read_fragments
from pore_c2.overlaps import FragmentOverlapper


@dataclass
class ConcatemerAlignData:
    read_name: str

    raw_alignments: List[AlignedSegment]

    def __repr__(self):
        return (
            f"<ConcatemerAligns read={self.read_name} "
            f"raw_aligns={len(self.raw_alignments)}>"
        )


def test_mappy(scenario):
    import mappy as mp

    a = mp.Aligner(str(scenario.reference_fasta), preset="map-ont")
    enzyme = _get_enzyme(scenario.enzyme)
    for rec in FastxFile(scenario.concatemer_fastq):
        print(rec.comment)
        for read_frag in sequence_to_read_fragments(enzyme, rec):
            seq, _ = read_frag.slice_fastq(rec)
            for hit in a.map(seq):
                if not hit.is_primary:
                    continue
                print(
                    read_frag.read_fragment_id,
                    hit.ctg,
                    hit.r_st,
                    hit.r_en,
                    hit.q_st,
                    hit.q_en,
                )


@pytest.mark.skip
def test_map(scenario):
    logger.debug(scenario)
    ref_fasta = scenario.reference_fasta
    idx_file = Path(str(ref_fasta) + ".bwt")
    if not idx_file.exists():
        sp.run(["bwa", "index", str(ref_fasta)], capture_output=True, check=True)
        # logger.debug(list(ref_fasta.parent.glob("*.*")))
    fastq = scenario.concatemer_fastq
    threads = 1
    bam_file = fastq.with_suffix(".sam")
    sp.run(
        f"bwa bwasw -t {threads} -b 5 -q 2 -r 1 -T 15 -z 10 "
        f"{ref_fasta} {fastq} | samtools sort -n -o {bam_file} -",
        shell=True,
        capture_output=True,
    )

    ol = FragmentOverlapper.from_fragments_df(scenario.fragments_df)
    af = AlignmentFile(str(bam_file))
    for query_name, aligns in groupby(af, lambda x: x.query_name):
        aligns = list(aligns)
        concatemer_id, fragments = query_name.split(":")
        fragments = [int(_) for _ in fragments.split("_")]
        d = ConcatemerAlignData(
            concatemer_id, sorted(aligns, key=lambda x: x.query_alignment_start)
        )
        for a in d.raw_alignments:
            logger.debug(
                f"{a.query_name}:{a.query_alignment_start}-{a.query_alignment_end} "
                f"{a.reference_start}-{a.reference_end}"
            )
            try:
                olap = ol.find_overlap(
                    a.reference_name, a.reference_start, a.reference_end
                )
            except ZeroDivisionError:
                print(a)
    # res = simulate_fasta_with_cut_sites(fasta, lengths, expected, enzyme='EcoRI')
    # logger.debug(res)
