from typing import List, Optional, Tuple

import mappy as mp
from attr import Factory, define

from .digest import Cutter, ReadFragment, sequence_to_read_fragments
from .overlaps import FragmentOverlapper, OverlapStats
from .reads import AlignData


@define(kw_only=True)
class ReadFragmentAligns:
    read_fragment: ReadFragment
    aligns: List[mp.Alignment]
    overlaps: List[OverlapStats]


@define(kw_only=True)
class ConcatemerAlignData:
    read_id: str
    seq: str
    qual: str
    tags: List[str] = Factory(list)
    aligns: List[ReadFragmentAligns]

    def to_sam(self) -> List[str]:
        res = []
        read_id = self.read_id
        for read_frag in self.aligns:
            frag_id = read_frag.read_fragment.read_id
            s, e = read_frag.read_fragment.read_start, read_frag.read_fragment.read_end
            seq, qual = self.seq[s:e], self.qual[s:e]
            rc_seq = None
            seq_len = len(seq)
            if len(read_frag.aligns) == 0:
                flag = 4
                sam = (
                    # QNAME,FLAG,RNAME,POS
                    f"{frag_id}\t{flag}\t*\t0\t"
                    # MAPQ,CIGAR,SEQ,QUAL
                    f"0\t*\t"
                    # RNEXT,PNEXT,TLEN,
                    "*\t0\t0\t"
                    # SEQ, QUAL
                    f"{seq}\t{qual}\t"
                    # Optional fields
                    f"NM:i:0\tMI:Z:{read_id}"
                )
                res.append(sam)
            else:
                for rec in read_frag.aligns:
                    # TODO: primary/secondary/supplementary
                    flag = 0
                    clip_end_bases = seq_len - rec.q_en
                    clip_st = f"{rec.q_st}S" if rec.q_st else ""
                    clip_en = f"{clip_end_bases}S" if clip_end_bases else ""
                    if rec.strand == +1:
                        flag = 0
                        cigar = f"{clip_st}{rec.cigar_str}{clip_en}"
                        _seq = seq
                    else:
                        flag = 16
                        cigar = f"{clip_en}{rec.cigar_str}{clip_st}"
                        if rc_seq is None:
                            rc_seq = mp.revcomp(seq)
                        _seq = rc_seq
                    sam = (
                        # QNAME,FLAG,RNAME,POS
                        f"{frag_id}\t{flag}\t{rec.ctg}\t{rec.r_st + 1}\t"
                        # MAPQ,CIGAR,SEQ,QUAL
                        f"{rec.mapq}\t{cigar}\t"
                        # RNEXT,PNEXT,TLEN,
                        "*\t0\t0\t"
                        # SEQ, QUAL
                        f"{_seq}\t{qual}\t"
                        # Optional fields
                        f"NM:i:{rec.NM}\tMD:Z:{rec.MD}\tMI:Z:{read_id}"
                    )
                    res.append(sam)
        return res


def split_and_map_read(
    *,
    aligner: mp.Aligner,
    read: AlignData,
    cutter: Cutter,
    overlapper: FragmentOverlapper,
) -> ConcatemerAlignData:
    aligns = []
    for read_frag in sequence_to_read_fragments(cutter, read):
        seq, _ = read_frag.slice_fastq(read)
        hits = list(aligner.map(seq))
        overlaps = [overlapper.overlaps(h.ctg, h.r_st, h.r_en) for h in hits]
        aligns.append(
            ReadFragmentAligns(read_fragment=read_frag, aligns=hits, overlaps=overlaps)
        )
    return ConcatemerAlignData(
        read_id=read.name, aligns=aligns, seq=read.seq, qual=read.qual
    )


def map_concatemer_read(
    *,
    aligner: mp.Aligner,
    read: AlignData,
    cutter: Cutter,
    thread_buf: Optional[mp.ThreadBuffer] = None,
) -> Tuple[AlignData, List[ReadFragment], List[mp.Alignment]]:
    read_frags = sequence_to_read_fragments(cutter, read, store_sequence=False)
    hits = [
        list(aligner.map(read.seq[_.read_start : _.read_end], buf=thread_buf))
        for _ in read_frags
    ]
    return (read, read_frags, hits)


def merge_colinear(c_align_data: ConcatemerAlignData) -> ConcatemerAlignData:
    return c_align_data


def pick_optimal(c_align_data: ConcatemerAlignData) -> ConcatemerAlignData:
    return c_align_data
