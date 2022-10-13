import array
import re
from abc import ABCMeta, abstractmethod
from typing import Any, Dict, List, Mapping, Optional, Tuple

from attrs import Factory, define
from Bio.Seq import Seq
from pysam import AlignedSegment, FastxRecord

from .settings import DEFAULT_ALIGN_HEADER, FASTQ_TAG_RE, MOD_TAGS

# copied from pysam.libcalignedsegment.pyx
SAM_TYPES = "iiiiiif"
HTSLIB_TYPES = "cCsSiIf"
PARRAY_TYPES = "bBhHiIf"

MM_TAG_SKIP_SCHEME_RE = re.compile(r"[.?]")

ARRAY_TO_HTSLIB_TRANS = str.maketrans(PARRAY_TYPES, HTSLIB_TYPES)
PYSAM_TO_SAM_TRANS = str.maketrans(HTSLIB_TYPES, SAM_TYPES)


TAG_MI_RE = re.compile(r"MI\:Z\:(\S+)")
# XC==concatemer metadata
TAG_XC_RE = re.compile(
    r"Xc:B:i,(?P<start>\d+),(?P<end>\d+),(?P<subread_idx>\d+),(?P<subread_total>\d+)"
)
# XW==walk metadata
TAG_XW_RE = re.compile(r"Xw\:Z\:(.+)")  # TODO fill this out


def tag_tuple_to_str(key: str, val: Any, value_type: str):
    # TODO not sure if this works for H
    value_type = value_type.translate(PYSAM_TO_SAM_TRANS)
    if value_type == "B":
        assert isinstance(val, array.array)
        elemtype = val.typecode.translate(ARRAY_TO_HTSLIB_TRANS)
        val = ",".join(map(str, val))
        res = f"{key}:{value_type}:{elemtype},{val}"
    elif value_type in "AfiZ":
        res = f"{key}:{value_type}:{val}"
    else:
        # hope that python string formatting is correct
        # TODO: warn about this?
        res = f"{key}:{value_type}:{val}"
    if ":S:" in res:
        raise ValueError(key, val, value_type)
    return res


def downgrade_mm_tag(align: AlignedSegment) -> AlignedSegment:
    orig_tag = None
    for tag in ["MM", "Mm"]:
        try:
            orig_tag = align.get_tag(tag)
        except KeyError:
            continue
        if orig_tag:
            break
    # alignment doesn't have tag, return
    if not orig_tag:
        return align
    m = MM_TAG_SKIP_SCHEME_RE.search(orig_tag)
    # tag does't contain skip_scheme character
    if not m:
        return align
    new_tag = f"{tag}:Z:" + MM_TAG_SKIP_SCHEME_RE.sub("", orig_tag)
    d = align.to_dict()
    d["tags"] = [t for t in d["tags"] if not t.startswith(tag)] + [new_tag]
    return AlignedSegment.from_dict(d, header=align.header)


class Cutter(metaclass=ABCMeta):
    @abstractmethod
    def get_cut_sites(self, seq: str) -> List[int]:
        ...

    def get_cut_intervals(self, seq: str) -> Tuple[List[int], List[int]]:
        sites = self.get_cut_sites(seq)
        seq_len = len(seq)
        if not sites:
            return [0], [seq_len]
        if sites[-1] == seq_len:
            sites = sites[:-1]
        starts = sites
        if starts[0] != 0:
            starts.insert(0, 0)
        ends = starts[1:] + [len(seq)]
        return starts, ends


class EnzymeCutter(Cutter):
    def __init__(self, enzyme: Any):
        self.enzyme = enzyme

    @classmethod
    def from_name(cls, enzyme_id: str):
        from Bio import Restriction

        enz = getattr(Restriction, enzyme_id, None)
        if enz is None:
            raise ValueError(f"Enzyme not found: {enzyme_id}")
        if enz.cut_twice():
            raise NotImplementedError(
                f"Enzyme cuts twice, not currently supported: {enzyme_id}"
            )
        return cls(enz)

    def get_cut_sites(self, seq: str) -> List[int]:
        return [_ - 1 for _ in self.enzyme.search(Seq(seq))]


def get_subread(
    sequence: str,
    quality: Optional[str],
    start: int,
    end: int,
    modified_bases: Optional[Mapping] = None,
) -> Tuple[str, Optional[str], Optional[str], Optional[str]]:
    seq = sequence[start:end]
    if quality:
        qual = quality[start:end]
    else:
        qual = None
    if modified_bases:
        mm_str, ml_str = "MM:Z:", "ML:B:"
        base_indices = {}
        for mod_key, mod_data in modified_bases.items():
            # find the modifications that overlap the subread
            idx = [x for x in range(len(mod_data)) if start <= mod_data[x][0] < end]
            if not idx:  # no mods in this subread
                continue
            try:
                canonical_base, strand, skip_scheme, mod_type = mod_key
            except ValueError:
                canonical_base, strand, mod_type = mod_key
                skip_scheme = ""

            canonical_base = mod_key[0]
            # find the posiitions of that canonical base in the subread
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


@define(kw_only=True)
class AlignInfo:
    ref_name: str = "*"
    ref_pos: int = 0
    flag: int = 4
    map_quality: int = 0
    cigar: str = "*"
    length: int = 0


@define(kw_only=True)
class ReadSeq:
    name: str
    sequence: str
    quality: Optional[str] = None
    mod_bases: Optional[Dict] = None
    align_info: Optional[AlignInfo] = None
    tags: Dict[str, str] = Factory(dict)

    @property
    def tag_str(self):
        return "\t".join(map(str, self.tags.values()))

    @classmethod
    def from_fastq(cls, rec: FastxRecord):
        if rec.comment:
            tags = {
                item.split(":", 1)[0]: item
                for item in rec.comment.split()
                if FASTQ_TAG_RE.match(item.strip())
            }
        else:
            tags = {}
        return ReadSeq(
            name=rec.name, sequence=rec.sequence, quality=rec.quality, tags=tags
        )

    def to_fastq_str(self) -> str:
        assert self.quality is not None
        return f"@{self.name}\t{self.tag_str}\n{self.sequence}\n+\n{self.quality}\n"

    def to_fastq(self) -> FastxRecord:
        return FastxRecord(
            name=self.name,
            sequence=self.sequence,
            quality=self.quality,
            comment=self.tag_str,
        )

    @classmethod
    def from_align(
        cls,
        rec: AlignedSegment,
        as_unaligned: bool = False,
        init_mod_bases: bool = False,
    ):
        tags = {}
        for (tag, tag_data, tag_type) in rec.get_tags(with_value_type=True):
            tags[tag] = tag_tuple_to_str(tag, tag_data, tag_type)
        if as_unaligned or rec.is_unmapped:
            align_info = None
        else:
            align_info = AlignInfo(
                ref_name=rec.reference_name,
                ref_pos=rec.reference_start,
                flag=rec.flag,
                map_quality=rec.mapping_quality,
                cigar=rec.cigarstring,
                length=rec.template_length,
            )
        if init_mod_bases:
            mod_bases = rec.modified_bases
            if not mod_bases:
                mod_bases = None
        else:
            mod_bases = None
        return cls(
            name=rec.query_name,
            sequence=rec.query_sequence,
            quality=rec.qual,
            tags=tags,
            mod_bases=mod_bases,
            align_info=align_info,
        )

    def to_align(
        self, header=DEFAULT_ALIGN_HEADER, as_unaligned: bool = False
    ) -> AlignedSegment:
        tags = self.tag_str
        if as_unaligned or self.align_info is None:
            flag = 4
            sam_str = (
                # QNAME,FLAG,RNAME,POS
                f"{self.name}\t{flag}\t*\t0\t"
                # MAPQ,CIGAR
                f"0\t*\t"
                # RNEXT,PNEXT,TLEN,
                "*\t0\t0\t"
                # SEQ, QUAL
                f"{self.sequence}\t{self.quality}\t"
                # Optional fields
                f"{tags}"
            )
        else:
            sam_str = (
                # QNAME,FLAG,RNAME,POS
                f"{self.name}\t{self.align_info.flag}\t"
                f"{self.align_info.ref_name}\t{self.align_info.ref_pos + 1}\t"
                # MAPQ,CIGAR
                f"{self.align_info.map_quality}\t{self.align_info.cigar}\t"
                # RNEXT,PNEXT,TLEN,
                f"*\t0\t{self.align_info.length}\t"
                # SEQ, QUAL
                f"{self.sequence}\t{self.quality}\t"
                # Optional fields
                f"{tags}"
            )
        return AlignedSegment.fromstring(sam_str, header=header)


@define(kw_only=True, frozen=True)
class ConcatemerCoords:
    """Relative position of monomer on concatemer"""

    start: int
    # start coordinate on the concatemer
    end: int
    # end coordinate on the concatemer
    subread_idx: int
    # the index of the monomer on the read
    subread_total: int
    # total number of monomers in the read

    @classmethod
    def from_tag(cls, tag: str):
        try:
            tag, _, tag_data = tag.split(":", 2)
            _, coords = tag_data.split(",", 1)
            start, end, subread_idx, subread_total = map(int, coords.split(","))
        except Exception:
            raise ValueError(f"Error parsing concatmer coords from {tag}")
        return cls(
            start=start, end=end, subread_idx=subread_idx, subread_total=subread_total
        )

    def to_tag(self) -> str:
        return f"Xc:B:i,{self.start},{self.end},{self.subread_idx},{self.subread_total}"


@define(kw_only=True)
class MonomerReadSeq:
    monomer_id: str
    # unique id for the monomer. Format {read_name}:{read_start}:{read_end}
    concatemer_id: str
    # the id of the read this monomer is derived from
    read_seq: ReadSeq
    # a ReadSeq consisting of just the sequence/quality/mod_bases of the subread
    coords: ConcatemerCoords
    # the coordinates of the monomer within the concatemer

    def _update_tags(self):
        self.read_seq.tags["MI"] = f"MI:Z:{self.concatemer_id}"
        self.read_seq.tags["Xc"] = self.coords.to_tag()

    @classmethod
    def from_readseq(cls, rec: ReadSeq):
        if "MI" not in rec.tags:
            raise ValueError(f"No MI tag in {rec.name}, can't create monomer")
        if "Xc" not in rec.tags:
            raise ValueError(
                f"No Xc tags in {rec.name}, can't reconstruct monomer coords"
            )
        coords = ConcatemerCoords.from_tag(rec.tags["Xc"])
        monomer_id = rec.name
        concatemer_id = rec.tags["MI"].rsplit(":", 1)[-1]
        return cls(
            monomer_id=monomer_id,
            concatemer_id=concatemer_id,
            read_seq=rec,
            coords=coords,
        )

    @classmethod
    def from_fastq(cls, rec: FastxRecord):
        return cls.from_readseq(ReadSeq.from_fastq(rec))

    @classmethod
    def from_align(
        cls,
        rec: AlignedSegment,
        as_unaligned: bool = False,
        init_mod_bases: bool = False,
    ):
        return cls.from_readseq(
            ReadSeq.from_align(
                rec, as_unaligned=as_unaligned, init_mod_bases=init_mod_bases
            )
        )

    def to_fastq(self) -> FastxRecord:
        self._update_tags()
        return self.read_seq.to_fastq()

    def to_align(
        self, header=DEFAULT_ALIGN_HEADER, as_unaligned: bool = False
    ) -> AlignedSegment:
        self._update_tags()
        return self.read_seq.to_align(header=header, as_unaligned=as_unaligned)


def splits_to_intervals(positions: List[int], length: int) -> List[Tuple[int, int]]:
    if len(positions) == 0:
        return [(0, length)]
    prefix, suffix = [], []
    if positions[0] != 0:
        prefix = [0]
    if positions[-1] != length:
        suffix = [length]
    breaks = prefix + positions + suffix
    return [(start, end) for start, end in zip(breaks[:-1], breaks[1:])]


@define(kw_only=True)
class ConcatemerReadSeq:
    concatemer_id: str
    read_seq: ReadSeq

    @classmethod
    def from_readseq(cls, rec: ReadSeq):
        return cls(concatemer_id=rec.name, read_seq=rec)

    @classmethod
    def from_fastq(cls, rec: FastxRecord):
        read_seq = ReadSeq.from_fastq(rec)
        return cls.from_readseq(read_seq)

    @classmethod
    def from_align(
        cls,
        rec: AlignedSegment,
        as_unaligned: bool = False,
        init_mod_bases: bool = False,
    ):
        read_seq = ReadSeq.from_align(
            rec, as_unaligned=as_unaligned, init_mod_bases=init_mod_bases
        )
        return cls.from_readseq(read_seq)

    def cut(self, cutter: Cutter) -> List[MonomerReadSeq]:
        positions = cutter.get_cut_sites(self.read_seq.sequence)
        return self.split(positions)

    def split(self, positions: List[int]) -> List[MonomerReadSeq]:
        intervals = splits_to_intervals(positions, len(self.read_seq.sequence))
        num_intervals = len(intervals)
        res = []
        for idx, (start, end) in enumerate(intervals):
            monomer_id = f"{self.concatemer_id}:{start}:{end}"
            coords = ConcatemerCoords(
                start=start,
                end=end,
                subread_idx=idx,
                subread_total=num_intervals,
            )
            seq, qual, mm_str, ml_str = get_subread(
                sequence=self.read_seq.sequence,
                quality=self.read_seq.quality,
                start=start,
                end=end,
                modified_bases=self.read_seq.mod_bases,
            )
            tags = {k: v for k, v in self.read_seq.tags.items() if k not in MOD_TAGS}
            tags["Xc"] = coords.to_tag()
            tags["MI"] = f"MI:Z:{self.concatemer_id}"
            if mm_str and ml_str:
                tags["MM"] = mm_str
                tags["ML"] = ml_str
            elif mm_str or ml_str:
                raise ValueError(mm_str, ml_str)

            _read_seq = ReadSeq(name=monomer_id, sequence=seq, quality=qual, tags=tags)
            res.append(
                MonomerReadSeq(
                    concatemer_id=self.concatemer_id,
                    monomer_id=monomer_id,
                    read_seq=_read_seq,
                    coords=coords,
                )
            )
        return res
