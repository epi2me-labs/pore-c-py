"""Model."""
from abc import ABCMeta, abstractmethod
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any, Dict, List, Literal, Mapping, Optional, Tuple

from Bio.Seq import Seq
from pysam import AlignedSegment, FastxRecord

from pore_c_py.sam_utils import (
    CONCATEMER_TAG,
    FASTQ_TAG_RE,
    MOD_TAGS,
    MOLECULE_TAG,
    SamFlags,
    tag_tuple_to_str,
    WALK_SEGMENT_RE,
    WALK_TAG
)
from pore_c_py.settings import DEFAULT_ALIGN_HEADER


class Cutter(metaclass=ABCMeta):
    """Cutter."""

    @abstractmethod
    def get_cut_sites(self, seq: str) -> List[int]:
        """Get cut sites."""
        ...

    def get_cut_intervals(self, seq: str) -> Tuple[List[int], List[int]]:
        """Get cut intervals."""
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
    """Enzyme cutter."""

    def __init__(self, enzyme: Any):
        """Init."""
        self.enzyme = enzyme

    @classmethod
    def from_name(cls, enzyme_id: str):
        """From name."""
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
        """Get cut sites."""
        return [_ - 1 for _ in self.enzyme.search(Seq(seq))]


def get_subread(
    sequence: str,
    quality: Optional[str],
    start: int,
    end: int,
    modified_bases: Optional[Mapping] = None,
) -> Tuple[str, Optional[str], Optional[str], Optional[str]]:
    """Get sub read."""
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


@dataclass()
class AlignInfo:
    """Align Info."""

    ref_name: str = "*"
    ref_pos: int = 0
    flag: int = 4
    map_quality: int = 0
    cigar: str = "*"
    template_length: int = 0
    length: int = 0

    @property
    def ref_end(self) -> int:
        """Ref end."""
        return self.ref_pos + self.length

    @property
    def strand(self) -> Literal["+", "-", "."]:
        """Strand."""
        return SamFlags.int_to_strand(self.flag)

    @cached_property
    def sam_flags(self) -> SamFlags:
        """Sam flags."""
        return SamFlags.from_int(self.flag)


@dataclass()
class ReadSeq:
    """Read seq."""

    name: str
    sequence: str
    flags: SamFlags = field(default_factory=lambda: SamFlags(unmap=True))
    quality: Optional[str] = None
    mod_bases: Optional[Dict] = None
    align_info: Optional[AlignInfo] = None
    next_align: Optional[Tuple[str, int]] = None
    tags: Dict[str, str] = field(default_factory=dict)

    @property
    def tag_str(self):
        """Tag str."""
        return "\t".join(map(str, self.tags.values()))

    @classmethod
    def from_fastq(
            cls, rec: FastxRecord, remove_tags: Optional[List[str]] = None):
        """From fastq."""
        if rec.comment:
            tags = {
                item.split(":", 1)[0]: item
                for item in rec.comment.split()
                if FASTQ_TAG_RE.match(item.strip())
            }
            if remove_tags is not None:
                for t in remove_tags:
                    if t in tags:
                        tags.pop(t)
        else:
            tags = {}
        return ReadSeq(
            name=rec.name, sequence=rec.sequence,
            quality=rec.quality, tags=tags
        )

    def to_fastq_str(self) -> str:
        """To fastq str."""
        assert self.quality is not None
        return f"@{self.name}\t{self.tag_str}\n{self.sequence}\n+\n{self.quality}\n" # noqa

    def to_fastq(self) -> FastxRecord:
        """To fastq."""
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
        init_mod_bases: bool = True,
        remove_tags: Optional[List[str]] = None,
    ):
        """From align."""
        tags = {}
        for (tag, tag_data, tag_type) in rec.get_tags(with_value_type=True):
            if remove_tags and tag in remove_tags:
                continue
            tags[tag] = tag_tuple_to_str(tag, tag_data, tag_type)
        if as_unaligned or rec.is_unmapped:
            align_info = None
            flags = SamFlags(unmap=True)
        else:
            flags = SamFlags.from_int(rec.flag)
            align_info = AlignInfo(
                ref_name=rec.reference_name,  # type: ignore
                ref_pos=rec.reference_start,
                flag=rec.flag,
                map_quality=rec.mapping_quality,
                cigar=rec.cigarstring,  # type: ignore
                template_length=rec.template_length,
                length=rec.reference_end - rec.reference_start,  # type: ignore
            )
        if init_mod_bases:
            mod_bases = rec.modified_bases  # type: ignore
            if not mod_bases:
                mod_bases = None
        else:
            mod_bases = None

        if rec.next_reference_name:
            next_align = (
                rec.next_reference_name, int(rec.next_reference_start))
        else:
            next_align = None
        return cls(
            name=rec.query_name,  # type: ignore
            sequence=rec.query_sequence,  # type: ignore
            quality=rec.qual,  # type: ignore
            flags=flags,
            tags=tags,
            mod_bases=mod_bases,
            align_info=align_info,
            next_align=next_align,
        )

    def to_align(
        self,
        header=DEFAULT_ALIGN_HEADER,
        **kwds,
    ) -> AlignedSegment:
        """To align."""
        return AlignedSegment.fromstring(self.to_sam(**kwds), header=header)

    #  TODO: create a dataclass for the to_sam args
    def to_sam(
        self,
        as_unaligned: bool = False,
        read_name: Optional[str] = None,
        flag: Optional[SamFlags] = None,
        strip_tags: bool = False,
        next_reference_name: str = "*",
        next_reference_start: int = 0,
        template_length: Optional[int] = None,
    ) -> str:
        """To sam."""
        if read_name is None:
            read_name = self.name
        if strip_tags:
            tags = ""
        else:
            tags = self.tag_str.strip()

        is_unmapped = as_unaligned or self.align_info is None

        if flag is None:
            flag = self.flags.copy()  # no side-effects

        if is_unmapped and not flag.unmap:
            # TODO: should we check other flags here too
            flag.unmap = True

        if is_unmapped or flag.paired is False:
            template_length = 0
        elif template_length is None:
            template_length = self.align_info.template_length  # type: ignore
        if is_unmapped:
            sam_str = (
                # QNAME,FLAG,RNAME,POS
                f"{read_name}\t{flag.to_int()}\t*\t0\t"
                # MAPQ,CIGAR
                f"0\t*\t"
                # RNEXT,PNEXT,TLEN,
                f"{next_reference_name}\t{next_reference_start}\t{template_length}\t" # noqa
                # SEQ, QUAL
                f"{self.sequence}\t{self.quality}\t"
                # Optional fields
                f"{tags}"
            )
        else:
            assert self.align_info is not None
            sam_str = (
                # QNAME,FLAG,RNAME,POS
                f"{read_name}\t{flag.to_int()}\t"
                f"{self.align_info.ref_name}\t"
                f"{self.align_info.ref_pos + 1}\t"
                # MAPQ,CIGAR
                f"{self.align_info.map_quality}\t"
                f"{self.align_info.cigar}\t"
                # RNEXT,PNEXT
                f"{next_reference_name}\t{next_reference_start}\t"
                # TLEN
                f"{template_length}\t"
                # SEQ, QUAL
                f"{self.sequence}\t{self.quality}\t"
                # Optional fields
                f"{tags}"
            )
        return sam_str


@dataclass()
class ConcatemerCoords:
    """Relative position of monomer on concatemer."""

    start: int
    # start coordinate on the concatemer
    end: int
    # end coordinate on the concatemer
    read_length: int
    # total concatemer length
    subread_idx: int
    # the index of the monomer on the read
    subread_total: int
    # total number of monomers in the read

    @classmethod
    def from_tag(cls, tag: str):
        """From tag."""
        try:
            tag, _, tag_data = tag.split(":", 2)
            _, coords = tag_data.split(",", 1)
            start, end, read_length, subread_idx, subread_total = map(
                int, coords.split(",")
            )
        except Exception:
            raise ValueError(f"Error parsing concatmer coords from {tag}")
        return cls(
            start=start,
            end=end,
            read_length=read_length,
            subread_idx=subread_idx,
            subread_total=subread_total,
        )

    def to_tag(self) -> str:
        """To tag."""
        return (
            f"{CONCATEMER_TAG}:B:i,{self.start},{self.end},{self.read_length},"
            f"{self.subread_idx},{self.subread_total}"
        )


@dataclass()
class MonomerReadSeq:
    """Monomer read seq."""

    monomer_id: str
    # unique id for the monomer. Format {read_name}:{read_start}:{read_end}
    concatemer_id: str
    # the id of the read this monomer is derived from
    read_seq: ReadSeq
    # a ReadSeq consisting of just the sequence/quality/mod_bases
    # of the subread
    coords: ConcatemerCoords
    # the coordinates of the monomer within the concatemer

    @property
    def is_aligned(self) -> bool:
        """Is aligned."""
        return self.read_seq.align_info is not None

    @staticmethod
    def generate_id(concatemer_id: str, coords: ConcatemerCoords) -> str:
        """Create a unique, lexographically sortable monomer_id.

        It's useful to have a monomer id that when lexographically sorted
        recapitulates the order of the monomers within the original read.
        To do this we append zero-padded interval of the monomer within
        the concatemer: <read_id>:<zero-padded start>:<zero-padded end>
        """
        num_digits = len(str(coords.read_length))
        return (
            f"{concatemer_id}:"
            f"{coords.start:0{num_digits}d}:{coords.end:0{num_digits}d}"
        )

    def _update_tags(self):
        """Update tags."""
        self.read_seq.tags[MOLECULE_TAG] = \
            f"{MOLECULE_TAG}:Z:{self.concatemer_id}"
        self.read_seq.tags[CONCATEMER_TAG] = self.coords.to_tag()

    @classmethod
    def from_readseq(cls, rec: ReadSeq):
        """From readseq."""
        if MOLECULE_TAG not in rec.tags:
            raise ValueError(
                f"No {MOLECULE_TAG} tag in {rec.name}, can't create monomer"
            )
        if CONCATEMER_TAG not in rec.tags:
            raise ValueError(
                f"No {CONCATEMER_TAG} tags in {rec.name}, "
                "can't reconstruct monomer coords"
            )
        coords = ConcatemerCoords.from_tag(rec.tags[CONCATEMER_TAG])
        monomer_id = rec.name
        concatemer_id = rec.tags[MOLECULE_TAG].rsplit(":", 1)[-1]
        return cls(
            monomer_id=monomer_id,
            concatemer_id=concatemer_id,
            read_seq=rec,
            coords=coords,
        )

    @classmethod
    def from_fastq(
            cls, rec: FastxRecord, remove_tags: Optional[List[str]] = None):
        """From fastq."""
        return cls.from_readseq(ReadSeq.from_fastq(
            rec, remove_tags=remove_tags))

    @classmethod
    def from_align(
        cls,
        rec: AlignedSegment,
        as_unaligned: bool = False,
        init_mod_bases: bool = False,
        remove_tags: Optional[List[str]] = None,
    ):
        """From align."""
        return cls.from_readseq(
            ReadSeq.from_align(
                rec,
                as_unaligned=as_unaligned,
                init_mod_bases=init_mod_bases,
                remove_tags=remove_tags,
            )
        )

    def to_fastq(self) -> FastxRecord:
        """To fastq."""
        self._update_tags()
        return self.read_seq.to_fastq()

    def to_fastq_str(self) -> str:
        """To fastq str."""
        self._update_tags()
        return self.read_seq.to_fastq_str()

    def to_align(
        self, header=DEFAULT_ALIGN_HEADER, as_unaligned: bool = False
    ) -> AlignedSegment:
        """To align."""
        self._update_tags()
        return self.read_seq.to_align(header=header, as_unaligned=as_unaligned)


def splits_to_intervals(
        positions: List[int], length: int) -> List[Tuple[int, int]]:
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


@dataclass()
class ConcatemerReadSeq:
    """Concatemer read seq."""

    concatemer_id: str
    read_seq: ReadSeq

    @classmethod
    def from_readseq(cls, rec: ReadSeq):
        """From read seq."""
        return cls(concatemer_id=rec.name, read_seq=rec)

    @classmethod
    def from_fastq(
            cls, rec: FastxRecord, remove_tags: Optional[List[str]] = None):
        """From fastq."""
        read_seq = ReadSeq.from_fastq(rec, remove_tags=remove_tags)
        return cls.from_readseq(read_seq)

    @classmethod
    def from_align(
        cls,
        rec: AlignedSegment,
        as_unaligned: bool = False,
        init_mod_bases: bool = False,
        remove_tags: Optional[List[str]] = None,
    ):
        """From align."""
        read_seq = ReadSeq.from_align(
            rec,
            as_unaligned=as_unaligned,
            init_mod_bases=init_mod_bases,
            remove_tags=remove_tags,
        )
        return cls.from_readseq(read_seq)

    def to_fastq_str(self, walk: Optional["Walk"]):
        """To fastq str."""
        if walk:
            self.read_seq.tags[WALK_TAG] = walk.to_tag()
        return self.read_seq.to_fastq_str()

    def cut(self, cutter: Cutter) -> List[MonomerReadSeq]:
        """Cut."""
        positions = cutter.get_cut_sites(self.read_seq.sequence)
        return self.split(positions)

    def split(self, positions: List[int]) -> List[MonomerReadSeq]:
        """Split."""
        read_length = len(self.read_seq.sequence)
        intervals = splits_to_intervals(positions, read_length)
        num_intervals = len(intervals)
        res = []
        for idx, (start, end) in enumerate(intervals):
            coords = ConcatemerCoords(
                start=start,
                end=end,
                read_length=read_length,
                subread_idx=idx,
                subread_total=num_intervals,
            )
            monomer_id = MonomerReadSeq.generate_id(self.concatemer_id, coords)
            seq, qual, mm_str, ml_str = get_subread(
                sequence=self.read_seq.sequence,
                quality=self.read_seq.quality,
                start=start,
                end=end,
                modified_bases=self.read_seq.mod_bases,
            )
            tags = {
                k: v for k, v in self.read_seq.tags.items()
                if k not in MOD_TAGS}
            tags[CONCATEMER_TAG] = coords.to_tag()
            tags[MOLECULE_TAG] = f"{MOLECULE_TAG}:Z:{self.concatemer_id}"
            if mm_str and ml_str:
                tags["MM"] = mm_str
                tags["ML"] = ml_str
            elif mm_str or ml_str:
                raise ValueError(mm_str, ml_str)

            _read_seq = ReadSeq(
                name=monomer_id, sequence=seq, quality=qual, tags=tags)
            res.append(
                MonomerReadSeq(
                    concatemer_id=self.concatemer_id,
                    monomer_id=monomer_id,
                    read_seq=_read_seq,
                    coords=coords,
                )
            )
        if len(res) != num_intervals:
            raise ValueError(intervals)
        return res


@dataclass
class WalkSegment:
    """Walk segment."""

    read_start: int
    read_end: int
    chrom: Optional[str] = None
    genome_start: Optional[int] = None
    genome_end: Optional[int] = None
    orientation: Optional[Literal["+", "-", "."]] = None

    @classmethod
    def from_monomer(cls, monomer: MonomerReadSeq):
        """From monomer."""
        if monomer.read_seq.align_info is not None:
            return cls(
                read_start=monomer.coords.start,
                read_end=monomer.coords.end,
                chrom=monomer.read_seq.align_info.ref_name,
                genome_start=monomer.read_seq.align_info.ref_pos,
                genome_end=monomer.read_seq.align_info.ref_end,
                orientation=monomer.read_seq.align_info.strand,
            )
        else:
            return cls(
                read_start=monomer.coords.start,
                read_end=monomer.coords.end,
            )

    @classmethod
    def from_string(cls, val: str):
        """From string."""
        if val.startswith("*"):
            try:
                _, coords = val.split(":")
                read_start, read_end = map(int, coords.split("-"))
                return cls(read_start, read_end)
            except Exception:
                raise ValueError(
                    f"Error parsing unaligned walk segment from {val}")
        else:
            m = WALK_SEGMENT_RE.match(val)
            if not m:
                raise ValueError(
                    f"Error parsing aligned walk segment from {val}")
            else:
                _ = m.groupdict()
                return cls(
                    read_start=int(_["read_start"]),
                    read_end=int(_["read_end"]),
                    chrom=_["chrom"],
                    genome_start=int(_["genome_start"]),
                    genome_end=int(_["genome_end"]),
                    orientation=_["orientation"],  # type: ignore
                )

    def to_string(self) -> str:
        """To string."""
        if self.chrom is None:
            assert (
                all(
                    [
                        _ is None
                        for _ in (
                            self.genome_start, self.genome_end,
                            self.orientation)
                    ]
                )
                is True
            )
            return f"*:{self.read_start}-{self.read_end}"
        elif ";" in self.chrom:
            raise ValueError(
                f"Chromosome name must not include ';' character {self.chrom}"
            )
        else:
            assert (
                any(
                    [
                        _ is None
                        for _ in (
                            self.genome_start, self.genome_end,
                            self.orientation)
                    ]
                )
                is False
            )
            return (
                f"{self.chrom}:{self.orientation}:"
                f"{self.genome_start}-{self.genome_end}:"
                f"{self.read_start}-{self.read_end}"
            )


@dataclass
class Walk:
    """Walk."""

    segments: List[WalkSegment]

    @classmethod
    def from_aligns(cls, aligns: List[MonomerReadSeq]):
        """From aligns."""
        return cls([WalkSegment.from_monomer(_) for _ in aligns])

    def __len__(self):
        """Length."""
        return len(self.segments)

    @classmethod
    def from_tag(cls, tag: str):
        """From tag."""
        # TODO should we check the tag here?
        _, _, tag_data = tag.split(":", 2)
        segments = [WalkSegment.from_string(_) for _ in tag_data.split(";")]
        return cls(segments)

    def to_tag(self) -> str:
        """To tag."""
        return f"{WALK_TAG}:Z:" + ";".join(
            [_.to_string() for _ in self.segments])
