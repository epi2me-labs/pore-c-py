import re
from typing import List, Mapping, Optional, Tuple

from attrs import Factory, asdict, define
from pysam import AlignedSegment, AlignmentHeader, FastxRecord

from .settings import DEFAULT_ALIGN_HEADER, FASTQ_TAG_RE, MOD_TAGS

TAG_MI_RE = re.compile(r"MI\:Z\:(\S+)")
TAG_XC_RE = re.compile(
    r"Xc:B:i,(?P<start>\d+),(?P<end>\d+),(?P<subread_idx>\d+),(?P<subread_total>\d+)"
)


@define(kw_only=True, frozen=True)
class ConcatemerMetaData:
    concatemer: str
    # The read_id of the concatemer
    start: int
    # start coordinate on the read
    end: int
    # end coordinate on the read
    subread_idx: int
    # the index of the monomer on the read
    subread_total: int
    # total number of monomers in the read

    @classmethod
    def from_tags(cls, data: List[str]) -> "ConcatemerMetaData":
        d = {_.split(":", 1)[0]: _ for _ in data}
        if "MI" not in d:
            raise ValueError(f"No concatemer id  in tags: {data}")
        m = TAG_MI_RE.match(d["MI"].strip())
        if m:
            concatemer = m.group(1)
        else:
            raise ValueError(f"Error parsing concatemer id: {d['MI']}")
        if "Xc" not in d:
            raise ValueError(f"No concatemer coords  in tags: {data}")
        m = TAG_XC_RE.match(d["Xc"].strip())
        if m:
            coord_data = m.groupdict()
            return cls(
                concatemer=concatemer,
                start=int(coord_data["start"]),
                end=int(coord_data["end"]),
                subread_idx=int(coord_data["subread_idx"]),
                subread_total=int(coord_data["subread_total"]),
            )
        else:
            raise ValueError(f"Error parsing concatemer coords: {d['Xc']}")


@define(kw_only=True, frozen=True)
class AlignData:
    """Utility class for building pysam AlignedSegments"""

    name: str
    seq: str
    flag: int = 4  # unaligned
    ref_name: str = "*"
    ref_pos: int = 0
    map_quality: int = 0
    cigar: str = "*"
    next_ref_name: str = "*"
    next_ref_pos: int = 0
    length: int = 0
    qual: str = "*"
    tags: List[str] = Factory(list)

    @property
    def concatemer_metadata(self) -> Optional[ConcatemerMetaData]:
        return ConcatemerMetaData.from_tags(self.tags)

    @property
    def has_mods(self) -> bool:
        for t in self.tags:
            if t.split(":", 1)[0] in MOD_TAGS:
                return True
        return False

    def split(self, positions: List[int]) -> List["AlignData"]:
        return split_align_data(self, positions)

    def get_subread(
        self,
        start: int,
        end: int,
        name: str,
        add_mi_tag: bool = True,
        add_xc_tag: bool = True,
        subread_index: Optional[Tuple[int, int]] = None,
        modified_bases: Optional[Mapping] = None,
    ):
        return get_subread(
            self,
            start=start,
            end=end,
            name=name,
            add_mi_tag=add_mi_tag,
            add_xc_tag=add_xc_tag,
            subread_index=subread_index,
            modified_bases=modified_bases,
        )

    def to_pysam(
        self, header: AlignmentHeader = DEFAULT_ALIGN_HEADER
    ) -> AlignedSegment:
        return AlignedSegment.from_dict(
            {k: str(v) if k != "tags" else v for k, v in asdict(self).items()}, header
        )

    @classmethod
    def from_fastq(cls, record: FastxRecord) -> "AlignData":
        if record.comment:
            tags = [
                item
                for item in record.comment.split()
                if FASTQ_TAG_RE.match(item.strip())
            ]
        else:
            tags = []
        return cls(
            name=record.name, seq=record.sequence, qual=record.quality, tags=tags
        )

    @classmethod
    def from_alignment(cls, align: AlignedSegment) -> "AlignData":
        d = align.to_dict()
        for k in ["flag", "map_quality", "ref_pos", "next_ref_pos", "length"]:
            d[k] = int(d[k])
        return cls(**d)

    def to_fastq(self, with_tags: bool = True) -> str:
        if with_tags:
            tag_str = "\t".join(self.tags)
        else:
            tag_str = ""
        if self.qual == "*":
            raise ValueError(f"No quality, can't write fastq for {self.name}")
        return f"@{self.name} {tag_str}\n{self.seq}\n+\n{self.qual}\n"


def split_align_data(record: AlignData, positions: List[int]) -> List["AlignData"]:
    read_length = len(record.seq)
    if not positions:  # no cut sites, return whole read
        positions = [0, read_length]
    else:
        positions = sorted(positions)
        if positions[0] != 0:
            positions = [0] + positions
        if positions[-1] != read_length:
            positions = positions + [read_length]
    intervals = list(zip(positions[:-1], positions[1:]))
    num_subreads = len(intervals)
    if record.has_mods:
        modified_bases = record.to_pysam().modified_bases  # pyright: ignore
    else:
        modified_bases = None
    subreads = []
    for x, (start, end) in enumerate(intervals):
        name = f"{record.name}:{x+1}_{num_subreads}"
        subreads.append(
            record.get_subread(
                start,
                end,
                name,
                modified_bases=modified_bases,
                subread_index=(x, num_subreads),
            )
        )
    return subreads


def get_subread(
    record,
    start: int,
    end: int,
    name: str,
    add_mi_tag: bool = True,
    add_xc_tag: bool = True,
    subread_index: Optional[Tuple[int, int]] = None,
    modified_bases: Optional[Mapping] = None,
):
    seq = record.seq[start:end]
    qual = record.qual[start:end]
    tags = []
    for t in record.tags:
        key = t.split(":", 1)[0]
        if key in MOD_TAGS:
            continue
        if key == "MI" and add_mi_tag:
            continue
        tags.append(t)
    if add_mi_tag:
        tags.append(f"MI:Z:{record.name}")
    if add_xc_tag and subread_index:
        tags.append(f"Xc:B:i,{start},{end},{subread_index[0]},{subread_index[1]}")

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
        tags.append(mm_str)
        tags.append(ml_str)
    return AlignData(name=name, seq=seq, qual=qual, tags=tags)
