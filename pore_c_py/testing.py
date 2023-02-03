"""Testing."""
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from functools import cached_property
from itertools import count
import json
from pathlib import Path
from shutil import which
import subprocess as sp
from tempfile import mkdtemp
from typing import Any, Dict, Iterable, List, Literal, Mapping, Optional, Tuple

import numpy as np
from numpy.random import default_rng, Generator
import numpy.typing as npt
import polars as pl
from pysam import faidx
from pysam import (
    FastaFile, tabix_compress,
    tabix_index, VariantFile, VariantHeader
)

from pore_c_py.io import fastq_to_ubam, FileCollection
from pore_c_py.log import get_logger
from pore_c_py.model import (
    AlignInfo,
    ConcatemerCoords,
    ConcatemerReadSeq,
    MonomerReadSeq,
    ReadSeq,
    SamFlags,
    Walk,
    WalkSegment,
)
from pore_c_py.monomers import GenomicFragment
from pore_c_py.sam_utils import MOLECULE_TAG

logger = get_logger()


def simulate_walk(
    num_monomers: int,
    contact_probs: npt.NDArray[np.float_],
    random_state: Generator,
    fragments_df: pl.DataFrame,
    max_attempt_factor: float = 2.0,
    stranded: bool = False,
) -> Tuple[Walk, List[int]]:
    """Simulate walk."""
    if stranded:
        strands = random_state.choice([True, False], size=num_monomers)
    else:
        strands = [True] * num_monomers
    walk_fragments = simulate_walk_path(
        num_monomers, contact_probs,
        random_state, max_attempt_factor=max_attempt_factor
    )
    walk_genome_coords = fragments_df[
        pl.Series(walk_fragments), ["chrom", "start", "end"]
    ]
    concatemer_start = 0
    walk = Walk([])
    for _, ((chrom, start, end), strand) in enumerate(
        zip(walk_genome_coords.rows(), strands)
    ):
        length = end - start
        concatemer_end = concatemer_start + length
        walk.segments.append(
            WalkSegment(
                read_start=concatemer_start,
                read_end=concatemer_end,
                chrom=chrom,
                genome_start=start,
                genome_end=end,
                orientation="+" if strand else "-",
            )
        )
        concatemer_start = concatemer_end
    return walk, walk_fragments


def build_monomer(
    segment: WalkSegment,
    ff: FastaFile,
    read_length: int,
    subread_idx: int,
    subread_total: int,
    concatemer_id: str,
    snps: Optional[List["SNPHaplotypes"]] = None,
    haplotype: int = 0,
) -> MonomerReadSeq:
    """Build monomer."""
    coords = ConcatemerCoords(
        start=segment.read_start,
        end=segment.read_end,
        read_length=read_length,
        subread_idx=subread_idx,
        subread_total=subread_total,
    )
    length = segment.read_end - segment.read_start
    monomer_id = MonomerReadSeq.generate_id(concatemer_id, coords)
    if snps:
        offsets = []
        allele_values = []
        for s in snps:
            allele_idx = s.haplotypes[haplotype]
            if allele_idx == 0:
                allele_values.append(s.ref)
            else:
                allele_values.append(s.alt)
            offsets.append(s.offset)
    else:
        allele_values, offsets = None, None

    # make type checker happy
    assert segment.chrom is not None
    assert segment.genome_start is not None
    assert segment.genome_end is not None

    _seq, _, _ = simulate_read_sequence(
        ff=ff,
        chrom=segment.chrom,
        start=segment.genome_start,
        end=segment.genome_end,
        strand=True if segment.orientation == "+" else False,
        offsets=offsets,
        allele_values=allele_values,
    )
    read_seq = ReadSeq(
        name=monomer_id,
        sequence=_seq,
        quality="5" * length,
        align_info=AlignInfo(
            ref_name=segment.chrom,
            ref_pos=segment.genome_start,
            flag=SamFlags(reverse=segment.orientation == "-").to_int(),
            length=length,
        ),
    )
    return MonomerReadSeq(
        concatemer_id=concatemer_id,
        monomer_id=monomer_id,
        read_seq=read_seq,
        coords=coords,
    )


def simulate_walk_path(
    num_monomers: int,
    contact_probs: npt.NDArray[np.float_],
    random_state: Generator,
    max_attempt_factor: float = 2.0,
) -> List[int]:
    """Simulate walk path."""
    n_frags = contact_probs.shape[0]
    # pick a random fragement to start walk
    start_idx = random_state.integers(0, n_frags)
    path = [start_idx]
    n_iter = 0
    max_iter = int(num_monomers * max_attempt_factor)
    while (len(path) < num_monomers) and n_iter < max_iter:
        # choose next fragment in proportion to the column for that fragment
        next_frag = random_state.choice(n_frags, p=contact_probs[:, start_idx])
        if next_frag not in path:
            path.append(next_frag)
        start_idx = next_frag
        n_iter += 1
    return path


def random_concatemer_generator(
    *,
    reference_fasta: FastaFile,
    fragments_df: pl.DataFrame,
    contact_probs: npt.NDArray[np.float_],
    random_state: Generator,
    mean_frags_per_concatemer: int = 5,
    min_frags_per_concatemer: int = 1,
    max_frags_per_concatemer: int = 50,
    fragment_to_snps: Mapping[int, List["SNPHaplotypes"]],
    num_haplotypes: int = 0,
    max_concatemers: int = 0,
) -> Iterable[Tuple[ConcatemerReadSeq, List[MonomerReadSeq], Walk, int]]:
    """Random concatemer generator."""
    concatemer_idx = 0
    while True:
        concatemer_id = f"CONCAT{concatemer_idx}"
        num_monomers = random_state.poisson(
            mean_frags_per_concatemer, size=1).clip(
            min_frags_per_concatemer, max_frags_per_concatemer
        )[0]
        if num_haplotypes > 0:
            haplotype = (
                0  # TODO:reinstate random_state.choice(num_haplotypes, size=1)[0] # noqa
            )
        else:
            haplotype = 0
        walk, fragments = simulate_walk(
            num_monomers=num_monomers,
            contact_probs=contact_probs,
            random_state=random_state,
            fragments_df=fragments_df,
        )
        # actual walk may be shorter than requested
        num_monomers = len(walk.segments)
        read_length = sum([(_.read_end - _.read_start) for _ in walk.segments])
        monomers = [
            build_monomer(
                segment,
                reference_fasta,
                read_length,
                monomer_idx,
                num_monomers,
                concatemer_id,
                fragment_to_snps.get(fragments[monomer_idx], None),
                haplotype,
            )
            for monomer_idx, segment in enumerate(walk.segments)
        ]

        concatemer = ConcatemerReadSeq(
            concatemer_id=concatemer_id,
            read_seq=ReadSeq(
                name=concatemer_id,
                sequence="".join([m.read_seq.sequence for m in monomers]),
                quality="".join([m.read_seq.quality for m in monomers]),
            ),
        )
        yield (concatemer, monomers, walk, haplotype)
        concatemer_idx += 1
        if max_concatemers and concatemer_idx > max_concatemers:
            break


def simulate_sequence_with_cut_sites(
    enzyme: str = "EcoRI",
    seq_length: int = 1000,
    cut_sites: List[int] = [50, 300, 750],
    random_state: Optional[Generator] = None,
    correct_ends: bool = True,
) -> Tuple[List[int], str]:
    """Simulate sequence with cut sites."""
    from Bio import Restriction
    from Bio.Seq import Seq

    if random_state is None:
        random_state = default_rng()
    enz = getattr(Restriction, enzyme)
    offset = enz.fst5
    length = len(enz.site)
    logger.debug(
        f"Simulating sequence with \
            {len(cut_sites)} cut sites in {seq_length} bases"
    )
    seq = random_state.choice(list("ACGT"), size=seq_length)
    niter = 0
    # remove any recognition sites we've introduced
    while True:
        matches = enz.search(Seq("".join(seq)))
        if len(matches) == 0:
            break
        for pos in matches:
            start = pos - offset - 1
            end = start + length
            new_site = random_state.choice(list("ACGT"), size=length)
            seq[start:end] = new_site
        niter += 1
        if niter > 10:
            raise ValueError(seq, matches)
    # TODO: handle case were inserting a cut site actually creates two sites
    for p in cut_sites:
        start = p - offset
        seq[start: start + length] = list(enz.site)
    seq = "".join(seq)
    if correct_ends:
        cut_left = enz.site[:offset]
        if not seq.endswith(cut_left):
            seq = seq[: -len(cut_left)] + cut_left
        cut_right = enz.site[offset:]
        if not seq.startswith(cut_right):
            seq = cut_right + seq[len(cut_right):]
    return cut_sites, seq


def simulate_fasta_with_cut_sites(
    fasta: Path,
    seq_length: Dict[str, int],
    cut_sites: Dict[str, List[int]],
    enzyme: str = "EcoRI",
    random_state: Optional[Generator] = None,
) -> Path:
    """Simulate fasta with cut sites."""
    if random_state is None:
        random_state = default_rng()
    fh = fasta.open("w")
    for chrom, length in seq_length.items():
        positions = cut_sites[chrom]
        _, seq = simulate_sequence_with_cut_sites(
            enzyme, cut_sites=positions, seq_length=length,
            random_state=random_state
        )
        fh.write(f">{chrom}\n{seq}\n")
    fh.close()
    faidx(str(fasta))
    return fasta


@dataclass
class SNPHaplotypes:
    """SNP Haplotypes."""

    chrom: str
    pos: int
    ref: str
    alt: str
    haplotypes: List[Literal[0, 1]]
    fragment_id: Optional[int] = None
    offset: Optional[int] = None

    def as_dict(self) -> Dict[str, Any]:
        """As Dict."""
        d = {
            "chrom": self.chrom, "pos": self.pos,
            "ref": self.ref, "alt": self.alt}
        for x, h in enumerate(self.haplotypes):
            d[f"HT{x}"] = h
        return d

    def to_vcf_str(self) -> str:
        """As VCF Str."""
        ht_str = "|".join(map(str, self.haplotypes))
        return (
            f"{self.chrom}\t{self.pos}\t.\t{self.ref}\t"
            f"{self.alt}\t.\tPASS\t.\tGT:PS\t{ht_str}:0\n"
        )


def assign_snps_to_fragments(
    chrom: str,
    snp_positions: List[int],
    fragment_df: pl.DataFrame,
    fragment_end_buffer: int,
) -> Tuple[List[int], List[int], List[int]]:
    """Assign snps to fragments."""
    # match snp positions to fragments
    # remove snps within tol of fragment ends
    #  don't want to destory RE site
    df = (
        (
            pl.DataFrame({"pos": snp_positions}, orient="col").with_row_count(
                name="snp_index"
            )
        )
        .join_asof(
            fragment_df.filter(pl.col("chrom") == chrom)
            .select(["start", "end", "fragment_id"])
            .with_column(pl.col("start").alias("asof_end")),
            left_on="pos",
            right_on="asof_end",
        )
        .filter(pl.col("pos") >= pl.col("start"))
        .with_columns(
            [
                (pl.col("pos") - pl.col("start")).alias("offset"),
                (pl.col("end") - pl.col("pos")).alias("end_offset"),
            ]
        )
        .filter(
            (pl.col("offset") >= fragment_end_buffer)
            & (pl.col("end_offset") >= fragment_end_buffer)
        )
    )
    return (
        df["pos"].to_list(),  # type: ignore
        df["offset"].to_list(),  # type: ignore
        df["fragment_id"].to_list(),  # type: ignore
    )


def simulate_sequence_haplotypes(
    chrom: str,
    ff: FastaFile,
    var_density: float,
    random_state: Generator,
    num_haplotypes: int = 2,
    fragment_df: Optional[pl.DataFrame] = None,
    fragment_end_buffer: int = 0,
) -> List[SNPHaplotypes]:
    """Simulate sequence haplotypes."""
    mutations = {}
    bases = list("ATGC")
    for b in bases:
        mutations[b] = [_ for _ in bases if _ != "b"]

    length = ff.get_reference_length(chrom)
    num_vars = max(2, int(np.rint(length * var_density)))
    snp_positions = sorted(
        random_state.choice(int(length), size=num_vars, replace=False)
    )
    if fragment_df:
        snp_positions, fragment_offsets, fragment_ids = assign_snps_to_fragments( # noqa
            chrom, snp_positions, fragment_df, fragment_end_buffer
            )
        num_vars = len(snp_positions)
    else:
        fragment_offsets, fragment_ids = None, None
    # assume bi-allelic
    ref_allele = [ff.fetch(chrom, _, _ + 1) for _ in snp_positions]
    alt_allele = [random_state.choice(mutations[_], 1)[0] for _ in ref_allele]
    # matrix of snps x haplotypes, 0==ref allele; 1==alt allele
    haplotypes = random_state.choice([0, 1], size=(num_vars, num_haplotypes))

    return [
        SNPHaplotypes(
            chrom,
            snp_positions[x],
            ref_allele[x],
            alt_allele[x],
            list(haplotypes[x, :]),
            None if fragment_ids is None else fragment_ids[x],
            None if fragment_offsets is None else fragment_offsets[x],
        )
        for x in range(num_vars)
    ]


def simulate_haplotypes(
    reference_fasta: Path,
    num_haplotypes: int = 2,
    var_density: float = 0.01,
    fragment_df: Optional[pl.DataFrame] = None,
    random_state: Optional[Generator] = None,
) -> List[SNPHaplotypes]:
    """Simulate haplotypes."""
    if random_state is None:
        random_state = default_rng()
    try:
        ff = FastaFile(str(reference_fasta))
    except Exception:
        logger.error(f"Error opeing {reference_fasta}")
        raise
    snp_haplotypes = []
    for chrom in ff.references:
        snp_haplotypes.extend(
            simulate_sequence_haplotypes(
                chrom,
                ff,
                var_density=var_density,
                num_haplotypes=num_haplotypes,
                random_state=random_state,
                fragment_df=fragment_df,
            )
        )
    assert len(snp_haplotypes) > 0
    return snp_haplotypes


def haplotypes_to_df(snp_haplotypes: List[SNPHaplotypes]) -> pl.DataFrame:
    """Haplotypes to df."""
    return (
        pl.DataFrame([_.as_dict() for _ in snp_haplotypes], orient="row")
        .melt(
            id_vars=["chrom", "pos", "ref", "alt"],
            variable_name="haplotype_id",
            value_name="allele_idx",
        )
        .with_column(
            pl.when(pl.col("allele_idx") == 0)
            .then(pl.col("ref"))
            .otherwise(pl.col("alt"))
            .alias("allele_value")
        )
    )


def create_phased_vcf(
    snp_haplotypes: List[SNPHaplotypes], vcf: Path, reference_fasta: Path
):
    """Create phased vcf."""
    if vcf.suffixes[-1] == ".gz":
        uncomp_vcf = vcf.with_suffix("".join(vcf.suffixes[:-1]))
        comp_vcf = vcf
    else:
        comp_vcf = Path(str(vcf) + ".gz")
        uncomp_vcf = vcf

    ff = FastaFile(str(reference_fasta))
    header = VariantHeader()
    header.add_line(f"##reference=file:///{reference_fasta.resolve()}")
    for chrom in ff.references:
        header.add_line(
            f"##contig=<ID={chrom},length={ff.get_reference_length(chrom)}>"
        )
    header.add_line(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line(
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">')
    header.add_line(
        '##FORMAT=<ID=PQ,Number=1,Type=Integer,Description="Phasing Quality">'
    )
    header.add_sample("SAMPLE_01")
    vcf_out = VariantFile(str(uncomp_vcf), "w", header=header)
    vcf_out.close()
    with uncomp_vcf.open("a") as fh:
        for snp in snp_haplotypes:
            fh.write(snp.to_vcf_str())
    tabix_compress(str(uncomp_vcf), str(comp_vcf))
    tabix_index(str(comp_vcf), preset="vcf")
    uncomp_vcf.unlink()
    return comp_vcf


def simulate_contact_prob_matrix(
    fragments_df, p_cis: float = 0.8, allow_self: bool = False
) -> npt.NDArray[np.float_]:
    """Simulate contact probability matrix."""
    n = fragments_df.shape[0]
    p_trans = 1.0 - p_cis

    # set default value to probabilty of a trans conctact
    p = np.ones((n, n)) * p_trans
    for row in (
        fragments_df.groupby("chrom")
        .agg(
            [
                pl.col("fragment_idx").min().alias("min"),
                pl.col("fragment_idx").max().alias("max"),
            ]
        )
        .rows()
    ):
        _, min_idx, max_idx = row
        # set each cis block to probability cis (no length dependence)
        p[min_idx - 1: max_idx, min_idx - 1: max_idx] = p_cis
    # no self-contacts
    if not allow_self:
        np.fill_diagonal(p, 0)
    # np.nan_to_num(p, copy=False, nan=0.0)
    # renormalize
    p = p / p.sum(axis=1)

    return p


def simulate_concatemers(
    prob_matrix: npt.NDArray[np.float_],
    sizes: npt.NDArray[np.int_],
    random_state: Optional[Generator] = None,
) -> List[List[int]]:
    """Simulate concatemers."""
    if random_state is None:
        random_state = default_rng()
    n_frags = prob_matrix.shape[0]
    concatemers = []
    for size in sizes:
        # pick a random fragement to start walk
        start_idx = random_state.integers(0, n_frags)
        path = [start_idx]
        n_iter = 0
        max_iter = int(size * 1.5)
        while (len(path) < size) and n_iter < max_iter:
            # choose next fragment in proportion
            # to the column for that fragment
            next_frag = random_state.choice(
                n_frags, p=prob_matrix[:, start_idx])
            if next_frag not in path:
                path.append(next_frag)
            start_idx = next_frag
            n_iter += 1
        concatemers.append([_ for _ in path])  # store the dataframe index
    return concatemers


def simulate_read_sequence(
    *,
    ff: FastaFile,
    chrom: str,
    start: int,
    end: int,
    offsets: Optional[List[int]] = None,
    allele_values: Optional[List[str]] = None,
    strand: bool = True,
) -> Tuple[str, Tuple[int], Tuple[str]]:
    """Simulate read sequence."""
    seq = ff.fetch(chrom, start, end)
    snps = []
    if offsets is not None:
        assert allele_values is not None
        s = list(seq)
        assert len(offsets) == len(allele_values)
        for x, a in zip(offsets, allele_values):
            if seq[x] != a:
                snps.append((x, a))
                s[x] = a
        seq = "".join(s)
    if strand is False:
        seq = seq[::-1].translate(str.maketrans("ATCG", "TAGC"))
    else:
        assert strand is True
    if snps:
        _offsets, _alleles = zip(*snps)
    else:
        _offsets, _alleles = tuple(), tuple()
    return seq, _offsets, _alleles


def simulate_concatemer_fastqs(
    fastq: Path,
    reference_fasta: Path,
    fragments_df: pl.DataFrame,
    monomer_fastq: Optional[Path],
    p_cis: float = 0.8,
    num_concatemers: int = 100,
    mean_frags_per_concatemer: int = 5,
    max_frags_per_concatemer: int = 50,
    snp_haplotypes: Optional[List[SNPHaplotypes]] = None,
    random_state: Optional[Generator] = None,
    contact_probs: Optional[npt.NDArray[np.float_]] = None,
) -> Tuple[
    Path,
    Optional[Path],
    List[ConcatemerReadSeq],
    List[List[MonomerReadSeq]],
    List[Walk],
    List[int],
]:
    """Simulate concatemer fastq."""
    ff = FastaFile(str(reference_fasta))

    if random_state is None:
        random_state = default_rng()

    if contact_probs is None:
        # pairwise probabilities of contacts between fragments
        contact_probs = simulate_contact_prob_matrix(fragments_df, p_cis=p_cis)

    fragment_to_snps = defaultdict(list)
    if snp_haplotypes:
        for snp in snp_haplotypes:
            if snp.fragment_id is not None:
                fragment_to_snps[snp.fragment_id].append(snp)
    concat_gen = random_concatemer_generator(
        fragments_df=fragments_df,
        reference_fasta=ff,
        contact_probs=contact_probs,
        random_state=random_state,
        max_concatemers=num_concatemers,
        fragment_to_snps=fragment_to_snps,
        mean_frags_per_concatemer=mean_frags_per_concatemer,
        max_frags_per_concatemer=max_frags_per_concatemer,
    )
    outfh = fastq.open("w")
    if monomer_fastq:
        monomer_fh = monomer_fastq.open("w")
    else:
        monomer_fh = None

    monomers = []
    concatemers = []
    walks = []
    haplotypes = []
    for (concatemer, _monomers, walk, haplotype) in concat_gen:
        outfh.write(concatemer.to_fastq_str(walk))
        concatemers.append(concatemer)
        monomers.append(_monomers)
        walks.append(walk)
        haplotypes.append(haplotype)
        if monomer_fh:
            for monomer in _monomers:
                monomer_fh.write(monomer.to_fastq_str())

    outfh.close()
    if monomer_fh:
        monomer_fh.close()
    return fastq, monomer_fastq, concatemers, monomers, walks, haplotypes


def monomers_to_dataframe(
    monomers: List[List[MonomerReadSeq]], haplotypes: List[int]
) -> pl.DataFrame:
    """Monomers to dataframe."""
    data = []

    for x, _m in enumerate(monomers):
        for m in _m:
            assert m.read_seq.align_info is not None
            data.append(
                {
                    "monomer_id": m.monomer_id,
                    "concatemer_id": m.concatemer_id,
                    "align.chrom": m.read_seq.align_info.ref_name,
                    "align.start": m.read_seq.align_info.ref_pos,
                    "align.end": m.read_seq.align_info.ref_end,
                    "align.strand": m.read_seq.align_info.strand,
                    "haplotype": haplotypes[x],
                }
            )

    return pl.DataFrame(
        data,
        orient="row",
    )


def concatemers_to_dataframe(
    concatemers: List[ConcatemerReadSeq], walks: List[Walk]
) -> pl.DataFrame:
    """Concatemer to dataframe."""
    return pl.DataFrame(
        [
            {"concatemer_id": c.concatemer_id, "num_monomers": len(walks[x])}
            for x, c in enumerate(concatemers)
        ],
        orient="row",
    )


@dataclass
class TestScenarioFileCollection(FileCollection):
    """Test scenario file collection."""

    ns_bam: Path = Path("{prefix}.name_sorted.bam")
    reference_fasta: Path = Path("{prefix}.fasta")
    concatemer_fastq: Path = Path("{prefix}.concatemers.fastq")
    concatemer_ubam: Path = Path("{prefix}.concatemers.bam")
    monomer_fastq: Path = Path("{prefix}.monomer.fastq")
    monomer_parquet: Path = Path("{prefix}.monomer.pq")
    phased_vcf: Path = Path("{prefix}.phased_variants.vcf.gz")
    contact_prob: Path = Path("{prefix}.contact_probs.npy")
    params_json: Path = Path("{prefix}.params.json")

    def truth_files(self):
        """Truth files."""
        for k, v in self.items():
            if k in [
                "monomer_fastq", "monomer_marquet",
                    "contact_prob", "params_json"]:
                yield (k, v)


@dataclass
class Scenario:
    """Scenario."""

    seed: int = 421
    genome_size: int = 5_000
    num_chroms: int = 2
    cut_rate: float = 0.005
    enzyme: str = "NlaIII"
    num_concatemers: int = 100
    num_haplotypes: int = 0
    variant_density: float = 0.05
    p_cis: float = 0.8
    mean_frags_per_concatemer: int = 5
    max_frags_per_concatemer: int = 10
    temp_path: Path = field(default_factory=lambda: Path(mkdtemp()))

    fc: TestScenarioFileCollection = field(init=False)

    def __post_init__(self):
        """Post init."""
        self.fc = TestScenarioFileCollection.with_prefix(self.temp_path)
        self.to_json()
        self.random_state = default_rng(self.seed)
        self.chrom_lengths = {
            f"chr{x+1}": v
            for x, v in enumerate(
                sorted(
                    self.random_state.choice(
                        self.genome_size, size=self.num_chroms, replace=False
                    )
                )
            )
        }

    @classmethod
    def from_json(cls, p: Path, temp_path: Optional[Path] = None):
        """From json."""
        d = json.loads(p.read_text())
        if temp_path:
            d["temp_path"] = temp_path
        return cls(**d)

    def to_json(self, p: Optional[Path] = None):
        """To json."""
        p = self.fc.params_json if p is None else p
        with p.open("w") as fh:
            d = {
                k: v for k, v in asdict(self).items()
                if k not in ("temp_path", "fc")
                }
            fh.write(json.dumps(d))

    @cached_property
    def cut_sites(self) -> Dict[str, List[int]]:
        """Cut sites."""
        cut_sites = {}
        for chrom, length in self.chrom_lengths.items():
            cut_sites[chrom] = sorted(
                list(
                    self.random_state.integers(
                        # TODO: fix this arbitrary padding
                        20,
                        length - 20,
                        size=max(1, int(length * self.cut_rate)),
                    )
                )
            )
        return cut_sites

    @cached_property
    def fragments(self) -> List[GenomicFragment]:
        """Fragments."""
        id_iter = count(1)
        res = []
        for chrom, length in self.chrom_lengths.items():
            res.extend(
                GenomicFragment.from_cuts(
                    chrom, length, self.cut_sites[chrom], id_iter)
            )
        return res

    @cached_property
    def fragments_df(self) -> pl.DataFrame:
        """Fragments df."""
        return GenomicFragment.to_dataframe(self.fragments)

    @cached_property
    def snp_haplotypes(self) -> Optional[List[SNPHaplotypes]]:
        """Snp haplotypes."""
        if self.num_haplotypes > 0:
            return simulate_haplotypes(
                self.reference_fasta,
                self.num_haplotypes,
                self.variant_density,
                random_state=self.random_state,
            )
        else:
            return None

    @cached_property
    def fragment_to_snps(self):
        """Fragments to snps."""
        res = defaultdict(list)
        if self.snp_haplotypes:
            for snp in self.snp_haplotypes:
                if snp.fragment_id is not None:
                    res[snp.fragment_id].append(snp)
        return res

    @cached_property
    def haplotype_df(self) -> Optional[pl.DataFrame]:
        """Haplotype df."""
        if self.num_haplotypes > 0 and self.snp_haplotypes:
            return haplotypes_to_df(self.snp_haplotypes)
        else:
            return None

    @property
    def phased_vcf(self):
        """Phased vcf."""
        vcf = self.fc.phased_vcf
        if not vcf.exists():
            if self.snp_haplotypes is None:
                raise ValueError("Can't created phased vcf without snps")
            create_phased_vcf(self.snp_haplotypes, vcf, self.reference_fasta)
        return vcf

    @property
    def reference_fasta(self):
        """Ref fasta."""
        fasta = self.fc.reference_fasta
        if not fasta.exists():
            simulate_fasta_with_cut_sites(
                fasta,
                self.chrom_lengths,
                self.cut_sites,
                self.enzyme,
                random_state=self.random_state,
            )
        return fasta

    @property
    def ff(self):
        """Fasta file."""
        return FastaFile(str(self.reference_fasta))

    @cached_property
    def contact_prob_matrix(self) -> npt.NDArray[np.float_]:
        """Contact probability matrix."""
        return simulate_contact_prob_matrix(
            self.fragments_df, p_cis=self.p_cis)

    @cached_property
    def contact_prob_npy(self) -> Path:
        """Contact prob npy."""
        np.save(self.fc.contact_prob, self.contact_prob_matrix)
        return self.fc.contact_prob

    @cached_property
    def _concatemers_res(self):
        """Concatemers res."""
        return simulate_concatemer_fastqs(
            self.fc.concatemer_fastq,
            self.reference_fasta,
            self.fragments_df,
            monomer_fastq=self.fc.monomer_fastq,
            snp_haplotypes=self.snp_haplotypes,
            random_state=self.random_state,
            num_concatemers=self.num_concatemers,
            p_cis=self.p_cis,
            contact_probs=self.contact_prob_matrix,
            mean_frags_per_concatemer=self.mean_frags_per_concatemer,
            max_frags_per_concatemer=self.max_frags_per_concatemer,
        )

    @cached_property
    def concatemer_fastq(self) -> Path:
        """Concatemer fastq."""
        return self._concatemers_res[0]

    @cached_property
    def concatemer_ubam(self) -> Path:
        """Concatemer ubam."""
        return fastq_to_ubam([self.concatemer_fastq], self.fc.concatemer_ubam)

    @cached_property
    def monomer_fastq(self) -> Optional[Path]:
        """Monomer fastq."""
        return self._concatemers_res[1]

    @cached_property
    def concatemer_metadata(self) -> pl.DataFrame:
        """Concatemer metadata."""
        return concatemers_to_dataframe(self.concatemers, self.walks)

    @cached_property
    def concatemers(self) -> List[ConcatemerReadSeq]:
        """Concatemers."""
        return self._concatemers_res[2]

    @cached_property
    def monomers(self) -> List[List[MonomerReadSeq]]:
        """Monomers."""
        return self._concatemers_res[3]

    @cached_property
    def walks(self) -> List[Walk]:
        """Walk."""
        return self._concatemers_res[4]

    @cached_property
    def haplotypes(self) -> List[int]:
        """Haplotypes."""
        return self._concatemers_res[5]

    @cached_property
    def monomer_metadata(self) -> pl.DataFrame:
        """Monomer meta data."""
        return monomers_to_dataframe(self.monomers, self.haplotypes)

    @cached_property
    def monomer_parquet(self) -> Path:
        """Monomer parquet."""
        pq = self.fc.monomer_parquet
        if not pq.exists():
            self.monomer_metadata.write_parquet(pq)
        return pq

    @property
    def namesorted_bam(self):
        """Name sorted bam."""
        ns_bam = self.fc.ns_bam
        if not ns_bam.exists():
            minimap_exe = which("minimap2")
            if not minimap_exe:
                # logger.error("No minimap executable found")
                return None
            try:
                sp.check_call(
                    f"minimap2 -ax map-ont {self.reference_fasta} "
                    f"{self.monomer_fastq} | samtools sort -t "
                    f"{MOLECULE_TAG} -o {ns_bam}"
                )
            except Exception:
                raise

            return ns_bam

    def __str__(self):
        """Str."""
        return f"<Scenario {self.temp_path}>"

    def save_files(self):
        """Save files."""
        pass

    #    n = self.fragments_df.shape[0]
    #    p_cis = 1.0
    #    p_trans = 1.0 - p_cis

    #    # set default value to probabilty of a trans conctact
    #    p = np.ones((n, n)) * p_trans
    #    for row in (
    #        self.fragments_df.groupby("chrom")
    #        .fragment_id.agg(["min", "max"])
    #        .itertuples()
    #    ):
    #        # set each cis block to probability cis (no length dependence)
    #        p[row.min - 1 : row.max, row.min - 1 : row.max] = p_cis
    #    # no self-contacts
    #    np.fill_diagonal(p, 0)
    #    # renormalize
    #    p = p / p.sum(axis=1)

    #    num_concatemers = 100
    #    concatemer_size = self.random_state.poisson(5, num_concatemers)
    #    concatemers = []
    #    for x, size in enumerate(concatemer_size):
    #        start_idx = self.random_state.integers(0, n)
    #        path = [start_idx]
    #        counter = 0
    #        while (len(path) < size) and counter < 100:
    #            next_frag = np.random.choice(n, p=p[:, start_idx])
    #            if next_frag not in path:
    #                path.append(next_frag)
    #            counter += 1
    #        concatemers.append([_ + 1 for _ in path])  # convert back to
    #        fragement ids

    #    df = self.fragments_df.set_index(["fragment_id"]).sort_index()
    #    logger.debug(df.index[-2:])
    #    ff = FastaFile(self.reference_fasta)
    #    for idx, path in enumerate(concatemers):
    #        segments = []
    #        for row in df.loc[path, [
    #           "chrom", "start", "end"]].itertuples(index=False):
    #            logger.debug(row)
    #            segments.append(ff.fetch(row.chrom, row.start, row.end))
    #        seq = "".join(segments)
