from dataclasses import dataclass, field
from itertools import count
from pathlib import Path
from tempfile import mkdtemp
from typing import Dict, List, Optional, Tuple

import numpy as np
import numpy.typing as npt
import polars as pl
from numpy.random import Generator, default_rng
from pysam import FastaFile, VariantFile, VariantHeader, tabix_compress, tabix_index

from pore_c2.monomers import GenomicFragment


def simulate_sequence_with_cut_sites(
    enzyme: str = "EcoRI",
    seq_length: int = 1000,
    cut_sites: List[int] = [50, 300, 750],
    random_state: Optional[Generator] = None,
) -> Tuple[List[int], str]:

    from Bio import Restriction
    from Bio.Seq import Seq

    if random_state is None:
        random_state = default_rng()
    enz = getattr(Restriction, enzyme)
    offset = enz.fst5
    length = len(enz.site)
    while True:
        seq = random_state.choice(list("ACGT"), size=seq_length)
        if len(enz.search(Seq("".join(seq)))) == 0:
            break
    for p in cut_sites:
        start = p - offset
        seq[start : start + length] = list(enz.site)
    return cut_sites, "".join(seq)


def simulate_fasta_with_cut_sites(
    fasta: Path,
    seq_length: Dict[str, int],
    cut_sites: Dict[str, List[int]],
    enzyme: str = "EcoRI",
    random_state: Optional[Generator] = None,
) -> Path:
    if random_state is None:
        random_state = default_rng()
    fh = fasta.open("w")
    for chrom, length in seq_length.items():
        positions = cut_sites[chrom]
        _, seq = simulate_sequence_with_cut_sites(
            enzyme, cut_sites=positions, seq_length=length, random_state=random_state
        )
        fh.write(f">{chrom}\n{seq}\n")
    fh.close()
    return fasta


def simulate_haplotypes(
    reference_fasta: Path,
    num_haplotypes: int = 3,
    var_density: float = 0.01,
    random_state: Optional[Generator] = None,
) -> pl.DataFrame:
    if random_state is None:
        random_state = default_rng()

    ff = FastaFile(str(reference_fasta))
    mutations = {}
    bases = list("ATGC")
    for b in bases:
        mutations[b] = [_ for _ in bases if _ != "b"]

    df = None
    for chrom in ff.references:
        length = ff.get_reference_length(chrom)
        num_vars = max(2, int(np.rint(length * var_density)))
        snp_positions = sorted(
            random_state.choice(int(length), size=num_vars, replace=False)
        )
        # assume bi-allelic
        ref_allele = [ff.fetch(chrom, _, _ + 1) for _ in snp_positions]
        alt_allele = [random_state.choice(mutations[_], 1)[0] for _ in ref_allele]
        # matrix of
        haplotypes = random_state.choice([0, 1], size=(num_vars, num_haplotypes))
        _df = (
            pl.DataFrame(
                {
                    **{
                        "chrom": [chrom] * num_vars,
                        "pos": snp_positions,
                        "ref": ref_allele,
                        "alt": alt_allele,
                    },
                    **{f"HT{x}": haplotypes[:, x] for x in range(num_haplotypes)},
                }
            )
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
        if df is None:
            df = _df
        else:
            df = df.vstack(_df)
    assert df is not None
    return df


def create_phased_vcf(haplotype_df: pl.DataFrame, vcf: Path, reference_fasta: Path):

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
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">')
    header.add_line(
        '##FORMAT=<ID=PQ,Number=1,Type=Integer,Description="Phasing Quality">'
    )
    header.add_sample("SAMPLE_01")
    vcf_out = VariantFile(str(uncomp_vcf), "w", header=header)
    vcf_out.close()

    with uncomp_vcf.open("a") as fh:
        for (chrom, pos, ref, alt, gt0, gt1) in haplotype_df.rows():
            # record = vcf_out.new_record(
            #    contig=chrom,
            #    start=pos,
            #    stop=pos,
            #    alleles=[ref, alt],
            #    samples={"SAMPLE_01": ""
            # )
            # vcf_out.write(record)
            fh.write(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:PS\t{gt0}|{gt1}:0\n"
            )
    tabix_compress(str(uncomp_vcf), str(comp_vcf))
    tabix_index(str(comp_vcf), preset="vcf")
    uncomp_vcf.unlink()
    return comp_vcf


def simulate_contact_prob_matrix(
    fragments_df, p_cis: float = 0.8
) -> npt.NDArray[np.float_]:

    p_trans = 1 - p_cis
    n = fragments_df.shape[0]
    p_cis = 1.0
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
        p[min_idx - 1 : max_idx, min_idx - 1 : max_idx] = p_cis
    # no self-contacts
    np.fill_diagonal(p, 0)
    # renormalize
    p = p / p.sum(axis=1)

    return p


def simulate_concatemers(
    prob_matrix: npt.NDArray[np.float_],
    sizes: npt.NDArray[np.int_],
    random_state: Optional[Generator] = None,
) -> List[List[int]]:
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
            # choose next fragment in proportion to the row for that fragment
            next_frag = np.random.choice(n_frags, p=prob_matrix[:, start_idx])
            if next_frag not in path:
                path.append(next_frag)
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
):
    seq = ff.fetch(chrom, start, end)
    snps = []
    if offsets is not None:
        assert allele_values is not None
        s = list(seq)
        assert len(offsets) == len(allele_values)
        for x, a in zip(offsets, allele_values):
            if seq[x] != a:
                snps.append(x)
                s[x] = a
        seq = "".join(s)
    return seq, snps


def simulate_concatemer_fastqs(
    fastq: Path,
    reference_fasta: Path,
    fragments_df: pl.DataFrame,
    p_cis: float = 0.8,
    num_concatemers: int = 100,
    mean_frags_per_concatemer: int = 5,
    max_frags_per_concatemer: int = 50,
    haplotype_df: Optional[pl.DataFrame] = None,
    num_haplotypes: int = 0,
    random_state: Optional[Generator] = None,
):
    if random_state is None:
        random_state = default_rng()

    # pairwise probabilities of contacts between fragments
    p = simulate_contact_prob_matrix(fragments_df, p_cis=p_cis)

    # generate distribution of fragments per concatemer
    sizes = random_state.poisson(mean_frags_per_concatemer, num_concatemers).clip(
        1, max_frags_per_concatemer
    )
    concatemers = simulate_concatemers(p, sizes, random_state=random_state)
    if num_haplotypes > 0:
        assert haplotype_df is not None
        snps = (
            haplotype_df.join_asof(
                (
                    fragments_df.select(
                        ["chrom", "start", "end", "fragment_id"]
                    ).with_column(pl.col("start").alias("asof_end"))
                ),
                by="chrom",
                left_on="pos",
                right_on="asof_end",
            )
            .filter(pl.col("pos") >= pl.col("start"))
            .with_column((pl.col("pos") - pl.col("start")).alias("offset"))
            .with_column(pl.concat_list(["ref", "alt"]).alias("alleles"))
            # .with_column(pl.col("alleles").take(pl.col("HT.0")).alias("test"))
            .groupby(
                ["chrom", "start", "end", "fragment_id", "haplotype_id"],
                maintain_order=True,
            )
            .agg(
                [
                    # pl.count().alias("num_vars"),
                    pl.col("allele_value"),
                    pl.col("offset"),
                ]
            )
        )
        haplotypes = random_state.integers(0, num_haplotypes, size=len(concatemers))
    else:
        snps = None
        haplotypes = [None] * len(concatemers)

    ff = FastaFile(str(reference_fasta))
    outfh = fastq.open("w")
    data = []
    for idx, (path, haplotype) in enumerate(zip(concatemers, haplotypes)):
        segments = []
        id = []

        frag_locs = fragments_df[
            pl.Series(path), ["chrom", "start", "end", "fragment_id"]
        ]
        if snps:
            frag_locs = frag_locs.join(
                snps.filter(pl.col("haplotype_id") == f"HT{haplotype}"),
                on=["chrom", "start", "end", "fragment_id"],
                how="left",
            ).sort(["fragment_id"])
        else:
            frag_locs = frag_locs.with_columns(
                [
                    pl.lit(None).alias("haplotype_id"),
                    pl.lit(None).alias("allele_value"),
                    pl.lit(None).alias("offset"),
                ]
            )
            # snps = frag_locs.join_asof(
            #    haplotype_df.with_column(pl.col("pos").alias("loc")),
            #    by="chrom",
            #    left_on="start",
            #    right_on="pos",
            # ).filter(pl.col("loc") <= pl.col("end"))

            # raise ValueError(snps)
        for (
            chrom,
            start,
            end,
            fragment_id,
            _,  # haplotype_id
            allele_value,
            offset,
        ) in frag_locs.rows():
            id.append(fragment_id)
            _seq, _ = simulate_read_sequence(
                ff=ff,
                chrom=chrom,
                start=start,
                end=end,
                offsets=offset,
                allele_values=allele_value,
            )
            segments.append(_seq)

        seq = "".join(segments)
        fragment_str = f"fragments={','.join(id)}"
        qual = "5" * len(seq)
        outfh.write(f"@CMER{idx} {fragment_str}\n{seq}\n+\n{qual}\n")
        data.append({"concatemer_id": f"CMER{idx}", "num_segments": len(segments)})
        # logger.debug(seq)
    return pl.DataFrame(data, orient="row")


@dataclass
class Scenario:
    chrom_lengths: Dict[str, int]
    cut_rate: float = 0.005
    enzyme: str = "EcoRI"
    random_state: np.random.Generator = field(
        default_factory=lambda: np.random.default_rng()
    )
    cut_sites: Dict[str, List[int]] = field(init=False)
    fragments: List[GenomicFragment] = field(init=False)
    fragments_df: pl.DataFrame = field(init=False)
    temp_path: Path = field(default_factory=lambda: Path(mkdtemp()))
    num_concatemers: int = 100
    num_haplotypes: int = 0
    variant_density: float = 0.0

    def __post_init__(self):
        self.cut_sites = {}
        self.fragments = []
        id_iter = count(1)
        for chrom, length in self.chrom_lengths.items():
            self.cut_sites[chrom] = sorted(
                list(
                    self.random_state.integers(
                        20, length - 20, size=int(length * self.cut_rate)
                    )
                )
            )
            self.fragments.extend(
                GenomicFragment.from_cuts(chrom, length, self.cut_sites[chrom], id_iter)
            )
        self.fragments_df = GenomicFragment.to_dataframe(self.fragments)
        if self.num_haplotypes > 0:
            self.haplotype_df = simulate_haplotypes(
                self.reference_fasta, self.num_haplotypes, self.variant_density
            )
        else:
            self.haplotype_df = None

    @property
    def phased_vcf(self):
        vcf = self.temp_path / "genotypes.vcf.gz"
        if self.haplotype_df is None:
            raise ValueError
        if not vcf.exists():
            create_phased_vcf(self.haplotype_df, vcf, self.reference_fasta)
        return vcf

    @property
    def reference_fasta(self):
        fasta = self.temp_path / "genome.fasta"
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
        return FastaFile(str(self.reference_fasta))

    @property
    def concatemer_fastq(self):
        fastq = self.temp_path / "concatemers.fastq"
        if not fastq.exists():
            self.concatemer_metadata = simulate_concatemer_fastqs(
                fastq,
                self.reference_fasta,
                self.fragments_df,
                haplotype_df=self.haplotype_df,
                num_haplotypes=self.num_haplotypes,
                random_state=self.random_state,
                num_concatemers=self.num_concatemers,
            )
        return fastq

    # @property
    # def contact_prob_matrix(self):

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
    #        concatemers.append([_ + 1 for _ in path])  # convert back to fragement ids

    #    df = self.fragments_df.set_index(["fragment_id"]).sort_index()
    #    logger.debug(df.index[-2:])
    #    ff = FastaFile(self.reference_fasta)
    #    for idx, path in enumerate(concatemers):
    #        segments = []
    #        for row in df.loc[path, ["chrom", "start", "end"]].itertuples(index=False):
    #            logger.debug(row)
    #            segments.append(ff.fetch(row.chrom, row.start, row.end))
    #        seq = "".join(segments)
