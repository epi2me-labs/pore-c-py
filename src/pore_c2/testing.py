from dataclasses import dataclass, field
from itertools import count
from pathlib import Path
from tempfile import mkdtemp
from typing import Dict, List, Optional, Tuple

import numpy as np
import numpy.typing as npt
from loguru import logger
from numpy.random import Generator, default_rng
from pandera.typing import DataFrame
from pysam import FastaFile

from pore_c2.digest import DigestFragment, VirtualDigestSchema


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
    enz = getattr(Restriction, enzyme, None)
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


def simulate_contact_prob_matrix(
    fragments_df: DataFrame[VirtualDigestSchema], p_cis: float = 0.8
) -> npt.NDArray[np.float_]:

    p_trans = 1 - p_cis
    n = fragments_df.shape[0]
    p_cis = 1.0
    p_trans = 1.0 - p_cis

    # set default value to probabilty of a trans conctact
    p = np.ones((n, n)) * p_trans
    for row in (
        fragments_df.groupby("chrom").fragment_id.agg(["min", "max"]).itertuples()
    ):
        # set each cis block to probability cis (no length dependence)
        p[row.min - 1 : row.max, row.min - 1 : row.max] = p_cis
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
        concatemers.append([_ + 1 for _ in path])  # convert back to fragement ids
    return concatemers


def simulate_concatemer_fastqs(
    fastq: Path,
    reference_fasta: Path,
    fragments_df: DataFrame[VirtualDigestSchema],
    p_cis: float = 0.8,
    num_concatemers: int = 100,
    mean_frags_per_concatemer: int = 5,
    max_frags_per_concatemer: int = 50,
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
    df = fragments_df.set_index("fragment_id")
    ff = FastaFile(reference_fasta)
    outfh = fastq.open("w")
    for idx, path in enumerate(concatemers):
        segments = []
        for row in df.loc[path, ["chrom", "start", "end"]].itertuples(index=False):
            segments.append(ff.fetch(row.chrom, row.start, row.end))
        seq = "".join(segments)
        id = "_".join(map(str, path))
        qual = "5" * len(seq)
        outfh.write(f"@CMER{idx}:{id}\n{seq}\n+\n{qual}\n")
        logger.debug(seq)


@dataclass
class Scenario:
    chrom_lengths: Dict[str, int]
    cut_rate: float = 0.005
    enzyme: str = "EcoRI"
    random_state: np.random.Generator = field(
        default_factory=lambda: np.random.default_rng()
    )

    cut_sites: Dict[str, List[int]] = field(init=False)
    fragments: List[DigestFragment] = field(init=False)
    fragments_df: DataFrame[VirtualDigestSchema] = field(init=False)
    temp_path: Path = field(default_factory=lambda: Path(mkdtemp()))

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
                DigestFragment.from_cuts(chrom, length, self.cut_sites[chrom], id_iter)
            )
        self.fragments_df = DataFrame[VirtualDigestSchema](self.fragments)

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
    def concatemer_fastq(self):
        fastq = self.temp_path / "concatemers.fastq"
        if not fastq.exists():
            simulate_concatemer_fastqs(
                fastq,
                self.reference_fasta,
                self.fragments_df,
                random_state=self.random_state,
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
