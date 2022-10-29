import json
from pathlib import Path

import polars as pl
import pytest
from attrs import define

from pore_c2.io import FileCollection


@define
class WorkflowFC(FileCollection):
    concat_fq: Path = Path("{prefix}/input/scenario.concatemers.fastq")
    monomer_pq: Path = Path("{prefix}/input/scenario.monomer.pq")
    params_json: Path = Path("{prefix}/input/scenario.params.json")
    pairs_stats: Path = Path("{prefix}/output/scenario.processed.pairs.stats.txt")

    def monomer_df(self) -> pl.DataFrame:
        return pl.read_parquet(self.monomer_pq)

    def pairs_stats_data(self):
        d = dict([tuple(_.strip().split("\t")) for _ in self.pairs_stats.open()])
        return d

    def params(self):
        return json.loads(self.params_json.read_text())


@pytest.mark.workflow("full")
def test_cis_trans(workflow_dir, tol=0.1):
    fc = WorkflowFC.with_prefix(workflow_dir)
    expected = fc.params()["p_cis"]

    pairs_stats = fc.pairs_stats_data()
    cis, trans = map(int, [pairs_stats[_] for _ in ["cis", "trans"]])
    total_contacts = cis + trans
    observed = float(cis) / total_contacts
    delta = expected - observed
    assert delta < tol
