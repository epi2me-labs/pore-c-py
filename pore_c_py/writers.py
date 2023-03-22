"""Collection of file writers."""
import collections
import json
from typing import List

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from pore_c_py import align_tools, utils


class StatsWriter:
    """Stats writer."""

    def __init__(self, path):
        """Init."""
        self.path = path
        self.concatemer_count = 0
        self.cardinality_count = collections.Counter()
        self.pair_count = collections.Counter()
        self.cis_trans = collections.Counter()

    def append(self, pair):
        """Write records."""
        self.concatemer_count += 1
        cardinality = utils.MonomerData.from_pysam(pair.left).subread_idx
        if cardinality in self.cardinality_count.keys():
            self.cardinality_count[cardinality] += 1
        else:
            self.cardinality_count[cardinality] = 1
        self.pair_count[pair.state.name] += 1
        if pair.state == align_tools.PairState.both:
            if pair.left.reference_name == pair.right.reference_name:
                self.cis_trans["cis"] += 1
            else:
                self.cis_trans["trans"] += 1

    def close(self):
        """Close."""
        d = {
            "cardinality": self.cardinality_count,
            "pair_count": self.pair_count,
            "cis_trans": self.cis_trans
        }
        with self.path.open("w") as fh:
            fh.write(json.dumps(d))

    def __enter__(self):
        """Enter context."""
        pass

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """Exit context, closing json file."""
        self.close()


class ChromunityWriter:
    """Chromunity writer."""

    def __init__(self, path, merge_distance=None):
        """Its the initializer init."""
        # TODO: see if chromunity can just parse this from the BAM
        self.path = path
        self.merge_distance = merge_distance
        self.schema = pa.schema([
            ("cid", pa.string()),
            ("chrom", pa.string()),
            ("start", pa.uint32()),
            ("end", pa.uint32()),
            ("num_fragments", pa.uint32())])
        self._writer = pq.ParquetWriter(str(self.path), self.schema)
        self.counter = 0

    def _aln_to_record(self, aln):
        return {
            "cid": aln.get_tag("MI"),
            "chrom": aln.reference_name,
            "start": aln.reference_start,
            "end": aln.reference_end,
            "num_fragments": 1}

    def get_pylist(self, aligns):
        """Get pylist."""
        pylist = list()
        for block in align_tools.group_colinear(
                [r for r in aligns],
                tol=self.merge_distance):
            if len(block) == 1:
                if not block[0].is_unmapped:
                    pylist.append(self._aln_to_record(block[0]))
            else:
                rec = self._aln_to_record(block[0])
                rec["start"] = min(
                    block[0].reference_start, block[-1].reference_start)
                rec["end"] = max(
                    block[0].reference_end, block[-1].reference_end)
                rec["num_fragments"] = len(block)
                pylist.append(rec)
        return pylist

    def write(self, alignments: List[pysam.AlignedSegment]):
        """Write records."""
        pylist = list()
        if self.merge_distance is None:
            pylist = [
                self._aln_to_record(aln) for aln in alignments
                if not aln.is_unmapped]
        else:
            pylist = self.get_pylist(alignments)
        if len(pylist) > 0:
            batch = pa.RecordBatch.from_pylist(pylist, schema=self.schema)
            self._writer.write_batch(batch)
            self.counter += len(pylist)

    def close(self):
        """Close the parquet writer."""
        if self._writer is not None:
            self._writer.close()

    def __enter__(self):
        """Enter context."""
        # parquet writer is opened on creation
        pass

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """Exit context, closing parquet file."""
        self.close()
