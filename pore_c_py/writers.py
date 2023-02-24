"""Collection of file writers."""
from typing import List

import pyarrow as pa
import pyarrow.parquet as pq
import pysam


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

    def write(self, alignments: List[pysam.AlignedSegment]):
        """Write records."""

        def _aln_to_record(aln):
            return {
                "cid": aln.get_tag("MI"),
                "chrom": aln.reference_name,
                "start": aln.reference_start,
                "end": aln.reference_end,
                "num_fragments": 1}

        def is_colinear(
                align1: pysam.AlignedSegment,
                align2: pysam.AlignedSegment, tol: int = 1):
            """Check if two alignments are co-linear in the genome."""
            # if either is is unmapped then they can't be co-linear
            if align1.is_unmapped or align2.is_unmapped:
                return False
            if align1.reference_name != align2.reference_name:
                return False
            if align1.is_reverse != align2.is_reverse:
                return False
            delta = min(align1.reference_end, align2.reference_end) \
                - max(align1.reference_start, align2.reference_start)
            # overlaps or within a distance
            return delta > 0 or delta < tol

        def group_colinear(
                aligns: List[pysam.AlignedSegment], tol: int = 1
                ) -> List[List[pysam.AlignedSegment]]:
            """Group alignments into co-linear blocks."""
            if len(aligns) < 2:
                return [list(range(len(aligns)))]
            res = []
            block = []
            last = aligns[0]
            block = [last]
            for aln in aligns[1:]:
                if is_colinear(last, aln, tol=tol):
                    block.append(aln)
                else:
                    res.append(block)
                    block = [aln]
                last = aln
            res.append(block)
            return res

        pylist = list()
        if self.merge_distance is None:
            pylist = [
                _aln_to_record(aln) for aln in alignments
                if not aln.is_unmapped]
        else:
            pylist = list()
            for block in group_colinear(
                    [r.read_seq.align_info for r in alignments],
                    tol=self.merge_distance):
                if len(block) == 1:
                    if not block[0].is_unmapped:
                        pylist.append(_aln_to_record(block[0]))
                else:
                    rec = _aln_to_record(block[0])
                    rec["start"] = min(
                        block[0].reference_start, block[-1].reference_start)
                    rec["end"] = max(
                        block[0].reference_end, block[-1].reference_end)
                    rec["num_fragments"] = len(block)
                    pylist.append(rec)
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
