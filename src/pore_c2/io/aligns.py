from pathlib import Path
from threading import Thread
from typing import Any, Iterator, List

from pysam import AlignedSegment, AlignmentFile, AlignmentHeader
from pysam import sort as sort_bam  # pyright: ignore

from ..align import MappingFileCollection
from ..log import get_logger

logger = get_logger()


class MapWriter(Thread):
    def __init__(
        self,
        result_iter: Iterator[Any],
        fc: MappingFileCollection,
        sam_header: AlignmentHeader,
        reference_filename: Path,
        batch_size: int = 1000,
    ):
        super().__init__()
        self.result_iter = result_iter
        self.fc = fc
        self.sam_header = sam_header
        self.reference_filename = reference_filename
        self.batch_size = batch_size
        self._cache: List[Any] = []
        self.af = AlignmentFile(
            str(self.fc.temp_bam),
            "wb",
            header=self.sam_header,
            reference_filename=str(self.reference_filename),
        )

    def run(self):
        for item in self.result_iter:
            self(item)

    def __call__(self, rec):
        self._cache.append(rec)
        if len(self._cache) > self.batch_size:
            self.flush()

    def flush(self):
        logger.info("Flushing cache")
        for rec in self._cache:
            for sam_str in rec.to_sam():
                self.af.write(AlignedSegment.fromstring(sam_str, self.af.header))
        self._cache = []

    def close(self):
        # TODO: need to add some exception handline here
        self.flush()
        self.af.close()
        logger.debug(f"Sorting bam {self.fc.temp_bam}")
        sort_bam("--write-index", "-o", str(self.fc.bam), str(self.fc.temp_bam))
        self.fc.temp_bam.unlink()
