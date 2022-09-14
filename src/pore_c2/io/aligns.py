from pathlib import Path
from threading import Thread
from typing import Any, Iterable, Iterator, List, Optional

from pysam import AlignedSegment, AlignmentFile, AlignmentHeader
from pysam import sort as sort_bam  # pyright: ignore

from ..aligns import MappingFileCollection
from ..log import get_logger
from ..model import AlignData

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


class AlignIter:
    def __init__(self, bam: Path, header: Optional[AlignmentHeader] = None):
        self.bam = bam
        self.header = header

    def __iter__(self) -> Iterator[AlignData]:
        reader = AlignmentFile(str(self.bam))
        if self.header is None:
            self.header = reader.header
        for align in reader.fetch(until_eof=True):
            yield AlignData.from_alignment(align)

    def iter_pysam(self) -> Iterator[AlignedSegment]:
        for align in self:
            yield (align.to_pysam(self.header))


def get_aligns(paths: Iterable[Path]) -> Iterable[AlignData]:
    for p in paths:
        reader = AlignIter(p)
        for align in reader:
            yield align
