from pathlib import Path
from typing import TYPE_CHECKING, List

from attrs import Factory, define, field
from pysam import AlignedSegment, AlignmentFile, AlignmentHeader
from pysam import sort as sort_bam  # pyright: ignore

from .log import get_logger
from .utils import FileCollection

logger = get_logger()

if TYPE_CHECKING:
    from .map import ConcatemerAlignData


@define
class MappingFileCollection(FileCollection):
    bam: Path = Path("{prefix}.bam")
    temp_bam: Path = Path("{prefix}.temp.bam")
    concatemers: Path = Path("{prefix}.concatemers.parquet")
    contacts: Path = Path("{prefix}.contacts.parquet")
    unmapped: Path = Path("{prefix}.unmapped.fastq")


@define
class MapWriter:
    fc: MappingFileCollection
    sam_header: AlignmentHeader
    reference_filename: Path
    _cache: List["ConcatemerAlignData"] = Factory(list)
    batch_size: int = 1000

    af: AlignmentFile = field(init=False)

    def __attrs_post_init__(self):

        self.af = AlignmentFile(
            str(self.fc.temp_bam),
            "wb",
            header=self.sam_header,
            reference_filename=str(self.reference_filename),
        )

    def __call__(self, rec: "ConcatemerAlignData"):
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
        sort_bam("--write-index", "-o", str(self.fc.bam), str(self.fc.temp_bam))
        self.fc.temp_bam.unlink()
