from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Iterator, List, Optional

from pysam import AlignedSegment, AlignmentFile, AlignmentHeader, FastqFile

from ..model import AlignData
from ..settings import DEFAULT_ALIGN_HEADER

if TYPE_CHECKING:
    pass


class ReadIter(metaclass=ABCMeta):
    def __init__(self, path: Path, header: AlignmentHeader = DEFAULT_ALIGN_HEADER):
        self.path = path
        self.header = header

    @staticmethod
    def load(p: Path, **kwds) -> "ReadIter":
        map = {
            FastqReadIter: [".fq", ".fastq", ".fq.gz", ".fastq.gz"],
            SamReadIter: [".sam", ".bam", ".cram"],
        }
        for cls, exts in map.items():
            for ext in exts:
                if p.name.endswith(ext):
                    return cls(p, **kwds)
        raise ValueError(f"No handler for file extension: {p.name}")

    @abstractmethod
    def __iter__(self) -> Iterator[AlignData]:
        ...

    def iter_pysam(self) -> Iterator[AlignedSegment]:
        for align in self:
            yield (align.to_pysam(self.header))


class FastqReadIter(ReadIter):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

    def __iter__(self) -> Iterator[AlignData]:
        reader = FastqFile(str(self.path))
        for record in reader:
            yield AlignData.from_fastq(record)
        return


class SamReadIter(ReadIter):
    def __iter__(self) -> Iterator[AlignData]:
        reader = AlignmentFile(
            str(self.path),
            check_sq=False,
        )  # will warn about missing index for ubams
        for align in reader.fetch(until_eof=True):
            if align.is_unmapped or not (align.is_supplementary | align.is_secondary):
                yield AlignData.from_alignment(align)
        return


class ReadWriter(metaclass=ABCMeta):
    def __init__(self, path: Path, header: AlignmentHeader = DEFAULT_ALIGN_HEADER):
        self.path = path
        self.header = header
        self.read_counter = 0
        self.base_counter = 0

    @classmethod
    def load(cls, path: Path, header: Optional[AlignmentHeader] = None):
        if path.suffix in [".fq", ".fastq"]:
            return FastQReadWriter(path)
        elif path.suffix in [".cram", ".sam", ".bam"]:
            if header is None:
                header = DEFAULT_ALIGN_HEADER
            return UBamReadWriter(path, header=header)
        else:
            raise ValueError(f"Unsupported suffix {path}")

    @abstractmethod
    def write_item(self, align: AlignData):
        ...

    def consume(self, align_iter: Iterable[List[AlignData]]):
        for aligns in align_iter:
            for a in aligns:
                self(a)

    def __call__(self, align: AlignData):
        self.write_item(align)
        self.read_counter += 1
        self.base_counter += len(align.seq)


class UBamReadWriter(ReadWriter):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        if self.path.suffix == ".bam":
            mode = "wb"
        elif self.path.suffix == ".sam":
            mode = "w"
        elif self.path.suffix == ".cram":
            mode = "wc"
        else:
            raise NotImplementedError(self.path.suffix)
        self.af = AlignmentFile(str(self.path), mode=mode, header=self.header)

    def write_item(self, align: AlignData):
        self.af.write(align.to_pysam(header=self.header))


class FastQReadWriter(ReadWriter):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.fh = self.path.open("w")

    def write_item(self, align: AlignData):
        self.fh.write(align.to_fastq())
