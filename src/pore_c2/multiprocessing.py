import queue
from itertools import count
from threading import Thread
from typing import Iterator

from mappy import Aligner, ThreadBuffer

from pore_c2.digest import Cutter
from pore_c2.map import map_concatemer_read
from pore_c2.reads import Read


class MappyWorkerThread(Thread):
    def __init__(
        self,
        *,
        aligner: Aligner,
        cutter: Cutter,
        input_queue: queue.Queue,
        output_queue: queue.Queue,
    ):
        super().__init__()
        self.aligner = aligner
        self.cutter = cutter
        self.input_queue = input_queue
        self.output_queue = output_queue

    def run(self):
        thrbuf = ThreadBuffer()
        counter = 0
        while True:
            read: Read = self.input_queue.get()
            if read is StopIteration:
                self.output_queue.put(read)
                break
            counter += 1
            if counter % 20 == 0:
                del thrbuf
                thrbuf = ThreadBuffer()
            results = map_concatemer_read(
                aligner=self.aligner, read=read, cutter=self.cutter, thread_buf=thrbuf
            )
            self.output_queue.put(results)


class MappyThreadPool(Thread):
    def __init__(
        self,
        *,
        read_iter: Iterator[Read],
        aligner: Aligner,
        cutter: Cutter,
        n_threads: int,
        maxsize: int = 2,
    ):
        super().__init__()
        self.read_iter = read_iter
        self.n_threads = n_threads
        self.work_queues, self.output_queues, self.workers = [], [], []
        for _ in range(self.n_threads):
            in_q, out_q = queue.Queue(maxsize), queue.Queue(maxsize)
            worker = MappyWorkerThread(
                aligner=aligner, cutter=cutter, input_queue=in_q, output_queue=out_q
            )
            self.work_queues.append(in_q)
            self.output_queues.append(out_q)
            self.workers.append(worker)

    def start(self):
        for worker in self.workers:
            worker.start()
        super().start()

    def __iter__(self):
        self.start()
        for i in count():
            item = self.output_queues[i % self.n_threads].get()
            if item is StopIteration:
                # drain other output queues
                for j in range(i + 1, i + self.n_threads):
                    self.output_queues[j % self.n_threads].get()
                break
            yield item

    def run(self):
        for i, read in enumerate(self.read_iter):
            self.work_queues[i % self.n_threads].put(read)
        for q in self.work_queues:
            q.put(StopIteration)
        for worker in self.workers:
            worker.join()
