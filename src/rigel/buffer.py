"""rigel.buffer — the columnar fragment buffer, with Arrow IPC spill to disk.

The native scanner hands over finalized chunks of compact NumPy arrays (``inject_chunk``), variable-length transcript
sets in CSR layout. When the in-memory chunks exceed *max_memory_bytes*, a background writer spills the oldest to disk
as Arrow IPC (Feather v2) with LZ4. Scoring consumes the chunks once, in order (``iter_chunks_consuming``).
"""

import logging
import queue
import shutil
import tempfile
import threading
import weakref
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np

from .config import CONSTANTS, BamScanConfig

# The fragment classes ``fragment_classes`` assigns: the codes the native scorer routes on, read from
# it rather than restated. Other modules import them from here.
from .native import (
    FRAG_AMBIG_OPP_STRAND,
    FRAG_AMBIG_SAME_STRAND,
    FRAG_MULTIMAPPER,
    FRAG_UNAMBIG,
)

logger = logging.getLogger(__name__)

__all__ = ["FragmentBuffer"]


# ---------------------------------------------------------------------------
# FinalizedChunk — immutable chunk in compact NumPy arrays
# ---------------------------------------------------------------------------


#: A chunk's per-fragment arrays and their dtypes, the contract with the native side: one entry per fragment.
_FRAGMENT_ARRAYS = {
    "splice_type": np.uint8,
    "align_strand": np.uint8,
    "num_hits": np.int32,
    "ambig_strand": np.uint8,
    "frag_id": np.int64,
    "read_length": np.uint16,
    "genomic_footprint": np.int32,
    "genomic_start": np.int32,
    "nm": np.uint16,
    "clip": np.uint16,
}

#: A chunk's per-candidate arrays, parallel to ``t_indices`` over ``t_offsets``' CSR, by the column each is spilled
#: as, with their dtypes.
_CANDIDATE_ARRAYS = {
    "t_inds": ("t_indices", np.int32),
    "frag_lengths": ("frag_lengths", np.int32),
    "exon_bp": ("exon_bp", np.uint16),
}

#: A multimapper's place masses, CSR over rows (``place_offsets``), by the column each is spilled as, with their
#: dtypes; a chunk carries all of them or none.
_PLACE_ARRAYS = {
    "place_kind": np.int8,
    "place_index": np.int64,
    "place_mass": np.float32,
}


@dataclass(slots=True)
class FinalizedChunk:
    """An immutable chunk of fragment data in compact columnar arrays.

    Fixed-width fields are typed vectors; transcript index sets are CSR (``t_offsets[i]:t_offsets[i+1]`` indexes
    ``t_indices``). The dtypes are a contract with the native side:

    * ``frag_id`` is int64, for BAMs with more than 2 billion fragments;
    * ``t_offsets`` and ``t_indices`` are int32, so a chunk holds at most ~2.1 billion candidate entries;
    * per-candidate exon overlap bases and ``read_length`` are uint16, guarded at native append;
    * fragment lengths are int32: a transcript-space length on a long intron-spanning candidate can exceed 65535.
    """

    splice_type: np.ndarray  # uint8[N]
    align_strand: np.ndarray  # uint8[N]
    num_hits: np.ndarray  # int32[N]
    t_offsets: np.ndarray  # int32[N+1]
    t_indices: np.ndarray  # int32[M_t]
    frag_lengths: np.ndarray  # int32[M_t]  (parallel to t_indices, -1 = missing)
    exon_bp: np.ndarray  # uint16[M_t]  (parallel to t_indices)
    ambig_strand: np.ndarray  # uint8[N]
    frag_id: np.ndarray  # int64[N]
    read_length: np.ndarray  # uint16[N]
    genomic_footprint: np.ndarray  # int32[N]
    genomic_start: np.ndarray  # int32[N]
    nm: np.ndarray  # uint16[N]
    clip: np.ndarray  # uint16[N]  the hit's soft-clipped read bases
    size: int
    _fragment_classes: np.ndarray | None = None  # cached uint8[N]
    # Each row's place masses, CSR over rows: the deposit rule's mass at each object (region, boundary or sj by
    # ``rigel.native.PLACE_*``, flat ids) a multimapper's genic place deposits on, read by its capture weight. Empty on
    # every other row; ``None`` when the chunk carries none.
    place_offsets: np.ndarray | None = None  # int64[N+1]
    place_kind: np.ndarray | None = None  # int8[M_p]
    place_index: np.ndarray | None = None  # int64[M_p]
    place_mass: np.ndarray | None = None  # float32[M_p]

    @classmethod
    def from_raw(cls, raw: dict) -> "FinalizedChunk":
        """Build a chunk from the dict C++ ``FragmentChunk.finalize`` returns. The arrays arrive C-contiguous in
        the declared dtypes, so nothing is copied or cast."""
        return cls(
            **{name: raw[name] for name in _FRAGMENT_ARRAYS},
            **{name: raw[name] for name, _dtype in _CANDIDATE_ARRAYS.values()},
            t_offsets=raw["t_offsets"],
            size=int(raw["size"]),
            **{name: raw.get(name) for name in ("place_offsets", *_PLACE_ARRAYS)},
        )

    @property
    def memory_bytes(self) -> int:
        """Total bytes consumed by the underlying NumPy arrays."""
        names = (
            *_FRAGMENT_ARRAYS,
            "t_offsets",
            *(name for name, _dtype in _CANDIDATE_ARRAYS.values()),
            "place_offsets",
            *_PLACE_ARRAYS,
        )
        return sum(a.nbytes for a in (getattr(self, name) for name in names) if a is not None)

    @property
    def fragment_classes(self) -> np.ndarray:
        """Each fragment's class, uint8[N], computed once:

        * ``FRAG_UNAMBIG``: one strand, one transcript, one hit;
        * ``FRAG_AMBIG_SAME_STRAND``: one strand, several transcripts, one hit;
        * ``FRAG_AMBIG_OPP_STRAND``: transcripts on both strands, one hit;
        * ``FRAG_MULTIMAPPER``: more than one hit, whatever else holds.
        """
        if self._fragment_classes is not None:
            return self._fragment_classes
        n_transcripts = np.diff(self.t_offsets).astype(np.intp)
        classes = np.full(self.size, FRAG_UNAMBIG, dtype=np.uint8)
        classes[(self.ambig_strand == 0) & (n_transcripts > 1) & (self.num_hits == 1)] = (
            FRAG_AMBIG_SAME_STRAND
        )
        classes[(self.ambig_strand > 0) & (self.num_hits == 1)] = FRAG_AMBIG_OPP_STRAND
        classes[self.num_hits > 1] = FRAG_MULTIMAPPER
        object.__setattr__(self, "_fragment_classes", classes)
        return classes

    def to_scoring_arrays(self) -> dict:
        """The arrays ``StreamingScorer.score_chunk`` reads, by name.

        ``from_raw`` and ``load_chunk`` both store every array C-contiguous in the dtype the C++ kernel expects, so
        these are the chunk's own arrays.
        """
        return {
            "t_offsets": self.t_offsets,
            "t_indices": self.t_indices,
            "frag_lengths": self.frag_lengths,
            "exon_bp": self.exon_bp,
            "splice_type": self.splice_type,
            "align_strand": self.align_strand,
            "fragment_classes": self.fragment_classes,
            "frag_id": self.frag_id,
            "read_length": self.read_length,
            "genomic_footprint": self.genomic_footprint,
            "genomic_start": self.genomic_start,
            "nm": self.nm,
            "clip": self.clip,
        }


# ---------------------------------------------------------------------------
# Arrow IPC (Feather v2) spill / load
# ---------------------------------------------------------------------------


def spill_chunk(chunk: FinalizedChunk, path: Path) -> None:
    """Write a finalized chunk to disk as Arrow IPC with LZ4: each per-fragment array a column, each per-candidate
    and place-mass array a list column over its CSR."""
    import pyarrow as pa
    import pyarrow.feather as pf

    columns = {name: getattr(chunk, name) for name in _FRAGMENT_ARRAYS}
    for column, (name, _dtype) in _CANDIDATE_ARRAYS.items():
        columns[column] = pa.ListArray.from_arrays(chunk.t_offsets, getattr(chunk, name))
    if chunk.place_offsets is not None:
        place_offsets = chunk.place_offsets.astype(np.int32)
        for name in _PLACE_ARRAYS:
            columns[name] = pa.ListArray.from_arrays(place_offsets, getattr(chunk, name))
    pf.write_feather(pa.table(columns), str(path), compression="lz4")


def load_chunk(path: Path) -> FinalizedChunk:
    """Load a spilled chunk from an Arrow IPC file, every array in its declared dtype and owned."""
    import pyarrow.feather as pf

    table = pf.read_table(str(path))
    arrays = {
        name: table.column(name).to_numpy().astype(dtype)
        for name, dtype in _FRAGMENT_ARRAYS.items()
    }
    for column, (name, dtype) in _CANDIDATE_ARRAYS.items():
        values = table.column(column).combine_chunks()
        arrays[name] = values.values.to_numpy().astype(dtype)
        if column == "t_inds":
            arrays["t_offsets"] = values.offsets.to_numpy().astype(np.int32)
    if "place_kind" in table.column_names:
        kind = table.column("place_kind").combine_chunks()
        arrays["place_offsets"] = kind.offsets.to_numpy().astype(np.int64)
        for name, dtype in _PLACE_ARRAYS.items():
            arrays[name] = table.column(name).combine_chunks().values.to_numpy().astype(dtype)
    return FinalizedChunk(**arrays, size=len(table))


@dataclass(slots=True)
class PendingSpill:
    """A chunk handed to the background spill writer."""

    path: Path
    done: threading.Event
    memory_bytes: int
    size: int
    error: BaseException | None = None


class _SpillWriter:
    """One background thread that writes spilled chunks, in submission order, into a temporary
    directory it creates and removes. At most ``capacity`` submitted chunks are unwritten at a time.

    The buffer's finalizer holds this object and not the buffer, so it does not keep the buffer alive.
    """

    _STOP = object()

    def __init__(self, spill_dir: Path | None, capacity: int):
        if spill_dir is not None:
            spill_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir = Path(tempfile.mkdtemp(dir=spill_dir, prefix="rigel_buf_"))
        # SimpleQueue.put is reentrant, so a finalizer that interrupts the writer thread inside get()
        # can still queue the stop.
        self._queue: queue.SimpleQueue[object] = queue.SimpleQueue()
        self._slots = threading.BoundedSemaphore(capacity)
        self._closed = False
        self._thread = threading.Thread(
            target=self._run,
            name="rigel-spill-writer",
            daemon=True,
        )
        self._thread.start()

    def submit(self, chunk: FinalizedChunk, pending: PendingSpill) -> None:
        if self._closed:
            raise RuntimeError("Cannot submit to a closed spill writer.")
        self._slots.acquire()
        try:
            self._queue.put((chunk, pending))
        except BaseException:
            self._slots.release()
            raise

    def stop(self, *, wait: bool = True) -> None:
        """Queue the stop; the thread writes what was submitted before it, removes the directory and
        exits. With *wait*, joins it."""
        if self._closed:
            return
        self._closed = True
        self._queue.put(self._STOP)
        if wait:
            self._thread.join()

    def _run(self) -> None:
        while True:
            item = self._queue.get()
            if item is self._STOP:
                shutil.rmtree(self.temp_dir, ignore_errors=True)
                return

            chunk, pending = item
            del item
            try:
                logger.info(
                    "Writing spilled chunk (%s fragments, %.1f MB) -> %s",
                    f"{pending.size:,}",
                    pending.memory_bytes / 1024**2,
                    pending.path,
                )
                spill_chunk(chunk, pending.path)
                logger.info("Finished spilled chunk -> %s", pending.path)
            except BaseException as exc:
                pending.error = exc
                logger.exception("Failed to write spilled chunk -> %s", pending.path)
            finally:
                # Released before the chunk is reported written: kept while the writer blocks on its queue,
                # a written chunk is one more than the memory budget allows.
                del chunk
                pending.done.set()
                self._slots.release()


# ---------------------------------------------------------------------------
# FragmentBuffer — public API
# ---------------------------------------------------------------------------


class FragmentBuffer:
    """The columnar buffer of resolved fragments, spilling to disk past a memory budget.

    Holds the chunks the native scanner hands over (:meth:`inject_chunk`) and yields them in order
    (:meth:`iter_chunks_consuming`). Call :meth:`cleanup`, or use it as a context manager, to remove spilled files.

    Parameters
    ----------
    max_memory_bytes : int
        Memory for in-memory chunks before the oldest is spilled; 0 disables spilling.
    spill_dir : Path or None
        Directory for spilled chunk files; the system temp dir by default.
    """

    def __init__(
        self,
        max_memory_bytes: int = BamScanConfig.buffer_size_bytes,
        spill_dir: Path | None = None,
    ):
        self.max_memory_bytes = max_memory_bytes
        self._spill_dir = spill_dir

        self._chunks: deque[FinalizedChunk | PendingSpill] = deque()
        self._total_size = 0
        self._memory_bytes = 0
        self._n_spilled = 0
        self._spill_writer: _SpillWriter | None = None
        self._spill_finalizer: weakref.finalize | None = None

    # -- Properties -----------------------------------------------------------

    @property
    def total_fragments(self) -> int:
        """Total fragments handed to the buffer."""
        return self._total_size

    @property
    def memory_bytes(self) -> int:
        """Bytes consumed by in-memory finalized chunks."""
        return self._memory_bytes

    @property
    def n_spilled(self) -> int:
        """Number of chunks spilled to disk."""
        return self._n_spilled

    # -- Accepting chunks -----------------------------------------------------

    def inject_chunk(self, chunk: FinalizedChunk) -> None:
        """Append a finalized chunk and spill if over memory budget."""
        self._raise_completed_spill_errors()
        self._total_size += chunk.size
        self._memory_bytes += chunk.memory_bytes
        self._chunks.append(chunk)

        if self.max_memory_bytes > 0:
            while self._memory_bytes > self.max_memory_bytes:
                if not self._spill_oldest():
                    break

    def _spill_oldest(self) -> bool:
        """Spill the oldest in-memory chunk to disk.  Return True if spilled."""
        for i, chunk in enumerate(self._chunks):
            if isinstance(chunk, FinalizedChunk):
                writer = self._ensure_spill_writer()
                path = writer.temp_dir / f"chunk_{self._n_spilled:04d}.arrow"
                freed = chunk.memory_bytes
                pending = PendingSpill(
                    path=path,
                    done=threading.Event(),
                    memory_bytes=freed,
                    size=chunk.size,
                )
                writer.submit(chunk, pending)
                self._chunks[i] = pending
                self._memory_bytes -= freed
                self._n_spilled += 1
                logger.info(
                    "Queued spill for chunk %d (%s fragments, %.1f MB) -> %s",
                    i,
                    f"{chunk.size:,}",
                    freed / 1024**2,
                    path,
                )
                return True
        return False

    def _ensure_spill_writer(self) -> _SpillWriter:
        if self._spill_writer is None:
            self._spill_writer = _SpillWriter(
                self._spill_dir, capacity=CONSTANTS.resources.spill_queue_chunks
            )
            # Stops the writer, which then removes its directory, if cleanup() is never called. It never waits:
            # garbage collection can run it on any thread, the writer's own included, and on one holding a lock the
            # writer needs (a logging handler's, when a collection runs inside a log call), where waiting would
            # deadlock.
            self._spill_finalizer = weakref.finalize(self, self._spill_writer.stop, wait=False)
        return self._spill_writer

    def _stop_spill_writer(self) -> None:
        if self._spill_finalizer is not None:
            self._spill_finalizer.detach()  # the stop below, which waits, replaces it
            self._spill_writer.stop()
            self._spill_finalizer = None
            self._spill_writer = None

    def _spill_error(self, pending: PendingSpill) -> RuntimeError:
        return RuntimeError(f"Failed to spill buffer chunk to {pending.path}")

    def _wait_pending_spill(self, pending: PendingSpill) -> Path:
        pending.done.wait()
        if pending.error is not None:
            raise self._spill_error(pending) from pending.error
        return pending.path

    def _wait_all_pending_spills(self) -> None:
        failed: PendingSpill | None = None
        for chunk_ref in list(self._chunks):
            if isinstance(chunk_ref, PendingSpill):
                chunk_ref.done.wait()
                if chunk_ref.error is not None and failed is None:
                    failed = chunk_ref
        if failed is not None:
            raise self._spill_error(failed) from failed.error

    def _raise_completed_spill_errors(self) -> None:
        for chunk_ref in self._chunks:
            if (
                isinstance(chunk_ref, PendingSpill)
                and chunk_ref.done.is_set()
                and chunk_ref.error is not None
            ):
                raise self._spill_error(chunk_ref) from chunk_ref.error

    # -- Iteration ------------------------------------------------------------

    def iter_chunks_consuming(self) -> Iterator[FinalizedChunk]:
        """Yield the chunks in order, releasing each once the caller advances; spilled files are deleted as they are
        read, and the buffer is empty afterwards."""
        while self._chunks:
            chunk_ref = self._chunks.popleft()
            if isinstance(chunk_ref, PendingSpill):
                path = self._wait_pending_spill(chunk_ref)
                chunk = load_chunk(path)
                path.unlink(missing_ok=True)
            else:
                chunk = chunk_ref
                self._memory_bytes -= chunk.memory_bytes
            yield chunk
            # The caller has it; held here as well, it would stay in memory while the next one loads.
            del chunk

    # -- Cleanup --------------------------------------------------------------

    def cleanup(self) -> None:
        """Wait for pending spills, then stop the writer, which removes the spill directory."""
        error: RuntimeError | None = None
        try:
            self._wait_all_pending_spills()
        except RuntimeError as exc:
            error = exc
        finally:
            self._stop_spill_writer()
        if error is not None:
            raise error

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.cleanup()
