"""
rigel.buffer — Columnar fragment buffer with Arrow IPC disk spill.

Stores resolved fragment data in memory-efficient columnar NumPy arrays
using CSR (Compressed Sparse Row) layout for variable-length transcript
and gene index sets.  When memory usage exceeds a configurable threshold,
completed chunks are spilled to disk as Arrow IPC (Feather v2) files
with LZ4 compression.

Architecture
------------
- The native scanner hands over finalized chunks of compact NumPy arrays
  (``inject_chunk``).
- When total in-memory chunk size exceeds *max_memory_bytes*,
  the oldest chunk is spilled to disk as Arrow IPC with LZ4.
- Scoring consumes the chunks once, in order (``iter_chunks_consuming``).
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

logger = logging.getLogger(__name__)

__all__ = ["FragmentBuffer"]

# Fragment classification constants used by fragment_classes property
FRAG_UNAMBIG: int = 0  # same-strand, 1 transcript, NH=1
FRAG_AMBIG_SAME_STRAND: int = 1  # same-strand, >1 transcript, NH=1
FRAG_AMBIG_OPP_STRAND: int = 2  # ambig-strand transcripts, NH=1
FRAG_MULTIMAPPER: int = 3  # NH > 1 (multimapped molecule)
FRAG_CHIMERIC: int = 4  # chimeric fragment (disjoint transcript sets)


# ---------------------------------------------------------------------------
# _FinalizedChunk — immutable chunk in compact NumPy arrays
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class _FinalizedChunk:
    """Immutable chunk of fragment data in compact columnar arrays.

    Fixed-width fields are stored as typed NumPy vectors.
    Variable-width transcript index sets use CSR (Compressed Sparse Row)
    layout: ``offsets[i]:offsets[i+1]`` indexes the flat ``indices`` array.

    Per-fragment strand mixing (``ambig_strand``) is cached as a uint8
    array emitted by the native resolver.

        Buffer dtype contract:
        * ``frag_id`` stays int64 to accommodate BAM files with more than
            2 billion fragments.
        * ``t_offsets`` and ``t_indices`` stay int32, limiting a chunk to
            ~2.1 billion candidate entries.
        * Per-candidate exon overlap bases are uint16; native append
            guards reject values above 65535.
        * Fragment lengths stay int32 because real transcript-space
            fragments can exceed 65535 on long intron-spanning candidates.
        * ``read_length`` is uint16 and guarded at native append.
    """

    splice_type: np.ndarray  # uint8[N]
    align_strand: np.ndarray  # uint8[N]
    num_hits: np.ndarray  # uint16[N]
    chimera_type: np.ndarray  # uint8[N]
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
    size: int
    _fragment_classes: np.ndarray | None = None  # cached uint8[N]

    @classmethod
    def from_raw(cls, raw: dict) -> "_FinalizedChunk":
        """Build a chunk from the dict C++ ``FragmentAccumulator.finalize`` returns. The native side
        hands over its vectors as C-contiguous arrays already carrying the declared dtypes, so nothing
        is copied or cast."""
        return cls(
            splice_type=raw["splice_type"],
            align_strand=raw["align_strand"],
            num_hits=raw["num_hits"],
            chimera_type=raw["chimera_type"],
            t_offsets=raw["t_offsets"],
            t_indices=raw["t_indices"],
            frag_lengths=raw["frag_lengths"],
            exon_bp=raw["exon_bp"],
            ambig_strand=raw["ambig_strand"],
            frag_id=raw["frag_id"],
            read_length=raw["read_length"],
            genomic_footprint=raw["genomic_footprint"],
            genomic_start=raw["genomic_start"],
            nm=raw["nm"],
            size=int(raw["size"]),
        )

    @property
    def memory_bytes(self) -> int:
        """Total bytes consumed by the underlying NumPy arrays."""
        return sum(
            a.nbytes
            for a in (
                self.splice_type,
                self.align_strand,
                self.num_hits,
                self.chimera_type,
                self.t_offsets,
                self.t_indices,
                self.frag_lengths,
                self.exon_bp,
                self.ambig_strand,
                self.frag_id,
                self.read_length,
                self.genomic_footprint,
                self.genomic_start,
                self.nm,
            )
        )

    @property
    def fragment_classes(self) -> np.ndarray:
        """Vectorized uint8[N] classification of each fragment.

        Returns
        -------
        np.ndarray
            ``FRAG_UNAMBIG`` (0): same-strand, 1 transcript, NH=1.
            ``FRAG_AMBIG_SAME_STRAND`` (1): same-strand, >1 transcript, NH=1.
            ``FRAG_AMBIG_OPP_STRAND`` (2): ambig-strand transcripts, NH=1.
            ``FRAG_MULTIMAPPER`` (3): NH > 1.
            ``FRAG_CHIMERIC`` (4): chimeric fragment.
        """
        if self._fragment_classes is not None:
            return self._fragment_classes
        n_transcripts = np.diff(self.t_offsets).astype(np.intp)
        classes = np.full(self.size, FRAG_UNAMBIG, dtype=np.uint8)
        # Same strand, multiple transcripts, single mapper → ambig same-strand
        classes[(self.ambig_strand == 0) & (n_transcripts > 1) & (self.num_hits == 1)] = (
            FRAG_AMBIG_SAME_STRAND
        )
        # Mixed strand, single mapper → ambig opposite-strand
        classes[(self.ambig_strand > 0) & (self.num_hits == 1)] = FRAG_AMBIG_OPP_STRAND
        # Multimapper → always FRAG_MULTIMAPPER (highest priority)
        classes[self.num_hits > 1] = FRAG_MULTIMAPPER
        # Chimeric → highest priority (overrides all others)
        classes[self.chimera_type > 0] = FRAG_CHIMERIC
        object.__setattr__(self, "_fragment_classes", classes)
        return classes

    def to_scoring_arrays(self) -> tuple:
        """Return the contiguous array tuple expected by ``StreamingScorer.score_chunk``.

        ``from_raw`` and ``_load_chunk`` both guarantee that every
        stored array is C-contiguous with the exact dtype the C++
        scoring kernel expects, so this method is a zero-copy view
        constructor (apart from the lazy ``fragment_classes`` compute
        on first access).
        """
        return (
            self.t_offsets,
            self.t_indices,
            self.frag_lengths,
            self.exon_bp,
            self.splice_type,
            self.align_strand,
            self.fragment_classes,
            self.frag_id,
            self.read_length,
            self.genomic_footprint,
            self.genomic_start,
            self.nm,
        )


# ---------------------------------------------------------------------------
# Arrow IPC (Feather v2) spill / load
# ---------------------------------------------------------------------------


def _spill_chunk(chunk: _FinalizedChunk, path: Path) -> None:
    """Write a finalized chunk to disk as Arrow IPC with LZ4."""
    import pyarrow as pa
    import pyarrow.feather as pf

    t_list = pa.ListArray.from_arrays(
        chunk.t_offsets,
        chunk.t_indices,
    )
    frag_lengths_list = pa.ListArray.from_arrays(
        chunk.t_offsets,
        chunk.frag_lengths,
    )
    exon_bp_list = pa.ListArray.from_arrays(
        chunk.t_offsets,
        chunk.exon_bp,
    )
    table = pa.table(
        {
            "splice_type": chunk.splice_type,
            "align_strand": chunk.align_strand,
            "num_hits": chunk.num_hits,
            "chimera_type": chunk.chimera_type,
            "t_inds": t_list,
            "frag_lengths": frag_lengths_list,
            "exon_bp": exon_bp_list,
            "ambig_strand": chunk.ambig_strand,
            "frag_id": chunk.frag_id,
            "read_length": chunk.read_length,
            "genomic_footprint": chunk.genomic_footprint,
            "genomic_start": chunk.genomic_start,
            "nm": chunk.nm,
        }
    )

    pf.write_feather(table, str(path), compression="lz4")


def _load_chunk(path: Path) -> _FinalizedChunk:
    """Load a spilled chunk from an Arrow IPC file."""
    import pyarrow.feather as pf

    table = pf.read_table(str(path))

    t_col = table.column("t_inds").combine_chunks()
    t_offsets = t_col.offsets.to_numpy().astype(np.int32)
    t_indices = t_col.values.to_numpy().copy()

    # Per-candidate CSR arrays (parallel to t_indices)
    frag_lengths_col = table.column("frag_lengths").combine_chunks()
    frag_lengths_arr = frag_lengths_col.values.to_numpy().astype(np.int32, copy=False)

    exon_bp_col = table.column("exon_bp").combine_chunks()
    exon_bp_arr = exon_bp_col.values.to_numpy().astype(np.uint16, copy=False)

    return _FinalizedChunk(
        splice_type=table.column("splice_type").to_numpy().copy(),
        align_strand=table.column("align_strand").to_numpy().copy(),
        num_hits=table.column("num_hits").to_numpy().copy(),
        chimera_type=table.column("chimera_type").to_numpy().copy(),
        t_offsets=t_offsets,
        t_indices=t_indices,
        frag_lengths=frag_lengths_arr,
        exon_bp=exon_bp_arr,
        ambig_strand=table.column("ambig_strand").to_numpy().copy(),
        frag_id=table.column("frag_id").to_numpy().astype(np.int64),
        read_length=table.column("read_length").to_numpy().copy().astype(np.uint16),
        genomic_footprint=table.column("genomic_footprint").to_numpy().copy(),
        genomic_start=table.column("genomic_start").to_numpy().copy(),
        nm=table.column("nm").to_numpy().copy().astype(np.uint16),
        size=len(table),
    )


@dataclass(slots=True)
class _PendingSpill:
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

    def __init__(self, spill_dir: Path | None, capacity: int = 2):
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

    def submit(self, chunk: _FinalizedChunk, pending: _PendingSpill) -> None:
        if self._closed:
            raise RuntimeError("Cannot submit to a closed spill writer.")
        self._slots.acquire()
        try:
            self._queue.put((chunk, pending))
        except BaseException:
            self._slots.release()
            raise

    def stop(self) -> None:
        """Queue the stop; the thread writes what was submitted before it, removes the directory and
        exits. Joins the thread, except on the thread itself, where garbage collection can run the
        buffer's finalizer."""
        if self._closed:
            return
        self._closed = True
        self._queue.put(self._STOP)
        if threading.current_thread() is not self._thread:
            self._thread.join()

    def _run(self) -> None:
        while True:
            item = self._queue.get()
            if item is self._STOP:
                shutil.rmtree(self.temp_dir, ignore_errors=True)
                return

            chunk, pending = item
            try:
                logger.info(
                    "Writing spilled chunk (%s fragments, %.1f MB) -> %s",
                    f"{pending.size:,}",
                    pending.memory_bytes / 1024**2,
                    pending.path,
                )
                _spill_chunk(chunk, pending.path)
                logger.info("Finished spilled chunk -> %s", pending.path)
            except BaseException as exc:
                pending.error = exc
                logger.exception("Failed to write spilled chunk -> %s", pending.path)
            finally:
                pending.done.set()
                self._slots.release()


# ---------------------------------------------------------------------------
# FragmentBuffer — public API
# ---------------------------------------------------------------------------


class FragmentBuffer:
    """Columnar buffer for resolved fragments with disk-spill support.

    Holds the chunks the native scanner hands over (:meth:`inject_chunk`).
    When total in-memory size exceeds *max_memory_bytes*, the oldest chunk
    is spilled to disk as an Arrow IPC (Feather v2) file with LZ4
    compression. :meth:`iter_chunks_consuming` yields the chunks in order.

    After quantification is complete, call :meth:`cleanup` (or use the
    context-manager protocol) to remove spilled files.

    Parameters
    ----------
    max_memory_bytes : int
        Maximum memory for in-memory chunks before spilling to disk.
        Default 2 GiB.  Set to 0 to disable spilling.
    spill_dir : Path or None
        Directory for spilled chunk files.  Default: system temp dir.
    """

    def __init__(
        self,
        max_memory_bytes: int = 2 * 1024**3,
        spill_dir: Path | None = None,
    ):
        self.max_memory_bytes = max_memory_bytes
        self._spill_dir = spill_dir

        self._chunks: deque[_FinalizedChunk | _PendingSpill] = deque()
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
    def n_chunks(self) -> int:
        """Number of finalized chunks (in-memory + spilled)."""
        return len(self._chunks)

    @property
    def n_spilled(self) -> int:
        """Number of chunks spilled to disk."""
        return self._n_spilled

    # -- Accepting chunks -----------------------------------------------------

    def inject_chunk(self, chunk: _FinalizedChunk) -> None:
        """Append a finalized chunk and spill if over memory budget."""
        self._raise_completed_spill_errors()
        self._total_size += chunk.size
        self._memory_bytes += chunk.memory_bytes
        self._chunks.append(chunk)

        # Spill if over memory budget
        if self.max_memory_bytes > 0:
            while self._memory_bytes > self.max_memory_bytes:
                if not self._spill_oldest():
                    break

    def _spill_oldest(self) -> bool:
        """Spill the oldest in-memory chunk to disk.  Return True if spilled."""
        for i, chunk in enumerate(self._chunks):
            if isinstance(chunk, _FinalizedChunk):
                writer = self._ensure_spill_writer()
                path = writer.temp_dir / f"chunk_{self._n_spilled:04d}.arrow"
                freed = chunk.memory_bytes
                pending = _PendingSpill(
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
            self._spill_writer = _SpillWriter(self._spill_dir, capacity=2)
            # Stops the writer and removes its directory if cleanup() is never called.
            self._spill_finalizer = weakref.finalize(self, self._spill_writer.stop)
        return self._spill_writer

    def _stop_spill_writer(self) -> None:
        if self._spill_finalizer is not None:
            self._spill_finalizer()  # runs writer.stop() now; a finalizer runs at most once
            self._spill_finalizer = None
            self._spill_writer = None

    def _spill_error(self, pending: _PendingSpill) -> RuntimeError:
        return RuntimeError(f"Failed to spill buffer chunk to {pending.path}")

    def _wait_pending_spill(self, pending: _PendingSpill) -> Path:
        pending.done.wait()
        if pending.error is not None:
            raise self._spill_error(pending) from pending.error
        return pending.path

    def _wait_all_pending_spills(self) -> None:
        failed: _PendingSpill | None = None
        for chunk_ref in list(self._chunks):
            if isinstance(chunk_ref, _PendingSpill):
                chunk_ref.done.wait()
                if chunk_ref.error is not None and failed is None:
                    failed = chunk_ref
        if failed is not None:
            raise self._spill_error(failed) from failed.error

    def _raise_completed_spill_errors(self) -> None:
        for chunk_ref in self._chunks:
            if (
                isinstance(chunk_ref, _PendingSpill)
                and chunk_ref.done.is_set()
                and chunk_ref.error is not None
            ):
                raise self._spill_error(chunk_ref) from chunk_ref.error

    # -- Iteration ------------------------------------------------------------

    def iter_chunks_consuming(self) -> Iterator[_FinalizedChunk]:
        """Yield chunks one at a time, releasing each after the caller advances.

        After this method returns, the buffer is empty.  Spilled chunk
        files are deleted as they are consumed.
        """
        while self._chunks:
            chunk_ref = self._chunks.popleft()
            if isinstance(chunk_ref, _PendingSpill):
                path = self._wait_pending_spill(chunk_ref)
                chunk = _load_chunk(path)
                path.unlink(missing_ok=True)
            else:
                chunk = chunk_ref
                self._memory_bytes -= chunk.memory_bytes
            yield chunk

    # -- Cleanup --------------------------------------------------------------

    def cleanup(self) -> None:
        """Remove any spilled chunk files from disk."""
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
