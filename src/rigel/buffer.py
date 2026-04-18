"""
rigel.buffer — Columnar fragment buffer with Arrow IPC disk spill.

Stores resolved fragment data in memory-efficient columnar NumPy arrays
using CSR (Compressed Sparse Row) layout for variable-length transcript
and gene index sets.  When memory usage exceeds a configurable threshold,
completed chunks are spilled to disk as Arrow IPC (Feather v2) files
with LZ4 compression.

Memory efficiency vs. Python objects::

    ResolvedFragment objects:   ~580-640 bytes per fragment
    Columnar buffer:           ~40-50 bytes per fragment  (15× reduction)

Architecture
------------
- Fragments are appended one at a time during the BAM scan.
- Accumulated into chunks of configurable size (default 1M fragments).
- Each chunk is finalized to compact NumPy arrays.
- When total in-memory chunk size exceeds *max_memory_bytes*,
  the oldest chunk is spilled to disk as Arrow IPC with LZ4.
- Phase 2/3 iterate the buffer, yielding lightweight
  ``BufferedFragment`` views backed by array slices.
"""

import logging
import shutil
import tempfile
import weakref
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np

from .native import FragmentAccumulator
from .splice import SpliceType
from .types import ChimeraType, Strand

logger = logging.getLogger(__name__)

__all__ = ["FragmentBuffer", "BufferedFragment"]

# Fragment classification constants used by fragment_classes property
FRAG_UNAMBIG: int = 0  # same-strand, 1 transcript, NH=1
FRAG_AMBIG_SAME_STRAND: int = 1  # same-strand, >1 transcript, NH=1
FRAG_AMBIG_OPP_STRAND: int = 2  # ambig-strand transcripts, NH=1
FRAG_MULTIMAPPER: int = 3  # NH > 1 (multimapped molecule)
FRAG_CHIMERIC: int = 4  # chimeric fragment (disjoint transcript sets)


# ---------------------------------------------------------------------------
# BufferedFragment — lightweight view into a finalized chunk
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class BufferedFragment:
    """Lightweight view into a columnar buffer chunk.

    Provides the same duck-typed interface as ``ResolvedFragment`` so
    that the C++ scoring functions work without modification.

    ``t_inds`` is a NumPy array slice (supports iteration, ``len()``,
    indexing) rather than a frozenset.  Strand mixing is represented
    by ``ambig_strand`` (derived from transcript strand array
    at buffer finalization time).
    """

    t_inds: np.ndarray
    ambig_strand: int
    splice_type: int
    exon_strand: int
    sj_strand: int
    frag_lengths: np.ndarray
    num_hits: int
    merge_criteria: int
    chimera_type: int = ChimeraType.NONE
    frag_id: int = 0
    read_length: int = 0
    genomic_footprint: int = -1
    genomic_start: int = -1
    nm: int = 0
    exon_bp: np.ndarray | None = None
    intron_bp: np.ndarray | None = None

    @property
    def is_chimeric(self) -> bool:
        return self.chimera_type != ChimeraType.NONE

    @property
    def is_same_strand(self) -> bool:
        return not self.ambig_strand

    @property
    def is_strand_qualified(self) -> bool:
        return (
            self.splice_type == SpliceType.SPLICED_ANNOT
            and self.is_same_strand
            and self.exon_strand in (Strand.POS, Strand.NEG)
            and self.sj_strand in (Strand.POS, Strand.NEG)
        )


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
    array, derived from transcript strand array at finalization time.

    ``frag_id`` is stored as **int64** to accommodate BAM files with
    more than 2 billion fragments.  ``t_offsets`` is **int32**, which
    limits the maximum number of candidates per chunk to ~2.1 billion
    — well above practical chunk sizes.
    """

    splice_type: np.ndarray  # uint8[N]
    exon_strand: np.ndarray  # uint8[N]
    sj_strand: np.ndarray  # uint8[N]
    num_hits: np.ndarray  # uint16[N]
    merge_criteria: np.ndarray  # uint8[N]
    chimera_type: np.ndarray  # uint8[N]
    t_offsets: np.ndarray  # int32[N+1]
    t_indices: np.ndarray  # int32[M_t]
    frag_lengths: np.ndarray  # int32[M_t]  (parallel to t_indices)
    exon_bp: np.ndarray  # int32[M_t]  (parallel to t_indices)
    intron_bp: np.ndarray  # int32[M_t]  (parallel to t_indices)
    ambig_strand: np.ndarray  # uint8[N]
    frag_id: np.ndarray  # int64[N]
    read_length: np.ndarray  # uint32[N]
    genomic_footprint: np.ndarray  # int32[N]
    genomic_start: np.ndarray  # int32[N]
    nm: np.ndarray  # uint16[N]
    size: int
    _fragment_classes: np.ndarray | None = None  # cached uint8[N]

    @classmethod
    def from_raw(cls, raw: dict) -> "_FinalizedChunk":
        """Build a chunk from the raw dict returned by C++ FragmentAccumulator.

        Accepts both legacy bytes (from ``finalize()``) and zero-copy
        numpy arrays (from ``finalize_zero_copy()``).  When the value
        is already an ndarray the data is used directly; ``t_offsets``
        is stored as int32 (native from C++) and ``frag_id`` is kept
        as int64.
        """

        def _arr(val, dtype, src_dtype=None):
            """Convert bytes or ndarray to a NumPy array of *dtype*.

            *src_dtype* overrides the dtype used when interpreting raw
            bytes (needed when the C++ source type differs from the
            desired Python storage type, e.g. int64 → int32).
            """
            if isinstance(val, np.ndarray):
                if val.dtype == dtype:
                    return np.ascontiguousarray(val)
                return val.astype(dtype, copy=False)
            # Legacy path: raw bytes from finalize()
            read_dtype = src_dtype if src_dtype is not None else dtype
            arr = np.frombuffer(val, dtype=read_dtype).copy()
            if read_dtype != dtype:
                return arr.astype(dtype)
            return arr

        return cls(
            splice_type=_arr(raw["splice_type"], np.uint8),
            exon_strand=_arr(raw["exon_strand"], np.uint8),
            sj_strand=_arr(raw["sj_strand"], np.uint8),
            num_hits=_arr(raw["num_hits"], np.uint16),
            merge_criteria=_arr(raw["merge_criteria"], np.uint8),
            chimera_type=_arr(raw["chimera_type"], np.uint8),
            t_offsets=_arr(raw["t_offsets"], np.int32),
            t_indices=_arr(raw["t_indices"], np.int32),
            frag_lengths=_arr(raw["frag_lengths"], np.int32),
            exon_bp=_arr(raw["exon_bp"], np.int32),
            intron_bp=_arr(raw["intron_bp"], np.int32),
            ambig_strand=_arr(raw["ambig_strand"], np.uint8),
            frag_id=_arr(raw["frag_id"], np.int64),
            read_length=_arr(raw["read_length"], np.uint32),
            genomic_footprint=_arr(raw["genomic_footprint"], np.int32),
            genomic_start=_arr(raw["genomic_start"], np.int32),
            nm=_arr(raw["nm"], np.uint16),
            size=raw["size"] if isinstance(raw["size"], int) else int(raw["size"]),
        )

    @property
    def memory_bytes(self) -> int:
        """Total bytes consumed by the underlying NumPy arrays."""
        return sum(
            a.nbytes
            for a in (
                self.splice_type,
                self.exon_strand,
                self.sj_strand,
                self.num_hits,
                self.merge_criteria,
                self.chimera_type,
                self.t_offsets,
                self.t_indices,
                self.frag_lengths,
                self.exon_bp,
                self.intron_bp,
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
        object.__setattr__(self, '_fragment_classes', classes)
        return classes

    def to_scoring_arrays(self) -> tuple:
        """Return the contiguous array tuple expected by ``StreamingScorer.score_chunk``.

        ``t_offsets`` is int32 (native) and ``frag_id`` is int64 (native),
        matching the C++ scoring kernel types directly without copying.
        Computes the derived ``fragment_classes`` column.
        All returned arrays are C-contiguous.
        """
        return (
            np.ascontiguousarray(self.t_offsets, dtype=np.int32),
            np.ascontiguousarray(self.t_indices, dtype=np.int32),
            np.ascontiguousarray(self.frag_lengths, dtype=np.int32),
            np.ascontiguousarray(self.exon_bp, dtype=np.int32),
            np.ascontiguousarray(self.intron_bp, dtype=np.int32),
            np.ascontiguousarray(self.splice_type, dtype=np.uint8),
            np.ascontiguousarray(self.exon_strand, dtype=np.uint8),
            np.ascontiguousarray(self.fragment_classes, dtype=np.uint8),
            np.ascontiguousarray(self.frag_id, dtype=np.int64),
            np.ascontiguousarray(self.read_length, dtype=np.uint32),
            np.ascontiguousarray(self.genomic_footprint, dtype=np.int32),
            np.ascontiguousarray(self.genomic_start, dtype=np.int32),
            np.ascontiguousarray(self.nm, dtype=np.uint16),
        )

    def __len__(self) -> int:
        return self.size

    def __getitem__(self, i: int) -> BufferedFragment:
        """Return a lightweight view for fragment *i*."""
        start = self.t_offsets[i]
        end = self.t_offsets[i + 1]
        return BufferedFragment(
            t_inds=self.t_indices[start:end],
            frag_lengths=self.frag_lengths[start:end],
            exon_bp=self.exon_bp[start:end],
            intron_bp=self.intron_bp[start:end],
            ambig_strand=int(self.ambig_strand[i]),
            splice_type=int(self.splice_type[i]),
            exon_strand=int(self.exon_strand[i]),
            sj_strand=int(self.sj_strand[i]),
            num_hits=int(self.num_hits[i]),
            merge_criteria=int(self.merge_criteria[i]),
            chimera_type=int(self.chimera_type[i]),
            frag_id=int(self.frag_id[i]),
            read_length=int(self.read_length[i]),
            genomic_footprint=int(self.genomic_footprint[i]),
            genomic_start=int(self.genomic_start[i]),
            nm=int(self.nm[i]),
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
    intron_bp_list = pa.ListArray.from_arrays(
        chunk.t_offsets,
        chunk.intron_bp,
    )

    table = pa.table(
        {
            "splice_type": chunk.splice_type,
            "exon_strand": chunk.exon_strand,
            "sj_strand": chunk.sj_strand,
            "num_hits": chunk.num_hits,
            "merge_criteria": chunk.merge_criteria,
            "chimera_type": chunk.chimera_type,
            "t_inds": t_list,
            "frag_lengths": frag_lengths_list,
            "exon_bp": exon_bp_list,
            "intron_bp": intron_bp_list,
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
    frag_lengths_arr = frag_lengths_col.values.to_numpy().astype(np.int32)

    exon_bp_col = table.column("exon_bp").combine_chunks()
    exon_bp_arr = exon_bp_col.values.to_numpy().astype(np.int32)

    intron_bp_col = table.column("intron_bp").combine_chunks()
    intron_bp_arr = intron_bp_col.values.to_numpy().astype(np.int32)

    return _FinalizedChunk(
        splice_type=table.column("splice_type").to_numpy().copy(),
        exon_strand=table.column("exon_strand").to_numpy().copy(),
        sj_strand=table.column("sj_strand").to_numpy().copy(),
        num_hits=table.column("num_hits").to_numpy().copy(),
        merge_criteria=table.column("merge_criteria").to_numpy().copy(),
        chimera_type=table.column("chimera_type").to_numpy().copy(),
        t_offsets=t_offsets,
        t_indices=t_indices,
        frag_lengths=frag_lengths_arr,
        exon_bp=exon_bp_arr,
        intron_bp=intron_bp_arr,
        ambig_strand=table.column("ambig_strand").to_numpy().copy(),
        frag_id=table.column("frag_id").to_numpy().astype(np.int64),
        read_length=table.column("read_length").to_numpy().copy(),
        genomic_footprint=table.column("genomic_footprint").to_numpy().copy(),
        genomic_start=(
            table.column("genomic_start").to_numpy().copy()
            if "genomic_start" in table.column_names
            else np.full(len(table), -1, dtype=np.int32)
        ),
        nm=table.column("nm").to_numpy().copy().astype(np.uint16),
        size=len(table),
    )


# ---------------------------------------------------------------------------
# FragmentBuffer — public API
# ---------------------------------------------------------------------------


class FragmentBuffer:
    """Columnar buffer for resolved fragments with disk-spill support.

    Accumulates ``ResolvedFragment`` data during the BAM scan into
    compact NumPy chunks.  When total in-memory size exceeds
    *max_memory_bytes*, the oldest chunk is spilled to disk as an
    Arrow IPC (Feather v2) file with LZ4 compression.

    After all fragments are appended, call :meth:`finalize` to flush
    any remaining data.  Then iterate with :meth:`__iter__` (yields
    ``BufferedFragment``) or :meth:`iter_chunks` (yields
    ``_FinalizedChunk``).

    After quantification is complete, call :meth:`cleanup` (or use the
    context-manager protocol) to remove spilled files.

    Parameters
    ----------
    t_strand_arr : np.ndarray
        Per-transcript strand array ``int8[n_transcripts]``, used to
        compute per-fragment ``ambig_strand`` when finalizing chunks.
    chunk_size : int
        Number of fragments per chunk (default 1,000,000).
    max_memory_bytes : int
        Maximum memory for in-memory chunks before spilling to disk.
        Default 2 GiB.  Set to 0 to disable spilling.
    spill_dir : Path or None
        Directory for spilled chunk files.  Default: system temp dir.

    Examples
    --------
    >>> buf = FragmentBuffer(t_strand_arr, chunk_size=500_000)
    >>> for resolved in resolved_fragments:
    ...     buf.append(resolved)
    >>> buf.finalize()
    >>> for frag in buf:
    ...     print(frag.splice_type, len(frag.t_inds))
    >>> buf.cleanup()
    """

    def __init__(
        self,
        t_strand_arr: np.ndarray,
        chunk_size: int = 1_000_000,
        max_memory_bytes: int = 2 * 1024**3,
        spill_dir: Path | None = None,
    ):
        self._t_strand_arr = t_strand_arr
        self.chunk_size = chunk_size
        self.max_memory_bytes = max_memory_bytes
        self._spill_dir = spill_dir
        self._temp_dir: str | None = None

        self._native_acc = FragmentAccumulator()
        self._chunks: deque[_FinalizedChunk | Path] = deque()
        self._total_size = 0
        self._memory_bytes = 0
        self._n_spilled = 0

        # Safety net: ensure spill directory is cleaned up even if
        # cleanup() is never called.  weakref.finalize is reliable
        # unlike __del__ (guaranteed ordering, runs during GC).
        self._ensure_cleanup = weakref.finalize(self, FragmentBuffer._weak_cleanup, self)

    # -- Properties -----------------------------------------------------------

    @property
    def total_fragments(self) -> int:
        """Total buffered fragments (finalized + pending)."""
        return self._total_size + self._native_acc.size

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

    # -- Append ---------------------------------------------------------------

    def append(self, resolved, frag_id: int = 0) -> None:
        """Append a resolved fragment to the buffer.

        Parameters
        ----------
        resolved : ResolvedFragment
            The resolved fragment to buffer (C++ native result).
        frag_id : int
            Fragment identifier for grouping multimapper alignments.
            All alignments of the same molecule should share a frag_id.
        """
        self._native_acc.append(resolved, frag_id)
        if self._native_acc.size >= self.chunk_size:
            self._finalize_native()

    # -- Finalization ---------------------------------------------------------

    def finalize(self) -> None:
        """Finalize any remaining accumulated data.

        Must be called after all fragments have been appended and
        before iterating.
        """
        self._finalize_native()

    def _finalize_native(self) -> None:
        """Finalize the C++ native accumulator into a chunk."""
        if self._native_acc.size == 0:
            return

        t_strand_list = self._t_strand_arr.tolist()
        raw = self._native_acc.finalize(t_strand_list)
        chunk = _FinalizedChunk.from_raw(raw)

        self._accept_chunk(chunk)
        self._native_acc = FragmentAccumulator()

    def _accept_chunk(self, chunk: _FinalizedChunk) -> None:
        """Append a finalized chunk and spill if over memory budget."""
        self._total_size += chunk.size
        self._memory_bytes += chunk.memory_bytes
        self._chunks.append(chunk)

        # Spill if over memory budget
        if self.max_memory_bytes > 0:
            while self._memory_bytes > self.max_memory_bytes:
                if not self._spill_oldest():
                    break

    def inject_chunk(self, chunk: _FinalizedChunk) -> None:
        """Inject an externally-built chunk into the buffer."""
        self._accept_chunk(chunk)

    def release(self) -> None:
        """Release all in-memory chunks and spilled files."""
        self.cleanup()
        self._chunks.clear()
        self._memory_bytes = 0

    def _spill_oldest(self) -> bool:
        """Spill the oldest in-memory chunk to disk.  Return True if spilled."""
        for i, chunk in enumerate(self._chunks):
            if isinstance(chunk, _FinalizedChunk):
                path = self._get_spill_path()
                _spill_chunk(chunk, path)
                freed = chunk.memory_bytes
                self._chunks[i] = path
                self._memory_bytes -= freed
                self._n_spilled += 1
                logger.info(
                    f"Spilled chunk {i} ({chunk.size:,} fragments, "
                    f"{freed / 1024**2:.1f} MB) → {path}"
                )
                return True
        return False

    def _get_spill_path(self) -> Path:
        """Return a unique path for a spilled chunk file."""
        if self._temp_dir is None:
            if self._spill_dir is not None:
                self._spill_dir.mkdir(parents=True, exist_ok=True)
                self._temp_dir = tempfile.mkdtemp(dir=self._spill_dir, prefix="rigel_buf_")
            else:
                self._temp_dir = tempfile.mkdtemp(prefix="rigel_buf_")
        idx = self._n_spilled
        return Path(self._temp_dir) / f"chunk_{idx:04d}.arrow"

    # -- Iteration ------------------------------------------------------------

    def __len__(self) -> int:
        return self.total_fragments

    def __iter__(self) -> Iterator[BufferedFragment]:
        """Iterate over all buffered fragments.

        Yields lightweight ``BufferedFragment`` views.  Spilled chunks
        are loaded from disk one at a time so only one extra chunk is
        in memory at any moment.
        """
        for chunk in self.iter_chunks():
            for i in range(chunk.size):
                yield chunk[i]

    def iter_chunks(self) -> Iterator[_FinalizedChunk]:
        """Iterate over finalized chunks.

        In-memory chunks are yielded directly.  Spilled chunks are
        loaded from disk on demand.
        """
        for chunk_ref in self._chunks:
            if isinstance(chunk_ref, Path):
                yield _load_chunk(chunk_ref)
            else:
                yield chunk_ref

    def iter_chunks_consuming(self) -> Iterator[_FinalizedChunk]:
        """Yield chunks one at a time, releasing each after the caller advances.

        After this method returns, the buffer is empty.  Spilled chunk
        files are deleted as they are consumed.

        This is the streaming counterpart to :meth:`iter_chunks`.  Use
        it when chunks will not be revisited — the typical case during
        scoring, where each chunk is scored exactly once.
        """
        while self._chunks:
            chunk_ref = self._chunks.popleft()
            if isinstance(chunk_ref, Path):
                chunk = _load_chunk(chunk_ref)
                chunk_ref.unlink(missing_ok=True)
            else:
                chunk = chunk_ref
                self._memory_bytes -= chunk.memory_bytes
            yield chunk

    # -- Cleanup --------------------------------------------------------------

    def cleanup(self) -> None:
        """Remove any spilled chunk files from disk."""
        if self._temp_dir is not None:
            shutil.rmtree(self._temp_dir, ignore_errors=True)
            self._temp_dir = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.cleanup()

    @staticmethod
    def _weak_cleanup(buf: "FragmentBuffer") -> None:
        """Release callback for weakref.finalize."""
        if buf._temp_dir is not None:
            shutil.rmtree(buf._temp_dir, ignore_errors=True)
            buf._temp_dir = None

    # -- Diagnostics ----------------------------------------------------------

    def summary(self) -> dict:
        """Return a JSON-serializable summary of buffer state."""
        n_mem = sum(1 for c in self._chunks if isinstance(c, _FinalizedChunk))
        n_disk = sum(1 for c in self._chunks if isinstance(c, Path))
        return {
            "total_fragments": self.total_fragments,
            "n_chunks": len(self._chunks),
            "in_memory_chunks": n_mem,
            "on_disk_chunks": n_disk,
            "memory_bytes": self._memory_bytes,
            "memory_mb": round(self._memory_bytes / 1024**2, 1),
        }
