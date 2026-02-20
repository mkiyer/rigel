"""
hulkrna.buffer — Columnar fragment buffer with Arrow IPC disk spill.

Stores resolved fragment data in memory-efficient columnar NumPy arrays
using CSR (Compressed Sparse Row) layout for variable-length transcript
and gene index sets.  When memory usage exceeds a configurable threshold,
completed chunks are spilled to disk as Arrow IPC (Feather v2) files
with LZ4 compression.

Memory efficiency vs. Python objects::

    ResolvedFragment objects:  ~580-640 bytes per fragment
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
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np

from .categories import SpliceType
from .types import ChimeraType, Strand

logger = logging.getLogger(__name__)

__all__ = ["FragmentBuffer", "BufferedFragment"]

# Fragment classification constants used by fragment_classes property
FRAG_UNIQUE: int = 0          # 1 gene, 1 transcript, NH=1
FRAG_ISOFORM_AMBIG: int = 1   # 1 gene, >1 transcript, NH=1
FRAG_GENE_AMBIG: int = 2      # >1 gene, NH=1
FRAG_MULTIMAPPER: int = 3     # NH > 1 (multimapped molecule)
FRAG_CHIMERIC: int = 4        # chimeric fragment (disjoint transcript sets)


# ---------------------------------------------------------------------------
# BufferedFragment — lightweight view into a finalized chunk
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class BufferedFragment:
    """Lightweight view into a columnar buffer chunk.

    Provides the same duck-typed interface as ``ResolvedFragment`` so
    that ``ReadCounter.assign_unique / assign_isoform_bayesian / assign_gene_bayesian``
    work without modification.

    ``t_inds`` is a NumPy array slice (supports iteration, ``len()``,
    indexing) rather than a frozenset.  Gene information is represented
    by ``n_genes`` (derived from transcript indices via ``t_to_g_arr``
    at buffer finalization time).
    """

    t_inds: np.ndarray
    n_genes: int
    splice_type: int
    exon_strand: int
    sj_strand: int
    insert_size: int
    num_hits: int
    merge_criteria: int
    chimera_type: int = ChimeraType.NONE
    frag_id: int = 0
    frag_length: int = 0
    frag_start: int = 0
    frag_end: int = 0
    exon_bp: np.ndarray | None = None
    intron_bp: np.ndarray | None = None

    @property
    def is_chimeric(self) -> bool:
        return self.chimera_type != ChimeraType.NONE

    @property
    def is_unique_gene(self) -> bool:
        return self.n_genes == 1

    @property
    def is_strand_qualified(self) -> bool:
        return (
            self.splice_type == SpliceType.SPLICED_ANNOT
            and self.is_unique_gene
            and self.exon_strand in (Strand.POS, Strand.NEG)
            and self.sj_strand in (Strand.POS, Strand.NEG)
        )


# ---------------------------------------------------------------------------
# _AccumulatorChunk — mutable append buffer (Python lists)
# ---------------------------------------------------------------------------

class _AccumulatorChunk:
    """Mutable accumulator for one chunk of fragment data.

    Uses plain Python lists for O(1) amortized append.  Converted to
    compact NumPy arrays by ``_finalize()``.
    """

    __slots__ = (
        "splice_type", "exon_strand", "sj_strand", "insert_size",
        "num_hits", "merge_criteria", "chimera_type",
        "t_indices", "t_offsets",
        "exon_bp_flat", "intron_bp_flat",
        "frag_id", "frag_length",
        "frag_start", "frag_end",
        "size",
    )

    def __init__(self):
        self.splice_type: list[int] = []
        self.exon_strand: list[int] = []
        self.sj_strand: list[int] = []
        self.insert_size: list[int] = []
        self.num_hits: list[int] = []
        self.merge_criteria: list[int] = []
        self.chimera_type: list[int] = []
        self.t_indices: list[int] = []
        self.t_offsets: list[int] = [0]
        self.exon_bp_flat: list[int] = []
        self.intron_bp_flat: list[int] = []
        self.frag_id: list[int] = []
        self.frag_length: list[int] = []
        self.frag_start: list[int] = []
        self.frag_end: list[int] = []
        self.size: int = 0

    def append(self, resolved, frag_id: int = 0) -> None:
        """Append a ResolvedFragment's data to the accumulator."""
        self.splice_type.append(int(resolved.splice_type))
        self.exon_strand.append(int(resolved.exon_strand))
        self.sj_strand.append(int(resolved.sj_strand))
        self.insert_size.append(resolved.insert_size)
        self.num_hits.append(resolved.num_hits)
        self.merge_criteria.append(int(resolved.merge_criteria))
        self.chimera_type.append(int(getattr(resolved, 'chimera_type', 0)))

        # Overlap profiles: stored parallel to t_indices in CSR layout.
        # Lookup each t_idx in the resolved overlap_bp dict;
        # default (frag_length, 0) for exon_bp and intron_bp.
        overlap_bp = getattr(resolved, 'overlap_bp', None) or {}
        fl = getattr(resolved, 'frag_length', 0)
        for t_idx in resolved.t_inds:
            self.t_indices.append(t_idx)
            t_key = int(t_idx) if not isinstance(t_idx, int) else t_idx
            profile = overlap_bp.get(t_key)
            if profile is not None:
                self.exon_bp_flat.append(profile[0])
                self.intron_bp_flat.append(profile[1])
            else:
                self.exon_bp_flat.append(fl)
                self.intron_bp_flat.append(0)
        self.t_offsets.append(len(self.t_indices))

        self.frag_id.append(frag_id)
        self.frag_length.append(fl)
        self.frag_start.append(getattr(resolved, 'frag_start', 0))
        self.frag_end.append(getattr(resolved, 'frag_end', 0))
        self.size += 1


# ---------------------------------------------------------------------------
# _FinalizedChunk — immutable chunk in compact NumPy arrays
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class _FinalizedChunk:
    """Immutable chunk of fragment data in compact columnar arrays.

    Fixed-width fields are stored as typed NumPy vectors.
    Variable-width transcript index sets use CSR (Compressed Sparse Row)
    layout: ``offsets[i]:offsets[i+1]`` indexes the flat ``indices`` array.

    Gene count per fragment (``n_genes``) is cached as a uint8 array,
    derived from transcript indices via ``t_to_g_arr`` at finalization
    time.  This replaces the former gene CSR storage (~12 bytes/fragment
    savings, ~30% reduction).
    """

    splice_type: np.ndarray       # uint8[N]
    exon_strand: np.ndarray     # uint8[N]
    sj_strand: np.ndarray       # uint8[N]
    insert_size: np.ndarray     # int32[N]
    num_hits: np.ndarray        # uint16[N]
    merge_criteria: np.ndarray  # uint8[N]
    chimera_type: np.ndarray    # uint8[N]
    t_offsets: np.ndarray       # int64[N+1]
    t_indices: np.ndarray       # int32[M_t]
    exon_bp: np.ndarray         # int32[M_t]  (parallel to t_indices)
    intron_bp: np.ndarray       # int32[M_t]  (parallel to t_indices)
    n_genes: np.ndarray         # uint8[N]
    frag_id: np.ndarray         # int64[N]
    frag_length: np.ndarray     # uint32[N]
    frag_start: np.ndarray      # int32[N]
    frag_end: np.ndarray        # int32[N]
    size: int

    @property
    def memory_bytes(self) -> int:
        """Total bytes consumed by the underlying NumPy arrays."""
        return sum(
            a.nbytes for a in (
                self.splice_type, self.exon_strand, self.sj_strand,
                self.insert_size, self.num_hits, self.merge_criteria,
                self.chimera_type,
                self.t_offsets, self.t_indices,
                self.exon_bp, self.intron_bp,
                self.n_genes, self.frag_id, self.frag_length,
                self.frag_start, self.frag_end,
            )
        )

    @property
    def ambiguous_mask(self) -> np.ndarray:
        """Vectorized bool[N] mask: True where fragment is ambiguous."""
        return (self.n_genes > 1) | (self.num_hits > 1)

    @property
    def fragment_classes(self) -> np.ndarray:
        """Vectorized uint8[N] classification of each fragment.

        Returns
        -------
        np.ndarray
            ``FRAG_UNIQUE`` (0): 1 gene, 1 transcript, NH=1.
            ``FRAG_ISOFORM_AMBIG`` (1): 1 gene, >1 transcript, NH=1.
            ``FRAG_GENE_AMBIG`` (2): >1 gene, NH=1.
            ``FRAG_MULTIMAPPER`` (3): NH > 1.
            ``FRAG_CHIMERIC`` (4): chimeric fragment.
        """
        n_transcripts = np.diff(self.t_offsets).astype(np.intp)
        classes = np.full(self.size, FRAG_UNIQUE, dtype=np.uint8)
        # Single gene, multiple transcripts, single mapper → isoform-ambiguous
        classes[
            (self.n_genes == 1) & (n_transcripts > 1) & (self.num_hits == 1)
        ] = FRAG_ISOFORM_AMBIG
        # Inter-gene, single mapper → gene-ambiguous
        classes[(self.n_genes > 1) & (self.num_hits == 1)] = FRAG_GENE_AMBIG
        # Multimapper → always FRAG_MULTIMAPPER (highest priority)
        classes[self.num_hits > 1] = FRAG_MULTIMAPPER
        # Chimeric → highest priority (overrides all others)
        classes[self.chimera_type > 0] = FRAG_CHIMERIC
        return classes

    def __len__(self) -> int:
        return self.size

    def __getitem__(self, i: int) -> BufferedFragment:
        """Return a lightweight view for fragment *i*."""
        start = self.t_offsets[i]
        end = self.t_offsets[i + 1]
        return BufferedFragment(
            t_inds=self.t_indices[start:end],
            exon_bp=self.exon_bp[start:end],
            intron_bp=self.intron_bp[start:end],
            n_genes=int(self.n_genes[i]),
            splice_type=int(self.splice_type[i]),
            exon_strand=int(self.exon_strand[i]),
            sj_strand=int(self.sj_strand[i]),
            insert_size=int(self.insert_size[i]),
            num_hits=int(self.num_hits[i]),
            merge_criteria=int(self.merge_criteria[i]),
            chimera_type=int(self.chimera_type[i]),
            frag_id=int(self.frag_id[i]),
            frag_length=int(self.frag_length[i]),
            frag_start=int(self.frag_start[i]),
            frag_end=int(self.frag_end[i]),
        )


# ---------------------------------------------------------------------------
# Chunk finalization (lists → numpy)
# ---------------------------------------------------------------------------

def _finalize(acc: _AccumulatorChunk, t_to_g_arr: np.ndarray) -> _FinalizedChunk:
    """Convert a mutable accumulator to immutable NumPy arrays.

    Parameters
    ----------
    acc : _AccumulatorChunk
        The accumulator with appended fragment data.
    t_to_g_arr : np.ndarray
        Mapping array ``t_index → g_index``, used to compute
        the per-fragment gene count (``n_genes``).
    """
    t_offsets = np.array(acc.t_offsets, dtype=np.int64)
    t_indices = np.array(acc.t_indices, dtype=np.int32)
    exon_bp_arr = np.array(acc.exon_bp_flat, dtype=np.int32)
    intron_bp_arr = np.array(acc.intron_bp_flat, dtype=np.int32)

    # Compute n_genes per fragment from transcript indices
    n_genes_list = np.empty(acc.size, dtype=np.uint8)
    for i in range(acc.size):
        start = t_offsets[i]
        end = t_offsets[i + 1]
        if start == end:
            n_genes_list[i] = 0
        else:
            t_slice = t_indices[start:end]
            n_genes_list[i] = len(np.unique(t_to_g_arr[t_slice]))

    return _FinalizedChunk(
        splice_type=np.array(acc.splice_type, dtype=np.uint8),
        exon_strand=np.array(acc.exon_strand, dtype=np.uint8),
        sj_strand=np.array(acc.sj_strand, dtype=np.uint8),
        insert_size=np.array(acc.insert_size, dtype=np.int32),
        num_hits=np.array(acc.num_hits, dtype=np.uint16),
        merge_criteria=np.array(acc.merge_criteria, dtype=np.uint8),
        chimera_type=np.array(acc.chimera_type, dtype=np.uint8),
        t_offsets=t_offsets,
        t_indices=t_indices,
        exon_bp=exon_bp_arr,
        intron_bp=intron_bp_arr,
        n_genes=n_genes_list,
        frag_id=np.array(acc.frag_id, dtype=np.int64),
        frag_length=np.array(acc.frag_length, dtype=np.uint32),
        frag_start=np.array(acc.frag_start, dtype=np.int32),
        frag_end=np.array(acc.frag_end, dtype=np.int32),
        size=acc.size,
    )


# ---------------------------------------------------------------------------
# Arrow IPC (Feather v2) spill / load
# ---------------------------------------------------------------------------

def _spill_chunk(chunk: _FinalizedChunk, path: Path) -> None:
    """Write a finalized chunk to disk as Arrow IPC with LZ4."""
    import pyarrow as pa
    import pyarrow.feather as pf

    t_list = pa.ListArray.from_arrays(
        chunk.t_offsets.astype(np.int32), chunk.t_indices,
    )
    exon_bp_list = pa.ListArray.from_arrays(
        chunk.t_offsets.astype(np.int32), chunk.exon_bp,
    )
    intron_bp_list = pa.ListArray.from_arrays(
        chunk.t_offsets.astype(np.int32), chunk.intron_bp,
    )

    table = pa.table({
        "splice_type": chunk.splice_type,
        "exon_strand": chunk.exon_strand,
        "sj_strand": chunk.sj_strand,
        "insert_size": chunk.insert_size,
        "num_hits": chunk.num_hits,
        "merge_criteria": chunk.merge_criteria,
        "chimera_type": chunk.chimera_type,
        "t_inds": t_list,
        "exon_bp": exon_bp_list,
        "intron_bp": intron_bp_list,
        "n_genes": chunk.n_genes,
        "frag_id": chunk.frag_id,
        "frag_length": chunk.frag_length,
        "frag_start": chunk.frag_start,
        "frag_end": chunk.frag_end,
    })

    pf.write_feather(table, str(path), compression="lz4")


def _load_chunk(path: Path) -> _FinalizedChunk:
    """Load a spilled chunk from an Arrow IPC file."""
    import pyarrow.feather as pf

    table = pf.read_table(str(path))

    t_col = table.column("t_inds").combine_chunks()

    # Overlap profiles: parallel CSR arrays of BP counts.
    # Legacy files without these columns default to fill values.
    n_vals = len(t_col.values)
    if "exon_bp" in table.column_names:
        exon_col = table.column("exon_bp").combine_chunks()
        exon_bp_arr = exon_col.values.to_numpy().astype(np.int32)
    elif "exon_fracs" in table.column_names:
        # Legacy format: exon_fracs float32 → approximate exon_bp
        # Can't recover exact BP without frag_length; use 0 as safe default
        exon_bp_arr = np.zeros(n_vals, dtype=np.int32)
    else:
        exon_bp_arr = np.zeros(n_vals, dtype=np.int32)

    if "intron_bp" in table.column_names:
        intron_col = table.column("intron_bp").combine_chunks()
        intron_bp_arr = intron_col.values.to_numpy().astype(np.int32)
    elif "genic_fracs" in table.column_names:
        intron_bp_arr = np.zeros(n_vals, dtype=np.int32)
    else:
        intron_bp_arr = np.zeros(n_vals, dtype=np.int32)

    # frag_length: per-fragment field (not per-candidate)
    if "frag_length" in table.column_names:
        frag_length_arr = table.column("frag_length").to_numpy().copy()
    else:
        frag_length_arr = np.zeros(len(table), dtype=np.uint32)

    # frag_start / frag_end: genomic coordinates for coverage islands
    if "frag_start" in table.column_names:
        frag_start_arr = table.column("frag_start").to_numpy().copy()
        frag_end_arr = table.column("frag_end").to_numpy().copy()
    else:
        frag_start_arr = np.zeros(len(table), dtype=np.int32)
        frag_end_arr = np.zeros(len(table), dtype=np.int32)

    return _FinalizedChunk(
        splice_type=table.column("splice_type").to_numpy().copy(),
        exon_strand=table.column("exon_strand").to_numpy().copy(),
        sj_strand=table.column("sj_strand").to_numpy().copy(),
        insert_size=table.column("insert_size").to_numpy().copy(),
        num_hits=table.column("num_hits").to_numpy().copy(),
        merge_criteria=table.column("merge_criteria").to_numpy().copy(),
        chimera_type=table.column("chimera_type").to_numpy().copy(),
        t_offsets=t_col.offsets.to_numpy().astype(np.int64),
        t_indices=t_col.values.to_numpy().copy(),
        exon_bp=exon_bp_arr,
        intron_bp=intron_bp_arr,
        n_genes=table.column("n_genes").to_numpy().copy(),
        frag_id=table.column("frag_id").to_numpy().copy(),
        frag_length=frag_length_arr,
        frag_start=frag_start_arr,
        frag_end=frag_end_arr,
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

    After counting is complete, call :meth:`cleanup` (or use the
    context-manager protocol) to remove spilled files.

    Parameters
    ----------
    t_to_g_arr : np.ndarray
        Mapping array ``t_index → g_index``, used to compute
        per-fragment gene counts when finalizing chunks.
    chunk_size : int
        Number of fragments per chunk (default 1,000,000).
    max_memory_bytes : int
        Maximum memory for in-memory chunks before spilling to disk.
        Default 2 GiB.  Set to 0 to disable spilling.
    spill_dir : Path or None
        Directory for spilled chunk files.  Default: system temp dir.

    Examples
    --------
    >>> buf = FragmentBuffer(t_to_g_arr, chunk_size=500_000)
    >>> for resolved in resolved_fragments:
    ...     buf.append(resolved)
    >>> buf.finalize()
    >>> for frag in buf:
    ...     print(frag.splice_type, len(frag.t_inds))
    >>> buf.cleanup()
    """

    def __init__(
        self,
        t_to_g_arr: np.ndarray,
        chunk_size: int = 1_000_000,
        max_memory_bytes: int = 2 * 1024**3,
        spill_dir: Path | None = None,
    ):
        self._t_to_g_arr = t_to_g_arr
        self.chunk_size = chunk_size
        self.max_memory_bytes = max_memory_bytes
        self._spill_dir = spill_dir
        self._temp_dir: str | None = None

        self._current = _AccumulatorChunk()
        self._chunks: list[_FinalizedChunk | Path] = []
        self._total_size = 0
        self._memory_bytes = 0
        self._n_spilled = 0

    # -- Properties -----------------------------------------------------------

    @property
    def total_fragments(self) -> int:
        """Total buffered fragments (finalized + pending)."""
        return self._total_size + self._current.size

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
            The resolved fragment to buffer.
        frag_id : int
            Fragment identifier for grouping multimapper alignments.
            All alignments of the same molecule should share a frag_id.
        """
        self._current.append(resolved, frag_id=frag_id)
        if self._current.size >= self.chunk_size:
            self._finalize_current()

    # -- Finalization ---------------------------------------------------------

    def finalize(self) -> None:
        """Finalize any remaining accumulated data.

        Must be called after all fragments have been appended and
        before iterating.
        """
        self._finalize_current()

    def _finalize_current(self) -> None:
        """Finalize the current accumulator into a chunk."""
        if self._current.size == 0:
            return

        chunk = _finalize(self._current, self._t_to_g_arr)
        self._total_size += chunk.size
        self._memory_bytes += chunk.memory_bytes
        self._chunks.append(chunk)
        self._current = _AccumulatorChunk()

        # Spill if over memory budget
        if self.max_memory_bytes > 0:
            while self._memory_bytes > self.max_memory_bytes:
                if not self._spill_oldest():
                    break  # nothing left to spill

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
                logger.debug(
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
                self._temp_dir = tempfile.mkdtemp(
                    dir=self._spill_dir, prefix="hulkrna_buf_"
                )
            else:
                self._temp_dir = tempfile.mkdtemp(prefix="hulkrna_buf_")
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

    def __del__(self):
        self.cleanup()

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
