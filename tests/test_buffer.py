"""`rigel.buffer` — the fragment buffer the EM reads, over the native accumulator.

The accumulator behind it, the chunks the buffer takes and hands back, `frag_id` assignment, the
fragment classes, spilling to disk and reading back, the finalizer, and the resolved-fragment surface.
Every chunk is built the way the scanner builds it — C++ `ResolvedFragment`s from
`FragmentResolver.resolve_fragment`, appended to the native `FragmentAccumulator`, finalized and handed
to `FragmentBuffer.inject_chunk` — so these exercise the real native path rather than a Python stand-in
that could diverge from it.
"""

import gc
import threading
import weakref

import numpy as np
import pytest

import rigel.buffer as buffer_mod
from rigel.native import FragmentAccumulator
from rigel.types import Strand, GenomicInterval
from rigel.splice import SpliceType
from _resolution_reference import make_fragment, resolve_fragment
from rigel.buffer import (
    FragmentBuffer,
    FRAG_AMBIG_SAME_STRAND,
    FRAG_MULTIMAPPER,
    _FinalizedChunk,
)


# =====================================================================
# Fixtures -- use shared mini_index from conftest.py
# =====================================================================


def _resolve(index, exons, introns=()):
    """Create a Fragment and resolve it via C++ FragmentResolver.

    Returns the C++ ResolvedFragment (or None for intergenic).
    """
    frag = make_fragment(exons=tuple(exons), introns=tuple(introns))
    return resolve_fragment(frag, index)


def _exon(ref, start, end, strand=Strand.POS):
    """Shorthand for creating a GenomicInterval."""
    return GenomicInterval(ref, start, end, strand)


def _inject(buf, fragments, chunk_size=100, frag_ids=None):
    """Hand ``fragments`` to ``buf`` the way the scanner does: the native accumulator finalizes every
    ``chunk_size`` of them into a chunk. ``frag_ids`` defaults to each fragment's position."""
    if frag_ids is None:
        frag_ids = range(len(fragments))
    frag_ids = list(frag_ids)
    for lo in range(0, len(fragments), chunk_size):
        acc = FragmentAccumulator()
        for resolved, frag_id in zip(
            fragments[lo : lo + chunk_size], frag_ids[lo : lo + chunk_size]
        ):
            acc.append(resolved, frag_id)
        buf.inject_chunk(_FinalizedChunk.from_raw(acc.finalize()))


def _consume(buf):
    return list(buf.iter_chunks_consuming())


def _only_chunk(buf):
    (chunk,) = _consume(buf)
    return chunk


def _spill_dirs(path):
    return [p for p in path.iterdir() if p.is_dir()]


_CHUNK_ARRAY_FIELDS = (
    "splice_type",
    "align_strand",
    "num_hits",
    "chimera_type",
    "t_offsets",
    "t_indices",
    "frag_lengths",
    "exon_bp",
    "ambig_strand",
    "frag_id",
    "read_length",
    "genomic_footprint",
    "genomic_start",
    "nm",
)


def _chunk_payloads(buf):
    return [
        (chunk.size, {name: getattr(chunk, name).copy() for name in _CHUNK_ARRAY_FIELDS})
        for chunk in buf.iter_chunks_consuming()
    ]


def _assert_chunk_payloads_equal(left, right):
    assert len(left) == len(right)
    for (left_size, left_arrays), (right_size, right_arrays) in zip(left, right):
        assert left_size == right_size
        for name in _CHUNK_ARRAY_FIELDS:
            assert np.array_equal(left_arrays[name], right_arrays[name]), name


# =====================================================================
# FragmentAccumulator -- direct C++ tests
# =====================================================================


class TestFragmentAccumulator:
    """Test FragmentAccumulator directly (no FragmentBuffer)."""

    def test_append_and_size(self, mini_index):
        from rigel._resolve_impl import FragmentAccumulator

        acc = FragmentAccumulator()
        assert acc.size == 0

        result = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert result is not None
        acc.append(result, 0)
        assert acc.size == 1

    def test_finalize_returns_dict(self, mini_index):
        from rigel._resolve_impl import FragmentAccumulator

        acc = FragmentAccumulator()
        result = _resolve(mini_index, [_exon("chr1", 120, 180)])
        acc.append(result, 0)

        raw = acc.finalize()
        assert isinstance(raw, dict)
        assert raw["size"] == 1
        assert "splice_type" in raw
        assert "t_offsets" in raw
        assert "t_indices" in raw

    def test_finalize_carries_the_stored_ambig_strand(self, mini_index):
        from rigel._resolve_impl import FragmentAccumulator

        acc = FragmentAccumulator()
        result = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert result.ambig_strand == 0
        assert len(result.t_inds) >= 2

        acc.append(result, 0)
        raw = acc.finalize()

        assert raw["ambig_strand"].tolist() == [0]

    def test_finalize_multiple(self, mini_index):
        from rigel._resolve_impl import FragmentAccumulator

        acc = FragmentAccumulator()

        # Exon in g1 region -> hits t1, t2
        r1 = _resolve(mini_index, [_exon("chr1", 120, 180)])
        # Exon in g2 region -> hits t3
        r2 = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])
        acc.append(r1, 0)
        acc.append(r2, 1)
        assert acc.size == 2

        raw = acc.finalize()
        assert raw["size"] == 2

    def test_finalize_buffer_payload_dtypes(self, mini_index):
        from rigel._resolve_impl import FragmentAccumulator

        acc = FragmentAccumulator()
        result = _resolve(mini_index, [_exon("chr1", 120, 180)])
        acc.append(result, 0)

        raw = acc.finalize()
        assert raw["frag_lengths"].dtype == np.int32
        assert raw["exon_bp"].dtype == np.uint16
        assert raw["read_length"].dtype == np.uint16
        assert "intron_bp" not in raw

    def test_append_read_length_overflow_names_column(self, mini_index):
        from rigel._resolve_impl import FragmentAccumulator

        result = _resolve(mini_index, [_exon("chr1", 10_000, 10_000 + 65_536)])
        assert result is not None
        assert result.read_length == 65_536

        acc = FragmentAccumulator()
        with pytest.raises(RuntimeError, match="read_length"):
            acc.append(result, 0)

    def test_append_read_length_uint16_boundary(self, mini_index):
        from rigel._resolve_impl import FragmentAccumulator

        result = _resolve(mini_index, [_exon("chr1", 10_000, 10_000 + 65_535)])
        assert result is not None
        assert result.read_length == 65_535

        acc = FragmentAccumulator()
        acc.append(result, 0)
        raw = acc.finalize()
        assert raw["read_length"].tolist() == [65_535]


# =====================================================================
# FragmentBuffer -- the chunks it takes and hands back
# =====================================================================


class TestFragmentBufferBasic:
    def test_empty_buffer(self, mini_index):
        buf = FragmentBuffer()
        assert buf.total_fragments == 0
        assert buf.n_chunks == 0
        assert _consume(buf) == []

    def test_single_fragment(self, mini_index):
        buf = FragmentBuffer()
        result = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert result is not None
        _inject(buf, [result])

        assert buf.total_fragments == 1
        assert buf.n_chunks == 1
        chunk = _only_chunk(buf)
        assert chunk.size == 1
        # t1 and t2 both have exon (99,200), so this should hit both
        assert int(np.diff(chunk.t_offsets)[0]) >= 1

    def test_roundtrip_preserves_splice_type(self, mini_index):
        """Splice type should survive through the C++ accumulator."""
        r_unspliced = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r_unspliced is not None

        buf = FragmentBuffer()
        _inject(buf, [r_unspliced])

        assert _only_chunk(buf).splice_type[0] == int(SpliceType.UNSPLICED)

    def test_roundtrip_preserves_strand(self, mini_index):
        """Exon strand should survive through the C++ accumulator."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180, Strand.POS)])
        buf = FragmentBuffer()
        _inject(buf, [r])

        assert _only_chunk(buf).align_strand[0] == int(Strand.POS)

    def test_multiple_fragments_different_genes(self, mini_index):
        """Fragments from different genes should be buffered correctly."""
        r1 = _resolve(mini_index, [_exon("chr1", 120, 180)])
        r2 = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])

        buf = FragmentBuffer()
        _inject(buf, [r1, r2])

        assert _only_chunk(buf).size == 2

    def test_chunks_are_consumed_in_order(self, mini_index):
        buf = FragmentBuffer()
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 25, chunk_size=10)

        assert buf.total_fragments == 25
        assert buf.n_chunks == 3  # 10, 10, 5

        chunks = _consume(buf)
        assert [chunk.size for chunk in chunks] == [10, 10, 5]
        assert np.concatenate([chunk.frag_id for chunk in chunks]).tolist() == list(range(25))
        assert buf.n_chunks == 0
        assert buf.memory_bytes == 0

    def test_num_hits_preserved(self, mini_index):
        """num_hits set on ResolvedFragment should survive buffer round-trip."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        r.num_hits = 3

        buf = FragmentBuffer()
        _inject(buf, [r])

        assert _only_chunk(buf).num_hits[0] == 3

    def test_nm_preserved(self, mini_index):
        """NM edit distance should survive buffer round-trip."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        r.nm = 5

        buf = FragmentBuffer()
        _inject(buf, [r])

        assert _only_chunk(buf).nm[0] == 5

    def test_nm_default_zero(self, mini_index):
        """Default NM should be 0."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])

        buf = FragmentBuffer()
        _inject(buf, [r])

        assert _only_chunk(buf).nm[0] == 0

    def test_chunk_payload_dtypes(self, mini_index):
        """Buffer chunks store bounded hot-path payloads as uint16."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])

        buf = FragmentBuffer()
        _inject(buf, [r])

        chunk = _only_chunk(buf)
        assert chunk.frag_lengths.dtype == np.int32
        assert chunk.exon_bp.dtype == np.uint16
        assert chunk.read_length.dtype == np.uint16
        assert not hasattr(chunk, "intron_bp")

    def test_memory_bytes_positive(self, mini_index):
        buf = FragmentBuffer()
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 50)

        assert buf.memory_bytes > 0
        assert _only_chunk(buf).memory_bytes > 0


# =====================================================================
# FragmentBuffer -- frag_id round-trip
# =====================================================================


class TestFragId:
    def test_frag_id_preserves_value(self, mini_index):
        """Explicit frag_id should survive append -> finalize -> inject -> consume."""
        buf = FragmentBuffer()
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 3, frag_ids=[42, 42, 99])

        assert _only_chunk(buf).frag_id.tolist() == [42, 42, 99]

    def test_frag_id_chunk_array(self, mini_index):
        """frag_id should be accessible as chunk array."""
        buf = FragmentBuffer()
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 5, frag_ids=[i // 2 for i in range(5)])

        assert list(_only_chunk(buf).frag_id) == [0, 0, 1, 1, 2]

    def test_frag_id_survives_spill(self, mini_index, tmp_path):
        """frag_id should survive Arrow IPC spill and reload."""
        buf = FragmentBuffer(
            max_memory_bytes=1,  # force spill
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        frag_ids = [i // 3 for i in range(100)]
        _inject(buf, [r] * 100, chunk_size=50, frag_ids=frag_ids)

        assert buf.n_spilled > 0
        result_ids = np.concatenate([chunk.frag_id for chunk in _consume(buf)]).tolist()
        assert result_ids == frag_ids


# =====================================================================
# FragmentBuffer -- fragment_classes
# =====================================================================


class TestFragmentClasses:
    def test_unique_gene_single_transcript(self, mini_index):
        """g2 exon region -> t3 + synthetic nRNA -> FRAG_AMBIG_SAME_STRAND."""
        r = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])
        assert r is not None
        assert r.ambig_strand == 0
        t_inds = list(r.t_inds)
        assert len(t_inds) == 2  # t3 + synthetic nRNA

        buf = FragmentBuffer()
        _inject(buf, [r])

        assert _only_chunk(buf).fragment_classes[0] == FRAG_AMBIG_SAME_STRAND

    def test_isoform_ambiguous(self, mini_index):
        """g1 shared exon region -> t1 + t2 (same strand) -> FRAG_AMBIG_SAME_STRAND."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r is not None
        assert r.ambig_strand == 0
        assert len(list(r.t_inds)) == 3  # t1, t2 + synthetic nRNA

        buf = FragmentBuffer()
        _inject(buf, [r])

        assert _only_chunk(buf).fragment_classes[0] == FRAG_AMBIG_SAME_STRAND

    def test_multimapper(self, mini_index):
        """NH > 1 -> FRAG_MULTIMAPPER regardless of gene count."""
        r = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])
        r.num_hits = 3

        buf = FragmentBuffer()
        _inject(buf, [r])

        assert _only_chunk(buf).fragment_classes[0] == FRAG_MULTIMAPPER

    def test_mixed_classes(self, mini_index):
        """Multiple fragment classes in one chunk."""
        r_unambig = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])
        r_iso = _resolve(mini_index, [_exon("chr1", 120, 180)])
        r_mm = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])
        r_mm.num_hits = 2

        buf = FragmentBuffer()
        _inject(buf, [r_unambig, r_iso, r_mm])

        fc = _only_chunk(buf).fragment_classes
        assert fc[0] == FRAG_AMBIG_SAME_STRAND  # t3 + synthetic nRNA
        assert fc[1] == FRAG_AMBIG_SAME_STRAND
        assert fc[2] == FRAG_MULTIMAPPER


# =====================================================================
# FragmentBuffer -- disk spill
# =====================================================================


class TestDiskSpill:
    def test_spill_triggers_above_threshold(self, mini_index, tmp_path):
        """When in-memory chunks exceed max_memory_bytes, spill to disk."""
        buf = FragmentBuffer(
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 150, chunk_size=50)

        assert buf.n_spilled > 0
        assert buf.total_fragments == 150
        buf.cleanup()

    def test_spill_preserves_data(self, mini_index, tmp_path):
        """Data roundtrips correctly through Arrow IPC spill."""
        buf = FragmentBuffer(
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 100, chunk_size=50)

        assert buf.n_spilled > 0
        chunks = _consume(buf)
        assert sum(chunk.size for chunk in chunks) == 100
        for chunk in chunks:
            assert (chunk.splice_type == int(SpliceType.UNSPLICED)).all()
        buf.cleanup()

    def test_forced_spill_matches_in_memory_chunks(self, mini_index, tmp_path):
        """Forced-spill chunks match the in-memory path exactly."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])

        in_memory = FragmentBuffer(
            max_memory_bytes=0,
        )
        spilled = FragmentBuffer(
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        _inject(in_memory, [r] * 75, chunk_size=25)
        _inject(spilled, [r] * 75, chunk_size=25)

        try:
            assert spilled.n_spilled > 0
            _assert_chunk_payloads_equal(_chunk_payloads(in_memory), _chunk_payloads(spilled))
        finally:
            in_memory.cleanup()
            spilled.cleanup()

    def test_iter_chunks_consuming_waits_and_deletes_spill(self, mini_index, tmp_path):
        buf = FragmentBuffer(
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 40, chunk_size=20)

        spill_paths = [c.path for c in buf._chunks if isinstance(c, buffer_mod._PendingSpill)]
        assert spill_paths
        chunks = _consume(buf)

        assert sum(chunk.size for chunk in chunks) == 40
        assert len(buf._chunks) == 0
        assert all(not path.exists() for path in spill_paths)
        buf.cleanup()

    def test_writer_exception_reaches_next_consumer(self, mini_index, tmp_path, monkeypatch):
        def fail_spill(chunk, path):
            raise OSError("synthetic spill failure")

        monkeypatch.setattr(buffer_mod, "_spill_chunk", fail_spill)
        buf = FragmentBuffer(
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 10, chunk_size=10)

        with pytest.raises(RuntimeError, match="Failed to spill buffer chunk"):
            buf.cleanup()
        with pytest.raises(RuntimeError, match="Failed to spill buffer chunk"):
            _consume(buf)
        assert _spill_dirs(tmp_path) == []

    def test_cleanup_waits_for_pending_spill(self, mini_index, tmp_path, monkeypatch):
        real_spill = buffer_mod._spill_chunk
        started = threading.Event()
        allow_write = threading.Event()

        def slow_spill(chunk, path):
            started.set()
            assert allow_write.wait(timeout=5)
            real_spill(chunk, path)

        monkeypatch.setattr(buffer_mod, "_spill_chunk", slow_spill)
        buf = FragmentBuffer(
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 10, chunk_size=10)

        assert started.wait(timeout=5)
        (pending,) = [c for c in buf._chunks if isinstance(c, buffer_mod._PendingSpill)]
        assert not pending.done.is_set()

        cleanup_done = threading.Event()
        cleanup_errors = []

        def run_cleanup():
            try:
                buf.cleanup()
            except BaseException as exc:
                cleanup_errors.append(exc)
            finally:
                cleanup_done.set()

        cleanup_thread = threading.Thread(target=run_cleanup)
        cleanup_thread.start()
        assert not cleanup_done.wait(timeout=0.05)
        allow_write.set()
        cleanup_thread.join(timeout=5)

        assert cleanup_done.is_set()
        assert cleanup_errors == []
        assert _spill_dirs(tmp_path) == []

    def test_cleanup_idempotent_with_spills(self, mini_index, tmp_path):
        buf = FragmentBuffer(
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 40, chunk_size=20)

        buf.cleanup()
        buf.cleanup()

        assert _spill_dirs(tmp_path) == []

    def test_cleanup_removes_files(self, mini_index, tmp_path):
        buf = FragmentBuffer(
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 100, chunk_size=50)

        assert buf.n_spilled > 0
        assert len(_spill_dirs(tmp_path)) > 0

        buf.cleanup()
        assert _spill_dirs(tmp_path) == []

    def test_context_manager_cleanup(self, mini_index, tmp_path):
        with FragmentBuffer(
            max_memory_bytes=1,
            spill_dir=tmp_path,
        ) as buf:
            r = _resolve(mini_index, [_exon("chr1", 120, 180)])
            _inject(buf, [r] * 100, chunk_size=50)
            assert buf.n_spilled > 0

        assert _spill_dirs(tmp_path) == []

    def test_no_spill_when_disabled(self, mini_index):
        """max_memory_bytes=0 disables spilling."""
        buf = FragmentBuffer(
            max_memory_bytes=0,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 25, chunk_size=10)

        assert buf.n_spilled == 0
        assert buf.n_chunks == 3

    def test_no_spill_under_threshold(self, mini_index):
        """Small buffer should not spill."""
        buf = FragmentBuffer(
            max_memory_bytes=100 * 1024**2,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 25)

        assert buf.n_spilled == 0


# =====================================================================
# FragmentBuffer -- the finalizer, for a buffer dropped without cleanup()
# =====================================================================


class TestFinalizer:
    def test_dropped_buffer_is_collected_and_its_spill_dir_removed(self, mini_index, tmp_path):
        """The finalizer must not hold the buffer: one that does keeps it, its chunks and its writer
        thread alive until the interpreter exits."""
        buf = FragmentBuffer(max_memory_bytes=1, spill_dir=tmp_path)
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 20, chunk_size=10)
        assert buf.n_spilled > 0
        assert _spill_dirs(tmp_path)

        alive = weakref.ref(buf)
        del buf
        gc.collect()

        assert alive() is None
        assert _spill_dirs(tmp_path) == []

    def test_finalizer_on_the_writer_thread_does_not_wait_on_itself(
        self, mini_index, tmp_path, monkeypatch
    ):
        """Garbage collection can run the finalizer on the spill writer's own thread. There it must
        neither join nor wait on that thread; the thread still writes what it holds, removes the
        directory and exits."""
        real_spill = buffer_mod._spill_chunk
        ready = threading.Event()
        finalized = threading.Event()
        errors = []
        box = {}

        def spill_then_finalize(chunk, path):
            assert ready.wait(timeout=5)
            real_spill(chunk, path)
            try:
                box["finalizer"]()  # what garbage collection would call, on this thread
            except BaseException as exc:
                errors.append(exc)
            finally:
                finalized.set()

        monkeypatch.setattr(buffer_mod, "_spill_chunk", spill_then_finalize)
        buf = FragmentBuffer(max_memory_bytes=1, spill_dir=tmp_path)
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        _inject(buf, [r] * 10, chunk_size=10)
        box["finalizer"] = buf._spill_finalizer
        thread = buf._spill_writer._thread
        ready.set()

        assert finalized.wait(timeout=5)
        assert errors == []
        thread.join(timeout=5)
        assert not thread.is_alive()
        assert _spill_dirs(tmp_path) == []


# =====================================================================
# ResolvedFragment -- C++ object properties
# =====================================================================


class TestResolvedFragment:
    """Test the C++ ResolvedFragment object returned by resolve_fragment."""

    def test_t_inds_type(self, mini_index):
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r is not None
        t_list = list(r.t_inds)
        assert len(t_list) >= 1

    def test_unspliced_fragment(self, mini_index):
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r.splice_type == int(SpliceType.UNSPLICED)

    def test_spliced_fragment(self, mini_index):
        """Fragment with annotated intron -> SPLICED_ANNOT."""
        r = _resolve(
            mini_index,
            exons=[
                _exon("chr1", 120, 200),
                _exon("chr1", 299, 380),
            ],
            introns=[
                GenomicInterval("chr1", 200, 299, Strand.POS),
            ],
        )
        assert r is not None
        assert r.splice_type == int(SpliceType.SPLICED_ANNOT)

    def test_intergenic_resolves_with_no_candidates(self, mini_index):
        """Fragment in intergenic region -> resolves with empty t_inds.

        Zero-candidate fragments are not dropped at the resolver boundary; they
        flow through with empty ``t_inds`` so calibration can deposit them as
        INTERGENIC.
        """
        r = _resolve(mini_index, [_exon("chr1", 1500, 1600)])
        assert r is not None
        assert len(r.t_inds) == 0

    def test_unique_gene_property(self, mini_index):
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r.is_same_strand is True

    def test_num_hits_default(self, mini_index):
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r.num_hits == 1

    def test_num_hits_mutable(self, mini_index):
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        r.num_hits = 5
        assert r.num_hits == 5

    def test_genomic_footprint(self, mini_index):
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r.genomic_footprint == 60  # 180 - 120

    def test_genomic_start(self, mini_index):
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r.genomic_start == 120
