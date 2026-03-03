"""Tests for hulkrna.buffer -- FragmentBuffer with native C++ accumulator.

All buffer append tests use C++ ResolvedFragment objects produced by
FragmentResolver.resolve_fragment(), exercising the real native code path.
"""

import numpy as np
import pytest

from hulkrna.types import Strand, MergeOutcome, GenomicInterval, ChimeraType
from hulkrna.splice import SpliceType
from hulkrna.resolution import make_fragment, resolve_fragment
from hulkrna.buffer import (
    FragmentBuffer,
    BufferedFragment,
    FRAG_UNIQUE,
    FRAG_AMBIG_SAME_STRAND,
    FRAG_AMBIG_OPP_STRAND,
    FRAG_MULTIMAPPER,
    FRAG_CHIMERIC,
)
from hulkrna.index import TranscriptIndex


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


# =====================================================================
# BufferedFragment duck-typing properties (no buffer needed)
# =====================================================================


class TestBufferedFragment:
    def test_is_same_strand_single(self):
        bf = BufferedFragment(
            t_inds=np.array([0, 1], dtype=np.int32),
            ambig_strand=0,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250, 250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeOutcome.INTERSECTION),
        )
        assert bf.is_same_strand is True

    def test_is_ambiguous_multi_gene(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            ambig_strand=1,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeOutcome.INTERSECTION),
        )
        assert bf.ambig_strand > 0

    def test_is_ambiguous_multimapped(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            ambig_strand=0,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=3,
            merge_criteria=int(MergeOutcome.INTERSECTION),
        )
        assert bf.num_hits > 1

    def test_has_annotated_sj(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            ambig_strand=0,
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NEG),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeOutcome.INTERSECTION),
        )
        assert bf.splice_type == int(SpliceType.SPLICED_ANNOT)

    def test_is_strand_qualified(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            ambig_strand=0,
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NEG),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeOutcome.INTERSECTION),
        )
        assert bf.is_strand_qualified is True

    def test_not_strand_qualified_wrong_cat(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            ambig_strand=0,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NEG),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeOutcome.INTERSECTION),
        )
        assert bf.is_strand_qualified is False

    def test_t_inds_iterable_and_len(self):
        bf = BufferedFragment(
            t_inds=np.array([10, 20, 30], dtype=np.int32),
            ambig_strand=0,
            splice_type=0,
            exon_strand=1,
            sj_strand=0,
            frag_lengths=np.array([100, 100, 100], dtype=np.int32),
            num_hits=1,
            merge_criteria=0,
        )
        assert len(bf.t_inds) == 3
        assert list(bf.t_inds) == [10, 20, 30]


# =====================================================================
# FragmentAccumulator -- direct C++ tests
# =====================================================================


class TestFragmentAccumulator:
    """Test FragmentAccumulator directly (no FragmentBuffer)."""

    def test_append_and_size(self, mini_index):
        from hulkrna._resolve_impl import FragmentAccumulator

        acc = FragmentAccumulator()
        assert acc.size == 0

        result = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert result is not None
        acc.append(result, 0)
        assert acc.size == 1

    def test_finalize_returns_dict(self, mini_index):
        from hulkrna._resolve_impl import FragmentAccumulator

        acc = FragmentAccumulator()
        result = _resolve(mini_index, [_exon("chr1", 120, 180)])
        acc.append(result, 0)

        raw = acc.finalize(mini_index.t_to_strand_arr.tolist())
        assert isinstance(raw, dict)
        assert raw["size"] == 1
        assert "splice_type" in raw
        assert "t_offsets" in raw
        assert "t_indices" in raw

    def test_finalize_multiple(self, mini_index):
        from hulkrna._resolve_impl import FragmentAccumulator

        acc = FragmentAccumulator()

        # Exon in g1 region -> hits t1, t2
        r1 = _resolve(mini_index, [_exon("chr1", 120, 180)])
        # Exon in g2 region -> hits t3
        r2 = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])
        acc.append(r1, 0)
        acc.append(r2, 1)
        assert acc.size == 2

        raw = acc.finalize(mini_index.t_to_strand_arr.tolist())
        assert raw["size"] == 2


# =====================================================================
# FragmentBuffer -- basic accumulation via native path
# =====================================================================


class TestFragmentBufferBasic:
    def test_empty_buffer(self, mini_index):
        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.finalize()
        assert buf.total_fragments == 0
        assert buf.n_chunks == 0
        assert list(buf) == []

    def test_single_fragment(self, mini_index):
        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        result = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert result is not None
        buf.append(result)
        buf.finalize()

        assert buf.total_fragments == 1
        assert buf.n_chunks == 1
        frags = list(buf)
        assert len(frags) == 1
        # t1 and t2 both have exon (99,200), so this should hit both
        assert len(frags[0].t_inds) >= 1

    def test_roundtrip_preserves_splice_type(self, mini_index):
        """Splice type should survive through the C++ accumulator."""
        r_unspliced = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r_unspliced is not None

        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r_unspliced)
        buf.finalize()

        bf = list(buf)[0]
        assert bf.splice_type == int(SpliceType.UNSPLICED)

    def test_roundtrip_preserves_strand(self, mini_index):
        """Exon strand should survive through the C++ accumulator."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180, Strand.POS)])
        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r)
        buf.finalize()

        bf = list(buf)[0]
        assert bf.exon_strand == int(Strand.POS)

    def test_multiple_fragments_different_genes(self, mini_index):
        """Fragments from different genes should be buffered correctly."""
        r1 = _resolve(mini_index, [_exon("chr1", 120, 180)])
        r2 = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])

        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r1, frag_id=0)
        buf.append(r2, frag_id=1)
        buf.finalize()

        frags = list(buf)
        assert len(frags) == 2

    def test_many_fragments_chunking(self, mini_index):
        """Buffer should create multiple chunks when chunk_size is exceeded."""
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr, chunk_size=10
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(25):
            buf.append(r, frag_id=i)
        buf.finalize()

        assert buf.total_fragments == 25
        assert buf.n_chunks == 3  # 10, 10, 5

        chunks = list(buf.iter_chunks())
        assert len(chunks) == 3
        assert chunks[0].size == 10
        assert chunks[1].size == 10
        assert chunks[2].size == 5

    def test_num_hits_preserved(self, mini_index):
        """num_hits set on ResolvedFragment should survive buffer round-trip."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        r.num_hits = 3

        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r)
        buf.finalize()

        bf = list(buf)[0]
        assert bf.num_hits == 3

    def test_nm_preserved(self, mini_index):
        """NM edit distance should survive buffer round-trip."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        r.nm = 5

        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r)
        buf.finalize()

        bf = list(buf)[0]
        assert bf.nm == 5

    def test_nm_default_zero(self, mini_index):
        """Default NM should be 0."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])

        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r)
        buf.finalize()

        bf = list(buf)[0]
        assert bf.nm == 0


# =====================================================================
# FragmentBuffer -- frag_id round-trip
# =====================================================================


class TestFragId:
    def test_frag_id_default_zero(self, mini_index):
        """Default frag_id should be 0."""
        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        buf.append(r)
        buf.finalize()

        bf = list(buf)[0]
        assert bf.frag_id == 0

    def test_frag_id_preserves_value(self, mini_index):
        """Explicit frag_id should survive append -> finalize -> iterate."""
        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        buf.append(r, frag_id=42)
        buf.append(r, frag_id=42)
        buf.append(r, frag_id=99)
        buf.finalize()

        frags = list(buf)
        assert frags[0].frag_id == 42
        assert frags[1].frag_id == 42
        assert frags[2].frag_id == 99

    def test_frag_id_chunk_array(self, mini_index):
        """frag_id should be accessible as chunk array."""
        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(5):
            buf.append(r, frag_id=i // 2)
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert list(chunk.frag_id) == [0, 0, 1, 1, 2]

    def test_frag_id_survives_spill(self, mini_index, tmp_path):
        """frag_id should survive Arrow IPC spill and reload."""
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr,
            chunk_size=50,
            max_memory_bytes=1,  # force spill
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        frag_ids = []
        for i in range(100):
            fid = i // 3
            buf.append(r, frag_id=fid)
            frag_ids.append(fid)
        buf.finalize()

        assert buf.n_spilled > 0
        result_ids = [bf.frag_id for bf in buf]
        assert result_ids == frag_ids


# =====================================================================
# FragmentBuffer -- fragment_classes
# =====================================================================


class TestFragmentClasses:
    def test_unique_gene_single_transcript(self, mini_index):
        """Single-transcript hit, NH=1 -> FRAG_UNIQUE."""
        r = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])
        assert r is not None
        assert r.ambig_strand == 0
        t_inds = list(r.t_inds)
        assert len(t_inds) == 1

        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r)
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.fragment_classes[0] == FRAG_UNIQUE

    def test_isoform_ambiguous(self, mini_index):
        """g1 shared exon region -> t1 + t2 (same strand) -> FRAG_AMBIG_SAME_STRAND."""
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        assert r is not None
        assert r.ambig_strand == 0
        assert len(list(r.t_inds)) == 2

        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r)
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.fragment_classes[0] == FRAG_AMBIG_SAME_STRAND

    def test_multimapper(self, mini_index):
        """NH > 1 -> FRAG_MULTIMAPPER regardless of gene count."""
        r = _resolve(mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)])
        r.num_hits = 3

        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r)
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.fragment_classes[0] == FRAG_MULTIMAPPER

    def test_mixed_classes(self, mini_index):
        """Multiple fragment classes in one chunk."""
        r_unique = _resolve(
            mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)]
        )
        r_iso = _resolve(mini_index, [_exon("chr1", 120, 180)])
        r_mm = _resolve(
            mini_index, [_exon("chr1", 1020, 1080, Strand.NEG)]
        )
        r_mm.num_hits = 2

        buf = FragmentBuffer(t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100)
        buf.append(r_unique)
        buf.append(r_iso)
        buf.append(r_mm)
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        fc = chunk.fragment_classes
        assert fc[0] == FRAG_UNIQUE
        assert fc[1] == FRAG_AMBIG_SAME_STRAND
        assert fc[2] == FRAG_MULTIMAPPER


# =====================================================================
# FragmentBuffer -- disk spill
# =====================================================================


class TestDiskSpill:
    def test_spill_triggers_above_threshold(self, mini_index, tmp_path):
        """When in-memory chunks exceed max_memory_bytes, spill to disk."""
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr,
            chunk_size=50,
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(150):
            buf.append(r, frag_id=i)
        buf.finalize()

        assert buf.n_spilled > 0
        assert buf.total_fragments == 150

    def test_spill_preserves_data(self, mini_index, tmp_path):
        """Data roundtrips correctly through Arrow IPC spill."""
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr,
            chunk_size=50,
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(100):
            buf.append(r, frag_id=i)
        buf.finalize()

        assert buf.n_spilled > 0
        result = list(buf)
        assert len(result) == 100
        for bf in result:
            assert bf.splice_type == int(SpliceType.UNSPLICED)

    def test_cleanup_removes_files(self, mini_index, tmp_path):
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr,
            chunk_size=50,
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(100):
            buf.append(r, frag_id=i)
        buf.finalize()

        assert buf.n_spilled > 0
        spill_dirs = list(tmp_path.iterdir())
        assert len(spill_dirs) > 0

        buf.cleanup()
        remaining = [p for p in tmp_path.iterdir() if p.is_dir()]
        assert len(remaining) == 0

    def test_context_manager_cleanup(self, mini_index, tmp_path):
        with FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr,
            chunk_size=50,
            max_memory_bytes=1,
            spill_dir=tmp_path,
        ) as buf:
            r = _resolve(mini_index, [_exon("chr1", 120, 180)])
            for i in range(100):
                buf.append(r, frag_id=i)
            buf.finalize()
            assert buf.n_spilled > 0

        remaining = [p for p in tmp_path.iterdir() if p.is_dir()]
        assert len(remaining) == 0

    def test_no_spill_when_disabled(self, mini_index):
        """max_memory_bytes=0 disables spilling."""
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr,
            chunk_size=10,
            max_memory_bytes=0,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(25):
            buf.append(r, frag_id=i)
        buf.finalize()

        assert buf.n_spilled == 0
        assert buf.n_chunks == 3

    def test_no_spill_under_threshold(self, mini_index):
        """Small buffer should not spill."""
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr,
            chunk_size=100,
            max_memory_bytes=100 * 1024**2,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(25):
            buf.append(r, frag_id=i)
        buf.finalize()

        assert buf.n_spilled == 0


# =====================================================================
# FragmentBuffer -- summary
# =====================================================================


class TestBufferSummary:
    def test_summary_structure(self, mini_index):
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(150):
            buf.append(r, frag_id=i)
        buf.finalize()

        s = buf.summary()
        assert s["total_fragments"] == 150
        assert s["n_chunks"] == 2
        assert s["in_memory_chunks"] == 2
        assert s["on_disk_chunks"] == 0
        assert s["memory_bytes"] > 0
        assert isinstance(s["memory_mb"], float)

    def test_summary_with_spill(self, mini_index, tmp_path):
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr,
            chunk_size=50,
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(100):
            buf.append(r, frag_id=i)
        buf.finalize()

        s = buf.summary()
        assert s["total_fragments"] == 100
        assert s["on_disk_chunks"] > 0


# =====================================================================
# FragmentBuffer -- iter_chunks
# =====================================================================


class TestIterChunks:
    def test_iter_chunks_yields_correct_count(self, mini_index):
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr, chunk_size=10
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(25):
            buf.append(r, frag_id=i)
        buf.finalize()

        chunks = list(buf.iter_chunks())
        assert len(chunks) == 3
        total = sum(c.size for c in chunks)
        assert total == 25

    def test_iter_chunks_getitem(self, mini_index):
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(5):
            buf.append(r, frag_id=i)
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.size == 5
        for i in range(5):
            bf = chunk[i]
            assert len(bf.t_inds) >= 1
            assert bf.frag_id == i

    def test_memory_bytes_positive(self, mini_index):
        buf = FragmentBuffer(
            t_strand_arr=mini_index.t_to_strand_arr, chunk_size=100
        )
        r = _resolve(mini_index, [_exon("chr1", 120, 180)])
        for i in range(50):
            buf.append(r, frag_id=i)
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.memory_bytes > 0


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

    def test_intergenic_returns_none(self, mini_index):
        """Fragment in intergenic region -> None."""
        r = _resolve(mini_index, [_exon("chr1", 1500, 1600)])
        assert r is None

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
