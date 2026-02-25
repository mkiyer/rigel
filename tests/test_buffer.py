"""Tests for hulkrna.buffer — FragmentBuffer with columnar storage and disk spill."""

import numpy as np
import pytest

from hulkrna.types import Strand, MergeCriteria
from hulkrna.categories import SpliceType
from hulkrna.resolution import ResolvedFragment
from hulkrna.buffer import FragmentBuffer, BufferedFragment, FRAG_UNIQUE, FRAG_ISOFORM_AMBIG, FRAG_GENE_AMBIG, FRAG_MULTIMAPPER, FRAG_CHIMERIC
from hulkrna.types import ChimeraType


# =====================================================================
# Helpers
# =====================================================================


def _make_resolved(**kwargs):
    """Create a ResolvedFragment with sensible defaults."""
    defaults = dict(
        t_inds=frozenset({0, 1}),
        n_genes=1,
        splice_type=SpliceType.UNSPLICED,
        exon_strand=Strand.POS,
        sj_strand=Strand.NONE,
        frag_lengths={0: 250, 1: 250},
        merge_criteria=MergeCriteria.INTERSECTION,
        num_hits=1,
        genomic_footprint=250,
    )
    defaults.update(kwargs)
    return ResolvedFragment(**defaults)


def _make_t_to_g_arr(n=5):
    """Create a simple t_to_g_arr for testing. t0-t4 → g0-g2 (round-robin)."""
    import numpy as np
    return np.array([i % 3 for i in range(n)], dtype=np.int64)


def _fill_buffer(n, *, chunk_size=100, **buffer_kwargs):
    """Create and fill a buffer with n fragments, returning (buffer, fragments)."""
    t_to_g = _make_t_to_g_arr()
    buf = FragmentBuffer(t_to_g_arr=t_to_g, chunk_size=chunk_size, **buffer_kwargs)
    fragments = []
    for i in range(n):
        t_set = frozenset({i % 5, (i + 1) % 5})
        # Derive n_genes from t_set using t_to_g
        n_genes = len(set(int(t_to_g[t]) for t in t_set))
        flen = 200 + i
        r = _make_resolved(
            t_inds=t_set,
            n_genes=n_genes,
            splice_type=SpliceType(i % len(SpliceType)),
            exon_strand=Strand(1 + i % 2),  # POS or NEG
            sj_strand=Strand(1 + (i + 1) % 2),
            frag_lengths={t: flen for t in t_set},
            num_hits=1 + i % 2,
        )
        fragments.append(r)
        buf.append(r)
    buf.finalize()
    return buf, fragments


# =====================================================================
# BufferedFragment duck-typing properties
# =====================================================================


class TestBufferedFragment:
    def test_is_unique_gene_single(self):
        bf = BufferedFragment(
            t_inds=np.array([0, 1], dtype=np.int32),
            n_genes=1,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250, 250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert bf.is_unique_gene is True
        assert (bf.n_genes > 1 or bf.num_hits > 1) is False

    def test_is_ambiguous_multi_gene(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            n_genes=2,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert (bf.n_genes > 1 or bf.num_hits > 1) is True

    def test_is_ambiguous_multimapped(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            n_genes=1,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=3,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert (bf.n_genes > 1 or bf.num_hits > 1) is True

    def test_has_annotated_sj(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            n_genes=1,
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NEG),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert bf.splice_type == int(SpliceType.SPLICED_ANNOT)

    def test_is_strand_qualified(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            n_genes=1,
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NEG),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert bf.is_strand_qualified is True

    def test_not_strand_qualified_wrong_cat(self):
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            n_genes=1,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NEG),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert bf.is_strand_qualified is False

    def test_t_inds_iterable_and_len(self):
        bf = BufferedFragment(
            t_inds=np.array([10, 20, 30], dtype=np.int32),
            n_genes=1,
            splice_type=0,
            exon_strand=1,
            sj_strand=0,
            frag_lengths=np.array([100, 100, 100], dtype=np.int32),
            num_hits=1,
            merge_criteria=0,
        )
        assert len(bf.t_inds) == 3
        assert list(bf.t_inds) == [10, 20, 30]

    def test_is_isoform_ambiguous_true(self):
        """1 gene, 2 transcripts, NH=1 → isoform-ambiguous."""
        bf = BufferedFragment(
            t_inds=np.array([0, 1], dtype=np.int32),
            n_genes=1,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250, 250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert (bf.n_genes == 1 and len(bf.t_inds) > 1 and bf.num_hits == 1) is True

    def test_is_isoform_ambiguous_false_unique(self):
        """1 gene, 1 transcript, NH=1 → truly unique, not isoform-ambiguous."""
        bf = BufferedFragment(
            t_inds=np.array([0], dtype=np.int32),
            n_genes=1,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert (bf.n_genes == 1 and len(bf.t_inds) > 1 and bf.num_hits == 1) is False

    def test_is_isoform_ambiguous_false_multi_gene(self):
        """2 genes → gene-ambiguous, not isoform-ambiguous."""
        bf = BufferedFragment(
            t_inds=np.array([0, 1], dtype=np.int32),
            n_genes=2,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250, 250], dtype=np.int32),
            num_hits=1,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert (bf.n_genes == 1 and len(bf.t_inds) > 1 and bf.num_hits == 1) is False

    def test_is_isoform_ambiguous_false_multimapped(self):
        """1 gene, 2 transcripts, NH=3 → multimapper, not isoform-ambiguous."""
        bf = BufferedFragment(
            t_inds=np.array([0, 1], dtype=np.int32),
            n_genes=1,
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            frag_lengths=np.array([250, 250], dtype=np.int32),
            num_hits=3,
            merge_criteria=int(MergeCriteria.INTERSECTION),
        )
        assert (bf.n_genes == 1 and len(bf.t_inds) > 1 and bf.num_hits == 1) is False


# =====================================================================
# FragmentBuffer — frag_id round-trip
# =====================================================================


class TestFragId:
    def test_frag_id_default_zero(self):
        """Default frag_id should be 0."""
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.append(_make_resolved())
        buf.finalize()

        bf = list(buf)[0]
        assert bf.frag_id == 0

    def test_frag_id_preserves_value(self):
        """Explicit frag_id should survive append → finalize → iterate."""
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.append(_make_resolved(), frag_id=42)
        buf.append(_make_resolved(), frag_id=42)
        buf.append(_make_resolved(), frag_id=99)
        buf.finalize()

        frags = list(buf)
        assert frags[0].frag_id == 42
        assert frags[1].frag_id == 42
        assert frags[2].frag_id == 99

    def test_frag_id_chunk_array(self):
        """frag_id should be accessible as chunk array."""
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        for i in range(5):
            buf.append(_make_resolved(), frag_id=i // 2)
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert list(chunk.frag_id) == [0, 0, 1, 1, 2]

    def test_frag_id_survives_spill(self, tmp_path):
        """frag_id should survive Arrow IPC spill and reload."""
        buf = FragmentBuffer(
            t_to_g_arr=_make_t_to_g_arr(),
            chunk_size=50,
            max_memory_bytes=1,  # force spill
            spill_dir=tmp_path,
        )
        frag_ids = []
        for i in range(100):
            fid = i // 3  # groups of 3
            buf.append(_make_resolved(), frag_id=fid)
            frag_ids.append(fid)
        buf.finalize()

        assert buf.n_spilled > 0
        result_ids = [bf.frag_id for bf in buf]
        assert result_ids == frag_ids


# =====================================================================
# FragmentBuffer — basic accumulation
# =====================================================================


class TestFragmentBufferBasic:
    def test_empty_buffer(self):
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.finalize()
        assert buf.total_fragments == 0
        assert buf.n_chunks == 0
        assert list(buf) == []

    def test_single_fragment(self):
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        r = _make_resolved()
        buf.append(r)
        buf.finalize()

        assert buf.total_fragments == 1
        assert buf.n_chunks == 1
        frags = list(buf)
        assert len(frags) == 1
        assert frags[0].frag_lengths[0] == 250
        assert frags[0].splice_type == int(SpliceType.UNSPLICED)

    def test_roundtrip_preserves_data(self):
        """All fields should survive append → finalize → iterate."""
        r = _make_resolved(
            t_inds=frozenset({3, 7, 11}),
            n_genes=2,
            splice_type=SpliceType.SPLICED_ANNOT,
            exon_strand=Strand.NEG,
            sj_strand=Strand.POS,
            frag_lengths={3: 350, 7: 350, 11: 350},
            merge_criteria=MergeCriteria.UNION,
            num_hits=4,
        )
        # Need a t_to_g that maps indices 3,7,11 → two different genes
        import numpy as np
        t_to_g = np.zeros(12, dtype=np.int64)
        t_to_g[3] = 0
        t_to_g[7] = 0
        t_to_g[11] = 1
        buf = FragmentBuffer(t_to_g_arr=t_to_g, chunk_size=100)
        buf.append(r)
        buf.finalize()

        bf = list(buf)[0]
        assert sorted(bf.t_inds) == [3, 7, 11]
        assert bf.n_genes == 2
        assert bf.splice_type == int(SpliceType.SPLICED_ANNOT)
        assert bf.exon_strand == int(Strand.NEG)
        assert bf.sj_strand == int(Strand.POS)
        assert bf.frag_lengths[0] == 350
        assert bf.merge_criteria == int(MergeCriteria.UNION)
        assert bf.num_hits == 4

    def test_multiple_fragments(self):
        buf, fragments = _fill_buffer(50, chunk_size=100)
        assert buf.total_fragments == 50
        assert buf.n_chunks == 1  # 50 < chunk_size of 100

        result = list(buf)
        assert len(result) == 50
        for i, bf in enumerate(result):
            assert bf.frag_lengths[0] == 200 + i

    def test_chunking(self):
        buf, fragments = _fill_buffer(250, chunk_size=100)
        assert buf.total_fragments == 250
        assert buf.n_chunks == 3  # 100, 100, 50

        chunks = list(buf.iter_chunks())
        assert len(chunks) == 3
        assert chunks[0].size == 100
        assert chunks[1].size == 100
        assert chunks[2].size == 50

    def test_large_overlap_and_frag_length_preserved(self):
        """Large per-candidate overlap bp and fragment length do not overflow."""
        t_to_g = np.array([0, 0, 1], dtype=np.int64)
        buf = FragmentBuffer(t_to_g_arr=t_to_g, chunk_size=10)

        big_exon_bp = 120_000
        big_intron_bp = 80_000
        big_frag_len = 200_000

        r = _make_resolved(
            t_inds=frozenset({0, 2}),
            n_genes=2,
            overlap_bp={
                0: (big_exon_bp, big_intron_bp, big_intron_bp),
                2: (big_exon_bp - 1, big_intron_bp - 1, big_intron_bp - 1),
            },
            read_length=big_frag_len,
        )
        buf.append(r, frag_id=7)
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        bf = chunk[0]

        assert bf.read_length == big_frag_len
        assert int(bf.exon_bp[0]) == big_exon_bp
        assert int(bf.intron_bp[0]) == big_intron_bp
        assert int(bf.exon_bp[1]) == big_exon_bp - 1
        assert int(bf.intron_bp[1]) == big_intron_bp - 1
        assert int(bf.unambig_intron_bp[0]) == big_intron_bp
        assert int(bf.unambig_intron_bp[1]) == big_intron_bp - 1

    def test_nm_field_roundtrip(self):
        """NM edit distance is preserved through buffer append → read cycle."""
        t_to_g = np.array([0, 0, 1], dtype=np.int64)
        buf = FragmentBuffer(t_to_g_arr=t_to_g, chunk_size=10)

        r = _make_resolved(
            t_inds=frozenset({0, 2}), n_genes=2,
            overlap_bp={0: (100, 0, 0), 2: (100, 0, 0)},
            read_length=100, nm=7,
        )
        buf.append(r)
        buf.finalize()

        bf = list(buf)[0]
        assert bf.nm == 7

    def test_nm_default_zero(self):
        """When nm not set on ResolvedFragment, buffered value is 0."""
        t_to_g = np.array([0, 0], dtype=np.int64)
        buf = FragmentBuffer(t_to_g_arr=t_to_g, chunk_size=10)

        r = _make_resolved(
            t_inds=frozenset({0}), n_genes=1,
            overlap_bp={0: (50, 0, 0)}, read_length=50,
        )
        buf.append(r)
        buf.finalize()

        bf = list(buf)[0]
        assert bf.nm == 0


# =====================================================================
# FragmentBuffer — fragment_classes (tri-state classification)
# =====================================================================


class TestFragmentClasses:
    def test_truly_unique(self):
        """1 gene, 1 transcript, NH=1 → FRAG_UNIQUE."""
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.append(_make_resolved(
            t_inds=frozenset({0}), n_genes=1, num_hits=1
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.fragment_classes[0] == FRAG_UNIQUE

    def test_isoform_ambiguous(self):
        """1 gene, 2 transcripts, NH=1 → FRAG_ISOFORM_AMBIG."""
        # t0 and t1 both map to same gene (need consistent t_to_g)
        t_to_g = np.array([0, 0, 1, 1, 2], dtype=np.int64)
        buf = FragmentBuffer(t_to_g_arr=t_to_g, chunk_size=100)
        buf.append(_make_resolved(
            t_inds=frozenset({0, 1}), n_genes=1, num_hits=1
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.fragment_classes[0] == FRAG_ISOFORM_AMBIG

    def test_gene_ambiguous_multi_gene(self):
        """2 genes, NH=1 → FRAG_GENE_AMBIG."""
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.append(_make_resolved(
            t_inds=frozenset({0, 1}), n_genes=2, num_hits=1
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.fragment_classes[0] == FRAG_GENE_AMBIG

    def test_gene_ambiguous_multimapped(self):
        """1 gene, 1 transcript, NH=3 → FRAG_MULTIMAPPER."""
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.append(_make_resolved(
            t_inds=frozenset({0}), n_genes=1, num_hits=3
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.fragment_classes[0] == FRAG_MULTIMAPPER

    def test_mixed_classes(self):
        """All three classes in one chunk."""
        t_to_g = np.array([0, 0, 1, 1, 2], dtype=np.int64)
        buf = FragmentBuffer(t_to_g_arr=t_to_g, chunk_size=100)
        buf.append(_make_resolved(  # unique
            t_inds=frozenset({0}), n_genes=1, num_hits=1
        ))
        buf.append(_make_resolved(  # isoform-ambig (t0,t1 same gene)
            t_inds=frozenset({0, 1}), n_genes=1, num_hits=1
        ))
        buf.append(_make_resolved(  # gene-ambig (t0,t2 different genes)
            t_inds=frozenset({0, 2}), n_genes=2, num_hits=1
        ))
        buf.append(_make_resolved(  # gene-ambig (multimapped)
            t_inds=frozenset({0}), n_genes=1, num_hits=2
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        fc = chunk.fragment_classes
        assert fc[0] == FRAG_UNIQUE
        assert fc[1] == FRAG_ISOFORM_AMBIG
        assert fc[2] == FRAG_GENE_AMBIG
        assert fc[3] == FRAG_MULTIMAPPER

    def test_consistent_with_ambiguous_mask(self):
        """fragment_classes and ambiguous_mask should agree."""
        t_to_g = np.array([0, 0, 1, 1, 2], dtype=np.int64)
        buf = FragmentBuffer(t_to_g_arr=t_to_g, chunk_size=100)
        buf.append(_make_resolved(
            t_inds=frozenset({0}), n_genes=1, num_hits=1
        ))
        buf.append(_make_resolved(
            t_inds=frozenset({0, 1}), n_genes=1, num_hits=1
        ))
        buf.append(_make_resolved(
            t_inds=frozenset({0, 2}), n_genes=2, num_hits=1
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        fc = chunk.fragment_classes
        mask = chunk.ambiguous_mask
        # FRAG_UNIQUE → not ambiguous
        assert not mask[0]
        assert fc[0] == FRAG_UNIQUE
        # FRAG_ISOFORM_AMBIG → not ambiguous (1 gene, model considers unique)
        assert not mask[1]
        assert fc[1] == FRAG_ISOFORM_AMBIG
        # FRAG_GENE_AMBIG → ambiguous
        assert mask[2]
        assert fc[2] == FRAG_GENE_AMBIG

    def test_chimeric_overrides_other_classes(self):
        """Chimeric fragments should be classified as FRAG_CHIMERIC regardless of gene/hit count."""
        t_to_g = np.array([0, 0, 1, 1, 2], dtype=np.int64)
        buf = FragmentBuffer(t_to_g_arr=t_to_g, chunk_size=100)
        # Chimeric + gene-ambig: chimeric should win
        buf.append(_make_resolved(
            t_inds=frozenset({0, 2}), n_genes=2, num_hits=1,
            chimera_type=ChimeraType.CIS_STRAND_SAME,
        ))
        # Chimeric + multimapper: chimeric should win
        buf.append(_make_resolved(
            t_inds=frozenset({0}), n_genes=1, num_hits=3,
            chimera_type=ChimeraType.TRANS,
        ))
        # Non-chimeric gene-ambig: should remain FRAG_GENE_AMBIG
        buf.append(_make_resolved(
            t_inds=frozenset({0, 2}), n_genes=2, num_hits=1,
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        fc = chunk.fragment_classes
        assert fc[0] == FRAG_CHIMERIC
        assert fc[1] == FRAG_CHIMERIC
        assert fc[2] == FRAG_GENE_AMBIG

    def test_chimera_type_round_trip(self):
        """chimera_type should survive append → finalize → iterate."""
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.append(_make_resolved(chimera_type=ChimeraType.NONE))
        buf.append(_make_resolved(chimera_type=ChimeraType.TRANS))
        buf.append(_make_resolved(chimera_type=ChimeraType.CIS_STRAND_SAME))
        buf.append(_make_resolved(chimera_type=ChimeraType.CIS_STRAND_DIFF))
        buf.finalize()

        frags = list(buf)
        assert frags[0].chimera_type == ChimeraType.NONE
        assert frags[1].chimera_type == ChimeraType.TRANS
        assert frags[2].chimera_type == ChimeraType.CIS_STRAND_SAME
        assert frags[3].chimera_type == ChimeraType.CIS_STRAND_DIFF


# =====================================================================
# FragmentBuffer — ambiguous mask
# =====================================================================


class TestAmbiguousMask:
    def test_unique_fragments(self):
        """Fragments with single gene and num_hits=1 are not ambiguous."""
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        for _ in range(5):
            buf.append(_make_resolved(
                t_inds=frozenset({0, 3}),  # both → g0
                n_genes=1, num_hits=1
            ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        mask = chunk.ambiguous_mask
        assert mask.sum() == 0

    def test_multi_gene_is_ambiguous(self):
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.append(_make_resolved(
            t_inds=frozenset({0, 1}),  # g0, g1
            n_genes=2, num_hits=1,
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.ambiguous_mask[0] is np.True_

    def test_multimapped_is_ambiguous(self):
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.append(_make_resolved(
            t_inds=frozenset({0, 3}),  # both → g0
            n_genes=1, num_hits=3,
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        assert chunk.ambiguous_mask[0] is np.True_

    def test_mixed(self):
        buf = FragmentBuffer(t_to_g_arr=_make_t_to_g_arr(), chunk_size=100)
        buf.append(_make_resolved(
            t_inds=frozenset({0, 3}), n_genes=1, num_hits=1,   # unique
        ))
        buf.append(_make_resolved(
            t_inds=frozenset({0, 1}), n_genes=2, num_hits=1,   # multi-gene
        ))
        buf.append(_make_resolved(
            t_inds=frozenset({0, 3}), n_genes=1, num_hits=2,   # multimapped
        ))
        buf.finalize()

        chunk = list(buf.iter_chunks())[0]
        mask = chunk.ambiguous_mask
        assert not mask[0]
        assert mask[1]
        assert mask[2]


# =====================================================================
# FragmentBuffer — disk spill
# =====================================================================


class TestDiskSpill:
    def test_spill_triggers_above_threshold(self, tmp_path):
        """When in-memory chunks exceed max_memory_bytes, spill to disk."""
        # Use a very small threshold to force spilling
        buf, _ = _fill_buffer(
            300,
            chunk_size=100,
            max_memory_bytes=1,  # 1 byte → always spill
            spill_dir=tmp_path,
        )

        assert buf.n_spilled > 0
        assert buf.memory_bytes == 0 or buf.n_spilled >= buf.n_chunks - 1

    def test_spill_preserves_data(self, tmp_path):
        """Data roundtrips correctly through Arrow IPC spill."""
        buf, fragments = _fill_buffer(
            200,
            chunk_size=100,
            max_memory_bytes=1,  # force all chunks to spill
            spill_dir=tmp_path,
        )

        result = list(buf)
        assert len(result) == 200
        for i, bf in enumerate(result):
            assert bf.frag_lengths[0] == 200 + i
            assert bf.splice_type == int(SpliceType(i % len(SpliceType)))

    def test_cleanup_removes_files(self, tmp_path):
        buf, _ = _fill_buffer(
            200,
            chunk_size=100,
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )

        # Files should exist before cleanup
        assert buf.n_spilled > 0
        spill_dirs = list(tmp_path.iterdir())
        assert len(spill_dirs) > 0

        buf.cleanup()

        # Temp dir should be removed
        remaining = [p for p in tmp_path.iterdir() if p.is_dir()]
        assert len(remaining) == 0

    def test_context_manager_cleanup(self, tmp_path):
        with FragmentBuffer(
            t_to_g_arr=_make_t_to_g_arr(),
            chunk_size=100,
            max_memory_bytes=1,
            spill_dir=tmp_path,
        ) as buf:
            for i in range(200):
                buf.append(_make_resolved(frag_lengths={0: i, 1: i}))
            buf.finalize()
            assert buf.n_spilled > 0

        # After exiting context, spill files should be removed
        remaining = [p for p in tmp_path.iterdir() if p.is_dir()]
        assert len(remaining) == 0

    def test_no_spill_when_disabled(self):
        """max_memory_bytes=0 disables spilling."""
        buf, _ = _fill_buffer(
            300,
            chunk_size=100,
            max_memory_bytes=0,
        )
        assert buf.n_spilled == 0
        assert buf.n_chunks == 3

    def test_no_spill_under_threshold(self):
        """Small buffers should not spill."""
        buf, _ = _fill_buffer(
            50,
            chunk_size=100,
            max_memory_bytes=100 * 1024**2,  # 100 MB
        )
        assert buf.n_spilled == 0


# =====================================================================
# FragmentBuffer — summary
# =====================================================================


class TestBufferSummary:
    def test_summary_structure(self):
        buf, _ = _fill_buffer(150, chunk_size=100)
        s = buf.summary()
        assert s["total_fragments"] == 150
        assert s["n_chunks"] == 2
        assert s["in_memory_chunks"] == 2
        assert s["on_disk_chunks"] == 0
        assert s["memory_bytes"] > 0
        assert isinstance(s["memory_mb"], float)

    def test_summary_with_spill(self, tmp_path):
        buf, _ = _fill_buffer(
            200,
            chunk_size=100,
            max_memory_bytes=1,
            spill_dir=tmp_path,
        )
        s = buf.summary()
        assert s["total_fragments"] == 200
        assert s["on_disk_chunks"] > 0


# =====================================================================
# FragmentBuffer — iter_chunks
# =====================================================================


class TestIterChunks:
    def test_iter_chunks_yields_correct_count(self):
        buf, _ = _fill_buffer(250, chunk_size=100)
        chunks = list(buf.iter_chunks())
        assert len(chunks) == 3
        total = sum(c.size for c in chunks)
        assert total == 250

    def test_iter_chunks_getitem(self):
        buf, _ = _fill_buffer(5, chunk_size=100)
        chunk = list(buf.iter_chunks())[0]
        assert chunk.size == 5
        for i in range(5):
            bf = chunk[i]
            assert bf.frag_lengths[0] == 200 + i

    def test_memory_bytes_positive(self):
        buf, _ = _fill_buffer(100, chunk_size=100)
        chunk = list(buf.iter_chunks())[0]
        assert chunk.memory_bytes > 0
