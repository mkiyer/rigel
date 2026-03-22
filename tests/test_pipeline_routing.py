"""Tests for pipeline routing counters and multimapper shadow behavior."""

from dataclasses import dataclass

import numpy as np
import pytest

from rigel.buffer import (
    FRAG_AMBIG_OPP_STRAND,
    FRAG_AMBIG_SAME_STRAND,
    FRAG_MULTIMAPPER,
    FRAG_UNAMBIG,
)
from rigel.splice import SpliceType
from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator
from rigel.frag_length_model import FragmentLengthModels
from rigel.scoring import FragmentScorer
from rigel.scan import FragmentRouter
from rigel.stats import PipelineStats
from rigel.strand_model import StrandModels
from rigel.types import Strand


@dataclass
class _BF:
    t_inds: np.ndarray
    splice_type: int
    exon_strand: int
    frag_lengths: np.ndarray | None
    exon_bp: np.ndarray
    intron_bp: np.ndarray
    read_length: int
    nm: int = 0
    genomic_footprint: int = 200
    genomic_start: int = -1


class _Chunk:
    """Columnar mock chunk compatible with the C++ native scan path.

    Has ``t_offsets`` and other columnar arrays so ``_scan_native``
    fires.  Also supports ``__getitem__`` for multimapper handling.
    """

    def __init__(self, bfs, fragment_classes, frag_ids):
        self._bfs = bfs
        n = len(bfs)

        # Per-fragment columnar arrays
        self.splice_type = np.array(
            [bf.splice_type for bf in bfs], dtype=np.uint8)
        self.exon_strand = np.array(
            [bf.exon_strand for bf in bfs], dtype=np.uint8)
        self.fragment_classes = np.array(fragment_classes, dtype=np.uint8)
        self.frag_id = np.array(frag_ids, dtype=np.int64)
        self.read_length = np.array(
            [bf.read_length for bf in bfs], dtype=np.uint32)
        self.genomic_footprint = np.array(
            [bf.genomic_footprint for bf in bfs], dtype=np.int32)
        self.genomic_start = np.array(
            [bf.genomic_start for bf in bfs], dtype=np.int32)
        self.nm = np.array([bf.nm for bf in bfs], dtype=np.uint16)
        self.size = n

        # CSR transcript indices (variable-length per fragment)
        offsets = [0]
        flat_t = []
        flat_fl = []
        flat_exon = []
        flat_intron = []
        for bf in bfs:
            n_cand = len(bf.t_inds)
            flat_t.extend(bf.t_inds)
            if bf.frag_lengths is not None:
                flat_fl.extend(bf.frag_lengths)
            else:
                flat_fl.extend([200] * n_cand)
            flat_exon.extend(bf.exon_bp)
            flat_intron.extend(bf.intron_bp)
            offsets.append(len(flat_t))

        self.t_offsets = np.array(offsets, dtype=np.int64)
        self.t_indices = np.array(flat_t, dtype=np.int32)
        self.frag_lengths = np.array(flat_fl, dtype=np.int32)
        self.exon_bp = np.array(flat_exon, dtype=np.int32)
        self.intron_bp = np.array(flat_intron, dtype=np.int32)

    def __getitem__(self, idx):
        return self._bfs[idx]

    def to_scoring_arrays(self):
        return (
            np.ascontiguousarray(self.t_offsets, dtype=np.int64),
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


class _Buffer:
    def __init__(self, chunks):
        self._chunks = chunks
        self.total_fragments = sum(ch.size for ch in chunks)

    def iter_chunks(self):
        yield from self._chunks


class _Index:
    def __init__(self, t_to_g, t_to_strand, g_to_strand):
        import pandas as pd
        self.t_to_g_arr = np.array(t_to_g, dtype=np.int64)
        self.t_to_strand_arr = np.array(t_to_strand, dtype=np.int8)
        self.g_to_strand_arr = np.array(g_to_strand, dtype=np.int8)
        self.num_transcripts = len(t_to_g)
        self.num_genes = len(g_to_strand)
        self.t_df = pd.DataFrame({
            "t_id": [f"t{i}" for i in range(self.num_transcripts)],
            "ref": ["chr1"] * self.num_transcripts,
            "start": np.zeros(self.num_transcripts, dtype=np.int32),
            "end": np.full(self.num_transcripts, 10000, dtype=np.int32),
            "length": np.full(self.num_transcripts, 1000, dtype=np.int32),
            "is_synthetic_nrna": np.zeros(self.num_transcripts, dtype=bool),
        })

    def get_exon_intervals(self, t_idx):
        return None


def _make_env(index):
    strand_models = StrandModels()
    for _ in range(20):
        strand_models.exonic_spliced.observe(Strand.POS, Strand.POS)
    strand_models.finalize()
    frag_length_models = FragmentLengthModels(max_size=1000)
    frag_length_models.observe(200, SpliceType.UNSPLICED)
    estimator = AbundanceEstimator(index.num_transcripts,
                                  em_config=EMConfig(seed=1))
    stats = PipelineStats()
    return strand_models, frag_length_models, estimator, stats


def _scan_em_data(
    buffer, index, strand_models, frag_length_models, estimator, stats,
    log_every=1_000_000, *, gdna_splice_penalties=None,
    overhang_log_penalty=None, mismatch_log_penalty=None, annotations=None,
):
    """Build EM data from a buffer using FragmentScorer + FragmentRouter."""
    ctx = FragmentScorer.from_models(
        strand_models, frag_length_models, index, estimator,
        overhang_log_penalty=overhang_log_penalty,
        mismatch_log_penalty=mismatch_log_penalty,
        gdna_splice_penalties=gdna_splice_penalties,
    )
    builder = FragmentRouter(
        ctx, estimator, stats, index, strand_models,
        annotations=annotations,
    )
    return builder.scan(buffer, log_every)


def test_multimapper_spliced_annot_skips_shadows():
    index = _Index(
        t_to_g=[0, 0],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, frag_length_models, estimator, stats = _make_env(index)

    bfs = [
        _BF(
            t_inds=np.array([0], dtype=np.int32),
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
            intron_bp=np.array([0], dtype=np.int16),
            read_length=100,
        ),
        _BF(
            t_inds=np.array([1], dtype=np.int32),
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
            intron_bp=np.array([0], dtype=np.int16),
            read_length=100,
        ),
    ]
    chunk = _Chunk(
        bfs=bfs,
        fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
        frag_ids=[10, 10],
    )
    buffer = _Buffer([chunk])

    em = _scan_em_data(
        buffer,
        index,
        strand_models,
        frag_length_models,
        estimator,
        stats,
        log_every=1_000_000,
    )

    assert em.n_units == 1
    assert stats.em_routed_multimapper_units == 1
    # Spliced multimappers: only transcript candidates
    assert np.all(em.t_indices < index.num_transcripts)


def test_multimapper_unspliced_adds_nrna_shadows():
    """Unspliced multimappers get transcript candidates in global CSR.

    gDNA candidates are no longer in the global CSR — they are added
    per-locus during locus EM construction. nRNA shadows no longer exist
    (synthetic nRNA transcripts are regular transcripts).
    """
    index = _Index(
        t_to_g=[0, 0],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, frag_length_models, estimator, stats = _make_env(index)

    bfs = [
        _BF(
            t_inds=np.array([0], dtype=np.int32),
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([90], dtype=np.int16),
            intron_bp=np.array([10], dtype=np.int16),
            read_length=100,
        ),
        _BF(
            t_inds=np.array([1], dtype=np.int32),
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([95], dtype=np.int16),
            intron_bp=np.array([5], dtype=np.int16),
            read_length=100,
        ),
    ]
    chunk = _Chunk(
        bfs=bfs,
        fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
        frag_ids=[11, 11],
    )
    buffer = _Buffer([chunk])

    em = _scan_em_data(
        buffer,
        index,
        strand_models,
        frag_length_models,
        estimator,
        stats,
        log_every=1_000_000,
    )

    assert em.n_units == 1
    assert stats.em_routed_multimapper_units == 1
    # All candidates are transcript-level (no nRNA shadows)
    assert np.all(em.t_indices < index.num_transcripts)
    # gDNA candidates NOT in global CSR (added per-locus)
    # gDNA log-likelihood is pre-computed per unit for unspliced
    assert em.gdna_log_liks[0] > -np.inf


def test_route_counters_are_exclusive_per_unit():
    index = _Index(
        t_to_g=[0, 0, 1],
        t_to_strand=[int(Strand.POS), int(Strand.POS), int(Strand.NEG)],
        g_to_strand=[int(Strand.POS), int(Strand.NEG)],
    )
    strand_models, frag_length_models, estimator, stats = _make_env(index)

    bfs = [
        # Deterministic unique: FRAG_UNAMBIG + SPLICED_ANNOT
        _BF(np.array([0], dtype=np.int32), int(SpliceType.SPLICED_ANNOT), int(Strand.POS), np.array([200], dtype=np.int32),
            np.array([100], dtype=np.int16), np.array([0], dtype=np.int16), 100),
        # Unique routed to EM: FRAG_UNAMBIG + UNSPLICED
        _BF(np.array([0], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200], dtype=np.int32),
            np.array([100], dtype=np.int16), np.array([0], dtype=np.int16), 100),
        # Isoform ambiguous
        _BF(np.array([0, 1], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200, 200], dtype=np.int32),
            np.array([95, 90], dtype=np.int16), np.array([5, 10], dtype=np.int16), 100),
        # Gene ambiguous
        _BF(np.array([0, 2], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200, 200], dtype=np.int32),
            np.array([90, 85], dtype=np.int16), np.array([10, 15], dtype=np.int16), 100),
        # Multimapper group (2 hits, same frag_id)
        _BF(np.array([0], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200], dtype=np.int32),
            np.array([95], dtype=np.int16), np.array([5], dtype=np.int16), 100),
        _BF(np.array([2], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200], dtype=np.int32),
            np.array([90], dtype=np.int16), np.array([10], dtype=np.int16), 100),
    ]

    chunk = _Chunk(
        bfs=bfs,
        fragment_classes=[
            FRAG_UNAMBIG,
            FRAG_UNAMBIG,
            FRAG_AMBIG_SAME_STRAND,
            FRAG_AMBIG_OPP_STRAND,
            FRAG_MULTIMAPPER,
            FRAG_MULTIMAPPER,
        ],
        frag_ids=[1, 2, 3, 4, 5, 5],
    )
    buffer = _Buffer([chunk])

    em = _scan_em_data(
        buffer,
        index,
        strand_models,
        frag_length_models,
        estimator,
        stats,
        log_every=1_000_000,
    )

    assert stats.deterministic_unambig_units == 1
    assert stats.em_routed_unambig_units == 1
    assert stats.em_routed_ambig_same_strand_units == 1
    assert stats.em_routed_ambig_opp_strand_units == 1
    assert stats.em_routed_multimapper_units == 1

    total_units = (
        stats.deterministic_unambig_units
        + stats.em_routed_unambig_units
        + stats.em_routed_ambig_same_strand_units
        + stats.em_routed_ambig_opp_strand_units
        + stats.em_routed_multimapper_units
    )
    assert total_units == 5
    assert em.n_units == 4


def test_nm_penalty_discriminates_multimapper_hits():
    """Multimapper with different NM values → NM penalty produces different log-liks.

    Hit A: maps to t0 with NM=0 (perfect match)
    Hit B: maps to t1 with NM=4 (4 mismatches)
    Both have identical geometry (oh=0), strand, fragment length.
    The NM penalty should make t0 score significantly higher than t1.
    """
    import math
    index = _Index(
        t_to_g=[0, 1],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS), int(Strand.POS)],
    )
    strand_models, frag_length_models, estimator, stats = _make_env(index)

    # Hit A: maps to t0 with NM=0
    bf_a = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        exon_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=0,
    )
    # Hit B: maps to t1 with NM=4
    bf_b = _BF(
        t_inds=np.array([1], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        exon_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=4,
    )

    chunk = _Chunk(
        bfs=[bf_a, bf_b],
        fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
        frag_ids=[1, 1],
    )
    buffer = _Buffer([chunk])

    em = _scan_em_data(
        buffer, index, strand_models, frag_length_models,
        estimator, stats, log_every=1_000_000,
    )

    assert em.n_units == 1
    # Both transcripts should be candidates
    mRNA_t = em.t_indices
    mRNA_ll = em.log_liks
    assert len(mRNA_t) == 2
    # t0 (NM=0) should have higher log-lik than t1 (NM=4)
    t0_idx = np.where(mRNA_t == 0)[0][0]
    t1_idx = np.where(mRNA_t == 1)[0][0]
    # Default mismatch_alpha = 0.1 → log(0.1) ≈ -2.3 per mismatch
    # 4 mismatches → Δ ≈ -9.2
    assert mRNA_ll[t0_idx] > mRNA_ll[t1_idx]
    assert mRNA_ll[t0_idx] - mRNA_ll[t1_idx] == pytest.approx(
        -4 * math.log(0.1), abs=0.01,
    )


def test_nm_penalty_zero_when_disabled():
    """With mismatch_alpha=1.0 (no penalty), NM has no effect on log-liks.

    Two identical fragments differing only in NM should score identically.
    """
    index = _Index(
        t_to_g=[0],
        t_to_strand=[int(Strand.POS)],
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, frag_length_models, estimator_a, stats_a = _make_env(index)
    _, _, estimator_b, stats_b = _make_env(index)

    bf_nm0 = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        exon_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=0,
    )
    bf_nm5 = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        exon_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=5,
    )

    # Run with mismatch_alpha=1.0 (disabled)
    import math
    mismatch_lp = math.log(1.0)  # = 0.0

    chunk_a = _Chunk(bfs=[bf_nm0], fragment_classes=[FRAG_UNAMBIG], frag_ids=[1])
    chunk_b = _Chunk(bfs=[bf_nm5], fragment_classes=[FRAG_UNAMBIG], frag_ids=[1])

    em_a = _scan_em_data(
        _Buffer([chunk_a]), index, strand_models, frag_length_models,
        estimator_a, stats_a, log_every=1_000_000,
        mismatch_log_penalty=mismatch_lp,
    )
    em_b = _scan_em_data(
        _Buffer([chunk_b]), index, strand_models, frag_length_models,
        estimator_b, stats_b, log_every=1_000_000,
        mismatch_log_penalty=mismatch_lp,
    )

    # Log-liks should be identical since penalty is 0
    assert em_a.log_liks[0] == pytest.approx(em_b.log_liks[0])
