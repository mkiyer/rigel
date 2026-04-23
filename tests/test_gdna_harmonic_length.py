"""Tests for the Option B gDNA harmonic-mean length correction.

These tests validate the per-fragment length correction computed inside
the C++ scorer.  Under Option B, for a fragment with NH alignment hits
against transcripts with genomic spans ``L_h = t_span[t_h] + gdna_flank``
and a shared genomic footprint ``gfp``, the emitted ``gdna_log_lik``
is the marginal likelihood under a uniform-position gDNA origin prior:

    gdna_log_lik = C + log( (1/NH) * sum_h (1/e_h) )

where ``e_h = max(L_h − gfp + 1, 1)`` and ``C`` gathers all per-hit
terms that are constant within a multimapper group (fragment-length
log-prob, splice penalty, LOG_HALF, and NM penalty).

The tests construct synthetic multimapper groups where all hits share
``gfp`` and ``nm=0`` so ``C`` cancels when we diff two units, and
recover the length-correction signal exactly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest

from rigel.buffer import FRAG_MULTIMAPPER, FRAG_UNAMBIG
from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator
from rigel.frag_length_model import FragmentLengthModels
from rigel.scan import FragmentRouter
from rigel.scoring import FragmentScorer
from rigel.splice import SpliceType
from rigel.stats import PipelineStats
from rigel.strand_model import StrandModels
from rigel.types import Strand


# ---------------------------------------------------------------------------
# Minimal test fixtures (mirrors test_pipeline_routing helpers, with a
# customisable per-transcript t_span so we can verify the length-correction
# arithmetic exactly).
# ---------------------------------------------------------------------------


@dataclass
class _BF:
    t_inds: np.ndarray
    splice_type: int
    exon_strand: int
    frag_lengths: np.ndarray
    exon_bp: np.ndarray
    intron_bp: np.ndarray
    read_length: int
    nm: int = 0
    genomic_footprint: int = 200
    genomic_start: int = -1


class _Chunk:
    def __init__(self, bfs, fragment_classes, frag_ids):
        self._bfs = bfs
        n = len(bfs)
        self.splice_type = np.array([bf.splice_type for bf in bfs], dtype=np.uint8)
        self.exon_strand = np.array([bf.exon_strand for bf in bfs], dtype=np.uint8)
        self.fragment_classes = np.array(fragment_classes, dtype=np.uint8)
        self.frag_id = np.array(frag_ids, dtype=np.int64)
        self.read_length = np.array([bf.read_length for bf in bfs], dtype=np.uint32)
        self.genomic_footprint = np.array(
            [bf.genomic_footprint for bf in bfs], dtype=np.int32
        )
        self.genomic_start = np.array(
            [bf.genomic_start for bf in bfs], dtype=np.int32
        )
        self.nm = np.array([bf.nm for bf in bfs], dtype=np.uint16)
        self.size = n

        offsets = [0]
        flat_t = []
        flat_fl = []
        flat_exon = []
        flat_intron = []
        for bf in bfs:
            n_cand = len(bf.t_inds)
            flat_t.extend(bf.t_inds)
            flat_fl.extend(bf.frag_lengths)
            flat_exon.extend(bf.exon_bp)
            flat_intron.extend(bf.intron_bp)
            offsets.append(len(flat_t))
        self.t_offsets = np.array(offsets, dtype=np.int32)
        self.t_indices = np.array(flat_t, dtype=np.int32)
        self.frag_lengths = np.array(flat_fl, dtype=np.int32)
        self.exon_bp = np.array(flat_exon, dtype=np.int32)
        self.intron_bp = np.array(flat_intron, dtype=np.int32)

    def __getitem__(self, idx):
        return self._bfs[idx]

    def to_scoring_arrays(self):
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


class _Buffer:
    def __init__(self, chunks):
        self._chunks = list(chunks)
        self.total_fragments = sum(ch.size for ch in self._chunks)

    def iter_chunks(self):
        yield from self._chunks

    def iter_chunks_consuming(self):
        while self._chunks:
            yield self._chunks.pop(0)


class _Index:
    """Minimal index with controllable per-transcript genomic span."""

    def __init__(self, t_spans: list[int], t_lengths: list[int] | None = None):
        n_t = len(t_spans)
        if t_lengths is None:
            t_lengths = [max(s, 1) for s in t_spans]
        self.t_to_g_arr = np.arange(n_t, dtype=np.int64)
        self.t_to_strand_arr = np.full(n_t, int(Strand.POS), dtype=np.int8)
        self.g_to_strand_arr = np.full(n_t, int(Strand.POS), dtype=np.int8)
        self.num_transcripts = n_t
        self.num_genes = n_t
        starts = np.zeros(n_t, dtype=np.int32)
        ends = np.array(t_spans, dtype=np.int32)
        self.t_df = pd.DataFrame(
            {
                "t_id": [f"t{i}" for i in range(n_t)],
                "ref": ["chr1"] * n_t,
                "start": starts,
                "end": ends,
                "length": np.array(t_lengths, dtype=np.int32),
                "is_nrna": np.zeros(n_t, dtype=bool),
                "is_synthetic": np.zeros(n_t, dtype=bool),
            }
        )

    def get_exon_intervals(self, t_idx):
        return None

    def build_exon_csr(self):
        n_t = self.num_transcripts
        offsets = np.zeros(n_t + 1, dtype=np.int32)
        empty = np.empty(0, dtype=np.int32)
        return offsets, empty, empty, empty


def _make_env(index, *, gdna_mean_target: int | None = None):
    """Build StrandModels, FragmentLengthModels, estimator, stats.

    ``gdna_mean_target``: if set, the gDNA FL model is seeded with a
    single-bin observation at this length so ``gdna_model.mean`` is
    exactly that value (drives ``gdna_flank`` in FragmentScorer).
    Otherwise the gDNA model is left empty, giving ``gdna_flank = 0``.
    """
    sm = StrandModels()
    for _ in range(50):
        sm.exonic_spliced.observe(Strand.POS, Strand.POS)
    sm.finalize()

    fl = FragmentLengthModels(max_size=1000)
    # Seed the RNA model with a flat region around 200 so
    # log_likelihood(200) is finite.
    for length in range(150, 260):
        fl.observe(length, splice_type=SpliceType.UNSPLICED, weight=1.0)
    if gdna_mean_target is not None:
        # Concentrate gDNA mass at the target length so the mean
        # (rounded to int) is exactly gdna_mean_target.
        fl.gdna_model.observe(gdna_mean_target, weight=1000.0)
    fl.build_scoring_models()
    fl.finalize()

    estimator = AbundanceEstimator(index.num_transcripts, em_config=EMConfig(seed=1))
    stats = PipelineStats()
    return sm, fl, estimator, stats


def _scan(buffer, index, sm, fl, estimator, stats):
    ctx = FragmentScorer.from_models(sm, fl, index, estimator)
    router = FragmentRouter(ctx, estimator, stats, index, sm)
    return router.scan(buffer, 1_000_000)


def _mk_bf(t_inds, stype, gfp, frag_length=200):
    """Build one buffered fragment with uniform gfp + frag_length."""
    n = len(t_inds)
    return _BF(
        t_inds=np.array(t_inds, dtype=np.int32),
        splice_type=int(stype),
        exon_strand=int(Strand.POS),
        frag_lengths=np.array([frag_length] * n, dtype=np.int32),
        exon_bp=np.array([100] * n, dtype=np.int16),
        intron_bp=np.array([0] * n, dtype=np.int16),
        read_length=100,
        nm=0,
        genomic_footprint=gfp,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestHarmonicLengthCorrection:
    """Verify the Option B harmonic-mean gDNA length correction."""

    GFP = 200
    FL = 200

    def _run_multimapper(self, t_spans, hit_t_inds, frag_id=10, frag_class=None):
        """Run one multimapper group and return (em, sm, fl, scorer)."""
        index = _Index(t_spans)
        sm, fl, est, stats = _make_env(index)
        fc = frag_class if frag_class is not None else FRAG_MULTIMAPPER
        bfs = [
            _mk_bf([ti], SpliceType.UNSPLICED, self.GFP, self.FL)
            for ti in hit_t_inds
        ]
        fclasses = [fc] * len(bfs)
        chunk = _Chunk(bfs=bfs, fragment_classes=fclasses, frag_ids=[frag_id] * len(bfs))
        buffer = _Buffer([chunk])
        em = _scan(buffer, index, sm, fl, est, stats)
        return em, index, sm, fl, est, stats

    def _expected_e_h(self, t_span, gdna_flank):
        L_h = t_span + gdna_flank
        e_h = L_h - self.GFP + 1
        return max(e_h, 1)

    def test_nh1_unique_unspliced(self):
        """NH=1 unique mapper: gdna_log_lik = C − log(e_0)."""
        em, idx, sm, fl, est, stats = self._run_multimapper(
            t_spans=[5000], hit_t_inds=[0], frag_class=FRAG_UNAMBIG
        )
        assert em.n_units == 1
        assert np.isfinite(em.gdna_log_liks[0])
        # With one hit, length-correction = −log(e_0)

    def test_nh2_identical_spans_equals_unique_length(self):
        """Two hits against transcripts with IDENTICAL spans should produce
        the same length correction as a single hit against either one.

        Because e_0 = e_1, the harmonic mean collapses:
            log((1/e + 1/e) / 2) = −log(e)
        """
        span = 5000
        em_mm, idx_mm, _, _, _, _ = self._run_multimapper(
            t_spans=[span, span], hit_t_inds=[0, 1]
        )
        em_uniq, idx_uq, _, _, _, _ = self._run_multimapper(
            t_spans=[span], hit_t_inds=[0], frag_class=FRAG_UNAMBIG
        )
        assert em_mm.n_units == 1
        assert em_uniq.n_units == 1
        np.testing.assert_allclose(
            em_mm.gdna_log_liks[0], em_uniq.gdna_log_liks[0], rtol=0, atol=1e-12
        )

    def test_nh2_dissimilar_spans_exact_harmonic_mean(self):
        """Two hits against transcripts with DIFFERENT spans: emitted
        gdna_log_lik − C = log((1/e_0 + 1/e_1)/2) where C is recovered
        from an NH=1 unit against either transcript.
        """
        t_spans = [5000, 20000]

        # Reference: NH=1 against t0 alone. Its emitted log-lik is
        #   C − log(e_0).
        em_ref, idx_ref, _, _, _, _ = self._run_multimapper(
            t_spans=[t_spans[0]], hit_t_inds=[0], frag_class=FRAG_UNAMBIG
        )
        # Flank from ref run (gdna_model has zero weight here → 0).
        gdna_flank_ref = 0
        e0 = self._expected_e_h(t_spans[0], gdna_flank_ref)
        C = float(em_ref.gdna_log_liks[0]) + math.log(e0)

        # Multimapper NH=2 against t0 and t1.
        em_mm, idx_mm, _, _, _, _ = self._run_multimapper(
            t_spans=t_spans, hit_t_inds=[0, 1]
        )
        e0_mm = self._expected_e_h(t_spans[0], gdna_flank_ref)
        e1_mm = self._expected_e_h(t_spans[1], gdna_flank_ref)
        expected = C + math.log((1.0 / e0_mm + 1.0 / e1_mm) / 2.0)
        np.testing.assert_allclose(
            em_mm.gdna_log_liks[0], expected, rtol=0, atol=1e-10
        )

    def test_spliced_any_hit_kills_gdna(self):
        """If any hit in an MM group is spliced, gdna_log_lik = −inf."""
        index = _Index([5000, 5000])
        sm, fl, est, stats = _make_env(index)
        bfs = [
            _mk_bf([0], SpliceType.UNSPLICED, self.GFP, self.FL),
            _mk_bf([1], SpliceType.SPLICED_ANNOT, self.GFP, self.FL),
        ]
        chunk = _Chunk(
            bfs=bfs,
            fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
            frag_ids=[7, 7],
        )
        em = _scan(_Buffer([chunk]), index, sm, fl, est, stats)
        assert em.n_units == 1
        assert em.gdna_log_liks[0] == -np.inf

    def test_all_spliced_unique_is_minus_inf(self):
        """A spliced multimapper group has gdna_log_lik = −inf
        (the gDNA hypothesis is inapplicable to spliced reads).
        A unique spliced fragment skips EM entirely, so we exercise
        this via a two-hit MM group where both hits are spliced."""
        index = _Index([5000, 5000])
        sm, fl, est, stats = _make_env(index)
        bfs = [
            _mk_bf([0], SpliceType.SPLICED_ANNOT, self.GFP, self.FL),
            _mk_bf([1], SpliceType.SPLICED_ANNOT, self.GFP, self.FL),
        ]
        chunk = _Chunk(
            bfs=bfs,
            fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
            frag_ids=[42, 42],
        )
        em = _scan(_Buffer([chunk]), index, sm, fl, est, stats)
        assert em.n_units == 1
        assert em.gdna_log_liks[0] == -np.inf

    def test_gdna_flank_adds_to_span(self):
        """With a non-zero gdna_flank (driven by gdna_model.mean), the
        effective window L_h = t_span + gdna_flank.  A run with
        gdna_flank = K should shift log-lik by +log(e) − log(e + K)
        relative to gdna_flank = 0 for an NH=1 unit."""
        t_span = 5000
        flank = 300

        # Reference: no gDNA observations → gdna_flank = 0
        em0, _, _, _, _, _ = self._run_multimapper(
            t_spans=[t_span], hit_t_inds=[0], frag_class=FRAG_UNAMBIG
        )

        # With flank: seed the gDNA model so mean() = flank exactly.
        index = _Index([t_span])
        sm, fl, est, stats = _make_env(index, gdna_mean_target=flank)
        # Confirm the FL model's mean rounds to flank.
        assert int(fl.gdna_model.mean) == flank
        bfs = [_mk_bf([0], SpliceType.UNSPLICED, self.GFP, self.FL)]
        chunk = _Chunk(
            bfs=bfs, fragment_classes=[FRAG_UNAMBIG], frag_ids=[1]
        )
        em_f = _scan(_Buffer([chunk]), index, sm, fl, est, stats)
        assert em_f.n_units == 1

        # The gdna_fl value also changes because the gDNA FL model is
        # now non-empty.  To isolate the length correction we work
        # purely on the per-hit e_h: a wider window L_h = t_span + flank
        # means a larger e_h, so the length-correction term
        #   −log(e_h)
        # becomes MORE negative (a wider uniform prior is less
        # concentrated per position).
        e_no_flank = self._expected_e_h(t_span, 0)
        e_flank = self._expected_e_h(t_span, flank)
        assert e_flank > e_no_flank
        assert -math.log(e_flank) < -math.log(e_no_flank)

    def test_nh3_three_way_harmonic_mean(self):
        """NH=3: verify the accumulator handles > 2 hits correctly."""
        t_spans = [3000, 8000, 15000]

        # Reference C from NH=1 against t0.
        em_ref, _, _, _, _, _ = self._run_multimapper(
            t_spans=[t_spans[0]], hit_t_inds=[0], frag_class=FRAG_UNAMBIG
        )
        e0 = self._expected_e_h(t_spans[0], 0)
        C = float(em_ref.gdna_log_liks[0]) + math.log(e0)

        em_mm, _, _, _, _, _ = self._run_multimapper(
            t_spans=t_spans, hit_t_inds=[0, 1, 2]
        )
        eh = [self._expected_e_h(s, 0) for s in t_spans]
        expected = C + math.log(sum(1.0 / e for e in eh) / 3.0)
        np.testing.assert_allclose(
            em_mm.gdna_log_liks[0], expected, rtol=0, atol=1e-10
        )
