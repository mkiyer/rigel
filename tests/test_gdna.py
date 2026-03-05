"""Tests for gDNA contamination modeling in locus-level EM.

Verifies that:
- _score_gdna_candidate computes correct log-likelihoods
- gDNA shadow component competes with mRNA in locus EM
- gDNA assignments propagate to correct counters
- _compute_gdna_rate_from_strand returns correct rates
- _compute_nrna_init returns correct nRNA priors
- StrandModels container works correctly
- PipelineStats gDNA fields are correct
"""

import numpy as np
import pytest

from hulkrna.splice import SpliceType, SpliceStrandCol
from hulkrna.config import EMConfig
from hulkrna.estimator import AbundanceEstimator
from hulkrna.scoring import score_gdna_standalone as _score_gdna_candidate, GDNA_SPLICE_PENALTIES as _GDNA_SPLICE_PENALTIES
from hulkrna.frag_length_model import FragmentLengthModels
from hulkrna.strand_model import StrandModel, StrandModels
from hulkrna.types import Strand

from conftest import _UNSPLICED_SENSE, _make_locus_em_data, _run_and_assign


# =====================================================================
# Helpers
# =====================================================================


def _make_frag_length_models(n_obs=200, size=250):
    """FragmentLengthModels with observations for all categories + intergenic."""
    im = FragmentLengthModels()
    for _ in range(n_obs):
        im.observe(size, SpliceType.SPLICED_ANNOT)
        im.observe(size, SpliceType.UNSPLICED)
        im.observe(size, SpliceType.SPLICED_UNANNOT)
        im.observe(size, None)  # intergenic
    return im


def _make_strand_models_default():
    """StrandModels with untrained intergenic model (P_strand=0.5)."""
    return StrandModels()


def _make_strand_models_with_ss(ss: float) -> StrandModels:
    """Create StrandModels with a target strand specificity.

    SS = max(p_r1_sense, 1 - p_r1_sense).  For SS >= 0.5, we
    observe ``n_sense`` + ``n_anti`` reads that produce the desired SS.
    """
    sm = StrandModels()
    total = 1000
    n_sense = int(ss * total)
    n_anti = total - n_sense
    for _ in range(n_sense):
        sm.exonic_spliced.observe(Strand.POS, Strand.POS)
    for _ in range(n_anti):
        sm.exonic_spliced.observe(Strand.POS, Strand.NEG)
    return sm


# =====================================================================
# _score_gdna_candidate
# =====================================================================


class TestScoreGDNA:
    """Tests for pipeline._score_gdna_candidate."""

    def test_unspliced_no_penalty(self):
        """UNSPLICED gets splice_penalty=1.0 → log(1)=0."""
        im = _make_frag_length_models()
        cat = int(SpliceType.UNSPLICED)
        score = _score_gdna_candidate(
            int(Strand.POS), cat, 250,
            im,
        )
        # Uses global insert model
        global_model = im.global_model
        assert score < 0
        assert score == pytest.approx(
            np.log(0.5) + global_model.log_likelihood(250)
        )

    def test_spliced_unannot_penalty(self):
        """SPLICED_UNANNOT gets heavy penalty (default 0.01)."""
        im = _make_frag_length_models()
        score_unannot = _score_gdna_candidate(
            int(Strand.POS), int(SpliceType.SPLICED_UNANNOT), 250,
            im,
        )
        score_unspliced = _score_gdna_candidate(
            int(Strand.POS), int(SpliceType.UNSPLICED), 250,
            im,
        )
        assert score_unannot < score_unspliced
        penalty_diff = score_unspliced - score_unannot
        assert penalty_diff == pytest.approx(-np.log(0.01), rel=0.01)

    def test_custom_penalties(self):
        """Custom penalty dict overrides defaults."""
        im = _make_frag_length_models()
        cat = int(SpliceType.UNSPLICED)
        custom = {SpliceType.UNSPLICED: 0.5}
        score = _score_gdna_candidate(
            int(Strand.POS), cat, 250,
            im,
            gdna_splice_penalties=custom,
        )
        global_model = im.global_model
        expected = np.log(0.5) + global_model.log_likelihood(250) + np.log(0.5)
        assert score == pytest.approx(expected)

    def test_zero_frag_length(self):
        """frag_length=0 → P_insert=1.0 (log=0)."""
        im = _make_frag_length_models()
        score = _score_gdna_candidate(
            int(Strand.POS), int(SpliceType.UNSPLICED), 0,
            im,
        )
        assert score == pytest.approx(np.log(0.5))


# =====================================================================
# gDNA component in locus EM
# =====================================================================


class TestLocusGDNATheta:
    """Verify gDNA competes with mRNA in locus EM theta."""

    def test_theta_has_gdna_component(self):
        """Theta has 2*n_t + 1 elements (mRNA + nRNA + 1 gDNA)."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0]] * 10,
            num_transcripts=3,
            include_gdna=True,
            gdna_init=10.0,
        )
        theta, _alpha = rc.run_locus_em(lem, em_iterations=1)
        # 2*3 + 1 = 7
        assert theta.shape == (7,)
        assert theta.sum() == pytest.approx(1.0)

    def test_high_gdna_init_biases_theta(self):
        """Large gdna_init biases gDNA theta upward."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0]] * 100,
            num_transcripts=2,
            include_gdna=True,
            gdna_init=500.0,
        )
        theta, _alpha = rc.run_locus_em(lem, em_iterations=10)

        gdna_idx = 2 * 2  # n_t = 2
        assert theta[gdna_idx] > theta[0], (
            f"gDNA theta ({theta[gdna_idx]:.4f}) should dominate "
            f"over t0 ({theta[0]:.4f})"
        )

    def test_no_gdna_init_minimal_absorption(self):
        """With zero gdna_init and weak gDNA likelihood, gDNA theta is small."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 1000.0
        lem = _make_locus_em_data(
            [[0]] * 100,
            num_transcripts=2,
            rc=rc,
            include_gdna=True,
            gdna_init=0.0,
            gdna_log_lik=-5.0,
        )
        theta, _alpha = rc.run_locus_em(lem, em_iterations=10)

        gdna_idx = 2 * 2
        assert theta[gdna_idx] < 0.01, f"gDNA theta too high: {theta[gdna_idx]:.4f}"


# =====================================================================
# gDNA assignment routing in locus EM
# =====================================================================


class TestLocusGDNAAssignment:
    """Verify gDNA assignments go to pool_counts, not em_counts."""

    def test_gdna_assignment_not_in_em_counts(self):
        """Fragments assigned to gDNA are counted in pool_counts return."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0]] * 100,
            log_liks_per_unit=[[-10.0]] * 100,
            num_transcripts=2,
            include_gdna=True,
            gdna_init=1000.0,
            gdna_log_lik=0.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert gdna_count > 0
        total = rc.em_counts.sum() + rc.nrna_em_counts.sum() + gdna_count
        assert total == pytest.approx(100.0)

    def test_total_counts_preserved(self):
        """em_counts + nrna_em + gdna == n_units."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=3,
            include_nrna=True,
            include_gdna=True,
            gdna_init=25.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        total = rc.em_counts.sum() + rc.nrna_em_counts.sum() + gdna_count
        assert total == pytest.approx(200.0)

    def test_strong_transcript_signal_beats_gdna(self):
        """When transcript likelihood >> gDNA, most go to transcript."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 500.0
        lem = _make_locus_em_data(
            [[0]] * 500,
            log_liks_per_unit=[[0.0]] * 500,
            num_transcripts=2,
            rc=rc,
            include_gdna=True,
            gdna_init=1.0,
            gdna_log_lik=-20.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert rc.em_counts[0].sum() > 490
        assert gdna_count < 10

    def test_high_gdna_init_absorbs_fragments(self):
        """With large gdna_init, gDNA takes larger share."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0]] * 500,
            num_transcripts=2,
            include_gdna=True,
            gdna_init=500.0,
            gdna_log_lik=0.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert gdna_count > 100, (
            f"Expected gDNA to absorb many fragments, got {gdna_count}"
        )


# =====================================================================
# gDNA locus attribution
# =====================================================================


class TestGDNALocusAttribution:
    """When a fragment goes to gDNA, locus_counts tracks it."""

    def test_locus_counts_populated(self):
        """gDNA assignments populate gdna_locus_counts."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        cc = _UNSPLICED_SENSE
        lem = _make_locus_em_data(
            [[0]] * 100,
            log_liks_per_unit=[[-10.0]] * 100,
            count_cols_per_unit=[[cc]] * 100,
            num_transcripts=2,
            include_gdna=True,
            gdna_init=1000.0,
            gdna_log_lik=0.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert gdna_count > 0
        assert rc.gdna_locus_counts[0].sum() > 0

    def test_locus_counts_zero_when_no_gdna(self):
        """No gDNA assignments → locus counts stay zero."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 1000.0
        lem = _make_locus_em_data(
            [[0]] * 50,
            log_liks_per_unit=[[0.0]] * 50,
            num_transcripts=2,
            rc=rc,
            include_gdna=True,
            gdna_init=0.0,
            gdna_log_lik=-50.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        # With fractional assignment, a vanishingly small posterior may
        # leak to gDNA even when it's extremely unlikely (log_lik=-50).
        assert gdna_count < 1e-10
        assert rc.gdna_locus_counts.sum() < 1e-10


# =====================================================================
# gDNA properties on AbundanceEstimator
# =====================================================================


class TestGDNAProperties:
    def test_gdna_total_from_em_total(self):
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc._gdna_em_total = 80.0
        assert rc.gdna_total == 80.0

    def test_gdna_contamination_rate(self):
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc._gdna_em_total = 10.0
        rc.unambig_counts[0, 0] = 80.0
        rc.em_counts[0, 0] = 10.0
        assert rc.gdna_contamination_rate == pytest.approx(0.1)

    def test_contamination_rate_zero_when_empty(self):
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        assert rc.gdna_contamination_rate == 0.0


# =====================================================================
# Stats integration
# =====================================================================


class TestStatsGDNA:
    """Verify PipelineStats gDNA field changes."""

    def test_intergenic_split(self):
        from hulkrna.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 100
        stats.n_intergenic_spliced = 10
        assert stats.n_intergenic == 110

    def test_gdna_unambig_property(self):
        from hulkrna.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 80
        stats.n_intergenic_spliced = 20
        assert stats.n_gdna_unambig == 100

    def test_gdna_total_property(self):
        from hulkrna.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 80
        stats.n_intergenic_spliced = 20
        stats.n_gdna_em = 50
        assert stats.n_gdna_total == 150

    def test_to_dict_includes_gdna(self):
        from hulkrna.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 10
        stats.n_intergenic_spliced = 5
        stats.n_gdna_em = 3

        d = stats.to_dict()
        assert d["n_intergenic"] == 15
        assert d["n_gdna_unambig"] == 15
        assert d["n_gdna_total"] == 18


# =====================================================================
# _compute_gdna_rate_from_strand
# =====================================================================


class TestComputeGdnaRate:
    """Tests for pipeline._compute_gdna_rate_from_strand."""

    def test_zero_counts_returns_zero(self):
        from hulkrna.locus import compute_gdna_rate_from_strand as _compute_gdna_rate_from_strand
        rate = _compute_gdna_rate_from_strand(0.0, 0.0, 1.0)
        assert rate == 0.0

    def test_no_antisense_gives_zero(self):
        from hulkrna.locus import compute_gdna_rate_from_strand as _compute_gdna_rate_from_strand
        rate = _compute_gdna_rate_from_strand(100.0, 0.0, 1.0)
        assert rate == 0.0

    def test_antisense_at_perfect_ss(self):
        """At SS=1.0, G = 2 × antisense, rate = G / (S + A)."""
        from hulkrna.locus import compute_gdna_rate_from_strand as _compute_gdna_rate_from_strand
        # S=100, A=50, SS=1.0 → G=2*50=100, rate=100/150
        rate = _compute_gdna_rate_from_strand(100.0, 50.0, 1.0)
        expected = (2.0 * 50.0) / (100.0 + 50.0)
        assert rate == pytest.approx(expected, abs=0.01)

    def test_strand_correction_at_imperfect_ss(self):
        """At SS < 1.0, exact formula with strand correction."""
        from hulkrna.locus import compute_gdna_rate_from_strand as _compute_gdna_rate_from_strand
        # S=200, A=50, SS=0.9
        ss = 0.9
        denom_ss = 2.0 * ss - 1.0
        numerator = 50.0 * ss - 200.0 * (1.0 - ss)
        g = 2.0 * max(0.0, numerator / denom_ss)
        expected_rate = g / (200.0 + 50.0)
        rate = _compute_gdna_rate_from_strand(200.0, 50.0, ss)
        assert rate == pytest.approx(expected_rate, abs=0.01)

    def test_returns_zero_at_low_ss(self):
        """At SS ≤ 0.6, strand info too weak → returns 0."""
        from hulkrna.locus import compute_gdna_rate_from_strand as _compute_gdna_rate_from_strand
        # denom_ss = 2*0.5-1 = 0 ≤ 0.2 → returns 0
        rate = _compute_gdna_rate_from_strand(500.0, 500.0, 0.5)
        assert rate == 0.0
        # SS=0.6 → denom_ss=0.2 ≤ 0.2 → returns 0
        rate2 = _compute_gdna_rate_from_strand(500.0, 500.0, 0.6)
        assert rate2 == 0.0

    def test_clamps_at_zero(self):
        """Corrected G is clamped at zero (no negative rate)."""
        from hulkrna.locus import compute_gdna_rate_from_strand as _compute_gdna_rate_from_strand
        # S=1000, A=10, SS=0.8 → G = 2(10*0.8 - 1000*0.2)/0.6 < 0 → 0
        rate = _compute_gdna_rate_from_strand(1000.0, 10.0, 0.8)
        assert rate == 0.0

    def test_rate_clamped_at_one(self):
        """Rate is clamped to [0, 1]."""
        from hulkrna.locus import compute_gdna_rate_from_strand as _compute_gdna_rate_from_strand
        # Extreme antisense dominance
        rate = _compute_gdna_rate_from_strand(1.0, 1000.0, 0.99)
        assert rate <= 1.0


# =====================================================================
# _compute_nrna_init (transcript-level intronic evidence)
# =====================================================================


class TestComputeNrnaInit:
    """Tests for pipeline._compute_nrna_init."""

    def test_zero_counts_returns_zeros(self):
        """No intronic counts → all inits = 0."""
        from hulkrna.locus import compute_nrna_init as _compute_nrna_init

        sense = np.zeros(3, dtype=np.float64)
        anti = np.zeros(3, dtype=np.float64)
        spans = np.array([10000.0, 5000.0, 1000.0])
        exonic_len = np.array([699.0, 499.0, 800.0])
        sm = _make_strand_models_with_ss(1.0)
        nrna = _compute_nrna_init(sense, anti, spans, exonic_len, 200.0, sm)
        np.testing.assert_array_equal(nrna, [0.0, 0.0, 0.0])

    def test_sense_excess_produces_nrna(self):
        """Intronic sense > antisense → nrna_init > 0."""
        from hulkrna.locus import compute_nrna_init as _compute_nrna_init

        sense = np.array([100.0, 50.0])
        anti = np.array([20.0, 10.0])
        spans = np.array([10000.0, 5000.0])
        exonic_len = np.array([699.0, 499.0])
        sm = _make_strand_models_with_ss(1.0)
        nrna = _compute_nrna_init(sense, anti, spans, exonic_len, 200.0, sm)
        # At SS≈1.0, formula ≈ (sense - anti) / 1.0 ≈ sense - anti
        assert nrna[0] == pytest.approx(80.0, abs=2.0)  # 100 - 20
        assert nrna[1] == pytest.approx(40.0, abs=2.0)  # 50 - 10

    def test_antisense_excess_clamped_to_zero(self):
        """Intronic antisense > sense → nrna_init = 0 (clamped)."""
        from hulkrna.locus import compute_nrna_init as _compute_nrna_init

        sense = np.array([10.0])
        anti = np.array([50.0])
        spans = np.array([10000.0])
        exonic_len = np.array([699.0])
        sm = _make_strand_models_with_ss(1.0)
        nrna = _compute_nrna_init(sense, anti, spans, exonic_len, 200.0, sm)
        assert nrna[0] == 0.0

    def test_single_exon_zeroed(self):
        """Single-exon transcripts (span ≈ exonic length) → nrna = 0."""
        from hulkrna.locus import compute_nrna_init as _compute_nrna_init

        sense = np.array([100.0])
        anti = np.array([10.0])
        # transcript_span = exonic length (single exon, no introns)
        exonic_len = np.array([699.0])
        mean_frag = 200.0
        spans = np.array([699.0])  # no introns
        sm = _make_strand_models_with_ss(1.0)
        nrna = _compute_nrna_init(sense, anti, spans, exonic_len, mean_frag, sm)
        assert nrna[0] == 0.0

    def test_multi_exon_not_zeroed(self):
        """Multi-exon transcripts (span >> exonic length) → nrna preserved."""
        from hulkrna.locus import compute_nrna_init as _compute_nrna_init

        sense = np.array([100.0])
        anti = np.array([10.0])
        exonic_len = np.array([699.0])
        mean_frag = 200.0
        spans = np.array([5000.0])  # large introns
        sm = _make_strand_models_with_ss(1.0)
        nrna = _compute_nrna_init(sense, anti, spans, exonic_len, mean_frag, sm)
        assert nrna[0] == pytest.approx(90.0, abs=2.0)

    def test_output_shape(self):
        """Returns array of shape (num_transcripts,)."""
        from hulkrna.locus import compute_nrna_init as _compute_nrna_init

        sense = np.zeros(5, dtype=np.float64)
        anti = np.zeros(5, dtype=np.float64)
        spans = np.full(5, 10000.0)
        exonic_len = np.full(5, 699.0)
        sm = _make_strand_models_with_ss(1.0)
        nrna = _compute_nrna_init(sense, anti, spans, exonic_len, 200.0, sm)
        assert nrna.shape == (5,)
        assert nrna.dtype == np.float64

    def test_ss_correction_moderate(self):
        """At SS=0.9, nRNA init is corrected upward by 1/(2SS-1)."""
        from hulkrna.locus import compute_nrna_init as _compute_nrna_init

        sense = np.array([100.0])
        anti = np.array([20.0])
        spans = np.array([10000.0])
        exonic_len = np.array([699.0])
        sm = _make_strand_models_with_ss(0.9)
        nrna = _compute_nrna_init(sense, anti, spans, exonic_len, 200.0, sm)
        ss = sm.exonic_spliced.strand_specificity
        denom = 2.0 * ss - 1.0
        expected = (100.0 - 20.0) / denom
        assert nrna[0] == pytest.approx(expected, abs=1.0)

    def test_ss_at_or_below_threshold_returns_zeros(self):
        """At SS ≤ 0.6, strand info too weak → returns zeros."""
        from hulkrna.locus import compute_nrna_init as _compute_nrna_init

        sense = np.array([100.0])
        anti = np.array([50.0])
        spans = np.array([10000.0])
        exonic_len = np.array([699.0])
        # SS=0.5
        sm = _make_strand_models_with_ss(0.5)
        nrna = _compute_nrna_init(sense, anti, spans, exonic_len, 200.0, sm)
        assert nrna[0] == 0.0
        # SS=0.6 (boundary)
        sm = _make_strand_models_with_ss(0.6)
        nrna = _compute_nrna_init(sense, anti, spans, exonic_len, 200.0, sm)
        assert nrna[0] == 0.0


# =====================================================================
# Locus gDNA behavior — integration tests
# =====================================================================


class TestLocusGDNABehavior:
    """Integration tests for gDNA behavior in locus EM."""

    def test_two_transcripts_share_one_gdna_shadow(self):
        """Two transcripts in one locus share a single gDNA component."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        lem = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=3,
            include_gdna=True,
            gdna_init=100.0,
        )
        theta, pool_counts = _run_and_assign(rc, lem, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        # n_components = 2*3 + 1 = 7, gDNA at index 6
        assert theta.shape == (7,)
        assert theta[6] > 0  # gDNA has non-zero theta
        total = rc.em_counts.sum() + rc.nrna_em_counts.sum() + gdna_count
        assert total == pytest.approx(200.0)

    def test_gdna_init_determines_absorption(self):
        """Higher gdna_init → more fragments absorbed by gDNA."""
        # Low gdna_init + strong RNA prior → no gDNA absorption
        rc_low = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc_low.unambig_counts[0, _UNSPLICED_SENSE] = 500.0
        lem_low = _make_locus_em_data(
            [[0]] * 200,
            num_transcripts=2,
            rc=rc_low,
            include_gdna=True,
            gdna_init=1.0,
            gdna_log_lik=-5.0,
        )
        _, pc_low = _run_and_assign(rc_low, lem_low, em_iterations=10)
        gc_low = pc_low["gdna"]

        # High gdna_init → gDNA absorbs fragments
        rc_high = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        lem_high = _make_locus_em_data(
            [[0]] * 200,
            num_transcripts=2,
            include_gdna=True,
            gdna_init=1000.0,
            gdna_log_lik=0.0,
        )
        _, pc_high = _run_and_assign(rc_high, lem_high, em_iterations=10)
        gc_high = pc_high["gdna"]

        assert gc_high > gc_low
        assert gc_high > 100  # most fragments go to gDNA
        assert gc_low < 0.1   # RNA dominates (tiny fractional leak ok)

    def test_gdna_em_count_property(self):
        """AbundanceEstimator.gdna_em_count accumulates across loci."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        assert rc.gdna_em_count == 0.0

        rc._gdna_em_total = 42.0
        assert rc.gdna_em_count == 42.0


# =====================================================================
# StrandModels container
# =====================================================================


class TestStrandModelsContainer:
    """Tests for the simplified StrandModels container."""

    def test_default_construction(self):
        from hulkrna.strand_model import StrandModels

        sm = StrandModels()
        assert sm.exonic_spliced.n_observations == 0
        assert sm.exonic.n_observations == 0
        assert sm.intergenic.n_observations == 0

    def test_delegation_to_exonic_spliced(self):
        from hulkrna.strand_model import StrandModels
        from hulkrna.types import Strand

        sm = StrandModels()
        for _ in range(100):
            sm.exonic_spliced.observe(Strand.POS, Strand.POS)
        for _ in range(2):
            sm.exonic_spliced.observe(Strand.POS, Strand.NEG)

        assert sm.p_r1_sense == sm.exonic_spliced.p_r1_sense
        assert sm.strand_specificity == sm.exonic_spliced.strand_specificity
        assert sm.read1_sense == sm.exonic_spliced.read1_sense
        assert sm.n_observations == sm.exonic_spliced.n_observations
        assert sm.strand_likelihood(Strand.POS, Strand.POS) == sm.exonic_spliced.strand_likelihood(Strand.POS, Strand.POS)

    def test_to_dict_structure(self):
        from hulkrna.strand_model import StrandModels

        sm = StrandModels()
        d = sm.to_dict()
        assert "exonic_spliced" in d
        assert "diagnostics" in d
        assert "exonic" in d["diagnostics"]
        assert "intergenic" in d["diagnostics"]
        assert "observations" in d["exonic_spliced"]

    def test_write_json(self, tmp_path):
        from hulkrna.strand_model import StrandModels
        from hulkrna.types import Strand
        import json

        sm = StrandModels()
        for _ in range(50):
            sm.exonic_spliced.observe(Strand.POS, Strand.POS)

        out = tmp_path / "strand_models.json"
        sm.write_json(out)
        with open(out) as f:
            data = json.load(f)
        assert "strand_models" in data
        assert "exonic_spliced" in data["strand_models"]
        assert "diagnostics" in data["strand_models"]

    def test_independent_models(self):
        """Each region's model is independent."""
        from hulkrna.strand_model import StrandModels
        from hulkrna.types import Strand

        sm = StrandModels()
        sm.exonic_spliced.observe(Strand.POS, Strand.POS)
        sm.exonic.observe(Strand.POS, Strand.NEG)
        sm.intergenic.observe(Strand.NEG, Strand.NEG)

        assert sm.exonic_spliced.n_observations == 1
        assert sm.exonic.n_observations == 1
        assert sm.intergenic.n_observations == 1

    def test_kappa_prior_on_finalize(self):
        """finalize() applies κ prior to exonic_spliced."""
        from hulkrna.strand_model import StrandModels
        from hulkrna.types import Strand

        sm = StrandModels()
        sm.strand_prior_kappa = 10.0
        for _ in range(10):
            sm.exonic_spliced.observe(Strand.POS, Strand.POS)
        sm.finalize()
        # κ=10 → alpha_prior=5, beta_prior=5
        # p_r1_sense = (5 + 10) / (5 + 10 + 5) = 15/20 = 0.75
        assert sm.p_r1_sense == pytest.approx(0.75)


# =====================================================================
# compute_gdna_rate_hybrid
# =====================================================================


class TestComputeGdnaRateHybrid:
    """Tests for the hybrid density+strand gDNA rate estimator."""

    def test_zero_counts_returns_zero(self):
        from hulkrna.locus import compute_gdna_rate_hybrid
        rate, ev = compute_gdna_rate_hybrid(0.0, 0.0, 1000.0, 1.0, 0.001)
        assert rate == 0.0
        assert ev == 0.0

    def test_strand_dominated_at_high_ss(self):
        """At SS=1.0, W=1 → pure strand estimate."""
        from hulkrna.locus import compute_gdna_rate_hybrid
        rate, ev = compute_gdna_rate_hybrid(100.0, 50.0, 1000.0, 1.0, 0.001)
        # Strand rate: G = 2*50 = 100, rate = 100/150
        expected_strand = (2.0 * 50.0) / 150.0
        assert rate == pytest.approx(expected_strand, abs=0.01)
        assert ev == 150.0

    def test_density_dominated_at_low_ss(self):
        """At SS=0.5, W≈0 → pure density estimate."""
        from hulkrna.locus import compute_gdna_rate_hybrid
        # SS=0.5 → denom_ss=0 → W=0
        # density: intergenic_density * exonic_bp / total = 0.001 * 1000 / 200
        rate, ev = compute_gdna_rate_hybrid(100.0, 100.0, 1000.0, 0.5, 0.001)
        expected = 0.001 * 1000.0 / 200.0
        assert rate == pytest.approx(expected, abs=0.01)
        assert ev == 200.0

    def test_density_zero_falls_back_to_strand(self):
        """With intergenic_density=0, only strand contributes."""
        from hulkrna.locus import compute_gdna_rate_hybrid
        rate, ev = compute_gdna_rate_hybrid(100.0, 50.0, 1000.0, 0.9, 0.0)
        # Should match strand-only
        from hulkrna.locus import compute_gdna_rate_from_strand
        strand_rate = compute_gdna_rate_from_strand(100.0, 50.0, 0.9)
        # W = (2*0.9-1)^2 = 0.64 → rate = 0.64*strand + 0.36*0
        W = (2.0 * 0.9 - 1.0) ** 2
        expected = W * strand_rate
        assert rate == pytest.approx(expected, abs=0.01)

    def test_unstranded_with_density_nonzero(self):
        """Unstranded (SS=0.5) with intergenic density gives nonzero rate."""
        from hulkrna.locus import compute_gdna_rate_hybrid
        rate, ev = compute_gdna_rate_hybrid(500.0, 500.0, 10000.0, 0.5, 0.01)
        # density component: 0.01 * 10000 / 1000 = 0.1
        expected = 0.01 * 10000.0 / 1000.0
        assert rate == pytest.approx(expected, abs=0.01)
        assert rate > 0.0  # Critical: nonzero for unstranded

    def test_unstranded_no_density_returns_zero(self):
        """Unstranded with no density → both components zero."""
        from hulkrna.locus import compute_gdna_rate_hybrid
        rate, _ = compute_gdna_rate_hybrid(500.0, 500.0, 10000.0, 0.5, 0.0)
        assert rate == 0.0

    def test_rate_clamped_to_unit_interval(self):
        """Rate is always in [0, 1]."""
        from hulkrna.locus import compute_gdna_rate_hybrid
        # Huge density → would exceed 1.0 without clamping
        rate, _ = compute_gdna_rate_hybrid(1.0, 1.0, 100000.0, 0.5, 1.0)
        assert 0.0 <= rate <= 1.0

    def test_hybrid_blend_at_moderate_ss(self):
        """At intermediate SS, both components contribute."""
        from hulkrna.locus import compute_gdna_rate_hybrid
        ss = 0.8  # W = (2*0.8-1)^2 = 0.36
        rate, _ = compute_gdna_rate_hybrid(100.0, 50.0, 1000.0, ss, 0.01)
        # Both strand and density components should be nonzero
        assert rate > 0.0

    def test_returns_tuple(self):
        """Returns (rate, evidence) tuple."""
        from hulkrna.locus import compute_gdna_rate_hybrid
        result = compute_gdna_rate_hybrid(100.0, 50.0, 1000.0, 1.0, 0.0)
        assert isinstance(result, tuple)
        assert len(result) == 2


# =====================================================================
# MoM κ auto-estimation in gDNA EB system
# =====================================================================


class TestGdnaMoMKappa:
    """Tests for MoM κ auto-estimation in compute_eb_gdna_priors."""

    def test_auto_kappa_matches_estimate_kappa(self):
        """Auto-estimated κ should match estimate_kappa on same data."""
        from hulkrna.estimator import estimate_kappa

        # Create synthetic gDNA rates and evidence
        rates = np.array([0.1, 0.15, 0.12, 0.08, 0.2, 0.11] * 5)
        evidence = np.full(30, 100.0)
        kappa = estimate_kappa(rates, evidence, 30.0)
        assert 2.0 <= kappa <= 200.0

    def test_fallback_when_few_loci(self):
        """With few loci, MoM returns fallback κ."""
        from hulkrna.estimator import estimate_kappa

        rates = np.array([0.1, 0.2])
        evidence = np.array([100.0, 100.0])
        kappa = estimate_kappa(
            rates, evidence, 30.0,
            kappa_fallback=5.0, kappa_min_obs=20,
        )
        assert kappa == 5.0  # Fallback: only 2 loci < 20


# =====================================================================
# Unstranded gDNA behavior
# =====================================================================


class TestUnstrandedGdna:
    """Verify gDNA rate estimation works for unstranded libraries."""

    def test_hybrid_nonzero_for_unstranded_with_density(self):
        """With intergenic density, unstranded gives nonzero gDNA rate."""
        from hulkrna.locus import compute_gdna_rate_hybrid

        # SS=0.5 (unstranded), with intergenic density
        rate, _ = compute_gdna_rate_hybrid(
            1000.0, 1000.0, 50000.0, 0.5, 0.005,
        )
        assert rate > 0.0
        # Expected: density_rate = 0.005 * 50000 / 2000 = 0.125
        expected = 0.005 * 50000.0 / 2000.0
        assert rate == pytest.approx(expected, abs=0.01)

    def test_strand_only_zero_for_unstranded(self):
        """Without density, unstranded gives zero gDNA rate."""
        from hulkrna.locus import compute_gdna_rate_from_strand

        rate = compute_gdna_rate_from_strand(1000.0, 1000.0, 0.5)
        assert rate == 0.0

    def test_hybrid_zero_for_unstranded_no_density(self):
        """Unstranded with no density falls back to zero (backward compat)."""
        from hulkrna.locus import compute_gdna_rate_hybrid

        rate, _ = compute_gdna_rate_hybrid(
            1000.0, 1000.0, 50000.0, 0.5, 0.0,
        )
        assert rate == 0.0

    def test_moderate_ss_blends_components(self):
        """At SS=0.7 (moderate), both strand and density contribute."""
        from hulkrna.locus import compute_gdna_rate_hybrid

        ss = 0.7
        rate, _ = compute_gdna_rate_hybrid(
            100.0, 40.0, 5000.0, ss, 0.001,
        )
        # Should be between pure-strand and pure-density
        assert rate > 0.0


# =====================================================================
# gDNA CLI/Config integration
# =====================================================================


class TestGdnaConfig:
    """Verify gDNA config fields and CLI args."""

    def test_emconfig_defaults(self):
        cfg = EMConfig()
        assert cfg.gdna_kappa_chrom is None
        assert cfg.gdna_kappa_locus is None
        assert cfg.gdna_mom_min_evidence_chrom == 50.0
        assert cfg.gdna_mom_min_evidence_locus == 30.0

    def test_emconfig_explicit(self):
        cfg = EMConfig(
            gdna_kappa_chrom=25.0,
            gdna_kappa_locus=10.0,
            gdna_mom_min_evidence_chrom=40.0,
            gdna_mom_min_evidence_locus=20.0,
        )
        assert cfg.gdna_kappa_chrom == 25.0
        assert cfg.gdna_kappa_locus == 10.0
        assert cfg.gdna_mom_min_evidence_chrom == 40.0
        assert cfg.gdna_mom_min_evidence_locus == 20.0

    def test_cli_parser_has_gdna_args(self):
        from hulkrna.cli import build_parser
        parser = build_parser()
        args = parser.parse_args([
            "quant",
            "--bam", "test.bam",
            "--index", "test_idx",
            "-o", "out",
            "--gdna-kappa-chrom", "30.0",
            "--gdna-kappa-locus", "15.0",
            "--gdna-mom-min-evidence-chrom", "60.0",
            "--gdna-mom-min-evidence-locus", "25.0",
        ])
        assert args.gdna_kappa_chrom == 30.0
        assert args.gdna_kappa_locus == 15.0
        assert args.gdna_mom_min_evidence_chrom == 60.0
        assert args.gdna_mom_min_evidence_locus == 25.0

    def test_cli_defaults_none(self):
        from hulkrna.cli import build_parser
        parser = build_parser()
        args = parser.parse_args([
            "quant",
            "--bam", "test.bam",
            "--index", "test_idx",
            "-o", "out",
        ])
        # Before _resolve_quant_args, CLI defaults are None
        assert args.gdna_kappa_chrom is None
        assert args.gdna_kappa_locus is None
