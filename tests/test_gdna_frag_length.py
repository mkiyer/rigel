"""Tests for gDNA/RNA fragment length distribution framework.

Covers:
- FragmentLengthModel.from_counts() factory for pre-built distributions
- Strand-weighted histogram mixing math (mix_models)
- Proof that correct distributions are recovered from synthetic data
- Graceful degradation when strand info is absent
- Integration with scoring (gdna_model vs rna_model likelihood separation)
- Serialization of new model fields
"""

import math

import numpy as np
import pytest

from rigel.frag_length_model import (
    FragmentLengthModel,
    FragmentLengthModels,
    DEFAULT_MAX_FRAG_SIZE,
)
from rigel.splice import SpliceType


# =====================================================================
# Helpers
# =====================================================================

def _build_models(
    *,
    spliced_lengths=None,
    intergenic_lengths=None,
    same_strand_lengths=None,
    opp_strand_lengths=None,
    max_size=500,
):
    """Build FragmentLengthModels with given observations."""
    m = FragmentLengthModels(max_size=max_size)
    if spliced_lengths is not None:
        for l in spliced_lengths:
            m.observe(l, splice_type=SpliceType.SPLICED_ANNOT)
    if intergenic_lengths is not None:
        for l in intergenic_lengths:
            m.observe(l, splice_type=None)  # intergenic path
    if same_strand_lengths is not None:
        m.unspliced_same_strand.observe_batch(
            np.asarray(same_strand_lengths, dtype=np.intp)
        )
    if opp_strand_lengths is not None:
        m.unspliced_opp_strand.observe_batch(
            np.asarray(opp_strand_lengths, dtype=np.intp)
        )
    return m


def _normal_counts(mean, std, n, max_size=500):
    """Discretized normal distribution as a histogram array."""
    rng = np.random.default_rng(42)
    samples = rng.normal(mean, std, size=n).astype(int)
    samples = np.clip(samples, 0, max_size)
    counts = np.bincount(samples, minlength=max_size + 1).astype(np.float64)
    return counts


# =====================================================================
# FragmentLengthModel.from_counts() factory
# =====================================================================

class TestFromCounts:
    """Tests for the from_counts() factory method."""

    def test_basic_creation(self):
        """Create a model from a simple histogram."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        counts[250] = 50.0
        model = FragmentLengthModel.from_counts(counts)

        assert model._finalized
        assert model.max_size == 500
        assert model.total_weight == pytest.approx(150.0)
        assert model.n_observations == 150
        assert model._log_prob is not None

    def test_max_size_inferred(self):
        """max_size is inferred as len(counts) - 1."""
        counts = np.ones(101)
        model = FragmentLengthModel.from_counts(counts)
        assert model.max_size == 100

    def test_explicit_max_size(self):
        """Explicit max_size with shorter counts array zero-fills."""
        counts = np.array([0, 0, 10, 20, 10])
        model = FragmentLengthModel.from_counts(counts, max_size=1000)
        assert model.max_size == 1000
        assert model.total_weight == pytest.approx(40.0)
        assert model.counts[2] == pytest.approx(10.0)
        assert model.counts[999] == pytest.approx(0.0)

    def test_log_likelihood_matches_trained(self):
        """from_counts() produces same log_likelihood as observe+finalize."""
        # Train a model via observations
        trained = FragmentLengthModel(max_size=500)
        for _ in range(100):
            trained.observe(200)
        for _ in range(50):
            trained.observe(300)
        trained.finalize()

        # Create same model via from_counts
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        counts[300] = 50.0
        factory = FragmentLengthModel.from_counts(counts)

        # Log-likelihoods must be identical
        for length in [0, 100, 200, 250, 300, 500]:
            assert factory.log_likelihood(length) == pytest.approx(
                trained.log_likelihood(length)
            )

    def test_tail_decay_works(self):
        """Queries beyond max_size use exponential tail decay."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        model = FragmentLengthModel.from_counts(counts)

        # Beyond max_size: each extra bp adds ~log(0.99)
        ll_500 = model.log_likelihood(500)
        ll_600 = model.log_likelihood(600)
        assert ll_600 < ll_500
        assert ll_600 == pytest.approx(
            model._tail_base + (600 - 500) * math.log(0.99)
        )

    def test_empty_counts(self):
        """Zero-count histogram yields uniform distribution."""
        counts = np.zeros(101, dtype=np.float64)
        model = FragmentLengthModel.from_counts(counts)
        assert model._finalized
        # All bins should have equal probability (uniform)
        ll_0 = model.log_likelihood(0)
        ll_50 = model.log_likelihood(50)
        assert ll_0 == pytest.approx(ll_50)
        assert ll_0 == pytest.approx(-np.log(101))

    def test_single_peak(self):
        """Single-peak histogram: mode has highest likelihood."""
        counts = np.zeros(501, dtype=np.float64)
        counts[250] = 1000.0
        model = FragmentLengthModel.from_counts(counts)
        assert model.log_likelihood(250) > model.log_likelihood(100)
        assert model.log_likelihood(250) > model.log_likelihood(400)
        assert model.mode == 250

    def test_statistics(self):
        """Mean, std, median, mode are computed correctly."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        counts[300] = 100.0
        model = FragmentLengthModel.from_counts(counts)
        assert model.mean == pytest.approx(250.0)
        # Median: cumulative 50% is at bin 200 (first peak has exactly 50%)
        assert model.median == pytest.approx(200.0)
        assert model.mode in (200, 300)  # either peak


# =====================================================================
# Mix model math correctness
# =====================================================================

class TestMixModelsBasic:
    """Core mixing math correctness."""

    def test_perfect_stranding(self):
        """s_rna=1.0: all antisense → gDNA, all sense → RNA."""
        m = _build_models(
            spliced_lengths=[200] * 100,
            intergenic_lengths=[300] * 50,
            same_strand_lengths=[250] * 80,
            opp_strand_lengths=[350] * 40,
            max_size=500,
        )
        m.mix_models(s_rna=1.0, p_r1_sense=0.8)

        S, A = 80.0, 40.0
        R = (S - A) / 1.0  # 40
        G = (S + A) - R    # 80
        W_sense = G * 0.5 / S  # 0.5
        W_anti = G * 0.5 / A   # 1.0
        # Directional certainty: gDNA uses max(2W-1, 0), RNA uses max(1-2W, 0)
        Cg_sense = max(2 * W_sense - 1, 0)  # 0.0
        Cg_anti = max(2 * W_anti - 1, 0)    # 1.0
        Cr_sense = max(1 - 2 * W_sense, 0)  # 0.0
        Cr_anti = max(1 - 2 * W_anti, 0)    # 0.0

        assert m.gdna_model.total_weight == pytest.approx(
            50 + W_sense * Cg_sense * S + W_anti * Cg_anti * A  # 90
        )
        assert m.rna_model.total_weight == pytest.approx(
            100 + (1 - W_sense) * Cr_sense * S + (1 - W_anti) * Cr_anti * A  # 100
        )

    def test_no_strand_info(self):
        """s_rna=0.5: no separation possible → gDNA=intergenic, RNA=spliced+all unspliced."""
        m = _build_models(
            spliced_lengths=[200] * 100,
            intergenic_lengths=[300] * 50,
            same_strand_lengths=[250] * 80,
            opp_strand_lengths=[350] * 40,
        )
        m.mix_models(s_rna=0.5, p_r1_sense=0.5)

        assert m.gdna_model.total_weight == pytest.approx(50.0)
        assert m.rna_model.total_weight == pytest.approx(220.0)

    def test_typical_stranding_075(self):
        """s_rna=0.75: intermediate weights are computed correctly."""
        m = _build_models(
            spliced_lengths=[200] * 100,
            intergenic_lengths=[300] * 50,
            same_strand_lengths=[250] * 80,
            opp_strand_lengths=[350] * 40,
        )
        m.mix_models(s_rna=0.75, p_r1_sense=0.75)

        S, A = 80.0, 40.0
        denom = 0.5
        R = (S - A) / denom  # 80
        G = (S + A) - R      # 40
        W_sense = G * 0.5 / S  # 0.25
        W_anti = G * 0.5 / A   # 0.5
        Cg_sense = max(2 * W_sense - 1, 0)  # 0.0
        Cg_anti = max(2 * W_anti - 1, 0)    # 0.0
        Cr_sense = max(1 - 2 * W_sense, 0)  # 0.5
        Cr_anti = max(1 - 2 * W_anti, 0)    # 0.0

        assert m.gdna_model.total_weight == pytest.approx(
            50 + W_sense * Cg_sense * S + W_anti * Cg_anti * A  # 50
        )
        assert m.rna_model.total_weight == pytest.approx(
            100 + (1 - W_sense) * Cr_sense * S + (1 - W_anti) * Cr_anti * A  # 130
        )

    def test_certainty_penalty_discards_uncertain(self):
        """Certainty penalty means gDNA + RNA <= total input.

        Uncertain fragments (W ≈ 0.5) are excluded from both models,
        so the sum of output weights is less than or equal to the input.
        """
        for s_rna in [0.55, 0.7, 0.85, 1.0]:
            m2 = _build_models(
                spliced_lengths=[200] * 130,
                intergenic_lengths=[300] * 70,
                same_strand_lengths=[250] * 90,
                opp_strand_lengths=[350] * 45,
            )
            m2.mix_models(s_rna=s_rna, p_r1_sense=0.8)
            total_in = 130 + 70 + 90 + 45
            total_out = m2.gdna_model.total_weight + m2.rna_model.total_weight
            # Anchors (intergenic + spliced) always included
            assert total_out >= 130 + 70, f"s_rna={s_rna}"
            # Total cannot exceed input
            assert total_out <= total_in + 0.01, f"s_rna={s_rna}"


class TestMixModelsStrandSwap:
    """Verify p_r1_sense < 0.5 swaps same/opp → anti/sense."""

    def test_r1_antisense_protocol(self):
        """p_r1_sense < 0.5: same_strand → antisense pool."""
        m = _build_models(
            spliced_lengths=[200] * 100,
            intergenic_lengths=[300] * 50,
            same_strand_lengths=[250] * 30,
            opp_strand_lengths=[350] * 90,
        )
        m.mix_models(s_rna=1.0, p_r1_sense=0.3)

        # sense = opp (90), anti = same (30)
        S, A = 90.0, 30.0
        R = (S - A) / 1.0  # 60
        G = 120 - R         # 60
        W_sense = G * 0.5 / S  # 1/3
        W_anti = G * 0.5 / A   # 1.0
        Cg_sense = max(2 * W_sense - 1, 0)  # 0
        Cg_anti = max(2 * W_anti - 1, 0)    # 1.0
        Cr_sense = max(1 - 2 * W_sense, 0)  # 1/3
        Cr_anti = max(1 - 2 * W_anti, 0)    # 0

        assert m.gdna_model.total_weight == pytest.approx(
            50 + W_sense * Cg_sense * S + W_anti * Cg_anti * A
        )
        assert m.rna_model.total_weight == pytest.approx(
            100 + (1 - W_sense) * Cr_sense * S + (1 - W_anti) * Cr_anti * A
        )

    def test_symmetry_of_swap(self):
        """Swapping same/opp counts AND flipping p_r1_sense gives same result."""
        m1 = _build_models(
            spliced_lengths=[200] * 50,
            intergenic_lengths=[300] * 20,
            same_strand_lengths=[250] * 60,
            opp_strand_lengths=[350] * 30,
        )
        m1.mix_models(s_rna=0.9, p_r1_sense=0.8)

        m2 = _build_models(
            spliced_lengths=[200] * 50,
            intergenic_lengths=[300] * 20,
            same_strand_lengths=[350] * 30,  # swapped
            opp_strand_lengths=[250] * 60,   # swapped
        )
        m2.mix_models(s_rna=0.9, p_r1_sense=0.2)  # flipped

        assert m1.gdna_model.total_weight == pytest.approx(
            m2.gdna_model.total_weight
        )
        assert m1.rna_model.total_weight == pytest.approx(
            m2.rna_model.total_weight
        )


class TestMixModelsEdgeCases:
    """Edge cases: empty data, single-strand, etc."""

    def test_empty_unspliced(self):
        """No unspliced observations: gDNA = intergenic, RNA = spliced."""
        m = _build_models(
            spliced_lengths=[200] * 100,
            intergenic_lengths=[300] * 50,
        )
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)
        assert m.gdna_model.total_weight == pytest.approx(50.0)
        assert m.rna_model.total_weight == pytest.approx(100.0)

    def test_only_same_strand(self):
        """Only same-strand observations: opp=0 → R depends on S/denom."""
        m = _build_models(
            spliced_lengths=[200] * 50,
            intergenic_lengths=[300] * 20,
            same_strand_lengths=[250] * 60,
        )
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)

        S, A = 60.0, 0.0
        denom = 0.8
        R = max(0.0, (S - A) / denom)  # 75
        G = (S + A) - R  # -15 → W_sense clamped to 0
        W_sense = min(max((G * 0.5 / S) if S > 0 else 0.0, 0.0), 1.0)
        Cg = max(2 * W_sense - 1, 0)
        Cr = max(1 - 2 * W_sense, 0)

        assert m.gdna_model.total_weight == pytest.approx(
            20 + W_sense * Cg * S
        )
        assert m.rna_model.total_weight == pytest.approx(
            50 + (1 - W_sense) * Cr * S
        )

    def test_only_opp_strand(self):
        """Only opp-strand observations: S=0 → R=0 → W_anti=0.5 → certainty=0."""
        m = _build_models(
            spliced_lengths=[200] * 50,
            intergenic_lengths=[300] * 20,
            opp_strand_lengths=[350] * 40,
        )
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)

        # S=0, A=40 → R=0, G=40, W_anti=0.5, C_anti=0.0
        # Entire antisense pool excluded by certainty penalty
        assert m.gdna_model.total_weight == pytest.approx(20.0)
        assert m.rna_model.total_weight == pytest.approx(50.0)

    def test_no_data_at_all(self):
        """Completely empty models don't crash."""
        m = FragmentLengthModels(max_size=100)
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)
        assert m.gdna_model.total_weight == pytest.approx(0.0)
        assert m.rna_model.total_weight == pytest.approx(0.0)

    def test_empty_finalize_gives_uniform(self):
        """Empty models finalize to uniform distributions."""
        m = FragmentLengthModels(max_size=100)
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)
        m.finalize()

        # Both should give uniform (since zero observations)
        assert m.gdna_model.log_likelihood(50) == pytest.approx(-np.log(101))
        assert m.rna_model.log_likelihood(50) == pytest.approx(-np.log(101))


# =====================================================================
# Proving correct distributions are recovered from data
# =====================================================================

class TestDistributionRecovery:
    """Verify that mix_models recovers the true gDNA and RNA distributions
    from synthetic observations with known ground truth."""

    def _simulate_observations(
        self, rng, n_rna, n_gdna, rna_mean, gdna_mean,
        std, s_rna, max_size=500,
    ):
        """Generate synthetic unspliced strand observations.

        Simulates a genic locus with known RNA and gDNA fragment pools:
        - RNA fragments: normal(rna_mean, std), oriented by s_rna
        - gDNA fragments: normal(gdna_mean, std), 50/50 strand orientation
        """
        # RNA fragment lengths
        rna_lengths = np.clip(
            rng.normal(rna_mean, std, size=n_rna).astype(int), 1, max_size
        )
        # gDNA fragment lengths
        gdna_lengths = np.clip(
            rng.normal(gdna_mean, std, size=n_gdna).astype(int), 1, max_size
        )

        # Assign strand orientation
        # RNA: s_rna fraction go to sense, (1-s_rna) to antisense
        rna_sense_mask = rng.random(n_rna) < s_rna
        # gDNA: 50/50
        gdna_sense_mask = rng.random(n_gdna) < 0.5

        sense_pool = np.concatenate([
            rna_lengths[rna_sense_mask],
            gdna_lengths[gdna_sense_mask],
        ])
        anti_pool = np.concatenate([
            rna_lengths[~rna_sense_mask],
            gdna_lengths[~gdna_sense_mask],
        ])

        return sense_pool, anti_pool, rna_lengths, gdna_lengths

    def test_recover_distinct_peaks(self):
        """With well-separated peaks (RNA=200, gDNA=400), mixing recovers both."""
        rng = np.random.default_rng(12345)
        n_rna, n_gdna = 2000, 1000
        rna_mean, gdna_mean, std = 200, 400, 20
        s_rna = 0.9
        max_size = 500

        sense, anti, true_rna, true_gdna = self._simulate_observations(
            rng, n_rna, n_gdna, rna_mean, gdna_mean, std, s_rna, max_size,
        )

        # Also generate intergenic (pure gDNA) and spliced (pure RNA)
        intergenic = np.clip(
            rng.normal(gdna_mean, std, size=500).astype(int), 1, max_size
        )
        spliced = np.clip(
            rng.normal(rna_mean, std, size=500).astype(int), 1, max_size
        )

        m = _build_models(
            spliced_lengths=spliced.tolist(),
            intergenic_lengths=intergenic.tolist(),
            same_strand_lengths=sense.tolist(),
            opp_strand_lengths=anti.tolist(),
            max_size=max_size,
        )
        m.mix_models(s_rna=s_rna, p_r1_sense=s_rna)
        m.finalize()

        # RNA model should peak near rna_mean
        rna_mode = m.rna_model.mode
        assert abs(rna_mode - rna_mean) < 3 * std

        # gDNA model should peak near gdna_mean
        gdna_mode = m.gdna_model.mode
        assert abs(gdna_mode - gdna_mean) < 3 * std

        # At the RNA peak, RNA model gives higher probability
        assert m.rna_model.log_likelihood(rna_mean) > (
            m.gdna_model.log_likelihood(rna_mean)
        )
        # At the gDNA peak, gDNA model gives higher probability
        assert m.gdna_model.log_likelihood(gdna_mean) > (
            m.rna_model.log_likelihood(gdna_mean)
        )

    def test_gdna_mean_recovery(self):
        """Certainty penalty yields a gDNA model mean close to truth."""
        rng = np.random.default_rng(99)
        n_rna, n_gdna = 3000, 1000
        s_rna = 0.85
        max_size = 500

        sense, anti, _, _ = self._simulate_observations(
            rng, n_rna, n_gdna, 200, 350, 30, s_rna, max_size,
        )
        intergenic = np.clip(
            rng.normal(350, 30, size=200).astype(int), 1, max_size
        )
        spliced = np.clip(
            rng.normal(200, 30, size=800).astype(int), 1, max_size
        )

        m = _build_models(
            spliced_lengths=spliced.tolist(),
            intergenic_lengths=intergenic.tolist(),
            same_strand_lengths=sense.tolist(),
            opp_strand_lengths=anti.tolist(),
            max_size=max_size,
        )
        m.mix_models(s_rna=s_rna, p_r1_sense=s_rna)

        # gDNA mean should be close to the true value (350)
        gdna_mean = m.gdna_model.mean
        assert abs(gdna_mean - 350) < 30, f"gDNA mean {gdna_mean:.1f} too far from 350"

        # RNA mean should be close to the true value (200)
        rna_mean = m.rna_model.mean
        assert abs(rna_mean - 200) < 20, f"RNA mean {rna_mean:.1f} too far from 200"

        # Anchors always included: min weight >= intergenic + spliced
        assert m.gdna_model.total_weight >= 200
        assert m.rna_model.total_weight >= 800

    def test_recovery_degrades_at_low_strand_specificity(self):
        """At s_rna=0.55, separation is poor but doesn't crash."""
        rng = np.random.default_rng(77)
        n_rna, n_gdna = 2000, 1000
        s_rna = 0.55  # near-unstranded
        max_size = 500

        sense, anti, _, _ = self._simulate_observations(
            rng, n_rna, n_gdna, 200, 400, 25, s_rna, max_size,
        )

        m = _build_models(
            spliced_lengths=[200] * 500,
            intergenic_lengths=[400] * 200,
            same_strand_lengths=sense.tolist(),
            opp_strand_lengths=anti.tolist(),
            max_size=max_size,
        )
        m.mix_models(s_rna=s_rna, p_r1_sense=s_rna)
        m.finalize()

        # At very low SS, most unspliced weight goes to RNA (W_sense, W_anti ≈ 0)
        # gDNA model dominated by intergenic
        assert m.gdna_model.total_weight > 0
        assert m.rna_model.total_weight > 0
        # Should still be able to query likelihoods
        assert np.isfinite(m.gdna_model.log_likelihood(300))
        assert np.isfinite(m.rna_model.log_likelihood(300))


# =====================================================================
# from_counts for test scenario construction
# =====================================================================

class TestFromCountsForScenarios:
    """Demonstrate using from_counts() to build test scenarios."""

    def test_custom_gdna_model_in_container(self):
        """Install a from_counts() model as the gdna_model on a container."""
        # Build a container with default training
        m = FragmentLengthModels(max_size=500)
        for _ in range(200):
            m.observe(200, SpliceType.SPLICED_ANNOT)
            m.observe(200, None)  # intergenic

        # Replace gdna_model with a custom distribution
        gdna_counts = np.zeros(501, dtype=np.float64)
        gdna_counts[350] = 500.0
        gdna_counts[400] = 300.0
        m.gdna_model = FragmentLengthModel.from_counts(gdna_counts)

        # Replace rna_model similarly
        rna_counts = np.zeros(501, dtype=np.float64)
        rna_counts[200] = 800.0
        rna_counts[250] = 200.0
        m.rna_model = FragmentLengthModel.from_counts(rna_counts)

        m.finalize()

        # gdna_model peaks at 350, rna_model peaks at 200
        assert m.gdna_model.log_likelihood(350) > m.gdna_model.log_likelihood(200)
        assert m.rna_model.log_likelihood(200) > m.rna_model.log_likelihood(350)

    def test_from_counts_matches_normal_training(self):
        """Manually constructed histogram via from_counts matches
        the same histogram built through observe() calls."""
        # Build via observe()
        trained = FragmentLengthModel(max_size=500)
        for _ in range(100):
            trained.observe(180)
        for _ in range(200):
            trained.observe(220)
        for _ in range(50):
            trained.observe(260)
        trained.finalize()

        # Build via from_counts
        counts = np.zeros(501, dtype=np.float64)
        counts[180] = 100
        counts[220] = 200
        counts[260] = 50
        factory = FragmentLengthModel.from_counts(counts)

        # Every bin must match
        for length in range(501):
            assert factory.log_likelihood(length) == pytest.approx(
                trained.log_likelihood(length)
            ), f"Mismatch at length={length}"


# =====================================================================
# Finalize + scoring integration
# =====================================================================

class TestMixModelsFinalizeScoring:
    """Integration: mix → finalize → scoring produces consistent results."""

    def test_finalize_after_mix_creates_luts(self):
        """After mix+finalize, both models have valid _log_prob arrays."""
        m = _build_models(
            spliced_lengths=[200] * 200,
            intergenic_lengths=[300] * 100,
            same_strand_lengths=[250] * 80,
            opp_strand_lengths=[350] * 40,
        )
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)
        m.finalize()

        assert m.rna_model._log_prob is not None
        assert m.gdna_model._log_prob is not None
        assert m.rna_model._finalized
        assert m.gdna_model._finalized

    def test_rna_gdna_discriminate(self):
        """RNA and gDNA models assign higher probability to their own peak."""
        m = _build_models(
            spliced_lengths=[150] * 300,
            intergenic_lengths=[400] * 200,
            same_strand_lengths=[200] * 100,
            opp_strand_lengths=[350] * 50,
        )
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)
        m.finalize()

        # At RNA peak (150): RNA model should win
        assert m.rna_model.log_likelihood(150) > m.gdna_model.log_likelihood(150)
        # At gDNA peak (400): gDNA model should win
        assert m.gdna_model.log_likelihood(400) > m.rna_model.log_likelihood(400)

    def test_log_likelihood_is_proper_distribution(self):
        """Exponentiated log-likelihoods approximately sum to 1."""
        m = _build_models(
            spliced_lengths=[200] * 500,
            intergenic_lengths=[350] * 300,
            same_strand_lengths=[220] * 200,
            opp_strand_lengths=[340] * 100,
            max_size=500,
        )
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)
        m.finalize()

        for model in [m.rna_model, m.gdna_model]:
            probs = np.exp([model.log_likelihood(k) for k in range(501)])
            # With Laplace smoothing, not exactly 1.0, but close
            assert abs(probs.sum() - 1.0) < 0.01


# =====================================================================
# Serialization
# =====================================================================

class TestToDict:
    """Verify serialization includes new models."""

    def test_to_dict_has_all_keys(self):
        m = _build_models(
            spliced_lengths=[200] * 50,
            intergenic_lengths=[300] * 20,
            same_strand_lengths=[250] * 30,
            opp_strand_lengths=[350] * 10,
        )
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)
        m.finalize()

        d = m.to_dict()
        assert "rna" in d
        assert "gdna" in d
        assert "unspliced_same_strand" in d
        assert "unspliced_opp_strand" in d
        assert "global" in d
        assert "intergenic" in d

    def test_to_dict_summary_accurate(self):
        """Serialized summary stats match live model."""
        m = _build_models(
            spliced_lengths=[200] * 100,
            intergenic_lengths=[350] * 50,
            same_strand_lengths=[200] * 60,
            opp_strand_lengths=[300] * 30,
        )
        m.mix_models(s_rna=0.9, p_r1_sense=0.9)
        m.finalize()

        d = m.to_dict()
        rna_summary = d["rna"]["summary"]
        assert rna_summary["total_weight"] == pytest.approx(
            m.rna_model.total_weight, rel=0.01
        )
        gdna_summary = d["gdna"]["summary"]
        assert gdna_summary["total_weight"] == pytest.approx(
            m.gdna_model.total_weight, rel=0.01
        )
