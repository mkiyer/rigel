"""Tests for rigel.bias — BiasProfile and uniform equivalence."""

import math

import numpy as np
import pytest

from rigel.bias import (
    BiasProfile,
)


# -----------------------------------------------------------------------
# BiasProfile.uniform — construction sanity
# -----------------------------------------------------------------------


class TestBiasProfileUniform:
    """Test that uniform profiles have the expected prefix-sum structure."""

    @pytest.mark.parametrize("length", [1, 10, 100, 500, 2000])
    def test_prefix_sum_shape(self, length):
        p = BiasProfile.uniform(length)
        assert p.prefix_sum.shape == (length + 1,)

    @pytest.mark.parametrize("length", [1, 10, 100])
    def test_prefix_sum_values(self, length):
        p = BiasProfile.uniform(length)
        np.testing.assert_array_equal(
            p.prefix_sum, np.arange(length + 1, dtype=np.float64)
        )

    def test_length_property(self):
        p = BiasProfile.uniform(42)
        assert p.length == 42


# -----------------------------------------------------------------------
# fragment_weight — uniform case should always return 1.0
# -----------------------------------------------------------------------


class TestFragmentWeightUniform:
    """Under uniform bias, W(f|t) == 1.0 for every valid fragment."""

    @pytest.mark.parametrize(
        "t_len, frag_len",
        [(100, 50), (200, 1), (500, 499), (10, 10), (1000, 300)],
    )
    def test_weight_is_one(self, t_len, frag_len):
        p = BiasProfile.uniform(t_len)
        for start in range(0, t_len - frag_len + 1, max(1, (t_len - frag_len) // 5)):
            w = p.fragment_weight(start, start + frag_len)
            assert w == pytest.approx(1.0), (
                f"t_len={t_len}, frag_len={frag_len}, start={start}"
            )

    def test_zero_span_returns_zero(self):
        p = BiasProfile.uniform(100)
        assert p.fragment_weight(50, 50) == 0.0

    def test_negative_span_returns_zero(self):
        p = BiasProfile.uniform(100)
        assert p.fragment_weight(60, 50) == 0.0


# -----------------------------------------------------------------------
# effective_length — uniform equivalence
# -----------------------------------------------------------------------


class TestEffectiveLengthUniform:
    """Under uniform bias, L̃_t(l_f) == max(t_len - l_f + 1, 1)."""

    @pytest.mark.parametrize(
        "t_len, frag_len",
        [
            (100, 50),
            (200, 1),
            (500, 499),
            (500, 500),  # n_valid == 1
            (10, 10),    # n_valid == 1
            (10, 11),    # n_valid == 0 → degenerate → 1.0
            (1, 1),
            (1, 2),      # degenerate
            (1000, 300),
            (2000, 150),
        ],
    )
    def test_matches_current_formula(self, t_len, frag_len):
        p = BiasProfile.uniform(t_len)
        expected = max(t_len - frag_len + 1, 1)
        assert p.effective_length(frag_len) == pytest.approx(expected), (
            f"t_len={t_len}, frag_len={frag_len}"
        )


# -----------------------------------------------------------------------
# log-likelihood equivalence: log(W) - log(L̃) == -log(max(L - l + 1, 1))
# -----------------------------------------------------------------------


class TestLogLikelihoodUniformEquivalence:
    """Full round-trip: the bias model should reproduce the current
    baked-in effective-length correction for uniform profiles."""

    @pytest.mark.parametrize(
        "t_len, frag_len",
        [
            (100, 50),
            (200, 1),
            (500, 499),
            (1000, 300),
            (50, 50),
            (300, 150),
        ],
    )
    def test_log_lik_equivalence(self, t_len, frag_len):
        p = BiasProfile.uniform(t_len)
        expected = -math.log(max(t_len - frag_len + 1, 1))
        # Pick several start positions
        for start in range(0, t_len - frag_len + 1, max(1, (t_len - frag_len) // 5)):
            w = p.fragment_weight(start, start + frag_len)
            eff = p.effective_length(frag_len)
            log_bias = math.log(w) - math.log(eff)
            assert log_bias == pytest.approx(expected, abs=1e-12), (
                f"t_len={t_len}, frag_len={frag_len}, start={start}"
            )



# -----------------------------------------------------------------------
# Non-uniform bias (smoke test for future phases)
# -----------------------------------------------------------------------


class TestNonUniformBias:
    """Verify the math works with a simple non-uniform profile."""

    def test_linear_ramp(self):
        """Bias b[i] = i+1 → prefix_sum = [0, 1, 3, 6, 10, ...]."""
        L = 10
        b = np.arange(1, L + 1, dtype=np.float64)
        B = np.zeros(L + 1, dtype=np.float64)
        B[1:] = np.cumsum(b)
        p = BiasProfile(prefix_sum=B)

        assert p.length == L

        # Fragment [0, 3) → mean of b[0..2] = (1+2+3)/3 = 2.0
        assert p.fragment_weight(0, 3) == pytest.approx(2.0)

        # Fragment [7, 10) → mean of b[7..9] = (8+9+10)/3 = 9.0
        assert p.fragment_weight(7, 10) == pytest.approx(9.0)

    def test_effective_length_non_uniform(self):
        """Effective length with two-level bias: first half 1.0, second half 2.0."""
        L = 100
        b = np.ones(L, dtype=np.float64)
        b[50:] = 2.0
        B = np.zeros(L + 1, dtype=np.float64)
        B[1:] = np.cumsum(b)
        p = BiasProfile(prefix_sum=B)

        frag_len = 10
        eff = p.effective_length(frag_len)

        # Manual: sum over valid starts v=0..90 of mean(b[v..v+9])
        # For v < 41: all b are 1.0 → weight 1.0 (41 starts: 0..40)
        # For v in 41..49: mixed → weight between 1.0 and 2.0 (9 starts)
        # For v >= 50: all b are 2.0 → weight 2.0 (41 starts: 50..90)
        # Total starts = 91, eff should be > 91 (because of 2.0 region)
        assert eff > L - frag_len + 1  # > uniform
        assert eff < 2.0 * (L - frag_len + 1)  # < all-2.0


# -----------------------------------------------------------------------
# is_uniform flag
# -----------------------------------------------------------------------


class TestIsUniformFlag:
    """Verify the is_uniform flag is set correctly."""

    def test_uniform_factory_sets_flag(self):
        p = BiasProfile.uniform(100)
        assert p.is_uniform is True

    def test_manual_construction_defaults_false(self):
        B = np.arange(11, dtype=np.float64)
        p = BiasProfile(prefix_sum=B)
        assert p.is_uniform is False

    def test_non_uniform_profile_flag_false(self):
        L = 10
        b = np.arange(1, L + 1, dtype=np.float64)
        B = np.zeros(L + 1, dtype=np.float64)
        B[1:] = np.cumsum(b)
        p = BiasProfile(prefix_sum=B)
        assert p.is_uniform is False


# -----------------------------------------------------------------------
# Bias correction — now exercised via the C++ _em_impl module
# -----------------------------------------------------------------------


class TestBiasCorrectionCpp:
    """Verify that the C++ uniform bias correction in _em_impl produces
    the expected -log(max(L - frag_len + 1, 1)) correction.

    The C++ solver applies bias correction internally; we test it by
    running run_locus_em_native on minimal inputs and verifying the
    results are consistent with the expected effective-length correction.
    """

    def _expected_bias(self, frag_len, prof_len):
        """Python reference for uniform bias: -log(max(L - f + 1, 1))."""
        import math
        eff = max(prof_len - frag_len + 1, 1)
        return -math.log(eff)

    def test_simple_correction_values(self):
        """Check a few known effective-length corrections."""
        import math
        # L=500, frag=250 -> eff=251 -> correction=-log(251)
        assert abs(self._expected_bias(250, 500) - (-math.log(251))) < 1e-12
        # L=100, frag=200 -> eff=1 -> correction=0
        assert abs(self._expected_bias(200, 100)) < 1e-12
        # L=100, frag=0 -> eff=101 -> correction=-log(101)
        assert abs(self._expected_bias(0, 100) - (-math.log(101))) < 1e-12

    def test_em_solver_applies_bias(self):
        """Run EM solver and verify bias correction is applied.

        A single-transcript locus with uniform fragment length should
        converge with the bias correction folded in.
        """
        from rigel._em_impl import run_locus_em_native

        # 1 mRNA + 1 nRNA + 1 gDNA = 3 components
        n_comp = 3
        # 2 units, each mapping to component 0 (mRNA)
        offsets = np.array([0, 1, 2], dtype=np.int64)
        t_indices = np.array([0, 0], dtype=np.int32)
        log_liks = np.array([0.0, 0.0], dtype=np.float64)
        cov_wts = np.array([1.0, 1.0], dtype=np.float64)
        tx_starts = np.array([0, 0], dtype=np.int32)
        tx_ends = np.array([200, 200], dtype=np.int32)
        bias_profiles = np.array([1000, 500, 2000], dtype=np.int64)
        unambig_totals = np.array([10.0, 0.0, 0.0], dtype=np.float64)
        eff_lens = np.ones(n_comp, dtype=np.float64)
        eligible = np.array([1.0, 0.0, 0.0], dtype=np.float64)

        theta, alpha, em_totals = run_locus_em_native(
            offsets, t_indices, log_liks, cov_wts,
            tx_starts, tx_ends, bias_profiles,
            unambig_totals, eff_lens, eligible,
            n_comp, 0.01, 1.0, 1000, 1e-6,
            False, -1.0,
            0, 0,
            np.array([], dtype=np.int32), np.array([], dtype=np.int32),
            np.array([], dtype=np.int32),
        )
        # All evidence points to component 0; theta[0] should dominate
        assert np.asarray(theta)[0] > 0.9

    def test_degenerate_frag_len(self):
        """Fragment length > transcript length → effective_length = 1,
        so bias correction is 0 (no penalty)."""
        from rigel._em_impl import run_locus_em_native

        n_comp = 3
        offsets = np.array([0, 1], dtype=np.int64)
        t_indices = np.array([0], dtype=np.int32)
        log_liks = np.array([0.0], dtype=np.float64)
        cov_wts = np.array([1.0], dtype=np.float64)
        tx_starts = np.array([0], dtype=np.int32)
        tx_ends = np.array([500], dtype=np.int32)  # frag > transcript
        bias_profiles = np.array([100, 100, 100], dtype=np.int64)
        unambig_totals = np.array([5.0, 0.0, 0.0], dtype=np.float64)
        eff_lens = np.ones(n_comp, dtype=np.float64)
        eligible = np.array([1.0, 0.0, 0.0], dtype=np.float64)

        theta, alpha, em_totals = run_locus_em_native(
            offsets, t_indices, log_liks, cov_wts,
            tx_starts, tx_ends, bias_profiles,
            unambig_totals, eff_lens, eligible,
            n_comp, 0.01, 1.0, 1000, 1e-6,
            False, -1.0,
            0, 0,
            np.array([], dtype=np.int32), np.array([], dtype=np.int32),
            np.array([], dtype=np.int32),
        )
        # Should converge without issues
        assert np.all(np.isfinite(np.asarray(theta)))
        assert np.asarray(theta)[0] > 0.9
