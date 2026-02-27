"""Tests for hulkrna.bias — BiasProfile and uniform equivalence."""

import math

import numpy as np
import pytest

from hulkrna.bias import (
    BiasProfile,
    build_uniform_profiles,
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
# build_uniform_profiles
# -----------------------------------------------------------------------


class TestBuildUniformProfiles:
    def test_correct_count(self):
        lengths = np.array([100, 200, 500])
        profiles = build_uniform_profiles(lengths)
        assert len(profiles) == 3

    def test_each_is_uniform(self):
        lengths = np.array([100, 200, 500])
        profiles = build_uniform_profiles(lengths)
        for prof, tl in zip(profiles, lengths):
            assert prof.length == tl
            np.testing.assert_array_equal(
                prof.prefix_sum, np.arange(tl + 1, dtype=np.float64)
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
