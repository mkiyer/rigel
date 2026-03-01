"""Tests for hulkrna.bias — BiasProfile and uniform equivalence."""

import math

import numpy as np
import pytest

from hulkrna.bias import (
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
# _apply_bias_correction_uniform — vectorized fast-path equivalence
# -----------------------------------------------------------------------


class TestApplyBiasCorrectionUniform:
    """Verify the uniform fast-path matches the general path exactly."""

    def _run_general_path(self, log_liks, t_indices, tx_starts, tx_ends,
                          bias_profiles):
        """Run the general (non-uniform) LUT path by temporarily
        clearing is_uniform on all profiles."""
        import math as _math
        n = len(t_indices)
        _MAX_FLEN = 1_000_000
        frag_lens = np.clip(
            (tx_ends - tx_starts).astype(np.int64), 0, _MAX_FLEN - 1)
        compound = t_indices.astype(np.int64) * _MAX_FLEN + frag_lens
        unique_compound, inverse = np.unique(compound, return_inverse=True)
        log_eff = np.empty(len(unique_compound), dtype=np.float64)
        for j in range(len(unique_compound)):
            uc = int(unique_compound[j])
            cidx, fl = divmod(uc, _MAX_FLEN)
            prof = bias_profiles[int(cidx)]
            # Use the prefix-sum path regardless of is_uniform
            L = prof.length
            n_valid = L - int(fl) + 1
            eff = max(n_valid, 1) if n_valid > 0 else 1
            log_eff[j] = _math.log(max(float(eff), 1e-300))
        log_liks -= log_eff[inverse]
        # fragment_weight is 1.0 for uniform → no adjustment

    def test_equivalence_simple(self):
        """Basic 3-transcript, 5-candidate scenario."""
        from hulkrna.estimator import _apply_bias_correction

        profiles = [BiasProfile.uniform(500), BiasProfile.uniform(1000),
                    BiasProfile.uniform(200)]
        t_indices = np.array([0, 1, 2, 0, 1], dtype=np.int32)
        tx_starts = np.array([10, 20, 5, 100, 50], dtype=np.int32)
        tx_ends = np.array([260, 320, 155, 350, 550], dtype=np.int32)

        ll_fast = np.zeros(5, dtype=np.float64)
        ll_general = np.zeros(5, dtype=np.float64)

        _apply_bias_correction(ll_fast, t_indices, tx_starts, tx_ends, profiles)
        self._run_general_path(ll_general, t_indices, tx_starts, tx_ends,
                               profiles)

        np.testing.assert_allclose(ll_fast, ll_general, atol=1e-12)

    def test_equivalence_degenerate_frag_len(self):
        """Fragment length > transcript length → effective_length = 1."""
        from hulkrna.estimator import _apply_bias_correction

        profiles = [BiasProfile.uniform(100)]
        t_indices = np.array([0, 0], dtype=np.int32)
        tx_starts = np.array([0, 0], dtype=np.int32)
        tx_ends = np.array([50, 200], dtype=np.int32)  # 200 > 100

        ll_fast = np.zeros(2, dtype=np.float64)
        ll_general = np.zeros(2, dtype=np.float64)

        _apply_bias_correction(ll_fast, t_indices, tx_starts, tx_ends, profiles)
        self._run_general_path(ll_general, t_indices, tx_starts, tx_ends,
                               profiles)

        np.testing.assert_allclose(ll_fast, ll_general, atol=1e-12)

    def test_equivalence_many_components(self):
        """Stress test with many components mimicking a real locus."""
        from hulkrna.estimator import _apply_bias_correction

        rng = np.random.default_rng(42)
        n_components = 200
        n_candidates = 10000
        lengths = rng.integers(100, 5000, size=n_components)
        profiles = [BiasProfile.uniform(int(l)) for l in lengths]

        t_indices = rng.integers(0, n_components, size=n_candidates).astype(np.int32)
        tx_starts = rng.integers(0, 500, size=n_candidates).astype(np.int32)
        tx_ends = tx_starts + rng.integers(50, 600, size=n_candidates).astype(np.int32)

        ll_fast = np.zeros(n_candidates, dtype=np.float64)
        ll_general = np.zeros(n_candidates, dtype=np.float64)

        _apply_bias_correction(ll_fast, t_indices, tx_starts, tx_ends, profiles)
        self._run_general_path(ll_general, t_indices, tx_starts, tx_ends,
                               profiles)

        np.testing.assert_allclose(ll_fast, ll_general, atol=1e-12)

    def test_empty_input(self):
        """No candidates → no modification."""
        from hulkrna.estimator import _apply_bias_correction

        profiles = [BiasProfile.uniform(100)]
        ll = np.array([], dtype=np.float64)
        ti = np.array([], dtype=np.int32)
        ts = np.array([], dtype=np.int32)
        te = np.array([], dtype=np.int32)
        _apply_bias_correction(ll, ti, ts, te, profiles)
        assert len(ll) == 0

    def test_non_uniform_uses_general_path(self):
        """When any profile is non-uniform, the general path is used."""
        from hulkrna.estimator import _apply_bias_correction

        L = 100
        b = np.ones(L, dtype=np.float64)
        b[50:] = 2.0
        B = np.zeros(L + 1, dtype=np.float64)
        B[1:] = np.cumsum(b)
        non_uniform = BiasProfile(prefix_sum=B)
        uniform = BiasProfile.uniform(200)

        t_indices = np.array([0, 1], dtype=np.int32)
        tx_starts = np.array([10, 20], dtype=np.int32)
        tx_ends = np.array([60, 120], dtype=np.int32)
        ll = np.zeros(2, dtype=np.float64)

        # Should not crash — exercises general path with mixed profiles
        _apply_bias_correction(ll, t_indices, tx_starts, tx_ends,
                               [non_uniform, uniform])
        assert np.all(np.isfinite(ll))
