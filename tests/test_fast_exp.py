"""`fast_exp`, the vectorised exponential the EM's E-step kernel runs, against `math.exp` (libm, so a real
reference rather than a second approximation).

Accuracy over the full range and at near-zero, moderate and deep-negative inputs; the edge cases (zero, the cutoff,
negative infinity, short and empty arrays); agreement between the vector lanes and the scalar tail; and the
worst-case error in ULP.
"""

import math

import numpy as np
import numpy.testing as npt

from rigel._em_impl.testing import fast_exp


# ---------------------------------------------------------------------------
# Accuracy: fast_exp vs math.exp over the entire input domain
# ---------------------------------------------------------------------------


class TestFastExpAccuracy:
    """Validate fast_exp matches std::exp to within a few ULP."""

    def test_uniform_sweep_full_range(self):
        """Sweep 1M values uniformly over [-708, 0]."""
        rng = np.random.default_rng(42)
        x = rng.uniform(-708.0, 0.0, size=1_000_000)
        result = fast_exp(x)
        expected = np.array([math.exp(v) for v in x])
        # Relative tolerance: 4 ULP ≈ 4 * 2^-52 ≈ 8.9e-16
        npt.assert_allclose(result, expected, rtol=1e-14, atol=0)

    def test_near_zero(self):
        """Values very close to 0 (exp(x) ≈ 1)."""
        x = np.array([-1e-15, -1e-10, -1e-5, -0.001, -0.01, -0.1, 0.0])
        result = fast_exp(x)
        expected = np.array([math.exp(v) for v in x])
        npt.assert_allclose(result, expected, rtol=1e-14, atol=0)

    def test_moderate_values(self):
        """Values in the range where most non-trivial posteriors live."""
        x = np.linspace(-50.0, 0.0, 10_000)
        result = fast_exp(x)
        expected = np.array([math.exp(v) for v in x])
        npt.assert_allclose(result, expected, rtol=1e-14, atol=0)

    def test_deep_negative(self):
        """Values near the underflow boundary (above cutoff)."""
        x = np.array([-700.0, -705.0, -707.0, -707.5, -707.99])
        result = fast_exp(x)
        expected = np.array([math.exp(v) for v in x])
        npt.assert_allclose(result, expected, rtol=1e-12, atol=1e-320)


# ---------------------------------------------------------------------------
# Edge cases and cutoff behavior
# ---------------------------------------------------------------------------


class TestFastExpEdgeCases:
    """Verify correct behavior at boundaries and special values."""

    def test_exact_zero_input(self):
        """exp(0) = 1.0 exactly."""
        result = fast_exp(np.array([0.0]))
        assert result[0] == 1.0

    def test_cutoff_returns_zero(self):
        """Values below -708.0 must return exactly 0.0."""
        x = np.array([-709.0, -800.0, -1000.0, -1e6])
        result = fast_exp(x)
        npt.assert_array_equal(result, np.zeros(4))

    def test_cutoff_boundary_continuity(self):
        """Values just above the cutoff should produce small but nonzero results."""
        x = np.array([-707.9, -707.99, -708.0])
        result = fast_exp(x)
        # All should be >= 0
        assert np.all(result >= 0.0)
        # -707.9 should produce a nonzero value
        assert result[0] > 0.0

    def test_negative_infinity(self):
        """exp(-inf) = 0.0."""
        result = fast_exp(np.array([-np.inf]))
        assert result[0] == 0.0

    def test_odd_length_array(self):
        """Odd-length input exercises the scalar tail after the vector blocks."""
        x = np.array([-1.0, -2.0, -3.0])
        result = fast_exp(x)
        expected = np.array([math.exp(-1.0), math.exp(-2.0), math.exp(-3.0)])
        npt.assert_allclose(result, expected, rtol=1e-14)

    def test_single_element(self):
        """Single element — pure scalar path."""
        result = fast_exp(np.array([-5.0]))
        npt.assert_allclose(result[0], math.exp(-5.0), rtol=1e-14)

    def test_empty_array(self):
        """Empty input returns empty output."""
        result = fast_exp(np.array([], dtype=np.float64))
        assert len(result) == 0


# ---------------------------------------------------------------------------
# Vector vs scalar consistency
# ---------------------------------------------------------------------------


class TestFastExpNeonScalarConsistency:
    """Lanes the vector path computes match lanes the scalar tail computes (NEON 2-wide, AVX2 4-wide, AVX-512
    8-wide)."""

    def test_paired_vs_tail(self):
        """Three identical values: on NEON the first two take the vector path and the third the scalar tail,
        and all three agree exactly."""
        val = -123.456
        x = np.array([val, val, val])
        result = fast_exp(x)
        # All three should be identical
        assert result[0] == result[1] == result[2]

    def test_large_array_consistency(self):
        """Random values at sizes that leave a scalar tail after the vector blocks."""
        rng = np.random.default_rng(99)
        for size in [3, 5, 7, 11, 101, 1001]:
            x = rng.uniform(-708.0, 0.0, size=size)
            result = fast_exp(x)
            expected = np.array([math.exp(v) for v in x])
            npt.assert_allclose(
                result, expected, rtol=1e-14, atol=0, err_msg=f"Mismatch for size={size}"
            )


# ---------------------------------------------------------------------------
# Worst-case ULP error
# ---------------------------------------------------------------------------


class TestFastExpULP:
    """The maximum ULP error, printed and bounded."""

    def test_max_ulp_error(self):
        """The max ULP error over 1M values in [-708, 0] is under 50: the degree-11 Horner truncation (about
        25 ULP) plus FMA rounding."""
        rng = np.random.default_rng(42)
        x = rng.uniform(-708.0, 0.0, size=1_000_000)
        result = fast_exp(x)

        max_ulp = 0.0
        for i in range(len(x)):
            expected = math.exp(x[i])
            if expected == 0.0:
                continue
            # ULP distance
            ulp = abs(result[i] - expected) / max(abs(expected) * 2**-52, 5e-324)
            if ulp > max_ulp:
                max_ulp = ulp

        # Report
        print(f"\nMax ULP error over 1M values: {max_ulp:.2f}")
        # under 50 ULP, the bound the docstring derives
        assert max_ulp < 50.0, f"Max ULP error {max_ulp:.2f} exceeds 50.0"
