"""gDNA FL effective lengths: boundary = μ_FL, region = E_f[max(0, L−ℓ)]."""

from __future__ import annotations

import numpy as np

from rigel.calibration.effective_length import (
    boundary_eff_length,
    fl_mean,
    region_eff_length,
)


def _brute_region(L: float, pmf: np.ndarray) -> float:
    """Σ_{ℓ≤L} (L−ℓ) f(ℓ), the definition region_eff_length implements."""
    p = pmf / pmf.sum()
    return float(sum((L - ell) * p[ell] for ell in range(p.shape[0]) if ell <= L))


def _spike(at: int, n: int = 801) -> np.ndarray:
    p = np.zeros(n, dtype=np.float64)
    p[at] = 1.0
    return p


def test_fl_mean_normalizes_unnormalized_counts():
    counts = np.zeros(801)
    counts[100] = 3.0  # raw counts, not a pmf
    counts[300] = 1.0
    # mean = (3*100 + 1*300) / 4 = 150
    assert fl_mean(counts) == 150.0


def test_boundary_eff_length_is_mu_fl_independent_of_regions():
    # THE regression against the old "tent" error: a boundary's gDNA exposure is
    # the FL mean, NOT min(L_left, L_right, …). With μ_FL = 300, the boundary
    # exposure is 300 regardless of the 400 bp / 200 bp neighbours.
    pmf = _spike(300)
    assert boundary_eff_length(pmf) == 300.0  # not capped at min(400, 200) = 200


def test_region_eff_length_spike():
    pmf = _spike(300)
    # L=400: only fragments ≤400 fit; the 300 bp spike → 400−300 = 100.
    # L=200: a 300 bp fragment cannot be contained in 200 bp → 0.
    np.testing.assert_allclose(
        region_eff_length(np.array([400.0, 200.0, 300.0]), pmf), [100.0, 0.0, 0.0]
    )


def test_region_eff_length_matches_brute_force():
    rng_pmf = np.zeros(801, dtype=np.float64)
    rng_pmf[50:801] = 1.0  # uniform over [50, 800]; mean = 425
    for L in (0.0, 100.0, 300.0, 425.0, 800.0, 1000.0):
        np.testing.assert_allclose(
            region_eff_length(np.array([L]), rng_pmf)[0], _brute_region(L, rng_pmf), rtol=1e-12
        )


def test_region_eff_length_large_L_is_L_minus_mu():
    pmf = np.zeros(801, dtype=np.float64)
    pmf[200:401] = 1.0  # mean = 300
    mu = fl_mean(pmf)
    # For L well beyond the support, contained exposure → L − μ_FL.
    np.testing.assert_allclose(region_eff_length(np.array([5000.0]), pmf)[0], 5000.0 - mu)


def test_region_eff_length_zero_and_nonnegative():
    pmf = _spike(150)
    eff = region_eff_length(np.array([0.0, 50.0, 150.0]), pmf)
    assert np.all(eff >= 0.0)
    assert eff[0] == 0.0  # zero-length region: no contained exposure
    assert eff[1] == 0.0  # 50 bp region cannot contain a 150 bp fragment
