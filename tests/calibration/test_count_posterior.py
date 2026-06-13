"""Phase 1 — count-module posterior variance (`_count_fraction_variance`, LOESS variance~mean).

Hybrid capture is bimodal (on/off-target chasm), so the imputation variance is a robust LOESS curve in
log-log space rather than a parametric `α·μ²` law. These pin the LOESS smoother (linear recovery,
robustness to outliers), the on/off-target chasm (non-constant relative variance, no gap contamination),
the Poisson count-noise floor, the no-anchor/observable cases, and the Beta cap.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.density_model import (
    _count_fraction_variance, _loess, density_variance_curve,
)


def test_loess_recovers_a_line():
    # Local-linear LOESS is exact on noise-free linear data (interior and at the boundary).
    x = np.linspace(0.0, 10.0, 50)
    y = 2.0 * x + 1.0
    xq = np.array([0.0, 1.0, 5.0, 9.0, 10.0])
    out = _loess(x, y, xq, frac=0.5)
    np.testing.assert_allclose(out, 2.0 * xq + 1.0, atol=1e-6)


def test_loess_is_robust_to_an_outlier():
    # A single wild outlier must not drag the local fit (bisquare robustness downweights it). Realistic
    # noisy data (MAD well-defined) — the no-noise case is degenerate for any MAD-based robustness.
    rng = np.random.default_rng(0)
    x = np.linspace(0.0, 10.0, 51)
    y = x + rng.normal(0.0, 0.3, 51)
    y[25] = 100.0  # outlier at x≈5
    out = _loess(x, y, np.array([5.0]), frac=0.4, robust_iters=2)
    assert abs(out[0] - 5.0) < 1.0  # ≈ the line, not dragged toward 100


def test_nonparametric_curve_captures_the_on_off_target_chasm():
    # Off-target (~1) and on-target (~100) clusters with DIFFERENT relative variance (0.01 vs 0.001).
    # A parametric α·μ² would force one v_rel; LOESS must keep them distinct — and distance weighting
    # must stop the far cluster from contaminating across the gap.
    n = 30
    low = np.linspace(0.8, 1.2, n)
    high = np.linspace(80.0, 120.0, n)
    density = np.concatenate([low, high])
    sl, sh = np.sqrt(0.01), np.sqrt(0.001)  # v_rel = ¼(d_L−d_R)²/μ² = (relative half-spread)²
    dl = np.concatenate([low * (1 + sl), high * (1 + sh)])
    dr = np.concatenate([low * (1 - sl), high * (1 - sh)])
    cg = np.full(2 * n, 0.1)
    var = _count_fraction_variance(
        cg, density, own=np.zeros(2 * n, bool), own_count=np.zeros(2 * n),
        d_left=dl, d_right=dr, n_anchor=np.full(2 * n, 1e6),
    )
    assert var[:n].mean() > var[n:].mean()  # the chasm: low-density nodes relatively noisier
    # distance weighting → no gap contamination, so the cluster MEANS recover the planted v_rel
    np.testing.assert_allclose(var[:n].mean(), 0.1**2 * 0.01, rtol=0.25)
    np.testing.assert_allclose(var[n:].mean(), 0.1**2 * 0.001, rtol=0.25)


def test_poisson_floor_binds_for_sparse_anchors():
    # Near-agreeing anchors ⇒ LOESS curve ≈ 0 ⇒ the Poisson floor 1/N dominates.
    n = 30
    mu = np.linspace(1.0, 5.0, n)
    dl, dr = mu * 1.0001, mu * 0.9999
    cg = np.full(n, 0.1)
    var = _count_fraction_variance(
        cg, 0.5 * (dl + dr), own=np.zeros(n, bool), own_count=np.zeros(n),
        d_left=dl, d_right=dr, n_anchor=np.full(n, 4.0),  # floor 1/4 = 0.25
    )
    np.testing.assert_allclose(var, cg**2 * 0.25, rtol=1e-6)


def test_no_anchor_node_is_uninformative():
    cg = np.array([0.1, 0.4])
    var = _count_fraction_variance(
        cg, np.zeros(2), own=np.array([False, False]), own_count=np.array([0.0, 0.0]),
        d_left=np.array([np.nan, np.nan]), d_right=np.array([np.nan, np.nan]),
        n_anchor=np.array([0.0, 0.0]),
    )
    np.testing.assert_allclose(var, np.minimum(cg**2, cg * (1.0 - cg)), rtol=1e-9)


def test_observable_node_uses_own_count_poisson():
    cg = np.array([0.2])
    var = _count_fraction_variance(
        cg, np.zeros(1), own=np.array([True]), own_count=np.array([50.0]),
        d_left=np.array([np.nan]), d_right=np.array([np.nan]), n_anchor=np.array([0.0]),
    )
    np.testing.assert_allclose(var, cg**2 * (1.0 / 50.0), rtol=1e-9)


def test_variance_capped_at_bernoulli_max():
    cg = np.array([0.9])
    var = _count_fraction_variance(
        cg, np.zeros(1), own=np.array([False]), own_count=np.array([0.0]),
        d_left=np.array([np.nan]), d_right=np.array([np.nan]), n_anchor=np.array([0.0]),
    )
    assert var[0] == pytest.approx(0.9 * 0.1)  # capped, not 0.9² · 1.0


def test_density_variance_curve_two_point_reduces_to_quarter_sqdiff():
    # The shared curve with only the two boundaries must recover ¼(d_L−d_R)² at the fit nodes:
    # plant a clean var∝mean² (raw_var = ¼(d_L−d_R)² = c·μ²) and check the curve returns ≈ c·density².
    n = 40
    mu = np.linspace(1.0, 50.0, n)
    c = 0.02
    half = np.sqrt(c) * mu  # so ¼(d_L−d_R)² = (half)² = c·μ²
    dl, dr = mu + half, mu - half
    ok = np.ones(n, bool)
    s2 = density_variance_curve(mu, d_left=dl, d_right=dr, left_ok=ok, right_ok=ok)
    np.testing.assert_allclose(s2, c * mu**2, rtol=0.2)


def test_density_variance_curve_triplet_uses_three_points():
    # Adding the count-observable contained observation makes each fit node a 3-point sample. A node
    # whose contained sits between agreeing boundaries has LOWER disagreement than the pair alone would
    # imply only if contained agrees; here contained disagrees, so the triplet variance is well-defined
    # and finite (the curve fits on 3-point samples, not just pairs).
    n = 40
    mu = np.linspace(1.0, 50.0, n)
    dl, dr = mu * 1.1, mu * 0.9
    contained = mu * 1.0  # third observation at the mean
    ok = np.ones(n, bool)
    s2_pair = density_variance_curve(mu, d_left=dl, d_right=dr, left_ok=ok, right_ok=ok)
    s2_trip = density_variance_curve(
        mu, d_left=dl, d_right=dr, left_ok=ok, right_ok=ok, contained=contained, contained_ok=ok
    )
    assert np.all(np.isfinite(s2_trip))
    # the triplet's variance-of-the-mean (k=3) is smaller than the pair's (k=2) when the extra point
    # sits at the mean (more observations → tighter mean estimate).
    assert np.nanmean(s2_trip) < np.nanmean(s2_pair)


def test_too_few_two_anchor_falls_back_to_floor():
    # Fewer than _LOESS_MIN_FIT 2-anchor points ⇒ LOESS skipped ⇒ imputed nodes use the Poisson floor.
    n = 3
    var = _count_fraction_variance(
        np.full(n, 0.1), np.full(n, 2.0), own=np.zeros(n, bool), own_count=np.zeros(n),
        d_left=np.full(n, np.nan), d_right=np.full(n, 5.0), n_anchor=np.full(n, 10.0),
    )
    np.testing.assert_allclose(var, 0.1**2 * (1.0 / 10.0), rtol=1e-9)
