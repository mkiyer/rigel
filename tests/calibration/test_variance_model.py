"""Monotone-increasing P-spline (SCAM) var~mean fitter + the node-PAIR IMPUTATION builder.

Pins: monotonicity by construction, power-law recovery, robustness to an outlier, flat extrapolation
outside the fit range, the too-few-points fallback, the diagnostics dataframe, and the node-PAIR
imputation reliability (CALIBRATION_PLAN_v5 §3 — the Step-2 imputation substrate). The DIRECT own-count
builder and the RNA-density assembly were removed in the Step-1 precision rebuild (count→fraction-precision).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.variance_model import (
    MonotoneVarMean,
    fit_pair_imputation_varmean,
)


def _powerlaw(n=400, a=0.01, b=2.0, seed=0):
    rng = np.random.default_rng(seed)
    mean = np.exp(rng.uniform(np.log(0.01), np.log(10.0), n))
    var = a * mean**b * np.exp(rng.normal(0.0, 0.3, n))  # log-normal noise around the law
    return mean, var


def test_fit_is_monotone_on_a_grid():
    mean, var = _powerlaw()
    fit = MonotoneVarMean.fit(mean, var)
    grid = np.logspace(-2, 1, 200)
    pred = fit.predict(grid)
    assert np.all(np.diff(pred) >= -1e-9)  # non-decreasing everywhere


def test_recovers_power_law_trend():
    mean, var = _powerlaw(a=0.02, b=2.0, n=600)
    fit = MonotoneVarMean.fit(mean, var)
    xs = np.array([0.05, 0.5, 5.0])
    pred = fit.predict(xs)
    truth = 0.02 * xs**2.0
    # within a factor of ~2 of the planted trend across 2 decades
    assert np.all(pred > truth / 2.5) and np.all(pred < truth * 2.5)


def test_robust_to_a_wild_outlier():
    mean, var = _powerlaw(n=400, seed=1)
    var = var.copy()
    var[mean.argmax()] = 1e6  # a wild high-variance outlier at the top mean
    robust = MonotoneVarMean.fit(mean, var, robust_iters=2)
    naive = MonotoneVarMean.fit(mean, var, robust_iters=0)
    top = float(mean.max())
    # robust fit at the top is far less dragged toward the 1e6 outlier than the naive fit
    assert robust.predict(np.array([top]))[0] < 0.5 * naive.predict(np.array([top]))[0]


def test_predict_clips_to_fit_range():
    mean, var = _powerlaw()
    fit = MonotoneVarMean.fit(mean, var)
    lo_in = fit.predict(np.array([float(np.exp(fit.x_lo))]))[0]
    lo_out = fit.predict(np.array([1e-9]))[0]  # far below the fit range
    hi_in = fit.predict(np.array([float(np.exp(fit.x_hi))]))[0]
    hi_out = fit.predict(np.array([1e9]))[0]  # far above
    assert np.isclose(lo_out, lo_in, rtol=1e-6)  # flat extrapolation, not a runaway
    assert np.isclose(hi_out, hi_in, rtol=1e-6)


def test_too_few_points_falls_back_monotone():
    mean = np.array([0.1, 1.0, 5.0])
    var = np.array([1e-3, 0.1, 2.0])
    fit = MonotoneVarMean.fit(mean, var)  # < k points ⇒ power-law fallback
    pred = fit.predict(np.logspace(-1, 1, 50))
    assert np.all(np.diff(pred) >= -1e-9)


def test_to_dataframe_has_points_and_curve():
    mean, var = _powerlaw(n=200)
    fit = MonotoneVarMean.fit(mean, var)
    df = fit.to_dataframe()
    assert set(df["kind"].unique()) == {"point", "curve"}
    assert (df["var"] > 0).all()


def test_pair_imputation_builder_densifies_and_is_monotone():
    # The node-PAIR builder: one point per (observable side → eligible region) adjacency, at
    # mean = region density, raw_var = (d_region − d_side)². A both-flanks region → 2 points; a
    # one-flank region → 1 point (the densification the both-sides triplet missed).
    rng = np.random.default_rng(7)
    n = 300
    rd = np.exp(rng.uniform(np.log(0.05), np.log(8.0), n))
    # planted: each side density is the region density plus log-normal multiplicative error.
    ld = rd * np.exp(rng.normal(0.0, 0.4, n))
    rrd = rd * np.exp(rng.normal(0.0, 0.4, n))
    ref = np.zeros(n, dtype=np.int64)
    elig = rng.random(n) < 0.7
    lok = rng.random(n) < 0.8
    rok = rng.random(n) < 0.8
    fit = fit_pair_imputation_varmean(
        rd, ld, rrd, region_eligible=elig, left_ok=lok, right_ok=rok, ref_id=ref
    )
    # point count == number of (eligible region, observable flank) adjacencies (both-flanks → 2)
    expected = int((elig & lok).sum() + (elig & rok).sum())
    assert fit.fit_mean.size == expected
    # monotone over the fit range
    g = np.logspace(np.log10(np.exp(fit.x_lo)), np.log10(np.exp(fit.x_hi)), 100)
    assert np.all(np.diff(fit.predict(g)) >= -1e-9)
    # the queried axis is the REGION density (means are drawn from rd[eligible & flank-ok])
    assert fit.fit_mean.min() >= rd[elig].min() * (1 - 1e-9)


def test_pair_imputation_single_flank_contributes():
    # A region eligible with only ONE observable flank still contributes exactly one point.
    rd = np.array([1.0, 2.0, 3.0])
    ld = np.array([0.5, 0.0, 0.0])
    rrd = np.array([0.0, 0.0, 0.0])
    elig = np.array([True, True, True])
    lok = np.array([True, False, False])  # only region 0 has a left flank
    rok = np.array([False, False, False])
    ref = np.zeros(3, dtype=np.int64)
    fit = fit_pair_imputation_varmean(
        rd, ld, rrd, region_eligible=elig, left_ok=lok, right_ok=rok, ref_id=ref
    )
    assert fit.fit_mean.size == 1
    assert np.isclose(fit.fit_mean[0], 1.0)  # mean = region density of region 0


def test_jensen_offset_values():
    # The Jensen df-offset Δ = log(dof/2) − ψ(dof/2): positive, decreasing in dof, → 0 as dof → ∞.
    from rigel.calibration.variance_model import _jensen_offset

    off = _jensen_offset(np.array([1.0, 2.0, 1e6]))
    assert off[0] == pytest.approx(1.2703628454614782, rel=1e-9)  # dof=1 (k=2 disagreement)
    assert off[1] == pytest.approx(0.5772156649015329, rel=1e-9)  # dof=2 (Euler–Mascheroni)
    assert off[2] == pytest.approx(0.0, abs=1e-5)  # dof → ∞ (no correction)
    assert off[0] > off[1] > off[2] >= 0.0


def test_jensen_offset_inflates_recovered_variance():
    # With dof passed, the fit target log(var) is shifted UP by Δ_k>0 ⇒ the recovered variance is
    # uniformly larger than the un-corrected fit (removes the small-dof over-confidence).
    mean, var = _powerlaw(n=400, seed=3)
    dof = np.full(mean.shape[0], 1.0)  # every point a k=2 disagreement
    corrected = MonotoneVarMean.fit(mean, var, dof=dof)
    plain = MonotoneVarMean.fit(mean, var)  # dof=None ⇒ back-compat, no offset
    g = np.logspace(-1.5, 0.5, 40)
    pc, pp = corrected.predict(g), plain.predict(g)
    assert np.all(pc > pp)  # inflated everywhere
    # the inflation ≈ exp(Δ_1) = exp(1.2704) ≈ 3.56× (the verified k=2 factor)
    ratio = np.median(pc / pp)
    assert ratio == pytest.approx(np.exp(1.2703628454614782), rel=0.15)


# --- the Poisson-offset fit (imputation_variance_model.md): learn σ²_bio, not σ²_bio + sampling ---


def _poisson_pair_points(lam_src, lam_dst, l_src, l_dst, rng):
    """Vectorized node-pair samples → (μ=ρ_dst, R²=(ρ_dst−ρ_src)², V_p=ρ_src/L_src+ρ_dst/L_dst).

    Each density ρ = C/L with C ~ Poisson(λ·L), so Var_samp(ρ) = C/L² = ρ/L (the computed offset)."""
    c_src = rng.poisson(np.maximum(np.asarray(lam_src) * l_src, 0.0))
    c_dst = rng.poisson(np.maximum(np.asarray(lam_dst) * l_dst, 0.0))
    rho_src, rho_dst = c_src / l_src, c_dst / l_dst
    r2 = (rho_dst - rho_src) ** 2
    v_p = rho_src / l_src + rho_dst / l_dst
    return rho_dst, r2, v_p


def test_offset_fit_uniform_field_recovers_zero_dispersion():
    # Uniform field (λ=5 everywhere ⇒ true σ²_bio=0). The squared residual is PURE Poisson sampling noise.
    # The Poisson-offset fit must learn σ²_bio ≈ 0 — whereas the naive log-log TOTAL fit reads ~the sampling
    # floor and over-states the biological dispersion. This is the contamination the decomposition removes.
    rng = np.random.default_rng(0)
    n, lam = 800, 5.0
    l_src = rng.uniform(2.0, 200.0, n)  # varied lengths ⇒ varied certainty at the SAME density
    l_dst = rng.uniform(2.0, 200.0, n)
    mu, r2, v_p = _poisson_pair_points(lam, lam, l_src, l_dst, rng)
    offset_fit = MonotoneVarMean.fit_offset(mu, r2, v_p)
    naive_fit = MonotoneVarMean.fit(mu, r2, dof=np.ones(n))  # the contaminated density-only total
    sigma_bio = float(offset_fit.predict(np.array([lam]))[0])
    naive = float(naive_fit.predict(np.array([lam]))[0])
    # σ²_bio ≈ 0, and far below both the sampling floor and the naive total (the de-confounding).
    assert sigma_bio < 0.25 * float(np.mean(v_p))
    assert sigma_bio < 0.25 * naive


def test_offset_fit_recovers_planted_dispersion():
    # Plant a real σ²_bio (adjacent true rates differ) on top of Poisson sampling; the offset fit recovers
    # σ²_bio (≈ the planted variance), not σ²_bio + sampling.
    rng = np.random.default_rng(1)
    n, lam0, sd_bio = 1000, 20.0, 4.0  # σ²_bio = 16
    l_src = rng.uniform(50.0, 200.0, n)
    l_dst = rng.uniform(50.0, 200.0, n)
    lam_dst = np.maximum(lam0 + rng.normal(0.0, sd_bio, n), 0.1)  # dst true rate differs by σ_bio
    mu, r2, v_p = _poisson_pair_points(np.full(n, lam0), lam_dst, l_src, l_dst, rng)
    fit = MonotoneVarMean.fit_offset(mu, r2, v_p)
    sigma_bio = float(fit.predict(np.array([lam0]))[0])
    assert 0.4 * sd_bio**2 < sigma_bio < 2.5 * sd_bio**2  # ≈ 16, within ~2×


def test_offset_fit_is_nonneg_monotone_finite():
    rng = np.random.default_rng(2)
    n = 500
    mu = np.exp(rng.uniform(np.log(0.1), np.log(50.0), n))
    v_p = mu * rng.uniform(0.05, 0.3, n)
    r2 = np.maximum(0.05 * mu**2 + v_p + rng.normal(0.0, 0.02, n), 0.0)  # σ²_bio = 0.05·μ² (rising)
    fit = MonotoneVarMean.fit_offset(mu, r2, v_p)
    grid = np.logspace(-1, np.log10(50.0), 120)
    pred = fit.predict(grid)
    assert np.all(np.isfinite(pred))
    assert np.all(pred >= 0.0)  # σ²_bio is a variance, floored at 0
    assert np.all(np.diff(pred) >= -1e-9)  # monotone non-decreasing


def test_offset_fit_too_few_points_flat_fallback():
    # < k points ⇒ a FLAT σ²_bio = the reliability-weighted mean excess max(mean(R²−V_p), 0).
    mu = np.array([1.0, 2.0, 3.0])
    r2 = np.array([0.5, 0.6, 0.55])
    v_p = np.array([0.4, 0.4, 0.4])
    fit = MonotoneVarMean.fit_offset(mu, r2, v_p)
    pred = fit.predict(np.array([0.5, 1.5, 3.0, 10.0]))
    assert np.all(np.isfinite(pred)) and np.all(pred >= 0.0)
    assert np.allclose(pred, pred[0])  # flat
    assert pred[0] == pytest.approx(0.15, abs=0.05)  # mean([0.1, 0.2, 0.15])

