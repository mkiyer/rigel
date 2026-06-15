"""Monotone-increasing P-spline (SCAM) var~mean fitter + the DIRECT/IMPUTATION split.

Pins: monotonicity by construction, power-law recovery, robustness to an outlier, flat extrapolation
outside the fit range, the too-few-points fallback, the diagnostics dataframe, and that the builders
partition the points by region-observability.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.variance_model import (
    MonotoneVarMean,
    VarMeanPoints,
    fit_direct_varmean,
    fit_imputation_varmean,
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


def test_builders_partition_by_region_observability():
    rng = np.random.default_rng(2)
    n = 500
    mean = np.exp(rng.uniform(np.log(0.01), np.log(10.0), n))
    var = 0.01 * mean**2 * np.exp(rng.normal(0.0, 0.3, n))
    obs = rng.random(n) < 0.5
    pts = VarMeanPoints(mean=mean, raw_var=var, region_observable=obs)
    direct = fit_direct_varmean(pts)
    imp = fit_imputation_varmean(pts)
    # each builder fit only its own subset
    assert direct.fit_mean.size == int(obs.sum())
    assert imp.fit_mean.size == int((~obs).sum())
    # both monotone
    g = np.logspace(-2, 1, 100)
    assert np.all(np.diff(direct.predict(g)) >= -1e-9)
    assert np.all(np.diff(imp.predict(g)) >= -1e-9)
