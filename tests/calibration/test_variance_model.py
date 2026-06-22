"""The message-reliability var~mean fitter (``MonotoneVarMean.fit_offset``).

Pins the Poisson-offset decomposition (``imputation_variance_model.md``): the fit learns the BIOLOGICAL
dispersion σ²_bio (the EXCESS over the computed Poisson sampling floor V_p), not σ²_bio + sampling. Also
pins monotonicity / non-negativity / finiteness, flat extrapolation outside the fit range, the
too-few-points flat fallback, and the diagnostics dataframe.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.variance_model import MonotoneMean, MonotoneVarMean


def _poisson_pair_points(lam_src, lam_dst, l_src, l_dst, rng):
    """Vectorized node-pair samples → (μ=ρ_dst, R²=(ρ_dst−ρ_src)², V_p=ρ_src/L_src+ρ_dst/L_dst).

    Each density ρ = C/L with C ~ Poisson(λ·L), so Var_samp(ρ) = C/L² = ρ/L (the computed offset)."""
    c_src = rng.poisson(np.maximum(np.asarray(lam_src) * l_src, 0.0))
    c_dst = rng.poisson(np.maximum(np.asarray(lam_dst) * l_dst, 0.0))
    rho_src, rho_dst = c_src / l_src, c_dst / l_dst
    r2 = (rho_dst - rho_src) ** 2
    v_p = rho_src / l_src + rho_dst / l_dst
    return rho_dst, r2, v_p


def _rising_offset_points(rng, n=600):
    """A rising σ²_bio = 0.05·μ² over a Poisson floor V_p — generic monotone-fit test data."""
    mu = np.exp(rng.uniform(np.log(0.2), np.log(20.0), n))
    v_p = mu * rng.uniform(0.05, 0.3, n)
    r2 = np.maximum(0.05 * mu**2 + v_p + rng.normal(0.0, 0.02, n), 0.0)
    return mu, r2, v_p


def test_offset_fit_uniform_field_recovers_zero_dispersion():
    # Uniform field (λ=5 everywhere ⇒ true σ²_bio=0). The squared residual is PURE Poisson sampling noise.
    # The Poisson-offset fit must learn σ²_bio ≈ 0 — far below the raw total (σ²_bio + sampling), which is
    # the contamination the decomposition removes.
    rng = np.random.default_rng(0)
    n, lam = 800, 5.0
    l_src = rng.uniform(2.0, 200.0, n)  # varied lengths ⇒ varied certainty at the SAME density
    l_dst = rng.uniform(2.0, 200.0, n)
    mu, r2, v_p = _poisson_pair_points(lam, lam, l_src, l_dst, rng)
    sigma_bio = float(MonotoneVarMean.fit_offset(mu, r2, v_p).predict(np.array([lam]))[0])
    assert sigma_bio < 0.25 * float(np.mean(v_p))  # ≈ 0, far below the sampling floor
    assert sigma_bio < 0.25 * float(np.mean(r2))  # and far below the raw total a naive fit would read


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


def test_offset_fit_predict_clips_to_fit_range():
    # Flat extrapolation outside the fitted log-mean range (not a runaway).
    mu, r2, v_p = _rising_offset_points(np.random.default_rng(4))
    fit = MonotoneVarMean.fit_offset(mu, r2, v_p)
    lo_in = fit.predict(np.array([float(np.exp(fit.x_lo))]))[0]
    lo_out = fit.predict(np.array([1e-9]))[0]  # far below the fit range
    hi_in = fit.predict(np.array([float(np.exp(fit.x_hi))]))[0]
    hi_out = fit.predict(np.array([1e9]))[0]  # far above
    assert np.isclose(lo_out, lo_in, rtol=1e-6)
    assert np.isclose(hi_out, hi_in, rtol=1e-6)


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


def test_offset_fit_to_dataframe_has_points_and_curve():
    df = MonotoneVarMean.fit_offset(*_rising_offset_points(np.random.default_rng(6))).to_dataframe()
    assert set(df["kind"].unique()) == {"point", "curve"}
    assert (df["var"] >= 0.0).all()  # σ²_bio ≥ 0


# --- MonotoneMean: the enrichment transfer ê(z)=E[ρ_g|z] (mean fit, count-recalibrated) -------------------

def _power_law_points(rng, n=500, C=3.0, b=1.08, sigma=0.5):
    """z ~ loguniform, ρ_g = C·z^b·lognormal(σ) (the validated near-log-linear enrichment law)."""
    z = np.exp(rng.uniform(np.log(0.05), np.log(7.0), n))
    rho = C * z**b * np.exp(rng.normal(0.0, sigma, n))
    return z, rho


def test_mean_fit_is_monotone_nonneg_finite():
    z, rho = _power_law_points(np.random.default_rng(10))
    fit = MonotoneMean.fit(z, rho)
    grid = np.logspace(-2, 1, 200)
    pred = fit.predict(grid)
    assert np.all(np.isfinite(pred))
    assert np.all(pred >= 0.0)
    assert np.all(np.diff(pred) >= -1e-9)  # ê monotone non-decreasing in z


def test_mean_fit_net_unbiased_on_recal_weight():
    # The count-recal scale forces Σ ê(z)·rw = Σ ρ·rw on the fit set, EXACTLY (the §F net fix).
    rng = np.random.default_rng(11)
    z, rho = _power_law_points(rng)
    rw = rng.uniform(50.0, 400.0, z.size)  # E_dst-like recal weight
    fit = MonotoneMean.fit(z, rho, weight=rw, recal_weight=rw)
    assert np.sum(rw * fit.predict(z)) == pytest.approx(np.sum(rw * rho), rel=1e-9)


def test_mean_fit_recovers_power_law_conditional_mean():
    # ê at a query z must ≈ the conditional MEAN C·z^b·exp(σ²/2) (the recal scale supplies the lognormal-mean
    # correction that a raw log-space fit's exp would miss — the Duan role, as one global scale).
    rng = np.random.default_rng(12)
    C, b, sigma = 3.0, 1.08, 0.5
    z, rho = _power_law_points(rng, n=4000, C=C, b=b, sigma=sigma)
    fit = MonotoneMean.fit(z, rho)
    for zq in (0.5, 2.0, 5.0):
        truth = C * zq**b * np.exp(sigma**2 / 2)
        assert fit.predict(np.array([zq]))[0] == pytest.approx(truth, rel=0.25)


def test_mean_fit_flat_input_gives_constant():
    # ρ_g independent of z ⇒ ê ≈ a constant (the recal-weighted mean), no spurious slope.
    rng = np.random.default_rng(13)
    z = np.exp(rng.uniform(np.log(0.05), np.log(7.0), 500))
    rho = 4.0 * np.exp(rng.normal(0.0, 0.4, z.size))  # no z-dependence
    fit = MonotoneMean.fit(z, rho)
    pred = fit.predict(np.logspace(-2, 1, 50))
    assert pred.max() / pred.min() < 1.6  # nearly flat
    assert float(np.mean(pred)) == pytest.approx(np.mean(rho), rel=0.25)


def test_mean_fit_predict_clips_to_fit_range():
    z, rho = _power_law_points(np.random.default_rng(14))
    fit = MonotoneMean.fit(z, rho)
    lo_in = fit.predict(np.array([float(np.exp(fit.x_lo))]))[0]
    hi_in = fit.predict(np.array([float(np.exp(fit.x_hi))]))[0]
    assert np.isclose(fit.predict(np.array([1e-9]))[0], lo_in, rtol=1e-6)  # flat below
    assert np.isclose(fit.predict(np.array([1e9]))[0], hi_in, rtol=1e-6)   # flat above


def test_mean_fit_too_few_points_flat_fallback():
    z = np.array([1.0, 2.0, 3.0])
    rho = np.array([5.0, 6.0, 5.5])
    fit = MonotoneMean.fit(z, rho)
    pred = fit.predict(np.array([0.5, 2.0, 10.0]))
    assert np.allclose(pred, pred[0])  # flat
    assert pred[0] == pytest.approx(np.mean(rho), rel=0.05)


def test_mean_constant_factory_is_flat_at_value():
    fit = MonotoneMean.constant(0.46)  # the off-capture / z-unusable ρ_global fallback
    pred = fit.predict(np.logspace(-3, 3, 40))
    assert np.allclose(pred, 0.46)
