"""M-step hyperparameter updates (doc 03 §5): closed forms + bounded 1-D fits."""

from __future__ import annotations

import numpy as np

from rigel.calibration.estep import _PI_CLIP
from rigel.calibration.mstep import (
    _RHO_D_BB_FALLBACK,
    fit_rho_d_bb,
    update_exposure_dispersion,
    update_pi_g_prior,
    update_rho_0,
)


def test_update_rho_0_closed_form():
    # Σ M_g_tot / Σ(ω·L) = (10+20) / (1·100 + 2·100) = 30/300 = 0.1.
    rho_0 = update_rho_0(np.array([10.0, 20.0]), np.array([1.0, 2.0]), np.array([100.0, 100.0]))
    np.testing.assert_allclose(rho_0, 0.1)


def test_update_rho_0_empty_exposure_floor():
    rho_0 = update_rho_0(np.array([0.0]), np.array([0.0]), np.array([100.0]))
    assert rho_0 > 0.0  # 1 / max(ΣL, 1)


def test_exposure_dispersion_low_variance_is_small():
    # Tightly concentrated exposure posteriors (ω≈1, α large) ⇒ small dispersion.
    omega = np.ones(6)
    log_omega_var = np.full(6, 1e-3)  # α = 1/log_omega_var = 1000 (very concentrated)
    disp = update_exposure_dispersion(omega, log_omega_var, floor=1e-8)
    assert disp < 0.1


def test_exposure_dispersion_high_variance_is_large():
    # Exposure posteriors spread far from 1 ⇒ large dispersion.
    omega = np.array([0.01, 0.02, 5.0, 8.0])
    log_omega_var = np.full(4, 1.0)  # α = 1 (diffuse posteriors)
    disp = update_exposure_dispersion(omega, log_omega_var, floor=1e-8)
    assert disp > 1.0


def test_exposure_dispersion_clips_to_ceil():
    # Extreme heterogeneity ⇒ dispersion pinned at the ceiling.
    omega = np.array([1e-4, 1e-4, 1e4, 1e4])
    log_omega_var = np.full(4, 1.0)
    disp = update_exposure_dispersion(omega, log_omega_var, floor=1e-8, ceil=50.0)
    assert disp == 50.0


def test_exposure_dispersion_empty_is_finite():
    disp = update_exposure_dispersion(np.array([]), np.array([]), floor=1e-8)
    assert disp > 0.0


def test_fit_rho_d_bb_balanced_is_small():
    # Sense counts at κ_d = 0.5 (5/10) → little overdispersion → ρ_d_bb small.
    rho = fit_rho_d_bb(np.array([5.0, 5.0, 5.0]), np.array([10.0, 10.0, 10.0]))
    assert rho < 0.1


def test_fit_rho_d_bb_all_or_nothing_is_large():
    # All-sense / all-antisense at n=10 → maximal overdispersion → ρ_d_bb large.
    rho = fit_rho_d_bb(np.array([0.0, 10.0, 0.0, 10.0]), np.array([10.0, 10.0, 10.0, 10.0]))
    assert rho > 0.5


def test_fit_rho_d_bb_too_few_obs_fallback():
    rho = fit_rho_d_bb(np.array([5.0]), np.array([10.0]))
    assert rho == _RHO_D_BB_FALLBACK


def test_update_pi_g_prior_balanced():
    # μ_g = ω ρ_0 L = 1·1·50 = 50; μ_d = n_u − μ_g = 50 ⇒ prior = 0.5.
    prior = update_pi_g_prior(np.array([1.0]), 1.0, np.array([50.0]), np.array([100.0]))
    np.testing.assert_allclose(prior, [0.5])


def test_update_pi_g_prior_gdna_dominated_clips():
    # μ_g (50) ≫ n_u (10) ⇒ all-gDNA ⇒ prior → 1, clipped to 1 − _PI_CLIP.
    prior = update_pi_g_prior(np.array([1.0]), 1.0, np.array([50.0]), np.array([10.0]))
    np.testing.assert_allclose(prior, [1.0 - _PI_CLIP])
