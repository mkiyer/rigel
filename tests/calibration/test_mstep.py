"""M-step hyperparameter updates (doc 03 §5): closed forms + bounded 1-D fits."""

from __future__ import annotations

import numpy as np

from rigel.calibration.estep import _PI_CLIP
from rigel.calibration.mstep import (
    _RHO_D_BB_FALLBACK,
    fit_phi,
    fit_rho_d_bb,
    update_eps_s,
    update_pi_g_prior,
    update_rho_0,
)
from rigel.calibration.strand_balance import _BB_FLOOR


def test_update_rho_0_closed_form():
    # Σ M_g_tot / Σ(ω·L) = (10+20) / (1·100 + 2·100) = 30/300 = 0.1.
    rho_0 = update_rho_0(np.array([10.0, 20.0]), np.array([1.0, 2.0]), np.array([100.0, 100.0]))
    np.testing.assert_allclose(rho_0, 0.1)


def test_update_rho_0_empty_exposure_floor():
    rho_0 = update_rho_0(np.array([0.0]), np.array([0.0]), np.array([100.0]))
    assert rho_0 > 0.0  # 1 / max(ΣL, 1)


def test_update_eps_s_beta11():
    # π_g = 0 everywhere ⇒ no spliced attributed to gDNA ⇒ ε_s = 1/(1+Σn_s).
    eps = update_eps_s(np.array([0.0, 0.0]), np.array([10.0, 10.0]))
    np.testing.assert_allclose(eps, 1.0 / 21.0)
    # π_g = 1 everywhere ⇒ all spliced gDNA-attributed ⇒ ε_s → 1, clipped.
    eps_hi = update_eps_s(np.array([1.0, 1.0]), np.array([10.0, 10.0]))
    assert eps_hi == 1.0 - _BB_FLOOR


def test_fit_phi_underdispersed_floors():
    # All counts exactly at the mean (variance 0) → no overdispersion → φ → floor.
    n_u = np.full(8, 50, dtype=np.int64)
    mu_g = np.full(8, 50.0)
    phi = fit_phi(n_u, mu_g, np.zeros(8), phi_floor=1e-8)
    assert phi < 1e-3


def test_fit_phi_overdispersed_is_large():
    # Bimodal counts far from the mean → large overdispersion.
    n_u = np.array([0, 0, 100, 100], dtype=np.int64)
    mu_g = np.full(4, 50.0)
    phi = fit_phi(n_u, mu_g, np.zeros(4), phi_floor=1e-8)
    assert phi > 0.5


def test_fit_phi_no_data_returns_floor():
    assert fit_phi(np.zeros(3, dtype=np.int64), np.zeros(3), np.zeros(3), phi_floor=1e-6) == 1e-6


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
