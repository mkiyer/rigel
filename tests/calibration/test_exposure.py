"""Closed-form Gamma exposure posterior (G4) with D1 side-attribution."""

from __future__ import annotations

import numpy as np

from rigel.calibration.exposure import Exposure, exposure_posterior


def test_closed_form_matches_hand_calc():
    # M_g_tot = 10 + 2 + 3 = 15; inv_phi = 1/0.5 = 2.
    # α = 2 + 15 = 17; β = 2 + 0.1*100 = 12; ω = 17/12; Var(log ω) = 1/17.
    exp = exposure_posterior(
        np.array([10.0]),
        np.array([2.0]),
        np.array([3.0]),
        rho_0=0.1,
        L_eff=np.array([100.0]),
        phi=0.5,
    )
    assert isinstance(exp, Exposure)
    np.testing.assert_allclose(exp.m_g_tot, [15.0])
    np.testing.assert_allclose(exp.omega, [17.0 / 12.0])
    np.testing.assert_allclose(exp.log_omega_var, [1.0 / 17.0])


def test_d1_aggregation_has_no_half_split():
    # D1: M_g_tot = cont + left + right (NOT cont + ½(left + right)).
    exp = exposure_posterior(
        np.array([10.0]),
        np.array([2.0]),
        np.array([4.0]),
        rho_0=1.0,
        L_eff=np.array([1.0]),
        phi=1.0,
    )
    np.testing.assert_allclose(exp.m_g_tot, [16.0])  # 10+2+4, not 10+0.5*6=13


def test_empty_region_gets_prior_exposure():
    # Zero gDNA mass → ω = (1/φ)/(1/φ + ρ_0 L) ; with ρ_0 L = 0 → ω = 1 (prior mean).
    exp = exposure_posterior(
        np.array([0.0]),
        np.array([0.0]),
        np.array([0.0]),
        rho_0=0.01,
        L_eff=np.array([0.0]),
        phi=0.2,
    )
    np.testing.assert_allclose(exp.omega, [1.0])
    np.testing.assert_allclose(exp.log_omega_var, [0.2])  # 1/α = 1/(1/φ) = φ


def test_vectorized_over_regions():
    exp = exposure_posterior(
        np.array([10.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        rho_0=0.1,
        L_eff=np.array([100.0, 50.0]),
        phi=1.0,
    )
    assert exp.omega.shape == (2,)
    assert np.all(exp.omega > 0.0)
    assert np.all(exp.log_omega_var > 0.0)
