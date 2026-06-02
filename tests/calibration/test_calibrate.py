"""calibrate(): the PR 5 outer-loop EM — convergence + structural invariants.

The loop runs E-step (count channel live) → exposure → AMBIG sweep → M-step until
the mass-change converges. These tests pin the *mechanics* (convergence, mass
conservation, parameter ranges, monotone diagnostic); converged *biology*
(paralog rescue, exon→RNA) needs realistic data and is covered by the scenario
suite (PR 7) — the shared 3-region synthetic is too sparse to fit a sensible RNA
strand dispersion.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from _synthetic import make_gdna_fl_pmf, make_synthetic_payload

from rigel.calibration import calibrate
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig


def _run(config=None):
    payload, ra = make_synthetic_payload()
    strand_model = SimpleNamespace(p_r1_sense=0.95, n_observations=40)
    return calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=make_gdna_fl_pmf(),
        config=config or CalibrationConfig(),
    )


def test_converges_with_monotone_mass_change():
    result = _run()
    assert isinstance(result, CalibrationResult)
    assert result.n_regions == 3
    assert result.converged is True
    assert 2 <= result.n_iterations <= CalibrationConfig().max_outer_iterations
    hist = result.mass_change_history
    assert hist.shape == (result.n_iterations,)
    assert np.all(np.isfinite(hist))
    # Non-increasing (the result also enforces this) and the last step is below tol.
    assert np.all(np.diff(hist) <= 1e-9 * (1.0 + np.abs(hist[:-1])))
    assert hist[-1] < CalibrationConfig().mass_rel_tol


def test_mass_conserved_per_view():
    result = _run()
    np.testing.assert_allclose(
        result.mass_g_contained + result.mass_d_contained, [15.0, 26.0, 15.0]
    )
    np.testing.assert_allclose(result.mass_g_left + result.mass_d_left, [0.0, 3.0, 1.5])
    np.testing.assert_allclose(result.mass_g_right + result.mass_d_right, [2.0, 4.5, 0.0])


def test_masses_bounded_and_exposure_positive():
    result = _run()
    for g, tot in (
        (result.mass_g_contained, np.array([15.0, 26.0, 15.0])),
        (result.mass_g_left, np.array([0.0, 3.0, 1.5])),
        (result.mass_g_right, np.array([2.0, 4.5, 0.0])),
    ):
        assert np.all(g >= -1e-9)
        assert np.all(g <= tot + 1e-9)
    assert np.all(result.omega > 0.0)
    assert np.all(result.log_omega_var > 0.0)


def test_fitted_hyperparameters_in_range():
    result = _run()
    assert result.rho_0 > 0.0
    assert result.exposure_dispersion > 0.0
    assert 0.0 < result.rho_d_bb < 1.0
    # κ_rna is the posterior mean (n_same+1)/(n_obs+2), pulled off the 0.95 MLE by
    # the Beta(1,1) prior (PR 9); exact value checked in test_strand_params_fixed_not_refit.
    assert 0.0 < result.kappa_rna < 1.0
    assert 0.0 < result.rho_r_bb < 1.0  # posterior-predictive overdispersion 1/(n_obs+3)


def test_strand_params_fixed_not_refit():
    # κ_rna and ρ_r_bb are the posterior-predictive strand fit (PR 9) and must
    # equal a one-shot fit on the same StrandModel — the M-step never touches them.
    from rigel.calibration.strand_balance import fit_strand_balance

    sb = fit_strand_balance(SimpleNamespace(p_r1_sense=0.95, n_observations=40))
    result = _run()
    assert result.kappa_rna == sb.kappa_rna
    assert result.rho_r_bb == sb.rho_r_bb
