"""calibrate(): the PR 4a single E-step pass — real deconvolution invariants.

The single pass uses ω=1, π_g_prior=0.5, μ_d=0 (count channel silent), so the
allocation is strand-driven; the count channel and the outer loop arrive in
PR 5. These tests pin the structural invariants, not converged biology.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from _synthetic import make_gdna_fl_pmf, make_synthetic_payload

from rigel.calibration import calibrate
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig


def _run():
    payload, ra = make_synthetic_payload()
    strand_model = SimpleNamespace(p_r1_sense=0.95)  # κ_rna source
    # No AMBIG region in the synthetic, so the PR 4b sweep is a no-op here.
    result = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=make_gdna_fl_pmf(),
        config=CalibrationConfig(),
    )
    return result


def test_returns_schema_valid_single_pass_result():
    result = _run()
    assert isinstance(result, CalibrationResult)
    assert result.n_regions == 3
    assert result.n_iterations == 1
    assert result.converged is False
    assert result.mass_change_history.shape == (1,)
    assert np.isfinite(result.mass_change_history).all()
    assert result.mass_change_history[0] >= 0.0


def test_real_hyperparameters():
    result = _run()
    # κ_rna is the StrandModel's p_r1_sense; ρ_r_bb is fit (PR 3).
    assert result.kappa_rna == 0.95
    assert 0.0 < result.rho_r_bb < 1.0
    # ρ_0 is the density seed: only r2 (intergenic) seeds → (1 + 16.5)/(1 + 100).
    np.testing.assert_allclose(result.rho_0, 17.5 / 101.0)
    assert result.phi > 0.0
    assert 0.0 < result.rho_d_bb < 1.0
    assert 0.0 < result.eps_s < 1.0


def test_mass_is_conserved_per_view():
    # Each view: M_g + M_d == that view's total mass (unspliced + spliced).
    result = _run()
    np.testing.assert_allclose(
        result.mass_g_contained + result.mass_d_contained, [15.0, 26.0, 15.0]
    )
    np.testing.assert_allclose(result.mass_g_left + result.mass_d_left, [0.0, 3.0, 1.5])
    np.testing.assert_allclose(result.mass_g_right + result.mass_d_right, [2.0, 4.5, 0.0])


def test_masses_are_bounded_and_nonnegative():
    result = _run()
    for g, tot in (
        (result.mass_g_contained, np.array([15.0, 26.0, 15.0])),
        (result.mass_g_left, np.array([0.0, 3.0, 1.5])),
        (result.mass_g_right, np.array([2.0, 4.5, 0.0])),
    ):
        assert np.all(g >= 0.0)
        assert np.all(g <= tot + 1e-9)


def test_exposure_is_positive_and_finite():
    result = _run()
    assert np.all(result.omega > 0.0)
    assert np.all(np.isfinite(result.omega))
    assert np.all(result.log_omega_var > 0.0)
