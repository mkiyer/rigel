"""calibrate(): real strand model + placeholder (no-gDNA) deconvolution."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from _synthetic import make_synthetic_payload

from rigel.calibration import calibrate
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig


def test_calibrate_strand_fit_with_placeholder_gdna():
    payload, ra = make_synthetic_payload()
    strand_model = SimpleNamespace(p_r1_sense=0.95)  # κ_rna source
    result = calibrate(
        payload=payload, region_arrays=ra, strand_model=strand_model, config=CalibrationConfig()
    )

    assert isinstance(result, CalibrationResult)
    assert result.n_regions == payload.r_total == 3

    # Strand model is real now: κ_rna comes from the StrandModel.
    assert result.kappa_rna == 0.95
    assert 0.0 < result.rho_r_bb < 1.0

    # gDNA deconvolution is still a placeholder: no gDNA, all mass → RNA.
    np.testing.assert_array_equal(result.mass_g_contained, [0.0, 0.0, 0.0])
    np.testing.assert_array_equal(result.mass_g_left, [0.0, 0.0, 0.0])
    np.testing.assert_array_equal(result.mass_g_right, [0.0, 0.0, 0.0])
    np.testing.assert_allclose(result.mass_d_contained, [15.0, 26.0, 15.0])
    np.testing.assert_allclose(result.mass_d_left, [0.0, 3.0, 1.5])
    np.testing.assert_allclose(result.mass_d_right, [2.0, 4.5, 0.0])
    np.testing.assert_array_equal(result.omega, [1.0, 1.0, 1.0])
    assert result.converged is True
    assert result.n_iterations == 0
    assert result.mass_change_history.shape == (0,)
    assert result.rho_0 > 0.0 and result.phi > 0.0
    assert 0.0 < result.rho_d_bb < 1.0
    assert 0.0 < result.eps_s < 1.0
    np.testing.assert_allclose(result.log_omega_var, np.full(3, result.phi))
