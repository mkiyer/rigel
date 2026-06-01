"""CalibrationResult.__post_init__ intrinsic invariants."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.errors import CalibrationConvergenceError
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig


def _valid_kwargs(n_regions: int = 2) -> dict:
    z = np.zeros(n_regions, dtype=np.float64)
    o = np.ones(n_regions, dtype=np.float64)
    return dict(
        mass_g_contained=z.copy(),
        mass_d_contained=o.copy(),
        mass_g_left=z.copy(),
        mass_d_left=z.copy(),
        mass_g_right=z.copy(),
        mass_d_right=z.copy(),
        omega=o.copy(),
        log_omega_var=o.copy(),
        rho_0=1e-3,
        exposure_dispersion=0.2,
        rho_d_bb=0.01,
        kappa_rna=0.9,
        rho_r_bb=0.01,
        eps_s=1e-3,
        n_iterations=0,
        converged=True,
        mass_change_history=np.empty(0, dtype=np.float64),
        n_regions=n_regions,
        config=CalibrationConfig(),
    )


def test_valid_result_constructs():
    CalibrationResult(**_valid_kwargs())


@pytest.mark.parametrize(
    "field,value",
    [
        ("mass_g_contained", np.array([-1.0, 0.0])),  # negative mass
        ("omega", np.array([0.0, 1.0])),  # non-positive exposure
        ("log_omega_var", np.array([1.0, -1.0])),  # non-positive variance
    ],
)
def test_rejects_bad_region_arrays(field, value):
    kw = _valid_kwargs()
    kw[field] = value
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


def test_rejects_non_float64_array():
    kw = _valid_kwargs()
    kw["mass_d_contained"] = np.array([1, 1], dtype=np.int64)
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


def test_rejects_length_mismatch():
    kw = _valid_kwargs()
    kw["mass_g_contained"] = np.zeros(3, dtype=np.float64)  # n_regions is 2
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


@pytest.mark.parametrize(
    "field,value",
    [
        ("rho_0", 0.0),
        ("exposure_dispersion", -1.0),
        ("rho_d_bb", 1.5),
        ("kappa_rna", 0.0),
        ("rho_r_bb", 0.0),
        ("eps_s", 1.0),
    ],
)
def test_rejects_bad_hyperparameters(field, value):
    kw = _valid_kwargs()
    kw[field] = value
    with pytest.raises(ValueError):
        CalibrationResult(**kw)


def test_mass_change_history_non_increasing_ok():
    kw = _valid_kwargs()
    kw["n_iterations"] = 3
    kw["mass_change_history"] = np.array([0.3, 0.1, 0.05], dtype=np.float64)
    CalibrationResult(**kw)  # decreasing → fine


def test_mass_change_history_increase_allowed():
    # Strict monotonicity is NOT required: the mass change legitimately spikes at
    # iteration 2 when the count channel activates (doc 03 §3.1), before
    # converging. A finite, non-monotone history must construct fine.
    kw = _valid_kwargs()
    kw["n_iterations"] = 3
    kw["mass_change_history"] = np.array([0.14, 57.3, 8e-5], dtype=np.float64)
    CalibrationResult(**kw)


def test_mass_change_history_non_finite_raises():
    # The divergence sentinel: a non-finite mass change means the EM blew up.
    kw = _valid_kwargs()
    kw["n_iterations"] = 2
    kw["mass_change_history"] = np.array([0.1, np.inf], dtype=np.float64)
    with pytest.raises(CalibrationConvergenceError):
        CalibrationResult(**kw)


def test_mass_change_history_length_must_match_iterations():
    kw = _valid_kwargs()
    kw["n_iterations"] = 2
    kw["mass_change_history"] = np.array([0.1], dtype=np.float64)
    with pytest.raises(ValueError):
        CalibrationResult(**kw)
