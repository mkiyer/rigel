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
        phi=0.2,
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
        ("phi", -1.0),
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


def test_mass_change_history_increasing_raises():
    kw = _valid_kwargs()
    kw["n_iterations"] = 2
    kw["mass_change_history"] = np.array([0.1, 0.2], dtype=np.float64)  # increases
    with pytest.raises(CalibrationConvergenceError):
        CalibrationResult(**kw)


def test_mass_change_history_length_must_match_iterations():
    kw = _valid_kwargs()
    kw["n_iterations"] = 2
    kw["mass_change_history"] = np.array([0.1], dtype=np.float64)
    with pytest.raises(ValueError):
        CalibrationResult(**kw)
