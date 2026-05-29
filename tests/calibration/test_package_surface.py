"""Public surface of rigel.calibration — names present, skeletons raise."""

from __future__ import annotations

import pytest

import rigel.calibration as cal

EXPECTED_NAMES = [
    "CalibrationConfig",
    "CalibrationResult",
    "CalibrationSubstrate",
    "CalibrationSubstrateError",
    "CalibrationConvergenceError",
    "calibrate",
    "assemble_priors",
]


def test_public_names_present():
    for name in EXPECTED_NAMES:
        assert hasattr(cal, name), f"rigel.calibration is missing {name!r}"
    assert set(cal.__all__) == set(EXPECTED_NAMES)


def test_error_types():
    assert issubclass(cal.CalibrationSubstrateError, ValueError)
    assert issubclass(cal.CalibrationConvergenceError, RuntimeError)


def test_calibration_result_instantiable():
    # Empty placeholder dataclass; the real schema lands in PR 2.
    assert cal.CalibrationResult() is not None


def test_skeletons_raise_not_implemented():
    with pytest.raises(NotImplementedError):
        cal.calibrate(payload=None, region_arrays=None, strand_model=None, config=None)
    with pytest.raises(NotImplementedError):
        cal.assemble_priors()
    with pytest.raises(NotImplementedError):
        cal.CalibrationSubstrate.from_payload(None, None)
