"""Public surface of rigel.calibration."""

from __future__ import annotations

import pytest

import rigel.calibration as cal

EXPECTED_NAMES = [
    "CalibrationConfig",
    "CalibrationResult",
    "CalibrationSubstrate",
    "SubstrateView",
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


def test_assemble_priors_still_stubbed():
    # The locus-prior bridge lands in PR 6; everything else is implemented.
    with pytest.raises(NotImplementedError):
        cal.assemble_priors()
