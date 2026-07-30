"""Public surface of rigel.calibration."""

from __future__ import annotations

import rigel.calibration as cal

EXPECTED_NAMES = [
    "CalibrationConfig",
    "CalibrationResult",
    "CalibrationSubstrate",
    "PopulationView",
    "StrandBalance",
    "CalibrationSubstrateError",
    "calibrate",
    "assemble_priors",
    "LocusPriors",
]


def test_public_names_present():
    for name in EXPECTED_NAMES:
        assert hasattr(cal, name), f"rigel.calibration is missing {name!r}"
    assert set(cal.__all__) == set(EXPECTED_NAMES)


def test_error_types():
    assert issubclass(cal.CalibrationSubstrateError, ValueError)


def test_assemble_priors_implemented():
    # The locus-prior bridge landed in PR 6: assemble_priors is a real callable
    # and LocusPriors is its dataclass result (behavior tested in test_priors.py).
    assert callable(cal.assemble_priors)
    assert hasattr(cal.LocusPriors, "__dataclass_fields__")
