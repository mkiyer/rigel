"""The public surface of ``rigel.calibration`` is exactly ``__all__``, and every name in it resolves.

The package is imported by the pipeline, by the instruments and by the rest of the test suite through
these names alone, so a name that disappears or one that is exported without existing is a break that
would otherwise only surface at a caller's import.
"""

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
    # The locus-prior bridge is real: assemble_priors is a callable and LocusPriors is its dataclass
    # result. Its behaviour is gated in test_priors.py; this only holds the surface.
    assert callable(cal.assemble_priors)
    assert hasattr(cal.LocusPriors, "__dataclass_fields__")
