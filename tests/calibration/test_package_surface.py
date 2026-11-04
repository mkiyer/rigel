"""The public surface of ``rigel.calibration`` is exactly ``__all__``, and every name in it resolves.

The pipeline, the instruments and the rest of the test suite import the package through these names alone, so a name
that disappears, or one exported without existing, would otherwise surface only at a caller's import.
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
    # The locus-prior bridge is real: assemble_priors is a callable and LocusPriors its dataclass result. test_priors.py
    # covers its behaviour; this holds only the surface.
    assert callable(cal.assemble_priors)
    assert hasattr(cal.LocusPriors, "__dataclass_fields__")
