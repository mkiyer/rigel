"""CalibrationResult — the calibration output schema.

This is a placeholder during the calibration-v6 rebuild. The full schema
(per-region deconvoluted mass, exposure posteriors, the library
hyperparameters ρ_0 / φ / ρ_d / ρ_r / ε_s, and the gDNA / RNA models) is
defined in ``docs/caljointmodel/04_interface_contract.md`` and is filled in
in PR 2.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CalibrationResult:
    """Placeholder for the calibration output (real schema lands in PR 2)."""


__all__ = ["CalibrationResult"]
