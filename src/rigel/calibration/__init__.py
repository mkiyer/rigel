"""Rigel calibration — Phase A burndown stub.

The v5 calibration surface was burned down on 2026-05-29 to make room for
the joint fractional-accumulator + calibration-v6 rewrite. See
``docs/acc_caljointmodel/00_implementation_plan.md`` for the phase plan.

Only three names remain public until the v6 calibrator lands:

- ``CalibrationConfig`` — sourced from :mod:`rigel.config`.
- ``CalibrationResult`` — empty placeholder; full schema lands in PR 2
  (see ``docs/caljointmodel/04_interface_contract.md``).
- ``calibrate`` — stub that raises :class:`NotImplementedError`.

``rigel quant`` therefore aborts at the calibration step until the v6
calibrator is wired end-to-end (PR 6).
"""

from __future__ import annotations

from dataclasses import dataclass

from ..config import CalibrationConfig

__all__ = ["CalibrationConfig", "CalibrationResult", "calibrate"]


@dataclass(frozen=True)
class CalibrationResult:
    """Placeholder for the Phase-D calibration output.

    The full schema (per-region deconvoluted mass, exposure posteriors,
    library hyperparameters) is defined in
    ``docs/caljointmodel/04_interface_contract.md`` and is implemented in PR 2.
    """


def calibrate(*args, **kwargs) -> CalibrationResult:  # pragma: no cover
    raise NotImplementedError(
        "rigel.calibration.calibrate is stubbed during the Phase A burndown. "
        "The replacement (joint fractional-accumulator + calibration-v6) "
        "lands in Phase D — see docs/acc_caljointmodel/00_implementation_plan.md."
    )
