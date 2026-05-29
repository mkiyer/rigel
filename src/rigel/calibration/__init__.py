"""Rigel calibration — Phase A burndown stub.

The v5 calibration surface was burned down on 2026-05-29 to make room for
the joint fractional-accumulator + calibration-v6 rewrite. See
``docs/acc_caljointmodel/00_implementation_plan.md`` for the phase plan.

Only three names remain public until Phase D lands the new calibrator:

- ``CalibrationConfig`` — sourced from :mod:`rigel.config`.
- ``CalibrationResult`` — empty placeholder; full schema lands in Phase D
  (see ``docs/caljointmodel/04_outputs.md``).
- ``calibrate`` — stub that raises :class:`NotImplementedError`.

``rigel quant`` therefore aborts at the calibration step until Phase D.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..config import CalibrationConfig

__all__ = ["CalibrationConfig", "CalibrationResult", "calibrate"]


@dataclass(frozen=True)
class CalibrationResult:
    """Placeholder for the Phase-D calibration output.

    The full schema (state tensors, log-evidence ledger, exposure tables,
    FL models, etc.) is defined in ``docs/caljointmodel/04_outputs.md``
    and will be implemented in Phase D.
    """


def calibrate(*args, **kwargs) -> CalibrationResult:  # pragma: no cover
    raise NotImplementedError(
        "rigel.calibration.calibrate is stubbed during the Phase A burndown. "
        "The replacement (joint fractional-accumulator + calibration-v6) "
        "lands in Phase D — see docs/acc_caljointmodel/00_implementation_plan.md."
    )
