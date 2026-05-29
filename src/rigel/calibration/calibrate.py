"""calibrate() — the joint fractional-accumulator + calibration-v6 orchestrator.

Consumes the accumulator payload, the region geometry, and the trained strand
model and returns a :class:`~rigel.calibration.result.CalibrationResult` with
the recovered library hyperparameters and per-region deconvoluted mass.

Skeleton during PR 1: the signature is wired with its real inputs (so
``run_pipeline`` reaches it with live arguments) but the body raises. The
implementation is filled in across PR 2–PR 5; see
``docs/acc_caljointmodel/00_implementation_plan.md``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .result import CalibrationResult

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays


def calibrate(
    *,
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    config: "CalibrationConfig",
) -> CalibrationResult:
    raise NotImplementedError(
        "rigel.calibration.calibrate is stubbed during the calibration-v6 rebuild. "
        "The joint deconvolution (substrate → strand balance → E-step/exposure → "
        "M-step/outer loop) lands in PR 2–PR 5 — see "
        "docs/acc_caljointmodel/00_implementation_plan.md."
    )


__all__ = ["calibrate"]
