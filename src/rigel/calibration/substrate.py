"""CalibrationSubstrate — the calibrator-facing view of the accumulator payload.

The substrate adapts the raw :class:`~rigel.scan_payload.AccumulatorPayload`
(per-region contained counts + per-boundary mass/flux) plus the region
geometry into the count/mass sufficient statistics the calibrator consumes:
the channel reductions and the D1 boundary→region mass attribution
(``region r ← boundary_mass_right[left_boundary(r)] +
boundary_mass_left[right_boundary(r)]``), built on the index mapping in
:mod:`rigel.calibration.region_arrays`.

Skeleton during PR 1; the real adapter lands in PR 2.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..scan_payload import AccumulatorPayload
    from .region_arrays import RegionArrays


@dataclass(frozen=True, slots=True)
class CalibrationSubstrate:
    """Placeholder for the calibrator-facing substrate (real adapter in PR 2)."""

    @classmethod
    def from_payload(
        cls,
        payload: "AccumulatorPayload",
        region_arrays: "RegionArrays",
    ) -> "CalibrationSubstrate":
        raise NotImplementedError(
            "CalibrationSubstrate is a PR 1 skeleton; the payload→substrate adapter "
            "(channel reductions + boundary→region mass attribution) lands in PR 2 — "
            "see docs/acc_caljointmodel/00_implementation_plan.md."
        )


__all__ = ["CalibrationSubstrate"]
