"""Per-region gDNA exposure surface.

In v4 calibration, ``RegionExposure`` is the consumer-facing view of the
gDNA density evidence. The uniform constructor remains for tests and
degenerate paths; the orchestrator uses :meth:`RegionExposure.from_density`
to lift a :class:`DensityEvidence` into the consumer-side exposure surface.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np


if TYPE_CHECKING:  # pragma: no cover - import for annotations only.
    from .density_model import DensityEvidence


__all__ = ["RegionExposure"]


@dataclass(frozen=True, slots=True)
class RegionExposure:
    """Per-region capture/gDNA exposure weights.

    ``A_r`` is the per-region relative exposure: 1.0 means the region matches
    the reference density ``rho_ref``. Values may exceed 1.0 when constructed
    from density evidence on regions with above-reference gDNA density.
    """

    mode: Literal["uniform", "density"]
    A_r: np.ndarray  # float32[R], non-negative; may exceed 1 when constructed from density evidence
    rho_r: np.ndarray  # float32[R] per-region gDNA density, frag/bp
    rho_ref: float
    reference_quantile: float
    eligible: np.ndarray  # bool[R]
    flags: np.ndarray  # uint8[R]

    @classmethod
    def uniform(cls, R: int) -> "RegionExposure":
        """Return an all-ones exposure surface for ``R`` fine regions."""
        n_regions = int(R)
        if n_regions < 0:
            raise ValueError(f"RegionExposure.uniform: R must be >= 0; got {R!r}")
        return cls(
            mode="uniform",
            A_r=np.ones(n_regions, dtype=np.float32),
            rho_r=np.zeros(n_regions, dtype=np.float32),
            rho_ref=0.0,
            reference_quantile=0.0,
            eligible=np.ones(n_regions, dtype=bool),
            flags=np.zeros(n_regions, dtype=np.uint8),
        )

    @classmethod
    def from_density(
        cls,
        density_evidence: "DensityEvidence",
        *,
        max_exposure: float | None = None,
    ) -> "RegionExposure":
        """Lift a :class:`DensityEvidence` into the consumer-side exposure surface.

        Parameters
        ----------
        density_evidence:
            Compact per-region evidence produced by ``fit_density_evidence``.
            Its ``relative_exposure`` array provides the per-region ``A_r``
            and its ``rho_post`` array provides the per-region ``rho_r``.
        max_exposure:
            Optional upper clip on ``A_r``. When provided, ``A_r`` is clipped
            in-place with ``np.minimum``. Phase 4 leaves this unset by
            default; the CLI/config knob lands in Phase 6.

        Notes
        -----
        Regions flagged with ``FLAG_LOW_BOUNDARY_OPPORTUNITY`` are marked
        ``eligible=False`` but keep their posterior ``A_r`` and ``rho_r``
        values (no overwrite). Locus EM is responsible for deciding how to
        treat ineligible regions.
        """
        from .density_model import FLAG_LOW_BOUNDARY_OPPORTUNITY

        A_r = np.asarray(density_evidence.relative_exposure, dtype=np.float32)
        if max_exposure is not None:
            cap = float(max_exposure)
            if not np.isfinite(cap) or cap <= 0.0:
                raise ValueError(
                    f"RegionExposure.from_density: max_exposure must be positive; got {cap!r}."
                )
            A_r = np.minimum(A_r, np.float32(cap))
        rho_r = np.asarray(density_evidence.rho_post, dtype=np.float32)
        flags = np.asarray(density_evidence.flags, dtype=np.uint8)
        eligible = (flags & FLAG_LOW_BOUNDARY_OPPORTUNITY) == 0
        return cls(
            mode="density",
            A_r=A_r,
            rho_r=rho_r,
            rho_ref=float(density_evidence.rho_ref),
            reference_quantile=0.0,
            eligible=eligible,
            flags=flags,
        )

    def to_summary_dict(self) -> dict[str, object]:
        """Return a compact JSON-safe exposure summary."""
        weights = np.asarray(self.A_r, dtype=np.float64)
        return {
            "mode": self.mode,
            "n_regions": int(weights.size),
            "n_regions_eligible": int(np.asarray(self.eligible, dtype=bool).sum()),
            "rho_ref": float(self.rho_ref),
            "reference_quantile": float(self.reference_quantile),
            "A_min": float(weights.min()) if weights.size else 1.0,
            "A_mean": float(weights.mean()) if weights.size else 1.0,
            "A_p99": float(np.quantile(weights, 0.99)) if weights.size else 1.0,
            "A_max": float(weights.max()) if weights.size else 1.0,
        }
