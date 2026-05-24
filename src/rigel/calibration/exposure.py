"""Per-region gDNA exposure surface.

Phase 4 only needs the uniform exposure object so calibration can expose a
stable ``RegionExposure`` contract. Phase 7 replaces the uniform constructor
call with an unsupervised density-ratio estimate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np


__all__ = ["RegionExposure"]


@dataclass(frozen=True, slots=True)
class RegionExposure:
    """Per-region capture/gDNA exposure weights.

    ``A_r = 1.0`` means a region is fully exposed. The uniform constructor is
    the only producer until Phase 7 adds unsupervised exposure inference.
    """

    mode: Literal["uniform", "unsupervised"]
    A_r: np.ndarray  # float32[R] in (0, 1]
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
            "A_max": float(weights.max()) if weights.size else 1.0,
        }
