"""Calibration strand-orientation helpers.

Numeric constants mirror ``src/rigel/native/calibration/orient.h``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

from ..types import Strand
from .regions import RegionStrand

if TYPE_CHECKING:
    from ..strand_model import StrandModel

ORIENT_SAME = 0
ORIENT_OPP = 1
ORIENT_UNINF = 2
ORIENT_N = 3


@dataclass(frozen=True, slots=True)
class StrandSummary:
    """Minimal strand-model summary used by calibration density correction."""

    p_r1_sense: float = 0.5
    n_observations: int = 0

    def __post_init__(self) -> None:
        p_r1_sense = float(self.p_r1_sense)
        if not math.isfinite(p_r1_sense) or not 0.0 <= p_r1_sense <= 1.0:
            raise ValueError(
                f"StrandSummary.p_r1_sense must be finite and in [0, 1]; got {self.p_r1_sense!r}."
            )
        n_observations = int(self.n_observations)
        if n_observations < 0:
            raise ValueError(
                f"StrandSummary.n_observations must be >= 0; got {self.n_observations!r}."
            )
        object.__setattr__(self, "p_r1_sense", p_r1_sense)
        object.__setattr__(self, "n_observations", n_observations)

    @property
    def p_r1_antisense(self) -> float:
        """Complement probability for read 1 aligning opposite to gene strand."""
        return 1.0 - self.p_r1_sense

    @property
    def signed_strand_contrast(self) -> float:
        """Signed contrast between sense and antisense read-1 orientation."""
        return 2.0 * self.p_r1_sense - 1.0

    @property
    def signed_strand_contrast_se(self) -> float:
        """Standard error of the signed strand contrast estimate."""
        if self.n_observations <= 0:
            return math.inf
        variance_p = self.p_r1_sense * self.p_r1_antisense / self.n_observations
        return 2.0 * math.sqrt(variance_p)

    def signed_strand_contrast_margin(self, *, confidence: float = 0.99) -> float:
        """Normal-approximation margin for detecting nonzero strand contrast."""
        z_by_confidence = {
            0.90: 1.644854,
            0.95: 1.959964,
            0.99: 2.575829,
        }
        try:
            z = z_by_confidence[float(confidence)]
        except KeyError as exc:
            raise ValueError("confidence must be one of 0.90, 0.95, or 0.99") from exc
        return z * self.signed_strand_contrast_se

    @property
    def strand_specificity(self) -> float:
        """Strand specificity in [0.5, 1.0]."""
        return max(self.p_r1_sense, self.p_r1_antisense)

    @property
    def read1_sense(self) -> bool:
        """True when the protocol is predominantly read-1 sense."""
        return self.p_r1_sense >= 0.5

    @classmethod
    def from_model(cls, model: "StrandModel") -> "StrandSummary":
        """Create a summary from a trained :class:`StrandModel`."""
        return cls(
            p_r1_sense=float(model.p_r1_sense),
            n_observations=int(model.n_observations),
        )

    @classmethod
    def uninformative(cls) -> "StrandSummary":
        """Return the no-strand-information summary."""
        return cls(p_r1_sense=0.5, n_observations=0)


def classify_orient(region_strand: int, fragment_strand: int) -> int:
    """Classify fragment strand relative to an unambiguous region strand."""
    if region_strand not in (int(RegionStrand.POS), int(RegionStrand.NEG)):
        return ORIENT_UNINF
    if fragment_strand not in (int(Strand.POS), int(Strand.NEG)):
        return ORIENT_UNINF
    return ORIENT_SAME if fragment_strand == region_strand else ORIENT_OPP
