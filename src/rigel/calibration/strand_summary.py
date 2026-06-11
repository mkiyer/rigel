"""StrandSummary — a public summary of :class:`StrandModel`.

Independent of region orient / per-fragment orient classification. Lives
in its own module so it can be imported without pulling in calibration
internals.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..strand_model import StrandModel


# Numerical floor for the signed strand contrast |2 p_r1_sense - 1| below
# which the spliced strand model is considered unidentifiable.
STRAND_CONTRAST_NUMERICAL_FLOOR: float = 1e-3

_Z_BY_CONFIDENCE = {0.90: 1.644854, 0.95: 1.959964, 0.99: 2.575829}


def strand_contrast_identifiable(
    p_r1_sense: float, n_observations: int, *, confidence: float = 0.99
) -> bool:
    """True when the strand contrast ``|2·p_r1_sense − 1|`` is distinguishable from 0 at ``confidence``.

    The strand channel can separate gDNA (sense ½) from RNA (sense κ) only when the observed sense
    rate is detectably off ½ given its sample size: ``|2p−1|`` must exceed both a numerical floor and
    its normal-approximation margin ``z·SE`` (``SE = 2·sqrt(p(1−p)/n)``). Works off the two universal
    strand-model scalars (``p_r1_sense``, ``n_observations``) so it applies to the ``StrandModels``
    container, a single ``StrandModel``, and :class:`StrandSummary` alike.
    """
    p = float(p_r1_sense)
    n = int(n_observations)
    if n <= 0:
        return False
    try:
        z = _Z_BY_CONFIDENCE[float(confidence)]
    except KeyError as exc:
        raise ValueError("confidence must be one of 0.90, 0.95, or 0.99") from exc
    se = 2.0 * math.sqrt(p * (1.0 - p) / n)
    return abs(2.0 * p - 1.0) >= max(STRAND_CONTRAST_NUMERICAL_FLOOR, z * se)


@dataclass(frozen=True, slots=True)
class StrandSummary:
    """Minimal strand-model summary used by calibration density correction."""

    p_r1_sense: float = 0.5
    n_observations: int = 0
    n_same: int = 0
    n_opposite: int = 0

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
        n_same = int(self.n_same)
        n_opposite = int(self.n_opposite)
        if n_same < 0:
            raise ValueError(f"StrandSummary.n_same must be >= 0; got {self.n_same!r}.")
        if n_opposite < 0:
            raise ValueError(f"StrandSummary.n_opposite must be >= 0; got {self.n_opposite!r}.")
        if n_same + n_opposite != n_observations:
            raise ValueError(
                "StrandSummary.n_same + n_opposite must equal n_observations; "
                f"got {n_same} + {n_opposite} != {n_observations}."
            )
        if n_observations > 0:
            observed_p = n_same / n_observations
            if not math.isclose(p_r1_sense, observed_p, rel_tol=0.0, abs_tol=1e-12):
                raise ValueError(
                    "StrandSummary.p_r1_sense must equal n_same / n_observations; "
                    f"got {p_r1_sense!r}, expected {observed_p!r}."
                )
        object.__setattr__(self, "p_r1_sense", p_r1_sense)
        object.__setattr__(self, "n_observations", n_observations)
        object.__setattr__(self, "n_same", n_same)
        object.__setattr__(self, "n_opposite", n_opposite)

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
        try:
            z = _Z_BY_CONFIDENCE[float(confidence)]
        except KeyError as exc:
            raise ValueError("confidence must be one of 0.90, 0.95, or 0.99") from exc
        return z * self.signed_strand_contrast_se

    def is_identifiable(self, *, confidence: float = 0.99) -> bool:
        """True when the strand contrast is distinguishable from 0 at ``confidence`` (see
        :func:`strand_contrast_identifiable`)."""
        return strand_contrast_identifiable(
            self.p_r1_sense, self.n_observations, confidence=confidence
        )

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
            n_same=int(model.n_same),
            n_opposite=int(model.n_opposite),
        )

    @classmethod
    def uninformative(cls) -> "StrandSummary":
        """Return the no-strand-information summary."""
        return cls(p_r1_sense=0.5, n_observations=0)


__all__ = ["StrandSummary", "STRAND_CONTRAST_NUMERICAL_FLOOR", "strand_contrast_identifiable"]
