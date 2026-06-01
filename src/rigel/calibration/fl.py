"""rigel.calibration.fl — gDNA / RNA fragment-length distributions (PR 4c).

Produces the library-wide FL distributions the calibrator's **effective lengths**
need (boundary ``μ_FL``, region ``E_f[max(0, L−ℓ)]``): the **gDNA FL** (from
gDNA-dominated regions/boundaries) and the **RNA FL** (spliced fragments), each
**smoothly empirical-Bayes-shrunk** toward the global FL.

This is **not** a per-fragment FL likelihood — that channel is deliberately
excluded from calibration (``docs/caljointmodel/03_inference.md`` §2). Only gDNA
consumes a modelled length today (RNA cannot be FL-corrected per transcript); the
RNA FL is produced for PR 5's RNA effective length.

The gDNA pool aggregates the accumulator's **intergenic + intronic** FL pools
(both compartments); exonic pools (mature mRNA) are excluded. The FL pools and
this pool index match the C++ ``rigel::accumulator::fl_pool_idx`` convention.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..scan_payload import N_FL_POOLS

if TYPE_CHECKING:
    from ..scan_payload import AccumulatorPayload

__all__ = [
    "FLModels",
    "POOL_EB_PRIOR_ESS",
    "N_FL_POOLS",
    "gdna_fl_mass",
    "build_fl_models",
]

# gDNA FL-pool indices = region_type * 2 + compartment, matching the C++
# fl_pool_idx (region_type 0=intergenic, 1=intronic, 2=exonic; compartment
# 0=contained, 1=boundary).
FL_POOL_INTERGENIC_CONTAINED = 0
FL_POOL_INTERGENIC_BOUNDARY = 1
FL_POOL_INTRONIC_CONTAINED = 2
FL_POOL_INTRONIC_BOUNDARY = 3
FL_POOL_EXONIC_CONTAINED = 4
FL_POOL_EXONIC_BOUNDARY = 5

#: The gDNA pool = intergenic + intronic, both compartments (PR04c §1).
_GDNA_POOLS = (
    FL_POOL_INTERGENIC_CONTAINED,
    FL_POOL_INTERGENIC_BOUNDARY,
    FL_POOL_INTRONIC_CONTAINED,
    FL_POOL_INTRONIC_BOUNDARY,
)

#: Dirichlet pseudo-count for the smooth EB shrink toward the global FL. Not a
#: cliff: ``pool_total ≫ ρ_ess`` → empirical; ``≪`` → the global anchor; ``= 0``
#: → global. Revisit the default on real-data pool sizes (PR04c decision 5).
POOL_EB_PRIOR_ESS: float = 1000.0


@dataclass(frozen=True, slots=True)
class FLModels:
    """Library-wide FL pmfs (float64[max_size + 1]) + their pool totals."""

    global_pmf: np.ndarray  # unconditional anchor (no prior)
    rna_pmf: np.ndarray  # spliced, EB-shrunk toward global
    gdna_pmf: np.ndarray  # gDNA pool, EB-shrunk toward global
    n_global: float
    n_rna: float
    n_gdna: float
    max_size: int


def gdna_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """Aggregate gDNA FL mass from the accumulator pools (intergenic + intronic).

    Returns float64[fl_max_size + 1], or an empty vector when the payload carries
    no FL pools (FL pooling disabled) — the caller then falls back to global FL.
    """
    fp = payload.fl_pool_mass
    if fp is None:
        return np.zeros(0, dtype=np.float64)
    return np.asarray(fp, dtype=np.float64)[list(_GDNA_POOLS)].sum(axis=0)


def _aligned(counts: np.ndarray, max_size: int) -> np.ndarray:
    """Align a raw FL count vector to ``max_size + 1`` bins (overflow → last bin)."""
    out = np.zeros(max_size + 1, dtype=np.float64)
    c = np.asarray(counts, dtype=np.float64)
    n = min(out.size, c.size)
    out[:n] = c[:n]
    if c.size > out.size:
        out[max_size] += float(np.sum(c[out.size :], dtype=np.float64))
    return out


def _normalized(v: np.ndarray) -> np.ndarray:
    s = float(v.sum())
    return v / s if s > 0.0 else np.full(v.size, 1.0 / v.size, dtype=np.float64)


def _smooth_eb(counts: np.ndarray, global_pmf: np.ndarray, max_size: int, prior_ess: float):
    """Smooth EB pmf: ``(counts + ρ_ess·global_pmf) / (total + ρ_ess)``.

    Continuous in the pool total — no quality threshold (PR04c decision 5).
    Returns ``(pmf, pool_total)``.
    """
    aligned = _aligned(counts, max_size)
    total = float(aligned.sum())
    denom = total + prior_ess
    pmf = (aligned + prior_ess * global_pmf) / denom if denom > 0.0 else global_pmf.copy()
    return _normalized(pmf), total


def build_fl_models(
    *,
    global_counts: np.ndarray,
    rna_counts: np.ndarray,
    gdna_counts: np.ndarray,
    max_size: int,
    prior_ess: float = POOL_EB_PRIOR_ESS,
) -> FLModels:
    """Build the global / RNA / gDNA FL pmfs under one smooth-EB policy.

    ``global_counts`` is the scanner's unconditional histogram (the EB anchor),
    ``rna_counts`` the spliced (SPLICED-ANNOT) histogram, ``gdna_counts`` the
    aggregated gDNA pool (:func:`gdna_fl_mass`). RNA and gDNA are EB-shrunk toward
    the global FL with a single Dirichlet ``prior_ess``.
    """
    global_pmf = _normalized(_aligned(global_counts, max_size))
    rna_pmf, n_rna = _smooth_eb(rna_counts, global_pmf, max_size, prior_ess)
    gdna_pmf, n_gdna = _smooth_eb(gdna_counts, global_pmf, max_size, prior_ess)
    return FLModels(
        global_pmf=global_pmf,
        rna_pmf=rna_pmf,
        gdna_pmf=gdna_pmf,
        n_global=float(np.asarray(global_counts, dtype=np.float64).sum()),
        n_rna=n_rna,
        n_gdna=n_gdna,
        max_size=int(max_size),
    )
