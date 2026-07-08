"""rigel.calibration.fl — gDNA / RNA fragment-length distributions (PR 4c).

Produces the library-wide FL distributions the calibrator's **effective lengths**
need (boundary ``fl_mean``, region ``E_f[max(0, L−ℓ)]``): the **gDNA FL** (from
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
    from ..frag_length_model import FragmentLengthModel
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
#: cliff: ``pool_total ≫ prior_ess`` → empirical; ``≪`` → the global anchor; ``= 0``
#: → global. Revisit the default on real-data pool sizes (PR04c decision 5).
POOL_EB_PRIOR_ESS: float = 1000.0


@dataclass(frozen=True, slots=True)
class FLModels:
    """Library-wide FL distributions (float64[max_size + 1]) + their pool totals.

    Two views of each distribution are carried:

    * ``*_pmf`` — the EB-smoothed pmf used for **scoring / calibration** (RNA and
      gDNA are shrunk toward ``global_pmf``; see :func:`build_fl_models`).
    * ``*_counts`` — the **raw, unsmoothed** histograms (aligned to
      ``max_size + 1`` bins) the pmfs were built from. These are the honest
      empirical distributions surfaced as QC (:meth:`rna_model` / :meth:`gdna_model`
      / :meth:`global_model`); the EB smoothing is a <1% perturbation at real
      library scale and stays internal to scoring.

    ``gdna_counts`` is the **structural-pool proxy** (intergenic + intronic FL
    mass; :func:`gdna_fl_mass`) — i.e. gDNA *plus* intronic/nascent RNA — not a
    deconvolved pure-gDNA distribution. ``rna_counts`` is the spliced-annotated
    (SPLICED-ANNOT) histogram.
    """

    global_pmf: np.ndarray  # unconditional anchor (no prior)
    rna_pmf: np.ndarray  # spliced, EB-shrunk toward global
    gdna_pmf: np.ndarray  # gDNA pool, EB-shrunk toward global
    global_counts: np.ndarray  # raw scanner histogram (all fragments)
    rna_counts: np.ndarray  # raw spliced-annotated histogram
    gdna_counts: np.ndarray  # raw structural pool (intergenic + intronic)
    n_global: float
    n_rna: float
    n_gdna: float
    max_size: int

    def _empirical(self, counts: np.ndarray) -> "FragmentLengthModel":
        """Wrap a raw count vector as an unfinalized ``FragmentLengthModel``.

        Unfinalized on purpose: summary stats (mean/std/median/mode) and
        ``to_dict()`` then read the *raw* histogram, matching how the scanner's
        global + per-category models are reported (no scoring-side pseudo-count
        smoothing enters the QC numbers).
        """
        from ..frag_length_model import FragmentLengthModel

        return FragmentLengthModel(max_size=self.max_size, counts=counts.copy())

    def global_model(self) -> "FragmentLengthModel":
        """Empirical global FL distribution as a ``FragmentLengthModel`` (QC)."""
        return self._empirical(self.global_counts)

    def rna_model(self) -> "FragmentLengthModel":
        """Empirical RNA (spliced-annotated) FL distribution (QC)."""
        return self._empirical(self.rna_counts)

    def gdna_model(self) -> "FragmentLengthModel":
        """Empirical gDNA structural-pool FL distribution (QC).

        The pool is intergenic + intronic FL mass (:func:`gdna_fl_mass`): a
        gDNA-dominated proxy that also includes intronic/nascent RNA, not a
        deconvolved pure-gDNA distribution.
        """
        return self._empirical(self.gdna_counts)


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


def _smooth_eb(aligned: np.ndarray, global_pmf: np.ndarray, prior_ess: float):
    """Smooth EB pmf: ``(counts + prior_ess·global_pmf) / (total + prior_ess)``.

    ``aligned`` is an already ``max_size``-aligned count vector. Continuous in the
    pool total — no quality threshold (PR04c decision 5). Returns ``(pmf, pool_total)``.
    """
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
    global_aligned = _aligned(global_counts, max_size)
    rna_aligned = _aligned(rna_counts, max_size)
    gdna_aligned = _aligned(gdna_counts, max_size)
    global_pmf = _normalized(global_aligned)
    rna_pmf, n_rna = _smooth_eb(rna_aligned, global_pmf, prior_ess)
    gdna_pmf, n_gdna = _smooth_eb(gdna_aligned, global_pmf, prior_ess)
    return FLModels(
        global_pmf=global_pmf,
        rna_pmf=rna_pmf,
        gdna_pmf=gdna_pmf,
        global_counts=global_aligned,
        rna_counts=rna_aligned,
        gdna_counts=gdna_aligned,
        n_global=float(global_aligned.sum()),
        n_rna=n_rna,
        n_gdna=n_gdna,
        max_size=int(max_size),
    )
