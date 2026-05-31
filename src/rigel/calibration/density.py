"""Global gDNA density seed — the ρ_0 Gamma posterior (doc 03 §5.1 / §7).

``ρ_0`` is the library-wide gDNA fragment density (fragments per bp of
accessible gDNA). It is the count channel's gDNA-only rate
(:mod:`rigel.calibration.estep`) and the exposure denominator
(:mod:`rigel.calibration.exposure`). This module seeds it.

**Seed set (PR04a §I.1 / III.2).** Pool the regions whose unspliced mass is
gDNA-dominated *by annotation* — **intergenic** (no transcript) **and intronic**
(no exon) regions — together with the boundary mass crossing into them. Exonic
regions are excluded (their unspliced mass mixes mature/nascent RNA with gDNA).

Seeding from intergenic regions *alone* fails on hybrid-capture libraries, where
intergenic gDNA is depleted relative to the captured (genic) gDNA; intronic
regions and boundary-crossing fragments near captured exons carry the
representative gDNA signal. Intronic unspliced mass also contains nascent RNA,
so this seed is upper-biased where introns are transcribed — an explicit,
benchmarked approximation, with the expressed/unexpressed latent as the planned
remedy (future PR). **No density clip** is applied: the burned
``top_t_fraction`` percentile cliff is gone.

The estimate is the ρ_0 M-step (doc 03 §5.1) restricted to the seed regions with
exposure ω ≡ 1, regularized by a unit-strength Gamma prior::

    ρ_0 = (α₀ + Σ_seed M_g_tot) / (β₀ + Σ_seed L_eff),    α₀ = β₀ = 1,

where ``M_g_tot = contained + left + right`` unspliced mass (the same D1
side-attribution the exposure uses; all of it called gDNA in the seed regions).
With no seed region (no non-exonic annotation — degenerate), it falls back to
the doc 03 §7 half-and-half density over all regions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .signature import BIT_EXON_NEG, BIT_EXON_POS

if TYPE_CHECKING:
    from .region_arrays import RegionArrays
    from .substrate import CalibrationSubstrate

# Unit-strength Gamma prior on the gDNA density: one pseudo-fragment over one
# pseudo-base. Weakly informative (keeps ρ_0 > 0 with no data); not a tunable.
# Spec: doc 03 §8 — "each numeric is a unit-strength Bayesian prior".
_RHO_PRIOR_ALPHA = 1.0
_RHO_PRIOR_BETA = 1.0


@dataclass(frozen=True, slots=True)
class GdnaDensity:
    """Gamma posterior for the library-wide gDNA density ρ_0."""

    rho_0: float  # posterior mean α/β (fragments per bp), > 0
    alpha: float  # posterior shape
    beta: float  # posterior rate
    n_seed_regions: int  # intergenic + intronic regions pooled
    seed_mass: float  # Σ_seed gDNA mass (numerator contribution, excl. prior)
    seed_length: float  # Σ_seed L_eff (denominator contribution, excl. prior)
    fallback_used: bool  # True → half-and-half over all regions (no seed region)


def _non_exonic(signature: np.ndarray) -> np.ndarray:
    """Boolean mask of regions with no exon bit — intergenic or intronic."""
    has_exon = (np.asarray(signature) & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    return ~has_exon


def estimate_gdna_density(
    substrate: "CalibrationSubstrate",
    region_arrays: "RegionArrays",
) -> GdnaDensity:
    """Seed ρ_0 from the gDNA-dominated (non-exonic) regions + their boundaries."""
    L_eff = np.asarray(substrate.L_eff, dtype=np.float64)
    # Total unspliced mass attributed to each region (D1: contained + both
    # boundary sides) — treated as gDNA in the seed regions.
    gdna_mass = (
        substrate.contained.mass_unspliced
        + substrate.left.mass_unspliced
        + substrate.right.mass_unspliced
    )

    seed = _non_exonic(region_arrays.signature)
    seed_length = float(L_eff[seed].sum())

    if seed_length > 0.0:
        seed_mass = float(gdna_mass[seed].sum())
        alpha = _RHO_PRIOR_ALPHA + seed_mass
        beta = _RHO_PRIOR_BETA + seed_length
        return GdnaDensity(
            rho_0=alpha / beta,
            alpha=alpha,
            beta=beta,
            n_seed_regions=int(seed.sum()),
            seed_mass=seed_mass,
            seed_length=seed_length,
            fallback_used=False,
        )

    # Degenerate: no non-exonic region anywhere. Fall back to the doc 03 §7
    # half-and-half density over all regions, expressed as the same Gamma
    # posterior (half the total mass is assumed gDNA). The unit prior floors it.
    half_mass = 0.5 * float(gdna_mass.sum())
    total_length = float(L_eff.sum())
    alpha = _RHO_PRIOR_ALPHA + half_mass
    beta = _RHO_PRIOR_BETA + total_length
    return GdnaDensity(
        rho_0=alpha / beta,
        alpha=alpha,
        beta=beta,
        n_seed_regions=0,
        seed_mass=0.0,
        seed_length=0.0,
        fallback_used=True,
    )


__all__ = ["GdnaDensity", "estimate_gdna_density"]
