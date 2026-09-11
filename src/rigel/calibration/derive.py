"""Derive the global gDNA density — a library-average QC scalar.

Aggregates the converged deconvolution over both axes into one scalar:

    gdna_density_global = (Σ_regions gdna_mass + Σ_boundaries gdna_mass) / (Σ_regions E_g + Σ_boundaries E_g)

Every entry on both axes is a real object with a real support: a contiguous boundary is the boundary
BETWEEN two adjacent regions, so ``E = N − n_refs`` and no slot has to be masked out.

It is a ratio of SUMS, never a mean of ratios (``ρ_bg = Σg/ΣE``) — a rate pooled over unequal supports
is not the average of the per-object rates.

The per-locus contraction of the gDNA component's effective length under capture is a different
quantity, computed downstream in ``priors.assemble_priors``; this aggregate produces only the
library-average density scalar.
"""

from __future__ import annotations

import numpy as np


def gdna_density_global(
    region_deconv,
    boundary_deconv,
    gdna_region_eff_len: np.ndarray,
    gdna_boundary_eff_len: np.ndarray,
) -> float:
    """Library-average gDNA density (QC scalar) = Σ gDNA mass / Σ gDNA effective length, over both axes.

    ``0.0`` when there is no support anywhere — an empty library has no density, and a floored
    division would report one.
    """
    total_g = float(
        np.asarray(region_deconv.gdna_mass, dtype=np.float64).sum()
        + np.asarray(boundary_deconv.gdna_mass, dtype=np.float64).sum()
    )
    total_l = float(
        np.asarray(gdna_region_eff_len, dtype=np.float64).sum()
        + np.asarray(gdna_boundary_eff_len, dtype=np.float64).sum()
    )
    return total_g / total_l if total_l > 0.0 else 0.0


__all__ = ["gdna_density_global"]
