"""Derive the global gDNA density — a library-average QC scalar.

An aggregate of the converged deconvolution into the library-average gDNA density scalar:

    gdna_density_global = (Σ_regions gdna_mass + Σ_edges gdna_mass) / (Σ_regions E_g + Σ_edges E_g)

⭐ **Two axes, and every object on them exists.** The predecessor summed a region's contained mass plus
its two boundary SIDES, and had to mask each side with ``same_ref_left_right`` because a reference
terminal's outer boundary has nothing on the far side — "a side that doesn't exist contributes no
length". A contiguous edge is the line BETWEEN two adjacent regions, so there is no such object to
exclude: ``E = N − n_refs`` and every entry is real. The mask goes with the terminal slots.

⚠ **It is a ratio of SUMS, never a mean of ratios** (``ρ_bg = Σg/ΣE``) — a rate
pooled over unequal supports is not the average of the per-object rates.

The per-locus **contraction** of the gDNA component's effective length under capture is the inverse
participation ratio of the deconvolved gDNA mass over the two supports, computed downstream in
``priors.assemble_priors``; this aggregate only produces the library-average density scalar.
"""

from __future__ import annotations

import numpy as np


def gdna_density_global(
    region_deconv,
    edge_deconv,
    gdna_region_eff_len: np.ndarray,
    gdna_edge_eff_len: np.ndarray,
) -> float:
    """Library-average gDNA density (QC scalar) = Σ gDNA mass / Σ gDNA effective length, over both axes.

    ``0.0`` when there is no support anywhere — an empty library has no density, and a floored
    division would report one.
    """
    total_g = float(
        np.asarray(region_deconv.gdna_mass, dtype=np.float64).sum()
        + np.asarray(edge_deconv.gdna_mass, dtype=np.float64).sum()
    )
    total_l = float(
        np.asarray(gdna_region_eff_len, dtype=np.float64).sum()
        + np.asarray(gdna_edge_eff_len, dtype=np.float64).sum()
    )
    return total_g / total_l if total_l > 0.0 else 0.0


__all__ = ["gdna_density_global"]
