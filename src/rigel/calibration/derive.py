"""Phase 4 — derive the global gDNA density (a QC scalar).

A one-shot, feed-forward aggregate of the Phase-3 per-node deconvolution (no loop, no feedback):

    gdna_density_global = Σ_nodes gdna_mass / Σ_nodes eff_len   (library-average density — a QC scalar)

Nodes are a region's contained mass (length ``region_eff_len``) and its two boundary **sides**
(length ``boundary_side_eff_len`` each); a side enters the totals only where its boundary exists
(same reference). The per-region **contraction** of the gDNA component's effective length under
capture is the inverse participation ratio (IPR) of the per-region gDNA *mass*, computed downstream
in ``priors.assemble_priors`` from the per-region effective supports (``gdna_region_eff_len`` +
``gdna_boundary_len``) directly — this aggregate only produces the library-average density scalar.
"""

from __future__ import annotations

import numpy as np

from .run_fill import same_ref_left_right


def gdna_density_global(
    region_deconv,
    left_deconv,
    right_deconv,
    region_eff_len: np.ndarray,
    boundary_side_eff_len: np.ndarray,
    ref_id: np.ndarray,
) -> float:
    """Library-average gDNA density (QC scalar) = Σ node gDNA mass / Σ node effective length."""
    region_eff_len = np.asarray(region_eff_len, dtype=np.float64)
    bside = np.asarray(boundary_side_eff_len, dtype=np.float64)
    left_exists, right_exists = same_ref_left_right(np.asarray(ref_id))
    # A side that doesn't exist contributes no length.
    left_eff = np.where(left_exists, bside, 0.0)
    right_eff = np.where(right_exists, bside, 0.0)

    total_g = float(
        region_deconv.gdna_mass.sum() + left_deconv.gdna_mass.sum() + right_deconv.gdna_mass.sum()
    )
    total_l = float(region_eff_len.sum() + left_eff.sum() + right_eff.sum())
    return total_g / total_l if total_l > 0.0 else 0.0


__all__ = ["gdna_density_global"]
