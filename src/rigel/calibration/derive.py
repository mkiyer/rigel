"""Phase 4 — derive the global gDNA density and the geometric gDNA length.

A one-shot, feed-forward aggregate of the Phase-3 per-node deconvolution (no loop, no feedback):

    gdna_density_global = Σ_nodes gdna_mass / Σ_nodes eff_len   (library-average density — a QC scalar)
    gdna_geom_len       = region_eff_len + left_eff + right_eff (the GEOMETRIC gDNA support per region)

Nodes are a region's contained mass (length ``region_eff_len``) and its two boundary **sides**
(length ``boundary_side_eff_len`` each); a side enters the totals only where its boundary exists
(same reference). The per-region **contraction** of the gDNA component's effective length under
capture is the inverse participation ratio (IPR) of the per-region gDNA *mass*, computed downstream
in ``priors.assemble_priors`` from ``gdna_geom_len`` — the length itself stays geometric here, so it
is well-defined even where calibration saw no reads.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class DerivedDensity:
    """The global gDNA density (QC scalar) and the per-region geometric gDNA length."""

    gdna_density_global: float
    gdna_geom_len: np.ndarray  # float64[R]  Σ_node L_node (geometric, capture-independent)


def _side_exists(ref_id: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Boolean masks: does region r have a left / right internal (same-ref) boundary?"""
    r = ref_id.shape[0]
    left = np.zeros(r, dtype=bool)
    right = np.zeros(r, dtype=bool)
    if r > 1:
        same = ref_id[:-1] == ref_id[1:]
        left[1:] = same
        right[:-1] = same
    return left, right


def derive(
    region_deconv,
    left_deconv,
    right_deconv,
    region_eff_len: np.ndarray,
    boundary_side_eff_len: np.ndarray,
    ref_id: np.ndarray,
) -> DerivedDensity:
    """Aggregate the Phase-3 deconvolution into the global density + geometric gDNA length."""
    region_eff_len = np.asarray(region_eff_len, dtype=np.float64)
    bside = np.asarray(boundary_side_eff_len, dtype=np.float64)
    left_exists, right_exists = _side_exists(np.asarray(ref_id))
    # A side that doesn't exist contributes no length.
    left_eff = np.where(left_exists, bside, 0.0)
    right_eff = np.where(right_exists, bside, 0.0)

    total_g = float(
        region_deconv.gdna_mass.sum() + left_deconv.gdna_mass.sum() + right_deconv.gdna_mass.sum()
    )
    total_l = float(region_eff_len.sum() + left_eff.sum() + right_eff.sum())
    gdna_density_global = total_g / total_l if total_l > 0.0 else 0.0

    gdna_geom_len = region_eff_len + left_eff + right_eff
    return DerivedDensity(gdna_density_global=gdna_density_global, gdna_geom_len=gdna_geom_len)


__all__ = ["DerivedDensity", "derive"]
