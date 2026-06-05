"""Phase 4 — derive the global gDNA density (``gdna_density_global``) + per-node exposure.

A one-shot, feed-forward aggregate of the Phase-3 per-node deconvolution (no loop, no feedback):

    gdna_density_global = Σ_nodes gdna_mass / Σ_nodes eff_len   (the GLOBAL average density)
    exposure            = (gdna_mass / eff_len) / gdna_density_global   (local density / global)

Exposure is the **pure ratio** of a node's local gDNA density to the library average — `exposure<1`
depleted, `=1` uniform, `>1` enriched (capture) — with **no shrinkage** (decision E: the mass
is already regularized by Phase 3's count prior; shrinking exposure too would double-count).

Nodes are a region's contained mass (length ``region_eff_len``) and its two boundary **sides**
(length ``boundary_side_eff_len`` each), per the Phase-3.1 unit. A side enters gdna_density_global only where its
boundary exists (same reference); otherwise it carries no exposure capacity.

Graceful across the whole spectrum (decision F): zero gDNA ⇒ gdna_density_global=0 ⇒ exposure:=1 (neutral); a pure-RNA
node ⇒ exposure≈0 (correctly depleted, not a collapse — other nodes carry the locus's gDNA length).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class DerivedExposure:
    """gdna_density_global and the per-node exposure weights (regions + the two boundary sides)."""

    gdna_density_global: float
    region_exposure: np.ndarray  # float64[R]
    left_exposure: np.ndarray  # float64[R]  (right side of each region's left boundary)
    right_exposure: np.ndarray  # float64[R]  (left side of each region's right boundary)
    gdna_geom_len: (
        np.ndarray
    )  # float64[R]  Σ_node L_node (geometric gDNA length, exposure-free — Option A)


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
) -> DerivedExposure:
    """Derive gdna_density_global (global density) and per-node exposure from the Phase-3 per-node deconvolution."""
    region_eff_len = np.asarray(region_eff_len, dtype=np.float64)
    bside = np.asarray(boundary_side_eff_len, dtype=np.float64)
    left_exists, right_exists = _side_exists(np.asarray(ref_id))
    # Exposure capacity per node (a side that doesn't exist contributes none).
    left_eff = np.where(left_exists, bside, 0.0)
    right_eff = np.where(right_exists, bside, 0.0)

    total_g = float(
        region_deconv.gdna_mass.sum() + left_deconv.gdna_mass.sum() + right_deconv.gdna_mass.sum()
    )
    total_l = float(region_eff_len.sum() + left_eff.sum() + right_eff.sum())
    gdna_density_global = total_g / total_l if total_l > 0.0 else 0.0

    def _exposure(gdna_mass: np.ndarray, eff: np.ndarray) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            local_density = np.where(eff > 0.0, gdna_mass / np.maximum(eff, 1e-12), 0.0)
        # exposure = local/global; neutral (1) where there is no capacity or no global density.
        return np.where(
            (eff > 0.0) & (gdna_density_global > 0.0),
            local_density / max(gdna_density_global, 1e-12),
            1.0,
        )

    region_exposure = _exposure(region_deconv.gdna_mass, region_eff_len)
    left_exposure = _exposure(left_deconv.gdna_mass, left_eff)
    right_exposure = _exposure(right_deconv.gdna_mass, right_eff)
    # Option A: the gDNA component's effective length is the GEOMETRIC genomic span of
    # the region's nodes (exposure-free). Exposure exposure lives only in the deconvolved gDNA mass,
    # never the length — so the eff-len is knowable even where calibration is blind
    # (multimap / silent loci, exposure=0) and can never collapse to the floor. A non-existent
    # boundary side has left_eff/right_eff = 0 and drops out. See
    # docs/futureprs/phase6_multimap_regression_analysis.md.
    gdna_geom_len = region_eff_len + left_eff + right_eff
    return DerivedExposure(
        gdna_density_global=gdna_density_global,
        region_exposure=region_exposure,
        left_exposure=left_exposure,
        right_exposure=right_exposure,
        gdna_geom_len=gdna_geom_len,
    )


__all__ = ["DerivedExposure", "derive"]
