"""Phase 4 — derive the global gDNA density ρ₀ and per-node exposure ω.

A one-shot, feed-forward aggregate of the Phase-3 per-node decode (no loop, no feedback):

    ρ₀ = Σ_nodes gdna_mass / Σ_nodes L_eff          (the GLOBAL average density, decision A)
    ω_node = (gdna_mass / L_eff) / ρ₀               (local density ÷ global; decision A/E — pure ratio)

Exposure is the **pure ratio** of a node's local gDNA density to the library average — `ω<1`
depleted, `=1` uniform, `>1` enriched (capture) — with **no shrinkage** (decision E: the mass
is already regularized by Phase 3's count prior; shrinking ω too would double-count).

Nodes are a region's contained mass (length ``region_eff_len``) and its two boundary **sides**
(length ``boundary_side_eff_len`` each), per the Phase-3.1 unit. A side enters ρ₀ only where its
boundary exists (same reference); otherwise it carries no exposure capacity.

Graceful across the whole spectrum (decision F): zero gDNA ⇒ ρ₀=0 ⇒ ω:=1 (neutral); a pure-RNA
node ⇒ ω≈0 (correctly depleted, not a collapse — other nodes carry the locus's gDNA length).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class DerivedExposure:
    """ρ₀ and the per-node exposure weights (regions + the two boundary sides)."""

    gdna_density_global: float
    region_exposure: np.ndarray  # float64[R]
    left_exposure: np.ndarray  # float64[R]  (right side of each region's left boundary)
    right_exposure: np.ndarray  # float64[R]  (left side of each region's right boundary)
    gdna_geom_len: (
        np.ndarray
    )  # float64[R]  Σ_node L_node (geometric gDNA length, ω-free — Option A)


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
    region_decode,
    left_decode,
    right_decode,
    region_eff_len: np.ndarray,
    boundary_side_eff_len: np.ndarray,
    ref_id: np.ndarray,
) -> DerivedExposure:
    """Derive ρ₀ (global density) and per-node exposure ω from the Phase-3 per-node decode."""
    region_eff_len = np.asarray(region_eff_len, dtype=np.float64)
    bside = np.asarray(boundary_side_eff_len, dtype=np.float64)
    left_exists, right_exists = _side_exists(np.asarray(ref_id))
    # Exposure capacity per node (a side that doesn't exist contributes none).
    left_eff = np.where(left_exists, bside, 0.0)
    right_eff = np.where(right_exists, bside, 0.0)

    total_g = float(
        region_decode.gdna_mass.sum() + left_decode.gdna_mass.sum() + right_decode.gdna_mass.sum()
    )
    total_l = float(region_eff_len.sum() + left_eff.sum() + right_eff.sum())
    gdna_density_global = total_g / total_l if total_l > 0.0 else 0.0

    def _exposure(gdna_mass: np.ndarray, eff: np.ndarray) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            local_density = np.where(eff > 0.0, gdna_mass / np.maximum(eff, 1e-12), 0.0)
        # ω = local/global; neutral (1) where there is no capacity or no global density.
        return np.where(
            (eff > 0.0) & (gdna_density_global > 0.0),
            local_density / max(gdna_density_global, 1e-12),
            1.0,
        )

    region_exposure = _exposure(region_decode.gdna_mass, region_eff_len)
    left_exposure = _exposure(left_decode.gdna_mass, left_eff)
    right_exposure = _exposure(right_decode.gdna_mass, right_eff)
    # Option A: the gDNA component's effective length is the GEOMETRIC genomic span of
    # the region's nodes (ω-free). Exposure ω lives only in the deconvolved gDNA mass,
    # never the length — so the eff-len is knowable even where calibration is blind
    # (multimap / silent loci, ω=0) and can never collapse to the floor. A non-existent
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
