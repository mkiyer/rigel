"""CalibrationResult — the acyclic calibrator's output schema.

Per-region deconvolved gDNA / RNA mass across the region's three nodes (contained
plus the two boundary sides), the per-node exposure ``ω`` stored explicitly for
QC/diagnostics, the gDNA component's exposure-weighted effective length, and the two
library scalars (``ρ_0``, ``κ_rna``). The calibrator is a single feed-forward pass,
so there are **no** convergence diagnostics. ``__post_init__`` enforces the intrinsic
invariants (shapes, finiteness, sign); mass conservation against the raw fragment
counts is checked by the calibrator / tests (it needs the substrate, which the result
does not carry).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import CalibrationConfig


def _check_region_array(arr: np.ndarray, name: str, n_regions: int) -> None:
    """Validate a per-region float64 array: shape, dtype, finite, non-negative."""
    if not isinstance(arr, np.ndarray):
        raise ValueError(f"CalibrationResult.{name} must be a numpy array.")
    if arr.dtype != np.float64:
        raise ValueError(f"CalibrationResult.{name} must be float64; got {arr.dtype}.")
    if arr.shape != (n_regions,):
        raise ValueError(
            f"CalibrationResult.{name} has shape {arr.shape}; expected ({n_regions},)."
        )
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"CalibrationResult.{name} contains non-finite values.")
    if np.any(arr < 0.0):
        raise ValueError(f"CalibrationResult.{name} must be non-negative.")


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Per-region deconvolved mass + per-node exposure + library scalars (acyclic)."""

    # --- deconvolved mass across the region's 3 nodes (float64[R]) ---
    mass_gdna_contained: np.ndarray
    mass_rna_contained: np.ndarray
    mass_gdna_left: np.ndarray  # right side of the left boundary
    mass_rna_left: np.ndarray
    mass_gdna_right: np.ndarray  # left side of the right boundary
    mass_rna_right: np.ndarray

    # --- per-node exposure ω, stored explicitly for QC / diagnostics (float64[R]) ---
    # ω < 1 depleted, = 1 uniform, > 1 enriched; 0 where a node carries no gDNA.
    exposure_contained: np.ndarray
    exposure_left: np.ndarray
    exposure_right: np.ndarray

    # --- gDNA component geometric effective length: Σ_node L_node (ω-free, Option A) ---
    # Exposure ω is carried in the deconvolved gDNA mass, NOT here — so this length is
    # well-defined even where calibration saw no reads, and never collapses.
    gdna_geom_len: np.ndarray  # float64[R]

    # --- directional boundary effective length 𝓔(L)=E[min(ℓ,L)] per region (geometry) ---
    # The per-region capacity to supply boundary-crossing fragments; consumed by the
    # boundary-flux transport in assemble_priors.
    gdna_side_len: np.ndarray  # float64[R]

    # --- library scalars ---
    gdna_density_global: float  # >= 0, global gDNA density (mass/bp); 0 in a zero-gDNA library
    rna_sense_frac: float  # in [0, 1], RNA sense fraction used by the strand clue

    # --- provenance ---
    n_regions: int
    config: CalibrationConfig

    def __post_init__(self) -> None:
        n = self.n_regions
        if n < 0:
            raise ValueError(f"CalibrationResult.n_regions must be >= 0; got {n}.")

        for name in (
            "mass_gdna_contained",
            "mass_rna_contained",
            "mass_gdna_left",
            "mass_rna_left",
            "mass_gdna_right",
            "mass_rna_right",
            "exposure_contained",
            "exposure_left",
            "exposure_right",
            "gdna_geom_len",
            "gdna_side_len",
        ):
            _check_region_array(getattr(self, name), name, n)

        if not np.isfinite(self.gdna_density_global) or self.gdna_density_global < 0.0:
            raise ValueError(
                f"CalibrationResult.gdna_density_global must be finite and >= 0; got {self.gdna_density_global}."
            )
        sense_frac = float(self.rna_sense_frac)
        if not 0.0 <= sense_frac <= 1.0:
            raise ValueError(
                f"CalibrationResult.rna_sense_frac must be in [0, 1]; got {sense_frac}."
            )


__all__ = ["CalibrationResult"]
