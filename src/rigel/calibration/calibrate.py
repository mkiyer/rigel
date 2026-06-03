"""calibrate() — the acyclic fractional-accumulator calibrator.

A single feed-forward pass: no EM loop, no ``ρ_0 → decode → ρ_0`` feedback. It
deconvolves each node (every region and the two boundary sides) into gDNA / RNA mass
from two **orthogonal** clues — a count-based gDNA density and the strand split — then
**derives** the global gDNA density ``ρ_0`` and the per-node exposure ``ω`` as outputs::

    substrate → node_gdna_density → decode_regions / decode_sides → derive
                  (count clue)        (joint count × strand)        (ρ_0, ω, eff-len)

Decoding each node from LOCAL evidence only — and deriving ``ρ_0`` / ``ω`` from the
aggregate afterward — is what dissolved the old EM loop's sparse-data collapse (the
loop re-estimated ``ρ_0`` from decisions it had itself driven). See
``docs/futureprs/acyclic_deconvolution_design.md``.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from .density_model import node_gdna_density
from .derive import derive
from .effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from .joint_decode import decode_regions, decode_sides
from .result import CalibrationResult
from .strand_balance import fit_strand_balance
from .substrate import CalibrationSubstrate

if TYPE_CHECKING:
    import numpy as np

    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays

logger = logging.getLogger(__name__)


def calibrate(
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    gdna_fl_pmf: "np.ndarray",
    config: "CalibrationConfig",
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per node, then derive ρ_0 and exposure.

    Single feed-forward pass; see the module docstring for the data flow. ``ρ_0`` may
    be ``0`` (a zero-gDNA library) and per-node ``ω`` may be ``0`` (a pure-RNA node);
    both are valid, graceful outputs — not failures.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)

    # gDNA fragment-length effective lengths (PR 4c geometry): the region-contained
    # length, the per-side boundary density length, and the region-free crossing mean.
    region_eff = region_eff_length(region_arrays.region_size_bp, gdna_fl_pmf)
    bside_eff = boundary_side_eff_length(gdna_fl_pmf, region_arrays.region_size_bp)
    mu_fl = boundary_eff_length(gdna_fl_pmf)

    # Count clue (per-node gDNA density via the region↔boundary sweep) + RNA strand
    # balance (κ_rna posterior mean from the spliced channel).
    node_density = node_gdna_density(substrate, region_arrays, region_eff, mu_fl)
    kappa_rna = float(fit_strand_balance(strand_model).kappa_rna)

    # Joint per-node deconvolution: count prior × Beta-Binomial strand likelihood.
    regions = decode_regions(
        substrate,
        region_arrays,
        node_density,
        region_eff,
        kappa_rna=kappa_rna,
        confidence=config.confidence,
        n_grid=config.n_grid,
    )
    left, right = decode_sides(
        substrate,
        region_arrays,
        node_density,
        bside_eff,
        kappa_rna=kappa_rna,
        confidence=config.confidence,
        n_grid=config.n_grid,
    )

    # Derive ρ_0 and per-node exposure ω (+ the exposure-weighted gDNA length).
    derived = derive(regions, left, right, region_eff, bside_eff, region_arrays.ref_id)

    result = CalibrationResult(
        mass_g_contained=regions.gdna_mass,
        mass_d_contained=regions.rna_mass,
        mass_g_left=left.gdna_mass,
        mass_d_left=left.rna_mass,
        mass_g_right=right.gdna_mass,
        mass_d_right=right.rna_mass,
        omega_contained=derived.region_omega,
        omega_left=derived.left_omega,
        omega_right=derived.right_omega,
        gdna_exposure_len=derived.gdna_exposure_len,
        rho_0=derived.rho_0,
        kappa_rna=kappa_rna,
        n_regions=region_arrays.n_regions,
        config=config,
    )
    logger.debug(
        "calibration: R=%d rho_0=%.4g kappa_rna=%.3f", result.n_regions, result.rho_0, kappa_rna
    )
    return result


__all__ = ["calibrate"]
