"""calibrate() — the acyclic fractional-accumulator calibrator.

A single feed-forward pass: no EM loop, no ``ρ_0 → decode → ρ_0`` feedback. The strand clue
(``κ_rna``) is fit first and used to **strand-clean** the count density, so the global gDNA
density ``ρ_0`` is read from clean gDNA, not gDNA + nascent RNA. Each node (every region and the
two boundary sides) is then deconvolved into gDNA / RNA by the **joint** count × strand
posterior; ``ρ_0`` and the per-node exposure ``ω`` are **derived** from the aggregate::

    substrate
      → strand balance (κ_rna) → strand-clean the count density (closed-form gDNA fraction)
      → node_gdna_density (count clue; density strand-cleaned by κ_rna)
      → decode_regions / decode_sides (joint count × strand)
      → derive (ρ_0, ω, eff-len)

The strand clue cleans ``ρ_0`` (a global hyperparameter) and also enters each node's joint
decode; the two are conditionally independent given ``(ρ_0, κ_rna)`` and a node's ~1/N
self-contribution to ``ρ_0`` is negligible (empirical Bayes). Deriving ``ρ_0`` / ``ω`` from the
aggregate — rather than looping — is what dissolved the old EM loop's sparse-data collapse. See
``docs/futureprs/acyclic_deconvolution_design.md``.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from .density_model import node_gdna_density
from .derive import derive
from .errors import CalibrationStrandError
from .effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from .joint_decode import decode_regions, decode_sides
from .result import CalibrationResult
from .signature import TS_NEG, TS_NONE
from .strand_balance import fit_strand_balance
from .substrate import CalibrationSubstrate

if TYPE_CHECKING:
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

    # RNA strand balance first (κ_rna posterior mean from the spliced channel) — the strand
    # clue is what cleans the count density, so it must precede Phase 1.
    balance = fit_strand_balance(strand_model)
    if balance.fallback_used:
        # No spliced reads ⇒ κ_rna degenerates to 0.5 (symmetric), indistinguishable
        # from gDNA. Fail loudly rather than mis-deconvolve (see CalibrationStrandError).
        raise CalibrationStrandError(
            "the library has zero spliced unique-mapper observations; the RNA strand "
            "orientation cannot be identified, so gDNA and sense RNA cannot be separated. "
            "A real RNA-seq library always carries spliced reads."
        )
    kappa_rna = float(balance.kappa_rna)

    # Strand-deconvolved gDNA fraction of each region's contained mass, used to clean the
    # nascent-RNA bias from the count density (ρ_0). Closed-form (unbiased) MLE from the oriented
    # sense fraction: π_g = (sense_frac − κ_rna)/(½ − κ_rna), clamped to [0, 1] — gDNA is
    # symmetric (½), RNA skewed (κ_rna), so a node at the RNA rate reads π_g = 0 (no inflation),
    # and any symmetric excess reads as gDNA. (The grid posterior MEAN would bias π_g toward ½
    # and under-clean.) Intergenic (NONE, no transcript) is pure gDNA; unstranded (κ_rna = ½)
    # cannot strand-clean, so ρ_0 falls back to the raw count density. Exonic / AMBIG regions are
    # not count-decodable, so their fraction never enters ρ_0.
    ts = np.asarray(region_arrays.ts_class)
    c = substrate.contained
    n_unspl = (c.n_unspliced_pos + c.n_unspliced_neg).astype(np.float64)
    sense = np.where(ts == TS_NEG, c.n_unspliced_neg, c.n_unspliced_pos).astype(np.float64)
    sense_frac = np.where(n_unspl > 0.0, sense / np.maximum(n_unspl, 1e-9), 0.5)
    denom = 0.5 - kappa_rna
    if abs(denom) < 1.0e-6:  # unstranded — strand cannot clean; keep the raw count density
        gdna_frac = np.ones(region_arrays.n_regions, dtype=np.float64)
    else:
        gdna_frac = np.clip((sense_frac - kappa_rna) / denom, 0.0, 1.0)
    gdna_frac = np.where(ts == TS_NONE, 1.0, gdna_frac)  # intergenic = pure gDNA

    # Count clue: per-node gDNA density via the region↔boundary sweep, on the strand-cleaned
    # density (ρ_0 is now clean gDNA, not gDNA+nascent). The count-prior precision is the
    # expected gDNA count μ_g = density·eff_len (NodeDensity.count_evidence), so the count clue
    # defers to strand where RNA dominates.
    node_density = node_gdna_density(
        substrate, region_arrays, region_eff, mu_fl, gdna_frac=gdna_frac
    )

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
        gdna_geom_len=derived.gdna_geom_len,
        gdna_side_len=bside_eff,
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
