"""calibrate() — the acyclic fractional-accumulator calibrator.

A single feed-forward pass: no EM loop, no density->deconv->density feedback. The strand clue
(``rna_sense_frac``) is fit first and used to **strand-clean** the count density, so the global
gDNA density (``gdna_density_global``) is read from clean gDNA, not gDNA + nascent RNA. Each node
(every region and the two boundary sides) is then deconvolved into gDNA / RNA by the **joint**
count x strand posterior; the global density and the per-node ``exposure`` are **derived** from
the aggregate::

    substrate
      -> strand balance (rna_sense_frac) -> strand-clean the count density (closed-form gDNA frac)
      -> node_gdna_density (count clue; density strand-cleaned by rna_sense_frac)
      -> deconv_regions / deconv_sides (joint count x strand)
      -> derive (gdna_density_global, exposure, eff-len)

The strand clue cleans the global density (a hyperparameter) and also enters each node's joint
deconv; the two are conditionally independent given ``(gdna_density_global, rna_sense_frac)`` and a
node's ~1/N self-contribution to the global density is negligible (empirical Bayes). Deriving the
global density and ``exposure`` from the aggregate -- not looping -- is what dissolved the old EM
loop's sparse-data collapse. See ``docs/futureprs/acyclic_deconvolution_design.md``.
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
from .joint_deconv import deconv_regions, deconv_sides
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
    """Deconvolve the library into gDNA / RNA per node, then derive gdna_density_global and exposure.

    Single feed-forward pass; see the module docstring for the data flow. ``gdna_density_global`` may
    be ``0`` (a zero-gDNA library) and per-node ``exposure`` may be ``0`` (a pure-RNA node);
    both are valid, graceful outputs — not failures.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)

    # gDNA fragment-length effective lengths (PR 4c geometry): the region-contained
    # length, the per-side boundary density length, and the region-free crossing mean.
    region_eff_len = region_eff_length(region_arrays.region_size_bp, gdna_fl_pmf)
    boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, region_arrays.region_size_bp)
    fl_mean = boundary_eff_length(gdna_fl_pmf)

    # RNA strand balance first (rna_sense_frac posterior mean from the spliced channel) — the strand
    # clue is what cleans the count density, so it must precede Phase 1.
    balance = fit_strand_balance(strand_model)
    if balance.fallback_used:
        # No spliced reads ⇒ rna_sense_frac degenerates to 0.5 (symmetric), indistinguishable
        # from gDNA. Fail loudly rather than mis-deconvolve (see CalibrationStrandError).
        raise CalibrationStrandError(
            "the library has zero spliced unique-mapper observations; the RNA strand "
            "orientation cannot be identified, so gDNA and sense RNA cannot be separated. "
            "A real RNA-seq library always carries spliced reads."
        )
    rna_sense_frac = float(balance.rna_sense_frac)

    # Strand-deconvolved gDNA fraction of each region's contained mass, used to clean the
    # nascent-RNA bias from the count density (gdna_density_global). Closed-form (unbiased) MLE from the oriented
    # sense fraction: gdna_frac = (sense_frac − rna_sense_frac)/(½ − rna_sense_frac), clamped to [0, 1] — gDNA is
    # symmetric (½), RNA skewed (rna_sense_frac), so a node at the RNA rate reads gdna_frac = 0 (no inflation),
    # and any symmetric excess reads as gDNA. (The grid posterior MEAN would bias gdna_frac toward ½
    # and under-clean.) Intergenic (NONE, no transcript) is pure gDNA; unstranded (rna_sense_frac = ½)
    # cannot strand-clean, so gdna_density_global falls back to the raw count density. Exonic / AMBIG regions are
    # not count-observable, so their fraction never enters gdna_density_global.
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    n_unspl = (c.n_unspliced_pos + c.n_unspliced_neg).astype(np.float64)
    sense = np.where(ts == TS_NEG, c.n_unspliced_neg, c.n_unspliced_pos).astype(np.float64)
    sense_frac = np.where(n_unspl > 0.0, sense / np.maximum(n_unspl, 1e-9), 0.5)
    denom = 0.5 - rna_sense_frac
    if abs(denom) < 1.0e-6:  # unstranded — strand cannot clean; keep the raw count density
        gdna_frac = np.ones(region_arrays.n_regions, dtype=np.float64)
    else:
        gdna_frac = np.clip((sense_frac - rna_sense_frac) / denom, 0.0, 1.0)
    gdna_frac = np.where(ts == TS_NONE, 1.0, gdna_frac)  # intergenic = pure gDNA

    # Count clue: per-node gDNA density via the region↔boundary sweep, on the strand-cleaned
    # density (gdna_density_global is now clean gDNA, not gDNA+nascent). The count-prior precision is the
    # expected gDNA count count_evidence = density·eff_len (NodeDensity.count_evidence), so the count clue
    # defers to strand where RNA dominates.
    node_density = node_gdna_density(
        substrate, region_arrays, region_eff_len, fl_mean, gdna_frac=gdna_frac
    )

    # Joint per-node deconvolution: count prior × Beta-Binomial strand likelihood.
    regions = deconv_regions(
        substrate,
        region_arrays,
        node_density,
        region_eff_len,
        rna_sense_frac=rna_sense_frac,
        confidence=config.confidence,
        n_grid=config.n_grid,
    )
    left, right = deconv_sides(
        substrate,
        region_arrays,
        node_density,
        boundary_eff_len,
        rna_sense_frac=rna_sense_frac,
        confidence=config.confidence,
        n_grid=config.n_grid,
    )

    # Derive gdna_density_global and per-node exposure (+ the exposure-weighted gDNA length).
    derived = derive(regions, left, right, region_eff_len, boundary_eff_len, region_arrays.ref_id)

    result = CalibrationResult(
        mass_gdna_contained=regions.gdna_mass,
        mass_rna_contained=regions.rna_mass,
        mass_gdna_left=left.gdna_mass,
        mass_rna_left=left.rna_mass,
        mass_gdna_right=right.gdna_mass,
        mass_rna_right=right.rna_mass,
        exposure_contained=derived.region_exposure,
        exposure_left=derived.left_exposure,
        exposure_right=derived.right_exposure,
        gdna_geom_len=derived.gdna_geom_len,
        gdna_boundary_len=boundary_eff_len,
        gdna_density_global=derived.gdna_density_global,
        rna_sense_frac=rna_sense_frac,
        n_regions=region_arrays.n_regions,
        config=config,
    )
    logger.debug(
        "calibration: R=%d gdna_density_global=%.4g rna_sense_frac=%.3f",
        result.n_regions,
        result.gdna_density_global,
        rna_sense_frac,
    )
    return result


__all__ = ["calibrate"]
