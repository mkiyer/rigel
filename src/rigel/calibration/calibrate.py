"""calibrate() — the acyclic fractional-accumulator calibrator (decoupled strand/count).

A single feed-forward pass: no EM loop, no density->deconv->density feedback. Each node (every region
and the two boundary sides) is deconvolved into gDNA / RNA by a **strand-vs-count handoff**, not a
joint product: a strand-observable node in a strand-identifiable library uses the **strand module**
(Beta-Binomial posterior), every other node uses the **count module** (raw density ratio). The two
estimators act on disjoint nodes and never multiply — the decoupling that fixed the joint model's
bias-mixing (see ``docs/calibration/decoupled_calibration_design.md``; joint archived in
``docs/calibration/archive/joint_deconvolution.md``)::

    substrate
      -> strand balance: rna_sense_frac (κ) + strand_identifiable (the global routing gate)
      -> node_gdna_density (count module; RAW density, intergenic-anchored global, no strand clean)
      -> gdna/rna strand Beta-Binomial overdispersions (strand module parameters)
      -> deconv_regions / deconv_sides (per-node handoff: strand posterior OR count fraction)
      -> derive (gdna_density_global, gdna_geom_len)

The global density and geometric gDNA length are **derived** from the aggregate deconvolved mass --
not looped -- which dissolved the old EM loop's sparse-data collapse. The per-region gDNA length
contraction under capture is the IPR of the deconvolved gDNA mass, applied later in
``priors.assemble_priors``.
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
from .gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_substrate,
    overdispersion_for_beta,
)
from .result import CalibrationResult
from .strand_balance import fit_strand_balance
from .strand_deconv import deconv_regions, deconv_sides
from .strand_summary import strand_contrast_identifiable
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
    """Deconvolve the library into gDNA / RNA per node, then derive gdna_density_global.

    Single feed-forward pass; see the module docstring for the data flow. ``gdna_density_global``
    may be ``0`` (a zero-gDNA library) and a node's deconvolved gDNA mass may be ``0`` (a pure-RNA
    node); both are valid, graceful outputs — not failures.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)

    # gDNA fragment-length effective lengths (PR 4c geometry): the region-contained
    # length, the per-side boundary density length, and the region-free crossing mean.
    region_eff_len = region_eff_length(region_arrays.region_size_bp, gdna_fl_pmf)
    boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, region_arrays.region_size_bp)
    fl_mean = boundary_eff_length(gdna_fl_pmf)

    # RNA strand balance: rna_sense_frac (κ) = posterior-mean spliced sense fraction, and whether the
    # strand channel is IDENTIFIABLE (κ distinguishable from ½ at 99% given its sample size). The
    # identifiability is the GLOBAL routing gate: a strand-identifiable library routes its
    # strand-observable nodes to the strand module, an unstranded one routes everything to count.
    balance = fit_strand_balance(strand_model)
    if balance.fallback_used:
        # No spliced reads at all — not a usable RNA-seq library. Fail loudly (a real RNA-seq library
        # always carries spliced reads); see CalibrationStrandError.
        raise CalibrationStrandError(
            "the library has zero spliced unique-mapper observations; this does not look like an "
            "RNA-seq library. A real RNA-seq library always carries spliced reads."
        )
    rna_sense_frac = float(balance.rna_sense_frac)
    strand_identifiable = strand_contrast_identifiable(
        strand_model.p_r1_sense,
        strand_model.n_observations,
        confidence=config.strand_identifiability_confidence,
    )

    # Count clue (the count module): per-region gDNA density by LOCAL boundary-anchored imputation on
    # RAW unspliced counts (no strand cleaning — the strand module owns the strand channel). An
    # observable region uses its own contained density; a non-observable (exon/AMBIG) region is
    # anchored from its observable boundary sides; the global fallback comes from intergenic regions.
    # node_gdna_density returns the count module's gDNA fraction (count_gdna_frac) per node, consumed
    # for count-routed nodes and as the gDNA strand-fit seed weight.
    node_density = node_gdna_density(substrate, region_arrays, region_eff_len, fl_mean)

    # gDNA strand Beta-Binomial overdispersion: fitted from the count-observable seed regions
    # (intergenic + intronic), using the count-clue gDNA weight (overdispersion-free, since the
    # density is cleaned by the strand MEAN ½, not the dispersion) to break the circularity.
    # gDNA is unstranded (mean ½) but overdispersed in real libraries; this is what keeps a
    # noise-skewed pure-gDNA node reading as gDNA rather than mis-called RNA. The RNA strand is
    # given a symmetric overdispersion just below (no longer Binomial). See docs/em_strand/03+05.
    gdna_strand = fit_gdna_strand_from_substrate(
        substrate,
        region_arrays,
        node_density,
        boundary_eff_len,
        rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(config.gdna_strand_prior_alpha_beta),
        prior_weight=config.gdna_strand_prior_weight,
    )
    gdna_strand_overdispersion = gdna_strand.gdna_strand_overdispersion

    # RNA strand Beta-Binomial overdispersion: fitted from boundary-side SPLICED counts (spliced ⇒
    # pure RNA, so node mean = κ and the whole node is RNA). Symmetric with the gDNA fit and shrunk
    # toward the *same* default prior, so under sparse data both components collapse to one
    # distribution and an unstranded node (κ = ½) is uninformative — the symmetry an earlier
    # gDNA-only overdispersion broke. See docs/em_strand/05 + gdna_strand.py.
    rna_strand = fit_rna_strand_from_substrate(
        substrate,
        rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(config.rna_strand_prior_alpha_beta),
        prior_weight=config.rna_strand_prior_weight,
    )
    rna_strand_overdispersion = rna_strand.rna_strand_overdispersion

    # Per-node deconvolution by the strand/count HANDOFF: each node routes to the strand module
    # (Beta-Binomial posterior) when the library is strand-identifiable AND the node is
    # strand-observable, else to the count module (count_gdna_frac). The two never multiply.
    regions = deconv_regions(
        substrate,
        region_arrays,
        node_density,
        strand_identifiable=strand_identifiable,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        gdna_strand_llr_bias=config.gdna_strand_llr_bias,
        n_grid=config.n_grid,
    )
    left, right = deconv_sides(
        substrate,
        region_arrays,
        node_density,
        boundary_eff_len,
        strand_identifiable=strand_identifiable,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        gdna_strand_llr_bias=config.gdna_strand_llr_bias,
        n_grid=config.n_grid,
    )

    # Derive gdna_density_global (QC scalar) and the geometric gDNA length.
    derived = derive(regions, left, right, region_eff_len, boundary_eff_len, region_arrays.ref_id)

    result = CalibrationResult(
        mass_gdna_contained=regions.gdna_mass,
        mass_rna_contained=regions.rna_mass,
        mass_gdna_left=left.gdna_mass,
        mass_rna_left=left.rna_mass,
        mass_gdna_right=right.gdna_mass,
        mass_rna_right=right.rna_mass,
        gdna_geom_len=derived.gdna_geom_len,
        gdna_boundary_len=boundary_eff_len,
        gdna_density_global=derived.gdna_density_global,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_regions=region_arrays.n_regions,
        config=config,
    )
    # Diagnostic: the boundary-spliced sense fraction should agree with the StrandModel κ (the
    # decode mean). A large gap flags a strand-model / accumulator mismatch (we do NOT refit the
    # mean from boundary spliced — κ stays the StrandModel posterior; this is QC only).
    spl_sense = float(
        np.sum(substrate.left.n_spliced_sense) + np.sum(substrate.right.n_spliced_sense)
    )
    spl_total = spl_sense + float(
        np.sum(substrate.left.n_spliced_antisense) + np.sum(substrate.right.n_spliced_antisense)
    )
    boundary_sense_frac = spl_sense / spl_total if spl_total > 0.0 else float("nan")
    logger.debug(
        "calibration: R=%d gdna_density_global=%.4g rna_sense_frac=%.3f "
        "gdna_strand_overdispersion=%.4g (%d seed nodes, %d frags%s) "
        "rna_strand_overdispersion=%.4g (%d seed sides, %d frags%s) "
        "[boundary-spliced sense_frac=%.3f vs κ=%.3f]",
        result.n_regions,
        result.gdna_density_global,
        rna_sense_frac,
        gdna_strand_overdispersion,
        gdna_strand.n_seed_nodes,
        gdna_strand.n_seed_fragments,
        ", FALLBACK" if gdna_strand.fallback_used else "",
        rna_strand_overdispersion,
        rna_strand.n_seed_nodes,
        rna_strand.n_seed_fragments,
        ", FALLBACK" if rna_strand.fallback_used else "",
        boundary_sense_frac,
        rna_sense_frac,
    )
    return result


__all__ = ["calibrate"]
