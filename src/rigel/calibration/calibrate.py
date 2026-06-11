"""calibrate() — the acyclic fractional-accumulator calibrator (decoupled strand/count).

A single feed-forward pass: no EM loop, no density->deconv->density feedback. Each node (every region
and the two boundary sides) is deconvolved into gDNA / RNA by **precision-weighted strand→count
deference**: ``g = w·g_strand + (1−w)·g_count`` with ``w=(2κ−1)²`` (the strand discriminability). At
high strand specificity ``w→1`` (the unbiased strand posterior governs); at unstranded ``w→0`` (the
count module governs) — a smooth blend with no gate, the decoupling that fixed the joint model's
bias-mixing (see ``docs/calibration/decoupled_calibration_design.md`` +
``count_channel_capture_design.md``; joint archived in ``archive/joint_deconvolution.md``)::

    substrate
      -> strand balance: rna_sense_frac (κ) -> strand weight w=(2κ−1)²
      -> node_gdna_density (count module; RAW density, no strand clean)
      -> gdna/rna strand Beta-Binomial overdispersions (strand module parameters)
      -> deconv_regions / deconv_sides (per-node blend: w·strand posterior + (1−w)·count fraction)
      -> derive (gdna_density_global, gdna_geom_len)

The global density and geometric gDNA length are **derived** from the aggregate deconvolved mass --
not looped -- which dissolved the old EM loop's sparse-data collapse. The per-region gDNA length
contraction under capture is the IPR of the deconvolved gDNA mass, applied later in
``priors.assemble_priors``.
"""

from __future__ import annotations

import dataclasses
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
from .splice_junction import region_splice_gdna_frac
from .strand_balance import fit_strand_balance
from .strand_deconv import deconv_regions, deconv_sides
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
    rna_fl_pmf: "np.ndarray",
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
    rna_fl_mean = boundary_eff_length(rna_fl_pmf)  # RNA boundary-crossing eff length (splice fraction)

    # RNA strand balance: rna_sense_frac (κ) = posterior-mean spliced sense fraction. The strand
    # channel's discriminability w=(2κ−1)² (set inside the deconv) is the smooth strand→count
    # deference weight — there is no hard identifiability gate (an unstranded library has κ≈½ ⇒ w≈0 ⇒
    # count governs, regardless of depth).
    balance = fit_strand_balance(strand_model)
    if balance.fallback_used:
        # No spliced reads at all — not a usable RNA-seq library. Fail loudly (a real RNA-seq library
        # always carries spliced reads); see CalibrationStrandError.
        raise CalibrationStrandError(
            "the library has zero spliced unique-mapper observations; this does not look like an "
            "RNA-seq library. A real RNA-seq library always carries spliced reads."
        )
    rna_sense_frac = float(balance.rna_sense_frac)

    # Count clue (the count module): per-region gDNA density by LOCAL boundary-anchored imputation on
    # RAW unspliced counts (no strand cleaning — the strand module owns the strand channel). An
    # observable region uses its own contained density; a non-observable (exon/AMBIG) region is
    # anchored from its observable boundary sides; the global fallback comes from intergenic regions.
    # node_gdna_density returns the count module's gDNA fraction (count_gdna_frac) per node, consumed
    # for count-routed nodes and as the gDNA strand-fit seed weight.
    node_density = node_gdna_density(substrate, region_arrays, region_eff_len, fl_mean)

    # Splice-junction gDNA-fraction (Phase 4-mean): for exon regions with an eligible splice-junction
    # boundary, replace the absolute-density count fraction with the boundary gDNA-fraction — the
    # crossing-spliced reads are a clean mature reference, so the fraction partitions the region's own
    # (capture-enriched) total and avoids the absolute density's boundary-depletion under-count.
    # Ineligible regions keep the absolute-density fallback. This upgrades only the REGION count
    # fraction consumed by deconv_regions; the seed fit and side deconvolution keep node_density.
    region_count_frac, n_splice_upgraded = region_splice_gdna_frac(
        substrate, region_arrays, node_density.count_gdna_frac, eff_gdna=fl_mean, eff_rna=rna_fl_mean
    )
    region_node_density = dataclasses.replace(node_density, count_gdna_frac=region_count_frac)

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

    # Per-node deconvolution by precision-weighted strand→count deference: g = w·g_strand +
    # (1−w)·g_count with w=(2κ−1)² (set inside the deconv). w→1 at high strand specificity, →0 at
    # unstranded (count governs) — a smooth blend, no gate. Strand-unobservable nodes are count-only.
    regions = deconv_regions(
        substrate,
        region_arrays,
        region_node_density,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        deconv_quantile=config.gdna_deconv_quantile,
        n_grid=config.n_grid,
    )
    left, right = deconv_sides(
        substrate,
        region_arrays,
        node_density,
        boundary_eff_len,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        deconv_quantile=config.gdna_deconv_quantile,
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
        "splice_gdna_frac_upgraded=%d "
        "gdna_strand_overdispersion=%.4g (%d seed nodes, %d frags%s) "
        "rna_strand_overdispersion=%.4g (%d seed sides, %d frags%s) "
        "[boundary-spliced sense_frac=%.3f vs κ=%.3f]",
        result.n_regions,
        result.gdna_density_global,
        rna_sense_frac,
        n_splice_upgraded,
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
