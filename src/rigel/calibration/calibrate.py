"""calibrate() — the fractional-accumulator calibrator (iterative odds-propagation simplex sweep).

The contained-region split is solved by the **odds-propagation grid sum-product** on the per-node
2-simplex (gDNA / mature-RNA / nascent-RNA), iterated over an **all-gDNA bootstrap**. Each node combines
three Bayesian terms — a **strand likelihood** (the Beta-Binomial tilt of the per-strand counts, the
capture-INVARIANT *direction* signal), a **count prior** (the local gDNA-density imputation weighted by its
``var~mean`` precision ``τ_count``, the *magnitude* signal), and a **global gDNA prior** (the foundation,
``ρ_global`` at precision ``τ_global``) — and propagates the per-strand RNA:gDNA log-odds along same-strand
exon stretches so that AMBIG (overlapping opposite-strand) nodes are resolved from their stranded
neighbours::

    substrate
      -> strand balance: rna_sense_frac (κ)
      -> node_gdna_density (RAW) -> fit gdna/rna strand Beta-Binomial overdispersions (seed)
      -> strand_deconvolve -> cleaned_gdna_count: clean the BOUNDARY crossings (nascent removal)
      -> node_gdna_density (RAW contained + CLEANED boundaries) -> count prior g_count
         (+ splice-junction gDNA-FRACTION upgrade, 3-term, where eligible)
      -> deconv_sides ONCE: the fixed boundary gDNA anchors for the var~mean's boundary->region imputation
      -> for each pass (all-gDNA init):
           rho_global + var~mean (direct / imputation) re-fit on the PREVIOUS pass's gDNA estimate
           -> deconv_regions_sweep: the per-node 3-term simplex solve + odds propagation
           -> converge on per-node f_g
      -> derive (gdna_density_global, gdna_geom_len)

Pass 0 is the all-gDNA init (every unspliced fragment is gDNA ⇒ ``ρ_global`` = the count-observable total
density, a deliberate over-estimate the iteration drives down as RNA is removed). The ``var~mean`` reports
genuinely-high variance where the boundaries cannot predict the region (RNA-rich exons) ⇒ the count yields
there ⇒ the strand + propagation + global solve ⇒ the next pass's estimate is better. The contained count
stays RAW (a region's strand enters once, via the strand likelihood); only the boundary crossings are
cleaned, because they feed the spatial imputation of exon / AMBIG regions. The global density and geometric
gDNA length are **derived** from the aggregate deconvolved mass. The per-region gDNA length contraction
under capture is the IPR of the deconvolved gDNA mass, applied later in ``priors.assemble_priors``.
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
from .simplex_sweep import deconv_regions_sweep
from .splice_junction import region_splice_gdna_frac
from .strand_balance import fit_strand_balance
from .strand_deconv import cleaned_gdna_count, deconv_sides, strand_deconvolve
from .substrate import CalibrationSubstrate
from .variance_model import fit_direct_varmean, fit_imputation_varmean, varmean_points

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
    # RNA region (contained) eff length — pairs with region_eff_len (gDNA) to convert the splice
    # boundary *density* fraction into the region *count* fraction (FL-consistency; see
    # docs/calibration/fl_consistency_diagnostic.md).
    region_eff_len_rna = region_eff_length(region_arrays.region_size_bp, rna_fl_pmf)

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

    # Count clue on RAW counts (the count module, pre-cleaning): per-region gDNA density by LOCAL
    # boundary-anchored imputation. Needed here only to fit the gDNA strand overdispersion (its seed
    # identification is pre-cleaning) — the cleaning below depends on that overdispersion, so the raw
    # pass must come first. This is NOT the region answer (the cleaned pass below is). The count
    # posterior variance feeds only the FP-rate quantile of the boundary-side deconvolution (the sweep's
    # per-node τ_count comes from the var~mean, not from here), so compute it only for a non-½ quantile;
    # else skip its O(R²) LOESS.
    need_count_variance = float(config.gdna_deconv_quantile) != 0.5
    node_density_raw = node_gdna_density(
        substrate, region_arrays, region_eff_len, fl_mean, need_count_variance=need_count_variance
    )

    # Strand-module parameters — the two Beta-Binomial overdispersions. gDNA (mean ½) fitted from the
    # count-observable seed regions/sides using the raw count-clue gDNA weight (breaks the circularity:
    # the seed weight is the strand MEAN ½, not the dispersion). RNA (mean κ) fitted from boundary-side
    # spliced counts. Both shrunk toward the SAME default prior, so under sparse data they collapse to
    # one distribution and an unstranded node (κ=½) is uninformative. See docs/em_strand/03+05.
    gdna_strand = fit_gdna_strand_from_substrate(
        substrate,
        region_arrays,
        node_density_raw,
        boundary_eff_len,
        rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(config.gdna_strand_prior_alpha_beta),
        prior_weight=config.gdna_strand_prior_weight,
    )
    gdna_strand_overdispersion = gdna_strand.gdna_strand_overdispersion
    rna_strand = fit_rna_strand_from_substrate(
        substrate,
        rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(config.rna_strand_prior_alpha_beta),
        prior_weight=config.rna_strand_prior_weight,
    )
    rna_strand_overdispersion = rna_strand.rna_strand_overdispersion

    # Strand deconvolution → CLEAN the BOUNDARY counts for the count module's spatial imputation. The
    # strand emits per node the gDNA fraction g_strand + its information I=N·(2κ−1)²; cleaned_gdna_count
    # removes the strand-identified RNA from a count by w·g_strand+(1−w), w=I/(I+I₀). Cleaning the
    # boundary crossings (left/right) makes the imputed density at exon / AMBIG regions drop the nascent
    # the count clue can't see (the Phase-2 AMBIG fix). The CONTAINED count stays RAW — a region's strand
    # enters the solve once, via the sweep's strand likelihood, so g_count carries count MAGNITUDE only
    # (orthogonal to strand direction ⇒ no double-count). Only the boundary splits are used here; the
    # contained split is discarded (the sweep recomputes the contained strand likelihood internally).
    _, left_split, right_split = strand_deconvolve(
        substrate,
        region_arrays,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        deconv_quantile=config.gdna_deconv_quantile,
        n_grid=config.n_grid,
    )
    i0 = config.gdna_strand_info_scale

    def _raw_count(view):
        return view.n_unspliced_pos.astype(np.float64) + view.n_unspliced_neg.astype(np.float64)

    cleaned_left = cleaned_gdna_count(left_split, _raw_count(substrate.left), i0)
    cleaned_right = cleaned_gdna_count(right_split, _raw_count(substrate.right), i0)

    # Count module g_count: per-region gDNA density by LOCAL imputation — RAW contained own-density for
    # signature count-observable regions (intergenic / intron); the CLEANED boundary crossings impute
    # exon / AMBIG regions (nascent removed → the Phase-2 AMBIG fix). This is the count-only spatial
    # estimate the combine blends with the strand (g_count = magnitude; g_strand = direction).
    node_density = node_gdna_density(
        substrate,
        region_arrays,
        region_eff_len,
        fl_mean,
        need_count_variance=need_count_variance,
        gdna_counts=(_raw_count(substrate.contained), cleaned_left, cleaned_right),
    )

    # gDNA-FRACTION imputation (splice-junction method, 3-term): for exon regions with an eligible
    # splice-junction boundary, replace the absolute gDNA-DENSITY count fraction with the boundary
    # gDNA-FRACTION — the crossing-spliced reads are a clean mature reference, so the fraction partitions
    # the region's own (capture-enriched) total and avoids the gDNA-density method's boundary-depletion
    # under-count. The **3-term** form uses the strand-CLEANED gDNA crossings (cleaned_left/right) so
    # nascent moves to the RNA side. This upgrades the count prior (g_count) for eligible regions; the
    # sweep's count term weighs it by τ_count (it is not an override — where the strand is informative its
    # likelihood governs anyway).
    region_count_frac, n_splice_upgraded = region_splice_gdna_frac(
        substrate,
        region_arrays,
        node_density.count_gdna_frac,
        eff_gdna=fl_mean,
        eff_rna=rna_fl_mean,
        eff_gdna_region=region_eff_len,
        eff_rna_region=region_eff_len_rna,
        left_gdna_unspl=cleaned_left,
        right_gdna_unspl=cleaned_right,
    )

    # THE SOLVE — the iterative odds-propagation simplex sweep (CALIBRATION_PLAN_v2 §2/§5). Each pass
    # re-fits the gDNA var~mean + the global density ρ_global on the PREVIOUS pass's gDNA estimate, then
    # re-solves. Pass 0 is the **all-gDNA init** (every unspliced fragment is gDNA ⇒ ρ_global = the
    # count-observable total density, a deliberate over-estimate that iteration drives down as RNA is
    # removed). The var~mean reports genuinely-high variance where the boundaries cannot predict the region
    # (RNA-rich exons) ⇒ the count yields there ⇒ the strand + propagation + global solve ⇒ the next pass's
    # estimate is better. The boundary SIDES are deconvolved ONCE (the per-node strand/count combine,
    # deconv_sides) as the fixed boundary gDNA anchors for the var~mean's boundary→region imputation; their
    # gDNA mass also feeds the per-locus prior via the boundary-flux transport (priors), post-loop.
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
        info_scale=i0,
    )
    c = substrate.contained
    u_tot = c.n_unspliced_pos.astype(np.float64) + c.n_unspliced_neg.astype(np.float64)
    eff_len = np.maximum(np.asarray(region_eff_len, dtype=np.float64), 1e-9)
    mass_u = np.maximum(np.asarray(c.mass_unspliced, dtype=np.float64), 1e-9)
    geom2 = (eff_len / mass_u) ** 2  # density-var → fraction-var
    obs = np.asarray(node_density.region_count_observable, dtype=bool)
    gdna_left = np.asarray(left.gdna_mass, dtype=np.float64)
    gdna_right = np.asarray(right.gdna_mass, dtype=np.float64)
    gdna_c = u_tot.copy()  # Pass-0 all-gDNA: every unspliced fragment is gDNA
    regions = None
    prev_fg = None
    for _pass in range(int(config.sweep_max_passes)):
        # ρ_global = the gDNA-baseline density, from the COUNT-OBSERVABLE nodes only (intergenic/intron —
        # where gDNA is known), using their current gDNA estimate (iterating). NOT all nodes: an all-node
        # average is irreducibly inflated by AMBIG/unstranded RNA the strand can't remove, so it never
        # drops and the tightening global locks in a phantom (observed). The var~mean DOES use all nodes
        # (it needs the enriched range); only the baseline density restricts to the observable gDNA.
        rho_global = (
            float(gdna_c[obs].sum() / max(float(eff_len[obs].sum()), 1e-9)) if obs.any() else 0.0
        )
        # var~mean on the CURRENT gDNA estimate (region contained + the fixed boundary anchors). DIRECT
        # (count-observable own variance) and IMPUTATION (boundary→region prediction error); read each
        # node in its own class at its current gDNA density μ. Density-var → fraction-var via (eff/mass)²;
        # τ = 1/σ². Count τ ≤ own-count Poisson ceiling (mass); global τ ≥ 1-pseudo-obs foundation,
        # ≤ own data. σ²_global = the between-node SPREAD (robust MAD) — wide/uncertain on Pass 0
        # (all-gDNA spans depleted→enriched), tightening as iteration converges (the zero-DNA pin).
        pts = varmean_points(
            substrate, region_arrays, region_eff_len, fl_mean,
            gdna_views=(gdna_c, gdna_left, gdna_right),
        )
        direct = fit_direct_varmean(pts)
        imputation = fit_imputation_varmean(pts)
        mu = gdna_c / eff_len
        var_d = np.where(obs, direct.predict(mu), imputation.predict(mu))
        sig2_frac = np.maximum(var_d * geom2, 1e-12)
        active_mu = mu[u_tot > 0.0]
        sigma_d_global = (
            1.4826 * float(np.median(np.abs(active_mu - np.median(active_mu))))
            if active_mu.size else rho_global
        )
        sig2_glob = np.maximum(sigma_d_global**2 * geom2, 1e-12)
        tau_count = np.minimum(1.0 / sig2_frac, mass_u)
        tau_global = np.clip(1.0 / sig2_glob, 1.0, mass_u)
        regions = deconv_regions_sweep(
            substrate,
            region_arrays,
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            count_gdna_frac=region_count_frac,
            count_precision=tau_count,
            n_grid=config.sweep_n_grid,
            rho_global=rho_global,
            region_eff_len=region_eff_len,
            info_scale=i0,
            global_tau=tau_global,
        )
        if prev_fg is not None and float(np.mean(np.abs(regions.gdna_frac - prev_fg))) < 1e-3:
            break
        prev_fg = regions.gdna_frac
        gdna_c = np.asarray(regions.gdna_mass, dtype=np.float64)

    # Derive gdna_density_global (QC scalar) and the geometric gDNA length.
    derived = derive(regions, left, right, region_eff_len, boundary_eff_len, region_arrays.ref_id)

    # Spliced RNA mass per region (Σ over the 3 nodes). The deconv adds the full per-node
    # spliced mass to rna_mass (rna = (1−g)·M_unspliced + M_spliced), so this is exactly the
    # spliced component of mass_rna_*. assemble_priors withholds it from rna_prior_count — a spliced
    # fragment is guaranteed-RNA in the EM (no gDNA candidate), so it must not load the RNA side of
    # the gDNA-vs-RNA *unspliced* split (the prior's a_g+a_r should equal the unspliced competing
    # mass). mass_rna_* stays spliced-inclusive so node conservation gdna+rna = total holds.
    mass_rna_spliced = (
        np.asarray(substrate.contained.mass_spliced, dtype=np.float64)
        + np.asarray(substrate.left.mass_spliced, dtype=np.float64)
        + np.asarray(substrate.right.mass_spliced, dtype=np.float64)
    )

    result = CalibrationResult(
        mass_gdna_contained=regions.gdna_mass,
        mass_rna_contained=regions.rna_mass,
        mass_gdna_left=left.gdna_mass,
        mass_rna_left=left.rna_mass,
        mass_gdna_right=right.gdna_mass,
        mass_rna_right=right.rna_mass,
        mass_rna_spliced=mass_rna_spliced,
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
    # deconv mean). A large gap flags a strand-model / accumulator mismatch (we do NOT refit the
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
