"""calibrate() — the fractional-accumulator calibrator (per-node simplex solve, iterated all-gDNA bootstrap).

The contained-region split is solved per node on the **2-simplex** ``(sense-RNA / antisense-RNA / gDNA)``
(calibration models RNA-vs-gDNA only; the per-locus EM separates nascent from mature downstream), iterated
over an **all-gDNA bootstrap**. Per the count-zero-information architecture (CALIBRATION_ARCHITECTURE.md §0),
a node's gDNA/RNA composition is determined by exactly two terms here — a **strand likelihood** (the
Beta-Binomial tilt of the per-strand counts, the only INTRINSIC per-node signal; the count enters ONLY as its
overdispersed Fisher information) and a **global gDNA prior** (the foundation, ``ρ_global`` at the MAD-spread
precision ``τ_global``, governing nodes the strand leaves undetermined)::

    substrate
      -> strand balance: rna_sense_frac (κ)
      -> node_gdna_density (RAW) -> fit gdna/rna strand Beta-Binomial overdispersions (seed)
      -> strand_deconvolve -> cleaned_gdna_count: clean the BOUNDARY crossings (nascent removal)
      -> node_gdna_density (RAW contained + CLEANED boundaries): the boundary anchors' count input
      -> deconv_sides ONCE: the fixed boundary gDNA anchors (the boundary-flux transport into priors)
      -> for each pass (signature-binary all-gDNA init):
           rho_global re-fit on the PREVIOUS pass's gDNA estimate (the iterating foundation)
           -> deconv_regions_sweep: the per-node strand + global simplex solve
           -> converge on per-node f_g
      -> gdna_density_global (the library-average density QC scalar)

Pass 0 is the **signature-binary all-gDNA init**: every node starts at ``f_gdna=1`` carrying its full unspliced
MASS (count-free; CALIBRATION_ARCHITECTURE §3) ⇒ ``ρ_global`` = the observable-node total mass density, a
deliberate over-estimate the iteration drives down as the strand removes RNA. Cross-node imputation (the
reliability-weighted neighbour propagation) is deferred to Step 2; until then a strand-weak node falls to the
global foundation (``ρ_global ≈ 0`` in a pure-RNA library ⇒ no phantom gDNA). The boundary SIDES are
deconvolved ONCE (``deconv_sides``) as the fixed gDNA anchors whose mass feeds the per-locus prior via the
boundary-flux transport (``priors``), post-loop. The library-average gDNA density (a QC scalar) is **derived**
from the aggregate deconvolved mass. The per-region gDNA length contraction under capture is the IPR of the
deconvolved gDNA mass over the per-region effective supports (``gdna_region_eff_len`` + ``gdna_boundary_len``),
applied later in ``priors.assemble_priors``.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from .density_model import node_gdna_density
from .derive import gdna_density_global
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
from .imputation import gdna_imputation_prior
from .result import CalibrationResult
from .run_fill import same_ref_left_right
from .simplex_sweep import deconv_regions_sweep
from .strand_balance import fit_strand_balance
from .strand_deconv import cleaned_gdna_count, deconv_sides, strand_deconvolve
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

    Iterative all-gDNA bootstrap (``config.sweep_max_passes`` passes, re-fitting ``ρ_global`` + the
    gDNA var~mean each pass and converging on per-node ``f_g``); see the module docstring for the data
    flow. ``gdna_density_global``
    may be ``0`` (a zero-gDNA library) and a node's deconvolved gDNA mass may be ``0`` (a pure-RNA
    node); both are valid, graceful outputs — not failures.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)

    # gDNA fragment-length effective lengths (PR 4c geometry): the region-contained
    # length, the per-side boundary density length, and the region-free crossing mean.
    region_eff_len = region_eff_length(region_arrays.region_size_bp, gdna_fl_pmf)
    boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, region_arrays.region_size_bp)
    fl_mean = boundary_eff_length(gdna_fl_pmf)
    # ``rna_fl_pmf`` feeds the RNA-side effective lengths for the cross-node RNA imputation, which is
    # deferred to Step 2 (CALIBRATION_ARCHITECTURE §8); it is unused by the Step-1 strand+global solve.

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
    # pass must come first. This is NOT the region answer (the cleaned pass below is).
    node_density_raw = node_gdna_density(substrate, region_arrays, region_eff_len, fl_mean)

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
        gdna_counts=(_raw_count(substrate.contained), cleaned_left, cleaned_right),
    )

    # THE SOLVE — the iterative per-node simplex solve (CALIBRATION_ARCHITECTURE §4). Each pass re-fits the
    # global density ρ_global on the PREVIOUS pass's gDNA estimate, then re-solves. Pass 0 is the
    # signature-binary **all-gDNA init** (every node ``f_gdna=1`` carrying its full unspliced MASS, count-free
    # ⇒ ρ_global = the observable-node mass density, a deliberate over-estimate iteration drives down as the
    # strand removes RNA). The boundary SIDES are deconvolved ONCE (deconv_sides) as the fixed gDNA anchors
    # whose mass feeds the per-locus prior via the boundary-flux transport (priors), post-loop.
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
    eff_len = np.maximum(np.asarray(region_eff_len, dtype=np.float64), 1e-9)
    mass_unspl = np.asarray(c.mass_unspliced, dtype=np.float64)
    mass_u = np.maximum(mass_unspl, 1e-9)
    # geom2 = the density→fraction jacobian (eff_len/mass_u)² — the SAME mass_u that DEFINES the fraction
    # f_g = M_gdna/mass_u (the legitimate normalizer, NOT a count-confidence; CALIBRATION_ARCHITECTURE §2).
    # It converts the global prior's MAD density-spread σ²_global to a fraction-variance.
    geom2 = (eff_len / mass_u) ** 2
    obs = np.asarray(node_density.region_count_observable, dtype=bool)  # signature partition (intergenic/intron)
    has_data = mass_unspl > 0.0  # data-presence mask (which nodes carry unspliced fragments)
    # Signature-binary all-gDNA init (CALIBRATION_ARCHITECTURE §3): every node f_gdna=1 ⇒ its gDNA MASS is the
    # full unspliced MASS (count-FREE; not the raw count u_tot). The locked strand axes are realized by the
    # sweep's allow_pos/allow_neg forbid mask; here the seed is only the gDNA mass the loop re-derives ρ_global
    # from. Iteration drives it down as the strand separates RNA.
    gdna_c = mass_u.copy()
    # gDNA imputation adjacency (CALIBRATION_ARCHITECTURE §1.2): the imputed dests are the non-count-
    # observable (exon / AMBIG) regions; each flanking boundary side that exists (same-ref), is
    # count-observable, and carries deconv'd gDNA mass is a predictor. The side gDNA densities use the FIXED
    # deconv_sides gDNA mass / the per-side density length (boundary_eff_len) — boundaries remain fixed
    # anchors until Step 3 promotes them to solved nodes, so these predictors are stable across passes.
    bco = np.asarray(node_density.boundary_count_observable, dtype=bool)
    ls_same, rs_same = same_ref_left_right(np.asarray(region_arrays.ref_id))
    left_obs = np.zeros(region_arrays.n_regions, dtype=bool)
    right_obs = np.zeros(region_arrays.n_regions, dtype=bool)
    if region_arrays.n_regions > 1:
        left_obs[1:] = bco[:-1] & ls_same[1:]  # region r's LEFT side = boundary (r−1, r)
        right_obs[:-1] = bco[:-1] & rs_same[:-1]  # region r's RIGHT side = boundary (r, r+1)
    region_eligible_g = ~obs
    inv_side_len = 1.0 / np.maximum(boundary_eff_len, 1e-9)
    gdna_left = np.asarray(left.gdna_mass, dtype=np.float64)
    gdna_right = np.asarray(right.gdna_mass, dtype=np.float64)
    d_left = gdna_left * inv_side_len  # fixed boundary-side gDNA crossing densities (the predictors)
    d_right = gdna_right * inv_side_len
    left_ok_g = left_obs & (gdna_left > 0.0)
    right_ok_g = right_obs & (gdna_right > 0.0)
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
        # σ²_global = the robust between-node SPREAD of the per-node densities (MAD) — the population spread an
        # unanchored node faces (NOT a single node's sampling variance). The global gDNA prior is the
        # foundation (CALIBRATION_ARCHITECTURE §1.3): it governs the nodes the strand leaves undetermined
        # (κ=½ / AMBIG / thin). 1.4826 = the normal-consistency MAD constant. The spread → 0 at zero-DNA where
        # all observable nodes agree density ≈ 0 (the zero-DNA pin ⇒ no phantom gDNA). τ_global = 1/σ²_global
        # via the density→fraction jacobian geom2, floored at the 1-pseudo-observation foundation; there is NO
        # count cap (the mass_u ceiling was a §0 violation — CALIBRATION_ARCHITECTURE §6.2).
        mu = gdna_c / eff_len
        active_mu = mu[has_data]
        sigma_d_global = (
            1.4826 * float(np.median(np.abs(active_mu - np.median(active_mu))))
            if active_mu.size
            else rho_global
        )
        sig2_glob = np.maximum(sigma_d_global**2 * geom2, 1e-12)
        tau_global = np.maximum(1.0 / sig2_glob, 1.0)
        # gDNA imputation prior on f_g (CALIBRATION_ARCHITECTURE §1.2): the exon/AMBIG regions get a prior
        # from their observable boundary-side gDNA crossings — mean = the side density (identity), precision
        # = 1/(σ²_bio(ρ_g) + the predictor's Poisson sampling noise). This is the strand-weak rescue (the
        # propagation removed in Step 1) and reinforces f_g=0 at zero-gDNA (the flanks carry ~0 gDNA there).
        gdna_imp = gdna_imputation_prior(
            mu, d_left, d_right,
            region_eff_len=eff_len, side_eff_len=boundary_eff_len, mass_u=mass_u,
            region_eligible=region_eligible_g, left_ok=left_ok_g, right_ok=right_ok_g,
            ref_id=region_arrays.ref_id,
        )
        regions = deconv_regions_sweep(
            substrate,
            region_arrays,
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            n_grid=config.sweep_n_grid,
            rho_global=rho_global,
            region_eff_len=region_eff_len,
            global_tau=tau_global,
            gdna_imp_frac=gdna_imp.frac,
            gdna_imp_precision=gdna_imp.precision,
        )
        if (
            prev_fg is not None
            and float(np.mean(np.abs(regions.gdna_frac - prev_fg)))
            < config.sweep_convergence_delta
        ):
            break
        prev_fg = regions.gdna_frac
        gdna_c = np.asarray(regions.gdna_mass, dtype=np.float64)

    # Derive gdna_density_global (the library-average density QC scalar).
    density_global = gdna_density_global(
        regions, left, right, region_eff_len, boundary_eff_len, region_arrays.ref_id
    )

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
        gdna_boundary_len=boundary_eff_len,
        gdna_region_eff_len=np.asarray(region_eff_len, dtype=np.float64),
        gdna_density_global=density_global,
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
