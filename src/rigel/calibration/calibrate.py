"""calibrate() — the calibrator: deconvolve each node's UNSPLICED mass into ``(RNA₊, RNA₋, gDNA)``.

Calibration models **RNA vs gDNA only** (the per-locus EM separates nascent from mature downstream). Per the
count-zero-information architecture (`CALIBRATION_ARCHITECTURE.md` §0) a node's composition is set by three
sources — the **strand likelihood** (the Beta-Binomial tilt of the per-strand counts, the only INTRINSIC
signal; the count enters ONLY as its overdispersed Fisher information), the **cross-node imputation** (the
neighbour density messages, at the belief-free Poisson disagreement-variance precision
``pr = n_src/(n_src·σ²_imp + 1)``), and the **population gDNA prior** (the conservative intergenic+intron floor +
the Phase-2 density KDE). The solver is the belief-propagation SWEEP over the unified region↔boundary chain
(:mod:`rigel.calibration.bp_solver`)::

    substrate (+ boundary_substrate)
      -> strand balance: rna_sense_frac (κ)
      -> node_gdna_density (RAW) -> fit gDNA / RNA strand Beta-Binomial overdispersions (seed)
      -> build chain + geometry + statics -> signature-binary init (G1/G2/G3)
      -> compute the scalar total-density σ²_imp once (adjacent_disagreement_variance)
      -> PASS 1 node_sweep (no KDE prior): ONE forward + ONE backward pass — each node integrates its strand
           likelihood + the conservative gDNA floor + the belief-free gDNA & per-strand RNA density messages
      -> train the Phase-2 gDNA-density KDE on the pass-1 belief
      -> PASS 2 node_sweep (the KDE prior added per node) -> the converged per-node pie
      -> chain_region_deconv  -> per-region gDNA / RNA contained mass
      -> chain_boundary_side_deconv -> per-region boundary-side flux (gDNA / RNA), for the per-locus prior
      -> gdna_density_global (the library-average density QC scalar)

The boundary nodes are first-class (they carry the one-sided, motif-stranded spliced RNA as a fixed floor +
density term); their per-side gDNA/RNA mass feeds the per-locus prior via :mod:`rigel.calibration.priors`. The
per-region gDNA length contraction under capture is the IPR of the deconvolved gDNA mass over the per-region
effective supports (``gdna_region_eff_len`` + ``gdna_boundary_len``), applied later in ``priors``. A
zero-gDNA library (``gdna_density_global == 0``, per-node gDNA mass ``0``) is a valid, graceful output.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from .bp_solver import (
    build_node_geometry,
    build_node_statics,
    adjacent_disagreement_variance,
    adjacent_disagreement_shrinkage_weight,
    chain_boundary_side_deconv,
    chain_region_deconv,
    init_beliefs,
    node_sweep,
)
from .density_model import node_gdna_density
from .derive import gdna_density_global
from .errors import CalibrationStrandError
from .gdna_density_prior import GdnaDensityPrior, build_training_substrate
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
from .node_chain import REGION, build_node_chain
from .result import CalibrationResult, RnaWarmStart
from .strand_balance import fit_strand_balance
from .substrate import BoundarySubstrate, CalibrationSubstrate

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
    _debug: dict | None = None,
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per node, then derive gdna_density_global.

    Runs the belief-propagation sweep (``config.sweep_max_passes`` passes, converging on the per-node
    pie); see the module docstring for the data flow. ``gdna_density_global`` may be ``0`` (a zero-gDNA
    library) and a node's deconvolved gDNA mass may be ``0`` (a pure-RNA node); both are valid, graceful
    outputs — not failures.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)

    # gDNA fragment-length effective lengths: the region-contained length, the per-side boundary
    # density length, and the region-free crossing mean. ``rna_fl_pmf`` feeds the RNA-side effective
    # lengths the sweep's per-strand RNA messages use (geometry built in bp_solver).
    region_eff_len = region_eff_length(region_arrays.region_size_bp, gdna_fl_pmf)
    boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, region_arrays.region_size_bp)
    fl_mean = boundary_eff_length(gdna_fl_pmf)

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

    # THE SOLVE — the bipartite belief-propagation sweep (bp_solver): build the unified region↔boundary chain
    # + its per-node geometry / statics, the signature-binary init, then a single forward-backward pass
    # (gDNA + per-strand RNA identity-density messages at the belief-free σ²_imp precision fit once, plus the
    # anchored global gDNA prior). The region nodes give the per-region gDNA fraction; the boundary
    # nodes give the per-side boundary flux feeding the per-locus prior (chain_boundary_side_deconv) — the
    # first-class boundary nodes the sweep solves, projected to per-region sides for the prior.
    boundary_substrate = BoundarySubstrate.from_payload(payload)
    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geometry = build_node_geometry(
        chain, substrate, boundary_substrate, region_arrays, gdna_fl_pmf, rna_fl_pmf
    )
    statics = build_node_statics(chain, substrate, boundary_substrate, region_arrays)
    belief = init_beliefs(
        chain, substrate, boundary_substrate, region_arrays,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_grid=config.sweep_n_grid, n_grid_ss=config.sweep_n_grid_single_strand,
        logodds_window=config.sweep_logodds_window, statics=statics,
    )
    # Belief-free Poisson message precision (`disagreement_shrinkage_prior_design_v2.md`): the scalar
    # total-density σ²_imp — the empirical adjacent-node imputation floor. σ²_msg = σ²_imp + 1/n_src; the
    # single production basis for every message channel, both passes.
    sig_total = adjacent_disagreement_variance(chain, geometry)
    # EB shrinkage weight w (data-fit signal fraction; NO tunable — see adjacent_disagreement_shrinkage_weight).
    # Used only by RIGEL_MSG_MODE=eb; the prod/off/max paths ignore it. Cheap (one pass over adjacent edges).
    w_eb = adjacent_disagreement_shrinkage_weight(chain, geometry)
    logger.debug("calibration: total-density σ²_imp=%.4f  EB shrinkage w=%.4f", sig_total, w_eb)

    def _sweep(prior):
        return node_sweep(
            chain, statics, geometry, belief, region_arrays, boundary_substrate,
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            n_grid=config.sweep_n_grid, max_passes=config.sweep_max_passes,
            convergence_delta=config.sweep_convergence_delta,
            logodds_window=config.sweep_logodds_window,
            n_tilt=config.sweep_n_tilt, n_grid_ss=config.sweep_n_grid_single_strand, gdna_prior=prior,
            disagreement_sigma2=sig_total, disagreement_weight=w_eb,
        )

    # PASS 1 — single-strand solve with the total-density floor.
    belief = _sweep(None)
    # PHASE 2 — train the nonparametric gDNA-density mixture prior on the solved single-strand nodes, then
    # PASS 2 — re-solve ALL nodes with that mixture as the per-node prior (self-scaling; fills the tilt's
    # null space on AMBIG). Falls back to the pass-1 belief if the substrate is too small to fit.
    gdna_prior = None
    train_sub = build_training_substrate(
        chain, belief, geometry, statics, region_arrays, boundary_substrate,
        min_eff_length=fl_mean,  # exclude regions too short to contain a gDNA fragment (§8e)
    )
    if train_sub.n >= config.calib_kde_min_training_nodes:
        gdna_prior = GdnaDensityPrior.fit(
            train_sub, bandwidth=config.gdna_prior_bandwidth,
            mixture_bridge=config.gdna_prior_mixture_bridge,
            bridge_trim_pct=config.calib_kde_bridge_trim_pct,
        )
        belief = _sweep(gdna_prior)
    regions = chain_region_deconv(chain, belief, substrate)
    left, right = chain_boundary_side_deconv(chain, belief, substrate)
    rna_warm_start = _build_rna_warm_start(chain, belief, geometry, substrate)
    logger.debug("calibration: %s", "two-pass (Phase-2 mixture prior)" if gdna_prior else "single pass")

    if _debug is not None:  # inert diagnostic hook — the solved chain internals (Phase-2 substrate + plots)
        _debug.update(
            chain=chain, belief=belief, geometry=geometry, statics=statics, substrate=substrate,
            boundary_substrate=boundary_substrate, region_arrays=region_arrays,
            rna_sense_frac=rna_sense_frac, region_eff_len=region_eff_len, boundary_eff_len=boundary_eff_len,
        )

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
        rna_warm_start=rna_warm_start,
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


def _build_rna_warm_start(chain, belief, geometry, substrate) -> RnaWarmStart:
    """Project the solved chain into the per-region per-strand RNA densities the EM warm-start builder needs.

    A pure additive read of the final ``belief`` (``f_pos`` / ``f_neg`` read directly — raw, non-zeroed for
    ALL nodes) over the static per-face ``geometry`` (masses + RNA / spliced effective lengths, already
    floored). The three ROLES map onto the chain by adjacency (verified: a region node's ``chain.right`` is
    the seam boundary to region ``r+1`` when that boundary itself has a right neighbour; its ``chain.left`` /
    ``chain.right`` are always valid boundary nodes):

    * CONTAINED: the region node's own face — ``f_s · mass / eff_rna``.
    * CROSSING (seam ``r``↔``r+1``, left-keyed): the seam boundary node's belief over its POOLED two-side
      crossing mass and averaged RNA support — mirroring the gDNA pooled seam so the builder's RNA and
      reconstructed-gDNA densities share one node basis.
    * SPLICED (per region SIDE): the flanking boundary's one-sided motif-stranded spliced mass over the
      half-triangle spliced support — kept SEPARATE from crossing, single value + its fixed strand.
    """
    kind = np.asarray(chain.kind)
    order = np.asarray(chain.order)
    reg_nodes = order[kind == REGION]
    r = np.asarray(chain.ref_idx, dtype=np.int64)[reg_nodes]  # region index per region node
    lb = np.asarray(chain.left, dtype=np.int64)[reg_nodes]  # region r's left-flank boundary node
    rb = np.asarray(chain.right, dtype=np.int64)[reg_nodes]  # region r's right-flank (seam) boundary node
    seam = np.asarray(chain.right, dtype=np.int64)[rb] != -1  # rb is an internal seam (has region r+1)
    R = int(substrate.n_regions)
    fp, fn, g = belief.f_pos, belief.f_neg, geometry

    rho_c_pos, rho_c_neg = np.zeros(R), np.zeros(R)  # contained, per strand
    rho_x_pos, rho_x_neg = np.zeros(R), np.zeros(R)  # crossing (seam), per strand
    rho_sp_pos_l, rho_sp_neg_l = np.zeros(R), np.zeros(R)  # spliced, region r's LEFT side, per strand
    rho_sp_pos_r, rho_sp_neg_r = np.zeros(R), np.zeros(R)  # spliced, region r's RIGHT side, per strand

    # CONTAINED: region node's own face (mass_left == mass_right, eff_rna_left == eff_rna_right for regions).
    rho_c_pos[r] = fp[reg_nodes] * g.mass_left[reg_nodes] / g.eff_rna_left[reg_nodes]
    rho_c_neg[r] = fn[reg_nodes] * g.mass_left[reg_nodes] / g.eff_rna_left[reg_nodes]

    # CROSSING: pooled seam mass (both boundary faces) / averaged RNA support, tilted by the seam belief.
    m_seam = g.mass_left[rb] + g.mass_right[rb]
    s_seam = 0.5 * (g.eff_rna_left[rb] + g.eff_rna_right[rb])
    rho_x_pos[r[seam]] = fp[rb[seam]] * m_seam[seam] / s_seam[seam]
    rho_x_neg[r[seam]] = fn[rb[seam]] * m_seam[seam] / s_seam[seam]

    # SPLICED (mature, per strand): region r's RIGHT side (donor) rides the seam boundary rb's LEFT (exon)
    # face; its LEFT side (acceptor) rides the left boundary lb's RIGHT (exon) face. Each junction is
    # single-stranded (one of pos/neg nonzero per face); both are nonzero only where a +/− pair share the
    # exact splice coordinate, kept SEPARATE so the antisense isoform reads on its own strand.
    rho_sp_pos_r[r] = g.spliced_pos_left[rb] / g.eff_spl_left[rb]
    rho_sp_neg_r[r] = g.spliced_neg_left[rb] / g.eff_spl_left[rb]
    rho_sp_pos_l[r] = g.spliced_pos_right[lb] / g.eff_spl_right[lb]
    rho_sp_neg_l[r] = g.spliced_neg_right[lb] / g.eff_spl_right[lb]

    return RnaWarmStart(
        rho_contained_pos=rho_c_pos,
        rho_contained_neg=rho_c_neg,
        rho_crossing_pos=rho_x_pos,
        rho_crossing_neg=rho_x_neg,
        rho_spliced_pos_left=rho_sp_pos_l,
        rho_spliced_neg_left=rho_sp_neg_l,
        rho_spliced_pos_right=rho_sp_pos_r,
        rho_spliced_neg_right=rho_sp_neg_r,
    )


__all__ = ["calibrate"]
