"""calibrate() — the calibrator: deconvolve each node's UNSPLICED mass into ``(RNA₊, RNA₋, gDNA)``.

Calibration models **RNA vs gDNA only** (the per-locus EM separates nascent from mature downstream). Per the
count-zero-information architecture (`CALIBRATION_ARCHITECTURE.md` §0) a node's composition is set by three
sources — the **strand likelihood** (the Beta-Binomial tilt of the per-strand counts, the only INTRINSIC
signal; the count enters ONLY as its overdispersed Fisher information), the **cross-node imputation** (the
neighbour density messages, at the modeled var~mean reliability), and the **global gDNA prior** (the
population baseline ``ρ_global`` at the MAD-spread precision). The solver is the bipartite belief-propagation
SWEEP over the unified region↔boundary chain (:mod:`rigel.calibration.bp_solver`)::

    substrate (+ boundary_substrate)
      -> strand balance: rna_sense_frac (κ)
      -> node_gdna_density (RAW) -> fit gDNA / RNA strand Beta-Binomial overdispersions (seed)
      -> build chain + geometry + statics -> signature-binary init (G1/G2/G3)
      -> node_sweep: directional (L→R, R→L) Gauss-Seidel passes — each node integrates its strand
           likelihood + the node-class prior + the gDNA & per-strand RNA identity-density messages from
           its sweep-direction neighbour (precision = the per-pass frozen-snapshot var~mean reliability);
           ρ_global re-fit each pass. Converge on the per-node pie.
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
import os
from dataclasses import replace
from typing import TYPE_CHECKING

import numpy as np

from .bp_solver import (
    build_node_geometry,
    build_node_statics,
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
from .node_chain import build_node_chain
from .result import CalibrationResult
from .strand_balance import fit_strand_balance
from .substrate import BoundarySubstrate, CalibrationSubstrate

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays

logger = logging.getLogger(__name__)

#: Below this many single-strand teacher nodes the KDE is degenerate (a density from a handful of points has
#: no meaningful modes), so calibration falls back to the pass-1 prior-free belief rather than train + apply
#: a garbage mixture. Real genomes have millions of teachers; only tiny single-locus scenarios trip this.
_MIN_KDE_TEACHERS = 10


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

    # THE SOLVE — the bipartite belief-propagation sweep (plan v3 / bp_solver): build the unified
    # region↔boundary chain + its per-node geometry / statics, the signature-binary init, then the directional
    # Gauss-Seidel sweep (gDNA + per-strand RNA identity-density messages, per-pass frozen-snapshot var~mean
    # reliability, the global gDNA prior). The region nodes give the per-region gDNA fraction; the boundary
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
        n_grid=config.sweep_n_grid, logodds_window=config.sweep_logodds_window, statics=statics,
    )
    def _sweep(prior):
        return node_sweep(
            chain, statics, geometry, belief, region_arrays, boundary_substrate,
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            n_grid=config.sweep_n_grid, max_passes=config.sweep_max_passes,
            convergence_delta=config.sweep_convergence_delta,
            logodds_window=config.sweep_logodds_window,
            n_tilt=config.sweep_n_tilt, gdna_prior=prior,
        )

    # PASS 1 — the single-strand solve with the extremely-weak stability floor (Phase 1).
    belief = _sweep(None)
    # PHASE 2 — train the nonparametric gDNA-density mixture prior on the solved single-strand nodes, then
    # PASS 2 — re-solve ALL nodes with that mixture as the per-node prior (self-scaling; fills the tilt's
    # null space on AMBIG). Falls back to the pass-1 belief if the substrate is too small to fit.
    gdna_prior = None
    train_sub = build_training_substrate(
        chain, belief, geometry, statics, region_arrays, boundary_substrate,
        min_eff_length=fl_mean,  # exclude regions too short to contain a gDNA fragment (§8e)
    )
    if train_sub.n >= _MIN_KDE_TEACHERS:
        gdna_prior = GdnaDensityPrior.fit(
            train_sub, bandwidth=config.gdna_prior_bandwidth,
            mixture_bridge=config.gdna_prior_mixture_bridge,
        )
        belief = _sweep(gdna_prior)
    regions = chain_region_deconv(chain, belief, substrate)
    left, right = chain_boundary_side_deconv(chain, belief, substrate)
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
    )
    # EXPERIMENTAL (eff-len study, env-gated): replace the DECONVOLVED per-node masses with the TRUE
    # per-node masses (oracle calibration) so any residual EM leak is purely the eff-len formula, not
    # deconvolution error. Keeps the geometry (eff-lens) + hyperparameters. See scripts/debug/oracle_calibration.py.
    _oracle = os.environ.get("RIGEL_ORACLE_CALIB")
    if _oracle:
        z = np.load(_oracle)
        result = replace(
            result,
            mass_gdna_contained=z["mass_gdna_contained"],
            mass_rna_contained=z["mass_rna_contained"],
            mass_gdna_left=z["mass_gdna_left"],
            mass_rna_left=z["mass_rna_left"],
            mass_gdna_right=z["mass_gdna_right"],
            mass_rna_right=z["mass_rna_right"],
            mass_rna_spliced=z["mass_rna_spliced"],
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
