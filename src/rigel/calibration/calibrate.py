"""calibrate() — the calibrator: deconvolve each node's UNSPLICED mass into ``(RNA₊, RNA₋, gDNA)``.

Calibration models **RNA vs gDNA only** (the per-locus EM separates nascent from mature downstream). Per the
count-zero-information architecture (`CALIBRATION_ARCHITECTURE.md` §0) a node's composition is set by three
sources — the **strand likelihood** (the Beta-Binomial tilt of the per-strand counts, the only INTRINSIC
signal; the count enters ONLY as its overdispersed Fisher information), the **cross-node imputation** (the
neighbour density messages, at the source's own honest belief precision (strand + count)
``pr = n_src/(n_src·vb_src + 1)``), and the **population gDNA prior** (the conservative intergenic+intron floor +
the Phase-2 density KDE). The solver is the belief-propagation SWEEP over the unified region↔boundary chain
(:mod:`rigel.calibration.bp_solver`)::

    substrate (+ boundary_substrate)
      -> strand balance: rna_sense_frac (κ)
      -> node_gdna_density (RAW) -> fit gDNA / RNA strand Beta-Binomial overdispersions (seed)
      -> build chain + geometry + statics -> signature-binary init (G1/G2/G3)
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
    chain_boundary_side_deconv,
    chain_region_deconv,
    init_beliefs,
    node_global_geometry,
    node_sweep,
)
from .density_model import node_gdna_density
from .derive import gdna_density_global
from .errors import CalibrationStrandError
from .gdna_rate_prior import GdnaRatePrior
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


def calibrate(
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    gdna_fl_pmf: "np.ndarray",
    rna_fl_pmf: "np.ndarray",
    config: "CalibrationConfig",
    _debug: dict | None = None,
    diagnostics_out: dict | None = None,
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per node, then derive gdna_density_global.

    Runs the belief-propagation sweep (a single forward-backward pass per phase, resolving the per-node
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
    # (gDNA + per-strand RNA identity-density messages at each source's own honest belief precision, plus the
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
        chain,
        substrate,
        boundary_substrate,
        region_arrays,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_grid=config.sweep_n_grid,
        n_grid_ss=config.sweep_n_grid_single_strand,
        logodds_window=config.sweep_logodds_window,
        statics=statics,
    )
    # Message precision is now the source's OWN honest belief precision (strand + count), computed inside
    # the sweep — there is nothing to fit here. The adjacent-pair overdispersion σ²_transfer that used to be
    # estimated at this point is a PRIOR (and, fit on total density, one that weakened every message
    # dramatically); it is ZERO while the prior-free solve is developed. See `bp_solver._scan`.

    # When ``_debug`` is on, the LAST sweep also fills ``_debug["capture"]`` with the per-node message
    # internals (local vs final belief, each channel's message mode/precision) — the substrate for the
    # message-corruption trace (`scripts/debug/msg_trace.py`). Inert in production.
    def _sweep(prior):
        capture = {} if _debug is not None else None
        out = node_sweep(
            chain,
            statics,
            geometry,
            belief,
            region_arrays,
            boundary_substrate,
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            n_grid=config.sweep_n_grid,
            logodds_window=config.sweep_logodds_window,
            n_tilt=config.sweep_n_tilt,
            n_grid_ss=config.sweep_n_grid_single_strand,
            gdna_prior=prior,
            fold_coarse_k=config.fold_coarse_k,
            fold_fine_k=config.fold_fine_k,
            fold_sigma_coverage=config.fold_sigma_coverage,
            fold_refine_iters=config.fold_refine_iters,
            _capture=capture,
        )
        if capture is not None:
            _debug["capture"] = capture
        return out

    # THE PASS-0 gDNA-RATE PRIOR — fit ONCE on ALL nodes' total unspliced density (f_g=1, belief-free) via the
    # Fixed-Kernel Poisson-lognormal Mixture NPMLE (`GdnaRatePrior`). It is EXTREMELY WEAK (n_eff≈0.15
    # pseudo-obs, never >1), so the strand likelihood and the boundary messages dominate the sweep and peel
    # RNA out of the f_g=1 start. Replaces the retired seed/floor/global + the Phase-2 density KDE with one
    # count-space model (docs/calibration/npmle_struggles.md). Belief-free ⇒ ONE sweep suffices (no train-then-
    # re-solve); a refit on the solved belief is a future enhancement.
    mass_global, eff_global = node_global_geometry(chain, geometry)
    gdna_prior = GdnaRatePrior.fit(
        mass_global, eff_global, bandwidth=config.gdna_rate_prior_bandwidth
    )
    # PASS 0 — the gentle bootstrap. The prior comes from TOTAL density, so it is deliberately conservative
    # (extremely weak, n_eff≈0.15). Messages are gentle for a DERIVED reason rather than a fitted one: a
    # source can send no more precision than its own belief supports (`pr` saturates at `1/vb_src`), so a
    # node the strand cannot peel stays quiet on its own. The strand peels what it can here.
    belief = _sweep(gdna_prior)
    # REFIT ITERATIONS — now that the belief is (partly) peeled, re-estimate the population model on the
    # SOLVED gDNA rates: P(ρ) on the believed gDNA counts, with the belief width τ² as the observation WIDTH.
    # Message precision needs no refit: it tracks the belief automatically (`vb_src` IS the belief variance),
    # so messages strengthen exactly where the data has earned it. Deterministic (fixed kernels; no spline).
    for _ in range(int(config.calib_refit_iters)):
        g_hat = np.asarray(belief.f_g, dtype=np.float64) * mass_global
        gdna_prior = GdnaRatePrior.fit(
            g_hat,
            eff_global,
            var_g=np.asarray(belief.var_gdna, dtype=np.float64),
            bandwidth=config.gdna_rate_prior_bandwidth,
        )
        # ⚠ HISTORICAL, and it still bites: an adjacent-pair σ²_transfer(ρ) refit on the SOLVED rates was
        # built and A/B'd (`message_precision.adjacent_imputation_variance`) and it BACKFIRES — measured on a
        # not-yet-peeled belief, adjacent nodes agree on their *wrong* gDNA rates ⇒ σ² small ⇒ messages become
        # confident and propagate the error (mean mwae 0.151→0.189; the zero-gDNA/unstranded/nascent/capture
        # corner 0.52→0.95). Honesty measured against a wrong belief is not honesty. σ²_transfer is now ZERO
        # (it is a prior); when it returns it must be an NPMLE projection gated on belief QUALITY. Keeping
        # messages weak lets the PRIOR refit alone do the work (mean mwae 0.151→0.118). See npmle_roadmap.md.
        belief = _sweep(gdna_prior)
    regions = chain_region_deconv(chain, belief, substrate)
    left, right = chain_boundary_side_deconv(chain, belief, substrate)
    logger.debug(
        "calibration: pass-0 gDNA-rate NPMLE prior (%d cells) + %d refit iteration(s)",
        gdna_prior.n_cells,
        int(config.calib_refit_iters),
    )

    if (
        _debug is not None
    ):  # inert diagnostic hook — the solved chain internals (Phase-2 substrate + plots)
        _debug.update(
            chain=chain,
            belief=belief,
            geometry=geometry,
            statics=statics,
            substrate=substrate,
            boundary_substrate=boundary_substrate,
            region_arrays=region_arrays,
            gdna_prior=gdna_prior,
            rna_sense_frac=rna_sense_frac,
            region_eff_len=region_eff_len,
            boundary_eff_len=boundary_eff_len,
        )

    # Report-facing diagnostics: the fitted gDNA-rate prior P(ρ) (bimodal ⇒ capture enrichment). Consumed by
    # the QC report, never by the EM.
    if diagnostics_out is not None:
        from .diagnostics import CalibrationDiagnostics

        diagnostics_out["calibration"] = CalibrationDiagnostics.from_prior(gdna_prior)

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
