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
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .background_reference import BackgroundReference, measure_background
from .bp_solver import (
    NODE,
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
from .npmle import DensityNPMLE
from .density_deconv import GdnaBackground, density_lambda_factor, fit_intron_background
from .gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_sj_table,
)
from .node_chain import build_node_chain
from .result import CalibrationResult
from .gdna_landscape import GdnaLandscape, fit_gdna_landscape
from .signature import RegionType, coarse_type_array
from .simplex_logodds import _logodds_grid
from .strand_balance import fit_strand_balance
from .substrate import CalibrationSubstrate

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class InjectedCalibrationPriors:
    """Population-scale calibration priors — the objects that require genome-scale (or many-gene) data to fit and
    are physically **directly observable** (no deconvolution / no solving): the RNA strand balance, the strand
    Beta-Binomial overdispersions, the strand-Fisher noise-floor sample sizes, the enrichment-density NPMLE (the
    σ²_transfer landscape), the intergenic intron-factory background, and the aggregate ρ_bg background.

    A tiny (single-transcript) toy CANNOT fit these — so :func:`calibrate` accepts them pre-fit from a
    population scenario and injects them, letting the toy provide only the controlled per-node GEOMETRY. Every
    field is optional; ``None`` ⇒ fit that prior internally (the default, byte-identical). ``calibrate`` also
    stashes the fitted-or-injected bundle in ``_debug["calibration_priors"]`` so a population scenario's fitted
    priors can be extracted and re-injected into a toy (`scripts/debug/toy_inject.py`)."""

    rna_sense_frac: float | None = None
    n_rna_obs: float | None = None
    n_gdna_obs: float | None = None
    gdna_strand_overdispersion: float | None = None
    rna_strand_overdispersion: float | None = None
    enrichment_prior: DensityNPMLE | None = None
    intron_background: GdnaBackground | None = None
    background: BackgroundReference | None = None


def _build_intron_prior(chain, substrate, region_arrays, region_eff_len, config, bg=None):
    """The gDNA intron factory λ-factor per chain node (`docs/CARRY_FORWARD.md`).

    Fits the intergenic-background NegBinom (`fit_intron_background`) and tabulates, for each INTRON NODE node,
    ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` over the σ(λ) solve grid → a ``(n_nodes, K)`` array, ZERO on every
    non-intron node (a no-op there). Returns ``None`` when the factory is disabled, the background pool is
    uninformative, or there are no intron nodes — in which case the sweep is byte-identical to the pre-factory
    path. gDNA is strand-symmetric, so this factor lives purely on ``λ`` (peels gDNA; the residual RNA's tilt is
    left to the solver), and is consumed identically by the single-strand and AMBIG per-node solves.

    ``bg`` (an injected population :class:`GdnaBackground`) overrides the internal fit — a tiny toy's own
    intergenic pool is too sparse to fit the background the introns are peeled against."""
    if bg is None:
        bg = fit_intron_background(substrate, region_arrays, region_eff_len, include_introns=False)
    if not bg.informative:
        return None
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(
        np.int64
    )  # 0/1/2 per NODE
    R = rtype.shape[0]
    is_intron = (kind == NODE) & (rtype[np.clip(idx, 0, R - 1)] == 1)  # INTRON == 1
    if not bool(is_intron.any()):
        return None
    _, fg = _logodds_grid(int(config.sweep_n_grid), float(config.sweep_logodds_window))
    prior = np.zeros((kind.shape[0], fg.shape[0]), dtype=np.float64)
    ridx = idx[is_intron]
    count = np.asarray(substrate.contained.n_unspliced, dtype=np.float64)[ridx]
    eff_g = np.asarray(region_eff_len, dtype=np.float64)[ridx]
    prior[is_intron] = density_lambda_factor(bg, count, eff_g, fg)
    return prior


#: Minimum training nodes for a hyperprior fit — below this the population is not a population.
_MIN_TRAIN = 5


def _fit_gdna_hyperprior(
    chain, belief, statics, region_arrays, mass_global, eff_global, *, strength
):
    """Select the training substrate from the chain and fit the :class:`GdnaLandscape` on the initial solve's
    peeled gDNA — the composition (gDNA) arm of ψ for the Phase-2 refit. ``None`` if it cannot be fit.

    **This affects only the PRIOR fit — never the solve's gDNA messages** (the G1/TSS/TES boundary emissions
    in ``node_sweep`` are a separate mechanism and are untouched).

    The substrate is the whole of this function's job; the estimator itself is in
    :mod:`.gdna_landscape`. Four axes decide membership, and conflating them is the mistake this replaces
    (production plan §2.2):

    * **circularity → structural exclusion.** AMBIG nodes are out. They are the two-root ambiguity the prior
      exists to resolve, so training on them would have it predict itself. ⚠ Note the *empirical* case
      alongside this is stratum-dependent — it holds over all conditions but reverses on the ones where an
      enriched mode exists — so the argument carrying the decision is the circularity, not the EMD.
    * **identifiability → structural INCLUSION.** A live region with no unspliced mass has density ``0`` for
      every ``f_g``: "gDNA is absent here" is the strongest depletion evidence there is, and it anchors the
      depleted mode. Dropping it costs +0.26 / +0.61 EMD (+1.04 on zero-gDNA libraries).
    * **precision → a continuous weight**, never admission (`gdna_landscape._reliability`).
    * **geometry → boundaries are EXCLUDED** (owner, 2026-07-27). They cross rather than contain, are ~as
      numerous as regions, and only 5.1 % of them are truly enriched against the regions' 12.1 %; their
      two-flank mixture lands between the two true modes and supplies 74 % of all the mass in the valley.

    ⛔ **Admitting AMBIG into the FINAL fit was implemented, measured and REFUTED (2026-07-28) — do not
    re-propose it without new evidence.** The non-circular form (``pass-0 → fit #1 (out) → re-solve → … →
    final fit (IN) → FINAL solve``, so the admitting fit trains on AMBIG estimates produced by a prior that
    never saw AMBIG) was A/B'd over four paired arms on all 32 conditions. It is **worse on 25/32** at the
    shipped depth: ALL-32 mwae ``0.046675 → 0.048393``, VSTRONG ``+0.005568``, and the zero-gDNA
    false-positive guard regresses ``0.010766 → 0.011065``. ``n_train`` grew 1179 → 1395 (+18.3 %), so it was
    active, not inert. The production plan's "wins on every axis" was scored against the **retired refit=1
    baseline**; at equal depth its own table already showed fabrication regressing 4.7× → 6.7×.
    `SESSION_2026_07_28_HANDOFF_18.md` §2 R3.
    """
    isr = np.asarray(chain.kind) == NODE
    fp = np.asarray(statics.free_pos, dtype=bool)
    fn = np.asarray(statics.free_neg, dtype=bool)
    rtype = coarse_type_array(np.asarray(region_arrays.signature))
    ridx = np.clip(np.asarray(chain.obj_idx, dtype=np.int64), 0, rtype.shape[0] - 1)
    expressed = isr & (eff_global > 1.0e-9) & (mass_global > 1.0e-12)
    # the zero-count structural anchor: an intergenic or intronic region that sequenced no unspliced mass
    anchor = (
        isr & (eff_global > 1.0e-9) & (mass_global <= 1.0e-12) & (rtype[ridx] != RegionType.EXON)
    )
    sel = (expressed & ((fp ^ fn) | (~fp & ~fn))) | anchor
    if int(sel.sum()) < _MIN_TRAIN:
        return None
    mass = np.asarray(mass_global, dtype=np.float64)[sel]
    return fit_gdna_landscape(
        np.asarray(belief.f_g, dtype=np.float64)[sel] * mass,
        mass,
        np.asarray(eff_global, dtype=np.float64)[sel],
        np.asarray(belief.var_gdna, dtype=np.float64)[sel],
        anchor=anchor[sel],
        strength=strength,
    )


def calibrate(
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    gdna_fl_pmf: "np.ndarray",
    rna_fl_pmf: "np.ndarray",
    config: "CalibrationConfig",
    _debug: dict | None = None,
    diagnostics_out: dict | None = None,
    injected_priors: "InjectedCalibrationPriors | None" = None,
    boundary_flags: "np.ndarray | None" = None,
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per node, then derive gdna_density_global.

    Runs the belief-propagation sweep (a single forward-backward pass per phase, resolving the per-node
    pie); see the module docstring for the data flow. ``gdna_density_global`` may be ``0`` (a zero-gDNA
    library) and a node's deconvolved gDNA mass may be ``0`` (a pure-RNA node); both are valid, graceful
    outputs — not failures.

    ``boundary_flags`` is the splice graph's per-boundary structural bits
    (:func:`~rigel.calibration.splice_graph.build_boundary_flags_array`), carried onto the chain as
    ``NodeStatics.boundary_flags``. ``None`` (the default) leaves them zero.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    inj = injected_priors  # population-scale priors to inject in place of the internal (toy-untrustworthy) fits

    # gDNA fragment-length effective lengths: the region-contained length, the per-side boundary
    # density length, and the region-free crossing mean. ``rna_fl_pmf`` feeds the RNA-side effective
    # lengths the sweep's per-strand RNA messages use (geometry built in bp_solver).
    # ⛔ S5.f. Three things below this line no longer exist, and they are named here rather than in a
    # commit message because the body underneath is the specification of what S5.f has to rewire:
    #   * `CalibrationSubstrate` reads payload fields S4 deleted;
    #   * `region_eff_length` became `effective_length.contained_eff_length` (same quantity);
    #   * `boundary_side_eff_length` and `boundary_eff_length` have NO successor — they divided a
    #     per-FACE mass, and a contiguous edge has no faces. Its divisor is `crossing_eff_length`,
    #     one number per edge, and whether that carries the RNA reach taper is the open A7 decision.
    raise NotImplementedError(
        "calibrate() is S5.f. Its substrate still reads payload fields S4 deleted, and its effective "
        "lengths still assume per-face boundary geometry that S5.c removed. "
        "See docs/S5_DESIGN_LOG.md §2."
    )
    # ⚠ Unreachable, and bound only so the pre-S5 body below stays readable AND lintable — `ruff check`'s
    # undefined-name list is S6's authoritative delete list, so it must not fill with noise now.
    region_eff_len = boundary_eff_len = fl_mean = boundary_substrate = None

    # RNA strand balance: rna_sense_frac (κ) = posterior-mean spliced sense fraction. The strand
    # channel's discriminability w=(2κ−1)² (set inside the deconv) is the smooth strand→count
    # deference weight — there is no hard identifiability gate (an unstranded library has κ≈½ ⇒ w≈0 ⇒
    # count governs, regardless of depth).
    if inj is not None and inj.rna_sense_frac is not None:
        # INJECTED population κ + spliced sample size — a tiny toy cannot fit either (scarce spliced).
        rna_sense_frac = float(inj.rna_sense_frac)
        n_rna_obs = float(inj.n_rna_obs) if inj.n_rna_obs is not None else 0.0
    else:
        balance = fit_strand_balance(strand_model)
        if balance.fallback_used:
            # No spliced reads at all — not a usable RNA-seq library. Fail loudly (a real RNA-seq library
            # always carries spliced reads); see CalibrationStrandError.
            raise CalibrationStrandError(
                "the library has zero spliced unique-mapper observations; this does not look like an "
                "RNA-seq library. A real RNA-seq library always carries spliced reads."
            )
        rna_sense_frac = float(balance.rna_sense_frac)
        n_rna_obs = float(balance.n_observations)

    # Count clue on RAW counts (the count module, pre-cleaning): per-region gDNA density by LOCAL
    # boundary-anchored imputation. Needed here only to fit the gDNA strand overdispersion (its seed
    # identification is pre-cleaning) — the cleaning below depends on that overdispersion, so the raw
    # pass must come first. This is NOT the region answer (the cleaned pass below is).
    node_density_raw = node_gdna_density(substrate, region_arrays, region_eff_len, fl_mean)

    # Strand-module parameters — the two Beta-Binomial overdispersions. gDNA (mean ½) fitted from the
    # count-observable seed regions/sides using the raw count-clue gDNA weight (breaks the circularity:
    # the seed weight is the strand MEAN ½, not the dispersion). RNA (mean κ) fitted from the PER-JUNCTION
    # SJ strand table — the same strand-qualified population κ itself is the marginal of, so both halves of
    # the RNA Beta-Binomial come from one source (docs/CARRY_FORWARD.md). Both shrunk
    # toward the SAME default prior, so under sparse data they collapse to one distribution and an
    # unstranded node (κ=½) is uninformative. See docs/em_strand/03+05.
    _gd_seed = _rna_seed = (
        -1,
        -1,
        False,
    )  # (n_seed_nodes, n_seed_frags, fallback) — QC log only; -1 = injected
    if inj is not None and inj.gdna_strand_overdispersion is not None:
        gdna_strand_overdispersion = float(inj.gdna_strand_overdispersion)
    else:
        gdna_strand = fit_gdna_strand_from_substrate(
            substrate,
            region_arrays,
            node_density_raw,
            boundary_eff_len,
            rna_sense_frac=rna_sense_frac,
        )
        gdna_strand_overdispersion = gdna_strand.gdna_strand_overdispersion
        _gd_seed = (
            gdna_strand.n_seed_nodes,
            gdna_strand.n_seed_fragments,
            gdna_strand.fallback_used,
        )
    if inj is not None and inj.rna_strand_overdispersion is not None:
        rna_strand_overdispersion = float(inj.rna_strand_overdispersion)
    else:
        rna_strand = fit_rna_strand_from_sj_table(
            strand_model.sj_table,
            rna_sense_frac=rna_sense_frac,
        )
        rna_strand_overdispersion = rna_strand.rna_strand_overdispersion
        _rna_seed = (rna_strand.n_seed_nodes, rna_strand.n_seed_fragments, rna_strand.fallback_used)

    # THE SOLVE — the bipartite belief-propagation sweep (bp_solver): build the unified region↔boundary chain
    # + its per-node geometry / statics, the signature-binary init, then a single forward-backward pass
    # (gDNA + per-strand RNA identity-density messages at each source's own honest belief precision, plus the
    # anchored global gDNA prior). The region nodes give the per-region gDNA fraction; the boundary
    # nodes give the per-side boundary flux feeding the per-locus prior (chain_boundary_side_deconv) — the
    # first-class boundary nodes the sweep solves, projected to per-region sides for the prior.
    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geometry = build_node_geometry(
        chain, substrate, boundary_substrate, region_arrays, gdna_fl_pmf, rna_fl_pmf
    )
    statics = build_node_statics(
        chain, substrate, boundary_substrate, region_arrays, boundary_flags
    )

    # Strand-Fisher noise-floor SAMPLE SIZES (bp_solver τ seed): N_gdna (gDNA-eligible unspliced fragments in
    # the structurally pure-gDNA intergenic regions, coarse type 0) and N_spliced (the pure-RNA count κ_RNA was
    # fit from). gDNA's sense mean is ½ by biology (dsDNA symmetry — not fitted); the seed only needs the sample
    # sizes to size the sampling part of the floor ¼·(1/N + ω). N_gdna=0 (a gDNA-free library) ⇒ 1/N_gdna → ∞ ⇒
    # the strand seed is gated off (nothing to distinguish RNA from).
    if inj is not None and inj.n_gdna_obs is not None:
        n_gdna_obs = float(
            inj.n_gdna_obs
        )  # INJECTED intergenic gDNA sample size (toy intergenic is sparse)
    else:
        _inter = coarse_type_array(np.asarray(region_arrays.signature)) == 0
        _gpos = float(
            np.sum(np.asarray(substrate.contained.n_unspliced_pos, dtype=np.float64)[_inter])
        )
        _gneg = float(
            np.sum(np.asarray(substrate.contained.n_unspliced_neg, dtype=np.float64)[_inter])
        )
        n_gdna_obs = _gpos + _gneg
    # n_rna_obs is set above (injected or from the strand-balance fit).

    def _init_belief():
        return init_beliefs(
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

    belief = _init_belief()
    # The gDNA INTRON FACTORY λ-factor (`docs/CARRY_FORWARD.md`): peel confident gDNA
    # from intron nodes against the intergenic background, BEFORE the pass-0 solve. Built ONCE (belief-free —
    # only the intron count vs the background), applied in every sweep below. ``None`` (disabled / no
    # informative background / no introns) ⇒ byte-identical to the pre-factory pass-0.
    # INJECTED intergenic intron-factory background overrides the internal (toy-sparse) fit.
    intron_background = (
        inj.intron_background
        if (inj is not None and inj.intron_background is not None)
        else (
            fit_intron_background(substrate, region_arrays, region_eff_len, include_introns=False)
            if config.intron_factory
            else None
        )
    )
    intron_prior = (
        _build_intron_prior(
            chain, substrate, region_arrays, region_eff_len, config, bg=intron_background
        )
        if (config.intron_factory and intron_background is not None)
        else None
    )
    # Message precision is entirely SELF-CONTAINED in the sweep: the source's own honest belief precision
    # (strand + count, from `node_init.build_node_init`), degraded by the reframe's scale variance
    # (σ²_transfer = Var(log r)) and the DerSimonian–Laird composition-mismatch b̂² — all derived from counts and
    # effective lengths inside the pass (`message_variance_derivation.md`). There is nothing to fit here.

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
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            n_gdna_obs=n_gdna_obs,
            n_rna_obs=n_rna_obs,
            n_grid=config.sweep_n_grid,
            logodds_window=config.sweep_logodds_window,
            n_tilt=config.sweep_n_tilt,
            n_grid_ss=config.sweep_n_grid_single_strand,
            gdna_prior=prior,
            intron_prior=intron_prior,
            _capture=capture,
        )
        if capture is not None:
            _debug["capture"] = capture
        return out

    # THE ENRICHMENT NPMLE — fit ONCE on ALL nodes' TOTAL unspliced density (belief-free). It models the
    # hybrid-capture ENRICHMENT/DEPLETION landscape, NOT composition: a total-density prior is
    # composition-vacuous (count-zero-information — docs/CARRY_FORWARD.md §2/§5), so it is NEVER fed to the
    # composition (gDNA) arm. Its old second role — supplying the message σ²_transfer by projection — is
    # RETIRED (that was a density-uniformity proxy, invalid under capture, and identically 0 in pass-0); the
    # solver now derives σ²_transfer itself. What remains is the QC report's P(ρ) landscape + the toy-injection
    # substrate, so it is fit here and consumed below, never inside the sweep.
    mass_global, eff_global = node_global_geometry(chain, geometry)
    if inj is not None and inj.enrichment_prior is not None:
        # INJECTED population enrichment landscape — a toy has too few nodes to resolve enriched vs depleted
        # modes. (Scale note: the toy's densities must be generated at the reference library's depth so its
        # nodes project onto the right cells of this absolute log-density landscape.)
        enrichment_prior = inj.enrichment_prior
    else:
        enrichment_prior = DensityNPMLE.fit(
            mass_global, eff_global, bandwidth=config.npmle_bandwidth
        )
    # PHASE 1 — the INITIAL solve is PRIOR-FREE of the DNA composition prior: the inert Beta(½,½) reference
    # alone (``gdna_prior=None``) + the strand likelihood + the belief-free forward-backward messages.
    # Single-strand nodes self-solve from strand; unstranded
    # AMBIG nodes are grounded only by the messages here (the two-root DNA ambiguity, docs/CARRY_FORWARD.md §4,
    # is resolved by the DECONVOLVED-gDNA hyperprior in Phase 2 — fit on this solve's peeled DNA, then a refit).
    belief = _sweep(None)
    belief_pass0 = (
        belief  # the initial (prior-free) solve — kept for the refit before/after (movie / debug)
    )
    logger.debug(
        "calibration: PHASE 1 prior-free initial solve (enrichment NPMLE %d cells → QC only)",
        enrichment_prior.n_cells,
    )

    # PHASE 2 — the DECONVOLVED-gDNA hyperprior REFIT (docs/CARRY_FORWARD.md §4/§5). Fit the gDNA-rate NPMLE on
    # the initial solve's peeled gDNA, then RE-SOLVE with it as the composition arm — resolving the two-root DNA
    # ambiguity the prior-free pass leaves at unstranded AMBIG nodes. Repeated ``calib_refit_iters`` times.
    # ANCHORED, EXTREMELY WEAK. The aggregate DNA-background reference (`ρ_bg`, pooled pure intergenic/intron —
    # belief-free) is the refit floor; ``None`` when disabled.
    if inj is not None and inj.background is not None:
        background = (
            inj.background
        )  # INJECTED aggregate ρ_bg (pooled pure intergenic/intron — population-scale)
    elif config.background_floor:
        background = measure_background(
            substrate,
            region_arrays,
            region_eff_len,
            include_introns=config.background_include_introns,
            robust_trim_mad=config.background_robust_trim_mad,
        )
    else:
        background = None
    gdna_hyperprior: GdnaLandscape | None = None
    for it in range(int(config.calib_refit_iters)):
        gdna_hyperprior = _fit_gdna_hyperprior(
            chain,
            belief,
            statics,
            region_arrays,
            mass_global,
            eff_global,
            strength=config.gdna_prior_strength,
        )
        if gdna_hyperprior is None:
            break
        # FULL reset, then re-solve WITH the prior: nothing from pass-0 survives into the re-solve except
        # the fitted landscape itself, so an over-confident node cannot refuse to budge when the prior lands.
        belief = _init_belief()
        belief = _sweep(gdna_hyperprior)
        logger.debug(
            "calibration: PHASE 2 gDNA-hyperprior refit %d/%d (%d training nodes)",
            it + 1,
            config.calib_refit_iters,
            gdna_hyperprior.n_train,
        )

    regions = chain_region_deconv(chain, belief, substrate)
    left, right = chain_boundary_side_deconv(chain, belief, substrate)

    if (
        _debug is not None
    ):  # inert diagnostic hook — the solved chain internals (Phase-2 substrate + plots)
        _debug.update(
            chain=chain,
            belief=belief,  # the FINAL belief (refit if calib_refit_iters>0, else the initial solve)
            belief_pass0=belief_pass0,  # the initial prior-free solve (the refit before/after frame)
            geometry=geometry,
            statics=statics,
            substrate=substrate,
            boundary_substrate=boundary_substrate,
            region_arrays=region_arrays,
            gdna_prior=enrichment_prior,  # the ENRICHMENT NPMLE (QC landscape / injection substrate)
            gdna_hyperprior=gdna_hyperprior,  # the DECONVOLVED-gDNA composition hyperprior (None if no refit)
            rna_sense_frac=rna_sense_frac,
            region_eff_len=region_eff_len,
            boundary_eff_len=boundary_eff_len,
            # the fitted-or-injected population priors — extract from a population scenario, inject into a toy
            calibration_priors=InjectedCalibrationPriors(
                rna_sense_frac=rna_sense_frac,
                n_rna_obs=n_rna_obs,
                n_gdna_obs=n_gdna_obs,
                gdna_strand_overdispersion=gdna_strand_overdispersion,
                rna_strand_overdispersion=rna_strand_overdispersion,
                enrichment_prior=enrichment_prior,
                intron_background=intron_background,
                background=background,
            ),
        )

    # Report-facing diagnostics: the fitted gDNA hyperprior P(ρ) (bimodal ⇒ capture enrichment). Consumed by
    # the QC report, never by the EM.
    if diagnostics_out is not None:
        from .diagnostics import CalibrationDiagnostics

        diagnostics_out["calibration"] = CalibrationDiagnostics.from_prior(enrichment_prior)

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
        "rna_strand_overdispersion=%.4g (%d junctions, %d frags%s) "
        "[boundary-spliced sense_frac=%.3f vs κ=%.3f]",
        result.n_regions,
        result.gdna_density_global,
        rna_sense_frac,
        gdna_strand_overdispersion,
        _gd_seed[0],
        _gd_seed[1],
        ", FALLBACK" if _gd_seed[2] else ("" if _gd_seed[0] >= 0 else ", INJECTED"),
        rna_strand_overdispersion,
        _rna_seed[0],
        _rna_seed[1],
        ", FALLBACK" if _rna_seed[2] else ("" if _rna_seed[0] >= 0 else ", INJECTED"),
        boundary_sense_frac,
        rna_sense_frac,
    )
    return result


__all__ = ["calibrate"]
