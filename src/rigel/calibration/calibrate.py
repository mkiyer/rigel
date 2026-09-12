"""calibrate() — the calibrator: deconvolve each object's UNSPLICED count into ``(RNA₊, RNA₋, gDNA)``.

This is LAYER 7's assembly stage. It owns the ORDER the pieces run in and the wiring between them;
every number it reports is produced by a module below it, and it computes no model of its own.

Calibration models RNA vs gDNA only (the per-locus EM separates transcripts downstream). An object's
composition has three sources: the STRAND LIKELIHOOD (the Beta-Binomial tilt of the per-strand
counts, the only intrinsic signal, entering as its overdispersed Fisher information rather than as a
raw count); the MESSAGES from its chain neighbours (:mod:`.messages`, each hop priced inside the
sweep from the two nodes' own counts); and the POPULATION gDNA PRIOR (the intergenic-only background
pool plus the phase-2 density landscape — both `fit_intron_background` call sites pass
``include_introns=False``, because an intron-inclusive pool is inflated by nascent RNA worst exactly
where gDNA is scarce). The solver is the belief-propagation SWEEP over the ``N E N E … N`` chain
(:mod:`rigel.calibration.sweep`)::

    substrate  (five populations on three axes)
      -> build chain + geometry + statics      (the geometry owns EVERY divisor)
      -> strand balance: rna_sense_frac (κ)
      -> count_observable_masks -> fit gDNA / RNA strand Beta-Binomial overdispersions (seed)
      -> signature-binary init (G1/G2/G3)
      -> PASS 1 solve_chain (no fitted prior): ONE forward + ONE backward pass — each object integrates
           its strand likelihood, the intron factory and its neighbours' messages
      -> fit the phase-2 gDNA-density landscape on the pass-1 belief
      -> PASS 2 solve_chain (the landscape added per object) -> the converged per-object pie
      -> chain_region_deconv  -> per-REGION gDNA / RNA contained mass
      -> chain_boundary_deconv  -> per-BOUNDARY gDNA / RNA crossing mass, for the per-locus prior
      -> gdna_density_global (the library-average density QC scalar)

⛔ THE GEOMETRY IS BUILT FIRST, AND IT OWNS EVERY DIVISOR. ``build_region_geometry`` produces the
per-slot ``eff_gdna``/``eff_rna`` before anything reads a count, and everything downstream — the
density clue, the background pool, the intron factory, the result's two supports — reads THAT array.
Computing a second length model anywhere below is how a consumer comes to divide by a length the
solver never used.

A contiguous boundary is a 0-bp boundary with one count and one divisor, ``crossing_eff_length``;
there is no per-face machinery and no ½. A zero-gDNA library (``gdna_density_global == 0``, per-object
gDNA mass ``0``) is a valid, graceful output.

Known bias (`TRAPS: prove-the-substrate`): the RNA half of an unspliced crossing takes
``UNBOUNDED_REACH`` rather than its transcript's real remaining length, which over-calls gDNA
genome-wide and worst in the last region before a polyA site. SpliceJunction boundaries DO take their
real exonic reach.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

import numpy as np

from .messages.silent import SilentPolicy
from .messages.transfer import TransferPolicy
from .region_chain import BOUNDARY, REGION
from .region_geometry import (
    build_region_geometry,
    build_region_statics,
    init_beliefs,
    region_gdna_geometry,
)
from .sweep import MessageMemo, chain_boundary_deconv, chain_region_deconv, solve_chain
from .density_model import count_observable_masks
from .derive import gdna_density_global
from .errors import CalibrationStrandError
from .density_deconv import (
    GdnaBackground,
    density_lambda_factor,
    fit_intron_background,
)
from .abundance_landscape import AbundanceLandscape, fit_abundance_landscape
from .total_abundance import (
    build_region_wall_mask,
    region_counts_and_exposure,
    w_max_from_deposited_lengths,
)
from .gdna_strand import (
    _MAX_OVERDISPERSION,
    reconcile_overdispersions,
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_sj_table,
)
from .region_chain import build_region_chain
from .result import CalibrationResult
from .landscape import DensityLandscape, fit_landscape
from .signature import RegionType, coarse_type_array
from .simplex_logodds import _logodds_grid
from .strand_balance import fit_strand_balance
from .substrate import CalibrationSubstrate
from ..types import Strand

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays
    from .splice_graph import SpliceJunctionGeometry

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class InjectedCalibrationPriors:
    """Population-scale calibration priors — the objects that require genome-scale (or many-gene) data to fit and
    are physically **directly observable** (no deconvolution / no solving): the RNA strand balance, the strand
    Beta-Binomial overdispersions, the strand-Fisher noise-floor sample sizes, the intergenic
    intron-factory background, and the pre-solve TOTAL-density landscape.

    A tiny (single-transcript) toy CANNOT fit these — so :func:`calibrate` accepts them pre-fit from a
    population scenario and injects them, letting the toy provide only the controlled per-region GEOMETRY. Every
    field is optional; ``None`` ⇒ fit that prior internally (the default, byte-identical). ``calibrate`` also
    stashes the fitted-or-injected bundle in ``_debug["calibration_priors"]`` so a population scenario's fitted
    priors can be extracted and re-injected into a toy (`scripts/debug/toy_inject.py`)."""

    rna_sense_frac: float | None = None
    n_rna_obs: float | None = None
    n_gdna_obs: float | None = None
    gdna_strand_overdispersion: float | None = None
    rna_strand_overdispersion: float | None = None
    intron_background: GdnaBackground | None = None
    #: the pre-pass-0 TOTAL-density field + mode census — population-scale (a toy cannot fit a
    #: landscape from a handful of regions), injectable exactly like the enrichment prior it is
    #: planned to replace. ``None`` ⇒ fit internally when the config asks, else absent.
    abundance_landscape: AbundanceLandscape | None = None


def _empty_sj_geometry() -> "SpliceJunctionGeometry":
    """A sj axis with no rows — a graph whose references are all single-exon.

    Legal, and not the same as "no sj flux": a payload scanned against such a graph has
    ``n_sj == 0``, so an empty axis is the only consistent input and the alignment check below will say
    so if it is not.
    """
    from .splice_graph import SpliceJunctionGeometry

    e_i = np.zeros(0, dtype=np.int64)
    return SpliceJunctionGeometry(
        src_region=e_i,
        dst_region=e_i.copy(),
        strand=np.zeros(0, dtype=np.int8),
        reach_lo=np.zeros(0, dtype=np.float64),
        reach_hi=np.zeros(0, dtype=np.float64),
    )


def _project_eff(chain, eff_slots, payload) -> tuple[np.ndarray, np.ndarray]:
    """Split a per-SLOT divisor back onto the region and contiguous-boundary axes.

    ⛔ A projection, not a recomputation. ``build_region_geometry`` already applied
    ``contained_eff_length`` at REGION slots and ``crossing_eff_length`` at BOUNDARY slots; calling
    those again here would put two implementations of one quantity in the tree, and two
    implementations of one divisor are how a factor of ½ survives unnoticed. Whatever the solver
    divided by is what the result reports.

    It takes the ARRAY, not the geometry, so ONE function serves both populations rather than a
    second copy that drifts. Only the gDNA projection has a live consumer — the prior divides by
    nothing on the mass path — and the RNA one is retained deliberately (`result.py` carries the
    reasoning).
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    eff = np.asarray(eff_slots, dtype=np.float64)
    is_region, is_boundary = kind == REGION, kind == BOUNDARY
    region_eff = np.zeros(int(payload.n_regions), dtype=np.float64)
    boundary_eff = np.zeros(int(payload.n_boundaries), dtype=np.float64)
    region_eff[obj[is_region]] = eff[is_region]
    boundary_eff[obj[is_boundary]] = eff[is_boundary]
    return region_eff, boundary_eff


def _build_intron_prior(chain, substrate, region_arrays, region_eff_len, config, bg=None):
    """The gDNA intron factory λ-factor per chain slot.

    Fits the intergenic-background NegBinom (`fit_intron_background`) and tabulates, for each INTRON REGION
    slot, ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` over the σ(λ) solve grid → an ``(n_slots, K)`` array,
    ZERO on every other slot (a no-op there). Returns ``None`` when the factory is disabled, the
    background pool is uninformative, or there are no intron regions — in which case the sweep is
    byte-identical to the pre-factory path. gDNA is strand-symmetric, so this factor lives purely on
    ``λ`` (deconvolves gDNA; the residual RNA's tilt is left to the solver), and is consumed identically by the
    single-strand and AMBIG per-region solves.

    BOUNDARY slots are zero, structurally. The factor scores a CONTAINED count against a CONTAINED
    support; a boundary's count is a crossing and its divisor a different formula, so applying the same
    NegBinom there would be scoring one frame's evidence against another frame's support.

    ``bg`` (an injected population :class:`GdnaBackground`) overrides the internal fit — a tiny toy's own
    intergenic pool is too sparse to fit the background the introns are deconvolved against."""
    if bg is None:
        bg = fit_intron_background(substrate, region_arrays, region_eff_len, include_introns=False)
    if not bg.informative:
        return None
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(
        np.int64
    )  # 0/1/2 per REGION
    R = rtype.shape[0]
    is_intron = (kind == REGION) & (rtype[np.clip(idx, 0, R - 1)] == 1)  # INTRON == 1
    if not bool(is_intron.any()):
        return None
    _, fg = _logodds_grid(int(config.sweep_n_grid), float(config.sweep_logodds_window))
    prior = np.zeros((kind.shape[0], fg.shape[0]), dtype=np.float64)
    ridx = idx[is_intron]
    # GENOME-strand columns summed: gDNA is strand-symmetric, so the deconvolution is against a total rate.
    count = np.asarray(substrate.region_contained.count, dtype=np.float64).sum(axis=1)[ridx]
    eff_g = np.asarray(region_eff_len, dtype=np.float64)[ridx]
    prior[is_intron] = density_lambda_factor(bg, count, eff_g, fg)
    return prior


#: Minimum training regions for a hyperprior fit — below this the population is not a population.
_MIN_TRAIN = 5


def _fit_gdna_hyperprior(
    chain, belief, statics, region_arrays, mass_global, eff_global, *, strength, prev=None
):
    """Select the training substrate from the chain and fit the :class:`DensityLandscape` on the
    initial solve's deconvolved gDNA — the composition (gDNA) arm of ψ for the phase-2 refit.
    ``None`` if it cannot be fit.

    This affects only the PRIOR fit, never the solve's own messages.

    The substrate is the whole of this function's job; the estimator in :mod:`.landscape` is
    component-agnostic and knows nothing about gDNA. Four axes decide membership, and conflating them
    is the mistake to avoid:

    * circularity → structural EXCLUSION. AMBIG regions are out: they are the two-root ambiguity the
      prior exists to resolve, so training on them would have it predict itself. The empirical case
      alongside this is stratum-dependent, so the argument carrying the decision is the circularity.
    * identifiability → structural INCLUSION. A live region that sequenced no unspliced mass has gDNA
      density ``0`` for every ``f_g``: "gDNA is absent here" is the strongest depletion evidence there
      is, and this zero-count anchor is what grounds the depleted mode. Dropping it is a large,
      measured loss, worst on zero-gDNA libraries.
    * precision → a continuous WEIGHT, never admission (`landscape._reliability`).
    * geometry → BOUNDARIES ARE EXCLUDED. They cross rather than contain, are about as numerous as
      regions but far less often truly enriched, and their two-flank mixture fills the valley between
      the two true modes.

    ⛔ Admitting AMBIG into the FINAL fit — even in the non-circular form where the admitting fit
    trains on AMBIG estimates from a prior that never saw AMBIG — was implemented, A/B'd over every
    condition and REFUTED: worse on most of them, with the zero-gDNA false-positive guard regressing
    too. It grew the training set, so it was active rather than inert. Do not re-propose it without
    new evidence.
    """
    isr = np.asarray(chain.kind) == REGION
    fp = np.asarray(statics.free_pos, dtype=bool)
    fn = np.asarray(statics.free_neg, dtype=bool)
    rtype = coarse_type_array(np.asarray(region_arrays.signature))
    ridx = np.clip(np.asarray(chain.obj_idx, dtype=np.int64), 0, rtype.shape[0] - 1)
    expressed = isr & (eff_global > 1.0e-9) & (mass_global > 1.0e-12)
    # the zero-count structural anchor: an intergenic or intronic region that sequenced no unspliced mass
    anchor = (
        isr & (eff_global > 1.0e-9) & (mass_global <= 1.0e-12) & (rtype[ridx] != RegionType.EXON)
    )
    sel = expressed & ((fp ^ fn) | (~fp & ~fn))
    # ⛔ A SLOT WHOSE ONLY EVIDENCE IS A BOUND, OR WHICH HAS NONE, DOES NOT TRAIN THE PRIOR
    # (`RegionBelief.informed`). Its value is where the prior put it — outright, or inside the
    # half-line a level or ceiling admits — so re-fitting on it re-seeds the landscape's tail at the
    # density its own total implies. ⛔ Do not soften the cut by DISCOUNTING the delivered rows
    # instead (a one-sided row as a bound, the own-evidence variance as a weight): those readings lose
    # the deferred stratum badly, because those rows are the enriched mode's witness under capture.
    # The anchor trains regardless: it is a structural statement, not a solve.
    # THE GRID IS THE CONSUMERS' DOMAIN — EVERY slot the prior is read at, regions and boundaries
    # alike — so the composition cut above changes which kernels are summed and never the axis
    # (`fit_landscape`'s ``domain``). ⛔ Without that, a gDNA-free library whose exons are all blind
    # trains on its anchors alone, the grid collapses to the floor, and every exon reads a flat prior.
    domain_sel = np.asarray(eff_global, dtype=np.float64) > 1.0e-9
    # ⛔ The substrate guard measures the population the ANNOTATION admits, not the population left
    # after the composition cut: "is there enough of a population to fit a prior for" is not the
    # training set's question. Guarding the cut population instead refuses the refit on a gDNA-free
    # library whose only survivors are a handful of anchors, which then leaves its AMBIG regions at
    # their prior-free share. `fit_landscape` still refuses a training set under two.
    if int((sel | anchor).sum()) < _MIN_TRAIN:
        return None
    if belief.informed is not None:
        sel &= np.asarray(belief.informed, dtype=bool)
    sel |= anchor
    mass = np.asarray(mass_global, dtype=np.float64)[sel]
    return fit_landscape(
        np.asarray(belief.f_g, dtype=np.float64)[sel] * mass,
        mass,
        np.asarray(eff_global, dtype=np.float64)[sel],
        np.asarray(belief.var_gdna, dtype=np.float64)[sel],
        anchor=anchor[sel],
        strength=strength,
        # the previous refit's landscape: the E-step on the location-free kernels (`landscape._estep_kernels`)
        prev=prev,
        domain=(
            np.asarray(mass_global, dtype=np.float64)[domain_sel],
            np.asarray(eff_global, dtype=np.float64)[domain_sel],
        ),
    )


def _scaled_grid(n: int, window: float, required: float) -> int:
    """Grid points at the widened bracket, holding the lattice spacing ``dlam = 2L/(n−1)`` FIXED.

    ⛔ The bracket and the RESOLUTION are two different knobs, and an arm that moves both at once is
    uninterpretable: widening ``L`` at a FIXED ``n`` coarsens ``dlam``, and the result then reverses
    with the bracket instead of saturating. Saturation is what distinguishes a truncated bracket from
    an improper ψ, so the spacing must not move and ``n`` is linear in ``L``, exactly.

    The rounding is the only slack: ``dlam`` is preserved to within half a grid point, which is the
    best a discrete lattice can do.
    """
    return int(round(1.0 + (n - 1) * (required / window)))


def calibrate(
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    gdna_fl_pmf: "np.ndarray",
    rna_fl_pmf: "np.ndarray",
    config: "CalibrationConfig",
    sj: "SpliceJunctionGeometry | None" = None,
    boundary_rna_reach=None,
    _debug: dict | None = None,
    diagnostics_out: dict | None = None,
    injected_priors: "InjectedCalibrationPriors | None" = None,
    boundary_flags: "np.ndarray | None" = None,
    mature_walls=None,
    boundary_reach=None,
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per object, then derive gdna_density_global.

    Runs the belief-propagation sweep (a single forward-backward pass per phase, resolving the
    per-object pie); see the module docstring for the data flow. ``gdna_density_global`` may be ``0``
    (a zero-gDNA library) and an object's deconvolved gDNA mass may be ``0`` (a pure-RNA object); both
    are valid, graceful outputs — not failures.

    ``sj`` is the splice graph's sj axis
    (:func:`~rigel.calibration.splice_graph.build_sj_geometry_arrays`), in the accumulator's own
    sj slot order — where each sj attaches, its transcript strand, and its exonic reach.
    ``None`` means "this library's graph has no sj boundaries", which is legal (a single-exon-only
    reference) and is NOT the same as "no sj flux".

    ``boundary_flags`` is the splice graph's per-contiguous-boundary structural bits
    (:func:`~rigel.calibration.splice_graph.build_boundary_flags_array`), carried onto the chain as
    ``RegionStatics.boundary_flags``. ``None`` (the default) leaves them zero.

    ``mature_walls`` / ``boundary_reach`` are the two annotation-only WALL inputs the MEASURED-TOTAL
    exposure needs (:func:`~rigel.calibration.splice_graph.build_mature_wall_distances` and
    :func:`~rigel.calibration.splice_graph.build_contiguous_boundary_reach_arrays`). ⛔ They are
    consulted ONLY when ``config.background_abundance == "measured_total"``, and that setting REFUSES
    to run without them rather than silently falling back — a background rate that quietly changed
    estimator would be the worst of both.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    inj = injected_priors  # population-scale priors to inject in place of the internal (toy-untrustworthy) fits
    sj = _empty_sj_geometry() if sj is None else sj
    if int(sj.n_sj) != int(substrate.n_sj):
        raise ValueError(
            f"the sj axis has {int(sj.n_sj)} boundaries but the payload has "
            f"{int(substrate.n_sj)}. Build it with "
            "splice_graph.build_sj_geometry_arrays(index) against the SAME index the payload "
            "was scanned on — a sj axis addressing a different graph would place every splice "
            "on the wrong boundary."
        )

    # THE CHAIN AND ITS GEOMETRY COME FIRST, and the geometry owns every divisor from here down:
    # `eff_gdna`/`eff_rna` are the CONTAINED placements at a REGION slot and the CROSSING placements at
    # a BOUNDARY slot, one rule, one array. Nothing below computes a second length model.
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geometry = build_region_geometry(
        chain, substrate, region_arrays, sj, gdna_fl_pmf, rna_fl_pmf, boundary_rna_reach
    )
    statics = build_region_statics(chain, region_arrays, boundary_flags)

    # The result's two gDNA supports, PROJECTED off the geometry rather than recomputed — so the number
    # `priors` divides by is byte-identically the one the solver divided by. Two implementations
    # of one quantity is how they come to disagree.
    region_eff_gdna, boundary_eff_gdna = _project_eff(chain, geometry.eff_gdna, payload)

    # THE MEASURED-TOTAL (counts, exposure) PAIR for the pooled gDNA background estimators, built
    # ONCE here and handed to both. ``None`` under the shipped default, and then every consumer takes
    # its own contained pair instead.
    # ⛔ It REFUSES rather than falling back: a background rate that silently changed estimator because
    # an argument was missing is worse than either estimator.
    background_pair = None
    if config.background_abundance == "measured_total":
        if mature_walls is None or boundary_reach is None:
            raise ValueError(
                "CalibrationConfig.background_abundance = 'measured_total' needs the wall inputs: "
                "pass mature_walls=build_mature_wall_distances(index, region_arrays) and "
                "boundary_reach=build_contiguous_boundary_reach_arrays(index) (both are in "
                "scan_cache.index_derived_inputs). Refusing rather than falling back to the contained "
                "pair, which would change the background rate without saying so."
            )
        _wall_mask = build_region_wall_mask(
            region_arrays,
            mature_walls,
            boundary_reach[0],
            boundary_reach[1],
            w_max=w_max_from_deposited_lengths(payload.deposited_lengths),
        )
        _bg_counts, _bg_exposure, _ = region_counts_and_exposure(
            substrate, region_arrays, _wall_mask
        )
        background_pair = (_bg_counts, _bg_exposure)

    # THE ABUNDANCE LANDSCAPE — the pre-pass-0 TOTAL-density field + mode census, fitted at INIT
    # from counts and lengths only, so it is not circular with anything solved below. It is a QC and
    # injection surface; nothing in the solve reads it.
    abundance_landscape = None
    if config.abundance_landscape:
        if inj is not None and inj.abundance_landscape is not None:
            abundance_landscape = inj.abundance_landscape
        elif mature_walls is None or boundary_reach is None:
            # ⛔ SKIPPED, LOUDLY — it must not raise. The flag is ON by default, this object is the
            # sole source of the QC report's density panel, and many unit and toy callers legitimately
            # have no wall arrays and never wanted a panel. The alternative to skipping was never to
            # fit on unmasked totals — nothing fits unmasked either way — so the choice is between a
            # missing panel and a refusal, and a missing panel is right here. It is NOT a silent
            # fallback: the object stays None, never a quietly-different estimate, and the skip is
            # logged. `background_abundance` above KEEPS its refusal, because that pair feeds ψ while
            # this one is read only by the report and the debug bundle.
            logger.warning(
                "calibration: abundance_landscape is enabled but the wall inputs are missing "
                "(mature_walls / boundary_reach, both in scan_cache.index_derived_inputs) — skipping "
                "the total-density landscape, so the QC density panel will be omitted. Nothing in the "
                "solve reads it, so no solved number changes."
            )
        else:
            _al_mask = build_region_wall_mask(
                region_arrays,
                mature_walls,
                boundary_reach[0],
                boundary_reach[1],
                w_max=w_max_from_deposited_lengths(payload.deposited_lengths),
            )
            abundance_landscape = fit_abundance_landscape(substrate, region_arrays, _al_mask)
    # And the RNA twin, on the same two axes — the RNA population's own opportunity per object.
    # It has no consumer in the solve: the prior is a conserved FRAGMENT COUNT and divides by nothing
    # on the mass path. Kept because it is the RNA divisor any density-based prior needs and because
    # it is byte-identically the opportunity the solver used — see `result.py`.
    region_eff_rna, boundary_eff_rna = _project_eff(chain, geometry.eff_rna, payload)

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

    # The gDNA strand fit's SEED SELECTOR: count-observability, straight off the signature. The
    # away-half moment (`gdna_strand`) needs no seed weight, so two masks are the whole input — there
    # is no local density imputation behind them.
    _region_obs, _boundary_obs = count_observable_masks(
        np.asarray(region_arrays.signature), np.asarray(region_arrays.ref_id)
    )

    # Strand-module parameters — the two Beta-Binomial overdispersions.
    # RNA FIRST (mean κ; fitted from the PER-SJ SJ strand table, certified pure RNA): it is the gDNA
    # fit's fallback. THEN gDNA (mean ½) by the AWAY-HALF moment over every genic count- and
    # strand-observable object — intron regions, exon|intron and gene-edge boundaries — which is
    # unbiased under any RNA content of the seeds (`gdna_strand`'s lemma), so no seed is weighted and
    # no class is asserted pure: the fit must hold regardless of the reference transcriptome.
    # Intergenic and AMBIG objects cannot be oriented and are out. The gDNA fit is the raw pooled
    # moment clipped to the physical support — no location prior (gate
    # `tests/calibration/test_gdna_strand_fit.py`).
    # (n_seed_regions, n_seed_frags, fallback, effective_seeds, raw_od) — QC log only; -1 = injected
    _gd_seed = (-1, -1, False, float("nan"), float("nan"))
    _rna_seed = (-1, -1, False, float("nan"))
    if inj is not None and inj.rna_strand_overdispersion is not None:
        rna_strand_overdispersion = float(inj.rna_strand_overdispersion)
    else:
        rna_strand = fit_rna_strand_from_sj_table(
            strand_model.sj_table,
            rna_sense_frac=rna_sense_frac,
        )
        rna_strand_overdispersion = rna_strand.rna_strand_overdispersion
        _rna_seed = (
            rna_strand.n_seed_regions,
            rna_strand.n_seed_fragments,
            rna_strand.fallback_used,
            rna_strand.raw_overdispersion,
        )
    if inj is not None and inj.gdna_strand_overdispersion is not None:
        gdna_strand_overdispersion = float(inj.gdna_strand_overdispersion)
    else:
        gdna_strand = fit_gdna_strand_from_substrate(
            substrate,
            region_arrays,
            region_count_observable=_region_obs,
            boundary_count_observable=_boundary_obs,
            rna_sense_frac=rna_sense_frac,
        )
        gdna_strand_overdispersion = gdna_strand.gdna_strand_overdispersion
        _gd_seed = (
            gdna_strand.n_seed_regions,
            gdna_strand.n_seed_fragments,
            gdna_strand.fallback_used,
            gdna_strand.effective_seeds,
            gdna_strand.raw_overdispersion,
        )

    # RECONCILE THE TWO COMPONENTS, with no conjured shrinkage target. The weaker-measured dispersion
    # shrinks toward the better-measured one, weighted by their
    # own null informations; neither is pulled toward a conjured number, and with neither measured the two
    # coincide at the ceiling, which is what leaves the strand channel uninformative rather than confident.
    # ⛔ ONLY when NEITHER is injected: an injected value is the arm's whole point, and letting the other
    # component shrink toward it would silently change what every od-injection instrument measures.
    if inj is None or (
        inj.rna_strand_overdispersion is None and inj.gdna_strand_overdispersion is None
    ):
        rna_strand_overdispersion, gdna_strand_overdispersion = reconcile_overdispersions(
            rna_strand.raw_overdispersion,
            rna_strand.information,
            gdna_strand.raw_overdispersion,
            gdna_strand.information,
        )

    # Strand-Fisher noise-floor SAMPLE SIZES (the sweep's τ seed): N_gdna (gDNA-eligible unspliced fragments in
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
        n_gdna_obs = float(
            np.asarray(substrate.region_contained.count, dtype=np.float64)[_inter].sum()
        )
    # n_rna_obs is set above (injected or from the strand-balance fit).

    def _init_belief():
        return init_beliefs(
            chain,
            geometry,
            statics,
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            n_grid=config.sweep_n_grid,
            n_grid_ss=config.sweep_n_grid_single_strand,
            logodds_window=config.sweep_logodds_window,
        )

    belief = _init_belief()
    # The gDNA INTRON FACTORY λ-factor: deconvolve confident gDNA
    # from intron regions against the intergenic background, BEFORE the pass-0 solve. Built ONCE (belief-free —
    # only the intron count vs the background), applied in every sweep below. ``None`` (disabled / no
    # informative background / no introns) ⇒ byte-identical to the pre-factory pass-0.
    # INJECTED intergenic intron-factory background overrides the internal (toy-sparse) fit.
    intron_background = (
        inj.intron_background
        if (inj is not None and inj.intron_background is not None)
        else (
            fit_intron_background(
                substrate,
                region_arrays,
                region_eff_gdna,
                include_introns=False,
                counts_exposure=background_pair,
            )
            if config.intron_factory
            else None
        )
    )
    intron_prior = (
        _build_intron_prior(
            chain, substrate, region_arrays, region_eff_gdna, config, bg=intron_background
        )
        if (config.intron_factory and intron_background is not None)
        else None
    )
    # ⛔⛔ **THE INTRON FACTORY'S λ-FACTOR IS EVALUATED *ON* THE SOLVE GRID, so it is a function of
    # (n_grid, L) and must be REBUILT when the bracket moves — not regridded.** `_regrid_global`
    # interpolates a prior BETWEEN two grids of the same L; there is no such map onto a WIDER domain,
    # because the factor was never evaluated out there; regridding onto a wider domain raises
    # `IndexError`, which is the honest failure. Memoised, because only two brackets ever occur in one
    # call (the configured one and the derived one) and `_build_intron_prior` is a per-slot NegBinom
    # over the whole grid.
    _intron_priors = {(int(config.sweep_n_grid), float(config.sweep_logodds_window)): intron_prior}

    def _intron_prior_at(n_grid: int, window: float):
        key = (int(n_grid), float(window))
        if key not in _intron_priors:
            _intron_priors[key] = _build_intron_prior(
                chain,
                substrate,
                region_arrays,
                region_eff_gdna,
                replace(config, sweep_n_grid=int(n_grid), sweep_logodds_window=float(window)),
                bg=intron_background,
            )
        return _intron_priors[key]

    # Every message's price is SELF-CONTAINED in the sweep: each hop charges the two nodes' counting
    # and the pair's own disagreement, all derived from counts and opportunities inside the pass.
    # There is nothing to fit here.

    # When ``_debug`` is on, the LAST sweep also fills ``_debug["capture"]`` with the per-region
    # message internals. Inert in production.
    # ONE policy instance for every phase; the intron factory's rows reach it on each sweep's context
    # (``factory_rows``, the same memoised array the sweep adds as ψ's λ-factor).
    # ⛔ THE POLICY NAME MUST SELECT THE POLICY. An unreadable knob is worse than no knob: an
    # arm that silently runs a different policy than it names is a benchmark that cannot be
    # trusted, so an unknown name RAISES here.
    if not config.message_propagation or config.message_policy == "silent":
        policy = SilentPolicy()
    elif config.message_policy == "transfer":
        # the composition-transfer policy on the two-phase backbone. An intron's own claim is the
        # factory's row for it, which the sweep hands over on the context; `None` there (factory off /
        # uninformative background) makes it byte-identical to silence.
        policy = TransferPolicy(
            strand=(rna_sense_frac, gdna_strand_overdispersion, rna_strand_overdispersion),
        )
    else:
        raise ValueError(
            f"unknown message_policy {config.message_policy!r} — expected 'silent' or 'transfer'"
        )

    def _sweep(prior, memo=None):
        capture = {} if _debug is not None else None
        # THE λ BRACKET IS `max(the reference's floor, the fitted prior's own demand)`. ψ evaluates
        # the landscape at `log ρ = log f + log M − log E` and can only offer `f ∈ [σ(−L), σ(L)]`, so a
        # bracket narrower than the prior's own support leaves ψ with no coordinate for what the prior is
        # telling it — and the answer then depends on `L`, which `simplex_logodds` calls its own
        # acceptance test. The demand is DERIVED (`required_logodds_window`), never chosen.
        # ⛔ The prior-free Phase-1 solve has no landscape to widen for, and is L-invariant to seven
        # digits, so it keeps the floor. That is the whole of the conditional.
        # ⛔ `dlam` is held FIXED, so the grids scale with the bracket: widening `L` at a fixed grid
        # size coarsens the lattice and confounds two knobs, and the answer then reverses with `L`
        # instead of saturating.
        # The TILT axis does NOT scale. `θ` is the RNA-internal SHARE and has no bracket problem, so it
        # stays at the configured resolution, which keeps the AMBIG cube's growth roughly linear in the
        # bracket rather than quadratic.
        window = float(config.sweep_logodds_window)
        n_grid, n_grid_ss = int(config.sweep_n_grid), int(config.sweep_n_grid_single_strand)
        n_tilt = (
            config.sweep_n_tilt if config.sweep_n_tilt is not None else int(config.sweep_n_grid)
        )
        if prior is not None:
            required = prior.required_logodds_window(mass_global, eff_global)
            if required > window:
                n_grid = _scaled_grid(n_grid, window, required)
                n_grid_ss = _scaled_grid(n_grid_ss, window, required)
                window = required
                logger.debug(
                    "calibration: λ bracket %.4f (the landscape's support), n_grid %d, n_grid_ss %d, "
                    "n_tilt %d (unscaled)",
                    window,
                    n_grid,
                    n_grid_ss,
                    n_tilt,
                )
        lam_factor = _intron_prior_at(n_grid, window)
        # ⛔ ψ has NO reference location. A located reference is a prior ASSERTION at fixed strength,
        # and it becomes the entire answer wherever the strand channel is dead. The reference is the
        # symmetric Jeffreys measure; background information enters as the intron-factory λ-factor
        # above, a likelihood whose precision scales with counts.
        out = solve_chain(
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
            n_grid=n_grid,
            logodds_window=window,
            n_tilt=n_tilt,
            n_grid_ss=n_grid_ss,
            gdna_prior=prior,
            intron_prior=lam_factor,
            # Message propagation ships ON; `message_policy` selects which policy that installs —
            # see the dispatch above. The two halves of the panel are judged against different bars
            # and are never pooled: messages must WIN on unstranded data, where kappa = 1/2 zeroes the
            # strand lambda-term so a slot has no own composition evidence and a message is its only
            # source, and do minimal HARM against silence on stranded data, where a sighted slot's own
            # solve is already good.
            policy=policy,
            block_slots=config.sweep_block_slots,
            message_memo=memo,
            _capture=capture,
        )
        if capture is not None:
            _debug["capture"] = capture
        return out

    # ⛔ A TOTAL density over ONE component's opportunity model is not a composition estimate: an
    # estimator fitted on `mass / eff_gdna` over all slots is answering the wrong question, however
    # well it fits. The total-density field this module does use is the `AbundanceLandscape` above,
    # fitted on the wall-exact measured totals over each region's own LENGTH — a geometry rather than
    # a model in the divisor — and it reaches the report, never the solve.
    mass_global, eff_global = region_gdna_geometry(geometry)
    # PHASE 1 — the INITIAL solve carries no fitted composition prior: the inert Beta(½,½) reference
    # alone (``gdna_prior=None``) plus the strand likelihood and the messages. Single-strand regions
    # self-solve from strand; AMBIG regions on unstranded data are grounded only by the messages here,
    # and their two-root ambiguity is what the phase-2 hyperprior resolves.
    belief = _sweep(None)
    belief_pass0 = (
        belief  # the initial (prior-free) solve — kept for the refit before/after (movie / debug)
    )
    logger.debug(
        "calibration: PHASE 1 prior-free initial solve (abundance landscape: %s)",
        "none"
        if abundance_landscape is None
        else f"{abundance_landscape.n_train} training regions, {len(abundance_landscape.modes)} modes",
    )

    # PHASE 2 — the DECONVOLVED-gDNA hyperprior REFIT. Fit `landscape.DensityLandscape` on the
    # initial solve's deconvolved gDNA, then RE-SOLVE with it as ψ's composition arm, resolving the
    # two-root ambiguity the prior-free pass leaves at unstranded AMBIG regions. Repeated
    # ``calib_refit_iters`` times. It is anchored and extremely weak: `fit_landscape`'s own `anchor`
    # argument is what grounds it, and the only intergenic background that reaches ψ is
    # `fit_intron_background`'s.
    gdna_hyperprior: DensityLandscape | None = None
    # THE REFIT SWEEPS SHARE THEIR MESSAGE LAYER. The belief is reset before each, and the messages
    # never read the prior, so for one grid every input the layer reads is identical from refit to
    # refit; a refit pays its two ψ solves and is served the rest (`sweep.MessageMemo`, content-keyed:
    # a refit whose bracket widens changes the grid and misses). Pass 0's grid is never reused.
    memo = MessageMemo()
    for it in range(int(config.calib_refit_iters)):
        gdna_hyperprior = _fit_gdna_hyperprior(
            chain,
            belief,
            statics,
            region_arrays,
            mass_global,
            eff_global,
            strength=config.gdna_prior_strength,
            prev=gdna_hyperprior,  # None at the first fit; the fit before, after
        )
        if gdna_hyperprior is None:
            break
        # FULL reset, then re-solve WITH the prior: nothing from pass-0 survives into the re-solve except
        # the fitted landscape itself, so an over-confident region cannot refuse to budge when the prior lands.
        belief = _init_belief()
        belief = _sweep(gdna_hyperprior, memo)
        logger.debug(
            "calibration: PHASE 2 gDNA-hyperprior refit %d/%d (%d training regions)",
            it + 1,
            config.calib_refit_iters,
            gdna_hyperprior.n_train,
        )

    regions = chain_region_deconv(chain, belief, substrate)
    boundaries = chain_boundary_deconv(chain, belief, substrate)

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
            sj=sj,
            region_arrays=region_arrays,
            gdna_prior=abundance_landscape,  # the TOTAL-density landscape (QC / injection substrate)
            gdna_hyperprior=gdna_hyperprior,  # the DECONVOLVED-gDNA composition hyperprior (None if no refit)
            rna_sense_frac=rna_sense_frac,
            region_eff_gdna=region_eff_gdna,
            boundary_eff_gdna=boundary_eff_gdna,
            # the fitted-or-injected population priors — extract from a population scenario, inject into a toy
            calibration_priors=InjectedCalibrationPriors(
                rna_sense_frac=rna_sense_frac,
                n_rna_obs=n_rna_obs,
                n_gdna_obs=n_gdna_obs,
                gdna_strand_overdispersion=gdna_strand_overdispersion,
                rna_strand_overdispersion=rna_strand_overdispersion,
                intron_background=intron_background,
                abundance_landscape=abundance_landscape,
            ),
        )
        _debug["abundance_landscape"] = abundance_landscape

    # Report-facing diagnostics: the fitted gDNA hyperprior P(ρ) (bimodal ⇒ capture enrichment). Consumed by
    # the QC report, never by the EM.
    if diagnostics_out is not None:
        from .diagnostics import CalibrationDiagnostics

        # The QC density panel comes from the total-density landscape's CENSUS. ``None`` when the
        # landscape was not fit (no wall inputs) — the report then omits the panel.
        if abundance_landscape is not None:
            diagnostics_out["calibration"] = CalibrationDiagnostics.from_abundance_landscape(
                abundance_landscape
            )

    # Derive gdna_density_global (the library-average density QC scalar).
    density_global = gdna_density_global(regions, boundaries, region_eff_gdna, boundary_eff_gdna)

    # The certified-RNA crossings per BOUNDARY: molecules that crossed contiguously having spliced
    # elsewhere. ``chain_boundary_deconv`` adds the whole of this to ``rna_mass`` (rna = (1−g)·unspliced +
    # spliced), so it is exactly the spliced component of ``mass_rna_boundary``. ``assemble_priors``
    # withholds it from ``rna_prior_count`` — a spliced fragment is guaranteed-RNA in the EM (no gDNA
    # candidate), so it must not load the RNA side of the gDNA-vs-RNA *unspliced* split. ``mass_rna_boundary``
    # stays spliced-inclusive so per-boundary conservation gdna + rna = unspliced + spliced holds.
    #
    # There is no REGION twin, and that is structural: ``region_contained`` is credited only when the
    # fragment used no sj, so a region's contained population cannot hold a spliced molecule.
    mass_rna_spliced_boundary = np.asarray(substrate.boundary_spliced.count, dtype=np.float64).sum(
        axis=1
    )
    # GEOMETRY, not a split: the mean conserved fragment-mass one crossing at this boundary carries.
    # ``assemble_priors`` needs it to turn a per-boundary object-incidence total into a fragment count.
    boundary_mass_per_crossing = substrate.boundary_unspliced.mass_per_crossing

    # The JUMPING population, exported verbatim. A sj boundary is pure RNA by construction, so there
    # is nothing to deconvolve: this is ``sj_count`` summed over the genome-strand columns.
    # ``assemble_priors`` does not read it — it is certified RNA in exactly the sense the spliced
    # crossings are withheld for — but the calibration's output should not be silent about the
    # population that at a donor boundary IS the gene's whole spliced output.
    count_rna_sj = np.asarray(substrate.sj.count, dtype=np.float64).sum(axis=1)

    # The two remaining INCIDENCE→FRAGMENT conversions, alongside `boundary_mass_per_crossing`. Each is
    # its own population's `mass / count`: applying one population's ratio to another is
    # `TRAPS: a-pooled-conversion-applied-per-component`. `CalibrationResult.library_rna_fragments`
    # derives the library count from them — a property, never a stored scalar, so an oracle arm that
    # swaps the mass arrays cannot inherit a count describing the arrays it replaced.
    boundary_spliced_mass_per_crossing = substrate.boundary_spliced.mass_per_crossing
    sj_mass_per_crossing = substrate.sj.mass_per_crossing

    result = CalibrationResult(
        mass_gdna_region=regions.gdna_mass,
        mass_rna_region=regions.rna_mass,
        mass_gdna_boundary=boundaries.gdna_mass,
        mass_rna_boundary=boundaries.rna_mass,
        mass_rna_spliced_boundary=mass_rna_spliced_boundary,
        boundary_mass_per_crossing=boundary_mass_per_crossing,
        count_rna_sj=count_rna_sj,
        boundary_spliced_mass_per_crossing=boundary_spliced_mass_per_crossing,
        sj_mass_per_crossing=sj_mass_per_crossing,
        gdna_region_eff_len=region_eff_gdna,
        gdna_boundary_eff_len=boundary_eff_gdna,
        rna_region_eff_len=region_eff_rna,
        rna_boundary_eff_len=boundary_eff_rna,
        # The simplex ψ solved, published per object rather than summed away. `mass_*` above is this
        # same answer with the two RNA strands added together.
        gdna_frac_region=regions.gdna_frac,
        rna_pos_frac_region=regions.rna_pos_frac,
        rna_neg_frac_region=regions.rna_neg_frac,
        gdna_frac_boundary=boundaries.gdna_frac,
        rna_pos_frac_boundary=boundaries.rna_pos_frac,
        rna_neg_frac_boundary=boundaries.rna_neg_frac,
        gdna_density_global=density_global,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_regions=int(substrate.n_regions),
        n_boundaries=int(substrate.n_boundaries),
        n_sj=int(substrate.n_sj),
        config=config,
    )
    # Diagnostic: the SJ sense fraction should agree with the StrandModel κ. A large gap flags a
    # strand-model / accumulator mismatch (κ stays the StrandModel posterior — this is QC only).
    # It is also "sense derived, never stored" in one line: the accumulator's columns are GENOME
    # strand, and which of them is *sense* is read off each sj's own annotated transcript strand.
    _flux = np.asarray(substrate.sj.count, dtype=np.float64)
    _is_pos = np.asarray(sj.strand) == np.int8(Strand.POS)
    spl_sense = float(np.where(_is_pos, _flux[:, 0], _flux[:, 1]).sum())
    spl_total = float(_flux.sum())
    sj_sense_frac = spl_sense / spl_total if spl_total > 0.0 else float("nan")
    logger.debug(
        "calibration: N=%d E=%d J=%d gdna_density_global=%.4g rna_sense_frac=%.3f "
        "gdna_strand_overdispersion=%.4g (%d seed regions, %d frags, %.1f effective%s%s) "
        "rna_strand_overdispersion=%.4g (%d sj, %d frags%s) "
        "[own-evidence od: rna=%.4g gdna=%.4g] "
        "[sj sense_frac=%.3f vs κ=%.3f]",
        result.n_regions,
        result.n_boundaries,
        result.n_sj,
        result.gdna_density_global,
        rna_sense_frac,
        gdna_strand_overdispersion,
        _gd_seed[0],
        _gd_seed[1],
        _gd_seed[3],
        ", FALLBACK" if _gd_seed[2] else ("" if _gd_seed[0] >= 0 else ", INJECTED"),
        (
            f", CLAMPED at the ceiling from a raw {_gd_seed[4]:.3f} - NOT a measurement"
            if (not _gd_seed[2]) and _gd_seed[4] > _MAX_OVERDISPERSION
            else ""
        ),
        rna_strand_overdispersion,
        _rna_seed[0],
        _rna_seed[1],
        ", FALLBACK" if _rna_seed[2] else ("" if _rna_seed[0] >= 0 else ", INJECTED"),
        _rna_seed[3],
        _gd_seed[4],
        sj_sense_frac,
        rna_sense_frac,
    )
    return result


__all__ = ["calibrate"]
