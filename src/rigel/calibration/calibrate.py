"""calibrate() — the calibrator: deconvolve each object's UNSPLICED count into ``(RNA₊, RNA₋, gDNA)``.

This is LAYER 7's assembly stage. It owns the ORDER the pieces run in and the wiring between them;
every number it reports is produced by a module below it, and it computes no model of its own.

Calibration models RNA vs gDNA only (the per-locus EM separates transcripts downstream). An object's
composition has three sources: the STRAND LIKELIHOOD (the Beta-Binomial tilt of the per-strand
counts, the only intrinsic signal, entering as its overdispersed Fisher information rather than as a
raw count); the MESSAGES from its chain neighbours (:mod:`.messages`, each hop priced inside the
sweep from the two nodes' own counts); and the POPULATION gDNA PRIOR (the intergenic-only background
pool plus the phase-2 density landscape). The solver is the belief-propagation SWEEP over the
``N E N E … N`` chain (:mod:`rigel.calibration.sweep`)::

    substrate  (five populations on three axes)
      -> build chain + geometry + statics      (the geometry owns EVERY divisor)
      -> strand balance: rna_sense_frac (κ)
      -> strand overdispersion: 0 for both components (binomial, by policy)
      -> signature-binary init (G1/G2/G3)
      -> PASS 1 solve_chain (no fitted prior): ONE forward + ONE backward pass — each object integrates
           its strand likelihood, the intron factory and its neighbours' messages
      -> REFITS (``calib_refit_iters``): fit the gDNA-density landscape on the current belief, reset
           the belief, solve_chain again with the landscape added per object -> the per-object pie
      -> chain_region_deconv  -> per-REGION gDNA / RNA contained mass
      -> chain_boundary_deconv  -> per-BOUNDARY gDNA / RNA crossing mass, for the per-locus prior
      -> gdna_density_global (the library-average density QC scalar) and gdna_reference_density (the
         located enriched mode of the last refit's landscape — the ruler's reference; None ⇒ no contraction)

⛔ THE GEOMETRY IS BUILT FIRST, AND IT OWNS EVERY DIVISOR. ``build_region_geometry`` produces the
per-slot ``eff_gdna``/``eff_rna`` before anything reads a count, and everything downstream — the
density clue, the background pool, the intron factory, the result's two supports — reads THAT array.
Computing a second length model anywhere below is how a consumer comes to divide by a length the
solver never used.

A contiguous boundary is a 0-bp boundary with one count and one divisor, ``crossing_eff_length``;
there is no per-face machinery and no ½. A zero-gDNA library (``gdna_density_global == 0``, per-object
gDNA mass ``0``) is a valid, graceful output.

Known bias: the RNA half of an unspliced crossing takes
``UNBOUNDED_REACH`` rather than its transcript's real remaining length, which over-calls gDNA
genome-wide and worst in the last region before a polyA site. SJ boundaries DO take their
real exonic reach.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
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
from .sweep import chain_boundary_deconv, chain_region_deconv, solve_chain
from .derive import gdna_density_global
from .density_deconv import (
    GdnaBackground,
    fit_intron_background,
)
from .capture_efficiency import capture_efficiencies
from .effective_length import UNBOUNDED_REACH, conserved_cut_shares
from .blocks import SweepCapture
from .region_arrays import boundary_region_indices
from .region_chain import build_region_chain
from .result import CalibrationResult
from .landscape import DensityLandscape, fit_landscape, located_enriched_mode
from .signature import RegionType, coarse_type_array
from .strand_balance import fit_strand_balance
from .substrate import CalibrationSubstrate
from ..config import CONSTANTS
from ..types import Strand

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays
    from .splice_graph import SpliceJunctionGeometry

logger = logging.getLogger(__name__)


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


class FactoryRows:
    """The gDNA intron factory's λ-factor rows as their INPUTS: the background, the intron mask, each
    intron's count and opportunity — the solve's kernel builds the rows per block from these
    (`native/solve_kernel.cpp`, ``factory_row``), so no ``(n_slots, K)`` array ever exists.

    For each INTRON REGION slot the row is ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` over the σ(λ) solve
    grid — the factory deconvolves confident gDNA from introns against the intergenic background; NONE on
    every other slot, a no-op there. BOUNDARY slots carry none structurally: the factor scores a CONTAINED
    count against a CONTAINED support, and a boundary's count is a crossing with a different divisor. gDNA
    is strand-symmetric, so the factor lives purely on ``λ`` and is consumed identically by every slot
    class. ``kernel()`` is the tuple the sweep hands the kernel.
    """

    def __init__(self, background: GdnaBackground, chain, substrate, region_arrays, region_eff_len):
        self.background = background
        kind = np.asarray(chain.kind)
        idx = np.asarray(chain.obj_idx, dtype=np.int64)
        rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
        self.is_intron = (kind == REGION) & (
            rtype[np.clip(idx, 0, rtype.shape[0] - 1)] == RegionType.INTRON
        )
        ridx = idx[self.is_intron]
        n = kind.shape[0]
        # GENOME-strand columns summed: gDNA is strand-symmetric, so the deconvolution is against a total
        self.count = np.zeros(n)
        self.count[self.is_intron] = np.asarray(
            substrate.region_contained.count, dtype=np.float64
        ).sum(axis=1)[ridx]
        self.eff = np.zeros(n)
        self.eff[self.is_intron] = np.asarray(region_eff_len, dtype=np.float64)[ridx]

    def kernel(self) -> tuple:
        """The factory as the kernel takes it: its inputs — the intron mask, every slot's contained count and
        gDNA opportunity, the background's location, over-dispersion, size and whether it is informative."""
        bg = self.background
        return (
            "inputs",
            np.ascontiguousarray(self.is_intron, bool),
            np.ascontiguousarray(self.count, np.float64),
            np.ascontiguousarray(self.eff, np.float64),
            float(bg.log_mu_bg),
            float(bg.alpha),
            float(bg.size),
            bool(bg.informative),
        )


def _fit_gdna_hyperprior(
    chain, belief, statics, region_arrays, mass_global, eff_global, *, prev=None
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
    * location → ADMISSION AT THE FLOOR, precision above it → a continuous WEIGHT
      (`landscape._reliability`). A slot whose solve is wider than one nat² in ``log f_g``
      (``CONSTANTS.landscape.located_var``, the count rule's one-fragment wall through ``Var(log c) = 1/c``) has no
      location, whatever produced its solve: a strand term at a pure-RNA vertex, a factory row on an empty
      intron, a one-sided delivered row. Its median is where the reference measure sits under its bound,
      and training on it re-seeds the landscape at that resolution.
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
    constants = CONSTANTS.calibration
    expressed = (
        isr & (eff_global > constants.opportunity_floor) & (mass_global > constants.mass_floor)
    )
    # the zero-count structural anchor: an intergenic or intronic region that sequenced no unspliced mass
    anchor = (
        isr
        & (eff_global > constants.opportunity_floor)
        & (mass_global <= constants.mass_floor)
        & (rtype[ridx] != RegionType.EXON)
    )
    sel = expressed & ~(fp & fn)
    # ⛔ A SLOT WHOSE ONLY EVIDENCE IS A BOUND, OR WHICH HAS NONE, DOES NOT TRAIN THE PRIOR
    # (`RegionBelief.has_composition`). Its value is where the prior put it — outright, or inside the
    # half-line a level or ceiling admits — so re-fitting on it re-seeds the landscape's tail at the
    # density its own total implies. ⛔ Do not soften the cut by DISCOUNTING the delivered rows
    # instead (a one-sided row as a bound, the own-evidence variance as a weight): those readings lose
    # the deferred stratum badly, because those rows are the enriched mode's witness under capture.
    # The anchor trains regardless: it is a structural statement, not a solve.
    # THE GRID IS THE CONSUMERS' DOMAIN — EVERY slot the prior is read at, regions and boundaries
    # alike — so the composition cut above changes which kernels are summed and never the axis
    # (`fit_landscape`'s ``domain``). ⛔ Without that, a gDNA-free library whose exons are all blind
    # trains on its anchors alone, the grid collapses to the floor, and every exon reads a flat prior.
    domain_sel = np.asarray(eff_global, dtype=np.float64) > constants.opportunity_floor
    # ⛔ The substrate guard measures the population the ANNOTATION admits, not the population left
    # after the composition cut: "is there enough of a population to fit a prior for" is not the
    # training set's question. Guarding the cut population instead refuses the refit on a gDNA-free
    # library whose only survivors are a handful of anchors, which then leaves its AMBIG regions at
    # their prior-free share. `fit_landscape` still refuses a training set under two.
    if int((sel | anchor).sum()) < constants.min_training_regions:
        return None
    if belief.has_composition is not None:
        # a composition, AND a solve that locates it: both come from a sweep, so an initial belief
        # (no predicate, no solve) reads the annotation alone
        sel &= np.asarray(belief.has_composition, dtype=bool)
        sel &= np.asarray(belief.var_gdna, dtype=np.float64) <= CONSTANTS.landscape.located_var
    sel |= anchor
    mass = np.asarray(mass_global, dtype=np.float64)[sel]
    return fit_landscape(
        np.asarray(belief.f_g, dtype=np.float64)[sel] * mass,
        mass,
        np.asarray(eff_global, dtype=np.float64)[sel],
        np.asarray(belief.var_gdna, dtype=np.float64)[sel],
        anchor=anchor[sel],
        # the previous refit's landscape: the E-step on the location-free kernels (`landscape._estep_kernels`)
        prev=prev,
        domain=(
            np.asarray(mass_global, dtype=np.float64)[domain_sel],
            np.asarray(eff_global, dtype=np.float64)[domain_sel],
        ),
    )


def lattice_points(window: float, step: float) -> int:
    """The λ lattice's point count at bracket ``window``: ``λ ∈ [−L, L]`` at
    ``CalibrationConfig.sweep_logodds_step``, ``K = round(2L/step) + 1``.

    ⛔ The bracket and the RESOLUTION are two different knobs, and an arm that moves both at once is
    uninterpretable: widening ``L`` at a fixed point count coarsens the lattice, and the result then
    reverses with the bracket instead of saturating — and saturation is what distinguishes a truncated
    bracket from an improper ψ. So the STEP is what stays fixed when the landscape prior's derived demand
    widens the bracket (`_sweep`), and the point count follows. The rounding is the only slack: the step
    is preserved to within half a point, the best a discrete lattice can do.
    """
    return int(round(2.0 * float(window) / float(step))) + 1


@dataclass(frozen=True, slots=True)
class _Strand:
    """The library's strand model as the solve reads it: the RNA sense fraction ``κ``, ``N_rna`` — the spliced count
    κ was fit from, the sample the strand channel's protocol decision reads (`region_init.strand_discriminability`)
    — and the two Beta-Binomial overdispersions, both 0 (binomial) by policy."""

    rna_sense_frac: float
    n_rna_obs: float
    gdna_strand_overdispersion: float
    rna_strand_overdispersion: float

    @property
    def model(self) -> tuple[float, float, float]:
        """``(κ, od_gdna, od_rna)`` — the triple the transfer policy's own strand claims need."""
        return (
            self.rna_sense_frac,
            self.gdna_strand_overdispersion,
            self.rna_strand_overdispersion,
        )


def _fit_strand(strand_models) -> _Strand:
    """The strand model: ``κ``, the posterior-mean spliced sense fraction (`fit_strand_balance`), and the
    spliced count behind it; the strand channel's discriminability is ``(2κ−1)²`` where the protocol preserves
    strand and 0 where it does not (`region_init.strand_discriminability`, a decision on that same 2×2).

    Both overdispersions are 0 — the strand channel's width is binomial — by policy: no population measures an od
    cleanly on real data (antisense RNA reaches the gDNA seeds, splice artifacts the junctions, and at n ≤ 3 a
    seed's spread is composition as much as dispersion), and an od set too high mutes the channel. A library with
    no spliced read is unstranded: κ = ½ and a dead channel, never an error. Gates:
    `tests/calibration/test_calibrate.py`'s ``test_the_strand_overdispersion_is_binomial_by_policy`` and
    ``test_an_empty_spliced_census_calibrates_as_unstranded``."""
    balance = fit_strand_balance(strand_models)
    return _Strand(float(balance.rna_sense_frac), float(balance.n_observations), 0.0, 0.0)


class _IntronFactory:
    """The gDNA INTRON FACTORY: the intergenic background fitted here, and its λ-factor
    rows as a :class:`FactoryRows` (the inputs the kernel builds each block's rows from). ``rows`` is
    ``None`` when there is nothing to factor: the background uninformative, or no intron region."""

    def __init__(self, chain, substrate, region_arrays, region_eff_gdna):
        self.background = fit_intron_background(substrate, region_arrays, region_eff_gdna)
        self.rows = None
        if self.background.informative:
            rows = FactoryRows(self.background, chain, substrate, region_arrays, region_eff_gdna)
            if rows.is_intron.any():
                self.rows = rows


def _policy(config, strand: _Strand):
    """The message-composition policy the config names. ⛔ THE NAME MUST SELECT THE POLICY — an arm that
    silently runs a different policy than it names is a benchmark that cannot be trusted, so an unknown
    name raises. The transfer policy's own strand claims read the library's strand model; an intron's
    own claim is the factory's row for it, which the sweep hands over on the context."""
    if config.message_policy == "silent":
        return SilentPolicy()
    if config.message_policy == "transfer":
        return TransferPolicy(strand=strand.model)
    raise ValueError(
        f"unknown message_policy {config.message_policy!r} — expected 'silent' or 'transfer'"
    )


@dataclass(frozen=True, slots=True)
class _Solve:
    """Everything a sweep reads, fixed for the whole calibration: the chain, its statics and geometry,
    the region arrays, the strand model, the intron factory, the policy, the config, and the per-slot
    gDNA support ``(mass_global, eff_global)`` the landscape prior is fit and read on."""

    chain: object
    statics: object
    geometry: object
    region_arrays: object
    strand: _Strand
    factory: _IntronFactory
    policy: object
    config: object
    mass_global: np.ndarray
    eff_global: np.ndarray


def _init_belief(s: _Solve):
    """The signature-binary G1/G2/G3 belief on the chain, before any sweep."""
    return init_beliefs(
        s.geometry,
        s.statics,
        rna_sense_frac=s.strand.rna_sense_frac,
        gdna_strand_overdispersion=s.strand.gdna_strand_overdispersion,
        rna_strand_overdispersion=s.strand.rna_strand_overdispersion,
        n_grid=lattice_points(s.config.sweep_logodds_window, s.config.sweep_logodds_step),
        logodds_window=s.config.sweep_logodds_window,
        n_threads=int(s.config.n_threads),
    )


def _sweep(s: _Solve, belief, prior, capture=None):
    """One sweep of the chain from ``belief``, with the composition prior ``prior`` (``None``: the
    prior-free pass) — the λ bracket first, then `solve_chain`.

    THE λ BRACKET IS ``max(the reference's floor, the fitted prior's own demand)``. ψ evaluates the
    landscape at ``log ρ = log f + log M − log E`` and can only offer ``f ∈ [σ(−L), σ(L)]``, so a bracket
    narrower than the prior's support leaves ψ no coordinate for what the prior says — and the answer
    then depends on ``L``. The demand is DERIVED (`required_logodds_window`), never chosen; the
    prior-free pass has nothing to widen for and keeps the floor. The lattice STEP is held fixed, so the
    point count scales with the bracket (`lattice_points`): widening ``L`` at a fixed count would coarsen
    the lattice and confound two knobs. The TILT axis does not scale — θ is a share with no bracket
    problem — which keeps the AMBIG cube linear in the bracket.

    ⛔ ψ has NO reference location. A located reference is a prior assertion at fixed strength and
    becomes the whole answer wherever the strand channel is dead; background information enters as the
    factory's λ-factor, a likelihood whose precision scales with counts. Every message's price is
    self-contained in the sweep (each hop charges the two nodes' counting and the pair's disagreement)
    — there is nothing to fit here."""
    cfg = s.config
    window = float(cfg.sweep_logodds_window)
    if prior is not None:
        required = prior.required_logodds_window(s.mass_global, s.eff_global)
        if required > window:
            window = required
            logger.debug(
                "calibration: λ bracket %.4f (the landscape's support), %d points at step %g",
                window,
                lattice_points(window, cfg.sweep_logodds_step),
                cfg.sweep_logodds_step,
            )
    n_grid = lattice_points(window, cfg.sweep_logodds_step)
    return solve_chain(
        s.chain,
        s.statics,
        s.geometry,
        belief,
        s.region_arrays,
        rna_sense_frac=s.strand.rna_sense_frac,
        gdna_strand_overdispersion=s.strand.gdna_strand_overdispersion,
        rna_strand_overdispersion=s.strand.rna_strand_overdispersion,
        n_rna_obs=s.strand.n_rna_obs,
        n_grid=n_grid,
        logodds_window=window,
        gdna_prior=prior,
        intron_prior=s.factory.rows,
        policy=s.policy,
        block_slots=cfg.sweep_block_slots,
        n_threads=int(cfg.n_threads),
        _capture=capture,
    )


def _solve(s: _Solve, _debug):
    """The two phases. PHASE 1, the INITIAL solve, carries no fitted composition prior: the inert
    Beta(½,½) reference alone plus the strand likelihood and the messages — single-strand regions
    self-solve from strand, AMBIG regions on unstranded data are grounded only by the messages, and
    their two-root ambiguity is what phase 2 resolves. PHASE 2, the DECONVOLVED-gDNA hyperprior REFIT:
    fit `landscape.DensityLandscape` on the previous solve's deconvolved gDNA (`_fit_gdna_hyperprior`
    selects the training substrate), reset the belief in FULL — nothing from pass 0 survives but the
    fitted landscape, so an over-confident region cannot refuse to budge when the prior lands — and
    re-solve with it as ψ's composition arm; ``calib_refit_iters`` times, each refit's landscape the
    E-step's start for the next. Every sweep runs the whole message layer: the messages never read the
    prior, so on one grid a refit's messages equal the previous refit's.

    Returns ``(belief, hyperprior)`` — the final belief and the last fitted landscape (``None`` if no
    refit ran). With ``_debug`` the last sweep fills ``_debug["capture"]`` (a
    :class:`~.blocks.SweepCapture`)."""
    capture = SweepCapture() if _debug is not None else None
    belief = _sweep(s, _init_belief(s), None, capture=capture)
    hyperprior: DensityLandscape | None = None
    for it in range(int(s.config.calib_refit_iters)):
        hyperprior = _fit_gdna_hyperprior(
            s.chain,
            belief,
            s.statics,
            s.region_arrays,
            s.mass_global,
            s.eff_global,
            prev=hyperprior,
        )
        if hyperprior is None:
            break
        capture = SweepCapture() if _debug is not None else None
        belief = _sweep(s, _init_belief(s), hyperprior, capture=capture)
        logger.debug(
            "calibration: PHASE 2 gDNA-hyperprior refit %d/%d (%d training regions)",
            it + 1,
            s.config.calib_refit_iters,
            hyperprior.n_train,
        )
    if _debug is not None:
        _debug["capture"] = capture
    return belief, hyperprior


def _result(
    substrate,
    regions,
    boundaries,
    strand: _Strand,
    region_eff_gdna,
    boundary_eff_gdna,
    gdna_reference_density: float | None,
    gdna_reference_members: int,
    efficiency: np.ndarray,
    efficiency_boundary: np.ndarray,
    boundary_conserved_gdna: np.ndarray,
) -> CalibrationResult:
    """The solved chain projected onto the two payload axes and published as the
    :class:`CalibrationResult`, with the library-average gDNA density QC scalar.

    ``count_rna_spliced_boundary`` is the certified-RNA crossings per BOUNDARY — molecules that crossed
    contiguously having spliced elsewhere; `chain_boundary_deconv` adds the whole of it to ``rna_mass``
    (``rna = (1−g)·unspliced + spliced``); the EM's prior reads only the gDNA mass
    (`priors.assemble_priors`). There is no REGION twin, structurally: a
    region's contained population cannot hold a spliced molecule. ``count_rna_sj`` is the JUMPING
    population, exported verbatim — pure RNA by construction, nothing to deconvolve. The three
    ``mass_per_crossing`` are each their own population's incidence→fragment conversion, never applied
    to another population."""
    return CalibrationResult(
        count_gdna_region=regions.gdna_mass,
        count_rna_region=regions.rna_mass,
        count_gdna_boundary=boundaries.gdna_mass,
        count_rna_boundary=boundaries.rna_mass,
        count_rna_spliced_boundary=np.asarray(
            substrate.boundary_spliced.count, dtype=np.float64
        ).sum(axis=1),
        boundary_mass_per_crossing=substrate.boundary_unspliced.mass_per_crossing,
        count_rna_sj=np.asarray(substrate.sj.count, dtype=np.float64).sum(axis=1),
        boundary_spliced_mass_per_crossing=substrate.boundary_spliced.mass_per_crossing,
        sj_mass_per_crossing=substrate.sj.mass_per_crossing,
        gdna_region_eff_len=region_eff_gdna,
        gdna_boundary_eff_len=boundary_eff_gdna,
        gdna_boundary_conserved_len=boundary_conserved_gdna,
        # the simplex ψ solved, published per object; the masses above are the same answer with the
        # two RNA strands added together
        gdna_frac_region=regions.gdna_frac,
        rna_pos_frac_region=regions.rna_pos_frac,
        rna_neg_frac_region=regions.rna_neg_frac,
        gdna_frac_boundary=boundaries.gdna_frac,
        rna_pos_frac_boundary=boundaries.rna_pos_frac,
        rna_neg_frac_boundary=boundaries.rna_neg_frac,
        gdna_density_global=gdna_density_global(
            regions, boundaries, region_eff_gdna, boundary_eff_gdna
        ),
        gdna_reference_density=gdna_reference_density,
        gdna_reference_members=gdna_reference_members,
        gdna_capture_efficiency_region=efficiency,
        gdna_capture_efficiency_boundary=efficiency_boundary,
        rna_sense_frac=strand.rna_sense_frac,
        gdna_strand_overdispersion=strand.gdna_strand_overdispersion,
        rna_strand_overdispersion=strand.rna_strand_overdispersion,
        n_regions=int(substrate.n_regions),
        n_boundaries=int(substrate.n_boundaries),
        n_sj=int(substrate.n_sj),
    )


def _gdna_boundary_conserved_len(region_arrays, gdna_fl_pmf: np.ndarray) -> np.ndarray:
    """gDNA's conserved share at every boundary (`effective_length.conserved_cut_shares`): the boundary's
    two flanking regions are the pieces beside the cut, and gDNA's template does not end."""
    lo, hi = boundary_region_indices(np.asarray(region_arrays.ref_id))
    length = np.asarray(region_arrays.end, dtype=np.float64) - np.asarray(
        region_arrays.start, dtype=np.float64
    )
    left, right = conserved_cut_shares(
        gdna_fl_pmf, length[lo], length[hi], UNBOUNDED_REACH, UNBOUNDED_REACH
    )
    return left + right


def _log_summary(result: CalibrationResult, strand: _Strand, substrate, sj) -> None:
    """The debug line: the sizes, the density scalar, the strand model, and
    the sj sense fraction against κ — a large gap flags a strand-model / accumulator mismatch (κ stays
    the StrandModel posterior; this is QC only). Sense is derived, never stored: the accumulator's
    columns are GENOME strand, and which is sense is read off each sj's annotated transcript strand."""
    flux = np.asarray(substrate.sj.count, dtype=np.float64)
    is_pos = np.asarray(sj.strand) == np.int8(Strand.POS)
    spl_sense = float(np.where(is_pos, flux[:, 0], flux[:, 1]).sum())
    spl_total = float(flux.sum())
    sj_sense_frac = spl_sense / spl_total if spl_total > 0.0 else float("nan")
    if result.gdna_reference_density is None:
        logger.info(
            "calibration: no located enriched gDNA mode — effective lengths are not corrected for "
            "capture (the reference needs about √n located probed pieces at one gDNA fragment or more)"
        )
    else:
        logger.info(
            "calibration: capture reference %.3e gDNA fragments/bp, the located enriched mode of %d "
            "kernels at one fragment or more",
            result.gdna_reference_density,
            result.gdna_reference_members,
        )
    logger.debug(
        "calibration: N=%d E=%d J=%d gdna_density_global=%.4g rna_sense_frac=%.3f "
        "strand overdispersion 0 (binomial, by policy) [sj sense_frac=%.3f vs κ=%.3f]",
        result.n_regions,
        result.n_boundaries,
        result.n_sj,
        result.gdna_density_global,
        strand.rna_sense_frac,
        sj_sense_frac,
        strand.rna_sense_frac,
    )


def calibrate(
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    gdna_fl_pmf: "np.ndarray",
    rna_fl_pmf: "np.ndarray",
    config: "CalibrationConfig",
    sj: "SpliceJunctionGeometry | None" = None,
    _debug: dict | None = None,
    boundary_flags: "np.ndarray | None" = None,
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per object, then derive gdna_density_global — the
    stages in the module docstring's order, each a function above with one job.

    ``gdna_density_global`` may be ``0`` (a zero-gDNA library) and an object's deconvolved gDNA mass
    may be ``0`` (a pure-RNA object); both are valid, graceful outputs — not failures.

    ``sj`` is the splice graph's sj axis (:func:`~rigel.calibration.splice_graph.build_sj_geometry_arrays`),
    in the accumulator's own sj slot order; ``None`` means the graph has no sj boundaries, which is legal
    (a single-exon-only reference) and is NOT the same as "no sj flux". ``boundary_flags`` is the graph's
    per-contiguous-boundary structural bits, carried onto the chain as ``RegionStatics.boundary_flags``.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
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
    # `eff_gdna` / `eff_rna` are the CONTAINED placements at a REGION slot and the CROSSING placements at
    # a BOUNDARY slot, one rule, one array. Nothing below computes a second length model; the result's
    # supports are PROJECTED off the geometry (`_project_eff`), so the number `priors` divides by is
    # byte-identically the one the solver divided by.
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geometry = build_region_geometry(chain, substrate, region_arrays, sj, gdna_fl_pmf, rna_fl_pmf)
    statics = build_region_statics(chain, region_arrays, boundary_flags)
    region_eff_gdna, boundary_eff_gdna = _project_eff(chain, geometry.eff_gdna, payload)

    strand = _fit_strand(strand_model)
    factory = _IntronFactory(chain, substrate, region_arrays, region_eff_gdna)
    # ⛔ A TOTAL density over ONE component's opportunity model is not a composition estimate; the
    # per-slot gDNA support below is the basis the landscape prior is fit and read on.
    mass_global, eff_global = region_gdna_geometry(geometry)
    solve = _Solve(
        chain,
        statics,
        geometry,
        region_arrays,
        strand,
        factory,
        _policy(config, strand),
        config,
        mass_global,
        eff_global,
    )
    belief, gdna_hyperprior = _solve(solve, _debug)
    # THE RULER'S REFERENCE: the fully-captured gDNA level is the located enriched mode of the fitted
    # landscape, or nothing — capture-OFF and gDNA-free libraries carry no enriched mode and contract
    # nothing. One definition, read by `capture_eff_length` and `priors`.
    enriched = located_enriched_mode(gdna_hyperprior) if gdna_hyperprior is not None else None
    gdna_reference_density = float(np.exp(enriched.mode.log_rho)) if enriched is not None else None
    gdna_reference_members = enriched.n_members if enriched is not None else 0
    # THE CAPTURE EFFICIENCIES: every object's clipped gDNA density against the reference, the posterior
    # mean under the landscape from its own count — a region's contained, a boundary's crossing
    # (`capture_efficiency`); published on the result for the ruler and the locus prior. No reference ⇒
    # every efficiency is exactly 1.
    regions = chain_region_deconv(chain, belief, substrate)
    boundaries = chain_boundary_deconv(chain, belief, substrate)
    if gdna_reference_density is None:
        efficiency = np.ones(int(substrate.n_regions))
        efficiency_boundary = np.ones(int(substrate.n_boundaries))
    else:
        efficiency, efficiency_boundary = capture_efficiencies(
            gdna_hyperprior,
            gdna_reference_density,
            regions.gdna_mass,
            region_eff_gdna,
            boundaries.gdna_mass,
            boundary_eff_gdna,
        )
    result = _result(
        substrate,
        regions,
        boundaries,
        strand,
        region_eff_gdna,
        boundary_eff_gdna,
        gdna_reference_density,
        gdna_reference_members,
        efficiency,
        efficiency_boundary,
        _gdna_boundary_conserved_len(region_arrays, gdna_fl_pmf),
    )

    if _debug is not None:  # inert diagnostic hook — the solved chain internals
        _debug.update(
            chain=chain,
            belief=belief,  # the FINAL belief (refit if calib_refit_iters>0, else the initial solve)
            geometry=geometry,
            statics=statics,
            substrate=substrate,
            sj=sj,
            region_arrays=region_arrays,
            gdna_hyperprior=gdna_hyperprior,  # the DECONVOLVED-gDNA hyperprior (None if no refit)
            rna_sense_frac=strand.rna_sense_frac,
        )
    _log_summary(result, strand, substrate, sj)
    return result


__all__ = ["calibrate"]
