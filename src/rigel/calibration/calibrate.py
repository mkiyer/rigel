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
from .message_cache import MessageCache
from .sweep import chain_boundary_deconv, chain_region_deconv, solve_chain
from .density_model import count_observable_masks
from .derive import gdna_density_global
from .errors import CalibrationStrandError
from .density_deconv import (
    GdnaBackground,
    density_lambda_factor,
    fit_intron_background,
)
from .abundance_landscape import AbundanceLandscape, fit_abundance_landscape
from .blocks import SweepCapture
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
    are physically **directly observable** (no deconvolution / no solving): the RNA strand balance and the
    spliced sample size it was fit from, the strand Beta-Binomial overdispersions, the intergenic
    intron-factory background, and the pre-solve TOTAL-density landscape.

    A tiny (single-transcript) toy CANNOT fit these — so :func:`calibrate` accepts them pre-fit from a
    population scenario and injects them, letting the toy provide only the controlled per-region GEOMETRY. Every
    field is optional; ``None`` ⇒ fit that prior internally (the default, byte-identical). ``calibrate`` also
    stashes the fitted-or-injected bundle in ``_debug["calibration_priors"]`` so a population scenario's fitted
    priors can be extracted and re-injected into a toy (`scripts/debug/toy_inject.py`)."""

    rna_sense_frac: float | None = None
    n_rna_obs: float | None = None
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


class FactoryRows:
    """The gDNA intron factory's λ-factor rows, built per block on demand and never held for the chain.

    For each INTRON REGION slot, ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` over the σ(λ) solve grid
    (`density_lambda_factor`) — the factory deconvolves confident gDNA from introns against the
    intergenic background; ZERO on every other slot, a no-op there. BOUNDARY slots are zero structurally:
    the factor scores a CONTAINED count against a CONTAINED support, and a boundary's count is a
    crossing with a different divisor. gDNA is strand-symmetric, so the factor lives purely on ``λ`` and
    is consumed identically by every slot class.

    ``rows[sl]`` is the ``(len(sl), K)`` array for the slots ``sl`` — the sweep asks for each block's
    rows as it solves the block, so the chain-wide ``(n_slots, K)`` array (2–5 GB on the human chain,
    one per bracket) never exists; ``np.asarray(rows)`` materialises the whole chain and is for the
    diagnostic capture only.
    """

    def __init__(
        self, background: GdnaBackground, chain, substrate, region_arrays, region_eff_len, config
    ):
        self.background = background
        kind = np.asarray(chain.kind)
        idx = np.asarray(chain.obj_idx, dtype=np.int64)
        rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
        self.is_intron = (kind == REGION) & (rtype[np.clip(idx, 0, rtype.shape[0] - 1)] == 1)
        ridx = idx[self.is_intron]
        n = kind.shape[0]
        # GENOME-strand columns summed: gDNA is strand-symmetric, so the deconvolution is against a total
        self.count = np.zeros(n)
        self.count[self.is_intron] = np.asarray(
            substrate.region_contained.count, dtype=np.float64
        ).sum(axis=1)[ridx]
        self.eff = np.zeros(n)
        self.eff[self.is_intron] = np.asarray(region_eff_len, dtype=np.float64)[ridx]
        _, self.fg = _logodds_grid(
            lattice_points(config.sweep_logodds_window, config.sweep_logodds_step),
            float(config.sweep_logodds_window),
        )
        self.shape = (n, int(self.fg.shape[0]))

    def __getitem__(self, sl) -> np.ndarray:
        sel = self.is_intron[sl]
        out = np.zeros((sel.shape[0], self.shape[1]), dtype=np.float64)
        if sel.any():
            out[sel] = density_lambda_factor(
                self.background, self.count[sl][sel], self.eff[sl][sel], self.fg
            )
        return out

    def __array__(self, dtype=None, copy=None):
        return self[:] if dtype is None else self[:].astype(dtype)


#: Minimum training regions for a hyperprior fit — below this the population is not a population.
_MIN_TRAIN = 5


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
    domain_sel = np.asarray(eff_global, dtype=np.float64) > 1.0e-9
    # ⛔ The substrate guard measures the population the ANNOTATION admits, not the population left
    # after the composition cut: "is there enough of a population to fit a prior for" is not the
    # training set's question. Guarding the cut population instead refuses the refit on a gDNA-free
    # library whose only survivors are a handful of anchors, which then leaves its AMBIG regions at
    # their prior-free share. `fit_landscape` still refuses a training set under two.
    if int((sel | anchor).sum()) < _MIN_TRAIN:
        return None
    if belief.has_composition is not None:
        sel &= np.asarray(belief.has_composition, dtype=bool)
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


#: the QC seeds of an INJECTED strand value: no fit ran, so no seed regions, no fragments
_NO_GDNA_SEED = (-1, -1, False, float("nan"), float("nan"))
_NO_RNA_SEED = (-1, -1, False, float("nan"))


@dataclass(frozen=True, slots=True)
class _Strand:
    """The library's strand model as the solve reads it, fitted or injected: the RNA sense fraction
    ``κ``, ``N_rna`` — the spliced count κ was fit from, the sample the strand channel's protocol
    decision reads (`region_init.strand_discriminability`) — and the two Beta-Binomial overdispersions.
    The two seeds are QC for the log only — ``(n_seed_regions, n_seed_fragments, fallback, …)``,
    ``-1`` where the value was injected."""

    rna_sense_frac: float
    n_rna_obs: float
    gdna_strand_overdispersion: float
    rna_strand_overdispersion: float
    gdna_seed: tuple = _NO_GDNA_SEED
    rna_seed: tuple = _NO_RNA_SEED

    @property
    def model(self) -> tuple[float, float, float]:
        """``(κ, od_gdna, od_rna)`` — the triple the transfer policy's own strand claims need."""
        return (
            self.rna_sense_frac,
            self.gdna_strand_overdispersion,
            self.rna_strand_overdispersion,
        )


def _fit_strand(substrate, region_arrays, strand_models, inj) -> _Strand:
    """The strand model, in the order the fits lean on each other.

    ``κ`` first — the posterior-mean spliced sense fraction (`fit_strand_balance`) and the spliced count
    behind it; the strand channel's discriminability is ``(2κ−1)²`` where the protocol preserves strand
    and 0 where it does not (`region_init.strand_discriminability`, a decision on that same 2×2, so an
    unstranded library's count governs at any depth). A library with no spliced reads is not an RNA-seq
    library and raises. Then the RNA overdispersion
    (mean κ, from the per-sj strand table, certified pure RNA), which is the gDNA fit's fallback; then
    the gDNA overdispersion (mean ½ by dsDNA symmetry) by the AWAY-HALF moment over every genic count-
    and strand-observable object — unbiased under any RNA content of the seeds (`gdna_strand`'s
    lemma), so no seed is weighted and no class asserted pure; intergenic and AMBIG objects cannot be
    oriented and are out. The two are RECONCILED with no conjured target: the weaker-measured
    dispersion shrinks toward the better-measured one by their own null informations — ⛔ only when
    neither is injected, since an injected value is an arm's whole point."""
    if inj is not None and inj.rna_sense_frac is not None:
        rna_sense_frac = float(inj.rna_sense_frac)
        n_rna_obs = float(inj.n_rna_obs) if inj.n_rna_obs is not None else 0.0
    else:
        balance = fit_strand_balance(strand_models)
        if balance.fallback_used:
            raise CalibrationStrandError(
                "the library has zero spliced unique-mapper observations; this does not look like an "
                "RNA-seq library. A real RNA-seq library always carries spliced reads."
            )
        rna_sense_frac = float(balance.rna_sense_frac)
        n_rna_obs = float(balance.n_observations)

    gdna_seed, rna_seed = _NO_GDNA_SEED, _NO_RNA_SEED
    if inj is not None and inj.rna_strand_overdispersion is not None:
        rna_od = float(inj.rna_strand_overdispersion)
    else:
        rna_strand = fit_rna_strand_from_sj_table(
            strand_models.sj_table, rna_sense_frac=rna_sense_frac
        )
        rna_od = rna_strand.rna_strand_overdispersion
        rna_seed = (
            rna_strand.n_seed_regions,
            rna_strand.n_seed_fragments,
            rna_strand.fallback_used,
            rna_strand.raw_overdispersion,
        )
    if inj is not None and inj.gdna_strand_overdispersion is not None:
        gdna_od = float(inj.gdna_strand_overdispersion)
    else:
        # the seed selector is count-observability, straight off the signature
        region_obs, boundary_obs = count_observable_masks(
            np.asarray(region_arrays.signature), np.asarray(region_arrays.ref_id)
        )
        gdna_strand = fit_gdna_strand_from_substrate(
            substrate,
            region_arrays,
            region_count_observable=region_obs,
            boundary_count_observable=boundary_obs,
            rna_sense_frac=rna_sense_frac,
        )
        gdna_od = gdna_strand.gdna_strand_overdispersion
        gdna_seed = (
            gdna_strand.n_seed_regions,
            gdna_strand.n_seed_fragments,
            gdna_strand.fallback_used,
            gdna_strand.effective_seeds,
            gdna_strand.raw_overdispersion,
        )
    if inj is None or (
        inj.rna_strand_overdispersion is None and inj.gdna_strand_overdispersion is None
    ):
        rna_od, gdna_od = reconcile_overdispersions(
            rna_strand.raw_overdispersion,
            rna_strand.information,
            gdna_strand.raw_overdispersion,
            gdna_strand.information,
        )

    return _Strand(rna_sense_frac, n_rna_obs, gdna_od, rna_od, gdna_seed, rna_seed)


class _IntronFactory:
    """The gDNA INTRON FACTORY: the intergenic background — fitted here with ``include_introns=False``
    (an intron-inclusive pool is inflated by nascent RNA worst exactly where gDNA is scarce), or
    injected — and its λ-factor rows on a solve grid (:class:`FactoryRows`, built per block as the sweep
    asks). ``background`` is ``None`` when the factory is off, and ``rows`` is then ``None`` too, which
    leaves every sweep byte-identical to the pre-factory path.

    ⛔ The rows are evaluated ON the solve grid, so they are a function of ``(n_grid, L)`` and are
    REBUILT when the bracket widens, never regridded: `_regrid_global` maps between two grids of the
    same ``L``, and there is no map onto a wider domain the factor was never evaluated on."""

    def __init__(
        self, chain, substrate, region_arrays, region_eff_gdna, config, inj, background_pair
    ):
        self._site = (chain, substrate, region_arrays, region_eff_gdna, config)
        if inj is not None and inj.intron_background is not None:
            self.background = inj.intron_background
        else:
            self.background = fit_intron_background(
                substrate,
                region_arrays,
                region_eff_gdna,
                include_introns=False,
                counts_exposure=background_pair,
            )
        self._rows: dict = {}

    def rows(self, n_grid: int, window: float):
        """The λ-factor rows on the grid ``(n_grid, window)`` as a :class:`FactoryRows` — sliced per block
        by the sweep — or ``None`` when there is nothing to factor (the background uninformative, no
        intron regions), which leaves every sweep byte-identical to the pre-factory
        path."""
        if self.background is None or not self.background.informative:
            return None
        key = (int(n_grid), float(window))
        if key not in self._rows:
            chain, substrate, region_arrays, region_eff_gdna, config = self._site
            rows = FactoryRows(
                self.background,
                chain,
                substrate,
                region_arrays,
                region_eff_gdna,
                replace(config, sweep_logodds_window=float(window)),
            )
            self._rows[key] = rows if bool(rows.is_intron.any()) else None
        return self._rows[key]


def _wall_mask(payload, region_arrays, mature_walls, boundary_reach):
    return build_region_wall_mask(
        region_arrays,
        mature_walls,
        boundary_reach[0],
        boundary_reach[1],
        w_max=w_max_from_deposited_lengths(payload.deposited_lengths),
    )


def _background_pair(payload, substrate, region_arrays, config, mature_walls, boundary_reach):
    """The MEASURED-TOTAL ``(counts, exposure)`` pair for the pooled gDNA background estimators, built
    once and handed to each. ``None`` under the shipped default (``"contained"``), and every consumer
    then takes its own contained pair. ⛔ ``"measured_total"`` REFUSES to run without the wall inputs
    rather than falling back: a background rate that silently changed estimator because an argument
    was missing is worse than either estimator."""
    if config.background_abundance != "measured_total":
        return None
    if mature_walls is None or boundary_reach is None:
        raise ValueError(
            "CalibrationConfig.background_abundance = 'measured_total' needs the wall inputs: "
            "pass mature_walls=build_mature_wall_distances(index, region_arrays) and "
            "boundary_reach=build_contiguous_boundary_reach_arrays(index) (both are in "
            "scan_cache.index_derived_inputs). Refusing rather than falling back to the contained "
            "pair, which would change the background rate without saying so."
        )
    mask = _wall_mask(payload, region_arrays, mature_walls, boundary_reach)
    counts, exposure, _ = region_counts_and_exposure(substrate, region_arrays, mask)
    return (counts, exposure)


def _abundance_landscape(
    payload, substrate, region_arrays, config, inj, mature_walls, boundary_reach
):
    """THE ABUNDANCE LANDSCAPE — the pre-pass-0 TOTAL-density field + mode census, fitted at INIT from
    counts and lengths only (the wall-exact measured totals), so it is circular with nothing solved. A QC
    and injection surface: it is the sole source of the QC report's gDNA-density panel
    (`CalibrationDiagnostics.from_abundance_landscape`) and nothing in the solve reads it. Without the
    wall inputs (``mature_walls``, ``boundary_reach``) it is SKIPPED, LOUDLY, never raised for: many unit
    and toy callers have no wall arrays and want no panel, so the object stays ``None`` and the report
    omits the panel rather than carrying a quietly different estimate — deliberately not the policy of
    ``background_abundance``'s pair, which feeds ψ and REFUSES. (An on/off switch on this fit was a
    tunable nothing ever turned off; retired 2026-09-13.)"""
    if inj is not None and inj.abundance_landscape is not None:
        return inj.abundance_landscape
    if mature_walls is None or boundary_reach is None:
        logger.warning(
            "calibration: the wall inputs are missing (mature_walls / boundary_reach, both in "
            "scan_cache.index_derived_inputs) — skipping "
            "the total-density landscape, so the QC density panel will be omitted. Nothing in the "
            "solve reads it, so no solved number changes."
        )
        return None
    mask = _wall_mask(payload, region_arrays, mature_walls, boundary_reach)
    return fit_abundance_landscape(substrate, region_arrays, mask)


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
        s.chain,
        s.geometry,
        s.statics,
        rna_sense_frac=s.strand.rna_sense_frac,
        gdna_strand_overdispersion=s.strand.gdna_strand_overdispersion,
        rna_strand_overdispersion=s.strand.rna_strand_overdispersion,
        n_grid=lattice_points(s.config.sweep_logodds_window, s.config.sweep_logodds_step),
        logodds_window=s.config.sweep_logodds_window,
    )


def _sweep(s: _Solve, belief, prior, cache=None, capture=None):
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
        intron_prior=s.factory.rows(n_grid, window),
        policy=s.policy,
        block_slots=cfg.sweep_block_slots,
        message_cache=cache,
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
    E-step's start for the next. THE REFIT SWEEPS SHARE THEIR MESSAGE LAYER: the belief is reset before
    each and the messages never read the prior, so for one grid every input the layer reads is
    identical from refit to refit, and a refit pays only its two ψ solves (`MessageCache`,
    content-keyed; a widened bracket misses). Pass 0's grid is never reused.

    Returns ``(belief, belief_pass0, hyperprior)`` — the final belief, the prior-free one, and the last
    fitted landscape (``None`` if no refit ran). With ``_debug`` the last sweep fills
    ``_debug["capture"]`` (a :class:`~.blocks.SweepCapture`)."""
    capture = SweepCapture() if _debug is not None else None
    belief = _sweep(s, _init_belief(s), None, capture=capture)
    belief_pass0 = belief
    hyperprior: DensityLandscape | None = None
    cache = MessageCache()
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
        belief = _sweep(s, _init_belief(s), hyperprior, cache, capture=capture)
        logger.debug(
            "calibration: PHASE 2 gDNA-hyperprior refit %d/%d (%d training regions)",
            it + 1,
            s.config.calib_refit_iters,
            hyperprior.n_train,
        )
    if _debug is not None:
        _debug["capture"] = capture
    return belief, belief_pass0, hyperprior


def _result(
    substrate, sj, chain, belief, strand: _Strand, region_eff, boundary_eff, config
) -> CalibrationResult:
    """The solved chain projected onto the two payload axes and published as the
    :class:`CalibrationResult`, with the library-average gDNA density QC scalar.

    ``count_rna_spliced_boundary`` is the certified-RNA crossings per BOUNDARY — molecules that crossed
    contiguously having spliced elsewhere; `chain_boundary_deconv` adds the whole of it to ``rna_mass``
    (``rna = (1−g)·unspliced + spliced``), and `assemble_priors` withholds it from the RNA prior count,
    since a spliced fragment is guaranteed-RNA in the EM. There is no REGION twin, structurally: a
    region's contained population cannot hold a spliced molecule. ``count_rna_sj`` is the JUMPING
    population, exported verbatim — pure RNA by construction, nothing to deconvolve. The three
    ``mass_per_crossing`` are each their own population's incidence→fragment conversion: applying one
    population's ratio to another is `TRAPS: a-pooled-conversion-applied-per-component`."""
    regions = chain_region_deconv(chain, belief, substrate)
    boundaries = chain_boundary_deconv(chain, belief, substrate)
    region_eff_gdna, region_eff_rna = region_eff
    boundary_eff_gdna, boundary_eff_rna = boundary_eff
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
        rna_region_eff_len=region_eff_rna,
        rna_boundary_eff_len=boundary_eff_rna,
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
        rna_sense_frac=strand.rna_sense_frac,
        gdna_strand_overdispersion=strand.gdna_strand_overdispersion,
        rna_strand_overdispersion=strand.rna_strand_overdispersion,
        n_regions=int(substrate.n_regions),
        n_boundaries=int(substrate.n_boundaries),
        n_sj=int(substrate.n_sj),
        config=config,
    )


def _log_summary(result: CalibrationResult, strand: _Strand, substrate, sj) -> None:
    """The debug line: the sizes, the density scalar, the strand model with its seeds and clamps, and
    the sj sense fraction against κ — a large gap flags a strand-model / accumulator mismatch (κ stays
    the StrandModel posterior; this is QC only). Sense is derived, never stored: the accumulator's
    columns are GENOME strand, and which is sense is read off each sj's annotated transcript strand."""
    flux = np.asarray(substrate.sj.count, dtype=np.float64)
    is_pos = np.asarray(sj.strand) == np.int8(Strand.POS)
    spl_sense = float(np.where(is_pos, flux[:, 0], flux[:, 1]).sum())
    spl_total = float(flux.sum())
    sj_sense_frac = spl_sense / spl_total if spl_total > 0.0 else float("nan")
    gd, rn = strand.gdna_seed, strand.rna_seed
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
        strand.rna_sense_frac,
        strand.gdna_strand_overdispersion,
        gd[0],
        gd[1],
        gd[3],
        ", FALLBACK" if gd[2] else ("" if gd[0] >= 0 else ", INJECTED"),
        (
            f", CLAMPED at the ceiling from a raw {gd[4]:.3f} - NOT a measurement"
            if (not gd[2]) and gd[4] > _MAX_OVERDISPERSION
            else ""
        ),
        strand.rna_strand_overdispersion,
        rn[0],
        rn[1],
        ", FALLBACK" if rn[2] else ("" if rn[0] >= 0 else ", INJECTED"),
        rn[3],
        gd[4],
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
    diagnostics_out: dict | None = None,
    injected_priors: "InjectedCalibrationPriors | None" = None,
    boundary_flags: "np.ndarray | None" = None,
    mature_walls=None,
    boundary_reach=None,
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per object, then derive gdna_density_global — the
    stages in the module docstring's order, each a function above with one job.

    ``gdna_density_global`` may be ``0`` (a zero-gDNA library) and an object's deconvolved gDNA mass
    may be ``0`` (a pure-RNA object); both are valid, graceful outputs — not failures.

    ``sj`` is the splice graph's sj axis (:func:`~rigel.calibration.splice_graph.build_sj_geometry_arrays`),
    in the accumulator's own sj slot order; ``None`` means the graph has no sj boundaries, which is legal
    (a single-exon-only reference) and is NOT the same as "no sj flux". ``boundary_flags`` is the graph's
    per-contiguous-boundary structural bits, carried onto the chain as ``RegionStatics.boundary_flags``.
    ``mature_walls`` / ``boundary_reach`` are the two annotation-only WALL inputs the measured-total
    exposure and the abundance landscape need. ``injected_priors`` are population-scale priors a tiny toy
    cannot fit, injected in place of the internal fits (:class:`InjectedCalibrationPriors`).
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    inj = injected_priors
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
    background_pair = _background_pair(
        payload, substrate, region_arrays, config, mature_walls, boundary_reach
    )
    abundance_landscape = _abundance_landscape(
        payload, substrate, region_arrays, config, inj, mature_walls, boundary_reach
    )
    # the RNA twin: no consumer in the solve (the prior is a conserved FRAGMENT COUNT and divides by
    # nothing on the mass path), kept because it is byte-identically the opportunity the solver used
    region_eff_rna, boundary_eff_rna = _project_eff(chain, geometry.eff_rna, payload)

    strand = _fit_strand(substrate, region_arrays, strand_model, inj)
    factory = _IntronFactory(
        chain, substrate, region_arrays, region_eff_gdna, config, inj, background_pair
    )
    # ⛔ A TOTAL density over ONE component's opportunity model is not a composition estimate; the
    # per-slot gDNA support below is the basis the landscape prior is fit and read on, and the
    # total-density field this module does use is the abundance landscape above, which reaches the
    # report and never the solve.
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
    belief, belief_pass0, gdna_hyperprior = _solve(solve, _debug)
    logger.debug(
        "calibration: PHASE 1 prior-free initial solve (abundance landscape: %s)",
        "none"
        if abundance_landscape is None
        else f"{abundance_landscape.n_train} training regions, {len(abundance_landscape.modes)} modes",
    )

    result = _result(
        substrate,
        sj,
        chain,
        belief,
        strand,
        (region_eff_gdna, region_eff_rna),
        (boundary_eff_gdna, boundary_eff_rna),
        config,
    )

    if _debug is not None:  # inert diagnostic hook — the solved chain internals
        _debug.update(
            chain=chain,
            belief=belief,  # the FINAL belief (refit if calib_refit_iters>0, else the initial solve)
            belief_pass0=belief_pass0,  # the prior-free solve (the refit before/after frame)
            geometry=geometry,
            statics=statics,
            substrate=substrate,
            sj=sj,
            region_arrays=region_arrays,
            gdna_prior=abundance_landscape,  # the TOTAL-density landscape (QC / injection substrate)
            gdna_hyperprior=gdna_hyperprior,  # the DECONVOLVED-gDNA hyperprior (None if no refit)
            rna_sense_frac=strand.rna_sense_frac,
            region_eff_gdna=region_eff_gdna,
            boundary_eff_gdna=boundary_eff_gdna,
            # the fitted-or-injected population priors — extract from a population scenario, inject
            # into a toy
            calibration_priors=InjectedCalibrationPriors(
                rna_sense_frac=strand.rna_sense_frac,
                n_rna_obs=strand.n_rna_obs,
                gdna_strand_overdispersion=strand.gdna_strand_overdispersion,
                rna_strand_overdispersion=strand.rna_strand_overdispersion,
                intron_background=factory.background,
                abundance_landscape=abundance_landscape,
            ),
            abundance_landscape=abundance_landscape,
        )
    if diagnostics_out is not None and abundance_landscape is not None:
        # the QC density panel comes from the total-density landscape's census; the report omits it
        # when the landscape was not fit
        from .diagnostics import CalibrationDiagnostics

        diagnostics_out["calibration"] = CalibrationDiagnostics.from_abundance_landscape(
            abundance_landscape
        )
    _log_summary(result, strand, substrate, sj)
    return result


__all__ = ["calibrate"]
