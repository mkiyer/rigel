"""calibrate() — the calibrator: deconvolve each object's UNSPLICED count into ``(RNA₊, RNA₋, gDNA)``.

Calibration models **RNA vs gDNA only** (the per-locus EM separates nascent from mature downstream). Per the
count-zero-information architecture a region's composition is set by three sources — the **strand likelihood**
(the Beta-Binomial tilt of the per-strand counts, the only INTRINSIC signal; the count enters ONLY as its
overdispersed Fisher information), the **cross-region imputation** (the neighbour density messages, at the
source's own honest belief precision (strand + count) ``pr = n_src/(n_src·vb_src + 1)``), and the
**population gDNA prior** (the conservative intergenic+intron floor + the Phase-2 density KDE). The solver
is the belief-propagation SWEEP over the ``N E N E … N`` chain (:mod:`rigel.calibration.sweep`)::

    substrate  (five populations on three axes)
      -> build chain + geometry + statics      (the geometry owns EVERY divisor)
      -> strand balance: rna_sense_frac (κ)
      -> region_gdna_density (RAW) -> fit gDNA / RNA strand Beta-Binomial overdispersions (seed)
      -> signature-binary init (G1/G2/G3)
      -> PASS 1 solve_chain (no KDE prior): ONE forward + ONE backward pass — each object integrates its
           strand likelihood + the conservative gDNA floor + the belief-free density messages
      -> train the Phase-2 gDNA-density KDE on the pass-1 belief
      -> PASS 2 solve_chain (the KDE prior added per object) -> the converged per-object pie
      -> chain_region_deconv  -> per-REGION gDNA / RNA contained mass
      -> chain_boundary_deconv  -> per-BOUNDARY gDNA / RNA crossing mass, for the per-locus prior
      -> gdna_density_global (the library-average density QC scalar)

⭐ **THE GEOMETRY IS BUILT FIRST, AND IT OWNS EVERY DIVISOR.** The predecessor computed three effective
lengths up front (``region_eff_length``, ``boundary_side_eff_length``, ``boundary_eff_length``) and
handed them around; the density model took a fourth, a scalar ``fl_mean``. That is four chances for a
consumer to divide by a length model the solver never used. ``build_region_geometry`` now produces the
per-slot ``eff_gdna``/``eff_rna`` before anything reads a count, and everything downstream — the density
clue, the background pool, the intron factory, the result's two supports — reads THAT array


⛔ **The per-face boundary machinery is gone.** ``boundary_substrate``, the ``left``/``right`` views and
``boundary_side_eff_length``'s ``E[min(ℓ,L)]/2`` had no successor: a contiguous boundary is a 0-bp boundary with
one count and one divisor, ``crossing_eff_length``. A zero-gDNA library
(``gdna_density_global == 0``, per-object gDNA mass ``0``) remains a valid, graceful output.

⚠ **TRAPS: prove-the-substrate IS OFF AT CONTIGUOUS BOUNDARIES, AND THE FIRST BASELINE CARRIES ITS BIAS.** The RNA half of an
unspliced crossing takes ``UNBOUNDED_REACH`` rather than its transcript's real remaining length
(owner-ruled). Cost, already measured: an **11.0 %** genome-wide gDNA
over-call and **+0.36** in the last region before a polyA site. It is
sequenced after this step so that turning it on can be A/B'd against the baseline this step produces.
SpliceJunction boundaries DO take their real exonic reach — they are a new population with no predecessor divisor.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

import numpy as np

from .background_reference import BackgroundReference, measure_background
from .messages.currency import CurrencyPolicy
from .messages.relay import RelayPolicy
from .messages.silent import SilentPolicy
from .region_chain import BOUNDARY, REGION
from .region_geometry import (
    build_region_geometry,
    build_region_statics,
    init_beliefs,
    region_gdna_geometry,
)
from .sweep import chain_boundary_deconv, chain_region_deconv, solve_chain
from .density_model import region_gdna_density
from .derive import gdna_density_global
from .errors import CalibrationStrandError
from .npmle import DensityNPMLE
from .density_deconv import GdnaBackground, density_lambda_factor, fit_intron_background
from .gdna_strand import (
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
    Beta-Binomial overdispersions, the strand-Fisher noise-floor sample sizes, the enrichment-density NPMLE (the
    σ²_transfer landscape), the intergenic intron-factory background, and the aggregate ρ_bg background.

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
    enrichment_prior: DensityNPMLE | None = None
    intron_background: GdnaBackground | None = None
    background: BackgroundReference | None = None


def _empty_sj_geometry() -> "SpliceJunctionGeometry":
    """A sj axis with no rows — a graph whose references are all single-exon.

    ⚠ Legal, and **not** the same as "no sj flux": a payload scanned against such a graph has
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

    ⭐ **A projection, not a recomputation.** ``build_region_geometry`` already applied
    ``contained_eff_length`` at REGION slots and ``crossing_eff_length`` at BOUNDARY slots; calling those
    again here would put two implementations of one quantity in the tree, which is precisely how the
    prose and the code came to disagree about a ½ for months ( and 27).
    Whatever the solver divided by is what the result reports.

    ⚠ It takes the ARRAY, not the geometry, so ONE function serves both populations rather than a
    second copy that drifts. ⛔ It was written because ``assemble_priors`` needed the RNA divisor as
    well as the gDNA one; the prior no longer divides by anything, so only the gDNA projection has a
    live consumer today. The RNA one is retained deliberately — `result.py` carries the reasoning.
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

    ⚠ **BOUNDARY slots are zero, structurally.** The factor scores a CONTAINED count against a CONTAINED
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
    # ⚠ GENOME-strand columns summed: gDNA is strand-symmetric, so the deconvolution is against a total rate.
    count = np.asarray(substrate.region_contained.count, dtype=np.float64).sum(axis=1)[ridx]
    eff_g = np.asarray(region_eff_len, dtype=np.float64)[ridx]
    prior[is_intron] = density_lambda_factor(bg, count, eff_g, fg)
    return prior


#: Minimum training regions for a hyperprior fit — below this the population is not a population.
_MIN_TRAIN = 5


def _fit_gdna_hyperprior(
    chain, belief, statics, region_arrays, mass_global, eff_global, *, strength
):
    """Select the training substrate from the chain and fit the :class:`DensityLandscape` on the initial solve's
    deconvolved gDNA — the composition (gDNA) arm of ψ for the Phase-2 refit. ``None`` if it cannot be fit.

    **This affects only the PRIOR fit — never the solve's gDNA messages** (the G1/TSS/TES boundary emissions
    in ``sweep.solve_chain`` are a separate mechanism and are untouched).

    ⭐⭐ **The substrate is the whole of this function's job; the estimator in :mod:`.landscape` is
    component-agnostic and knows nothing about gDNA.** Two facts that used to live in that module's
    docstring are stated here instead, because they are decisions about WHICH OBJECTS TRAIN IT rather
    than about how it fits:

    * **BOUNDARIES ARE EXCLUDED** (owner, 2026-07-27). They are ~as numerous as regions but only 5.1 %
      of them are truly enriched against the regions' 12.1 %, so including them nearly halves the
      enriched component's census, and their two-flank mixture fills the valley between the two true
      modes — 74 % of all valley mass.
    * **THE ZERO-COUNT ANCHOR IS CRITICALLY LOAD-BEARING**, and what it is evidence *of* is specific to
      gDNA: a live region that sequenced no unspliced mass has gDNA density ``0`` for every ``f_g``.
      Dropping it costs **+0.26 / +0.61 EMD, and +1.04 on zero-gDNA libraries**.

    Four axes decide membership, and conflating them is the mistake this replaces
    (production plan §2.2):

    * **circularity → structural exclusion.** AMBIG regions are out. They are the two-root ambiguity the prior
      exists to resolve, so training on them would have it predict itself. ⚠ Note the *empirical* case
      alongside this is stratum-dependent — it holds over all conditions but reverses on the ones where an
      enriched mode exists — so the argument carrying the decision is the circularity, not the EMD.
    * **identifiability → structural INCLUSION.** A live region with no unspliced mass has density ``0`` for
      every ``f_g``: "gDNA is absent here" is the strongest depletion evidence there is, and it anchors the
      depleted mode. Dropping it costs +0.26 / +0.61 EMD (+1.04 on zero-gDNA libraries).
    * **precision → a continuous weight**, never admission (`landscape._reliability`).
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
    sel = (expressed & ((fp ^ fn) | (~fp & ~fn))) | anchor
    if int(sel.sum()) < _MIN_TRAIN:
        return None
    mass = np.asarray(mass_global, dtype=np.float64)[sel]
    return fit_landscape(
        np.asarray(belief.f_g, dtype=np.float64)[sel] * mass,
        mass,
        np.asarray(eff_global, dtype=np.float64)[sel],
        np.asarray(belief.var_gdna, dtype=np.float64)[sel],
        anchor=anchor[sel],
        strength=strength,
    )


def _scaled_grid(n: int, window: float, required: float) -> int:
    """Grid points at the widened bracket, holding the lattice spacing ``dlam = 2L/(n−1)`` FIXED.

    ⛔ **The bracket and the RESOLUTION are two different knobs, and an arm that moves both at once is
    uninterpretable.** Measured on the zero-gDNA control: widening ``L`` at a FIXED ``n`` coarsens
    ``dlam`` and reads better at 20 but **REVERSES at 40** (`g50 ss_0.50 off` 1.0383×), while holding
    ``dlam`` fixed SATURATES (0.9843× → 0.9842×). Saturation is what distinguishes a truncated bracket
    from an improper ψ, so the spacing must not move. ``n`` is therefore linear in ``L``, exactly.

    ⚠ The rounding is the only slack: ``dlam`` is preserved to within half a grid point, which is the
    best a discrete lattice can do and is far below the effect being measured.
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

    # ⭐ THE CHAIN AND ITS GEOMETRY COME FIRST, and the geometry owns every divisor from here down:
    # `eff_gdna`/`eff_rna` are the CONTAINED placements at a REGION slot and the CROSSING placements at an
    # BOUNDARY slot, one rule, one array. Nothing below computes a second length model.
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geometry = build_region_geometry(
        chain,
        substrate,
        region_arrays,
        sj,
        gdna_fl_pmf,
        rna_fl_pmf,
        boundary_rna_reach,
        region_abundance_bank=config.region_abundance_bank,
    )
    statics = build_region_statics(chain, region_arrays, boundary_flags)

    # The result's two gDNA supports, PROJECTED off the geometry rather than recomputed — so the number
    # `priors` divides by is byte-identically the one the solver divided by. Two implementations
    # of one quantity is how they come to disagree.
    region_eff_gdna, boundary_eff_gdna = _project_eff(chain, geometry.eff_gdna, payload)
    # ⭐ And the RNA twin, on the same two axes — the RNA population's own opportunity per object.
    # ⛔ It USED to be what turned ``mass_rna_*`` into a density in ``assemble_priors``. That consumer
    # is gone: the prior became a conserved FRAGMENT COUNT and divides by nothing on the mass path
    # (`d045d820`, `5591cc01`). Kept because it is the RNA divisor any density-based prior needs and
    # because it is byte-identically the opportunity the solver used — see `result.py`.
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

    # Count clue on RAW counts (the count module, pre-cleaning): per-region gDNA density by LOCAL
    # boundary-anchored imputation. Needed here only to fit the gDNA strand overdispersion (its seed
    # identification is pre-cleaning) — the cleaning below depends on that overdispersion, so the raw
    # pass must come first. This is NOT the region answer (the sweep below is).
    region_density_raw = region_gdna_density(chain, geometry, region_arrays)

    # Strand-module parameters — the two Beta-Binomial overdispersions. gDNA (mean ½) fitted from the
    # count-observable seed regions/sides using the raw count-clue gDNA weight (breaks the circularity:
    # the seed weight is the strand MEAN ½, not the dispersion). RNA (mean κ) fitted from the PER-SJ
    # SJ strand table — the same strand-qualified population κ itself is the marginal of, so both halves of
    # the RNA Beta-Binomial come from one source. Both shrunk
    # toward the SAME default prior, so under sparse data they collapse to one distribution and an
    # unstranded region (κ=½) is uninformative. See docs/em_strand/03+05.
    _gd_seed = _rna_seed = (
        -1,
        -1,
        False,
    )  # (n_seed_regions, n_seed_frags, fallback) — QC log only; -1 = injected
    if inj is not None and inj.gdna_strand_overdispersion is not None:
        gdna_strand_overdispersion = float(inj.gdna_strand_overdispersion)
    else:
        gdna_strand = fit_gdna_strand_from_substrate(
            substrate,
            region_arrays,
            region_density_raw,
            rna_sense_frac=rna_sense_frac,
        )
        gdna_strand_overdispersion = gdna_strand.gdna_strand_overdispersion
        _gd_seed = (
            gdna_strand.n_seed_regions,
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
        _rna_seed = (
            rna_strand.n_seed_regions,
            rna_strand.n_seed_fragments,
            rna_strand.fallback_used,
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
            fit_intron_background(substrate, region_arrays, region_eff_gdna, include_introns=False)
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
    # because the factor was never evaluated out there. ⚠ Found by the arm raising `IndexError` in
    # `_regrid_global`, which is the honest failure: the array simply was not the shape the solve
    # claimed. ⭐ Memoised, because only two brackets ever occur in one call (the configured one and the
    # derived one) and `_build_intron_prior` is a per-slot NegBinom over the whole grid.
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

    # Message precision is entirely SELF-CONTAINED in the sweep: the source's own honest belief precision
    # (strand + count, from `region_init.build_region_init`), degraded by the reframe's scale variance
    # (σ²_transfer = Var(log r)) and the DerSimonian–Laird composition-mismatch b̂² — all derived from counts and
    # effective lengths inside the pass. There is nothing to fit here.

    # When ``_debug`` is on, the LAST sweep also fills ``_debug["capture"]`` with the per-region message
    # internals (local vs final belief, each channel's message mode/precision) — the substrate for the
    # message-corruption trace (`scripts/debug/msg_trace.py`). Inert in production.
    def _sweep(prior):
        capture = {} if _debug is not None else None
        # ⭐⭐⭐ THE λ BRACKET IS `max(the reference's floor, the fitted prior's own demand)`. ψ evaluates
        # the landscape at `log ρ = log f + log M − log E` and can only offer `f ∈ [σ(−L), σ(L)]`, so a
        # bracket narrower than the prior's own support leaves ψ with no coordinate for what the prior is
        # telling it — and the answer then depends on `L`, which `simplex_logodds` calls its own
        # acceptance test. The demand is DERIVED (`required_logodds_window`), never chosen.
        # ⛔ The prior-free Phase-1 solve has no landscape to widen for, and is L-invariant to seven
        # digits, so it keeps the floor. That is the whole of the conditional.
        # ⛔ `dlam` is held FIXED, so the grids scale with the bracket: widening `L` at a fixed grid size
        # coarsens the lattice and confounds two knobs — measured, that variant REVERSES at L = 40 while
        # this one saturates.
        # ⭐ The TILT axis does NOT scale. `θ` is the RNA-internal SHARE and has no bracket problem, so it
        # stays at the configured resolution; that is what keeps the AMBIG cube ~2.3× rather than ~5.2×.
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
            # ⭐ ψ's reference MEAN from the annotation. OFF ⇒ `location=None` ⇒ no term ⇒ bit-identical.
            structural_reference=config.structural_reference,
            # ⛔⛔⛔ **MESSAGE PROPAGATION IS OFF (owner, 2026-08-07), AND A MEASUREMENT PUT IT THERE.**
            # `SilentPolicy` sends nothing: psi carries each slot's OWN evidence alone — its two strand
            # counts, its spliced count, the derived reference, the fitted gDNA prior and the intron
            # factory. Measured on the 36-condition ladder — RETIRED, rebuilt at 16 conditions on
            # 2026-08-13, so these are as-recorded and not reproducible as written — muting the message
            # layer is a net IMPROVEMENT on THREE OF THE FOUR STRATA:
            #
            #     stranded   x capture ON    -58.3 %   16/16 conditions better
            #     stranded   x capture OFF   -43.7 %   16/16 better
            #     unstranded x capture OFF   -32.1 %   14/16 better
            #     unstranded x capture ON   +154.8 %    0/16 better
            #
            # ⚠ The panel TOTAL is +99.9 % worse, because the one stratum the messages help carries 73 % of
            # the panel's error — so this is a deliberate trade and not a free win. `ROADMAP.md` carries it.
            # ⭐ That stratum is exactly where kappa = 1/2 makes the strand lambda-term exactly 0, so a
            # slot has no own composition evidence at all and a message is the only source there is.
            # ⛔ The planned way out — a theta-independent FRAGMENT-LENGTH channel — was built, measured
            # on the drained arm and DELETED (2026-08-10): its answer is not a function of the length gap
            # at all. `TRAPS.md` carries the mechanism.
            #
            # ⛔ Switching back is ONE WORD, `RelayPolicy()`, and every operator is still there behind its
            # own named switch, so this is reversible and each operator stays individually priceable.
            # ⭐ `message_policy` picks WHICH policy the flag installs: the shipped relay, or the
            # Stage-3 CurrencyPolicy under development. One config value IS the A/B.
            policy=(
                (CurrencyPolicy() if config.message_policy == "currency" else RelayPolicy())
                if config.message_propagation
                else SilentPolicy()
            ),
            _capture=capture,
        )
        if capture is not None:
            _debug["capture"] = capture
        return out

    # THE ENRICHMENT NPMLE — fit ONCE on ALL regions' TOTAL unspliced density (belief-free). It models the
    # hybrid-capture ENRICHMENT/DEPLETION landscape, NOT composition: a total-density prior is
    # composition-vacuous (count-zero-information — /§5), so it is NEVER fed to the
    # composition (gDNA) arm. Its old second role — supplying the message σ²_transfer by projection — is
    # RETIRED (that was a density-uniformity proxy, invalid under capture, and identically 0 in pass-0); the
    # solver now derives σ²_transfer itself. What remains is the QC report's P(ρ) landscape + the toy-injection
    # substrate, so it is fit here and consumed below, never inside the sweep.
    mass_global, eff_global = region_gdna_geometry(geometry)
    if inj is not None and inj.enrichment_prior is not None:
        # INJECTED population enrichment landscape — a toy has too few regions to resolve enriched vs depleted
        # modes. (Scale note: the toy's densities must be generated at the reference library's depth so its
        # regions project onto the right cells of this absolute log-density landscape.)
        enrichment_prior = inj.enrichment_prior
    else:
        enrichment_prior = DensityNPMLE.fit(
            mass_global, eff_global, bandwidth=config.npmle_bandwidth
        )
    # PHASE 1 — the INITIAL solve is PRIOR-FREE of the DNA composition prior: the inert Beta(½,½) reference
    # alone (``gdna_prior=None``) + the strand likelihood + the belief-free forward-backward messages.
    # Single-strand regions self-solve from strand; unstranded
    # AMBIG regions are grounded only by the messages here (the two-root DNA ambiguity,
    # is resolved by the DECONVOLVED-gDNA hyperprior in Phase 2 — fit on this solve's deconvolved DNA, then a refit).
    belief = _sweep(None)
    belief_pass0 = (
        belief  # the initial (prior-free) solve — kept for the refit before/after (movie / debug)
    )
    logger.debug(
        "calibration: PHASE 1 prior-free initial solve (enrichment NPMLE %d cells → QC only)",
        enrichment_prior.n_cells,
    )

    # PHASE 2 — the DECONVOLVED-gDNA hyperprior REFIT. Fit the gDNA-rate NPMLE on
    # the initial solve's deconvolved gDNA, then RE-SOLVE with it as the composition arm — resolving the two-root DNA
    # ambiguity the prior-free pass leaves at unstranded AMBIG regions. Repeated ``calib_refit_iters`` times.
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
            region_eff_gdna,
            include_introns=config.background_include_introns,
            robust_trim_mad=config.background_robust_trim_mad,
        )
    else:
        background = None
    gdna_hyperprior: DensityLandscape | None = None
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
        # the fitted landscape itself, so an over-confident region cannot refuse to budge when the prior lands.
        belief = _init_belief()
        belief = _sweep(gdna_hyperprior)
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
            gdna_prior=enrichment_prior,  # the ENRICHMENT NPMLE (QC landscape / injection substrate)
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
    density_global = gdna_density_global(regions, boundaries, region_eff_gdna, boundary_eff_gdna)

    # The certified-RNA crossings per BOUNDARY: molecules that crossed contiguously having spliced
    # elsewhere. ``chain_boundary_deconv`` adds the whole of this to ``rna_mass`` (rna = (1−g)·unspliced +
    # spliced), so it is exactly the spliced component of ``mass_rna_boundary``. ``assemble_priors``
    # withholds it from ``rna_prior_count`` — a spliced fragment is guaranteed-RNA in the EM (no gDNA
    # candidate), so it must not load the RNA side of the gDNA-vs-RNA *unspliced* split. ``mass_rna_boundary``
    # stays spliced-inclusive so per-boundary conservation gdna + rna = unspliced + spliced holds.
    #
    # ⚠ There is no REGION twin, and that is structural: ``region_contained`` is credited only when the
    # fragment used no sj, so a region's contained population cannot hold a spliced molecule.
    mass_rna_spliced_boundary = np.asarray(substrate.boundary_spliced.count, dtype=np.float64).sum(
        axis=1
    )
    # ⭐ GEOMETRY, not a split: the mean conserved fragment-mass one crossing at this boundary carries.
    # ``assemble_priors`` needs it to turn a per-boundary object-incidence total into a fragment count.
    boundary_mass_per_crossing = substrate.boundary_unspliced.mass_per_crossing

    # ⭐ The JUMPING population, exported verbatim (owner ruling, 2026-07-30). A sj boundary is pure
    # mature RNA by construction, so there is nothing to deconvolve: this is ``sj_count`` summed over
    # the genome-strand columns. ``assemble_priors`` does not read it — it is certified RNA in exactly
    # the sense the spliced crossings are withheld for — but the calibration's output should not be
    # silent about the population that at a donor boundary IS the gene's whole mature output.
    count_rna_sj = np.asarray(substrate.sj.count, dtype=np.float64).sum(axis=1)

    # ⭐ The two remaining INCIDENCE→FRAGMENT conversions, alongside `boundary_mass_per_crossing`. Each is
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
        # ⭐ The simplex ψ solved, published per object rather than summed away. `mass_*` above is this
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
    # ⭐ It is also "sense derived, never stored" in one boundary: the accumulator's columns are GENOME
    # strand, and which of them is *sense* is read off each sj's own annotated transcript strand.
    _flux = np.asarray(substrate.sj.count, dtype=np.float64)
    _is_pos = np.asarray(sj.strand) == np.int8(Strand.POS)
    spl_sense = float(np.where(_is_pos, _flux[:, 0], _flux[:, 1]).sum())
    spl_total = float(_flux.sum())
    sj_sense_frac = spl_sense / spl_total if spl_total > 0.0 else float("nan")
    logger.debug(
        "calibration: N=%d E=%d J=%d gdna_density_global=%.4g rna_sense_frac=%.3f "
        "gdna_strand_overdispersion=%.4g (%d seed regions, %d frags%s) "
        "rna_strand_overdispersion=%.4g (%d sj, %d frags%s) "
        "[sj sense_frac=%.3f vs κ=%.3f]",
        result.n_regions,
        result.n_boundaries,
        result.n_sj,
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
        sj_sense_frac,
        rna_sense_frac,
    )
    return result


__all__ = ["calibrate"]
