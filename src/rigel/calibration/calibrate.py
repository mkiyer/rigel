"""calibrate() — the calibrator: deconvolve each object's UNSPLICED count into ``(RNA₊, RNA₋, gDNA)``.

Calibration models **RNA vs gDNA only** (the per-locus EM separates nascent from mature downstream). Per the
count-zero-information architecture a node's composition is set by three sources — the **strand likelihood**
(the Beta-Binomial tilt of the per-strand counts, the only INTRINSIC signal; the count enters ONLY as its
overdispersed Fisher information), the **cross-node imputation** (the neighbour density messages, at the
source's own honest belief precision (strand + count) ``pr = n_src/(n_src·vb_src + 1)``), and the
**population gDNA prior** (the conservative intergenic+intron floor + the Phase-2 density KDE). The solver
is the belief-propagation SWEEP over the ``N E N E … N`` chain (:mod:`rigel.calibration.bp_solver`)::

    substrate  (five populations on three axes)
      -> build chain + geometry + statics      (the geometry owns EVERY divisor)
      -> strand balance: rna_sense_frac (κ)
      -> node_gdna_density (RAW) -> fit gDNA / RNA strand Beta-Binomial overdispersions (seed)
      -> signature-binary init (G1/G2/G3)
      -> PASS 1 node_sweep (no KDE prior): ONE forward + ONE backward pass — each object integrates its
           strand likelihood + the conservative gDNA floor + the belief-free density messages
      -> train the Phase-2 gDNA-density KDE on the pass-1 belief
      -> PASS 2 node_sweep (the KDE prior added per object) -> the converged per-object pie
      -> chain_node_deconv  -> per-NODE gDNA / RNA contained mass
      -> chain_edge_deconv  -> per-EDGE gDNA / RNA crossing mass, for the per-locus prior
      -> gdna_density_global (the library-average density QC scalar)

⭐ **THE GEOMETRY IS BUILT FIRST, AND IT OWNS EVERY DIVISOR.** The predecessor computed three effective
lengths up front (``region_eff_length``, ``boundary_side_eff_length``, ``boundary_eff_length``) and
handed them around; the density model took a fourth, a scalar ``fl_mean``. That is four chances for a
consumer to divide by a length model the solver never used. ``build_node_geometry`` now produces the
per-slot ``eff_gdna``/``eff_rna`` before anything reads a count, and everything downstream — the density
clue, the background pool, the intron factory, the result's two supports — reads THAT array


⛔ **The per-face boundary machinery is gone.** ``boundary_substrate``, the ``left``/``right`` views and
``boundary_side_eff_length``'s ``E[min(ℓ,L)]/2`` had no successor: a contiguous edge is a 0-bp line with
one count and one divisor, ``crossing_eff_length``. A zero-gDNA library
(``gdna_density_global == 0``, per-object gDNA mass ``0``) remains a valid, graceful output.

⚠ **A7 IS OFF AT CONTIGUOUS EDGES, AND THE FIRST BASELINE CARRIES ITS BIAS.** The RNA half of an
unspliced crossing takes ``UNBOUNDED_REACH`` rather than its transcript's real remaining length
(owner-ruled). Cost, already measured: an **11.0 %** genome-wide gDNA
over-call and **+0.36** in the last node before a polyA site. It is
sequenced after this step so that turning it on can be A/B'd against the baseline this step produces.
Junction edges DO take their real exonic reach — they are a new population with no predecessor divisor.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .background_reference import BackgroundReference, measure_background
from .bp_solver import (
    EDGE,
    NODE,
    build_node_geometry,
    build_node_statics,
    chain_edge_deconv,
    chain_node_deconv,
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
from ..types import Strand

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays
    from .splice_graph import JunctionGeometry

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


def _empty_junction_geometry() -> "JunctionGeometry":
    """A junction axis with no rows — a graph whose references are all single-exon.

    ⚠ Legal, and **not** the same as "no junction flux": a payload scanned against such a graph has
    ``n_sj == 0``, so an empty axis is the only consistent input and the alignment check below will say
    so if it is not.
    """
    from .splice_graph import JunctionGeometry

    e_i = np.zeros(0, dtype=np.int64)
    return JunctionGeometry(
        src_node=e_i,
        dst_node=e_i.copy(),
        strand=np.zeros(0, dtype=np.int8),
        reach_lo=np.zeros(0, dtype=np.float64),
        reach_hi=np.zeros(0, dtype=np.float64),
    )


def _project_eff(chain, eff_slots, payload) -> tuple[np.ndarray, np.ndarray]:
    """Split a per-SLOT divisor back onto the node and contiguous-edge axes.

    ⭐ **A projection, not a recomputation.** ``build_node_geometry`` already applied
    ``contained_eff_length`` at NODE slots and ``crossing_eff_length`` at EDGE slots; calling those
    again here would put two implementations of one quantity in the tree, which is precisely how the
    prose and the code came to disagree about a ½ for months ( and 27).
    Whatever the solver divided by is what the result reports.

    ⚠ It takes the ARRAY, not the geometry: ``assemble_priors`` needs the RNA divisor as well as the
    gDNA one (a mass becomes a density only against its own component's opportunity), and one function
    called twice is the alternative to a second copy that drifts.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    eff = np.asarray(eff_slots, dtype=np.float64)
    is_node, is_edge = kind == NODE, kind == EDGE
    node_eff = np.zeros(int(payload.n_nodes), dtype=np.float64)
    edge_eff = np.zeros(int(payload.n_edges), dtype=np.float64)
    node_eff[obj[is_node]] = eff[is_node]
    edge_eff[obj[is_edge]] = eff[is_edge]
    return node_eff, edge_eff


def _build_length_loglik(chain, geometry, region_arrays, gdna_fl_pmf, rna_fl_pmf, config):
    """The per-slot FRAGMENT-LENGTH log-likelihood on the ψ solve grid → ``(n_slots, K)``, or ``None``.

    ⭐ **The fourth information source, finally wired**. The accumulator
    has measured ``inv_length_sum`` and ``length_sum`` on every population since S5.a and nothing read
    them; `length_likelihood` turns them into evidence about ``f_g`` on the same ``λ`` grid the strand
    term lives on, and `node_init` registers its curvature as ``I_length``.

    ⚠ **Built here, once, for the same reason ``intron_prior`` is**: this is the only layer holding the
    chain, the geometry AND both fitted pmfs at the same time, and building it twice would put two
    implementations of one quantity in the tree.

    ⚠ **Strand-agnostic**: the two genome-strand columns are summed. Which strand a read aligned to says
    nothing about whether its molecule was gDNA or RNA; the strand Beta-Binomial keeps its own columns.

    ⛔ Returns ``None`` when the switch is off, which is byte-identical to the pre-P2 path.
    """
    if not getattr(config, "length_likelihood", False):
        return None
    from .length_likelihood import build_slot_moments, length_loglik

    _, fg_grid = _logodds_grid(int(config.sweep_n_grid), float(config.sweep_logodds_window))
    return length_loglik(
        build_slot_moments(chain, region_arrays, gdna_fl_pmf),
        build_slot_moments(chain, region_arrays, rna_fl_pmf),
        np.asarray(geometry.unspliced_count, np.float64).sum(axis=1),
        np.asarray(geometry.unspliced_inv_length_sum, np.float64).sum(axis=1),
        np.asarray(geometry.unspliced_length_sum, np.float64).sum(axis=1),
        fg_grid,
    )


def _build_intron_prior(chain, substrate, region_arrays, node_eff_len, config, bg=None):
    """The gDNA intron factory λ-factor per chain slot.

    Fits the intergenic-background NegBinom (`fit_intron_background`) and tabulates, for each INTRON NODE
    slot, ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` over the σ(λ) solve grid → an ``(n_slots, K)`` array,
    ZERO on every other slot (a no-op there). Returns ``None`` when the factory is disabled, the
    background pool is uninformative, or there are no intron nodes — in which case the sweep is
    byte-identical to the pre-factory path. gDNA is strand-symmetric, so this factor lives purely on
    ``λ`` (peels gDNA; the residual RNA's tilt is left to the solver), and is consumed identically by the
    single-strand and AMBIG per-node solves.

    ⚠ **EDGE slots are zero, structurally.** The factor scores a CONTAINED count against a CONTAINED
    support; a line's count is a crossing and its divisor a different formula, so applying the same
    NegBinom there would be scoring one frame's evidence against another frame's support.

    ``bg`` (an injected population :class:`GdnaBackground`) overrides the internal fit — a tiny toy's own
    intergenic pool is too sparse to fit the background the introns are peeled against."""
    if bg is None:
        bg = fit_intron_background(substrate, region_arrays, node_eff_len, include_introns=False)
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
    # ⚠ GENOME-strand columns summed: gDNA is strand-symmetric, so the peel is against a total rate.
    count = np.asarray(substrate.node_contained.count, dtype=np.float64).sum(axis=1)[ridx]
    eff_g = np.asarray(node_eff_len, dtype=np.float64)[ridx]
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
    junctions: "JunctionGeometry | None" = None,
    edge_rna_reach=None,
    _debug: dict | None = None,
    diagnostics_out: dict | None = None,
    injected_priors: "InjectedCalibrationPriors | None" = None,
    edge_flags: "np.ndarray | None" = None,
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per object, then derive gdna_density_global.

    Runs the belief-propagation sweep (a single forward-backward pass per phase, resolving the
    per-object pie); see the module docstring for the data flow. ``gdna_density_global`` may be ``0``
    (a zero-gDNA library) and an object's deconvolved gDNA mass may be ``0`` (a pure-RNA object); both
    are valid, graceful outputs — not failures.

    ``junctions`` is the splice graph's junction axis
    (:func:`~rigel.calibration.splice_graph.build_junction_geometry_arrays`), in the accumulator's own
    junction slot order — where each junction attaches, its transcript strand, and its exonic reach.
    ``None`` means "this library's graph has no junction edges", which is legal (a single-exon-only
    reference) and is NOT the same as "no junction flux".

    ``edge_flags`` is the splice graph's per-contiguous-edge structural bits
    (:func:`~rigel.calibration.splice_graph.build_edge_flags_array`), carried onto the chain as
    ``NodeStatics.edge_flags``. ``None`` (the default) leaves them zero.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    inj = injected_priors  # population-scale priors to inject in place of the internal (toy-untrustworthy) fits
    junctions = _empty_junction_geometry() if junctions is None else junctions
    if int(junctions.n_junctions) != int(substrate.n_junctions):
        raise ValueError(
            f"the junction axis has {int(junctions.n_junctions)} edges but the payload has "
            f"{int(substrate.n_junctions)}. Build it with "
            "splice_graph.build_junction_geometry_arrays(index) against the SAME index the payload "
            "was scanned on — a junction axis addressing a different graph would place every splice "
            "on the wrong line."
        )

    # ⭐ THE CHAIN AND ITS GEOMETRY COME FIRST, and the geometry owns every divisor from here down:
    # `eff_gdna`/`eff_rna` are the CONTAINED placements at a NODE slot and the CROSSING placements at an
    # EDGE slot, one rule, one array. Nothing below computes a second length model.
    chain = build_node_chain(payload.ref_node_offsets, payload.ref_edge_offsets)
    geometry = build_node_geometry(
        chain, substrate, region_arrays, junctions, gdna_fl_pmf, rna_fl_pmf, edge_rna_reach
    )
    statics = build_node_statics(chain, region_arrays, edge_flags)

    # The result's two gDNA supports, PROJECTED off the geometry rather than recomputed — so the number
    # `priors` divides by is byte-identically the one the solver divided by. Two implementations
    # of one quantity is how they come to disagree.
    node_eff_gdna, edge_eff_gdna = _project_eff(chain, geometry.eff_gdna, payload)
    # ⭐ And the RNA twin, on the same two axes. It is what turns ``mass_rna_*`` into a DENSITY in
    # ``assemble_priors``; without it the RNA side had no divisor at all and the prior's g:r ratio
    # carried the ``Σ A_g / Σ A_r`` length tilt.
    node_eff_rna, edge_eff_rna = _project_eff(chain, geometry.eff_rna, payload)

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

    # Count clue on RAW counts (the count module, pre-cleaning): per-node gDNA density by LOCAL
    # edge-anchored imputation. Needed here only to fit the gDNA strand overdispersion (its seed
    # identification is pre-cleaning) — the cleaning below depends on that overdispersion, so the raw
    # pass must come first. This is NOT the node answer (the sweep below is).
    node_density_raw = node_gdna_density(chain, geometry, region_arrays)

    # Strand-module parameters — the two Beta-Binomial overdispersions. gDNA (mean ½) fitted from the
    # count-observable seed regions/sides using the raw count-clue gDNA weight (breaks the circularity:
    # the seed weight is the strand MEAN ½, not the dispersion). RNA (mean κ) fitted from the PER-JUNCTION
    # SJ strand table — the same strand-qualified population κ itself is the marginal of, so both halves of
    # the RNA Beta-Binomial come from one source. Both shrunk
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

    # Strand-Fisher noise-floor SAMPLE SIZES (bp_solver τ seed): N_gdna (gDNA-eligible unspliced fragments in
    # the structurally pure-gDNA intergenic nodes, coarse type 0) and N_spliced (the pure-RNA count κ_RNA was
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
            np.asarray(substrate.node_contained.count, dtype=np.float64)[_inter].sum()
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
    # The gDNA INTRON FACTORY λ-factor: peel confident gDNA
    # from intron nodes against the intergenic background, BEFORE the pass-0 solve. Built ONCE (belief-free —
    # only the intron count vs the background), applied in every sweep below. ``None`` (disabled / no
    # informative background / no introns) ⇒ byte-identical to the pre-factory pass-0.
    # INJECTED intergenic intron-factory background overrides the internal (toy-sparse) fit.
    intron_background = (
        inj.intron_background
        if (inj is not None and inj.intron_background is not None)
        else (
            fit_intron_background(substrate, region_arrays, node_eff_gdna, include_introns=False)
            if config.intron_factory
            else None
        )
    )
    intron_prior = (
        _build_intron_prior(
            chain, substrate, region_arrays, node_eff_gdna, config, bg=intron_background
        )
        if (config.intron_factory and intron_background is not None)
        else None
    )
    length_loglik_arr = _build_length_loglik(
        chain, geometry, region_arrays, gdna_fl_pmf, rna_fl_pmf, config
    )
    # Message precision is entirely SELF-CONTAINED in the sweep: the source's own honest belief precision
    # (strand + count, from `node_init.build_node_init`), degraded by the reframe's scale variance
    # (σ²_transfer = Var(log r)) and the DerSimonian–Laird composition-mismatch b̂² — all derived from counts and
    # effective lengths inside the pass. There is nothing to fit here.

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
            length_loglik=length_loglik_arr,
            _capture=capture,
        )
        if capture is not None:
            _debug["capture"] = capture
        return out

    # THE ENRICHMENT NPMLE — fit ONCE on ALL nodes' TOTAL unspliced density (belief-free). It models the
    # hybrid-capture ENRICHMENT/DEPLETION landscape, NOT composition: a total-density prior is
    # composition-vacuous (count-zero-information — /§5), so it is NEVER fed to the
    # composition (gDNA) arm. Its old second role — supplying the message σ²_transfer by projection — is
    # RETIRED (that was a density-uniformity proxy, invalid under capture, and identically 0 in pass-0); the
    # solver now derives σ²_transfer itself. What remains is the QC report's P(ρ) landscape + the toy-injection
    # substrate, so it is fit here and consumed below, never inside the sweep.
    mass_global, eff_global = node_global_geometry(geometry)
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
    # AMBIG nodes are grounded only by the messages here (the two-root DNA ambiguity,
    # is resolved by the DECONVOLVED-gDNA hyperprior in Phase 2 — fit on this solve's peeled DNA, then a refit).
    belief = _sweep(None)
    belief_pass0 = (
        belief  # the initial (prior-free) solve — kept for the refit before/after (movie / debug)
    )
    logger.debug(
        "calibration: PHASE 1 prior-free initial solve (enrichment NPMLE %d cells → QC only)",
        enrichment_prior.n_cells,
    )

    # PHASE 2 — the DECONVOLVED-gDNA hyperprior REFIT. Fit the gDNA-rate NPMLE on
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
            node_eff_gdna,
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

    nodes = chain_node_deconv(chain, belief, substrate)
    edges = chain_edge_deconv(chain, belief, substrate)

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
            junctions=junctions,
            region_arrays=region_arrays,
            gdna_prior=enrichment_prior,  # the ENRICHMENT NPMLE (QC landscape / injection substrate)
            gdna_hyperprior=gdna_hyperprior,  # the DECONVOLVED-gDNA composition hyperprior (None if no refit)
            rna_sense_frac=rna_sense_frac,
            node_eff_gdna=node_eff_gdna,
            edge_eff_gdna=edge_eff_gdna,
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
    density_global = gdna_density_global(nodes, edges, node_eff_gdna, edge_eff_gdna)

    # The certified-RNA crossings per LINE: molecules that crossed contiguously having spliced
    # elsewhere. ``chain_edge_deconv`` adds the whole of this to ``rna_mass`` (rna = (1−g)·unspliced +
    # spliced), so it is exactly the spliced component of ``mass_rna_edge``. ``assemble_priors``
    # withholds it from ``rna_prior_count`` — a spliced fragment is guaranteed-RNA in the EM (no gDNA
    # candidate), so it must not load the RNA side of the gDNA-vs-RNA *unspliced* split. ``mass_rna_edge``
    # stays spliced-inclusive so per-line conservation gdna + rna = unspliced + spliced holds.
    #
    # ⚠ There is no NODE twin, and that is structural: ``node_contained`` is credited only when the
    # fragment used no junction, so a node's contained population cannot hold a spliced molecule.
    mass_rna_spliced_edge = np.asarray(substrate.edge_spliced.count, dtype=np.float64).sum(axis=1)

    # ⭐ The JUMPING population, exported verbatim (owner ruling, 2026-07-30). A junction edge is pure
    # mature RNA by construction, so there is nothing to deconvolve: this is ``sj_count`` summed over
    # the genome-strand columns. ``assemble_priors`` does not read it — it is certified RNA in exactly
    # the sense the spliced crossings are withheld for — but the calibration's output should not be
    # silent about the population that at a donor seam IS the gene's whole mature output.
    mass_rna_junction = np.asarray(substrate.junction.count, dtype=np.float64).sum(axis=1)

    result = CalibrationResult(
        mass_gdna_node=nodes.gdna_mass,
        mass_rna_node=nodes.rna_mass,
        mass_gdna_edge=edges.gdna_mass,
        mass_rna_edge=edges.rna_mass,
        mass_rna_spliced_edge=mass_rna_spliced_edge,
        mass_rna_junction=mass_rna_junction,
        gdna_node_eff_len=node_eff_gdna,
        gdna_edge_eff_len=edge_eff_gdna,
        rna_node_eff_len=node_eff_rna,
        rna_edge_eff_len=edge_eff_rna,
        gdna_density_global=density_global,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_nodes=int(substrate.n_nodes),
        n_edges=int(substrate.n_edges),
        n_junctions=int(substrate.n_junctions),
        config=config,
    )
    # Diagnostic: the JUNCTION sense fraction should agree with the StrandModel κ. A large gap flags a
    # strand-model / accumulator mismatch (κ stays the StrandModel posterior — this is QC only).
    # ⭐ It is also "sense derived, never stored" in one line: the accumulator's columns are GENOME
    # strand, and which of them is *sense* is read off each junction's own annotated transcript strand.
    _flux = np.asarray(substrate.junction.count, dtype=np.float64)
    _is_pos = np.asarray(junctions.strand) == np.int8(Strand.POS)
    spl_sense = float(np.where(_is_pos, _flux[:, 0], _flux[:, 1]).sum())
    spl_total = float(_flux.sum())
    junction_sense_frac = spl_sense / spl_total if spl_total > 0.0 else float("nan")
    logger.debug(
        "calibration: N=%d E=%d J=%d gdna_density_global=%.4g rna_sense_frac=%.3f "
        "gdna_strand_overdispersion=%.4g (%d seed nodes, %d frags%s) "
        "rna_strand_overdispersion=%.4g (%d junctions, %d frags%s) "
        "[junction sense_frac=%.3f vs κ=%.3f]",
        result.n_nodes,
        result.n_edges,
        result.n_junctions,
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
        junction_sense_frac,
        rna_sense_frac,
    )
    return result


__all__ = ["calibrate"]
