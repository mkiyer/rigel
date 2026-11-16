"""Calibration's entry point: deconvolve each object's unspliced count into ``(RNA+, RNA−, gDNA)``.

Layer 7. This module owns the order the stages run in and the wiring between them; every number it reports comes from a
module below it. An object's composition has three sources: the binomial strand likelihood of its per-strand counts
(the only intrinsic signal, entering as its Fisher information), the messages from its chain neighbours
(:mod:`.messages`), and the population gDNA prior (the density landscape). The solver is the sweep over the
``N E N E … N`` chain (:mod:`.sweep`)::

    substrate (observation banks on three axes)
      -> chain + geometry + statics               the geometry owns every divisor
      -> strand balance: κ
      -> signature-binary init (G1/G2/G3)
      -> pass 1: solve_chain with no fitted prior, one forward and one backward pass
      -> refits (calib_refit_iters): fit the gDNA landscape on the belief, reset the belief, solve again
      -> chain_region_deconv, chain_boundary_deconv: per-object gDNA / RNA mass
      -> gdna_density_global: the library-average gDNA density (QC)
      -> capture_efficiencies: each object's gDNA density against the last landscape (1 everywhere without one)
      -> the capture weight: the visible share of the gDNA the final solve booked; the reading is raised to it

Gate: ``tests/calibration/test_calibrate.py``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .capture_efficiency import CaptureReading, capture_efficiencies, visible_share, weighted
from .diagnostics import CalibrationDiagnostics, SweepDiagnostics
from .effective_length import observable_cut_shares
from .landscape import DensityLandscape, fit_landscape
from .messages.silent import SilentPolicy
from .messages.transfer import TransferPolicy
from .region_chain import (
    BOUNDARY,
    REGION,
    RegionChain,
    RegionDeconv,
    build_region_chain,
    on_object_axis,
)
from .region_geometry import (
    RegionGeometry,
    RegionStatics,
    build_region_geometry,
    build_region_statics,
    init_beliefs,
    region_gdna_geometry,
    slot_signature,
)
from .result import CalibrationResult
from .signature import RegionType, coarse_type_array
from .strand_balance import StrandBalance, fit_strand_balance
from .substrate import CalibrationSubstrate
from .sweep import chain_boundary_deconv, chain_region_deconv, solve_chain
from ..config import CONSTANTS
from ..types import Strand

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays
    from .splice_graph import SpliceJunctionGeometry

logger = logging.getLogger(__name__)


def _fit_gdna_hyperprior(
    chain, belief, statics, region_arrays, mass_global, eff_global, *, prev=None
):
    """Select the training substrate from the chain and fit the :class:`DensityLandscape` on the deconvolved gDNA.

    The landscape is ψ's gDNA arm in the next refit; it changes the prior only, never the messages. Returns ``None``
    if there is too little substrate to fit. The estimator in :mod:`.landscape` is component-agnostic, so the
    selection is this function's whole job, decided on four separate axes:

    * circularity: ambiguous (both-strand) regions are excluded. They are the two-root ambiguity the prior exists to
      resolve, so training on them would have the prior predict itself.
    * identifiability: a live intronic or intergenic region with no unspliced mass is included as a zero-count
      anchor. "No gDNA here" is the strongest depletion evidence there is, and it grounds the depleted mode.
    * location: a slot whose solve is wider than ``CONSTANTS.landscape.located_var`` in ``log f_g`` (one fragment's
      ``Var(log c) = 1/c``) has no location and is excluded; above that floor, precision is a continuous weight
      (``landscape._reliability``).
    * geometry: boundaries are excluded. They cross rather than contain and are rarely truly enriched, so their
      two-flank mixture would fill the valley between the two modes.
    """
    isr = np.asarray(chain.kind) == REGION
    fp = np.asarray(statics.free_pos, dtype=bool)
    fn = np.asarray(statics.free_neg, dtype=bool)
    slot_type = coarse_type_array(slot_signature(chain, region_arrays.signature))
    constants = CONSTANTS.calibration
    expressed = (
        isr & (eff_global > constants.opportunity_floor) & (mass_global > constants.mass_floor)
    )
    # the zero-count structural anchor: an intergenic or intronic region that sequenced no unspliced mass
    anchor = (
        isr
        & (eff_global > constants.opportunity_floor)
        & (mass_global <= constants.mass_floor)
        & (slot_type != RegionType.EXON)
    )
    sel = expressed & ~(fp & fn)
    # The grid must span every slot the prior is read at, regions and boundaries alike (`fit_landscape`'s
    # ``domain``), so the training cut changes which kernels are summed and never the axis. Otherwise a gDNA-free
    # library whose exons are all blind would train on its anchors alone and collapse the grid to the floor.
    domain_sel = np.asarray(eff_global, dtype=np.float64) > constants.opportunity_floor
    # The guard counts the population the annotation admits, not the one left after the composition cut: guarding
    # the cut population would refuse the refit on a gDNA-free library whose only survivors are a few anchors,
    # leaving its ambiguous regions at their prior-free share. `fit_landscape` still refuses fewer than two.
    if int((sel | anchor).sum()) < constants.min_training_regions:
        return None
    # A slot whose only evidence is a bound, or which has none, does not train the prior (`RegionBelief.has_composition`,
    # which the sweep publishes): its value is where the prior put it, so training on it would re-seed the landscape's
    # tail at the density its own total implies. The anchor trains regardless: it is a structural statement, not a
    # solve.
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
        # the previous refit's landscape: the E-step on the location-free kernels (`landscape.estep_kernels`)
        prev=prev,
        domain=(
            np.asarray(mass_global, dtype=np.float64)[domain_sel],
            np.asarray(eff_global, dtype=np.float64)[domain_sel],
        ),
    )


def lattice_points(window: float, step: float) -> int:
    """The λ lattice's point count for the bracket ``λ ∈ [−L, L]`` at ``step``: ``K = round(2L/step) + 1``.

    The step stays fixed when the landscape widens the bracket (:func:`_sweep`) and the point count follows. The
    bracket and the resolution are separate knobs: widening ``L`` at a fixed count would coarsen the lattice, and a
    result would then move with the bracket instead of saturating, which hides a truncated bracket. Rounding keeps the
    step to within half a point.
    """
    return int(round(2.0 * float(window) / float(step))) + 1


@dataclass(frozen=True, slots=True)
class _SweepInputs:
    """Everything a sweep reads, fixed for the whole calibration: the chain, its statics and geometry, the region
    arrays, the strand model, the policy, the config, and the per-slot gDNA support ``(mass_global, eff_global)`` the
    landscape prior is fit and read on."""

    chain: RegionChain
    statics: RegionStatics
    geometry: RegionGeometry
    region_arrays: RegionArrays
    strand: StrandBalance
    policy: SilentPolicy | TransferPolicy
    config: CalibrationConfig
    mass_global: np.ndarray
    eff_global: np.ndarray


def _sweep_inputs(
    payload,
    substrate,
    region_arrays,
    sj,
    gdna_fl_pmf,
    rna_fl_pmf,
    boundary_flags,
    strand_model,
    config,
) -> _SweepInputs:
    """The chain and its geometry, statics and strand model. The geometry owns every divisor: ``eff_gdna`` /
    ``eff_rna`` are the contained placements at a region slot and the crossing placements at a boundary slot.

    The strand likelihood is binomial for both components, with no overdispersion term: no population measures one
    cleanly on real data (antisense RNA reaches the gDNA seeds, splice artifacts the junctions), and one set too high
    mutes the channel. A library with no spliced read is unstranded, κ = ½ with a dead channel, never an error."""
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geometry = build_region_geometry(chain, substrate, region_arrays, sj, gdna_fl_pmf, rna_fl_pmf)
    strand = fit_strand_balance(strand_model)
    mass_global, eff_global = region_gdna_geometry(geometry)
    return _SweepInputs(
        chain=chain,
        statics=build_region_statics(chain, region_arrays, boundary_flags),
        geometry=geometry,
        region_arrays=region_arrays,
        strand=strand,
        policy=_policy(config, strand),
        config=config,
        mass_global=mass_global,
        eff_global=eff_global,
    )


def _policy(config, strand: StrandBalance):
    """The message-composition policy the config names. An unknown name raises, because an arm that silently ran a
    different policy than it names could not be trusted. The transfer policy reads the library's κ."""
    if config.message_policy == "silent":
        return SilentPolicy()
    if config.message_policy == "transfer":
        return TransferPolicy(kappa=strand.rna_sense_frac)
    raise ValueError(
        f"unknown message_policy {config.message_policy!r} — expected 'silent' or 'transfer'"
    )


def _init_belief(s: _SweepInputs):
    """The signature-binary G1/G2/G3 belief on the chain, before any sweep."""
    return init_beliefs(
        s.geometry,
        s.statics,
        rna_sense_frac=s.strand.rna_sense_frac,
        n_grid=lattice_points(s.config.sweep_logodds_window, s.config.sweep_logodds_step),
        logodds_window=s.config.sweep_logodds_window,
        n_threads=int(s.config.n_threads),
    )


def _sweep(s: _SweepInputs, belief, prior, diagnostics: SweepDiagnostics | None):
    """One sweep of the chain from ``belief`` under the composition prior ``prior`` (``None``: the prior-free pass).

    The λ bracket is ``max(the configured window, the prior's own demand)``. ψ evaluates the landscape at
    ``log ρ = log f + log M − log E`` and can only offer ``f ∈ [σ(−L), σ(L)]``, so a bracket narrower than the prior's
    support would make the answer depend on ``L``. The demand is derived (`required_logodds_window`), never chosen.
    The step stays fixed, so the point count scales with the bracket (:func:`lattice_points`); the tilt axis θ is a
    share with no bracket and does not scale, which keeps the ambiguous cube linear in the bracket. Each message is
    priced inside the sweep from the two nodes' own counts, so nothing is fitted here."""
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
    return solve_chain(
        s.chain,
        s.statics,
        s.geometry,
        belief,
        s.region_arrays,
        rna_sense_frac=s.strand.rna_sense_frac,
        n_rna_obs=float(s.strand.n_observations),
        n_grid=lattice_points(window, cfg.sweep_logodds_step),
        logodds_window=window,
        gdna_prior=prior,
        policy=s.policy,
        block_slots=cfg.sweep_block_slots,
        n_threads=int(cfg.n_threads),
        diagnostics=diagnostics,
    )


def _solve(s: _SweepInputs, diagnostics: SweepDiagnostics | None):
    """The two phases; returns ``(belief, landscape)``.

    Phase 1 solves with no fitted prior: the Beta(½,½) reference, the strand likelihood and the messages.
    Single-strand regions solve from strand; ambiguous regions on unstranded data are grounded only by the messages,
    and their two-root ambiguity is what phase 2 resolves. Phase 2 runs ``calib_refit_iters`` times: fit the
    landscape on the previous solve's deconvolved gDNA (:func:`_fit_gdna_hyperprior`), reset the belief in full so an
    over-confident region cannot resist the prior, and solve again with the landscape as ψ's gDNA arm. Each refit's
    landscape starts the next one's E-step. The messages never read the prior, so on the same grid every refit
    recomputes the same messages; caching them would hold memory through calibration behind a key that must name
    every input the layer reads.

    The landscape is ``None`` if no refit fitted one. Every sweep fills ``diagnostics`` when given one, so it ends
    holding the last."""
    belief = _sweep(s, _init_belief(s), None, diagnostics)
    landscape: DensityLandscape | None = None
    for it in range(int(s.config.calib_refit_iters)):
        fitted = _fit_gdna_hyperprior(
            s.chain,
            belief,
            s.statics,
            s.region_arrays,
            s.mass_global,
            s.eff_global,
            prev=landscape,
        )
        if fitted is None:
            # the belief stands as solved under the previous landscape, so that landscape is the answer's
            break
        landscape = fitted
        belief = _sweep(s, _init_belief(s), landscape, diagnostics)
        logger.debug(
            "calibration: PHASE 2 gDNA-landscape refit %d/%d (%d training regions)",
            it + 1,
            s.config.calib_refit_iters,
            landscape.n_train,
        )
    return belief, landscape


def _result(
    s: _SweepInputs, substrate, sj, belief, landscape, gdna_fl_pmf: np.ndarray
) -> CalibrationResult:
    """Project the solved chain onto the two payload axes and publish the :class:`CalibrationResult`.

    The published gDNA effective lengths are projected off the geometry, never recomputed, so the published divisor
    is the one the solver used. gDNA's conserved share at a boundary (`effective_length.observable_cut_shares`) is not a
    solver divisor: it partitions gDNA's starts over the objects for the EM's gDNA length, and gDNA's template does not
    end, so a side with nothing masked has unbounded reach.

    ``count_rna_spliced_boundary`` is the certified-RNA crossings per boundary (molecules that crossed contiguously
    having spliced elsewhere); `chain_boundary_deconv` adds all of it to ``rna_mass``. A region has no such count:
    its contained population cannot hold a spliced molecule. ``count_rna_sj`` is the jumping population, pure RNA,
    exported verbatim, beside the two regions each sj joins (``sj``'s own, so a consumer can place the sj mass). Each
    ``mass_per_crossing`` converts only its own population's incidences to fragments."""
    n_regions, n_boundaries = int(substrate.n_regions), int(substrate.n_boundaries)
    regions = chain_region_deconv(s.chain, belief, substrate)
    boundaries = chain_boundary_deconv(s.chain, belief, substrate)
    region_eff = on_object_axis(s.chain, s.geometry.eff_gdna, REGION, n_regions)
    boundary_eff = on_object_axis(s.chain, s.geometry.eff_gdna, BOUNDARY, n_boundaries)
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
        sj_src_region=np.asarray(sj.src_region, dtype=np.int64),
        sj_dst_region=np.asarray(sj.dst_region, dtype=np.int64),
        gdna_region_eff_len=region_eff,
        gdna_boundary_eff_len=boundary_eff,
        gdna_boundary_conserved_len=observable_cut_shares(s.region_arrays, gdna_fl_pmf),
        # the simplex ψ solved, published per object; the masses above are the same answer with the
        # two RNA strands added together
        gdna_frac_region=regions.gdna_frac,
        rna_pos_frac_region=regions.rna_pos_frac,
        rna_neg_frac_region=regions.rna_neg_frac,
        gdna_frac_boundary=boundaries.gdna_frac,
        rna_pos_frac_boundary=boundaries.rna_pos_frac,
        rna_neg_frac_boundary=boundaries.rna_neg_frac,
        gdna_density_global=_gdna_density_global(regions, boundaries, region_eff, boundary_eff),
        **_capture_fields(
            landscape, regions, boundaries, region_eff, boundary_eff, _capture_weight(s, belief)
        ),
        rna_sense_frac=s.strand.rna_sense_frac,
        n_regions=n_regions,
        n_boundaries=n_boundaries,
        n_sj=int(substrate.n_sj),
    )


def _gdna_density_global(
    regions: RegionDeconv,
    boundaries: RegionDeconv,
    region_eff: np.ndarray,
    boundary_eff: np.ndarray,
) -> float:
    """The library-average gDNA density, a QC scalar: ``Σ gDNA mass / Σ gDNA effective length`` over both axes.

    A ratio of sums, never a mean of ratios: a rate pooled over unequal supports is not the average of the per-object
    rates. ``0.0`` when there is no support anywhere: an empty library has no density, and a floored division would
    report one.
    """
    total_g = float(
        np.asarray(regions.gdna_mass, dtype=np.float64).sum()
        + np.asarray(boundaries.gdna_mass, dtype=np.float64).sum()
    )
    total_l = float(
        np.asarray(region_eff, dtype=np.float64).sum()
        + np.asarray(boundary_eff, dtype=np.float64).sum()
    )
    return total_g / total_l if total_l > 0.0 else 0.0


def _capture_weight(s: _SweepInputs, belief) -> float:
    """The weight the capture reading is raised to: the visible share (`capture_efficiency.visible_share`) of the gDNA
    the final solve booked, over the chain's slots. A slot's mass is its unspliced count, the support the landscape is
    fitted on, and its booking is ``f_g`` of that; a slot whose annotation admits both RNA strands does not vote."""
    both = np.asarray(s.statics.free_pos, dtype=bool) & np.asarray(s.statics.free_neg, dtype=bool)
    return visible_share(
        np.asarray(belief.f_g, dtype=np.float64) * s.mass_global, s.mass_global, ~both
    )


def _capture_fields(
    landscape: DensityLandscape | None,
    regions: RegionDeconv,
    boundaries: RegionDeconv,
    region_eff: np.ndarray,
    boundary_eff: np.ndarray,
    weight: float,
) -> dict:
    """The result's capture reading: each object's own gDNA density read against the last landscape's bulk, with the
    landscape as a two-group prior and the object's count as the likelihood (`capture_efficiency`), raised to the
    library's capture weight. With no landscape every efficiency is exactly 1, there is no bulk density and the
    saturation level is 1; the weight is published either way."""
    if landscape is None:
        return dict(
            gdna_bulk_density=None,
            gdna_capture_weight=weight,
            gdna_capture_saturation=1.0,
            gdna_capture_efficiency_region=np.ones(region_eff.shape[0]),
            gdna_capture_efficiency_boundary=np.ones(boundary_eff.shape[0]),
        )
    reading: CaptureReading = weighted(
        capture_efficiencies(
            landscape, regions.gdna_mass, region_eff, boundaries.gdna_mass, boundary_eff
        ),
        weight,
    )
    return dict(
        gdna_bulk_density=reading.bulk_density,
        gdna_capture_weight=weight,
        gdna_capture_saturation=reading.saturation,
        gdna_capture_efficiency_region=reading.efficiency_region,
        gdna_capture_efficiency_boundary=reading.efficiency_boundary,
    )


def _log_summary(result: CalibrationResult, strand: StrandBalance, substrate, sj) -> None:
    """Log the sizes, the density scalar, the strand model and the capture reading.

    The sj sense fraction is logged beside κ as QC: a large gap flags a strand-model / accumulator mismatch. Sense is
    derived per sj from its annotated transcript strand, because the accumulator's columns are genome strand."""
    flux = np.asarray(substrate.sj.count, dtype=np.float64)
    is_pos = np.asarray(sj.strand) == np.int8(Strand.POS)
    spl_sense = float(np.where(is_pos, flux[:, 0], flux[:, 1]).sum())
    spl_total = float(flux.sum())
    sj_sense_frac = spl_sense / spl_total if spl_total > 0.0 else float("nan")
    if result.gdna_bulk_density is None:
        logger.info(
            "calibration: no gDNA landscape — every capture efficiency is 1, and effective lengths are not "
            "corrected for capture"
        )
    else:
        logger.info(
            "calibration: capture read against the gDNA bulk at %.3e fragments/bp at weight %.4f (the share of the "
            "reads whose own gDNA would show); the saturated level is %.4g× the bulk at that weight (1: no level "
            "above it)",
            result.gdna_bulk_density,
            result.gdna_capture_weight,
            result.gdna_capture_saturation,
        )
    logger.debug(
        "calibration: N=%d E=%d J=%d gdna_density_global=%.4g rna_sense_frac=%.3f "
        "[sj sense_frac=%.3f vs κ=%.3f]",
        result.n_regions,
        result.n_boundaries,
        result.n_sj,
        result.gdna_density_global,
        strand.rna_sense_frac,
        sj_sense_frac,
        strand.rna_sense_frac,
    )


def _check_sj_axis(sj: "SpliceJunctionGeometry", substrate: CalibrationSubstrate) -> None:
    """Refuse an sj axis whose length disagrees with the payload's: one addressing a different graph would place every
    splice on the wrong boundary."""
    if int(sj.n_sj) != int(substrate.n_sj):
        raise ValueError(
            f"the sj axis has {int(sj.n_sj)} boundaries but the payload has "
            f"{int(substrate.n_sj)}. Build it with "
            "splice_graph.build_sj_geometry_arrays(index) against the SAME index the payload "
            "was scanned on — a sj axis addressing a different graph would place every splice "
            "on the wrong boundary."
        )


def calibrate(
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    gdna_fl_pmf: "np.ndarray",
    rna_fl_pmf: "np.ndarray",
    config: "CalibrationConfig",
    *,
    sj: "SpliceJunctionGeometry",
    boundary_flags: "np.ndarray",
    diagnostics: CalibrationDiagnostics | None = None,
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per object and publish the :class:`CalibrationResult`.

    The stages run in the module docstring's order. ``sj`` is the splice graph's sj axis
    (:func:`~rigel.calibration.splice_graph.build_sj_geometry_arrays`) in the accumulator's sj slot order, and
    ``boundary_flags`` the graph's per-contiguous-boundary structural bits
    (:func:`~rigel.calibration.splice_graph.build_boundary_flags_array`), carried onto the chain as
    ``RegionStatics.boundary_flags``. A zero ``gdna_density_global`` and a zero per-object gDNA mass are valid outputs.

    ``diagnostics``, a record instruments and tests pass, is filled with the solve's intermediates
    (:class:`~.diagnostics.CalibrationDiagnostics`); production passes none.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    _check_sj_axis(sj, substrate)
    inputs = _sweep_inputs(
        payload,
        substrate,
        region_arrays,
        sj,
        gdna_fl_pmf,
        rna_fl_pmf,
        boundary_flags,
        strand_model,
        config,
    )
    belief, landscape = _solve(inputs, None if diagnostics is None else diagnostics.sweep)
    result = _result(inputs, substrate, sj, belief, landscape, gdna_fl_pmf)
    if diagnostics is not None:
        diagnostics.chain = inputs.chain
        diagnostics.geometry = inputs.geometry
        diagnostics.statics = inputs.statics
        diagnostics.substrate = substrate
        diagnostics.sj = sj
        diagnostics.region_arrays = region_arrays
        diagnostics.belief = belief
        diagnostics.landscape = landscape
        diagnostics.rna_sense_frac = inputs.strand.rna_sense_frac
    _log_summary(result, inputs.strand, substrate, sj)
    return result


__all__ = ["calibrate"]
