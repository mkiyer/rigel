"""rigel.calibration.region_init — the pass-0 per-slot INITIALIZATION (the message-free self-solve).

Calibration's prior-free first pass ("pass-0") deconvolves each slot's unspliced fragment mass into
``(f_pos, f_neg, f_g)`` before any message passing, and records what the slot's OWN data says about that
split, as two things the backbone and the message policy read:

* the composition MODE ``f_*`` — the strand deconvolution (`simplex_logodds`) at every solvable slot; a
  slot that does not deconvolve its own split (no admissible RNA strand, or no counts) keeps its
  signature-binary belief;
* ``tau_lam`` — the slot's own λ-axis composition evidence, the DATA's Fisher information and never a
  prior's, from two sources: the STRAND deconvolution (a Beta-Binomial that is RANK-1, so it informs
  only ``p``; at a single-strand slot the tilt is structurally locked and the strand PINS ``f_g``, at
  an AMBIG slot the tilt is free and the strand cancels out of ``f_g`` — the Schur marginal — so the
  strand term is gated to single-strand slots; identically zero on unstranded data by a derived
  noise-floor deadband) and the INTRON FACTORY (the curvature of the density deconvolution's per-slot
  λ-factor, `density_deconv.density_factor_precision`). The message policy reads ``tau_lam`` as the
  liveness of a node's strand channel; `has_own_composition_evidence` is the instruments' one
  predicate on it.

Structural certainty is not recorded here: the one predicate for a pure-gDNA object is
`region_geometry.g1_locked`, which the instruments import.

Layer: LAYER 6. It imports DOWN to `region_chain` (0), `region_geometry` and `simplex_logodds` (3) and
`density_deconv` (5) — never `sweep`, so it sits cleanly beneath the backbone that consumes
:func:`build_region_init`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .density_deconv import density_factor_precision
from .region_chain import RegionChain
from .region_geometry import RegionGeometry, RegionStatics
from .simplex_logodds import _logodds_grid, _solve_regions_logodds_all

__all__ = [
    "RegionInit",
    "has_own_composition_evidence",
    "strand_discriminability",
    "strand_evidence",
    "build_region_init",
]

_EPS = 1.0e-9


def has_own_composition_evidence(tau_lam) -> np.ndarray:
    """THE ONE DEFINITION of "this slot has own composition evidence", and it lives here so every
    consumer imports it instead of restating the number.

    ``tau_lam`` is the λ-axis Fisher precision summed over the sources (:func:`build_region_init`);
    anything above the divide-by-zero guard is a live channel. That is the whole content — the
    predicate is read off the solver's behaviour, not chosen. The transfer policy's own liveness test
    is ``tau_lam > 0`` (`messages.transfer`); the two agree wherever a positive ``tau_lam`` exceeds the
    guard.

    ⛔ It is NOT a resolving-power test and must not become one. ``τ`` is continuous across the
    interesting region, so a floor on it is a tuned constant. On an unstranded library the strand arm
    carries ``I ≈ Var(κ̂)·N_eff/(p(1−p))`` — roughly the region's depth over the library's spliced
    depth — which is genuinely nonzero and physically nil, and no derivation makes it exactly zero.
    The consumer's defence is a FIXED-DENOMINATOR score, not a tighter bound here
    (``solvability_audit.summarise``'s ``all_mwae`` / ``abs_err``, gated in
    ``test_solvability_audit.py``).

    The home is production, rather than each instrument restating the constant beside a comment saying
    it must match the solver, because the predicate is a production concept and ``scripts/`` is
    deliberately not importable.
    """
    return np.asarray(tau_lam, np.float64) > _EPS


@dataclass(frozen=True, slots=True)
class RegionInit:
    """The per-slot message-free self-solve (length ``n_slots``) — the sweep's starting beliefs.

    ⛔ The first axis is the unified region+boundary CHAIN, not the region axis:
    :func:`build_region_init` reads `RegionGeometry`, whose every array is ``n_slots``, so a reader
    sizing a new array off "n_regions" builds the wrong shape.

    ``f_*`` are the composition MODE (strand / intron / structural-default fractions), read by the
    diagnostics capture; ``tau_lam`` is the combined ``λ``-axis evidence (I_strand + I_factory) — the
    DATA's information, never a prior's, see :func:`build_region_init` — the liveness of a node's
    strand channel, which the message policy reads."""

    f_g: np.ndarray
    f_pos: np.ndarray
    f_neg: np.ndarray
    tau_lam: np.ndarray


# ── source 3: the strand composition evidence (I_strand) + the structural lock ─────────────────────────────


def strand_discriminability(kappa, od_g, od_r, n_gdna_obs, n_rna_obs) -> float:
    """The library's strand DISCRIMINABILITY ``disc = 4·max(0, (κ−½)² − σ²_d)`` — the strand Fisher
    information's library-level factor (:func:`strand_evidence`), with the DERIVED noise floor
    ``σ²_d = ¼·(1/N_rna + ω_r) + ¼·(1/N_gdna + ω_g)``: a κ within √σ²_d of ½ is not composition signal,
    the deadband that kills the unstranded phantom, and ``1/N_gdna`` gates a gDNA-free library
    (N_gdna=0 ⇒ σ²_d→∞ ⇒ disc=0). Exactly zero or strictly positive, and it is the same for every slot:
    a positive ``disc`` is the one condition under which any counted single-strand slot has a live
    strand channel, so it is what the message layer reads as "is the strand split a witness at all"
    (`messages.ChainView.strand_live`). One definition, used by both."""
    sig2_d = 0.25 * (1.0 / max(float(n_rna_obs), _EPS) + od_r) + 0.25 * (
        1.0 / max(float(n_gdna_obs), _EPS) + od_g
    )
    return 4.0 * max(0.0, (kappa - 0.5) ** 2 - sig2_d)


def strand_evidence(u_pos, u_neg, fg_loc, *, kappa, od_g, od_r, n_gdna_obs, n_rna_obs):
    """The reference-free strand composition evidence ``τ₀_λ`` (**I_strand**), evaluated at the
    message-free local ``fg_loc``. Pure; no cross-region coupling.

    ``I_strand(λ) = N_eff·disc·[f_g(1−f_g)]² / (4 p(1−p))``, ``p = κ + f_g(½−κ)`` — the strand Fisher
    information, IDENTICALLY 0 at κ=½ (unstranded). The count enters as the OVERDISPERSED effective count
    ``N_eff = N/(1+(N−1)ω_r)`` (power saturates at ~1/ω, not the raw depth), and the discriminability
    ``disc = 4·max(0, (κ−½)² − σ²_d)`` carries the DERIVED noise floor
    ``σ²_d = ¼·(1/N_rna + ω_r) + ¼·(1/N_gdna + ω_g)`` — a κ within √σ²_d of ½ is not composition signal (the
    deadband that kills the unstranded phantom). ``1/N_gdna`` gates a gDNA-free library (N_gdna=0 ⇒ σ²_d→∞ ⇒
    disc=0).

    Structural composition CERTAINTY is not this function's to declare: the one predicate is
    `region_geometry.g1_locked` (neither RNA strand admissible), and the instruments that classify
    objects by it import that."""
    n_raw = np.asarray(u_pos, np.float64) + np.asarray(u_neg, np.float64)
    n_str = n_raw / (1.0 + np.maximum(n_raw - 1.0, 0.0) * od_r)
    fgl = np.clip(np.asarray(fg_loc, np.float64), _EPS, 1.0 - _EPS)
    pmix = np.clip(kappa + fgl * (0.5 - kappa), _EPS, 1.0 - _EPS)
    disc = strand_discriminability(kappa, od_g, od_r, n_gdna_obs, n_rna_obs)
    return n_str * disc * (fgl * (1.0 - fgl)) ** 2 / (4.0 * pmix * (1.0 - pmix))


# ── source 2 (the intron-factory density-deconv precision `I_density`) lives in `density_deconv.py`
#    (`density_factor_precision`) — a factor's curvature is the density deconvolution's own precision. ──


# ── the assembly ───────────────────────────────────────────────────────────────────────────────────────────


def build_region_init(
    chain: RegionChain,
    statics: RegionStatics,
    geometry: RegionGeometry,
    *,
    kappa: float,
    od_g: float,
    od_r: float,
    n_gdna_obs: float,
    n_rna_obs: float,
    n_grid: int,
    logodds_window: float,
    n_tilt: int | None,
    n_grid_ss: int | None,
    belief,
    priors=None,
    intron_prior=None,
) -> RegionInit:
    """The pass-0 per-slot self-solve → :class:`RegionInit`. Runs the message-free strand deconvolution
    (`simplex_logodds`) and compiles the strand + intron-factory composition evidence ``tau_lam`` (see
    the module docstring).

    The strand deconvolution reference (`fg_ref`/`fpos_ref`/`fneg_ref`) is the incoming ``belief`` — the
    count-zero-information variance freeze evaluates the composition variance near the truth, not at a flat ½.
    ``priors`` (ψ's two composition arms, each ``(m, K)`` or ``None``) and ``intron_prior`` (the intron
    factory ``λ``-factor, ``(m, K)``) enter ψ; ``intron_prior`` additionally seeds I_factory."""
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    # The counts come from the GEOMETRY, their single source: the unspliced ``count`` is both the
    # density numerator and the Poisson n, so there is no second copy to drift.
    count = np.asarray(geometry.unspliced_count, np.float64)
    u_pos, u_neg = count[:, 0], count[:, 1]
    n_region = count.sum(axis=1)
    spliced = np.asarray(geometry.spliced_count, np.float64).sum(axis=1)

    # ── source 3: the message-free strand deconvolution (1-DOF solves; AMBIG partial) ──
    dc = _solve_regions_logodds_all(
        u_pos,
        u_neg,
        fp,
        fn,
        n_region,
        spliced,
        kappa=kappa,
        od_g=od_g,
        od_r=od_r,
        n_grid=int(n_grid),
        L=float(logodds_window),
        n_tilt=n_tilt,
        n_grid_ss=n_grid_ss,
        priors=priors,
        lam_logprior=intron_prior,
        fg_ref=np.asarray(belief.f_g, np.float64),
        fpos_ref=np.asarray(belief.f_pos, np.float64),
        fneg_ref=np.asarray(belief.f_neg, np.float64),
    )
    fg_loc = np.asarray(dc.gdna_frac, np.float64)
    fp_loc = np.asarray(dc.rna_pos_frac, np.float64)
    fn_loc = np.asarray(dc.rna_neg_frac, np.float64)

    # a region that does not deconvolve its own split (G1 sink / empty) keeps the signature-binary init.
    solvable = (fp | fn) & (n_region > 0.0)
    locked = ~solvable
    fg_loc = np.where(locked, np.asarray(belief.f_g, np.float64), fg_loc)
    fp_loc = np.where(locked, np.asarray(belief.f_pos, np.float64), fp_loc)
    fn_loc = np.where(locked, np.asarray(belief.f_neg, np.float64), fn_loc)

    # ── the strand's composition evidence τ_λ (the Schur-marginal gDNA-level precision) ──────────────
    # `strand_evidence` returns the SINGLE-STRAND strand λ-Fisher I_strand = c·a² (the value at the locked tilt).
    i_strand = strand_evidence(
        u_pos,
        u_neg,
        fg_loc,
        kappa=kappa,
        od_g=od_g,
        od_r=od_r,
        n_gdna_obs=n_gdna_obs,
        n_rna_obs=n_rna_obs,
    )
    # The strand Beta-Binomial is
    # RANK-1: it depends on (λ,θ) only through p = ½+(κ−½)(1−f_g)sinθ. So the honest MARGINAL gDNA-level
    # precision (the Schur complement of the 2×2 composition Fisher) is:
    #   * SINGLE-STRAND (1-DOF): θ is STRUCTURALLY locked ⇒ τ_λ gets the full strand λ-term c·a² (strand pins f_g);
    #   * AMBIG (2-DOF): θ is a FREE nuisance ⇒ the strand CANCELS out of f_g (Schur ⇒ 0) — it constrains only
    #     the tilt, never the gDNA level. Crediting c·a² to an AMBIG region is a (bounded) phantom precision on
    #     exactly the regions calibration exists to resolve. Gate the strand λ-term to single-strand regions.
    single_strand = np.asarray(fp, bool) ^ np.asarray(fn, bool)
    tau_lam = np.where(single_strand, i_strand, 0.0)
    lam_grid, _ = _logodds_grid(int(n_grid), float(logodds_window))
    tau_fac = density_factor_precision(
        intron_prior, lam_grid
    )  # I_density (NB curvature) on the λ axis
    if tau_fac is not None:
        tau_lam = tau_lam + tau_fac
    # ⛔ ψ's composition reference (its location term) contributes NOTHING here, and that is deliberate
    # rather than an omission: ``tau_lam`` is the DATA's Fisher information on λ, and the location term
    # carries no count, so crediting it a curvature would flip the own-evidence predicate on slots that
    # gained no measurement.

    return RegionInit(f_g=fg_loc, f_pos=fp_loc, f_neg=fn_loc, tau_lam=tau_lam)
