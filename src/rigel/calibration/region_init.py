"""rigel.calibration.region_init — the pass-0 per-region INITIALIZATION (the message-free self-solve).

Calibration's prior-free first pass ("pass-0") deconvolves each region's unspliced fragment mass into
``(f_pos, f_neg, f_g)``. Before any message passing, every region is given its OWN belief — a per-component
density ``ρ_c`` (the message currency) and a per-component precision ``p_c`` — from the four information
sources below:

1. **MEASURED** counts get **Poisson** precision. Intergenic / intergenic-exon regions are structurally pure
   gDNA (``f_g = 1``, composition CERTAIN); their gDNA density carries only the count precision ``1/n`` (so
   the own precision is the raw count ``n``). This is the anchor the whole prior-free pass leans on.
2. **DENSITY DECONVOLUTION** — a region's gDNA is peeled against a gDNA density prior via NegBinom
   (`density_deconv`); the curvature of that per-region ``λ``-factor is honest, count-derived composition
   evidence, registered here as ``τ_λ`` (`density_factor_precision`). The **intron factory** is its special
   case (the gDNA prior = the intergenic region distribution).
3. **STRAND DECONVOLUTION** — the strand Beta-Binomial is RANK-1 (it informs only ``p``), so its
   gDNA-level (``λ``) precision is the Schur-marginal (approach E): a 1-DOF
   (single-strand) region has its tilt structurally locked, so the strand PINS ``f_g`` (``τ_λ`` gets the strand
   λ-term ``c·a²``); a 2-DOF (AMBIG) region's tilt is free, so the strand CANCELS out of ``f_g`` and contributes
   ZERO (it constrains only the tilt). `strand_evidence` returns the single-strand λ-term; `build_region_init`
   gates it to single-strand regions. Identically zero on unstranded data (κ=½) by a derived noise-floor deadband.
4. **UNSOLVED** regions default to **100 % gDNA at ZERO precision** — the honest "no information" state
   (``τ_λ = 0 ⇒ Var(log f) = ∞ ⇒ p = 0``), left for the sweep + population prior to resolve.

The precision arithmetic (sources → per-component ``Var(log f_c)`` → precision) is pure and unit-tested here.
Layer: imports only lower layers (`region_geometry`, `simplex_logodds`, `density_deconv`), never
`sweep`, so it sits cleanly beneath the backbone that consumes :func:`build_region_init`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .density_deconv import density_factor_precision
from .messages.variance import count_logvar
from .region_chain import REGION, RegionChain
from .region_geometry import RegionGeometry, RegionStatics, region_gdna_geometry
from .simplex_logodds import _logodds_grid, _solve_regions_logodds_all

__all__ = [
    "RegionInit",
    "count_logvar",
    "has_own_composition_evidence",
    "own_composition_logvar",
    "own_precision",
    "strand_evidence",
    "build_region_init",
]

_EPS = 1.0e-9


def has_own_composition_evidence(tau_lam) -> np.ndarray:
    """⭐⭐ **THE ONE DEFINITION of "this slot has own composition evidence", and it lives here so
    every consumer imports it instead of restating the number.**

    ``tau_lam`` is the λ-axis Fisher precision summed over the sources
    (:func:`build_region_init`); anything above the divide-by-zero guard is a live channel, and
    :func:`own_composition_logvar` returns a finite variance there and ``∞`` below it. That is the
    whole content — the predicate is read off the solver's behaviour, not chosen.

    ⛔ **IT IS NOT A RESOLVING-POWER TEST, AND MUST NOT BECOME ONE.** ``τ`` is continuous across the
    interesting region, so a floor on it is a tuned constant: one at ``1/(2L)²`` was derived,
    implemented and refuted by its own insensitivity gate. On an unstranded library the strand arm
    carries ``I ≈ Var(κ̂)·N_eff/(p(1−p))`` — roughly the region's depth over the library's spliced
    depth — which is genuinely nonzero and physically nil, and no derivation makes it exactly zero.
    ⭐ **The consumer's defence is a FIXED-DENOMINATOR score, not a better region_bound**
    (``solvability_audit.summarise``'s ``all_mwae`` / ``abs_err``, gated in
    ``test_solvability_denominator.py``).

    ⚠ Three instruments each restated this as ``_EPS = 1.0e-9`` beside a comment saying it must match
    the solver; changing the solver would have moved none of them. The home is production because the
    predicate is a production concept and ``scripts/`` is deliberately not importable.
    """
    return np.asarray(tau_lam, np.float64) > _EPS


@dataclass(frozen=True, slots=True)
class RegionInit:
    """The per-region message-free self-solve (length ``n_regions``) — the pass-0 relay's starting beliefs.

    ``f_*`` are the composition MODE (strand / intron / structural-default fractions). ``rho_*`` are the own
    per-component DENSITIES ``ρ_c = f_c·M_region/E_c`` (the message currency; a killed / count-free component is
    0). ``prec_*`` are the own per-component PRECISIONS from the four sources — the strand + intron-factory
    composition evidence combined with the Poisson count power (0 where a region is genuinely uninformed).
    ``struct_lock`` marks the structurally composition-certain (pure-gDNA) regions; ``tau_lam`` is the combined
    ``λ``-axis evidence (I_strand + I_factory), retained for tests/diagnostics."""

    f_g: np.ndarray
    f_pos: np.ndarray
    f_neg: np.ndarray
    rho_g: np.ndarray
    rho_pos: np.ndarray
    rho_neg: np.ndarray
    prec_g: np.ndarray
    prec_pos: np.ndarray
    prec_neg: np.ndarray
    struct_lock: np.ndarray
    tau_lam: np.ndarray


# ── the precision arithmetic (pure) ───────────────────────────────────────────────────────────────────────


def own_composition_logvar(f_g, tau_lam, struct_lock):
    """The message-free composition variance of a region's own belief, on the two log-fraction arms — the
    honest "how well does this region know its gDNA/RNA split?".

    Three states, nothing between:

    * **structural lock** → composition CERTAIN → variance ``0`` (an intergenic pure-gDNA region);
    * **real evidence ``τ_λ > 0``** → ``Var(log f_g) = (1−f_g)²/τ_λ`` and ``Var(log f_r) = f_g²/τ_λ`` (the
      Jacobians of ``log f_g`` / ``log f_r`` w.r.t. ``λ = logit f_g``);
    * **no evidence ``τ_λ = 0``** → composition UNSEEN → variance ``∞`` (⇒ own precision 0 — the unsolved
      default). ``∞`` is set DIRECTLY (never ``fr²·∞``, which would be ``0·∞ = nan`` at the window boundary).

    Returns ``(v_log_fg, v_log_fr)``, same shape as ``f_g``."""
    fg = np.clip(np.asarray(f_g, np.float64), 0.0, 1.0)
    fr = 1.0 - fg
    tau = np.asarray(tau_lam, np.float64)
    lock = np.asarray(struct_lock, bool)
    # ⭐ the predicate is :func:`has_own_composition_evidence` — ONE definition, imported by the
    #   instruments rather than restated, so the solver and everything that reports on it cannot drift.
    seen = has_own_composition_evidence(tau)
    with np.errstate(divide="ignore", invalid="ignore"):
        v_fg = np.where(lock, 0.0, np.where(seen, fr * fr / np.maximum(tau, _EPS), np.inf))
        v_fr = np.where(lock, 0.0, np.where(seen, fg * fg / np.maximum(tau, _EPS), np.inf))
    return v_fg, v_fr


def own_precision(n, v_log, live):
    """The own-belief precision of one component: ``p = 1/(Var(log f) + Var(log ρ_count))`` — the
    composition variance combined with the count's own.

    ⭐ The count term is :func:`count_logvar`, the EXACT Poisson log-rate variance, in place of the
    ``1/n`` asymptote it replaces. That deletes the ``n > 0`` branch outright: a zero count is a
    measurement with finite precision, so there is nothing to special-case.

    It is still 0 where there is genuinely nothing to say — no evidence (``Var(log f) = ∞``) or a
    component that is not ``live`` (structurally inadmissible). ⭐ **Those two remain and they are the
    point of the distinction**: ignorance and impossibility both silence a source, a zero count does
    not. No ``0·∞`` nan — the ``∞`` is masked to a finite value BEFORE the sum, matching ``np.where``'s
    both-branch evaluation."""
    v = np.asarray(v_log, np.float64)
    ok = np.isfinite(v)
    v_fin = np.where(ok, v, 0.0)
    return np.where(
        ok & np.asarray(live, bool),
        1.0 / (np.maximum(v_fin, 0.0) + count_logvar(n)),
        0.0,
    )


# ── source 3: the strand composition evidence (I_strand) + the structural lock ─────────────────────────────


def strand_evidence(
    u_pos, u_neg, fg_loc, *, kappa, od_g, od_r, n_gdna_obs, n_rna_obs, is_region, locked
):
    """The reference-free strand composition evidence ``τ₀_λ`` (**I_strand**) + the structural-certainty mask,
    evaluated at the message-free local ``fg_loc``. Pure; no cross-region coupling.

    ``I_strand(λ) = N_eff·disc·[f_g(1−f_g)]² / (4 p(1−p))``, ``p = κ + f_g(½−κ)`` — the strand Fisher
    information, IDENTICALLY 0 at κ=½ (unstranded). The count enters as the OVERDISPERSED effective count
    ``N_eff = N/(1+(N−1)ω_r)`` (power saturates at ~1/ω, not the raw depth), and the discriminability
    ``disc = 4·max(0, (κ−½)² − σ²_d)`` carries the DERIVED noise floor
    ``σ²_d = ¼·(1/N_rna + ω_r) + ¼·(1/N_gdna + ω_g)`` — a κ within √σ²_d of ½ is not composition signal (the
    deadband that kills the unstranded phantom). ``1/N_gdna`` gates a gDNA-free library (N_gdna=0 ⇒ σ²_d→∞ ⇒
    disc=0).

    ``struct_lock`` (**I_struct**) — composition CERTAIN — is scoped to true intergenic REGION regions, never a
    G1 BOUNDARY boundary (TSS/TES): a boundary is structurally gDNA but sits between RNA-carrying exons, so its
    crossing mass is RNA-contaminated and a certainty there compounds into a phantom-gDNA emitter; a true
    intergenic region carries ~0 mass in a zero-gDNA library, so it is safe."""
    n_raw = np.asarray(u_pos, np.float64) + np.asarray(u_neg, np.float64)
    n_str = n_raw / (1.0 + np.maximum(n_raw - 1.0, 0.0) * od_r)
    fgl = np.clip(np.asarray(fg_loc, np.float64), _EPS, 1.0 - _EPS)
    pmix = np.clip(kappa + fgl * (0.5 - kappa), _EPS, 1.0 - _EPS)
    sig2_d = 0.25 * (1.0 / max(float(n_rna_obs), _EPS) + od_r) + 0.25 * (
        1.0 / max(float(n_gdna_obs), _EPS) + od_g
    )
    disc = 4.0 * max(0.0, (kappa - 0.5) ** 2 - sig2_d)
    i_strand = n_str * disc * (fgl * (1.0 - fgl)) ** 2 / (4.0 * pmix * (1.0 - pmix))
    struct_lock = np.asarray(locked, bool) & np.asarray(is_region, bool)
    return i_strand, struct_lock


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
    global_logprior=None,
    intron_prior=None,
) -> RegionInit:
    """The pass-0 per-region self-solve → :class:`RegionInit`. Runs the message-free strand deconvolution
    (`simplex_logodds`), compiles the strand + intron-factory composition evidence, and assembles each region's
    own per-component density + precision from the four sources (see the module docstring).

    The strand deconvolution reference (`fg_ref`/`fpos_ref`/`fneg_ref`) is the incoming ``belief`` — the
    count-zero-information variance freeze evaluates the composition variance near the truth, not at a flat ½.
    ``global_logprior`` (the anchored population gDNA prior, ``(m, K)``) and ``intron_prior`` (the intron
    factory ``λ``-factor, ``(m, K)``) enter ψ; ``intron_prior`` additionally seeds I_factory."""
    is_reg = np.asarray(chain.kind) == REGION
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    # ⭐ The counts come from the GEOMETRY, which is their single source since S5.e: the unspliced
    # ``count`` is both the density numerator and the Poisson n, so there is no second copy to drift.
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
        global_logprior=global_logprior,
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

    # ── sources 1 & 3: the composition evidence τ_λ (the Schur-marginal gDNA-level precision) + struct lock ──
    # `strand_evidence` returns the SINGLE-STRAND strand λ-Fisher I_strand = c·a² (the value at the locked tilt).
    i_strand, struct_lock = strand_evidence(
        u_pos,
        u_neg,
        fg_loc,
        kappa=kappa,
        od_g=od_g,
        od_r=od_r,
        n_gdna_obs=n_gdna_obs,
        n_rna_obs=n_rna_obs,
        is_region=is_reg,
        locked=locked,
    )
    # APPROACH E (verified). The strand Beta-Binomial is
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

    # ── the own per-component densities + precisions — ONE set of numbers per slot, no faces to pool ──
    mass_global, eff_global = region_gdna_geometry(geometry)
    eff_r = np.asarray(geometry.eff_rna, np.float64)
    v_log_fg, v_log_fr = own_composition_logvar(fg_loc, tau_lam, struct_lock)

    # ⭐⭐ THE LOCATION IS UNCHANGED AND THAT IS DELIBERATE — only the PRECISION was ever broken.
    #    ``rho_c = f_c*M/E_c`` keeps the composition identity ``sum_c rho_c*E_c = M`` exactly, which is
    #    what the relay's mass pin exists to enforce; and the relay fuses in LINEAR density space
    #    (the scan's inverse-variance fuse), so ``rho = 0`` is perfectly expressible. A zero density was never the
    #    problem — an INFINITE precision on it was (TRAPS: a-zero-count-is-a-measurement).
    #    ⛔ A first version of this fix also moved the location to the ``Gamma(a+½, E)`` posterior mean
    #    ``(a+½)/E``. That is right for one rate in isolation and WRONG here: three components each
    #    gaining ``+½`` breaks ``sum_c rho_c*E_c = M`` by exactly 3/2, which `test_relay_mass_pin` caught
    #    as ``R_own = 0.5 + 1/M``. The half belongs to the rate's VARIANCE, not to a share of a total.
    #
    # ⛔ What DID change is the ``live`` predicate: STRUCTURAL (opportunity, strand admissibility), never
    #    the count. A zero count is a measurement; only a zero opportunity is an absence of data.
    #
    # gDNA (source 1 measured / source 3 strand / source 4 default all flow through here). It is
    # genomically continuous, so it is admissible wherever it has opportunity — there is no gDNA analogue
    # of a forbidden strand.
    rho_g = np.maximum(
        np.where(
            eff_global > 0.0,
            fg_loc * mass_global / np.where(eff_global > 0.0, eff_global, 1.0),
            0.0,
        ),
        0.0,
    )
    prec_g = own_precision(fg_loc * mass_global, v_log_fg, eff_global > 0.0)

    # RNA per strand — admissible iff the annotation allows that strand here AND there is opportunity.
    # ⚠ ``np.isfinite(var_loc)`` used to gate this too and is GONE as redundant, not as a change of
    # meaning: ``own_precision`` already returns exactly 0 on a non-finite composition variance, so the
    # unresolvable region emitted nothing before and emits nothing now.
    def _rna(frac_loc, free_s):
        admissible = np.asarray(free_s, bool)
        count = np.asarray(frac_loc, np.float64) * mass_global
        live = admissible & (eff_r > 0.0)
        return (
            np.where(live, count / np.where(live, eff_r, 1.0), 0.0),
            own_precision(count, v_log_fr, live),
        )

    rho_pos, prec_pos = _rna(fp_loc, fp)
    rho_neg, prec_neg = _rna(fn_loc, fn)

    return RegionInit(
        f_g=fg_loc,
        f_pos=fp_loc,
        f_neg=fn_loc,
        rho_g=rho_g,
        rho_pos=rho_pos,
        rho_neg=rho_neg,
        prec_g=prec_g,
        prec_pos=prec_pos,
        prec_neg=prec_neg,
        struct_lock=struct_lock,
        tau_lam=tau_lam,
    )
