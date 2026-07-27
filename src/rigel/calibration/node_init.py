"""rigel.calibration.node_init — the pass-0 per-node INITIALIZATION (the message-free self-solve).

Calibration's prior-free first pass ("pass-0") deconvolves each node's unspliced fragment mass into
``(f_pos, f_neg, f_g)``. Before any message passing, every node is given its OWN belief — a per-component
density ``ρ_c`` (the message currency) and a per-component precision ``p_c`` — from the four information
sources of ``docs/calibration/variance_model_concepts.md``:

1. **MEASURED** counts get **Poisson** precision. Intergenic / intergenic-exon nodes are structurally pure
   gDNA (``f_g = 1``, composition CERTAIN); their gDNA density carries only the count precision ``1/n`` (so
   the own precision is the raw count ``n``). This is the anchor the whole prior-free pass leans on.
2. **DENSITY DECONVOLUTION** — a node's gDNA is peeled against a gDNA density prior via NegBinom
   (`density_deconv`); the curvature of that per-node ``λ``-factor is honest, count-derived composition
   evidence, registered here as ``τ_λ`` (`density_factor_precision`). The **intron factory** is its special
   case (the gDNA prior = the intergenic node distribution).
3. **STRAND DECONVOLUTION** — the strand Beta-Binomial is RANK-1 (it informs only ``p``), so its
   gDNA-level (``λ``) precision is the Schur-marginal (`variance_foundation_proposal.md`, approach E): a 1-DOF
   (single-strand) node has its tilt structurally locked, so the strand PINS ``f_g`` (``τ_λ`` gets the strand
   λ-term ``c·a²``); a 2-DOF (AMBIG) node's tilt is free, so the strand CANCELS out of ``f_g`` and contributes
   ZERO (it constrains only the tilt). `strand_evidence` returns the single-strand λ-term; `build_node_init`
   gates it to single-strand nodes. Identically zero on unstranded data (κ=½) by a derived noise-floor deadband.
4. **UNSOLVED** nodes default to **100 % gDNA at ZERO precision** — the honest "no information" state
   (``τ_λ = 0 ⇒ Var(log f) = ∞ ⇒ p = 0``), left for the sweep + population prior to resolve.

The precision arithmetic (sources → per-component ``Var(log f_c)`` → precision) is pure and unit-tested here.
Layer: imports only lower layers (`node_geometry`, `simplex_logodds`, `density_deconv`), never
`bp_solver`, so it sits cleanly beneath the solver that consumes :func:`build_node_init`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .density_deconv import density_factor_precision
from .node_chain import REGION, NodeChain
from .node_geometry import NodeGeometry, NodeStatics, node_global_geometry
from .simplex_logodds import _logodds_grid, _solve_nodes_logodds_all

__all__ = [
    "NodeInit",
    "own_composition_logvar",
    "own_precision",
    "strand_evidence",
    "build_node_init",
]

_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class NodeInit:
    """The per-node message-free self-solve (length ``n_nodes``) — the pass-0 relay's starting beliefs.

    ``f_*`` are the composition MODE (strand / intron / structural-default fractions). ``rho_*`` are the own
    per-component DENSITIES ``ρ_c = f_c·M_node/E_c`` (the message currency; a killed / count-free component is
    0). ``prec_*`` are the own per-component PRECISIONS from the four sources — the strand + intron-factory
    composition evidence combined with the Poisson count power (0 where a node is genuinely uninformed).
    ``struct_lock`` marks the structurally composition-certain (pure-gDNA) nodes; ``tau_lam`` is the combined
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
    """The message-free composition variance of a node's own belief, on the two log-fraction arms — the
    honest "how well does this node know its gDNA/RNA split?" (`variance_model_concepts.md`).

    Three states, nothing between (`docs/calibration/archive/emission_and_precision_derivation.md` §2):

    * **structural lock** → composition CERTAIN → variance ``0`` (an intergenic pure-gDNA node);
    * **real evidence ``τ_λ > 0``** → ``Var(log f_g) = (1−f_g)²/τ_λ`` and ``Var(log f_r) = f_g²/τ_λ`` (the
      Jacobians of ``log f_g`` / ``log f_r`` w.r.t. ``λ = logit f_g``);
    * **no evidence ``τ_λ = 0``** → composition UNSEEN → variance ``∞`` (⇒ own precision 0 — the unsolved
      default). ``∞`` is set DIRECTLY (never ``fr²·∞``, which would be ``0·∞ = nan`` at the window edge).

    Returns ``(v_log_fg, v_log_fr)``, same shape as ``f_g``."""
    fg = np.clip(np.asarray(f_g, np.float64), 0.0, 1.0)
    fr = 1.0 - fg
    tau = np.asarray(tau_lam, np.float64)
    lock = np.asarray(struct_lock, bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        v_fg = np.where(lock, 0.0, np.where(tau > _EPS, fr * fr / np.maximum(tau, _EPS), np.inf))
        v_fr = np.where(lock, 0.0, np.where(tau > _EPS, fg * fg / np.maximum(tau, _EPS), np.inf))
    return v_fg, v_fr


def own_precision(n, v_log, live):
    """The own-belief precision of one component: ``p = n / (n·Var(log f) + 1) = 1/(Var(log f) + 1/n)`` — the
    composition variance combined with the Poisson count power (`CALIBRATION_ARCHITECTURE.md` §1.2). It is 0
    when there is no count (``n = 0``), no evidence (``Var(log f) = ∞``), or the component is not ``live``
    (density ≤ 0) — a composition-vacuous source emits nothing, with no ``0·∞`` nan (the ``∞`` is masked to a
    finite value BEFORE the product, matching ``np.where``'s both-branch evaluation)."""
    n = np.asarray(n, np.float64)
    v = np.asarray(v_log, np.float64)
    ok = np.isfinite(v)
    v_fin = np.where(ok, v, 0.0)
    return np.where(
        (n > 0.0) & ok & np.asarray(live, bool),
        n / (n * np.maximum(v_fin, 0.0) + 1.0),
        0.0,
    )


# ── source 3: the strand composition evidence (I_strand) + the structural lock ─────────────────────────────


def strand_evidence(
    u_pos, u_neg, fg_loc, *, kappa, od_g, od_r, n_gdna_obs, n_rna_obs, is_region, locked
):
    """The reference-free strand composition evidence ``τ₀_λ`` (**I_strand**) + the structural-certainty mask,
    evaluated at the message-free local ``fg_loc``. Pure; no cross-node coupling.

    ``I_strand(λ) = N_eff·disc·[f_g(1−f_g)]² / (4 p(1−p))``, ``p = κ + f_g(½−κ)`` — the strand Fisher
    information, IDENTICALLY 0 at κ=½ (unstranded). The count enters as the OVERDISPERSED effective count
    ``N_eff = N/(1+(N−1)ω_r)`` (power saturates at ~1/ω, not the raw depth), and the discriminability
    ``disc = 4·max(0, (κ−½)² − σ²_d)`` carries the DERIVED noise floor
    ``σ²_d = ¼·(1/N_rna + ω_r) + ¼·(1/N_gdna + ω_g)`` — a κ within √σ²_d of ½ is not composition signal (the
    deadband that kills the unstranded phantom). ``1/N_gdna`` gates a gDNA-free library (N_gdna=0 ⇒ σ²_d→∞ ⇒
    disc=0).

    ``struct_lock`` (**I_struct**) — composition CERTAIN — is scoped to true intergenic REGION nodes, never a
    G1 BOUNDARY seam (TSS/TES): a seam is structurally gDNA but sits between RNA-carrying exons, so its
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


def build_node_init(
    chain: NodeChain,
    statics: NodeStatics,
    geometry: NodeGeometry,
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
) -> NodeInit:
    """The pass-0 per-node self-solve → :class:`NodeInit`. Runs the message-free strand deconvolution
    (`simplex_logodds`), compiles the strand + intron-factory composition evidence, and assembles each node's
    own per-component density + precision from the four sources (see the module docstring).

    The strand deconvolution reference (`fg_ref`/`fpos_ref`/`fneg_ref`) is the incoming ``belief`` — the
    count-zero-information variance freeze evaluates the composition variance near the truth, not at a flat ½.
    ``global_logprior`` (the anchored population gDNA prior, ``(m, K)``) and ``intron_prior`` (the intron
    factory ``λ``-factor, ``(m, K)``) enter ψ; ``intron_prior`` additionally seeds I_factory."""
    is_reg = np.asarray(chain.kind) == REGION
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)

    # ── source 3: the message-free strand deconvolution (1-DOF solves; AMBIG partial) ──
    dc = _solve_nodes_logodds_all(
        statics.u_pos,
        statics.u_neg,
        fp,
        fn,
        statics.mass_unspliced,
        statics.mass_spliced,
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
    vp_loc = np.asarray(dc.rna_pos_frac_var, np.float64)
    vn_loc = np.asarray(dc.rna_neg_frac_var, np.float64)

    # a node that does not deconvolve its own split (G1 sink / empty) keeps the signature-binary init.
    solvable = (fp | fn) & (np.asarray(statics.mass_unspliced, np.float64) > 0.0)
    locked = ~solvable
    fg_loc = np.where(locked, np.asarray(belief.f_g, np.float64), fg_loc)
    fp_loc = np.where(locked, np.asarray(belief.f_pos, np.float64), fp_loc)
    fn_loc = np.where(locked, np.asarray(belief.f_neg, np.float64), fn_loc)

    # ── sources 1 & 3: the composition evidence τ_λ (the Schur-marginal gDNA-level precision) + struct lock ──
    # `strand_evidence` returns the SINGLE-STRAND strand λ-Fisher I_strand = c·a² (the value at the locked tilt).
    i_strand, struct_lock = strand_evidence(
        statics.u_pos,
        statics.u_neg,
        fg_loc,
        kappa=kappa,
        od_g=od_g,
        od_r=od_r,
        n_gdna_obs=n_gdna_obs,
        n_rna_obs=n_rna_obs,
        is_region=is_reg,
        locked=locked,
    )
    # APPROACH E (docs/calibration/variance_foundation_proposal.md, verified). The strand Beta-Binomial is
    # RANK-1: it depends on (λ,θ) only through p = ½+(κ−½)(1−f_g)sinθ. So the honest MARGINAL gDNA-level
    # precision (the Schur complement of the 2×2 composition Fisher) is:
    #   * SINGLE-STRAND (1-DOF): θ is STRUCTURALLY locked ⇒ τ_λ gets the full strand λ-term c·a² (strand pins f_g);
    #   * AMBIG (2-DOF): θ is a FREE nuisance ⇒ the strand CANCELS out of f_g (Schur ⇒ 0) — it constrains only
    #     the tilt, never the gDNA level. Crediting c·a² to an AMBIG node is a (bounded) phantom precision on
    #     exactly the nodes calibration exists to resolve. Gate the strand λ-term to single-strand nodes.
    single_strand = np.asarray(fp, bool) ^ np.asarray(fn, bool)
    tau_lam = np.where(single_strand, i_strand, 0.0)
    lam_grid, _ = _logodds_grid(int(n_grid), float(logodds_window))
    tau_fac = density_factor_precision(intron_prior, lam_grid)  # I_density (NB curvature) on the λ axis
    if tau_fac is not None:
        tau_lam = tau_lam + tau_fac

    # ── the own per-component densities + precisions (node-level, both faces pooled) ──
    mass_global, eff_global = node_global_geometry(chain, geometry)
    n_node = np.where(
        is_reg,
        np.asarray(geometry.n_unspl_left, np.float64),
        np.asarray(geometry.n_unspl_left, np.float64)
        + np.asarray(geometry.n_unspl_right, np.float64),
    )
    eff_r = np.where(
        is_reg,
        np.asarray(geometry.eff_rna_left, np.float64),
        np.asarray(geometry.eff_rna_left, np.float64)
        + np.asarray(geometry.eff_rna_right, np.float64),
    )
    v_log_fg, v_log_fr = own_composition_logvar(fg_loc, tau_lam, struct_lock)

    # gDNA (source 1 measured / source 3 strand / source 4 default all flow through here):
    rho_g = np.where(
        (mass_global > _EPS) & (eff_global > _EPS),
        fg_loc * mass_global / np.maximum(eff_global, _EPS),
        0.0,
    )
    rho_g = np.maximum(rho_g, 0.0)
    prec_g = own_precision(n_node, v_log_fg, rho_g > _EPS)

    # RNA per strand — the density is live iff the strand is free, there is mass, and the local solve gave a
    # finite posterior variance (a node the strand could actually resolve on this axis):
    def _rna(frac_loc, var_loc, free_s):
        rho_raw = np.where(
            (mass_global > _EPS) & (eff_r > _EPS) & np.asarray(free_s, bool),
            np.asarray(frac_loc, np.float64) * mass_global / np.maximum(eff_r, _EPS),
            0.0,
        )
        live = (n_node > 0.0) & np.isfinite(np.asarray(var_loc, np.float64)) & (rho_raw > _EPS)
        rho = np.where(live, rho_raw, 0.0)
        prec = own_precision(n_node, v_log_fr, rho > _EPS)
        return rho, prec

    rho_pos, prec_pos = _rna(fp_loc, vp_loc, fp)
    rho_neg, prec_neg = _rna(fn_loc, vn_loc, fn)

    return NodeInit(
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
