"""The log-density 1-D/2-D per-node solver — the single production per-node solve driving
``bp_solver.node_sweep`` (the memory-prohibitive 2-simplex lattice it replaced is retired).

Design: ``docs/calibration/log_density_1d_solver_design.md``. The latent magnitude dof is the
gDNA-vs-RNA **log-odds** ``λ = logit(f_g) = log ρ_g − log ρ_rna`` (log-odds bounds the 5–6-decade ρ_g
range and resolves both ``f_g→0`` and ``f_g→1`` vertices, which the uniform linear lattice cannot). We
grid ``λ`` on a FIXED ``[−L, L]`` window (no node-adaptivity) and read out the linear fraction
``f_g = σ(λ)``. ``O(m·K)`` per node (vs the lattice's ``O(m·K²)`` 2-simplex), so it is genome-scale
tractable.

The ``ψ`` integrand sums the strand mixture, the sided spliced floor, the Jeffreys reference, the
count-space global, and the gDNA + per-strand RNA imputation messages, PLUS the change-of-variable
Jacobian ``log σ'(λ) = log(f_g(1−f_g))`` so the uniform-``λ`` Riemann sum is an unbiased posterior
integral. Single-strand nodes (exactly one of ``allow_pos`` / ``allow_neg``) are an exact 1-D solve
over ``λ``; AMBIG nodes (both set) marginalize the RNA tilt ``τ`` on a 2-D ``(λ, τ)`` grid
(``_solve_ambig_logodds``). ``_solve_nodes_logodds_all`` dispatches between the two.
"""

from __future__ import annotations

import numpy as np
from scipy.special import expit, log_expit, logsumexp

from .simplex import _mixture_strand_loglik
from .strand_deconv import NodeDeconv

__all__ = ["_logodds_grid", "_tilt_grid", "_single_strand_mask", "_ambig_mask",
           "_log_fg", "_log1m_fg", "_floor_log_density",
           "_local_loglik_logodds", "_solve_nodes_logodds", "_solve_ambig_logodds",
           "_solve_nodes_logodds_all"]

_EPS = 1.0e-9
_PRIOR_EPS = 1.0e-3  # the Jeffreys edge floor
_STRAND_PRIOR = 0.5  # Beta(½,½) Jeffreys strand prior
_DEFAULT_L = 10.0    # f_g ∈ [σ(−10), σ(10)] = [4.5e-5, 1−4.5e-5] — brackets Phase-0's vertex mass; a
#                      data-driven tuning knob (the (K,L) sweep), NOT a numerical ceiling: log_expit +
#                      the 1/E density floor keep the grid exact at any L (no underflow cap needed).
_AMBIG_BATCH = 8192  # memory-tiling batch for the AMBIG 2-D (λ,τ) cube — bounds the (B,K,K_t) materialized
#                      array at genome scale. NOT a model parameter: AMBIG nodes solve INDEPENDENTLY (no
#                      cross-node coupling), so any batch size yields bit-identical results — it only caps
#                      peak memory (the full subset at once is ~O(m·K²), same as the retired lattice).


def _log_fg(lam):
    """``log f_g = log σ(λ)`` — computed via ``scipy.special.log_expit`` (stable: it never forms
    ``1−σ(λ)``, so it stays exact in the depleted tail where the naive ``log(clip(σ(λ), ε))`` underflows)."""
    return log_expit(np.asarray(lam, dtype=np.float64))


def _log1m_fg(lam):
    """``log(1 − f_g) = log σ(−λ) = log f_rna`` — exact via ``log_expit(−λ)``."""
    return log_expit(-np.asarray(lam, dtype=np.float64))


def _floor_log_density(rho, eff):
    """``log(max(ρ, 1/E))`` — the consistent density floor (D5): a density from a finite count can never
    be exactly 0, the smallest being one ``_MSG_PSEUDOCOUNT=1`` pseudo-fragment over the eff-len ``E``.
    Used identically in the message mean ``μ_log``, ``pois_log``, and the σ²_bio offset so they agree. No
    new constant (ties to the existing pseudocount)."""
    r = np.asarray(rho, dtype=np.float64)
    e = np.maximum(np.asarray(eff, dtype=np.float64), _EPS)
    return np.log(np.maximum(r, 1.0 / e))


def _logodds_grid(n_grid: int, L: float = _DEFAULT_L):
    """The fixed log-odds lattice: ``λ`` uniform on ``[−L, L]`` (``K = n_grid`` points, ascending) and
    the matching ``f_g = σ(λ)`` (also ascending). Returns ``(lam, fg)``, each length ``K``."""
    lam = np.linspace(-float(L), float(L), int(n_grid))
    return lam, expit(lam)


def _tilt_grid(n_tilt: int) -> np.ndarray:
    """The RNA-internal tilt grid ``τ ∈ [−1, 1]`` (``K_t`` points). ``τ = ±1`` ⇒ all RNA on one strand
    (the single-strand edges); ``τ = 0`` ⇒ balanced. Only AMBIG nodes integrate over it."""
    return np.linspace(-1.0, 1.0, int(n_tilt))


def _single_strand_mask(allow_pos, allow_neg) -> np.ndarray:
    """The nodes the 1-D (Phase-1) solver is valid for: exactly one strand live (tilt determined)."""
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    return ap ^ an


def _ambig_mask(allow_pos, allow_neg) -> np.ndarray:
    """AMBIG nodes (both strands live) — the Phase-2 2-D ``(λ, τ)`` path."""
    return np.asarray(allow_pos, bool) & np.asarray(allow_neg, bool)


def _local_loglik_logodds(
    u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa, od_g, od_r, lam, fg,
    strand_obs=None, global_logprior=None,
    gdna_imp_mode=None, gdna_imp_prec=None, rna_imp_mode=None, rna_imp_prec=None,
):
    """ψ over the log-odds grid for single-strand nodes (strand mixture, sided spliced floor, Jeffreys,
    global, imputation), evaluated at ``f_g = σ(λ)`` with the live strand carrying ``f_active = 1 − f_g``,
    PLUS the ``log σ'(λ)`` Jacobian. Returns ``(m, K)``.

    ``global_logprior`` must already be evaluated on THIS ``fg`` grid → ``(m, K)`` (the caller builds
    it with ``bp_solver._global_logprior`` passing ``fg`` as the ``fgg`` arg)."""
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    pos_live = (ap & ~an)[:, None]  # (m,1)
    neg_live = (an & ~ap)[:, None]
    fg2 = fg[None, :]               # (1,K)
    f_act = 1.0 - fg2
    f_pos = np.where(pos_live, f_act, 0.0)  # (m,K)
    f_neg = np.where(neg_live, f_act, 0.0)
    n = (u_pos + u_neg)[:, None]
    # ── strand mixture (identical to the lattice term; broadcasts (m,1)×(m,K)→(m,K)) ──
    psi = _mixture_strand_loglik(u_pos[:, None], n, fg2, f_pos, f_neg, kappa, od_g, od_r)
    # ── sided spliced lower bound (identical) ──
    inv_n = np.where(n > 0.0, 1.0 / np.maximum(n, _EPS), 0.0)
    short_p = np.maximum(spliced_pos[:, None] * inv_n - f_pos, 0.0)
    short_n = np.maximum(spliced_neg[:, None] * inv_n - f_neg, 0.0)
    psi = psi - 0.5 * spliced_pos[:, None] * short_p**2 - 0.5 * spliced_neg[:, None] * short_n**2
    # ── Jeffreys Beta(½,½) reference at strand-observable nodes (identical) ──
    if strand_obs is not None:
        jeff = (_STRAND_PRIOR - 1.0) * (
            np.log(np.clip(fg2, _PRIOR_EPS, 1.0))
            + np.log(np.clip(f_pos, _PRIOR_EPS, 1.0))
            + np.log(np.clip(f_neg, _PRIOR_EPS, 1.0))
        )
        psi = psi + np.asarray(strand_obs, np.float64)[:, None] * jeff
    # ── count-space global (identical; pre-evaluated on this fg grid by the caller) ──
    if global_logprior is not None:
        psi = psi + np.asarray(global_logprior, np.float64)
    # ── imputation messages: LOG-FRACTION Gaussians (the overhaul). The mode is a log-FRACTION target
    #    (``log`` of the imputed fraction, built in ``_scan``); evaluated against ``log f_c(λ)``. No clip —
    #    an off-grid target (source denser than the dst can hold) is a bounded monotone pull toward the
    #    edge, governed by precision (D-plan P6, verify-don't-clip). ──
    log_fg = _log_fg(lam)[None, :]      # log f_g = log σ(λ)
    log_fact = _log1m_fg(lam)[None, :]  # log(1−f_g) = log f_active (the single live strand)
    if gdna_imp_mode is not None and gdna_imp_prec is not None:
        m_ = np.asarray(gdna_imp_mode, np.float64)[:, None]
        p_ = np.asarray(gdna_imp_prec, np.float64)[:, None]
        psi = psi - 0.5 * p_ * (log_fg - m_) ** 2
    if rna_imp_mode is not None and rna_imp_prec is not None:
        # single-strand: the live strand carries f_active = 1−f_g; the per-strand precision gates which
        # message applies (the dead strand's prec is 0 → no-op). Both evaluate against log f_active.
        for ms, ps in ((rna_imp_mode[0], rna_imp_prec[0]), (rna_imp_mode[1], rna_imp_prec[1])):
            psi = psi - 0.5 * np.asarray(ps, np.float64)[:, None] * (
                log_fact - np.asarray(ms, np.float64)[:, None]
            ) ** 2
    # ── change of variable: uniform-λ Riemann sum → uniform-f integral (so the median matches the
    #    linear lattice's uniform-f measure). log σ'(λ) = log σ(λ) + log(1−σ(λ)) = log f_g + log(1−f_g),
    #    exact via log_expit (D6). This is the ONE reference Jacobian (D3); the log-density empirical
    #    priors are added bare (no per-term Jacobian). ──
    psi = psi + (_log_fg(lam) + _log1m_fg(lam))[None, :]
    return psi, f_pos, f_neg


def _solve_nodes_logodds(
    u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, strand_obs,
    mass_unspl, mass_spliced, *, kappa, od_g, od_r, n_grid, L: float = _DEFAULT_L,
    global_logprior=None, gdna_imp_mode=None, gdna_imp_prec=None,
    rna_imp_mode=None, rna_imp_prec=None,
) -> NodeDeconv:
    """The log-odds 1-D per-node solve for SINGLE-STRAND nodes.

    ``f_g`` = posterior median over the ``λ`` grid (= median of ``f_g = σ(λ)``, monotone so the median
    is preserved); ``f_pos``/``f_neg`` = posterior MEANS (the current-state fractions). ``*_frac_var`` =
    the precision state ``Var(log f_c)`` (the LOG-fraction variance, the message currency — D2), NOT the
    fraction variance. Zero-mass nodes report 0; the dead strand is locked-certain (var 0). AMBIG nodes
    (both strands live) are out of contract — mask them out."""
    lam, fg = _logodds_grid(int(n_grid), L)
    psi, f_pos_g, f_neg_g = _local_loglik_logodds(
        u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa, od_g, od_r, lam, fg,
        strand_obs=strand_obs, global_logprior=global_logprior,
        gdna_imp_mode=gdna_imp_mode, gdna_imp_prec=gdna_imp_prec,
        rna_imp_mode=rna_imp_mode, rna_imp_prec=rna_imp_prec,
    )
    post = np.exp(psi - logsumexp(psi, axis=1, keepdims=True))  # (m,K)
    # f_g posterior median (fg ascending ⇒ cumulative CDF directly)
    cw = np.cumsum(post, axis=1)
    idx = np.clip((cw < 0.5).sum(axis=1), 0, fg.size - 1)
    f_g = fg[idx]
    # composition: f_g median + f_pos/f_neg posterior MEANS (the current-state fractions → node_densities).
    f_pos = np.sum(post * f_pos_g, axis=1)
    f_neg = np.sum(post * f_neg_g, axis=1)
    # precision state = Var(log f_c) (the message currency, D2), moment-matched on the grid. log f_g = log σ(λ);
    # the single live RNA strand carries f_active = 1−f_g → log f_active = log σ(−λ). The dead strand is
    # locked-certain (f=0) → var 0 (and emits nothing — the _scan gate). Capping is AUTOMATIC: the send
    # prec_log = 1/(var+σ²_bio+pois) ≤ 1/(σ²_bio+pois) since var ≥ 0 (D2a), so a window-truncated var only
    # approaches, never exceeds, the noise-model ceiling.
    Lg = _log_fg(lam)
    mLg = post @ Lg
    var_g = np.maximum(post @ (Lg * Lg) - mLg * mLg, 0.0)
    La = _log1m_fg(lam)
    mLa = post @ La
    var_act = np.maximum(post @ (La * La) - mLa * mLa, 0.0)
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    var_pos = np.where(ap & ~an, var_act, 0.0)
    var_neg = np.where(an & ~ap, var_act, 0.0)
    active = (u_pos + u_neg) > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_pos = np.where(active, np.clip(f_pos, 0.0, 1.0), 0.0)
    f_neg = np.where(active, np.clip(f_neg, 0.0, 1.0), 0.0)
    var_g = np.where(active, var_g, 0.0)
    var_pos = np.where(active, var_pos, 0.0)
    var_neg = np.where(active, var_neg, 0.0)
    return NodeDeconv(
        gdna_mass=f_g * mass_unspl, rna_mass=(1.0 - f_g) * mass_unspl + mass_spliced,
        gdna_frac=f_g, rna_pos_frac=f_pos, rna_neg_frac=f_neg,
        gdna_frac_var=var_g, rna_pos_frac_var=var_pos, rna_neg_frac_var=var_neg,
    )


def _solve_ambig_logodds(
    u_pos, u_neg, spliced_pos, spliced_neg, mass_unspl, mass_spliced, *,
    kappa, od_g, od_r, n_grid, L: float = _DEFAULT_L, n_tilt: int | None = None,
    global_logprior=None, gdna_imp_mode=None, gdna_imp_prec=None,
    rna_imp_mode=None, rna_imp_prec=None,
) -> NodeDeconv:
    """Phase 2 — the 2-D ``(λ, τ)`` solve for AMBIG nodes (both strands live; NO Jeffreys, matching the
    lattice). Grids the gDNA-vs-RNA log-odds ``λ`` (outer, ``K = n_grid``) and the RNA-internal tilt
    ``τ`` (inner, ``K_t = n_tilt`` or ``n_grid``), evaluates the SAME ψ terms as the lattice on the
    ``(m, K, K_t)`` cube, **marginalizes τ** (``logsumexp``), then adds the τ-INDEPENDENT 2-D
    change-of-variable Jacobian ``log|J| = log f_g + 2·log(1−f_g)`` (uniform ``(f_pos, f_neg)`` →
    ``(λ, τ)``; derivation in ``log_density_1d_solver_design.md`` §5.2). ``f_g`` = posterior median over
    the τ-marginal λ-posterior; ``f_pos``/``f_neg`` = means over the full 2-D posterior.

    ``global_logprior`` is ``(m, K)`` evaluated on the σ(λ) grid (broadcast over τ). The cube is only
    materialized for the AMBIG subset (the caller masks); ``K·K_t`` is the per-node cost."""
    lam, fg = _logodds_grid(int(n_grid), L)            # (K,)
    K = fg.size
    Kt = int(n_tilt) if n_tilt else int(n_grid)
    tau = _tilt_grid(Kt)                               # (Kt,)
    f_act = (1.0 - fg)[:, None]                         # (K,1)
    f_pos_kt = f_act * (1.0 + tau[None, :]) / 2.0       # (K,Kt)
    f_neg_kt = f_act * (1.0 - tau[None, :]) / 2.0       # (K,Kt)
    n = u_pos + u_neg
    # ── strand mixture over the cube: (m,1,1)×(1,K,1)×(1,K,Kt) → (m,K,Kt) ──
    psi = _mixture_strand_loglik(
        u_pos[:, None, None], n[:, None, None],
        fg[None, :, None], f_pos_kt[None, :, :], f_neg_kt[None, :, :], kappa, od_g, od_r,
    )
    # ── sided spliced floor ──
    inv_n = np.where(n > 0.0, 1.0 / np.maximum(n, _EPS), 0.0)[:, None, None]
    short_p = np.maximum(spliced_pos[:, None, None] * inv_n - f_pos_kt[None, :, :], 0.0)
    short_n = np.maximum(spliced_neg[:, None, None] * inv_n - f_neg_kt[None, :, :], 0.0)
    psi = psi - 0.5 * spliced_pos[:, None, None] * short_p**2 - 0.5 * spliced_neg[:, None, None] * short_n**2
    # ── LOG-fraction grids (the overhaul): log f_g (τ-independent) + log f_pos/f_neg over the cube,
    #    floored at one pseudo-fragment 1/(n+1) (D5: the τ=±1 edges have f_s=0 → log(0); the count floor
    #    keeps it finite + consistent with pois_log). ──
    log_fg_grid = _log_fg(lam)                                   # (K,) = log f_g
    frac_floor = (1.0 / (n + 1.0))[:, None, None]                # (m,1,1)
    log_fpos = np.log(np.maximum(f_pos_kt[None, :, :], frac_floor))  # (m,K,Kt)
    log_fneg = np.log(np.maximum(f_neg_kt[None, :, :], frac_floor))
    # ── log-density global on log f_g (τ-independent; pre-evaluated on the σ(λ) grid by the caller) ──
    if global_logprior is not None:
        psi = psi + np.asarray(global_logprior, np.float64)[:, :, None]
    # ── gDNA LOG-fraction message on log f_g (τ-independent) ──
    if gdna_imp_mode is not None and gdna_imp_prec is not None:
        mo = np.asarray(gdna_imp_mode, np.float64)[:, None, None]
        pr = np.asarray(gdna_imp_prec, np.float64)[:, None, None]
        psi = psi - 0.5 * pr * (log_fg_grid[None, :, None] - mo) ** 2
    # ── per-strand RNA LOG-fraction messages on log f_pos/log f_neg (τ-dependent — inside the cube) ──
    if rna_imp_mode is not None and rna_imp_prec is not None:
        for log_f, ms, ps in (
            (log_fpos, rna_imp_mode[0], rna_imp_prec[0]),
            (log_fneg, rna_imp_mode[1], rna_imp_prec[1]),
        ):
            psi = psi - 0.5 * np.asarray(ps, np.float64)[:, None, None] * (
                log_f - np.asarray(ms, np.float64)[:, None, None]
            ) ** 2
    # ── Reference measure (τ-independent → add once on the λ axis). NEUTRAL on the gDNA fraction:
    #    uniform in (f_g, τ) ⇒ Jacobian log σ'(λ) = log f_g + log(1−f_g) — the SAME measure the
    #    single-strand solver uses. The old uniform-(f_pos, f_neg) simplex measure (log f_g + 2 log(1−f_g))
    #    put an implicit ANTI-gDNA prior on AMBIG nodes (it peaks at f_g = ⅓ and prefers balanced ±RNA), so
    #    a gDNA-dominant AMBIG node — a near-50/50 count — was pulled to f_g ≈ 0.46 no matter what the strand
    #    tilt or the gDNA prior said. Neutral lets the STRAND resolve it: a strong tilt ⇒ RNA (unchanged),
    #    a balanced count ⇒ parsimoniously gDNA (the τ-marginal favours high f_g, where more τ fit the
    #    balance) with the gDNA prior setting the level. (Derivation: log_density_1d_solver_design.md §5.2.) ──
    log_jac = _log_fg(lam) + _log1m_fg(lam)  # (K,)
    psi_full = psi + log_jac[None, :, None]                 # (m,K,Kt) — the full 2-D log-posterior
    # τ-marginal λ-posterior (m,K)
    psi_lam = logsumexp(psi_full, axis=2)
    post_lam = np.exp(psi_lam - logsumexp(psi_lam, axis=1, keepdims=True))
    cw = np.cumsum(post_lam, axis=1)
    idx = np.clip((cw < 0.5).sum(axis=1), 0, K - 1)
    f_g = fg[idx]
    # precision state = Var(log f_g) over the τ-marginal λ-posterior (D2).
    mLg = post_lam @ log_fg_grid
    var_g = np.maximum(post_lam @ (log_fg_grid * log_fg_grid) - mLg * mLg, 0.0)
    # f_pos/f_neg MEANS + Var(log f_pos/neg) over the FULL 2-D posterior.
    flat = psi_full.reshape(psi_full.shape[0], -1)
    post2d = np.exp(flat - logsumexp(flat, axis=1, keepdims=True)).reshape(psi_full.shape)  # (m,K,Kt)
    fp_grid = f_pos_kt[None, :, :]
    fn_grid = f_neg_kt[None, :, :]
    f_pos = np.sum(post2d * fp_grid, axis=(1, 2))
    f_neg = np.sum(post2d * fn_grid, axis=(1, 2))
    mLp = np.sum(post2d * log_fpos, axis=(1, 2))
    mLn = np.sum(post2d * log_fneg, axis=(1, 2))
    var_pos = np.maximum(np.sum(post2d * log_fpos * log_fpos, axis=(1, 2)) - mLp * mLp, 0.0)
    var_neg = np.maximum(np.sum(post2d * log_fneg * log_fneg, axis=(1, 2)) - mLn * mLn, 0.0)
    active = n > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_pos = np.where(active, np.clip(f_pos, 0.0, 1.0), 0.0)
    f_neg = np.where(active, np.clip(f_neg, 0.0, 1.0), 0.0)
    var_g = np.where(active, var_g, 0.0)
    var_pos = np.where(active, var_pos, 0.0)
    var_neg = np.where(active, var_neg, 0.0)
    return NodeDeconv(
        gdna_mass=f_g * mass_unspl, rna_mass=(1.0 - f_g) * mass_unspl + mass_spliced,
        gdna_frac=f_g, rna_pos_frac=f_pos, rna_neg_frac=f_neg,
        gdna_frac_var=var_g, rna_pos_frac_var=var_pos, rna_neg_frac_var=var_neg,
    )


def _solve_nodes_logodds_all(
    u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, strand_obs,
    mass_unspl, mass_spliced, *, kappa, od_g, od_r, n_grid, L: float = _DEFAULT_L,
    n_tilt: int | None = None, global_logprior=None,
    gdna_imp_mode=None, gdna_imp_prec=None, rna_imp_mode=None, rna_imp_prec=None,
) -> NodeDeconv:
    """The full per-node log-odds dispatcher (Phase 3 #1): routes single-strand nodes to the 1-D
    ``λ`` solve (:func:`_solve_nodes_logodds`) and AMBIG nodes to the 2-D ``(λ, τ)`` solve
    (:func:`_solve_ambig_logodds`), scattering both into full-length arrays. G1 / zero-mass nodes
    report 0 (``node_sweep`` keeps their signature-binary init via the ``solvable`` write-back). A
    drop-in for the lattice ``_local_loglik``+``_node_marginals`` pair: same ψ terms (still
    linear-fraction messages / global — the log-density migration is Phase 3 #2), only the grid differs.

    All array inputs are full length ``m``; ``global_logprior`` is ``(m, K)`` on the σ(λ) grid;
    ``gdna_imp_*`` are ``(m,)``; ``rna_imp_*`` are 2-tuples of ``(m,)``. Each is sub-indexed per class."""
    m = int(np.asarray(u_pos).shape[0])
    out = {k: np.zeros(m, dtype=np.float64) for k in
           ("fg", "fp", "fn", "vg", "vp", "vn", "gmass", "rmass")}
    # Skip EMPTY nodes — no per-strand counts AND no unspliced/spliced mass. Both per-class solvers zero
    # every output for an inactive node (gdna/rna_mass = f_g·M = (1−f_g)·M + S = 0 when all are 0), so an
    # empty node's solve is identical to the zero-initialized `out` — skipping is BIT-IDENTICAL. At genome
    # scale most region/boundary nodes carry no fragments (unexpressed genes, intergenic deserts), so this
    # is the dominant cost saver, not a slice artifact. (A spliced-only node has signal ⇒ still solved.)
    signal = (np.asarray(u_pos, np.float64) + np.asarray(u_neg, np.float64)
              + np.asarray(mass_unspl, np.float64) + np.asarray(mass_spliced, np.float64)) > 0.0
    ss = _single_strand_mask(allow_pos, allow_neg) & signal
    amb = _ambig_mask(allow_pos, allow_neg) & signal

    def _s(a, msk):
        return None if a is None else np.asarray(a)[msk]

    def _sp(pair, msk):
        return None if pair is None else (np.asarray(pair[0])[msk], np.asarray(pair[1])[msk])

    def _scatter(msk, dc):
        out["fg"][msk] = dc.gdna_frac
        out["fp"][msk] = dc.rna_pos_frac
        out["fn"][msk] = dc.rna_neg_frac
        out["vg"][msk] = dc.gdna_frac_var
        out["vp"][msk] = dc.rna_pos_frac_var
        out["vn"][msk] = dc.rna_neg_frac_var
        out["gmass"][msk] = dc.gdna_mass
        out["rmass"][msk] = dc.rna_mass

    if bool(ss.any()):
        _scatter(ss, _solve_nodes_logodds(
            u_pos[ss], u_neg[ss], spliced_pos[ss], spliced_neg[ss], allow_pos[ss], allow_neg[ss],
            _s(strand_obs, ss), mass_unspl[ss], mass_spliced[ss], kappa=kappa, od_g=od_g, od_r=od_r,
            n_grid=n_grid, L=L, global_logprior=_s(global_logprior, ss),
            gdna_imp_mode=_s(gdna_imp_mode, ss), gdna_imp_prec=_s(gdna_imp_prec, ss),
            rna_imp_mode=_sp(rna_imp_mode, ss), rna_imp_prec=_sp(rna_imp_prec, ss)))
    if bool(amb.any()):
        # The 2-D (λ,τ) cube is (B,K,K_t); materialized for ALL ambig nodes at once it is ~O(m·K²) (the
        # memory the lattice OOM'd on). AMBIG nodes solve independently, so tile the subset into row-batches
        # — bit-identical results, peak memory bounded to one (≤_AMBIG_BATCH, K, K_t) cube.
        amb_idx = np.where(amb)[0]
        for s0 in range(0, amb_idx.size, _AMBIG_BATCH):
            bidx = amb_idx[s0:s0 + _AMBIG_BATCH]
            _scatter(bidx, _solve_ambig_logodds(
                u_pos[bidx], u_neg[bidx], spliced_pos[bidx], spliced_neg[bidx], mass_unspl[bidx],
                mass_spliced[bidx], kappa=kappa, od_g=od_g, od_r=od_r, n_grid=n_grid, L=L, n_tilt=n_tilt,
                global_logprior=_s(global_logprior, bidx),
                gdna_imp_mode=_s(gdna_imp_mode, bidx), gdna_imp_prec=_s(gdna_imp_prec, bidx),
                rna_imp_mode=_sp(rna_imp_mode, bidx), rna_imp_prec=_sp(rna_imp_prec, bidx)))
    return NodeDeconv(
        gdna_mass=out["gmass"], rna_mass=out["rmass"], gdna_frac=out["fg"],
        rna_pos_frac=out["fp"], rna_neg_frac=out["fn"], gdna_frac_var=out["vg"],
        rna_pos_frac_var=out["vp"], rna_neg_frac_var=out["vn"],
    )
