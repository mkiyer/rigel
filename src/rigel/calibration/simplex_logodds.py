"""The log-density 1-D/2-D per-node solver — the single production per-node solve driving
``bp_solver.node_sweep`` (the memory-prohibitive 2-simplex lattice it replaced is retired).

Design: ``docs/calibration/archive/log_density_1d_solver_design.md``. The latent magnitude dof is the
gDNA-vs-RNA **log-odds** ``λ = logit(f_g) = log ρ_g − log ρ_rna`` (log-odds bounds the 5–6-decade ρ_g
range and resolves both ``f_g→0`` and ``f_g→1`` vertices, which the uniform linear lattice cannot). We
grid ``λ`` on a FIXED ``[−L, L]`` window (no node-adaptivity) and read out the linear fraction
``f_g = σ(λ)``. ``O(m·K)`` per node (vs the lattice's ``O(m·K²)`` 2-simplex), so it is genome-scale
tractable.

The ``ψ`` integrand is ``strand + (gDNA arm) + (RNA arm) + the imputation messages``, where each **arm** is
that component group's fitted log-rate prior when we have one, else the **Jeffreys reference**
``+½·log f`` (``_JEFFREYS_REF``). Derivation, review, and the resolved design:
``docs/calibration/archive/reference_prior_derivation.md`` (§10 is authoritative).

Three facts that determine this file's shape:

1. **Omitting a component's term is not "no prior"** — the grid's own measure supplies one. A bare uniform-λ
   grid IS Haldane per component ⇒ Beta(0,0) on the composition: improper at BOTH vertices, a vertex
   amplifier. There is no third option. (This is what shipped before; it hid while ψ stayed symmetric.)
2. **The composition is a TWO-GROUP split on the λ axis** — gDNA vs RNA-**total** — which is what calibration
   models. The per-strand tilt is a NUISANCE parameter. On the two-group axis the measure residual is
   **exactly zero**: each ``logP`` is a density in LOG-rate, so its linear-rate conversion ``−log ρ`` cancels
   ``log σ'(λ)`` exactly, once per group. **No Jacobian is written.**
3. **The tilt is gridded as ``θ = arcsin(τ)``, NOT ``τ``.** The Berger–Bernardo reference prior for this model
   (``f_g`` of interest, ``τ`` nuisance — the two are information-ORTHOGONAL, ``I_{f_g,τ} = 0`` exactly) has a
   ``(1−τ²)^{−½}`` tilt conditional. Under ``θ = arcsin(τ)`` the Jacobian ``|dτ/dθ| = cos θ = (1−τ²)^{½}``
   cancels it **identically**, so the tilt term is **exactly 0** and the reference collapses to ONE expression
   for both node classes: ``ψ_ref = ½·log f_g + ½·log(1−f_g)``. No class branch, no endpoint singularity, no
   quadrature weights. θ is to the tilt what λ is to ``f_g``: the coordinate the geometry asks for.
   *(This vanishing is a property of the BB reference specifically — a Dirichlet(½,¼,¼) reference would leave
   a residual ``−¼·log(1−τ²)``.)*

There is NO spliced term: ``mass_spliced`` is consumed only by the returned ``rna_mass``, never by ψ. That is
correct — at a junction mature RNA *splices*, so the unspliced crossing mass is gDNA + nascent, a channel
genuinely disjoint from the (directly observed, already-pure-RNA) spliced mass.

Single-strand nodes (exactly one of ``allow_pos`` / ``allow_neg``) are an exact 1-D solve over ``λ``; AMBIG
nodes (both set) marginalize the tilt on a 2-D ``(λ, θ)`` grid (``_solve_ambig_logodds``).
``_solve_nodes_logodds_all`` dispatches between the two. Structurally RNA-free nodes (neither strand live —
intergenic / TSS / TES) have no composition DOF and never reach either solver: ``bp_solver.node_sweep`` gates
them out via ``solvable``, so no reference is applied to a node whose composition is known structurally.
"""

from __future__ import annotations

import numpy as np
from scipy.special import expit, log_expit

from .simplex import _mixture_strand_loglik
from .strand_deconv import NodeDeconv

# Public surface consumed by bp_solver / node_geometry. The remaining private helpers stay importable
# for tests but are not part of the module's external API.
__all__ = ["_logodds_grid", "_solve_nodes_logodds_all"]

_EPS = 1.0e-9

# The reference exponent for an UNFITTED component group, as a density in LOG-rate.
#
# DERIVED, not tuned, by two agreeing routes: (a) Jeffreys for a Poisson rate — `g ~ Poisson(ρE)` ⇒
# `I(ρ) ∝ 1/ρ` ⇒ `p(ρ) ∝ ρ^(−½)`, which as a LOG-rate density is `ρ^(+½)` ⇒ `+½·log f`; (b) the
# Berger–Bernardo reference prior for the composition with `f_g` of interest and the tilt as nuisance, whose
# `f_g` marginal is Beta(½,½) — the SAME `+½·log f_g + ½·log(1−f_g)`.
#
# Its ONLY job is to make ψ PROPER (Beta(½,½) integrates; Beta(0,0) does not). It carries no information, so
# it must be REPLACED by `logP_g`/`logP_r` as those are fitted — never added on top, which double-counts that
# arm and is what broke the first attempt at this.
#
# ⚠ A DECLARED CHOICE, not forced by the likelihood: the observed-data Fisher information for f_g is
# `∝ n(½−κ)²` = EXACTLY 0 on unstranded libraries, where the strand term is bit-flat and the posterior simply
# IS this reference. Licensed as the "structural Jeffreys" prior (docs/calibration/archive/reference_prior_derivation.md §10.6); §10.5
# records the known cost — it forbids the simplex vertices, where some truth genuinely lives.
_JEFFREYS_REF = 0.5

# f_g ∈ [σ(−10), σ(10)] = [4.5e-5, 1−4.5e-5]. A pure STATE-SPACE bracket: the widest f_g the grid can
# represent, NOT an accuracy knob — but that is a PROPERTY OF A PROPER ψ, not of this constant. It holds
# because both arms are now always written (`_JEFFREYS_REF`): under Beta(½,½) ~0.9% of the reference's mass
# lies outside L=10, and the answer is L-invariant. An improper ψ (either arm omitted) has plateau mass
# growing linearly in L, and then L silently sets the prior strength — which is what the `+0.5·λ` ramp was.
# **L-invariance is the acceptance test for this file.** NB: production does not read this default —
# `bp_solver.node_sweep` threads `logodds_window` (=10.0) explicitly.
_DEFAULT_L = 10.0

# Memory-tiling batch for the AMBIG 2-D (λ,θ) cube — bounds the (B,K,K_t) materialized array at genome scale.
# NOT a model parameter: AMBIG nodes solve INDEPENDENTLY (no cross-node coupling), so any batch size yields
# bit-identical results — it only caps peak memory (the full subset at once is ~O(m·K²), like the old lattice).
_AMBIG_BATCH = 8192


def _lse(a, axis, keepdims=False):
    """Lean numpy log-sum-exp — a drop-in for ``scipy.special._lse(a, axis, keepdims)`` without the
    scipy wrapper overhead (arg validation, ``b``/``return_sign`` handling), which the profiler flagged as
    ~9 s of pure per-call cost across the AMBIG solve. Same max-shift stabilisation; the ``m→0`` guard makes
    an all-``-inf`` slice give ``log(0) = -inf`` and avoids ``(-inf)-(-inf) = nan``."""
    m = np.max(a, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    with np.errstate(
        divide="ignore"
    ):  # all-(-inf) slice ⇒ log(0) = -inf (correct); suppress the warning
        # Accumulate the sum in float64 even for a float32 cube (a 60-/3600-wide reduction loses too much in
        # f32), then return the input dtype — keeps the f64 path byte-identical, gives the f32 cube an
        # accurate normalizer.
        s = np.sum(np.exp(a - m), axis=axis, keepdims=True, dtype=np.float64)
        r = (m + np.log(s)).astype(m.dtype, copy=False)
    return r if keepdims else np.squeeze(r, axis=axis)


def _log_fg(lam):
    """``log f_g = log σ(λ)`` — computed via ``scipy.special.log_expit`` (stable: it never forms
    ``1−σ(λ)``, so it stays exact in the depleted tail where the naive ``log(clip(σ(λ), ε))`` underflows)."""
    return log_expit(np.asarray(lam, dtype=np.float64))


def _log1m_fg(lam):
    """``log(1 − f_g) = log σ(−λ) = log f_rna`` — exact via ``log_expit(−λ)``."""
    return log_expit(-np.asarray(lam, dtype=np.float64))


def _logodds_grid(n_grid: int, L: float = _DEFAULT_L):
    """The fixed log-odds lattice: ``λ`` uniform on ``[−L, L]`` (``K = n_grid`` points, ascending) and
    the matching ``f_g = σ(λ)`` (also ascending). Returns ``(lam, fg)``, each length ``K``."""
    lam = np.linspace(-float(L), float(L), int(n_grid))
    return lam, expit(lam)


def _posterior_median_fg(post, fg):
    """Per-node point estimate of ``f_g``: the grid MEDIAN — ``fg`` at the grid point where the CDF first
    reaches 0.5. Transform-invariant and robust to the SKEW of the f_g posterior; the log-odds quantization
    (Δf_g≈0.085 at n_grid=60) is de-quantized by the finer single-strand grid (``sweep_n_grid_single_strand``),
    NOT by a different estimator (a sub-grid mode would silently under-call skewed/vertex-near posteriors).
    ``post``: (m,K) normalized posterior; ``fg``=σ(λ). Returns (m,)."""
    K = fg.shape[0]
    cw = np.cumsum(post, axis=1)
    idx = np.clip((cw < 0.5).sum(axis=1), 0, K - 1)
    return fg[idx]


def _regrid_global(glp, n_from, n_to, L):
    """Interpolate a ``(m, n_from)`` global-logprior (evaluated on the σ(λ) grid at ``n_from``) onto the
    ``n_to`` σ(λ) grid (Fix 3). The single-strand solve runs on a finer grid than the AMBIG cube, so the
    coarse-grid global prior is regridded to feed it; the global is smooth in ``f_g``, so linear
    interpolation is exact to interpolation accuracy. ``None`` / equal grids ⇒ passthrough (bit-identical)."""
    if glp is None or int(n_from) == int(n_to):
        return glp
    _, fc = _logodds_grid(int(n_from), L)
    _, ff = _logodds_grid(int(n_to), L)
    j = np.clip(np.searchsorted(fc, ff), 1, int(n_from) - 1)
    x0, x1 = fc[j - 1], fc[j]
    t = np.clip((ff - x0) / np.maximum(x1 - x0, _EPS), 0.0, 1.0)  # (n_to,)
    g = np.asarray(glp, np.float64)
    return g[:, j - 1] + t[None, :] * (g[:, j] - g[:, j - 1])


def _gdna_arm(lam, global_logprior):
    """The gDNA arm of ψ over the λ grid → broadcastable to ``(m, K)``.

    ``_JEFFREYS_REF`` reference ``+½·log f_g`` → ``(1, K)`` when no fitted prior. When a prior IS fitted, the
    fitted ``logP_g(log ρ_g)`` (pre-evaluated on THIS ``f_g`` grid → ``(m, K)``) REPLACES the reference — it
    carries no information and must not be double-counted.

    ``None`` means "not fitted", **not** "no term": writing nothing would silently select Haldane at this arm
    and leave ψ improper at ``f_g → 0``."""
    ref = _JEFFREYS_REF * _log_fg(lam)[None, :]
    if global_logprior is None:
        return ref
    return np.asarray(global_logprior, np.float64)


def _rna_arm(lam, rna_logprior=None):
    """The RNA-**total** arm of ψ over the λ grid → broadcastable to ``(m, K)``.

    ``logP_r(log ρ_r)`` when fitted, else the ``_JEFFREYS_REF`` reference ``+½·log(1 − f_g)`` → ``(1, K)``.
    This is the **two-group** arm (gDNA vs RNA-total): the per-strand split is the nuisance tilt, integrated
    out on the θ axis, and needs no prior of its own. `logP_r` is not fitted yet — the reference is what bounds
    the ``f_g → 1`` vertex today, and it is the ONLY thing doing so."""
    if rna_logprior is not None:
        return np.asarray(rna_logprior, np.float64)
    return _JEFFREYS_REF * _log1m_fg(lam)[None, :]


def _tilt_grid(n_tilt: int) -> np.ndarray:
    """The RNA-internal tilt grid as the ANGLE ``θ ∈ [−π/2, π/2]`` (``K_t`` points), with ``τ = sin θ``.

    Gridding θ (not τ) is what makes the Berger–Bernardo tilt conditional ``(1−τ²)^{−½}`` vanish identically:
    ``|dτ/dθ| = cos θ = (1−τ²)^{½}`` cancels it exactly, so **no tilt term is written at all** and the ψ
    reference is the same expression for AMBIG as for single-strand nodes (module docstring §3). It also
    removes the endpoint singularity outright — no clipping, no Gauss–Jacobi weights, no constant.

    Resolution follows the reference measure rather than being uniform in τ: at ``K_t=60`` the τ-spacing is
    ~0.053 near balanced tilt and ~0.0014 near strand-purity (vs a flat 0.034). That is the intended trade —
    it spends grid on the strand-pure edges, where distinguishing a pure strand from a small antisense leak is
    the high-stakes call, and economizes on the balanced middle, where the distinction rarely matters.

    ``θ = ±π/2`` ⇒ ``τ = ±1`` ⇒ all RNA on one strand; ``θ = 0`` ⇒ balanced. Only AMBIG nodes integrate it."""
    return np.linspace(-0.5 * np.pi, 0.5 * np.pi, int(n_tilt))


def _single_strand_mask(allow_pos, allow_neg) -> np.ndarray:
    """The nodes the 1-D (Phase-1) solver is valid for: exactly one strand live (tilt determined)."""
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    return ap ^ an


def _ambig_mask(allow_pos, allow_neg) -> np.ndarray:
    """AMBIG nodes (both strands live) — the Phase-2 2-D ``(λ, τ)`` path."""
    return np.asarray(allow_pos, bool) & np.asarray(allow_neg, bool)


def _local_loglik_logodds(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    kappa,
    od_g,
    od_r,
    lam,
    fg,
    f_g_ref,
    f_pos_ref,
    f_neg_ref,
    global_logprior=None,
    gdna_imp_mode=None,
    gdna_imp_prec=None,
    rna_imp_mode=None,
    rna_imp_prec=None,
    lam_logprior=None,
    lam_imp_mode=None,
    lam_imp_prec=None,
):
    """ψ over the log-odds grid for single-strand nodes (strand mixture, the two arms, imputation), evaluated
    at ``f_g = σ(λ)`` with the live strand carrying ``f_active = 1 − f_g``. Returns ``(m, K)``.

    ψ = strand + ``_gdna_arm`` + ``_rna_arm`` + messages. **No Jacobian** — on the two-group axis the log-rate
    conversions cancel ``log σ'(λ)`` exactly (module docstring §2). Both arms are ALWAYS written: a fitted
    ``logP`` if we have one, else the ``_JEFFREYS_REF`` reference. Omitting one is not neutral.

    ``global_logprior`` must already be evaluated on THIS ``fg`` grid → ``(m, K)``; ``None`` ⇒ the gDNA arm
    takes its reference (a PRIOR-FREE solve is not a REFERENCE-FREE solve).

    ``f_g_ref`` / ``f_pos_ref`` / ``f_neg_ref`` (per-node ``(m,)``) are the count-zero-information freeze
    reference (`docs/calibration/archive/node_prior_design.md` §2): the strand mixture's variance is evaluated at THIS fixed
    composition — not the grid ``f_g`` being integrated — so the count sets precision, not composition."""
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    pos_live = (ap & ~an)[:, None]  # (m,1)
    neg_live = (an & ~ap)[:, None]
    fg2 = fg[None, :]  # (1,K)
    f_act = 1.0 - fg2
    f_pos = np.where(pos_live, f_act, 0.0)  # (m,K)
    f_neg = np.where(neg_live, f_act, 0.0)
    n = (u_pos + u_neg)[:, None]
    # ── strand mixture — variance FROZEN at the reference composition (count-zero-info §2). Broadcasts
    #    the (m,1) reference against the (m,K) grid. ──
    psi = _mixture_strand_loglik(
        u_pos[:, None],
        n,
        fg2,
        f_pos,
        f_neg,
        kappa,
        od_g,
        od_r,
        np.asarray(f_g_ref, np.float64)[:, None],
        np.asarray(f_pos_ref, np.float64)[:, None],
        np.asarray(f_neg_ref, np.float64)[:, None],
    )
    # ── the two composition arms: gDNA and RNA-total. Each is its fitted logP if we have one, else the
    #    derived Jeffreys reference. ALWAYS both — see `_gdna_arm` / `_rna_arm`. Together they make ψ proper
    #    (Beta(½,½) when neither is fitted); the gDNA arm alone would leave f_g→1 unbounded, and the RNA arm
    #    alone would leave f_g→0 unbounded. ──
    psi = psi + _gdna_arm(lam, global_logprior) + _rna_arm(lam)
    # ── the gDNA INTRON FACTORY λ-factor (gdna_intron_factory_design.md): a per-node (m,K) log-likelihood on
    #    the λ axis, ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)``, ADDED (not folded into the gDNA arm — that arm
    #    REPLACES the Jeffreys reference; folding would drop the f_g→1 bound). It peels confident gDNA from
    #    introns against the intergenic background; zero on non-intron nodes ⇒ a no-op there. ──
    if lam_logprior is not None:
        psi = psi + np.asarray(lam_logprior, np.float64)
    # ── imputation messages: LOG-FRACTION Gaussians (the overhaul). The mode is a log-FRACTION target
    #    (``log`` of the imputed fraction, built in ``_scan``); evaluated against ``log f_c(λ)``. No clip —
    #    an off-grid target (source denser than the dst can hold) is a bounded monotone pull toward the
    #    edge, governed by precision (D-plan P6, verify-don't-clip). ──
    log_fg = _log_fg(lam)[None, :]  # log f_g = log σ(λ)
    log_fact = _log1m_fg(lam)[None, :]  # log(1−f_g) = log f_active (the single live strand)
    if gdna_imp_mode is not None and gdna_imp_prec is not None:
        m_ = np.asarray(gdna_imp_mode, np.float64)[:, None]
        p_ = np.asarray(gdna_imp_prec, np.float64)[:, None]
        psi = psi - 0.5 * p_ * (log_fg - m_) ** 2
    if rna_imp_mode is not None and rna_imp_prec is not None:
        # single-strand: the live strand carries f_active = 1−f_g; the per-strand precision gates which
        # message applies (the dead strand's prec is 0 → no-op). Both evaluate against log f_active.
        for ms, ps in ((rna_imp_mode[0], rna_imp_prec[0]), (rna_imp_mode[1], rna_imp_prec[1])):
            psi = (
                psi
                - 0.5
                * np.asarray(ps, np.float64)[:, None]
                * (log_fact - np.asarray(ms, np.float64)[:, None]) ** 2
            )
    # ── the SINGLE-λ composition message (the M6 rank-1 fix): ONE Gaussian on the log-odds grid variable λ
    #    DIRECTLY (not on log f_c) — the one gDNA-vs-RNA-total DOF, so ψ counts it ONCE, not twice
    #    (`message_variance_derivation.md` §4). Enrichment-invariant: λ carries no reframe. ──
    if lam_imp_mode is not None and lam_imp_prec is not None:
        lm_ = np.asarray(lam_imp_mode, np.float64)[:, None]
        lp_ = np.asarray(lam_imp_prec, np.float64)[:, None]
        psi = psi - 0.5 * lp_ * (lam[None, :] - lm_) ** 2
    # ── NO change-of-variable Jacobian, and NO reference prior. Both are deliberate, and they are the SAME
    #    fact: `DensityNPMLE.logP` is a density in LOG-rate, so its conversion to a linear-rate density
    #    (−log f_g, up to a constant) cancels log σ'(λ) = log f_g + log(1−f_g) exactly, once per component.
    #    Writing either alone is what produced the improper +0.5·λ ramp (docs/calibration/archive/prior_ramp_and_bp_roadmap.md §2).
    #    ⇒ ψ_λ = strand + logP_g + logP_r, bare. ──
    return psi, f_pos, f_neg


def _solve_nodes_logodds(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    mass_unspl,
    mass_spliced,
    f_g_ref,
    f_pos_ref,
    f_neg_ref,
    *,
    kappa,
    od_g,
    od_r,
    n_grid,
    L: float = _DEFAULT_L,
    global_logprior=None,
    gdna_imp_mode=None,
    gdna_imp_prec=None,
    rna_imp_mode=None,
    rna_imp_prec=None,
    lam_logprior=None,
    lam_imp_mode=None,
    lam_imp_prec=None,
) -> NodeDeconv:
    """The log-odds 1-D per-node solve for SINGLE-STRAND nodes.

    Read-out: ``f_g`` = posterior median over the ``λ`` grid; ``f_pos``/``f_neg`` = posterior MEANS;
    ``*_frac_var`` = the grid-moment ``Var(log f_c)``. The dead strand is locked-certain (var 0); zero-mass
    nodes report 0. ``f_g_ref``/``f_pos_ref``/``f_neg_ref`` (per-node) are the count-zero-info variance
    freeze reference (§2). AMBIG nodes are out of contract — masked out."""
    lam, fg = _logodds_grid(int(n_grid), L)
    psi, f_pos_g, f_neg_g = _local_loglik_logodds(
        u_pos,
        u_neg,
        allow_pos,
        allow_neg,
        kappa,
        od_g,
        od_r,
        lam,
        fg,
        f_g_ref,
        f_pos_ref,
        f_neg_ref,
        global_logprior=global_logprior,
        gdna_imp_mode=gdna_imp_mode,
        gdna_imp_prec=gdna_imp_prec,
        rna_imp_mode=rna_imp_mode,
        rna_imp_prec=rna_imp_prec,
        lam_logprior=lam_logprior,
        lam_imp_mode=lam_imp_mode,
        lam_imp_prec=lam_imp_prec,
    )
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    post = np.exp(psi - _lse(psi, axis=1, keepdims=True))  # (m,K)
    # f_g posterior median (fg ascending ⇒ cumulative CDF directly)
    f_g = _posterior_median_fg(post, fg)
    # composition: f_g median + f_pos/f_neg posterior MEANS (the current-state fractions).
    f_pos = np.sum(post * f_pos_g, axis=1)
    f_neg = np.sum(post * f_neg_g, axis=1)
    # precision state = Var(log f_c), moment-matched on the grid. The dead strand is locked-certain
    # (f=0) → var 0. Capping is AUTOMATIC: the send prec_log = 1/(var+σ²+pois) ≤ 1/(σ²+pois).
    Lg = _log_fg(lam)
    mLg = post @ Lg
    var_g = np.maximum(post @ (Lg * Lg) - mLg * mLg, 0.0)
    La = _log1m_fg(lam)
    mLa = post @ La
    var_act = np.maximum(post @ (La * La) - mLa * mLa, 0.0)
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
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
        gdna_frac_var=var_g,
        rna_pos_frac_var=var_pos,
        rna_neg_frac_var=var_neg,
    )


def _solve_ambig_logodds(
    u_pos,
    u_neg,
    mass_unspl,
    mass_spliced,
    f_g_ref,
    f_pos_ref,
    f_neg_ref,
    *,
    kappa,
    od_g,
    od_r,
    n_grid,
    L: float = _DEFAULT_L,
    n_tilt: int | None = None,
    global_logprior=None,
    gdna_imp_mode=None,
    gdna_imp_prec=None,
    rna_imp_mode=None,
    rna_imp_prec=None,
    lam_logprior=None,
    lam_imp_mode=None,
    lam_imp_prec=None,
    theta_imp_mode=None,
    theta_imp_prec=None,
) -> NodeDeconv:
    """The 2-D ``(λ, θ)`` solve for AMBIG nodes (both strands live). Grids the gDNA-vs-RNA-total log-odds
    ``λ`` (outer, ``K = n_grid``) and the tilt ANGLE ``θ = arcsin(τ)`` (inner, ``K_t = n_tilt`` or
    ``n_grid``), evaluates ψ on the ``(m, K, K_t)`` cube, and **marginalizes θ** (``logsumexp``) for the
    ``f_g`` read-out. ``f_g`` = posterior median over the θ-marginal λ-posterior; ``f_pos``/``f_neg`` = means
    over the full 2-D posterior.

    **No Jacobian and no tilt term are written** — and that is the point of the θ coordinate. The
    Berger–Bernardo tilt conditional ``(1−τ²)^{−½}`` is cancelled identically by ``|dτ/dθ| = cos θ``
    (``_tilt_grid``), and on the two-group λ axis the log-rate conversions cancel ``log σ'(λ)``. So ψ here is
    the SAME expression as the 1-D path: strand + ``_gdna_arm`` + ``_rna_arm`` + messages. The 1-DOF/AMBIG
    reference asymmetry is closed **identically**, not approximately.

    ``global_logprior`` is ``(m, K)`` evaluated on the σ(λ) grid (broadcast over θ). The cube is only
    materialized for the AMBIG subset (the caller masks); ``K·K_t`` is the per-node cost."""
    lam, fg = _logodds_grid(int(n_grid), L)  # (K,)
    Kt = int(n_tilt) if n_tilt else int(n_grid)
    theta = _tilt_grid(Kt)  # (Kt,) the ANGLE; τ = sin θ
    tau = np.sin(theta)  # (Kt,) the tilt itself — smooth and exact across the whole θ domain
    f_act = (1.0 - fg)[:, None]  # (K,1)
    f_pos_kt = f_act * (1.0 + tau[None, :]) / 2.0  # (K,Kt) f64 (reused f64 in the moment sums)
    f_neg_kt = f_act * (1.0 - tau[None, :]) / 2.0  # (K,Kt) f64
    n = u_pos + u_neg
    # ── float32 cube (authorized small-tolerance perf): the (m,K,Kt) log-posterior and its per-cell inputs
    #    are f32 — halves the cube's memory, ~2x the elementwise exp/log. Every REDUCTION below accumulates
    #    in f64 and the (m,K) marginals are lifted to f64, so medians / means / variances keep precision. ──
    F = np.float32
    fg32 = fg.astype(F)
    fpk = f_pos_kt.astype(F)  # (K,Kt) f32 cube grid
    fnk = f_neg_kt.astype(F)
    # ── strand mixture over the cube: (m,1,1)×(1,K,1)×(1,K,Kt) → (m,K,Kt) f32 ──
    psi = _mixture_strand_loglik(
        np.asarray(u_pos, F)[:, None, None],
        np.asarray(n, F)[:, None, None],
        fg32[None, :, None],
        fpk[None, :, :],
        fnk[None, :, :],
        kappa,
        od_g,
        od_r,
        np.asarray(f_g_ref, F)[:, None, None],
        np.asarray(f_pos_ref, F)[:, None, None],
        np.asarray(f_neg_ref, F)[:, None, None],
    )
    # ── LOG-fraction grids (the overhaul): log f_g (τ-independent) + log f_pos/f_neg over the cube,
    #    floored at one pseudo-fragment 1/(n+1) (D5: the τ=±1 edges have f_s=0 → log(0); the count floor
    #    keeps it finite + consistent with pois_log). ──
    log_fg_grid = _log_fg(lam)  # (K,) f64 = log f_g (moments use f64)
    log_fg32 = log_fg_grid.astype(F)  # (K,) f32 for the cube message
    frac_floor = (1.0 / (n + 1.0)).astype(F)[:, None, None]  # (m,1,1) f32
    log_fpos = np.log(np.maximum(fpk[None, :, :], frac_floor))  # (m,K,Kt) f32
    log_fneg = np.log(np.maximum(fnk[None, :, :], frac_floor))
    # ── the two composition arms (θ-independent — they live on the λ axis, which is exactly what makes the
    #    tilt a nuisance). Identical call to the 1-D path; see `_gdna_arm` / `_rna_arm`. ──
    psi += np.asarray(_gdna_arm(lam, global_logprior) + _rna_arm(lam), F)[:, :, None]
    # ── the gDNA INTRON FACTORY λ-factor (θ-independent — it lives on the λ axis), ADDED like the arms; the
    #    [:, :, None] broadcast makes it constant across the tilt, so θ is integrated out cleanly. ──
    if lam_logprior is not None:
        psi += np.asarray(lam_logprior, F)[:, :, None]
    # ── gDNA LOG-fraction message on log f_g (τ-independent) ──
    if gdna_imp_mode is not None and gdna_imp_prec is not None:
        mo = np.asarray(gdna_imp_mode, F)[:, None, None]
        pr = np.asarray(gdna_imp_prec, F)[:, None, None]
        psi -= F(0.5) * pr * (log_fg32[None, :, None] - mo) ** 2
    # ── per-strand RNA LOG-fraction messages on log f_pos/log f_neg (τ-dependent — inside the cube) ──
    if rna_imp_mode is not None and rna_imp_prec is not None:
        for log_f, ms, ps in (
            (log_fpos, rna_imp_mode[0], rna_imp_prec[0]),
            (log_fneg, rna_imp_mode[1], rna_imp_prec[1]),
        ):
            psi -= (
                F(0.5)
                * np.asarray(ps, F)[:, None, None]
                * (log_f - np.asarray(ms, F)[:, None, None]) ** 2
            )
    # ── the SINGLE-λ composition message on λ DIRECTLY (θ-INDEPENDENT — it lives on the λ axis, which is
    #    exactly what makes the tilt a nuisance): one Gaussian, ψ counts the g-vs-R DOF ONCE. ──
    if lam_imp_mode is not None and lam_imp_prec is not None:
        lm_ = np.asarray(lam_imp_mode, F)[:, None, None]
        lp_ = np.asarray(lam_imp_prec, F)[:, None, None]
        psi -= F(0.5) * lp_ * (lam.astype(F)[None, :, None] - lm_) ** 2
    # ── the TILT message on θ (λ-INDEPENDENT — the separate strand-tilt DOF an AMBIG node needs; not part of
    #    the g-vs-R double-count): a Gaussian on the θ = arcsin(τ) grid. ──
    if theta_imp_mode is not None and theta_imp_prec is not None:
        tm_ = np.asarray(theta_imp_mode, F)[:, None, None]
        tp_ = np.asarray(theta_imp_prec, F)[:, None, None]
        psi -= F(0.5) * tp_ * (theta.astype(F)[None, None, :] - tm_) ** 2
    psi_full = psi  # (m,K,Kt) f32
    # θ-marginal λ-posterior (m,K) — lift to f64 so the posterior median + moments are full-precision.
    psi_lam = _lse(psi_full, axis=2).astype(np.float64)
    post_lam = np.exp(psi_lam - _lse(psi_lam, axis=1, keepdims=True))
    f_g = _posterior_median_fg(post_lam, fg)
    # precision state = Var(log f_g) over the θ-marginal λ-posterior (D2).
    mLg = post_lam @ log_fg_grid
    var_g = np.maximum(post_lam @ (log_fg_grid * log_fg_grid) - mLg * mLg, 0.0)
    # f_pos/f_neg MEANS + Var(log f_pos/neg) over the FULL 2-D posterior (f32 cube; sums accumulate in f64).
    flat = psi_full.reshape(psi_full.shape[0], -1)
    post2d = np.exp(flat - _lse(flat, axis=1, keepdims=True)).reshape(
        psi_full.shape
    )  # (m,K,Kt) f32
    fp_grid = fpk[None, :, :]
    fn_grid = fnk[None, :, :]
    f_pos = np.sum(post2d * fp_grid, axis=(1, 2), dtype=np.float64)
    f_neg = np.sum(post2d * fn_grid, axis=(1, 2), dtype=np.float64)
    mLp = np.sum(post2d * log_fpos, axis=(1, 2), dtype=np.float64)
    mLn = np.sum(post2d * log_fneg, axis=(1, 2), dtype=np.float64)
    var_pos = np.maximum(
        np.sum(post2d * log_fpos * log_fpos, axis=(1, 2), dtype=np.float64) - mLp * mLp, 0.0
    )
    var_neg = np.maximum(
        np.sum(post2d * log_fneg * log_fneg, axis=(1, 2), dtype=np.float64) - mLn * mLn, 0.0
    )
    active = n > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_pos = np.where(active, np.clip(f_pos, 0.0, 1.0), 0.0)
    f_neg = np.where(active, np.clip(f_neg, 0.0, 1.0), 0.0)
    var_g = np.where(active, var_g, 0.0)
    var_pos = np.where(active, var_pos, 0.0)
    var_neg = np.where(active, var_neg, 0.0)
    return NodeDeconv(
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
        gdna_frac_var=var_g,
        rna_pos_frac_var=var_pos,
        rna_neg_frac_var=var_neg,
    )


def _solve_nodes_logodds_all(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    mass_unspl,
    mass_spliced,
    *,
    kappa,
    od_g,
    od_r,
    n_grid,
    L: float = _DEFAULT_L,
    n_tilt: int | None = None,
    n_grid_ss: int | None = None,
    global_logprior=None,
    gdna_imp_mode=None,
    gdna_imp_prec=None,
    rna_imp_mode=None,
    rna_imp_prec=None,
    lam_logprior=None,
    lam_imp_mode=None,
    lam_imp_prec=None,
    theta_imp_mode=None,
    theta_imp_prec=None,
    fg_ref=None,
    fpos_ref=None,
    fneg_ref=None,
) -> NodeDeconv:
    """The full per-node log-odds dispatcher (Phase 3 #1): routes single-strand nodes to the 1-D
    ``λ`` solve (:func:`_solve_nodes_logodds`) and AMBIG nodes to the 2-D ``(λ, τ)`` solve
    (:func:`_solve_ambig_logodds`), scattering both into full-length arrays. G1 / zero-mass nodes
    report 0 (``node_sweep`` keeps their signature-binary init via the ``solvable`` write-back). A
    drop-in for the lattice ``_local_loglik``+``_node_marginals`` pair: same ψ terms — the log-density
    log-fraction Gaussian messages + the global prior — evaluated on the ``σ(λ)`` log-odds grid.

    All array inputs are full length ``m``; ``global_logprior`` is ``(m, K)`` on the σ(λ) grid;
    ``gdna_imp_*`` are ``(m,)``; ``rna_imp_*`` are 2-tuples of ``(m,)``. Each is sub-indexed per class."""
    m = int(np.asarray(u_pos).shape[0])
    ap_all = np.asarray(allow_pos, bool)
    an_all = np.asarray(allow_neg, bool)
    # Count-zero-info variance-freeze reference (§2, B1). Supplied by the sweep as the incoming belief; at
    # init (None) use the structural-neutral default: f_g=½ with the remaining ½ split among the LIVE
    # strands (single-strand → ½ on its strand; AMBIG → ¼ each). Location is prior/likelihood-set; the
    # reference only fixes the variance/precision.
    if fg_ref is None or fpos_ref is None or fneg_ref is None:
        nlive = ap_all.astype(np.float64) + an_all.astype(np.float64)
        half = np.where(nlive > 0.0, 0.5 / np.maximum(nlive, 1.0), 0.0)
        fg_ref = np.full(m, 0.5)
        fpos_ref = np.where(ap_all, half, 0.0)
        fneg_ref = np.where(an_all, half, 0.0)
    else:
        fg_ref = np.asarray(fg_ref, np.float64)
        fpos_ref = np.asarray(fpos_ref, np.float64)
        fneg_ref = np.asarray(fneg_ref, np.float64)
    out = {
        k: np.zeros(m, dtype=np.float64)
        for k in ("fg", "fp", "fn", "vg", "vp", "vn")
    }
    # Skip EMPTY nodes — no per-strand counts AND no unspliced/spliced mass. Both per-class solvers zero
    # every output for an inactive node (gdna/rna_mass = f_g·M = (1−f_g)·M + S = 0 when all are 0), so an
    # empty node's solve is identical to the zero-initialized `out` — skipping is BIT-IDENTICAL. At genome
    # scale most region/boundary nodes carry no fragments (unexpressed genes, intergenic deserts), so this
    # is the dominant cost saver, not a slice artifact. (A spliced-only node has signal ⇒ still solved.)
    signal = (
        np.asarray(u_pos, np.float64)
        + np.asarray(u_neg, np.float64)
        + np.asarray(mass_unspl, np.float64)
        + np.asarray(mass_spliced, np.float64)
    ) > 0.0
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

    if bool(ss.any()):
        # Single-strand nodes solve on the FINE 1-D grid (Fix 3, n_grid_ss); the coarse-grid global prior is
        # regridded onto it. AMBIG keeps the coarse n_grid (the expensive 2-D cube). n_grid_ss=None ⇒ n_grid.
        k_ss = int(n_grid_ss) if n_grid_ss else int(n_grid)
        _scatter(
            ss,
            _solve_nodes_logodds(
                u_pos[ss],
                u_neg[ss],
                allow_pos[ss],
                allow_neg[ss],
                mass_unspl[ss],
                mass_spliced[ss],
                fg_ref[ss],
                fpos_ref[ss],
                fneg_ref[ss],
                kappa=kappa,
                od_g=od_g,
                od_r=od_r,
                n_grid=k_ss,
                L=L,
                global_logprior=_regrid_global(_s(global_logprior, ss), n_grid, k_ss, L),
                gdna_imp_mode=_s(gdna_imp_mode, ss),
                gdna_imp_prec=_s(gdna_imp_prec, ss),
                rna_imp_mode=_sp(rna_imp_mode, ss),
                rna_imp_prec=_sp(rna_imp_prec, ss),
                lam_logprior=_regrid_global(_s(lam_logprior, ss), n_grid, k_ss, L),
                lam_imp_mode=_s(lam_imp_mode, ss),
                lam_imp_prec=_s(lam_imp_prec, ss),
            ),
        )
    if bool(amb.any()):
        # The 2-D (λ,τ) cube is (B,K,K_t); materialized for ALL ambig nodes at once it is ~O(m·K²) (the
        # memory the lattice OOM'd on). AMBIG nodes solve independently, so tile the subset into row-batches
        # — bit-identical results, peak memory bounded to one (≤_AMBIG_BATCH, K, K_t) cube.
        amb_idx = np.where(amb)[0]
        for s0 in range(0, amb_idx.size, _AMBIG_BATCH):
            bidx = amb_idx[s0 : s0 + _AMBIG_BATCH]
            _scatter(
                bidx,
                _solve_ambig_logodds(
                    u_pos[bidx],
                    u_neg[bidx],
                    mass_unspl[bidx],
                    mass_spliced[bidx],
                    fg_ref[bidx],
                    fpos_ref[bidx],
                    fneg_ref[bidx],
                    kappa=kappa,
                    od_g=od_g,
                    od_r=od_r,
                    n_grid=n_grid,
                    L=L,
                    n_tilt=n_tilt,
                    global_logprior=_s(global_logprior, bidx),
                    gdna_imp_mode=_s(gdna_imp_mode, bidx),
                    gdna_imp_prec=_s(gdna_imp_prec, bidx),
                    rna_imp_mode=_sp(rna_imp_mode, bidx),
                    rna_imp_prec=_sp(rna_imp_prec, bidx),
                    lam_logprior=_s(lam_logprior, bidx),
                    lam_imp_mode=_s(lam_imp_mode, bidx),
                    lam_imp_prec=_s(lam_imp_prec, bidx),
                    theta_imp_mode=_s(theta_imp_mode, bidx),
                    theta_imp_prec=_s(theta_imp_prec, bidx),
                ),
            )
    return NodeDeconv(
        gdna_frac=out["fg"],
        rna_pos_frac=out["fp"],
        rna_neg_frac=out["fn"],
        gdna_frac_var=out["vg"],
        rna_pos_frac_var=out["vp"],
        rna_neg_frac_var=out["vn"],
    )
