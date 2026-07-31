"""Generic DENSITY DECONVOLUTION — deconvolve a node's counts into gDNA + RNA against a gDNA prior.

A node's unspliced count ``C`` (over its gDNA effective length ``E_g``) is ``gDNA + RNA``. If we know the
**gDNA density prior** — how dense gDNA is at this node, as a distribution ``π_bg`` — we can peel the gDNA
(``g ≈ ρ_bg·E_g``) and read the residual as RNA, with honest, count-derived precision. This is the generic
count-deconvolution primitive; the **intron factory** is its special case (`fit_intron_background`), where the
gDNA prior is the intergenic node distribution (introns are off-target, at the same capture depletion as
intergenic — `docs/CARRY_FORWARD.md` §2).

Model (owner-ratified 2026-07-20): the gDNA count ``g ~ NegBinom(mean = ρ_bg·E_g, size = α_eff)``, a
Gamma-Poisson — the per-region background rate ``ρ ~ Gamma`` (over-dispersion ``α`` = per-region CNV /
mappability spread, fitted from the pool) mixed by the per-node Poisson sampling. The observed ``C = g + r``
with a **flat one-sided RNA prior** ``r ≥ 0`` (RNA carries no informative density prior, and this dodges any
cutoff), so the posterior on the gDNA fraction is ``P(g | C) ∝ NegBinom(g; ρ_bg·E_g, α_eff)·1[g ≤ C]``. On the
``f_g`` solve grid this is the factor ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` (the truncation is automatic
since ``f_g ≤ 1 ⇒ g ≤ C``). It peels **gDNA only** (gDNA is strand-symmetric, so no strand is assigned; the
residual RNA's strand is the tilt ``θ``, left to the strand solver — the SYNERGY of the two deconvolutions).

Precision is honest and count-over-length, with NO tuned constant: three variance sources add as
``Var(g) = μ + μ²/α_eff`` with ``1/α_eff = 1/α + 1/(Σg + n0)`` — the per-node Poisson (``μ``), the per-region
over-dispersion (``α``, fitted), and the pooled-mean/resolution uncertainty (``Σg + n0``, the floor-blend
Fisher info from ``gdna_background_floor_derivation``). Graceful ``Σg → 0``: ``μ`` falls back to the resolution
wall and ``α_eff`` widens, so the peel becomes imprecise (never confidently wrong) as capture depletes the pool.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from .background_reference import measure_background
from .signature import coarse_type_array

__all__ = [
    "GdnaBackground",
    "fit_gdna_background",
    "fit_intron_background",
    "density_lambda_factor",
    "density_factor_precision",
]

_EPS = 1.0e-12


@dataclass(frozen=True, slots=True)
class GdnaBackground:
    """The fitted gDNA density prior the density deconvolution scores each node against.

    ``log_mu_bg`` is the natural-log background gDNA DENSITY location: ``log(Σg/ΣE)`` (the pooled rate) when the
    pool has counts, else the resolution wall ``log ρ_res`` (``Σg = 0``). ``alpha`` is the NegBinom size = the
    per-region over-dispersion (``+inf`` = Poisson-clean, no CNV/mappability spread). ``sg``/``n0`` are the
    pooled counts and the zero-count-region tally — the pooled-mean/resolution Fisher info ``Σg + n0`` that
    combines with ``alpha`` into the effective size (§ module docstring). ``informative`` is False when the pool
    is empty (no background at all) ⇒ the factor is flat."""

    log_mu_bg: float
    alpha: float
    sg: float
    n0: float
    n_regions: int
    informative: bool


def fit_gdna_background(g_counts, eff_g, *, log_rho_floor: float) -> GdnaBackground:
    """The GENERIC fit: a gDNA NegBinom background from a pool of **pure-gDNA** nodes.

    ``g_counts``/``eff_g`` are the pooled nodes' gDNA counts and gDNA effective lengths; ``log_rho_floor`` is
    the resolution-wall location for the ``Σg = 0`` case (from :func:`background_reference.measure_background`).
    ``α`` by method-of-moments on the pool given the pooled mean ``ρ_bg = Σg/ΣE`` (so ``μ_i = ρ_bg·E_i`` and
    the residuals sum to zero at the MLE): from ``E[Σ(g_i − μ_i)²] = Σμ_i + (1/α)·Σμ_i²``,
    ``α = Σμ_i² / max(Σ(g_i − μ_i)² − Σμ_i, 0⁺)`` — ``+inf`` (Poisson) when the pool is not over-dispersed. No
    tuned constant. Callers select the pool (introns → intergenic; post-pass-0 exons → clean-gDNA nodes)."""
    g = np.asarray(g_counts, dtype=np.float64)
    E = np.asarray(eff_g, dtype=np.float64)
    sg = float(g.sum())
    se = float(E.sum())
    n0 = float(np.count_nonzero(g <= _EPS))
    n_reg = int(g.shape[0])

    if sg > _EPS and se > _EPS:
        log_mu_bg = float(np.log(sg) - np.log(se))
    else:
        log_mu_bg = float(log_rho_floor)  # Σg=0 ⇒ the wall (−inf iff no support at all)
    informative = bool(np.isfinite(log_mu_bg) and n_reg > 0)

    alpha = np.inf
    if informative and sg > _EPS:
        rho_bg = sg / se
        mu = rho_bg * E
        S = float(np.sum((g - mu) ** 2))  # Σ residual²
        sum_mu = float(np.sum(mu))
        sum_mu2 = float(np.sum(mu * mu))
        excess = S - sum_mu  # over-dispersion beyond Poisson
        if excess > _EPS and sum_mu2 > _EPS:
            alpha = sum_mu2 / excess  # finite ⇒ over-dispersed; else Poisson (α=∞)

    return GdnaBackground(
        log_mu_bg=log_mu_bg,
        alpha=float(alpha),
        sg=sg,
        n0=n0,
        n_regions=n_reg,
        informative=informative,
    )


def fit_intron_background(
    substrate, region_arrays, node_eff_g, *, include_introns: bool = False
) -> GdnaBackground:
    """The INTRON special case of :func:`fit_gdna_background`: the gDNA prior IS the intergenic node
    distribution (introns are off-target at the same capture depletion as intergenic —
    `docs/CARRY_FORWARD.md` §2).

    ``include_introns=False`` (default) pools **intergenic only** — the clean, non-circular reference (introns
    are what we deconvolve). The real-data path may add RNA-free introns for resolution
    (`background_reference` docstring), traded for a little circularity. Delegates the NB fit to
    :func:`fit_gdna_background`; the intergenic pool + resolution floor come from
    :func:`background_reference.measure_background`.

    ``node_eff_g`` is the gDNA CONTAINED effective length per node
    (:func:`effective_length.contained_eff_length`) — the same array ``measure_background`` pools over,
    passed through so the fit and the floor sit on one support."""
    bg = measure_background(substrate, region_arrays, node_eff_g, include_introns=include_introns)
    sig = np.asarray(region_arrays.signature)
    eff = np.asarray(node_eff_g, dtype=np.float64)
    ctype = coarse_type_array(sig)  # 0 intergenic / 1 intron / 2 exon
    # the pure-gDNA pool: intergenic only (or non-exonic when include_introns), matching measure_background.
    pool = ((ctype != 2) if include_introns else (ctype == 0)) & (eff > _EPS)
    # ⚠ GENOME-strand columns, summed: gDNA is strand-symmetric, so the background is a total rate.
    counts = np.asarray(substrate.node_contained.count, dtype=np.float64).sum(axis=1)
    return fit_gdna_background(counts[pool], eff[pool], log_rho_floor=bg.log_rho_floor)


def _log_negbinom(g: np.ndarray, mu: np.ndarray, size: float) -> np.ndarray:
    """``log NegBinom(g; mean=mu, size)`` for CONTINUOUS ``g ≥ 0`` (the ``f_g``-grid gives fractional counts).

    Mean-``μ`` / size-``r`` parameterization: ``p = r/(r+μ)``,
    ``log NB = Γln(g+r) − Γln(r) − Γln(g+1) + r·log(r/(r+μ)) + g·log(μ/(r+μ))``. ``r → ∞`` is the exact Poisson
    limit ``g·log μ − μ − Γln(g+1)`` (used directly to avoid the ``Γln(∞)`` overflow)."""
    g = np.asarray(g, dtype=np.float64)
    mu = np.maximum(np.asarray(mu, dtype=np.float64), _EPS)
    if not np.isfinite(size):  # α = ∞ ⇒ Poisson
        return g * np.log(mu) - mu - gammaln(g + 1.0)
    r = max(float(size), _EPS)
    rpm = r + mu
    return (
        gammaln(g + r)
        - gammaln(r)
        - gammaln(g + 1.0)
        + r * (np.log(r) - np.log(rpm))
        + g * (np.log(mu) - np.log(rpm))
    )


def density_lambda_factor(
    background: GdnaBackground, count: np.ndarray, eff_g: np.ndarray, fg_grid: np.ndarray
) -> np.ndarray:
    """The per-node factor on ``λ`` over the ``f_g = σ(λ)`` solve grid: ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)``.

    ``count`` (``C``), ``eff_g`` (``E_g``) are per-node arrays ``(n,)``; ``fg_grid`` is the ``(K,)`` grid.
    Returns ``(n, K)``, peaked at ``f_g = ρ_bg/ρ_obs`` (the confident gDNA peel), with the honest
    ``μ + μ²/α_eff`` curvature. Each row is offset so its max is 0 (an ``f_g``-independent constant is
    irrelevant to ψ; this only keeps the numbers well-scaled). A non-informative background returns all-zero."""
    C = np.asarray(count, dtype=np.float64)
    Eg = np.maximum(np.asarray(eff_g, dtype=np.float64), _EPS)
    fg = np.clip(np.asarray(fg_grid, dtype=np.float64), _EPS, 1.0 - _EPS)
    n, K = C.shape[0], fg.shape[0]
    if not background.informative:
        return np.zeros((n, K), dtype=np.float64)

    # effective size: 1/α_eff = 1/α + 1/(Σg + n0). Per-region over-dispersion ⊕ pooled-mean/resolution info.
    denom = background.sg + background.n0
    inv_alpha = 0.0 if not np.isfinite(background.alpha) else 1.0 / max(background.alpha, _EPS)
    inv_alpha += (1.0 / denom) if denom > _EPS else 0.0
    alpha_eff = np.inf if inv_alpha <= _EPS else 1.0 / inv_alpha

    mu = np.exp(background.log_mu_bg) * Eg  # (n,) background gDNA count location
    g = fg[None, :] * C[:, None]  # (n, K) gDNA count at each grid point
    logf = _log_negbinom(g, mu[:, None], alpha_eff)  # (n, K)
    logf = logf - np.max(logf, axis=1, keepdims=True)  # per-node offset (ψ-irrelevant constant)
    return logf


def density_factor_precision(lam_logprior, lam_grid):
    """``I_density`` — the composition evidence a density-deconvolution ``λ``-factor carries, read off its own
    CURVATURE (`gdna_intron_factory_design.md` §4).

    The factor :func:`density_lambda_factor` is a genuine, reference-FREE likelihood on ``λ`` (external
    ``ρ_bg`` information about this node's composition), so its precision belongs on the composition ``λ``-axis
    alongside the strand evidence. ``τ_λ = 1/Var_λ`` under the normalized factor; the NegBinom
    ``Var(g) = μ + μ²/α_eff`` makes a low-count / high-overdispersion node peel IMPRECISELY and a dense one
    confidently, so it self-limits on thin data. A FLAT row (a non-deconvolved node, or an uninformative
    background) carries NO information ⇒ ``τ = 0``. Returns ``(m,)`` (all-zero when ``lam_logprior`` is None)."""
    if lam_logprior is None:
        return None
    lp = np.asarray(lam_logprior, dtype=np.float64)
    lam = np.asarray(lam_grid, dtype=np.float64)
    live = (
        np.ptp(lp, axis=1) > _EPS
    )  # a flat factor carries NO information (τ=0), never the grid's own width
    tau = np.zeros(lp.shape[0], dtype=np.float64)
    if not bool(live.any()):
        return tau
    w = np.exp(lp[live] - np.max(lp[live], axis=1, keepdims=True))
    w /= np.maximum(w.sum(axis=1, keepdims=True), _EPS)
    mu = w @ lam
    var = w @ (lam * lam) - mu * mu
    tau[live] = np.where(var > _EPS, 1.0 / np.maximum(var, _EPS), 0.0)
    return tau
