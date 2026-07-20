"""The gDNA intron factory — peel confident gDNA out of introns against the intergenic background.

Derivation + design: ``docs/calibration/gdna_intron_factory_design.md``. Introns are **off-target** (the same
capture depletion as intergenic regions), so the pooled intergenic gDNA rate ``ρ_bg`` (``background_reference``)
is the intron's *actual* gDNA density — a **two-sided** estimate, not the ``§8`` one-sided floor that only
bounds on-target exons. This module deconvolves each intron's unspliced count ``C`` (over its gDNA effective
length ``E_g``) into gDNA (explained by the background) + nascent RNA (the excess), producing a per-intron factor
on ``λ = logit(f_g)`` — the gDNA-vs-RNA-total log-odds — with **honest, count-derived precision**. It **peels
gDNA, not RNA** (gDNA is strand-symmetric, so no strand is assigned; the residual RNA's strand is the tilt ``θ``,
left to the solver).

Model (owner-ratified 2026-07-20): the intron gDNA count ``g ~ NegBinom(mean = ρ_bg·E_g, size = α_eff)``, a
Gamma-Poisson — the per-region background rate ``ρ ~ Gamma`` (over-dispersion ``α`` = per-region CNV/mappability
spread, fitted from the intergenic pool) mixed by the per-intron Poisson sampling. The observed
``C = g + r`` with a **flat one-sided nascent prior** ``r ≥ 0`` (nascent is scarce and region-specific — no
informative prior, and this dodges any cutoff), so the posterior on the gDNA fraction is
``P(g | C) ∝ NegBinom(g; ρ_bg·E_g, α_eff)·1[g ≤ C]``. On the ``f_g`` solve grid this is the factor
``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` (the truncation is automatic since ``f_g ≤ 1 ⇒ g ≤ C``).

Precision is honest and count-over-length, with NO tuned constant: three variance sources add as
``Var(g) = μ + μ²/α_eff`` with ``1/α_eff = 1/α + 1/(Σg + n0)`` — the per-intron Poisson (``μ``), the per-region
over-dispersion (``α``, fitted), and the pooled-mean/resolution uncertainty (``Σg + n0``, the floor-blend Fisher
info from ``gdna_background_floor_derivation``). Graceful ``Σg → 0``: ``μ`` falls back to the resolution wall and
``α_eff`` widens, so the peel becomes imprecise (never confidently wrong) as capture depletes the off-target.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from .background_reference import measure_background
from .signature import coarse_type_array

__all__ = ["IntronBackground", "fit_intron_background", "intron_lambda_factor"]

_EPS = 1.0e-12


@dataclass(frozen=True, slots=True)
class IntronBackground:
    """The fitted background model the intron factory scores each intron against.

    ``log_mu_bg`` is the natural-log background gDNA DENSITY location: ``log(Σg/ΣE)`` (the pooled off-target rate)
    when the pool has counts, else the resolution wall ``log ρ_res`` (``Σg = 0``). ``alpha`` is the NegBinom size
    = the per-region over-dispersion (``+inf`` = Poisson-clean, no CNV/mappability spread). ``sg``/``n0`` are the
    pooled counts and the zero-count-region tally — the pooled-mean/resolution Fisher info ``Σg + n0`` that
    combines with ``alpha`` into the effective size (§ module docstring). ``informative`` is False when the pool
    is empty (no background at all) ⇒ the factor is flat."""

    log_mu_bg: float
    alpha: float
    sg: float
    n0: float
    n_regions: int
    informative: bool


def fit_intron_background(
    substrate,
    region_arrays,
    region_eff_g,
    *,
    include_introns: bool = False,
) -> IntronBackground:
    """Fit the background NegBinom from the pure-gDNA pool: ``ρ_bg`` (via :func:`measure_background`) + the
    per-region over-dispersion ``α`` (method-of-moments on the SAME pool).

    ``include_introns=False`` (default) pools **intergenic only** — the clean, non-circular reference (introns
    are what we deconvolve; nascent-bearing introns would inflate the over-dispersion). The real-data path may
    add RNA-free introns for resolution (``background_reference`` docstring), traded for a little circularity.

    ``α`` by method-of-moments on the intergenic counts ``g_i`` over gDNA eff-length ``E_i``, given the pooled
    mean ``ρ_bg = Σg/ΣE`` (so ``μ_i = ρ_bg·E_i`` and the residuals sum to zero at the MLE): from
    ``E[Σ(g_i − μ_i)²] = Σμ_i + (1/α)·Σμ_i²``,
    ``α = Σμ_i² / max(Σ(g_i − μ_i)² − Σμ_i, 0⁺)`` — ``+inf`` (Poisson) when the pool is not over-dispersed. No
    tuned constant: ``α`` is fitted, its Poisson limit is the natural ``α → ∞``.
    """
    bg = measure_background(
        substrate, region_arrays, region_eff_g, include_introns=include_introns
    )
    sig = np.asarray(region_arrays.signature)
    eff = np.asarray(region_eff_g, dtype=np.float64)
    ctype = coarse_type_array(sig)  # 0 intergenic / 1 intron / 2 exon
    # the pure-gDNA pool: intergenic only (or non-exonic when include_introns), matching measure_background.
    pool = ((ctype != 2) if include_introns else (ctype == 0)) & (eff > _EPS)
    counts = (
        np.asarray(substrate.contained.n_unspliced_pos, dtype=np.float64)
        + np.asarray(substrate.contained.n_unspliced_neg, dtype=np.float64)
    )
    g = counts[pool]
    E = eff[pool]
    sg = float(g.sum())
    se = float(E.sum())
    n0 = float(np.count_nonzero(g <= _EPS))
    n_reg = int(pool.sum())

    # log background gDNA DENSITY location μ_bg: the pooled MEAN Σg/ΣE where it has counts (the intron's expected
    # off-target gDNA density), the resolution wall at Σg=0 (measure_background.log_rho_floor).
    if sg > _EPS and se > _EPS:
        log_mu_bg = float(np.log(sg) - np.log(se))
    else:
        log_mu_bg = float(bg.log_rho_floor)  # Σg=0 ⇒ the wall (−inf iff no support at all)
    informative = np.isfinite(log_mu_bg) and n_reg > 0

    # α by method-of-moments (given ρ_bg). μ_i = ρ_bg·E_i; residuals sum to 0 at the pooled MLE.
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

    return IntronBackground(
        log_mu_bg=log_mu_bg,
        alpha=float(alpha),
        sg=sg,
        n0=n0,
        n_regions=n_reg,
        informative=bool(informative),
    )


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


def intron_lambda_factor(
    background: IntronBackground,
    count: np.ndarray,
    eff_g: np.ndarray,
    fg_grid: np.ndarray,
) -> np.ndarray:
    """The per-intron factor on ``λ`` over the ``f_g = σ(λ)`` solve grid: ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)``.

    ``count`` (``C``), ``eff_g`` (``E_g``) are per-intron arrays ``(n,)``; ``fg_grid`` is the ``(K,)`` grid.
    Returns ``(n, K)``, peaked at ``f_g = ρ_bg/ρ_obs`` (the confident gDNA peel), with the honest ``μ + μ²/α_eff``
    curvature. Each row is offset so its max is 0 (an ``f_g``-independent constant is irrelevant to ψ; this only
    keeps the numbers well-scaled). A non-informative background (empty pool) returns an all-zero (flat) factor.
    """
    C = np.asarray(count, dtype=np.float64)
    Eg = np.maximum(np.asarray(eff_g, dtype=np.float64), _EPS)
    fg = np.clip(np.asarray(fg_grid, dtype=np.float64), _EPS, 1.0 - _EPS)
    n, K = C.shape[0], fg.shape[0]
    if not background.informative:
        return np.zeros((n, K), dtype=np.float64)

    # effective size: 1/α_eff = 1/α + 1/(Σg + n0). Per-region over-dispersion ⊕ pooled-mean/resolution info.
    denom = background.sg + background.n0
    inv_alpha = (0.0 if not np.isfinite(background.alpha) else 1.0 / max(background.alpha, _EPS))
    inv_alpha += (1.0 / denom) if denom > _EPS else 0.0
    alpha_eff = np.inf if inv_alpha <= _EPS else 1.0 / inv_alpha

    mu = np.exp(background.log_mu_bg) * Eg  # (n,) background gDNA count location
    g = fg[None, :] * C[:, None]  # (n, K) gDNA count at each grid point
    logf = _log_negbinom(g, mu[:, None], alpha_eff)  # (n, K)
    logf = logf - np.max(logf, axis=1, keepdims=True)  # per-node offset (ψ-irrelevant constant)
    return logf
