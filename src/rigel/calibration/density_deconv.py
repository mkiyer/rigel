"""Generic DENSITY DECONVOLUTION — deconvolve a region's counts into gDNA + RNA against a gDNA prior.

A region's unspliced count ``C`` (over its gDNA effective length ``E_g``) is ``gDNA + RNA``. If we know the
**gDNA density prior** — how dense gDNA is at this region, as a distribution ``π_bg`` — we can deconvolve the gDNA
(``g ≈ ρ_bg·E_g``) and read the residual as RNA, with honest, count-derived precision. This is the generic
count-deconvolution primitive; the **intron factory** is its special case (`fit_intron_background`), where the
gDNA prior is the intergenic region distribution (introns are off-target, at the same capture depletion as
intergenic —).

Model (owner-ratified 2026-07-20): the gDNA count ``g ~ NegBinom(mean = ρ_bg·E_g, size = α_eff)``, a
Gamma-Poisson — the per-region background rate ``ρ ~ Gamma`` (over-dispersion ``α`` = per-region CNV /
mappability spread, fitted from the pool) mixed by the per-region Poisson sampling. The observed ``C = g + r``
with a **flat one-sided RNA prior** ``r ≥ 0`` (RNA carries no informative density prior, and this dodges any
cutoff), so the posterior on the gDNA fraction is ``P(g | C) ∝ NegBinom(g; ρ_bg·E_g, α_eff)·1[g ≤ C]``. On the
``f_g`` solve grid this is the factor ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` (the truncation is automatic
since ``f_g ≤ 1 ⇒ g ≤ C``). It deconvolves **gDNA only** (gDNA is strand-symmetric, so no strand is assigned; the
residual RNA's strand is the tilt ``θ``, left to the strand solver — the SYNERGY of the two deconvolutions).

**The background rate is a POSTERIOR, not an MLE with a fallback** (owner, 2026-08-18): with a Jeffreys
prior — the same ``c = ½`` convention ψ's reference already uses — the pooled Poisson observation gives

    rho_bg  ~  Gamma(Σg + ½, ΣE)        location (Σg+½)/ΣE ;  posterior shape (``size``) Σg + ½

ONE formula for every ``Σg``, smooth through zero: a first observed fragment moves the location by at most
3×, and the ``Σg ≫ 1`` limit is the pooled MLE ``Σg/ΣE`` exactly. Precision stays honest and
count-over-length with NO tuned constant: ``Var(g) = μ + μ²/α_eff`` with ``1/α_eff = 1/α + 1/size`` — the
per-region Poisson (``μ``), the fitted per-region over-dispersion (``α``), and the posterior's own width
(``size``). At ``Σg = 0`` the size is ½, so the factory says "around ½/ΣE, and I genuinely do not know" and
the region's own count decides.

⛔ **This replaces a branch that detonated** (2026-08-18): at ``Σg = 0`` the location used to fall back to a
"resolution wall" ``mean(1/E)`` — a mean of reciprocals, owned by the smallest regions of the partition
(TRAPS: a-mean-of-ratios-inherits-the-partition). One fragment-length intergenic sliver put the wall at
0.2985/bp on a library whose true background was exactly 0, at claimed precision ``Σg + n0 ≈ 1,300`` — an
empty region counted as a full unit of Fisher information — and the factory confidently called **80 % of all
nascent intron mass gDNA**. The falsification set is ``test_density_deconv``'s sliver-invariance /
smooth-through-zero / empty-pool-honesty gates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from .gdna_density import pooled_log_rate
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
    """The fitted gDNA density posterior the density deconvolution scores each region against.

    ``log_mu_bg`` is the Gamma-posterior location ``ln((Σg+½)/ΣE)`` — one formula at every ``Σg``, no
    branch. ``alpha`` is the fitted per-region over-dispersion (``+inf`` = Poisson-clean). ``size`` is the
    posterior shape ``Σg + ½`` — the information the pooled mean actually carries, in observed fragments
    (an EMPTY region is not a unit of information; counting the zero-count tally here is the 2026-08-18
    detonation, § module docstring). ``informative`` is False only when the pool has NO support at all
    (``ΣE = 0`` — a toy chain without intergenic regions) ⇒ the factor is flat."""

    log_mu_bg: float
    alpha: float
    size: float
    n_regions: int
    informative: bool


#: Jeffreys shape for a Poisson rate — the SAME ½ convention ψ's reference strength already uses.
_JEFFREYS_SHAPE = 0.5


def fit_gdna_background(g_counts, eff_g) -> GdnaBackground:
    """The GENERIC fit: the gDNA background posterior from a pool of **pure-gDNA** regions.

    ``g_counts``/``eff_g`` are the pooled regions' gDNA counts and gDNA effective lengths. The rate posterior
    is conjugate — ``rho_bg ~ Gamma(Σg + ½, ΣE)`` (§ module docstring) — so location and size are each one
    line and hold at every ``Σg``, including zero. ``α`` by method-of-moments against that location: from
    ``E[Σ(g_i − μ_i)²] = Σμ_i + (1/α)·Σμ_i²``, ``α = Σμ_i² / max(Σ(g_i − μ_i)² − Σμ_i, 0⁺)`` — ``+inf``
    (Poisson) when the pool is not over-dispersed. No tuned constant. Callers select the pool
    (introns → intergenic; post-pass-0 exons → clean-gDNA regions)."""
    g = np.asarray(g_counts, dtype=np.float64)
    E = np.asarray(eff_g, dtype=np.float64)
    sg = float(g.sum())
    se = float(E.sum())
    informative = se > _EPS
    # ⭐ The location comes from `gdna_density.pooled_log_rate` — the SAME formula, moved so that the two
    # estimators of the gDNA background rate live in one module (the other being the contamination-robust
    # `one_sided_rate` the fragment-length model uses). ⚠ This is a byte-identical refactor and nothing
    # about this fit changed; ⛔ but note that `pooled_log_rate` assumes the pool it is handed is PURE,
    # which this function's own pool is not — measured 53 % mature RNA at `g05` once unannotated
    # transcription exists (`TRAPS: purity-is-a-property-of-the-annotation`). Adopting the one-sided rate
    # here is a SEPARATE change to the composition reference and must be priced on its own.
    log_mu_bg = pooled_log_rate(g, E, shape=_JEFFREYS_SHAPE) if informative else -np.inf

    alpha = np.inf
    if informative and sg > _EPS:
        mu = np.exp(log_mu_bg) * E
        excess = float(np.sum((g - mu) ** 2)) - float(np.sum(mu))  # beyond-Poisson spread
        sum_mu2 = float(np.sum(mu * mu))
        if excess > _EPS and sum_mu2 > _EPS:
            alpha = sum_mu2 / excess  # finite ⇒ over-dispersed; else Poisson (α=∞)

    return GdnaBackground(
        log_mu_bg=log_mu_bg,
        alpha=float(alpha),
        size=sg + _JEFFREYS_SHAPE,
        n_regions=int(g.shape[0]),
        informative=informative,
    )


def fit_intron_background(
    substrate,
    region_arrays,
    region_eff_g,
    *,
    include_introns: bool = False,
    counts_exposure: "tuple[np.ndarray, np.ndarray] | None" = None,
) -> GdnaBackground:
    """The INTRON special case of :func:`fit_gdna_background`: the gDNA prior IS the intergenic region
    distribution (introns are off-target at the same capture depletion as intergenic —


    ``include_introns=False`` (default) pools **intergenic only** — the clean, non-circular reference (introns
    are what we deconvolve). The real-data path may add RNA-free introns for resolution
    traded for a little circularity. Delegates the fit to
    :func:`fit_gdna_background`; the whole of this function's own job is selecting the pool.

    ``region_eff_g`` is the gDNA CONTAINED effective length per region
    (:func:`effective_length.contained_eff_length`) — the support the pooled counts are a rate over.

    ⭐⭐ ``counts_exposure`` REPLACES the ``(contained count, contained effective length)`` pair with any
    other per-region ``(counts, exposure)`` — today the MEASURED TOTAL's side-selected START/END banks
    over the region's own length (`total_abundance.region_counts_and_exposure`), selected by
    ``CalibrationConfig.background_abundance``. ⛔ **It changes the PAIR, never the estimator**: the
    pool predicate, the Gamma conjugacy and the method-of-moments ``α`` are all untouched, and a zero
    exposure drops a region from the pool exactly as a zero effective length already does — which is
    how a double-walled region excludes itself without a second predicate. ⭐ The pmf leaves the divisor
    entirely: ``E[S] = ρ·ℓ`` for every fragment length, where ``E[count] = ρ·E_contained`` needs the
    length model that capture distorts (measured 1.8–4.3× under capture; `config.py` carries the
    table)."""
    sig = np.asarray(region_arrays.signature)
    if counts_exposure is None:
        eff = np.asarray(region_eff_g, dtype=np.float64)
        # ⚠ GENOME-strand columns, summed: gDNA is strand-symmetric, so the background is a total rate.
        counts = np.asarray(substrate.region_contained.count, dtype=np.float64).sum(axis=1)
    else:
        counts, eff = (np.asarray(a, dtype=np.float64) for a in counts_exposure)
    ctype = coarse_type_array(sig)  # 0 intergenic / 1 intron / 2 exon
    # the pure-gDNA pool: intergenic only, or non-exonic when ``include_introns``.
    pool = ((ctype != 2) if include_introns else (ctype == 0)) & (eff > _EPS)
    return fit_gdna_background(counts[pool], eff[pool])


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
    """The per-region factor on ``λ`` over the ``f_g = σ(λ)`` solve grid: ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)``.

    ``count`` (``C``), ``eff_g`` (``E_g``) are per-region arrays ``(n,)``; ``fg_grid`` is the ``(K,)`` grid.
    Returns ``(n, K)``, peaked at ``f_g = ρ_bg/ρ_obs`` (the confident gDNA deconvolve), with the honest
    ``μ + μ²/α_eff`` curvature. Each row is offset so its max is 0 (an ``f_g``-independent constant is
    irrelevant to ψ; this only keeps the numbers well-scaled). A non-informative background returns all-zero."""
    C = np.asarray(count, dtype=np.float64)
    Eg = np.maximum(np.asarray(eff_g, dtype=np.float64), _EPS)
    fg = np.clip(np.asarray(fg_grid, dtype=np.float64), _EPS, 1.0 - _EPS)
    n, K = C.shape[0], fg.shape[0]
    if not background.informative:
        return np.zeros((n, K), dtype=np.float64)

    # effective size: 1/α_eff = 1/α + 1/size — per-region over-dispersion ⊕ the posterior's own width.
    # ``size ≥ ½`` whenever informative, so an EMPTY pool yields a wide factor, never a confident one.
    inv_alpha = 0.0 if not np.isfinite(background.alpha) else 1.0 / max(background.alpha, _EPS)
    inv_alpha += 1.0 / max(background.size, _EPS)
    alpha_eff = np.inf if inv_alpha <= _EPS else 1.0 / inv_alpha

    mu = np.exp(background.log_mu_bg) * Eg  # (n,) background gDNA count location
    g = fg[None, :] * C[:, None]  # (n, K) gDNA count at each grid point
    logf = _log_negbinom(g, mu[:, None], alpha_eff)  # (n, K)
    logf = logf - np.max(logf, axis=1, keepdims=True)  # per-region offset (ψ-irrelevant constant)
    return logf


def density_factor_precision(lam_logprior, lam_grid):
    """``I_density`` — the composition evidence a density-deconvolution ``λ``-factor carries, read off its own
    CURVATURE.

    The factor :func:`density_lambda_factor` is a genuine, reference-FREE likelihood on ``λ`` (external
    ``ρ_bg`` information about this region's composition), so its precision belongs on the composition ``λ``-axis
    alongside the strand evidence. ``τ_λ = 1/Var_λ`` under the normalized factor; the NegBinom
    ``Var(g) = μ + μ²/α_eff`` makes a low-count / high-overdispersion region deconvolve IMPRECISELY and a dense one
    confidently, so it self-limits on thin data. A FLAT row (a non-deconvolved region, or an uninformative
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
