"""Generic density deconvolution — deconvolve a region's counts into gDNA + RNA against a gDNA prior.

A region's unspliced count ``C`` (over its gDNA effective length ``E_g``) is ``gDNA + RNA``. Given a gDNA
density prior — how dense gDNA is at this region, as a distribution ``π_bg`` — the gDNA part is
``g ≈ ρ_bg·E_g`` and the residual reads as RNA, with an honest, count-derived precision. This is the
generic count-deconvolution primitive; the intron factory (:func:`fit_intron_background`) is its special
case, in which the gDNA prior is the intergenic region distribution, introns being off-target at the same
capture depletion as intergenic.

The model: the gDNA count ``g ~ NegBinom(mean = ρ_bg·E_g, size = α_eff)``, a Gamma-Poisson — the
per-region background rate ``ρ ~ Gamma`` (over-dispersion ``α`` = per-region CNV / mappability spread,
fitted from the pool) mixed by the per-region Poisson sampling. The observed ``C = g + r`` carries a flat
one-sided RNA prior ``r ≥ 0`` (RNA has no informative density prior, and flat dodges any cutoff), so the
posterior on the gDNA fraction is ``P(g | C) ∝ NegBinom(g; ρ_bg·E_g, α_eff)·1[g ≤ C]``. On the ``f_g``
solve grid that is the factor ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)``, the truncation automatic since
``f_g ≤ 1 ⇒ g ≤ C``. The factor's rows are built per block inside the solve's kernel from the background
fitted here and the intron's count and opportunity (`native/solve_kernel.cpp`, ``factory_row``; the inputs
are `calibrate.FactoryRows`), and their curvature — the composition evidence a factor carries, ``1/Var_λ``
under the normalised row — is read there too (``factor_precision_row``); both are bound for the gates as
`native.transfer_rows.factory_rows` / `factor_precision`. It deconvolves gDNA only: gDNA is strand-symmetric, so no strand is assigned here and
the residual RNA's strand is the tilt ``θ``, left to the strand solver — the two deconvolutions are
complementary.

The background rate is a posterior, not an MLE with a fallback. Under a Jeffreys prior — the same ``c = ½``
convention ψ's reference uses — the pooled Poisson observation gives

    rho_bg  ~  Gamma(Σg + ½, ΣE)        location (Σg+½)/ΣE ;  posterior shape (``size``) Σg + ½

one formula for every ``Σg``, smooth through zero: a first observed fragment moves the location by at most
3×, and the ``Σg ≫ 1`` limit is the pooled MLE ``Σg/ΣE`` exactly. Precision stays honest and
count-over-length with no tuned constant: ``Var(g) = μ + μ²/α_eff`` with ``1/α_eff = 1/α + 1/size`` — the
per-region Poisson (``μ``), the fitted per-region over-dispersion (``α``), and the posterior's own width
(``size``). At ``Σg = 0`` the size is ½, so the factory says "around ½/ΣE, and I genuinely do not know" and
the region's own count decides.

⛔ There must be no fallback branch at ``Σg = 0``. A fallback location built from a mean of reciprocals is
owned by the smallest regions of the partition, so a
single fragment-length sliver can assert a large background on a library whose true background is exactly
zero, at a claimed precision that counts every empty region as a unit of Fisher information. The
falsification set is ``test_density_deconv``'s sliver-invariance / smooth-through-zero /
empty-pool-honesty gates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .gdna_density import pooled_log_rate
from .signature import coarse_type_array

__all__ = [
    "GdnaBackground",
    "fit_gdna_background",
    "fit_intron_background",
]

_EPS = 1.0e-12


@dataclass(frozen=True, slots=True)
class GdnaBackground:
    """The fitted gDNA density posterior the density deconvolution scores each region against.

    ``log_mu_bg`` is the Gamma-posterior location ``ln((Σg+½)/ΣE)`` — one formula at every ``Σg``, no
    branch. ``alpha`` is the fitted per-region over-dispersion (``+inf`` = Poisson-clean). ``size`` is the
    posterior shape ``Σg + ½`` — the information the pooled mean actually carries, in observed fragments.
    An empty region is not a unit of information, so the zero-count tally must never be added here.
    ``informative`` is False only when the pool has no support at all (``ΣE = 0`` — a toy chain without
    intergenic regions), and the factor is then flat."""

    log_mu_bg: float
    alpha: float
    size: float
    n_regions: int
    informative: bool


#: Jeffreys shape for a Poisson rate — the same ½ convention ψ's reference strength uses.
_JEFFREYS_SHAPE = 0.5


def fit_gdna_background(g_counts, eff_g) -> GdnaBackground:
    """The generic fit: the gDNA background posterior from a pool of pure-gDNA regions.

    ``g_counts``/``eff_g`` are the pooled regions' gDNA counts and gDNA effective lengths. The rate posterior
    is conjugate — ``rho_bg ~ Gamma(Σg + ½, ΣE)`` — so location and size are each one
    line and hold at every ``Σg``, including zero. ``α`` by method-of-moments against that location: from
    ``E[Σ(g_i − μ_i)²] = Σμ_i + (1/α)·Σμ_i²``, ``α = Σμ_i² / max(Σ(g_i − μ_i)² − Σμ_i, 0⁺)`` — ``+inf``
    (Poisson) when the pool is not over-dispersed. No tuned constant. Callers select the pool
    (introns → intergenic; post-pass-0 exons → clean-gDNA regions)."""
    g = np.asarray(g_counts, dtype=np.float64)
    E = np.asarray(eff_g, dtype=np.float64)
    sg = float(g.sum())
    se = float(E.sum())
    informative = se > _EPS
    # The location comes from `gdna_density.pooled_log_rate`, so that both estimators of the gDNA
    # background rate live in one module (the other being the contamination-robust `one_sided_rate` the
    # fragment-length model uses). Note that `pooled_log_rate` assumes the pool it is handed is pure,
    # which this function's pool is not once unannotated transcription exists. Adopting the one-sided
    # rate here would be a separate change to the composition reference and must be priced on its own.
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
) -> GdnaBackground:
    """The intron special case of :func:`fit_gdna_background`: the gDNA prior is the intergenic region
    distribution, introns being off-target at the same capture depletion as intergenic.

    ``include_introns=False`` (default) pools intergenic only — the clean, non-circular reference, since
    introns are what is being deconvolved. The real-data path may add RNA-free introns, buying resolution
    at the price of a little circularity. The fit itself is delegated to :func:`fit_gdna_background`; the
    whole of this function's own job is selecting the pool.

    ``region_eff_g`` is the gDNA contained effective length per region
    (:func:`effective_length.contained_eff_length`) — the support the pooled counts are a rate over."""
    sig = np.asarray(region_arrays.signature)
    eff = np.asarray(region_eff_g, dtype=np.float64)
    # Genome-strand columns, summed: gDNA is strand-symmetric, so the background is a total rate.
    counts = np.asarray(substrate.region_contained.count, dtype=np.float64).sum(axis=1)
    ctype = coarse_type_array(sig)  # 0 intergenic / 1 intron / 2 exon
    # the pure-gDNA pool: intergenic only, or non-exonic when ``include_introns``.
    pool = ((ctype != 2) if include_introns else (ctype == 0)) & (eff > _EPS)
    return fit_gdna_background(counts[pool], eff[pool])
