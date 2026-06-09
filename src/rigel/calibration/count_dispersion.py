"""Count overdispersion — the count-side twin of the strand Beta-Binomial overdispersion.

The count-prior concentration must be the count clue's *honest* precision. The raw expected gDNA
count ``density·eff_len`` (thousands) is the **Poisson** precision; RNA-seq / gDNA counts are
**NB-overdispersed**, so the honest precision is the overdispersion-limited **effective count**::

    N_eff = N / (1 + α·N)        (→ 1/α for large N — a ceiling on how much any count is trusted)

``α`` is the NB count dispersion (``Var = μ + α·μ²``), fit by a **common-dispersion method of
moments** from the count-observable gDNA seeds — the same seeds the density estimator and the gDNA
strand fit trust. Two count *types* are fit separately because they carry different heterogeneity:

* **contained** — count-observable regions (intron / intergenic). Uniformly depleted under capture,
  so low dispersion even there.
* **crossing** — count-observable boundary sides (exon–intron / exon–intergenic seams). These are
  the imputation source for exonic regions and carry the capture enrichment heterogeneity, so their
  dispersion is large under capture (``N_eff → 1/α`` ⇒ imputed exons fade and the strand clue
  governs — the leak-closing behaviour).

Each seed's count ``N`` is measured against the expectation ``μ = ρ̄·eff_len`` under the **global
pooled gDNA density** ``ρ̄`` (not the seed's own count → no circularity; the analog of DESeq2's
mean-trend fit). The per-type MoM is shrunk toward the **global pooled-seed trend** (all seeds
pooled) by seed count — abundant seeds follow their own type, sparse ones borrow from the pool —
with the conservative geometric ``α₀`` (``N_eff → 1``) used only when the pool itself is degenerate.
Mirrors :mod:`gdna_strand` (pooled MoM + prior shrinkage in seed-node units).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_ALPHA_FLOOR: float = 0.0  # an NB dispersion cannot be negative


@dataclass(frozen=True, slots=True)
class CountDispersionModel:
    """Fitted gDNA count overdispersion (one per count-type) + fit provenance."""

    alpha_contained: float  # NB dispersion of the count-observable contained-region counts
    alpha_crossing: float  # NB dispersion of the count-observable crossing boundary-side counts
    alpha_pooled: float  # the global pooled-seed trend α the per-type estimates were shrunk toward
    n_contained_seeds: int
    n_crossing_seeds: int
    fallback_used: bool  # pooled trend undefined (no usable seeds) ⇒ shrank toward the prior α₀


def effective_count(count: np.ndarray, alpha: float) -> np.ndarray:
    """Overdispersion-limited effective count ``N / (1 + α·N)`` (→ ``1/α`` as ``N → ∞``)."""
    count = np.asarray(count, dtype=np.float64)
    a = max(float(alpha), _ALPHA_FLOOR)
    return count / (1.0 + a * count)


def _nb_mom(count: np.ndarray, mu: np.ndarray) -> float | None:
    """Common-dispersion NB method of moments ``α̂ = Σ[(N−μ)² − μ] / Σ μ²``.

    ``Var = μ + α·μ²`` ⇒ the Poisson term ``μ`` is subtracted per seed before pooling. Returns
    ``None`` when the denominator is non-positive (no usable seeds) — the estimate is undefined.
    """
    count = np.asarray(count, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    valid = mu > 0.0
    if not np.any(valid):
        return None
    num = float(np.sum((count[valid] - mu[valid]) ** 2 - mu[valid]))
    den = float(np.sum(mu[valid] ** 2))
    if den <= 0.0 or not np.isfinite(num):
        return None
    return max(num / den, _ALPHA_FLOOR)


def _shrink(a_hat: float | None, n: int, target: float, prior_weight: float) -> float:
    """Precision-weighted shrinkage of ``a_hat`` toward ``target`` by seed count ``n``.

    ``(n·a_hat + w·target)/(n + w)``; an undefined per-type estimate (``None``) takes the target.
    """
    if a_hat is None:
        return max(float(target), _ALPHA_FLOOR)
    w = max(float(prior_weight), 0.0)
    total = n + w
    blended = (n * a_hat + w * float(target)) / total if total > 0.0 else a_hat
    return max(blended, _ALPHA_FLOOR)


def fit_gdna_count_overdispersion(
    contained_count: np.ndarray,
    contained_eff_len: np.ndarray,
    crossing_count: np.ndarray,
    crossing_eff_len: np.ndarray,
    *,
    prior_alpha: float,
    prior_weight: float,
) -> CountDispersionModel:
    """Fit the gDNA count overdispersion per count-type from the calibration seeds.

    Parameters
    ----------
    contained_count, contained_eff_len : np.ndarray
        Per-seed clean gDNA count ``N`` and effective length ``L`` for the count-observable
        contained regions (intron / intergenic).
    crossing_count, crossing_eff_len : np.ndarray
        The same for the count-observable boundary sides (exon-adjacent seams); ``L`` is the
        per-side gDNA-FL crossing mean (``fl_mean``).
    prior_alpha : float
        Conservative fallback ``α₀`` (geometric ``1.0`` ⇒ ``N_eff → 1``) used only when the pooled
        trend is degenerate.
    prior_weight : float
        Shrinkage strength in seed-node units (toward the pooled trend).
    """
    cc = np.asarray(contained_count, dtype=np.float64)
    cl = np.asarray(contained_eff_len, dtype=np.float64)
    xc = np.asarray(crossing_count, dtype=np.float64)
    xl = np.asarray(crossing_eff_len, dtype=np.float64)

    all_n = np.concatenate([cc, xc])
    all_l = np.concatenate([cl, xl])
    # Global pooled gDNA density trend ρ̄ (frags per effective bp) — the expectation each seed's
    # count is measured against, so the excess variance is around the GLOBAL density (no circularity
    # with the seed's own count). The analog of DESeq2's fitted mean-trend.
    tot_l = float(np.sum(all_l))
    rho = float(np.sum(all_n)) / tot_l if tot_l > 0.0 else 0.0

    pooled = _nb_mom(all_n, rho * all_l)  # the global trend (None ⇒ degenerate)
    target = pooled if pooled is not None else float(prior_alpha)
    a_contained = _shrink(_nb_mom(cc, rho * cl), cc.size, target, prior_weight)
    a_crossing = _shrink(_nb_mom(xc, rho * xl), xc.size, target, prior_weight)
    return CountDispersionModel(
        alpha_contained=a_contained,
        alpha_crossing=a_crossing,
        alpha_pooled=max(float(target), _ALPHA_FLOOR),
        n_contained_seeds=int(cc.size),
        n_crossing_seeds=int(xc.size),
        fallback_used=pooled is None,
    )


__all__ = [
    "CountDispersionModel",
    "effective_count",
    "fit_gdna_count_overdispersion",
]
