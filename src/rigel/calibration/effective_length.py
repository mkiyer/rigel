"""gDNA fragment-length effective lengths (PR 4b §I.1).

The gDNA count exposure of a region or boundary — the genomic measure of fragment
start-positions that produce the event — is a **pure geometric function of a
fragment-length (FL) pmf**, identical for any species; only the FL distribution
differs (the PR04b fairness note). Calibration applies these with the **gDNA** FL
because gDNA is the density × exposure component (`ρ_0`, `ω`); RNA is the observed
residual mass (doc 01 §9), so it consumes no modelled length.

* **Region** of physical length ``L`` — a fragment is *contained* iff it fits, so
  the exposure is the FL-corrected ``E_f[max(0, L − ℓ)] = Σ_{ℓ≤L} (L − ℓ) f(ℓ)``.
  Region size limits which fragments fit (a contained fragment is shorter than
  the region); `→ L` for `L ≫ μ_FL`, small for `L ≲ μ_FL`.
* **Boundary** — a fragment *crosses* a point iff its start lies in the ``ℓ`` bp
  upstream, **independent of region sizes** (a longer fragment simply spills onto
  further regions and still crosses), so the exposure is just ``μ_FL = Σ_ℓ ℓ
  f(ℓ)`` — the mean gDNA fragment length.

Region size constrains only the fractional *mass* the accumulator splits per side
(overlap in bases) — a separate quantity from this count exposure; the molecule is
still drawn from the unconstrained gDNA FL.
"""

from __future__ import annotations

import numpy as np

__all__ = ["fl_mean", "region_eff_length", "boundary_eff_length", "boundary_side_eff_length"]


def _as_pmf(fl_pmf: np.ndarray) -> np.ndarray:
    """Return the FL pmf as a 1-D float64 array (normalised if it has mass)."""
    p = np.asarray(fl_pmf, dtype=np.float64)
    if p.ndim != 1 or p.shape[0] == 0:
        raise ValueError("fl_pmf must be a non-empty 1-D array indexed by fragment length.")
    total = float(p.sum())
    return p / total if total > 0.0 else p


def fl_mean(fl_pmf: np.ndarray) -> float:
    """``μ_FL = Σ_ℓ ℓ·f(ℓ)`` — the gDNA boundary crossing exposure (region-free)."""
    p = _as_pmf(fl_pmf)
    ell = np.arange(p.shape[0], dtype=np.float64)
    return float(np.dot(ell, p))


def boundary_eff_length(fl_pmf: np.ndarray) -> float:
    """gDNA boundary crossing effective length = ``μ_FL`` (independent of regions)."""
    return fl_mean(fl_pmf)


def region_eff_length(region_len_bp: np.ndarray, fl_pmf: np.ndarray) -> np.ndarray:
    """``E_f[max(0, L − ℓ)]`` per region, vectorised over region lengths.

    Uses ``Σ_{ℓ≤L}(L − ℓ) f(ℓ) = L·F(L) − S(L)`` with cumulative sums
    ``F(L)=Σ_{ℓ≤L} f(ℓ)`` and ``S(L)=Σ_{ℓ≤L} ℓ f(ℓ)``. For ``L`` beyond the pmf
    support the full sums apply (`F=1`, `S=μ_FL`) so the result is ``L − μ_FL``.
    """
    p = _as_pmf(fl_pmf)
    n = p.shape[0]
    ell = np.arange(n, dtype=np.float64)
    cum_f = np.cumsum(p)  # F(ℓ)
    cum_lf = np.cumsum(ell * p)  # S(ℓ)

    L = np.asarray(region_len_bp, dtype=np.float64)
    idx = np.clip(np.floor(L).astype(np.int64), 0, n - 1)
    eff = L * cum_f[idx] - cum_lf[idx]
    return np.maximum(eff, 0.0)


def boundary_side_eff_length(fl_pmf: np.ndarray, region_side_len_bp: np.ndarray) -> np.ndarray:
    """Per-side boundary **density** effective length ``E_f[min(ℓ, R_side)]`` (R2/R3, Phase 3.1).

    A crossing fragment of length ``ℓ`` contributes at most ``min(ℓ, R_side)`` bases of mass
    to the side bounded by its region of length ``R_side`` — so the length the fractional
    crossing mass is divided by, integrated over the FL pmf, is

        ``E_f[min(ℓ, R)] = Σ_{ℓ≤R} ℓ f(ℓ)  +  R·Σ_{ℓ>R} f(ℓ) = S(R) + R·(1 − F(R))``.

    For ``R ≫ support`` this → ``μ_FL`` (the region never binds); for ``R`` small it → ``R``
    (every crossing fragment spills past the short side). This is the **density** length —
    distinct from the region-free **count exposure** ``μ_FL`` (statistical power), which keeps
    the full FL because long fragments still cross (doc above; review R5).
    """
    p = _as_pmf(fl_pmf)
    n = p.shape[0]
    ell = np.arange(n, dtype=np.float64)
    cum_f = np.cumsum(p)  # F(ℓ)
    cum_lf = np.cumsum(ell * p)  # S(ℓ)

    R = np.asarray(region_side_len_bp, dtype=np.float64)
    idx = np.clip(np.floor(R).astype(np.int64), 0, n - 1)
    eff = cum_lf[idx] + R * (1.0 - cum_f[idx])
    return np.maximum(eff, 0.0)
