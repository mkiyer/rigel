"""gDNA fragment-length effective lengths (PR 4b §I.1).

The count effective length of a region or boundary — the genomic measure of fragment
start-positions that produce the event — is a **pure geometric function of a
fragment-length (FL) pmf**, identical for any species; only the FL distribution
differs. Calibration applies these PER COMPONENT: ``bp_solver.build_node_geometry`` computes
each node's gDNA eff-lengths from the **gDNA** FL and its RNA (nascent unspliced + spliced)
eff-lengths from the **RNA** FL.

* **Region** of physical length ``L`` — a fragment is *contained* iff it fits, so
  the effective length is the FL-corrected ``E_f[max(0, L − ℓ)] = Σ_{ℓ≤L} (L − ℓ) f(ℓ)``.
  Region size limits which fragments fit (a contained fragment is shorter than
  the region); `→ L` for `L ≫ fl_mean`, small for `L ≲ fl_mean`.
* **Boundary** — a fragment *crosses* a point iff its start lies in the ``ℓ`` bp
  upstream, **independent of region sizes** (a longer fragment simply spills onto
  further regions and still crosses), so the effective length is just ``fl_mean = Σ_ℓ ℓ
  f(ℓ)`` — the mean gDNA fragment length.

Region size constrains only the fractional *mass* the accumulator splits per side
(overlap in bases) — a separate quantity from this count effective length; the molecule is
still drawn from the unconstrained gDNA FL.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "region_eff_length",
    "boundary_eff_length",
    "boundary_side_eff_length",
    "spliced_side_eff_length",
]


def _as_pmf(fl_pmf: np.ndarray) -> np.ndarray:
    """Return the FL pmf as a 1-D float64 array (normalised if it has mass)."""
    p = np.asarray(fl_pmf, dtype=np.float64)
    if p.ndim != 1 or p.shape[0] == 0:
        raise ValueError("fl_pmf must be a non-empty 1-D array indexed by fragment length.")
    total = float(p.sum())
    return p / total if total > 0.0 else p


def fl_mean(fl_pmf: np.ndarray) -> float:
    """``fl_mean = Σ_ℓ ℓ·f(ℓ)`` — the gDNA boundary crossing effective length (region-free)."""
    p = _as_pmf(fl_pmf)
    ell = np.arange(p.shape[0], dtype=np.float64)
    return float(np.dot(ell, p))


def boundary_eff_length(fl_pmf: np.ndarray) -> float:
    """gDNA boundary crossing effective length = ``fl_mean`` (independent of regions)."""
    return fl_mean(fl_pmf)


def region_eff_length(region_len_bp: np.ndarray, fl_pmf: np.ndarray) -> np.ndarray:
    """``E_f[max(0, L − ℓ)]`` per region, vectorised over region lengths.

    Uses ``Σ_{ℓ≤L}(L − ℓ) f(ℓ) = L·F(L) − S(L)`` with cumulative sums
    ``F(L)=Σ_{ℓ≤L} f(ℓ)`` and ``S(L)=Σ_{ℓ≤L} ℓ f(ℓ)``. For ``L`` beyond the pmf
    support the full sums apply (`F=1`, `S=fl_mean`) so the result is ``L − fl_mean``.
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
    """Per-side boundary **density** effective length ``E_f[min(ℓ, R_side)]``.

    A crossing fragment of length ``ℓ`` contributes at most ``min(ℓ, R_side)`` bases of mass
    to the side bounded by its region of length ``R_side`` — so the length the fractional
    crossing mass is divided by, integrated over the FL pmf, is

        ``E_f[min(ℓ, R)] = Σ_{ℓ≤R} ℓ f(ℓ)  +  R·Σ_{ℓ>R} f(ℓ) = S(R) + R·(1 − F(R))``.

    For ``R ≫ support`` this → ``fl_mean`` (the region never binds); for ``R`` small it → ``R``
    (every crossing fragment spills past the short side). This is the **density** length —
    distinct from the region-free **count effective length** ``fl_mean`` (statistical power), which
    keeps the full FL because long fragments still cross.
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


def spliced_side_eff_length(fl_pmf: np.ndarray, region_side_len_bp: np.ndarray) -> np.ndarray:
    """Per-side **one-sided crossing** density effective length ``E_f[min(ℓ, R)² / (2ℓ)]``.

    A *spliced* fragment crosses a junction on **one side only** — its exon flank is credited; the intron
    flank is never touched (an intron region carries no slice). Under a uniform mature density ``ρ`` the
    one-sided deposited mass is ``Σ_{a=1}^{min(ℓ,R)} (a/ℓ) ≈ min(ℓ,R)²/(2ℓ)`` (a **half-triangle**: the exon
    coverage ``a`` ranges ``0→min(ℓ,R)`` and its mass share ``a/ℓ`` is linear in ``a``), so the divisor that
    recovers ``ρ`` is ``E_f[min(ℓ,R)²/(2ℓ)]``. For ``R ≫ support`` this → ``fl_mean/2`` — exactly HALF the
    two-sided crossing measure ``boundary_side_eff_length`` (``E[min(ℓ,R)]``). Dividing a one-sided spliced
    mass by the two-sided measure (as a contiguous gDNA/nascent crossing is) understates the mature density
    ~2×; this is the correct one-sided divisor.

    Split ``min(ℓ,R)²/(2ℓ) = ℓ/2`` for ``ℓ ≤ R`` and ``= R²/(2ℓ)`` for ``ℓ > R``; the ``ℓ=0`` term is 0.
    """
    p = _as_pmf(fl_pmf)
    n = p.shape[0]
    ell = np.arange(n, dtype=np.float64)
    cum_lf = np.cumsum(ell * p)  # S(R) = Σ_{ℓ≤R} ℓ f(ℓ)
    inv_ell = np.zeros(n, dtype=np.float64)
    np.divide(p, ell, out=inv_ell, where=ell > 0)  # f(ℓ)/ℓ, 0 at ℓ=0
    cum_g = np.cumsum(inv_ell)  # H(R) = Σ_{1≤ℓ≤R} f(ℓ)/ℓ
    h_total = cum_g[-1]

    R = np.asarray(region_side_len_bp, dtype=np.float64)
    idx = np.clip(np.floor(R).astype(np.int64), 0, n - 1)
    eff = 0.5 * cum_lf[idx] + 0.5 * R * R * (h_total - cum_g[idx])
    return np.maximum(eff, 0.0)
