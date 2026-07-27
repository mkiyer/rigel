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
    """``E_f[max(0, L − ℓ + 1)]`` per region — the DISCRETE count of start positions at which a
    length-``ℓ`` fragment lies wholly inside a length-``L`` region.

    **The ``+1`` is the discrete count, not a correction factor.** A fragment ``[s, s+ℓ)`` is contained in
    ``[a, a+L)`` iff ``a ≤ s`` and ``s+ℓ ≤ a+L``, i.e. ``s ∈ [a, a+L−ℓ]`` — that is ``L − ℓ + 1`` positions,
    verified by brute-force enumeration (L=100,ℓ=100 → 1; L=150,ℓ=100 → 51; L=500,ℓ=100 → 401). It matches
    the accumulator's own rule, which is discrete: a fragment with all bp in one region increments that
    region by exactly ``+1`` (`accumulator/00_design.md` §4.1-4.2). The divisor of an integer count must be
    the integer opportunity count.

    ⚠ **This previously read ``E_f[max(0, L − ℓ)]``, which is 0 whenever ``L`` equals the shortest fragment**
    — not merely inaccurate but a DIVISION BY ZERO, floored to ``_EPS`` by the callers and yielding densities
    of ~1e9. Measured on `ambig_dense_10mb`: 211 / 1698 regions (12.4 %) had ``eff ≤ 1e-9``; one 100 bp region
    with 2 fragments (gDNA FL concentrated at ℓ=100) produced ρ = 2e9 and consumed 84 % of the NPMLE prior's
    fitted grid. Away from ``L ≈ ℓ`` the change is negligible (L=2000 → 1900 vs 1901).

    **Short regions still relay rather than project — that is the point, and this fix is what makes it
    honest.** A short exon gets a SMALL eff-length and therefore a genuinely IMPRECISE density (few or zero
    contained fragments over little opportunity: ρ = 2/1 with ``Var(log ρ) = 1/n = 0.5``, versus ~0.001 for a
    long region), so its weak local belief defers to its neighbours' messages; at zero count it emits nothing
    at all. The old ``eff = 0`` did not make it imprecise — it made it *confidently absurd*.

    Uses ``Σ_{ℓ≤L}(L − ℓ + 1) f(ℓ) = (L+1)·F(L) − S(L)`` with ``F(L)=Σ_{ℓ≤L} f(ℓ)`` and
    ``S(L)=Σ_{ℓ≤L} ℓ f(ℓ)``. Beyond the pmf support the full sums apply (``F=1``, ``S=fl_mean``) so the
    result is ``L + 1 − fl_mean``.

    (Checked at the same time: :func:`boundary_eff_length` and
    :func:`boundary_side_crossing_count_eff_length` are EXACT against enumeration. :func:`spliced_side_eff_length`
    is NOT merely "0.5-2 % low / cosmetic" — that holds only when the mature transcript TERMINATES at the flank's
    far edge. When the exon CONTINUES past ``R`` (or that far edge is itself a junction), the correct one-sided
    opportunity is ``E[min(ℓ,R)]/2`` = :func:`boundary_side_eff_length`, and the half-triangle is **2–199× low**
    on such faces (≈26 % of junction faces on ``ambig_dense_10mb``). It is functionally **inert today** — an A/B
    through the sweep moves ``f_g`` by ≤1e-4 because the spliced channel is heavily down-weighted — but it is a
    real geometry error that becomes load-bearing once the spliced-measurement channel is un-gagged (priority #3).
    See ``docs/calibration/archive/mature_crossing_gate.md`` §6.2. **Deferred, not cosmetic.**)
    """
    p = _as_pmf(fl_pmf)
    n = p.shape[0]
    ell = np.arange(n, dtype=np.float64)
    cum_f = np.cumsum(p)  # F(ℓ)
    cum_lf = np.cumsum(ell * p)  # S(ℓ)

    L = np.asarray(region_len_bp, dtype=np.float64)
    idx = np.clip(np.floor(L).astype(np.int64), 0, n - 1)
    eff = (L + 1.0) * cum_f[idx] - cum_lf[idx]
    return np.maximum(eff, 0.0)


def boundary_side_crossing_count_eff_length(
    fl_pmf: np.ndarray, region_side_len_bp: np.ndarray
) -> np.ndarray:
    """Per-side boundary **COUNT** effective length ``E_f[min(ℓ, R_side)]`` — the expected NUMBER of
    fragments crossing the boundary whose near part lies in a flank region of length ``R_side``:

        ``E_f[min(ℓ, R)] = Σ_{ℓ≤R} ℓ f(ℓ)  +  R·Σ_{ℓ>R} f(ℓ) = S(R) + R·(1 − F(R))``.

    For ``R ≫ support`` this → ``fl_mean``; for small ``R`` → ``R``. This is a **count**, so it is the
    right scale for statistical power (Fisher information) — **not** for a density. The accumulator does
    not deposit a whole fragment onto a face; it deposits that face's *share*. Use
    :func:`boundary_side_eff_length` to divide a per-face MASS.
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


def boundary_side_eff_length(fl_pmf: np.ndarray, region_side_len_bp: np.ndarray) -> np.ndarray:
    """Per-side boundary **DENSITY** effective length ``E_f[min(ℓ, R_side)] / 2``.

    The divisor that recovers ``ρ`` from a boundary face's deposited **mass**. The factor ½ is
    **derived from the accumulator's deposit rule**, not a tuning constant
    (`accumulator/00_design.md` §4.3; reference: ``tests/native/_accumulator_reference.py``):
    a slice's deposited share is ``(slice_len/ℓ) / n_cross``, where ``n_cross`` is the number of
    boundaries that slice crosses. For a fragment straddling one boundary with near-part length ``x``:

    * ``x ≤ R`` — the near slice is an END slice (``n_cross=1``) ⇒ deposits ``x/ℓ``;
    * ``x > R`` — the fragment overruns the flank, so that slice is INTERIOR (``n_cross=2``) ⇒
      deposits ``R/(2ℓ)``.

    Integrating over the uniform start position, the per-face expected mass per unit ``ρ`` is

        ``∫₀^min(ℓ,R) (x/ℓ) dx  +  (ℓ−R)·R/(2ℓ)  =  min(ℓ, R) / 2``   — **exactly, for every R**.

    Hence ``E_face = E_f[min(ℓ,R)]/2``. Verified against the reference accumulator at fixed ``ℓ=200``:
    measured face mass 99.5 / 50.0 / 25.0 at ``R`` = 600 / 100 / 50 vs this function's 100 / 50 / 25
    (the short-``R`` cells discriminate it from every other candidate).

    **This ½ was previously missing**, so every boundary↔region message carried a systematic 2× frame
    error (a boundary told its neighbour "you are half gDNA"; a region told a boundary "you are 2×",
    clipped to 1). The retired total-density σ²_imp moment estimator measured the resulting −0.69-nat offset
    and removed it from the *variance* while leaving it in the message *mode*.

    Distinct from :func:`boundary_side_crossing_count_eff_length` (the un-halved **count** length —
    statistical power, where the whole fragment does cross) and from
    :func:`spliced_side_eff_length` (``E[min(ℓ,R)²/(2ℓ)]`` — a *one-sided* spliced deposit, where the
    far flank is never credited so the slice is always an END slice).
    """
    return 0.5 * boundary_side_crossing_count_eff_length(fl_pmf, region_side_len_bp)


def spliced_side_eff_length(fl_pmf: np.ndarray, region_side_len_bp: np.ndarray) -> np.ndarray:
    """Per-side **one-sided crossing** density effective length ``E_f[min(ℓ, R)² / (2ℓ)]``.

    A *spliced* fragment crosses a junction on **one side only** — its exon flank is credited; the intron
    flank is never touched (an intron region carries no slice). Under a uniform mature density ``ρ`` the
    one-sided deposited mass is ``Σ_{a=1}^{min(ℓ,R)} (a/ℓ) ≈ min(ℓ,R)²/(2ℓ)`` (a **half-triangle**: the exon
    coverage ``a`` ranges ``0→min(ℓ,R)`` and its mass share ``a/ℓ`` is linear in ``a``), so the divisor that
    recovers ``ρ`` is ``E_f[min(ℓ,R)²/(2ℓ)]``.

    ⚠ **This is NOT "half of** :func:`boundary_side_eff_length` **"** — that claim was true before the ½ moved
    into that function (see its body) and is now stale by exactly 2×. At ``R ≫ FL support`` **the two are
    EQUAL**: measured on a Gaussian FL (μ=200, σ=50), ``R = 400 / 1000 / 5000`` give ratio **1.0000**. They
    diverge only at ``R ≲ FL support`` (0.531 at R=100, 0.269 at R=50), where the one-sided half-triangle and
    the two-sided ``E[min(ℓ,R)]`` genuinely measure different geometry. **Keeping the two divisors separate is
    still correct — but for the short-flank reason, not a 2× one.** Anyone deriving a message variance from
    the old wording lands 2× off.

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
