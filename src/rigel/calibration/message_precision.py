"""σ²_imp(ρ) — the honest, NON-CONSTANT adjacent-node gDNA-rate imputation variance (the message precision).

A message from node A to its neighbour B says "your gDNA rate is like mine". Its reliability is the genuine
spatial variation of the gDNA rate between adjacent loci — which is **not** a constant:

* a **depleted** pair (both at the off-target floor) agrees closely → σ² LOW → precision HIGH;
* an **enriched** pair (both on-target) varies by the enriched population's own spread → σ² MODERATE;
* a pair straddling a **probe edge** (one enriched, one depleted) genuinely disagrees → σ² HIGH → precision
  very low — the message *should* be ignored there.

Production previously used ONE scalar (`bp_solver.adjacent_disagreement_variance`) fit on the **total**
density ("RNA included" — its own docstring). Because adjacent total densities differ wildly (an expressed
exon beside an intron), that scalar is inflated: measured σ²_imp ≈ 0.96–3.54 ⇒ the precision
`pr = n_src/(n_src·σ²_imp + 1)` saturates at `1/σ²_imp` ≈ **0.28–1.05 pseudo-observations**, so no message can
ever outweigh ~1 read and RNA cannot be peeled on unstranded data (the gDNA over-call).

THE LOOP (why this is honest, not circular): a node's gDNA rate is exactly the unknown, so no *belief-free*
rate exists — the raw density is RNA-contaminated, and a prior projection collapses high-density nodes onto
the population mean (manufacturing false agreement). So:

* **pass-0** keeps the scalar total-density σ²_imp — deliberately over-stated ⇒ **gentle** messages that nudge
  rather than ruin the solve;
* **each refit** re-estimates σ²_imp(ρ) from the *solved* gDNA rates, so the messages become honest and
  stronger exactly where the data has earned it, and the solve converges.

Estimator: a **fixed-bandwidth kernel regression** of the Poisson-corrected excess disagreement on the
source's log-rate — deterministic weighted sums, NO spline / λ / GCV / argmin (the machinery whose `argmin`
bistability was the old cross-process nondeterminism root cause). Binned first, so it is O(pairs) + O(grid²).
"""

from __future__ import annotations

import numpy as np

from .bp_solver import _adjacent_edges, node_global_geometry

__all__ = ["adjacent_imputation_variance"]

_EPS = 1.0e-12


def adjacent_imputation_variance(
    chain,
    geometry,
    belief,
    *,
    fallback: float,
    bandwidth: float = 0.5,
    n_grid: int = 96,
    min_pairs: float = 32.0,
):
    """Per-node σ²_imp from the SOLVED belief's gDNA rates. Returns a ``(n_nodes,)`` array to be consumed by
    ``node_sweep(disagreement_sigma2=…)``.

    ``ρ_i = f_g,i·M_i/E_i`` is node i's solved gDNA rate. For every adjacent pair the Poisson-corrected excess
    disagreement is ``r² = (log ρ_src − log ρ_dst)² − (1/g_src + 1/g_dst)`` (the same moment correction as the
    retired scalar). ``σ²(ρ)`` is then the kernel-weighted local mean of ``r²`` against the SOURCE's log-rate —
    we cannot condition on the destination's rate, since that is precisely what the message imputes (the
    destination's own count corrects it inside the solve; proper Bayes).

    ``fallback`` (the pass-0 scalar) is used for nodes with no gDNA rate (``g=0`` ⇒ no message is emitted
    anyway) and where the local pair support is too thin (< ``min_pairs`` effective). ``bandwidth`` is in
    decades — a fixed KDE-style knob, never optimised."""
    mass, eff = node_global_geometry(chain, geometry)
    eff = np.maximum(np.asarray(eff, dtype=np.float64), _EPS)
    g = np.asarray(belief.f_g, dtype=np.float64) * np.asarray(mass, dtype=np.float64)  # gDNA count
    has_rate = g > _EPS
    lr = np.full(g.shape, np.nan, dtype=np.float64)
    lr[has_rate] = np.log(g[has_rate] / eff[has_rate])  # solved log gDNA rate

    src, dst, _ = _adjacent_edges(chain)
    ok = has_rate[src] & has_rate[dst]
    src, dst = src[ok], dst[ok]
    out = np.full(g.shape, float(fallback), dtype=np.float64)
    if src.size < min_pairs:
        return out  # too few solved pairs to learn anything — stay gentle

    x = lr[src]  # the source's log-rate (the conditioning variable)
    r2 = (lr[src] - lr[dst]) ** 2 - (1.0 / g[src] + 1.0 / g[dst])  # Poisson-corrected excess

    # ---- bin, then smooth with a FIXED kernel: O(pairs) + O(grid²), deterministic ----
    h = float(bandwidth) * np.log(10.0)
    lo, hi = float(np.min(x)), float(np.max(x))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < _EPS:
        return out
    grid = np.linspace(lo - h, hi + h, int(n_grid))
    step = grid[1] - grid[0]
    idx = np.clip(((x - grid[0]) / step).astype(np.int64), 0, int(n_grid) - 1)
    sum_r2 = np.bincount(idx, weights=r2, minlength=int(n_grid))
    cnt = np.bincount(idx, minlength=int(n_grid)).astype(np.float64)

    d = grid[:, None] - grid[None, :]
    kern = np.exp(-0.5 * (d / h) ** 2)  # (G,G) fixed Gaussian kernel
    num = kern @ sum_r2
    den = kern @ cnt
    sigma2_grid = np.where(
        den >= min_pairs, np.maximum(num / np.maximum(den, _EPS), 0.0), float(fallback)
    )

    out[has_rate] = np.interp(
        lr[has_rate], grid, sigma2_grid, left=sigma2_grid[0], right=sigma2_grid[-1]
    )
    return out
