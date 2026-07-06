# Fix 3 — single-strand log-odds grid quantization

**Status:** implemented 2026-07-06 (full-suite A/B pending at time of writing). Scenario:
`gdna300 / ss0.99 / capture-on`. Follows Fix 1 (mixture bridge). Companion:
`boundary_kde_valley_collapse_and_simplex_precision.md` (Fix 1/2), `log_density_1d_solver_design.md`.

## Summary

After Fix 1, the #1 remaining calibration error is **not** the AMBIG class it fixed — it is **grid
quantization of the single-strand gDNA fraction**. The per-node solve reads `f_g` as the posterior **median
snapped to the σ(λ) grid**; a deep-count single-strand exon has a **sharp** posterior that concentrates at
the nearest of the `n_grid=60` grid points, so `f_g` is quantized to the lattice {0.300, 0.376, 0.458, 0.542,
0.624} (Δf_g ≈ 0.085 near 0.5). The true off-grid `f_g` differs by up to ±0.04, and × high exon mass this
was **61 % of the remaining error** (top-40: 18,920 single-strand vs 2,180 AMBIG; 36/40 pinned at the grid).

**Fix:** solve single-strand nodes on a **finer 1-D grid** (`sweep_n_grid_single_strand=256`), decoupled from
the AMBIG 2-D `(λ,τ)` cube (which stays `sweep_n_grid=60` — a fine grid there is `O(K·K_t)` and a genome-scale
memory risk). Keep the **snap median** estimator. Target result: calibration Σ|err| **30,977 → 16,842
(−45.6 %)**; the pure-gDNA / precision unit tests still pass; goldens regenerated (tiny shift, max rel 2.9e-6).

## Estimator study (why snap-median, not interpolation / mean / mode)

The instinct was to *interpolate* the median instead of finer grid. The study (target condition, Σ|err|,
`RIGEL_MEDIAN_MODE` × `RIGEL_N_GRID`):

| grid | snap | cdf-interp | mean | parabola(mode) |
|---|---|---|---|---|
| 60  | 30,977 | 59,920 | 28,396 | 20,061 |
| 256 | 19,478 | 26,249 | 18,897 | 15,242 |
| 512 | 18,150 | 20,699 | 19,192 | 15,087 |

Lessons:

* **CDF-interpolation is WRONG for a sharp posterior** (worst at every grid). It models the density as
  piecewise-uniform (linear CDF between grid points), so for a near-point-mass at grid `k` it places the
  median at the *middle of the CDF jump* ≈ half a step **below** the mass → a large systematic under-call
  (net −55k at n_grid=60). Finer grid shrinks but never removes the bias.
* **`mean`** (E[f_g]) is continuous but only mildly better than snap: it is **not transform-invariant**
  (E[σ(λ)] ≠ σ(E[λ])), so the skew of the f_g posterior biases it.
* **`parabola`** (sub-grid MODE via a 3-point log-posterior fit) de-quantizes best on *symmetric* posteriors
  (parabola@60 ≈ snap@256, ~free), BUT it silently changes the estimator from **median → mode**. For a
  **skewed** posterior — a confident pure-gDNA node whose posterior piles up near the f_g→1 vertex — the mode
  is below the median, so it **under-calls** (a pure-gDNA 50/50-count node reads ~0.75 instead of >0.8; unit
  test `test_pure_gdna_node_confident_at_near_binomial_od`). The median was chosen deliberately
  (transform-invariance + robustness to the f_g skew), so replacing it with the mode is a regression at the
  vertices even though it wins on the aggregate.
* **`snap` on a finer grid** keeps the median (correct everywhere, passes the vertex unit tests) and removes
  the quantization directly. So Fix 3 = **finer grid + snap**, NOT a new estimator. (The `RIGEL_MEDIAN_MODE`
  env is retained for study only; production default is `snap`.)

The finer grid also exposes a small **secondary** effect: the coarse snap uses the upper CDF bracket
(`fg[idx]`), biasing f_g slightly *high* (net +2,030 at n_grid=60); at n_grid=256 that bias is gone and a
small residual strand-posterior under-call surfaces (net −5,848). This is a separate, smaller lever (the
snap direction / a strand-likelihood centring question), not addressed here.

## Why the 1-D path only

Single-strand solve is `O(m·K)` (a 1-D λ grid) — quadrupling K is cheap. The AMBIG solve is the 2-D `(λ,τ)`
cube `O(m·K·K_t)`; at K=256 a genome-scale run would materialise ~2 GB per batch (18× the K=60 cube), a real
memory risk. AMBIG also carries little of the residual (2,180 of 31,020) and has its own issues (the κ-driven
strand-attribution flip on region 236 — the next dissection target). So the grids are decoupled:
`_solve_nodes_logodds_all` runs single-strand on `n_grid_ss` and AMBIG on `n_grid`, regridding the
coarse-grid global prior onto the fine grid for the single-strand solve (`_regrid_global`; the global is
smooth in f_g so linear interpolation is exact).

## Config / CLI

`CalibrationConfig.sweep_n_grid_single_strand` (default 256), CLI `--sweep-n-grid-single-strand`, `--config`
YAML. 256 is where the pair saturates (512 adds <1 %). Decoupled from `sweep_n_grid` (AMBIG, default 60).

## Open / next

* Full-suite A/B: confirm mature-transcript accuracy stays flat and no FP regression (Fix 3's effect is at
  per-node calibration; the gene-pool deliverable moves little).
* The AMBIG residual (region 236 strand-flip; AMBIG grid quantization) — the next dissection.
* The secondary snap-direction / strand-centring bias (net −5.8k at fine grid).
