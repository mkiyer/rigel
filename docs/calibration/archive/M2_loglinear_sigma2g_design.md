# M2 — deterministic log-linear σ²_g(μ) (reference design & implementation plan)

**Status:** design reference, 2026-07-05. The deterministic replacement for the bistable `MonotoneVarMean`
σ²_g P-spline (the calibration-nondeterminism root cause — `calibrate_cross_process_nondeterminism.md`).
Kept as a reference; **may be made moot** by the pass-1-prior-free restructure (see §7). Companion:
`background_gdna_density_prior_model.md` (why a single scalar fails — σ²_g must be density-local).

## 1. What σ²_g is and why it must be density-local

`σ²_g` is the between-node spread of log gDNA density used by `_global_logprior` as the genome-wide baseline
precision `N = 1/(var_mean + σ²_g)` (capped at one pseudo-observation). It is queried at **one point**,
`ρ_global`. It is **NOT** a single constant: measured (background_gdna_density_prior_model.md §7) the
load-bearing value is the CONDITIONAL variance at ρ_global on a monotone-increasing σ²_g(μ) curve —
`~0` at ρ_global≈0 (no-gDNA → strong prior → suppress FP), `~31` at ρ_global≈0.24 (gdna300 → weak prior →
spare enriched). Five single-scalar estimators all regress. The (retired) spline captured this monotone μ
dependence but its GCV-λ `argmin` was machine-ε bistable.

## 2. The model

Fit a **log-linear** variance–mean law and evaluate it at ρ_global:

```
σ²_g(μ)  =  a + b · log μ ,      report  σ̂²_g = max( a + b · log ρ_global , 0 )
```

Justification: the spline was ≈ linear in log μ over the observed range — points (log μ, σ²_g) ≈
(−4.3, 24.5), (−1.45, 31), (0, 36.5), (2.3, 42.9), slope ≈ 2.7/decade. A 2-parameter log-linear law captures
this with **no** smoothing hyperparameter, no `argmin`, no IRLS — hence deterministic. (If a curved trend later
proves necessary, the same closed-form machinery extends to `a + b·log μ + c·(log μ)²`; do NOT reintroduce a
spline.)

## 3. The fit — closed-form weighted least squares (deterministic)

Reuse the seed-edge sufficient statistics the spline used (`_fit_seed_varmean` inputs): for each adjacent
seed–seed edge, the Poisson-corrected excess response and the midpoint log-density predictor,

```
yₑ = (Δlog ρ)²ₑ − offsetₑ        # offsetₑ = 1/(ρ_i E_i + 1) + 1/(ρ_j E_j + 1)   (per-edge Poisson log-var)
xₑ = log μₑ ,  μₑ = ½(ρ_i + ρ_j)  # edge midpoint density
wₑ = min(seed_w_i, seed_w_j)      # the weaker endpoint limits reliability
```

Weighted OLS (closed form, `Σw`-weighted; all sums are continuous in the data ⇒ ε-robust):

```
x̄ = Σwx/Σw ,  ȳ = Σwy/Σw
b  = Σw(x−x̄)(y−ȳ) / max(Σw(x−x̄)², ε)
a  = ȳ − b·x̄
σ̂²_g = max(a + b·(log ρ_global − x̄) + b·x̄ , 0)   # = max(a + b·log ρ_global, 0)
```

**Degenerate guards (deterministic):** `Σw(x−x̄)² ≤ ε` (no x spread — e.g. all seeds at one density, as in
no-gDNA) ⇒ set `b=0`, `σ̂²_g = max(ȳ, 0)` (the flat weighted-mean excess). `< 2` edges ⇒ `σ̂²_g = 0`.

**Robustness (optional, only if the benchmark shows outlier sensitivity):** one round of fixed
bisquare/Huber down-weighting with a **data-derived** scale (the weighted MAD of residuals) — a *single*
deterministic reweight (NOT iterated to convergence, which would reintroduce nondeterminism). Start WITHOUT
it (plain WLS); add only if needed. No tuned constant beyond the standard bisquare 4.685·MAD (a fixed
statistical constant, not a fit knob) — pause and confirm before adding.

## 4. Why this is ε-robust (the whole point)

Every quantity is a continuous algebraic function of the (yₑ, xₑ, wₑ) — sums, ratios, a `max(·,0)`. A
machine-ε change in the FL pmf (→ eff-lengths → densities) moves `a`, `b`, `σ̂²_g` by machine-ε. **No discrete
selection** (the spline's `argmin` GCV / IRLS iteration count / monotone active-set are all gone). Amplification
→ ≤1× (verified requirement).

## 5. Implementation plan (surgical, in `bp_solver.py`)

1. **Replace** `_fit_seed_varmean(chain, dens, eff, is_seed, seed_w)` (returns `MonotoneVarMean`) with
   `_fit_loglinear_varmean(...)` → returns `(a, b, x̄)` (or a tiny frozen dataclass), computed by §3. Keep the
   exact same edge-gathering loop (the sufficient statistics are unchanged).
2. `_gdna_seed_estimate`: return `(rho_global, (a, b), var_mean)` — replacing the spline object with the 2 (or
   3) fit coefficients. Update its docstring.
3. `_global_logprior`: signature `sigma2_g` becomes the coeffs; compute
   `s2_between = max(a + b·log(max(ρ_global, ε)), 0.0)` (evaluated at the passed `rho_global`). Update docstring.
4. `node_sweep`: the `_capture` diagnostic stores the coeffs instead of the spline (or drop the field).
5. **Delete** `src/rigel/calibration/variance_model.py` (432 lines: `MonotoneVarMean`, `MonotoneMean`,
   `_bspline_design`, `_select_lambda`, `_slope_evidence_weight`, `_fit_monotone`) and
   `tests/calibration/test_variance_model.py` — confirmed dead once the spline is gone (`grep` shows only the
   bp_solver usage). Remove the `from .variance_model import MonotoneVarMean` import.

## 6. Validation (in order)
1. Controlled-ε amplification (`scripts/debug/isolate_amplification.py 1e-10`) — pass-1 & pass-2 Δ ≤ machine
   level; amplification ≤1×.
2. Cross-process determinism (multithreaded scan, ≥3 processes) — bit-identical (unstranded/stranded/no-gDNA).
3. σ̂²_g magnitude check — ≈31 at gdna300 ρ_global, ≈0 at no-gDNA (matches the spline ⇒ no-regression-by-value).
4. 16-scenario net-flow + abundance benchmark vs the committed baseline — **no substantial regression**.
5. `pytest tests/ --update-golden` (goldens shift: σ²_g scalar ≠ spline value) → `pytest tests/` green.

## 7. M2 IS REQUIRED — the "eliminate the spline" restructure was tested and FAILS (2026-07-05)

σ²_g lives ONLY in `_global_logprior`, which runs in BOTH passes. Two toggled experiments settled whether the
pass-1-prior-free / restructure framing could eliminate it:

- **Pass-1 prior-free (RIGEL_PASS1_NO_GLOBAL):** VIABLE — no catastrophe; stranded unchanged; unstranded-capon
  flips under→over (±12%, symmetric, and the KDE is the real lever there); **no-gDNA stays tiny (32)** because
  pass-2's floor+KDE does the FP suppression, not pass-1. Good simplification (clean training data) — but it
  does NOT fix determinism, because the σ²_g spline is used by **pass-2**'s `_global_logprior`.
- **Drop pass-2's genome-wide `ρ_global`+σ²_g baseline (keep floor + KDE):** FAILS — no-gDNA capture-on exon
  false-positive gDNA **explodes 32 → 110,684**. The depleted floor covers only intergenic+intron REGIONS; the
  **genome-wide baseline is load-bearing for no-gDNA EXON suppression.** So it cannot be removed.

**Conclusion: the density-dependent σ²_g baseline is required in pass-2 → M2 is the fix (the spline cannot be
eliminated).** This *validates* M2's shape: σ²_g(ρ_global)≈0 at low ρ_global suppresses no-gDNA exons; ≈31 at
high ρ_global spares enriched exons — exactly the monotone law M2 fits deterministically.

**Independently**, pass-1-prior-free is a clean, viable simplification (training pass ⇒ raw data for the KDE)
worth adopting on its own merits — but it is orthogonal to the determinism fix. Sequence: ship M2 (determinism)
first; consider pass-1-prior-free with the deferred floor/KDE redesign.
