# The background gDNA density prior — count model & the σ²_g estimator (M1)

**Status:** design + implementation (2026-07-05, branch `calib-gdna-accuracy`). Replaces the bistable
`MonotoneVarMean` σ²_g P-spline (the calibration-nondeterminism root cause, see
`calibrate_cross_process_nondeterminism.md`) with a theory-grounded, deterministic, single-scalar estimator.

## 1. The estimand

The **background gDNA density** `ρ_bg` — a single scalar rate (gDNA fragments per unit effective length) —
estimated from **intergenic + intron REGIONS** (contained-fragment counts). These are many independent
observations of the same background level: under hybrid capture, the *depleted off-target* level; without
capture, the flat *contamination* level. The prior built from it pulls a node's gDNA density toward `ρ_bg`
with a precision set by the background's spread.

## 2. The count model (grounded in RNA-seq practice — DESeq2 / edgeR use NB)

Node `i`: gDNA count `G_i` over gDNA eff-length `E_i`, mean `μ_i = ρ_bg·E_i`.

| model | Var(G_i) | verdict |
|---|---|---|
| **M0** Poisson | `μ_i` | too rigid — coverage overdisperses (GC, mappability, CNV) |
| **M1** Negative Binomial | `μ_i + α·μ_i²` | the standard (DESeq2/edgeR core); one overdispersion `α` |
| **M2** NB mean–dispersion trend | `α(μ)=a₀+a₁/μ` | 2-param deterministic curve (DESeq2 `parametric`) — fallback if a real trend survives on background-only nodes |
| **M3** Zero-inflated NB | + zero-inflation π | for the sparse/no-gDNA regime (excess zeros) |

## 3. The variance the prior needs — the key decomposition

The prior is a Gaussian on `log ρ_g`. By the delta method under NB, the spread of a background node's
log-density separates into two INDEPENDENT pieces:

```
Var(log ρ̂_i) ≈ Var(G_i)/μ_i²  =  1/(ρ_bg·E_i)  +  α
                                  └─ Poisson ──┘   └ overdispersion ┘
```

1. **`1/(ρ_bg·E_i)` — per-node Poisson/sampling noise.** Density/count-dependent, large for low-count nodes.
   Already handled per node by the count precision (the Poisson–Gamma rate posterior `var_mean = 1/(1+G)`,
   which also handles `G=0` gracefully — the zero-inflation guard).
2. **`α` — the overdispersion.** A **single constant scalar**: the irreducible between-background-node spread.
   **This is `σ²_g`.**

**Therefore σ²_g is ONE scalar (`α`), not a curve.** Two errors in the retired spline follow directly:
- it fit *total* variance vs μ nonparametrically, **conflating** the falling Poisson term with the constant α;
- its *increasing* σ²_g(μ) `[24.5→42.9]` was an **artifact of including the enriched single-strand-exon seeds**
  (high μ, high real variance). On the **background regions only**, α is a clean constant.
- and its GCV-λ `argmin` was bistable under machine-ε → the ~2.6% cross-run nondeterminism.

## 4. The prior's role (why weak, and how it serves each regime)

Precision `N = 1/(var_mean + σ²_g)`, **capped at one pseudo-observation** so it can never overrule a node's own
strand evidence. σ²_g = α sets the weakness, and that is exactly right:

- **Pass-1 (capture unknown → must be weak):** an enriched node sits far above `ρ_bg`; a weak (capped-N) prior
  centered at `ρ_bg` cannot drag it down, so enriched nodes are safe. A *sparse / zero-inflated* node with no
  strand evidence is gently shrunk toward `ρ_bg`. α = the true background overdispersion is this weak-but-real
  level — larger α ⇒ weaker (safer for enriched, less robustifying); smaller α ⇒ stronger.
- **Pass-2 (discriminate):** **bimodal** under capture (a node shrinks to the nearest of the depleted / enriched
  modes — the KDE), **unimodal** without capture (robust shrinkage to the single background). Pass-1 stays weak
  and agnostic; pass-2 does the capture-aware discrimination.

## 5. M1 — the estimator (implemented)

`σ²_g = α`, a single **NB overdispersion scalar** from the **background (intergenic+intron) regions only**, by
**method-of-moments** in the robust, deterministic, continuous form:

```
α̂ = weighted_mean_i[ max( (log ρ̂_i − log ρ_bg)² − 1/(G_i + ½), 0 ) ]     over background regions
```

- `(log ρ̂_i − log ρ_bg)²` = the observed squared log-density deviation; `1/(G_i+½)` = the per-node Poisson
  log-variance (subtracted so only the *excess* = overdispersion remains); `max(·,0)` = a variance is ≥0.
- **Deterministic & continuous** in the data (a weighted mean, no `argmin`/GCV/IRLS/spline) ⇒ ε-robust.
- **Background-only** (no enriched-exon seeds) ⇒ the clean constant α, not the spline's inflated/artifactual
  trend. Homogeneous background ⇒ a plain weighted mean is robust enough; escalate to a Huberized mean only if
  the benchmark shows outlier sensitivity.
- This is essentially what `_floor_estimate`'s `s2_floor` already computes for the depleted floor — M1 unifies:
  **one background level `ρ_bg` + one background overdispersion `α`**, used for both the genome-wide baseline and
  the floor.

**M2 fallback:** if a real mean–dispersion trend survives on background-only nodes, fit `α(μ)=a₀+a₁/μ` by
deterministic weighted least squares (still not a spline), evaluate at `ρ_bg`.

## 6. Validation plan
Controlled-ε amplification test (must be ≤1×) + cross-process determinism (bit-identical) + the 16-scenario
net-flow & abundance benchmark (must show **no substantial regression** vs the committed baseline; ideally an
improvement). Then delete `variance_model.py` + `test_variance_model.py` (dead), regen goldens, commit.

---

## 7. RESULT — M1 as a single scalar CANNOT work; σ²_g must be density-local (evaluated at ρ_global)

Implemented + measured M1 (single scalar, floored between-seed log-density variance). Determinism: perfect
(controlled ε → 0× amplification; unstranded + no-gDNA bit-identical cross-process). **But it regresses
accuracy**: capon-unstranded total gDNA 1.09M vs baseline ~1.47M (true 1.66M) — it makes the known leak worse.

Cause, measured directly: the load-bearing σ²_g is the CONDITIONAL variance AT ρ_global (a point on a
monotone-increasing σ²_g(μ) curve):

| condition | ρ_global | spline σ²_g@ρ_global | M1 scalar (unconditional) | s2_floor |
|---|---|---|---|---|
| gdna300 ss0.50 capon | 0.235 | **31** | 9.7 (3× too small → crush) | 0.0 |
| gdna300 ss0.99 capon | — | ~31 | 10.1 | — |
| no-gDNA ss0.50 capon | 5e-7 | ~0 | 1.7 | 0.0 |

A single (unconditional) scalar is weighted toward the many low-density background seeds and lands BELOW the
value at ρ_global. Five estimators tried (all fail): spline (accurate, bistable); per-edge mean (too large →
FP); interp-median (too small → crush); background-overdispersion s2_floor (≈0 → crush); floored between-seed
variance (≈10 → crush). **Density-dependence is essential and load-bearing — the μ-dependence I twice dismissed
is real.**

## 8. The two clean deterministic paths

**(A) M2 — deterministic parametric monotone σ²_g(μ), evaluated at ρ_global.** Fit σ²_g = a + b·log μ (or a
power law) by CLOSED-FORM weighted least squares on the seed-edge Poisson-corrected excess vs log μ — no GCV,
no IRLS, no active-set, so no bistable argmin. Replicates the spline's value (~31 at ρ_global) deterministically.
Minimal, no-regression-by-construction (it reproduces the load-bearing value). Keeps the current architecture.

**(B) Architectural restructure (the count-model framing).** Split the roles so the genome-wide baseline no
longer needs a density-dependent σ²_g: a UNIMODAL background prior (center ρ_floor + the background-overdispersion
scalar) for depleted/no-capture, and the Phase-2 KDE (BIMODAL) for the enriched/capture case; pass-1 stays weak.
Then σ²_g becomes the single background overdispersion scalar (M1) and the enrichment spread is handled by the
KDE modes, not by a wide σ²_g. Cleaner and aligned with the theory, but a larger change entangled with the
deferred pass-1-floor / pass-2-KDE work.

**Recommendation:** (A) M2 now — it is the surgical determinism fix that provably does not regress (it
reproduces the spline value deterministically); pursue (B) as part of the deferred floor/KDE redesign.
