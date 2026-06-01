# PR 5 — M-step + outer loop (the working calibrator)

**Parent plan:** [`../00_implementation_plan.md`](../00_implementation_plan.md) §7 (PR 5).
**Spec:** [`../../caljointmodel/03_inference.md`](../../caljointmodel/03_inference.md) §1, §5, §9.
**Builds on:** PR 4a (E-step/exposure), 4b (AMBIG sweep), 4c (gDNA FL).
**Type:** Python-only. **Build required:** no.
**Status:** **DESIGN FINAL — awaiting go-ahead.** Decisions III.1–.5 resolved;
constants enumerated (§II.1, Q6).

PR 5 turns the PR 4a **single pass** into the full EM. Three things change:

1. The **count channel goes live** — the E-step's NB log-Bayes-factor uses the
   previous iteration's RNA mean `μ_d` (silent on iteration 1, `μ_d = 0`), so the
   density signal (intergenic→gDNA, the paralog rescue of doc 01 §10) finally
   engages.
2. The **five library hyperparameters are fitted** each iteration (M-step), not
   held at their inits.
3. The loop **iterates to convergence** (mass-change tolerance), with the
   monotonicity sentinel already in `CalibrationResult`.

This is the PR that makes the calibrator *work*.

---

# Part I — Theory (doc 03 §1/§5/§9)

## I.1 The outer loop

Per reference-wide iteration (initialise `π_g=0.5`, `ω=1`, hyperparameters at
their PR 3 / 4a / 4c values):

```
for iter in 1..max_outer:                                  # config.max_outer_iterations = 25
    E-step over the 3 views (count LLR now uses μ_d = M_d_unspl from iter-1)   # §I.3
    exposure   = Gamma posterior, D1 aggregation (no ½)                        # PR 4a
    sweep      = re-impute AMBIG exposure from neighbours                      # PR 4b
    M-step     = fit ρ_0, ε_s (closed) ; φ, ρ_d_bb (1-D Newton)               # §I.2; ρ_r_bb fixed
    π_g_prior  = ω ρ_0 L_eff / N   (data-driven, doc 03 §5.6)
    δ = max_r |M_g_tot − M_g_tot_prev| / (M_g_tot_prev + 1)
    if δ < mass_rel_tol: converged; break                                     # = 1e-4
```

## I.2 The M-step (doc 03 §5)

| Param | Update | Note |
|---|---|---|
| `ρ_0` | `Σ M_g_tot / Σ(ω·L_eff)` | closed form (§5.1) |
| `ε_s` | `(1 + Σ π_g·n_s) / (1 + Σ n_s)` | closed, Beta(1,1) failsafe (§5.5) |
| `φ` | 1-D Newton on the NB NLL, MoM warm start, bounds `(phi_floor, _PHI_MAX)` | §5.4 |
| `ρ_d_bb` | 1-D Newton on the gDNA BB NLL (`κ_d=0.5` fixed), MoM warm start | §5.2 |
| `ρ_r_bb` | **fixed** at PR 3's spliced fit (not re-fit) | III.1 |

So the M-step fits **four** scalars (`ρ_0, ε_s, φ, ρ_d_bb`); `κ_rna` and `ρ_r_bb`
stay fixed from the clean spliced channel (PR 3), and `κ_d=0.5` is biology.

The two strand Newtons consume **soft-allocated sense counts** (doc 03 §3.5):
`k⁺_g = max(k⁺ − κ_rna·M_d_unspl, 0)` (gDNA-attributable) and
`k⁺_d = max(k⁺ − 0.5·M_g_unspl, 0)` (RNA-attributable), so the E-step must also
surface the oriented `k_sense` and `M_d_unspl` (a small `Allocation` extension).
Each Newton is a single global scalar via `scipy.optimize.minimize_scalar`
(bounded), warm-started by the moment estimator — no per-region optimisation.

## I.3 Count channel live

The E-step already takes `m_d_unspl_prev` (PR 4a). PR 5 feeds it the previous
iteration's `M_d_unspl` per view (iteration 1: zeros → the channel is silent, as
in PR 4a). So no E-step change beyond threading the per-iteration `M_d_unspl`,
`ω`, and `π_g_prior`.

## I.4 Convergence (doc 03 §9)

Standard EM: each step weakly increases the observed-data log-likelihood. The
mass-change diagnostic `δ` must be non-increasing; `CalibrationResult.__post_init__`
already raises `CalibrationConvergenceError` if it increases (the runtime
sentinel). `n_iterations` / `converged` / the full `mass_change_history` are
reported. Typical convergence 5–10 iterations.

---

# Part II — Implementation plan

```
src/rigel/calibration/
  mstep.py     # NEW: ρ_0, ε_s (closed); φ, ρ_d_bb[, ρ_r_bb] (1-D Newton + MoM); π_g_prior
  estep.py     # EDIT: Allocation also surfaces k_sense + m_d_unspl (for the BB Newtons)
  calibrate.py # EDIT: replace the single pass with the outer loop (E → exposure → sweep → M)
```

- `calibrate()` keeps its signature (PR 4c's `gdna_fl_pmf`); the gDNA FL stays a
  **fixed input** to the loop (one-shot from PR 4c; iterative refinement is
  §III.2). Returns the converged result (`n_iterations ≥ 1`, real
  `mass_change_history`).
### II.1 Constants & heuristics — full enumeration (Q6)

**A. Reused (already in code/config; PR 5 adds no new value):**

| Name | Value | Source | Role in PR 5 |
|---|---|---|---|
| `max_outer_iterations` | 25 | `CalibrationConfig` (PR 2) | outer-loop cap |
| `mass_rel_tol` | 1e-4 | `CalibrationConfig` (PR 2) | convergence `δ` threshold |
| `phi_floor` | 1e-8 | `CalibrationConfig` (PR 2, empirically validated) | φ Newton **lower** bound + exposure floor |
| `_BB_FLOOR` | 1e-6 | `strand_balance` (PR 3) | ρ_d_bb Newton bounds `(_BB_FLOOR, 1−_BB_FLOOR)` |
| `_PI_CLIP` | 1e-6 | `estep` (PR 4a) | π_g + π_g_prior clip |
| `_RHO_D_BB_INIT` | 0.01 | `calibrate` (PR 4a) | ρ_d_bb warm start, iter 1 |
| `_EPS_S_INIT` | 1e-3 | `calibrate` (PR 4a) | ε_s value, iter 1 |
| `_PI_G_PRIOR_INIT` | 0.5 | `calibrate` (PR 4a) | π_g_prior, iter 1 |

**B. NEW in PR 5 — three only (all machine floors / unit priors / range caps;
none flip qualitative behaviour):**

| Name | Value | Role | Kind |
|---|---|---|---|
| `_PHI_MAX` | 100.0 | φ Newton **upper** bound (NB-overdispersion ceiling) | range guard (doc 03 §8 `1e2`) |
| ε_s Beta(1,1) pseudocount | 1.0 | the `1+` in `ε_s = (1 + Σ π_g n_s)/(1 + Σ n_s)` | unit-strength Bayesian prior (doc 03 §5.5) |
| `_PI_FLOOR` | 1e-9 | `max(N_r − μ_g, _PI_FLOOR)` in the π_g_prior update | machine floor / no-div-0 (doc 03 §5.6) |

*(If a `log()` guard is needed inside an NLL, `np.finfo(float64).tiny` — machine
epsilon, not a tunable.)*

**C. Derived (not free parameters):** the MoM warm starts — `φ_init =
max(Var(nᵤ)/n̄ᵤ² − 1/n̄ᵤ, phi_floor)`; `ρ_d_bb_init` from the §5.2 moment formula
using `0.25 = κ_d(1−κ_d)` (κ_d=0.5 fixed) — and the `+1` unit-smoothing in
`δ = max|ΔM_g_tot|/(M_g_tot_prev + 1)` (doc 03 §1; already in PR 4a).

**φ-floor reconciliation:** doc 03 §8 lists `_PHI_FLOOR = 1e-6`, but `CalibrationConfig.phi_floor = 1e-8`
is the empirically-validated value (PR 2). PR 5 uses **`config.phi_floor` (1e-8)**
as the φ Newton lower bound — no duplicate 1e-6 floor.
- **Tests:** M-step closed forms vs hand calc; each Newton recovers a known
  dispersion from synthetic counts; the loop converges + is monotone on a
  synthetic locus; the paralog case (doc 01 §10) — count channel pushes a
  126/26 gDNA region to mostly-gDNA over iterations; `converged=True`,
  `n_iterations` small. Optional: `scripts/debug/dump_calibration_state.py`
  (master plan §6) for per-iteration diagnostics.
- **Acceptance gate:** new tests pass; PR 1–4 green; `ruff` clean; full-suite
  failure mode stays the post-calibrate `NotImplementedError` (PR 6).

---

# Part III — Decisions (resolved)

- **III.1 — `ρ_r_bb` fixed from PR 3 ✓.** Only unspliced fragments are
  deconvolved; spliced are fixed pure RNA, so `κ_rna` (mean, from the scan) and
  `ρ_r_bb` (overdispersion, from the spliced channel) do not change with the EM.
  The M-step fits `ρ_0, ε_s, φ, ρ_d_bb` only. *(Future: deconvolved unspliced RNA
  could iteratively update `ρ_r_bb`; not now.)*
- **III.2 — One-shot gDNA FL ✓.** Iterative refinement would require maintaining
  a per-region/boundary FL map; deferred. The one-shot FL from the gDNA-dominated
  init regions/boundaries (INTERGENIC, INTRONIC + their EXON crossings) is the
  fixed loop input; benchmarks vs truth decide whether refinement is warranted.
- **III.3 — FL distributions are calibration OUTPUTS, not calibration
  likelihood ✓.** The calibration EM uses count + strand balance only; the gDNA
  and RNA FL distributions are **outputs** the downstream primary EM consumes.
  So the live count channel keeps `μ_d` as the observed RNA residual — no RNA
  `E_{f_rna}` term in calibration.
- **III.4 — scipy solver ✓ (provisional).** `scipy.optimize.minimize_scalar`
  (bounded) + MoM warm start to start; we keep evaluating whether a hand-rolled
  Newton or different bounds are warranted.
- **III.5 — Diagnostics ✓.** Add `scripts/debug/dump_calibration_state.py`
  (per-iteration `ρ_0, φ, ρ_d_bb, ε_s, δ`, + the fixed `κ_rna, ρ_r_bb`) to verify
  convergence/behaviour. Counted as a debug script, not calibrator surface.

## Rollback

Revert the sub-PR. `calibrate()` reverts to the PR 4a single pass (n_iterations=1,
converged=False). No on-disk artifacts change.
