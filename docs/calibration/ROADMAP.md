# Calibration ROADMAP — status, architecture, and the path to production

**This is the single entry point for calibration work. Read it first.** Last updated: 2026-07-24.

> **Status in one line:** we are mid-migration from a density-transfer solver to a composition
> (enrichment-ratio) solver for the **prior-free first pass (“pass-0”)**. The new solver is not yet correct,
> and — the current blocker — **its variance model is wrong** (the old one assumed genome-wide density
> uniformity, which hybrid capture breaks). **NOT ready to ship.**

The only other docs that are live (everything else is in `archive/`, kept for history, NOT to be referenced):
* `CALIBRATION_ARCHITECTURE.md` — the authoritative theory (count-zero-information; the three information
  sources). Still correct; read second.
* `unified_solver_design.md` — the target solver's architecture (the reframe + ÷M_dst mode). Its **precision /
  variance sections (§8 R1–R4) are SUPERSEDED** by `variance_model_handoff.md`; the mode design stands.
* `gdna_intron_factory_design.md` — a shipped feature (the intron gDNA factory). Live.
* `variance_model_handoff.md` — the variance-model derivation work, to be **redone** next session (handoff).
* `SESSION_2026_07_24_HANDOFF.md` — what this session did, what remains, the next-session prompt.

---

## 1. What calibration is, and what “pass-0” means

Calibration deconvolves each genomic node’s **unspliced** fragment mass into a composition
**(f_rna₊, f_rna₋, f_g)** — sense-RNA / antisense-RNA / gDNA. It runs in stages:

```
   PASS-0 (prior-free)  →  fit the gDNA HYPERPRIOR  →  RE-SOLVE (with the hyperprior)
   an APPROXIMATION        (from the pass-0 result)     the actual answer
```

* **Pass-0 is not a solution — it is an approximation.** It solves each node from only the two *intrinsic*
  information sources (the strand likelihood + cross-node imputation), with **no population prior**. On
  unstranded data the strand likelihood is flat (`CALIBRATION_ARCHITECTURE.md`), so pass-0 leans almost
  entirely on cross-node message propagation.
* **The pass-0 result trains a gDNA hyperprior** — the population baseline gDNA density landscape.
* **The hyperprior is then required to re-solve**, in particular to resolve **AMBIG** (opposite-strand
  overlap) nodes, which have no intrinsic gDNA/RNA signal at all (see
  `strand_likelihood_constrains_tilt_not_fg` in memory: on AMBIG nodes the strand likelihood constrains only
  the tilt, never f_g).

Everything below is **pass-0** unless stated. The hyperprior fit is a separate, also-WIP workstream (§4).

## 2. The two-solver problem (the core mess)

`bp_solver.py` (~1870 lines) currently contains **two** solve paths:

| path | flag | status | what it is |
|---|---|---|---|
| **`_scan`** (legacy) | default (`RIGEL_UNIFIED=0`) | **production path today** | density-transfer relay; ~470 lines + many experiment flags (`_B1B`, `_N4A`, `_N4B`, `_E2`, …) |
| **`_unified_solve`** | `RIGEL_UNIFIED=1` | **the target, default OFF** | composition (enrichment-ratio) relay; ~300-line nested closure |

**Goal: productionize `_unified_solve`, put it on the default path, and delete `_scan` and all its flags.**
This has not happened because (a) the unified solver still has bugs / is not fully correct, and (b) there is
no correct variance model for its message propagation (§3). Until both are fixed, the unified solver *loses*
the A/B (measured pass-0-vs-oracle mwae: unified **0.1280** vs legacy **0.0949**), so it cannot be flipped on.

## 3. ⛔ THE BLOCKER — the variance model

The conceptual shift is: the old solver compared **absolute densities** between nodes; the new solver compares
**compositions**, normalizing between nodes by an **enrichment ratio** that cancels hybrid-capture
enrichment/depletion. This works for the **mode** in almost all cases (measured: exon message f_g 0.682 vs
oracle 0.677). But the **variance model** was built for density transfer under a **genome-wide uniformity
assumption**, and hybrid capture breaks uniformity (it enriches targeted regions, depletes off-target ones).
So the variance model is now wrong for the solver we are building.

The derivation work done this session, and the reason it must be redone, is in **`variance_model_handoff.md`**.
The short version: we derived a share-weighted precision for the RNA graft and a difference-variance for the
peel and validated them by Monte-Carlo — but discovered late that they were written in the **ratio (k)
parameterization**, which is singular at the pure-gDNA anchor, whereas the code uses the **÷M_dst density
mode**, for which the right form is a **per-component density variance**. The unifying framework
(`Var(log ρ_c) = Σ_k (ρ_k/ρ_c)²·v_k + σ²_transfer`) is identified but **not yet derived clean or validated**,
and the **transfer-variance term (`σ²_transfer`) needs a ground-up redo** for the composition solver — the
current `σ²_transfer` is a density-uniformity proxy.

## 4. The gDNA hyperprior (the second WIP workstream)

Fitting the gDNA hyperprior from the pass-0 result was working well and is much of the way there, but
**pass-0 development this session exposed solver bugs**, which is why effort pivoted back to the solver. The
hyperprior is what makes the re-solve possible (esp. AMBIG). It is not done. (Production config ships
`calib_refit_iters=1` — one hyperprior-refit pass — but that refit currently *regresses* unstranded
capture-ON pass-0 badly; see the session handoff §“open problems”.)

## 5. What SHIPPED (is correct and on by default) vs WIP

**Shipped / correct (default path):**
* The **gDNA intron factory** — `intron_factory=True` (default, changed this session). Deconvolves introns
  against the intergenic background NegBinom and now **carries its derived precision** as composition evidence
  (`_lambda_factor_precision`), so a factory-solved intron can actually propagate. Pass-0 vs oracle over the
  32-scenario suite: **mwae 0.1361 → 0.0949, corr 0.688 → 0.736**, 20 better / 1 worse / 11 flat; intron mwae
  0.1781 → 0.0117; every stranded scenario better-or-flat.

**Work in progress (not shipped, flag-gated off):**
* The **unified solver** (`RIGEL_UNIFIED`) — mode largely correct, variance model wrong (§3), still behind.
* The **variance model** — being redone (§3, `variance_model_handoff.md`).
* The **gDNA hyperprior refit** (§4).

## 6. The path to production (ordered)

1. **Derive the correct variance model** for the composition solver — per-component density variance +
   a composition-appropriate transfer variance. Validate by MC. *(This is the next session’s first task.)*
2. **Implement it in `_unified_solve`**, as pure, tested arithmetic functions (mirroring `enrichment_frame.py`
   / `gdna_intron_factory.py`) so the closure shrinks. Re-run the loop (worst scenario → dissect → fix).
3. **Make the unified solver win the A/B** (≥ legacy 0.0949, no stranded regression).
4. **Converge:** flip `RIGEL_UNIFIED` on, delete `_scan` + its flags, collapse the `_UNIFIED_*` diagnostic
   gates, regenerate goldens. `bp_solver.py` roughly halves.
5. **Return to the gDNA hyperprior refit** (§4) on the clean solver, then the re-solve, then a ship candidate.

## 7. How we work (methodology — see memory `pass0_debug_iteration_loop`)

Debug by the loop: **run the full ambig_dense_10mb suite → find the worst scenario (by error MASS) → dissect
its worst nodes → trace to root cause → fix → repeat.** Compare **pass-0 vs oracle only** (`calib_refit_iters=0`)
unless testing the refit. Everything is cached (`_selfsolve_cache`), so a scenario solve is ~1 s. Standing
directives: **no magic numbers** (pause and discuss before any new constant); develop on toy/controlled
scenarios; keep the module/constant count small. **⚠ The synthetic suite is Poisson by construction** (memory
`synthetic_suite_is_poisson_omega_zero`), so it **cannot validate anything overdispersion-dependent.**

Key tools (in `scripts/debug/`): `pass0_error_concentration.py` (suite → where the error is),
`pass0_node_dissect.py` (exact-replay ψ channel ablation + per-node dumps), `unified_message_audit.py`
(Σ_c f_c invariant), `message_precision_mc.py` (variance-law MC validation).
