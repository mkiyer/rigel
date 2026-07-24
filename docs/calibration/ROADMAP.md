# Calibration ROADMAP — status, architecture, and the path to production

**This is the single entry point for calibration work. Read it first.** Last updated: 2026-07-24.

> **Status in one line:** the migration to the composition (enrichment-ratio) solver for the **prior-free
> first pass (“pass-0”)** has **CONVERGED to one solver** — the legacy density-transfer `_scan` path is
> deleted and the unified solver is the sole path. The **initialization phase is now a hardened, unit-tested
> module (`node_init.py`)**. The current blocker is unchanged: **the message variance model is wrong** (the old
> one assumed genome-wide density uniformity, which hybrid capture breaks) — that is the immediate next task.
> **NOT ready to ship.**
>
> **Update 2026-07-24 (post-handoff session).** Retired `_scan` + all its flags/helpers (`bp_solver.py` 1871 →
> ~730 lines); extracted the per-node self-solve into `node_init.build_node_init` (the four init sources of
> `variance_model_concepts.md`, one unit test each), behavior-preserving (byte-identical to the pre-refactor
> unified path across all 32 scenarios). Goldens regenerated to the unified default. §2 below is now historical.

The only other docs that are live (everything else is in `archive/`, kept for history, NOT to be referenced):
* `CALIBRATION_ARCHITECTURE.md` — the authoritative theory (count-zero-information; the three information
  sources). Still correct; read second.
* `unified_solver_design.md` — the target solver's architecture (the reframe + ÷M_dst mode). Its **precision /
  variance sections (§8 R1–R4) are SUPERSEDED** by `variance_model_handoff.md`; the mode design stands.
* `gdna_intron_factory_design.md` — a shipped feature (the intron gDNA factory). Live.
* `variance_model_concepts.md` — the owner's spec for the **initialization** phase (the four sources) that
  `node_init.py` implements. Read for the init model.
* `variance_foundation_plan.md` — **the current implementation plan**: isolate the composition precision
  `(τ_λ,τ_θ)` (foundation) from the sampling `1/n` (messages). Hardened against `variance_foundation_critique.md`
  (a 6-lens adversarial critique). **The next task's spec.** `SESSION_2026_07_24_HANDOFF_2.md` is the handoff.
* `variance_model_handoff.md` — the MESSAGE variance-model derivation, to be redone **after** the foundation
  (`variance_foundation_plan.md`) lands. Live.
* `SESSION_2026_07_24_HANDOFF.md` — the prior session's record. Its "two-solver / uncommitted" specifics are
  now superseded (the convergence landed); the variance-model derivation §§ still stand.

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

## 2. The two-solver problem — RESOLVED (historical)

`bp_solver.py` used to contain **two** solve paths — the legacy density-transfer `_scan` (default) and the
composition `_unified_solve` (flag-gated). **As of 2026-07-24 this is resolved: `_scan` and all its flags
(`RIGEL_B1B/N4A/N4B/E2`, the `_UNIFIED` gate) and helpers were deleted; the unified solver is the sole path.**
`bp_solver.py` roughly a third of its former size (1871 → ~730 lines); the per-node INITIALIZATION self-solve
now lives in `node_init.py` (`build_node_init` — the four sources, unit-tested). The unified solver still
**loses the A/B** (measured pass-0-vs-oracle mwae: unified ~**0.15** vs the legacy-with-factory baseline) — this
is **expected and accepted on this WIP branch**; the variance model (§3) is what recovers it. Nothing ships
until it does.

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

0. ✅ **DONE (2026-07-24):** converge to one solver — deleted `_scan` + flags + the `_UNIFIED_*` gates;
   extracted + hardened + unit-tested the **initialization** phase (`node_init.py`, the four sources);
   regenerated goldens. `bp_solver.py` 1871 → ~730 lines. Behavior-preserving (byte-identical A/B).
1. ✅ **DONE — the variance FOUNDATION** (`variance_foundation_plan.md` v4, `variance_foundation_proposal.md`).
   The honest local composition precision is a **single Schur-marginal scalar `τ_λ`** (approach E) — a diagonal
   `(τ_λ, τ_θ)` is prohibited (the strand Fisher is rank-1). Derived by a 5-approach workflow, numerically
   validated + independently re-verified, and independently critiqued (both converge on approach E + Option B).
   Phase 1 (the strand-gate bug fix — AMBIG gets zero strand f_g credit) is landed + A/B-validated (stranded arm
   improves, no regression), committed `c6df8c50`.
2. **The MESSAGE variance model** (`variance_model_handoff.md`) — **the NEXT task**. Builds the message precision
   on the `τ_λ` foundation; absorbs the deferred **composition/sampling separation** (move `1/n` out of the
   fusion weight into the transport transfer function + the `struct_lock` hard-override; plan v4 §4). Per-node
   density variance + a composition-appropriate transfer variance. Validate by MC + the per-condition A/B.
3. **Make the unified solver win the A/B** (≥ the legacy-with-factory baseline, no stranded regression).
4. **Return to the gDNA hyperprior refit** (§4) on the clean solver, then the re-solve, then a ship candidate.

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
