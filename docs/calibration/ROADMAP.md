# Calibration ROADMAP — status, architecture, and the path to production

**This is the single entry point for calibration work. Read it first.** Last updated: 2026-07-25.

> **Status in one line:** the message-variance model is **COMPLETE** — derived, MC-validated, independently
> verified, implemented, and A/B-won. A message's precision is now
> `1/(Var(log f_c^src) + 1/n_src + σ²_transfer + b̂²)`: the source's earned composition+count precision, the
> reframe's **scale** uncertainty (M5 `Var(log r)`), and the **DerSimonian–Laird composition-mismatch** `b̂²`
> (M7) — the honest replacement for the `(log r)²` proxy. **Best aggregate on record: 0.0969 (refit=0) /
> 0.0828 (refit=1)** vs the 0.1267/0.1234 pre-fix baseline. **The next task is Phase 2 — the gDNA hyperprior**
> (`SESSION_2026_07_25_HANDOFF_6.md`). **NOT ready to ship** (the hyperprior refit still regresses
> unstranded-capON, and that is now the single binding constraint).
>
> **Update 2026-07-25 (DL cliff-term session).** `(log r)²` charged the WHOLE enrichment cliff as composition
> drift, which recovered the stranded arm but over-damped extreme capture. The delivered message error splits
> exactly into a composition-SHARE mismatch plus the reframe's scale noise, so the two are now priced
> separately; `b̂²` is estimated prior-free against the node's own self-solve (a two-study random-effects
> meta-analysis) with **no tuned constant**. Its safety property is exact: **a message can out-weigh a node's own
> belief only if it agrees within `√2·σ_own`** — the governing principle as an inequality rather than a knob.
> Attribution is clean and the two effects are disjoint: deleting `(log r)²` recovers verystrong, `b̂²` recovers
> stranded, neither costs the other. Also retired the dead NPMLE-projection σ²_transfer plumbing: **no prior of
> any kind now enters message precision.**
>
> **GOVERNING PRINCIPLE (owner):** **pass-0 must be WEAK and correctable** — an over-confident message that pins
> a node wrong is worse than a weak one slightly off. Prefer the under-confident option when unsure.
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
* **`SESSION_2026_07_25_HANDOFF_6.md` — the LIVE handoff. START HERE for the next session** (Phase 2: the gDNA
  hyperprior — fix the refit's unstranded-capON regression, then feed the hyperprior into DL's `v_own`); it has
  the A/B state to beat, the measured Phase-2 experiment, the invariants, and a kickoff prompt.
  `SESSION_2026_07_24_HANDOFF_5.md` (the DL term's plan, now DONE) and `..._HANDOFF_4.md` (the full arc + the
  audit/design workflow findings) are the reference for how the variance model got here. (Handoffs 1–3 are
  historical.)
* `message_variance_derivation.md` — the derived + MC-validated + independently-verified message-variance laws
  (M1–M5), the M6 combine finding, and the empirical results (§6). Live.
* `variance_model_concepts.md` — the owner's spec for the **initialization** phase (the four sources).
* `variance_foundation_proposal.md` — the SETTLED foundation model (approach E, the single Schur scalar `τ_λ`),
  derived + numerically verified. `variance_foundation_plan.md` v4 — the invariants + the deferred-work spec;
  `variance_foundation_critique.md` — the adversarial critique ledger.
* `variance_model_handoff.md` — the MESSAGE variance-model derivation substrate (§3-4), the NEXT task's math,
  built on the `τ_λ` foundation.

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

## 3. ✅ RESOLVED — the variance model (historical; see `message_variance_derivation.md` for the final laws)

The conceptual shift was: the old solver compared **absolute densities** between nodes; the new solver compares
**compositions**, normalizing by an **enrichment ratio** `r` that cancels hybrid-capture enrichment/depletion.
The old variance model assumed genome-wide density **uniformity**, which capture breaks — so it had to be
rebuilt for a composition transport. It now is (`message_variance_derivation.md`, laws M1–M7, every one
MC-validated in `scripts/debug/message_variance_mc.py`, which runs 0 failures end-to-end):

```
    p_message = 1 / ( Var(log f_c^src) + 1/n_src  +  σ²_transfer  +  b̂² )
                     \__ strand ___/   \_count_/    \_ SCALE _/    \_ COMPOSITION _/
```

The last two are the cross-cliff cost, and the session's central result is that they are **different objects**
that must be priced separately. The delivered message error splits EXACTLY (to machine precision) into
`log(s_c^src/s_c^dst,true) + log(r̂/r_true)` — a composition-SHARE mismatch plus the reframe's own scale noise.
`σ²_transfer = Var(log r)` (M5) prices the scale; `b̂²` (M7) prices the imputation PREMISE ("my neighbour and I
share a composition") being false. The retired `(log r)²` proxy charged the whole cliff as mismatch, which is
why it fixed the stranded arm but over-damped extreme capture, where the composition really is preserved
across a 1000× enrichment step.

`b̂²` is a population quantity we lack prior-free — but the destination has an **independent** estimate of its
own composition: its message-free self-solve. Two estimators of one quantity is a two-study random-effects
meta-analysis, so the **DerSimonian–Laird** between-study estimator supplies it with **no tuned constant**:
`b̂² = max(0, G² − v_msg − v_own)`, giving `p_eff = 1/max(v_msg, G² − v_own)`. Its three regimes fall out of
`v_own` (from the `τ_λ` foundation) with no gate, and the middle one is the safety property, exact: **a message
can out-weigh a node's own belief only if it agrees within `√2·σ_own`.** Where a node has NO composition
evidence (`τ_own = 0` — every AMBIG node, all unstranded data) `v_own = ∞` and the term switches itself off, so
messages propagate untouched exactly where they are the only information. That inertness is deliberate — and it
is also why the remaining AMBIG error is Phase 2's problem, not this term's.

## 4. ⛔ THE BLOCKER — the gDNA hyperprior refit

This is now the single binding constraint. Fitting the hyperprior from the pass-0 result was working well and
is much of the way there, and it is what makes the re-solve possible (esp. AMBIG). But production config ships
`calib_refit_iters=1` and **that refit regresses unstranded capture-ON**, which is the largest error-mass arm.

That regression is no longer just a scoring problem — it is what blocks the AMBIG fix. Measured this session
(`SESSION_2026_07_25_HANDOFF_6.md` §3): feeding the hyperprior's own λ-curvature into DL's `v_own` — the
committed Phase-2 step, ~6 lines — improves exactly the arms it should (stranded 0.0376→0.0333, verystrong
0.1292→0.1196, capture-off 0.0354→0.0168) and regresses unstranded-capON 0.1702→0.2177, because where the
hyperprior is wrong DL now *enforces* it against the messages that would have corrected it. **Fix the
hyperprior first; the AMBIG fix then lands almost for free.**

## 5. What SHIPPED (is correct and on by default) vs WIP

**Shipped / correct (default path):**
* The **gDNA intron factory** — `intron_factory=True` (default, changed this session). Deconvolves introns
  against the intergenic background NegBinom and now **carries its derived precision** as composition evidence
  (`_lambda_factor_precision`), so a factory-solved intron can actually propagate. Pass-0 vs oracle over the
  32-scenario suite: **mwae 0.1361 → 0.0949, corr 0.688 → 0.736**, 20 better / 1 worse / 11 flat; intron mwae
  0.1781 → 0.0117; every stranded scenario better-or-flat.

* The **message-variance model** (§3, M1–M7) — the sole message-precision law, no flags, no prior input.
  Pass-0-vs-oracle over the 32-scenario suite: **0.1267 → 0.0969 (refit=0), 0.1234 → 0.0828 (refit=1)**.

**Work in progress (NOT ready to ship):**
* The **gDNA hyperprior refit** (§4) — the blocker.
* **AMBIG nodes** — prior-free they have no composition evidence at all (`τ_own = 0`), so they are carried by
  messages alone and the DL term does not protect them. The minimal reproduction is the factor-1 bedrock toy
  (`test_gdna_sweep_factor1_ambig_recovery`, xfail): on a uniform ρ=0.5 chain the AMBIG node between two exact
  anchors reads **0.3914**. This is the designed weakness, NOT a mode defect — the shortfall shrinks
  monotonically with depth (21.7% at ρ=0.5 → 0.8% at ρ=5000), so the transported mode is right and what is
  missing is WEIGHT: messages arrive at their honest `1/n`, and ψ's uninformative reference holds the node off
  the vertex until the data earn it. Fixed by §4 (a trained prior), not by more damping.

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
2. ✅ **DONE — the MESSAGE variance model** (§3; `message_variance_derivation.md`, laws M1–M7). The M1–M5
   per-component transport laws, the single-λ combine on a three-stream relay (the M6 rank-1 double-count), and
   the cross-cliff cost split into the M5 scale term + the M7 DerSimonian–Laird composition-mismatch `b̂²`. The
   composition/sampling separation came out of the three-stream relay (τ carries composition only, the
   measurement stream the counts). Every law MC-validated; `b̂²` additionally re-derived and adversarially
   audited by a 4-agent workflow. Commits `44f1ecc6`…`1a3e0a89`.
3. ✅ **DONE — the unified solver wins the A/B**: 0.0969 (refit=0) / 0.0828 (refit=1) vs the 0.0949
   legacy-with-factory target — and note the current suite gained the hard `verystrong`/`gdna1`/`gdna5`
   scenarios since that number was set, so gate on in-run A/B deltas, not the absolute.
4. **The gDNA hyperprior refit** (§4) — **the NEXT task**, on the now-clean solver. Fix the unstranded-capON
   refit regression, then re-apply the measured 6-line Phase-2 step (hyperprior → DL `v_own`, which fixes AMBIG),
   then the re-solve, then a ship candidate. Exact plan + numbers: `SESSION_2026_07_25_HANDOFF_6.md`.

## 7. How we work (methodology — see memory `pass0_debug_iteration_loop`)

Debug by the loop: **run the full ambig_dense_10mb suite → find the worst scenario (by error MASS) → dissect
its worst nodes → trace to root cause → fix → repeat.** Compare **pass-0 vs oracle only** (`calib_refit_iters=0`)
unless testing the refit. Everything is cached (`_selfsolve_cache`), so a scenario solve is ~1 s. Standing
directives: **no magic numbers** (pause and discuss before any new constant); develop on toy/controlled
scenarios; keep the module/constant count small. **⚠ The synthetic suite is Poisson by construction** (memory
`synthetic_suite_is_poisson_omega_zero`), so it **cannot validate anything overdispersion-dependent.**

Key tools (in `scripts/debug/`): `pass0_oracle_bench.py` (THE A/B — `--arm`, `P0_REFIT`, `--report`),
`pass0_error_concentration.py` (suite → where the error is), `pass0_node_dissect.py` (exact-replay ψ channel
ablation + per-node dumps), `message_variance_mc.py` (the variance-law MC arbiter, M1–M7),
`unified_message_audit.py` (Σ_c f_c invariant). In `scratchpad/`: `dl_dissect.py` (error mass attributed by
DL-protection state / strand DOF / node class, + per-node message provenance), `dump_node.py`.
Env `RIGEL_S2T_OFF=1` disables both cliff terms for isolation; `_capture["_dl"]` publishes the per-message
gaps and the τ-stream kill.
