> ## ⛔ SUPERSEDED — THIS IS NOT THE LIVE HANDOFF
> The live handoff is **`SESSION_2026_07_26_HANDOFF_14.md`**; the entry point is **`ROADMAP.md`**.
> This file is kept as HISTORY — what was tried, measured and refuted. Its numbers describe the
> code as it was on its own date and many are now superseded. **Do not act on it without checking
> the live handoff first.** Its own DO-NOT-RE-RUN findings (if it has any) remain HERE and still
> stand — HANDOFF_14 §3 indexes which files carry them.

# Session handoff — 2026-07-24 (session 2: convergence + the variance foundation)

Read `ROADMAP.md` first, then this. **Supersedes `SESSION_2026_07_24_HANDOFF.md`** (session 1).

> **UPDATE (later in session 2): the variance-foundation derivation is DONE and the fix is landing.** A
> 5-approach derivation workflow + independent verification (`scratchpad/verify_foundation.py`, all claims PASS)
> resolved the rank-1 blocker: the honest foundational gDNA-level precision is the **Schur-marginal scalar `τ_λ`**
> (**approach E**, `variance_foundation_proposal.md`) — the strand pins f_g only for single-strand nodes; for
> AMBIG it contributes ZERO (it constrains only the tilt). The plan of record is now `variance_foundation_plan.md`
> v3. **Phase 1 (the `τ_λ` strand-gate bug fix) is implemented + validated** (A/B: stranded arm improves —
> refit=1 aggregate 0.0910→0.0889, stranded/capON −0.0102, unstranded neutral, no regression). Phases 2 (the 1/n
> composition/sampling separation) and 3 (cleanup) remain. §3–§5 below (the "critique the plan / open blockers")
> are now historical — the blockers were resolved by the derivation.

---

## 1. What session 2 accomplished

1. **Converged to ONE solver.** Deleted the legacy `_scan` path + all its flags (`RIGEL_B1B/N4A/N4B/E2`, the
   `_UNIFIED` gate) and helpers (`_mode_shift/_mode_density/_pred_precision/_routing_precision/_fold_lambda/
   _compile_strand_evidence/_lambda_factor_precision/_boundary_spliced_mass_increment`). `bp_solver.py`
   1871 → 726 lines. **Committed `b224bcd8` (WIP).**
2. **Extracted + hardened INITIALIZATION into `node_init.py`** (`build_node_init`): each node's own
   `(density, precision)` from the four sources of `variance_model_concepts.md` — MEASURED (intergenic→Poisson),
   density-deconvolution, strand-deconv (I_strand deadband), unsolved-default. Pure precision arithmetic
   (`own_composition_logvar`, `own_precision`); one unit test per source. `node_global_geometry`/
   `node_total_density` moved to `node_geometry.py`. **Behavior-preserving: byte-identical to the pre-refactor
   unified path across all 32 `ambig_dense_10mb` scenarios** (A/B, max Δ=0). Goldens regenerated. In `b224bcd8`.
3. **Extracted the generic `density_deconv.py`** (counts + gDNA prior → gDNA/RNA deconvolution, NB precision);
   the intron factory is now `fit_intron_background` (the special case where the gDNA prior = the intergenic
   node distribution). `gdna_intron_factory.py` deleted; `factory_precision`→`density_factor_precision`.
   Byte-identical. **Committed `05b70516` (WIP).**
4. **Resolved the precision-model confusion** (owner comment #3): there are TWO precisions — the **composition**
   precision `(τ_λ, τ_θ)` (the node's belief, fraction space, foundational) and the **density/sampling**
   precision `1/n` (message propagation only). count→density is deterministic given composition; each atomic
   source lands its precision directly in composition space; `1/n` is a separate sampling quantity for transport.
   `node_init.own_precision = 1/(Var(log f)+1/n)` currently merges them.
5. **Wrote + adversarially hardened the variance-foundation implementation plan.** Draft → a 6-lens adversarial
   critique (27 findings: 4 blocker, 16 major) → v2 plan. Artifacts:
   * `docs/calibration/variance_foundation_plan.md` — the v2 plan (§0 open decisions, §9 finding→resolution ledger).
   * `docs/calibration/variance_foundation_critique.md` — the raw critique ledger (all 27 findings, file:line).

## 2. Exact code state

`git log --oneline -2`: `05b70516` (density_deconv), `b224bcd8` (convergence + node_init). Working tree holds
the **uncommitted** plan + critique + this handoff (docs only) + the ROADMAP/CLAUDE.md/memory updates.
Gates: `pytest tests/calibration tests/native` green (370 pass, 3 xfail, 1 xpass); `ruff check src/ tests/`
clean; full suite 1207 pass. Pass-0-vs-oracle (refit=0, unified default): mwae ~**0.150** (stranded arm
**0.026** clean; unstranded **0.225** = the variance-model gap). Loses the A/B vs legacy — **expected, WIP,
nothing ships** until the variance model lands.

## 3. The next task — the variance FOUNDATION (verify → critique → implement)

The goal (owner): **isolate setup + init + solver kernels from message propagation; the message-propagation
precision must be built on top of a solidified foundational variance model.** Separate the **composition**
precision `τ` (the foundation) from the **sampling** `1/n` (message layer). The plan is
`variance_foundation_plan.md`.

**⛔ Two BLOCKERS are OPEN design decisions — resolve BEFORE writing code (plan §0):**

* **D1 — the strand Fisher information is RANK-1.** It depends on `(λ,θ)` only through the scalar
  `p = ½+(κ−½)(1−f_g)sinθ`, so a **diagonal `(τ_λ,τ_θ)` cannot represent it** and adding a `τ_θ` alongside the
  existing `τ_λ=I_strand` **double-counts the strand**. This **refutes the τ_θ approach as drafted.**
  **RECOMMENDATION: DEFER τ_θ** — do the composition/sampling separation with the current RNA arm untouched
  (no double-count), and make τ_θ a separate future workstream that first resolves the rank-1 representation
  (single-precision-on-`t` + null ridge / full 2×2 / derive from the ψ joint posterior). The core separation
  does not need τ_θ.
* **D2 — Phase C is not byte-identical**, so it needs the falsifiable, per-condition, refit≥1 gate of plan §6.
  Ratify the noise-derived tolerance method.

**The sound core (implement, in order):** Phase B (`NodeState` + a **verbatim** byte-identical adapter — reuse
`own_composition_logvar+own_precision`, NOT `1/(1/τ+1/n)`) → Phase C (move `1/n` out of the **fusion weight
only**, with the **struct_lock HARD-OVERRIDE** to avoid the `∞→nan` at anchors — plan §5) → Phase D (honest
deletions: re-express the RNA liveness gate before dropping `rna_*_frac_var`; keep `strand_likelihood.py` as the
test oracle or port its tests; announce `lam_var` retirement to `message_precision_mc.py`) → Phase E (goldens +
the full gate). Every finding's resolution is in plan §9; the deep pitfalls are in §5 and the critique doc.

**Non-negotiables the critique surfaced:** the struct_lock nan (B1), the exact adapter (M1), the corrected
"dead" census — `rna_*_frac_var` and `strand_likelihood.py` are LIVE (M2/M3), the refit≥1 validation (B4), the
per-condition unstranded-capON guardrail (M8), and the debug-diagnostic contract (M6, plan §8) — the tools that
run the A/B must not be broken by the rename.

## 4. How the next session should proceed

1. **Verify the plan** against the code (spot-check the blocker claims: the rank-1 Fisher, the `∞→nan` in
   `_fuse`, the `own_precision` adapter algebra, the `_rna` liveness gate). Critique further — the critique doc
   is a starting point, not exhaustive.
2. **Ratify §0 D1/D2** (recommend: defer τ_θ; adopt the §6 gate). Then implement Phase B→E on the sound core.
3. Keep the debug loop (worst scenario → dissect → fix) and the standing directives (no magic numbers; develop
   on toys; keep modules/constants few).

## 5. ▶ PROMPT TO START THE NEXT SESSION (copy-paste)

> We are developing Rigel calibration — the pass-0 variance FOUNDATION. Read `docs/calibration/ROADMAP.md`,
> then `docs/calibration/variance_foundation_plan.md` (the implementation plan) and
> `docs/calibration/variance_foundation_critique.md` (the adversarial critique it was hardened against), then
> `docs/calibration/SESSION_2026_07_24_HANDOFF_2.md`. Do not read `docs/calibration/archive/`.
>
> The plan separates the node's **composition** precision `(τ_λ, τ_θ)` (the foundational belief) from the
> **density/sampling** precision `1/n` (message propagation only). **First: verify the plan and critique it
> further** — especially the two open blockers in §0: (D1) the strand Fisher is rank-1, which refutes a diagonal
> `τ_θ` (recommendation: defer τ_θ, do the sound core); (D2) the Phase-C gate. Confirm the blocker claims against
> the code (the `∞→nan` at struct_lock anchors, the exact `own_precision` adapter, the live `rna_*_frac_var`
> gate, the rank-1 Fisher). **Then implement the sound core: Phase B (byte-identical adapter) → C (move `1/n`
> out of the fusion weight, struct_lock hard-override) → D (honest deletions) → E (goldens + the falsifiable,
> per-condition, refit≥1 gate).** No magic numbers; pause and discuss before any new constant or any deviation
> from the plan's blocker fixes.
