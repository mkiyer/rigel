# Session handoff — 2026-07-24

Read `ROADMAP.md` first, then this. This records what the session did, the exact code state, what remains,
and the prompt to start the next session.

---

## 1. What this session accomplished

1. **Fixed the unified solver’s composition arithmetic.** Six bugs; the message’s implied composition
   `Σ_c f_c` went from **75.6 → 1.02** (it must be ~1). Root cause: the relay carried three *independent*
   densities with no tie to the node’s own mass — the same defect the legacy `(λ,θ)` relay had retired. Fix =
   a composition **pin** at each node + relay-frame + per-face mature routing + coherent own-precision.
2. **Removed an over-applied relay pin** (found by an adversarial audit): pinning *inside* the relay cancelled
   the reframe. After removal, unified message modes are essentially correct — **exon message f_g 0.682 vs
   oracle 0.677**, boundary 0.666 vs 0.646. The unified solver’s **mode is now good; its variance is not.**
3. **★ Resurrected the intron-factory precision → shipped as default.** The factory’s NegBinom curvature was
   never registered as composition evidence, so a factory-solved intron had `τ = 0` and could not propagate.
   Fix = `bp_solver._lambda_factor_precision` (curvature → `τ_λ`), and **`intron_factory` flipped to default
   `True`**. Pass-0 vs oracle, 32-scenario suite: **mwae 0.1361 → 0.0949, corr 0.688 → 0.736**, 20 better /
   1 worse / 11 flat; intron mwae 0.1781 → 0.0117; every stranded scenario better-or-flat.
4. **Derived + MC-validated the graft/peel variance laws**, then found they were in the wrong
   parameterization (ratio vs per-component). See `variance_model_handoff.md` — this is what next session
   redoes.
5. **Measured junction overdispersion ω:** synthetic suite is **Poisson by construction (ω = 0)**, real
   ω ≈ 0.02 and **not fittable**. Decision: **model counts as Poisson**, no ω term.
6. **Established a key architectural fact:** the strand likelihood constrains only the **tilt**, never `f_g`,
   so AMBIG nodes get **zero f_g information from strand at any κ** (memory
   `strand_likelihood_constrains_tilt_not_fg`).
7. **Cleanup:** archived 51 stale calibration docs → `archive/`; wrote `ROADMAP.md`, `variance_model_handoff.md`,
   this handoff; updated CLAUDE.md + memory.

## 2. Exact code state (UNCOMMITTED — the working tree holds everything)

Nothing is committed. `git status` (branch `calib-ambig-init-wip`):

```
 M src/rigel/calibration/bp_solver.py        # unified-solver fixes + _lambda_factor_precision + capture hooks
 M src/rigel/calibration/node_geometry.py     # per-face spliced/geometry accessors for the unified path
 M src/rigel/config.py                        # intron_factory: bool = True  (the default flip)
 M tests/calibration/test_bp_solver.py        # +2 tests: _lambda_factor_precision
?? src/rigel/calibration/enrichment_frame.py   # NEW: pure reframe/mode/variance primitives (the pure layer)
?? tests/calibration/test_enrichment_frame.py  # NEW
?? tests/calibration/test_message_frames.py    # NEW
```

**Gates now:** `pytest tests/calibration tests/native` = **373 pass**, ruff clean on `src/ tests/`. **Goldens:
7 `test_golden_output.py` failures** — 5 are pre-existing (reproduce at HEAD), 2 are from the intron-factory
default flip and are EXPECTED (regenerate with `pytest tests/ --update-golden` when ready to lock the new
default). Flag-off phantom guard is no longer 3,704,635 because the factory is now on by default.

**Recommend committing before the next session** so this work is durable (I did not commit — that is your
call). A reasonable message: `calib(pass-0): intron-factory precision → default on; unified-solver composition
fixes; variance-model handoff`.

## 3. What remains (ordered — this is §6 of the ROADMAP, tactical detail)

1. **Derive the correct variance model** (`variance_model_handoff.md` §3): the **per-component density
   variance** `Var(log ρ_c) = Σ_k (ρ_k/ρ_c)²·v_k + σ²_transfer`, which unifies the intergenic-anchor count
   term, the graft share-weighting, and a transfer hook. **MC-validate including the pure-gDNA anchor limit.**
2. **Derive the composition transfer variance** (`σ²_transfer` redo, §4): direction-dependent (`Var(log r)`
   on the peel, ~0 on the graft), covering **relay and combine** together.
3. **Implement both as pure, tested functions** in the pure layer (extend `enrichment_frame.py`, rename it
   `message_arithmetic.py` when it clearly owns all message math). The unified closure must SHRINK.
4. **Loop** (worst scenario → dissect → fix) until the unified solver ≥ legacy 0.0949, no stranded regression.
5. **Converge:** flip `RIGEL_UNIFIED` on, delete `_scan` + flags, collapse `_UNIFIED_*` gates, regen goldens.

## 4. Open problems the next session should know about

* **The phase-2 refit REGRESSES unstranded capture-ON badly** (pass-0 0.2232 → refit **0.6079** on
  gdna300_ss_0.50_nrna_none_capture_on). It overwrites the factory’s gain, so production (`calib_refit_iters=1`)
  nets only ~−0.0035. Once pass-0 is clean this becomes the **dominant production error** and is the gDNA
  hyperprior workstream’s problem.
* **Remaining pass-0 error is exon-dominated** (83.7% of suite error after the factory). The lead is the
  gDNA-message enrichment cliff at >200×-enriched exons — which is exactly what the *correct variance model*
  (weighting the message by how well the source knows its composition) should fix. `RIGEL_E2=1` (an
  unconditional reframe) fixes high-gDNA exons but phantoms low-gDNA ones — it is NOT the answer; the
  evidence-weighted variance model is.
* **`strand_deconv.py:50` docstring is wrong** — `gdna_frac_var` is documented as `Var(f_g)` but carries a
  **log-scale** variance (max 4.049 > ¼). Use `NodeDeconv.lam_var` for `Var(λ)`. Audit `composition_logvar`’s
  callers for this confusion.

## 5. ▶ PROMPT TO START THE NEXT SESSION (copy-paste)

> We are developing the **calibration** stage of Rigel — specifically the prior-free first pass (“pass-0”).
> Read `docs/calibration/ROADMAP.md`, then `docs/calibration/variance_model_handoff.md`, then
> `docs/calibration/SESSION_2026_07_24_HANDOFF.md`. Do not read anything in `docs/calibration/archive/` unless
> I point you there — those are stale.
>
> The working tree (branch `calib-ambig-init-wip`) holds uncommitted work: the intron factory is now on by
> default (a real win, keep it), and the new **unified solver** (`RIGEL_UNIFIED=1`, default off) has a correct
> *mode* but a **wrong variance model**. The old variance model assumed genome-wide density uniformity, which
> hybrid capture breaks; the unified solver compares *compositions* via enrichment ratios, and its message
> propagation needs a new variance model.
>
> **First task: derive the new variance model from scratch**, per `variance_model_handoff.md` §3–§4 — the
> per-component density variance for the graft, the difference-variance for the peel, and a
> direction-dependent transfer variance — model counts as **Poisson** (no overdispersion). Validate every law
> by Monte-Carlo (`scripts/debug/message_precision_mc.py`) **before** writing solver code, including the
> pure-gDNA anchor limit where the old ratio form was singular. Then plan a clean implementation that puts the
> math in the pure arithmetic layer and shrinks the `_unified_solve` closure. Keep the debug loop:
> full ambig_dense_10mb suite → worst scenario by error mass → dissect worst nodes → root cause → fix. No
> magic numbers; pause and discuss before any new constant.
