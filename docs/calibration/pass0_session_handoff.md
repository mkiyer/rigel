# Pass-0 Session Handoff (2026-07-21)

**Read this first on resume.** It captures the full state of the calibration **pass-0** work so a fresh session
continues seamlessly. Branch: `calib-ambig-init-wip` (22+ commits ahead of `main`; **do NOT push to main**). The
authoritative map is [`pass0_roadmap.md`](pass0_roadmap.md); this doc is the *session* handoff (state + the
immediate next task + the tools + the lessons).

---

## 0. THE IMMEDIATE NEXT TASK — re-derive the boundary rule (the enrichment cliff)

Follow the mantra: **(1) design & derive → (2) plan → (3) execute.** Careful, one step at a time.

**The finding that opened it** (`solve_gate_design.md` §3-4, `scripts/debug/dof_nodetype.py`): a node-type
CORRELATION test of the §6B DOF solvability criterion over the 20 unstranded `ambig_dense_10mb` scenarios shows
the criterion is **correct for regions** (single-strand + AMBIG: the nodes it withholds are genuine coin-flips,
corr ≈ 0; the ones it solves correlate 0.63–0.69) but **INVERTED for boundaries**: the boundaries it marks
"solvable" are coin-flips (corr **0.13**), the ones it "withholds" are meaningful (corr **0.68**).

**The owner's reframe (the key insight — this is NOT a technical solvability problem):** the exon–intron boundary
is the **hardest node** because it **straddles the enrichment cliff**. Under hybrid capture it is partially
enriched/depleted, and it receives **mixed messages** — from the **intron** (depleted) and the **exon** (enriched)
— which produce a **muddled gDNA estimate**. So the "boundary rule" is really about **correctly reconciling
cliff-discrepant messages at a partially-enriched crossing node**, not about a DOF gate. This is a standing TODO
area.

**Where this connects (do not re-derive these — build on them):**
- [`message_propagation_arithmetic.md`](message_propagation_arithmetic.md) — the boundary IS the cliff-crossing
  node in the `intron ⟷ boundary ⟷ exon` chain. §4a: the **cliff lives on `intron⟷boundary`** (handled by the
  FL-frame shift). §4b: `boundary⟷exon` is cliff-free-but-mature. So a boundary's belief must fuse a
  depleted-intron message and an enriched-exon message that are at *different* enrichment scales.
- [`transfer_variance_formal_derivation.md`](transfer_variance_formal_derivation.md) — `σ²_transfer` is the
  enrichment-crossing message variance (projection onto the enrichment NPMLE). **Check whether it is correctly
  damping the cliff-discrepant boundary messages** — a boundary receiving a depleted-intron and an enriched-exon
  message at high confidence would muddle if `σ²_transfer` under-damps the crossing.
- The boundary's **own** density is the *crossing* (partial) enrichment — neither the intron's depleted nor the
  exon's enriched. Its gDNA estimate should reflect that.

**Suggested first move (design/derive):** extend `scripts/debug/dof_nodetype.py` to *dissect boundaries* — split
the boundary group by (enrichment gap between its intron and exon flanks, message agreement/disagreement,
one-sided-spliced presence, `σ²_transfer`) and find WHAT distinguishes the meaningful boundaries (corr 0.68) from
the coin-flip ones (corr 0.13). The corrected rule should fall out of that. Likely it is a **cliff-aware** rule
(reconcile flank messages by their enrichment level), not a strand/mass DOF gate.

---

## 1. State of the pass-0 message system (the four concerns + the destination)

A message = per-component `(mode, precision, gate)`; the destination decides whether to solve. Status:

| concern | status | where |
|---|---|---|
| **MODE** (value) | derived + **MC-validated**; code proven exact (S1); helpers extracted+tested (S2); flip set characterized (S3); the flip (S4) **blocked** on the mature cliff | `message_propagation_arithmetic.md`, `message_mode_implementation_plan.md` |
| **EMISSION** (whether) | **LANDED** — always-emit, `ev_λ=∞`, `_pred_precision`; fixed the emission↔density coupling bug; retired the τ-gag bug class | `emission_and_precision_derivation.md`, commit `e53cbc8e` |
| **PRECISION** (confidence) | honest 3-term model `1/(Var+1/M+σ²_transfer)`, `Var` from reference-free `τ`; ceiling measured (80.5 % identifiability floor / 1.9 % precision-attackable) | `emission_and_precision_derivation.md`, `density_imputation_precision.md` |
| **GATE** (which strands) | structural `free_s`, stable | — |
| **SOLVE** (destination) | DOF criterion **validated for regions**, **inverted for boundaries** (the next task) | `solve_gate_design.md` |

**Hardening (F/G/H) is DONE** (commits `78bbf441`, `9b13f2a4`, `bcb254c5`): nan-robustness sweep, emission
invariant tests (density⊥evidence, `pr→0` as `τ→0`), and determinism (bit-reproducible within+cross-process at
`OMP_NUM_THREADS=1`). The solver is robust, principled, deterministic — a trustworthy base.

---

## 2. What landed this session (commits on `calib-ambig-init-wip`, newest first)

- `f9a12c81` — **corrected** the solve-gate record (retracted a flawed "refuted"; DOF valid for regions, boundary bug).
- `e062f3f5` — solve-gate first pass (the flawed "refuted" — superseded by `f9a12c81`; code reverted, byte-identical).
- `bcb254c5` — hardening H (determinism test).
- `9b13f2a4`, `78bbf441` — hardening F (nan sweep) + G (invariant tests).
- `d5b9dd44` — the pass-0 roadmap.
- `e53cbc8e` — **principled emission** + the emission↔density **coupling bug fix** (1 golden regen: antisense).
- `fecc8159` — the derivations + plans + MC harnesses (mode arithmetic, emission, mode-flip plan).
- earlier: the τ-gag fix (`243bd5ef`), A0/A1 refactors.

---

## 3. The tools (harnesses) — reuse these

- **`scripts/debug/dof_nodetype.py`** — the node-type CORRELATION analysis (THE tool for the boundary task; pools
  the cached unstranded ambig scenarios; prints per-type corr(f_g, oracle)). Extend it to dissect boundaries.
- **`scratchpad/pass0_bench.py`** — the 32-scenario `ambig_dense_10mb` benchmark; runs the real `calibrate` on a
  cached scan substrate (`_selfsolve_cache`), computes mass-weighted `f_g` error (mwae) + per-scenario detail.
  ⚠ **mwae counts UNSOLVED nodes at their arbitrary default** — do NOT use it to judge the solve-gate; use
  correlation / precision-weighted error on *solved* nodes (§5 lesson). Supports env-flag A/B (the `RIGEL_TAU_GAG`
  pattern). Copy + set `calib_refit_iters=1` for the hyperprior arm.
- **`scripts/debug/chain_mode_mc.py`** — MC validation of the message-mode arithmetic (binary + `capture_mode="proportional"`).
- **`scripts/debug/mode_verify.py`** — node-level `code_mode vs derived vs oracle` on toy scenarios (via `toy_prod`).
- **`scripts/debug/toy_prod.py`** — the production-faithful toy driver (real `calibrate` on hand-specified genes,
  per-region oracle f_g). Set `toy_prod.SCRATCH` to your scratch dir.
- **Cached suite:** `/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb` (32 conditions, `_selfsolve_cache` so
  the scan is cached and only `calibrate` re-runs). Skill: `calibration-benchmark`.
- **The `_capture` hook** in `bp_solver.node_sweep` (inert unless `_debug`/`_capture` passed) exposes per-node
  `f_g / var_g / prec_g,p,n / mode_g,p,n / free_pos,free_neg / solvable / mass_global / eff_global` + per-edge
  `_edge_modes`. Add fields temporarily as needed (guard with `if _capture is not None`; remove before commit).

---

## 4. The code (key locations in `src/rigel/calibration/`)

- **`bp_solver.py`** — the BP sweep. `node_sweep` (~285). The two message MODES: `_mode_shift` / `_mode_density`
  (pure, tested). `_pred_precision` (the honest-precision helper: 0 when unseen `v=∞` or countless; nan-safe).
  The evidence variance (three-state: lock=0 / evidence=jac²·1/τ / no-evidence=∞) ~596. **Structural emission**
  (no evidence gate) `emit_g/emit_p/emit_n` ~634. The `use_shift` mode predicate + the dormant `RIGEL_GATE_SHIFT`
  A/B flag (default off, mode-flip blocked) ~640. **`solvable = (fp|fn) & mass>0`** (structural; also drives
  `locked`/`struct_lock` — do NOT feed a DOF gate into it) line **416**. The `tau_lam` relay (line ~764,
  per-scan local; `tau_th` NOT relayed). The final **write-back** `f_g = where(solvable, solved, init)` line
  **~813** (where a solve-gate would go — as a SEPARATE `write_solvable`, not `solvable`).
- **`calibrate.py`** — the orchestrator. **Two passes:** PASS-1 prior-free (`gdna_prior=None`), then
  `for it in range(calib_refit_iters)`: **fit the hyperprior on the PASS-1 solve** (`_fit_gdna_hyperprior`, line
  108 — **solved single-strand + structural-gDNA regions only; NO AMBIG, NO boundary — non-circular**) → re-run
  `_sweep(gdna_hyperprior)` (line ~355, the refit RE-RUNS node_sweep). **A solve-gate must count the fitted prior
  as an identification source**, else the refit re-skips deferred nodes and the prior never resolves them.
- **`node_geometry.py`** — `init_beliefs` (G1 intergenic = f_g=1 locked var=0; **G2 single-strand = the strand
  solve**; **G3 AMBIG = f_g=1 at MAX var** — so a withheld AMBIG node reverts to all-gDNA, its default). The
  `NodeBelief` (f_pos/f_neg/f_g + vars; var=∞ ⇒ unsolved).

---

## 5. LESSONS from this session (do not repeat)

1. **Never measure error on zero-precision (unsolved) nodes.** Their `f_g` is an arbitrary init default
   (AMBIG → f_g=1), not a prediction. The mass-weighted mwae counts it → a false "regression." Use
   **correlation** (does the forced solve track the oracle?) or **precision-weighted error over solved nodes**.
   This is what falsely "refuted" the solve-gate.
2. **The hyperprior fit is already non-circular** (solved-only, excludes AMBIG/boundary). Deferred nodes do not
   pollute it. And a solve-gate's **refit** must count the prior as a source (or the deferred nodes never resolve).
3. **Do not jump to the hyperprior before finishing the solver tasks.** The hyperprior (refit=1) already halves
   the mwae on today's solve-everything pass-0 (0.2030 → 0.0998) — it is the eventual lever, but the solver tasks
   (the boundary rule, the correctness items) come first.
4. **A regression that only appears in a downstream metric may be an artifact of the metric or a bug in the new
   code — investigate node-by-node before concluding.** (The emission "golden shift" was a real bug found; the
   solve-gate "refutation" was a false one — both needed node-level dissection.)

---

## 6. Constraints & mantra (owner)

- **Mantra:** design & derive → plan → execute. One careful step at a time. Pause & re-evaluate after each.
- **Elegance/simplicity > raw accuracy.** No magic numbers / no lingering knobs (env-flag A/B is OK *temporarily*,
  removed on decision). ≤ ~25 calibration modules, ≤ ~8 constants.
- **Gates:** golden byte-identity for refactors; `gdna_none` phantom guard (hard) + node-level honesty for
  behavioural steps. Never land on the aggregate benchmark alone.
- **Develop on toys** (`toy_prod.py`), validate on the full suite. **Owner drives commits + sequencing.**
- **Env:** all build/test/lint inside the activated `rigel` conda env; `OMP_NUM_THREADS=1`. After C++ changes:
  `pip install --no-build-isolation -e .` (no C++ touched this session).

---

## 7. The rest of the roadmap (after the boundary rule)

Per [`pass0_roadmap.md`](pass0_roadmap.md) §4/§7: finish the solve-gate (boundary rule → prior-as-source →
correlation metric → re-run), then **I** correctness items (discretization-frame `+1` vs continuous; the AMBIG
two-root; the identifiability wall), then the feature levers **B** (mature capture-scale correction — the mode-flip
unblocker) → **C** (retry the mode flip) → **D** (TSS/TES seam lower anchor) → **E** (precision merge), then **J**
(single-exon holdout) + the **Phase-2 hyperprior** (the real under-call lever, fed an honestly-solved pass-0).
