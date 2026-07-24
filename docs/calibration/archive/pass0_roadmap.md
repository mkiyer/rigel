# Pass-0 Roadmap — finish the message system, then harden the solver

**Status (2026-07-21):** the consolidated view of what pass-0 is, what has landed, and what remains — organized so
we **tighten and harden confidence in the solver before moving on to the Phase-2 hyperprior**. Pass-0 is the
*prior-free* belief-propagation sweep that deconvolves each node's unspliced mass into `(f_g, f_+, f_−)` from two
sources only: the **strand likelihood** (intrinsic) and **cross-node messages** (imputation). The third source —
the **gDNA hyperprior** — enters in pass-2 and is out of scope here (it is the real under-call lever; §5).

---

## 1. The organizing frame — a message has four concerns + a destination decision

Every node-to-node message carries, per component `c ∈ {gDNA, +RNA, −RNA}`:

| concern | question | governed by | doc |
|---|---|---|---|
| **MODE** (value) | what composition target? | density currency + eff-length frame | `message_propagation_arithmetic.md` |
| **PRECISION** (confidence) | how sure? | reference-free evidence `τ` + count + `σ²_transfer` | `emission_and_precision_derivation.md`, `message_precision_derivation.md` |
| **EMISSION** (whether sent) | send it at all? | **always** — `pr=0` is ignored | `emission_and_precision_derivation.md` §3 |
| **GATE** (which components) | which strands cross? | **structural** (`free_s` continuity) | — |
| **SOLVE** (destination) | can I solve, or defer? | DOF criterion (§6B) | `message_system_derivation.md` §6B |

**The currency is DENSITY** (frame-invariant); the composition fraction is a frame-dependent projection of it. The
count enters only as statistical power, never as a composition vote (count-zero-info, `CALIBRATION_ARCHITECTURE.md`).

---

## 2. State of each concern

| concern | state | what remains |
|---|---|---|
| **MODE** | derived + **MC-validated**; code proven to match bit-for-bit (S1); helpers extracted + tested (S2); flip set characterized (S3) | the flip itself (S4-6) is **blocked** on the mature cliff (§4-B); the unequal-gate seam anchor (§4-D) |
| **EMISSION** | **landed** — always-emit, `ev_λ=∞`, one `_pred_precision` helper; fixed the emission↔density coupling bug; retired the τ-gag bug class | — (the solve-gate is its destination-side twin, §4-A) |
| **PRECISION** | honest 3-term model `1/(Var + 1/M + σ²_transfer)` derived; `Var` from reference-free `τ`; ceiling measured (§3) | the prediction⊕measurement **merge provenance** (§4-E); C4 conditioning; `σ²_transfer` at unequal gates |
| **GATE** | structural `free_s`, stable | — |
| **SOLVE** | structural `solvable=(fp\|fn)&mass>0` — over-commits unidentifiable nodes | upgrade to the **DOF gate** (§4-A) |

---

## 3. The honest error budget (why pass-0 alone can't finish the under-call)

Measured over all 32 `ambig_dense_10mb` scenarios (`density_imputation_precision.md`):
- **80.5 % under-determined** — an *identifiability floor* (unstranded, low-gDNA). **No message precision or mode
  can fix it.** The **gDNA hyperprior resolves it** (§5).
- **1.9 % confidently-wrong** — the only slice message precision can attack.

**⟹ Pass-0's job is HONESTY, not a big benchmark number:** hand pass-2 a belief that is *weak* where
unidentifiable and *never confident-wrong*. Everything below is judged against that — and against **not
regressing the `gdna_none` phantom guard**, the hard gate.

---

## 4. Remaining feature work (each derive → plan → stage, A/B-gated)

> **BOUNDARY MODE — LANDED (2026-07-22, `boundary_rule_rederivation.md`).** The whole boundary/σ²_transfer/mode
> thread below (items A boundary-rule, B, C) is **superseded** by two shipped fixes: **Fix #1** the directional
> spliced-density σ²_transfer (§12–13 there — the precision bug: spliced excluded from the transfer variance) and
> **Fix #2** the geo-mean cliff-interpolated mature-free crossing mode (§14 — count-legal, zero constants). The
> antisense-intronic leak is fixed; boundary corr 0.105→0.233. The mature-gate-dismantle's exon→intron *direct*
> nascent emission is superseded by the residual relay (change 1). **RETIRED:** `RIGEL_GATE_SHIFT` (the mode-flip
> A/B — the geo-mode owns boundary edges now). **NEW open item (item L below):** the **nascent≫mature** over-call
> (the unstranded identifiability floor, not the boundary mode) — the quintuple grid is its regression harness.

**A. The solve-gate (§6B). ACTIVE — validated for regions; boundary rule + metric to fix**
([`solve_gate_design.md`](solve_gate_design.md)). A node-type CORRELATION test (not the mass-weighted-error
metric, which wrongly counts withheld nodes' arbitrary default) confirms: the DOF criterion is **correct for
regions (77 % of mass)** — the single-strand + AMBIG regions it withholds are genuine coin-flips (corr ≈ 0),
those it solves correlate 0.63–0.69. A first "refutation" was a **flawed measurement** (mwae counted arbitrary
defaults; the gate omitted the prior as a source so the refit re-skipped instead of resolving). **To ship:** (1) a
correlation/precision-aware metric; (2) count the fitted prior as an identification source (so the refit resolves
deferred nodes); (3) **re-derive the boundary rule** (currently inverted — a real bug). Then re-run ON/OFF.

**B. The mature capture-scale correction (§6a #3) — the mode-flip UNBLOCKER.** Under proportional capture, the
mature (spliced, one-sided) channel is captured at a different rate than the unspliced channel, so the mature
add/subtract carries a bias ∝ mature-fraction. This is *why* the Stage-4 exon shift regressed. Needs its own
derivation (the capture-frame ratio) + implementation.

**C. Retry the mode flip (S4-6) — after B.** With the mature bias corrected, re-run the `RIGEL_GATE_SHIFT` A/B;
land `use_shift = gates_equal` if it clears Gate 3 (class A is empirically empty, so S6 is a confirm-and-land).

**D. The unequal-gate lower anchor (TSS/TES seams) — the high-impact seam fix.** The intergenic↔exon seam is
gate-unequal: the gDNA density transfers only as a **one-sided lower anchor** (§5/§6 of `message_system_derivation`),
never the crushing density point-estimate. Needs the ψ form (half-quadratic / hinge) — precision-adjacent.

**E. Precision refinements.** The prediction⊕measurement **merge**: the antisense bug surfaced that the RNA-total
*mode* (both-strand density) and its *precision* (a +spliced measurement) have different provenances — decoupling
was right, but the principled merge is open. Plus C4 conditioning (mature≫nascent) and `σ²_transfer` at unequal
gates.

**L. The nascent≫mature over-call (NEW — the next investigation).** On unstranded data, when nascent RNA far
outweighs mature, the solver **over-calls gDNA** (middle-exon f_g 0.81 vs oracle 0.31; 0.69 vs 0.05 near zero
gDNA — `boundary_rule_rederivation.md` §14 interrogation 2). Nascent is unspliced ⇒ locally indistinguishable
from gDNA; the only RNA proof is the scarce spliced/mature. This is the **unstranded identifiability floor**, NOT
the boundary mode (the geo-mode gives identical results). The lever: the **intron→exon nascent relay** (the
quintuple `intron ↔ IE ↔ exon ↔ EI ↔ intron` — the introns are transcribed and well-defined; the exon is a relay
whose gDNA is imputed) or the global gDNA prior. Regression harness: `scripts/debug/quintuple_grid.py` (sweep
gDNA × nascent × mature; the solver must "just work" across the grid).

---

## 5. The handoff — the Phase-2 hyperprior is THE lever (measured)

**The Phase-2 gDNA hyperprior** is the real under-call lever, and it is **already wired** (`calib_refit_iters`):
on `ambig_dense_10mb`, refit=1 (with today's solve-everything pass-0) takes the mwae **0.2030 → 0.0998 — it
halves the error**. *(Whether a properly-fixed solve-gate then does better is OPEN — the earlier "solving beats
deferring" claim came from a broken measurement, `solve_gate_design.md` §2; it must be re-tested with the
prior-as-source fix.)* Pass-0's contract to the prior: solve what it can, honestly **withhold** what it cannot
(the coin-flip nodes, §4-A), and hand the prior a **fit trained on solved nodes only** (already the case —
`_fit_gdna_hyperprior` excludes AMBIG/boundary). Plus withholding single-exon transcripts from its training (§6-J).

---

## 6. HARDENING — tighten confidence in the solver (the near-term priority)

Before investing in §5, make the solver *trustworthy*. These are ordered by confidence-risk.

- **F. Numerical robustness sweep. ✅ DONE** (`test_bp_solver.py::test_message_primitives_never_nan`,
  `::test_sweep_finite_over_extreme_configs`). Fuzzed `_pred_precision`/`_mode_shift`/`_mode_density` over extreme
  inputs (incl. `v=∞` AND `v=nan` — `_pred_precision` self-guards) + the real sweep over 20 configs
  (pure-gDNA / pure-RNA / empty / tiny / huge). **No further latent nan found** beyond the two the refactor
  already fixed; every emitted mode/precision + final belief is finite (∞-variance = honest "unsolved" allowed).
- **G. Invariant test coverage. ✅ DONE** (`::test_pred_precision_honest_semantics`,
  `::test_vacuous_unstranded_source_zero_precision_but_density_flows`; plus the existing `test_tau_gag_fix_*`,
  `test_compile_strand_evidence_deadband_*`). The principles are now pinned directly: no-evidence⇒pr=0 (ev_λ=∞);
  monotone-in-count; and the **density ⊥ evidence decoupling** — a vacuous source emits 0 composition precision
  yet a well-defined density mode (the exact principle the coupling bug violated). A future refactor now fails a
  *named principle*, not just "the golden moved."
- **H. Determinism / reproducibility. ✅ DONE.** Pass-0 is bit-reproducible at `OMP_NUM_THREADS=1`, both
  within-process (`test_bp_solver.py::test_node_sweep_deterministic` — belief + every emitted message
  bit-identical run-to-run) and **cross-process** (two independent `calibrate` runs → byte-identical across all
  32 `ambig_dense_10mb` scenarios). The sweep is sequential Python with no parallel reduction; the historical
  ~2.6 % wander was C++ scan/EM FP-reduction order, absent at single-thread — and the property survived the whole
  mode+emission refactor.
- **I. Known correctness items.**
  - **Discretization-frame consistency** (§6a #2): `region_eff_length` (+1 discrete) vs `spliced/boundary`
    (continuous) — a several-% frame error on short flanks (`L≲fl_mean`). Align the conventions; gate on golden.
  - **The AMBIG two-root problem** (`CALIBRATION_MASTER.md` §4): an unstranded AMBIG node's `λ` is bimodal under
    a gDNA-only message — the solve-gate (§4-A) guards the zero-precision axis, but the *informed-but-bimodal*
    case needs an explicit stance.
  - **The identifiability wall** (`L<fl_mean`, 0/0 density): confirm the honest-imprecision path (emit nothing,
    defer) is what actually fires — it is *by design*, but should be tested, not assumed.
- **J. Single-exon holdout.** Withhold single-exon transcripts from hyperprior training (a fundamental
  identifiability limit on unstranded data). Small, owner-flagged, sets up §5.
- **K. The map itself.** Keep the error-budget decomposition (§3) current as each item lands — it is how we know
  a change moved the *right* slice (confident-wrong ↓) rather than just perturbing the identifiability floor.

---

## 7. Sequencing — the remaining path to SHIP pass-0, then the hyperprior

**Ship criterion for pass-0 (the contract to pass-2, §3):** the solve is *HONEST* — **never confident-wrong**,
*weak* where unidentifiable, and the **`gdna_none` phantom guard held**. Pass-0 does NOT need a big benchmark
number (the hyperprior owns the 80.5 % identifiability floor); it needs to hand the prior a belief that is
trustworthy where it is confident and honestly weak where it is not.

1. **Harden — ✅ DONE.** F (nan sweep) + G (invariant tests) + H (determinism).
2. **Boundary transfer mode — ✅ LANDED (2026-07-22, commit `4f2f5511`).** Fix #1 (spliced-density σ²_transfer) +
   Fix #2 (geo-mean crossing mode). This subsumed the boundary parts of the old items A/B/C: the "boundary rule
   inversion" was a **Simpson metric artifact** (no inversion — the DOF criterion is fine for regions), the mature
   reconciliation at boundaries is handled by the geo-mean's nascent-as-residual, and `RIGEL_GATE_SHIFT` is retired.
3. **⟵ SHIP-BLOCKER — item L, the `nascent ≫ mature` over-call (§4-L). THE NEXT TASK.** This is a
   **confident-wrong** pattern (over-calling gDNA where the mass is nascent RNA), exactly what the ship criterion
   forbids — so it must be fixed OR made honestly weak before pass-0 ships. It is the unstranded
   nascent-vs-gDNA identifiability floor; the lever is the **intron→exon nascent relay** (the quintuple) or the
   global gDNA prior. Harness: `scripts/debug/quintuple_grid.py`. Two acceptable outcomes: (a) the relay resolves
   it (a real fix), or (b) the over-call is driven to LOW precision (honest-weak, deferred to the hyperprior).
4. **The solve-gate — A (§4-A), the honesty refinement.** Decide withhold-vs-solve for the coin-flip nodes; needs
   the correlation/precision metric + the fitted prior as an identification source. Reframed: no boundary
   inversion (that was the artifact); the open part is the region withhold decision + the prior-as-source.
5. **The correctness items — I (§6).** Discretization-frame consistency; the AMBIG two-root stance; the
   identifiability-wall path (test it fires). Smaller hardening.
6. **SHIP pass-0**, then **the hyperprior (§5)** — the real under-call lever, fed the honest pass-0.
   (Lower-priority polish D seam-anchor / E precision-merge can follow the hyperprior.)

**The one rule throughout:** golden byte-identity for refactors; `gdna_none` + node-level-honesty for behavioural
steps; the error-budget map (§3) as the scoreboard. Never land on the aggregate benchmark alone.

---

## 8. CONSOLIDATED ROADMAP (2026-07-22) — supersedes §7 sequencing

§7's ship path was written before the message-arithmetic thread. That thread was not a digression: it found
that pass-0's **mode** arithmetic was structurally wrong across the capture cliff, which is the root cause of
the unstranded+capture-ON failure §7 was trying to route around. This section is the complete open set.

### 8.1 Landed this session
* **Toolkit** — `injected_priors` hook in `calibrate` (population priors: κ, strand overdispersions, noise-floor
  sample sizes, enrichment NPMLE, intron background, ρ_bg) + `scripts/debug/toy_inject.py` + the FULL-TRANSCRIPT
  toy (intergenic-TSS-exon-intron-exon-intron-exon-TES-intergenic). Injection verified byte-identical.
  *A tiny toy cannot fit these priors — this is what makes toy results trustworthy at all.*
* **intron ↔ IE-boundary** → composition invariance (shift). Identical active components (gDNA+nascent).
* **exon ↔ IE-boundary** → shift ± `c_b`, `c_b = log(1+S_B/D_B)` (`exon_boundary_mature_dilution_plan.md`).
  Enrichment-invariant, zero constants, validated on 22 toy cells (premise `mature-cross=0.000` universal).

### 8.2 Message-arithmetic reconciliation (`message_arithmetic_reconciliation.md`) — ACTIVE
| # | item | status |
|---|---|---|
| R1 | **MASS→COUNT**: `n_unspl_left/right` plumbed into `NodeGeometry`; `_pred_precision` uses integer counts. | ✅ **LANDED** (phantom-neutral) |
| R2 | **Decouple the spliced MEASUREMENT channel from σ²_transfer**. ⛔ **BLOCKED on E** — measured +49 % `gdna_none` phantom: `pr += S` attaches MEASUREMENT confidence to the PREDICTION's mode, which mature-absorption clamps to ~0 ⇒ confident "no RNA". Needs the measurement to carry its OWN mode first. | blocked |
| R3 | **Retire `rho_g_cross`** (unweighted pre-scan geo-mean). | ✅ **LANDED** — measured INERT (byte-identical on the phantom guard, the toy, and all suites): once both exon↔boundary directions carry real messages its only remaining destinations were structurally-pinned seams. |
| R4 | **Neutralize the σ²_transfer cliff term** `(μ_dst−μ_src)²` — a proxy for a mode that is now correct; measured to throttle the corrected messages to prec≈0.03. **Behind a switch; per-condition A/B on STRANDED and UNSTRANDED arms** (it currently protects stranded data). | the big lever |
| R5 | **Measure residual disagreement** after R1–R4; decide whether ANY replacement damping is warranted. | after R4 |
| R6 | **Seam mode derivations**: TSS/TES (intergenic↔exon) and exon–exon AMBIG. Components differ by RNA *presence*, not mature. **Completing R6 is what allows the DENSITY MODE to be retired entirely.** | last mode item |

### 8.3 Remaining pass-0 items (carried forward)
* **A′. Solve-gate, reframed.** "Always emit" is implemented. The node-side gate is NOT: the refuted DOF gate
  kept the `f_g=1` init (an all-gDNA lock). The intended paradigm is an **honest precision-0 state** for
  unidentified nodes, which the gDNA hyperprior later resolves. Materially different from what was refuted;
  needs derivation + implementation + the hyperprior handoff contract.
* **E. Precision prediction⊕measurement merge** (§4-E). **PROMOTED — the NEXT task and a PREREQUISITE for R2.**
  **DERIVED** in `prediction_measurement_merge_derivation.md`: the bug is precision-ADDITION (replicate-measurement
  algebra) applied to an ADDITIVE decomposition `ρ_r = ρ_m + ρ_n`; the fix is a share-weighted delta-method
  log-variance whose mode is the SUM with ρ_n floored at 0 — which makes the measured mature a lower bound on
  RNA that the prediction cannot erase (structurally killing the confident-"no RNA" ⇒ phantom pathology).
  σ²_xfer (point-flux → region-containment) applies to BOTH components: nothing here measures the exon directly.
* **L. nascent≫mature over-call** — may be substantially addressed by the mode fixes; **re-measure** on the
  grounded toy before treating it as open.
* **Depletion bias** (probe-attenuated junctions): bounded (~0.053 f_g at 4× depletion), one-sided, contained by
  the intron-side anchor. NOT fixable by any precision term (it is a mode bias). Real fix = cross-junction
  consistency detector (a transcript's junctions share one mature abundance) — non-local; follow-up.

### 8.4 Standing gates (every behavioural stage)
`msg_audit` direction · mature-dilution identity test · **`gdna_none` phantom guard (hard)** · calibration
suite · **benchmark A/B per-condition (stranded AND unstranded)** · goldens regenerated LAST, only after A/B review.

### 8.5 Ship criterion (unchanged)
Pass-0 is HONEST: never confident-wrong, weak where unidentifiable, `gdna_none` guard held. The gDNA hyperprior
owns the residual identifiability floor.
