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

**A. The solve-gate (§6B). ⛔ REFUTED — do not re-attempt** ([`solve_gate_design.md`](solve_gate_design.md)).
Derived, implemented (flag-gated), and measured: the DOF skip (unidentified → keep `f_g=1` init, defer to the
prior) regresses **both** standalone (refit=0 +0.010, `gdna_none` rises) **and with the hyperprior** (refit=1
+0.025, 15/32 worse, 0 better). The premise is empirically false — the prior resolves an imperfectly-**solved**
node (which carries the observed density) far better than a deferred all-gDNA init. Reverted; `solvable` stays
structural. **"Destination decides" is retired as a pass-0 mechanism** — honest emission (pr→0 when vacuous)
already prevents over-commitment.

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

---

## 5. The handoff — the Phase-2 hyperprior is THE lever (measured)

**The Phase-2 gDNA hyperprior** is the real under-call lever, and it is **already wired** (`calib_refit_iters`):
on `ambig_dense_10mb`, refit=1 alone takes the pass-0 mwae **0.2030 → 0.0998 — it halves the error**
([`solve_gate_design.md`](solve_gate_design.md) §4). Crucially it works **best with pass-0 solving every node**
(the solve-gate deferral made it *worse*, §4-A). So pass-0's contract is simpler than §6B assumed: **solve
honestly (imperfect-but-informative), keep it weak where the evidence is weak, and hand the prior a solved belief
— not a deferral.** The remaining lever is improving the **hyperprior itself** (the live DNA-prior track,
`dna_prior_session_resume.md`), plus withholding single-exon transcripts from its training (§6-J).

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

## 7. Sequencing

1. **Harden first — ✅ DONE.** F (nan sweep) + G (invariant tests) + H (determinism) landed; the solver is now
   robust, principled, and bit-reproducible.
2. ~~Finish the clean architecture — A (the solve-gate)~~ **⛔ REFUTED** (§4-A): the DOF skip regresses even with
   the prior. "Destination decides" is retired; honest emission already prevents over-commitment.
3. **The lever is the hyperprior** (§5) — measured to halve the pass-0 error, and it wants pass-0 to *solve*, not
   defer. The substantive next step is improving the **hyperprior itself** (the DNA-prior track).
4. **The pass-0 correctness/feature items** remain useful but are second-order to the hyperprior: I (discretization;
   the AMBIG two-root — now WITHOUT a solve-gate crutch, so it needs its own stance; the wall), then B (mature
   correction) → C (mode flip) → D (seam anchor) → E (precision merge).
5. **Hand off** — J (single-exon holdout) + the Phase-2 hyperprior (§5), fed an honestly-*solved* pass-0.

**The one rule throughout:** golden byte-identity for refactors; `gdna_none` + node-level-honesty for behavioural
steps; the error-budget map (§3) as the scoreboard. Never land on the aggregate benchmark alone.
