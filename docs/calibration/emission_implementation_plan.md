# Emission Thread — cautious, stepwise implementation plan

**Status (2026-07-21):** the phased plan to realise "always emit; let the destination decide" in
`bp_solver` — the pivot after the mode flip proved blocked. Derivation:
[`emission_and_precision_derivation.md`](emission_and_precision_derivation.md). Same discipline as the mode plan
([`message_mode_implementation_plan.md`](message_mode_implementation_plan.md)): each step independently
verifiable, temporary env flag, `gdna_none` the hard guard.

---

## 0. The honest shape of this thread (read first)

Grounding the derivation in the code surfaced a **material fact** that reshapes expectations:

> **The emission de-gating is BIT-IDENTICAL in effect.** The boolean gate (`emit_g = sm>0 AND lam_ev`) and the
> `ev_λ=∞` fix make the **same** decision in every case: a composition-vacuous source (`τ=0`, unlocked) is *not
> emitted* by the gate, and emits `pr=0` under `ev_λ=∞` — and a `pr=0` message is **never appended to
> `lam_factors`** ([bp_solver.py:685,726](../../src/rigel/calibration/bp_solver.py)), so it is never folded.
> Both append-decisions agree for locked, real-evidence, and vacuous sources alike (§2 of the derivation).

So this thread splits into two very different pieces:

| piece | what it is | risk | direct impact |
|---|---|---|---|
| **`ev_λ=∞` + de-gate** (§2–§3 deriv.) | a **bit-identical refactor** — honest architecture, kills the τ-gag *bug class* (no more all-or-nothing gate that can silence a real channel) | **low** (golden byte-identical) | **none** (by construction) |
| **the solve-gate** (§5 deriv., §6B) | **behavioural** — an unidentifiable node *skips* its pass-0 solve (keeps the honest signature init) instead of committing to a weak/wrong value | medium | pays off **with the hyperprior** |

**And the part that is NOT in this thread but IS the real under-call lever:** the **Phase-2 gDNA hyperprior**
(the DNA-prior track). The enriched-exon under-call is **80.5 % an identifiability floor**
([`density_imputation_precision.md`](density_imputation_precision.md)) — a vacuous gDNA message is *genuinely*
zero-information (`pr=0` is *correct*, not a gate artifact), so no emission/precision change creates the missing
information. The hyperprior is the designed resolver; the solve-gate's job is to hand it an *honest weak* belief
(skip → wide posterior) instead of a confident-wrong one.

**⟹ Net:** this thread is a **safe architectural cleanup + a solve-gate that sets up the hyperprior.** It will
*not* move the benchmark much on its own (the ≤1.9 % confidently-wrong ceiling). Its value is honesty,
bug-class elimination, and unblocking the hyperprior — stated plainly so we don't over-expect.

---

## 1. Verification assets (reused)

- **`gdna_none` phantom guard** — the hard gate (zero-gDNA scenarios must not gain false gDNA).
- **`pass0_bench.py`** on the cached `ambig_dense_10mb` substrate (32 scenarios, mwae vs oracle) — the A/B, via a
  temporary `RIGEL_EMISSION` env flag (levels below), the same pattern as `RIGEL_GATE_SHIFT`.
- **`mode_verify.py`** `_capture` hook (`_edge_modes` carries per-edge `prec_g/prec_p/prec_n`) — to probe the
  emitted precision on vacuous nodes (§Stage 1) with no new code.
- **Golden byte-identity** — the gate on the refactor stages.
- Flag: `RIGEL_EMISSION ∈ {off, evinf, degate, dofgate}`, additive; `off` = today.

---

## 2. The stages

### Stage 1 + 2 — Implemented; NOT bit-identical → a real coupling bug found & fixed ✅ (pending golden regen)
Landed the principled form directly (no flag — a refactor should read clean): the `ev_λ=∞` three-state variance
(lock=0 / evidence=jac²/τ / **no-evidence=∞**, set *directly* — nan-safe at the λ-window edge), one
`_pred_precision(count, v_log, s2t)` helper (0 when unseen or countless — absorbs **two** latent nans the
`has_comp_ev` gate had masked), and **structural-only emission** (no evidence gate).

**The make-or-break came back NOT bit-identical — and the discrepancy is exactly the bug we were hunting.** One
golden test shifts (`antisense_overlap`, the sparse − / +-overlap case); 20/21 golden byte-identical. Traced
(git-stash A/B + a node-14 probe): the OLD emission gate was **coupled to the RNA-total density accumulation** —
`rho_r += rho_neg` lives *inside* `if emit_n`, so gating a **structurally-present but composition-vacuous** strand
(a spliced junction, τ=0, at node 12) silently **dropped that strand's density from the RNA-total message**. The
density (the *mode*) was being suppressed by an *evidence* gate — conflating mode and precision, the very thing
the count-zero-info principle forbids. De-gating decouples them: the RNA-total density now reflects all
structurally-present RNA; τ still governs only the *precision*.

**Validated** (`pass0_bench` on `ambig_dense_10mb`, new-emission vs stashed-old): mean mwae **0.2033 → 0.2030**,
**`gdna_none` phantom guard HELD** (every zero-gDNA scenario flat or improved toward 0), **worst regression
+0.0009** (< the 0.002 gate), 3 scenarios improve > 0.002. So the fix is *more principled AND slightly better*.
**Decision (owner):** accept the behavioural change + regenerate the one golden. Open subtlety for the precision
thread proper: the RNA-total *mode* (density, both strands) and its *precision* (here a +spliced measurement) have
different provenances — decoupling is correct, but their merge is §4's deferred question.

### Stage 2 — Land `ev_λ=∞` + de-gate as a refactor *(byte-identical; golden)*
Replace the `ev_λ=0`-at-`τ=0` quirk with the `∞` limit ([bp_solver.py:587](../../src/rigel/calibration/bp_solver.py)),
remove the boolean `emit_g/emit_p/emit_n` gates and the `has_comp_ev` prediction-precision gate
([:613-616](../../src/rigel/calibration/bp_solver.py)) — every source always emits; the precision self-zeros on
vacuous sources. Keep the per-strand **structural** `free_s` continuity (that gates *which strand physically
crosses*, not evidence) and the spliced-measurement credit. **Gate:** golden + `pass0_bench` byte-identical
(proven in Stage 1); the diff is pure architecture. **Value:** the τ-gag *bug class* is gone — no gate can ever
again silence a real measurement channel.

### Stage 3 — The solve-gate: structural → DOF (§6B) *(behavioural; flag=`dofgate`; measured)*
Replace `solvable = (fp|fn) & mass>0` ([bp_solver.py:398](../../src/rigel/calibration/bp_solver.py)) with the DOF
criterion: solvable iff **every free axis** (`λ`; `θ` for AMBIG) has ≥1 nonzero-precision source among {strand,
gDNA message, per-strand RNA message, prior}. An unidentifiable node **skips** the ψ solve and keeps its
signature init (`f_g=1`, max variance). Reuses the `StrandEvidence` compiler. **Gates:**
1. **`gdna_none` non-regression** (hard).
2. **`pass0_bench` A/B** — report the split: does the confidently-wrong slice shrink? Measure **both** pass-0-only
   **and with `calib_refit_iters=1`** (the hyperprior payoff — skip → the prior resolves it). Expect the pass-0-
   only number to be roughly neutral (skip vs weak-solve) and the *with-refit* number to improve where the
   hyperprior can lift the now-honest weak nodes.
3. **No confident-wrong regressions** — an unidentifiable node must not end more committed than before.

### Stage 4 — Land + clean up
Land the decided defaults (likely: `ev_λ=∞`+de-gate always on — it's byte-identical; the solve-gate per Stage-3
evidence), **delete `RIGEL_EMISSION`** (no lingering knob), regenerate goldens for any intended change, update the
derivation + [`message_system_derivation.md`](message_system_derivation.md) §7/§6B to "landed."

---

## 3. Deferred (explicitly not this thread)

- **The Phase-2 gDNA hyperprior** — the real under-call lever (the 80.5 % identifiability lifting); the solve-gate
  sets it up. The live DNA-prior track ([[dna_prior_projection_resume]]).
- **The mature capture-scale correction** (`message_propagation_arithmetic.md` §6a #3) — the *mode-flip* blocker,
  independent of emission.
- **The C4 mature≫nascent conditioning** (derivation §4).

---

## 4. Execution order & the one rule

**S1 → (STOP: bit-identical?) → S2 → S3 → (review) → S4.** S1–S2 are byte-identical and low-risk (do them back to
back once S1 proves bit-identity). S3 is the only behavioural change and its payoff is coupled to the hyperprior —
measure it **with and without refit** and review before landing.

**The one rule (unchanged):** every step gates on golden byte-identity (refactor) or `gdna_none` +
node-level-honesty (behavioural). Never land on the benchmark aggregate alone — here especially, since the
aggregate is ceiling-bounded and the real payoff routes through the hyperprior.
