# Message-System Implementation Plan — fix the τ-gate bug + harden the code

**Status (2026-07-21):** the meticulous, code-grounded plan to (1) fix the τ-gate bug, (2) restore `I_spliced`
as a first-class identifiability channel, (3) move the identifiability logic from **emission** to the **solve**
decision (the owner's reframe: *emission is never dangerous; an unidentifiable node solving itself without a
prior is*), and (4) refactor `bp_solver._scan` so this class of bug cannot hide again. Derivation +
message-format context: [`message_system_derivation.md`](message_system_derivation.md). Owner's mental model:
[`intergenic_boundary_behavior.md`](intergenic_boundary_behavior.md). Prior spliced derivation (do not
re-derive): [`spliced_precision_status.md`](spliced_precision_status.md).

---

## 0. Provenance (so we know what we are changing)

- The **τ-emission gate does not ship.** Added 2026-07-20 (`b7d0300e`, "τ-precision phantom fix") on this WIP
  branch; `git tag --contains b7d0300e` is empty. **v0.7.1 (2026-07-12) does not have it.** This is entirely a
  WIP experiment — no released behavior is at risk.
- The **spliced credit was strangled, not lost.** `b7d0300e` made `_ISPLICED_PRECISION` always-on **and** added
  the τ-emission gate in the same commit; the credit (`bp_solver.py:630`) is now unreachable when `τ=0`.
- **Measured impact of the bug** (`scratchpad/msg_gate_verify.py`, `gdna300_ss0.50_capON_nrna`): **52% of the 764
  spliced-carrying boundaries emit nothing** (`τ_src=0 ⇒ emit_g=emit_p=False`); junction 1054 (96 spliced reads)
  → exon 1055 emits `pr=0`.

---

## 1. The organizing concept — the identifiability model

A node's composition (`λ = logit f_g`, and the tilt `θ`) is identifiable from exactly three **evidence
channels**, plus relayed inflow. Making these explicit is the refactor's spine.

| channel | what it is | present when | precision it carries |
|---|---|---|---|
| **`I_strand`** | the unspliced strand-tilt Fisher info (differential-κ deadband, `bp_solver.py:439-443`) | κ≠½ **and** counts present ⇒ 0 on unstranded | `∝ N·((κ−½)²−σ²_d)` — the composition (λ,θ) precision |
| **`I_spliced`** *(restore)* | the motif-stranded **mature-RNA measurement** (a count that legitimately measures RNA — `spliced_precision_status.md` §3) | a boundary carries spliced mass `SPs+SNs > 0` | the RNA-message precision credit `S_eff/(1+S_eff·s2t)` (`bp_solver.py:630`), **and** marks the node non-vacuous |
| **`I_struct`** | structural composition-**certainty** by signature | a node is **pure gDNA by structure** (see §1.1) | a *certainty* (`ev=0` — zero composition variance), not a magnitude; the message precision is then the honest count + `s2t` only |

**Local identifiability:** `identifiable = (I_strand > 0) or (I_spliced > 0) or I_struct`. A node with none of
these, and no equal-gate composition inflow, is **unidentifiable** — it must NOT solve itself in pass-0.

### 1.1 `I_struct` answered (the owner's questions)
- **Purpose.** Some nodes' composition is *known from structure*, needing no strand/spliced evidence: an
  intergenic region is pure gDNA; a TSS/TES (intergenic↔exon) seam is pure gDNA (no RNA crosses a transcript
  start/end — `enrichment_sensitivity_worklog.md` §8c: 430 such seams, oracle `f_g=1.00` exactly). `I_struct`
  marks these as composition-certain so they emit a confident **gDNA** message.
- **When it fires — TODAY vs. TARGET.** Today: `struct_lock = locked & is_region` (`bp_solver.py:452`) — **only
  intergenic REGIONS**; TSS/TES boundary seams are deliberately excluded to avoid the exon↔exon phantom (a
  seam "between RNA-carrying exons" whose crossing is RNA-contaminated). **The exclusion is too coarse:**
  intergenic↔exon seams are genuinely pure gDNA (safe); exon↔exon seams are not. **Target:** `I_struct` fires at
  intergenic REGIONS **and intergenic↔exon (TSS/TES) seams**, never at exon↔exon seams.
- **How its precision is computed.** `I_struct` is not a magnitude — it is `ev_composition = 0` (zero variance on
  the composition axis, i.e. certainty). The emitted message precision is then governed **only** by the honest
  count `1/M_src` + `σ²_transfer` (the crossing is thin + capture-depleted ⇒ genuinely weak, as it should be).
  So a pure-gDNA seam sends a *correctly-signed* (gDNA) but *weak* message — §6 of the derivation.

---

## 2. What is broken, precisely (the two coupled defects)

1. **The emission gate discards real evidence.** `emit_g = (sm>0) and lam_ev`, `emit_p = … and (lam_ev or th_ev)`,
   `lam_ev = struct_lock or (τ>0)` (`bp_solver.py:528-548`). On unstranded data `I_strand=0`, and **`I_spliced`
   is not in τ**, so a spliced junction has `τ=0 ⇒ gated ⇒ the spliced credit (line 630, inside the gated block)
   never runs.** The measurement can't open the gate it would then fund.
2. **The `ev_lam=0`-at-`τ=0` quirk makes the gate load-bearing.** `ev_lam = 0.0 if lock_s else (1/τ if τ>0 else
   0.0)` (`bp_solver.py:529`) — when `τ=0` and unlocked, `ev_lam=0` *means certain*, which is backwards (no
   evidence ⇒ **∞** variance). The boolean gate exists to suppress that wrong value. **So we cannot simply delete
   the gate; we must first fix `τ=0 ⇒ ev_lam=∞`.**

---

## 3. The target design

- **Emission: always on.** Every node emits its per-component message. Precision is sourced from the evidence
  (`τ` + the spliced credit + `I_struct` certainty), so a genuinely vacuous source emits `pr≈0` — *bit-identical
  in effect to not emitting*, so nothing cascades. The **phantom fix stays** (it was the τ-*precision*, not the
  gate; §2.2 of the derivation).
- **The τ=0 fix:** `τ=0, unlocked ⇒ ev_lam = ∞ ⇒ pr → 0` (numerically a large variance / `pr=0`). Then the
  boolean emission gate is **removed**.
- **`I_spliced` restored:** the spliced credit fires whenever `SPs+SNs>0` (no gate), and `I_spliced>0` marks the
  junction identifiable. Mature stays a **component, not a vote** — mode carries `ρ_mature`, precision `S_eff`;
  it does **not** inflate the λ (gDNA-vs-nascent) composition precision (that would revive the refuted
  RNA-vote). (`spliced_precision_status.md` §3/§6.)
- **The solve gate (the owner's #2):** a node **solves iff identifiable**; else it **skips** — keeps its
  signature-binary init (`f_g=1`, max variance, ARCHITECTURE §3) and waits for the pass-1 gDNA prior. This is
  where the danger actually lives, and where the τ/identifiability logic belongs — *not* on emission.

---

## 4. Phased implementation (each phase test-gated; behavior changes measured)

**Gates for EVERY phase:** the 234 calibration unit tests + golden (`pytest tests/calibration`, `--update-golden`
only for *intended* output changes); the **`gdna_none` false-positive guard** (the phantom must not return); the
single-strand under-call (`scratchpad/pass0_ss_dissect.py`); net fragment-flow (`oracle_and_benchmarking.md`).
Toys first (owner directive). `OMP_NUM_THREADS=1`.

- **Phase A — behavior-preserving refactor (health first; byte-identical).**
  - **✅ Regression tests LANDED + FALSIFIER-VERIFIED (2026-07-21, `tests/calibration/test_bp_solver.py`).**
    `test_tau_gag_fix_spliced_junction_emits_when_unstranded` (a spliced junction delivers +RNA to its exon at
    κ=½; the no-spliced control delivers zero — no phantom) and
    `test_tau_gag_fix_deconvolution_prediction_stays_gated`. Proven falsifiers: under the old gag the first
    fails with `prec_p==0.0` (the exact bug); with the fix both pass. Suite 247→**249 pass**, lint clean. These
    lock the fix so the bug class cannot silently return.
  - **TODO (the structural extraction):** pull the identifiability compiler + per-edge message construction out
    of the `_scan` monolith into named, testable units — `Evidence(I_strand, I_spliced, I_struct)` and a
    `ComponentMessage(density, precision, gate)` — replacing the inline `tau0_lam`/`struct_lock` tangle and the
    six parallel `amg/apg/…` arrays. *Gate: golden + all tests byte-identical.* Best done as its own reviewable
    change, separate from the behavioral fix above.
- **Phase B — the τ-gate fix (the "huge win", measured).** **✅ STEP B1 LANDED + VALIDATED (2026-07-21).**
  The minimal, safe first cut: the **spliced RNA MEASUREMENT is independent evidence — always emitted; the
  deconvolution PREDICTION stays τ-gated** (`bp_solver._scan`: `emit_p/emit_n` open on `SPs/SNs>0`; the
  `n_eff` prediction precision is gated by `has_comp_ev = lam_ev or th_ev`, the spliced credit added
  unconditionally). This unblocks the 52% *without touching the vacuous-node path* — so the phantom mechanism is
  structurally untouched. **Measured** (`scratchpad/{msg_gate_verify,pass0_ss_dissect,gdna_none_fp}.py`):

  | metric | BEFORE (gagged) | AFTER (fix) |
  |---|---|---|
  | spliced junctions emitting nothing | 52% | **19%** |
  | gdna300 enriched single-strand mean f_g (oracle 0.909) | 0.564 | **0.691** |
  | **gdna_none ss0.50 false-gDNA** (zero-gDNA library — the phantom guard) | 31.6% | **10.2%** (~3× better) |
  | gdna_none ss0.99 (stranded control) | 0.21% | 0.02% (neutral) |

  The fix **removes** false gDNA (the spliced RNA messages correctly pull transcribed nodes toward RNA) AND
  recovers real gDNA (the relay reconnects). 247 calibration units pass; lint clean; 6 golden files shift (tiny
  false-gDNA values dropping — the intended effect; regen on land). **Remaining in Phase B:** the `ev_lam=∞`
  quirk fix + de-gating the gDNA side + `I_spliced` as an explicit τ/identifiability channel (folds into the
  Phase A refactor). **Add the regression tests** (spliced junction emits; vacuous emits pr≈0).

  **✅ FULL ambig_dense_10mb BENCHMARK (32 scenarios, pass-0 mass-weighted f_g error vs oracle, OLD vs NEW):**

  | group | n | OLD (gagged) | NEW (fix) | Δ | better / worse |
  |---|---|---|---|---|---|
  | ALL | 32 | 0.2444 | 0.2033 | **−17%** | **22 / 0** |
  | UNSTRANDED (ss0.50) | 20 | 0.3742 | 0.3097 | **−17%** | **20 / 0** |
  | STRANDED (ss0.99) | 12 | 0.0280 | 0.0260 | −7% | 2 / 0 (neutral) |
  | unstranded, gDNA-rich | 9 | 0.2311 | 0.2006 | −13% | 9 / 0 |
  | unstranded, gDNA-poor/none | 11 | 0.4913 | 0.3989 | −19% | 11 / 0 |

  **Zero regressions across all 32.** `scratchpad/pass0_bench.py` (A/B via `RIGEL_TAU_GAG`, now removed from
  `bp_solver`; the JSONs are cached).

  **✅ I_strand DEADBAND AUDIT (the real phantom fix, `bp_solver.py:439-443`, commit `b7d0300e`).** Verified the
  strand-tilt-noise correction is in place and working — `disc = 4·max(0, (κ−½)² − σ²_d)`,
  `σ²_d = ¼·(1/N_rna+ω_r) + ¼·(1/N_gdna+ω_g)`: **unstranded (κ≈0.500) ⇒ disc=0 on ALL 20 scenarios** (I_strand
  killed — the κ≈½ sampling whisper is neutralized), **stranded (κ≈0.010) ⇒ disc≈0.94–0.96** (strand informative),
  **gDNA-free ⇒ σ²_d→∞ ⇒ disc=0** (correctly gated). This is WHY de-gating is safe: the phantom is killed at its
  source (I_strand=0 on unstranded), so messages pass without manufacturing strand-derived confidence.
- **Phase C — the solve/skip gate (the owner's #2).** Route the identifiability quantity to the **solve**
  decision: an unidentifiable node keeps its init (skips the ψ solve), waits for the prior. *Gate: unidentifiable
  nodes no longer over-commit; the under-call and net-flow improve or hold; `gdna_none` holds.*
- **Phase D — `I_struct` at TSS/TES + the unequal-gate lower bound.** Extend `I_struct` to intergenic↔exon seams;
  make the gDNA-only (unequal-gate) message a one-sided **lower anchor**, never a downward crush (derivation §6).
  Fixes the 38% wrong-signed seam messages. *Gate: the TSS/TES-flanked enriched exons stop being crushed below
  the 0.5 floor; `gdna_none` holds.* (This is the larger derivation; sequence it last and carefully — prior
  composition-transfer attempts regressed, `cliff_message_derivation.md` §9.)

---

## 5. Code-health measures (so this bug class cannot recur)

- **Explicit typed evidence + message** (Phase A): `Evidence(I_strand, I_spliced, I_struct)` and
  `ComponentMessage(density, precision, gate)` make "what evidence exists" and "what a message carries" readable
  and unit-testable, instead of six parallel arrays and a boolean gate buried mid-loop.
- **One place decides identifiability**, consumed by both the message precision and the solve/skip gate — no
  duplicated or drifting gate logic.
- **Invariant regression tests** that pin the exact defects we just found:
  1. a spliced junction on unstranded data emits a non-zero RNA message (the τ-gag);
  2. a vacuous unstranded node emits `pr≈0` (the phantom guard, at the unit level);
  3. `gdna_none` produces no false gDNA→RNA (integration).
- **No new magic numbers** (owner directive): `I_spliced` reuses the settled `S_eff/(1+S_eff·s2t)`; `ev_lam=∞`
  is a limit, not a constant; `I_struct` reuses the signature. Any new constant pauses for discussion.

---

## 6. Open questions to resolve DURING implementation (not blockers for A–C)

1. **Does receiving a mature (spliced) measurement make the *recipient exon* identifiable?** The exon then has its
   mature pinned, but the gDNA-vs-nascent residual still needs the prior. Provisional: the exon is **partially**
   identifiable (solve the mature component; leave the gDNA amount to the prior) — but for A–C, the simple rule
   (identifiable if local `I_*>0` OR equal-gate inflow) is enough; the mature-vs-nascent transfer is Phase D+
   (`message_system_derivation.md` open issue #4, the known-hard one).
2. **`ev_lam=∞` numerics.** Represent as a capped large variance vs. an explicit `pr=0` branch — pick the one that
   is cleanest and keeps `_fold_lambda` stable.
3. **`I_struct` relay through unequal gates.** A relayed gDNA-only certainty must not mark a downstream exon
   composition-identifiable (it only bounds gDNA). Ensure the relay carries the **gate**, so an unequal-gate hop
   does not launder gDNA-certainty into composition-certainty.
4. **Belief storage `(λ,θ)` vs. 3-variance** — deferred (derivation open issue #6); not needed for A–D.

---

## 7. Execution order

**A → B → C** now (B is the headline fix; A makes B safe and permanent; C is the owner's #2 and is small once the
identifiability quantity is explicit). **D** after, as its own carefully-measured change. Proceed on A immediately
(byte-identical, zero-risk), then B behind the `gdna_none` gate.
