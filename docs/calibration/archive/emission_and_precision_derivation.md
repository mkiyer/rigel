# Emission & Message Precision — the "always emit; let the destination decide" derivation

**Status (2026-07-21):** the derivation step for the **emission thread** — the pivot after the mode flip proved
blocked ([`message_mode_implementation_plan.md`](message_mode_implementation_plan.md) Stage 4: the exon-edge
shift is emission-limited + mature-cliff-biased). This is the **confidence half** of the message system (the
mode derivation, [`message_propagation_arithmetic.md`](message_propagation_arithmetic.md), was the value half).
It builds on the settled precision work — it does **not** re-derive it:
- the reference-free evidence compiler `τ` ([`message_precision_derivation.md`](message_precision_derivation.md));
- the honest sampling variance `Var(log f_c) + 1/M` and the **error-budget ceiling** ([`density_imputation_precision.md`](density_imputation_precision.md));
- the enrichment transfer variance `σ²_transfer` ([`transfer_variance_formal_derivation.md`](transfer_variance_formal_derivation.md));
- the spliced measurement precision `S_eff/(1+S_eff·σ²_transfer)` ([`spliced_precision_status.md`](spliced_precision_status.md)).

**Goal.** Realise the owner's architecture — *a source always emits; a zero-information message carries zero
precision and is ignored; the destination decides whether it can solve* ([`message_system_derivation.md`](message_system_derivation.md)
§7 + §6B) — by (1) fixing the `ev_lam=∞` sign quirk so a vacuous source emits `pr→0` naturally, (2) removing the
boolean emission gates, (3) upgrading the structural solve-gate to the DOF criterion. **This unblocks the mode
flip** (messages flow, so their mode matters) and hands the Phase-2 hyperprior an *honest weak* belief instead of
a gated silence or a confident-wrong phantom.

---

## 0. What this thread can and cannot move (honesty first)

`density_imputation_precision.md` measured the pass-0 error budget over all 32 `ambig_dense_10mb` scenarios:
**80.5 % is under-determined** (an *identifiability floor* on unstranded low-gDNA nodes — no message precision can
fix it; the gDNA **hyperprior** is designed to resolve it) and only **1.9 % is confidently-wrong** (the sole
precision-attackable slice). **So this thread has a small *direct* benchmark ceiling.** Its value is elsewhere:

1. **Architecture** — "always emit, destination decides" replaces a holdover of ad-hoc gates (one of which was
   the τ-gag bug).
2. **Honesty** — an unidentifiable node ends pass-0 *weak* (high variance), not gated-silent and not
   confident-wrong. That is exactly the input the hyperprior needs (`message_system_derivation.md` §9).
3. **Unblocking** — the mode flip (and any future message work) can only act where a message *flows*; de-gating
   is the prerequisite (Stage-1/3 of the mode thread measured the gate suppressing most gDNA messages).

We hold `gdna_none` (the zero-gDNA phantom guard) as the **hard** gate: de-gating must not resurrect false gDNA.

---

## 1. The message precision (recap — the settled model)

A message on component `c` carries a **log-fraction target** (the mode, §mode-doc) and a **precision**

```
   pr_c = 1 / ( Var(log f_c)  +  1/M_src  +  σ²_transfer )
          └─ composition ──┘   └ count ┘   └ enrichment ┘
```

- **`Var(log f_c)`** — the *composition* uncertainty, sourced from the **reference-free evidence** `τ`, never from
  the belief variance (which contains the shared Beta(½,½) reference and pooled inflow ⇒ manufactured confidence,
  `message_precision_derivation.md` §2). In log-fraction coordinates `Var(log f_g) = (1−f_g)²·ev_λ`,
  `Var(log f_r) = f_g²·ev_λ`, with **`ev_λ ≡ 1/τ_λ`** the evidence variance ([bp_solver.py:587-591](../../src/rigel/calibration/bp_solver.py)).
- **`1/M_src`** — the count sampling term. **Effective length enters here, through the count** `n = ρ·E`: a short
  node collects few fragments ⇒ large `1/M` (there is no separate "length" precision axis —
  `density_imputation_precision.md` takeaway).
- **`σ²_transfer`** — the enrichment-crossing damping (the NPMLE projection variance, `transfer_variance…`).

Two **independent** evidence channels feed one message (owner's model, `spliced_precision_status.md` §3):
the **deconvolution PREDICTION** (the deconvolved unspliced RNA — confidence `τ`) and the **spliced MEASUREMENT**
(a splice junction's motif-stranded fragments — confidence its own count `S_eff`). §4 derives how they merge.

---

## 2. The `ev_λ = ∞` fix (the core of this thread)

**The quirk** ([bp_solver.py:587](../../src/rigel/calibration/bp_solver.py)):
```python
   ev_lam = 0.0 if lock_s else (1.0/tau_lam if tau_lam > 0 else 0.0)     # ← the 0-at-τ=0 branch is BACKWARDS
```

`ev_λ` is the evidence **variance** `1/τ_λ`. Take the two limits:

| source state | `τ_λ` | correct `ev_λ = 1/τ_λ` | meaning | current code |
|---|---|---|---|---|
| **structural lock** (intergenic = pure gDNA) | `∞` | **0** | composition CERTAIN | `0` ✅ |
| real evidence (strand tilt / spliced / relayed) | `>0` | `1/τ_λ` | finite | `1/τ_λ` ✅ |
| **composition-vacuous** (unstranded, unlocked, no relay) | **0** | **∞** | **NO information** | **`0` ❌ — says CERTAIN** |

The `τ=0` branch conflates *"no evidence"* with *"certain"* — the exact inversion. With `ev_λ=0`, a vacuous source
gets `Var(log f_c)=0`, so `pr_c = 1/(0 + 1/M + σ²_transfer) > 0` — a **confident** message about a composition it
cannot see. **That is the phantom.** The boolean emission gate (`lam_ev`) exists only to *suppress this wrong
value* — it is a patch over the sign error, which is why it became load-bearing (and why the τ-gag bug could hide
behind it).

**The fix — the honest limit:**
```python
   ev_lam = 0.0 if lock_s else (1.0/tau_lam if tau_lam > _EPS else INF)   # τ→0 ⇒ ev→∞ ⇒ Var→∞ ⇒ pr→0
```
(numerically `INF` = a large constant, or a `pr=0` short-circuit — chosen for `_fold_lambda` stability). Now a
vacuous source has `Var(log f_c)=∞ ⇒ pr_c = 1/(∞ + …) → 0` **by construction** — it emits a *zero-precision*
message, which the recipient's ψ ignores (`−½·pr·(…)² = 0`). The lock and real-evidence rows are unchanged, so
**the fix is behavioural only on the τ=0-unlocked rows** — precisely the rows the boolean gate was suppressing.

---

## 3. Why the emission gates then become removable (§7)

**Principle** ([`message_system_derivation.md`](message_system_derivation.md) §7): a source *always* emits; a
zero-precision message is *identically ignored* by the recipient, so "not emitting" ≡ "emitting `pr=0`." Emission
may be skipped **only** as a bit-identical performance optimisation.

Once `ev_λ=∞` (§2) makes a vacuous source emit `pr≈0`, the boolean gates `emit_g/emit_p/emit_n` are **redundant**:
they gate exactly the `pr≈0` messages that the recipient would ignore anyway. Removing them is therefore
**bit-identical in effect** — *provided* a genuinely vacuous source really does reach `τ≈0` (§3.1). What removing
them *gains*: a source with **partial** information (e.g. τ=0 but a spliced measurement, or a gDNA lock but no
strand) is no longer all-or-nothing silenced — it emits its real channel at its real precision and zeros the
rest. This is the τ-gag class of bug, structurally eliminated.

### 3.1 The make-or-break: does `τ` actually reach 0 on a vacuous chain?

**The risk** (`spliced_precision_status.md` §7.8, dated *before* the deadband): on a **fitted-κ** unstranded chain,
`τ` did **not** reach 0 — a fitted `κ≈0.500` but not *exactly* ½ seeds a tiny `I_strand`, and the **τ cavity**
relays absorbed precision onward, so a vacuous chain could pool a small `τ>0` into confidence. If `τ` stays `>0`,
de-gating re-opens the phantom.

**Why it is likely now resolved** — the **I_strand deadband** (`b7d0300e`, 2026-07-20, *after* §7.8):
```
   disc = 4·max(0, (κ−½)² − σ²_d),   σ²_d = ¼(1/N_rna + ω_r) + ¼(1/N_gdna + ω_g)
```
kills `I_strand` (`disc=0`) whenever `|κ−½|` is within the sampling noise `√σ²_d` — i.e. on **every** unstranded
library the seed `τ` is *exactly* 0, not "small." And the τ cavity only relays what is **absorbed**; a chain of
`τ=0` sources emits `pr=0`, so nothing is absorbed, so `τ` **stays** 0. The §7.8 precondition should now be
**met**. **This is the first thing the implementation plan measures** (a read-only probe: the emitted `pr` on
vacuous unstranded nodes must already be ≈0 today — the gate is then provably redundant).

---

## 4. The prediction ⊕ measurement merge

A boundary carrying both a deconvolved-nascent PREDICTION and a spliced MEASUREMENT sends **one** RNA message. The
two are **statistically independent** evidence about the same log-fraction target (the prediction from the
unspliced strand/relay evidence; the measurement from a disjoint set of spliced reads). Independent Gaussian
evidence ⇒ **precisions add**:
```
   pr_RNA = pr_pred + pr_meas,     pr_pred = n_eff/(n_eff·(Var+σ²_transfer)+1),  pr_meas = S_eff/(1+S_eff·σ²_transfer)
```
This is the current code ([bp_solver.py](../../src/rigel/calibration/bp_solver.py) `pr += SPs/…`) and it has the
right limits: `S=0 ⇒ pr_meas=0` (prediction only); vacuous prediction (`Var→∞ ⇒ pr_pred→0`, §2) ⇒ measurement
only; both present ⇒ the sum, dominated by whichever is more certain. The mode is the precision-weighted combine
of the two channels' targets (already folded via `ρ = ρ_nasc + ρ_mat`). **Open sub-point (deferred):** the C4
conditioning when mature ≫ nascent (`message_propagation_arithmetic.md` §4b) — the recovered nascent is a
difference of large near-equal densities; its precision should shrink accordingly. Flagged, not blocking.

---

## 5. The destination decides — the solve-gate (§6B)

Emission is the source side; the **solve-gate** is the destination side of "let the destination decide." Today it
is **structural** ([bp_solver.py:398](../../src/rigel/calibration/bp_solver.py)):
```python
   solvable = (fp | fn) & (mass_unspliced > 0)     # "has a strand + has mass"
```
This lets an **unstranded node with mass but no incoming information** solve anyway — landing on the reference /
whatever weak inputs exist (the over-commitment the deadband + gates currently paper over). The **DOF criterion**
([`message_system_derivation.md`](message_system_derivation.md) §6B) is the honest replacement: a node is solvable
iff **every free axis** (`λ`; and `θ` for AMBIG) has ≥1 nonzero-precision source among {strand, a gDNA message, a
per-strand RNA message, the prior}. Equivalently — now that emission is honest (§2) — **skip the pass-0 solve iff
the total incoming precision on a free axis is 0**, and keep the signature-binary init (`f_g=1`, max variance).
The node then defers to the Phase-2 hyperprior with an *honest wide* posterior — no pass-0 phantom to un-learn.

This composes cleanly with §2–§3: the same `τ`/evidence quantities that zero a vacuous source's *emitted*
precision also zero its *self-solve* precision — one identifiability compiler (`StrandEvidence`, already extracted),
two consumers (send-precision and the solve-gate).

---

## 6. The full picture & the deferred edge

```
 SOURCE  (always emits, per component)          RECIPIENT / SELF
   mode   = §mode-doc (value)                    solve iff a free axis has precision > 0  (§5, DOF)
   pr_c   = 1/(Var(log f_c) + 1/M + σ²_transfer) else: keep signature init, defer to hyperprior (honest weak)
     Var  = (jac)²·ev_λ,  ev_λ = 1/τ_λ           fold each message: −½·pr·(log f_c(λ) − mode)²  (pr=0 ⇒ ignored)
     τ=0 ⇒ ev_λ=∞ ⇒ pr→0   (§2, the fix)
   pr_RNA = pr_pred ⊕ pr_meas  (§4, add)
   NO boolean emission gate    (§3, removed)
```

**Deferred (not this thread):** the mature capture-scale correction (`message_propagation_arithmetic.md` §6a #3 —
the *other* mode-flip blocker); the C4 mature≫nascent conditioning (§4); the Phase-2 hyperprior itself (the DNA-
prior track — the 80.5 % identifiability lifting).

---

## 7. What the implementation must validate (the gates)

1. **The `τ≈0` probe (§3.1)** — read-only: today, on vacuous unstranded nodes, is the *emitted* `pr` already ≈0
   (the gate redundant)? If yes, de-gating is bit-identical-in-effect and safe. If `τ>0` leaks, that is the real
   work (fix the leak before de-gating). **Make-or-break; do it first.**
2. **`gdna_none` phantom guard** — the hard gate on every step. De-gating + `ev_λ=∞` must not increase false gDNA.
3. **The mode flip re-measured** — with messages flowing, re-run the Stage-4 `RIGEL_GATE_SHIFT` A/B: does the
   B-safe shift's regression change once the gDNA messages actually emit? (Diagnostic, not a gate.)
4. **The benchmark** — expect a *small* aggregate move (the ≤1.9 % ceiling, §0); the win is honesty + unblocking,
   measured as: fewer confident-wrong nodes, unidentifiable nodes ending *weak* not *committed*.
