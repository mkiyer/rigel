<!-- title: Node state vs message content — the self-confidence ⊥ predictive-authority separation -->
# Node state vs message content — separating self-confidence from predictive authority

**Status:** design note (2026-06-27). Pre-refactor. Prompted by the muddle exposed when the gDNA-emission
gate was decoupled from the RNA strand-gate (the TSS/TES bug fix): once an intergenic / gene-boundary node is
*allowed* to emit, we must be explicit about **what** it is allowed to say and **how loudly** — and that turns
out to be a different axis from how sure it is about itself. This note settles the conceptual model and scopes
the cleanup. Builds on [`precision_state_design.md`](precision_state_design.md) (own-uncertainty is additive in
the send precision), [`honest_precision_message_design.md`](honest_precision_message_design.md) (a message is
`(density, precision)`, count-honest), [`node_state_representation.md`](node_state_representation.md) (the
fraction/density/count/mass currencies), and the count-zero-info principle in
[`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md). It also subsumes the "robust message weight" that
[`forward_backward_plan.md`](forward_backward_plan.md) / [`forward_backward_state.md`](forward_backward_state.md)
filed as a *separate* cliff patch — here it is not a patch, it is the definition of predictive authority.

## 0. TL;DR

A node has a **state** and it **emits messages**. These are different objects with different precisions:

- **State** `(f_c, τ_own_c)` per component `c ∈ {gDNA, RNA₊, RNA₋}`: the fraction + **how sure the node is
  about its own composition**. `τ_own=∞` locked (structural certainty); `0` no idea; finite = solved posterior.
- **Message** `(mode_c, τ_pred_c)` per neighbour per component: the source's *prediction* of the neighbour's
  `f_c` + **how much authority the source has to set the neighbour's value**.

**The principle:** `τ_own` and `τ_pred` are independent. A node may be *supremely* confident about itself
(`τ_own=∞`) and have *zero* authority over a neighbour (`τ_pred=0`). The intergenic node is the canonical case:
it knows for certain it is all-gDNA, and it knows *nothing* about whether the adjacent exon is expressed.

The system realizes this with **three separately legible knobs**, currently tangled into the one `solvable`
flag and a regime-blind `σ²_bio`:

```
τ_pred_c(src→dst) =   GATE_c(edge)                 # structural transmissibility ∈ {0,1}   ← from the SIGNATURE
                    ·  W_regime(ρ_src, θ̂_dst)        # predictive authority across regimes   ← the robust weight
                    ·  1 / ( vo + σ²_bio + pois )    # predictive uncertainty (own-bounded)   ← var~mean + own var
```

| principle | knob | source | status |
|---|---|---|---|
| structural-zero vs observed-zero | **`GATE_c`** | region signature (`free_pos`/`free_neg`; `mass_unspliced>0`) | gDNA done; RNA implicit → make explicit |
| self-confidence can only *lower* authority | **`vo` (own var, additive)** | `Var(f_c^src)` | ✓ already additive |
| self-confidence ≠ authority *across regimes* | **`W_regime`** | robust weight `(ν+1)/(ν+r²)` | **MISSING** — the cliff drag |

## 1. The two objects

### 1.1 State — what the node believes about itself

Per component, `(f_c, τ_own_c)`. `τ_own` is the node's *own* posterior precision; it is set by the node's
*own* evidence (its strand likelihood, its signature lock, the global prior) and is moved only by messages the
node chooses to *receive*. Semantically:

- `τ_own = ∞` (locked): structural certainty. An **intergenic** node's RNA (`f_RNA=0`, can never be RNA); a
  single-strand node's forbidden strand. A locked component **ignores all incoming messages** for that
  component — no message can move it.
- `τ_own = 0`: no own information (an unsolved AMBIG node's gDNA before any message).
- finite: a solved posterior; movable by sufficiently authoritative messages.

### 1.2 Message — what the node tells a neighbour

Per neighbour per component, `(mode_c, τ_pred_c)`. `mode_c` is the source's belief *projected onto the
destination's fraction scale* via the density currency (`ρ_src = f_c·M_src/E_src`; `mode = ρ_src·E_dst/M_dst`,
spliced-adjusted). `τ_pred_c` is the **authority** — how hard this prediction pulls the destination.

The whole point: **`τ_pred` is not `τ_own`.** A node computes one message per neighbour, and its authority is
governed by its *ability to predict that neighbour*, not by how sure it is about itself.

## 2. The discrimination you asked for: structural-zero vs observed-zero

Two nodes can both have `f_RNA = 0`. They are *not* the same, and the difference is **structural** (it comes
from the region signature, never from the observed counts — the count-zero-info principle):

| | RNA state | RNA emission to a neighbour |
|---|---|---|
| **intergenic** (structural zero) | `f_RNA=0, τ_own=∞` (locked; ignores all incoming RNA) | `GATE_RNA = 0` on every edge → `τ_pred(RNA)=0`. RNA cannot cross a gene boundary, so it carries **zero predictive information** about a neighbour's RNA. It must *not* emit "there's no RNA here." |
| **unexpressed exon** (observed zero) | `f_RNA≈0, τ_own=` modest (it *measured* ≈0; still movable) | `GATE_RNA = 1` across a strand-continuous edge → `mode≈0, τ_pred>0`. It *can* tell its sibling exon, across the intron, "this gene looks unexpressed." |

The discriminator is `GATE_c`, read from the signature: intergenic has no strand bits (`free_pos=free_neg=
False`); an exon has its strand bit *regardless of expression*. So "I'm zero because I'm intergenic" gates RNA
emission **off**; "I'm zero because I'm not expressed" leaves it **on** (with a near-zero mode). This is already
latent in `free_pos`/`free_neg`; the cleanup makes it an explicit per-component gate rather than a side effect
of `solvable`.

The same split applies to **gDNA**, mirror-imaged: an intergenic / TSS / TES seam is *structurally* all-gDNA
(`τ_own(gDNA)=∞`, locked) and **may** emit gDNA (`GATE_gDNA = 1` wherever `mass_unspliced>0` — genomic
continuity). That is the bug fix already in the tree (`gdna_emit`). What it must **not** do is emit with high
authority across a regime boundary — see §4.

## 3. Why self-confidence cannot inflate authority (already correct)

The send precision is, in density space, `1/(vo + σ²_bio + pois)` with `vo = (M_src/E_src)²·(1/τ_own_src)` the
source's own density variance ([`precision_state_design.md`](precision_state_design.md) Q2). `vo` is **additive
in the denominator**, so:

- an *unsure* source (`τ_own→0` ⇒ `vo→∞`) ⇒ `τ_pred→0`: it speaks quietly. ✓
- a *locked* source (`τ_own=∞` ⇒ `vo→0`) ⇒ `τ_pred → 1/(σ²_bio+pois)`: capped at the predictive ceiling, **not**
  `∞`. ✓

So self-confidence can only *lower* authority (an unsure source predicts worse), never raise it past `σ²_bio`.
This half of the principle is structurally present. The gap is that `σ²_bio` is the *wrong ceiling at a regime
boundary* — §4.

## 4. The missing knob: predictive authority across regimes (`W_regime`)

`σ²_bio` is the var~mean between-node spread learned on the bulk; it does **not** know that a gene boundary is a
*molecular regime change*. So a depleted intergenic node predicting an *enriched* exon's gDNA gets `τ_pred ≈
1/(σ²_bio+pois)` — high — and **bosses the exon down** to the depleted level (the capture cliff). That is the
node being "supremely confident about itself" *and wrongly granted authority over a neighbour in a different
regime* — exactly the failure the principle forbids.

The fix is a **regime-awareness multiplier** on `τ_pred`, the Student-t robust weight (prototyped in
`scripts/debug/fb_cliff_toy.py`, proven cliff-correct + uniform-harmless):

```
W_regime = (ν+1) / (ν + r²),     r² = (ρ_src − θ̂_dst)² / v_msg
```

where `θ̂_dst` is the destination's *own local* density belief and `v_msg` the message variance. When the
source's prediction agrees with the destination's own read (`r²` small) → `W≈1` (full authority). When it
disagrees — a depleted source predicting an enriched destination, i.e. a regime boundary — `r²` is large →
`W→0` → **the self-confident-but-poor-predictor sends a near-zero-authority message.** This is the literal
mechanism for "confident about myself, zero authority over you," and it is therefore *part of the definition of
`τ_pred`*, not a separate cliff patch. `ν` is the one new constant (Student-t dof; the heavy-tail knob — discuss
per the no-magic-numbers rule, default ν≈4 from the toy).

## 5. Where the current code is muddled

1. **`solvable = (fp|fn) & mass>0` is overloaded.** One boolean does three jobs: (a) "can I solve my own
   state," (b) "can I emit a message," (c) historically *suppressed* intergenic gDNA emission on the wrong
   theory that a locked node "carries no evidence." (a) and (b) are different questions, and they are
   per-component (gDNA vs each RNA strand), not one flag. The gDNA half is split off (`gdna_emit`); finish the
   job: explicit `state_solvable_c` and `emit_c` per component.

2. **The message is built inline in `_scan`** with `vo`, `σ²_bio`, `pois`, the `(M_dst/E_dst)²` Jacobian, and
   the `[0,1]` mode clip all interleaved per component × per direction. The three knobs of §0 are not separable
   by eye, which is *why* it was unclear whether intergenic emits RNA (it does not — the `fp` gate — but you
   have to trace it). Route every message through **one** function so the separation is structural.

3. **No `W_regime`.** Predictive authority is currently `GATE · 1/(vo+σ²_bio+pois)` — missing the regime factor,
   so cross-regime self-confidence leaks into authority (the cliff).

## 6. Target shape

```python
# state: per component, per node
f_c, tau_own_c            # tau_own = inf locked / 0 unknown / finite solved

# emission gate: per component, per edge — STRUCTURAL, from the signature
GATE_gdna(edge) = mass_unspliced[src] > 0                  # genomic continuity (strand-agnostic)
GATE_rna_s(edge) = free_s[src] & free_s[dst]               # strand-s continuity (both endpoints)

# ONE message function
def message(state_src, geom, edge, c):
    if not GATE_c(edge):
        return (None, 0.0)                                 # zero authority — structurally silent
    rho_src = f_c_src * M_src / E_src
    mode    = clip((rho_src * E_dst - spliced_dst) / M_dst, 0, 1)
    vo      = (M_src / E_src)**2 * (1 / tau_own_c_src)      # own var — ADDITIVE only
    tau_pred = W_regime(rho_src, theta_hat_dst, v_msg) * (M_dst/E_dst)**2 / (vo + sigma2_bio + pois)
    return (mode, tau_pred)
```

- `GATE_c` = the structural-zero vs observed-zero discriminator (§2).
- `vo` additive = self-confidence can only lower authority (§3).
- `W_regime` = self-confidence does not transfer across regimes (§4).
- A locked component (`tau_own=∞`) **receives** nothing (state ignores messages); it still **emits** per its
  gate + authority. Self-confidence (receive side) and authority (send side) are now physically separate code.

## 7. Refactor plan (after sign-off)

1. **Per-component gates.** Replace `solvable` with `state_solvable_c` (the write/solve gate) and `emit_c` (the
   message gate), per component. gDNA: `mass_unspliced>0`. RNA_s: `free_s`. (Mostly a rename + split; the gДНК
   `gdna_emit` already lands here.) Confirm intergenic/TSS/TES: `emit_gdna=1`, `emit_rna=0`, `tau_own(RNA)=∞`.
2. **Extract `message(...)`** from the inline `_scan` arithmetic — one function, the §6 shape, called for
   forward α and backward β identically. Behaviour-preserving except as below; lock with an FB-vs-current
   agreement test on a small chain.
3. **Add `W_regime`** (the robust weight) as the *only* behavioural change in this pass. Validate on the
   benchmark (`evaluate_suite.py` + `precision_benchmark_report.py`): the cliff scenarios
   (unstranded+gDNA, capture) must improve; flagship + zero-gDNA must not regress.
4. **Regression tests** for the principle itself (not just the numbers):
   - an intergenic node emits **zero** RNA authority on every edge (`τ_pred(RNA)=0`), regardless of counts;
   - a TSS/TES seam with gDNA crossing **does** emit a gDNA message (the bug-fix lock);
   - a locked component is **not moved** by any incoming message (state ⊥ receive);
   - a depleted source predicting an enriched destination sends **low** authority (`W_regime` cliff lock).

## 8. Open questions

- **`ν` (Student-t dof).** One new constant. Default ν≈4 (the toy); revisit under the no-magic-numbers rule —
  can it be tied to the var~mean residual tail rather than fixed?
- **`θ̂_dst` for `W_regime`.** Use the destination's *local* (message-free) belief as the disagreement anchor;
  it is available in FB step (A) before the scans. Confirm it is the right anchor (vs the running forward
  belief, which would couple the weight to scan order).
- **Does `W_regime` subsume the `(M_dst/E_dst)²` Jacobian dishonesty** flagged in
  [`honest_precision_message_design.md`](honest_precision_message_design.md) §1, or are they orthogonal? The
  Jacobian over-credits high-mass destinations; `W_regime` down-credits regime-discordant ones. Likely
  orthogonal — keep the question open and measure.
