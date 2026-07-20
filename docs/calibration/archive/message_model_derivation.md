# The calibration message model — a derivation, and where the current relay goes wrong

**Status:** theory + design derivation, written 2026-07-16 (branch `calib-ambig-init-wip`), the session-close
record of the DOF-pie-relay work. **Audience:** a reviewer who does not know the codebase; everything needed is
stated here. **Companions:** [`dof_pie_relay_derivation.md`](dof_pie_relay_derivation.md) (the coherent `(λ,θ)`
relay — LANDED), [`dof_pie_relay_implementation_plan.md`](dof_pie_relay_implementation_plan.md) (S1–S4 + the A/B
and diagnostics), [`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) (count-zero-information — authoritative),
[`nascent_rna_sourcing_regression.md`](nascent_rna_sourcing_regression.md) (the nascent factory),
[`mature_crossing_gate.md`](mature_crossing_gate.md) (item 1 — the gate this document recommends dismantling).
**Evidence:** `scripts/debug/` + the session scratchpad (`message_model_experiments.py`,
`high_nascent_gate_test.py`, `intron_leak_diag.py`, `gate_ab_dof.py`).

> **The one-line result.** A message is **not** a pie. It is a product of independent **per-component density
> factors**, each a soft observation of one component with its own precision; the **recipient** combines them
> (with its own self-solve) into one coherent pie. The node's *state* is one degree of freedom per composition
> axis (one mode + one precision); the *messages* are free to speak about a single component. The DOF pie fix got
> the state right. What is still broken is the **precision** of the factors (`σ²_transfer = 0` ⇒ a dense node's
> over-confident message overrules a neighbour's correct self-solve) and one **band-aid** (the mature gate) that
> the coherent pie no longer needs and that inhibits it.

---

## 1. The problem in one paragraph

Calibration deconvolves each genomic node's fragment mass into a composition `(f_pos, f_neg, f_g)` and propagates
belief along a region↔boundary chain by belief propagation. The DOF pie fix (LANDED) made the per-node **state** a
coherent composition — coordinates `(λ = logit f_g[, θ])`, one precision `σ²_λ` per degree of freedom — so the
relayed pie is a valid composition by construction (`n_c ≤ M` everywhere; the 424–925 nodes/condition with a
component fraction > 1 → 0). That fixed the *representation*. The diagnostics below show the *residual* errors are
not in the representation but in the **message precision** and in one **heuristic** (item 1's mature gate). This
document derives what a message **is**, states the precision model, and recommends the next moves.

---

## 2. What a message is (the derivation)

### 2.1 Belief propagation on the chain

Node `i`'s posterior over its composition `x_i` is
```
    b_i(x_i)  ∝  ℓ_i(x_i)  ·  ∏_{j ∈ neighbours(i)}  m_{j→i}(x_i)
```
where `ℓ_i` is `i`'s **local likelihood** (its own strand data + reference + prior — its *self-solve*) and
`m_{j→i}` is the **message** from `j`. The message is the neighbour's evidence about `x_i`, obtained by
integrating the **edge potential** `ψ_{ij}(x_i, x_j)` against `j`'s belief:
```
    m_{j→i}(x_i)  =  ∫  ψ_{ij}(x_i, x_j)  b_j(x_j) / m_{i→j}(x_j)  dx_j .
```

### 2.2 The edge potential factorizes over components ⇒ the message is per-component

The edge potential encodes how adjacent nodes' compositions relate. Physically that relation is **per component,
in DENSITY** (`CALIBRATION_ARCHITECTURE.md` §5): gDNA density is smooth genomically; **nascent** RNA density is
continuous *within a transcript* (across exon↔intron seams); **mature** RNA is anchored and crosses a junction
only *as spliced*, never as an unspliced crossing. So the potential is a **product of per-component density
couplings**, each gated by that component's continuity across the edge:
```
    ψ_{ij}(x_i, x_j)  =  ∏_{c continuous across (i,j)}  ψ_c( ρ_c^i , ρ_c^j ; σ²_transfer,c )
```
with `ψ_c` a soft "`ρ_c^i ≈ ρ_c^j`" (a Gaussian on `log ρ_c^i − log ρ_c^j`, spread `σ²_transfer,c`). Substituting,
the message factorizes:
```
    m_{j→i}(x_i)  =  ∏_{c continuous}  f_c( ρ_c^i ; a_c , w_c )
```
— **a product of independent per-component density factors.** Each factor is a soft observation of **one
component's density** at the destination, with its **own** precision `w_c`. **The message carries no constraint
that the components form a pie.** It *can* say "I see gDNA at density X (precision P); I say nothing about RNA."

### 2.3 The two levels (this is the crux, and it was the source of the confusion)

| level | object | DOF / precision |
|---|---|---|
| **node STATE** (stored, relayed) | one coherent composition (a pie) | **one** mode + **one** precision per composition axis (`λ`; `θ` for AMBIG) — a single-strand node stores ONE precision, not separate gDNA/RNA precisions |
| **MESSAGE** (on an edge) | a **set of per-component density factors** | each factor its OWN mode + precision; they need not form a pie |

The recipient multiplies all incoming factors × its own local likelihood → the coherent pie (which the constrained
grid solve produces). **The pie lives at the recipient, never on the wire.** The DOF pie fix implemented exactly
this: `_scan` stores `(μ_λ, σ²_λ[, μ_θ, σ²_θ])` (state = one DOF) and emits per-component density messages
(`ρ_g`, `ρ_pos`, `ρ_neg`); the fold reconciles them onto the recipient's `λ`. *(A prior write-up mis-described the
recipient-frame reading of a density factor as "the message's f_g", and two same-component messages from two
neighbours as "gDNA vs RNA duelling" — both were description errors, not code: the wire carries `ρ_c = f_c·M/E`.)*

---

## 3. The precision model — count-zero-information, and the gap

A per-component density factor's precision is the honest inverse-variance of the source's density estimate,
transported:
```
    w_c  =  1 / ( Var(log f_c^src)   +   1/M_src   +   σ²_transfer,c )
             ╰── composition ──╯    ╰ magnitude ╯   ╰─ transport ─╯
                (STRAND Fisher info)  (total count)   (adjacent-node spread)
```

Reading it against the count-zero-information invariant (`CALIBRATION_ARCHITECTURE.md` §0 — *the count carries
ZERO composition information; it may enter ONLY as statistical power*):

* **`Var(log f_c^src)`** — the source's composition uncertainty. The count enters this ONLY through the strand
  Beta-Binomial Fisher information `∝ N(2κ−1)²` (a legitimate count→precision channel). At `κ=½` it is the
  reference width (count-independent) — correct: an unstranded node's composition is not sharpened by *any* count.
* **`1/M_src`** — the source's total-density (magnitude) sampling. Legitimate: it is the predictor's Poisson
  sampling noise (`CALIBRATION_ARCHITECTURE.md` §1.2), about the density's *scale*, not the composition.
* **`σ²_transfer,c`** — how much adjacent nodes' true densities genuinely differ (enrichment for gDNA;
  along-transcript continuity for RNA). **Count-independent.** **Currently ZERO.** This is the load-bearing gap.

**Why `σ²_transfer = 0` breaks things (measured, §7 E1/E3).** With the transport spread at zero, a **dense**
source's message precision is capped only by its own `M` (the toy exon: 1672 ≈ its `sm`), with **no per-hop
decay**, so it **overrules a sparse neighbour's correct self-solve**. The boundary in the toy self-solves to
`f_g = 0.963` (right — its own strand says gDNA) and the over-confident exon message drags it to 0.81. On the real
suite the message-free boundary self-solve is *better* than after messages on capOFF (e.g. 0.0387 vs 0.0519). The
DOF pie fix stopped the *unbounded* inflation (`n_src ≤ M`); it did **not** supply `σ²_transfer`. **This is the
next problem to solve, and it is the deepest one.** *(The count-term `1/M` is not a count-zero violation — it is
honest magnitude sampling; the problem is the missing transport term, not the count-term.)*

---

## 4. Where the current relay goes wrong (diagnosed, not guessed)

Two things, both at the precision / edge-potential layer, **neither** a representation bug:

1. **`σ²_transfer = 0` ⇒ over-confident transport overrides the self-solve** (§3, §7 E1/E3). The information to
   solve a boundary is usually already in its own self-solve; the message layer, un-damped, breaks it.
2. **The mature gate silences the RNA factor instead of defining it.** At an exon→intron edge the continuous
   components are gDNA **and nascent**; only unspliced *mature* is discontinuous. The correct nascent factor is
   `ρ_nascent^dst ≈ ρ_RNA^src − ρ_mature^src`, which for a pure-mature exon is **≈ 0** — a real "no nascent here"
   message. Because the 2-component pie cannot separate mature from nascent, item 1 blocked the **whole** RNA
   factor. Silence = "no information about RNA" ⇒ the recipient's coherent pie assigns the un-explained gDNA
   deficit to RNA (the toy junction's 18 % phantom RNA). A `nascent≈0` factor removes it (§7 E2).

---

## 5. The mature gate — DISMANTLED (landed 2026-07-16, commit `5e54fdc5`)

**Decision (owner, 2026-07-16): dismantle the mature gate temporarily and pursue an honest, gate-free solution;
re-introduce a gate only if we *cannot* solve accurately without it.** **LANDED** — the asymmetric
`send_s = mrna_active[dst] or not mrna_active[src]` gate is removed from `bp_solver._scan`; only the structural
per-strand `free_s` continuity gate remains. Revert `5e54fdc5` to restore the gate (the checkpoint `9b0f7419`
is the gate-intact state). The `mrna_active_*` mask stays computed in the statics for the nascent factory (§6.2).
The reasoning, with the evidence and its limits stated honestly:

* **The gate is a band-aid for problems the DOF pie fix now owns.** Item 1 was created because the incoherent
  three-fraction relay let an exon over-send its (mostly-mature) RNA into introns with no way to bound it. The
  coherent pie + an honest `σ²_transfer` bound that at the recipient's own self-solve — the mechanism the gate
  approximated. The **released Rigel has no mature gate and performs well**; the largest headroom is nascent-RNA
  resolution, which the DOF pie fix targets directly.
* **The gate inhibits the elegant solution.** It creates a gDNA/RNA **asymmetry** (it blocks the RNA factor but
  not the gDNA factor), and under the coherent pie a lone gDNA factor forces phantom RNA at the recipient (§4.2).
  It is **annotation-based** (keys on the exon *bit*, not the actual mature content), so it blocks *legitimate*
  nascent flow out of an actively-transcribed exon.
* **Honest limits of the evidence.** (a) On the current 7-condition suite the gate **helps 7/7** under the
  coherent pie (mean mwae 0.169 vs 0.193; `gate_ab_dof.py`) — because that suite is mature-heavy and the gate
  blocks a real mature leak. **So dismantling will regress the suite until `σ²_transfer` lands.** (b) The specific
  "very-high-nascent breaks the gate" hypothesis was **NOT confirmed** (`high_nascent_gate_test.py`: gate ON = gate
  OFF): an intron-boundary always has an intron flank, so the intron-side message (never gated) + the boundary's
  self-solve carry the nascent regardless of the exon-side gate. The gate is *neutral* there, not destructive.
* **The principle decides it.** We prefer a pristine architecture: solve the problem honestly (coherent pie +
  `σ²_transfer` + the nascent factor), and add a gate back only as a last resort. Two clean paths, and dismantling
  keeps the second open:
  * **keep the gate** ⇒ we are then *obligated* to manufacture a `nascent = RNA − mature ≈ 0` counter-factor to
    stop the pie manufacturing RNA (§4.2). A second band-aid to prop the first.
  * **dismantle the gate** ⇒ the exon→intron RNA factor flows again (mature absorbed at the junction, nascent
    imputed honestly); the mature-leak it prevented is instead bounded by `σ²_transfer` + the recipient's
    self-solve. **The pristine path.**

---

## 6. The two problems the next session must solve

### 6.1 The transfer-variance (`σ²_transfer`) model — the load-bearing one

> **DERIVED + PROVEN belief-free (2026-07-16): [`sigma2_transfer_derivation.md`](sigma2_transfer_derivation.md).**
> `σ²_transfer,g` is a **pair** quantity, **stratified on the belief-free total-density regime** of both
> endpoints (`[0, 1.6, 25]` on-capture: dep–dep / enr–enr / crossing; flat ≈0 off-capture), fit **belief-free
> on RNA-free anchor pairs**. Source-only is refuted (it blends the reliable same-mode edge with the crossing).
> The shipped total-density scalar is 12–28× too large off-capture (pure RNA noise). The NPMLE-on-total-density
> already shipped IS the regime classifier; the variances come from adjacent RNA-free pairs, not the NPMLE
> marginal width. Circularity cut: regime label + `σ²_transfer` are belief-free and fixed across passes.

Per-component, count-independent, honest, and **fit on raw observables, not relay-smoothed posteriors**
(`archive/count_space_relay_implementation_plan.md` §9): the gDNA leg from **gDNA-identifiable** edges
(RNA-free enriched exons — where the exon's gDNA density is directly observable) with an irreducible unobservable-
probe floor; the RNA leg from deep constitutive exon segments; **stratified** on observable probe-geometry proxies
(spliced mass at the boundary, density regime, boundary type). It must give the relay **per-hop decay** so a dense
node's reach is finite, and it must let a correct self-solve survive an incorrect transported message. Target: the
message-free self-solve is a floor the messages may only *improve* on, never *break* (§7 E1). It is the same
`σ²_bio` these earlier notes derive (`archive/imputation_variance_model.md`, `archive/honest_precision_unified_design.md`).

### 6.2 The nascent factor (`RNA − mature`) — replacing the gate's job, honestly

Nascent density = total RNA density − mature density. Mature is measurable (the motif-stranded spliced count).
Nascent is continuous, so the flanking **intron's RNA density** is the nascent baseline. So an exon→boundary
factor imputes the boundary's nascent as ≈ the flanking intron's RNA — which for a pure-mature locus is ≈ 0
(the counter-message the gate currently fakes by silencing). This is the deferred nascent-factory work
(`nascent_rna_sourcing_regression.md`); with the gate dismantled it becomes the *principled* RNA-side message.

---

## 7. Experiments (the evidence)

Toy boundary **B2** (exon→intron `+` junction; truth `f_g = 1.0`; `message_model_experiments.py`):

| experiment | result | reads |
|---|---|---|
| **E1** decompose B2 | self-solve **0.963** → +gDNA msg **0.81** | the message BREAKS a correct self-solve |
| **E2** nascent-counter | +`f_pos≈0` factor at prec 50 → **0.984** | a `nascent≈0` factor is the missing counterweight |
| **E3** `σ²_transfer` | weaken exon prec 2453→~20 → **0.95** | a transport floor lets the self-solve reassert |
| **E4** forward transient | exon fwd 0.235→B2 0.74; final 0.30→0.94; truth 0.32→0.99 | the source's forward-only belief understates |

Suite (`gate_ab_dof.py`, DOF-fix, 7 conds): the gate **helps 7/7** (mean 0.169 vs 0.193) — mature-heavy suite.
High-nascent toy (`high_nascent_gate_test.py`): gate ON = gate OFF (both 0.097) — the gate is neutral, not
destructive, in high nascent. Accuracy A/B (DOF-fix vs shipped, `ablate_replay`): flagship unstranded+capture leak
**−21…−25 %**, capOFF neutral/better, **stranded regresses** — the pie helps most in the count-zero-info regime.

---

## 8. The path forward

1. **Dismantle the mature gate** (§5) — ✅ **DONE** (commit `5e54fdc5`). Gate-free relay; the temporary suite
   regression is documented and measured (gate helped 7/7: ALL mwae 0.1690 → 0.1925, concentrated in introns —
   the mature leak; verified equivalent to the all-True mask). Recover it with (2), not with the gate.
2. **Build the `σ²_transfer` model** (§6.1) — the load-bearing fix; per-component, count-independent, fit on raw
   observables, giving honest per-hop decay. This is where the "precision is fundamentally broken" is repaired.
   **NEXT.**
3. **The nascent factor** (§6.2) — replace the gate's job honestly (`RNA − mature`, intron-baseline).
4. Re-measure the flagship + the whole suite; only if the pristine solution is *unacceptable* does a gate return
   as a documented last resort.

**Established vs assumed:** the message-is-per-component derivation (§2) — **derived** (BP factorization).
`σ²_transfer=0` overrides self-solves — **measured** (E1/E3). The gate silences instead of defining the RNA
factor — **measured** (E2). The gate helps the current suite — **measured**. The gate destroys high nascent —
**REFUTED for the tested construction** (neutral); the case against it is the pristine-architecture principle +
the coherent-pie asymmetry, not a demonstrated high-nascent catastrophe.
