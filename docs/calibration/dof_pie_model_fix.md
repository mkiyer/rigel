# Roadmap item 2 — degrees of freedom: relay the pie, not three fractions

**Status:** design + stepwise plan. Written 2026-07-16, branch `calib-ambig-init-wip`.

> **➡ THE THEORY IS DERIVED + REVIEWED; THE PLAN IS READY (2026-07-16).** Three companions supersede this
> document's §2–§4 as the live artefacts:
> * [`dof_pie_relay_derivation.md`](dof_pie_relay_derivation.md) — the rigorous, reviewer-facing derivation
>   (third-party reviewed: *proceed*): the coherent `(λ,θ)` EP relay, the two-level model (stored coordinates /
>   transmitted density), the count-zero-info count-term `1/M_src`, the preserved priority-#3 measurement, the
>   honest `σ²_transfer=0` liabilities.
> * [`dof_pie_relay_implementation_plan.md`](dof_pie_relay_implementation_plan.md) — the implementation-ready plan
>   (review's three landmines defused; the **two-stage self-correcting fold**; minimal footprint; self-critiqued;
>   sign-off gate).
> * [`scripts/debug/dof_pie_relay_check.py`](../../scripts/debug/dof_pie_relay_check.py) — the numpy prototype
>   (C1–C9 all pass; C8 disproved the first single-grid fold and drove the two-stage design).
>
> **This document's §1 measurement is re-derived POST-GATE below and by
> [`scripts/debug/pie_probe.py`](../../scripts/debug/pie_probe.py); it stands** (62–70 % of solvable nodes
> incoherent; 424–925 nodes/condition with a component fraction > 1; MAX 600× under capture). **Per the agreed
> sequencing, NO production code (S1–S6) lands until the plan's §10 sign-off gate is cleared.**
> **⚠ The framing below is stale; §1's measurements are not.** Item 1's §3.2 has been **retired** as an exact
> algebraic no-op and replaced by [`mature_crossing_gate.md`](mature_crossing_gate.md). This document's stated
> justification — *"it bounds item 1's residual"* — was **conditional on §3.2 shipping** (under §3.2 the residual
> is 0.1517 and the argument is sound). Under the gate the capOFF residual is largely gone and the capON
> remainder is the **gDNA** channel, not composition. **This document must re-derive its case on the post-gate
> residual** (`mature_crossing_gate.md` §6.3) — it very likely still earns it, since §1's defect is real,
> measured, and independent of item 1. Sequencing: **after** the gate, unconfounded.
>
> Two instruments here cost no behaviour and can land any time: `assert n_src <= sm` on every edge (B740
> violates it 52×) and the pie probe `|fbg+fbp+fbn − 1| < 1e-9` as `xfail(strict=True)`.

**Companion:** [`mature_crossing_gate.md`](mature_crossing_gate.md) — the live item 1;
[`splice_junction_absorption_fix.md`](archive/splice_junction_absorption_fix.md) — its superseded predecessor.
**Predicted long ago:** `calibration_roadmap.md` §3 Phase 5, and `mynotes.md` §2 ("Respect message
interdependence"). Both were right. This is the measurement and the prescription.

---

## 1. The defect, measured

A composition **cannot** sum to anything but 1. The relay's does:

```
fbg + fbp + fbn over 2538 solvable nodes:
    p1=0.822   p10=0.991   p50=1.086   p90=4.348   p99=45.773   MAX=93.351
    pie > 1.01 : 1819 nodes (71.7%)
    pie < 0.99 :  251 nodes ( 9.9%)
```

**71.7% of solvable nodes carry a relay belief that is not a composition.** The traced case:

```
B740 relay : fbg=0.7883  fbp=51.9132  fbn=0.0000   SUM = 52.70
B740 solve : f_g=0.0160  f_pos=0.9840 f_neg=0.0000  SUM =  1.0000
```

The node simultaneously "believes" it is **5191% RNA and 79% gDNA** — and *that* is the state it relays.

### 1.1 The split: the solve is constrained, the relay is not

This is the whole finding, and it is exactly what `calibration_roadmap.md` §3 Phase 5 predicted.

**The solve is correct.** `simplex_logodds._local_loglik_logodds`:

```python
fg2   = fg[None, :]          # f_g = sigma(lambda) on ONE grid
f_act = 1.0 - fg2            # the live strand carries 1 - f_g  <-- the constraint, BY CONSTRUCTION
f_pos = np.where(pos_live, f_act, 0.0)
f_neg = np.where(neg_live, f_act, 0.0)
```

One free parameter (`λ`); both messages are functions of it; `f_g + f_act ≡ 1` identically. **The 1-DOF
concept IS implemented** — in the final solve.

**The relay is not.** `bp_solver._scan`, three independent updates with nothing linking them:

```python
fbg[i] = math.exp((pg_loc[i]*lfg_loc[i] + pr*mo)/pt);  vbg[i] = 1.0/pt   # line 418
fbp[i] = math.exp((pp_loc[i]*lfp_loc[i] + pr*mo)/pt);  vbp[i] = 1.0/pt   # line 459
fbn[i] = math.exp((pn_loc[i]*lfn_loc[i] + pr*mo)/pt);  vbn[i] = 1.0/pt   # line 477
```

Three independent geometric means. **The quantity propagated between nodes is not a composition at all.**

`node_geometry._type_belief` seeds the same defect — it sets `var_g` / `var_p` / `var_n` **independently**
from ONE 1-DOF solve, so a node can be born with "`var_g=0` but `var_p=∞`", which on a 1-DOF node is
self-contradictory.

### 1.2 Why it compounds — the precision inflation

`pr = n_src/(n_src*vb_src + 1)` with **`n_src = fbp[src]*sm`** — the message precision is computed **from the
unconstrained relay belief**. So at B740's outgoing edge `fbp = 51.9` makes `n_src` **52× the node's actual
mass**.

**The broken pie does not merely corrupt the message content — it inflates the message precision by the same
factor.** The louder the wrong answer, the more confidently it is asserted, and the next hop inherits both.

Measured at B740:

```
B740's own unspliced count n         = 67
its strand Fisher information         = 16.1        <- all the node itself can muster
incoming RNA  message precision       = 27,208.6    <- 1,691x the node's own evidence
incoming gDNA message precision       =    447.9
RNA outguns gDNA by                   = 61x
```

Both channels describe **the same 67 fragments**. The RNA channel wins by 61× purely because the source exon
is RNA-dense. This is the user's *"higher densities should not necessarily mean higher precision — then the
highest-density piece will steal fragments from the other parts"*, quantified.

### 1.3 The consequence item 1 cannot fix

**9.4% of genuine junction faces demand `f_RNA > 1` even after a correct mature subtraction (max 41.5).** A
fraction cannot be 41.5. That the arithmetic can *ask* for it is this defect, and no absorption rule, no
allocation rule, and no `σ²_transfer` value bounds it — because nothing in the representation says a
composition is a composition.

---

## 2. The prescription

> **Relay the free coordinates, not the fractions. One number + one precision per degree of freedom.
> Derive the fractions.**

The coordinates are already chosen and already validated — `reference_prior_derivation.md` §10:

| node class | DOF | free coordinates | fractions |
|---|---|---|---|
| single-strand | **1** | `λ = logit(f_g)` | `f_g = σ(λ)`, `f_live = 1 − f_g`, dead strand `= 0` |
| AMBIG | **2** | `λ`, `θ = arcsin(τ)` | `f_g = σ(λ)`; `f_± = (1−f_g)(1 ± sin θ)/2` |
| G1 (intergenic/TSS/TES) | **0** | — | `{0,0,1}`, locked (already correct; `solvable` gates them out) |

**This is the same parameterisation the solve already uses.** The solve integrates ψ on exactly this
`(λ[,θ])` grid. Item 2 is not a new model — it is **making the relay speak the solve's language.**

`fbp = 51.9` becomes **unrepresentable by construction**: `f_pos = (1−σ(λ))(1+sin θ)/2 ∈ [0,1]` for every
`(λ, θ) ∈ ℝ × [−π/2, π/2]`. The pie sums to 1 identically, and `n_src = f_pos·sm ≤ sm` — the precision
inflation dies with it.

### 2.1 The precision architecture (the user's open question)

*"a two-stranded node has 2-DOF. The node's belief state should have the correct precision architecture.. I am
not sure how to model it."*

**One Gaussian per free coordinate.** 1-DOF: `N(λ; μ_λ, σ²_λ)`. 2-DOF: `N(λ; μ_λ, σ²_λ) × N(θ; μ_θ, σ²_θ)`.

Two facts make the product form legitimate rather than convenient:

* **`λ` and `θ` are information-ORTHOGONAL.** `I_{f_g,τ} = 0` **exactly** — verified numerically to ~1e-11
  (`reference_prior_derivation.md` §10.1, `scripts/debug/reference_prior_bb_check.py`). The Fisher matrix is
  diagonal in these coordinates, so a product of independents is the correct second-order belief, not an
  approximation.
* **The reference already factorises the same way** — Beta(½,½)_λ ⊗ Beta(½,½)_θ, and in `(λ, θ)` its θ term
  **vanishes identically** (§10.3). The relay and the reference then live in the same frame.

**A component's precision is no longer a thing that exists.** "The RNA precision" and "the gDNA precision"
were never two facts — they are two views of `σ²_λ`. Asking whether the solver "respects the DNA precision
when it steals RNA" dissolves: there is one number, and moving it moves both fractions. That IS the pie.

### 2.2 What a message becomes

Today: three `(mode, prec)` pairs on `log f_g`, `log f_pos`, `log f_neg`.
After: the message is still a **density** claim (the currency is settled and unchanged) — but the destination
folds it onto **`λ`**, where its own belief lives, instead of onto an unconstrained fraction.

The gDNA and RNA messages stop being competitors and become **two observations of one coordinate**. A gDNA
message saying "you are dense in gDNA" and an RNA message saying "you are dense in RNA" are then automatically
in tension on a single axis, precision-weighted — which is the interdependence the pie is supposed to encode.

---

## 3. What this does NOT fix (state it now)

* **It does not make a wrong message right.** B740's exon still sends mature it should not — that is item 1.
  Item 2 only guarantees the wrong answer stays a *composition* and cannot inflate its own precision 52×.
* **It does not choose `σ²_transfer`.** Still open, still to be decided from data (see the paradigm work).
* **It does not touch the reference or the prior.** Item 1 (#1 reference) and the NPMLE are unaffected.
* **It will move goldens.** Expected and correct.

---

## 4. Stepwise implementation plan

Each step small, independently verifiable, and ordered so nothing lands unmeasured.

**S1 — Instrument first, change nothing.** Land the pie-coherence probe as a permanent diagnostic:
`fbg+fbp+fbn` per node at the end of each scan, reported by node class. *Verify: reproduces §1's numbers
(p50=1.086, 71.7% > 1.01, max 93.35). This is the metric every later step is judged on.*

**S2 — A test that fails on today's code.** Assert the relay's pie sums to 1 (within float tolerance) on a
synthetic chain that includes an RNA-dense exon beside a boundary. **It must FAIL now** — verify by running it
before the fix (the `region_eff_length` lesson: a test that does not fail on the broken code is a
restatement, not a test). Mark `xfail(strict=True)` until S4 flips it.

**S3 — `_type_belief` seeds one coordinate, not three.** Derive `(μ_λ, σ²_λ)` (+ `(μ_θ, σ²_θ)` for AMBIG)
from the same 1-DOF solve it already runs, and construct `var_g`/`var_p`/`var_n` **from** them rather than
setting them independently. *Verify: the born-belief pie sums to 1 on every node; "var_g=0 with var_p=∞"
becomes unconstructible.*

**S4 — `_scan` relays `(λ, [θ])`.** Replace the three independent geometric means with one update per free
coordinate; recover `fbg`/`fbp`/`fbn` by the §2 formulas wherever the loop needs a fraction. *Verify: S2 flips
to pass; the pie is 1.000 everywhere; `n_src = f_pos·sm ≤ sm` on every edge (assert it).*

**S5 — Re-measure the message precisions.** `intron_message_trace` at B740: `prec_p` must fall from 27,208
toward O(n_src) with `n_src ≤ 67`, and the RNA/gDNA imbalance from 61× toward ~1 (they now share `σ²_λ`).
*This is the step that proves the compounding is gone.*

**S6 — `selfsolve_diag --stage both`, capOFF and capON.** Report against the item-1 predictions. **Sequencing
note:** if item 1 lands first, its residual tail (§4 there) should shrink here; if item 2 lands first, expect
the bulk intron damage to REMAIN (item 2 bounds it, it does not remove it — the message is still wrong). Land
them in either order but **measure each alone before both**, or the attribution is lost — as it was for the
`region_eff_length` fix, where we credited a real bug with everything nearby that looked wrong.

---

## 5. Falsification tests (theory space)

* **The pie is 1.** `|fbg+fbp+fbn − 1| < 1e-9` on every solvable node, every condition. Binary.
* **Fractions are fractions.** `f_g, f_pos, f_neg ∈ [0,1]` by construction — assert on random `(λ, θ)`,
  including `θ = ±π/2` and `|λ| = L`.
* **Precision is bounded by mass.** `n_src ≤ sm` on every edge. Today B740 violates it 52×.
* **1-DOF reduces.** An AMBIG node with one strand structurally dead must give bit-identical results to the
  single-strand path. Currently untestable (the classes carry different objects); after S4 it is a one-line
  assert.
* **Orthogonality holds in-code.** Re-verify `I_{f_g,θ} = 0` on the implemented grids
  (`reference_prior_bb_check.py` already does this analytically; make it a test).

---

## 6. Established vs assumed

| | |
|---|---|
| relay pie sums to 1.086 median, 93.35 max; 71.7% > 1.01 | **measured** |
| the solve IS constrained; the relay is NOT | **read from the code + measured** (B740: relay 52.70, solve 1.0000) |
| `n_src = fbp·sm` inflates precision by the pie violation | **measured** (B740: `fbp=51.9`, `prec_p=27,209` vs own evidence 16.1) |
| RNA outguns gDNA 61× on the same 67 fragments | **measured** |
| 9.4% of junction faces demand `f_RNA > 1` after a correct subtraction | **measured** |
| `λ` ⊥ `θ` ⇒ a product of independent Gaussians is the right belief | **verified** (`I_{f_g,τ} = 0` to 1e-11) |
| relaying `(λ,θ)` removes the precision inflation | **follows** from `f_pos ≤ 1` ⇒ `n_src ≤ sm`; verify at S5 |
| item 2 bounds item 1's residual tail | **assumed** — the 9.4% figure supports it; S6 tests it |

---

## 7. Why this is the model fix and item 1 is not

Item 1 is arithmetic: an exon sends a molecule that does not cross. Fix the taxonomy, the message is right.

Item 2 is representational: **the system has no way to say "this is a composition."** So a wrong message is
unbounded, its precision self-amplifies, and the two channels compete for the same fragments as if they were
independent. Item 1 makes the messages true; **item 2 makes them possible.**

The roadmap called this in Phase 5 — *"carry ONE number (+ its precision) on the free axis for 1-DOF, TWO for
2-DOF; derive the rest"* — and `mynotes` §2 called it before that. Neither was acted on because nothing
measured the pie. It is measured now: **71.7%, max 93.35.**
