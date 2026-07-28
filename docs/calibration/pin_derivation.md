# THE PIN — what it is, what it must be, and where the circularity comes from

**Owner-directed derivation, 2026-07-27.** Written because `bp_solver._pin_v` was measured to feed a node
its own guess back as an incoming message, and the question is whether that is an oversight with a quick fix
or a symptom of something structural. **Answer: one conceptual conflation, a two-part local fix, but with
suite-wide reach and a second term already built on top of the contamination.**

Evidence: `scratchpad/gdna_d{5,6,7,8,9}_*.py`; the measured trace is in
`gdna_hyperprior_production_plan.md` §"THE BUG, FOUND".

---

## 1. What belief propagation requires of a message

Sum–product:

```
    m_{j→i}(x_i)  ∝  ∫ ψ_j(x_j) · ψ_ij(x_i, x_j) · Π_{k ∈ N(j)\{i}} m_{k→j}(x_j) dx_j
    b_i(x_i)      ∝  ψ_i(x_i) · Π_{j ∈ N(i)} m_{j→i}(x_i)
```

The message is a **function of** `x_i` but is **built from** node `j`'s evidence and `j`'s *other*
neighbours. If `m_{j→i}` contains any factor of `ψ_i`, then `ψ_i` appears **twice** in `b_i` — the node
confirms itself, and its confidence is the square of what it earned.

### The sharp line: node `i`'s CONSTANTS are allowed, node `i`'s BELIEFS are not

This distinction is what the whole question turns on.

* ✅ **Allowed — node `i`'s structural constants and the change of variables.** The latent at node `i` is
  its composition `f^i` on the 2-simplex; densities are `ρ_c^i = f_c^i · M_i / E_c^i`. The effective lengths
  `E_c^i` are geometry, and `M_i` is the node's observed unspliced mass. Writing a claim in node `i`'s
  coordinates *requires* both. That is a **reparameterization of `x_i`**, which lives in the pairwise
  potential `ψ_ij` — and `ψ_ij` is allowed to involve both endpoints, by definition. It is also exactly the
  **count-zero-information** principle in `CALIBRATION_ARCHITECTURE.md`: `M_i` carries no gDNA/RNA
  information, so using it to express *what composition means at node i* smuggles nothing in.
* ⛔ **Forbidden — node `i`'s inferred quantities.** Its message-free self-solve `f^own`, its own densities
  `(og, op, on)`, its posterior. These are `ψ_i`'s **output**. Any of them inside `m_{j→i}` is
  self-confirmation.

**`_pin_v` uses `M_i` and `E_c^i` (fine) and `og, op, on` (not fine).** Everything below follows from that
one sentence.

---

## 2. What the pin actually is

`bp_solver._pin_v`, verbatim in substance:

```python
sg = np.where(pg_ > 0.0, g, og)          # og/op/on = the DESTINATION's own self-solve densities
sp = np.where(pp_ > 0.0, p, op)
sn = np.where(pn_ > 0.0, n, on)
S  = sg*E_g + (sp + sn)*E_r              # a "mass budget" for the destination
k  = M / S                               # M = the DESTINATION's observed mass
return g*k, p*k, n*k
```

It rescales an incoming density claim by a common factor so the claim's implied fragment count equals the
destination's observed mass.

---

## 3. What the pin is FOR — two jobs, both legitimate

**Job A — transfer COMPOSITION, not scale.** The imputation premise is *"my neighbour and I share a
composition"* — never *"we share a density"* (hybrid capture makes densities differ by 1000×). A neighbour's
density vector must therefore be mapped onto the destination's own scale, keeping only its direction on the
simplex. A common rescale to `M_i` does exactly that.

**Job B — cancel the reframe.** The enrichment ratio `r` is the least trustworthy quantity in the relay. A
common rescale to `M_i` removes it **exactly** (verified to 1.8e-15), which is why substituting the ORACLE
capture step for `r` buys ≈ 0. That robustness is real and must be preserved by any replacement.

Both jobs are right. The pin is not a hack; it is the operator those two jobs call for.

---

## 4. THE DERIVATION — what the operator must be

Let node `j` be the source. Its accumulator satisfies, **identically**,

```
    Σ_c ρ_c^j · E_c^j  =  M_j                                                     (conservation at the source)
```

Under premise P (`f^i = f^j`), the message is a likelihood over `f^i`. The source's own share of any
component it measured is therefore computable **from source data alone**:

```
    f̂_c^j  =  ρ̂_c^j · E_c^j / M_j                                                        (★)
```

and the delivered claim at node `i` is

```
    ρ_c^i  =  f̂_c^j · M_i / E_c^i                                                        (★★)
```

**(★)–(★★) has every property the pin was reaching for, and none of the defects:**

1. **No destination belief.** `f̂^j` uses only `j`'s densities, `j`'s effective lengths and `j`'s observed
   mass. `M_i, E_c^i` enter only as the coordinate change of §1. **BP-admissible.**
2. **Partial claims are handled natively.** A gDNA-only message delivers `f̂_g^j = ρ̂_g^j E_g^j / M_j < 1`
   automatically, because `M_j` already contains the source's RNA. This is precisely the property
   `_pin_v`'s docstring wants ("a message carrying gDNA only still delivers `f_g < 1`") — obtained without
   touching the destination. **The deficit `1 − Σ_{c∈A} f̂_c^j` is measured at the source, not guessed at
   the destination.**
3. ⭐ **The reframe is not needed at all.** `ρ̂_c^j` and `M_j` are both in the source's frame, so `r` cancels
   *before* it is ever formed. Job B is achieved structurally rather than by repair. **This explains why the
   pin "cancels `r`" — a composition transfer never needed `r` in the first place.**

**And for a COMPLETE claim the shipped pin already agrees with (★★).** With `A = {g,+,−}`,
`S = r·Σ_c ρ_c^j E_c^i` and `k = M_i/S`, so the delivered `ρ_c^i = ρ̂_c^j M_i / Σ_{c'} ρ̂_{c'}^j E_{c'}^i` —
identical to (★★) whenever `E_c^i = E_c^j`, and otherwise differing only in which frame the normalizer is
evaluated in (the pin's choice, the destination's `E`-frame, is arguably the better one). **No destination
belief appears. The pin is CORRECT for a full claim.**

> **The entire defect is the partial-claim branch.** `S` needs the unsupplied components; the shipped code
> takes them from the destination, the derivation takes them from the source's own mass.

---

## 5. THE SECOND CONFLATION — and it is the one that actually fires

(★) needs `M_j > 0`. **The nodes that emit the false gDNA have `M_j = 0`** — they are splice-junction
boundaries (measured: **96.4 % intron↔exon**) with zero unspliced mass and real spliced flux.

That is not a corner case to patch. It is the correct answer: **a node with no unspliced mass has no
unspliced composition.** It cannot make a composition claim at all.

What it *can* claim is an **absolute RNA density** measured from junction flux. That is a different object:

| claim | what transfers | needs `r`? | pin? |
|---|---|---|---|
| **COMPOSITION** (source has mass) | `f̂^j`, scale-free | **no** — cancels in (★) | **yes**, normalized by `M_j` |
| **DENSITY / measurement** (spliced flux) | `ρ̂_R`, an absolute rate | **yes**, with `Var(log r)` | **NO** |

**ψ already consumes a density claim correctly**: `mo_p = log(cp·E_r/M_i)` is exactly `log f̂_R` under the
coordinate change, so `−½p·(log[(1−f_g)] − mo_p)²` is a proper likelihood on the composition given an
absolute RNA rate. And the solver's own three-stream combine already declares the distinction —
*"(3) SPLICED RNA MEASUREMENT → rna_imp … INDEPENDENT of the composition, so fused separately"*.

**The architecture is right. The plumbing pins both streams through one operator that is valid for only
one of them.** `mo_p` is built from `cp`, and `cp` is the **pinned** density.

### The consequence, in closed form

On those nodes the message supplies RNA only (measured: gDNA-message precision is 0 on 433/433), so with
`E_g = E_r = E` and an incoming RNA density equal to the destination's own total `ρ_tot = M_i/E`:

```
    n_pinned = n·M_i / (og·E + n·E) = ρ_tot² / (ρ_tot + og) = ρ_tot · 1/(1 + f_g^own)
```

**The delivered RNA fraction is `1/(1 + f_g^own)` — a pure function of the destination's own strand-only
self-solve, containing no source information whatsoever.** Verified to 2.1e-16.

**And it is self-fulfilling.** On unstranded data the strand channel is flat, so `f_g^own` sits at the
uninformative ½; the pin reserves half the mass budget for that imaginary gDNA; the RNA message is capped at
`1/(1+½) = 0.667`; the solve reads back ~34 % gDNA — **the number it started from.** The control is exact:

| library | `f_g^own` | budget reserved for gDNA never claimed | false-positive rate |
|---|---|---|---|
| zero-gDNA, **unstranded** | 0.5065 | **33.6 %** | **29.3 %** |
| zero-gDNA, **stranded** | **0.0130** | **1.2 %** | **1.4 %** |

---

## 6. THE CIRCULARITY, RANKED — not all destination dependence is equal

Three places the destination's own state enters the message pipeline. **Only one crosses the line.**

| # | site | what it touches | verdict |
|---|---|---|---|
| **1** | **`_pin_v` partial branch** | the **MODE** | ⛔ **self-confirming. This is the bug.** |
| 2 | P1e `δ = log(M/S)` | a variance, but computed from the contaminated `S` | ⚠ inherited; already a standing debt |
| 3 | M7 DL `b̂²` | the **VARIANCE** only | ✅ deliberate, documented, admissible |

**The categorical distinction:** a destination-dependent **variance** mis-weights a message but cannot
invent a location — a message that agrees with the node moves the fused mode nowhere. A
destination-dependent **mode** manufactures the very answer it is asked to test. M7 is a two-study
random-effects meta-analysis and is a legitimate use of the destination's independent self-estimate;
`_pin_v` is not.

⭐ **Site 2 is already documented as contaminated, in the code, and was scoped around rather than fixed.**
P1e's own comment reads: *"`S` is a COMPLETE budget: `_pin_v`'s partial-claim semantics fill every component
the message does not supply from the node's OWN density. So a shortfall … can be the node's own density
being too low just as easily as the message being wrong — **it does not attribute**."* That is this bug,
seen from downstream, and the response was to charge only one direction of `δ`.

**`variance_ledger.md` §6 says P1e "must SHRINK when the bias strata are diagnosed". This is that
diagnosis.** With `S` built from the message alone, `δ = log(M/S)` becomes an attributable conservation
residual for the first time, and the one-sided scope may no longer be needed.

⚠ **Not a violation:** the reframe `r = dst_face_v / src_face_v[src]` uses the destination's *frame*, but at
`_RHO_ITERS = 1` that frame comes from `_init_belief()`, which is belief-free. This is what P4b established
and it still holds.

---

## 7. HOW BIG IS IT — the partial branch is not a corner case

Measured over all 32 conditions (`gdna_d9_pin_scope.py`):

* **7,211 of 66,150 message-carrying nodes (10.9 %) take the partial branch**, carrying **16.3 % of the
  destination mass**.
* On RNA-only claims the reservation is **33–35 % on every unstranded condition regardless of gDNA level** —
  a constant of the unstranded regime, not a zero-gDNA artefact.
* Stranded conditions: 1.2–16.4 % (the strand likelihood resolves `f_g^own`, so the distortion shrinks with
  it — the causal signature).
* ⭐ **`capture_verystrong` has the highest partial rate of any stratum, 26.6–32.8 %** — which is the *other*
  open regression (VSTRONG enriched-census recovery 0.29). **The two ends of the spectrum may share this
  cause**, and that should be tested rather than assumed.

---

## 8. VERDICT

**Not a typo. Not an architectural collapse. One conflation, with a two-part local fix and suite-wide reach.**

* The pin's **goal** is right and its **full-claim form is correct and BP-admissible**.
* The **partial-claim branch** is a genuine modelling error: it answers "what do the unsupplied components
  hold?" with the destination's guess instead of the source's measurement.
* The deeper issue is that **two kinds of claim are routed through one operator**. A composition claim must
  be pinned (normalized by the source's mass); a density measurement must not be pinned at all. The
  three-stream combine already names the distinction; the plumbing does not honour it.
* It is **worse than it looks** because a second term (P1e) was built on top of the contaminated budget and
  scoped around it, and **better than it looks** because fixing it should make that term attributable.

### The proposal — two changes, each measurable alone

* **P-1 — normalize the composition pin by the SOURCE's own mass.** Replace the destination-own substitution
  with (★): `f̂_c^j = ρ̂_c^j E_c^j / M_j`. Where `M_j = 0` the source has no composition and must emit no
  composition message (the λ-emission gate already exists for exactly this reason).
* **P-2 — do not pin the measurement stream.** Deliver the reframed spliced RNA density unpinned, with its
  `Var(log r)`; ψ already consumes it correctly.

### Pre-registered predictions (so the A/B cannot be read after the fact)

1. zero-gDNA unstranded false positives fall sharply — the 33.6 % reservation is the mechanism;
2. **stranded conditions barely move** (their reservation is already 1.2 %) — if they move a lot, the
   diagnosis is wrong;
3. **VSTRONG moves**, since it has 3× the partial-branch rate;
4. **capture-OFF is the risk**: `_pin_v`'s docstring records that rescaling all three components blindly
   regresses capture-OFF 3.6×. ⚠ That is a *different* operator from (★) — but it is the stratum to watch;
5. P1e's firing set and its bias share change; re-measure `δ` before deciding whether the one-sided scope
   is still needed.

**Gates:** full 32-condition A/B at refit 0 and 1 (the pin is on every message in every library, so
`gdna_none` — 9 of 32 — cannot decide it); held-fixed `z2`; then goldens.

---

## 9. ⭐ P-2 IS IMPLEMENTED AND MEASURED (2026-07-27, uncommitted). EVERY PREDICTION HELD.

**The change is three lines and one rename.** `_pin_v`'s result is no longer assigned back over the
delivered densities; it is bound to `pin_g/pin_p/pin_n` and used **only** by the DL mismatch gaps. The
delivered `tg/tp/tn` stay exactly as the source measured and the reframe delivered them.

**Why that is the whole fix, verified before editing:** on all 433 affected exons the λ (composition) stream
is **dead** (`c_tau = 0` on 433/433, the emission gate already firing) and the gDNA measurement is dead
(0/433) — **only the RNA measurement is live (433/433)**. And `tlam`/`tth` are scale-free, so the pin's
common factor cancels from them *identically*. So the pin's only live consumer was the measurement mode —
the one claim that must not be composition-normalised. Nothing else in the packet moves.

| stratum | refit=0 base → pin fix | Δ | refit=1 base → pin fix | Δ |
|---|---|---|---|---|
| **ALL 32** | 0.0841 → **0.0788** | **−0.0053** (11 b / 5 w) | 0.0565 → **0.0525** | **−0.0040** (8 b / 7 w) |
| stranded | 0.0311 → 0.0307 | −0.0005 | 0.0264 → 0.0258 | −0.0007 |
| unstranded × capON | 0.1658 → **0.1513** | −0.0145 | 0.1201 → **0.1047** | −0.0155 |
| VERYSTRONG | 0.1823 → **0.1705** | −0.0118 | 0.0969 → 0.0937 | −0.0033 |
| capture OFF | 0.0350 → 0.0339 | −0.0010 | 0.0195 → 0.0204 | **+0.0009** |
| **`gdna_none`** | 0.0931 → **0.0667** | **−0.0263** | 0.0312 → **0.0152** | **−0.0160** |

The two conditions the mechanism was diagnosed on, at refit=0: **0.2934 → 0.1715** and **0.3318 → 0.1891**.

**The pre-registered predictions, scored:**

1. ✅ zero-gDNA unstranded falls sharply — −42 % on both diagnosed conditions.
2. ✅ **stranded barely moves — the falsification test. −0.0005 / −0.0007**, and the individual stranded
   zero-gDNA conditions are unmoved to 4 decimals. Had they moved, the diagnosis would have been wrong.
3. ✅ VSTRONG moves (−0.0118 / −0.0033), as its 3× partial-branch rate predicted.
4. ⚠ **capture-OFF is the risk, and it materialised**: +0.0009 at refit=1 (1 b / 5 w), concentrated in
   **unstranded × capture-OFF × gDNA-bearing** (gdna100 nrna_none +0.0181, gdna300 nrna_none +0.0078).
   Stranded capture-OFF is untouched (+0.0000/+0.0001).

⚠ **Honest counter-reading: mass-weighted correlation falls** −0.0070 (r0) / −0.0089 (r1) while mwae
improves. Accuracy up, tracking slightly down.

⭐ **And pass-0 ITSELF improved 6.3 % (0.0841 → 0.0788)** — the harder gate, on a solver the ROADMAP had
recorded as near its prior-free ceiling.

## 10. ⭐ THE RESIDUAL IS EXACTLY P-1's TARGET — measured, not assumed

The capture-OFF regression is the mirror of the bug just fixed: with nothing to limit it, an RNA-only
message now asserts *too much* RNA. That is precisely what P-1 (normalize by the SOURCE's own mass) exists
to bound — and the populations separate perfectly:

| condition | RNA-only destinations | **sources with `M_src > 0`** |
|---|---|---|
| `none · ss0.50 · nrna_none · capOFF` (the FIXED case) | 433 | **64/863 = 7.4 %** |
| `gdna300 · ss0.50 · nrna_none · capOFF` (the REGRESSING case) | 18 | **36/36 = 100 %** |
| `gdna100 · ss0.50 · nrna_none · capOFF` | 18 | **36/36 = 100 %** |
| `gdna300 · ss0.50 · nrna_present · capOFF` | 22 | **44/44 = 100 %** |

**Where the fix helps, the sources have no mass and therefore no composition — an unbounded density claim is
the correct statement. Where it regresses, every source has mass and can state its own share, so P-1 bounds
it.** The two halves of the derivation are complementary, not redundant, and P-1 targets the residual and
nothing else.

## 11. ⛔⛔ P-1 WAS IMPLEMENTED AND **REVERTED** — and it CORRECTS §4 of this document

It hit **every** pre-registered target and still lost, decisively:

| stratum | refit=0 Δ | refit=1 Δ |
|---|---|---|
| ✅ **unstranded × capOFF × gDNA-bearing** (its target) | **−0.0022** | **−0.0020** |
| ✅ capture OFF overall | −0.0015 | −0.0011 |
| ✅ `gdna_none` (must not move) | −0.0009 | −0.0000 |
| ⛔ **capture ON** | **+0.0264** | **+0.0274** |
| ⛔ **unstranded × capON** | **+0.0378** | **+0.0415** |
| ⛔ stranded | +0.0062 | +0.0052 |
| ⛔ **ALL 32** | **+0.0119** (3 b / 11 w) | **+0.0121** (3 b / 10 w) |

corr −0.033. Reverted; the revert is **bit-identical 32/32** to the P-2 state at both refit settings.

### ⭐ THE LESSON, and it is worth more than the fix would have been

**"Composition transfers" is a WEAKER premise than "density transfers, reframed".** §4 above derives (★)
from the imputation premise `f^i = f^j` and is correct *as mathematics* — but that premise is
**empirically false across a capture step**. An exon and its flanking boundary do not share a composition:
the boundary sits on the capture SLOPE, measured at **0.125× the exon and 2113× the intron at verystrong**
(`ROADMAP` §11). Replacing the reframe with a mass ratio therefore discards the enrichment information `r`
carries, and the capture-ON regime pays for it immediately.

⚠ **And this re-reads an old result.** The record said substituting the ORACLE capture step for `r` "buys
≈ 0", which made `r` look inert. It was inert *only while `_pin_v` was cancelling it*. Once the pin stopped
rewriting the delivered density (P-2), **`r` became load-bearing again** — so that measurement does not
license removing the reframe, and §4's step 3 ("the reframe is not needed at all") is **withdrawn**: it is
true of a composition transfer and the composition transfer is not the better estimator.

**What survives of §4:** the diagnosis (the destination's belief must not set a message's mode), the
severity ranking of §6, and P-2, which is measured and landed. **What does not:** the claim that (★) is the
right replacement operator.

**→ The P-2 residual (capture-OFF +0.0009 at refit=1, concentrated in unstranded × capture-OFF ×
gDNA-bearing) is therefore OPEN, and (★) is not the answer.** Whatever bounds a partial RNA-only claim must
do so without discarding the reframe. ⛔ Do not re-derive (★) — it is implemented, measured and refuted.

⚠ **Status:** P-2 implemented, measured, uncommitted; **21 goldens pending (LAST)**. P-1 derived and
targeted, **not implemented**.

---

## 12. ⭐⭐⭐ THE RESIDUAL IS DIAGNOSED — AND §10's MECHANISM IS WRONG IN EVERY PARTICULAR

**Measured 2026-07-27, refit=0, on the tree that ships P-2.** Scripts `scratchpad/p2r_{a,b,c,d,e,f}_*.py`;
working notes `scratchpad/P2_RESIDUAL_NOTES.md`. §10 asserted the residual was *"an RNA-only message
asserting too much RNA"*, i.e. a **partial** claim **under**-calling `f_g`. That was an inference from
§10's source-mass table, never a measurement. It is false:

| §10 predicted | measured on the regressing stratum |
|---|---|
| partial (RNA-only) claims | **0.2 %** of the harmed mass is partial |
| `f_g` **UNDER**-called | **96–99 %** of the regression mass is **OVER**-calling |
| the message asserts too much **RNA** | it asserts too much **gDNA**: `e^moG` = 0.42–0.77 against an oracle `f_g` of 0.008–0.043 |

An off-simplex reading was pre-registered and also falsified: `mo_R > 0` carries **3.3 %** of the stratum's
error mass, on boundaries, in the *over*-calling direction.

### 12.1 The residual is on the branch where the pin was CORRECT — but there is nothing there to correct

Split every message by whether the pin's budget borrowed the destination's own density for a component the
message did not supply (a structurally dead strand has own density 0 and lends nothing):

| condition | Δ from P-2 | **CLEAN** | CONTAMINATED |
|---|---|---|---|
| gdna100 ss0.50 nrna_none capOFF | **+0.0148** | **+0.0129** | +0.0019 |
| gdna300 ss0.50 nrna_none capOFF | +0.0076 | **+0.0067** | +0.0009 |
| **none ss0.50 nrna_none capOFF** (P-2's win) | **−0.1219** | +0.0000 (**0 clean nodes**) | **−0.1219** |
| gdna100 ss0.50 nrna_none capON | −0.0012 | −0.0028 | +0.0016 |

P-2's entire `gdna_none` win is contaminated-branch — exactly as §5 diagnosed — and the residual is
clean-branch, **entirely EXONS** (net +0.0042 of the stratum; intron/boundary/intergenic ≈ 0).

**So the obvious repair is to restore the pin only where it is BP-legal. It was implemented and it does not
work**, and the reason is measurable: on the clean branch the conservation violation is
**|δ| = 0.073**, against **|δ| = 0.46** on the contaminated branch. There is almost nothing there to
correct. Two arms, both run:

| arm | suite r0 | suite r1 | unstr × capOFF × gDNA | unstr × capON r1 |
|---|---|---|---|---|
| **`w = σ_cm²/(αᵀΣα + 1/n_dst)`** — the conditional-mean shrinkage, derived, no new constant | −0.0000 | +0.0000 | −0.0002 / −0.0001 | +0.0003 |
| **`w ≡ 1`** — restore the pin outright wherever it is BP-legal | −0.0001 | **+0.0002** | −0.0004 / −0.0004 | **+0.0009** |

The full restore recovers **11 %** of the residual and costs more than it recovers. ⛔ **Do not re-try a
pin-shaped fix for this residual.** *(The derived form is worth recording even though it is inert: the pin
and P-2 are the `w = 1` and `w = 0` endpoints of one shrinkage, `E[ε_cm|δ] = δ·σ_cm²/(αᵀΣα + 1/n_dst)`,
built entirely from quantities P1e already computes. It is inert because the model attributes the violation
to component-specific error, not to the common scale — which is the correct answer here.)*

### 12.2 ⭐⭐ WHAT THE PIN WAS MASKING — and it is worth 20× the residual

On the **contaminated** branch the old pin applied a **×0.648** median shrink to every delivered component
of a message whose budget it had completed from the destination's own density. `pG` dominates `pP` on those
exons (5.9–13.9 vs 0.35–5.6), so a common down-shift of both level modes pulls `f_g` down — which is why
removing it reads as an over-call. **The shrink was masking a separate and much larger defect.**

The delivered **gDNA** level into exons is wrong by **5.4×**, and the reframe is why. `r =
ρ_tot(dst)/ρ_tot(src)` exists to carry the hybrid-capture step; between an INTRON or boundary (the source on
**100 %** of these edges) and an EXON in a gDNA-bearing library the total-density ratio is dominated by
**RNA**, not capture — so it multiplies the gDNA claim by the destination's RNA content. gDNA is uniform, so
across a capture-OFF hop its correct transfer is `r = 1`:

| condition | median `r` | err REFRAMED | err **UN**-REFRAMED | mass improved by dropping `r` |
|---|---|---|---|---|
| gdna100 ss0.50 nrna_none **capOFF** | 1.48 | 0.732 dec | **0.170** | **80.2 %** |
| gdna300 ss0.50 nrna_none **capOFF** | 1.44 | 0.556 | **0.110** | 74.2 % |
| gdna100 ss0.99 nrna_none **capOFF** | 1.26 | 0.731 | **0.154** | 72.5 % |
| gdna100 ss0.50 nrna_none **capON** | 2.82 | **0.639** | 1.608 ✗ | 10.9 % |

*(err = mass-weighted mean `|log10(claim/true ρ_g)|` over exon destinations.)* **The sign flips exactly at
capture**, which is the physics: at capture-ON `r` IS the enrichment step and dropping it is fatal.

**The prize, sized by a full A/B of `r_g ≡ 1` on the plain reframe** (prior-free — "gDNA density is
uniform" — but wrong under capture, so NOT shippable):

| stratum | r0 | r1 |
|---|---|---|
| ⭐ **unstranded × capOFF × gDNA** (the residual's own stratum) | **0.0495 → 0.0146 (−0.0349), 6/6 better** | **0.0313 → 0.0077 (−0.0237), 6/6** |
| capture OFF (all 14) | 0.0339 → **0.0161** (−0.0178), 9 b / 0 w | 0.0204 → **0.0080** (−0.0124), 9 b / 0 w |
| stranded × capOFF | −0.0028 | −0.0025 |
| `gdna_none` | −0.0031 | −0.0022 |
| ⛔ **unstranded × capON** | **+0.1728** | **+0.2069** |
| ⛔ ALL 32 | +0.0457 | +0.0591 |

> **The P-2 residual (+0.0009) is a 5 %-sized symptom of a defect worth −0.0349 on its own stratum.**
> Fixing the frame of the gDNA component is the single largest identified lever on capture-OFF, and it is
> *bounded by the same measurement* — any operator must reduce to `r` at capture-ON.

⚠ **This re-opens `ROADMAP` §11** ("REFUTED, do not build: per-channel enrichment ratios at the boundary
face"). That refutation rested on substituting the oracle capture step for `r` buying ≈ 0 — measured **while
`_pin_v` was cancelling `r`**, which §11 of this document already flags as invalidating. The refutation must
be re-taken on the P-2 tree.

**What is NOT yet answered, and is the next task's whole difficulty:** what sets `r_g` prior-free. The two
endpoints are `r` (right at capture-ON) and `1` (right at capture-OFF), and pass-0 cannot tell them apart
from one edge — the population gDNA-density landscape is exactly the object that can, which puts this at the
pass-0 / Phase-2 boundary rather than inside either. One belief-free candidate exists and is untested: the
frames already come from `_init_belief()` at `_RHO_ITERS = 1` (P4b), so a **gDNA-specific** frame
`ρ_g^init(dst)/ρ_g^init(src)` is available at no worse BP standing than the current `r`.

### 12.3 Verdict on the residual

**Not fixable at the pin, and it should not be chased there.** P-2 stands as measured (suite −0.0053 r0 /
−0.0040 r1); its +0.0009 capture-OFF cost at refit=1 is the honest price of removing a self-confirming mode,
and it is now a *localised, sized, and named* pointer to the gDNA frame defect rather than an open question
about the pin. The tree is unchanged by this investigation — **bit-identical 32/32** to the P-2 baseline.
