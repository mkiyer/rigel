# The gDNA factory, absolute density, and why the SHIFT is invalid into an exon

**Status:** derived. Supersedes the "vacuous-source variance" options in `message_arithmetic_reconciliation.md`
§E.4 — no such factor is needed. **Date:** 2026-07-22. **Scope:** pass-0 message propagation.

---

## 1. The owner's locus model

> A **locus** is a genomic span of overlapping transcripts, bookended on both sides by **intergenic** regions.
> An intergenic region holds a **measured** amount of pure DNA — guaranteed true, no solve, invariant. It is a
> **factory** of observed gDNA that propagates into the locus at counting precision, with no penalty. It is
> also the **sink**: no RNA crosses it and a message that reaches it dies — which is also what stops gDNA being
> averaged across a whole chromosome.
>
> The intergenic↔exon **boundaries** are pure DNA too, just like the intergenic regions; if such a boundary is
> partially enriched by a probe on the adjacent exon, that is recoverable through a **shift**, because both
> endpoints are pure DNA and there is **no solve** between them.
>
> **Density is absolute and can always be imputed across a boundary.** A TSS intergenic–exon boundary has no
> information about RNA, so it must send a message with **zero precision on the RNA component** onto the exon —
> but its **gDNA component has precision**.

That last sentence is the design, and §7 shows the current code violates it in a specific, fixable way.

---

## 2. Two currencies, two validity conditions

A node's per-component state has two readings, and they are **not** interchangeable:

| currency | definition | what it is |
|---|---|---|
| **density** `ρ_c = M_c/E_c` | mass per effective length | **ABSOLUTE.** Comparable between nodes. Physically continuous for gDNA (a genomic background) and for nascent RNA (a transcription-unit property). |
| **composition** `f_c = M_c/ΣM_c'` | share of the node's own unspliced pool | **RELATIVE.** Meaningless without the node's component set. |

The solver has exactly two message vehicles, and each is valid under a different, checkable condition:

**SHIFT (composition transfer).** `M_c^dst = ρ_c^src·E_c^dst·g_c`, `mode_c = log(M_c^dst/Σ M_c'^dst)`.

* ✅ **Enrichment cancels identically** — `e(src)` is common to every `ρ_c^src` and therefore to the normalizer;
  `e(dst)` never appears. Correct across an arbitrary capture cliff.
* ⛔ **Requires component-set equality**: the source must be able to supply **every component the destination
  possesses**. Any component the destination has and the source cannot supply is implicitly asserted to be
  **absent** — the normalizer omits it.

**DENSITY (absolute imputation).** `mode_c = log(ρ_c^src·E_c^dst / md_dst)`, normalized by the destination's own
**observed** total.

* ✅ **Needs no knowledge of the destination's other components** — `md_dst` is measured, so whatever else the
  destination holds is accounted for implicitly. This is what makes it the right vehicle for a *partial* source.
* ⛔ **Enrichment does NOT cancel**: it substitutes `ρ̂_c(src)` for `ρ̂_c(dst)`, so it is unbiased only when
  `e(dst) = e(src)`. §6.

> The two failure modes are complements. The shift cannot handle a **missing component**; the density mode
> cannot handle an **enrichment change**. Every pass-0 mode bug in this thread is one of these two.

---

## 3. The theorem: no unspliced crossing can ever carry an exon's MATURE

An exon region's unspliced pool is `{g, ν, μ}`. A boundary's unspliced **crossing** pool is `{g, ν}` — mature
is spliced, so it leaves the genome contiguously and cannot appear as an unspliced junction crossing (this is
the same fact that gives the `c_b` peel).

> **Therefore every `boundary → exon` edge is component-UNEQUAL, unconditionally.** Not sometimes, not when a
> strand gate closes — *always*. The destination possesses a component (μ) that no unspliced crossing can carry.

**Corollary: the SHIFT is invalid into an exon region, period.** Its normalizer omits μ, so it asserts the
exon's mature is absent — and the more mature the exon holds, the harder it asserts the exon is pure gDNA.

The current `use_shift_g` predicate tests **strand continuity** (`fp/fn`), which is the wrong question. Strand
continuity says nascent may flow; it says nothing about mature. This is why the `c_b` term exists at all — it
is a patch that restores μ *when the boundary happens to measure spliced fragments*. When `S_B = 0`, `c_b = 0`
and the shift asserts "no mature" at full confidence, with no term left to object.

### 3.1 Exact numerical demonstration

`tests/calibration/test_bp_solver.py::_mature_exon_chain` is a closed-form chain — `intron⁺ | exon⁺ | intron⁺`,
`ρ_g = 0.5`, `ρ_m = 1.0`, δ-function FLs (gDNA 300, RNA 200), regions 2000 bp:

```
  E_g(contained) = 1701.0     E_r(contained) = 1801.0     E_g(crossing) = 150.0     (computed, not assumed)
  exon: u_pos = ρ_g·E_g/2 + ρ_m·E_r = 425.2 + 1801.0 = 2226.2 ,  u_neg = 425.2  ⇒  md = 2651.5
  exon gDNA mass = 850.5                     ⇒  TRUE f_g(exon) = 850.5/2651.5 = 0.320762
```

The boundary flanking the exon carries a gDNA-only crossing at `ρ_g = 0.5` (exact, by construction). The two
vehicles, both fed that same correct source density:

| vehicle | computation | result | truth |
|---|---|---|---|
| **DENSITY** | `ρ_g·E_g(dst)/md = 0.5·1701.0/2651.5` | **0.320762** | 0.320762 ✅ **exact to 6 d.p.** |
| **SHIFT** | `M_g/(M_g+M_ν)`, `M_ν = 0` (no mature can cross, none measured) | **1.000000** | 0.320762 ❌ |

(Verified against `region_eff_length` / `boundary_side_eff_length` directly, not hand-derived.)

The density mode recovers the truth to the last digit; the shift is wrong by the entire mature content. Feeding
that shift into the solver at the certainty the factory unlock supplies pinned the exon at **f_g = 0.9616**
(`test_tau_gag_fix_deconvolution_prediction_stays_gated`, which requires 0.2 < f_g < 0.8).

**The revert in §E.8 was therefore misattributed.** The factory unlock was not the defect; it merely supplied
enough precision to expose an already-wrong mode. `struct_lock` was innocent.

---

## 4. Per-component precision — the owner's statement, formalized

> *"The intergenic–exon boundary at the TSS has no information about RNA and has to send a message with zero
> precision on the RNA component onto the exon. But the gDNA component has precision."*

This separates two things the code currently conflates:

| | claim about the component's **density** | claim's **precision** |
|---|---|---|
| gDNA at a pure-DNA seam | `ρ_g` = the node's own mass/eff — **measured, no solve** | counting: `1/(1/n + …)` — **real** |
| RNA at a pure-DNA seam | *none* — the seam observes no RNA | **exactly 0** |

**Zero precision is not a claim of zero density.** The shift conflates them: by omitting RNA from its
normalizer it converts *"I know nothing about RNA"* into *"there is no RNA"* at the gDNA factor's full
confidence. The density mode keeps them separate — it states only the gDNA density and lets the destination's
own `md` absorb everything else.

This is exactly the emission architecture already in place (`emit_p`/`emit_n` are gated by strand continuity,
so a seam emits no RNA factor). The single defect is that the surviving **gDNA** factor is computed by the
wrong vehicle.

---

## 5. The pure↔pure edge: intergenic region ↔ intergenic–exon boundary

Both endpoints are composition-certain pure gDNA (`f_g ≡ 1`). Consequences:

1. **Component-set equality holds trivially** (`{g}` on both sides) ⇒ the **shift is exact**, and it is also
   vacuous as a *composition* statement: 1 → 1. There is **no solve** on this edge.
2. **The density ratio is an IDENTIFIED enrichment measurement.** With composition known on both sides,
   ```
        log ρ̂_g(B) − log ρ̂_g(I)  =  log( e(B)/e(I) )        exactly
   ```
   Nothing else can contribute to the gap. This is the owner's "partial enrichment from a probe on the adjacent
   exon, determined through a shift" — and it is the **only** place in the chain where an enrichment ratio is
   identified from data alone (§6).
3. **No σ²_transfer penalty is warranted.** That term exists to damp *model error in a composition transfer*.
   Here there is no composition being transferred, hence no model error to damp. "No penalty for that."

So certainty must survive the intergenic→seam hop. `struct_lock = locked & is_region` blocks it purely because
the seam is a boundary, which §3.1 now shows was never the problem.

---

## 6. The enrichment obstruction — stated precisely, and its boundary

Write `ρ̂_c(x) = e(x)·ρ_c(x)`; capture is nucleic-acid-agnostic at pass-0, so `e` is component-independent.
Assume the gDNA background is genomically uniform (`ρ_g ≡ ρ_bg` — the basis of `measure_background`/ρ_bg). Then
for any imputation src→dst:

```
    log f_g^true(dst)  =  _mode_density(ρ̂_g(src), E_g(dst), md(dst))  +  log( e(dst)/e(src) )
```

> **The density mode's error is exactly the log enrichment ratio.** It is a **bias in the MODE**, not a
> variance. σ²_transfer's `(μ_proj[dst] − μ_proj[src])²` is that bias being *paid for as uncertainty* — the
> same pathology the density→shift fix already corrected once.

**Identifiability.** The observed total-density gap decomposes as

```
    μ_proj(dst) − μ_proj(src)  =  log(e(dst)/e(src))  +  log(composition/content ratio)
```

so using the gap to correct the mode is **circular** in general — the second term is what we are solving for.
It is identified only when src and dst have the **same composition**, which is exactly the pure↔pure edge of
§5. Hence:

| edge | ratio identified? | consequence |
|---|---|---|
| intergenic ↔ intergenic-exon boundary | ✅ **yes** — both pure gDNA | exact enrichment measurement; anchor |
| pure-gDNA seam → RNA-bearing exon | ❌ no | density mode carries an unknown bias `log(e(dst)/e(src))` |

**One-sidedness is what survives.** Even unidentified, the *sign* is often structural: a seam that straddles
intergenic territory is enrichment-dominated by its depleted half, so `e(dst) ≥ e(src)` and the density mode is
a genuine **lower bound** on `f_g(dst)`. In the reverse direction (a boundary flanking an *enriched* exon
imputing into a depleted intron) the sign flips and the mode over-claims — measured "fractions" of **186, 294,
324**. So a density floor may be trusted **only in the direction of non-decreasing enrichment**, and that is a
structural, checkable condition, not a fitted one.

---

## 7. What this implies for the code

The corrections are *removals*, and they compose:

1. **`use_shift_g` must test component-set equality, not strand continuity.** By §3 a destination exon is
   always component-unequal ⇒ **boundary→exon uses the DENSITY mode** for gDNA. (The old comment at
   `bp_solver.py` "an EXON endpoint keeps the DENSITY mode … the correct permanent design" was right; the
   §E-era shift work partly overrode it.)
2. **`c_b` and the mature graft both become unnecessary for the gDNA factor.** Both exist to reconcile μ inside
   a *composition* normalizer. The density mode never forms that normalizer — `md_dst` already contains the
   exon's mature. The peel `+c_b` remains meaningful only for exon→boundary, where the destination genuinely is
   mature-free. **This is the consolidation**: one vehicle chosen by one predicate, and two correction terms
   disappear rather than being maintained.
3. **`struct_lock` becomes structural** (`locked`, dropping `& is_region`) — the intergenic factory reaches the
   seam. Safe **only together with (1)**; alone it amplifies the §3 mode bug (measured: f_g 0.32 → 0.96).
4. **No σ²_transfer on a pure↔pure edge** (§5.3): no composition is transferred there.
5. **The RNA channel already emits nothing at a seam** — no change needed. The fix is entirely in the gDNA
   vehicle (§4).

**Staging** (each gated on the `gdna_none` delta + the suites):

| # | step | risk |
|---|---|---|
| F1 | `use_shift_g` ⇒ component-set equality (boundary→exon on the density mode) | behavioural; the mature graft/`c_b` interaction must be removed in the same step |
| F2 | `struct_lock` structural + no σ²_transfer on pure↔pure | inert without F1; unlocks the factory |
| F3 | one-sided density floor gated on non-decreasing enrichment (§6) | needs the structural sign test |
| F4 | migrate the identified enrichment ratio (§5.2) from σ²_transfer's variance into the **mode** | this is the corrected R4 |

**R4 is now specified rather than open:** the cliff term is not deleted, it is **moved from the variance into
the mode**, and only where the ratio is identified.

---

## 8. Summary

| question | answer |
|---|---|
| Does pass-0 need a "vacuous-source" variance (uniform 1/12, Popoviciu ¼, population prior)? | **No.** The infinity was never the obstacle — a locked node is composition-*certain*, and the real defect was the mode it was certain about. |
| Why did unlocking the factory pin an exon at 0.96? | The shift asserts an exon's **mature** is absent. Certainty amplified an already-wrong mode; `struct_lock` was innocent. |
| Which vehicle is right into an exon? | **DENSITY** — exact on the closed-form fixture (0.3208 vs 0.3208), where the shift gives 1.0000. |
| What does "zero precision on RNA" mean? | Emit **no** RNA factor — never a normalizer that omits RNA, which converts *no information* into *no RNA*. |
| Where is enrichment identified? | Only on **pure↔pure** edges (intergenic ↔ intergenic-exon boundary), where composition is known on both sides. |
| What is σ²_transfer really? | An uncorrected **mode bias** (`log e(dst)/e(src)`) being paid for as variance. |

---

## 9. What a +stranded TSS exon can solve with — the two-message question

**Setup** (owner). A `+`-stranded first exon `E`: left flank `B_L` = intergenic–exon seam (gDNA only), right
flank `B_R` = exon–intron junction (gDNA + nascent + a measured spliced channel). Enrichment rises stepwise:
`e(I) → e(B_L) → e(E)`. *"At solve time the exon sees RNA from at least one boundary and gDNA from both. Can't
we compute the compositional shift?"*

### 9.1 The composition closes on ONE scalar — verified

`md(E)` is observed, so `f_g(E) = M_g(E)/md(E)` has exactly one unknown. Writing `ρ̂ = π·ρ` (π = positional
probe coverage, §9.4) and using `ρ_g ≡ ρ_bg` (uniform background, `π(I) = 1` by definition of intergenic):

```
    f_g(E)  =  π(E) · ρ_bg · E_g(E) / md(E)
```

**Measured across 12 toy conditions** (gDNA ∈ {1,3,10} × capture {off,on} × nascent {50,200}): this reproduces
the oracle `f_g` to **3 decimal places in 12/12 cells** (0.184, 0.240, 0.505, 0.503, 0.384, 0.472, 0.758, 0.754,
0.674, 0.748, 0.911, 0.911). The parameterization is therefore **complete**: the intergenic factory supplies
`ρ_bg`, the destination supplies `md`, and everything else reduces to one scalar — **the exon's own probe
coverage π(E)**.

### 9.2 Yes — but the right operation is combining DENSITIES, not compositions

The owner's construction is correct and is a real architectural finding. `B_L` supplies gDNA; `B_R` supplies
gDNA **and** RNA; between them the exon's whole component set `{g, ν, μ}` is covered. So the "missing
component" that invalidates the shift for a *single* message (§3) is repaired by using **both flanks jointly**.

But the operation must be: **pool per-component densities from both flanks first, and normalize once.** Today
the solver normalizes *each message separately* into a λ-factor and then combines the factors (`_comb`, a
precision-weighted geometric mean of log-fractions) — i.e. it averages **compositions**. That is what lets
`B_L` contribute the claim "f_g = 1": a composition formed from a partial source. Pooling densities first and
normalizing once makes `B_L`'s zero-precision-on-RNA mean *"I add nothing to the RNA total"* instead of
*"the RNA total is zero"*.

### 9.3 What the two messages do NOT give: the flanks sit at different enrichments

Combining across flanks does **not** cancel π, because `B_L` and `B_R` are at different coverage levels. The
residual unknown is the ratio `π(E)/π(B_L)`. Measured (oracle-anchored):

| capture | π(E)/π(B_L) across gDNA ∈ {1,3,10} | uncorrected `f_g` vs oracle |
|---|---|---|
| off | 1.61, 1.39 / 1.33, 1.18 / 1.00, 1.06 | 0.115 vs 0.184 … 0.709 vs 0.748 |
| **on** | **2.99, 3.15 / 2.42, 2.55 / 2.12, 2.16** | **0.313 vs 0.758** (a 2.4× under-call) |

**This is the whole of the enriched-gDNA under-call**, quantified: the seam's density claim is right about
*gDNA* and wrong about *where it is being imputed to*.

**The encouraging structure:** the ratio is essentially **invariant to RNA content** — doubling nascent 50→200
moves it by <6 % in every cell (2.99→3.15, 2.42→2.55, 2.12→2.16, 1.33→1.18). It is a property of **capture
geometry, not of composition**, which is exactly what makes it learnable from pure-gDNA structures
(intergenic regions and seams) without a solve — the factory idea extended from a *level* to a *gradient*.

⚠️ It is **not** a universal constant (2.1–3.2 under capture), and its drift with gDNA level is unexplained by
geometry — the probe's `e(B_L)` estimate here is the seam's *total* density and may carry contamination
(it reads below the intergenic at capture-off, which pure geometry does not predict). **Resolve before
building on it.**

### 9.4 Capture is NOT nucleic-acid-agnostic — and the shift survives anyway

`bp_solver` states "capture is nucleic-acid-agnostic ⇒ the crossing is a discontinuity in gDNA AND RNA".
The simulator disagrees: `toy_inject.CAP` carries `gdna_split_penalty = 0.2`, so gDNA is captured less
efficiently than RNA. Factorize `e_c(x) = π(x)·γ_c` (positional coverage × a component constant). Then:

* **The SHIFT still cancels enrichment** — `e_c(dst)/e_c(src) = π(dst)/π(src)` is component-independent, so
  the γ's drop out of the composition ratio. The §2 validity claim is unharmed.
* **The FACTORY still works in absolute terms** — `ρ̂_g(E) = π(E)·γ_g·ρ_bg` and `ρ̂_g(I) = γ_g·ρ_bg`, so
  `ρ̂_g(E) = π(E)·ρ̂_g(I)`: **γ_g cancels against the intergenic reference.** The factory never needs to know
  the gDNA capture penalty.

So the assumption is false but **harmless for both vehicles** — provided γ_c is position-independent. Worth
correcting in the comment, and worth checking on real data, where a probe panel's gDNA/RNA efficiency split
may well vary by target.

### 9.5 Answer

> **Yes** — the exon's composition is fully determined, and the two messages jointly cover its component set.
> The fix is to pool **densities** per component across both flanks and normalize **once**, rather than
> averaging two separately-normalized compositions. What that does *not* deliver for free is the coverage
> gradient `π(E)/π(B_L)` — measured at **2.1–3.2× under capture**, and the entire residual under-call. It is
> RNA-invariant, hence learnable from pure-gDNA structure, which is where the factory should go next.
