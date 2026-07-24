# The junction enrichment-scaling solve — derivation + implementation plan

**Status:** derivation, not implemented. **Date:** 2026-07-23. **Owner's design**, formalized.
**Scope:** the splice-junction boundary and its two flanking messages in pass-0.
**Companion:** `message_layer_derivation.md` (§6.5 the δ taxonomy, §12 the relay + sink model, §13 N4).

---

## 0. The problem in one paragraph

Hybrid capture creates an **enrichment cliff**: a probe-covered exon can sit ~10³× above an off-target intron,
and the junction boundary between them is *partially* enriched depending on probe placement. A message carries
a **density**, so it is only transportable if we know the enrichment ratio between source and destination.
Total density measures that ratio — **but total density itself depends on composition**, because gDNA and RNA
have different fragment-length distributions and therefore different effective lengths. So the ratio we need
to transport a composition depends on the composition we are trying to find.

**The owner's resolution:** at a junction we have enough independent observations to break it — the exon
message, the boundary's own **measured spliced** density, and the intron message. Solve for the three total
densities, form the two ratios, scale both messages into the boundary's frame, then fuse and solve.

---

## 1. Setup and notation

One junction boundary `B`, exon flank `E`, intron flank `I`. All densities are **observed**
(enrichment-included); `e(x)` is the capture efficiency at `x`.

| node | unspliced pool | observed | effective lengths |
|---|---|---|---|
| `E` exon region | `{g, R}` — one RNA (§13.0) | mass `M_E`, strand counts | `E_g(E)`, `E_r(E)` contained |
| `B` junction boundary | `{g, R_cont}` + a **measured** spliced channel | crossing mass `M_B`, spliced mass `S_B` | `E_g(B)`, `E_r(B)` crossing; `E_spl(B)` |
| `I` intron region | `{g, R_cont}` | mass `M_I` | `E_g(I)`, `E_r(I)` contained |

**The mass identity** at every node — the observation constraint (§2 of the companion):

```
    M_x = ρ_g(x)·E_g(x) + ρ_R(x)·E_r(x)
```

**Total density** — nucleic acid per base, the quantity an enrichment ratio must be formed from:

```
    ρ_tot(x) = ρ_g(x) + ρ_R(x)  =  M_x·[ f_g(x)/E_g(x) + (1−f_g(x))/E_r(x) ]
```

The bracket is why `ρ_tot` **depends on composition**: it is a composition-weighted harmonic blend of the two
effective lengths, and it is *not* `M_x/E_anything`.

---

## 2. The key lemma — a density RATIO is enrichment-free, a density is not

```
    ρ_g(x) = e(x)·ρ̃_g(x),   ρ_R(x) = e(x)·ρ̃_R(x)          (capture is component-agnostic at pass-0)

    ⇒   k(x) ≡ ρ_g(x)/ρ_R(x) = ρ̃_g(x)/ρ̃_R(x)              e(x) CANCELS
```

**`k` is transportable where the underlying content ratio is shared; `ρ` alone is not.** This is the lever
that breaks the circularity, and it is why the intron can resolve the boundary's composition without knowing
either node's enrichment.

`f_g` is **not** transportable, even between component-matched nodes, because

```
    f_g(x) = k(x)·E_g(x) / ( k(x)·E_g(x) + E_r(x) )
```

depends on the node's own effective lengths — crossing vs contained. **Transport `k`, then re-form `f_g` in
the destination's frame.** (Copying `f_g` across a frame change is a known past defect.)

---

## 3. The solve — three total densities from the available observations

### 3.1 The boundary's composition, from the intron

`B`'s unspliced crossing and `I` carry the **same component set** `{g, R_cont}` — gDNA is genomically uniform
and the contiguous RNA is continuous into the intron. So by §2, `k(B) = k(I)`, and with the mass identity:

```
    ρ_R(B) = M_B / ( k(I)·E_g(B) + E_r(B) )
    ρ_g(B) = k(I)·ρ_R(B)
    f_g(B) = k(I)·E_g(B) / ( k(I)·E_g(B) + E_r(B) )
```

**No enrichment enters** — `M_B` is the boundary's own observation and `k(I)` is enrichment-free. This is the
formal content of the owner's *"if the boundary has no belief, use the intron's split."* It is not an
approximation of convenience: it is the component-matched condition of §6.5, and it transports the **ratio**,
not the fraction.

### 3.2 The three total densities

```
    ρ_tot(E)          = M_E·[ f_g(E)/E_g(E) + (1−f_g(E))/E_r(E) ]        f_g(E) from the exon's own solve
    ρ_tot(B)          = ρ_g(B) + ρ_R(B) + S_B/E_spl(B)                    ← INCLUDES the measured spliced
    ρ_tot,unspl(B)    = ρ_g(B) + ρ_R(B)                                   ← EXCLUDES it
    ρ_tot(I)          = M_I·[ f_g(I)/E_g(I) + (1−f_g(I))/E_r(I) ]
```

### 3.3 The two enrichment ratios — and why each is unbiased

```
    r₁ = ρ_tot(E) / ρ_tot(B)                    = e(E)/e(B)      TOTAL frame  — B includes spliced
    r₂ = ρ_tot,unspl(B) / ρ_tot(I)              = e(B)/e(I)      UNSPLICED frame — spliced excluded
```

Each isolates the efficiency ratio **only if the two sides share their underlying content**, which is exactly
why the frames differ:

* `r₁`: the exon holds all its RNA; the boundary's **total** holds the same RNA across both routes plus the
  same gDNA. Matched `{g, R}` ⇒ the content ratio is 1 ⇒ `r₁ = e(E)/e(B)`.
* `r₂`: both sides are `{g, R_cont}`. Including `S_B` here would break the match — the intron has no spliced.

These are precisely **families C and B** of `message_layer_derivation.md` §6.5 — **the taxonomy transfers**
(which edges are component-set-matched, and in which frame: TOTAL for C, unspliced-only for B).

> ### ⚠ CORRECTION (adversarial review, 2026-07-23). An earlier draft said "the estimator is already
> validated on all 32 scenarios." **That was false on two counts**, both verified in-repo:
>
> 1. **Wrong estimator.** §6.5's census computes `ρ̂` from `node_global_geometry`, which returns
>    `eff = eff_gdna_*` — i.e. `M/E_g`, a **gDNA-frame, belief-free, composition-free** density. §3.2's
>    `ρ_tot = M·[f_g/E_g + (1−f_g)/E_r]` is a **composition-blended** quantity taken from a *solve*. They
>    coincide only at `f_g = 1`. The validation does not transfer to the estimator this document proposes.
> 2. **Wrong scenario set.** `enrichment_ratio_census.py` filters `"0.50" in d.name` — **all 20 conditions are
>    UNSTRANDED**. The 12 `ss_0.99` scenarios were never censused, and stranded is precisely the arm that
>    regressed in both prior attempts (§13.6b +0.0047, §13.6h +0.0044).
>
> What IS measured: the component-set taxonomy and the effective-length frames, on 20 unstranded scenarios
> (§6.5 E1/E2 — capture-off medians **+0.008** (C) and **−0.022** (B) vs an unmeasurable family D at
> **−1.116**).
>
> ### ✅ GAP CLOSED — E0 executed 2026-07-23. Gate PASSED.
> `enrichment_ratio_census.py` rewritten for **all 32 conditions** with **both** estimators. Capture-off
> medians under the **blended** `ρ_tot`:
>
> | | A | B | C | D *(control)* |
> |---|---|---|---|---|
> | unstranded | −0.090 | −0.122 | −0.191 | **−2.684** |
> | **STRANDED** | −0.103 | −0.078 | −0.174 | **−2.418** |
>
> All six measurable cells centre on ~0; the unmeasurable control is far from 0 in both arms. **Stranded is
> the cleaner arm** (sd 0.236/0.403/1.514). Family A is bit-identical between estimators, as it must be
> (`f_g ≡ 1` there). ⚠ The blend is slightly worse-centred than `M/E_g` on B and C — it inherits the solve's
> composition error — so **keep `M/E_g` as a fallback and A/B the two when the solver is wired.**

### 3.4 Scale both messages into the boundary's frame, then fuse

```
    exon → B :   ρ_c(B) ← ρ_c(E) / r₁
    intron → B:  ρ_c(B) ← ρ_c(I) · r₂
```

Both messages now live in `B`'s frame. Fuse them with `B`'s current belief by the existing
precision-weighted rule (§12.4), then solve. **All of this happens at solve time**, when both messages, the
spliced density and the unspliced masses are simultaneously in hand.

---

## 4. What this fixes, stated against measurement

The corrected diagnosis (`message_layer_derivation.md` §13.6k) is that the damage is **not** in the
subtraction (condition number 1.1 median) but in **transport across a scale change**: an error at the boundary
lands in an intron whose RNA density is ~86× smaller, amplifying ×6.7. Scaling by a *measured* ratio is
exactly the correction for a scale change. The paired diff (§13.6e) says where it must show up:
`intron adj-junction` classes carry **46 % of the damage**.

---

## 5. Known risks — to be attacked, not assumed away

| # | risk |
|---|---|
| **K1** | **Circularity residue.** `r₂`'s numerator is built from `k(I)`, so `r₂` is not a fully independent comparison of `B` and `I`. Its independent content is `M_B` (observed). Does `r₂` degenerate toward 1 by construction? Needs checking algebraically. |
| **K2** | **Uncertainty propagation.** `r₁`, `r₂` are ratios of *solved* quantities. If applied as exact, this repeats §13.6h — an uncertain correction treated as certain. Each ratio's variance must enter the message precision. |
| **K3** | **`k(I)` when the intron is unsolved.** Precision 0 ⇒ no `k` ⇒ no `ρ_tot(B)` ⇒ no ratios. Fallback must be honest (emit nothing), not a default. |
| **K4** | **`f_g(E)` is itself uncertain**, so `ρ_tot(E)` is too, and it sits in `r₁`'s numerator. |
| **K5** | **Component-agnostic capture is false** — the simulator carries `gdna_split_penalty = 0.2`. §6.2 argues a *constant* `γ_c` cancels in a ratio; a junction-spanning probe makes it position-dependent (§6.3, deferred). |
| **K6** | **`S_B` depletion** (trap #3, ~300×) corrupts `ρ_tot(B)` and hence `r₁` — the same sensitivity as the xfailed depletion test. |
| **K7** | **Non-junction boundaries** have `S_B = 0`; then `r₁` and `r₂` compare the same quantity and the scheme must degrade gracefully to §12's plain relay. |
| **K8** | **Terminal/degenerate geometry:** `E` or `I` absent (a terminal), `M_B = 0`, `E_spl → 0`. |

---

## 5b. K1 RESOLVED — and it exposes an asymmetry between the two ratios

Substituting §3.1 into `r₂` and simplifying (verified numerically to 5.8e−16 over 20,000 random draws):

```
    r₂  =  (M_B / M_I) · ( k·E_g(I) + E_r(I) ) / ( k·E_g(B) + E_r(B) )
```

**The `(k+1)` factor cancels exactly.** `r₂` is therefore *not* degenerate — it is dominated by the two
**observed masses**, with composition entering only through a weak effective-length blend. K1 is closed.

How weak? Across the entire plausible range `k ∈ [0.01, 100]` the blend factor moves 17.9 → 11.4, i.e. a
**1.57× swing end-to-end** (not the order-of-magnitude I first assumed). So `r₂` tolerates a fairly poor `k(I)`
and remains useful — good news for the design, though K2's variance propagation still stands.

**But `r₁` does not share this property.** Because the spliced term is *added*,

```
    r₁ = ρ_tot(E) / [ ρ_tot,unspl(B) + S_B/E_spl(B) ]
```

the `(k+1)` does **not** cancel, and `r₁` depends on `k_E` — **the exon's own composition, which is precisely
what the exon→boundary message exists to establish.** That is a genuine circularity, not a residue.

> ### The asymmetry to act on
> **`r₂` (intron → boundary) is robust and should be applied with confidence.**
> **`r₁` (exon → boundary) is self-referential and must be treated as weak** — or derived differently.
>
> This lines up with everything measured independently: the intron→junction edge is the well-conditioned one
> (§6.5 family B, δ median −0.022), and the exon→junction edge has been the problem all along. The design
> should exploit `r₂` and be sceptical of `r₁`, rather than treating them as a symmetric pair.

**Candidate escape from `r₁`'s circularity (not yet derived):** once `r₂` has placed the intron's message in
`B`'s frame, `B` can be solved from *that* plus its own strand likelihood and spliced channel — and only then
is `ρ_tot(B)` known well enough to *measure* `r₁` and use it in the reverse direction (B → exon). That makes
the junction solve **sequential** (intron first, exon second) rather than simultaneous, which is a materially
different and possibly better plan than §3.4's joint scaling.

## 6. Implementation plan

| step | change | gate |
|---|---|---|
| **J1** | `junction_frame.py` — pure functions: `k_from_belief`, `boundary_composition_from_k`, `total_density`, `enrichment_ratios`. No solver coupling; unit-tested against closed forms. | new tests only; production byte-identical |
| **J2** | Wire `r₁`/`r₂` into `_scan` at junction boundaries, **behind a flag**, applying the scaling to both incoming messages. Retain the existing path when `S_B = 0` or either flank is unsolved (K3/K7). | `gdna_none` delta; suites |
| **J3** | Propagate each ratio's variance into the message precision (K2/K4). | paired diff: `intron adj-junction` damage share must fall |
| **J4** | Full battery per condition, stranded **and** unstranded; then decide the default. | 32-scenario A/B |

**Standing gates:** `gdna_none` as a delta (current default 3,704,635); `pytest tests/calibration/ tests/native/`
(326 passed, 1 xfailed, 2 xpassed); paired diff on stranded; goldens last. **No new constants.**
