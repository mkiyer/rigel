Yes, the sequential option is reasonable. It avoids the challenge of splitting the exon into spliced RNA + unspliced RNA + gDNA.

However, I don't think it will be the most accurate, and here's why:
- under hybrid capture, exons will be highly enriched. The signal in an exon will inherently be much more precise because of the high counts.
- the intron will be depleted. the discrete counts in a deplete region will be sparse and prone to sampling inaccuracy in a depleted low count region

Therefore:
- the intron message will be imprecise
- the exon message will be more precise

Consider "very strong" hybrid capture, where intron basically has zero counts, and exon has high counts. We will get no information from the intron.

Then, how would we compute the solution?

The exon message should still give us a message with gDNA + RNA. If it is highly enriched, we should get a gDNA estimate from that message.

Then we just have ONE message (intron message is empty).

Boundary is spliced (known fixed rna) + unspliced (unknown rna + gdna)
Exon message has RNA+gDNA belief.

- We CAN compute the exon message total density, because it has composition
- Computing the boundary total density requires knowing the composition of the boundary.

So we are STUCK and cannot solve. We would have to make an assumption to solve!

We could assume the boundary unspliced fragments are 100% gDNA. That would actually be a very reasonable initial assumption (remember, this is a pass-0 solve.. it's not FINAL and we get to re-solve after we fit our gDNA hyperprior).

Typically with real data we have almost no nascent RNA.

So I've talked myself into agreeing with you.

Step 1: Solve unspliced using intron message. If intron message is ZERO, the boundary uses its 100% gDNA default assumption (at 0 precision, because the current belief has no measured support).

This assumption + the intron (when it actually has data) allows a step-wise solve of the unspliced density first, then the boundary's TOTAL density.

Then we have total densities! We can compute enrichment ratios. Then we can solve properly.

So.. I am fairly confident in this plan.

The key is that a 100% gDNA assumption is the correct assumption for the unspliced fragments at an exon-intron boundary.

We have to be CAREFUL! 

This does not always hold true for all boundaries... when we jump into complex transcripts, we may need to reconsider or elaborate on this. We can have complex overlapping transcripts on opposite strands giving rise to challenging scenarios.

This is the starting point and the straightforward case. I believe it GENERALIZES TOO!!!

Generalization:
- we need total densities of left message, node being solved, and right message
- once we have total densities, we can compute enrichment ratios (left / node), (node / right)
- we use enrichment ratios to scale densities (left densities x left enrichment ratio), (right densities x right enrichment ratio)
- now we have everything in the NODE's frame of reference
- then solve

This procedure should be general.

The question is -- HOW do we compute this for every situation. Many are trivial, which is why we didn't derive this formally yet. Now we are getting to more complex situations where we need the formal derivation.




---
---

# FORMALIZATION (2026-07-23)

*The owner's design above, made precise. Companions: `junction_enrichment_scaling.md` (the junction
instance), `message_layer_derivation.md` (§6.5 the δ taxonomy, §12 the relay + sink model, §13 the corrected
diagnosis).*

## 1. Why TOTAL DENSITY is the right pivot

Enrichment is component-agnostic at pass-0, so it multiplies *all* components equally — it is visible in the
**total, and only in the total**. Composition *ratios* are enrichment-free (`k = ρ_g/ρ_R` cancels `e`); the
total is where `e` lives. Hence the procedure keys on totals.

The obstruction, exactly as stated above:

```
    M_x = ρ_g(x)·E_g(x) + ρ_R(x)·E_r(x)
    ρ_tot(x) = ρ_g(x) + ρ_R(x) = M_x·[ f_g(x)/E_g(x) + (1 − f_g(x))/E_r(x) ]
```

## 2. ⭐ THE BOUNDING LEMMA — the composition assumption costs almost nothing

The bracket is a **convex combination** of `1/E_g` and `1/E_r`, so **for any composition whatsoever**

```
    M_x / max(E_g, E_r)   ≤   ρ_tot(x)   ≤   M_x / min(E_g, E_r)
```

The entire composition uncertainty is bounded by the effective-length ratio. Measured on the real FL models
(gDNA ~N(300,60), RNA ~N(200,50)):

| frame | `E_g` | `E_r` | worst-case factor |
|---|---|---|---|
| REGION contained, L = 1000 | 701.0 | 801.0 | **1.143** |
| REGION contained, L = 3000 | 2701.0 | 2801.0 | **1.037** |
| BOUNDARY crossing, L ≥ 1000 | 300.0 | 200.0 | **1.500** |
| REGION contained, L = 300 (short) | 24.4 | 101.4 | 4.149 ⚠ |

> **A totally wrong composition still pins the total density to within ~1.04–1.5×.** Against enrichment cliffs
> of 10²–10³×, that is negligible. **This is what makes the whole framework work**, and it is the formal
> justification for the "100 % gDNA at zero precision" fallback: the assumption is not load-bearing for the
> *ratio*, only for a second-order correction to it.

Independently corroborated: `junction_enrichment_scaling.md` §5b measured `r₂`'s composition sensitivity as a
**1.57× end-to-end swing** across `k ∈ [0.01, 100]` — the same bound by a different route.

**The exception is SHORT regions** (`L ≲ fl_mean`), where `E_g` collapses faster than `E_r` and the factor
exceeds 4×. Those are the same nodes §12.2 flags as structurally data-free. **They must be excluded as
enrichment references** — a structural condition, not a tuned threshold.

## 3. Obtaining the node's own total density, case by case

| node being solved | recipe |
|---|---|
| **intergenic region** | `f_g ≡ 1` certain ⇒ `ρ_tot = M/E_g`, exact. The anchor. |
| **TSS/TES seam** | same — no RNA crosses. |
| **intron / exon region** | from its own solved `f_g`; §2 bounds the damage if that solve is poor. |
| **splice-junction boundary** | `ρ_tot = ρ_unspl + S_B/E_spl`; the spliced part is **measured**, the unspliced needs §4. |
| **non-junction boundary** | `S_B = 0` ⇒ the unspliced case alone. |

## 4. The junction boundary — the step-wise solve

**Step 1.** The boundary's unspliced crossing and the **intron** share `{g, R_cont}`, so their density ratio
`k` transfers (enrichment-free), and the boundary's own mass identity gives its composition **in its own
frame**:

```
    ρ_R(B) = M_B / ( k(I)·E_g(B) + E_r(B) ) ,   ρ_g(B) = k(I)·ρ_R(B)
```

`f_g` itself is **not** transferable — only `k` is — because `f_g` depends on the node's own `E_g, E_r`
(crossing vs contained). Copying `f_g` across a frame change is a known past defect.

**Step 2 — the silent-intron fallback.** Under very strong capture the intron holds zero counts and there is
no `k(I)`. Assume the unspliced crossing is **100 % gDNA at ZERO precision**. Justified because (a) real
libraries carry little contiguous unspliced RNA, (b) §2 bounds the cost at **1.5×**, and (c) zero precision is
the honest part — it resolves the *frame* so a ratio can be formed, without asserting a belief. Pass-2's gDNA
hyperprior re-solves.

**Step 3.** `ρ_tot(B) = ρ_g(B) + ρ_R(B) + S_B/E_spl(B)`, then ratios → scale → fuse → solve.

**Why step-wise, not simultaneous.** `r₂` (intron↔boundary) is robust — the `(k+1)` cancels, leaving observed
masses times a weak blend. `r₁` (exon↔boundary) is **circular**: it depends on the exon's own composition, the
very quantity the exon→boundary message exists to establish. Resolving `B` from the intron (or the fallback)
**first** breaks that circularity, after which `r₁` becomes a *measurement*.

## 5. ⚠ Where the 100 % gDNA default is WRONG — enforce structurally, do not assume

| situation | why it fails | detectable by |
|---|---|---|
| **opposite-strand overlap (AMBIG)** | the crossing carries another transcript's RNA on the other strand | `free_pos & free_neg` |
| **exon↔exon boundary** | contiguous mature crosses; the unspliced pool is RNA-rich | `mrna_active_s` |
| **retained-intron / nascent-rich loci** | the contiguous RNA is not scarce | the intron's own solve |
| **boundary with no junction** | no routing split; RNA simply continues | `S_B = 0` |

In each, the fallback must be **refused** (emit nothing) rather than applied. §2's bound does not cover a
component set that is qualitatively different — a wrong frame is worse than no message.

## 6. Open list — the per-case recipes still to derive

1. **exon↔intron junction** — §4. Ready.
2. **AMBIG boundaries** — which flank supplies `k`, and per strand. **Not covered by §4.**
3. **exon↔exon boundary that is also a junction** (§12.7 topology) — `mrna_active_s` and `S_B > 0` together.
4. **TSS/TES seam → exon** — the one edge where `δ` is not measurable at all (§6.5 family D). It should stay
   the honest exception rather than acquire a recipe.
5. **short regions** — excluded as references per §2; confirm excluded, not defaulted.

## 7. Implementation shape

One pure module, `enrichment_frame.py`, no solver coupling, closed-form unit tests:

```python
    total_density(mass, f_g, E_g, E_r)            -> rho_tot        # §1
    k_from_belief(f_g, E_g, E_r)                  -> rho_g/rho_R    # enrichment-free
    boundary_unspliced_from_k(k, mass, E_g, E_r)  -> (rho_g, rho_R) # §4 step 1
    enrichment_ratio(rho_tot_src, rho_tot_dst)    -> r (+ variance)
    is_valid_reference(node)                      -> bool           # §2 short-region + §5 exclusions
```

The solver then, at solve time only: resolve the node's own total → form `r_L`, `r_R` → scale → fuse → solve.
**Gates unchanged:** `gdna_none` as a delta (current 3,704,635), the 32-scenario battery per condition
stranded *and* unstranded, the paired diff on `intron adj-junction`, goldens last. **No new constants.**







