# Message propagation across the density cliff — composition imputation + geometric mature removal

**Status:** design note (2026-07-20). Derives the correct region↔boundary↔region message model for the
exon→boundary→intron chain, resolving two long-standing problems: (1) the **capture-enrichment density cliff**
(adjacent nodes differ in total density by ~10²–10³×) and (2) the **mature-RNA absorption** (the spliced siphon).
Grounded in worked examples and validated on the cached `ambig_dense_10mb` scenarios. **This note specifies the
target model; it also states the hard dependency on the phantom-confidence fix (§7) — the model must NOT be
shipped before that fix, or it makes the aggregate worse.**

Companions: [`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) §0–1 (count-zero-information),
[`message_precision_derivation.md`](message_precision_derivation.md) (the phantom + the τ evidence compiler),
[`message_absorption_fix.md`](message_absorption_fix.md) (the prior absorption diagnosis this supersedes on the
arithmetic).

---

## 1. The problem

Calibration deconvolves each node's **unspliced** fragment mass into `(f_pos, f_neg, f_g)`. Nodes form a bipartite
chain `region ↔ boundary ↔ region`. Messages propagate a source node's solution to a destination node so a
higher-precision node can inform a lower-precision one.

The canonical hard case is a **splice-junction chain**: `EXON (region) — BOUNDARY (splice junction) — INTRON
(region)`. Physically, three fragment populations coexist, each with different propagation behaviour:

| component | where it lives | crosses the junction? | genomic smoothness |
|---|---|---|---|
| **gDNA** | everywhere | yes (unspliced) | smooth (modulo capture) |
| **nascent RNA** (pre-mRNA) | whole gene body (exon+intron), contiguous | yes (unspliced) | smooth (modulo capture) |
| **mature RNA** | exons only | **no** — it *splices* across (counted as spliced reads) | exon-local; **absent in introns** |

The calibration solves the **unspliced** composition `f_g = gDNA / (gDNA + nascent + mature_unspliced)`. The
mature/nascent split is *not* calibration's job (it is the per-locus EM's, and is only identifiable via intronic
RNA — [`nascent_rna_identifiability_intron_required`]). But the mature **does** matter to the message model,
because it is present in the exon and absent at the boundary/intron, so it changes the composition along the chain.

Two distinct discontinuities are superimposed between an exon and its boundary/intron, and **they must be handled
separately**:

1. **The capture cliff (density).** Hybrid capture enriches on-target (exon) fragments and depletes off-target
   (intron) fragments. Capture is **nucleic-acid-agnostic** — it enriches gDNA and RNA by the *same* factor at a
   given position. So it scales all components equally: the **total density jumps** across the cliff, but the
   **composition (`f_g`) is invariant**.
2. **The mature discontinuity (composition).** The exon carries mature RNA that the intron does not (it spliced
   out). This is *not* capture — it is a genuine composition change: the exon's RNA denominator includes mature,
   the boundary/intron's does not.

The current message model conflates these: it imputes the source's **absolute density** and reconciles it against
the destination's total, which is wrong for **both** discontinuities.

---

## 2. The worked example (real numbers — `gdna5_ss_0.50_capture_on`, chain nodes 2931/2930/2929)

`EXON 2931 — BOUNDARY 2930 — INTRON 2929`, an expressed single-strand(+) gene. Oracle decomposition of the
**unspliced** mass (from the accumulator; `scripts/debug/worked_example.py`):

| node | unspliced total | gDNA | mature-unspl | nascent-unspl | spliced (mature) | true f_g |
|---|---|---|---|---|---|---|
| EXON 2931 | 119,021 | 493 | 73,639 | 44,889 | 0 | **0.0041** |
| BOUNDARY 2930 | 962 | 7 | 0 | 955 | 7,502 | **0.0073** |
| INTRON 2929 | 78 | 0 | 0 | 78 | 0 | **0.0000** |

Effective lengths: exon `eff_r ≈ 1074`; boundary crossing `eff_r ≈ 100`, `eff_spl ≈ 100`.

**The capture cliff is ~970×.** Nascent *intrinsic* density (count/bp): intron `78/1806 = 0.043`, exon
`44,889/1074 = 41.8`. The exon is capture-enriched ≈ **970×** for nascent (and, agnostically, for gDNA). Yet the
composition `f_g` is ≈ invariant along the chain (0.004 → 0.007 → 0) — exactly as a capture-agnostic cliff
predicts.

**The mature discontinuity.** The exon's unspliced RNA is 73,639 mature + 44,889 nascent (62% mature). The
boundary's unspliced RNA is 955 nascent, **0 mature** (the mature left as the 7,502 spliced reads). So the exon's
RNA composition is *not* the boundary's — the mature must be removed.

---

## 3. Why the current message model fails on this example

The shipped gDNA message mode (`bp_solver._scan`) is, for source `s` → destination `d`:

```
   mode = log( ρ_g^s / ρ_total^d ) ,     ρ_g^s = f_g^s · M_s / E_g^s ,   ρ_total^d = M_d / E^d
```

i.e. it imputes the source's **absolute gDNA density** and divides by the **destination's total density** to get
the destination's `f_g`. With *correct* beliefs on the worked example:

```
   ρ_g^exon = 0.0041 · 119021 / 1085.6 = 0.45          (exon gDNA density)
   ρ_total^bnd = 556.6 / 100 = 5.57                    (boundary total density, off the cliff)
   ⇒ imputed f_g^bnd = 0.45 / 5.57 = 0.081             vs truth 0.0073  →  11× over-call
```

The over-call is **purely the cliff**: an enriched-frame numerator (0.45) over a depleted-frame denominator
(5.57). σ²_transfer damps the *precision* of this message but does nothing to the biased *mode*. The RNA channel
has the mirror problem, and it *also* over-states the boundary's RNA because it carries the exon's mature.

---

## 4. The derivation

### 4.1 The capture cliff → impute COMPOSITION, not density

Let component density `ρ_c(x) = e(x) · d_c(x)` where `e(x)` is the position-dependent capture enrichment
(agnostic across components) and `d_c(x)` the intrinsic (pre-capture) density. The composition is

```
   f_c(x) = ρ_c(x) / Σ_k ρ_k(x) = e(x)·d_c(x) / (e(x)·Σ_k d_k(x)) = d_c(x) / Σ_k d_k(x)
```

— **independent of `e(x)`**. The composition is capture-invariant; the density is not. So the message should
impute the source **composition** (a ratio of counts) and never the absolute density. In count space the
enrichment ratio `e(s)/e(d)` **cancels identically** — there is no ratio to estimate, no `ρ_total^d` to divide by.
This is the count-space, scale-invariant imputation the "enrichment ratio" intuition points at, done exactly.

### 4.2 The mature discontinuity → remove the contained mature (count space, geometry)

Composition-imputation alone is still wrong across the exon→boundary transition, because the exon's composition
includes mature and the boundary's does not. `f_g^exon = 493/119021 = 0.0041`, but the boundary needs the
**mature-removed** composition `493/(493+44889) = 0.0109 ≈ 0.0073` (truth). So we must subtract the exon's
mature from its RNA before imputing.

The exon cannot separate its own mature from nascent (both are unspliced exon-body fragments). But the boundary
**directly measures** the spliced count `S` (the mature flux through the junction), and the exon's **contained**
mature is geometrically tied to `S`:

```
   mature_contained(exon)  ≈  S · ( E_r^exon / E_spl^boundary )
```

where `E_r^exon` is the exon's RNA contained eff-length (`region_eff_length`) and `E_spl^boundary` is the
boundary's one-sided spliced half-triangle eff-length (`spliced_side_eff_length`). Rationale: both the contained
mature and the spliced count are proportional to the same mature abundance, differing only by the geometric
opportunity to be *contained in the exon body* vs *span the junction*. This is **computed geometry** (from the
FL pmf and the region lengths), **not a tuned constant**.

**Validation (`scripts/debug/mature_geometry.py`, 292 exon–boundary pairs, gdna5):**
`corr(log mature_contained, log [S·E_r/E_spl]) = 0.939`; median predicted/actual = 0.87 (p25–p75 0.72–1.10).
On the worked example: `S·E_r/E_spl = 7502 · 1074/100 = 80,571` vs true contained mature 73,639 — the right pool,
~9% high.

> **Note — this vindicates the current absorption arithmetic.** The shipped code subtracts `absorb = S/E_spl`
> (`SPd/ESPd`) as a *density in the exon frame*: `mature_contained/E_r = S·(E_r/E_spl)/E_r = S/E_spl`. So `S/E_spl`
> is **exactly** the exon-frame contained-mature density — `E_spl` (the half-triangle) is the **correct** divisor.
> A prior diagnosis (this note's author, and [`message_absorption_fix.md`] §5) proposed replacing `E_spl` with the
> crossing eff-length `boundary_side_eff_length`; **that is wrong** — it would compute the ~4× smaller *flux*
> mature. The absorption *arithmetic* is right; its problem is that it is applied inside a density-imputation mode
> that then saturates on the cliff (§3), and that it is currently mixed into the destination's own `f_g` solve
> rather than the onward relay (§5).

### 4.3 The combined message

Impute the destination's `f_g` as the source's **mature-removed composition**, in count space:

```
   g_s     = f_g^s · M_s                                   (source gDNA count)
   RNA_s   = f_r^s · M_s                                   (source total RNA count: mature + nascent)
   mat_s   = S · ( E_r^s / E_spl^bnd )                     (contained mature, from the boundary spliced; ≥ 0)
   nasc_s  = max( RNA_s − mat_s , 0 )                      (the nascent that continues)

   f_g^{s→d}  =  g_s / ( g_s + nasc_s )                    (the imputed destination composition; mode = log f_g)
```

Worked example (correct beliefs): `493 / (493 + max(118528 − 80571, 0)) = 493/38450 = 0.0128 ≈ 0.0073`. ✓
(vs 0.081 from the shipped density mode — an 11× improvement.)

**Degenerate-safe:** `S = 0` (no junction, e.g. region→region) ⇒ `mat_s = 0` ⇒ ordinary composition-imputation;
pure intron source ⇒ `mat_s = 0` ⇒ nascent-only, unchanged. The `max(·,0)` handles over-subtraction on short
flanks without a clamp cliff (a fully-absorbed source relays ≈0 nascent — a weak zero, per §5/§7).

### 4.4 Why NEITHER piece alone works (the two must be combined)

| mode | worked-example f_g^bnd | truth 0.0073 | failure |
|---|---|---|---|
| density-imputation + absorption (shipped) | saturates (≈1) | ✗ | cliff over the depleted denominator |
| composition-imputation, **no** mature removal | 0.0041 | ✗ | mature inflates the RNA denominator → under-call |
| **composition-imputation + geometric mature removal** | **0.0128** | ✓ | — |

---

## 5. Where the message applies — own solve vs onward relay

The mature removal must shape only what **propagates onward as nascent**, not a node's *own* `f_g`:

- **Exon's own `f_g`:** total RNA (mature+nascent) vs gDNA. No mature removal — mature *is* RNA for the exon's
  gDNA-vs-RNA split.
- **Boundary's / intron's own `f_g`:** nascent vs gDNA (their unspliced has no mature). The message *into* them
  is the source's mature-removed composition (§4.3).
- **The relay onward:** the mature-removed (nascent) composition continues; the mature does not. This is the
  §4 formula applied on the boundary's **outgoing** side.

This is the resolution of the empirical finding that turning the absorption *fully off* currently improves
boundary accuracy: the shipped absorption injects the mature removal into the destination's own solve **and** in
a density mode that saturates, so it over-calls gDNA. In the composition frame, applied to the relay, it is
correct.

---

## 6. What changes in `bp_solver._scan`

1. **Mode → composition.** Replace the gDNA `log(ρ_g^s / ρ_total^d)` and the RNA-total `log(ρ_r^s / ρ_total^d)`
   modes with the count-space **composition** modes (§4.3). Per-strand modes analogously on `f_pos/f_neg`.
2. **Keep the absorption arithmetic** (`S/E_spl` = contained-mature density; §4.2). Do **not** change `E_spl` to
   the crossing eff-length.
3. **Move the mature removal to the relay** (the nascent that continues), not the destination's own `f_g`
   solve (§5).
4. **Precision unchanged in spirit** — the message still carries its honest precision (count + composition +
   σ²_transfer + the spliced measurement precision). The mode change is orthogonal to precision (no backdoor;
   composition-imputation is an ordinary message, still transfer-variance-damped).

Config/flags: this replaces the `_RNA_ABSORB` and the `_IMPUTE_MODE` experiment with a single production path;
no new magic numbers (`E_r/E_spl` is existing geometry).

---

## 7. THE HARD DEPENDENCY — the phantom must be fixed first

Composition-imputation is a **more faithful** copy of the source belief than density-imputation (it removes the
cliff's accidental dilution). That is only good if the **source belief is honest**. On unstranded chains it is
**not** — the phantom-confidence cascade ([`message_precision_derivation.md`]; verified 2026-07-20) manufactures
a confident-but-wrong `f_g` at nodes with zero intrinsic information. A/B (production-faithful pass-0, exon
|Δf_g|):

| condition | density mode | composition mode |
|---|---|---|
| unstranded capON gDNA=0 | 0.383 | **0.475** (worse) |
| unstranded capON gDNA=high | 0.230 | **0.187** (better) |

Composition-imputation *helps* where the source belief is meaningful (enriched, gDNA-present) and *hurts* where
it is phantom (zero/low gDNA), because it propagates the phantom faithfully. **Therefore this model must land
AFTER the phantom-confidence fix (τ-precision, `message_precision_derivation.md`), not before.** With honest
source beliefs, composition-imputation + geometric mature removal is the correct model; without them it is a
higher-fidelity copier of garbage.

**Sequencing:** (1) fix the phantom (message precision from reference-free evidence, not belief variance);
(2) land this composition + mature-removal message model; (3) re-evaluate the gDNA-hyperprior fit on the now-clean
pass-0 substrate.

---

## 8. Open items

- **Mature-estimator calibration.** `S·E_r/E_spl` predicts contained mature at corr 0.94 but ~10% high on the
  worked example; check whether the small bias needs a second geometric term (both flanking boundaries contribute
  spliced; the two-junction bookkeeping) — geometry, still no tuned constant.
- **Per-strand vs λ-axis.** §4.3 is written on the λ-axis (`f_g`); confirm the θ-axis (strand split) composition
  imputation is the mirror and that AMBIG nodes behave.
- **Real-data cliff magnitudes.** The sim's ~970× enrichment and abundant nascent are extreme; the sign of the
  fixes is architectural, but the final validation is on real caches (LBX0190, MO_3005).
- **Backdoor check (owner's concern), resolved:** composition-imputation is not a transfer-variance bypass — it
  is an ordinary message with the same precision machinery; only the *mode* changed from absolute density to a
  scale-invariant ratio. The mature removal is a **mode** operation (physical accounting of what crosses the
  junction), decoupled from precision.
