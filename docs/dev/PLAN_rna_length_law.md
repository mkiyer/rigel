# PLAN — the RNA fragment-length law: the sj-observability bias, and the fix already in the tree

    ⚠ **A DEV DOC and a SANDBOX.** Nothing here is authoritative and nothing may cite into it.
    ⛔ **NOTHING HERE IS BUILT.** Every number was measured 2026-08-31 by prototypes outside `src/`;
    the state of what IS built is `ROADMAP.md` §0.

## 1. The claim, in one paragraph

`rna_pmf` is fitted from the spliced pool and is **−4 bp too short as a law of mature RNA**, and a
further **−0.3 to −1.7 bp** short as a law of *all* RNA. It reaches the deliverable through
`effective_lengths_em`, where injecting the true law is worth an attributable **−18,208 transcripts** on
the ladder's `g05` capture-ON row. ⭐ **The fix needs no new input and no new model**: the two-pool
contrast already solves for the contaminant law and throws it away, and that discarded law is unspliced
RNA measured on pools with no observability selection — against truth it lands **inside 0.2 bp**.

## 2. Where it hurts — the consumer split, measured

Injecting the simulator's true RNA pmf into one consumer at a time, ladder `g05` capture-ON, transcript
`count_abs_err`:

| the true RNA law into… | Δ transcripts | Δ fp_mass |
|---|---|---|
| the EM **scorer** (`rna_fl`) | +156 | −46 |
| the second-pass **drain** (`score_held_fragments`) | −822 | −465 |
| **GEOMETRY** (`effective_lengths_em` → the prior) | **−18,208** | **−7,972** |

⭐⭐ **So RNA is NOT the gDNA problem repeated.** gDNA needed two estimands because its two consumers
wanted different populations; RNA's scorer is indifferent and essentially all of the value is in geometry.
**This is a plain bias fix at one consumer, not a routing question** — which is why it was measured rather
than assumed symmetric, and the measurement is what stopped a second realized-law field being built for
nothing.

## 3. The mechanism, decomposed by stage

Ladder, mean fragment length in bp:

| stage | capture-OFF `g50` | capture-ON `g05` |
|---|---|---|
| raw spliced pool | 218.35 | 232.56 |
| after the sj de-tilt | 204.08 (**−14.27**) | 221.52 (**−11.03**) |
| after EB shrinkage | 204.09 (+0.01) | 221.53 (+0.00) |
| **MATURE truth** | **208.18** | **225.94** |
| residual | **−4.10** | **−4.41** |

⭐ **The de-tilt OVER-corrects.** The raw pool is ~+10 bp long because longer fragments cross more
junctions — that part is real and the de-tilt is right to remove it — but it removes ~14 and overshoots by
~4. The EB shrinkage is exonerated (0.01 bp; the pools are 1.5–3.8 M against an ESS of 1,000).

⛔ **WHY, and it is written in `fl.py` already.** `sj_opportunity` is `pi(w)` — the probability a
uniformly placed length-`w` fragment **CROSSES** an annotated sj. But the pool is selected on the splice
being **OBSERVED**: `sj_implicit` fragments (a junction that fell in the unsequenced mate gap) are excluded
by the accumulator. Observability **FALLS** with length — the junction must land inside a sequenced read,
and the read pair covers a shrinking fraction of a growing fragment. So the true selection is
`pi(w) x P(observed | crosses, w)` with the second factor decreasing, the net selection is weaker than
`pi(w)` alone, and dividing by `pi(w)` alone removes too much. ⭐ The sign and the stability across all four
conditions (−3.97 … −4.43) are both what that predicts.

⚠ **A second, smaller term — the same ESTIMAND mismatch as the gDNA frame error, a third instance.**
`rna_pmf` estimates the **MATURE** law (the spliced pool is certified mature), while `region_geometry` uses
it for the effective lengths of **ALL** RNA — `fl.py`'s own docstring says "nascent unspliced + spliced".
Nascent runs longer (216.5 vs 208.2 off capture), so the mismatch is worth **1.7 bp off capture** (nascent
20.8 % of RNA) and **0.3 bp under it** (2.7 %).

## 4. ⛔ REFUTED — do not retry

**The sj bank's model-free identity.** `EQUATIONS.md` §3c gives `count/inv_length_sum + 1 = E[w]` with no
pmf and no divisor, and `payload.sj_count` / `sj_inv_length_sum` look like exactly the banks for it.
Measured, it is **equally biased**: **−4.62** and **−5.34** against mature truth, versus the de-tilt's
−4.10 and −4.41. ⭐ **The reason is structural and worth carrying**: the sj bank counts the same
OBSERVED splices, so it samples the same selected population. **No estimator built on spliced fragments
escapes this bias, however it is derived.** The escape has to be a different population.

## 5. ⭐⭐⭐ The fix: keep what the contrast already discards

`_deconvolved_gdna_counts` solves the two-component system and uses `g`; the complement-first solve
computes `r_hat` — the **contaminant** law — and drops it. That law is unspliced RNA (nascent in introns,
unannotated transcription in intergenic space) measured on CONTAINED pools, which have **no sj-observability
selection at all**. Measured against truth:

| condition | `r_hat` | nascent truth | error |
|---|---|---|---|
| ladder `g50` ss0.99 capture-OFF | **216.45** | 216.48 | **−0.03** |
| ladder `g05` ss0.99 capture-OFF | **216.83** | 216.64 | **+0.19** |

⭐ Inside 0.2 bp, with no new input, no read-length model, and no new bank. The machinery is already in
`fl.py`; the change is to return `r_hat` instead of discarding it.

**The proposed estimator.** Two laws for one quantity, with complementary blind spots — exactly the shape
`_couple_estimands` already implements:

* `r_spliced` — the de-tilted spliced pool. Certified mature, well measured (millions of fragments),
  biased ~−4 bp by observability.
* `r_unspliced` — the contrast's `r_hat`. Unbiased by observability, noisier, and it measures the
  *unspliced* population (nascent-dominated), not the mature one.

⛔ **They are NOT two estimates of one law, and that must not be fudged.** Mature and nascent genuinely
differ (208.2 vs 216.5). What geometry wants is their **mixture at the library's own composition**, and
the composition is something calibration already estimates. So the honest combination is a *mixture*, not
a shrinkage:

    r_all(L) = (1 - phi_nascent) * r_spliced_corrected(L) + phi_nascent * r_hat(L)

with `phi_nascent` the nascent share of RNA. ⚠ **And `r_spliced_corrected` is the open piece**: the −4 bp
observability bias is still there and the mixture does not remove it. Two candidates, neither derived:

1. **Model observability.** `P(observed | crosses, w)` is a geometry problem given the read length —
   the junction must fall in one of the two reads. The read length is knowable from the BAM. ⭐ Principled
   and complete. ⛔ A new input to `fl.py` and a new model; it is the honest but heavier route.
2. **Calibrate the overshoot against `r_hat`.** Where both laws are well measured, the difference between
   `r_spliced` and `r_hat` decomposes into (mature − nascent) plus the observability bias. ⛔ Those two are
   not separable from one comparison, so this needs a second contrast to break the degeneracy and may not
   be identifiable at all. **Check identifiability before building.**

⚠ **The mixture alone (fixing only the estimand mismatch, leaving the −4 bp) is worth roughly 1.7 bp off
capture and 0.3 bp under it** — a fraction of the −18,208 ceiling. It is the cheap half; do not present it
as the whole fix.

## 6. Sequencing, and how to judge it

0. ⭐ **Re-measure the consumer split first.** The three-way numbers above predate the gDNA routing and
   coupling that landed after them, and geometry's sensitivity to the pmf tail is documented as knife-edge.
   ⛔ Do not build against a stale −18,208.
1. **Return `r_hat` from the contrast** (a few lines) and gate it against truth on both fl-gap arms and the
   ladder — it is a new estimator and needs its own falsification, including the zero controls.
2. **Decide the observability route** (§5 candidates 1 vs 2). ⛔ Check candidate 2's identifiability on
   paper before writing any of it.
3. **A/B on `quant_accuracy.py` transcript `count_abs_err`**, against the −18,208 ceiling re-derived at
   step 0, on both sign arms plus the equal-length control, capture ON and OFF. ⛔ NOT on `gdna_frac_est`:
   that metric is what made this session close the gDNA shape question wrongly, and the fl model's consumer
   is the EM.
4. ⛔ **Watch geometry.** This changes an `effective_lengths_em` input, and a ~3 % tail change there was
   measured at ±188k/−15.5k transcripts. Every claim needs the reseed floor beside it.

## 7. Open issues

* ⛔ **The observability correction is genuinely underived** (§5). Everything else here is measured; this
  is the piece that could turn out to need the read-length model, and it carries most of the −4 bp.
* ⚠ **`r_hat`'s accuracy is measured at TWO capture-OFF conditions only.** Under capture the contrast
  degenerates to the intergenic pool alone (`a_0` clips to 1), and `r_hat` is then whatever that
  degeneracy leaves — **unmeasured, and probably not usable**. The nascent share under capture is 2.7 %,
  so the mixture term is small there; that is luck, not design.
* ⚠ **`phi_nascent` is a composition the solver estimates**, so using it here is a coupling from the length
  model to the solve. It is weak (a 20 % share error moves the mixture ~1.7 bp × 0.2), but it is a new
  dependency direction and `fl.py` is layer 2 while the composition lives above it. **Check `_layers.py`
  before wiring anything.**
* ⚠ The **drain** (−822) is a third consumer with its own reading of `rna_pmf`; it is small but it is not
  zero and it is not covered by a geometry-only fix.
