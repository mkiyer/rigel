# A robust strand-overdispersion estimator — the problem, and the design

**Owner-directed, 2026-07-28.** Written because the fitted strand overdispersion saturates its `Beta(2,2)`
ceiling on every real library, and because "raise the ceiling" is the wrong answer for a biological reason.
Evidence: `scratchpad/od_{outliers,screen,shape,design}.py`. Nothing here is implemented.

---

## 1. What a "seed" actually is — the concrete version

The strand overdispersion is fitted from **seed nodes**: places the model believes it can see gDNA without
RNA contaminating the view. Two kinds (`gdna_strand.fit_gdna_strand_from_substrate`):

* **contained regions** — `region_count_observable = (signature & EXON_BITS) == 0`, i.e. a region the
  **annotation** marks as intergenic or intron-only, excluding AMBIG;
* **boundary sides** — exon↔intron and exon↔intergenic seams (`strand_deconv.boundary_side_seeds`), needed
  because hybrid capture depletes off-target intergenic/intronic gDNA.

Each seed contributes a triple `(sense, total, gdna_weight)`:

* `total` = its raw unspliced fragment count (pos + neg),
* `sense` = how many of those are on the transcript-sense strand,
* `gdna_weight` = `count_gdna_frac`, "what fraction of this node is gDNA", supplied by the **count module**
  so that the strand fit is not circular (the weight uses no strand information).

The estimator then pools them:

```
    mean_s   = ½·w_s + κ·(1 − w_s)                                  the seed's mixture sense rate
    excess_s = (sense_s − N_s·mean_s)² − N_s·mean_s(1 − mean_s)     observed − Binomial
    scale_s  = n_c(n_c − 1)·¼,      n_c = w_s·N_s                   BetaBinom excess-variance scale
    od       = Σ_s excess_s / Σ_s scale_s                           ← the shipped estimator
```

then shrinks toward a prior by **seed-node count** and clips to `[0, 0.2]`.

**The anchor that makes all of this well-posed: gDNA's strand mean is ½ BY CONSTRUCTION.** It is genomic
DNA; there is no strand preference. We are estimating a *dispersion around a known mean*, which means
contamination and dispersion live in **different moments** — contamination *displaces* the mean, dispersion
*spreads* it. That asymmetry is the entire basis of the design in §3.

## 2. What is wrong — three defects, measured

### D1 — the gDNA weight is the ANNOTATION'S LABEL, not a measurement

For a count-observable region, the count module reads density *directly from the region's own unspliced
counts* — so `count_gdna_frac ≈ 1.0` **by construction** for exactly the regions selected as seeds. It is
not an independent estimate of gDNA-ness; it is a re-encoding of "the annotation called this intergenic".

Measured consequence — the top seeds on `vcap`, by share of the pooled numerator:

| excess_s | share of Σexcess | N | sense | **sense fraction** | **gdna_weight** |
|---|---|---|---|---|---|
| 571,912 | **12.3 %** | 1,523 | 5 | **0.003** | **1.000** |
| 466,141 | 10.0 % | 1,392 | 13 | **0.009** | **1.000** |
| 348,981 | 7.5 % | 1,200 | 9 | **0.007** | **1.000** |

gDNA cannot be 0.3 % sense. These are **transcripts** — in territory the annotation calls intergenic or
intronic — entering a mean-½ fit at **full weight**. That is what "the count-clue fix" means: *a seed's
gDNA weight must be a measurement, not a restatement of the annotation.*

⚠ And it is not only an *unannotated*-transcription problem. On the **stranded synthetic**, where the
simulator emits only annotated transcripts, the same failure occurs via *annotated* nascent RNA in introns
(13.8 % of seeds rejected by a consistency screen). The general statement is **RNA anywhere in
"gDNA-observable" territory**; unannotated transcription is merely the case the annotation cannot warn us
about.

### D2 — ⭐ the prior shrinkage counts SEED NODES, but the information is PAIRS

`od = (n_seed_nodes·od_mom + prior_weight·od₀)/(n_seed_nodes + prior_weight)` with `prior_weight = 30`.

| sample | n_seed_nodes | **prior's actual weight** | seeds with n_c ≥ 2 | Σ pairs `n_c(n_c−1)/2` | **median n_c** |
|---|---|---|---|---|---|
| LBX0190 | 160,366 | **0.019 %** | 107,534 | 699,894 | **2.0** |
| LBX0588 | 297,293 | 0.010 % | 131,105 | 4,914,517 | **1.0** |
| MO_3021 | 624,913 | 0.005 % | 380,713 | 10,745,654 | **2.0** |
| vcap | 863,246 | **0.003 %** | 695,509 | 100,998,511 | 8.0 |

**The median seed carries one or two fragments.** A seed of size 1 carries *zero* information about a
variance; a seed of size 2 carries one pair. Yet 160 k such seeds outvote a `prior_weight` of 30 by four
orders of magnitude, so **the prior never bites and an unidentifiable dispersion is reported as a confident
fit.** This is a **power** failure being reported as a measurement, and it is the single clearest defect.

### D3 — the pooled ratio has unbounded per-seed influence

`scale_s ∝ n_c²`, so the estimate is fragment-**squared**-weighted with no bound on any one seed. A single
`N = 1,523` seed contributes ~12 % of the numerator on its own, and in transcribed territory the largest
seeds are precisely the transcripts. The bulk of the population (n_c = 1–2) contributes signed noise that
nearly cancels — which is why the numerator is a near-cancelling sum whose sign is decided by a handful of
seeds (on LBX0588 the top 0.1 % of seeds is **298 %** of the net numerator).

### ⚠ Two conclusions I reached and then had to RETRACT — do not repeat them

Both were artifacts of tiny seeds, and both looked convincing:

| retracted claim | why it was wrong |
|---|---|
| "60–87 % of seeds are DISPLACED (\|sense frac − ½\| > 0.3), so the contamination is pervasive" | the median such seed has **N = 1–2**. At N = 2 the only outcomes are sf ∈ {0, ½, 1}, so ¾ of *fair coins* land there. That is counting, not contamination. |
| "the MEDIAN of per-seed moments still saturates, so >50 % of seeds are genuinely overdispersed" | at n_c = 2 the per-seed ratio `excess/scale` is **exactly ±1**, so any quantile of that population returns ±1 or 0 and measures seed SIZE, not dispersion. |

What survives, from `od_shape.py`: **median z = 0.000 on every real sample**, and `IQR(z)/1.349 = 1.03–1.48`
— entirely explainable by small-N discreteness. **The seed population is not broadly overdispersed.**

## 3. THE DESIGN

Three changes, each addressing one defect. ⚠ **P1 was REFUTED by the owner and by measurement after this
section was first written — read the block immediately below before reading P1 itself.** P2 stands alone
and is the change that needs no assumption.

### ⛔⛔ P1 AS FIRST DRAFTED IS REFUTED (owner objection, 2026-07-28) — READ THIS BEFORE §3.1

The owner's objection: *"we need to be rather generous, or else we will just filter out our overdispersion.
If we have 10 fragments, 9 sense and 1 antisense, the mean is 90 %. This seems to boil down to a magic
number. What level do we reject at?"* **That is correct, and it is correct more deeply than stated.**

Testing a seed against **Binomial(½)** is circular: overdispersion *is* seeds landing away from ½, so the
screen removes the signal. The obvious repair is to test against the model's own **ceiling**, `Beta(2,2)`
(`od = 0.2`, "the most overdispersion allowed") — a seed inconsistent with the *ceiling* cannot be gDNA
under any admissible dispersion, and that would be an existing constant rather than a new one.

**Measured: the ceiling rejects NOTHING.** Under intraclass correlation `od`, a seed of `n` fragments
carries only `n_eff = n/[1 + (n−1)·od]` independent observations:

| n | od = 0.2 | 0.05 | 0.02 | 0.01 | 0 |
|---|---|---|---|---|---|
| 10 | 3.6 | 6.9 | 8.5 | 9.2 | 10 |
| 100 | 4.8 | 16.8 | 33.6 | 50.3 | 100 |
| **1,523** | **5.0** | 19.8 | 48.4 | 93.9 | 1,523 |

**At the ceiling, a 1,523-fragment seed is worth FIVE coin flips.** Five flips landing 5–0 is ordinary, so
even `vcap`'s worst seed (N = 1,523, sense = 5) has `p = 1.1e-4` against `BetaBinom(n, 2, 2)` and survives
Bonferroni over 160 k seeds. And the owner's 9/10 example has `p = 0.22` — comfortably kept, as it must be.

**And the level is not removable, only relocatable.** At what assumed maximum dispersion does that vcap
seed become rejectable (Bonferroni, `p < 6.2e-6`)?

| assumed max od | 0.200 | **0.100** | 0.050 | 0.020 | 0.010 |
|---|---|---|---|---|---|
| p-value | 1.1e-4 | **5.8e-9** | 6.9e-17 | 1.9e-12 | 1.7e-12 |
| verdict | keep | **REJECT** | REJECT | REJECT | REJECT |

> ⭐ **The rejection level is EXACTLY equivalent to asserting a maximum gDNA dispersion.** There is no
> screen without that assertion. The magic number is not eliminated by any choice of test statistic — it is
> the bound itself, wearing a different hat.

**The identifiability statement, which is the real result:** a Beta-Binomial with mean ½ and free
overdispersion can explain **any** single seed's split. **Contamination and overdispersion are not
separable at the level of one seed.** They differ only in the *population* shape — `Beta(a,a)` vanishes at
`p = 0` and `p = 1`, while transcripts pile up exactly there — which is separable in principle by a mixture
fit, but not on libraries whose **median seed carries 1–2 fragments**.

**Consequences for the design:**

1. ⭐ **P2 is the whole fix that requires no assumption**, and it should land alone first. Where the seeds
   carry no information it sends the estimate to the prior, which is the correct answer to "this is not
   identifiable".
2. **A screen is still worth having for the information-RICH, contaminated case** — `vcap` has 101 M pairs,
   so P2 will not shrink it, and its pooled fit (0.0923) is demonstrably driven by transcripts. But it can
   only be run **after** the owner asserts a bound on gDNA dispersion.
3. ⭐ **The natural bound already exists and is not new: the PRIOR, `Beta(14,14)` ⇒ `od₀ = 0.0345`** — the
   model's own a-priori statement of what gDNA strand dispersion should be. At `od = 0.05` the vcap seed is
   rejected at `p = 7e-17`. Using the prior as the screening null therefore adds **no constant**.
   ⚠ **But be honest about what that buys:** a screen against the prior can only ever *confirm* the prior.
   If the truth were `od = 0.1`, the seeds carrying that evidence would be the ones rejected. **This is a
   trade the owner must make explicitly, not a derivation.**

**Recommendation: land P2 alone. Treat any screen as a separate, owner-gated decision that begins with
"how dispersed can gDNA's strand split actually be?" — a biological question, answered once, and then used
for both the ceiling and the null.**

### ⭐ OWNER RULING (2026-07-28) — the constant must exist, and `Beta(a,a)` is where it belongs

> *"We assume gDNA is symmetric around 0.5 (biological truth). We accept overdispersion due to unknown
> sequencing factors. At some point the likelihood shifts from 'this is overdispersed gDNA' to 'this is an
> annotation error'. You are correct that this is a magic number. It needs to exist, and I propose
> `Beta(2,2)` should be that ceiling because of the mathematical behaviour of `Beta ≥ 2` versus
> `Beta ≤ 2`: at `Beta(2,2)` the probability of 0.0 or 1.0 becomes zero. In fact I think `Beta(3,3)` might
> be better. I agree we should reject seeds that violate this ceiling."*

**Accepted, and the mathematical rationale is the right one** — `Beta(a,a)` has density `∝ [p(1−p)]^(a−1)`,
which vanishes at both ends for `a ≥ 2`. That is what makes "all fragments on one strand at large `n`"
impossible under the model rather than merely unlikely, and it is why the constant belongs on `a` rather
than on a p-value.

**And `Beta(3,3)` is measurably better, which vindicates the instinct.** `P(all n on one strand)` decays as
`6/((n+2)(n+3)) ~ n⁻²` at `a = 2` but `60/((n+3)(n+4)(n+5)) ~ n⁻³` at `a = 3`. Two-sided tail p-values
(Bonferroni over 160,366 seeds ⇒ reject at `p < 6.2e-6`):

| seed | `Beta(2,2)` | **`Beta(3,3)`** | `Beta(4,4)` |
|---|---|---|---|
| **owner's 9 sense / 1 antisense** | 0.217 keep | **0.154 keep** ✅ | 0.120 keep |
| 10/10 | 0.077 keep | 0.044 keep ✅ | 0.029 keep |
| 90/100 | 0.070 keep | 0.026 keep ✅ | 0.011 keep |
| **vcap top seed** (N=1523, sense=5) | 1.08e-4 **keep** ⛔ | **1.88e-6 REJECT** ✅ | 3.85e-8 REJECT |
| vcap #2 (N=1392, sense=13) | 6.4e-4 keep | 2.4e-5 keep ⚠ | 1.0e-6 REJECT |
| vcap #3 (N=1200, sense=9) | 4.5e-4 keep | 1.5e-5 keep ⚠ | 5.6e-7 REJECT |

**`Beta(2,2)` cannot reject even the worst contaminant** — the seed that is 12 % of the pooled numerator on
its own. **`Beta(3,3)` rejects it while comfortably keeping every legitimate-looking split**, including the
9/10 case the owner raised specifically. `Beta(4,4)` also catches #2 and #3, at the cost of a tighter
admissible band (`od ≤ 0.111`).

`od = 1/(2a+1)`: **a = 2 ⇒ 0.200, a = 3 ⇒ 0.143, a = 4 ⇒ 0.111.**

⚠ The choice is finely balanced between a = 3 and a = 4, and it is genuinely a judgement about how
dispersed gDNA's strand split can be from sequencing effects alone. **a = 3 is the conservative choice**
(more dispersion admitted ⇒ weaker strand likelihood ⇒ the direction the governing principle prefers) and
it clears the owner's own worked example by a factor of 25,000.

### P1 (superseded — kept for the record) — screen each seed against its OWN premise (addresses D1)

A seed entering a **mean-½** fit asserts its gDNA fragments are Binomial(n_c, ½). Test that assertion with
the statistic the estimator already forms:

```
    z_s = (sense_s − N_s·mean_s) / sqrt(N_s·mean_s(1 − mean_s))      and note  excess_s ≡ N_s·mean_s(1−mean_s)·(z_s² − 1)
```

Reject seeds whose `|z_s|` contradicts the premise. **This is not circular**: the fit's estimand is the
second moment, the screen tests the first, and they are orthogonal under the null. It is also *N-aware for
free* — a 2-fragment seed cannot produce a large `|z|`, so tiny seeds are never rejected for being tiny.

**Validation on ground truth** (the synthetic simulator is multinomial ⇒ **true od = 0**):

| condition | shipped `pooled` | **after the screen** |
|---|---|---|
| gdna100 ss0.99 nrna_present capON | **0.1345** ✗ | **0.0001** ✅ |
| gdna100 ss0.50 nrna_none capOFF | −0.0000 ✅ | −0.0000 ✅ |
| gdna300 ss0.99 nrna_none capOFF | −0.0000 ✅ | −0.0000 ✅ |

**The shipped estimator invents 0.13 of dispersion where the truth is exactly 0, and the screen removes it.**

⚠ `|z| ≤ 3` is a constant and needs owner sign-off. Two escapes: it is the classical Gaussian-outlier
convention rather than a fitted knob, and the alternative — an information-weighted quantile — is
cutoff-free but degenerates on small seeds (see the retraction table). **Recommend: the `|z|` screen, with
the value recorded in the constants ledger as asserted-not-derived.**

### P2 — ⭐ shrink toward the prior by INFORMATION, not by seed count (addresses D2)

#### What a "pair count" is, and why it is the right currency

**Overdispersion is a statement about CORRELATION BETWEEN FRAGMENTS inside one node** — "if this fragment
is sense, is the next one from the same node more likely to be sense too?" To see a correlation you need
**two fragments to compare**. One fragment tells you nothing: it has nothing to be correlated *with*.

So the unit of evidence about dispersion is a **PAIR of fragments in the same seed** — not a fragment, and
certainly not a seed:

| seed size `n` | pairs `n(n−1)/2` | what it can tell you |
|---|---|---|
| 1 | **0** | nothing at all |
| 2 | 1 | one comparison |
| 10 | 45 | |
| 1,523 | 1,158,903 | |

Note the estimator's own denominator is `scale_s = n_c(n_c−1)·¼` — **exactly twice the pair count times ¼**.
**The estimator already weights by pairs internally.** The bug is only that the *prior shrinkage* switched
currency to seed count on the way out.

The analogy: estimating a variance from 160,000 samples **of size 1** gives you nothing, no matter how many
samples you have. Counting the samples says "160,000 observations, the prior is irrelevant". Counting the
pairs says "zero information, keep the prior". The second is right.

#### The change

Replace `n_seed_nodes` in the shrinkage with the **effective number of dispersion-informative pairs**,
`Σ_s n_c(n_c − 1)/2` (equivalently: make `prior_weight` a weight in the same currency as `Σ scale`). Then:

* a library with 160 k singleton seeds correctly **falls back to the prior**, because it has no information;
* `vcap`, with 101 M pairs, correctly **keeps its fit**;
* it is not a new constant — it is a **units fix**. The prior's strength was always meant to be "how much
  data would outvote it", and seed *count* is the wrong measure of that for a second moment.

This is the change most likely to release the saturation honestly, because saturation on LBX0190/MO_3021 is
a power failure, not a dispersion measurement.

### P3 — keep the ceiling, and report saturation as QC (no change, but state it)

The `Beta(2,2)` ceiling is **the guard**, and it is currently the only thing preventing annotation error
from propagating into the strand likelihood — **the only intrinsic gDNA/RNA information source there is**
(`CALIBRATION_ARCHITECTURE.md`). Saturation is the **canary** that the seed set is contaminated or the
library is underpowered. Keep it, and surface it in QC rather than silently clipping.

⛔ **Do not raise it.** Inflating the overdispersion to accommodate contamination dilutes our strongest
evidence — and the human transcriptome is more complex than the annotation records, so this will *always*
be a live contaminant on real data. The model must respect that rather than absorb it.

### ⛔ Rejected candidates, measured

| candidate | result |
|---|---|
| percentile trim of the top x % by excess variance | fails to release LBX0190 / MO_3021 at 0.1/1/5 % — it is a knob and it treats the symptom |
| information-weighted median / IQR-trimmed mean | degenerates where the median seed has n_c = 2 (per-seed ratio is exactly ±1). Correct on `vcap` (−0.006), meaningless on LBX0190 (1.0000) |
| bounded-influence cap at `Σscale/√n` | **worse than pooled** on the synthetic (0.5583 vs 0.1345 against a truth of 0) |

### ⚠⚠ THE gDNA ARGUMENT DOES **NOT** TRANSFER TO THE RNA FIT — measured 2026-07-28

`fit_rna_strand_overdispersion` is the twin of the gDNA fit with two substitutions: `node_mean = κ` and
`component_mean = κ`. Both break the reasoning above, and the second breaks it fatally.

**1. `κ` is a FITTED library parameter, not a biological truth.** "gDNA is symmetric around ½" licenses
calling a displaced seed an error. There is no equivalent statement for RNA: a spliced seed far from `κ` may
be **real antisense biology** — an antisense gene, a bidirectional promoter, a readthrough — and rejecting
it would delete signal rather than contamination.

**2. The ceiling stops being a ceiling exactly where stranded libraries live.** `Beta(a,a)` is symmetric, so
it is only the right null at mean ½. For RNA it would have to become `Beta(2aκ, 2a(1−κ))`, and the
owner's "density vanishes at both ends" property requires **both** shape parameters `> 1`:

| a | κ | implied Beta | density at p = 1 | density at p = 0 |
|---|---|---|---|---|
| 3 | 0.50 | Beta(3.00, 3.00) | vanishes ✅ | vanishes ✅ |
| 3 | 0.90 | Beta(5.40, **0.60**) | vanishes | ⛔ **diverges** |
| 3 | **0.99** | Beta(5.94, **0.06**) | vanishes | ⛔ **diverges** |

At `κ = 0.99` the antisense shape parameter is **0.06** — the density blows up at `p = 0`, so an
all-antisense seed is not merely admissible, it is *favoured*. **There is no ceiling to violate.**

**⇒ The screen (P1) is licensed for the gDNA fit ONLY.** The RNA overdispersion needs its own treatment,
and it does not have the symmetric anchor that makes the gDNA case tractable. ⚠ Note `od_r` saturates on
**4/4** real samples — worse than `od_g` — so this is not a hypothetical gap.

### ⛔⛔⛔ AND `od_r` IS NOT AN ESTIMATE AT ALL — IT IS ABSORBING A MEAN MISFIT (2026-07-28)

**Owner's reasoning, and it is correct:** *"We are already constrained to ANNOTATED splice junctions. These
are highly reliable. The danger of modelling unannotated transcription does not exist for RNA. So I don't
think we need an overdispersion ceiling. I'm interested to see what the RNA overdispersion estimate is on
real data."*

Measured on the four real cfRNA libraries, over **exactly the seeds the fit uses**:

| sample | **κ handed to the fit** | **EMPIRICAL sense frac of those seeds** | ratio | **raw pooled od** | clipped |
|---|---|---|---|---|---|
| LBX0190 | 0.00231 | 0.30316 | **131×** | **10.73** | 0.2000 |
| LBX0588 | 0.01203 | 0.33107 | 27.5× | **13.55** | 0.2000 |
| MO_3021 | 0.00203 | 0.46607 | **229×** | **79.88** | 0.2000 |
| vcap | 0.00006 | 0.06811 | **1142×** | **12.61** | 0.2000 |

⛔ **`od` is an intraclass correlation. It is bounded by 1 by definition. The raw estimates are 10.7–79.9.**
That is not a dispersion at all — it is mathematically impossible, and the ceiling is the only reason the
output is finite.

**The cause is a MEAN misfit, not dispersion.** `excess_s = (sense_s − N_s·κ)² − N_s·κ(1−κ)`, so a κ that
is wrong by 131× puts a systematic offset straight into the numerator, squared. The estimator is measuring
the distance between κ and its own seeds and reporting it as overdispersion.

### ⛔ CORRECTION (same day) — "κ IS WRONG" WAS THE WRONG ATTRIBUTION. **κ IS RIGHT.**

The first draft of this section blamed κ and hypothesised a transcript-relative / motif-relative frame
mismatch. **Measured, and refuted: κ agrees with the deep junctions almost exactly.** What is wrong is the
**seed population**.

**The owner's framing is the correct physical one and it is what the measurement confirms.** During dUTP
prep the reverse-strand synthesis is degraded so a fragment can be traced to its strand; the process is
imperfect. The aligner reports only Forward/Reverse, which is **not** a strand call — but **iff a fragment
is SPLICED, the genomic GT/AG motif gives the true strand independently of prep**. Comparing motif strand
to aligner orientation over spliced fragments therefore **measures library-prep efficiency**, and κ is its
mean. ⚠ And the owner's power point is the crux: **at κ ≈ 0.99 (or ≈ 0.01) you need a junction with many
reads to see even one disagreeing read.**

**Agreement rate by junction depth — this is the whole story:**

| sample | κ | all seeds | depth ≥ 10 | **depth ≥ 100** | **depth ≥ 1000** |
|---|---|---|---|---|---|
| LBX0190 | **0.00231** | 0.3032 | 0.2682 | **0.0029** ✅ | **0.0038** ✅ |
| LBX0588 | **0.01203** | 0.3311 | 0.0319 | **0.0011** ✅ | — |
| MO_3021 | **0.00203** | 0.4661 | **0.4876** ⚠ | **0.0136** ✅ | **0.0024** ✅ |
| vcap | **0.00006** | 0.0681 | 0.0542 | 0.0083 | **0.0014** ✅ |

**The deep junctions reproduce κ.** LBX0190's κ = 0.0023 against 0.0029 / 0.0038 measured at depth. κ is a
sound estimate of prep efficiency.

**The seed population is not.** It is dominated by *millions* of shallow boundary sides — 818 k seeds for
2.4 M fragments, i.e. ~1–2 each — whose agreement rate sits at **0.30–0.49**. MO_3021 at depth ≥ 10 reads
**0.4876**: dead on ½, the signature of **effectively random orientation**, not of prep efficiency. Those
shallow sides disagree with κ by ~100×, and `excess = (sense − N·κ)²` turns that offset into an enormous
numerator.

⚠ **A related but MINOR defect, worth fixing anyway:** a seed with `n = 1` contributes
`excess ≈ (1 − κ)² ≈ 1` to the numerator but `scale = n(n−1)·κ(1−κ) = 0` to the denominator — it inflates
`od` without being able to inform it. Measured share of `Σexcess`: **0.1 – 9.1 %**, and excluding those
seeds moves `od` only 10.73 → 10.22. **Real, but not the cause.** (It is the same principle as P2: a seed
with one fragment has no pair to correlate. Here the asymmetry is a plain bug — it counts in one sum and
not the other.)

### ⭐⭐ THE HONEST RNA OVERDISPERSION, AND THE OWNER WAS RIGHT

Restricted to junctions deep enough to carry the information, and re-centred on their own mean so a mean
misfit cannot masquerade as dispersion:

| sample | depth ≥ 100 | depth ≥ 1000 | implied `Beta(a,a)` |
|---|---|---|---|
| LBX0190 | 0.0047 | **0.0035** | a ≈ 107–141 |
| LBX0588 | **0.0020** | — | a ≈ 250 |
| MO_3021 | 0.0158 | **0.0023** | a ≈ 31–218 |
| vcap | 0.0075 | **0.0011** | a ≈ 66–464 |

> **The true RNA strand overdispersion is 0.001–0.016 — 13 to 180× BELOW the ceiling of 0.2** — and it
> agrees with the synthetic suite's 0.0008–0.0017. **The owner's prediction was exactly right on both
> counts: there is no danger of overestimating RNA overdispersion, and it is hard to measure — only
> 16–3,011 junctions per library reach depth 1000.**

**⇒ THE FIX, and it needs no new constant:** the RNA overdispersion must be fitted from junctions that can
inform it. That is P2's principle again — weight by **pairs**, `n(n−1)/2`, which is zero for a singleton
and negligible for a shallow side — plus repairing the numerator/denominator asymmetry so a seed that
contributes no denominator contributes no numerator either. A depth *threshold* is then unnecessary: the
pair weighting demotes shallow sides automatically.

**⇒ AND THE CEILING QUESTION IS ANSWERED:** the owner proposed dropping the RNA ceiling on the grounds that
annotated junctions are reliable. **That is right in principle, and the measurement supports it** — the
honest value is 100× below the ceiling, so the ceiling would never bind. ⚠ **But do not drop it until the
seed weighting is fixed**, because today it is the only thing standing between the strand likelihood and an
`od` of 10–80.

**⇒ THREE CONSEQUENCES:**

1. ⛔ **Do NOT remove the RNA ceiling yet.** The owner's reasoning about contamination is right, but the
   saturation is not contamination — and without the ceiling `od_r` would be **10–80**, destroying the
   strand likelihood on every real library. **The ceiling is currently load-bearing for a reason nobody
   intended.**
2. ⭐ **`od_r` carries no information today.** Every real library gets the same clipped 0.2 regardless of
   its data, so the RNA strand channel is running at maximum dispersion — i.e. **maximally uninformative**
   — on all real data. That is a much bigger deal than the gDNA side.
3. **The fix is upstream of the estimator**: get the mean and the seeds into the same orientation frame.
   No amount of robustification helps a mean that is wrong by three orders of magnitude.

⚠ **The synthetic suite hides this** (`od_r` = 0.0008–0.0017 there), so it is another defect only real data
could show — and the reason to check the frames explicitly rather than infer them.

## 3b. ⛔⛔⛔ INDEPENDENT REVIEW (11-agent workflow, 2026-07-28) — THREE OF MY CLAIMS ARE REFUTED

Five independent derivations, each adversarially verified, plus a re-measurement that reproduces the
shipped `od_mom` exactly on all 7 datasets. **Read this before acting on anything above it.**

### ⛔ R1 — "`Beta(3,3)` is measurably better" is VOID. A ceiling screen rejects NOTHING at a = 2 OR a = 3.

**My §OWNER-RULING table applied LBX0190's multiplicity (m = 160,366) to a *vcap* seed.** With vcap's own
`m = 863,246` that seed **survives**. Measured directly by two derivations: **at a ∈ {2,3}, every
multiple-testing frame rejects ZERO seeds on all four libraries.** The reason is structural —
`p_min ∝ n^(−a)`, so a Bonferroni level `α = 1/S` needs `n ≳ S^(1/a)` (= 552 at vcap) while the
contaminants are mostly small.

**⇒ Owner decision D-c ("reject seeds that violate the ceiling") does not survive contact with the data as
specified.** The ceiling is a sound *admissibility clamp* and a sound *QC canary*; it is not a usable
screen at any a the owner is willing to assert.

⚠ **But the ceiling choice is still the ONE decision here that moves real data.** `a = 2 → 3` takes
`od_max` 0.200 → 0.1429 — a **1.4× sharpening of the strand likelihood on every node of every real
library** (it changes `od_g` on LBX0190 and MO_3021 and `od_r` on all four). Directionally that is *toward*
the measured truth (0.001–0.03). **It needs its own pre-registered A/B with stranded conditions as the
falsifier, and must NOT be bundled with the constants cleanup.**

### ⛔ R2 — the prior-currency fix (my "D2, the single clearest defect") changes NOTHING on real data

Measured: fitted `od_g` **identical to 4 decimals** at `W ∈ {30, 100, 300, 588, 2658}` pairs **and with the
shrinkage deleted entirely**. The prior binds only when `I ≲ 10⁴` pairs, which never happens on a real
library (they carry 0.7 M – 101 M).

**The saturation is BIAS, not power.** `od_mom = 0.9313` on LBX0190 against `SE|₀ = 0.0012` sits
**66–600 σ above the ceiling**. No shrinkage rule reaches that. **My "power failure reported as a
measurement" framing was wrong** — it is a *contaminated / non-exchangeable seed channel*.

✅ **Do the cleanup anyway** — it is correct, it is honest, and it is what makes the low-information corner
(unit fixtures, and the `denom ≤ 0` fallback that fires on 4/32 synthetic conditions) behave — but expect
**no change**, and treat any movement as a bug.

### ⛔ R3 — "more dispersion = weaker, so the governing principle prefers it" is FALSE for this parameter

`od` is **not** a confidence knob. It is the gDNA component's *specificity at its own mean ½*, so
inflating it degrades **both** directions at once. Measured on a balanced pure-gDNA `N = 200` node, raising
`od₀` to 0.1 collapses the strand log-evidence **305.4 → 113.9 nats** (and it *saturates in N* — depth does
not recover it), while a clean pure-RNA node's posterior `f_g` rises 0.034 → 0.080. This is a real
correction to how the "pass-0 must be WEAK" principle applies here.

**⇒ `od₀ = od_max/2 = 0.1` (derived independently by two agents as the max-entropy mean) is REJECTED.** It
also breaks a live tested bug fix (`test_bp_solver.py:809` asserts `od₀ < 0.05`) and pushes `node_init`'s
strand deadband from κ\* = 0.631 to 0.724, killing `I_strand` for partially-stranded protocols.

### ✅ WHAT SURVIVED, AND IT IS A REAL CLEANUP

**4 config fields + 2 duplicated literals + 1 wrong currency → 2 asserted constants + 1 structural zero,**
all in one module:

| | constant | value | jobs | status |
|---|---|---|---|---|
| K1 | `_CEIL_ALPHA_BETA` | 2.0 ⇒ `od_max` = 0.2 | admissibility clamp; QC canary | **ASSERTED** (owner's `a ≥ 2` rationale) |
| K2 | `_PRIOR_ALPHA_BETA` | 14.0 ⇒ `od₀` = 0.0345 | shrinkage target; `denom ≤ 0` fallback; shared gDNA+RNA | **ASSERTED, now bracketed by measurement** |
| K3 | floor | 0.0 | lower clamp | **STRUCTURAL** — an ICC cannot be negative |

* ⭐ **`W` stops being asserted.** `Var(od_mom)|₀ = 1/I` with `I = Σ N(N−1)/2` is **exact at μ = ½** (not
  asymptotic; MC-verified), so the shrinkage is literally a Normal–Normal posterior mean and `W = 1/τ²`
  with `τ²` the prior variance — **909 pairs**, derived from K1 and K2. ⚠ Name the step honestly: reducing
  a bounded prior to two moments is **Bühlmann linear credibility**, not exact Bayes.
* ⭐ **`Beta(14,14)` is no longer "completely arbitrary".** Two assumption-light measurements now bracket
  the truth: gDNA from the exact `od(n=2) = 2·P(same strand) − 1` readout (vcap plateau **0.007–0.028**),
  and RNA from deep junctions (**0.0011–0.0158**). `od₀ = 0.0345` sits 1.2–30× *above* the top of both —
  the conservative end of measured reality. **Keep 14.**
* ⭐ **`gdna_weight ≡ 1.0` on 100 % of fragment-carrying seeds** (min − 1 = −2.2e−16), so `node_mean ≡ ½`
  and the entire gDNA/RNA mixture in that fit is **dead code**.
* ⚠ **A floating-point trap that corrupted my own D2 table:** at `w = 1 − 2.2e−16`, a 2-fragment seed gives
  `n_c = 1.9999999999999996`, so `n_c ≥ 2` under-counts. True informative-seed count on LBX0190 is
  **109,875**, not the 107,534 I reported. **Use the integer `N`.**
* ⚠ **`I = pairs` holds ONLY at μ = ½ — it is WRONG for RNA.** Measured `I_r/pairs` = **0.140 / 0.094 /
  0.095 / 0.050**. The RNA information is 5–14 % of the pair count and needs its own closed form.

### ⛔ R4 — do NOT refit the RNA mean from its own seeds

Tempting given §3's finding, and refuted: the inflated `od_r` is currently **the only thing pricing the
mean error** in the applied likelihood (the χ² term is exactly 1.00 at `od = b²/κ(1−κ)`). Removing it while
`_mixture_strand_loglik` still uses κ̂ makes the likelihood **287–1088× over-confident**. And it does not
even work — `clip(od_mom(κ_seed))` releases the ceiling on **1 of 4** libraries. **Diagnose the population
split first** (§3's shallow random-orientation sides), then revisit.

### ⚠ R5 — the missing gate nobody has built

The synthetic suite has `od = 0` **by construction**, so it can only reward estimators that return zero —
it validates **one side** of a two-sided question. **Build an injection harness**: at each library's *real*
seed-size distribution, inject Beta-Binomial dispersion at `od_true ∈ {0.005 … 0.2}`, with and without
contamination, and pre-register both **(i)** recover `od_true` and **(ii)** do not manufacture dispersion at
`od_true = 0`. The pooled MoM passes (i) and fails (ii); a mixture EM does the reverse. **Neither is
shippable until both pass.**

## 4. VALIDATION PLAN AND GATES

1. **Ground truth first.** The synthetic suite has `od = 0` by construction, so every estimator must return
   ~0 on all 32 conditions. This validates the *estimator* and is the only place truth exists.
   ⚠ It **cannot** validate the downstream benefit of a *non-zero* od — the suite is Poisson by
   construction (memory `synthetic_suite_is_poisson_omega_zero`).
2. **Real data.** The gate is not "od goes down". It is: **od should fall where the seeds are tiny (the
   prior takes over) and be earned where they are not** — LBX0190/LBX0588/MO_3021 should move toward the
   prior; `vcap` (101 M pairs) should keep a fitted value.
3. **Full 32-condition A/B at both refit settings.** `od_g` and `od_r` enter `strand_loglik` on every node
   of every library, so this is not a local change. Pre-register: **stranded conditions are the falsification
   test** — they have real strand information, so a correct de-contamination should help or be neutral
   there; a large stranded regression means the screen is removing signal, not contamination.
4. Held-fixed trust view, then goldens LAST.

## 5. WHAT THIS DOES NOT FIX

D1's root cause is upstream: `count_gdna_frac` re-encodes the annotation. P1 *detects* the resulting
contamination and removes it from the fit, but the same wrong weight is used elsewhere (it is the seed
selector and the signature-partition mask). A genuine count-clue fix — a gDNA weight that is a measurement
— belongs with the accumulator/index rework, where the path-level FL contrast can identify gDNA vs RNA
prior-free and strand-free.

**Scope this work as: make the strand model robust to a contaminated seed set. Not: make the seed set
clean.** The second is upstream.
