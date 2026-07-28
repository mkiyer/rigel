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
