# nRNA Siphon: Complete Root Cause Diagnosis

## Executive Summary

The nRNA siphon is caused by **asymmetric spatial reach** between mRNA and nRNA
components in the EM mixture model. Intronic-only fragments inflate θ_nRNA, and
that inflated mixing weight then overwhelms the per-fragment mRNA bias advantage
when allocating shared exonic fragments. This is a model misspecification, not a
bug or convergence failure.

---

## Scenario

| Parameter | Value |
|---|---|
| Gene | TA1 (3 exons: [1000,1020], [5000,5500], [12000,13000]) |
| Exonic length | 1,520 bp |
| Genomic span | 12,000 bp |
| mRNA abundance (TA1) | 20 |
| nRNA abundance (NTA1) | 500 |
| n_fragments | 10,000 |
| Expected mRNA | 43 |
| **Observed mRNA (VBEM)** | **27.3 (36.6% error)** |
| **Observed mRNA (MAP-EM)** | **28.6 (33.5% error)** |

## Fragment Routing

Of 10,000 total fragments:

| Category | Count | Description |
|---|---|---|
| Deterministic mRNA | 5 | Splice-junction spanning → locked to mRNA |
| Both mRNA + nRNA | 1,499 | In exonic regions → EM must split |
| nRNA only | 8,496 | In intronic regions → unambiguously nRNA |

Of the 1,499 shared fragments (pre-bias-correction log-likelihood comparison):

| Category | Count | Effect |
|---|---|---|
| Exactly equal (mRNA LL == nRNA LL) | 877 | Fragments fully within exons; identical scores |
| nRNA favored (overhang penalty on mRNA) | 621 | Near exon boundaries; penalties up to 1377 nats |
| mRNA favored | 1 | Fragment that extends past nRNA span |

## Discriminative Signals

### 1. Per-fragment bias correction (effective length)

The C++ EM applies `ll[i] -= log(max(profile_length - frag_len + 1, 1))`:

- mRNA: profile_length = 1,520 → correction ≈ −7.15 (for frag_len ≈ 250)
- nRNA: profile_length = 12,000 → correction ≈ −9.37
- **mRNA advantage: +2.22 nats per fragment**

After bias correction, ALL 877 "equal" fragments become mRNA-favored.

### 2. Fragment length model

The FL model has only 5 splice-junction observations → degenerates to flat uniform
(all log_prob = −6.9137). Provides **zero discriminative power**.

### 3. Coverage weight

mRNA mean=1.164, nRNA mean=1.045. Slight mRNA advantage (~0.12 per fragment), but
this is only used in the warm-start, not in the E-step posterior.

### 4. Overhang penalty

log(0.01) ≈ −4.605 per base of overhang. Even 1 base of overhang is −4.605 nats,
far exceeding the 2.22 nat bias correction advantage. The 621 nRNA-favored
fragments remain locked to nRNA after bias correction.

## The Siphon Mechanism

### E-step posterior formula

```
P(mRNA|f) ∝ exp(ll_corrected_mRNA + log(θ_mRNA))
P(nRNA|f) ∝ exp(ll_corrected_nRNA + log(θ_nRNA))
```

### The mixing weight imbalance

At EM convergence:
- θ_mRNA ≈ 28.5 (from 5 det + 23.5 EM)
- θ_nRNA ≈ 9,971 (from 8,496 intronic-only + 1,475 from exonic)

For the 877 equal-scored fragments (after bias correction):
- log(θ_mRNA) = log(28.5) = 3.35
- log(θ_nRNA) = log(9,971) = 9.21
- **θ ratio disadvantage for mRNA: −5.86 nats**

Net for mRNA: +2.22 (bias correction) − 5.86 (θ ratio) = **−3.64 nats**

Result: mRNA gets only **1.57%** of shared fragments.

### Self-consistent fixed point

The EM converges to a unique fixed point where θ_mRNA = 28.5. This was verified by:

1. **Python 2-component MAP-EM**: converges from θ_m=754 down to 28.53 in 42 iterations
2. **C++ MAP-EM**: gives 28.6 (28.56 with gDNA accounting)
3. **C++ VBEM**: gives 27.3 (slightly lower due to digamma vs log)

The fixed point is BELOW truth because of negative feedback:
1. θ_mRNA starts lower than truth
2. → mRNA posterior decreases
3. → fewer fragments allocated to mRNA
4. → θ_mRNA decreases further
5. → stabilizes at 28.5

### Why effective length normalization is not enough

The current per-fragment bias correction `ll -= log(eff_len)` IS the effective length
normalization — it's mathematically equivalent to dividing each fragment's likelihood
by the positional rate. This is the standard RNA-seq EM formulation.

The problem is not missing normalization. The problem is that the EM's generative
model treats θ as a **global** mixing weight:

> P(fragment from component c) ∝ θ_c

This means intronic-only fragments that are 100% nRNA contribute to θ_nRNA, which
then applies **uniformly** to all fragments — including exonic ones where mRNA is
a legitimate candidate. The model has no mechanism for "local" mixing weights.

## Sensitivity Analysis

| Scenario | Predicted mRNA | Notes |
|---|---|---|
| Full EM (current) | 28.5 | Standard intronic dilution |
| No intronic dilution (θ_n excludes intronic-only) | 808 | Massive overestimate |
| At TRUE θ values (θ_m=43, θ_n=9957) | 39.7 | Even truth can't recover 43 |

Even injecting the true θ values only recovers ~40 mRNA (not 43), because 877
shared exonic fragments include both true mRNA AND true nRNA. The EM can only
separate them via the 2.22 nat bias advantage, which gives mRNA ~3.8% share
at the true θ ratio.

## Not the Cause

| Factor | Evidence |
|---|---|
| VBEM digamma effect | MAP gives 28.6 vs VBEM 27.3 — nearly identical |
| gDNA absorption | Only 0.64 fragments to gDNA despite gdna_init=0.882 |
| SQUAREM convergence | Python scalar EM matches C++ SQUAREM exactly |
| Prior influence | total_pseudocount=1.0 vs 10,000 data fragments |
| Coverage weights | Only used in warm-start, not E-step posterior |

## Core Issue (Generalized)

This is the same mechanism as the mega-locus gDNA problem, at a smaller scale:

- **nRNA siphon**: nRNA's large genomic span creates intronic-only evidence that
  inflates θ_nRNA, suppressing mRNA in exonic regions
- **gDNA siphon**: gDNA is a universal candidate for all unspliced fragments,
  inflating θ_gDNA, suppressing nRNA in gene-specific regions

Both are instances of **a component with broad spatial reach accumulating θ from
unambiguous regions, then using that θ to steal fragments from narrower components
in ambiguous regions**.
