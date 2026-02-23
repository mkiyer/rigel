# Pristine Benchmark Report — 10 Regions × 5 Seeds

**Date:** 2026-02-23
**Conditions:** No gDNA, no nRNA, perfect strand specificity (1.0), 50k fragments/region
**Aligner:** HISAT2
**Tools compared:** hulkrna_mm (multimap), salmon, kallisto, htseq-count (gene-level only)
**Seeds:** 101, 202, 303, 404, 505
**Abundance:** random in [0.01, 100000], log-uniform

## Executive Summary

In pristine conditions (no contamination, no nascent RNA, perfect strandedness),
hulkrna achieves **competitive but not leading** transcript-level accuracy.
Hulkrna's **mean absolute error (MAE) of 58.9** is 1.8× salmon's 33.4 and 2.5× kallisto's 24.0.
However, hulkrna **beats salmon in 6 of 10 regions** and has **zero dropout rate** (0/1665 true-positive transcripts)
compared to salmon's 14.1% and kallisto's 8.8%.
The overall gap is driven by **two catastrophic regions** where hulkrna's EM
misallocates counts among highly overlapping isoforms.

## 1. Overall Transcript-Level Metrics

| Tool | MAE | RMSE | Pearson | Spearman |
|------|----:|-----:|--------:|---------:|
| kallisto | 24.0 | 65.0 | 0.9997 | 0.9286 |
| salmon | 33.4 | 97.1 | 0.9992 | 0.9221 |
| hulkrna_mm | 58.9 | 168.5 | 0.9920 | 0.8752 |

## 2. Gene-Level Metrics

| Tool | MAE | RMSE | Pearson | Spearman |
|------|----:|-----:|--------:|---------:|
| kallisto | 1.9 | 4.4 | 0.900 | 0.897 |
| salmon | 3.7 | 9.4 | 0.900 | 0.895 |
| hulkrna_mm | 33.8 | 68.9 | 0.898 | 0.890 |
| htseq | 482.9 | 1142.8 | 0.868 | 0.851 |

At gene level, hulkrna's MAE is ~9× salmon and ~18× kallisto. This is concerning —
gene-level quantification should benefit from pooling isoform counts, yet hulkrna's
EM sometimes moves counts **between genes**, not just between isoforms.

## 3. Per-Region Transcript MAE

| Region | hulkrna_mm | salmon | kallisto | hulk/salmon |
|--------|----------:|-------:|---------:|------------:|
| chr11 (HBB) | **280.8** | 34.8 | 23.7 | **8.1×** |
| chr10 (FGFR2) | **130.3** | 82.8 | 52.2 | **1.6×** |
| chr17 (BRCA1) | 46.6 | **61.2** | 42.2 | **0.8×** ✓ |
| chr7 (EGFR) | 39.4 | **42.5** | 44.3 | **0.9×** ✓ |
| chr1 | 22.4 | 17.9 | 14.6 | 1.3× |
| chr12 | 21.3 | 12.6 | 9.8 | 1.7× |
| chr19 | 16.7 | **34.5** | 17.8 | **0.5×** ✓ |
| chr22 | 13.8 | **15.9** | 14.0 | **0.9×** ✓ |
| chr11 (FLI1) | 12.7 | **23.8** | 17.6 | **0.5×** ✓ |
| chr7 (HOXA) | 5.3 | **8.0** | 4.0 | **0.7×** ✓ |

**hulkrna beats salmon in 6/10 regions** (✓). Two outlier regions dominate the
overall MAE gap.

## 4. Pool-Level Error (Phantom Counts)

In pristine conditions (0 gDNA, 0 nRNA), hulkrna assigns small phantom counts
to nRNA and gDNA pools:

| Pool | hulkrna_mm | salmon | kallisto |
|------|----------:|-------:|---------:|
| mature_rna | −4.1 | −0.02 | +0.02 |
| nascent_rna | **+5.4** | 0.0 | 0.0 |
| genomic_dna | **+3.4** | 0.0 | 0.0 |

hulkrna attributes ~5 fragments/region to nRNA and ~3 to gDNA when both
should be zero. This is expected from the 3-pool EM: with prior-smoothed
nRNA/gDNA parameters, a few ambiguous fragments get distributed to non-RNA pools.
The magnitude (5 + 3 = 8 out of 50,000 = 0.016%) is negligible for practical use.

## 5. Dropout Analysis

| Tool | Dropouts | Rate |
|------|--------:|-----:|
| hulkrna_mm | 0/1665 | **0.00%** |
| kallisto | 147/1665 | 8.83% |
| salmon | 235/1665 | 14.11% |

**hulkrna has zero dropouts** — every transcript with true expression receives
a non-zero count estimate. Salmon and kallisto both drop 9-14% of true-positive
transcripts to zero. This is a significant advantage for hulkrna in detecting
low-abundance transcripts.

## 6. Error by Abundance Quartile

| Quartile | hulkrna_mm | salmon | kallisto |
|----------|----------:|-------:|---------:|
| Q1 (low) | 22.9 | 5.8 | 4.1 |
| Q2 | 11.3 | 18.6 | 13.1 |
| Q3 | 35.9 | 30.3 | 21.6 |
| Q4 (high) | 140.3 | 108.9 | 79.3 |

hulkrna's biggest weakness is at both extremes: low-abundance (where phantom counts
from EM misallocation inflate MAE) and high-abundance (where large isoform swap errors
create absolute errors of hundreds to thousands of counts).

## 7. Over- vs. Under-Estimation

| Tool | Over-count | Mean + | Under-count | Mean − |
|------|----------:|-------:|------------:|-------:|
| hulkrna_mm | 730 | +46 | 832 | −65 |
| salmon | 708 | +41 | 797 | −49 |
| kallisto | 674 | +33 | 813 | −33 |

hulkrna under-estimates more severely than it over-estimates (mean −65 vs +46),
suggesting the EM tends to pile counts onto a few dominant isoforms while
under-counting their paralogs/siblings.

## 8. Cross-Seed Consistency

| Region | hulkrna_mm σ | salmon σ | kallisto σ |
|--------|------------:|--------:|----------:|
| chr11 (HBB) | **591.4** | 26.2 | 23.1 |
| chr7 (EGFR) | 52.1 | 56.5 | 62.0 |
| chr1 | 24.5 | 16.0 | 10.8 |
| All others | 3.5–17.9 | 3.2–40.5 | 2.7–22.4 |

hulkrna has **high variance** in the HBB locus (σ=591 MAE), meaning some seeds
produce near-perfect results while others produce catastrophic misallocation.
This suggests the EM is converging to different local optima depending on the
random read arrangement.

---

## 9. Root Cause Analysis of Worst Cases

### 9a. chr11_5225000_5310000 (HBB locus) — MAE 280.8 vs salmon 34.8

The HBB locus contains the hemoglobin beta gene family with densely overlapping
transcripts. In seed 101:

| Transcript | Gene | Truth | hulkrna | salmon | kallisto |
|-----------|------|------:|--------:|-------:|---------:|
| ENST00000642908.1 | ENSG00000284931.1 | 16,103 | **190** | 16,082 | 16,072 |
| ENST00000647543.1 | ENSG00000284931.1 | 1 | **6,948** | 3 | 7 |
| ENST00000336906.6 | ENSG00000196565.15 | 229 | **5,153** | 243 | 257 |
| ENST00000330597.5 | ENSG00000213934.9 | 0 | **4,016** | 0 | 0 |

**Pattern:** hulkrna moves ~16,000 counts from ENST00000642908.1 to three other
transcripts. ENST00000642908.1 and ENST00000647543.1 share gene ENSG00000284931.1 —
a classic intra-gene isoform swap where the EM converges to the wrong isoform.
But ENST00000336906.6 (different gene) and ENST00000330597.5 (a third gene!)
also receive phantom counts, indicating **cross-gene leakage** in the EM.

Salmon and kallisto solve this locus near-perfectly.

### 9b. chr10_121476640_121601584 (FGFR2 locus) — MAE 130.3 vs salmon 82.8

The FGFR2 gene (ENSG00000066468.24) has **41 annotated transcripts** — extreme
isoform complexity. Key errors in seed 101:

| Transcript | Truth | hulkrna | salmon | kallisto |
|-----------|------:|--------:|-------:|---------:|
| ENST00000684153.1 | 0 | **1,841** | 0 | 0 |
| ENST00000369056.5 | 0 | **231** | 0 | 41 |
| ENST00000638709.2 | 377 | **98** | 314 | 321 |
| ENST00000682550.1 | 365 | **39** | 407 | 402 |
| ENST00000360144.7 | 443 | **249** | 419 | 391 |

**Pattern:** hulkrna creates phantom counts on zero-truth transcripts (ENST00000684153.1
gets 1,841 despite abundance=2.4 and truth=0) while severely under-counting
neighbouring isoforms. The EM distributes reads to low-abundance isoforms
that happen to be compatible with the observed read positions.

---

## 10. Key Areas for Improving hulkrna Performance

### Priority 1: EM Convergence / Isoform Resolution (Critical)

The dominant source of error is **EM misallocation within and across highly
overlapping isoforms**. Specific issues:

1. **Prior/effective length interaction:** Transcripts with very low abundance
   (near-zero truth) receive substantial EM counts because the effective length
   correction makes them appear viable. The prior should more aggressively
   shrink low-abundance transcripts toward zero. Consider abundance-weighted
   prior initialization or L1 regularization.

2. **Cross-gene leakage:** In the HBB locus, counts leak from one gene
   (ENSG00000284931.1) to unrelated genes (ENSG00000213934.9, ENSG00000196565.15).
   This should not happen if gene-level resolution is working correctly.
   Investigate whether the locus EM is correctly grouping genes and whether
   gene-boundary constraints are being enforced.

3. **Local optima sensitivity:** The HBB locus shows σ=591 across seeds,
   meaning the EM converges to wildly different solutions depending on
   random read ordering. Consider:
   - Multiple random restarts with best-likelihood selection
   - SQUAREM acceleration to reduce iterations and improve convergence
   - Warm-starting from salmon-style initial estimates

### Priority 2: Effective Length Model Refinement

salmon and kallisto use analytically correct effective-length corrections
derived from the fragment length distribution. hulkrna's effective length
computation should be verified against these established implementations:
- Ensure the eCDF-based effective length matches salmon's formulation exactly
- Check that the effective length is correctly applied as a denominator
  in the EM update (counts / eff_len in the E-step)

### Priority 3: Isoform Pruning for Complex Loci

The FGFR2 gene has 41 transcripts, most with near-zero abundance. The EM
must distribute reads among all 41, creating opportunities for misallocation.
Consider:
- **Pre-filtering:** Remove transcripts with very low initial abundance
  estimates (e.g., < 1% of gene total after N EM iterations) before continuing
- **Collapse near-identical isoforms:** Many GENCODE transcripts differ by
  only a few bases. Collapsing these during EM and splitting proportionally
  afterward would reduce the effective state space
- **Adaptive regularization:** Apply stronger shrinkage to transcripts
  receiving < K fragments in early EM iterations

### Priority 4: Reduce Phantom Pool Counts

The ~5 nRNA + ~3 gDNA phantom counts per region are small but avoidable.
The 3-pool EM should apply a minimum-evidence threshold: if the nRNA/gDNA
posterior probability is below a threshold after N iterations, collapse those
pools to zero and redistribute to mRNA.

### Priority 5: Leverage Zero-Dropout Advantage

hulkrna's **zero dropout rate** is a genuine strength that salmon and kallisto
lack. This means hulkrna detects all expressed transcripts, even at very low
abundance. The quantification pipeline should be tuned to preserve this
property while improving accuracy for medium/high-abundance transcripts where
the EM currently struggles.

---

## 11. Summary Scorecard

| Metric | hulkrna | salmon | kallisto | Winner |
|--------|--------:|-------:|---------:|--------|
| Transcript MAE | 58.9 | 33.4 | **24.0** | kallisto |
| Gene MAE | 33.8 | 3.7 | **1.9** | kallisto |
| Pearson r | 0.992 | 0.999 | **1.000** | kallisto |
| Spearman ρ | 0.875 | 0.922 | **0.929** | kallisto |
| Dropout rate | **0.0%** | 14.1% | 8.8% | **hulkrna** |
| Regions beating salmon | **6/10** | — | — | **hulkrna** |
| Cross-seed stability | Low | Med | **High** | kallisto |
| gDNA/nRNA separation | **Yes** | No | No | **hulkrna** |

hulkrna has unique capabilities (zero dropout, 3-pool separation, multimap
handling) but needs improved EM convergence to match salmon/kallisto in
pristine conditions. **Fixing the EM isoform resolution in complex loci
would likely cut hulkrna's MAE by 50%+**, bringing it to parity or better
than salmon.
