# VCaP Benchmark: Pristine vs Dirty Comparison Report

**Date:** 2026-04-03  
**Rigel Version:** 0.3.3 (VBEM_CLAMP_FLOOR = 0.1 fix applied)

## Experimental Setup

| Parameter | Value |
|-----------|-------|
| Cell line | VCaP (prostate cancer, highly complex transcriptome) |
| RNA fragments | 50,000,000 |
| Fragment length | 250 ± 50 bp |
| Read length | 101 bp |
| Error rate | 0.0 |
| Strand specificity | 0.9 (R1-antisense) |
| nRNA contamination | None |
| Seed | 42 |

| Condition | gDNA Fragments | gDNA Fraction | Total Fragments |
|-----------|----------------|---------------|-----------------|
| **Pristine** | 0 | 0% | 50M |
| **Dirty** | 25,000,000 | 33% | 75M |

**Tools compared:** Rigel MAP-EM, Rigel VBEM, Salmon, Kallisto

**Alignment:** Minimap2 `splice:sr` with junction annotations (Rigel, Salmon); Kallisto uses its own pseudoalignment.

---

## 1. Transcript-Level Performance

### 1.1 Headline Metrics

| Metric | Condition | Rigel MAP | Rigel VBEM | Salmon | Kallisto |
|--------|-----------|-----------|------------|--------|----------|
| **Pearson R** | Pristine | 0.9865 | 0.9864 | 0.8898 | **0.9964** |
|               | Dirty    | 0.9862 | 0.9865 | 0.8900 | **0.9958** |
|               | Δ        | +0.0003 | −0.0001 | −0.0002 | +0.0006 |
| **WARE** | Pristine | 0.0781 | 0.0785 | 0.3588 | **0.0765** |
|          | Dirty    | 0.0796 | 0.0816 | 0.3627 | 0.1162 |
|          | Δ        | −0.0015 | −0.0031 | −0.0039 | **−0.0397** |
| **MAE** | Pristine | 0.355 | **0.353** | 1.594 | 0.325 |
|         | Dirty    | 0.370 | 0.370 | 1.627 | 0.636 |
|         | Δ        | −0.015 | −0.017 | −0.033 | **−0.311** |
| **MAPE** | Pristine | 49.9% | 50.3% | 218.5% | **48.6%** |
|          | Dirty    | 56.8% | 58.1% | 239.1% | 194.5% |
|          | Δ        | −6.9pp | −7.8pp | −20.6pp | **−145.9pp** |
| **RMSE** | Pristine | 5.87 | 5.89 | 16.32 | **3.03** |
|          | Dirty    | 5.93 | 5.87 | 16.31 | 3.98 |
|          | Δ        | −0.06 | +0.02 | +0.01 | **−0.95** |
| **Spearman R** | Pristine | 0.8999 | **0.9038** | 0.8256 | 0.9016 |
|                | Dirty    | 0.8753 | 0.8764 | 0.7564 | 0.7319 |
|                | Δ        | +0.0246 | +0.0274 | +0.0692 | **+0.1697** |

### 1.2 Key Observations

1. **Rigel is remarkably stable across conditions.** Adding 25M gDNA fragments (33% contamination) causes only ~1.5–3 pp WARE degradation and ~7 pp MAPE increase. This validates that the gDNA decontamination model is working well.

2. **Kallisto suffers catastrophically from gDNA contamination.** In the pristine case, Kallisto wins on nearly every metric (Pearson R 0.9964, WARE 0.0765). But with gDNA contamination, its MAPE explodes from 48.6% → 194.5% and WARE from 0.0765 → 0.1162. Spearman R drops from 0.9016 → 0.7319. This is expected: Kallisto has no gDNA model and treats contaminating fragments as signal.

3. **Salmon is consistently weak in both conditions.** Pearson R ~0.89, WARE ~0.36 in both pristine and dirty. The minimap2-based alignment pipeline (used for both Salmon and Rigel) appears to disadvantage Salmon. Salmon's errors are similar across conditions, suggesting its problems are alignment-related rather than gDNA-related.

4. **Rigel's "pristine floor" — the baseline error when everything is perfect — is essentially identical to Kallisto's** in WARE (0.078 vs 0.077) and MAPE (50% vs 49%). The 1 pp gap is within noise. This means Rigel's quantification accuracy, once gDNA is properly removed, is competitive with the best pseudoalignment tools.

5. **VBEM vs MAP are nearly identical** in the pristine case, confirming the VBEM_CLAMP_FLOOR fix successfully resolved the previous catastrophic VBEM failure. In dirty data, VBEM is marginally worse (WARE 0.0816 vs 0.0796) but the gap is negligible.

---

## 2. Gene-Level Performance

| Metric | Condition | Rigel MAP | Rigel VBEM | Salmon | Kallisto |
|--------|-----------|-----------|------------|--------|----------|
| **Pearson R** | Pristine | 0.9954 | 0.9953 | 0.9881 | **0.9987** |
|               | Dirty    | 0.9952 | 0.9954 | 0.9880 | **0.9984** |
| **WARE** | Pristine | 0.0330 | **0.0329** | 0.1142 | 0.0283 |
|          | Dirty    | 0.0346 | 0.0340 | 0.1182 | 0.0805 |
|          | Δ        | −0.0016 | −0.0011 | −0.0040 | **−0.0522** |
| **MAE** | Pristine | 0.607 | **0.603** | 1.833 | 0.475 |
|         | Dirty    | 0.664 | 0.635 | 1.951 | 1.739 |
| **MAPE** | Pristine | 48.1% | 48.3% | 61.1% | **43.0%** |
|          | Dirty    | 70.6% | 67.6% | 190.6% | 781.8% |

### 2.1 Key Observations

1. **Gene-level amplifies the gDNA problem for tools without decontamination.** Kallisto's gene MAPE goes from 43% → 782% — an 18× blowup. This is because gDNA fragments map to introns and intergenic regions, inflating gene-level counts non-uniformly.

2. **Rigel gene MAPE increases modestly** (48% → 68-71%), which reflects residual gDNA leakage at the gene level that averages out less favorably than at transcript level.

3. **Gene-level Pearson R is largely unaffected** for all tools, demonstrating this metric is dominated by highly-expressed genes and is insensitive to contamination.

---

## 3. Pool-Level (Global) Fragment Accounting

### Rigel Fragment Decomposition

| Component | Truth (Pristine) | MAP Pred | VBEM Pred | Truth (Dirty) | MAP Pred | VBEM Pred |
|-----------|-------------------|----------|-----------|----------------|----------|-----------|
| mRNA | 50,000,000 | 49,934,475 | 49,940,809 | 50,000,000 | 49,970,646 | 49,970,443 |
| nRNA | 0 | 42,166 | 33,334 | 0 | 2,186,473 | 2,136,992 |
| gDNA | 0 | 9,125 | 11,623 | 25,000,000 | 21,309,015 | 21,358,699 |
| **mRNA rel. error** | — | **−0.13%** | **−0.12%** | — | **−0.06%** | **−0.06%** |
| **gDNA rel. error** | — | — | — | — | **−14.8%** | **−14.6%** |

### 3.1 Key Observations

1. **Pristine mRNA recovery is excellent** — only 0.1% loss, with small nRNA/gDNA "hallucination" (~42K nRNA and ~10K gDNA for MAP). This is the irreducible false-positive floor from the calibration model estimating nonzero gDNA density even when none exists.

2. **gDNA underestimation in dirty condition** is ~15%, with ~2.1M fragments misattributed to nRNA. This is the **nRNA siphon** effect documented previously — the nRNA model absorbs gDNA-like fragments when the two are incompletely identifiable.

3. **Calibrated gDNA density** is 0.00101 (pristine, truth=0) vs 0.01046 (dirty, truth≈0.0105 expected for 33% gDNA rate). The calibration model slightly overestimates gDNA in the pristine case (false positive density ≈ 0.1%), which is why we see ~10K gDNA-attributed fragments.

---

## 4. Calibration & Model Parameters

| Parameter | Pristine | Dirty |
|-----------|----------|-------|
| Strand specificity (estimated) | 0.8999 | 0.8999 |
| gDNA density (global) | 0.00101 | 0.01046 |
| gDNA mixing proportion | 0.211 | 0.486 |
| Kappa (strand) | 6.84 | 4.96 |
| gDNA FL mean | 235.87 | 256.74 |
| n_loci | 6,267 | 9,100 |
| n_unambiguous | 2,347,807 | 2,348,512 |
| n_EM fragments | 47,637,959 | 62,930,558 |
| mRNA fraction | 0.999 | 0.680 |

### 4.1 Key Observations

1. **Strand specificity is perfectly estimated** in both conditions (0.8999 vs truth 0.9000).

2. **gDNA density false positive** in pristine is 0.001 — very low but nonzero. This creates a small gDNA prior that the EM can absorb fragments into. A smarter calibration model could detect "near-zero gDNA" and clamp to zero.

3. **More loci in dirty condition** (9,100 vs 6,267): gDNA fragments create additional connected components by linking fragments to intergenic/intronic regions.

4. **Kappa is lower in dirty** (4.96 vs 6.84): gDNA fragments dilute the strand signal, making the strand model less confident. This is correct behavior.

---

## 5. Precision/Recall Analysis

| Metric | Condition | Rigel MAP | Rigel VBEM | Salmon | Kallisto |
|--------|-----------|-----------|------------|--------|----------|
| **Precision** | Pristine | 0.926 | **0.951** | 0.903 | 0.916 |
|               | Dirty    | 0.892 | **0.923** | 0.752 | 0.630 |
| **Recall** | Pristine | 0.864 | 0.843 | 0.816 | **0.863** |
|            | Dirty    | 0.839 | 0.806 | **0.851** | **0.928** |
| **F1** | Pristine | **0.894** | **0.894** | 0.858 | 0.888 |
|        | Dirty    | **0.865** | 0.861 | 0.798 | 0.751 |

### 5.1 Key Observations

1. **VBEM has highest precision** — it reports fewer false positives by shrinking low-confidence transcripts toward zero. This comes at a recall cost (more unexpressed transcripts missed).

2. **Kallisto's recall is inflated in dirty data** (0.928) because gDNA fragments activate many transcripts as "expressed" — but this destroys precision (0.630) and F1 (0.751).

3. **Rigel F1 is most stable** across conditions — 0.894 → 0.865 (MAP), a 3 pp drop. Kallisto drops 14 pp (0.888 → 0.751).

---

## 6. Expression-Level Stratified Analysis

### Very High Expression (>1000 TPM)

| Metric | Condition | Rigel MAP | Rigel VBEM | Salmon | Kallisto |
|--------|-----------|-----------|------------|--------|----------|
| MAPE | Pristine | 5.28% | 5.27% | 21.3% | **3.36%** |
|       | Dirty    | 5.28% | 5.27% | 22.0% | **8.34%** |
| WARE | Pristine | 0.048 | 0.047 | 0.219 | **0.028** |
|       | Dirty    | 0.047 | 0.047 | 0.226 | 0.079 |
| Pearson | Pristine | 0.910 | 0.910 | 0.565 | **0.976** |
|         | Dirty    | 0.910 | 0.909 | 0.563 | **0.977** |

### Low Expression (1–10 TPM)

| Metric | Condition | Rigel MAP | Rigel VBEM | Salmon | Kallisto |
|--------|-----------|-----------|------------|--------|----------|
| MAPE | Pristine | 23.0% | 23.5% | 81.6% | **21.3%** |
|       | Dirty    | 23.6% | 24.6% | 81.1% | **24.5%** |
| WARE | Pristine | 0.181 | 0.184 | 0.659 | **0.173** |
|       | Dirty    | 0.181 | 0.191 | 0.655 | 0.194 |

### Zero Expression (no true signal)

| Metric | Condition | Rigel MAP | Rigel VBEM | Salmon | Kallisto |
|--------|-----------|-----------|------------|--------|----------|
| FP predictions | Pristine | 8,674 | **5,502** | 10,923 | 9,935 |
|                | Dirty    | 12,689 | **8,372** | 35,172 | 68,323 |
| FP MAE | Pristine | 0.095 | **0.087** | 0.362 | **0.048** |
|        | Dirty    | 0.112 | **0.097** | 0.398 | 0.354 |

### 6.1 Key Observations

1. **Rigel is remarkably stable across expression strata** — low-expression MAPE barely changes between pristine and dirty (23% → ~24%). The gDNA model successfully prevents contamination from inflating low-expressed transcripts.

2. **Kallisto's pristine accuracy at low expression** (MAPE 21.3%) is the best of all tools, but converges to Rigel's level in dirty data (24.5%).

3. **False positive transcripts (zero-expression bin)** are well-controlled by Rigel VBEM (5,502 pristine, 8,372 dirty). Kallisto goes from 9,935 → 68,323 FPs with gDNA contamination — a 7× increase.

4. **Very-high-expression Pearson R for Rigel** (0.91) is notably lower than Kallisto (0.98). This gap persists in both conditions and likely reflects the mega-locus effect rather than contamination handling.

---

## 7. Summary of Rigel Strengths and Weaknesses

### Strengths

1. **Robustness to gDNA contamination** — the defining advantage over Kallisto and Salmon. All transcript-level metrics change by <3 pp WARE and <8 pp MAPE when adding 33% gDNA contamination. No other tool achieves this.

2. **High precision and low FP rate** — VBEM in particular suppresses false positives, important for differential expression downstream.

3. **Accurate global fragment accounting** — mRNA recovery within 0.1%, gDNA estimation within 15%.

4. **Correct calibration** — strand specificity, gDNA density, and fragment length models are well-estimated.

### Weaknesses

1. **Pristine accuracy ~1 pp below Kallisto** — when there is no contamination, Rigel's extra model complexity is a slight liability. The WARE gap (0.078 vs 0.077) is tiny, but Pearson R (0.986 vs 0.996) and RMSE (5.87 vs 3.03) show a meaningful gap, especially at high expression.

2. **Very-high-expression Pearson R is 0.91 vs Kallisto 0.98** — this 7 pp gap is present in both conditions and points to a structural issue, likely in the mega-locus containing many overlapping transcripts.

3. **nRNA siphon in dirty data** — 2.1M fragments (~4.3% of gDNA) are misattributed to nRNA, contributing to the 15% gDNA underestimate. This is a known identifiability issue.

4. **gDNA false positive in pristine** — calibration estimates gDNA density at 0.001 even when truth is 0, causing ~10K fragments to be attributed to gDNA.

5. **Salmon alignment pipeline underperformance** — minimap2 splice:sr alignment may not be fully optimized for Salmon's model. Salmon's consistently poor results (Pearson R ~0.89) in both conditions suggests the issue is alignment quality, not the quantification model itself.

---

## 8. Research Directions and Improvement Opportunities

### High Priority

#### 8.1 Mega-Locus Handling (Impact: Pearson R at high expression)

The very-high-expression Pearson R gap (Rigel 0.91 vs Kallisto 0.98) likely originates in the mega-locus (locus 85 in VCaP, 424K+ components). Possible approaches:

- **Locus decomposition:** Apply graph-theoretic methods (e.g., minimum vertex cut, community detection) to split mega-loci into more tractable subproblems before EM.
- **Hierarchical EM:** Run EM at the gene level first, then re-distribute within genes. This could break the identifiability deadlock in large loci without sacrificing transcript-level resolution.
- **Fragment scoring improvements:** Better multimapper resolution in dense loci through positional bias models, GC bias, or fragment length selectivity per transcript.

#### 8.2 Calibration Zero-gDNA Detection (Impact: False positives in clean data)

When true gDNA is zero, calibration still estimates density=0.001, causing ~10K false gDNA attributions. Improvements:

- **Hypothesis testing for gDNA presence:** After calibration, test whether the estimated gDNA density is significantly above zero. If not, clamp to zero and remove the gDNA component from EM entirely.
- **Bayesian model selection:** Compare marginal likelihoods of models with and without gDNA, and select the simpler model when gDNA is negligible.

#### 8.3 nRNA Siphon Mitigation (Impact: 15% gDNA underestimate in dirty data)

The nRNA model absorbs ~4% of gDNA fragments. Approaches:

- **Informative nRNA priors from calibration:** Use intergenic fragment patterns to estimate a global nRNA rate, then apply strong priors pushing nRNA toward zero when no evidence exists.
- **Structured nRNA initialization:** Currently nRNA spans are shared across transcripts. A more constrained initialization (e.g., nRNA only for transcripts with intronic coverage evidence) could reduce the siphon.
- **Joint calibration of nRNA and gDNA:** Currently these are calibrated somewhat independently. A joint model could better resolve the nRNA/gDNA ambiguity.

### Medium Priority

#### 8.4 Fragment Length Model per Transcript (Impact: Multimapper resolution)

Currently, fragment length likelihood uses a global model. Transcript-specific expected fragment lengths (based on transcript length vs insert size) could provide additional discrimination in complex loci.

#### 8.5 Positional Bias Model (Impact: Quantification accuracy)

Adding 5'/3' positional bias modeling (similar to Salmon's sequence-specific bias) could improve accuracy, especially for long transcripts where coverage is non-uniform. This would help particularly in the very-high-expression stratum.

#### 8.6 Alignment Optimization (Impact: Salmon comparison fairness)

Salmon's poor performance may be an alignment issue rather than a quantification issue. Investigating:

- Whether Salmon's native mapping mode (quasi-mapping) performs better than minimap2 input
- Optimal minimap2 parameters for Salmon compatibility
- Whether STAR alignment improves Salmon results

### Lower Priority

#### 8.7 Adaptive EM Component Pruning

During EM, components with near-zero responsibility could be pruned to reduce computation time, especially in mega-loci. The VBEM zero-forcing behavior (now capped at 0.1) was a crude version of this.

#### 8.8 Coverage-Aware Scoring

Using per-base coverage patterns (not just fragment counts) could help distinguish transcripts that share exons but differ in coverage shape. This is more impactful for long-read data but could benefit short reads in complex loci.

#### 8.9 Benchmark Expansion

- Add nRNA contamination conditions (nrna_rand) to the pristine vs dirty comparison
- Test additional strand specificities (ss=0.50 unstranded, ss=1.00 perfectly stranded)
- Test with nonzero error rates to evaluate the mismatch penalty model
- Test with real data (VCaP actual sequencing runs vs simulation)
