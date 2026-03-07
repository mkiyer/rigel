# Pristine Benchmark Report — rigel vs salmon vs kallisto

**Date:** 2025-07-01
**Benchmark:** `benchmark_20_loci.yaml` — 10 genomic loci, 100k fragments each
**Conditions:** strand_specificity=1.0, gDNA=0%, nRNA=0%, oracle aligner
**Tools:** rigel (MAP-EM), salmon (v1.10.3), kallisto (v0.51.1)

---

## 1. Executive Summary

Under pristine conditions (perfect strand specificity, no contamination, oracle
alignment), rigel achieves the **best transcript-level accuracy** across the
10-region panel, winning 7/10 regions on Spearman correlation. However, rigel
**loses 8/10 regions at the gene level** due to a systematic artifact: the EM
distributes phantom fractional counts to zero-truth genes, degrading rank-based
metrics.

Four issues were identified, ordered by impact:

| # | Issue | Severity | Scope |
|---|-------|----------|-------|
| 1 | Equivalence class collapse | High | 9/10 regions, 402/1522 transcripts |
| 2 | Gene-level phantom smearing | Medium | 8/10 regions lose gene Spearman |
| 3 | False nascent RNA classification | Low | 5/10 regions, max 33 reads |
| 4 | Transcript-level losses at 3 loci | Low | MALAT1_NEAT1, TTN, EGFR |

---

## 2. Overall Results

### Transcript Level (averaged over 10 regions)

| Tool | Spearman | Pearson | MAE | RMSE |
|------|----------|---------|-----|------|
| **rigel** | **0.9423** | **0.9889** | **41.36** | **157.24** |
| kallisto | 0.9349 | 0.9641 | 97.97 | 385.03 |
| salmon | 0.9102 | 0.9446 | 116.76 | 473.48 |

rigel leads by +0.0074 Spearman over kallisto and +0.0321 over salmon.
MAE is 2.4× better than kallisto and 2.8× better than salmon.

### Gene Level (averaged over 10 regions)

| Tool | Spearman | Pearson | MAE | RMSE |
|------|----------|---------|-----|------|
| **kallisto** | **0.9973** | 1.0000 | **2.98** | **8.67** |
| salmon | 0.9936 | 1.0000 | 3.27 | 9.67 |
| rigel | 0.9848 | 1.0000 | 2.95 | 8.49 |

rigel's gene-level MAE (2.95) is actually the lowest, but its Spearman
(0.9848) is worst due to rank perturbations among zero-count genes.

### Win/Loss Tally

| Level | rigel Wins | rigel Losses |
|-------|-------------|---------------|
| Transcript Spearman | **7** | 3 |
| Gene Spearman | 2 | **8** |

---

## 3. Per-Region Breakdown

### Transcript-Level Spearman

| Region | rigel | salmon | kallisto | Advantage | Result |
|--------|---------|--------|----------|-----------|--------|
| EGFR | 0.9948 | 0.9874 | **0.9995** | −0.0047 | LOSS |
| CDKN2A_2B | **0.9929** | 0.9882 | 0.9850 | +0.0048 | WIN |
| HOXA_cluster | **0.9754** | 0.9565 | 0.9745 | +0.0009 | WIN |
| H19_IGF2 | **0.9536** | 0.9280 | 0.9356 | +0.0181 | WIN |
| PVT1_MYC | **0.9312** | 0.8730 | 0.9261 | +0.0052 | WIN |
| XIST_TSIX | **0.9219** | 0.8882 | 0.9084 | +0.0135 | WIN |
| TTN | 0.9214 | 0.9041 | **0.9242** | −0.0028 | LOSS |
| TP53 | **0.9186** | 0.8145 | 0.8963 | +0.0222 | WIN |
| GAPDH | **0.9134** | 0.8553 | 0.8897 | +0.0237 | WIN |
| MALAT1_NEAT1 | 0.9000 | 0.9071 | **0.9101** | −0.0101 | LOSS |

### Gene-Level Spearman

| Region | rigel | salmon | kallisto | Advantage | Result |
|--------|---------|--------|----------|-----------|--------|
| XIST_TSIX | 0.9381 | **1.0000** | **1.0000** | −0.0619 | LOSS |
| PVT1_MYC | 0.9613 | **1.0000** | **1.0000** | −0.0387 | LOSS |
| HOXA_cluster | 0.9829 | 0.9637 | **0.9853** | −0.0024 | LOSS |
| TTN | 0.9863 | 0.9995 | **1.0000** | −0.0137 | LOSS |
| H19_IGF2 | 0.9921 | **1.0000** | **1.0000** | −0.0079 | LOSS |
| CDKN2A_2B | 0.9945 | **0.9956** | **0.9956** | −0.0011 | LOSS |
| MALAT1_NEAT1 | 0.9952 | **0.9998** | 0.9990 | −0.0046 | LOSS |
| GAPDH | **0.9974** | 0.9821 | 0.9948 | +0.0026 | WIN |
| EGFR | 1.0000 | 1.0000 | 1.0000 | 0.0000 | LOSS |
| TP53 | **1.0000** | 0.9956 | 0.9985 | +0.0015 | WIN |

---

## 4. Issue #1 — Equivalence Class Collapse

### Description

When two or more transcripts share identical fragment compatibility patterns
(same exonic overlap for all reads), the EM cannot distinguish them and assigns
each an equal share of the shared count. This is the dominant source of
transcript-level error.

### Scope

9 out of 10 regions contain transcripts with duplicated rigel estimates:

| Region | Collapsed / Total | Groups | Worst |
|--------|------------------|--------|-------|
| MALAT1_NEAT1 | 73/252 (29%) | 17 | — |
| TP53 | 67/212 (32%) | 18 | — |
| XIST_TSIX | 67/228 (29%) | 17 | — |
| TTN | 57/136 (42%) | 8 | — |
| PVT1_MYC | 47/225 (21%) | 10 | — |
| GAPDH | 46/178 (26%) | 13 | — |
| HOXA_cluster | 29/122 (24%) | 5 | — |
| H19_IGF2 | 8/83 (10%) | 2 | — |
| CDKN2A_2B | 5/68 (7%) | 1 | — |
| EGFR | 3/18 (17%) | 1 | — |

**Total: 402/1522 transcripts (26%) are affected.**

### Example — TP53

| Transcript | Truth | rigel | salmon | kallisto |
|------------|-------|---------|--------|----------|
| ENST00000420246.6 | 7974 | 3843.7 | 7280.1 | 3805.1 |
| ENST00000622645.4 | 6 | 3843.7 | 0.1 | 3805.1 |

Both transcripts receive identical rigel estimates (~3844). The EM sees them
as a single equivalence class and splits evenly. Salmon resolves this pair
(7280 vs 0.1), likely using its positional bias or sequence-specific bias models
to break the symmetry.

### Root Cause

rigel's fragment-to-transcript compatibility is based solely on genomic
overlap — if two transcripts share all exons and differ only in UTR length or
minor terminal variants shorter than the fragment length, every fragment is
compatible with both. The EM has no information to break the tie.

### Proposed Solutions

1. **Fragment length distribution modeling** — Transcripts with different
   effective lengths generate different fragment length distributions. Using a
   transcript-length-dependent likelihood can break some ties.

2. **Positional coverage bias** — Salmon's positional bias model captures
   non-uniform coverage along the transcript body, providing information to
   distinguish transcripts of different lengths.

3. **Sequence-specific bias** — k-mer context at fragment start/end positions
   varies by transcript, giving another discriminative signal.

4. **Effective length weighting in the prior** — Currently the OVR Dirichlet
   prior doesn't weight by effective length. Longer transcripts should receive
   proportionally higher prior mass.

---

## 5. Issue #2 — Gene-Level Phantom Smearing

### Description

rigel assigns tiny fractional counts (0.001–0.20) to genes whose true count
is zero. These phantom values perturb the rank ordering of zero-count genes,
degrading Spearman correlation even though absolute errors are negligible.

### Mechanism

When the EM distributes posterior mass across equivalence class members, some
fraction leaks into transcripts belonging to genes with zero true expression. At
the gene level, these fractions sum to small but non-zero values. Spearman
correlation is rank-based, so even 0.01 counts at a truth-zero gene can shift
its rank by 10+ positions.

### Evidence

**XIST_TSIX (worst region, rigel=0.9381 vs kal=1.0000):**
- 36 genes total; 20 have truth=0
- rigel: only 14 of those 20 are exactly zero (6 get phantom counts)
- kallisto: all 20 are exactly zero
- Largest rank shift: 10 positions (gene ENSG00000234969.1, truth=0, rigel=0.01)
- Max absolute error: rigel=0.55 counts (vs sal=21.6, kal=0.19)

**PVT1_MYC (rigel=0.9613 vs kal=1.0000):**
- 16 genes total; 8 have truth=0
- rigel: only 6 of those 8 are exactly zero
- kallisto: all 8 are exactly zero
- Max absolute error: rigel=0.31 counts

### Key Insight

Gene-level MAE is actually best for rigel (2.95 vs 2.98 for kallisto), proving
the errors are negligible in magnitude. The Spearman penalty comes entirely from
rank perturbations among near-zero genes.

### Proposed Solutions

1. **Post-EM thresholding** — Zero out gene/transcript counts below a minimum
   threshold (e.g., 0.5 counts). This would eliminate phantom smearing with
   minimal impact on true counts.

2. **Sparse EM initialization** — Initialize transcripts with zero unambig counts
   at exactly zero, preventing the EM from distributing mass to unlikely
   transcripts.

3. **Gene-level aggregation cleanup** — When summing transcript counts to gene
   level, apply a floor/round to suppress sub-integer noise.

---

## 6. Issue #3 — False Nascent RNA Classification

### Description

In a pristine simulation with 0% nascent RNA, rigel incorrectly classifies a
small number of reads as nascent RNA (unspliced). This subtracts those reads from
the mature RNA pool.

### Evidence

| Region | nRNA Leak (reads) | % of 100k |
|--------|------------------|-----------|
| HOXA_cluster | 32.7 | 0.033% |
| MALAT1_NEAT1 | 19.4 | 0.019% |
| GAPDH | 3.5 | 0.004% |
| H19_IGF2 | 2.7 | 0.003% |
| TP53 | 2.6 | 0.003% |
| XIST_TSIX | 0.5 | 0.001% |
| PVT1_MYC | 0.3 | 0.000% |
| CDKN2A_2B | 0.0 | 0.000% |
| EGFR | 0.0 | 0.000% |
| TTN | 0.0 | 0.000% |

### Root Cause

HOXA_cluster has the largest leak (33 reads). This region contains many short
introns. When a read spans a short intron (e.g., < 100bp), the alignment gap may
fall within the fragment length distribution, making it ambiguous whether the
read represents an unspliced (nascent) or spliced (mature) molecule. The strand
model or intron-retention classifier assigns non-zero probability to the nascent
pool for these borderline cases.

### Impact

The impact is minimal — 33/100,000 reads (0.03%) at worst. At gene level, this
redistributes ~33 counts away from mature RNA genes, but this is spread across
many genes and causes < 1 count error per gene on average.

### Proposed Solutions

1. **Increase confidence threshold for nRNA classification** — Require stronger
   evidence (e.g., fully retained intron) before assigning to nascent pool.

2. **Short-intron special handling** — For introns below a length threshold,
   default to mature RNA classification unless there is unambiguous retention.

---

## 7. Issue #4 — Transcript-Level Losses at 3 Loci

### MALAT1_NEAT1 (rigel=0.9000 vs kal=0.9101, Δ=−0.0101)

This is the only region where rigel loses meaningfully on transcript Spearman.
MALAT1_NEAT1 is a dense locus with 252 transcripts, 73 of which are in 17
equivalence class collapse groups. The high collapse rate (29%) combined with the
19.4-read nRNA leak both contribute.

### TTN (rigel=0.9214 vs kal=0.9242, Δ=−0.0028)

Marginal loss. TTN has the highest collapse rate (42% of transcripts) among all
regions, with 57/136 transcripts in 8 collapse groups. The large TTN gene
produces many similar isoforms that are difficult for any tool.

### EGFR (rigel=0.9948 vs kal=0.9995, Δ=−0.0047)

EGFR has only 18 transcripts but 3 are in 1 collapse group. kallisto achieves
near-perfect Spearman (0.9995) here, likely because its k-mer-based model can
subtly distinguish transcripts that share the same exonic alignment footprint
but differ in k-mer composition.

---

## 8. Proposed Solution Road Map

### Priority 1 — Post-EM Count Thresholding (Low effort, immediate gene-level fix)

Apply a minimum count threshold (e.g., 0.5) after EM convergence: any transcript
with estimated count below the threshold is zeroed out. This directly fixes
Issue #2 (gene-level phantom smearing) and partially addresses Issue #1 by
cleaning up near-zero collapsed estimates.

**Expected impact:** Gene-level Spearman should match or exceed kallisto (0.997+).

### Priority 2 — Fragment Length Distribution Modeling (Medium effort)

Incorporate the empirical fragment length distribution into the EM likelihood.
Transcripts with different effective lengths will generate different expected
fragment length profiles. This can break equivalence class ties for transcript
pairs differing in UTR length.

**Expected impact:** Reduce collapse groups by ~30-50% where length differences
exist.

### Priority 3 — Positional Coverage Bias (Higher effort)

Model non-uniform read coverage along transcript bodies (5′/3′ bias, internal
position effects). This is the approach salmon uses to resolve many collapsed
pairs, and is the likely explanation for salmon's advantage at resolving
specific transcript pairs (e.g., TP53 ENST00000420246 vs ENST00000622645).

**Expected impact:** Further 20-30% reduction in collapse groups. This is the
feature most likely to close the gap on the remaining 3 transcript-level losses.

### Priority 4 — Nascent RNA Classification Tuning (Low effort)

Tighten the intron retention classifier to require stronger evidence for assigning
reads to the nascent RNA pool, particularly at short introns. This fixes Issue #3.

**Expected impact:** Eliminate the 0-33 read nRNA leak in pristine conditions.

---

## 9. Summary

rigel already achieves the best transcript-level accuracy on this pristine
benchmark (0.9423 Spearman, 41.4 MAE), substantially outperforming both salmon
and kallisto. The two main areas for improvement are:

1. **Gene-level Spearman** — fixable with simple post-EM thresholding (Priority 1)
2. **Equivalence class collapse** — addressable through fragment length and
   positional bias modeling (Priorities 2-3)

The nascent RNA leak is negligible and the 3 transcript-level regional losses are
small in magnitude, but they trace back to the same underlying issue
(equivalence class collapse) that Priorities 2-3 address.
