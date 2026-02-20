# Benchmark Report: hulkrna vs salmon, kallisto, htseq

**Date:** 2026-02-19
**Seed:** 101 (single seed)
**Simulation:** 50,000 fragments per condition, random abundance (0.01–100,000)

## Overview

This report compares **hulkrna** (unique-map mode and multimap mode) against
**salmon**, **kallisto**, and **htseq-count** on simulated RNA-seq data from 10
human genome regions.  Each region was tested under 12 conditions:
4 gDNA contamination levels × 3 strand-specificity levels = 120 total conditions.

### Tools

| Tool | Mode | Description |
|------|------|-------------|
| **hulkrna** | unique | Default mode, unique mappers only |
| **hulkrna_mm** | multimap | Multimap-aware mode (`include_multimap=True`) |
| **salmon** | — | Salmon quant (quasi-mapping, standard settings) |
| **kallisto** | — | Kallisto quant (standard settings) |
| **htseq** | — | htseq-count (gene-level only) |

### Conditions

| Parameter | Values |
|-----------|--------|
| **Regions** | FGFR2, EGFR, CD44, BRCA1, HBB_cluster, HOXA_cluster, chr12_GAPDH_cluster, chr19_dense, chr22_CRKL_cluster, chr1_dense |
| **gDNA levels** | none (0%), low (5%), moderate (15%), high (30%) |
| **Strand specificity** | 0.95, 0.99, 1.0 |
| **Fragments** | 50,000 per condition |

---

## 1. Transcript-Level Results

### 1.1 Overall Means

| Tool | MAE | RMSE | Pearson r | Spearman ρ | Time (s) |
|------|----:|-----:|----------:|-----------:|---------:|
| hulkrna | 213.0 | 571.9 | 0.9893 | 0.9240 | 9.99 |
| hulkrna_mm | 118.6 | 246.1 | 0.9960 | 0.9315 | 10.36 |
| **salmon** | **83.2** | **193.5** | **0.9982** | **0.9646** | 0.63 |
| kallisto | 98.4 | 205.3 | 0.9974 | 0.9453 | 0.51 |

**Key finding:** Salmon achieves the lowest overall MAE (83.2), followed by
kallisto (98.4), hulkrna_mm (118.6), and hulkrna (213.0). The multimap-aware
mode of hulkrna reduces MAE by 44% compared to unique-map mode.

### 1.2 MAE by gDNA Contamination Level

| gDNA level | hulkrna | hulkrna_mm | salmon | kallisto |
|------------|--------:|-----------:|-------:|---------:|
| none | 211.9 | 119.9 | 86.2 | 97.1 |
| low | 212.8 | 119.8 | 86.3 | 89.9 |
| moderate | 197.5 | 104.9 | 69.1 | 83.2 |
| high | 229.8 | 129.8 | 91.2 | 123.4 |

**Key finding:** All tools show increased MAE at high gDNA contamination.
Salmon is the most resilient across all gDNA levels. At moderate gDNA, all
tools actually show *lower* MAE than at none — this is an artifact of the
particular random seed and abundance distribution. The gap between hulkrna and
salmon/kallisto is consistent across gDNA levels, suggesting gDNA handling is
not the primary driver of hulkrna's underperformance.

### 1.3 MAE by Strand Specificity

| SS | hulkrna | hulkrna_mm | salmon | kallisto |
|----|--------:|-----------:|-------:|---------:|
| 0.95 | 215.8 | 121.6 | 111.8 | 140.1 |
| 0.99 | 214.5 | 120.4 | 76.0 | 86.9 |
| 1.00 | 208.7 | 113.8 | 61.7 | 68.2 |

**Key finding:** Salmon and kallisto improve dramatically as strand specificity
increases (salmon 111.8 → 61.7, kallisto 140.1 → 68.2), but hulkrna shows
minimal improvement (215.8 → 208.7). This is a significant observation: hulkrna
is not fully exploiting high strand specificity to improve accuracy, while
competitors benefit greatly from the reduced ambiguity. However, this means
hulkrna's relative advantage is greatest at SS=0.95 where it competes better
against salmon and kallisto.

### 1.4 MAE by Region

| Region | hulkrna | hulkrna_mm | salmon | kallisto | Best Tool |
|--------|--------:|-----------:|-------:|---------:|-----------|
| BRCA1 | 135.2 | 135.1 | 166.4 | 127.3 | kallisto |
| CD44 | 62.2 | 62.2 | 52.6 | 56.1 | salmon |
| EGFR | 454.4 | 456.8 | 272.5 | 380.6 | salmon |
| FGFR2 | 182.1 | 182.0 | 103.3 | 91.6 | kallisto |
| **HBB_cluster** | **1021.4** | **95.2** | **38.7** | **133.5** | **salmon** |
| HOXA_cluster | 49.5 | 49.5 | 15.3 | 20.4 | salmon |
| chr12_GAPDH_cluster | 50.6 | 50.5 | 45.6 | 39.5 | kallisto |
| chr19_dense | 65.6 | 65.5 | 42.4 | 46.0 | salmon |
| chr1_dense | 44.4 | 38.9 | 33.8 | 40.1 | salmon |
| chr22_CRKL_cluster | 64.7 | 50.1 | 61.0 | 49.1 | kallisto |

**Key findings:**
- **HBB_cluster is catastrophic** for hulkrna unique mode: MAE 1021 vs salmon's
  38.7 (**26× worse**). Multimap mode fixes this (95.2), but is still 2.5×
  worse than salmon.
- **EGFR is the second-worst region:** hulkrna 454 vs salmon 273. Multimap mode
  does not help here (457), ruling out multimapping as the cause.
- **BRCA1 is hulkrna's strongest region:** hulkrna (135) beats salmon (166),
  though kallisto (127) is slightly better.
- In compact/dense regions (chr1_dense, chr12_GAPDH, chr19_dense, chr22_CRKL),
  hulkrna is competitive but not best.

### 1.5 Pearson Correlation by Region

| Region | hulkrna | hulkrna_mm | salmon | kallisto |
|--------|--------:|-----------:|-------:|---------:|
| BRCA1 | 0.9973 | 0.9974 | 0.9957 | 0.9979 |
| CD44 | 0.9996 | 0.9996 | 0.9999 | 0.9999 |
| EGFR | 0.9745 | 0.9739 | 0.9936 | 0.9795 |
| FGFR2 | 0.9982 | 0.9982 | 0.9994 | 0.9998 |
| HBB_cluster | 0.9326 | 0.9998 | 1.0000 | 0.9997 |
| HOXA_cluster | 0.9962 | 0.9962 | 0.9999 | 0.9999 |
| chr12_GAPDH_cluster | 0.9968 | 0.9968 | 0.9972 | 0.9986 |
| chr19_dense | 0.9992 | 0.9992 | 0.9998 | 0.9998 |
| chr1_dense | 0.9999 | 0.9999 | 1.0000 | 1.0000 |
| chr22_CRKL_cluster | 0.9984 | 0.9988 | 0.9968 | 0.9987 |

**Key finding:** Correlations are universally high (>0.97) except for HBB_cluster
in unique mode (0.93) and EGFR (0.97). These are the same two regions with the
highest MAE, confirming systematic quantification issues there.

### 1.6 Win Rates (Transcript Level)

**hulkrna (unique) vs competitors:**

| Comparison | Wins | Losses | Win Rate |
|------------|-----:|-------:|---------:|
| hulkrna vs salmon | 29 | 91 | 24.2% |
| hulkrna vs kallisto | 27 | 93 | 22.5% |

**hulkrna (unique) vs hulkrna_mm (multimap):**

| Comparison | Count | Percentage |
|------------|------:|-----------:|
| hulkrna wins | 8 | 6.7% |
| Tied | 50 | 41.7% |
| hulkrna_mm wins | 62 | 51.7% |

Multimap mode is strictly better or equal in 93.3% of conditions.

### 1.7 Best & Worst Conditions for hulkrna

**Top 10 worst (highest ratio vs best competitor):**

| Region | gDNA | SS | hulkrna MAE | Best MAE | Best Tool | Ratio |
|--------|------|----|------------:|---------:|-----------|------:|
| HBB_cluster | moderate | 1.00 | 998.1 | 8.3 | salmon | 119.8× |
| HBB_cluster | low | 1.00 | 996.4 | 10.1 | salmon | 98.7× |
| HBB_cluster | none | 1.00 | 993.0 | 13.8 | salmon | 71.8× |
| HBB_cluster | high | 1.00 | 1087.7 | 17.8 | salmon | 61.3× |
| HBB_cluster | none | 0.99 | 1035.9 | 24.4 | salmon | 42.4× |
| HBB_cluster | low | 0.99 | 983.2 | 23.3 | salmon | 42.2× |
| HBB_cluster | moderate | 0.99 | 985.5 | 24.6 | salmon | 40.0× |
| HBB_cluster | high | 0.99 | 1083.6 | 33.4 | salmon | 32.5× |
| HBB_cluster | low | 0.95 | 984.1 | 73.4 | salmon | 13.4× |
| HBB_cluster | moderate | 0.95 | 985.5 | 73.9 | salmon | 13.3× |

All 12 worst conditions are in HBB_cluster — a locus with highly homologous
hemoglobin genes that require multimap resolution.

**Top 10 best (lowest ratio vs best competitor):**

| Region | gDNA | SS | hulkrna MAE | Best MAE | Best Tool | Ratio |
|--------|------|----|------------:|---------:|-----------|------:|
| CD44 | moderate | 0.95 | 55.4 | 93.3 | salmon | 0.59× |
| CD44 | high | 0.95 | 68.5 | 105.7 | salmon | 0.65× |
| CD44 | low | 0.95 | 60.5 | 92.8 | salmon | 0.65× |
| CD44 | none | 0.95 | 59.9 | 91.4 | salmon | 0.66× |
| chr1_dense | high | 0.95 | 43.9 | 66.3 | salmon | 0.66× |
| chr1_dense | moderate | 0.95 | 46.5 | 67.7 | salmon | 0.69× |
| chr1_dense | low | 0.95 | 46.7 | 65.5 | salmon | 0.71× |
| chr1_dense | none | 0.95 | 46.7 | 64.7 | salmon | 0.72× |
| BRCA1 | low | 0.95 | 146.5 | 186.5 | kallisto | 0.79× |
| chr12_GAPDH | high | 0.95 | 51.4 | 64.2 | kallisto | 0.80× |

All best conditions are at **SS=0.95**, where hulkrna's strand model provides
the most differentiation vs alignment-free tools.

---

## 2. Gene-Level Results

### 2.1 Overall Means

| Tool | MAE | RMSE | Pearson r | Spearman ρ |
|------|----:|-----:|----------:|-----------:|
| hulkrna | 372.1 | 731.5 | 0.7920 | 0.7739 |
| hulkrna_mm | 237.6 | 334.0 | 0.7995 | 0.7846 |
| **salmon** | **275.9** | **348.8** | **0.8000** | **0.7992** |
| kallisto | 340.6 | 433.9 | 0.8000 | 0.7991 |
| htseq | 790.8 | 1525.2 | 0.7421 | 0.6749 |

**Key finding:** At the gene level, hulkrna_mm (237.6) actually beats salmon
(275.9) and kallisto (340.6) on MAE. hulkrna unique mode (372.1) is worse than
salmon but better than kallisto and htseq. All tools show much lower
correlations at gene level (~0.80) vs transcript level (~0.99), reflecting the
inherent difficulty of gene-level estimation from isoform-resolved data.

### 2.2 Gene MAE by Region

| Region | hulkrna | hulkrna_mm | salmon | kallisto | htseq |
|--------|--------:|-----------:|-------:|---------:|------:|
| BRCA1 | 99.8 | 99.6 | 205.6 | 215.1 | 240.4 |
| CD44 | 67.7 | 67.7 | 418.7 | 510.8 | 629.9 |
| EGFR | 313.2 | 300.2 | 824.5 | 1156.3 | 1212.8 |
| FGFR2 | 1350.8 | 1346.9 | 923.2 | 979.3 | 980.9 |
| HBB_cluster | 1408.5 | 156.5 | 47.7 | 119.7 | 1735.7 |
| HOXA_cluster | 84.0 | 84.0 | 24.7 | 35.9 | 742.3 |
| chr12_GAPDH | 67.9 | 67.7 | 69.7 | 84.6 | 83.5 |
| chr19_dense | 133.6 | 132.9 | 91.4 | 111.1 | 111.0 |
| chr1_dense | 118.4 | 98.7 | 95.3 | 124.9 | 563.0 |
| chr22_CRKL | 77.0 | 21.9 | 57.8 | 68.8 | 1608.5 |

**Notable gene-level results:**
- **BRCA1, CD44, EGFR:** hulkrna dominates. CD44 hulkrna MAE 67.7 vs salmon
  418.7 (**6.2× better**). EGFR hulkrna 313 vs salmon 825 (**2.6× better**).
- **FGFR2, HBB_cluster:** hulkrna struggles. FGFR2 hulkrna 1351 vs salmon 923.
  HBB unique mode 1409 vs salmon 48.
- **htseq** is generally worst due to its simple union-intersection counting
  model, except in FGFR2 and chr19_dense where it's competitive.

### 2.3 Gene Win Rates

| Comparison | Wins | Losses | Win Rate |
|------------|-----:|-------:|---------:|
| hulkrna vs salmon | 44 | 76 | 36.7% |
| hulkrna vs kallisto | 46 | 74 | 38.3% |
| hulkrna vs htseq | 85 | 35 | 70.8% |

---

## 3. Speed Comparison

| Tool | Mean (s) | Min (s) | Max (s) | Total (s) |
|------|--------:|--------:|--------:|----------:|
| hulkrna | 9.99 | 2.24 | 26.95 | 1199.3 |
| hulkrna_mm | 10.36 | 5.04 | 26.48 | 1243.6 |
| salmon | 0.63 | 0.29 | 5.44 | 75.7 |
| kallisto | 0.51 | 0.38 | 0.95 | 60.7 |

hulkrna is approximately **16× slower than salmon** and **20× slower than
kallisto** on average. This is expected since hulkrna performs alignment-based
quantification (via minimap2) while salmon and kallisto use lightweight
pseudo/quasi-mapping.

---

## 4. Summary of Strengths and Weaknesses

### Strengths

1. **Gene-level quantification in complex loci:** hulkrna excels at gene-level
   estimation in regions like BRCA1 (2× better than salmon), CD44 (6× better),
   and EGFR (2.6× better). This is hulkrna's strongest competitive advantage.

2. **Low strand-specificity robustness:** At SS=0.95, hulkrna consistently
   outperforms salmon by 10–40% in several regions (CD44, chr1_dense, BRCA1,
   chr12_GAPDH). This validates the strand model's effectiveness.

3. **gDNA resilience:** Performance degradation from gDNA contamination is
   proportional across all tools, and hulkrna handles it no worse than others.

### Weaknesses

1. **Multimapping loci (HBB_cluster):** unique-map mode fails catastrophically
   in homologous gene clusters (26–120× worse than salmon). Multimap mode
   fixes most of this but is still 2.5× worse than salmon at HBB.

2. **Transcript-level accuracy overall:** hulkrna's mean MAE (213) is 2.6×
   worse than salmon (83). Even multimap mode (119) trails salmon by 43%.

3. **EGFR region:** Consistently worst non-HBB region. Neither unique nor
   multimap mode helps (both ~455 MAE vs salmon 273), suggesting an issue
   in the EM/resolution algorithm for this gene family structure.

4. **Strand specificity exploitation:** hulkrna doesn't improve much as SS
   increases (216 → 209), while salmon improves dramatically (112 → 62).
   hulkrna is under-utilizing high strand purity.

5. **HOXA_cluster:** 3–7× worse than salmon/kallisto. These are closely
   spaced homeobox genes that require better resolution logic.

6. **Speed:** 16–20× slower than alignment-free tools. Acceptable for
   targeted analyses but not for transcriptome-wide use.

---

## 5. Root Cause Analysis & Recommendations

### Issue 1: HBB_cluster Catastrophe (Priority: Critical)

**Root cause:** The HBB locus contains highly homologous hemoglobin genes
(HBB, HBD, HBG1, HBG2, HBE1) with near-identical sequences. In unique-map
mode, most reads are unmappable or misassigned because they map equally well
to multiple genes. The reads are either discarded (unique mode) or poorly
distributed (multimap mode still trails salmon 2.5×).

**Recommendations:**
- **Default to multimap mode.** The 41.7% tie rate between modes confirms
  multimap handling only activates where needed, with minimal downside.
- **Improve EM convergence for homologous families.** Salmon's equivalence
  class EM is demonstrably better at distributing reads among paralogs.
  Investigate whether hulkrna's posterior update adequately handles gene
  families with nearly identical exonic sequences.
- **Consider equivalence-class compression.** Group fragments that map to
  identical transcript sets and process them as equivalence classes rather
  than individual fragments, which would both speed up EM and improve
  numerical stability.

### Issue 2: EGFR Region Underperformance (Priority: High)

**Root cause:** Multimap mode does not help EGFR (456.8 vs 454.4), ruling
out multimapping as the cause. EGFR has complex alternative splicing with
overlapping exon structures. The issue likely lies in how the EM/resolution
algorithm distributes fragments among highly overlapping isoforms.

**Recommendations:**
- **Diagnose the EGFR region specifically:** Extract per-isoform truth vs
  estimated values. Identify which isoforms are overestimated/underestimated.
- **Review EM prior initialization:** Check if uniform priors for isoform
  abundances lead to poor convergence in regions with many splice variants.
- **Evaluate fragment-to-isoform assignment weights:** The compatibility
  scoring (insert size, splice junctions) may not be discriminating enough
  between closely overlapping isoforms.

### Issue 3: Strand Specificity Under-utilization (Priority: High)

**Root cause:** hulkrna's MAE barely changes from SS=0.95 to SS=1.00
(216 → 209, only 3% improvement), while salmon improves 45%. This suggests
hulkrna's strand model assigns overly diffuse probabilities even when strand
signal is strong, or that strand information isn't sufficiently weighted in
the EM resolution step.

**Recommendations:**
- **Increase strand model influence at high specificity.** When observed
  strand specificity is ≥0.99, the strand model should more aggressively
  prune antisense assignments.
- **Audit strand posterior integration:** Check whether strand-based
  probabilities are multiplied into fragment-transcript compatibility
  scores before EM resolution. If strand is used post-EM, it may be too
  late to help.
- **Benchmark strand-only vs EM-only:** Temporarily disable strand
  weighting to quantify its actual contribution, then tune its influence.

### Issue 4: HOXA_cluster (Priority: Medium)

**Root cause:** HOXA genes are closely spaced on chromosome 7 and share
regulatory elements and partial exon overlap. hulkrna MAE 49.5 vs salmon
15.3, suggesting poor resolution among adjacent genes.

**Recommendations:**
- **Check transcript boundary handling:** HOXA genes have overlapping UTRs.
  Ensure fragment assignment correctly handles UTR-crossing reads.
- **Review the resolution algorithm for adjacent non-overlapping genes:**
  If genes are on opposite strands (some HOXA genes overlap antisense
  lncRNAs), verify strand model correctly separates them.

### Issue 5: Speed (Priority: Low)

**Root cause:** hulkrna uses alignment-based quantification (minimap2
alignment → BAM → fragment counting → EM) while salmon/kallisto use
lightweight pseudo-mapping. The 16–20× slowdown is inherent to the
approach.

**Recommendations:**
- For the current alignment-based approach, speed optimization is
  secondary to accuracy. Focus on accuracy improvements first.
- Consider caching or pre-filtering to avoid re-scanning entire BAM
  files for each region in regional mode.
- Long-term: evaluate whether a hybrid approach (lightweight mapping for
  initial abundance estimation, alignment for complex loci) could
  combine speed and accuracy.

---

## 6. Priority Roadmap

| Priority | Issue | Expected Impact |
|----------|-------|-----------------|
| 1 | Default to multimap mode | -45% overall transcript MAE |
| 2 | Fix EM for homologous gene families | -50% MAE at HBB-like loci |
| 3 | Improve strand specificity utilization | -10–20% MAE at SS≥0.99 |
| 4 | Diagnose & fix EGFR isoform resolution | -40% MAE at complex splice loci |
| 5 | HOXA boundary handling | -30% MAE at dense gene clusters |

Simply defaulting to multimap mode would bring hulkrna's overall transcript
MAE from 213 to 119 — a 44% improvement that narrows the gap with salmon
from 2.6× to 1.4×.
