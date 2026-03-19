# Benchmark v6 Analysis: Pruning Impact × Oracle/Minimap2 vs Salmon/Kallisto

**Date:** 2026-03-20  
**Config:** [scripts/benchmark_v6_prune.yaml](../../scripts/benchmark_v6_prune.yaml)  
**Data:** HeLa/Cervix simulation, 10M RNA fragments/condition, 9 conditions (3 gDNA × 3 SS)  
**Analysis script:** [scripts/debug/benchmark_v6_analysis.py](../../scripts/debug/benchmark_v6_analysis.py)  
**Runtime:** 23,824s (~6.6 hours)

**Tools compared (6 total):**

| Tool | Description |
|------|-------------|
| `rigel_default_oracle` | Rigel with pruning ON (threshold=0.1), oracle alignments |
| `rigel_no_prune_oracle` | Rigel with pruning OFF, oracle alignments |
| `rigel_default_minimap2` | Rigel with pruning ON, minimap2 alignments |
| `rigel_no_prune_minimap2` | Rigel with pruning OFF, minimap2 alignments |
| `salmon` | Salmon (pseudo-alignment) |
| `kallisto` | Kallisto (pseudo-alignment) |

---

## Executive Summary

### Pruning Has Negligible Impact on Accuracy

**Post-EM pruning does not meaningfully affect results.** Disabling pruning changes the average transcript-level MAE by just +0.34% for oracle alignments and −0.00% for minimap2. The false negative rate — which we hypothesized pruning would improve — drops from 1,390 to 1,388 FNs with oracle (2 fewer out of ~14,600 highly expressed transcripts). Only **21 transcripts** were "rescued" by disabling pruning across all 118K+ transcripts, all in multi-isoform genes with median 7 isoforms. The false negatives are overwhelmingly an **identifiability problem**, not a pruning artifact.

### Rigel Dominates at Both Alignment Levels

**With oracle alignments**, rigel achieves 10.7× lower MAE than salmon and 16.8× lower than kallisto. **With realistic minimap2 alignments**, rigel still achieves 4.2× lower MAE than salmon and 6.6× lower than kallisto — establishing clear real-world superiority. Gene-level accuracy shows even larger advantages: 7.4× over salmon and 12.1× over kallisto with minimap2.

### The Oracle-to-Minimap2 Gap Is the Main Opportunity

The gap from oracle (MAE=0.95) to minimap2 (MAE=2.41) — a 2.5× degradation — is entirely due to alignment quality, not the EM algorithm or pruning. Closing this gap through alignment improvements or MAPQ-based weighting represents the largest opportunity for real-world accuracy gains.

---

## 1. Pruning Mechanism and Investigation

### How Pruning Works

Rigel's post-EM pruning (in `em_solver.cpp`) fires after SQUAREM convergence. For each EM component (except gDNA, which is never pruned), it prunes the component if **all** of:

1. `prune_threshold >= 0` (pruning is enabled)
2. `unambig_totals[i] == 0` (no fragments uniquely assigned to this component)
3. `evidence_ratio = (unambig + em_totals) / alpha < prune_threshold` (low EM evidence relative to prior)

After pruning, a single redistribution step is performed (VBEM or MAP), but **not** a full re-EM. This is a "hard" decision that zeros out pruned components permanently.

### Hypothesis

We hypothesized that pruning was responsible for rigel's 9.0% false negative rate (vs salmon's 6.7%) discovered in benchmark v5. The reasoning: when two isoforms share >90% of their sequence and neither has unique fragments, the EM may converge with all mass on one isoform, and pruning then zeros the other permanently, producing a false negative.

### Result: Hypothesis Largely Refuted

Disabling pruning yields negligible FN improvement. The FNs are caused by the **EM convergence** itself — the pruning step merely confirms what the EM already decided.

---

## 2. Pruning Impact: Detailed Results

### Transcript-Level MAE: Pruning ON vs OFF

**Oracle alignments:**

| Condition | Prune ON | Prune OFF | Δ% |
|-----------|----------|-----------|-----|
| gdna_none, ss=0.50 | 0.9448 | 0.9482 | +0.36% |
| gdna_none, ss=0.90 | 0.9319 | 0.9322 | +0.03% |
| gdna_none, ss=1.00 | 0.9162 | 0.9205 | +0.47% |
| gdna_low, ss=0.50 | 0.9324 | 0.9310 | −0.15% |
| gdna_low, ss=0.90 | 0.9447 | 0.9365 | −0.86% |
| gdna_low, ss=1.00 | 0.9505 | 0.9427 | −0.82% |
| gdna_high, ss=0.50 | 0.9861 | 0.9956 | +0.97% |
| gdna_high, ss=0.90 | 0.9414 | 0.9500 | +0.91% |
| gdna_high, ss=1.00 | 0.9836 | 1.0036 | +2.03% |
| **AVERAGE** | **0.9479** | **0.9511** | **+0.34%** |

**Minimap2 alignments:**

| Metric | Prune ON | Prune OFF | Δ% |
|--------|----------|-----------|-----|
| Avg MAE | 2.4087 | 2.4086 | −0.00% |
| Pearson | 0.985587 | 0.985590 | negligible |

**Gene-level:**

| Aligner | Prune ON | Prune OFF | Δ% |
|---------|----------|-----------|-----|
| Oracle | 0.4063 | 0.4107 | +1.08% |
| Minimap2 | 2.0619 | 2.0619 | −0.00% |

**Pruning has near-zero effect on minimap2** because the alignment noise dominates any post-EM fine-tuning. With oracle alignments, pruning provides a slight benefit at high gDNA levels (+2% lower MAE) but slightly hurts at low gDNA (−0.8%).

### False Negative Rates

| Aligner | Tool | Avg FP | Avg FN (truth>10) |
|---------|------|--------|-------------------|
| Oracle | prune ON | 1,091 | **1,390** |
| Oracle | prune OFF | 1,092 | **1,388** |
| Minimap2 | prune ON | 3,978 | **1,878** |
| Minimap2 | prune OFF | 3,983 | **1,864** |

Disabling pruning rescues only **2 FNs** with oracle and **14 FNs** with minimap2 — out of ~14,600 highly expressed transcripts.

### Rescued Transcripts (Oracle, gdna_low ss=0.90)

Only 21 transcripts were rescued (FN with pruning → detected without). All are in multi-isoform genes (median 7 isoforms per gene, no single-isoform rescues):

| Transcript | Gene | Truth | Prune ON | Prune OFF |
|------------|------|-------|----------|-----------|
| ENST00000646319 | DDX3X | 531 | 0.0 | 549.8 |
| ENST00000261772 | AARS1 | 287 | 0.0 | 269.7 |
| ENST00000444746 | RIF1 | 226 | 0.0 | 209.7 |
| ENST00000389266 | GARS1 | 197 | 0.0 | 218.6 |
| ENST00000356401 | ARIH2 | 155 | 0.0 | 113.9 |
| ENST00000267973 | SKIC8 | 144 | 0.0 | 145.2 |
| ... (15 more) | | | | |

These rescues are real and meaningful (DDX3X: truth=531 recovered to 549.8), but they're rare: 21 out of 1,325 FNs (1.6%).

### Per-Transcript Win Rate (Oracle, prune ON vs OFF)

| | Count | % |
|---|-------|---|
| Prune ON wins | 5,469 | 17.2% |
| Prune OFF wins | 15,547 | 48.9% |
| Ties | 10,789 | 33.9% |

**Paradox:** Pruning OFF "wins" on nearly 3× more transcripts, yet its average MAE is slightly *worse*. Explanation: disabling pruning produces many tiny improvements (rounding-level changes on thousands of transcripts) but also occasional large regressions (particularly at high gDNA), resulting in a net wash or slight degradation.

---

## 3. Realistic Comparison: Minimap2 Rigel vs Salmon vs Kallisto

This section presents the head-to-head comparison using **minimap2 alignments** for rigel — the realistic scenario for real RNA-seq data where perfect alignments are unavailable.

### Overall Averages (9 conditions)

| Tool | MAE | RMSE | Pearson | Rigel advantage |
|------|-----|------|---------|-----------------|
| **rigel_minimap2** | **2.41** | **28.06** | **0.9856** | — |
| salmon | 10.19 | 76.73 | 0.9540 | **4.2× better** |
| kallisto | 15.89 | 97.47 | 0.9302 | **6.6× better** |

### By gDNA Contamination

| gDNA | rigel_mm2 MAE | salmon MAE | kallisto MAE | Rigel advantage |
|------|---------------|------------|--------------|-----------------|
| none (0%) | 2.26 | 8.89 | 14.41 | 3.9× |
| low (20%) | 2.42 | 10.00 | 15.68 | 4.1× |
| high (50%) | 2.54 | 11.69 | 17.59 | 4.6× |

Rigel's advantage **increases** with gDNA contamination because the tripartite model absorbs gDNA reads correctly while salmon/kallisto count them as expression.

### By Strand Specificity

| SS | rigel_mm2 MAE | salmon MAE | kallisto MAE |
|----|---------------|------------|--------------|
| 0.50 (unstranded) | 2.67 | 12.34 | 18.38 |
| 0.90 | 2.28 | 8.81 | 13.92 |
| 1.00 (fully stranded) | 2.28 | 9.43 | 15.38 |

Rigel benefits slightly from strand information (MAE drops from 2.67 to 2.28 going from unstranded to stranded), but the improvement is modest. The minimap2 alignment noise dominates.

### Gene-Level Accuracy (Realistic)

| Tool | MAE | RMSE | Pearson |
|------|-----|------|---------|
| **rigel_minimap2** | **2.06** | **29.20** | **0.9893** |
| salmon | 15.35 | 111.58 | 0.9568 |
| kallisto | 25.03 | 157.17 | 0.9337 |

Rigel is **7.4× better** than salmon and **12.1× better** than kallisto at gene level with minimap2 alignments.

### Per-Transcript Win Rate (Realistic)

| Comparison | rigel_mm2 wins | Other wins |
|-----------|----------------|------------|
| rigel_mm2 vs salmon | 50.0% | 31.2% |
| rigel_mm2 vs kallisto | 64.3% | 22.9% |

Even with realistic alignments, rigel wins the majority of per-transcript comparisons.

---

## 4. Oracle Results: Upper Bound Accuracy

### Overall Averages (9 conditions)

| Tool | MAE | RMSE | Pearson |
|------|-----|------|---------|
| **rigel_oracle** | **0.95** | **7.70** | **0.9989** |
| salmon | 10.19 | 76.73 | 0.9540 |
| kallisto | 15.89 | 97.47 | 0.9302 |

Oracle results are **10.7× better** than salmon, demonstrating the theoretical ceiling of rigel's statistical model.

### Condition Stability

The oracle MAE is remarkably stable across all 9 conditions: 0.92–0.99 (7% variation), compared to salmon (7.79–14.46, 86% variation) and kallisto (12.80–20.87, 63% variation). Rigel's generative model is well-calibrated regardless of gDNA contamination level or strand specificity.

### Gene-Level

| Tool | MAE | Rigel advantage |
|------|-----|-----------------|
| rigel_oracle | 0.41 | — |
| salmon | 15.35 | **37.4×** |
| kallisto | 25.03 | **61.0×** |

---

## 5. Expression-Stratified Accuracy

### Oracle (gdna_low, ss=0.90)

| Stratum | rigel_oracle | salmon | kallisto | Rigel advantage |
|---------|-------------|--------|----------|-----------------|
| Low (0, Q25] | 1.50 | 10.33 | 15.72 | 6.9× |
| Mid (Q25, Q75] | 6.06 | 21.24 | 29.47 | 3.5× |
| High (>Q75) | 14.24 | 48.65 | 81.02 | 3.4× |
| Zero truth | 0.18 | 2.62 | 4.74 | **14.6×** |

### Minimap2 (gdna_low, ss=0.90)

| Stratum | rigel_mm2 MAE |
|---------|--------------|
| Low (0, Q25] | 2.90 |
| Mid (Q25, Q75] | 7.43 |
| High (>Q75) | 33.32 |
| Zero truth | 0.70 |

Minimap2 rigel is consistently better than salmon across all strata (salmon MAE: 10.33/21.24/48.65/2.62).

---

## 6. Error Distribution

### Oracle (gdna_low, ss=0.90)

| Metric | rigel_oracle | salmon | kallisto |
|--------|-------------|--------|----------|
| P50 abs error | 2.55 | 5.00 | 9.31 |
| P90 abs error | 16.05 | 47.72 | 74.67 |
| P95 abs error | 26.17 | 88.00 | 134.01 |
| P99 abs error | 59.32 | 292.45 | 465.72 |
| Mean signed error | −0.69 | +13.18 | +25.25 |
| FP rate (truth=0) | **1.4%** | 22.3% | 39.8% |
| FN rate (truth>10) | 9.1% | **6.7%** | **5.5%** |

### Minimap2

| Metric | rigel_mm2 |
|--------|----------|
| P50 abs error | 3.06 |
| P90 abs error | 23.53 |
| P95 abs error | 40.66 |
| P99 abs error | 137.06 |
| Mean signed error | −6.28 |
| FP rate | 4.6% |
| FN rate | 10.8% |

**Key trade-off:** Rigel has the lowest false positive rate (1.4% oracle, 4.6% mm2) by a wide margin, but the highest false negative rate (9.1% oracle, 10.8% mm2). This is inherited from the EM's identifiability problem — isoforms sharing sequence collapse to one winner — not from pruning.

---

## 7. nRNA Impact on Accuracy

(Condition: gdna_none, ss=0.90, oracle)

| nRNA:mRNA Ratio | N | rigel_oracle MAE | salmon MAE | kallisto MAE |
|-----------------|---|------------------|------------|--------------|
| 0 (no nRNA) | 4,352 | 1.74 | 20.52 | 23.04 |
| 0 < ratio ≤ 1 | 6,385 | 7.53 | 28.45 | 37.58 |
| 1 < ratio ≤ 5 | 10,623 | 7.99 | 25.94 | 40.31 |
| ratio > 5 | 10,445 | 6.94 | 18.60 | 38.98 |

Rigel outperforms at all nRNA strata. With minimap2, the rigel MAE ranges from 6.35 (no nRNA) to 15.78 (low nRNA), still substantially better than salmon/kallisto.

---

## 8. Isoform Resolution

(Condition: gdna_none, ss=0.90)

### Oracle

| Isoform Count | rigel_oracle MAE | salmon MAE | kallisto MAE |
|---------------|------------------|------------|--------------|
| 1 isoform | 2.14 | 21.56 | 31.10 |
| 2 isoforms | 6.28 | 20.60 | 38.26 |
| 3-5 isoforms | 7.61 | 20.93 | 33.03 |
| 6-10 isoforms | 8.37 | 26.76 | 44.08 |
| >10 isoforms | 8.36 | 31.65 | 44.66 |

### Minimap2

| Isoform Count | rigel_mm2 MAE |
|---------------|--------------|
| 1 isoform | 8.30 |
| 2 isoforms | 11.87 |
| 3-5 isoforms | 11.83 |
| 6-10 isoforms | 15.37 |
| >10 isoforms | 16.80 |

Rigel minimap2 is still better than salmon for all isoform strata (salmon: 21.56–31.65), and notably better than kallisto (31.10–44.66). The minimap2 gap grows with isoform count because alignment ambiguity compounds with isoform ambiguity.

---

## 9. Pool-Level mRNA Totals

Rigel's pool-level mRNA estimates remain excellent:

| Aligner | Best Condition | Worst Condition |
|---------|---------------|-----------------|
| Oracle | +0.03% (gdna_high, ss=0.90) | −0.71% (gdna_high, ss=1.00) |
| Minimap2 | −1.00% (gdna_high, ss=0.90) | −3.02% (gdna_high, ss=0.50) |

The tripartite model correctly partitions fragments regardless of pruning. Pruning ON vs OFF produces identical pool-level totals (within 0.01%).

---

## 10. Performance

| Tool | Avg Time (s) | Avg RSS (MB) | Throughput |
|------|-------------|-------------|------------|
| rigel_oracle | 416 | 21,616 | 64K frags/s |
| rigel_minimap2 | 480 | 22,745 | 66K frags/s |
| salmon | 215 | — | 57K frags/s |
| kallisto | 88 | — | 141K frags/s |

Pruning has zero impact on runtime (±0.5%). Rigel is ~2× slower than salmon and ~5× slower than kallisto, with ~22 GB RSS.

---

## 11. Competitive Ranking

### Overall (all 6 tools)

| Rank | Tool | Avg Rank | Wins |
|------|------|----------|------|
| 1 | rigel_default_oracle | 1.33 | 6/9 conditions |
| 2 | rigel_no_prune_oracle | 1.67 | 3/9 conditions |
| 3 | rigel_default_minimap2 | 3.44 | 0 |
| 4 | rigel_no_prune_minimap2 | 3.56 | 0 |
| 5 | salmon | 5.00 | 0 |
| 6 | kallisto | 6.00 | 0 |

### Realistic (minimap2 + salmon + kallisto)

| Rank | Tool | Avg Rank |
|------|------|----------|
| 1 | rigel_default_minimap2 | 1.44 |
| 2 | rigel_no_prune_minimap2 | 1.56 |
| 3 | salmon | 3.00 |
| 4 | kallisto | 4.00 |

Rigel wins **every condition** even with minimap2 alignments.

---

## 12. Summary of Findings

### Finding 1: Pruning Is Not the Problem

The post-EM pruning mechanism has negligible impact on accuracy. The 9% false negative rate is driven by **EM convergence identifiability** — when isoforms share >90% of their sequence, the EM assigns all mass to one winner before pruning ever fires. Disabling pruning rescues only 21/1,325 FNs (1.6%) with oracle alignments.

**Recommendation:** Keep pruning enabled (default threshold=0.1). It's benign and provides a minor benefit at high gDNA levels. Future FN reduction efforts should focus on the EM initialization or regularization, not post-EM pruning.

### Finding 2: Rigel Dominates Across All Scenarios

| Comparison | MAE Ratio |
|-----------|-----------|
| rigel_oracle vs salmon | 10.7× better |
| rigel_oracle vs kallisto | 16.8× better |
| rigel_minimap2 vs salmon | 4.2× better |
| rigel_minimap2 vs kallisto | 6.6× better |
| rigel_oracle gene vs salmon gene | 37.4× better |
| rigel_minimap2 gene vs salmon gene | 7.4× better |

### Finding 3: The Oracle-to-Minimap2 Gap Is the Key Opportunity

| Metric | Oracle | Minimap2 | Degradation |
|--------|--------|----------|-------------|
| MAE | 0.95 | 2.41 | 2.5× |
| Pearson | 0.9989 | 0.9856 | — |
| FN rate | 9.1% | 10.8% | +1.7pp |
| FP rate | 1.3% | 4.6% | +3.3pp |

The 2.5× gap is entirely alignment-driven. Improving alignment quality, MAPQ filtering, or minimap2 parameter tuning would directly translate to accuracy gains. This is the single largest opportunity for improving rigel's real-world performance.

### Finding 4: Rigel's Trade-off Profile Is Favorable

Rigel trades higher FN rate (9%) for dramatically lower FP rate (1.4% vs 22% salmon vs 40% kallisto) and 10× lower MAE. This trade-off is favorable for most applications because:

- False positives (phantom expression) contaminate downstream differential expression analysis
- False negatives in rigel are concentrated in multi-isoform genes where the *gene-level* count is still correct (MAE 0.41)
- Users can detect FN risk by checking isoform count and unique fragment availability

---

## 13. Opportunities for Improvement (Updated from v5)

| # | Opportunity | Impact | Difficulty | Status |
|---|-----------|--------|------------|--------|
| 1 | ~~Disable pruning to reduce FNs~~ | ~~Moderate~~ | ~~Low~~ | **Refuted** — pruning has negligible effect |
| 2 | **Close minimap2 gap** via alignment filtering/MAPQ weighting | High (2.5× MAE improvement potential) | Medium-High |
| 3 | **EM isoform regularization** — prevent full zeroing of near-identical isoforms within the EM | Moderate (targets 1,325 FNs) | Medium |
| 4 | **Speed optimization** — reduce 416s runtime | High (user experience) | High |
| 5 | **Memory reduction** from 22 GB | Moderate (broader accessibility) | Medium |
| 6 | **Report isoform uncertainty** when near-identical isoforms can't be disambiguated | Low (transparency) | Low |
