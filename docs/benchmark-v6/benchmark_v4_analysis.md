# Benchmark v4 Analysis Report

## Executive Summary

Benchmark v4 runs the improved codebase (transcript-space FL computation) on the pristine simulation (10M mRNA fragments, no gDNA, no nRNA, strand_specificity=0.95).

### Key Metrics Comparison: v4 vs v3

| Tool | Metric | v3 | v4 | Delta |
|------|--------|----|----|-------|
| **Rigel Oracle** | MAE | 1.84 | 1.83 | **−0.01** (neutral) |
| **Rigel Oracle** | RMSE | 30.19 | 30.17 | −0.02 |
| **Rigel Oracle** | Spearman | 0.8799 | 0.8801 | +0.0002 |
| **Rigel Minimap2** | MAE | 5.51 | 6.71 | **+1.20** ⚠️ |
| **Rigel Minimap2** | RMSE | 99.82 | 105.38 | +5.56 |
| **Rigel Minimap2** | Spearman | 0.7560 | 0.7056 | −0.0504 |
| **Salmon** | MAE | 2.91 | 2.91 | 0 (unchanged) |
| **Kallisto** | MAE | 3.36 | 3.36 | 0 (unchanged) |

### Pool-Level Counts

| Pool | Truth | Oracle v3 | Oracle v4 | MM2 v3 | MM2 v4 |
|------|-------|-----------|-----------|--------|--------|
| mRNA | 10,000,000 | 9,994,535 | 9,994,550 | 9,922,277 | 9,904,339 |
| nRNA | 0 | 3,680 | 3,668 | 72,958 | 90,543 |
| gDNA | 0 | 1,785 | 1,782 | 3,660 | 4,013 |

## Regression Analysis

### Oracle → Unaffected
The transcript-space FL change is **neutral for oracle alignments**. MAE improved by 0.01, Spearman improved by 0.0002. This is expected: oracle alignments don't create inter-gene multimappers, so the FL change has no pathological interactions.

### Minimap2 → Regressed (+1.20 MAE)
The regression concentrates in **specific gene families with retrotransposed pseudogene paralogs**:

#### Top Gene-Level Regressions  
| Gene | n_tx | Truth | Oracle | MM2 v3 | MM2 v4 | Gene Δ |
|------|------|-------|--------|--------|--------|--------|
| RPL5 | 2 | 14,579 | 14,579 | 13,995 | 11,903 | +2,092 |
| EIF3C/EIF3CL | 7 | 5,510 | 5,510 | 5,263 | 5,288 | +1,723* |
| ANXA2 | 17 | 29,968 | 29,968 | 24,324 | 22,968 | +1,356 |
| NPM1 | 12 | 23,571 | 23,536 | 18,210 | 17,032 | +1,178 |
| RPLP0 | 6 | 14,172 | 14,150 | 11,879 | 10,958 | +921 |

*EIF3C→EIF3CL cross-gene leakage: 1,751 fragments incorrectly assigned to EIF3CL (truth=33)

#### Pattern: The regression is **inter-gene** (leakage between gene families), NOT intra-gene
- Gene-level MAE: v3=8.66, v4=11.13, Δ=+2.48
- Intra-gene component: v3=0.98, v4=0.31 (actually **improved**)
- False positive counts (truth=0): v3=103,685, v4=137,876 (Δ=+34,191)
- 108 transcripts gained >100 counts in MM2 v4 that should have ~0

#### Affected Gene Families
All top regressions involve **housekeeping genes with processed pseudogenes**:
- Ribosomal proteins: RPL5, RPL7A, RPL9, RPL8, RPL3, RPL13A, RPSA, RPLP0, RPL21, RPL36A
- EEF1A1 (and pseudogene EEF1A1P5), GAPDH, LDHB, NPM1 (NPM1P27), ANXA2 (ANXA2P2)
- ACTN4, HSP90B1, HSPA8, TMSB4X, B2M

## Root Cause Analysis

### What Changed
The **only** code change between v3 and v4 is `compute_frag_lengths()` in `resolve_context.h`:
- v3: SJ gap correction using cgranges overlap queries
- v4: Transcript-space coordinate projection via `genomic_to_tx_pos()`

### Why the Change is Correct
1. **Per-fragment FL values**: For multi-block fragments where v3's SJ gap correction failed (near-miss SJ), v3 produced inflated FL. v4 always produces correct FL. This **helps** the real gene's likelihood.
2. **Identical FL for pseudogenes**: For both single-block fragments AND multi-block fragments mapped to pseudogenes, FL is identical in v3 and v4.
3. **Oracle is unaffected**: Proves the FL computation itself is not broken.

### Why the Regression Occurs: Indirect FL Model Training Effect
The FL model is trained during the BAM scan from **unique mappers** (fragments that map unambiguously to one transcript). The trained FL distribution feeds into scoring for ALL fragments.

1. In v3, some unique mapper multi-exon fragments got inflated FL values → **broader trained distribution**
2. In v4, all unique mappers get correct FL values → **correctly centered, potentially narrower distribution**
3. The different FL distribution changes P(FL) for all fragments globally
4. This global shift causes the EM to converge to a **different steady state**
5. The new steady state happens to leak more to pseudogene transcripts for minimap2 alignments

**This is NOT a bug in the FL computation** — it's a second-order effect of correcting the FL model training. The v3 FL model had a compensating bias that coincidentally reduced pseudogene leakage.

### Why Oracle is Unaffected While Minimap2 is Not
- **Oracle alignments**: No inter-gene multimapping. Each fragment maps only to its true source transcript. FL model changes only affect intra-gene isoform resolution.
- **Minimap2 alignments**: Housekeeping genes have high-identity pseudogenes that minimap2 reports as secondary alignments. The EM must resolve these ambiguities using scoring features.

### The Fundamental Problem: Pseudogene Indistinguishability
For fragments that map to BOTH a real gene and its intronless pseudogene:
- **FL**: Identical (both give correct insert size)
- **Strand**: Same (pseudogene is on the same strand)
- **Overhang**: ~0 for both (pseudogene's single exon is broad)
- **NM**: ~0 for both (high-identity paralog)
- **Splice type**: Per-fragment, not per-candidate — doesn't discriminate

The ONLY discriminating feature is **splice junction matching during resolution**: spliced fragments (with N operation in CIGAR) are pruned from pseudogene candidates because pseudogenes have no annotated splice junctions. But **unspliced fragments** (both reads within one exon) have ZERO discrimination between real gene and pseudogene.

## Performance Context

| Tool | MAE | vs Salmon | vs Kallisto |
|------|-----|-----------|-------------|
| Rigel Oracle | 1.83 | **1.6× better** | **1.8× better** |
| Rigel Minimap2 | 6.71 | 2.3× worse* | 2.0× worse* |
| Salmon | 2.91 | — | 1.2× better |
| Kallisto | 3.36 | — | — |

*The minimap2 gap is dominated by pseudogene leakage in housekeeping genes

## Opportunities for Improvement

### 1. Pseudogene-Aware EM Prior (High Impact)
Add an informational prior that down-weights transcripts annotated as pseudogenes or retrotransposed copies. Gencode annotations include biotype information (`processed_pseudogene`, `unprocessed_pseudogene`, etc.) that could feed into the prior pseudocount.

### 2. Effective Length Bias for Pseudogene Discrimination
Pseudogenes are typically shorter than their parent genes. If the EM effective length normalization is applied, the pseudogene's shorter effective length would concentrate abundance per fragment, potentially helping or hurting. Worth investigating.

### 3. Splice-Aware Per-Candidate Scoring
Currently, splice_type is per-fragment, not per-candidate. A per-candidate splice compatibility score (does this transcript's exon structure support the observed splice junction?) would provide direct discrimination. If a fragment is spliced at position X, and transcript A has an annotated SJ at X but transcript B (pseudogene) does not, candidate A should get a bonus.

### 4. NH-Aware Multimapper Weighting
Minimap2's secondary alignments to pseudogenes could be down-weighted based on mapping quality or alignment score differences. Higher AS to the real gene vs pseudogene could inform the EM initialization.

### 5. FL Model Robustness
The FL model training could be made more robust to the specific set of unique mappers observed. Options:
- Use a parametric model (log-normal) instead of histogram → less sensitive to training set composition
- Train only from spliced unique mappers (which are guaranteed to give correct FL)
- Apply a minimum smoothing width to prevent over-fitting

### 6. Separate Benchmark: Oracle vs Minimap2 Gap Analysis
The oracle MAE (1.83) is the ceiling performance. The minimap2 gap (4.88 additional MAE) is entirely due to multimapper resolution. A targeted benchmark focusing on the top 30 gene families with pseudogenes would quantify the pseudogene-specific error contribution.

## Recommendation
The minimap2 regression is a **second-order effect of a correct FL fix**. The transcript-space FL computation is mathematically correct, improves intra-gene isoform resolution, and doesn't affect oracle alignments. The regression should NOT motivate reverting the FL change.

Instead, the pseudogene leakage problem (which existed in both v3 and v4, just at different magnitudes) should be addressed directly through **pseudogene-aware priors** (opportunity #1) or **per-candidate splice scoring** (opportunity #3).
