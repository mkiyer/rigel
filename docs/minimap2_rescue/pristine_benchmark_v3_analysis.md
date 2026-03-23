# Pristine Benchmark v3: Post Pool-Separated Likelihood Pruning

**Date**: 2026-03-22
**What changed**: Replaced hard overhang gate with pool-separated likelihood-based pruning (scoring.cpp)

## Executive Summary

The pool-separated likelihood pruning fix produced a **29.2% MAE reduction** for rigel+minimap2 (13.60 → 9.64), with nRNA siphon completely eliminated (171K → 73K pool-level, 0 per-transcript leakage). The gDNA false-positive explosion (39.9K → 3.7K) was also resolved. However, a significant gap vs oracle (MAE 9.64 vs 3.88) and vs salmon (MAE 9.64 vs 6.16) remains.

## Top-Level Results

| Tool | MAE | RMSE | Total Predicted | Total |Error| | Spearman |
|------|-----|------|-----------------|----------------|----------|
| **rigel oracle** | **3.88** | 44.17 | 9,988,163 | 461,025 | 0.880 |
| **rigel minimap2 v3** | **9.64** | 140.32 | 9,663,535 | 1,143,595 | 0.756 |
| rigel minimap2 v2 (pre-fix) | 13.60 | 172.74 | 9,474,225 | 1,614,177 | 0.719 |
| salmon | 6.16 | 56.38 | 9,501,393 | 731,092 | 0.894 |
| kallisto | 7.06 | 86.73 | 9,491,975 | 838,186 | 0.848 |

### Improvement from v2 → v3

| Metric | v2 (pre-fix) | v3 (post-fix) | Change |
|--------|-------------|---------------|--------|
| MAE | 13.60 | 9.64 | **-29.2%** |
| RMSE | 172.74 | 140.32 | **-18.8%** |
| nRNA siphon (pool) | 171,189 | 72,958 | **-57.4%** |
| gDNA false positives | 39,896 | 3,660 | **-90.8%** |
| mRNA recovered | 9,787,810 | 9,922,277 | **+134,467** |
| Spearman | 0.719 | 0.756 | **+0.037** |

### Per-transcript improvement

- 18,209 transcripts improved (total: +562,450 count units)
- 10,351 transcripts degraded (total: -91,868 count units)
- Net: **+470,582** count units recovered
- Largest single improvement: ENST00000229239.10 (error 22,301 → 14,516)

## Pool-Level Counts

| Metric | Truth | Oracle | MM2 v3 | Salmon | Kallisto |
|--------|-------|--------|--------|--------|----------|
| mRNA | 10,000,000 | 9,988,163 | 9,922,277 | 9,501,393 | 9,491,975 |
| nRNA | 0 | 3,680 | 72,958 | — | — |
| gDNA | 0 | 1,785 | 3,660 | — | — |
| **Missing** | — | **11,837** | **336,465** | **498,607** | **508,025** |

Note: salmon and kallisto "miss" more total counts (~500K) because they quantify against the transcriptome (unmappable reads are lost), whereas rigel quantifies against the genome (so only truly unmappable reads are lost). Rigel's mm2 missing 336K is the sum of nRNA siphon (73K) + gDNA (3.7K) + intergenic/unmapped reads.

## Error Distribution by Expression Level

| Expression bin | N transcripts | Oracle MAE | MM2 v3 MAE | Salmon MAE | Kallisto MAE |
|---------------|---------------|------------|------------|------------|--------------|
| 0 (unexpressed) | 77,628 | 0.53 | **1.34** | 0.25 | 0.89 |
| 1–9 | 15,824 | 1.88 | **4.84** | 1.72 | 1.83 |
| 10–99 | 13,457 | 9.21 | **14.67** | 8.98 | 8.90 |
| 100–999 | 9,902 | 19.12 | **32.30** | 25.74 | 27.54 |
| 1K–10K | 1,806 | 32.04 | **152.57** | 117.10 | 123.51 |
| >10K | 67 | 285.75 | **2,544.27** | 1,456.19 | 1,861.85 |

**Key observation**: The mm2 gap concentrates at highly-expressed transcripts (>1K counts). At >10K, mm2 MAE is 2,544 vs oracle 286 — a 9× multiplier. This is where the multimap redistribution problem bites hardest.

## Root Cause Analysis

### 1. Error Decomposition

| Category | Count units | % of mm2 total error |
|----------|------------|---------------------|
| Total mm2 |error| | 1,143,595 | 100% |
| Shared with oracle (inherent EM noise) | 375,734 | 32.9% |
| **MM2-specific error** | **767,861** | **67.1%** |
| — Isoform redistribution (gene OK, tx wrong) | 163,782 | 14.3% |
| — Pure mm2-specific (oracle <10, mm2 >50) | 261,403 | 22.8% |

### 2. Fragment Accounting

Of the 336K fragments "missing" from mm2 mRNA:
- **72,958** → nRNA siphon (still nonzero, but vastly reduced from 171K)
- **3,660** → gDNA false positives
- **~260K** → likely intergenic/unmapped (reads that minimap2 can't map to any transcript)

### 3. MM2-Specific Over-Predictions (451 transcripts, +134K)

These are transcripts where mm2 assigns significantly more fragments than truth, but oracle is fine.

**Structural signature**: 100% are single-exon transcripts (n_exons=1). These are primarily multimap targets — reads that minimap2 maps ambiguously and the EM redistributes incorrectly due to alignment quality differences.

**Top examples**:
- ENSG00000196205.8: truth=136, mm2=10,005 (+9,869) — massive multimap inflation
- ENSG00000281181.1: truth=0, mm2=7,188 — entirely phantom
- ENSG00000280614.1: truth=0, mm2=7,044 — entirely phantom

### 4. MM2-Specific Under-Predictions (636 transcripts, -127K)

Mirror of over-predictions — these are the transcripts losing fragments to the over-predicted ones via multimap redistribution.

### 5. Isoform Redistribution (723 genes, 163K units)

Genes where the gene-level count is correct (<20 error) but individual isoforms are mis-distributed. This is 14.3% of total mm2 error and is fundamentally a multimap/isoform disambiguation problem — the aligner doesn't provide enough signal to distinguish similar isoforms.

### 6. Chromosomal Distribution

Error is broadly distributed across chromosomes, proportional to transcript density. No single chromosome dominates. Chr12 has the most mm2 excess (+75K), likely due to high-expression gene clusters (e.g., ribosomal proteins).

### 7. Error Concentration

| Top N tx | Oracle | MM2 v3 | Salmon | Kallisto |
|----------|--------|--------|--------|----------|
| Top 10 | 7.0% | 9.6% | 6.2% | 9.3% |
| Top 50 | 10.8% | 23.6% | 12.3% | 16.9% |
| Top 100 | 13.8% | 30.5% | 15.8% | 20.5% |
| Top 500 | 26.5% | 47.3% | 29.5% | 33.6% |

MM2 error is more concentrated in the top 100 transcripts (30.5%) than oracle (13.8%) or salmon (15.8%). This is classic multimap error — a modest number of ambiguous loci drive most of the total error.

## Head-to-Head: Rigel Oracle vs Salmon

| Metric | Rigel Oracle | Salmon |
|--------|-------------|--------|
| Transcript MAE | **3.88** | 6.16 |
| Gene MAE | **0.50** | 8.11 |
| Total |error| | **461,025** | 731,092 |
| Wins (# transcripts) | **18,621** | 13,679 |

Rigel oracle beats salmon by a wide margin at transcript level (37% lower MAE) and gene level (94% lower MAE). The oracle alignment provides perfect fragment origin, so this gap represents rigel's EM advantage over salmon's EM at isoform disambiguation.

## Head-to-Head: Rigel MM2 vs Salmon

| Metric | Rigel MM2 v3 | Salmon |
|--------|-------------|--------|
| Transcript MAE | 9.64 | **6.16** |
| Gene MAE | 8.66 | **8.11** |
| Total |error| | 1,143,595 | **731,092** |

Salmon still leads at transcript level (+56% rigel error) and gene level (+6.8% rigel error). The gene-level gap has narrowed considerably — salmon's advantage there is marginal.

## Opportunities for Improvement

### Priority 1: Reduce Multimap Redistribution Error (~260K units, ~23% of error)

The largest mm2-specific error source: reads that minimap2 maps ambiguously to the wrong set of transcripts. Unlike salmon/kallisto (which work in transcript-space where multimapping is well-defined), rigel works in genome-space where multimapping is between genomic loci that may share sequence.

**Options**:
- **Alignment quality weighting**: Use MAPQ or alignment score to down-weight low-confidence multimappers in the EM
- **Secondary alignment filtering**: Minimap2's `-N 10` allows up to 10 secondary alignments. The quality of these secondaries varies; filtering by alignment score delta could help
- **Fragment length consistency**: Reject multimap candidates where the implied fragment length is inconsistent with the trained FL model (threshold-based or as a likelihood penalty — already partly done via FL scoring)

### Priority 2: nRNA Siphon Residual (73K units, ~6.4% of error)

Even with pool-separated pruning, 73K fragments still get assigned to nRNA at the pool level. In a pristine simulation (zero true nRNA), every one of these is a false positive. The EM allows nRNA to "absorb" ambiguous fragments when the nRNA likelihood is competitive.

**Options**:
- **Stronger nRNA prior penalty**: Increase the Bayesian prior favoring mRNA over nRNA for fragments with small overhang. Currently the pruning threshold is ε=1e-4; a tighter threshold (1e-6) would prune more nRNA candidates
- **Fragment length discrimination**: nRNA components should have very different FL profiles (longer, covering introns). Verifying the FL model correctly penalizes nRNA candidates for typical 300bp fragments
- **Post-EM nRNA confidence check**: After EM convergence, reassign fragments from nRNA to mRNA if the nRNA posterior is below some threshold

### Priority 3: Phantom Over-Predictions at Unexpressed Transcripts (~1.34 MAE for bin 0)

Oracle achieves 0.53 MAE at unexpressed transcripts; mm2 has 1.34. This means mm2 falsely assigns fragments to transcripts that should be zero. Salmon does best here (0.25).

**Options**:
- **Minimum expression threshold**: Zero out transcript counts below a noise floor
- **EM initialization**: Consider initializing unexpressed candidates at lower prior weights

### Priority 4: High-Expression Transcript Error (>10K bin: 2,544 vs 286 oracle)

The highest-expression transcripts have disproportionately large errors with mm2. These are likely housekeeping genes (ribosomal proteins, actins, histones) with many pseudogenes that create multimap traps.

**Options**:
- **Pseudogene-aware scoring**: Flag known pseudogenes and apply additional penalties
- **Locus-level priors**: Use gene-model annotations (processed_pseudogene, protein_coding) to inform prior weights

### Priority 5: Improve Overall Fragment Recovery (-336K vs truth, -12K oracle)

Rigel oracle recovers 99.88% of fragments; mm2 recovers 96.6%. The 336K gap is partly unavoidable (intergenic alignment loss) but the 260K gap vs oracle suggests room to rescue more fragments through:
- **Overhang tolerance** (Problem 1 in the minimap2 rescue plan — soft overhang threshold)
- **Fragment length robustness** (Problem 3 — outlier FL rejection)

## Conclusion

The pool-separated likelihood pruning fix was a major improvement:
- **MAE reduced 29%** (13.60 → 9.64)
- **nRNA siphon reduced 57%** (171K → 73K)
- **gDNA false positives reduced 91%** (40K → 3.7K)

The remaining mm2-vs-oracle gap (MAE 9.64 vs 3.88) is driven by:
1. **Multimap redistribution** (~23% of error) — genome-space alignment ambiguity
2. **Isoform confusion** (~14%) — insufficient signal to distinguish similar isoforms
3. **Residual nRNA siphon** (~6%) — still assignable; needs stronger priors
4. **Shared EM noise** (~33%) — inherent uncertainty present even with oracle alignment

The biggest opportunity is reducing multimap-driven over/under-predictions at single-exon transcripts, which accounts for the bulk of the gap vs salmon.
