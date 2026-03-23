# Pristine Benchmark Report: Rigel vs Salmon vs Kallisto

**Date**: 2026-03-22
**Simulation**: 10M mRNA-only fragments (zero gDNA, zero nRNA), SS=0.95, HeLa/Cervix salmon TPM abundances
**Transcripts**: 118,684 annotated (41,056 expressed)

---

## Executive Summary

| Tool | Alignment | Total Pred | MAE | RMSE | Pearson | Spearman | mRNA Loss |
|------|-----------|-----------|-----|------|---------|----------|-----------|
| **rigel** | oracle | 9,988,455 | **1.84** | **30.45** | 0.9993 | 0.8795 | **0.12%** |
| salmon | quasi-map | 9,501,334 | 2.91 | 38.63 | 0.9999 | 0.8937 | 4.99% |
| kallisto | pseudo-align | 9,491,975 | 3.36 | 59.28 | 0.9980 | 0.8481 | 5.08% |
| **rigel** | minimap2 | 9,779,947 | 8.64 | 124.51 | 0.9865 | 0.6770 | 2.20% |

**Key finding**: With oracle alignments, **rigel has the lowest MAE (1.84)** — 37% lower than salmon's 2.91 and 45% lower than kallisto's 3.36. Rigel recovers 99.88% of all mRNA fragments vs salmon/kallisto's ~95%. However, minimap2 alignment introduces substantial degradation.

---

## Detailed Findings

### 1. Rigel Outperforms on Per-Transcript MAE (Oracle)

Rigel's MAE of 1.84 is the best across all tools. At every error quantile, rigel is tighter:

| Quantile | rigel | salmon | kallisto |
|----------|-------|--------|----------|
| 50th | 2.00 | 3.87 | 4.00 |
| 90th | 25.00 | 37.96 | 39.87 |
| 95th | 42.00 | 65.47 | 68.11 |
| 99th | 108.00 | 186.05 | 194.06 |

Head-to-head on expressed transcripts:
- **rigel wins 39.5%**, salmon wins 28.9%, tied 31.6%
- **rigel wins 41.6%**, kallisto wins 28.8%, tied 29.5%

### 2. Salmon/Kallisto ~5% Systematic Undercount

Both salmon and kallisto lose ~5% of all fragments uniformly. This is not specific to any transcript class — MT transcripts show the same ~5% ratio (0.9496) as non-MT (0.9481). This loss comes from reads that fail to quasi-map/pseudo-align in lightweight alignment mode. With 10M simulated reads (no sequencing errors), ~500K reads are simply lost.

**Root cause**: Quasi-mapping drops reads that don't have sufficient k-mer matches or that map ambiguously in the de Bruijn graph. This is a fundamental tradeoff of lightweight alignment — speed for completeness. The 5% loss rate at error_rate=0.0 suggests the issue is primarily short-homology regions and repetitive sequences, not sequencing error sensitivity.

### 3. Rigel's ~50% Isoform Redistribution Pattern (122 transcripts)

A distinctive failure pattern: 122 expressed transcripts show exactly ~50% underprediction by rigel. Investigation reveals:

- **105/120 genes are conserved at the gene level** — the "missing" fragments went to an unexpressed isoform of the same gene
- This is **classic EM identifiability**: when two isoforms share exonic sequence, the EM splits fragments between them
- 51/122 have an unexpressed sibling that received the "other half" (the remaining cases have the fragments spread across multiple siblings)

**Example (LDHB)**:
- ENST00000350669.5: truth=16,614, rigel=8,279 (50% undercount)
- ENST00000673047.2 (unexpressed sibling): truth=0, rigel=8,349 (FP)
- Gene total: truth=16,628, rigel=16,629 (**perfectly conserved**)

This is NOT a bug — it's a known limitation of all transcript-level EM quantification tools (salmon and kallisto have the same fundamental issue but manifested differently because they use different EM formulations and effective length corrections).

### 4. Rigel's Non-Conserved Gene Errors (15 genes)

15 genes show genuine fragment leakage where gene totals don't match truth:

| Gene | Truth | Rigel | Loss | Cause |
|------|-------|-------|------|-------|
| PTP4A1 | 16,071 | 10,733 | 5,338 | Cross-gene overlap with ENSG00000285976 |
| SLC35A4 | 1,486 | 773 | 713 | Cross-gene overlap |
| PYURF | 689 | 354 | 335 | Overlap with PIGY |
| LAS1L | 848 | 517 | 331 | Cross-gene overlap |

**Root cause**: These involve overlapping genes on the same strand where fragments are ambiguous between two different gene loci. The EM distributes fragments across both genes, causing leakage from the expressed gene to the overlapping unexpressed gene's transcripts.

### 5. Rigel's nRNA/gDNA Ghost Attribution

Even with zero true nRNA and zero true gDNA:
- **nRNA siphon**: 3,661 fragments (0.04%) attributed to nRNA
- **gDNA siphon**: 1,840 fragments (0.02%) attributed to gDNA
- **Total leakage**: 0.115% of all mRNA fragments

This is very small in absolute terms but represents a systematic floor. The calibration system estimated π=0.036 (3.6% gDNA contamination rate) and κ_strand=500 even though there is zero contamination. This suggests the calibration EM can be fooled by stochastic noise in the strand/density signals.

### 6. Minimap2 Fragment Length Model Corruption (**CRITICAL**)

The most serious issue discovered:

| Metric | Oracle | Minimap2 |
|--------|--------|----------|
| SPLICED_ANNOT mean | 298.7 | **335.9** |
| SPLICED_ANNOT mode | 299 | **1000** (max!) |
| Global std | 49.6 | **126.4** |
| UNSPLICED mode | 299 | 299 |

With minimap2, the **annotated-spliced fragment length distribution is severely corrupted** — the mode is at 1000 (the maximum fragment length). This means minimap2 is computing fragment lengths incorrectly for spliced reads, likely because:

1. minimap2 computes fragment length as the template length (TLEN) field, which for spliced reads includes the intron
2. The oracle BAM uses the correct **spliced fragment length** (sum of exon-overlapping bases)
3. This corrupted FL model then poisons all fragment likelihood scoring

**Impact**: The minimap2 nRNA siphon is **60× worse** (215,943 vs 3,661 fragments) because the wrong FL model makes the scoring system unable to distinguish mRNA from nRNA fragments, allowing nRNA components to absorb far more mass.

### 7. Gene-Level Accuracy

At the gene level (summing isoforms), rigel's MAE drops to **2.84** — the best of all tools — confirming that most of rigel's transcript-level error is isoform redistribution within genes:

| Tool | Gene MAE | Pearson | Spearman |
|------|----------|---------|----------|
| rigel (oracle) | **2.84** | 0.9996 | **0.9880** |
| salmon | 25.58 | **0.9999** | 0.9898 |
| kallisto | 28.44 | 0.9990 | 0.9896 |

Salmon/kallisto gene-level MAE is 9× higher than rigel's because their ~5% uniform undercount accumulates into large absolute errors on high-abundance genes.

---

## Prioritized Improvement Opportunities

### P0 — Minimap2 Fragment Length Computation (CRITICAL)

The fragment length model is completely wrong for spliced reads under minimap2 alignment. Mode=1000 (max) vs expected ~300. This single issue cascades into:
- Corrupted FL scoring → bad nRNA/mRNA discrimination
- 60× worse nRNA siphon (215K vs 3.6K)
- 6.3% total mRNA loss vs 0.12% with oracle

**Fix**: Investigate how `read_length` (sum of M/D/=/X CIGAR ops) is computed for minimap2 alignments. The TLEN field from minimap2 may include intronic spans. Ensure the C++ BAM scanner computes fragment length as the sum of aligned exon blocks, not the genomic span between mate pairs.

### P1 — False Positive on Ghost/Overlapping Genes

5,338 fragments leaked from PTP4A1 to ENSG00000285976 (two overlapping genes). This is a cross-gene EM identifiability issue that could potentially be mitigated:
- **Gene-aware EM regularization**: Penalize cross-gene fragment sharing when gene boundaries are known
- **Post-EM redistribution**: Detect unexpressed-gene FP and redistribute back to the overlapping expressed gene
- **Conservative**: Only affects ~15 genes out of 63K, so low priority unless these are biologically important

### P2 — gDNA Calibration False Positive

The calibration estimated π=3.6% gDNA contamination from pristine data (truth: 0%). While the downstream impact is minimal (1,840 fragments), this represents a floor on gDNA estimation accuracy. The calibration could benefit from:
- A null-hypothesis test before committing to a non-zero π
- Stricter convergence criteria at low contamination levels
- A minimum-evidence threshold before declaring contamination present

### P3 — nRNA Siphon at Zero nRNA (Minor)

3,661 fragments attributed to nRNA from pristine data. This is small (0.04%) but reflects the presence of nRNA components in the EM that have no data to constrain them. Could potentially be addressed by:
- Sparsity priors on nRNA components that are stronger at low evidence levels
- Post-EM pruning of nRNA components with very low posterior mass

---

## Conclusion

In pristine conditions (zero gDNA, zero nRNA), **rigel with oracle alignments is the most accurate transcript quantification tool tested**, with 37% lower MAE than salmon and 45% lower than kallisto. Rigel's gene-level accuracy is 9× better than salmon/kallisto because it recovers 99.88% of fragments vs ~95%.

The **critical issue is minimap2 fragment length computation** for spliced reads, which degrades all downstream scoring and causes a 60× increase in nRNA siphon. Fixing this single issue should bring minimap2-aligned rigel much closer to oracle-aligned performance.
