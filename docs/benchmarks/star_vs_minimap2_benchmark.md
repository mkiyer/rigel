# STAR vs Minimap2 Aligner Benchmark

**Date**: 2026-04-04
**Dataset**: VCaP prostate cancer simulation (50M RNA fragments, Gencode v47)
**Conditions**: Pristine (0% gDNA) and Dirty (33% gDNA, ss=0.9)

## Motivation

Root cause analysis ([vh_pearson_r_root_cause.md](vh_pearson_r_root_cause.md)) established that the entire VH Pearson R gap between Rigel (0.910) and Kallisto (0.976) was caused by minimap2 alignment failures:
1. **RPS24**: 3bp micro-exon at chr10:78037438-78037441 undetectable by minimap2
2. **RPL21**: 168 processed pseudogenes creating a 299K-transcript mega-locus

STAR was tested as an alternative aligner to see if it resolves these issues.

## Aligner Configuration

### STAR 2.7.11b
- Index: genome_controls.fasta + genes_controls.gtf (sjdbOverhang=149)
- `outFilterMultimapNmax 50` (matching minimap2's `-N 20` with secondary=yes)
- Chimeric detection enabled (chimSegmentMin=10, chimMultimapNmax=50)
- Output: BAM Unsorted → `samtools sort -n` for Rigel

### Minimap2 (existing)
- `minimap2 -ax splice:sr --secondary=yes -N 20` with junction BED

## Alignment Statistics

| Metric | Pristine (STAR) | Pristine (mm2) | Dirty (STAR) | Dirty (mm2) |
|--------|-----------------|----------------|--------------|-------------|
| Input reads | 50M | 50M | 75M | 75M |
| Uniquely mapped | 98.34% | — | 95.93% | — |
| Multi-mapped | 1.63% | ~15%* | 1.99% | — |
| Runtime | 915s | ~500s | ~870s | ~500s |

*Estimated from locus structure; minimap2 doesn't report directly.

## Transcript-Level Pearson R

### Pristine (0% gDNA)

| Tool | Overall | VH (≥1K) | High | Mid | Low |
|------|---------|----------|------|-----|-----|
| oracle_vbem | **0.9985** | **0.9886** | **0.9987** | **0.9931** | **0.9278** |
| kallisto | 0.9964 | 0.9761 | 0.9923 | 0.9861 | 0.8776 |
| **STAR+VBEM** | **0.9919** | **0.9466** | 0.9784 | 0.9783 | 0.7411 |
| **STAR+MAP** | 0.9909 | 0.9401 | 0.9782 | 0.9790 | 0.7408 |
| mm2+VBEM | 0.9864 | 0.9095 | 0.9775 | 0.9822 | 0.8041 |
| mm2+MAP | 0.9865 | 0.9095 | 0.9777 | 0.9822 | 0.8043 |
| salmon | 0.8898 | 0.5646 | 0.8130 | 0.5487 | 0.1522 |

### Dirty (33% gDNA)

| Tool | Overall | VH (≥1K) | High | Mid | Low |
|------|---------|----------|------|-----|-----|
| kallisto | 0.9958 | 0.9765 | 0.9925 | 0.9856 | 0.8484 |
| **STAR+VBEM** | **0.9919** | **0.9466** | 0.9785 | 0.9789 | 0.7408 |
| **STAR+MAP** | 0.9914 | 0.9435 | 0.9783 | 0.9789 | 0.7402 |
| mm2+VBEM | 0.9865 | 0.9094 | 0.9773 | 0.9813 | 0.8025 |
| mm2+MAP | 0.9862 | 0.9096 | 0.9773 | 0.9823 | 0.8033 |
| salmon | 0.8900 | 0.5634 | 0.8127 | 0.5490 | 0.1520 |

## Gene-Level Metrics (Dirty)

| Tool | Gene R | Gene MAPE | Gene WARE | Spearman |
|------|--------|-----------|-----------|----------|
| kallisto | **0.9984** | 781.8% | 0.0805 | 0.7010 |
| **STAR+VBEM** | 0.9962 | **46.5%** | **0.0336** | **0.9007** |
| **STAR+MAP** | 0.9960 | 49.5% | 0.0340 | 0.8990 |
| mm2+VBEM | 0.9954 | 67.6% | 0.0340 | 0.8865 |
| mm2+MAP | 0.9952 | 70.6% | 0.0346 | 0.8825 |
| salmon | 0.9880 | 190.6% | 0.1182 | 0.7785 |

Note: Kallisto has highest gene-level Pearson R but worst MAPE (781.8%) because it cannot distinguish gDNA from mRNA, inflating all gene counts proportionally. Rigel+STAR achieves the best balance of correlation and accuracy.

## Pool-Level gDNA Estimation (Dirty)

| Tool | mRNA Error | False nRNA | gDNA Predicted | gDNA Error |
|------|-----------|------------|----------------|------------|
| mm2+MAP | -0.06% | 2.19M | 21.3M | -14.8% |
| mm2+VBEM | -0.06% | 2.14M | 21.4M | -14.6% |
| **STAR+MAP** | -0.13% | **0.46M** | **23.1M** | **-7.8%** |
| **STAR+VBEM** | -0.14% | **0.28M** | **23.2M** | **-7.0%** |

STAR reduces false nRNA by 7× and halves gDNA estimation error.

## Mega-Locus Comparison

| Metric | Minimap2 | STAR | Oracle |
|--------|----------|------|--------|
| Total loci | 6,268 | 15,574 | 18,000 |
| Max locus (transcripts) | 299,031 | 46,360 | 511 |
| Loci >1000 tx | 2 | 2 | 0 |

STAR's smaller mega-locus (6.5× fewer transcripts) enables much better per-transcript inference.

## Root Cause Resolution

### RPS24 (3bp micro-exon)
| Tool | ENST00000360830.9 (7-exon) | ENST00000372360.9 (6-exon) |
|------|---------------------------|---------------------------|
| Truth | 1469 | 134 |
| Oracle | 1488 | 152 |
| Kallisto | 1451 | 140 |
| **STAR+VBEM** | **1383** | **128** |
| mm2+VBEM | **0** | **1093** |

STAR correctly detects the 3bp micro-exon splice junction, resolving 3,684 reads to the 7-exon isoform.

### RPL21 (168 pseudogenes)
| Tool | ENST00000311549.11 |
|------|-------------------|
| Truth | 1671 |
| Oracle | 1688 |
| Kallisto | 1648 |
| **STAR+VBEM** | **1606** |
| mm2+VBEM | **870** |

STAR's more restrictive multimapping (1.63% vs ~15%) means fewer pseudogene alignments enter the locus, reducing the mega-locus from 299K to 46K transcripts.

## Remaining Gap: ENST00000690968.2 (ENSG00000293320)

After removing the two minimap2 root causes, the remaining STAR VH gap is dominated by a single transcript:

- **ENST00000690968.2**: single-exon lncRNA, 712 bp, chr15:82372204-82372916
- Truth: 1506 TPM, STAR: 415 TPM, Minimap2: 1752 TPM, Kallisto: 1541 TPM
- Removing this one transcript: STAR VH R = 0.975 ≈ Kallisto VH R = 0.976

The cause is under investigation but may involve STAR mapping reads from this single-exon lncRNA to a processed pseudogene elsewhere.

## STAR Low-Expression Trade-off

STAR has worse Pearson R at low expression (0.74 vs 0.80 for minimap2). Root cause:

STAR's better alignment of ribosomal protein genes means more fragments correctly enter their loci, but EM then distributes some to minor non-coding isoforms:

| Transcript | Gene | Truth | STAR | Minimap2 |
|-----------|------|-------|------|----------|
| ENST00000084795.9 | RPL18 (non-coding) | 1.5 | **228** | 0.0 |
| ENST00000401609.5 | RPL3 (non-coding) | 3.9 | **158** | 0.0 |
| ENST00000594493.1 | RPS11 (non-coding) | 2.9 | **76** | 0.8 |

This is the classic isoform ambiguity problem — better alignment → more fragments → harder disambiguation between highly similar isoforms.

## Conclusions

1. **STAR resolves both minimap2 root causes**: micro-exon detection and pseudogene mega-locus
2. **Overall TX Pearson R improves**: 0.986 → 0.992 (STAR+VBEM vs mm2+VBEM)
3. **VH Pearson R improves dramatically**: 0.910 → 0.947
4. **gDNA estimation improves**: error -14.6% → -7.0%, false nRNA 7× reduction
5. **Gene-level MAPE drops**: 67.6% → 46.5% in dirty condition
6. **Trade-off**: slightly worse low-expression accuracy due to isoform spill
7. **After removing 1 problematic lncRNA**: STAR+Rigel matches Kallisto at VH level

STAR is the recommended aligner for Rigel. The low-expression trade-off is minor compared to the dramatic VH and gDNA improvements.
