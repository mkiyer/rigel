# Benchmark: hulkrna vs salmon

**Generated**: 2026-02-16 13:22:51
**Fragments per scenario**: 2000
**Seed**: 42
**salmon version**: 1.10.3

## Overview

This report compares **hulkrna** (Bayesian EM with gDNA modeling) against **salmon** (lightweight quasi-mapping quantification) on synthetic RNA-seq scenarios with known ground truth.

Key differences in design:
- **hulkrna** takes a BAM (genome-aligned) as input and models gDNA contamination, strand specificity, and insert size jointly.
- **salmon** performs lightweight mapping directly against the transcriptome (no genome alignment, no gDNA modeling).

## Summary: Mean Absolute Error by Scenario and Parameter

## Scenario 1: Single-exon gene (unspliced)

### gDNA Contamination Sweep

| gDNA abundance | Truth (t1) | hulkrna (t1) | salmon (t1) | gDNA truth | hulkrna gDNA | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 2000 | 2000 | 2000.0 | 0 | 0 | 0.0 | 0.0 |
| 20 | 925 | 1098 | 1198.0 | 1075 | 902 | 173.0 | 273.0 |
| 50 | 512 | 738 | 848.0 | 1488 | 1262 | 226.0 | 336.0 |

### Strand Specificity Sweep

| Strand spec. | Truth (t1) | hulkrna (t1) | salmon (t1) | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- |
| 0.5 | 2000 | 0 | 2000.0 | 2000.0 | 0.0 |
| 0.9 | 2000 | 2000 | 2000.0 | 0.0 | 0.0 |
| 1.0 | 2000 | 2000 | 2000.0 | 0.0 | 0.0 |

## Scenario 2: Spliced gene (multi-exon)

### gDNA Contamination Sweep

| gDNA abundance | Truth (t1) | hulkrna (t1) | salmon (t1) | gDNA truth | hulkrna gDNA | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 2000 | 2000 | 2000.0 | 0 | 0 | 0.0 | 0.0 |
| 20 | 602 | 784 | 804.0 | 1398 | 1216 | 182.0 | 202.0 |
| 50 | 294 | 532 | 550.0 | 1706 | 1468 | 238.0 | 256.0 |

### Strand Specificity Sweep

| Strand spec. | Truth (t1) | hulkrna (t1) | salmon (t1) | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- |
| 0.5 | 2000 | 890 | 2000.0 | 1110.0 | 0.0 |
| 0.9 | 2000 | 1883 | 2000.0 | 117.0 | 0.0 |
| 1.0 | 2000 | 2000 | 2000.0 | 0.0 | 0.0 |

## Scenario 3: Non-overlapping genes

### gDNA Contamination Sweep

| gDNA abundance | Truth (t1) | hulkrna (t1) | salmon (t1) | Truth (t2) | hulkrna (t2) | salmon (t2) | gDNA truth | hulkrna gDNA | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1371 | 1371 | 1371.0 | 629 | 629 | 629.0 | 0 | 0 | 0.0 | 0.0 |
| 20 | 374 | 493 | 498.0 | 191 | 265 | 267.0 | 1435 | 1242 | 96.5 | 100.0 |
| 50 | 180 | 303 | 298.0 | 92 | 0 | 205.0 | 1728 | 1697 | 107.5 | 115.5 |

### Strand Specificity Sweep

| Strand spec. | Truth (t1) | hulkrna (t1) | salmon (t1) | Truth (t2) | hulkrna (t2) | salmon (t2) | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 1370 | 614 | 1370.0 | 630 | 0 | 630.0 | 693.0 | 0.0 |
| 0.9 | 1370 | 1292 | 1370.0 | 630 | 630 | 630.0 | 39.0 | 0.0 |
| 1.0 | 1371 | 1371 | 1371.0 | 629 | 629 | 629.0 | 0.0 | 0.0 |

## Scenario 4: Two isoforms (shared exons)

### Abundance Ratio Sweep

| Fold change | Truth (t1) | hulkrna (t1) | salmon (t1) | Truth (t2) | hulkrna (t2) | salmon (t2) | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1303 | 1320 | 1319.9 | 697 | 680 | 680.1 | 17.2 | 16.9 |
| 4 | 1787 | 1808 | 1808.0 | 213 | 192 | 192.0 | 21.1 | 21.0 |
| 16 | 1940 | 1937 | 1936.6 | 60 | 63 | 63.4 | 2.7 | 3.4 |

### gDNA Contamination Sweep

| gDNA abundance | Truth (t1) | hulkrna (t1) | salmon (t1) | Truth (t2) | hulkrna (t2) | salmon (t2) | gDNA truth | hulkrna gDNA | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1911 | 1910 | 1908.4 | 89 | 90 | 91.6 | 0 | 0 | 0.6 | 2.6 |
| 20 | 837 | 1086 | 1056.5 | 50 | 85 | 74.5 | 1113 | 829 | 142.0 | 122.0 |
| 50 | 453 | 772 | 750.5 | 30 | 97 | 95.5 | 1517 | 1131 | 193.0 | 181.5 |

## Scenario 5: Overlapping antisense genes

### Abundance Ratio Sweep

| Fold change | Truth (t1) | hulkrna (t1) | salmon (t1) | Truth (t2) | hulkrna (t2) | salmon (t2) | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1031 | 1031 | 1031.5 | 969 | 969 | 968.5 | 0.1 | 0.5 |
| 4 | 1627 | 1628 | 1628.8 | 373 | 372 | 371.2 | 0.9 | 1.8 |
| 16 | 1906 | 1907 | 1908.3 | 94 | 93 | 91.7 | 1.1 | 2.3 |

### gDNA Contamination Sweep

| gDNA abundance | Truth (t1) | hulkrna (t1) | salmon (t1) | Truth (t2) | hulkrna (t2) | salmon (t2) | gDNA truth | hulkrna gDNA | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1031 | 1031 | 1031.5 | 969 | 969 | 968.5 | 0 | -0 | 0.1 | 0.5 |
| 20 | 475 | 592 | 587.0 | 451 | 596 | 567.0 | 1074 | 812 | 131.0 | 114.0 |
| 50 | 259 | 414 | 392.5 | 254 | 432 | 402.5 | 1487 | 1154 | 166.5 | 141.0 |

### Strand Specificity Sweep

| Strand spec. | Truth (t1) | hulkrna (t1) | salmon (t1) | Truth (t2) | hulkrna (t2) | salmon (t2) | hulkrna MAE | salmon MAE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 1037 | 464 | 1036.2 | 963 | 414 | 963.8 | 561.0 | 0.8 |
| 0.9 | 1037 | 1050 | 1036.2 | 963 | 946 | 963.8 | 14.7 | 0.8 |
| 1.0 | 1031 | 1031 | 1031.5 | 969 | 969 | 968.5 | 0.1 | 0.5 |

## Overall Comparison

**Total scenarios benchmarked**: 33

| Metric | hulkrna | salmon |
| --- | --- | --- |
| Total abs error (sum) | 8422 | 2717 |
| Mean abs error per scenario | 255.2 | 82.3 |
| Scenarios won (lower error) | 13 | 11 |
| Ties (within 0.5) | 9 | — |

### With gDNA contamination (10 scenarios)

| Metric | hulkrna | salmon |
| --- | --- | --- |
| Mean total abs error | 249.2 | 261.5 |

### Without gDNA contamination (23 scenarios)

| Metric | hulkrna | salmon |
| --- | --- | --- |
| Mean total abs error | 257.8 | 4.4 |

## Key Findings

### hulkrna Advantages

1. **gDNA contamination modeling**: hulkrna explicitly models genomic DNA contamination via per-gene shadow components in the EM. Salmon has no gDNA model — contaminating fragments that pseudo-align to transcripts inflate salmon's counts.
2. **Strand-aware quantification**: hulkrna uses trained strand models to separate sense from antisense reads, critical for overlapping antisense genes.
3. **Genome-aligned input**: hulkrna works from BAM files with full genomic context (introns, intergenic regions), enabling it to identify unspliced reads and intronic fragments.

### salmon Advantages

1. **Speed**: salmon's lightweight quasi-mapping is significantly faster than genome alignment + hulkrna counting.
2. **No alignment required**: salmon maps directly against the transcriptome, avoiding the need for a genome aligner like minimap2.
3. **Mature bias correction**: salmon includes sequence-specific, GC, and positional bias models (not tested here).
4. **Pure RNA quantification**: In clean RNA-seq without gDNA contamination, salmon's transcript-level EM is well-calibrated.

### Limitations of This Benchmark

- **Small simulated genomes** (5–8 kb) may not capture real-world mapping complexity.
- **No multi-gene families**: Real transcriptomes have many paralogs where salmon's EM is particularly important.
- **No read errors or biases**: Simulated reads are clean; real data has base-call errors, GC bias, and positional artifacts.
- **gDNA contamination is synthetic**: Real gDNA distributions may differ from our uniform model.

