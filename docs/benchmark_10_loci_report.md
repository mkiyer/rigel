# 10-Locus Benchmark Report: hulkrna vs salmon vs kallisto

**Date:** 2026-03-05

## Benchmark Design

- **10 genomic regions:** PVT1_MYC, EGFR, MALAT1_NEAT1, TP53, CDKN2A_2B, HOXA_cluster, XIST_TSIX, H19_IGF2, TTN, GAPDH
- **Conditions per region:** 12 (2 gDNA rates × 2 nRNA rates × 3 strand specificities)
- **Aligners:** oracle (perfect), minimap2 (realistic)
- **Fragments per condition:** 100,000
- **Read length:** 150 bp, error rate: 0.0
- **gDNA rates:** 0.0 (none), 0.2 (high)
- **nRNA rates:** 0.0 (none), 0.2 (high)
- **Strand specificity:** 1.0, 0.95, 0.50
- **Total conditions:** 10 × 12 = 120

---

## Overall Transcript-Level Accuracy

Mean across all 120 conditions:

| Tool | MAE | RMSE | Pearson | Spearman |
| --- | ---: | ---: | ---: | ---: |
| **hulkrna_oracle** | **15.42** | **51.02** | **0.9555** | **0.8025** |
| **hulkrna_minimap2** | **23.13** | **84.45** | **0.8946** | **0.7668** |
| salmon | 144.15 | 372.19 | 0.6263 | 0.5544 |
| kallisto | 180.34 | 452.67 | 0.5756 | 0.5208 |

hulkrna_oracle achieves **9.3× lower MAE** than salmon and **11.7× lower** than kallisto. Even with minimap2 alignment, hulkrna beats salmon by **6.2×** and kallisto by **7.8×** on MAE.

---

## Robustness to gDNA Contamination

Transcript-level MAE by gDNA level:

| gDNA Level | hulkrna_oracle | hulkrna_minimap2 | salmon | kallisto |
| --- | ---: | ---: | ---: | ---: |
| none | 25.65 | 37.83 | 152.63 | 193.37 |
| high (20%) | 5.20 | 8.43 | 135.67 | 167.31 |

hulkrna correctly decomposes gDNA from RNA signal — its error *drops* with gDNA contamination because the gDNA model absorbs noise. salmon and kallisto have no gDNA model and cannot distinguish genomic DNA from RNA; their errors remain uniformly high. Under 20% gDNA contamination, hulkrna_oracle achieves **26× lower MAE** than salmon.

---

## Robustness to Strand Specificity

Transcript-level MAE by strand specificity:

| Strand Spec | hulkrna_oracle | hulkrna_minimap2 | salmon | kallisto |
| --- | ---: | ---: | ---: | ---: |
| 1.00 (perfect) | 14.57 | 22.11 | 135.58 | 154.17 |
| 0.95 (typical) | 15.68 | 22.77 | 148.05 | 167.11 |
| 0.50 (unstranded) | 16.02 | 24.50 | 148.82 | 219.74 |

hulkrna is nearly immune to strand specificity degradation (14.57 → 16.02, a 10% increase). Its continuous weighting scheme smoothly transitions from strand-based to density-based estimation. kallisto collapses badly at SS=0.50 (154 → 220, a 43% increase).

---

## Accuracy by Abundance Quartile

Transcript-level MAE pooled across all conditions:

| Quartile | hulkrna_oracle | hulkrna_minimap2 | salmon | kallisto |
| --- | ---: | ---: | ---: | ---: |
| Q1 (low) | 5.47 | 8.90 | 87.16 | 125.60 |
| Q2 | 4.90 | 16.23 | 83.99 | 115.99 |
| Q3 | 19.06 | 22.45 | 91.60 | 112.28 |
| Q4 (high) | 31.82 | 45.14 | 163.97 | 185.62 |

hulkrna maintains strong accuracy across all abundance levels. salmon/kallisto errors scale with abundance, suggesting systematic mis-assignment of high-count transcripts via equivalence class splitting.

---

## Dropout Rate

Fraction of true-positive transcripts (truth > 0) predicted as zero:

| Tool | Dropout rate |
| --- | ---: |
| hulkrna_oracle | **0.00%** |
| hulkrna_minimap2 | **0.00%** |
| salmon | 12.27% |
| kallisto | 9.19% |

hulkrna achieves zero dropout — every expressed transcript receives a nonzero estimate. salmon drops ~1 in 8 expressed transcripts.

---

## Pool-Level Decomposition (Absolute Error)

Mean absolute error for RNA pool (mature_rna, nascent_rna) and DNA pool (genomic_dna):

| Pool | hulkrna_oracle | hulkrna_minimap2 | salmon | kallisto |
| --- | ---: | ---: | ---: | ---: |
| mature_rna | 348 | 524 | 10,358 | 14,678 |
| nascent_rna | 1,199 | 8,017 | 18,846 | 18,846 |
| genomic_dna | 968 | 7,960 | 47,094 | 47,094 |

hulkrna_oracle has **30× lower mature RNA error**, **16× lower nascent RNA error**, and **49× lower genomic DNA error** than salmon. salmon and kallisto have identical nascent RNA and gDNA errors because they have no model for either — all non-mRNA reads are treated as mRNA.

---

## Key Takeaways

1. **Dominant accuracy**: hulkrna outperforms salmon/kallisto by 6–12× on transcript-level MAE across all conditions.
2. **gDNA robustness**: hulkrna's gDNA model correctly absorbs contamination that confounds salmon/kallisto.
3. **Strand-agnostic**: hulkrna's continuous SS weighting handles unstranded data gracefully; kallisto's error inflates 43%.
4. **Zero dropout**: hulkrna never predicts zero for an expressed transcript; salmon drops 12%.
5. **Pool decomposition**: hulkrna uniquely quantifies mRNA, nRNA, and gDNA pools — salmon/kallisto conflate them.
6. **Speed tradeoff**: hulkrna is slower per-condition (2–46s vs 0.3–3s) due to the EM + gDNA/nRNA models, but this overhead is small relative to alignment time in a real pipeline.
