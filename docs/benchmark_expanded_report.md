# Expanded Benchmark Report: 20-Loci, 120 Conditions

## Setup

- **10 regions**: CDKN2A_2B, EGFR, GAPDH, H19_IGF2, HOXA_cluster, MALAT1_NEAT1, PVT1_MYC, TP53, TTN, XIST_TSIX
- **Conditions**: gDNA ∈ {0.0, 0.2}, nRNA ∈ {0.0, 0.1}, SS ∈ {1.0, 0.95, 0.5}
- **Aligners**: oracle (perfect placement), minimap2
- **Competitors**: salmon, kallisto
- **100k reads per region per condition**
- **Total**: 10 × 4 × 3 × 2 = 240 tool evaluations (120 unique conditions)

---

## 1. Headline Results

### Mean Transcript Spearman by Condition Class

| Condition | rigel_oracle | rigel_mm2 | salmon | kallisto |
|-----------|---------------|-------------|--------|----------|
| **Pristine** (gDNA=0, nRNA=0) | **0.9444** | 0.9270 | 0.9211 | 0.9371 |
| **nRNA only** (gDNA=0, nRNA=0.1) | **0.8123** | 0.7738 | 0.7120 | 0.6332 |
| **gDNA only** (gDNA=0.2, nRNA=0) | **0.5804** | 0.5740 | 0.3231 | 0.2741 |
| **gDNA + nRNA** (gDNA=0.2, nRNA=0.1) | **0.5830** | 0.5459 | 0.3142 | 0.2577 |

**rigel leads in every condition class.** The advantage grows with contamination:
pristine +0.007 → nRNA +0.10 → gDNA +0.26 → both +0.27 (vs best competitor).

### Mean Gene Spearman

| Condition | rigel_oracle | rigel_mm2 | salmon | kallisto |
|-----------|---------------|-------------|--------|----------|
| **Pristine** | 0.9840 | 0.9824 | **0.9900** | **0.9900** |
| **nRNA only** | **0.9594** | 0.9447 | 0.9320 | 0.9090 |
| **gDNA only** | 0.8072 | **0.8280** | 0.7884 | 0.7638 |
| **gDNA + nRNA** | **0.7869** | 0.7841 | 0.7615 | 0.7343 |

Gene-level: rigel **trails by 0.006** in pristine (the gene aggregation smooths out errors); leads everywhere else.

### Win/Loss Record (rigel_oracle vs best{salmon,kallisto})

| Condition | Wins | Losses | Ties | Mean Δ |
|-----------|------|--------|------|--------|
| Pristine | 11 | 0 | 19 | +0.007 |
| nRNA only | 25 | 4 | 1 | +0.100 |
| gDNA only | 27 | 3 | 0 | +0.255 |
| gDNA + nRNA | 27 | 3 | 0 | +0.265 |

**All 7 losses are from only 2 regions**: EGFR (6, gDNA conditions) and TP53 (3, nRNA-only) + XIST_TSIX (1, nRNA-only).

---

## 2. Pristine Case (gDNA=0, nRNA=0)

### nRNA Leak: ZERO

All 10 regions × 3 SS values: **nRNA observed = 0.0** (truth = 0). The empirical
evidence-weighted nrna_frac prior fix completely eliminates phantom nRNA attribution.

### Per-Region Transcript Spearman (SS=1.0)

| Region | rigel_oracle | salmon | kallisto | Result |
|--------|---------------|--------|----------|--------|
| CDKN2A_2B | 0.9922 | 0.9927 | 0.9926 | TIE |
| EGFR | **0.9995** | 0.9897 | 0.9995 | WIN |
| GAPDH | **0.8996** | 0.8695 | 0.8980 | WIN |
| H19_IGF2 | **0.9665** | 0.9209 | 0.9616 | WIN |
| HOXA_cluster | **0.9536** | 0.9472 | 0.9531 | WIN |
| MALAT1_NEAT1 | 0.9378 | 0.9352 | **0.9381** | TIE |
| PVT1_MYC | **0.8964** | 0.8373 | 0.8794 | WIN |
| TP53 | 0.9419 | 0.8739 | **0.9451** | TIE |
| TTN | **0.9635** | 0.9562 | 0.9618 | WIN |
| XIST_TSIX | **0.9217** | 0.9115 | 0.9062 | WIN |

**Pristine verdict**: 8 wins, 0 losses, 2 marginal ties (within 0.003).
rigel is clearly competitive in clean conditions.

---

## 3. nRNA Contamination (gDNA=0, nRNA=0.1)

### nRNA Pool Accuracy (rigel_oracle, SS=1.0)

| Region | nRNA truth | nRNA observed | Error % |
|--------|-----------|--------------|---------|
| PVT1_MYC | 93,700 | 93,802 | +0.1% |
| CDKN2A_2B | 82,215 | 82,007 | -0.3% |
| TTN | 44,816 | 45,057 | +0.5% |
| XIST_TSIX | 74,392 | 74,834 | +0.6% |
| H19_IGF2 | 31,208 | 31,720 | +1.6% |
| EGFR | 75,963 | 77,564 | +2.1% |
| HOXA_cluster | 41,908 | 41,142 | -1.8% |
| TP53 | 38,118 | 39,106 | +2.6% |
| GAPDH | 42,368 | 43,458 | +2.6% |
| MALAT1_NEAT1 | 29,508 | 26,404 | **-10.5%** |

**Pool accuracy is excellent** (8/10 within ±3%). MALAT1 under-counts nRNA by 10%.

### Per-Read nRNA Classification (oracle, SS=1.0)

| Region | nRNA reads | Correct pool | Wrong pool (→ mRNA) |
|--------|-----------|-------------|---------------------|
| EGFR | 75,963 | 98.4% | 1.6% |
| PVT1_MYC | 93,700 | 98.4% | 1.6% |
| CDKN2A_2B | 82,215 | 97.4% | 2.6% |
| GAPDH | 42,368 | 95.3% | 4.7% |
| TTN | 44,816 | 92.0% | 8.0% |
| XIST_TSIX | 74,392 | 91.1% | 8.9% |
| TP53 | 38,118 | 89.3% | 10.7% |
| HOXA_cluster | 41,908 | 87.8% | 12.2% |
| H19_IGF2 | 31,208 | 80.6% | 19.4% |
| MALAT1_NEAT1 | 29,508 | 72.9% | **27.1%** |

Some regions have significant nRNA→mRNA leakage (MALAT1 27%, H19 19%). These are
regions where nascent RNA reads overlap heavily with exonic regions.

### Losses in nRNA-only

| Region | SS | rigel | competitor | Δ |
|--------|-----|---------|-----------|------|
| TP53 | 1.0 | 0.7114 | 0.7679 | -0.057 |
| TP53 | 0.50 | 0.7251 | 0.7727 | -0.048 |
| TP53 | 0.95 | 0.7295 | 0.7704 | -0.041 |
| XIST_TSIX | 1.0 | 0.6841 | 0.6967 | -0.013 |

TP53 is a consistent (but modest) loss: ~5% gap. XIST_TSIX is marginal.

---

## 4. gDNA Contamination (gDNA=0.2, nRNA=0)

### gDNA Pool Accuracy

Consistently excellent: all regions within ±3.2% error. Hulkrna's gDNA model
captures the correct pool size.

### Residual nRNA Leak (gDNA_only → nRNA should be 0)

| Region | nRNA observed (should be 0) |
|--------|---------------------------|
| EGFR | **4,602** |
| GAPDH | 1,753 |
| TP53 | 807 |
| CDKN2A_2B | 675 |
| PVT1_MYC | 493 |
| H19_IGF2 | 510 |
| TTN | 461 |
| MALAT1_NEAT1 | 377 |
| HOXA_cluster | 145 |
| XIST_TSIX | 131 |

EGFR has **10x the leak** of most regions. This phantom nRNA is gDNA reads that the
empirical nrna_frac system incorrectly attributes to nRNA.

### mRNA Absorption by gDNA (percentage of mRNA reads misclassified as gDNA)

| Region | mRNA reads | → gDNA (wrong pool) | Correctly assigned |
|--------|-----------|--------------------|--------------------|
| XIST_TSIX | 1,638 | **76.3%** | 23.7% |
| MALAT1_NEAT1 | 1,698 | **68.7%** | 31.3% |
| HOXA_cluster | 443 | **62.3%** | 37.7% |
| H19_IGF2 | 2,077 | **58.2%** | 40.9% |
| CDKN2A_2B | 1,393 | **51.5%** | 48.4% |
| TTN | 6,998 | **46.1%** | 53.9% |
| EGFR | 4,773 | **42.7%** | 57.3% |
| TP53 | 1,883 | **40.4%** | 56.2% |
| PVT1_MYC | 696 | **31.0%** | 69.0% |
| GAPDH | 3,572 | **22.7%** | 76.6% |

**This is the #1 systematic weakness**: 23–76% of mRNA reads are incorrectly absorbed by
the gDNA model. The gDNA model's flat coverage profile competes with and "steals" reads
from transcripts, especially reads in single-exon or low-complexity transcript regions.

### Spearman Degradation (SS=1.0) — rigel still beats salmon despite mRNA absorption

| Region | rigel pristine→gDNA drop | sal drop | rigel still wins? |
|--------|------------------------|----------|-----------------|
| EGFR | -0.369 | -0.250 | **NO** (0.63 vs 0.74) |
| PVT1_MYC | -0.449 | -0.747 | YES (0.45 vs 0.09) |
| XIST_TSIX | -0.551 | -0.553 | YES (0.37 vs 0.36) |
| CDKN2A_2B | -0.202 | -0.636 | YES (0.79 vs 0.36) |
| GAPDH | -0.277 | -0.584 | YES (0.62 vs 0.29) |

Despite massive mRNA absorption, rigel still BEATS salmon/kallisto in 7/10 regions
because salmon/kallisto can't model gDNA at all and just hallucinate counts everywhere.
**The only loss is EGFR.**

---

## 5. Combined gDNA + nRNA (gDNA=0.2, nRNA=0.1)

### Critical Finding: gDNA Absorbs nRNA

| Region | nRNA truth | nRNA observed | Error | gDNA over-est |
|--------|-----------|--------------|-------|---------------|
| PVT1_MYC | 9,385 | 678 | **-92.8%** | +9.8% |
| TTN | 5,378 | 489 | **-90.9%** | +7.3% |
| XIST_TSIX | 4,543 | 179 | **-96.1%** | +5.6% |
| CDKN2A_2B | 6,052 | 967 | **-84.0%** | +5.5% |
| EGFR | 13,105 | 6,739 | **-48.6%** | +9.5% |
| GAPDH | 2,559 | 1,800 | **-29.7%** | +0.4% |
| MALAT1_NEAT1 | 706 | 367 | **-48.0%** | +0.7% |
| TP53 | 1,147 | 860 | **-25.1%** | +0.4% |

The nRNA+gDNA signed errors balance almost perfectly
(e.g., PVT1: nRNA −8707 + gDNA +8807 = +100), confirming **direct transfer from nRNA to gDNA**.

The EM cannot distinguish nRNA from gDNA because both have uniform (non-transcript-structured)
coverage profiles. The gDNA pool is 10–100x larger, so it "wins" the competition for
ambiguous reads.

---

## 6. EGFR Deep Dive — The One Persistent Loss

### Why EGFR is Unique

EGFR is the **only region where rigel consistently loses** to salmon/kallisto in gDNA conditions.

| Condition | rigel | salmon | kallisto |
|-----------|---------|--------|----------|
| Pristine SS=1.0 | **0.9995** | 0.9897 | 0.9995 |
| gDNA SS=1.0 | 0.6305 | **0.7400** | 0.7422 |
| gDNA+nRNA SS=1.0 | 0.5616 | **0.7902** | 0.7276 |

### The Dominant Transcript Problem

EGFR transcript ENST00000275493.7 (the canonical EGFR mRNA, abundance=6879 TPM):

| Condition | Truth | rigel | salmon | kallisto |
|-----------|-------|---------|--------|----------|
| Pristine | 42,173 | 41,615 | 39,563 | 41,495 |
| gDNA only | 1,991 | **0** | 2,183 | 1,693 |
| gDNA+nRNA | 1,722 | **0** | 1,923 | 1,277 |
| nRNA only | 10,085 | **9,420** | 8,195 | 5,048 |

The #1 EGFR transcript goes to **exactly zero** in rigel when gDNA is present!
The gDNA model completely absorbs it. Salmon doesn't have this failure mode because
it has no gDNA concept.

### Gene-Level Collapse

EGFR gene-level Spearman with gDNA: rigel 0.55 vs salmon 0.96.
The gDNA model not only misattributes transcript-level reads but also destroys gene-level structure.

---

## 7. Strand Specificity Impact

| SS | rigel_oracle | salmon | kallisto |
|----|---------------|--------|----------|
| 1.00 | **0.7318** | 0.5795 | 0.5458 |
| 0.95 | **0.7292** | 0.5666 | 0.5293 |
| 0.50 | **0.7290** | 0.5567 | 0.5014 |

SS has **minimal effect** on rigel (0.003 drop from 1.0→0.5), while
kallisto drops 0.04. This suggests rigel's strand model is not heavily leveraged yet —
potential for improvement.

## 8. Minimap2 Alignment Penalty

| Condition | Oracle | Minimap2 | Delta |
|-----------|--------|----------|-------|
| Pristine | 0.9444 | 0.9270 | -0.017 |
| nRNA only | 0.8123 | 0.7738 | -0.039 |
| gDNA only | 0.5804 | 0.5740 | -0.006 |
| gDNA + nRNA | 0.5830 | 0.5459 | -0.037 |

Alignment quality costs ~2–4% Spearman. Modest and consistent.

---

## 9. Root Cause Summary

### Issue #1: gDNA Model Over-Absorption of mRNA Reads (HIGH IMPACT)

**Magnitude**: 23–76% of mRNA reads misclassified as gDNA across regions.
**Root cause**: The gDNA model uses a flat/uniform coverage profile. Any mRNA fragment
that doesn't have strong positional evidence to distinguish it from "uniform" gets
claimed by gDNA. Single-exon fragments, fragments from long exons, and fragments from
uniformly-covered transcripts are most affected.
**Impact**: Transcript Spearman drops 0.2–0.55 from pristine in gDNA conditions.

### Issue #2: nRNA/gDNA Confusion in Combined Conditions (HIGH IMPACT)

**Magnitude**: nRNA under-estimated by 25–96% when gDNA is also present.
**Root cause**: Both nRNA and gDNA produce uniform (non-transcript-structured) read
coverage. The EM has no structural feature to distinguish them. Since gDNA is the
larger pool (10–100x), it absorbs the smaller nRNA signal.
**Impact**: nRNA pool is nearly invisible to the model when gDNA is present.

### Issue #3: Residual nRNA Leak with gDNA Only (LOW-MEDIUM IMPACT)

**Magnitude**: 131–4,602 reads attributed to nRNA when truth=0.
**Root cause**: The empirical nrna_frac prior correctly sends nrna_frac→0 globally, but some
residual assignment leaks through, particularly for EGFR (4,602 reads).
**Impact**: Small false-positive nRNA pool in gDNA-only conditions.

### Issue #4: EGFR-Specific gDNA Failure (MEDIUM IMPACT)

**Magnitude**: Only region where rigel loses to competitors with gDNA present.
**Root cause**: EGFR's dominant transcript (ENST00000275493.7) has long exons
producing reads that look uniform, identical to the gDNA profile. The gDNA model
absorbs it completely, driving the count to 0.
**Impact**: Hulkrna loses to salmon by 0.11–0.23 Spearman on EGFR.

---

## 10. Actionable Improvement Ideas

### A. Reduce gDNA Over-Absorption

The gDNA model's flat profile is too greedy. Potential fixes:
1. **Fragment-level features**: Use splice junction presence, fragment length distribution,
   or positional bias as hard evidence features that make fragments more "mRNA-like"
   and less available to the gDNA model.
2. **Prior constraint on gDNA fraction**: Cap the per-locus gDNA fraction based on genome-wide
   estimates, preventing the gDNA pool from growing unbounded.
3. **Post-EM rebalancing**: After EM convergence, check if any transcript has been
   completely zeroed out that had significant initial support; rescue those counts
   from the gDNA pool.

### B. Separate nRNA from gDNA

Both produce uniform coverage, but they differ in:
1. **Strand signal**: nRNA transcription is stranded (same strand as gene); gDNA is
   unstranded. With SS > 0.5, this is a distinguishing feature. Currently under-utilized.
2. **Per-gene structure**: nRNA reads should cluster around active genes proportional to
   expression; gDNA should be proportional to genomic length. This structural difference
   could be leveraged as a prior.
3. **Dedicated nRNA/gDNA separation step**: First estimate gDNA from intergenic regions
   (where nRNA doesn't exist), then compute the expected gDNA presence in genic regions,
   and attribute the remainder to nRNA.

### C. Improve Strand Model Utilization

The benchmark shows SS has minimal impact on rigel performance (0.003 drop from SS=1.0
to SS=0.5), suggesting strand information isn't being maximally leveraged. Better strand
modeling could help distinguish mRNA (stranded) from gDNA (unstranded).

### D. EGFR-Type Rescue

For cases where the gDNA model completely zeros out a transcript:
1. Detect when a well-supported transcript (pre-EM) goes to near-zero post-EM
2. Apply a "rescue" step that redistributes some gDNA reads back to that transcript
3. Use the splice-junction evidence as hard classification (spliced reads cannot be gDNA)
