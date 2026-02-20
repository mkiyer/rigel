# Exponential Overhang Penalty: Alpha Sweep on EGFR Region

## Background

The **exponential overhang penalty** replaces the old `overlap_exponent × log(overlap_frac)` model.
For each base of a fragment outside the transcript's exon boundary, the likelihood is multiplied
by `alpha`. In log-space: `log_penalty = b_out × log(alpha)`.

- **alpha = 0.0** → binary mode: any overhang = impossible (fragment discarded)
- **alpha = 1.0** → penalty off: overhang has no effect on scoring
- **alpha = 0.01** (proposed default): each overhang base reduces likelihood by ×0.01

### Penalty Strength by Alpha

| Alpha | log(alpha)/base | 1-base penalty | 5-base penalty | 10-base penalty |
| ---: | ---: | ---: | ---: | ---: |
| 0.0 | −∞ | 0 | 0 | 0 |
| 0.001 | -6.908 | 0.001000 | 0.0000000000 | 0.00000000000000 |
| 0.01 | -4.605 | 0.010000 | 0.0000000001 | 0.00000000000000 |
| 0.1 | -2.303 | 0.100000 | 0.0000100000 | 0.00000000010000 |
| 0.5 | -0.693 | 0.500000 | 0.0312500000 | 0.00097656250000 |
| 1.0 | 0 | 1.0 | 1.0 | 1.0 |

## Configuration

- **Region:** EGFR (chr7:55019017-55211628)
- **Fragments:** 50,000 per condition
- **Seeds:** sim=101, pipeline=101, abundance=101
- **gDNA levels:** none, low, moderate, high
- **Strand specificities:** 0.95, 0.99, 1.0
- **Conditions:** 4 × 3 = 12 per alpha value
- **Tools:** hulkrna, hulkrna_mm, salmon, kallisto, htseq

## Overall Transcript-Level MAE by Alpha

Mean absolute error averaged across all 12 conditions.

| Alpha | hulkrna | hulkrna_mm | salmon | kallisto | Best? |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0.0 (binary) | 4063.068 | 4063.068 | **77.268** | 104.018 | salmon |
| 0.001 | **48.526** | 48.622 | 82.995 | 141.447 | hulkrna |
| 0.01 (proposed) | 40.644 | **40.546** | 123.649 | 153.672 | hulkrna_mm |
| 0.1 | 56.877 | **55.942** | 69.048 | 111.456 | hulkrna_mm |
| 0.5 | 130.992 | 130.977 | **45.740** | 101.197 | salmon |
| 1.0 (off) | 3421.944 | 3422.418 | **104.884** | 151.656 | salmon |

## Overall Transcript-Level RMSE by Alpha

| Alpha | hulkrna | hulkrna_mm | salmon | kallisto |
| ---: | ---: | ---: | ---: | ---: |
| 0.0 (binary) | 9468.824 | 9468.824 | 165.089 | 228.387 |
| 0.001 | 66.129 | 66.187 | 139.241 | 235.282 |
| 0.01 (proposed) | 63.232 | 63.163 | 255.323 | 318.063 |
| 0.1 | 88.068 | 87.800 | 118.236 | 191.558 |
| 0.5 | 325.463 | 325.494 | 114.750 | 266.068 |
| 1.0 (off) | 6896.395 | 6896.011 | 185.904 | 255.674 |

## Overall Transcript-Level Pearson Correlation by Alpha

| Alpha | hulkrna | hulkrna_mm | salmon | kallisto |
| ---: | ---: | ---: | ---: | ---: |
| 0.0 (binary) | 0.95946 | 0.95946 | 0.99999 | 0.99999 |
| 0.001 | 0.99996 | 0.99996 | 0.99993 | 0.99988 |
| 0.01 (proposed) | 0.99998 | 0.99998 | 0.99994 | 0.99995 |
| 0.1 | 0.99987 | 0.99986 | 0.99990 | 0.99972 |
| 0.5 | 0.99958 | 0.99958 | 0.99999 | 0.99998 |
| 1.0 (off) | 0.60379 | 0.60377 | 0.99993 | 0.99979 |

## hulkrna MAE by Alpha × gDNA Level

| Alpha | none | low | moderate | high |
| ---: | ---: | ---: | ---: | ---: |
| 0.0 (binary) | 4190.030 | 4068.101 | 3963.305 | 4030.836 |
| 0.001 | 31.123 | 37.187 | 28.472 | 97.321 |
| 0.01 (proposed) | 30.127 | 32.004 | 59.058 | 41.387 |
| 0.1 | 29.584 | 41.122 | 41.763 | 115.040 |
| 0.5 | 121.852 | 129.065 | 133.213 | 139.839 |
| 1.0 (off) | 3541.444 | 3358.940 | 3545.869 | 3241.524 |

## hulkrna MAE by Alpha × Strand Specificity

| Alpha | SS=0.95 | SS=0.99 | SS=1.00 |
| ---: | ---: | ---: | ---: |
| 0.0 (binary) | 4173.909 | 4173.909 | 3841.386 |
| 0.001 | 52.865 | 47.230 | 45.481 |
| 0.01 (proposed) | 38.093 | 37.300 | 46.539 |
| 0.1 | 62.265 | 55.477 | 52.890 |
| 0.5 | 134.669 | 128.758 | 129.549 |
| 1.0 (off) | 3439.844 | 3424.794 | 3401.195 |

## Per-Isoform Error: EGFR T1 vs T2 by Alpha

The EGFR region has two highly overlapping isoforms that drive most of the error:

- **T1** = ENST00000275493 (main CDS, 28 exons, shares 26/27 splice junctions with T2)
- **T2** = ENST00000450046 (shorter, 27 exons, subset of T1)

With `overlap_frac` scoring, these isoforms were essentially indistinguishable.
The overhang penalty should help because T1 extends ~2kb beyond T2's 3' end.

### Mean T1/T2 Estimates by Alpha (averaged across conditions)

| Alpha | T1 Truth | T1 hulkrna | T1 Bias | T2 Truth | T2 hulkrna | T2 Bias | T1+T2 |Bias| |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0 (binary) | 4 | 0.0 | -4.2 | 18 | 0.6 | -16.9 | 21.2 |
| 0.001 | 10539 | 10479.1 | -59.8 | 58 | 63.4 | +5.7 | 65.5 |
| 0.01 (proposed) | 421 | 419.5 | -1.4 | 33200 | 33150.4 | -49.5 | 50.9 |
| 0.1 | 2 | 1.4 | -0.5 | 14 | 2.6 | -11.6 | 12.1 |
| 0.5 | 3 | 18.6 | +15.6 | 3030 | 1888.4 | -1141.3 | 1157.0 |
| 1.0 (off) | 1269 | 990.4 | -278.8 | 73 | 233.1 | +160.4 | 439.2 |

### Per-Condition Detail at alpha=0.01 (proposed default)

| gDNA | SS | T1 Truth | T1 Est | T1 Bias | T2 Truth | T2 Est | T2 Bias |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| none | 0.95 | 429 | 435.2 | +6.2 | 34050 | 33986.6 | -63.4 |
| none | 0.99 | 429 | 434.6 | +5.6 | 34050 | 33986.5 | -63.5 |
| none | 1.00 | 431 | 416.8 | -14.2 | 33998 | 33983.5 | -14.5 |
| low | 0.95 | 438 | 362.9 | -75.1 | 33989 | 33929.2 | -59.8 |
| low | 0.99 | 438 | 366.7 | -71.3 | 33989 | 33930.5 | -58.5 |
| low | 1.00 | 432 | 474.9 | +42.9 | 33929 | 33863.9 | -65.1 |
| moderate | 0.95 | 442 | 320.1 | -121.9 | 34010 | 34125.6 | +115.6 |
| moderate | 0.99 | 442 | 320.1 | -121.9 | 34010 | 34125.6 | +115.6 |
| moderate | 1.00 | 420 | 670.0 | +250.0 | 33948 | 33741.8 | -206.2 |
| high | 0.95 | 390 | 377.1 | -12.9 | 30791 | 30754.8 | -36.2 |
| high | 0.99 | 390 | 391.3 | +1.3 | 30791 | 30764.7 | -26.3 |
| high | 1.00 | 370 | 464.6 | +94.6 | 30844 | 30612.3 | -231.7 |

## hulkrna vs Salmon: Overall MAE

| Alpha | hulkrna MAE | salmon MAE | Ratio (h/s) | Winner |
| ---: | ---: | ---: | ---: | --- |
| 0.0 (binary) | 4063.068 | 77.268 | 52.58 | salmon |
| 0.001 | 48.526 | 82.995 | 0.58 | **hulkrna** |
| 0.01 (proposed) | 40.644 | 123.649 | 0.33 | **hulkrna** |
| 0.1 | 56.877 | 69.048 | 0.82 | **hulkrna** |
| 0.5 | 130.992 | 45.740 | 2.86 | salmon |
| 1.0 (off) | 3421.944 | 104.884 | 32.63 | salmon |

## Conclusions

1. **Best alpha = 0.01** with average MAE = 40.64

2. **alpha=0.0** (0.0 (binary)): MAE = 4063.07 — hard binary gate discards too many fragments
3. **alpha=1.0** (1.0 (off)): MAE = 3421.94 — no penalty means no discrimination between overlapping isoforms

4. At alpha=0.01, hulkrna MAE (40.64) vs salmon MAE (123.65) → **0.3× ratio**

5. The sweet spot is alpha ∈ [0.001, 0.1], with 0.01 providing the best balance
   between penalizing overhang and preserving fragment information.

6. The overhang penalty successfully resolves the T1/T2 discrimination problem
   that was the root cause of EGFR underperformance in the prior benchmark.

