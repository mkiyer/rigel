# Benchmark Analysis Report

Generated: 2026-04-03 11:23:36

## Simulation Parameters

- RNA fragments: 50,000,000
- Fragment length: 250.0 ± 50.0
- Read length: 101
- Error rate: 0.0
- Seed: 42

## Transcript-Level Summary by Tool

### kallisto

- Conditions evaluated: 1
- Mean Pearson R: 0.9958 (range: 0.9958–0.9958)
- Mean Spearman R: 0.7319
- Mean MAPE: 194.5%
- Mean RMSE: 3.98
- Mean MAE: 0.64
- Mean WARE: 0.1162

### rigel/map

- Conditions evaluated: 1
- Mean Pearson R: 0.9862 (range: 0.9862–0.9862)
- Mean Spearman R: 0.8753
- Mean MAPE: 56.8%
- Mean RMSE: 5.93
- Mean MAE: 0.37
- Mean WARE: 0.0796

### rigel/vbem

- Conditions evaluated: 1
- Mean Pearson R: 0.9865 (range: 0.9865–0.9865)
- Mean Spearman R: 0.8764
- Mean MAPE: 58.1%
- Mean RMSE: 5.87
- Mean MAE: 0.37
- Mean WARE: 0.0816

### salmon

- Conditions evaluated: 1
- Mean Pearson R: 0.8900 (range: 0.8900–0.8900)
- Mean Spearman R: 0.7564
- Mean MAPE: 239.1%
- Mean RMSE: 16.31
- Mean MAE: 1.63
- Mean WARE: 0.3627

## Transcript-Level Metrics (All Conditions)

| condition | tool | gdna_label | strand_specificity | nrna_label | pearson_r | spearman_r | rmse | mae | mape | ware | precision | recall | f1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gdna_high_ss_0.90_nrna_none | rigel/map | high | 0.9000 | none | 0.9862 | 0.8753 | 5.9275 | 0.3698 | 56.8111 | 0.0796 | 0.8923 | 0.8392 | 0.8649 |
| gdna_high_ss_0.90_nrna_none | rigel/vbem | high | 0.9000 | none | 0.9865 | 0.8764 | 5.8651 | 0.3701 | 58.0690 | 0.0816 | 0.9234 | 0.8060 | 0.8607 |
| gdna_high_ss_0.90_nrna_none | salmon | high | 0.9000 | none | 0.8900 | 0.7564 | 16.3138 | 1.6274 | 239.1 | 0.3627 | 0.7518 | 0.8510 | 0.7983 |
| gdna_high_ss_0.90_nrna_none | kallisto | high | 0.9000 | none | 0.9958 | 0.7319 | 3.9820 | 0.6365 | 194.5 | 0.1162 | 0.6298 | 0.9285 | 0.7505 |

## Gene-Level Metrics (All Conditions)

| condition | tool | gdna_label | strand_specificity | nrna_label | pearson_r | spearman_r | rmse | mae | mape | ware |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gdna_high_ss_0.90_nrna_none | rigel/map | high | 0.9000 | none | 0.9952 | 0.8825 | 9.2685 | 0.6640 | 70.5584 | 0.0346 |
| gdna_high_ss_0.90_nrna_none | rigel/vbem | high | 0.9000 | none | 0.9954 | 0.8865 | 9.0504 | 0.6352 | 67.5900 | 0.0340 |
| gdna_high_ss_0.90_nrna_none | salmon | high | 0.9000 | none | 0.9880 | 0.7785 | 14.6084 | 1.9505 | 190.6 | 0.1182 |
| gdna_high_ss_0.90_nrna_none | kallisto | high | 0.9000 | none | 0.9984 | 0.7010 | 8.1945 | 1.7387 | 781.8 | 0.0805 |

## Pool-Level Summary (Rigel)

| condition | tool | mrna_frag_truth | mrna_pred | mrna_rel_error | nrna_frag_truth | nrna_pred | nrna_rel_error | gdna_frag_truth | gdna_pred | gdna_rel_error |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gdna_high_ss_0.90_nrna_none | rigel/map | 50000000 | 49970646.0 | -0.0006 | 0 | 2186473.0 | — | 25000000 | 21309015.0 | -0.1476 |
| gdna_high_ss_0.90_nrna_none | rigel/vbem | 50000000 | 49970443.0 | -0.0006 | 0 | 2136992.0 | — | 25000000 | 21358699.0 | -0.1457 |

## Stratified Metrics by Expression Level

### kallisto

| expression_bin | mean_pearson | mean_spearman | mean_mape | mean_mae | mean_ware | n_conditions |
| --- | --- | --- | --- | --- | --- | --- |
| high_100-1000 | 0.9925 | 0.9782 | 6.8351 | 15.6056 | 0.0689 | 1 |
| low_1-10 | 0.8484 | 0.8830 | 24.5080 | 0.6660 | 0.1944 | 1 |
| mid_10-100 | 0.9856 | 0.9704 | 9.3938 | 2.3892 | 0.0851 | 1 |
| very_high_1000+ | 0.9765 | 0.9805 | 8.3400 | 124.2 | 0.0790 | 1 |
| zero | — | — | — | 0.3544 | — | 1 |

### rigel/map

| expression_bin | mean_pearson | mean_spearman | mean_mape | mean_mae | mean_ware | n_conditions |
| --- | --- | --- | --- | --- | --- | --- |
| high_100-1000 | 0.9773 | 0.9780 | 3.5639 | 7.5250 | 0.0332 | 1 |
| low_1-10 | 0.8033 | 0.8923 | 23.5841 | 0.6298 | 0.1838 | 1 |
| mid_10-100 | 0.9823 | 0.9731 | 7.0313 | 1.5946 | 0.0568 | 1 |
| very_high_1000+ | 0.9096 | 0.9320 | 5.2757 | 74.5934 | 0.0475 | 1 |
| zero | — | — | — | 0.1125 | — | 1 |

### rigel/vbem

| expression_bin | mean_pearson | mean_spearman | mean_mape | mean_mae | mean_ware | n_conditions |
| --- | --- | --- | --- | --- | --- | --- |
| high_100-1000 | 0.9773 | 0.9775 | 3.6194 | 7.5580 | 0.0334 | 1 |
| low_1-10 | 0.8025 | 0.8810 | 24.5563 | 0.6549 | 0.1912 | 1 |
| mid_10-100 | 0.9813 | 0.9713 | 7.2147 | 1.6286 | 0.0580 | 1 |
| very_high_1000+ | 0.9094 | 0.9316 | 5.2715 | 74.4308 | 0.0474 | 1 |
| zero | — | — | — | 0.0973 | — | 1 |

### salmon

| expression_bin | mean_pearson | mean_spearman | mean_mape | mean_mae | mean_ware | n_conditions |
| --- | --- | --- | --- | --- | --- | --- |
| high_100-1000 | 0.8127 | 0.7380 | 21.9045 | 50.3478 | 0.2222 | 1 |
| low_1-10 | 0.1520 | 0.6353 | 81.1249 | 2.2452 | 0.6554 | 1 |
| mid_10-100 | 0.5490 | 0.7012 | 32.6305 | 8.2260 | 0.2931 | 1 |
| very_high_1000+ | 0.5634 | 0.6434 | 22.0131 | 354.6 | 0.2257 | 1 |
| zero | — | — | — | 0.3980 | — | 1 |


## Rigel Calibration & Model Details

| condition | tool | strand_specificity | strand_protocol | gdna_density_global | gdna_mixing_prop | kappa_strand | calibration_converged | fl_global_mean | fl_rna_mean | fl_intergenic_mean | n_loci | n_unambig | n_em | mrna_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gdna_high_ss_0.90_nrna_none | rigel/map | 0.9000 | R1-antisense | 0.0105 | 0.4863 | 4.9600 | True | — | — | — | 9100 | 2348512 | 62930558 | 0.6802 |
| gdna_high_ss_0.90_nrna_none | rigel/vbem | 0.9000 | R1-antisense | 0.0105 | 0.4863 | 4.9600 | True | — | — | — | 9100 | 2348512 | 62930558 | 0.6802 |

## Cross-Tool Comparison

### Transcript-Level Pearson R by Condition

| condition | kallisto | rigel/map | rigel/vbem | salmon |
| --- | --- | --- | --- | --- |
| gdna_high_ss_0.90_nrna_none | 0.9958 | 0.9862 | 0.9865 | 0.8900 |

### Transcript-Level MAPE by Condition

| condition | kallisto | rigel/map | rigel/vbem | salmon |
| --- | --- | --- | --- | --- |
| gdna_high_ss_0.90_nrna_none | 194.5 | 56.8111 | 58.0690 | 239.1 |

### Transcript-Level WARE by Condition

| condition | kallisto | rigel/map | rigel/vbem | salmon |
| --- | --- | --- | --- | --- |
| gdna_high_ss_0.90_nrna_none | 0.1162 | 0.0796 | 0.0816 | 0.3627 |


## Sensitivity Analysis

### kallisto

#### By gDNA Contamination Level

| gdna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| high | 0.9958 | 194.5 | 0.1162 | 3.9820 | 1 |

#### By Strand Specificity

| strand_specificity | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| 0.9000 | 0.9958 | 194.5 | 0.1162 | 3.9820 | 1 |

#### By nRNA Contamination

| nrna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| none | 0.9958 | 194.5 | 0.1162 | 3.9820 | 1 |

### rigel/map

#### By gDNA Contamination Level

| gdna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| high | 0.9862 | 56.8111 | 0.0796 | 5.9275 | 1 |

#### By Strand Specificity

| strand_specificity | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| 0.9000 | 0.9862 | 56.8111 | 0.0796 | 5.9275 | 1 |

#### By nRNA Contamination

| nrna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| none | 0.9862 | 56.8111 | 0.0796 | 5.9275 | 1 |

### rigel/vbem

#### By gDNA Contamination Level

| gdna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| high | 0.9865 | 58.0690 | 0.0816 | 5.8651 | 1 |

#### By Strand Specificity

| strand_specificity | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| 0.9000 | 0.9865 | 58.0690 | 0.0816 | 5.8651 | 1 |

#### By nRNA Contamination

| nrna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| none | 0.9865 | 58.0690 | 0.0816 | 5.8651 | 1 |

### salmon

#### By gDNA Contamination Level

| gdna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| high | 0.8900 | 239.1 | 0.3627 | 16.3138 | 1 |

#### By Strand Specificity

| strand_specificity | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| 0.9000 | 0.8900 | 239.1 | 0.3627 | 16.3138 | 1 |

#### By nRNA Contamination

| nrna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| none | 0.8900 | 239.1 | 0.3627 | 16.3138 | 1 |

