# Benchmark Analysis Report

Generated: 2026-04-19 17:34:38

## Simulation Parameters

- RNA fragments: 10,000,000
- Fragment length: 250.0 ± 50.0
- Read length: 150
- Error rate: 0.0
- Seed: 42

## Transcript-Level Summary by Tool

### kallisto

- Conditions evaluated: 3
- Mean Spearman R: 0.7084
- Mean Pearson R: 0.9860 (range: 0.9681–0.9965)
- Mean log₂ Spearman R: 0.7084
- Mean log₂ Pearson R: 0.8967
- Mean MAPE: 259.1%
- Mean RMSE: 14.10
- Mean MAE: 2.18
- Mean WARE: 0.1957
- Count Spearman R: 0.7529
- Count MAE: 64.4
- Count RMSE: 505.1

### rigel/oracle_vbem

- Conditions evaluated: 3
- Mean Spearman R: 0.8215
- Mean Pearson R: 0.9977 (range: 0.9977–0.9977)
- Mean log₂ Spearman R: 0.8215
- Mean log₂ Pearson R: 0.9753
- Mean MAPE: 74.3%
- Mean RMSE: 4.63
- Mean MAE: 0.71
- Mean WARE: 0.0769
- Count Spearman R: 0.8125
- Count MAE: 61.4
- Count RMSE: 501.8

### rigel/star_vbem

- Conditions evaluated: 3
- Mean Spearman R: 0.8162
- Mean Pearson R: 0.9891 (range: 0.9854–0.9910)
- Mean log₂ Spearman R: 0.8162
- Mean log₂ Pearson R: 0.9711
- Mean MAPE: 78.5%
- Mean RMSE: 10.11
- Mean MAE: 0.92
- Mean WARE: 0.0991
- Count Spearman R: 0.8083
- Count MAE: 61.6
- Count RMSE: 506.4

### salmon

- Conditions evaluated: 3
- Mean Spearman R: 0.7509
- Mean Pearson R: 0.9065 (range: 0.9045–0.9075)
- Mean log₂ Spearman R: 0.7509
- Mean log₂ Pearson R: 0.8913
- Mean MAPE: 204.7%
- Mean RMSE: 28.76
- Mean MAE: 2.92
- Mean WARE: 0.3260
- Count Spearman R: 0.7505
- Count MAE: 69.6
- Count RMSE: 550.5

## Transcript-Level Metrics (All Conditions)

| condition | tool | gdna_label | strand_specificity | nrna_label | spearman_r | pearson_r | log2_spearman_r | log2_pearson_r | rmse | mae | mape | ware | precision | recall | f1 | count_spearman_r | count_pearson_r | count_mae | count_rmse |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pristine | rigel/oracle_vbem | unknown | 1.0000 | none | 0.8382 | 0.9977 | 0.8382 | 0.9765 | 4.5491 | 0.6959 | 72.0690 | 0.0759 | 0.8562 | 0.7887 | 0.8211 | 0.8246 | 0.6672 | 60.9574 | 499.2 |
| pristine | rigel/star_vbem | unknown | 1.0000 | none | 0.8340 | 0.9910 | 0.8340 | 0.9725 | 9.2374 | 0.8866 | 75.5554 | 0.0958 | 0.8534 | 0.7873 | 0.8190 | 0.8215 | 0.6603 | 61.1365 | 503.3 |
| pristine | salmon | unknown | 1.0000 | none | 0.8002 | 0.9075 | 0.8002 | 0.8993 | 28.3943 | 2.7597 | 173.7 | 0.3105 | 0.9376 | 0.6712 | 0.7823 | 0.7966 | 0.5772 | 68.6196 | 547.6 |
| pristine | kallisto | unknown | 1.0000 | none | 0.8577 | 0.9965 | 0.8577 | 0.9789 | 5.5982 | 0.6419 | 64.4717 | 0.0706 | 0.9277 | 0.7405 | 0.8236 | 0.8501 | 0.6654 | 61.3781 | 500.3 |
| gdna | rigel/oracle_vbem | unknown | 1.0000 | none | 0.8129 | 0.9977 | 0.8129 | 0.9754 | 4.6286 | 0.7034 | 76.5162 | 0.0758 | 0.8135 | 0.7854 | 0.7992 | 0.8060 | 0.6668 | 61.1686 | 499.4 |
| gdna | rigel/star_vbem | unknown | 1.0000 | none | 0.8080 | 0.9908 | 0.8080 | 0.9712 | 9.3660 | 0.8947 | 81.3769 | 0.0964 | 0.8093 | 0.7841 | 0.7965 | 0.8022 | 0.6595 | 61.3892 | 503.7 |
| gdna | salmon | unknown | 1.0000 | none | 0.7274 | 0.9074 | 0.7274 | 0.8906 | 28.5412 | 2.8660 | 219.0 | 0.3188 | 0.7934 | 0.7271 | 0.7588 | 0.7284 | 0.5765 | 69.2469 | 547.9 |
| gdna | kallisto | unknown | 1.0000 | none | 0.6410 | 0.9934 | 0.6410 | 0.8752 | 12.0770 | 2.1241 | 352.0 | 0.1776 | 0.6255 | 0.8662 | 0.7265 | 0.7057 | 0.6639 | 64.3357 | 501.2 |
| gdna_nrna | rigel/oracle_vbem | nrna | 1.0000 | none | 0.8134 | 0.9977 | 0.8134 | 0.9740 | 4.7029 | 0.7288 | 74.1796 | 0.0789 | 0.8158 | 0.7865 | 0.8008 | 0.8070 | 0.6637 | 62.0153 | 506.8 |
| gdna_nrna | rigel/star_vbem | nrna | 1.0000 | none | 0.8067 | 0.9854 | 0.8067 | 0.9696 | 11.7372 | 0.9738 | 78.5915 | 0.1050 | 0.8101 | 0.7831 | 0.7964 | 0.8012 | 0.6543 | 62.2236 | 512.3 |
| gdna_nrna | salmon | nrna | 1.0000 | none | 0.7251 | 0.9045 | 0.7251 | 0.8839 | 29.3321 | 3.1218 | 221.3 | 0.3488 | 0.7888 | 0.7327 | 0.7597 | 0.7263 | 0.5716 | 70.8079 | 556.2 |
| gdna_nrna | kallisto | nrna | 1.0000 | none | 0.6264 | 0.9681 | 0.6264 | 0.8361 | 24.6350 | 3.7663 | 360.9 | 0.3390 | 0.6173 | 0.8777 | 0.7248 | 0.7028 | 0.6520 | 67.6218 | 513.7 |

## Gene-Level Metrics (All Conditions)

| condition | tool | gdna_label | strand_specificity | nrna_label | spearman_r | pearson_r | log2_spearman_r | log2_pearson_r | rmse | mae | mape | ware |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pristine | rigel/oracle_vbem | unknown | 1.0000 | none | 0.8992 | 0.9996 | 0.8992 | 0.9931 | 3.1049 | 0.6036 | 56.9947 | 0.0374 |
| pristine | rigel/star_vbem | unknown | 1.0000 | none | 0.8952 | 0.9949 | 0.8952 | 0.9908 | 10.8198 | 0.8929 | 58.0826 | 0.0538 |
| pristine | salmon | unknown | 1.0000 | none | 0.8559 | 0.9794 | 0.8559 | 0.9727 | 23.1487 | 2.4216 | 73.2329 | 0.1512 |
| pristine | kallisto | unknown | 1.0000 | none | 0.8992 | 0.9984 | 0.8992 | 0.9927 | 5.8667 | 0.6080 | 56.7213 | 0.0366 |
| gdna | rigel/oracle_vbem | unknown | 1.0000 | none | 0.8435 | 0.9996 | 0.8435 | 0.9902 | 3.3969 | 0.6508 | 77.5996 | 0.0384 |
| gdna | rigel/star_vbem | unknown | 1.0000 | none | 0.8391 | 0.9948 | 0.8391 | 0.9878 | 10.9329 | 0.9363 | 78.3068 | 0.0552 |
| gdna | salmon | unknown | 1.0000 | none | 0.7584 | 0.9793 | 0.7584 | 0.9629 | 23.7794 | 2.6553 | 168.5 | 0.1621 |
| gdna | kallisto | unknown | 1.0000 | none | 0.6310 | 0.9960 | 0.6310 | 0.8804 | 17.3591 | 3.5719 | 667.9 | 0.1584 |
| gdna_nrna | rigel/oracle_vbem | nrna | 1.0000 | none | 0.8478 | 0.9996 | 0.8478 | 0.9887 | 3.7183 | 0.7036 | 72.8472 | 0.0417 |
| gdna_nrna | rigel/star_vbem | nrna | 1.0000 | none | 0.8423 | 0.9906 | 0.8423 | 0.9862 | 14.6086 | 1.0833 | 72.9862 | 0.0635 |
| gdna_nrna | salmon | nrna | 1.0000 | none | 0.7590 | 0.9773 | 0.7590 | 0.9583 | 26.1852 | 3.1273 | 157.4 | 0.1913 |
| gdna_nrna | kallisto | nrna | 1.0000 | none | 0.6150 | 0.9784 | 0.6150 | 0.8478 | 34.7247 | 6.2964 | 591.8 | 0.2982 |

## Pool-Level Summary (Rigel)

| condition | tool | mrna_frag_truth | mrna_pred | mrna_rel_error | nrna_frag_truth | nrna_pred | nrna_rel_error | gdna_frag_truth | gdna_pred | gdna_rel_error |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pristine | rigel/oracle_vbem | 10000000 | 9986411.0 | -0.0014 | 0 | 13401.0 | — | 0 | 188.0 | — |
| pristine | rigel/star_vbem | 10000000 | 9978719.0 | -0.0021 | 0 | 18296.0 | — | 0 | 1529.0 | — |
| gdna | rigel/oracle_vbem | 10000000 | 10013655.0 | 0.0014 | 0 | 182271.0 | — | 10000000 | 9802855.0 | -0.0197 |
| gdna | rigel/star_vbem | 10000000 | 10005122.0 | 0.0005 | 0 | 195482.0 | — | 10000000 | 9214704.0 | -0.0785 |
| gdna_nrna | rigel/oracle_vbem | 10139172 | 10073762.0 | -0.0065 | 4860828 | 5014198.0 | 0.0316 | 10000000 | 9910794.0 | -0.0089 |
| gdna_nrna | rigel/star_vbem | 10139172 | 10068418.0 | -0.0070 | 4860828 | 5029997.0 | 0.0348 | 10000000 | 9314676.0 | -0.0685 |

## Stratified Metrics by Expression Level

### kallisto

| expression_bin | mean_spearman | mean_pearson | mean_mape | mean_mae | mean_ware | n_conditions |
| --- | --- | --- | --- | --- | --- | --- |
| high_100-1000 | 0.9749 | 0.9889 | 16.2563 | 38.1681 | 0.1627 | 3 |
| low_1-10 | 0.7535 | 0.6315 | 43.6430 | 1.2963 | 0.3462 | 3 |
| mid_10-100 | 0.9549 | 0.9626 | 19.0800 | 5.5745 | 0.1800 | 3 |
| very_high_1000+ | 0.9735 | 0.9803 | 16.8702 | 312.1 | 0.1676 | 3 |
| zero | — | — | — | 1.0296 | — | 3 |

### rigel/oracle_vbem

| expression_bin | mean_spearman | mean_pearson | mean_mape | mean_mae | mean_ware | n_conditions |
| --- | --- | --- | --- | --- | --- | --- |
| high_100-1000 | 0.9852 | 0.9946 | 4.1139 | 8.5619 | 0.0365 | 3 |
| low_1-10 | 0.8125 | 0.7559 | 37.1215 | 1.0948 | 0.2924 | 3 |
| mid_10-100 | 0.9645 | 0.9774 | 10.8942 | 2.6999 | 0.0871 | 3 |
| very_high_1000+ | 0.9858 | 0.9907 | 2.6867 | 47.1256 | 0.0253 | 3 |
| zero | — | — | — | 0.1205 | — | 3 |

### rigel/star_vbem

| expression_bin | mean_spearman | mean_pearson | mean_mape | mean_mae | mean_ware | n_conditions |
| --- | --- | --- | --- | --- | --- | --- |
| high_100-1000 | 0.9720 | 0.9615 | 5.6243 | 13.2974 | 0.0567 | 3 |
| low_1-10 | 0.8029 | 0.4523 | 41.1536 | 1.2012 | 0.3208 | 3 |
| mid_10-100 | 0.9583 | 0.9597 | 11.8841 | 3.0210 | 0.0975 | 3 |
| very_high_1000+ | 0.9326 | 0.9467 | 7.1823 | 131.3 | 0.0705 | 3 |
| zero | — | — | — | 0.1631 | — | 3 |

### salmon

| expression_bin | mean_spearman | mean_pearson | mean_mape | mean_mae | mean_ware | n_conditions |
| --- | --- | --- | --- | --- | --- | --- |
| high_100-1000 | 0.7352 | 0.8171 | 23.7880 | 56.2875 | 0.2402 | 3 |
| low_1-10 | 0.5657 | 0.0998 | 97.9163 | 2.8787 | 0.7688 | 3 |
| mid_10-100 | 0.6884 | 0.4548 | 35.6558 | 10.0155 | 0.3232 | 3 |
| very_high_1000+ | 0.6878 | 0.6671 | 25.0897 | 468.4 | 0.2513 | 3 |
| zero | — | — | — | 0.3290 | — | 3 |


## Rigel Calibration & Model Details

| condition | tool | strand_specificity | strand_protocol | lambda_gdna | mixing_pi | mixing_pi_soft | strand_used | strand_z | mu_R | sigma_R | gdna_fraction | gdna_fl_mean | gdna_fl_observations | n_eligible | em_n_iter | calibration_converged | fl_global_mean | fl_rna_mean | fl_intergenic_mean | n_loci | n_unambig | n_em | mrna_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pristine | rigel/oracle_vbem | 1.0000 | R1-antisense | 0.0000 | 0.8387 | 0.9673 | True | 2053.6 | -4.7763 | 1.8452 | 0.0000 | 187.0 | 4672 | 400018 | 5 | True | 248.7 | 251.1 | 0.0000 | 19028 | 669172 | 9330828 | 0.9986 |
| pristine | rigel/star_vbem | 1.0000 | R1-antisense | 0.0000 | 0.7681 | 0.9609 | True | 2051.6 | -5.5029 | 1.8179 | 0.0000 | 168.7 | 3297 | 400000 | 6 | True | 248.8 | 251.2 | 278.0 | 18037 | 649758 | 9348783 | 0.9980 |
| gdna | rigel/oracle_vbem | 1.0000 | R1-antisense | 0.0033 | 0.8651 | 0.9872 | True | 1366.7 | -4.7187 | 1.9146 | 0.4434 | 350.1 | 5215503 | 375329 | 6 | True | 299.9 | 251.1 | 351.0 | 34250 | 669172 | 15125027 | 0.5007 |
| gdna | rigel/star_vbem | 1.0000 | R1-antisense | 0.0033 | 0.7883 | 0.9837 | True | 1365.7 | -5.5402 | 1.8705 | 0.4540 | 350.2 | 4677081 | 375943 | 7 | True | 298.2 | 251.2 | 351.0 | 29453 | 649763 | 15153876 | 0.5153 |
| gdna_nrna | rigel/oracle_vbem | 1.0000 | R1-antisense | 0.0037 | 0.7194 | 0.8226 | True | 2328.2 | -4.2534 | 1.7649 | 0.3811 | 349.9 | 5064112 | 380016 | 6 | True | 289.5 | 251.1 | 351.1 | 33986 | 662457 | 20130888 | 0.4030 |
| gdna_nrna | rigel/star_vbem | 1.0000 | R1-antisense | 0.0036 | 0.6900 | 0.8629 | True | 2326.4 | -4.4885 | 1.6523 | 0.3866 | 349.9 | 4563138 | 380867 | 4 | True | 287.7 | 251.2 | 351.1 | 28483 | 643030 | 20158098 | 0.4124 |

## Cross-Tool Comparison

### Transcript-Level Spearman R by Condition

| condition | kallisto | rigel/oracle_vbem | rigel/star_vbem | salmon |
| --- | --- | --- | --- | --- |
| gdna | 0.6410 | 0.8129 | 0.8080 | 0.7274 |
| gdna_nrna | 0.6264 | 0.8134 | 0.8067 | 0.7251 |
| pristine | 0.8577 | 0.8382 | 0.8340 | 0.8002 |

### Transcript-Level Pearson R by Condition

| condition | kallisto | rigel/oracle_vbem | rigel/star_vbem | salmon |
| --- | --- | --- | --- | --- |
| gdna | 0.9934 | 0.9977 | 0.9908 | 0.9074 |
| gdna_nrna | 0.9681 | 0.9977 | 0.9854 | 0.9045 |
| pristine | 0.9965 | 0.9977 | 0.9910 | 0.9075 |

### Transcript-Level MAPE by Condition

| condition | kallisto | rigel/oracle_vbem | rigel/star_vbem | salmon |
| --- | --- | --- | --- | --- |
| gdna | 352.0 | 76.5162 | 81.3769 | 219.0 |
| gdna_nrna | 360.9 | 74.1796 | 78.5915 | 221.3 |
| pristine | 64.4717 | 72.0690 | 75.5554 | 173.7 |

### Transcript-Level WARE by Condition

| condition | kallisto | rigel/oracle_vbem | rigel/star_vbem | salmon |
| --- | --- | --- | --- | --- |
| gdna | 0.1776 | 0.0758 | 0.0964 | 0.3188 |
| gdna_nrna | 0.3390 | 0.0789 | 0.1050 | 0.3488 |
| pristine | 0.0706 | 0.0759 | 0.0958 | 0.3105 |


## Sensitivity Analysis

### kallisto

#### By gDNA Contamination Level

| gdna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| nrna | 0.9681 | 360.9 | 0.3390 | 24.6350 | 1 |
| unknown | 0.9949 | 208.2 | 0.1241 | 8.8376 | 2 |

#### By Strand Specificity

| strand_specificity | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| 1.0000 | 0.9860 | 259.1 | 0.1957 | 14.1034 | 3 |

#### By nRNA Contamination

| nrna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| none | 0.9860 | 259.1 | 0.1957 | 14.1034 | 3 |

### rigel/oracle_vbem

#### By gDNA Contamination Level

| gdna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| nrna | 0.9977 | 74.1796 | 0.0789 | 4.7029 | 1 |
| unknown | 0.9977 | 74.2926 | 0.0758 | 4.5888 | 2 |

#### By Strand Specificity

| strand_specificity | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| 1.0000 | 0.9977 | 74.2549 | 0.0769 | 4.6269 | 3 |

#### By nRNA Contamination

| nrna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| none | 0.9977 | 74.2549 | 0.0769 | 4.6269 | 3 |

### rigel/star_vbem

#### By gDNA Contamination Level

| gdna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| nrna | 0.9854 | 78.5915 | 0.1050 | 11.7372 | 1 |
| unknown | 0.9909 | 78.4662 | 0.0961 | 9.3017 | 2 |

#### By Strand Specificity

| strand_specificity | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| 1.0000 | 0.9891 | 78.5080 | 0.0991 | 10.1135 | 3 |

#### By nRNA Contamination

| nrna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| none | 0.9891 | 78.5080 | 0.0991 | 10.1135 | 3 |

### salmon

#### By gDNA Contamination Level

| gdna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| nrna | 0.9045 | 221.3 | 0.3488 | 29.3321 | 1 |
| unknown | 0.9074 | 196.4 | 0.3147 | 28.4677 | 2 |

#### By Strand Specificity

| strand_specificity | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| 1.0000 | 0.9065 | 204.7 | 0.3260 | 28.7559 | 3 |

#### By nRNA Contamination

| nrna_label | mean_pearson | mean_mape | mean_ware | mean_rmse | n |
| --- | --- | --- | --- | --- | --- |
| none | 0.9065 | 204.7 | 0.3260 | 28.7559 | 3 |

