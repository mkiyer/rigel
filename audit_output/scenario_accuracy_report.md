# Scenario Accuracy Audit Report

## Overall Summary

- **Total transcript observations**: 406
- **Mean abs error**: 21.0
- **Median abs error**: 9.0
- **Max abs error**: 338.0
- **90th pctile abs error**: 56.5
- **Mean rel error** (where expected>0): 0.256
- **Median rel error**: 0.113
- **Max rel error**: 1.951

## Negative Control False Positives

- **N conditions**: 278
- **Median FP**: 1.0
- **Mean FP**: 3.4
- **Max FP**: 26.0
- **Conditions with FP > 5**: 64

### Worst FP cases

| Scenario | Condition | FP count | gDNA_exp | nRNA_exp |
| --- | --- | ---: | ---: | ---: |
| spliced | gdna=100_ss=0.65_nrna=30 | 26 | 437 | 25 |
| spliced | gdna=50_ss=0.65_nrna=0 | 22 | 426 | 0 |
| spliced | gdna=100_ss=0.65_nrna=50 | 19 | 423 | 41 |
| single_exon | gdna=50_ss=0.9_nrna=50 | 17 | 332 | 63 |
| spliced | gdna=50_ss=0.8_nrna=0 | 17 | 426 | 0 |
| spliced | gdna=50_ss=0.9_nrna=0 | 17 | 426 | 0 |
| single_exon | gdna=50_ss=0.65_nrna=30 | 16 | 350 | 40 |
| spliced | gdna=100_ss=0.65_nrna=0 | 16 | 460 | 0 |
| spliced | gdna=100_ss=0.8_nrna=50 | 16 | 423 | 41 |
| spliced | gdna=100_ss=0.9_nrna=50 | 16 | 423 | 41 |
| spliced | gdna=100_ss=0.95_nrna=50 | 16 | 423 | 41 |
| spliced | gdna=50_ss=0.95_nrna=0 | 15 | 426 | 0 |
| spliced | gdna=100_ss=0.8_nrna=30 | 15 | 437 | 25 |
| spliced | gdna=50_ss=1.0_nrna=0 | 14 | 426 | 0 |
| spliced | gdna=50_ss=1.0_nrna=50 | 14 | 366 | 71 |
| spliced | gdna=100_ss=0.8_nrna=0 | 14 | 460 | 0 |
| spliced | gdna=100_ss=0.9_nrna=0 | 14 | 460 | 0 |
| spliced | gdna=100_ss=0.9_nrna=30 | 14 | 437 | 25 |
| spliced | gdna=100_ss=0.95_nrna=30 | 14 | 437 | 25 |
| nonoverlap | gdna=50_ss=0.65_nrna=0 | 14 | 445 | 0 |

## Per-Scenario Accuracy

### antisense

- N observations: 48
- Mean abs error: 5.7, Median: 3.5, Max: 27.6
- Mean rel error: 0.066, Median: 0.053, Max: 0.237

### dist_paralogs

- N observations: 24
- Mean abs error: 15.3, Median: 12.9, Max: 33.5
- Mean rel error: 0.337, Median: 0.207, Max: 1.232

### iso_10_1

- N observations: 24
- Mean abs error: 8.3, Median: 4.4, Max: 47.3
- Mean rel error: 0.082, Median: 0.055, Max: 0.413

### iso_equal

- N observations: 24
- Mean abs error: 15.4, Median: 17.7, Max: 33.8
- Mean rel error: 0.072, Median: 0.060, Max: 0.200

### nonoverlap

- N observations: 100
- Mean abs error: 12.1, Median: 3.9, Max: 164.0
- Mean rel error: 0.214, Median: 0.114, Max: 1.118

### paralogs_spliced

- N observations: 12
- Mean abs error: 9.0, Median: 8.3, Max: 15.1
- Mean rel error: 0.174, Median: 0.098, Max: 0.607

### paralogs_unspliced

- N observations: 24
- Mean abs error: 33.2, Median: 35.2, Max: 63.0
- Mean rel error: 0.764, Median: 0.715, Max: 1.951

### single_exon

- N observations: 75
- Mean abs error: 58.8, Median: 56.0, Max: 338.0
- Mean rel error: 0.547, Median: 0.490, Max: 1.341

### spliced

- N observations: 75
- Mean abs error: 10.8, Median: 6.1, Max: 41.5
- Mean rel error: 0.082, Median: 0.070, Max: 0.268

## gDNA Accuracy (where gDNA expected > 0)

| Scenario | Condition | gDNA exp | gDNA pipe | gDNA diff | gDNA rel err |
| --- | --- | ---: | ---: | ---: | ---: |
| single_exon | gdna=20_ss=0.65_nrna=30 | 242 | 492 | 250 | 1.03 |
| single_exon | gdna=50_ss=0.65_nrna=50 | 332 | 486 | 154 | 0.46 |
| iso_10_1 | gdna=50_ss=1.0_nrna=30 | 674 | 547 | 127 | 0.19 |
| single_exon | gdna=100_ss=0.65_nrna=50 | 399 | 282 | 117 | 0.29 |
| single_exon | gdna=50_ss=0.65_nrna=30 | 350 | 234 | 116 | 0.33 |
| single_exon | gdna=50_ss=0.8_nrna=0 | 380 | 265 | 115 | 0.30 |
| single_exon | gdna=50_ss=0.9_nrna=0 | 380 | 265 | 115 | 0.30 |
| single_exon | gdna=50_ss=1.0_nrna=50 | 332 | 222 | 110 | 0.33 |
| single_exon | gdna=50_ss=0.9_nrna=50 | 332 | 223 | 109 | 0.33 |
| single_exon | gdna=50_ss=0.65_nrna=0 | 380 | 277 | 103 | 0.27 |
| spliced | gdna=100_ss=0.65_nrna=30 | 437 | 340 | 97 | 0.22 |
| paralogs_spliced | gdna=50_ss=0.9_nrna=0 | 439 | 347 | 92 | 0.21 |
| paralogs_spliced | gdna=50_ss=1.0_nrna=0 | 439 | 350 | 89 | 0.20 |
| spliced | gdna=50_ss=0.65_nrna=50 | 366 | 278 | 88 | 0.24 |
| spliced | gdna=100_ss=0.8_nrna=30 | 437 | 350 | 87 | 0.20 |
| iso_10_1 | gdna=50_ss=0.9_nrna=30 | 674 | 587 | 87 | 0.13 |
| antisense | gdna=20_ss=0.8_nrna=30 | 266 | 351 | 85 | 0.32 |
| single_exon | gdna=100_ss=0.65_nrna=30 | 412 | 494 | 82 | 0.20 |
| paralogs_spliced | gdna=20_ss=1.0_nrna=0 | 372 | 290 | 82 | 0.22 |
| antisense | gdna=20_ss=0.9_nrna=30 | 266 | 347 | 81 | 0.30 |
| spliced | gdna=50_ss=0.8_nrna=0 | 426 | 348 | 78 | 0.18 |
| iso_equal | gdna=50_ss=0.9_nrna=30 | 558 | 480 | 78 | 0.14 |
| iso_equal | gdna=50_ss=1.0_nrna=30 | 558 | 480 | 78 | 0.14 |
| antisense | gdna=20_ss=0.95_nrna=30 | 266 | 344 | 78 | 0.29 |
| single_exon | gdna=50_ss=0.8_nrna=50 | 332 | 259 | 73 | 0.22 |
| iso_10_1 | gdna=20_ss=1.0_nrna=30 | 452 | 380 | 72 | 0.16 |
| single_exon | gdna=50_ss=1.0_nrna=30 | 350 | 280 | 70 | 0.20 |
| paralogs_unspliced | gdna=50_ss=0.9_nrna=30 | 441 | 371 | 70 | 0.16 |
| single_exon | gdna=100_ss=1.0_nrna=0 | 432 | 363 | 69 | 0.16 |
| single_exon | gdna=50_ss=0.95_nrna=0 | 380 | 312 | 68 | 0.18 |

## nRNA Accuracy (where nRNA expected > 0)

| Scenario | Condition | nRNA exp | nRNA pipe | nRNA diff | nRNA rel err |
| --- | --- | ---: | ---: | ---: | ---: |
| antisense | gdna=0_ss=1.0_nrna=30 | 201 | 2 | 199 | 0.99 |
| antisense | gdna=0_ss=0.95_nrna=30 | 201 | 4 | 197 | 0.98 |
| antisense | gdna=0_ss=0.9_nrna=30 | 201 | 6 | 195 | 0.97 |
| antisense | gdna=0_ss=0.8_nrna=30 | 201 | 10 | 191 | 0.95 |
| iso_10_1 | gdna=50_ss=1.0_nrna=30 | 150 | 280 | 130 | 0.86 |
| paralogs_unspliced | gdna=0_ss=0.9_nrna=30 | 115 | 0 | 115 | 1.00 |
| paralogs_unspliced | gdna=0_ss=1.0_nrna=30 | 115 | 0 | 115 | 1.00 |
| iso_10_1 | gdna=50_ss=0.9_nrna=30 | 150 | 256 | 106 | 0.71 |
| iso_equal | gdna=50_ss=0.9_nrna=30 | 225 | 321 | 96 | 0.43 |
| antisense | gdna=20_ss=1.0_nrna=30 | 94 | 1 | 93 | 0.99 |
| antisense | gdna=20_ss=0.95_nrna=30 | 94 | 2 | 92 | 0.98 |
| antisense | gdna=20_ss=0.9_nrna=30 | 94 | 3 | 91 | 0.97 |
| single_exon | gdna=0_ss=0.65_nrna=50 | 188 | 98 | 90 | 0.48 |
| antisense | gdna=20_ss=0.8_nrna=30 | 94 | 6 | 88 | 0.93 |
| single_exon | gdna=0_ss=0.8_nrna=50 | 188 | 102 | 86 | 0.46 |
| single_exon | gdna=0_ss=0.95_nrna=50 | 188 | 105 | 83 | 0.44 |
| spliced | gdna=50_ss=0.65_nrna=50 | 71 | 154 | 83 | 1.17 |
| spliced | gdna=100_ss=0.8_nrna=30 | 25 | 107 | 82 | 3.29 |
| single_exon | gdna=0_ss=0.9_nrna=50 | 188 | 106 | 82 | 0.44 |
| single_exon | gdna=5_ss=0.95_nrna=50 | 157 | 75 | 82 | 0.52 |
| single_exon | gdna=0_ss=1.0_nrna=50 | 188 | 106 | 82 | 0.43 |
| single_exon | gdna=5_ss=0.9_nrna=50 | 157 | 81 | 76 | 0.49 |
| iso_equal | gdna=50_ss=1.0_nrna=30 | 225 | 301 | 76 | 0.34 |
| iso_equal | gdna=20_ss=1.0_nrna=30 | 338 | 413 | 75 | 0.22 |
| iso_10_1 | gdna=20_ss=1.0_nrna=30 | 251 | 323 | 72 | 0.29 |
| spliced | gdna=100_ss=0.65_nrna=30 | 25 | 97 | 72 | 2.87 |
| iso_10_1 | gdna=20_ss=0.9_nrna=30 | 251 | 323 | 72 | 0.29 |
| spliced | gdna=20_ss=0.9_nrna=50 | 126 | 197 | 71 | 0.56 |
| single_exon | gdna=20_ss=0.65_nrna=30 | 69 | 0 | 69 | 1.00 |
| single_exon | gdna=5_ss=0.8_nrna=50 | 157 | 91 | 66 | 0.42 |

## Worst Per-Transcript Errors (top 30)

| Scenario | Condition | t_id | Expected | Observed | AbsDiff | RelErr |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| single_exon | gdna=0_ss=0.65_nrna=0 | t1 | 338 | 0 | 338 | 1.00 |
| single_exon | gdna=0_ss=0.65_nrna=50 | t1 | 204 | 0 | 204 | 1.00 |
| nonoverlap | gdna=0_ss=0.65_nrna=0 | t2 | 164 | 0 | 164 | 1.00 |
| nonoverlap | gdna=0_ss=0.8_nrna=0 | t2 | 164 | 0 | 164 | 1.00 |
| single_exon | gdna=20_ss=0.65_nrna=30 | t1 | 132 | 0 | 132 | 1.00 |
| single_exon | gdna=0_ss=1.0_nrna=50 | t1 | 209 | 317 | 108 | 0.52 |
| single_exon | gdna=0_ss=0.8_nrna=50 | t1 | 204 | 304 | 100 | 0.49 |
| single_exon | gdna=0_ss=0.9_nrna=50 | t1 | 204 | 304 | 100 | 0.49 |
| single_exon | gdna=0_ss=0.95_nrna=50 | t1 | 204 | 304 | 100 | 0.49 |
| single_exon | gdna=5_ss=1.0_nrna=50 | t1 | 182 | 278 | 96 | 0.53 |
| single_exon | gdna=5_ss=0.65_nrna=50 | t1 | 172 | 257 | 85 | 0.49 |
| single_exon | gdna=5_ss=0.8_nrna=50 | t1 | 172 | 257 | 85 | 0.49 |
| single_exon | gdna=5_ss=0.9_nrna=50 | t1 | 172 | 257 | 85 | 0.49 |
| single_exon | gdna=20_ss=1.0_nrna=50 | t1 | 111 | 195 | 84 | 0.76 |
| single_exon | gdna=0_ss=0.65_nrna=30 | t1 | 237 | 314 | 77 | 0.32 |
| single_exon | gdna=0_ss=0.8_nrna=30 | t1 | 237 | 314 | 77 | 0.32 |
| single_exon | gdna=0_ss=0.9_nrna=30 | t1 | 237 | 314 | 77 | 0.32 |
| single_exon | gdna=0_ss=0.95_nrna=30 | t1 | 237 | 314 | 77 | 0.32 |
| single_exon | gdna=0_ss=1.0_nrna=30 | t1 | 244 | 321 | 77 | 0.32 |
| single_exon | gdna=20_ss=0.65_nrna=50 | t1 | 113 | 190 | 77 | 0.68 |
| single_exon | gdna=20_ss=0.8_nrna=50 | t1 | 113 | 190 | 77 | 0.68 |
| single_exon | gdna=50_ss=1.0_nrna=50 | t1 | 72 | 147 | 75 | 1.04 |
| single_exon | gdna=20_ss=1.0_nrna=30 | t1 | 129 | 202 | 73 | 0.57 |
| single_exon | gdna=50_ss=0.9_nrna=50 | t1 | 67 | 140 | 73 | 1.09 |
| single_exon | gdna=50_ss=0.95_nrna=50 | t1 | 67 | 140 | 73 | 1.09 |
| single_exon | gdna=50_ss=0.8_nrna=50 | t1 | 67 | 139 | 72 | 1.07 |
| single_exon | gdna=5_ss=1.0_nrna=30 | t1 | 194 | 265 | 71 | 0.37 |
| single_exon | gdna=5_ss=0.65_nrna=30 | t1 | 198 | 265 | 67 | 0.34 |
| single_exon | gdna=5_ss=0.8_nrna=30 | t1 | 198 | 265 | 67 | 0.34 |
| single_exon | gdna=5_ss=0.9_nrna=30 | t1 | 198 | 265 | 67 | 0.34 |

## RNA Count Inflation by gDNA Level

When gDNA is present, do RNA counts inflate (observed > expected)?

### antisense

| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |
| --- | ---: | ---: | ---: |
| gdna=0 | +6.2 | 9.0 | 0.056 |
| gdna=20 | +2.9 | 4.8 | 0.066 |
| gdna=50 | -0.5 | 3.2 | 0.078 |

### dist_paralogs

| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |
| --- | ---: | ---: | ---: |
| gdna=0 | -0.5 | 17.9 | 0.113 |
| gdna=20 | +6.0 | 14.0 | 0.288 |
| gdna=50 | +7.2 | 13.9 | 0.612 |

### iso_10_1

| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |
| --- | ---: | ---: | ---: |
| gdna=0 | -11.5 | 12.8 | 0.074 |
| gdna=20 | -2.0 | 4.6 | 0.051 |
| gdna=50 | -5.5 | 9.7 | 0.140 |

### iso_equal

| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |
| --- | ---: | ---: | ---: |
| gdna=0 | -10.1 | 22.9 | 0.074 |
| gdna=20 | -7.5 | 15.9 | 0.082 |
| gdna=50 | -1.7 | 7.4 | 0.060 |

### nonoverlap

| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |
| --- | ---: | ---: | ---: |
| gdna=0 | -15.0 | 33.3 | 0.221 |
| gdna=100 | -2.3 | 2.3 | 0.116 |
| gdna=20 | +4.5 | 6.0 | 0.140 |
| gdna=5 | +3.5 | 12.6 | 0.128 |
| gdna=50 | +5.0 | 5.7 | 0.325 |

### paralogs_spliced

| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |
| --- | ---: | ---: | ---: |
| gdna=0 | -0.2 | 8.0 | 0.032 |
| gdna=20 | +6.4 | 9.4 | 0.149 |
| gdna=50 | +9.5 | 9.5 | 0.341 |

### paralogs_unspliced

| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |
| --- | ---: | ---: | ---: |
| gdna=0 | +28.8 | 32.8 | 0.166 |
| gdna=20 | +33.3 | 34.8 | 0.712 |
| gdna=50 | +32.1 | 32.1 | 1.416 |

### single_exon

| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |
| --- | ---: | ---: | ---: |
| gdna=0 | +16.7 | 89.0 | 0.373 |
| gdna=100 | +3.1 | 39.5 | 0.867 |
| gdna=20 | +41.2 | 58.8 | 0.470 |
| gdna=5 | +53.9 | 53.9 | 0.288 |
| gdna=50 | +43.9 | 52.9 | 0.736 |

### spliced

| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |
| --- | ---: | ---: | ---: |
| gdna=0 | -21.0 | 22.5 | 0.086 |
| gdna=100 | -3.9 | 4.8 | 0.127 |
| gdna=20 | -0.2 | 9.2 | 0.074 |
| gdna=5 | -11.4 | 13.2 | 0.064 |
| gdna=50 | -3.3 | 4.2 | 0.061 |

## Strand Specificity Impact on Accuracy

### antisense

| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |
| --- | ---: | ---: | ---: | ---: |
| ss=0.8 | -0.7 | 3.9 | 0.063 | 1.0 |
| ss=0.9 | +2.0 | 4.4 | 0.047 | 0.8 |
| ss=0.95 | +2.6 | 6.7 | 0.064 | 0.8 |
| ss=1.0 | +7.6 | 7.7 | 0.092 | 2.2 |

### dist_paralogs

| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |
| --- | ---: | ---: | ---: | ---: |
| ss=0.9 | +4.6 | 16.4 | 0.345 | 1.3 |
| ss=1.0 | +3.9 | 14.2 | 0.330 | 1.0 |

### iso_10_1

| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |
| --- | ---: | ---: | ---: | ---: |
| ss=0.9 | -9.1 | 9.5 | 0.085 | 1.3 |
| ss=1.0 | -3.7 | 8.5 | 0.082 | 1.8 |

### iso_equal

| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |
| --- | ---: | ---: | ---: | ---: |
| ss=0.9 | -8.6 | 14.5 | 0.062 | 0.7 |
| ss=1.0 | -4.3 | 16.4 | 0.081 | 1.2 |

### nonoverlap

| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |
| --- | ---: | ---: | ---: | ---: |
| ss=0.65 | -3.9 | 18.5 | 0.262 | 4.4 |
| ss=0.8 | -5.5 | 18.6 | 0.239 | 4.1 |
| ss=0.9 | +1.6 | 10.8 | 0.173 | 3.8 |
| ss=0.95 | +1.9 | 10.2 | 0.170 | 3.4 |
| ss=1.0 | +2.5 | 7.4 | 0.125 | 2.1 |

### paralogs_spliced

| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |
| --- | ---: | ---: | ---: | ---: |
| ss=0.9 | +5.8 | 8.9 | 0.173 | 1.0 |
| ss=1.0 | +4.7 | 9.0 | 0.175 | 1.7 |

### paralogs_unspliced

| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |
| --- | ---: | ---: | ---: | ---: |
| ss=0.9 | +32.0 | 33.5 | 0.778 | 1.3 |
| ss=1.0 | +30.8 | 33.0 | 0.751 | 1.2 |

### single_exon

| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |
| --- | ---: | ---: | ---: | ---: |
| ss=0.65 | -23.0 | 89.0 | 0.711 | 3.3 |
| ss=0.8 | +40.6 | 52.9 | 0.526 | 2.9 |
| ss=0.9 | +45.1 | 50.6 | 0.502 | 3.1 |
| ss=0.95 | +43.3 | 48.7 | 0.489 | 2.3 |
| ss=1.0 | +52.9 | 52.9 | 0.508 | 1.9 |

### spliced

| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |
| --- | ---: | ---: | ---: | ---: |
| ss=0.65 | +0.4 | 7.3 | 0.062 | 8.1 |
| ss=0.8 | -9.0 | 12.1 | 0.097 | 6.3 |
| ss=0.9 | -11.7 | 12.5 | 0.091 | 6.1 |
| ss=0.95 | -10.5 | 11.3 | 0.074 | 5.8 |
| ss=1.0 | -9.1 | 10.7 | 0.087 | 5.1 |
