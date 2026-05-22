# VCAP Locus 3 Cap-At-1 Exposure Diagnostics - 2026-05-21

Input BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam`

Reference run for locus outputs: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1`

This diagnostic reran the native scan plus calibration only, then evaluated
alternate mappings from regional density `rho_hat` to exposure weight `A_r`.
It did not rerun EM.

## Current Scalar Denominator Check

| Metric | Value |
| --- | ---: |
| n_fragments rescanned | 31,985,707 |
| reconstructed footprint bp | 353,318,962 |
| saved locus_span_bp | 353,317,609 |
| saved gdna_eff_len_unweighted | 353,762,179 |
| saved gdna_em_exposure_weight | 0.261212 |
| saved gdna_eff_len | 92,407,095 |
| recomputed current mean A | 0.261213 |
| recomputed current gdna_eff_len | 92,407,283 |

## Normalizer Sweep

| Normalizer | rho ref | footprint mean A | L_gDNA | reduction | footprint bp at A=1 % | weighted denominator from A=1 % | global exon opportunity at cap % |
| --- | --- | --- | --- | --- | --- | --- | --- |
| current_global_q95 | 3.609e-04 | 0.261213 | 92,407,283 | 3.83 | 11.10 | 42.51 | 100.00 |
| global_q99 | 1.346e-03 | 0.101606 | 35,944,243 | 9.84 | 2.09 | 20.58 | 19.03 |
| global_q99.5 | 2.383e-03 | 0.064509 | 22,820,980 | 15.50 | 1.29 | 20.02 | 15.22 |
| global_q99.9 | 6.698e-03 | 0.026853 | 9,499,710 | 37.24 | 0.26 | 9.86 | 4.61 |
| class_q95 | 6.281e-03 | 0.178649 | 63,199,426 | 5.60 | 4.01 | 22.45 | 5.00 |
| class_q99 | 4.101e-02 | 0.054953 | 19,440,390 | 18.20 | 0.72 | 13.18 | 1.00 |
| class_q99.9 | 2.652e-01 | 0.015386 | 5,442,962 | 64.99 | 0.08 | 5.52 | 0.10 |
| locus_q99 | 2.933e-03 | 0.054541 | 19,294,429 | 18.33 | 1.00 | 18.34 | 13.17 |
| locus_q99.9 | 1.276e-02 | 0.014838 | 5,249,150 | 67.39 | 0.10 | 6.74 | 2.52 |
| locus_max | 3.495e+00 | 0.000138 | 48,911 | 7232.77 | 0.00 | 0.16 | 0.00 |
| global_max | 6.123e+00 | 0.000118 | 41,737 | 8475.89 | 0.00 | 0.00 | 0.00 |
| soft_global_q95 | 3.609e-04 | 0.184119 | 65,134,493 | 5.43 | 0.00 | 0.00 | 0.00 |

Full numeric tables:

- `results/vcap_cap_at_1_diagnostics_2026-05-21/normalizer_summary.tsv`
- `results/vcap_cap_at_1_diagnostics_2026-05-21/footprint_bin_decomposition.tsv`

## Footprint Bin Decomposition

For each normalizer, the `L_gDNA` contribution is the current scalar production
denominator contribution: `gdna_eff_len_unweighted * weighted_bp_bin / footprint_bp`.

| Normalizer | A bin | footprint bp % | weighted denominator % | L_gDNA contribution |
| --- | --- | --- | --- | --- |
| current_global_q95 | [1e-4,1e-3) | 0.00 | 0.00 | 0 |
| current_global_q95 | [1e-3,1e-2) | 1.55 | 0.05 | 46,977 |
| current_global_q95 | [1e-2,1e-1) | 47.71 | 8.25 | 7,627,654 |
| current_global_q95 | [1e-1,1) | 39.64 | 49.19 | 45,453,921 |
| current_global_q95 | A==1 | 11.10 | 42.51 | 39,278,730 |
| global_max | [1e-4,1e-3) | 99.69 | 90.63 | 37,826 |
| global_max | [1e-3,1e-2) | 0.30 | 5.22 | 2,177 |
| global_max | [1e-2,1e-1) | 0.02 | 3.70 | 1,545 |
| global_max | [1e-1,1) | 0.00 | 0.45 | 189 |
| global_max | A==1 | 0.00 | 0.00 | 0 |
| class_q95 | [1e-4,1e-3) | 0.00 | 0.00 | 0 |
| class_q95 | [1e-3,1e-2) | 3.98 | 0.18 | 110,607 |
| class_q95 | [1e-2,1e-1) | 56.58 | 14.50 | 9,162,090 |
| class_q95 | [1e-1,1) | 35.43 | 62.88 | 39,739,185 |
| class_q95 | A==1 | 4.01 | 22.45 | 14,187,543 |
| class_q99.9 | [1e-4,1e-3) | 10.87 | 0.49 | 26,715 |
| class_q99.9 | [1e-3,1e-2) | 60.61 | 14.61 | 795,450 |
| class_q99.9 | [1e-2,1e-1) | 26.44 | 47.52 | 2,586,455 |
| class_q99.9 | [1e-1,1) | 2.00 | 31.85 | 1,733,835 |
| class_q99.9 | A==1 | 0.08 | 5.52 | 300,507 |
| locus_q99.9 | [1e-4,1e-3) | 22.10 | 0.92 | 48,240 |
| locus_q99.9 | [1e-3,1e-2) | 53.58 | 12.83 | 673,459 |
| locus_q99.9 | [1e-2,1e-1) | 22.13 | 42.88 | 2,250,841 |
| locus_q99.9 | [1e-1,1) | 2.10 | 36.63 | 1,922,841 |
| locus_q99.9 | A==1 | 0.10 | 6.74 | 353,769 |
| locus_max | [1e-4,1e-3) | 99.21 | 81.48 | 39,853 |
| locus_max | [1e-3,1e-2) | 0.76 | 10.88 | 5,319 |
| locus_max | [1e-2,1e-1) | 0.03 | 6.10 | 2,983 |
| locus_max | [1e-1,1) | 0.00 | 1.38 | 677 |
| locus_max | A==1 | 0.00 | 0.16 | 79 |

## Immediate Interpretation

Under the current global Q95 normalizer, the footprint mean exposure remains
near `0.2612`. Saturation is important, but it is not the
whole explanation: `A==1` covers only
`11.10%` of the reconstructed footprint
while contributing `42.51%` of
the scalar denominator. The remaining large contribution comes mostly from the
`[0.1,1)` bin, so simply replacing the hard cap with a soft cap at the same Q95
reference scale is insufficient.

The strongest useful diagnostic is the reference scale. Raising the global
normalizer to Q99.9 shrinks `L_gDNA` to 9.5M. Class-aware Q99.9 and within-locus
Q99.9 land at similar scales, 5.44M and 5.25M respectively, close to the FLG2
mRNA break-even scale observed in the hotspot audit. The raw max normalizers give
40-50k denominators, which are useful lower bounds but almost certainly
over-normalize because they push nearly the whole footprint to the exposure
floor. A production candidate should be chosen by an explicit identifying
assumption about the exposure scale, then validated on normal loci and full
confusion matrices before changing EM-facing denominators.
