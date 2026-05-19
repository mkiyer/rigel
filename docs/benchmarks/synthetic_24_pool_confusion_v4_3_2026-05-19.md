# Synthetic 24 Pool Confusion Report - 2026-05-19

Base directory: `/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24`

Rigel output namespace: `<condition>/rigel_out`

Annotated BAM namespace: `<condition>/annotated.bam`

Analysis artifacts: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19`

## Executive Summary

This rerun completed all 24 synthetic scenarios with the current effective-length
normalization code. The pool-level confusion matrices use oracle read names as truth
and Rigel's annotated-BAM `ZF` winner flags as the predicted pool. Counts are therefore
hard winner assignments; the fractional EM pool masses are also written to the metrics
table from `loci.feather`.

In strand-specific cases (`ss=0.99`), Rigel recovers true gDNA well when gDNA is present:
mean gDNA winner recall across nonzero-gDNA stranded scenarios is
99.06%, with mean gDNA-to-RNA leakage of
0.93%. Mature RNA is also rarely mistaken for gDNA in clean
stranded, no-nRNA settings: mean mRNA-to-gDNA is
0.63%.

The remaining hard problem is not ordinary mature mRNA versus gDNA. It is genic,
unspliced nRNA versus genic gDNA. Across stranded scenarios with true nRNA, mean nRNA
winner recall is 68.46% and mean nRNA-to-gDNA leakage is
14.68%. This is a major improvement over unstranded nRNA
scenarios, where mean nRNA recall is 29.05% and mean
nRNA-to-gDNA leakage is 53.91%, but it is still not the
dramatic separation we want.

The zero-gDNA, high-nRNA contrast is the cleanest diagnostic. In the unstranded run,
Rigel sends 84.06% of true nRNA fragments
to gDNA. In the stranded run, this falls to
5.18%, and nRNA recall rises from
0.00% to
78.82%. Strand information is being used, but
genic unspliced evidence remains partially gDNA-like.

## Aggregate 3x3 Pool Confusion

These are hard winner assignments from `ZF` tags. Rates are row-normalized by true pool.

| Group   | True | Pred mRNA  | Pred nRNA | Pred gDNA  | mRNA rate | nRNA rate | gDNA rate |
| ------- | ---- | ---------- | --------- | ---------- | --------- | --------- | --------- |
| all 24  | mrna | 23,089,365 | 494,331   | 416,304    | 96.21%    | 2.06%     | 1.73%     |
| all 24  | nrna | 443,650    | 1,336,104 | 900,918    | 16.55%    | 49.83%    | 33.60%    |
| all 24  | gdna | 100,516    | 56,211    | 19,341,133 | 0.52%     | 0.29%     | 99.19%    |
| ss=0.50 | mrna | 11,497,224 | 180,637   | 322,139    | 95.81%    | 1.51%     | 2.68%     |
| ss=0.50 | nrna | 222,250    | 380,413   | 737,673    | 16.58%    | 28.38%    | 55.03%    |
| ss=0.50 | gdna | 56,993     | 22,820    | 9,669,117  | 0.58%     | 0.23%     | 99.17%    |
| ss=0.99 | mrna | 11,592,141 | 313,694   | 94,165     | 96.60%    | 2.61%     | 0.78%     |
| ss=0.99 | nrna | 221,400    | 955,691   | 163,245    | 16.52%    | 71.29%    | 12.18%    |
| ss=0.99 | gdna | 43,523     | 33,391    | 9,672,016  | 0.45%     | 0.34%     | 99.20%    |

## Per-Scenario Pool Confusion

The table below reports row-normalized winner assignments for each true pool. `nRNA to
gDNA` means true nascent fragments assigned to the gDNA pool. `gDNA to RNA` combines
gDNA assigned to mature or nascent RNA.

| Condition                   | Truth gDNA | Pred gDNA | Truth nRNA | nRNA to nRNA | nRNA to gDNA | gDNA to gDNA | gDNA to RNA | mRNA to gDNA | Pool acc |
| --------------------------- | ---------- | --------- | ---------- | ------------ | ------------ | ------------ | ----------- | ------------ | -------- |
| gdna_zero_ss_0.50_nrna_high | 0.00%      | 23.88%    | 21.62%     | 0.00%        | 84.06%       | n/a          | n/a         | 7.28%        | 71.95%   |
| gdna_zero_ss_0.99_nrna_high | 0.00%      | 1.42%     | 21.62%     | 78.82%       | 5.18%        | n/a          | n/a         | 0.39%        | 89.67%   |
| gdna_low_ss_0.50_nrna_high  | 16.38%     | 28.80%    | 18.08%     | 30.86%       | 52.53%       | 99.02%       | 0.97%       | 4.70%        | 82.75%   |
| gdna_low_ss_0.99_nrna_high  | 16.38%     | 17.85%    | 18.08%     | 75.71%       | 7.65%        | 98.23%       | 1.76%       | 0.57%        | 91.82%   |
| gdna_med_ss_0.50_nrna_high  | 43.94%     | 50.42%    | 12.12%     | 41.69%       | 41.81%       | 98.98%       | 1.01%       | 4.24%        | 89.36%   |
| gdna_med_ss_0.99_nrna_high  | 43.94%     | 45.40%    | 12.12%     | 70.89%       | 12.67%       | 98.73%       | 1.26%       | 1.11%        | 93.31%   |
| gdna_high_ss_0.50_nrna_high | 61.05%     | 65.71%    | 8.42%      | 39.48%       | 44.16%       | 99.05%       | 0.93%       | 4.97%        | 91.96%   |
| gdna_high_ss_0.99_nrna_high | 61.05%     | 62.48%    | 8.42%      | 65.91%       | 17.75%       | 98.96%       | 1.02%       | 1.85%        | 94.58%   |
| gdna_zero_ss_0.50_nrna_low  | 0.00%      | 6.09%     | 5.59%      | 0.00%        | 82.04%       | n/a          | n/a         | 1.59%        | 91.89%   |
| gdna_zero_ss_0.99_nrna_low  | 0.00%      | 0.37%     | 5.59%      | 76.89%       | 5.38%        | n/a          | n/a         | 0.08%        | 97.06%   |
| gdna_low_ss_0.50_nrna_low   | 19.09%     | 20.99%    | 4.53%      | 50.11%       | 32.30%       | 98.81%       | 1.18%       | 0.86%        | 95.83%   |
| gdna_low_ss_0.99_nrna_low   | 19.09%     | 19.63%    | 4.53%      | 71.06%       | 11.80%       | 98.74%       | 1.25%       | 0.31%        | 97.00%   |
| gdna_med_ss_0.50_nrna_low   | 48.56%     | 50.17%    | 2.88%      | 40.17%       | 41.78%       | 99.07%       | 0.92%       | 1.76%        | 96.35%   |
| gdna_med_ss_0.99_nrna_low   | 48.56%     | 49.28%    | 2.88%      | 59.08%       | 23.35%       | 99.18%       | 0.81%       | 0.92%        | 97.27%   |
| gdna_high_ss_0.50_nrna_low  | 65.38%     | 66.89%    | 1.94%      | 30.11%       | 52.62%       | 99.25%       | 0.74%       | 3.01%        | 96.78%   |
| gdna_high_ss_0.99_nrna_low  | 65.38%     | 66.15%    | 1.94%      | 49.35%       | 33.63%       | 99.36%       | 0.63%       | 1.64%        | 97.63%   |
| gdna_zero_ss_0.50_nrna_zero | 0.00%      | 0.00%     | 0.00%      | n/a          | n/a          | n/a          | n/a         | 0.00%        | 98.89%   |
| gdna_zero_ss_0.99_nrna_zero | 0.00%      | 0.00%     | 0.00%      | n/a          | n/a          | n/a          | n/a         | 0.00%        | 98.89%   |
| gdna_low_ss_0.50_nrna_zero  | 20.00%     | 20.05%    | 0.00%      | n/a          | n/a          | 99.11%       | 0.89%       | 0.28%        | 98.71%   |
| gdna_low_ss_0.99_nrna_zero  | 20.00%     | 20.04%    | 0.00%      | n/a          | n/a          | 99.29%       | 0.71%       | 0.23%        | 98.79%   |
| gdna_med_ss_0.50_nrna_zero  | 50.00%     | 50.22%    | 0.00%      | n/a          | n/a          | 99.28%       | 0.71%       | 1.16%        | 98.53%   |
| gdna_med_ss_0.99_nrna_zero  | 50.00%     | 50.14%    | 0.00%      | n/a          | n/a          | 99.47%       | 0.52%       | 0.80%        | 98.79%   |
| gdna_high_ss_0.50_nrna_zero | 66.67%     | 67.03%    | 0.00%      | n/a          | n/a          | 99.38%       | 0.61%       | 2.34%        | 98.46%   |
| gdna_high_ss_0.99_nrna_zero | 66.67%     | 66.87%    | 0.00%      | n/a          | n/a          | 99.55%       | 0.44%       | 1.50%        | 98.85%   |

## Strand-Specific Delta

Each row compares `ss=0.99` minus `ss=0.50` for the same gDNA/nRNA setting. Positive
`d nRNA recall` and negative `d nRNA to gDNA` are improvements.

| gDNA | nRNA | d pool acc | d nRNA recall | d nRNA to gDNA | d gDNA recall | d gDNA to RNA | d assigned gDNA |
| ---- | ---- | ---------- | ------------- | -------------- | ------------- | ------------- | --------------- |
| high | high | 2.61%      | 26.43%        | -26.41%        | -0.09%        | 0.09%         | -3.23%          |
| low  | high | 9.07%      | 44.85%        | -44.88%        | -0.79%        | 0.79%         | -10.95%         |
| med  | high | 3.95%      | 29.20%        | -29.14%        | -0.25%        | 0.25%         | -5.02%          |
| zero | high | 17.72%     | 78.82%        | -78.88%        | n/a           | n/a           | -22.46%         |
| high | low  | 0.85%      | 19.25%        | -19.00%        | 0.12%         | -0.12%        | -0.74%          |
| low  | low  | 1.17%      | 20.95%        | -20.50%        | -0.06%        | 0.06%         | -1.36%          |
| med  | low  | 0.93%      | 18.91%        | -18.43%        | 0.11%         | -0.11%        | -0.88%          |
| zero | low  | 5.18%      | 76.89%        | -76.66%        | n/a           | n/a           | -5.72%          |
| high | zero | 0.39%      | n/a           | n/a            | 0.18%         | -0.18%        | -0.16%          |
| low  | zero | 0.07%      | n/a           | n/a            | 0.18%         | -0.18%        | -0.00%          |
| med  | zero | 0.27%      | n/a           | n/a            | 0.19%         | -0.19%        | -0.09%          |
| zero | zero | -0.00%     | n/a           | n/a            | n/a           | n/a           | 0.00%           |

## What The Strand Signal Fixes

Strand-specific data sharply reduces intronic nRNA masquerading as gDNA. nRNA is
unspliced and spans introns, while mature mRNA is confined to exons and exon-exon splice
junctions. In stranded data, opposite-strand genic-unspliced evidence can be discounted
as gDNA-like, so the calibration correction suppresses much of the intron-only false
gDNA signal. This is visible in the zero-gDNA high-nRNA pair above and in the large
negative `d nRNA to gDNA` values for low/high nRNA settings.

## What Remains Confused

The strand signal does not fully resolve genic unspliced fragments that lack splice
junction evidence. Some are exon-intron boundary fragments that can feed the
EXON-INTRON calibration channel; others are exon-confined mature-RNA fragments, which
are compatible with both a mature exon and the corresponding nascent transcript span.
In zero-gDNA high-nRNA stranded data, the top error categories are:

| True | Pred | ZC                | ZS        | Count  |
| ---- | ---- | ----------------- | --------- | ------ |
| mrna | nrna | ambig_same_strand | unspliced | 59,532 |
| nrna | mrna | ambig_same_strand | unspliced | 33,895 |
| nrna | gdna | ambig_same_strand | unspliced | 10,671 |
| nrna | mrna | ambig_opp_strand  | unspliced | 10,196 |
| mrna | nrna | ambig_opp_strand  | unspliced | 9,952  |
| nrna | gdna | ambig_opp_strand  | unspliced | 3,608  |
| mrna | gdna | ambig_same_strand | unspliced | 2,974  |
| mrna | gdna | ambig_opp_strand  | unspliced | 897    |

When true gDNA and high nRNA coexist in stranded data, the largest residual errors are:

| True | Pred | ZC                | ZS        | Count  |
| ---- | ---- | ----------------- | --------- | ------ |
| nrna | gdna | ambig_same_strand | unspliced | 37,374 |
| mrna | nrna | ambig_same_strand | unspliced | 37,025 |
| nrna | mrna | ambig_same_strand | unspliced | 34,859 |
| mrna | gdna | ambig_same_strand | unspliced | 14,367 |
| nrna | gdna | ambig_opp_strand  | unspliced | 11,596 |
| nrna | mrna | ambig_opp_strand  | unspliced | 10,156 |
| gdna | nrna | ambig_same_strand | unspliced | 9,220  |
| mrna | nrna | ambig_opp_strand  | unspliced | 7,399  |

## Calibration Diagnostics In High-nRNA Scenarios

When true gDNA is absent or low, the gDNA FL model and genic gDNA densities can be
pulled toward nRNA-like evidence. `rho_ig` is the intergenic anchor; `rho_in` and
`rho_ex_in` are genic unspliced channels.

| Condition                   | SS   | gDNA FL | rho ig | rho intron | rho exon-intron | Pred gDNA | nRNA to gDNA |
| --------------------------- | ---- | ------- | ------ | ---------- | --------------- | --------- | ------------ |
| gdna_zero_ss_0.50_nrna_high | 0.50 | 250.1   | 0      | 0.227      | 0.225           | 23.88%    | 84.06%       |
| gdna_zero_ss_0.99_nrna_high | 0.99 | 250.1   | 0      | 7.79e-05   | 0.0775          | 1.42%     | 5.18%        |
| gdna_low_ss_0.50_nrna_high  | 0.50 | 303.6   | 0.0249 | 0.255      | 0.212           | 28.80%    | 52.53%       |
| gdna_low_ss_0.99_nrna_high  | 0.99 | 303.6   | 0.0249 | 0.0253     | 0.0901          | 17.85%    | 7.65%        |
| gdna_med_ss_0.50_nrna_high  | 0.50 | 332.8   | 0.1    | 0.332      | 0.271           | 50.42%    | 41.81%       |
| gdna_med_ss_0.99_nrna_high  | 0.99 | 332.8   | 0.1    | 0.0997     | 0.156           | 45.40%    | 12.67%       |
| gdna_high_ss_0.50_nrna_high | 0.50 | 341.2   | 0.2    | 0.433      | 0.366           | 65.71%    | 44.16%       |
| gdna_high_ss_0.99_nrna_high | 0.99 | 341.2   | 0.2    | 0.201      | 0.256           | 62.48%    | 17.75%       |

## Predicted-gDNA / No-transcript Bucket

Current annotated BAMs stamp gDNA winners with `ZT=.` and `ZL=-1`, so false gDNA
assignments cannot yet be localized back to the EM locus from the BAM alone. The table
below separates that diagnostic bucket from true biological loci. Preserving the source
locus id for gDNA winners would make future root-cause analysis much sharper.

| Condition                   | True | Pred | Count   |
| --------------------------- | ---- | ---- | ------- |
| gdna_zero_ss_0.50_nrna_high | nrna | gdna | 231,897 |
| gdna_low_ss_0.50_nrna_high  | nrna | gdna | 144,910 |
| gdna_high_ss_0.50_nrna_high | nrna | gdna | 121,822 |
| gdna_med_ss_0.50_nrna_high  | nrna | gdna | 115,342 |
| gdna_zero_ss_0.50_nrna_high | mrna | gdna | 72,812  |
| gdna_high_ss_0.50_nrna_high | mrna | gdna | 49,680  |
| gdna_high_ss_0.99_nrna_high | nrna | gdna | 48,970  |
| gdna_zero_ss_0.50_nrna_low  | nrna | gdna | 48,616  |
| gdna_low_ss_0.50_nrna_high  | mrna | gdna | 47,025  |
| gdna_med_ss_0.50_nrna_high  | mrna | gdna | 42,388  |
| gdna_med_ss_0.99_nrna_high  | nrna | gdna | 34,953  |
| gdna_high_ss_0.50_nrna_low  | nrna | gdna | 31,185  |

## Worst EM Loci By Non-gDNA Pool Error

| Condition                   | Locus | Errors | Err frac | Truth mRNA | Truth nRNA | Truth gDNA | Pred gDNA rate | n tx |
| --------------------------- | ----- | ------ | -------- | ---------- | ---------- | ---------- | -------------- | ---- |
| gdna_zero_ss_0.99_nrna_high | 23    | 24,831 | 15.24%   | 111,238    | 51,733     | 0          | 0.78%          | 15   |
| gdna_high_ss_0.99_nrna_high | 24    | 21,075 | 13.08%   | 110,518    | 48,828     | 1,814      | 5.16%          | 15   |
| gdna_med_ss_0.99_nrna_high  | 24    | 20,824 | 12.84%   | 110,915    | 50,181     | 1,086      | 2.96%          | 15   |
| gdna_low_ss_0.99_nrna_high  | 24    | 20,412 | 12.55%   | 111,177    | 51,204     | 308        | 1.37%          | 15   |
| gdna_high_ss_0.50_nrna_high | 24    | 19,789 | 13.20%   | 107,547    | 40,179     | 2,247      | 11.74%         | 15   |
| gdna_med_ss_0.50_nrna_high  | 24    | 18,910 | 12.57%   | 108,114    | 41,091     | 1,208      | 10.00%         | 15   |
| gdna_low_ss_0.50_nrna_high  | 24    | 16,633 | 11.65%   | 106,752    | 35,803     | 276        | 13.41%         | 15   |
| gdna_high_ss_0.99_nrna_high | 4     | 15,109 | 13.27%   | 83,569     | 29,006     | 1,303      | 4.79%          | 15   |
| gdna_med_ss_0.99_nrna_high  | 4     | 15,037 | 13.15%   | 83,823     | 29,817     | 744        | 2.72%          | 15   |
| gdna_low_ss_0.99_nrna_high  | 4     | 14,545 | 12.67%   | 84,068     | 30,464     | 231        | 1.19%          | 15   |
| gdna_zero_ss_0.99_nrna_high | 4     | 14,366 | 12.53%   | 84,105     | 30,591     | 0          | 0.83%          | 15   |
| gdna_high_ss_0.50_nrna_high | 4     | 13,691 | 12.85%   | 81,430     | 23,608     | 1,487      | 10.94%         | 15   |

## Root-Cause Interpretation

1. Strand-specific data works as intended for mature RNA versus gDNA. Mature fragments
   carry exonic or splice-junction evidence, and true gDNA has intergenic support plus
   continuous genomic compatibility. In stranded no-nRNA cases, mRNA-to-gDNA leakage is
   very small.
2. The main residual ambiguity is genic unspliced signal. nRNA and genic gDNA are both
    continuous over the transcript span and both cover introns. Mature mRNA also produces
    unspliced exon-confined fragments that are compatible with nRNA unless a splice
    junction or intronic evidence separates them. Strand correction removes much of the
    intron-only nRNA-as-gDNA signal, but it does not create a positive nRNA state during
    calibration; remaining genic-unspliced signal can still seed gDNA priors.
3. The gDNA fragment-length model is vulnerable when nRNA is abundant. If genic
   unspliced nRNA is admitted into gDNA training sources, the fitted gDNA FL can become
   RNA-like, which then reinforces nRNA/gDNA confusion in EM.
4. Pool-level mRNA behavior is better than transcript-level isoform behavior. Many mRNA
   mistakes are same-gene or same-locus isoform allocation issues rather than wrong-pool
   errors. That is a separate identifiability problem from the nRNA/gDNA confusion.

## Further Exploration

1. Build per-locus nRNA/gDNA truth summaries for the worst loci above and inspect the
   local transcript structures, especially loci with many nRNA entities or shared
   transcript spans.
2. Stratify nRNA errors by read geometry: intron-only, exon-intron boundary,
   exon-confined unspliced, and splice-junction-supported. The BAM `ZC` and `ZS` tags are
   a start, but a dedicated geometry label would make this much sharper.
3. Compare gDNA FL estimates trained from intergenic-only fragments versus the current
   pooled gDNA source in high-nRNA stranded cases.
4. Add synthetic fixtures for zero-gDNA plus high-nRNA at `ss=0.99`; this should be a
   regression guard for false gDNA from nascent RNA.
5. Preserve the source locus id for gDNA EM winners in annotated BAM output so false
    gDNA calls can be traced directly to their competing transcript/nRNA components.

## Improvements To Consider

1. Add an explicit nRNA-aware calibration state before projecting genic-unspliced signal
   into gDNA priors. Intergenic evidence should anchor true gDNA; genic intron/exon
   signal without intergenic support should not automatically become gDNA.
2. Train or regularize the gDNA FL model primarily from intergenic-only fragments when
   enough evidence exists, then use genic channels only after strand/nRNA correction.
3. Gate genic gDNA priors by consistency with intergenic density. If INTERGENIC density
   is near zero but INTRON or EXON-INTRON density is high, prefer nRNA or downweight the
   gDNA prior.
4. Introduce a locus-level competition feature between nRNA span coverage and gDNA: nRNA
   should be favored when unspliced coverage is transcript-span-local and lacks flanking
   intergenic support.
5. Report pool confusion diagnostics from annotated BAMs as a first-class benchmark
   artifact, because aggregate quant tables can hide where fragments moved.
6. Preserve locus ids for gDNA winners in annotated BAMs; this is a diagnostics
    improvement rather than a model change, but it would substantially shorten future
    root-cause analysis loops.

## Artifacts

- Condition metrics: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19/condition_metrics.tsv`
- Pool confusion counts: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19/pool_confusion_counts.tsv`
- Pool confusion matrices: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19/pool_confusion_matrices.tsv`
- Strand-pair deltas: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19/strand_pair_delta.tsv`
- ZC/ZS breakdown: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19/zc_zs_breakdown.tsv`
- ZF breakdown: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19/zf_breakdown.tsv`
- Locus pool breakdown: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19/locus_pool_breakdown.tsv`
- Predicted-gDNA/no-transcript bucket: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19/no_locus_pool_errors.tsv`
- Worst EM locus errors: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_v4_3_2026-05-19/worst_locus_pool_errors.tsv`
