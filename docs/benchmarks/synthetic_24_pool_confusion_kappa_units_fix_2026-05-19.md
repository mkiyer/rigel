# Synthetic 24 Pool Confusion Report - 2026-05-19

Base directory: `/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24`

Rigel output namespace: `<condition>/kappa_units_fix_out`

Annotated BAM namespace: `<condition>/annotated_kappa_units_fix.bam`

Analysis artifacts: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19`

## Executive Summary

This rerun completed all 24 synthetic scenarios with the current effective-length
normalization code. The pool-level confusion matrices use oracle read names as truth
and Rigel's annotated-BAM `ZF` winner flags as the predicted pool. Counts are therefore
hard winner assignments; the fractional EM pool masses are also written to the metrics
table from `loci.feather`.

In strand-specific cases (`ss=0.99`), Rigel recovers true gDNA well when gDNA is present:
mean gDNA winner recall across nonzero-gDNA stranded scenarios is
99.05%, with mean gDNA-to-RNA leakage of
0.93%. Mature RNA is also rarely mistaken for gDNA in clean
stranded, no-nRNA settings: mean mRNA-to-gDNA is
0.63%.

The remaining hard problem is not ordinary mature mRNA versus gDNA. It is genic,
unspliced nRNA versus genic gDNA. Across stranded scenarios with true nRNA, mean nRNA
winner recall is 68.52% and mean nRNA-to-gDNA leakage is
14.61%. This is a major improvement over unstranded nRNA
scenarios, where mean nRNA recall is 29.02% and mean
nRNA-to-gDNA leakage is 53.96%, but it is still not the
dramatic separation we want.

The zero-gDNA, high-nRNA contrast is the cleanest diagnostic. In the unstranded run,
Rigel sends 84.06% of true nRNA fragments
to gDNA. In the stranded run, this falls to
5.30%, and nRNA recall rises from
0.00% to
78.54%. Strand information is being used, but
genic unspliced evidence remains partially gDNA-like.

## Aggregate 3x3 Pool Confusion

These are hard winner assignments from `ZF` tags. Rates are row-normalized by true pool.

| Group   | True | Pred mRNA  | Pred nRNA | Pred gDNA  | mRNA rate | nRNA rate | gDNA rate |
| ------- | ---- | ---------- | --------- | ---------- | --------- | --------- | --------- |
| all 24  | mrna | 23,086,101 | 489,627   | 424,272    | 96.19%    | 2.04%     | 1.77%     |
| all 24  | nrna | 443,918    | 1,335,356 | 901,398    | 16.56%    | 49.81%    | 33.62%    |
| all 24  | gdna | 100,722    | 56,150    | 19,340,988 | 0.52%     | 0.29%     | 99.18%    |
| ss=0.50 | mrna | 11,490,230 | 180,239   | 329,531    | 95.75%    | 1.50%     | 2.75%     |
| ss=0.50 | nrna | 222,152    | 379,895   | 738,289    | 16.57%    | 28.34%    | 55.07%    |
| ss=0.50 | gdna | 57,082     | 22,901    | 9,668,947  | 0.59%     | 0.23%     | 99.17%    |
| ss=0.99 | mrna | 11,595,871 | 309,388   | 94,741     | 96.63%    | 2.58%     | 0.79%     |
| ss=0.99 | nrna | 221,766    | 955,461   | 163,109    | 16.54%    | 71.27%    | 12.17%    |
| ss=0.99 | gdna | 43,640     | 33,249    | 9,672,041  | 0.45%     | 0.34%     | 99.20%    |

## Per-Scenario Pool Confusion

The table below reports row-normalized winner assignments for each true pool. `nRNA to
gDNA` means true nascent fragments assigned to the gDNA pool. `gDNA to RNA` combines
gDNA assigned to mature or nascent RNA.

| Condition                   | Truth gDNA | Pred gDNA | Truth nRNA | nRNA to nRNA | nRNA to gDNA | gDNA to gDNA | gDNA to RNA | mRNA to gDNA | Pool acc |
| --------------------------- | ---------- | --------- | ---------- | ------------ | ------------ | ------------ | ----------- | ------------ | -------- |
| gdna_zero_ss_0.50_nrna_high | 0.00%      | 24.11%    | 21.62%     | 0.00%        | 84.06%       | n/a          | n/a         | 7.57%        | 71.72%   |
| gdna_zero_ss_0.99_nrna_high | 0.00%      | 1.45%     | 21.62%     | 78.54%       | 5.30%        | n/a          | n/a         | 0.39%        | 90.17%   |
| gdna_low_ss_0.50_nrna_high  | 16.38%     | 28.85%    | 18.08%     | 30.73%       | 52.60%       | 99.02%       | 0.97%       | 4.75%        | 82.70%   |
| gdna_low_ss_0.99_nrna_high  | 16.38%     | 17.83%    | 18.08%     | 75.83%       | 7.55%        | 98.23%       | 1.76%       | 0.57%        | 91.64%   |
| gdna_med_ss_0.50_nrna_high  | 43.94%     | 50.44%    | 12.12%     | 41.74%       | 41.80%       | 98.97%       | 1.01%       | 4.29%        | 89.36%   |
| gdna_med_ss_0.99_nrna_high  | 43.94%     | 45.39%    | 12.12%     | 70.92%       | 12.59%       | 98.72%       | 1.26%       | 1.11%        | 93.43%   |
| gdna_high_ss_0.50_nrna_high | 61.05%     | 65.74%    | 8.42%      | 39.39%       | 44.25%       | 99.06%       | 0.93%       | 5.04%        | 91.94%   |
| gdna_high_ss_0.99_nrna_high | 61.05%     | 62.50%    | 8.42%      | 65.81%       | 17.88%       | 98.97%       | 1.01%       | 1.86%        | 94.54%   |
| gdna_zero_ss_0.50_nrna_low  | 0.00%      | 6.31%     | 5.59%      | 0.00%        | 82.09%       | n/a          | n/a         | 1.82%        | 91.68%   |
| gdna_zero_ss_0.99_nrna_low  | 0.00%      | 0.37%     | 5.59%      | 77.10%       | 5.27%        | n/a          | n/a         | 0.08%        | 96.95%   |
| gdna_low_ss_0.50_nrna_low   | 19.09%     | 21.02%    | 4.53%      | 49.87%       | 32.54%       | 98.83%       | 1.16%       | 0.88%        | 95.80%   |
| gdna_low_ss_0.99_nrna_low   | 19.09%     | 19.61%    | 4.53%      | 71.31%       | 11.55%       | 98.72%       | 1.27%       | 0.32%        | 97.00%   |
| gdna_med_ss_0.50_nrna_low   | 48.56%     | 50.17%    | 2.88%      | 40.20%       | 41.79%       | 99.06%       | 0.93%       | 1.77%        | 96.34%   |
| gdna_med_ss_0.99_nrna_low   | 48.56%     | 49.29%    | 2.88%      | 59.19%       | 23.26%       | 99.18%       | 0.81%       | 0.94%        | 97.27%   |
| gdna_high_ss_0.50_nrna_low  | 65.38%     | 66.90%    | 1.94%      | 30.21%       | 52.52%       | 99.24%       | 0.75%       | 3.06%        | 96.76%   |
| gdna_high_ss_0.99_nrna_low  | 65.38%     | 66.16%    | 1.94%      | 49.48%       | 33.48%       | 99.37%       | 0.62%       | 1.68%        | 97.62%   |
| gdna_zero_ss_0.50_nrna_zero | 0.00%      | 0.00%     | 0.00%      | n/a          | n/a          | n/a          | n/a         | 0.00%        | 98.89%   |
| gdna_zero_ss_0.99_nrna_zero | 0.00%      | 0.00%     | 0.00%      | n/a          | n/a          | n/a          | n/a         | 0.00%        | 98.89%   |
| gdna_low_ss_0.50_nrna_zero  | 20.00%     | 20.05%    | 0.00%      | n/a          | n/a          | 99.11%       | 0.88%       | 0.28%        | 98.72%   |
| gdna_low_ss_0.99_nrna_zero  | 20.00%     | 20.04%    | 0.00%      | n/a          | n/a          | 99.29%       | 0.70%       | 0.22%        | 98.79%   |
| gdna_med_ss_0.50_nrna_zero  | 50.00%     | 50.21%    | 0.00%      | n/a          | n/a          | 99.27%       | 0.72%       | 1.14%        | 98.54%   |
| gdna_med_ss_0.99_nrna_zero  | 50.00%     | 50.13%    | 0.00%      | n/a          | n/a          | 99.46%       | 0.53%       | 0.79%        | 98.79%   |
| gdna_high_ss_0.50_nrna_zero | 66.67%     | 67.03%    | 0.00%      | n/a          | n/a          | 99.37%       | 0.62%       | 2.35%        | 98.46%   |
| gdna_high_ss_0.99_nrna_zero | 66.67%     | 66.87%    | 0.00%      | n/a          | n/a          | 99.55%       | 0.44%       | 1.52%        | 98.84%   |

## Strand-Specific Delta

Each row compares `ss=0.99` minus `ss=0.50` for the same gDNA/nRNA setting. Positive
`d nRNA recall` and negative `d nRNA to gDNA` are improvements.

| gDNA | nRNA | d pool acc | d nRNA recall | d nRNA to gDNA | d gDNA recall | d gDNA to RNA | d assigned gDNA |
| ---- | ---- | ---------- | ------------- | -------------- | ------------- | ------------- | --------------- |
| high | high | 2.60%      | 26.42%        | -26.37%        | -0.08%        | 0.08%         | -3.24%          |
| low  | high | 8.94%      | 45.10%        | -45.05%        | -0.79%        | 0.79%         | -11.02%         |
| med  | high | 4.07%      | 29.18%        | -29.22%        | -0.25%        | 0.25%         | -5.05%          |
| zero | high | 18.45%     | 78.54%        | -78.77%        | n/a           | n/a           | -22.66%         |
| high | low  | 0.86%      | 19.28%        | -19.04%        | 0.12%         | -0.12%        | -0.74%          |
| low  | low  | 1.20%      | 21.44%        | -20.99%        | -0.11%        | 0.11%         | -1.40%          |
| med  | low  | 0.93%      | 18.99%        | -18.53%        | 0.11%         | -0.11%        | -0.88%          |
| zero | low  | 5.27%      | 77.10%        | -76.82%        | n/a           | n/a           | -5.94%          |
| high | zero | 0.39%      | n/a           | n/a            | 0.18%         | -0.18%        | -0.16%          |
| low  | zero | 0.08%      | n/a           | n/a            | 0.18%         | -0.18%        | -0.01%          |
| med  | zero | 0.26%      | n/a           | n/a            | 0.19%         | -0.19%        | -0.08%          |
| zero | zero | 0.00%      | n/a           | n/a            | n/a           | n/a           | 0.00%           |

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
| mrna | nrna | ambig_same_strand | unspliced | 53,066 |
| nrna | mrna | ambig_same_strand | unspliced | 34,280 |
| nrna | gdna | ambig_same_strand | unspliced | 11,079 |
| nrna | mrna | ambig_opp_strand  | unspliced | 10,249 |
| mrna | nrna | ambig_opp_strand  | unspliced | 9,242  |
| nrna | gdna | ambig_opp_strand  | unspliced | 3,540  |
| mrna | gdna | ambig_same_strand | unspliced | 3,019  |
| mrna | gdna | ambig_opp_strand  | unspliced | 903    |

When true gDNA and high nRNA coexist in stranded data, the largest residual errors are:

| True | Pred | ZC                | ZS        | Count  |
| ---- | ---- | ----------------- | --------- | ------ |
| nrna | gdna | ambig_same_strand | unspliced | 37,598 |
| mrna | nrna | ambig_same_strand | unspliced | 37,558 |
| nrna | mrna | ambig_same_strand | unspliced | 34,770 |
| mrna | gdna | ambig_same_strand | unspliced | 14,498 |
| nrna | gdna | ambig_opp_strand  | unspliced | 11,725 |
| nrna | mrna | ambig_opp_strand  | unspliced | 10,183 |
| gdna | nrna | ambig_same_strand | unspliced | 9,045  |
| mrna | nrna | ambig_opp_strand  | unspliced | 7,895  |

## Calibration Diagnostics In High-nRNA Scenarios

When true gDNA is absent or low, the gDNA FL model and genic gDNA densities can be
pulled toward nRNA-like evidence. `rho_ig` is the intergenic anchor; `rho_in` and
`rho_ex_in` are genic unspliced channels.

| Condition                   | SS   | gDNA FL | rho ig | rho intron | rho exon-intron | Pred gDNA | nRNA to gDNA |
| --------------------------- | ---- | ------- | ------ | ---------- | --------------- | --------- | ------------ |
| gdna_zero_ss_0.50_nrna_high | 0.50 | 250.1   | 0      | 0.227      | 0.225           | 24.11%    | 84.06%       |
| gdna_zero_ss_0.99_nrna_high | 0.99 | 250.1   | 0      | 7.79e-05   | 0.0775          | 1.45%     | 5.30%        |
| gdna_low_ss_0.50_nrna_high  | 0.50 | 303.6   | 0.0249 | 0.255      | 0.212           | 28.85%    | 52.60%       |
| gdna_low_ss_0.99_nrna_high  | 0.99 | 303.6   | 0.0249 | 0.0253     | 0.0901          | 17.83%    | 7.55%        |
| gdna_med_ss_0.50_nrna_high  | 0.50 | 332.8   | 0.1    | 0.332      | 0.271           | 50.44%    | 41.80%       |
| gdna_med_ss_0.99_nrna_high  | 0.99 | 332.8   | 0.1    | 0.0997     | 0.156           | 45.39%    | 12.59%       |
| gdna_high_ss_0.50_nrna_high | 0.50 | 341.2   | 0.2    | 0.433      | 0.366           | 65.74%    | 44.25%       |
| gdna_high_ss_0.99_nrna_high | 0.99 | 341.2   | 0.2    | 0.201      | 0.256           | 62.50%    | 17.88%       |

## Predicted-gDNA / No-transcript Bucket

Current annotated BAMs stamp gDNA winners with `ZT=.` and `ZL=-1`, so false gDNA
assignments cannot yet be localized back to the EM locus from the BAM alone. The table
below separates that diagnostic bucket from true biological loci. Preserving the source
locus id for gDNA winners would make future root-cause analysis much sharper.

| Condition                   | True | Pred | Count   |
| --------------------------- | ---- | ---- | ------- |
| gdna_zero_ss_0.50_nrna_high | nrna | gdna | 231,922 |
| gdna_low_ss_0.50_nrna_high  | nrna | gdna | 145,126 |
| gdna_high_ss_0.50_nrna_high | nrna | gdna | 122,084 |
| gdna_med_ss_0.50_nrna_high  | nrna | gdna | 115,333 |
| gdna_zero_ss_0.50_nrna_high | mrna | gdna | 75,723  |
| gdna_high_ss_0.50_nrna_high | mrna | gdna | 50,417  |
| gdna_high_ss_0.99_nrna_high | nrna | gdna | 49,323  |
| gdna_zero_ss_0.50_nrna_low  | nrna | gdna | 48,648  |
| gdna_low_ss_0.50_nrna_high  | mrna | gdna | 47,508  |
| gdna_med_ss_0.50_nrna_high  | mrna | gdna | 42,897  |
| gdna_med_ss_0.99_nrna_high  | nrna | gdna | 34,732  |
| gdna_high_ss_0.50_nrna_low  | nrna | gdna | 31,126  |

## Worst EM Loci By Non-gDNA Pool Error

| Condition                   | Locus | Errors | Err frac | Truth mRNA | Truth nRNA | Truth gDNA | Pred gDNA rate | n tx |
| --------------------------- | ----- | ------ | -------- | ---------- | ---------- | ---------- | -------------- | ---- |
| gdna_zero_ss_0.99_nrna_high | 23    | 22,618 | 13.88%   | 111,203    | 51,697     | 0          | 0.82%          | 15   |
| gdna_high_ss_0.99_nrna_high | 24    | 21,122 | 13.11%   | 110,490    | 48,851     | 1,803      | 5.17%          | 15   |
| gdna_med_ss_0.99_nrna_high  | 24    | 20,840 | 12.85%   | 110,944    | 50,196     | 1,084      | 2.93%          | 15   |
| gdna_low_ss_0.99_nrna_high  | 24    | 20,500 | 12.59%   | 111,209    | 51,268     | 313        | 1.30%          | 15   |
| gdna_high_ss_0.50_nrna_high | 24    | 19,696 | 13.14%   | 107,540    | 40,060     | 2,237      | 11.82%         | 15   |
| gdna_med_ss_0.50_nrna_high  | 24    | 18,897 | 12.56%   | 107,987    | 41,233     | 1,227      | 9.98%          | 15   |
| gdna_low_ss_0.50_nrna_high  | 24    | 16,731 | 11.71%   | 106,716    | 35,847     | 284        | 13.40%         | 15   |
| gdna_high_ss_0.99_nrna_high | 4     | 15,263 | 13.42%   | 83,546     | 28,881     | 1,301      | 4.92%          | 15   |
| gdna_low_ss_0.99_nrna_high  | 4     | 15,169 | 13.22%   | 84,072     | 30,460     | 233        | 1.19%          | 15   |
| gdna_med_ss_0.99_nrna_high  | 4     | 14,973 | 13.09%   | 83,878     | 29,722     | 747        | 2.75%          | 15   |
| gdna_zero_ss_0.99_nrna_high | 4     | 14,264 | 12.44%   | 84,123     | 30,545     | 0          | 0.86%          | 15   |
| gdna_high_ss_0.50_nrna_high | 4     | 13,757 | 12.92%   | 81,329     | 23,596     | 1,544      | 10.99%         | 15   |

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

- Condition metrics: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19/condition_metrics.tsv`
- Pool confusion counts: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19/pool_confusion_counts.tsv`
- Pool confusion matrices: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19/pool_confusion_matrices.tsv`
- Strand-pair deltas: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19/strand_pair_delta.tsv`
- ZC/ZS breakdown: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19/zc_zs_breakdown.tsv`
- ZF breakdown: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19/zf_breakdown.tsv`
- Locus pool breakdown: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19/locus_pool_breakdown.tsv`
- Predicted-gDNA/no-transcript bucket: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19/no_locus_pool_errors.tsv`
- Worst EM locus errors: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_kappa_units_fix_2026-05-19/worst_locus_pool_errors.tsv`
