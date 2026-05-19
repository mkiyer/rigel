# Synthetic 24 Pool Confusion Report - 2026-05-15

Base directory: `/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24`

Rigel output namespace: `<condition>/rigel_strand_v2_out`

Annotated BAM namespace: `<condition>/annotated_strand_v2.bam`

Analysis artifacts: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15`

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
winner recall is 68.35% and mean nRNA-to-gDNA leakage is
14.75%. This is a major improvement over unstranded nRNA
scenarios, where mean nRNA recall is 30.96% and mean
nRNA-to-gDNA leakage is 52.00%, but it is still not the
dramatic separation we want.

The zero-gDNA, high-nRNA contrast is the cleanest diagnostic. In the unstranded run,
Rigel sends 80.47% of true nRNA fragments
to gDNA. In the stranded run, this falls to
5.66%, and nRNA recall rises from
3.55% to
77.83%. Strand information is being used, but
genic unspliced evidence remains partially gDNA-like.

## Per-Scenario Pool Confusion

The table below reports row-normalized winner assignments for each true pool. `nRNA to
gDNA` means true nascent fragments assigned to the gDNA pool. `gDNA to RNA` combines
gDNA assigned to mature or nascent RNA.

| Condition                   | Truth gDNA | Pred gDNA | Truth nRNA | nRNA to nRNA | nRNA to gDNA | gDNA to gDNA | gDNA to RNA | mRNA to gDNA | Pool acc |
| --------------------------- | ---------- | --------- | ---------- | ------------ | ------------ | ------------ | ----------- | ------------ | -------- |
| gdna_zero_ss_0.50_nrna_high | 0.00%      | 23.37%    | 21.62%     | 3.55%        | 80.47%       | n/a          | n/a         | 7.62%        | 72.45%   |
| gdna_zero_ss_0.99_nrna_high | 0.00%      | 1.51%     | 21.62%     | 77.83%       | 5.66%        | n/a          | n/a         | 0.37%        | 90.97%   |
| gdna_low_ss_0.50_nrna_high  | 16.38%     | 28.67%    | 18.08%     | 31.64%       | 51.78%       | 99.02%       | 0.97%       | 4.71%        | 82.88%   |
| gdna_low_ss_0.99_nrna_high  | 16.38%     | 17.87%    | 18.08%     | 75.63%       | 7.70%        | 98.24%       | 1.75%       | 0.58%        | 91.80%   |
| gdna_med_ss_0.50_nrna_high  | 43.94%     | 50.39%    | 12.12%     | 41.97%       | 41.55%       | 98.97%       | 1.02%       | 4.25%        | 89.38%   |
| gdna_med_ss_0.99_nrna_high  | 43.94%     | 45.38%    | 12.12%     | 71.03%       | 12.56%       | 98.73%       | 1.26%       | 1.10%        | 93.33%   |
| gdna_high_ss_0.50_nrna_high | 61.05%     | 65.72%    | 8.42%      | 39.57%       | 44.12%       | 99.06%       | 0.93%       | 5.01%        | 91.95%   |
| gdna_high_ss_0.99_nrna_high | 61.05%     | 62.51%    | 8.42%      | 65.63%       | 17.98%       | 98.97%       | 1.02%       | 1.88%        | 94.52%   |
| gdna_zero_ss_0.50_nrna_low  | 0.00%      | 5.47%     | 5.59%      | 10.07%       | 71.87%       | n/a          | n/a         | 1.53%        | 92.51%   |
| gdna_zero_ss_0.99_nrna_low  | 0.00%      | 0.38%     | 5.59%      | 76.95%       | 5.41%        | n/a          | n/a         | 0.08%        | 97.08%   |
| gdna_low_ss_0.50_nrna_low   | 19.09%     | 20.98%    | 4.53%      | 50.42%       | 32.00%       | 98.82%       | 1.17%       | 0.87%        | 95.83%   |
| gdna_low_ss_0.99_nrna_low   | 19.09%     | 19.64%    | 4.53%      | 71.14%       | 11.75%       | 98.77%       | 1.22%       | 0.32%        | 97.01%   |
| gdna_med_ss_0.50_nrna_low   | 48.56%     | 50.17%    | 2.88%      | 40.20%       | 41.77%       | 99.08%       | 0.92%       | 1.77%        | 96.35%   |
| gdna_med_ss_0.99_nrna_low   | 48.56%     | 49.30%    | 2.88%      | 58.90%       | 23.56%       | 99.19%       | 0.81%       | 0.94%        | 97.25%   |
| gdna_high_ss_0.50_nrna_low  | 65.38%     | 66.88%    | 1.94%      | 30.24%       | 52.43%       | 99.24%       | 0.75%       | 3.02%        | 96.77%   |
| gdna_high_ss_0.99_nrna_low  | 65.38%     | 66.15%    | 1.94%      | 49.66%       | 33.35%       | 99.37%       | 0.62%       | 1.66%        | 97.63%   |
| gdna_zero_ss_0.50_nrna_zero | 0.00%      | 0.00%     | 0.00%      | n/a          | n/a          | n/a          | n/a         | 0.00%        | 98.89%   |
| gdna_zero_ss_0.99_nrna_zero | 0.00%      | 0.00%     | 0.00%      | n/a          | n/a          | n/a          | n/a         | 0.00%        | 98.89%   |
| gdna_low_ss_0.50_nrna_zero  | 20.00%     | 20.04%    | 0.00%      | n/a          | n/a          | 99.12%       | 0.87%       | 0.27%        | 98.72%   |
| gdna_low_ss_0.99_nrna_zero  | 20.00%     | 20.04%    | 0.00%      | n/a          | n/a          | 99.30%       | 0.69%       | 0.23%        | 98.79%   |
| gdna_med_ss_0.50_nrna_zero  | 50.00%     | 50.22%    | 0.00%      | n/a          | n/a          | 99.29%       | 0.70%       | 1.15%        | 98.53%   |
| gdna_med_ss_0.99_nrna_zero  | 50.00%     | 50.13%    | 0.00%      | n/a          | n/a          | 99.47%       | 0.52%       | 0.80%        | 98.79%   |
| gdna_high_ss_0.50_nrna_zero | 66.67%     | 67.02%    | 0.00%      | n/a          | n/a          | 99.37%       | 0.62%       | 2.32%        | 98.47%   |
| gdna_high_ss_0.99_nrna_zero | 66.67%     | 66.86%    | 0.00%      | n/a          | n/a          | 99.55%       | 0.44%       | 1.50%        | 98.85%   |

## Strand-Specific Delta

Each row compares `ss=0.99` minus `ss=0.50` for the same gDNA/nRNA setting. Positive
`d nRNA recall` and negative `d nRNA to gDNA` are improvements.

| gDNA | nRNA | d pool acc | d nRNA recall | d nRNA to gDNA | d gDNA recall | d gDNA to RNA | d assigned gDNA |
| ---- | ---- | ---------- | ------------- | -------------- | ------------- | ------------- | --------------- |
| high | high | 2.57%      | 26.06%        | -26.15%        | -0.09%        | 0.09%         | -3.21%          |
| low  | high | 8.93%      | 43.99%        | -44.08%        | -0.78%        | 0.78%         | -10.81%         |
| med  | high | 3.95%      | 29.06%        | -29.00%        | -0.24%        | 0.24%         | -5.01%          |
| zero | high | 18.52%     | 74.28%        | -74.81%        | n/a           | n/a           | -21.86%         |
| high | low  | 0.87%      | 19.43%        | -19.08%        | 0.13%         | -0.13%        | -0.73%          |
| low  | low  | 1.18%      | 20.72%        | -20.25%        | -0.05%        | 0.05%         | -1.34%          |
| med  | low  | 0.91%      | 18.70%        | -18.21%        | 0.11%         | -0.11%        | -0.87%          |
| zero | low  | 4.57%      | 66.88%        | -66.46%        | n/a           | n/a           | -5.09%          |
| high | zero | 0.38%      | n/a           | n/a            | 0.18%         | -0.18%        | -0.16%          |
| low  | zero | 0.07%      | n/a           | n/a            | 0.18%         | -0.18%        | -0.00%          |
| med  | zero | 0.25%      | n/a           | n/a            | 0.18%         | -0.18%        | -0.09%          |
| zero | zero | 0.00%      | n/a           | n/a            | n/a           | n/a           | -0.00%          |

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
| mrna | nrna | ambig_same_strand | unspliced | 41,299 |
| nrna | mrna | ambig_same_strand | unspliced | 35,210 |
| nrna | gdna | ambig_same_strand | unspliced | 12,036 |
| nrna | mrna | ambig_opp_strand  | unspliced | 10,269 |
| mrna | nrna | ambig_opp_strand  | unspliced | 9,071  |
| nrna | gdna | ambig_opp_strand  | unspliced | 3,581  |
| mrna | gdna | ambig_same_strand | unspliced | 2,758  |
| mrna | gdna | ambig_opp_strand  | unspliced | 922    |

When true gDNA and high nRNA coexist in stranded data, the largest residual errors are:

| True | Pred | ZC                | ZS        | Count  |
| ---- | ---- | ----------------- | --------- | ------ |
| nrna | gdna | ambig_same_strand | unspliced | 38,031 |
| mrna | nrna | ambig_same_strand | unspliced | 37,020 |
| nrna | mrna | ambig_same_strand | unspliced | 35,007 |
| mrna | gdna | ambig_same_strand | unspliced | 14,628 |
| nrna | gdna | ambig_opp_strand  | unspliced | 11,562 |
| nrna | mrna | ambig_opp_strand  | unspliced | 10,166 |
| gdna | nrna | ambig_same_strand | unspliced | 9,166  |
| mrna | nrna | ambig_opp_strand  | unspliced | 8,024  |

## Calibration Diagnostics In High-nRNA Scenarios

When true gDNA is absent or low, the gDNA FL model and genic gDNA densities can be
pulled toward nRNA-like evidence. `rho_ig` is the intergenic anchor; `rho_in` and
`rho_ex_in` are genic unspliced channels.

| Condition                   | SS   | gDNA FL | rho ig | rho intron | rho exon-intron | Pred gDNA | nRNA to gDNA |
| --------------------------- | ---- | ------- | ------ | ---------- | --------------- | --------- | ------------ |
| gdna_zero_ss_0.50_nrna_high | 0.50 | 250.1   | 0      | 0.227      | 0.225           | 23.37%    | 80.47%       |
| gdna_zero_ss_0.99_nrna_high | 0.99 | 250.1   | 0      | 7.79e-05   | 0.0775          | 1.51%     | 5.66%        |
| gdna_low_ss_0.50_nrna_high  | 0.50 | 303.6   | 0.0249 | 0.255      | 0.212           | 28.67%    | 51.78%       |
| gdna_low_ss_0.99_nrna_high  | 0.99 | 303.6   | 0.0249 | 0.0253     | 0.0901          | 17.87%    | 7.70%        |
| gdna_med_ss_0.50_nrna_high  | 0.50 | 332.8   | 0.1    | 0.332      | 0.271           | 50.39%    | 41.55%       |
| gdna_med_ss_0.99_nrna_high  | 0.99 | 332.8   | 0.1    | 0.0997     | 0.156           | 45.38%    | 12.56%       |
| gdna_high_ss_0.50_nrna_high | 0.50 | 341.2   | 0.2    | 0.433      | 0.366           | 65.72%    | 44.12%       |
| gdna_high_ss_0.99_nrna_high | 0.99 | 341.2   | 0.2    | 0.201      | 0.256           | 62.51%    | 17.98%       |

## Predicted-gDNA / No-transcript Bucket

Current annotated BAMs stamp gDNA winners with `ZT=.` and `ZL=-1`, so false gDNA
assignments cannot yet be localized back to the EM locus from the BAM alone. The table
below separates that diagnostic bucket from true biological loci. Preserving the source
locus id for gDNA winners would make future root-cause analysis much sharper.

| Condition                   | True | Pred | Count   |
| --------------------------- | ---- | ---- | ------- |
| gdna_zero_ss_0.50_nrna_high | nrna | gdna | 222,000 |
| gdna_low_ss_0.50_nrna_high  | nrna | gdna | 142,846 |
| gdna_high_ss_0.50_nrna_high | nrna | gdna | 121,724 |
| gdna_med_ss_0.50_nrna_high  | nrna | gdna | 114,640 |
| gdna_zero_ss_0.50_nrna_high | mrna | gdna | 76,169  |
| gdna_high_ss_0.50_nrna_high | mrna | gdna | 50,071  |
| gdna_high_ss_0.99_nrna_high | nrna | gdna | 49,593  |
| gdna_low_ss_0.50_nrna_high  | mrna | gdna | 47,133  |
| gdna_zero_ss_0.50_nrna_low  | nrna | gdna | 42,592  |
| gdna_med_ss_0.50_nrna_high  | mrna | gdna | 42,511  |
| gdna_med_ss_0.99_nrna_high  | nrna | gdna | 34,638  |
| gdna_high_ss_0.50_nrna_low  | nrna | gdna | 31,073  |

## Worst EM Loci By Non-gDNA Pool Error

| Condition                   | Locus | Errors | Err frac | Truth mRNA | Truth nRNA | Truth gDNA | Pred gDNA rate | n tx |
| --------------------------- | ----- | ------ | -------- | ---------- | ---------- | ---------- | -------------- | ---- |
| gdna_high_ss_0.99_nrna_high | 24    | 21,438 | 13.31%   | 110,468    | 48,771     | 1,835      | 5.21%          | 15   |
| gdna_med_ss_0.99_nrna_high  | 24    | 20,809 | 12.82%   | 110,932    | 50,219     | 1,114      | 2.91%          | 15   |
| gdna_low_ss_0.99_nrna_high  | 24    | 20,305 | 12.47%   | 111,200    | 51,268     | 327        | 1.30%          | 15   |
| gdna_zero_ss_0.99_nrna_high | 23    | 19,959 | 12.26%   | 111,236    | 51,530     | 0          | 0.90%          | 15   |
| gdna_high_ss_0.50_nrna_high | 24    | 19,761 | 13.19%   | 107,572    | 40,005     | 2,293      | 11.80%         | 15   |
| gdna_med_ss_0.50_nrna_high  | 24    | 19,113 | 12.70%   | 108,091    | 41,158     | 1,223      | 9.96%          | 15   |
| gdna_low_ss_0.50_nrna_high  | 24    | 16,690 | 11.70%   | 106,765    | 35,653     | 276        | 13.49%         | 15   |
| gdna_high_ss_0.99_nrna_high | 4     | 15,198 | 13.35%   | 83,568     | 28,968     | 1,311      | 4.82%          | 15   |
| gdna_med_ss_0.99_nrna_high  | 4     | 14,968 | 13.07%   | 83,867     | 29,864     | 748        | 2.64%          | 15   |
| gdna_low_ss_0.99_nrna_high  | 4     | 14,643 | 12.76%   | 84,080     | 30,447     | 226        | 1.20%          | 15   |
| gdna_zero_ss_0.99_nrna_high | 4     | 14,085 | 12.28%   | 84,130     | 30,590     | 0          | 0.81%          | 15   |
| gdna_high_ss_0.50_nrna_high | 4     | 13,904 | 13.04%   | 81,407     | 23,722     | 1,511      | 10.85%         | 15   |

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

- Condition metrics: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15/condition_metrics.tsv`
- Pool confusion counts: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15/pool_confusion_counts.tsv`
- Pool confusion matrices: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15/pool_confusion_matrices.tsv`
- Strand-pair deltas: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15/strand_pair_delta.tsv`
- ZC/ZS breakdown: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15/zc_zs_breakdown.tsv`
- ZF breakdown: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15/zf_breakdown.tsv`
- Locus pool breakdown: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15/locus_pool_breakdown.tsv`
- Predicted-gDNA/no-transcript bucket: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15/no_locus_pool_errors.tsv`
- Worst EM locus errors: `/Users/mkiyer/proj/rigel/results/synthetic_24_pool_confusion_2026-05-15/worst_locus_pool_errors.tsv`
