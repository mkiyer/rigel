# Rigel Synthetic 24 Performance Report - 2026-05-14

Base directory: `/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24`

Evaluator run: `conda activate rigel && python scripts/sim/evaluate_suite.py --sim-base /Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24`

Deep-analysis outputs: `/Users/mkiyer/proj/rigel/results/synthetic_24_deep_analysis`

## Executive Summary

All 24 scenarios completed. The no-nRNA cases are mostly healthy: across nonzero-gDNA,
zero-nRNA conditions, the maximum absolute gDNA fraction error was
0.49%, gDNA-to-RNA leak was at most
0.87%, and mean mRNA correct-gene assignment
was 97.89%. The main residual no-nRNA weakness is
isoform resolution: exact transcript assignment averaged only
70.00%, while gene-level WAPE averaged
0.47%.

The major failure mode is nRNA/gDNA confounding. In unstranded high-nRNA data with no
true gDNA, Rigel assigned 25.11% of all fragments
to gDNA and recovered only 0.00% of true nRNA as nRNA.
With strand-specificity 0.99 on the same scenario, assigned gDNA fell to
1.89% and nRNA recall rose to
76.38%, so the strand-aware correction is doing real work but
is not sufficient for exon-boundary nRNA/gDNA ambiguity.

## Key Metrics

| Condition                   | Truth gDNA | Assigned gDNA | Truth nRNA | nRNA recall | nRNA to gDNA | gDNA FL | Tx WAPE | Gene WAPE |
| --------------------------- | ---------- | ------------- | ---------- | ----------- | ------------ | ------- | ------- | --------- |
| gdna_zero_ss_0.50_nrna_low  | 0.00%      | 6.49%         | 5.59%      | 0.00%       | 82.12%       | 250.2   | 11.77%  | 0.96%     |
| gdna_zero_ss_0.99_nrna_low  | 0.00%      | 0.64%         | 5.59%      | 73.38%      | 8.98%        | 250.2   | 3.91%   | 0.48%     |
| gdna_zero_ss_0.50_nrna_high | 0.00%      | 25.11%        | 21.62%     | 0.00%       | 84.29%       | 250.1   | 21.55%  | 4.45%     |
| gdna_zero_ss_0.99_nrna_high | 0.00%      | 1.89%         | 21.62%     | 76.38%      | 6.99%        | 250.1   | 10.99%  | 0.97%     |

## Healthy Baseline: No nRNA

| Condition                   | Truth gDNA | Assigned gDNA | rho_ex/rho_ig | gDNA to RNA | mRNA correct gene |
| --------------------------- | ---------- | ------------- | ------------- | ----------- | ----------------- |
| gdna_high_ss_0.50_nrna_zero | 66.67%     | 67.16%        | 0.980         | 0.60%       | 96.07%            |
| gdna_high_ss_0.99_nrna_zero | 66.67%     | 66.93%        | 0.972         | 0.42%       | 97.26%            |
| gdna_low_ss_0.50_nrna_zero  | 20.00%     | 20.10%        | 0.978         | 0.87%       | 98.31%            |
| gdna_low_ss_0.99_nrna_zero  | 20.00%     | 20.08%        | 0.981         | 0.67%       | 98.61%            |
| gdna_med_ss_0.50_nrna_zero  | 50.00%     | 50.34%        | 0.986         | 0.69%       | 97.33%            |
| gdna_med_ss_0.99_nrna_zero  | 50.00%     | 50.20%        | 0.977         | 0.50%       | 98.00%            |

Interpretation: calibration is coherent when nRNA is absent. The exon/intergenic density
ratio stays near 1.0, and gDNA leakage remains below the 1.5% acceptance guardrail.
This argues against a broad gDNA calibration regression.

## Issue 1: nRNA Is Treated As gDNA In Unstranded Libraries

In the zero-gDNA, low-nRNA unstranded case, Rigel assigned
6.49% of fragments to gDNA even though true gDNA
is zero. In the matched stranded run, assigned gDNA dropped to
0.64%. The high-nRNA case amplifies the same
pattern: 84.29% of true nRNA fragments went to gDNA
when unstranded, versus 6.99% when stranded.

Root cause: the calibration gDNA channels are region-mask based. The gDNA FL histogram
explicitly sums unspliced intron-only, exon-intron, and intergenic-only fragments in
[src/rigel/calibration/_fl_sources.py](src/rigel/calibration/_fl_sources.py). Global
densities are then estimated for INTERGENIC, INTRON, and EXON-INTRON channels in
[src/rigel/calibration/density_global.py](src/rigel/calibration/density_global.py).
Unstranded nRNA produces intronic and exon-boundary unspliced evidence with the same
observable geometry as genic gDNA, but without intergenic support. With no usable strand
contrast, the intron/exon channels are left uncorrected and become gDNA priors.

## Issue 2: Strand Correction Helps, But Exonic Boundary nRNA Still Leaks

For high nRNA with no gDNA, the 0.99 stranded run has intronic density almost zeroed by
strand correction, but still assigns 1.89% of all
fragments to gDNA. The remaining signal is concentrated in EXON-INTRON density: this is
consistent with boundary-crossing unspliced nRNA that is harder to separate from gDNA
even with strand information.

Improvement opportunity: add an explicit nRNA-aware calibration state before global gDNA
projection, or gate genic gDNA priors on intergenic support when intron/exon density is
strongly inconsistent with intergenic density. A practical diagnostic is the ratio
rho_in/rho_ig and rho_ex/rho_ig stratified by strand activity; this suite shows ratios
above 8x in the worst unstranded high-nRNA low-gDNA condition.

## Issue 3: The gDNA Fragment-Length Model Is Contaminated By nRNA

When true gDNA is absent but nRNA is present, the reported gDNA FL is RNA-like
(250.1 bp in high nRNA, no gDNA). In nonzero-gDNA,
high-nRNA conditions, the gDNA FL estimate falls as low as
303.6 bp, which is
-13.26% relative error against the simulated
350 bp gDNA mean. This follows directly from the same source selection in
`extract_gdna_counts`: nRNA fragments enter the unspliced intron/exon-intron pools used
to build the gDNA FL model.

Improvement opportunity: estimate gDNA FL primarily from intergenic-only fragments when
intergenic evidence is adequate, then use intron/exon channels for density after strand
and nRNA correction. The current pooled gDNA FL source is too permissive under nRNA.

## Issue 4: Transcript-Level Isoform Resolution Remains Weak

In the cleanest no-contamination condition, structural ambiguity explains much of the
transcript error:

| Structure class   | n_tx | Truth count | WAPE   | Median RE | FP mass |
| ----------------- | ---- | ----------- | ------ | --------- | ------- |
| same_intron_chain | 139  | 497,867     | 19.49% | 18.59%    | 36,843  |
| shared_nrna_span  | 100  | 423,041     | 1.70%  | 8.05%     | 179     |
| singleton         | 14   | 79,092      | 0.02%  | 0.00%     | 0       |

Fragment-level truth agrees with this: mature RNA is assigned to the correct gene at
about 99% in clean stranded data, but exact transcript assignment is only about 70%.
This is not primarily a cross-gene mapping failure; it is within-gene isoform
identifiability and EM mass splitting.

Top false-positive transcript counts:

| Condition                   | Transcript | Gene     | Pred count | Class             | Gene n_tx |
| --------------------------- | ---------- | -------- | ---------- | ----------------- | --------- |
| gdna_zero_ss_0.50_nrna_high | GENE0018.6 | GENE0018 | 31,182     | same_intron_chain | 10        |
| gdna_med_ss_0.50_nrna_low   | GENE0018.6 | GENE0018 | 28,072     | same_intron_chain | 10        |
| gdna_zero_ss_0.99_nrna_high | GENE0018.6 | GENE0018 | 28,032     | same_intron_chain | 10        |
| gdna_low_ss_0.50_nrna_high  | GENE0018.6 | GENE0018 | 26,433     | same_intron_chain | 10        |
| gdna_high_ss_0.99_nrna_high | GENE0018.6 | GENE0018 | 25,933     | same_intron_chain | 10        |
| gdna_high_ss_0.50_nrna_high | GENE0018.6 | GENE0018 | 24,685     | same_intron_chain | 10        |
| gdna_med_ss_0.50_nrna_high  | GENE0018.6 | GENE0018 | 24,007     | same_intron_chain | 10        |
| gdna_high_ss_0.50_nrna_low  | GENE0018.6 | GENE0018 | 23,537     | same_intron_chain | 10        |
| gdna_low_ss_0.50_nrna_low   | GENE0018.6 | GENE0018 | 22,162     | same_intron_chain | 10        |
| gdna_zero_ss_0.99_nrna_zero | GENE0018.6 | GENE0018 | 21,801     | same_intron_chain | 10        |
| gdna_zero_ss_0.50_nrna_zero | GENE0018.6 | GENE0018 | 21,794     | same_intron_chain | 10        |
| gdna_zero_ss_0.50_nrna_low  | GENE0018.6 | GENE0018 | 21,507     | same_intron_chain | 10        |

Top absolute transcript errors:

| Condition                   | Transcript | Truth  | Pred   | Abs err | Class             |
| --------------------------- | ---------- | ------ | ------ | ------- | ----------------- |
| gdna_high_ss_0.99_nrna_high | GENE0018.2 | 69,220 | 37,350 | 31,870  | same_intron_chain |
| gdna_zero_ss_0.50_nrna_high | GENE0018.6 | 0      | 31,182 | 31,182  | same_intron_chain |
| gdna_med_ss_0.50_nrna_low   | GENE0018.2 | 69,220 | 38,991 | 30,229  | same_intron_chain |
| gdna_low_ss_0.50_nrna_high  | GENE0018.2 | 69,220 | 40,583 | 28,637  | same_intron_chain |
| gdna_zero_ss_0.99_nrna_high | GENE0018.2 | 69,220 | 41,125 | 28,095  | same_intron_chain |
| gdna_med_ss_0.50_nrna_low   | GENE0018.6 | 0      | 28,072 | 28,072  | same_intron_chain |
| gdna_zero_ss_0.99_nrna_high | GENE0018.6 | 0      | 28,032 | 28,032  | same_intron_chain |
| gdna_high_ss_0.50_nrna_high | GENE0018.2 | 69,220 | 41,299 | 27,921  | same_intron_chain |
| gdna_zero_ss_0.99_nrna_zero | GENE0018.2 | 69,220 | 41,636 | 27,584  | same_intron_chain |
| gdna_zero_ss_0.50_nrna_zero | GENE0018.2 | 69,220 | 41,712 | 27,508  | same_intron_chain |
| gdna_med_ss_0.50_nrna_high  | GENE0018.2 | 69,220 | 42,223 | 26,997  | same_intron_chain |
| gdna_low_ss_0.50_nrna_high  | GENE0018.6 | 0      | 26,433 | 26,433  | same_intron_chain |

Improvement opportunity: interrogate EM responsibilities within same-intron-chain and
shared-nRNA-span transcript groups. The likely fixes are not in BAM resolution; they are
in isoform priors, effective-length normalization for short/near-duplicate isoforms, and
posterior regularization when the data cannot identify a unique isoform.

## Recommended Follow-Up Work

1. Build an nRNA-aware SRD calibration model with separate genic-unspliced nRNA and gDNA
   states, and use intergenic evidence as the anchor for true gDNA.
2. Change gDNA FL training to prefer intergenic-only fragments, falling back to genic
   unspliced evidence only after nRNA/strand correction.
3. Add acceptance tests for zero-gDNA plus nRNA scenarios. The current acceptance checks
   pass no-nRNA gDNA cases but do not guard against nRNA-induced false gDNA.
4. Add per-locus reports comparing truth nRNA, predicted nRNA, and predicted gDNA. The
   worst loci in this suite should become regression fixtures.
5. Profile same-intron-chain transcript groups: record per-fragment posterior entropy,
   count_unambig vs count_em, and short-transcript enrichment among false positives.

## Artifacts

- Stock evaluator report: `/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24/analysis_report.txt`
- Condition metrics: `/Users/mkiyer/proj/rigel/results/synthetic_24_deep_analysis/condition_metrics.tsv`
- Transcript detail: `/Users/mkiyer/proj/rigel/results/synthetic_24_deep_analysis/transcript_detail.tsv`
- Gene detail: `/Users/mkiyer/proj/rigel/results/synthetic_24_deep_analysis/gene_detail.tsv`
- Structural summary: `/Users/mkiyer/proj/rigel/results/synthetic_24_deep_analysis/structural_summary.tsv`
- Top false positives: `/Users/mkiyer/proj/rigel/results/synthetic_24_deep_analysis/top_false_positive_transcripts.tsv`
- Top absolute errors: `/Users/mkiyer/proj/rigel/results/synthetic_24_deep_analysis/top_abs_error_transcripts.tsv`
