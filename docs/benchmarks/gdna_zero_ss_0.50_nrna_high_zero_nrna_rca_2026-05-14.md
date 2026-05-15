# RCA: zero synthetic nRNA assignment in `gdna_zero_ss_0.50_nrna_high`

Date: 2026-05-14

## Question

The condition `gdna_zero_ss_0.50_nrna_high` has true counts:

- mRNA: 1,000,000 fragments
- nRNA: 275,885 fragments
- gDNA: 0 fragments

The surprising result was that no true nRNA fragments were assigned to synthetic
nRNA in the default Rigel output. This looked like a possible unstranded nRNA
routing bug.

## Short answer

This is not caused by missing synthetic nRNA candidates. The synthetic nRNA
parent candidate is present for every true nRNA fragment that enters EM.

The exact zero in the default discrete output is caused by a real model collapse
plus the post-EM assignment floor:

1. EM gives synthetic nRNA nonzero but tiny posterior mass: 38.03 expected
   synthetic nRNA fragments out of 275,885 true nRNA fragments.
2. No true nRNA fragment has synthetic nRNA as its MAP component.
3. Reconstructed per-fragment posteriors show the largest origin synthetic nRNA
   posterior is 0.001225, and the largest summed synthetic nRNA posterior across
   all synthetic candidates in a fragment is 0.002923.
4. The default `assignment_mode=sample` and `assignment_min_posterior=0.01`
   zero all components below 1 percent before sampling. Synthetic nRNA is below
   that floor everywhere, so the sampled synthetic nRNA output is exactly zero.

## Key evidence

Default annotated-BAM truth assignment for true nRNA fragments:

| true origin | assigned pool | fragments | fraction |
|---|---:|---:|---:|
| nRNA | gDNA | 232,541 | 84.29% |
| nRNA | mRNA | 43,282 | 15.69% |
| nRNA | unresolved | 62 | 0.02% |
| nRNA | synthetic nRNA | 0 | 0.00% |

By region class, intronic and exon-intron true nRNA are almost always sent to
gDNA, while exon-contained true nRNA is mostly sent to mRNA:

| region class | dominant output | fragments |
|---|---:|---:|
| intronic | gDNA | 202,973 of 204,139 |
| exon-intron | gDNA | 22,601 of 22,701 |
| exon-contained | mRNA | 42,078 of 49,045 |

Candidate-set reconstruction ruled out candidate omission:

| assigned pool | region class | n | origin synthetic nRNA present | any synthetic candidate | gDNA candidate |
|---|---|---:|---:|---:|---:|
| gDNA | exon-contained | 6,967 | 100% | 100% | 100% |
| gDNA | exon-intron | 22,601 | 100% | 100% | 100% |
| gDNA | intronic | 202,973 | 100% | 100% | 100% |
| mRNA | exon-contained | 42,078 | 100% | 100% | 100% |
| mRNA | exon-intron | 52 | 100% | 100% | 100% |
| mRNA | intronic | 1,152 | 100% | 100% | 100% |

Assignment-mode controls:

| mode | posterior floor | synthetic nRNA total |
|---|---:|---:|
| fractional | ignored | 38.029 |
| sample, seed 12345 | 0.00 | 48 |
| sample, seed 12345 | 0.01 | 0 |
| map | 0.00 | 0 |

Posterior reconstruction from the fractional EM state:

| assigned pool | region class | max origin nRNA posterior | max summed synthetic nRNA posterior | synthetic MAP fragments |
|---|---|---:|---:|---:|
| gDNA | exon-contained | 0.000733 | 0.001965 | 0 |
| gDNA | exon-intron | 0.001150 | 0.002652 | 0 |
| gDNA | intronic | 0.001225 | 0.002923 | 0 |
| mRNA | exon-contained | 0.000356 | 0.001060 | 0 |
| mRNA | exon-intron | 0.000123 | 0.000470 | 0 |
| mRNA | intronic | 0.000414 | 0.002035 | 0 |

## Mechanism

This condition is unstranded. The trained strand signal is essentially 0.5:

- spliced strand model: p_r1_sense = 0.4975, specificity = 0.5025
- exonic diagnostic: p_r1_sense = 0.4994, specificity = 0.5006

With no strand information, most true nRNA fragments look like unspliced genomic
coverage. The calibration prior then infers substantial gDNA even though the
simulation truth has zero gDNA:

- calibration mean_pi_gdna = 0.2701
- fractional gDNA EM total = 320,048.16
- locus 15 gDNA prior count = 12,699.80
- locus 23 gDNA prior count = 6,342.95
- locus 41 gDNA prior count = 4,286.18

The native assignment math also makes this especially hard for synthetic nRNA.
For transcript-like RNA components, assignment uses:

```text
log_likelihood + log(theta) - log(effective_length)
```

The gDNA component is added as a locus-level component with its own gDNA
likelihood and no transcript effective-length term. Synthetic nRNA spans are
large. In locus 15, synthetic nRNA effective lengths are about 42,000 bp, so
`-log(effective_length)` is about -10.65 for each synthetic nRNA component.

Example true nRNA fragment:

```text
qname: nrna_GENE0018.3:37629-37711:f:0
assigned: gDNA
winner posterior tag: 0.999967
origin synthetic nRNA candidate: present
RNA log_lik for synthetic nRNA candidates: about -11.75
synthetic nRNA effective length: about 42,223
gDNA log_lik: about -21.36
```

Although the raw gDNA log likelihood is worse than the synthetic nRNA likelihood,
the synthetic nRNA effective-length penalty and the very large gDNA mixture
weight make gDNA dominate. For exon-contained true nRNA fragments, mRNA and
synthetic nRNA often have the same likelihood, but mRNAs have much shorter
effective lengths and much larger learned mixture weights, so mRNA dominates.

## Source locations

- `src/rigel/native/em_solver.cpp`: `assign_posteriors` thresholds components
  below `assignment_min_posterior` for `sample` and `map`, but not for
  `fractional`.
- `src/rigel/native/em_solver.cpp`: transcript-like candidates use
  `log(theta) - log_eff_len` in assignment.
- `src/rigel/native/em_solver.cpp`: VBEM baseline prior is 0.5, while gDNA gets
  0.5 plus the calibration-derived `gdna_prior_count`.

## Artifacts

Focused diagnostic scripts:

- `scripts/debug/interrogate_zero_nrna_high.py`
- `scripts/debug/inspect_zero_nrna_candidate_sets.py`
- `scripts/debug/gdna_prior_scale_sweep_zero_nrna.py`

Focused outputs:

- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/origin_pool_counts.tsv`
- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/true_nrna_by_region_pool_zc.tsv`
- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/true_nrna_candidate_presence_summary.tsv`
- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/true_nrna_reconstructed_posterior_summary.tsv`
- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/example_candidate_sets.tsv`
- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/gdna_prior_scale_sweep.tsv`
- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/rigel_fractional_debug/`
- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/rigel_sample_min0_debug/`
- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/rigel_sample_seed12345_default_threshold_debug/`
- `results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high/rigel_map_min0_debug/`

## Interpretation

The exact zero is an output discretization effect, but the underlying posterior
collapse is real. Rigel is not failing to route true nRNA fragments to synthetic
nRNA candidates. It is routing them correctly, then assigning almost all
posterior mass to gDNA or mRNA because unstranded high-nRNA coverage is
structurally confounded with gDNA and because synthetic nRNA components carry a
large effective-length penalty while gDNA receives a large calibration prior.

## gDNA Prior Scale Sweep

A counterfactual fractional-EM sweep scaled only the calibration-derived gDNA
prior while keeping the scored fragments, gDNA likelihoods, and RNA likelihoods
fixed. This tests whether the zero-nRNA behavior is simply caused by an oversized
Dirichlet pseudocount.

| gDNA prior scale | nRNA count | gDNA count | gDNA rate |
|---:|---:|---:|---:|
| 0.00 | 56,069.5 | 224,958.5 | 17.63% |
| 0.01 | 31,448.6 | 250,865.6 | 19.66% |
| 0.10 | 559.0 | 287,232.9 | 22.51% |
| 0.50 | 51.3 | 303,118.1 | 23.76% |
| 1.00 | 38.0 | 320,048.2 | 25.09% |

Conclusion: the gDNA prior is too strong for this no-gDNA truth case, but it is
not the sole cause. Even with the gDNA prior scaled to zero, the gDNA likelihood
and the existing EM state still siphon roughly 225k fragments, while synthetic
nRNA recovers only about 56k of 275,885 true nRNA fragments. Lowering the gDNA
prior alone would not solve unstranded nRNA recovery; a useful fix needs either
better calibration evidence, a credible nRNA prior, or both.

## Improvement opportunities

1. Add a diagnostic output separating fractional posterior mass from sampled
   discrete BAM assignments, especially for nRNA/gDNA calls.
2. Revisit whether synthetic nRNA and gDNA are normalized comparably in the
   assignment objective. The current synthetic nRNA effective-length penalty is
   decisive in unstranded high-nRNA loci.
3. Make the calibration prior nRNA-aware for unstranded data. In this scenario,
   the gDNA prior is trained on true nRNA signal because strand cannot separate
   the pools.
4. Consider reporting low-posterior nRNA evidence instead of dropping it under
   the default 1 percent discrete-assignment floor.