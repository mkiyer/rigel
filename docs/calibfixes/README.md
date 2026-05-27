# Calibration Prior-Noise Implementation PR Sequence

Source design: `fixes_prior_noise_v5.md`.

This directory breaks the v5 design into implementation-ready PR recipes. The goal is to make each code PR small enough to review and benchmark while preserving the statistical contracts that motivated v5.

## Non-negotiable contracts

| Contract | Implementation implication |
|---|---|
| Latent states are strata | `unexpressed_offtarget`, `unexpressed_capture`, `expressed_capture`, and `expressed_offtarget` must never be consumed as RNA/gDNA source labels. |
| Source mass is explicit | EM priors consume `PriorMassDeconvolution.gdna_unspliced_mean` and `.rna_unspliced_mean`. |
| Exposure is opportunity times enrichment | `A_r` is built from capture geometry/opportunity and learned enrichment, never from raw local abundance ratios. |
| RNA prior is grouped | Do not add transcript-level pseudocount floors in Python or native EM. |
| Uncertainty travels | Source reliability, capture exposure identifiability, and ESS strength must be emitted in diagnostics. |

## Merge order

1. [pr01_continuous_strand_reliability.md](pr01_continuous_strand_reliability.md)
2. [pr02_empirical_bayes_capture_exposure.md](pr02_empirical_bayes_capture_exposure.md)
3. [pr03_ess_policy_sweep.md](pr03_ess_policy_sweep.md)
4. [pr04_simulator_overdispersion.md](pr04_simulator_overdispersion.md)

PR 1 and PR 2 are the actual calibration repair. PR 3 must not change defaults until PR 1 and PR 2 have benchmark evidence. PR 4 improves stress-test realism and can land before or after PR 3, but its new overdispersed scenarios should be used to validate any final ESS policy.

## Deferred follow-up recipes

The v5 document deliberately deferred two larger topics. They are captured separately so they do not sneak into the calibration repair PRs:

- [deferred01_unstranded_capture_fl_boundary_source_split.md](deferred01_unstranded_capture_fl_boundary_source_split.md)
- [deferred02_isoform_identifiability_diagnostics.md](deferred02_isoform_identifiability_diagnostics.md)

## Required scenario coverage

Every behavior-changing PR must be evaluated across the four RNA-seq strata:

| Strand-specific | Capture | PR 1 expectation | PR 2 expectation | Deferred? |
|---|---:|---|---|---|
| yes | no | improve or maintain source split | `A_r = 1` non-regression | no |
| yes | yes | improve source reliability | local EB exposure active when panel/opportunity is available | no |
| no | no | strand reliability inactive | `A_r = 1` non-regression | no |
| no | yes | strand reliability inactive | exposure can be estimated, but source split is not fully solved | yes, source split deferred |

The existing eight-condition hybrid-capture synthetic suite remains the acceptance gate after PR 1 + PR 2. PR 3 has a stricter gate: no material regression in any non-target stratum and no new transcript-collapse pattern.

## General validation commands

Use focused tests for each PR, then run the relevant synthetic suite. For documentation-only changes in this directory:

```bash
conda activate rigel && git --no-pager diff --check -- docs/calibfixes
```

For Python implementation PRs, use the smallest command that covers the touched files, then broaden before merge:

```bash
conda activate rigel && ruff check src/rigel tests/
conda activate rigel && pytest tests/test_strand_deconv.py tests/test_calibration_iteration.py -v
conda activate rigel && pytest tests/test_capture_exposure.py tests/test_per_locus_gdna_mass.py -v
```

If any native C++ files under `src/rigel/native/` are edited in a future PR, rebuild before testing:

```bash
conda activate rigel && pip install --no-build-isolation -e .
```
