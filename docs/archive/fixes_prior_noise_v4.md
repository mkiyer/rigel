# Fixes for Prior Noise and gDNA Handoff, v4

Date: 2026-05-27

This document replaces `fixes_prior_noise_v3.md`. v3 correctly rejected the old
latent-state/source-label confusion, but it still drifted toward local patches:
per-region exposure shrinkage constants, a strand-noise cap, and a transcript
floor inside the native EM solver. v4 removes those patches and restates the
problem as a small set of generative modeling contracts.

The core invariant remains:

```text
Latent states are strata: expressed/unexpressed and capture/off-target.
They are not source labels. All strata can contain gDNA.
```

The source split must be represented explicitly by `PriorMassDeconvolution`.
The capture exposure used by EM must be a geometric or low-dimensional exposure
surface, not a noisy per-region count ratio. Transcript allocation remains the
job of EM likelihoods; the grouped RNA prior may constrain aggregate RNA-vs-gDNA
mass, but it must not create transcript-level pseudocounts.

## 1. Diagnostics supporting v4

I added and ran:

```text
scripts/debug/review_prior_noise_v4.py
```

It writes:

```text
/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb/diagnostics/prior_noise_v4_review/
```

Key files:

- `prior_noise_v4_diagnostics.md`
- `kappa_summary.tsv`
- `capture_exposure_geometry.tsv`
- `exact_binomial_tails.tsv`
- `ess_scale_comparison.tsv`

Important facts from the synthetic suite:

| question | diagnostic fact | implication |
|---|---:|---|
| Is `kappa_d = 1e6` accidental? | all high-gDNA conditions hit `strand_channel_kappa_d = 1e6` | the simulated gDNA strand balance is binomial or tighter than binomial |
| Does no-gDNA capture-on have false source mass? | SS=0.99 capture-on probe exons have 158.9 false prior-gDNA counts | source repair must happen before exposure estimation |
| Why are local exposure ratios unstable? | the same probe exons have only 0.0419 off-target expected counts | local ratios explode when the denominator is tiny |
| Is the true capture exposure large? | high-gDNA SS=0.99 capture-on probe exons have true gDNA/off-target ratio 919 | capture exposure must be allowed to be large when supported |
| Is a Normal strand tail adequate at low N? | for `n=10`, `q=0.01`, Normal 3-sigma gives 1.04, exact binomial 0.999 quantile gives 2 | exact predictive tails are required |
| Is `0.25 * N` a small ESS in mega-loci? | at `N=100000`, the cap is 25000 | this is a major behavioral change, not a theorem |

## 2. Decisions

| v3 item | v4 decision |
|---|---|
| Fix A: source-specific per-region `A_r` with `A_R_SOURCE_PRIOR_COUNT = 100` and `A_R_MAX = 1000` | Reject constants and local ratios. Replace with a capture-exposure field fit from probe geometry and source-reliable gDNA evidence. |
| Fix B: excess-over-noise shrinkage with a Normal approximation | Keep the source-noise problem, but use exact predictive RNA-strand tails and treat `kappa_d` as a diagnostic/modeling issue, not a fudge factor. |
| Fix C: locus-scaled ESS `0.25 * N` | Keep only as a benchmarked candidate. Do not treat raw fragment count as the final measure of prior information. |
| Fix D: native RNA floor `a_r / K_rna` | Reject completely. It resurrects transcript-level pseudocounts and violates the grouped-prior design. |

Recommended sequence:

1. Remove Fix D from the plan and add transcript-collapse diagnostics only.
2. Repair strand/source confidence with exact predictive RNA-noise tests.
3. Replace `A_r` with a low-dimensional capture-exposure field.
4. Evaluate ESS scaling as a controlled policy sweep, then merge only if all
   eight scenarios pass.
5. Add simulator overdispersion knobs and rerun the suite under less perfect
   strand balance.

## 3. Rejected Fix D: no transcript-level floor

### Decision

Do not add `floor_i = alpha_rna_add / K_rna` to transcript counts in
`apply_grouped_prior_update`.

The current native comment is exactly the desired contract:

```text
The calibration prior constrains only the aggregate gDNA-vs-RNA split; RNA
pseudocount mass is dynamically distributed in proportion to current RNA
evidence, so no transcript receives a fixed floor.
```

Adding `alpha_rna_add / K_rna` breaks that contract. A transcript with no read
support would receive macroscopic expression simply because it is annotated in a
locus with a large aggregate RNA prior. That is not a numerical guard. It is a
false transcript prior.

### What to do instead

Keep the native grouped-prior redistribution evidence-proportional:

```cpp
out_counts[i] = R * nonnegative_finite(raw_counts[i]) / n_rna;
```

Retain the existing `EM_LOG_EPSILON` guard for logs and normalization. If a true
division-by-zero path exists, add only a machine-scale numerical guard that does
not alter biological mass, for example:

```cpp
constexpr double RNA_NUMERICAL_FLOOR = 1.0e-12;
```

This value must never be multiplied by `alpha_rna_add`, `R`, or locus size. It is
only a finite-arithmetic guard.

### Transcript-collapse diagnostic plan

GENE0008.3 becoming zero is not evidence that EM needs a transcript floor. It is
evidence that a separate isoform-resolution question needs debugging.

Add diagnostics for any transcript that drops from positive baseline count to
zero after a source-prior change:

- raw RNA responsibility before grouped-prior redistribution;
- grouped-prior output count;
- effective length and exposure weight;
- top equivalence classes containing the transcript;
- per-equivalence-class log likelihoods for sibling isoforms;
- whether the transcript has unique spliced, unique unspliced, or only shared
  evidence.

Files:

- `src/rigel/native/em_solver.cpp`, diagnostics only if native visibility is
  needed;
- `src/rigel/pipeline.py`, for locus/transcript debug table emission;
- `scripts/debug/diagnose_transcript_collapse.py`, for post-run analysis.

Acceptance criterion:

```text
No implementation may add macroscopic RNA mass to unsupported transcripts.
```

## 4. Source repair: exact RNA-strand predictive tests

### Problem

The all-region prior-mass path exposes false gDNA source mass in expressed
strand-specific regions. The earlier hard cap was wrong because it removed real
gDNA. The v3 Normal approximation was also too brittle at small `n`.

The correct question is not "how much should we cap gDNA?" The correct question
is:

```text
Is the minor-strand count more extreme than an all-RNA model predicts?
```

If not, the strand channel should not create confident gDNA mass for that region.

### Exact predictive null

For a transcript-stranded region, define:

```text
n = k_sense + k_antisense
k_minor = min-orientation count relative to the fitted RNA protocol
q = minor-orientation probability for RNA fragments
```

The predictive all-RNA tail must be exact at all depths.

Minimum implementation:

```python
from scipy.stats import binom

error_upper = binom.ppf(confidence, n, q)
```

Preferred implementation, because `q` itself is estimated from strand-model
training data:

```python
from scipy.stats import betabinom

error_upper = betabinom.ppf(confidence, n, alpha_q, beta_q)
```

where `alpha_q` and `beta_q` are the posterior parameters for the RNA minor
orientation rate. `StrandSummary` currently stores only `p_r1_sense` and
`n_observations`; v4 should extend it to carry the same/opposite counts or the
minor-rate posterior parameters directly. Reconstructing them from rounded
`p * n` should be avoided except as a temporary diagnostic.

Use the existing calibration confidence level as `confidence`. Do not introduce a
new ad hoc tail probability. If a more stringent confidence is needed, make it a
named calibration config with the same visibility as other calibration settings.

### How this feeds source mass

For each eligible region or compartment:

```python
tail_p = P(K_minor >= k_minor_observed | all RNA)
```

Then:

- If `tail_p` is not significant, the strand data do not support a gDNA source
  call. Set the strand-derived gDNA contribution to zero or mark it low
  reliability so the density/fallback channel must carry the source mass.
- If `tail_p` is significant, keep the existing full strand deconvolution
  posterior. Do not replace it with a capped heuristic.

This is a model-selection or reliability decision around the existing posterior,
not a global `min(gdna, cap)` operation.

### What about `kappa_d = 1e6`?

`kappa_d` is the gDNA strand-balance overdispersion. It is estimated from regions
that look like gDNA. The current estimator returns the maximum when residual
variance is at or below binomial expectation. In the synthetic high-gDNA suite,
that is expected because the simulator appears nearly perfect. It is not by
itself proof of an algorithmic bug.

However, it exposes an important validation gap:

```text
The simulator is too clean to test robustness against real strand-balance
overdispersion.
```

Add simulator knobs for:

- gDNA strand-balance overdispersion, controlling the true beta-binomial
  concentration for DNA strand balance;
- RNA minor-orientation overdispersion, controlling library strand leakage beyond
  a pure binomial model;
- optional region-level strand-bias heterogeneity.

Then rerun the eight-condition suite across at least two overdispersion settings:
near-binomial and moderately overdispersed. This determines whether failures are
caused by an over-perfect simulation, an algorithmic overconfidence problem, or
both.

### Files

- `src/rigel/strand_model.py`: expose minor-orientation posterior parameters.
- `src/rigel/calibration/strand_summary.py`: carry counts or posterior params.
- `src/rigel/calibration/strand_deconv.py`: add exact predictive all-RNA tail
  helper and reliability flag.
- `src/rigel/calibration/calibration_iteration.py`: propagate source flags into
  `PriorMassDeconvolution` diagnostics if needed.
- `tests/test_strand_deconv.py`: exact low-N predictive-tail tests.
- simulator files under `src/rigel` or `scripts/sim`, depending on where the
  current strand simulation lives.

### Tests

1. Low N exactness: `n=10`, `q=0.01`, `confidence=0.999` must give upper tail 2,
   not the continuous 3-sigma value 1.04.
2. Pure RNA: observations inside the predictive all-RNA tail must not produce a
   confident strand-derived gDNA source call.
3. Mixed gDNA: observations far outside the all-RNA tail must preserve the
   existing strand posterior mass.
4. Near-unstranded: if `abs(p - 0.5)` is below the existing contrast floor, this
   reliability test is inactive and the region remains strand-ineligible.
5. Training uncertainty: small strand-training `n_observations` must widen the
   predictive all-RNA tail.

## 5. Exposure repair: fit a capture-exposure field, not local ratios

### Problem

`A_r` is currently derived from a latent-state mixture. That is wrong because
latent states are strata, not source labels. v3 tried to replace it with
source-specific per-region ratios. That is still wrong because tiny denominators
turn small source errors into huge exposure weights.

The diagnostic example is decisive:

```text
No-gDNA SS=0.99 capture-on probe exons:
  prior_gdna_sum = 158.9
  off_target_expected_sum = 0.0419
  ratio = 3795
```

This ratio is not biology. It is a tiny denominator amplifying residual source
noise.

### Principle

`A_r` is not a per-region abundance estimate. It is an exposure surface. For
hybrid capture, exposure should be determined primarily by the probe panel and a
small number of fitted enrichment parameters.

Define:

```text
e_r = rho_off * L_r
c_r = FL-aware capture opportunity in [0, 1]
lambda = panel-level gDNA capture enrichment, lambda >= 1
A_r = 1 + c_r * (lambda - 1)
```

For capture-off data, `c_r = 0` everywhere and `A_r = 1`. For current generated
panels, one `lambda` is enough. Future panels may use a small number of strata,
for example probe-overlap, near-probe shoulder, and off-panel.

This model has no per-region pseudocount and no maximum exposure constant. If the
data imply a large `lambda`, the model reports a large `lambda` with diagnostics.
If the exposure is not identifiable, `lambda` stays at 1 and is flagged.

### Estimation

Estimate `lambda` only after source repair from Section 4. Use source-reliable
gDNA mass, not raw observed RNA+gDNA counts.

For source-reliable regions, fit the one-parameter Poisson or quasi-Poisson
offset model:

```text
d_r ~ Poisson(e_r * (1 + c_r * (lambda - 1)))
```

where `d_r` is the source-deconvolved gDNA mass from `PriorMassDeconvolution` and
`e_r` is the off-target expectation. Regions marked source-unreliable by the
exact RNA-strand predictive test should have zero or very low fit weight.

Implementation should use a bounded one-dimensional optimizer over
`log(lambda - 1)` or an equivalent nonnegative parameterization. The only hard
constraints are mathematical:

```text
lambda >= 1
e_r >= 0
0 <= c_r <= 1
```

There is no `A_R_SOURCE_PRIOR_COUNT`. There is no `A_R_MAX`.

### Identifiability gate

Use a likelihood-ratio or profile-likelihood test against `lambda = 1`. The
confidence level should be the calibration confidence already used by the
calibration package, not a new hidden constant.

If `lambda > 1` is not identifiable:

```text
lambda = 1
A_r = 1
flag = exposure_unidentifiable
```

If it is identifiable:

```text
A_r = 1 + c_r * (lambda_hat - 1)
```

Emit `lambda_hat`, confidence interval, fit weights, and the likelihood-ratio
statistic. Do not clip high values. High values should be explainable from probe
geometry and source mass; if they are not, the source-repair diagnostics should
show why.

### Capture opportunity `c_r`

`c_r` should be deterministic and geometric:

- use the probe BED12 panel when capture is enabled;
- compute FL-aware overlap opportunity for gDNA fragments, not just base-pair
  overlap;
- set `c_r = 0` for capture-off runs;
- support fractional values for boundary and partial-overlap regions.

This aligns exposure with how hybrid capture physically changes sampling. It
also solves the scale-invariance critique: small exons and large introns do not
receive unrelated pseudocounts. They share the same panel-level enrichment and
differ only through their geometric opportunity and effective length.

### Files

- `src/rigel/calibration/calibration_iteration.py`: stop using latent-state
  `mu_gdna` for `A_r`.
- `src/rigel/calibration/exposure.py` or a new `capture_exposure.py`: implement
  `CaptureExposureFit` and `A_r` construction.
- `src/rigel/calibration/_exposure.py`: add FL-aware capture-opportunity helper
  if the existing block averaging is insufficient.
- `src/rigel/index.py` or probe-panel loading code: ensure region/probe overlap
  is available to calibration.
- `tests/test_region_exposure_from_density.py` or a new exposure test module.

### Tests

1. Capture-off: `c_r = 0` everywhere gives `A_r = 1` exactly.
2. No identifiable gDNA source mass: likelihood test fails, `lambda = 1` and
   `A_r = 1` even if local denominators are tiny.
3. Synthetic capture enrichment: regions with `c_r = 1` and source-reliable
   `d_r = 100 * e_r` recover `lambda` near 100.
4. Fractional opportunity: `c_r = 0.25`, `lambda = 101` gives `A_r = 26`.
5. No hidden caps: a valid high-enrichment synthetic case can return
   `lambda > 1000` if the likelihood supports it.

### PR ordering note

Do not merge this exposure change before Section 4 source repair, unless it is
diagnostic-only. The no-gDNA SS=0.99 capture-on case already shows that source
noise can become an absurd exposure ratio when denominator opportunity is tiny.

## 6. ESS scaling: controlled policy, not theorem

### Problem

The current `MAX_ESS = 3000` underweights the source prior in very large loci.
But `0.25 * N` can place 25k effective prior counts into a 100k-fragment locus.
That can be correct if the source prior is high quality and the likelihood has a
known bias. It can also bully local RNA partitioning if source or exposure is
wrong.

### v4 policy

Do not merge a new ESS scaling rule as a standalone fix. Evaluate it only after:

1. source repair is in place;
2. capture exposure is geometric and identifiable;
3. diagnostics report source reliability and prior-posterior drift.

Candidate policies to benchmark:

| policy | formula | purpose |
|---|---|---|
| current | `min(locus_unspliced, 3000)` | baseline behavior |
| linear | `min(locus_unspliced, 0.25 * n_em)` | v3 candidate, strongest intervention |
| information-weighted | `min(locus_unspliced, f(source_reliability, exposure_identifiability, n_em))` | preferred endpoint |

The information-weighted policy should use quantities that already have model
meaning, for example:

- source reliability from the exact RNA-strand predictive test;
- exposure identifiability from the capture-exposure likelihood-ratio test;
- number of fragments with finite competing RNA and gDNA likelihoods, not merely
  total locus fragments.

If the linear policy improves the target scenario but degrades RNA MARD elsewhere,
do not tune the fraction. Promote the information-weighted policy instead.

### Files

- `src/rigel/calibration/adaptive_prior.py`
- `src/rigel/calibration/prior.py`
- `src/rigel/pipeline.py` for diagnostics
- `tests/test_adaptive_prior.py`

### Diagnostics to add before policy change

- `prior_ess_cap_policy`
- `prior_ess_cap_data`
- `prior_ess_cap_source_reliability`
- `prior_ess_cap_exposure_identifiability`
- `prior_ess_final`
- `prior_gdna_share_in`
- `posterior_gdna_share_out`
- `prior_posterior_gdna_drift`

### Acceptance gate

Linear ESS scaling is acceptable only if all eight scenarios maintain or improve
the condition-level metrics and the transcript-collapse diagnostics show no new
class of RNA misallocation.

## 7. Simulator hardening

The current suite is valuable but too clean for strand-overdispersion questions.
The kappa diagnostics show this directly: high-gDNA conditions hit the 1e6
maximum, meaning the observed gDNA strand balance is binomial or tighter.

Add simulation parameters:

```text
gdna_strand_kappa
rna_strand_error_kappa
region_strand_bias_sd
```

Then run at least:

1. near-binomial, matching the current suite;
2. moderately overdispersed, approximating real library and alignment
   heterogeneity.

The goal is to answer two separate questions:

- Is the current failure caused by a simulator that is too perfect?
- Does the algorithm remain calibrated when the simulator becomes less perfect?

The likely answer is both: the simulator is too clean for kappa stress-testing,
and the algorithm still needs exact predictive source-reliability checks.

## 8. Acceptance criteria

Evaluate on the existing eight-condition suite first, then repeat with the new
overdispersed simulator settings.

### Invariants

| invariant | gate |
|---|---|
| No transcript floor | no code path adds `alpha_rna_add / K_rna` or any prior-scaled mass to every transcript |
| Latent states are strata | no adaptive prior or exposure code treats `p_state_*` as source mass |
| Exposure is geometric | `A_r` comes from capture opportunity and fitted enrichment, not per-region observed/raw source ratios |
| Source uncertainty is explicit | exact RNA-strand predictive flags are emitted and summarized |

### Condition-level gates

| scenario group | gate |
|---|---|
| all no-gDNA scenarios | estimated gDNA fraction remains at current tolerance; no capture-on false exposure inflation |
| high-gDNA SS=0.99 capture-on | gDNA fraction moves materially toward 0.5 and mRNA MARD improves |
| high-gDNA SS=0.99 capture-off | no mRNA MARD regression beyond 5 percent relative |
| high-gDNA SS=0.50 capture-off | no mRNA MARD regression beyond 5 percent relative |
| high-gDNA SS=0.50 capture-on | no material regression; improvement requires a separate near-unstranded source channel |

### Diagnostic gates

- Capture-on high-gDNA should have identifiable `lambda > 1` with sensible
  profile-likelihood support.
- Capture-on no-gDNA should have `lambda = 1` or unidentifiable exposure after
  source repair.
- Regions whose minor-strand counts are compatible with all-RNA should not make
  high-confidence strand-derived gDNA source calls.
- Any transcript that drops to zero must have a diagnostic explanation in raw
  likelihood/equivalence-class support, not a hidden prior floor.

## 9. Validation commands

For the diagnostics added in this document:

```bash
conda activate rigel && python scripts/debug/review_prior_noise_v4.py
conda activate rigel && ruff check scripts/debug/review_prior_noise_v4.py
```

For the Python calibration changes:

```bash
conda activate rigel && ruff check src/rigel/calibration tests scripts/debug
conda activate rigel && pytest tests/test_strand_deconv.py tests/test_calibration_iteration.py tests/test_adaptive_prior.py -v
```

If native diagnostics are added but no EM math changes are made, rebuild only if
`src/rigel/native/*.cpp` or headers change:

```bash
conda activate rigel && pip install --no-build-isolation -e .
conda activate rigel && pytest tests/test_batch_em_impl.py -v
```

Full proof requires rerunning the eight-condition hybrid-capture suite and then
the same suite with strand-overdispersion enabled.

## 10. What is explicitly deferred

Near-unstranded high-gDNA capture-on remains a source-identification problem. The
current diagnostics show only about 3.3k prior-global gDNA when truth is 100k.
No exposure surface, ESS cap, or EM update can recover gDNA mass that the source
split never identifies. That requires a separate source channel, likely using
fragment length, capture geometry, and density jointly when strand is
uninformative.

Isoform-level failures are also deferred from the prior-noise fix. They deserve a
separate investigation into equivalence classes, effective lengths, exposure
weights, and raw likelihoods. The grouped prior must remain transcript-neutral.
