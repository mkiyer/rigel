# Fragment-Length Likelihood Asymmetry in gDNA/RNA Competition

Date: 2026-05-26

## Executive Summary

The zero-gDNA single-exon and wide-intron sentinel failures are caused by a
fragment-length (FL) likelihood asymmetry between the RNA and gDNA components.

The core problem is not the adaptive prior, boundary flux, gDNA exposure, RNA
strand likelihood, overhang, or effective length. The native scorer compares two
different FL surfaces:

- RNA uses a weak annotated-spliced-only FL model.
- gDNA uses the fallback/global FL model when no gDNA FL evidence exists.

In low-splicing or single-exon scenarios, the RNA FL model can be trained from a
small and unrepresentative spliced subset, while the gDNA fallback model is the
full library distribution. When an unspliced exonic fragment has equal RNA and
genomic footprint length, the gDNA component can receive a better FL likelihood
solely because its fallback model is better supported. This turns absence of
gDNA evidence into a gDNA likelihood advantage.

The worked examples below show this exactly. In both cases, replacing the RNA FL
term with the same global FL surface used by fallback gDNA collapses false gDNA
posterior mass to nearly zero.

## Relevant Code Paths

### FL Source Construction

`src/rigel/calibration/_fl_sources.py`

This file defines the raw FL count sources used by v6 calibration:

```python
def extract_global_counts(scan_trained):
    return scan_trained.global_model.counts

def extract_rna_counts(scan_trained):
    return scan_trained.category_models[SpliceType.SPLICED_ANNOT].counts

def extract_gdna_counts(payload):
    return gdna_fl_mass(payload)
```

Important consequence:

- `global` is every observed unique-mapper fragment.
- `rna` is annotated-spliced fragments only.
- `gdna` is intergenic plus intronic calibration FL mass.

The spliced-only RNA source is biologically clean because gDNA cannot create an
annotated splice junction. But it is not always statistically representative of
the library's RNA insert-length distribution, especially for single-exon or
low-splicing data.

### FL Model Finalization

`src/rigel/calibration/fl.py`

This file builds three finalized FL models:

- `global_`
- `rna`
- `gdna`

The current weak-pool policy is:

```python
def _adaptive_prior_ess(prior_ess, pool_weight):
    if prior_ess <= 0.0 or pool_weight <= 0.0:
        return 0.0
    return min(prior_ess, np.log1p(pool_weight))
```

and:

```python
if n >= good_threshold:
    return empirical_model, "good"
if n > 0.0 and n >= weak_threshold:
    return EB_shrunk_to_global, "weak"
return global_model, "fallback"
```

This is symmetric in code, but not symmetric in evidence. In a zero-gDNA sample:

- gDNA often has `n_gdna = 0`, so it identity-shares the global model.
- RNA can have `n_rna > 0` but very small, so it gets a weak empirical model
  whose deviations from global are treated as real class evidence.

That means a no-evidence gDNA model can be more accurate than the weak RNA model.

### Native Fragment Scoring

`src/rigel/native/scoring.cpp`

For unspliced non-multimapper EM-routed fragments, RNA scoring is effectively:

```cpp
log_lik = log_strand + frag_len_log_lik(flen)
        + overhang * overhang_log_penalty + nm * mismatch_log_penalty;
```

gDNA scoring is effectively:

```cpp
gdna_ll = gdna_frag_len_log_lik(genomic_footprint)
        + gdna_log_splice_penalty + LOG_HALF
        + nm * mismatch_log_penalty;
```

For contained single-exon/unspliced fragments in the sentinel examples:

- `flen == genomic_footprint`
- overhang is 0
- mismatch is 0
- gDNA splice penalty for unspliced is 0 in log space
- RNA strand term favors RNA under the perfect-stranded setup

So the raw gDNA-minus-RNA score is mostly:

```text
log h_gdna(length) - log h_rna(length) + (log 0.5 - log p_rna_strand)
```

The strand term is negative, so it favors RNA. The only term large enough to
make gDNA competitive is the FL log ratio.

### EM Effective-Length Normalization

`src/rigel/native/em_solver.cpp`

The EM applies component effective-length normalization by subtracting
`log_eff_len` from component weights. This is the right contract: the scorer
emits `log h(length)`, and EM applies `-log L_tilde`.

The diagnostics show effective length does not cause the pathology. It penalizes
gDNA in both sentinel cases.

### Prior Assembly

`src/rigel/calibration/prior.py` and
`src/rigel/calibration/adaptive_prior.py`

The adaptive prior gives very little gDNA mass in the sentinel failures. In the
single-exon example, the t1 locus receives `alpha_gdna_add ~= 0.113` and
`alpha_rna_add ~= 338.41`. In the wide-intron example, `alpha_gdna_add` is only
`6.8e-6`.

The prior is not causing gDNA to win. The prior is overwhelmingly RNA.

## Worked Example 1: Single-Exon, Zero gDNA

Diagnostic script:

`scripts/debug/diagnose_single_exon_gdna_likelihood.py`

Scenario:

- `t1`: single exon, 500-1500, true RNA only.
- helper transcript: spliced, supplies a small number of clean spliced RNA FL
  observations.
- no simulated gDNA.

Default sampled assignment result:

```text
Simulated fragments: 500
t1 expected: 339
t1 observed: 237
gDNA expected: 0
gDNA pipeline: 102
```

Locus prior:

```text
n_em_fragments: 339
alpha_gdna_add: 0.113436
alpha_rna_add: 338.409759
prior_rna_share_final: 0.999665
enable_gdna: 1
gdna_eff_len: 1198.455842
gdna_em_exposure_weight: 1.0
```

FL model summary:

```text
n_global: 500
n_rna: 39
n_gdna: 0
rna_quality: weak
gdna_quality: fallback
global_fl_mean: 199.71
rna_fl_mean: 213.29
gdna_fl_mean: 199.71
```

Fragment-level score decomposition for the t1 locus:

```text
Raw gDNA score > best RNA score: 248 / 339 units
Median raw gDNA-minus-RNA score: +1.7436
Median FL delta alone: +2.4367
Strand delta: -0.6931
Overhang delta: 0.0
Effective-length delta: -0.4197
Final theta delta after EM: -1.7540
Mean gDNA posterior: 0.2953
Fractional gDNA posterior mass: 100.11
MAP winner gDNA count: 0
```

Interpretation:

- Raw likelihood strongly prefers gDNA for most fragments because the fallback
  global FL surface assigns higher probability at common insert lengths.
- Strand and effective length both push against gDNA.
- After EM, gDNA is not the MAP winner for any fragment. But the posterior is
  broad enough that sampled assignment yields about 102 gDNA calls.

Counterfactual:

```text
If RNA also uses the global FL surface:
sum p_gdna: 0.0187 / 339
mean p_gdna: 5.5e-5
```

This isolates FL asymmetry as the root cause.

Boundary/prior audit:

For the t1 exon region 500-1500:

```text
contained_unspliced: 339
boundary_left_unspliced: 0
boundary_right_unspliced: 0
boundary_left_spliced: 0
boundary_right_spliced: 0
prior_gdna_mean: 0
prior_rna_mean: 339
p_background: 0.000335
p_gdna_only_capture: ~3.4e-13
p_expressed_offtarget: 0.999665
```

There is no boundary flux at the exon boundaries. The tiny gDNA prior comes from
the small background-state tail in `p_states`, not from boundary evidence.

## Worked Example 2: Wide-Intron, Zero gDNA

Diagnostic script:

`scripts/debug/diagnose_wide_intron_gdna_likelihood.py`

Scenario:

- `t1`: two exons, 1000-2000 and 3000-4000.
- no simulated gDNA.

Default sampled assignment result:

```text
Simulated fragments: 500
t1 expected: 500
t1 observed: 453
gDNA expected: 0
gDNA pipeline: 47
```

Locus prior:

```text
n_em_fragments: 453
alpha_gdna_add: 6.838e-6
alpha_rna_add: 448.999892
prior_rna_share_final: ~1.0
gdna_eff_len: 3198.700607
gdna_em_exposure_weight: 1.0
```

FL model summary:

```text
n_global: 500
n_rna: 47
n_gdna: 0
rna_quality: weak
gdna_quality: fallback
global_fl_mean: 199.71
rna_fl_mean: 214.41
gdna_fl_mean: 199.71
```

Fragment-level score decomposition:

```text
Raw gDNA score > best RNA score: 291 / 449 finite-gDNA t1 units
Median raw gDNA-minus-RNA score: +1.8583
Median FL delta alone: +2.5514
Strand delta: -0.6931
Overhang delta: 0.0
Effective-length delta: -0.5825
Final theta delta after EM: -2.9199
Mean gDNA posterior: 0.1081
Fractional gDNA posterior mass: 48.97
MAP winner gDNA count: 0
```

Counterfactual:

```text
If RNA also uses the global FL surface:
sum p_gdna: 9.0e-7 / 453
mean p_gdna: 2.0e-9
```

Again, the failure disappears when the FL likelihood ratio is made symmetric.

Boundary/prior audit:

The t1 exon regions have contained unspliced RNA mass and zero gDNA prior mass.
Boundary channels at the exon/intron interfaces carry spliced evidence, but the
strand-calibrated gDNA boundary means are zero. The adaptive prior remains RNA.

## What Is Actually Failing?

The failure is a misuse of fallback distributions in a class likelihood ratio.

A fallback FL model is useful as a numerical prior predictive distribution. It
prevents undefined likelihoods when a pool has no direct observations. But it is
not class-specific evidence. If gDNA has no FL observations and falls back to the
global library model, then the gDNA-vs-RNA FL Bayes factor should not be allowed
to favor gDNA merely because the RNA model has a small, noisy clean-spliced
sample.

The current system treats these as comparable class-specific likelihoods:

```text
log h_gdna_fallback_global(length) - log h_rna_weak_spliced_only(length)
```

That ratio is not a reliable gDNA-vs-RNA biological signal. It is mostly a
model-support artifact.

This also explains why increasing `rna_call_bias` did not fix the issue. The
problem is upstream of the v6 split dial: the native likelihood surface already
creates a broad gDNA posterior for many fragments.

## Non-Solutions

### High-level adaptive gDNA gating

A gate that disables gDNA whenever the adaptive prior is mostly RNA would make
these tests pass, but it would conflate two different questions:

- Is gDNA structurally possible for this fragment/locus?
- Does the data support gDNA strongly enough to receive posterior mass?

The first is a candidate-construction question. The second is a likelihood and
prior question. A high-level gate hides the FL likelihood pathology instead of
fixing it.

There is also a practical concern: the native EM currently accepts
`locus_enable_gdna`, but the diagnostic code review found no use of that array
inside the native extraction path. That should be fixed separately if explicit
gDNA disablement remains part of the ABI. It is not the root cause of the FL
asymmetry.

### Switching sampled assignment to MAP

MAP assignment would remove many false sampled gDNA calls in these examples,
because gDNA is not the final MAP winner for any t1 fragment. But the posterior
would still contain large false gDNA mass. Changing the assignment policy would
hide the ambiguity in output without fixing the model.

### Tuning test tolerances or adding scenario-specific exceptions

The sentinel tests are exposing a real statistical problem. The solution should
make the FL likelihood ratio well-defined and uncertainty-aware, not silence the
symptom.

## Principled Solution Space

### Option 1: Differential FL Reliability Weighting

This is the most direct and elegant fix.

Keep the existing raw FL models for diagnostics:

- `global_`
- raw/empirical `rna`
- raw/empirical `gdna`
- quality labels and support counts

But introduce scoring FL models that separate the shared library FL shape from
class-specific FL evidence.

Conceptually:

```text
log h_rna_score(l)  = log h_global(l) + w_rg * [log h_rna_raw(l)  - log h_global(l)]
log h_gdna_score(l) = log h_global(l) + w_rg * [log h_gdna_raw(l) - log h_global(l)]
```

where `w_rg` is the reliability of the RNA-vs-gDNA FL contrast. A natural first
definition is:

```text
w_rg = min(reliability_rna, reliability_gdna)
```

with:

```text
reliability_pool = 0                    if pool quality is fallback
reliability_pool = increasing in n_pool if pool quality is weak
reliability_pool = 1                    if pool quality is good
```

The exact reliability curve should be derived from the same empirical-Bayes
uncertainty used to build the FL models. For example, a posterior-predictive
interpretation could use:

```text
reliability_pool = n_pool / (n_pool + prior_ess_pool)
```

or a related monotone function based on posterior variance. The important
property is not the exact formula; it is the invariant:

> A class-specific FL ratio is allowed to influence RNA-vs-gDNA only to the
> extent that both sides of the ratio have class-specific evidence.

In the sentinel cases, `gdna_quality == fallback`, so `reliability_gdna == 0` and
`w_rg == 0`. Both RNA and gDNA score with the shared global FL surface. The FL
ratio becomes zero, and the gDNA posterior collapses as observed in the
counterfactual diagnostics.

When both RNA and gDNA have strong direct support, `w_rg` approaches 1 and the
existing FL discrimination is preserved.

Implementation sketch:

1. Extend `FLModels` with explicit scoring models, for example:
   - `rna_scoring`
   - `gdna_scoring`
   - `fl_contrast_weight`
2. Build these in `src/rigel/calibration/fl.py` after raw `rna` and `gdna` are
   finalized.
3. Normalize any tempered log-probability vector before wrapping it in a
   `FragmentLengthModel`, so `log h(l)` and effective lengths are computed from
   the same distribution.
4. Pass `fl_models.rna_scoring` and `fl_models.gdna_scoring` into
   `FragmentScorer.from_models()`.
5. Keep raw `rna` and `gdna` in summaries, so diagnostics still expose the true
   source support and quality.

Why this is principled:

- It preserves the library-wide FL shape as a shared nuisance distribution.
- It only uses class-specific deviations from global when those deviations are
  supported on both sides of the class comparison.
- It avoids special-casing single-exon loci, zero-gDNA samples, or tests.
- It remains useful when true gDNA has a distinct FL distribution.

### Option 2: Stronger Posterior-Predictive Shrinkage for Weak Pools

The current `_adaptive_prior_ess()` uses `log1p(pool_weight)`, which allows small
weak pools to diverge substantially from global. That was helpful for low-support
FL stress cases where true RNA and gDNA FL distributions were intentionally
different. But in locus EM, a sparse spliced-only RNA model is being used as a
class likelihood against fallback gDNA.

A more Bayesian posterior-predictive policy would keep weak pools closer to
global until their class-specific evidence is strong enough to justify a sharp
likelihood ratio.

Possible forms:

```text
prior_ess_pool = fixed internal ESS, capped below only for very large n
prior_ess_pool = function of posterior uncertainty, not log1p(n)
prior_ess_pool = selected by predictive calibration on held-out spliced lengths
```

This option is simpler than adding explicit scoring models, but less precise. It
reduces the artifact by dampening the weak RNA model, but it does not explicitly
encode the core invariant that fallback gDNA should not supply class-specific FL
evidence.

It is best viewed as a supporting change or a conservative version of Option 1.

### Option 3: Posterior-Weighted RNA/gDNA FL Re-estimation

The most comprehensive long-term solution is to estimate RNA and gDNA FL models
from posterior responsibilities rather than from hard source buckets.

Instead of saying:

```text
RNA FL = annotated-spliced only
gDNA FL = intergenic + intronic only
```

the calibration loop would estimate:

```text
RNA FL += fragment_length * P(fragment is RNA)
gDNA FL += fragment_length * P(fragment is gDNA)
```

using strand deconvolution, regional state probabilities, density evidence, and
possibly an outer EM over FL and region states.

Benefits:

- Single-exon RNA fragments can train the RNA FL model when calibration strongly
  identifies them as RNA.
- Low-gDNA samples do not give gDNA a fake global advantage.
- Real gDNA contamination can still learn a distinct gDNA FL distribution.

Costs:

- This likely requires additional sufficient statistics, perhaps per-region or
  per-fragment FL histograms, not just the current coarse global FL pools.
- Care is needed to avoid circularity: FL should not become self-confirming in
  ambiguous regions.
- This is a larger calibration redesign than needed to fix the immediate
  sentinel failures.

This is the best long-term model if Rigel needs highly accurate FL separation in
messy real samples, but it is heavier than the immediate reliability fix.

### Option 4: Use a Common FL Surface for Fallback Comparisons

A minimal principled subset of Option 1 is:

> If either RNA or gDNA is fallback for a pairwise comparison, use `global_` for
> both sides of that comparison.

For the current sentinel failures, this exactly matches the counterfactual that
collapsed false gDNA. It is easy to reason about and easy to test.

However, it is less smooth than reliability weighting. It creates hard quality
boundaries and does not handle weak-vs-weak or weak-vs-good comparisons as
gracefully.

This could be an acceptable first implementation if it is framed as enforcing a
clear invariant:

> Fallback FL is neutral, never discriminative.

But the reliability-weighted version is preferable because it generalizes.

## Recommended Path

Implement Option 1: differential FL reliability weighting, with Option 4 as the
initial invariant and test oracle.

Recommended design constraints:

1. Keep raw FL models and support counts for diagnostics.
2. Add explicit scoring FL models; do not silently mutate raw diagnostic models.
3. Ensure scoring log-probabilities and EM effective lengths are computed from
   the same scoring distributions.
4. Define a pairwise RNA/gDNA contrast weight that is zero if either side is
   fallback.
5. Choose the weak-pool reliability from posterior-predictive uncertainty, not a
   test-specific threshold.
6. Add tests that assert likelihood-ratio invariants rather than scenario-only
   output numbers.

Concrete scoring invariant:

```text
If gdna_quality == fallback or rna_quality == fallback:
    log h_gdna_score(l) - log h_rna_score(l) == 0
    for all l where both components use the same physical length.
```

For single-exon contained fragments, this means FL is neutral and the remaining
terms favor RNA:

```text
raw gDNA - RNA = log 0.5 - log p_rna_strand <= 0
```

For high-support divergent RNA/gDNA FL data, the contrast weight approaches 1,
so the existing biological FL discrimination remains available.

## Implementation Status

Implemented on 2026-05-26.

The implementation keeps the existing diagnostic FL models intact and adds
separate scoring surfaces:

- `FLModels.rna` and `FLModels.gdna` remain the finalized diagnostic/calibration
   pool models.
- `FLModels.rna_scoring` and `FLModels.gdna_scoring` are the EM scoring models.
- `FLModels.rna_fl_reliability`, `gdna_fl_reliability`, and
   `fl_contrast_weight` expose the joint contrast diagnostics.
- `POOL_SCORING_PRIOR_ESS = 200.0` defines the default score-side prior
   effective sample size for the reliability calculation. It is exposed as
   `CalibrationConfig.fl_scoring_prior_ess`, CLI flag
   `--cal-fl-scoring-prior-ess`, and YAML key `cal_fl_scoring_prior_ess`.

The scoring reliability is:

```text
R_pool = 0                                  if quality == fallback or n_pool == 0
R_pool = n_pool / (n_pool + kappa_score)    otherwise
w_rg   = min(R_rna, R_gdna)
```

The scoring PMFs are then built as:

```text
P_score(l | RNA)  = w_rg * P_empirical(l | RNA)  + (1 - w_rg) * P_global(l)
P_score(l | gDNA) = w_rg * P_empirical(l | gDNA) + (1 - w_rg) * P_global(l)
```

This makes fallback FL neutral by construction. If either class lacks enough FL
evidence to leave fallback, `w_rg = 0`, both scoring PMFs are identical to the
global PMF, and FL contributes no RNA-vs-gDNA likelihood ratio.

Effective lengths for RNA quantification now use `rna_scoring`, so the scoring
`log h(l)` and EM effective-length normalization remain based on the same FL
distribution. Calibration density and prior-exposure calculations still use the
raw/finalized `gdna` model, not `gdna_scoring`, because those steps estimate
gDNA opportunity rather than the RNA-vs-gDNA fragment likelihood contrast.

The implementation deliberately leaves `locus_enable_gdna` untouched. That ABI
cleanup is separate technical debt and should be handled in a focused native
refactor.

## Test Plan

### Unit Tests

Add tests around `build_fl_models()` or a new FL scoring-model helper:

1. `gdna` fallback, `rna` weak:
   - raw `rna` and raw `gdna` summaries remain unchanged.
   - scoring FL ratio is neutral when gDNA has no support.
2. both pools fallback:
   - scoring RNA and scoring gDNA are both global.
3. both pools good and divergent:
   - scoring distributions retain a strong class-specific FL ratio.
4. one weak, one good:
   - contrast is partially damped, not hard disabled unless fallback.
5. effective lengths are recomputed from scoring distributions.

### Scenario Tests

Use the two sentinel scenarios:

- `tests/scenarios/test_single_exon.py` zero-gDNA variants.
- `tests/scenarios/test_wide_intron.py::TestWideIntronInsertPenalty::test_baseline_no_false_nrna`.

Expected behavior:

- zero-gDNA false gDNA posterior mass collapses.
- sampled assignment no longer produces large false gDNA counts.
- RNA transcript accuracy passes without high-level gDNA gates.

### Regression/Stress Tests

Re-run FL stress diagnostics, especially cases where true RNA and true gDNA FL
distributions differ. The fix must not erase valid FL discrimination when both
class-specific pools are supported.

Also run golden tests only after confirming the statistical behavior is correct.

## Open Design Questions

1. What exact reliability curve should map weak-pool support to contrast weight?
   A posterior-predictive uncertainty derivation is preferable to a hand-tuned
   threshold.
2. Should reliability be a single scalar per pool or length-bin dependent?
   A scalar is simpler and likely sufficient first. Bin-dependent uncertainty is
   more exact but adds complexity.
3. Should `FLModels` expose raw and scoring models separately, or should the
   scorer accept a separate `FLScoringSurface` object? A separate object may keep
   summaries cleaner.
4. Should the native `locus_enable_gdna` ABI be fixed at the same time? It is not
   the FL root cause, but leaving an unused enable array is misleading.

## Bottom Line

The failure is not that Rigel lacks a strong enough RNA prior. The prior is
already overwhelmingly RNA. The failure is that the likelihood comparison is not
fair: fallback/global gDNA FL is being compared against sparse spliced-only RNA
FL as if both were equally reliable class-specific evidence.

The clean fix is to make FL class discrimination uncertainty-aware. Shared global
FL should remain available as a nuisance distribution, but class-specific FL
deviations should only affect RNA-vs-gDNA posterior odds when both classes have
enough evidence to support that contrast.