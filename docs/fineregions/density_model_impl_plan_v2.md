# gDNA Density Model Implementation Plan v2

Date: 2026-05-24
Status: implementation plan
Companion: `calibration_roadmap_v3.md`
Supersedes/refines: `density_model_impl_plan_v1.md`

## 0. Executive Summary

v1 made the right architectural move: use contained DNA-dominant anchors to fit
the background density prior, and use boundary flux as local Bayesian evidence
for region-specific enrichment. v2 keeps that model but tightens two important
implementation details before code is written.

First, background density must be explicitly overdispersed. Real intergenic and
intronic anchors are not Poisson-clean. Mappability, GC, copy-number, repetitive
sequence, alignment artifacts, and hybrid-capture leakage all make regional
counts much more variable than `Poisson(rho * L)`. The density prior must learn
that extra variation or be capped so it cannot become too stubborn.

Second, low-opportunity regions need diagnostics. Very short regions may have
near-zero contained opportunity and only small boundary-side opportunity. This
is not a fatal flaw because short exons cannot support many fully contained
gDNA fragments anyway, and crossing fragments should create boundary mass. But
the model must report when local boundary evidence is too weak to override the
background prior.

The v2 implementation should therefore use:

```text
contained anchors              -> exposure-adjusted method-of-moments Gamma prior
prior opportunity cap          -> prevents overconfident/stubborn background priors
boundary counts/opportunity    -> local Gamma posterior update
contained NB predictive        -> mean, variance, upper count, tail/tension
low-opportunity diagnostics    -> flags for short/weakly informed regions
```

This is simpler than a full numerical Negative Binomial MLE, robust enough for
v1 code, and aligned with the roadmap's requirement that density and strand stay
independent until probabilistic fusion.

## 1. Answers to the Outside Critique

### 1.1 Do we need to model overdispersion?

Yes. We should model overdispersion explicitly.

The nuance is that v1's written likelihood was not a pure Poisson model if
implemented literally. It was a conditional Poisson with a Gamma-distributed
regional rate, and the marginal distribution is Negative Binomial. That already
allows overdispersion in principle.

The practical trap remains real: an unconstrained fit can still choose a very
large `alpha` and `beta`, especially in sparse or nearly Poisson-looking anchor
sets. Since local boundary evidence updates the prior by adding `B_tot` to
`alpha` and `L_b` to `beta`, a huge `beta` makes the prior immovable. That is
bad for hybrid capture, where a captured boundary must be allowed to override
the off-target background.

v2 makes two changes:

```text
1. Estimate overdispersion with a fast exposure-adjusted method of moments.
2. Cap the fitted prior opportunity beta while preserving the fitted mean.
```

The cap is not a biological claim. It is an engineering guardrail: local data
must have a finite path to update the prior.

### 1.2 MoM or numerical NB MLE?

Use method of moments first.

Reasons:

- It is fast, vectorized, deterministic, and easy to test.
- It works with fractional counts without pretending they are exact integers.
- It gives direct diagnostics: observed variance, Poisson variance, extra
  variance, overdispersion factor, and prior strength.
- It avoids optimizer failure modes at the Poisson boundary.
- It is good enough for the first independent density surface because the main
  downstream requirement is a calibrated uncertainty scale, not perfect
  asymptotic efficiency.

Numerical NB MLE is a reasonable later refinement, but only after MoM outputs
are validated against no-capture and hybrid-capture simulations. If MLE is added
later, it should use MoM as its starting point and keep the same beta cap.

### 1.3 Is the short-exon issue critical?

It is mostly a diagnostic issue for v1, not a blocker.

Contained opportunity truly collapses for regions shorter than the gDNA
fragment-length distribution, which is correct: those regions cannot contain
complete unspliced gDNA fragments. Boundary opportunity does not vanish as
aggressively. For a side of length `S`, the current helper computes:

```text
L_side(S) = sum_ell h(ell) * min((ell - 1) / 2, S / 2)
```

For a very short region and ordinary gDNA fragment lengths, each side is about
`S / 2`, so both sides together are about `S`. That is small, but it is not
mathematically zero.

The user's interpretation is right: tiny exons should not support contained
unspliced fragments. Unspliced fragments that overlap them should usually cross
into adjacent sequence and enter boundary channels. Spliced boundary flux is
also useful evidence of expression/capture context, but it should not be used
as gDNA count evidence in v2.

Still, v2 should add explicit low-opportunity diagnostics:

```text
FLAG_LOW_CONTAINED_OPPORTUNITY
FLAG_LOW_BOUNDARY_OPPORTUNITY
FLAG_PRIOR_DOMINATED_POSTERIOR
```

These flags let Phase 8 spatial borrowing or capture-exposure learning identify
short regions that need neighboring context.

## 2. Simplified v2 Statistical Contract

### 2.1 Observations

For each region `r`:

```text
C_r   = contained_unspliced_pos + contained_unspliced_neg
B_l   = boundary_left_unspliced_pos + boundary_left_unspliced_neg
B_r   = boundary_right_unspliced_pos + boundary_right_unspliced_neg
B_tot = B_l + B_r

L_c   = contained effective length under gDNA FL
L_l   = left boundary-side opportunity under gDNA FL
L_r   = right boundary-side opportunity under gDNA FL
L_b   = L_l + L_r
```

Only unspliced channels enter the density likelihood. Spliced channels may be
recorded as diagnostics but not as gDNA numerator mass.

### 2.2 Anchor classes

Fit background priors from contained anchors:

```text
INTERGENIC: signature == 0x0
INTRON:     intron bits present and no exon bits (0x4, 0x8, 0xC)
```

Eligibility:

```text
L_c >= density_min_effective_length
C_r is finite and nonnegative
```

Exon and mixed exon/intron rows receive predictions but do not train the
background prior in v2.

### 2.3 Primary estimator: exposure-adjusted MoM

For an anchor family with counts `C_i` and opportunities `L_i`, estimate the
mean density:

```text
rho = sum_i C_i / sum_i L_i
```

Under the Gamma-Poisson model:

```text
lambda_i ~ Gamma(alpha, beta)
C_i | lambda_i ~ Poisson(lambda_i * L_i)

E[C_i]   = rho * L_i
Var[C_i] = rho * L_i + phi * (rho * L_i)^2
phi      = 1 / alpha
beta     = alpha / rho
```

Estimate `phi` in count space:

```text
mu_i = rho * L_i
S    = sum_i (C_i - mu_i)^2
B    = sum_i mu_i
A    = sum_i mu_i^2

phi_hat = max(0, (S - B) / max(A, eps))
alpha_raw = 1 / max(phi_hat, phi_floor)
beta_raw  = alpha_raw / max(rho, eps)
```

This is the same spirit as the existing beta-binomial strand overdispersion
MoM: compare observed residual variance to the variance expected under the
non-overdispersed model, then attribute the excess to a concentration parameter.

### 2.4 Prior opportunity cap

After MoM, cap the prior opportunity:

```text
beta = min(beta_raw, beta_cap)
alpha = rho * beta
```

This preserves the fitted mean density while bounding how much local boundary
evidence is needed to update the posterior.

The first implementation should compute `beta_cap` from the data rather than
exposing a new CLI knob:

```text
positive_boundary = L_b over regions with L_b > 0
typical_boundary_opportunity = median(positive_boundary)
beta_cap = density_prior_equiv_boundaries * typical_boundary_opportunity
```

Initial internal default:

```text
density_prior_equiv_boundaries = 20
```

This means the background prior can be roughly as informative as 20 typical
region-boundary observations, but not thousands. The exact value should be
validated and can become a config field later if needed.

### 2.5 Fallbacks

Use deterministic fallback depth when an anchor family is sparse:

```text
0: class genome-wide prior
1: all DNA-dominant anchors genome-wide
2: broad weak fallback from total observed density
```

Reference-specific pooling is deferred. The dataclasses should still carry
`family` and `fallback_depth` arrays so reference-aware pooling can land later
without reshaping downstream consumers.

### 2.6 Local boundary update

For every region, choose the best available prior and update with local boundary
evidence:

```text
alpha_post = alpha_prior + B_tot
beta_post  = beta_prior  + L_b
rho_post   = alpha_post / beta_post
```

Diagnostic information weight:

```text
w_boundary_opportunity = L_b / (beta_prior + L_b)
w_boundary_count       = B_tot / (alpha_prior + B_tot)  if B_tot > 0 else 0
```

These weights are diagnostics, not fusion weights. They tell us whether the
posterior was actually informed by local boundary evidence.

### 2.7 Contained predictive distribution

Predict held-out contained background with the Negative Binomial predictive:

```text
N_c^gdna | boundary ~ NB(
    shape = alpha_post,
    p = beta_post / (beta_post + L_c),
)
```

Mean and variance:

```text
mean_c = rho_post * L_c
var_c  = mean_c * (1 + L_c / beta_post)
```

The region-level unbounded density evidence is:

```text
mean_unbounded  = B_tot + mean_c
upper_unbounded = B_tot + NB_quantile(confidence, alpha_post, beta_post, L_c)
var_unbounded   = var_c
```

Do not clip these values to the observed count. Fusion will bound them later.

### 2.8 Tail/tension diagnostics

Record at least:

```text
observed_compatible_count = B_tot + C_r
tail_probability = P(N_c^gdna > C_r | boundary)
expected_tail_count = E[(N_c^gdna - C_r)+ | boundary]
density_over_observed_ratio = mean_unbounded / max(observed_compatible_count, eps)
```

If exact `expected_tail_count` is expensive, v2 may start with a stable
truncated-sum implementation for small means and a normal approximation for
large means, but the field should exist from the start.

## 3. Short-Region Diagnostics

The density model should not treat short regions as errors. It should mark the
quality of the density prediction.

Add flags to `DensityEvidence`:

```text
FLAG_LOW_CONTAINED_OPPORTUNITY: L_c < density_min_effective_length
FLAG_LOW_BOUNDARY_OPPORTUNITY:  L_b < density_min_boundary_opportunity
FLAG_PRIOR_DOMINATED:           w_boundary_opportunity < min_boundary_info
FLAG_HIGH_TAIL_TENSION:         tail_probability > tail_probability_warn
```

Suggested internal defaults:

```text
density_min_effective_length = 1.0
density_min_boundary_opportunity = 1.0
min_boundary_info = 0.05
tail_probability_warn = 0.5
```

These flags are for summaries, fusion quality, and later exposure learning.
They should not disable density evidence by themselves.

Spliced contained and spliced boundary channels remain diagnostic. In v2 they
can help explain short targeted exons, but they should not increase gDNA count
or boundary likelihood.

## 4. Dataclass Adjustments

### 4.1 `DensityObservation`

Keep the v1 observation schema, with one addition for diagnostics:

```python
@dataclass(frozen=True, slots=True)
class DensityObservation:
    contained_count: np.ndarray
    boundary_left_count: np.ndarray
    boundary_right_count: np.ndarray
    boundary_count: np.ndarray
    observed_compatible_count: np.ndarray

    contained_leff: np.ndarray
    boundary_left_leff: np.ndarray
    boundary_right_leff: np.ndarray
    boundary_leff: np.ndarray
    total_leff: np.ndarray

    spliced_count: np.ndarray
    class_code: np.ndarray
    anchor_intergenic: np.ndarray
    anchor_intron: np.ndarray
```

### 4.2 `GammaRatePrior`

Make overdispersion and cap diagnostics first-class:

```python
@dataclass(frozen=True, slots=True)
class GammaRatePrior:
    family: str
    alpha: float
    beta: float
    mean_density: float
    phi: float
    beta_raw: float
    beta_cap: float
    cap_applied: bool
    n_regions: int
    n_fragments: float
    eff_length: float
    residual_sum: float
    poisson_variance_sum: float
    extra_variance_basis_sum: float
    fallback_depth: int
    fit_status: str
```

### 4.3 `DensityEvidence`

Add local-information and flags fields:

```python
@dataclass(frozen=True, slots=True)
class DensityEvidence:
    mean_unbounded: np.ndarray
    upper_unbounded: np.ndarray
    variance_unbounded: np.ndarray
    precision_unbounded: np.ndarray

    relative_exposure: np.ndarray
    rho_post: np.ndarray
    alpha_post: np.ndarray
    beta_post: np.ndarray
    w_boundary_opportunity: np.ndarray
    w_boundary_count: np.ndarray

    tail_probability: np.ndarray
    expected_tail_count: np.ndarray
    density_over_observed_ratio: np.ndarray

    family: np.ndarray
    fallback_depth: np.ndarray
    quality: np.ndarray
    flags: np.ndarray
    confidence: float
```

`precision_unbounded` should be a monotone quality proxy derived from posterior
variance and local information, not a strand or fusion weight. A simple initial
definition is acceptable:

```text
cv = sqrt(var_unbounded) / max(mean_unbounded, eps)
precision_unbounded = 1 / (1 + cv)
```

## 5. Implementation Scope Simplifications

v2 intentionally removes or defers several things from v1:

1. No reference-specific priors in the first implementation.
2. No hard capture-mode detector.
3. No numerical NB MLE in the first implementation.
4. No use of strand output inside density fitting.
5. No EM prior assembly in the density module.
6. No spliced-count contribution to gDNA likelihood.

The first density implementation should be three pieces:

```text
calibration/evidence.py             -> memory-light count ledger
calibration/density_observation.py  -> counts and opportunities
calibration/density_model.py        -> MoM priors, local update, predictions
```

If this feels like too many modules during implementation, `DensityObservation`
can live in `density_model.py` temporarily, but the API boundary should remain
clear: observations are data, evidence is model output.

## 6. Calibration Integration

The orchestrator should build density evidence after FL models and before
returning `CalibrationResult`:

```text
build FL models
build RegionArrays and PayloadArrays
build RegionCountLedger
compute contained global densities for backward-compatible summaries
build strand-only RegionGdnaEstimate
build DensityObservation
fit DensityEvidence with MoM priors and local boundary updates
return CalibrationResult with strand and density surfaces
```

`DensityEvidence` should be added to `CalibrationResult` alongside the current
strand-only `region_gdna`. Do not rename `region_gdna` yet; probabilistic
density-strand fusion will introduce a fused result later.

## 7. Config Changes

Add to `CalibrationConfig`:

```python
gdna_density_confidence: float = 0.95
density_min_effective_length: float = 1.0
```

Internal constants in `density_model.py`:

```python
DENSITY_PRIOR_EQUIV_BOUNDARIES = 20.0
DENSITY_PHI_FLOOR = 1.0e-6
DENSITY_MIN_BOUNDARY_OPPORTUNITY = 1.0
DENSITY_MIN_BOUNDARY_INFO = 0.05
DENSITY_TAIL_PROBABILITY_WARN = 0.5
```

Expose `--gdna-density-confidence` in the CLI. Keep the other constants
internal until validation shows they need to be user-tunable.

## 8. Tests

### 8.1 MoM overdispersion tests

- Simulated Gamma-Poisson anchors recover mean density and overdispersion within
  a loose tolerance.
- Pure Poisson-like anchors do not create an immovable prior because beta is
  capped.
- Larger residual variance yields smaller `alpha` and broader priors.
- Fractional counts produce finite parameters.
- Sparse anchors use deterministic fallback depth.

### 8.2 Prior cap tests

- `beta_raw > beta_cap` preserves `alpha / beta == rho` after capping.
- A high-boundary region can move `rho_post` above background when the cap is
  applied.
- Cap diagnostics are reported in the density summary.

### 8.3 Short-region tests

- A region shorter than the gDNA FL has near-zero contained opportunity.
- Boundary opportunity for a short region follows the current side formula.
- Short regions get low-opportunity flags but still receive finite predictions.
- Spliced counts do not change density evidence.

### 8.4 End-to-end calibration tests

- `calibrate()` returns `density_observation` and `density_evidence`.
- Density evidence is unchanged by strand-specificity changes when the payload
  is identical.
- Summary reports prior overdispersion, cap usage, low-opportunity flag counts,
  and tail/tension metrics.

## 9. Decision Summary

For the first implementation:

```text
Use MoM, not numerical NB MLE.
Model overdispersion explicitly.
Cap prior opportunity to keep boundary evidence useful.
Treat short-exon vulnerability as a flagged low-information case, not a blocker.
Keep density independent of strand until fusion.
Keep density independent of EM prior allocation until footprint-aware priors land.
```
