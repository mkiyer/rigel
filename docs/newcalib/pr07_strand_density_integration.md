# PR 07 - Strand and Density Evidence Integration

## Status

Draft implementation plan.

This PR supersedes the tiered regional gDNA mass handoff from PR 03. The old production logic treated
strand as Tier 1, boundary/density as Tier 2, and background as Tier 3. That framing is conceptually
wrong for the final model. Strand and density are not ordered fallbacks. They are two evidence
channels about the same latent regional gDNA mass.

PR 07 promotes the fused evidence model to the production calibration path:

```text
posterior(D_r) proportional to density_evidence(D_r) * strand_evidence(D_r) * constraints(D_r)
```

where `D_r` is the gDNA portion of the unspliced-compatible mass in region `r`.

The central rule is simple: evidence should influence the answer in proportion to its uncertainty.
Strong strand signal should matter. Weak or unstranded signal should become nearly flat and matter
little. Strong boundary/density signal should matter regardless of strand. Sparse density evidence
should be broad and matter gently. No estimator should own the answer just because it appears first
in an `if` chain.

## Trigger

PR 06 validation found a no-gDNA, low-strand-specificity, nRNA-heavy case where calibration inferred
about 19% gDNA mass before EM. The immediate root cause was the PR 03 production bridge:

```text
if strand reliability > epsilon:
    use strand slab as regional gDNA mass
elif boundary evidence exists:
    use boundary sweep
else:
    use background fallback
```

The local strand deconvolution produced noisy positive mixed-slab estimates. Because the production
handoff let strand own the mass whenever any reliability was present, those noisy estimates became
`RegionUnsplicedMass.gdna_mass`, even though density evidence for the same case said there was no
detectable gDNA.

The PR 06 tactical fix forces zero gDNA mass when density returns `rho_ref_source == "ZERO"`. That
fix is acceptable as a sentinel guard, but it is still gate-like. PR 07 replaces that bridge with a
smooth uncertainty-aware fusion model.

## Goals

1. Make strand evidence expose uncertainty explicitly.
2. Make density evidence expose uncertainty explicitly.
3. Fuse strand and density synergistically as evidence about one latent `D_r`.
4. Handle unstranded and weakly stranded libraries by making strand evidence broad or flat, not by
   pretending a confident strand call exists.
5. Handle transcript-strand-ambiguous regions by omitting strand evidence structurally and using the
   density channel.
6. Delete the production Tier 1 -> Tier 2 -> Tier 3 ownership model.
7. Remove hard behavior cliffs where a tiny change in estimated gDNA or strand contrast changes
   which model owns the answer.
8. Preserve mass conservation:

```text
0 <= D_r <= T_r
R_r = T_r - D_r
D_r + R_r = T_r
```

## Non-Goals

- Do not redesign transcript-level EM scoring in this PR.
- Do not introduce gene-level logic into calibration, priors, locus construction, scoring, or EM.
- Do not keep both the tiered production path and the fused production path as competing modes.
- Do not solve every downstream no-gDNA EM assignment issue in this PR. PR 07 fixes the calibration
  evidence handoff and regional prior/exposure inputs.

## Design Principles

### Evidence Sources Are Not Tiers

The old names `METHOD_STRAND`, `METHOD_BOUNDARY`, and `METHOD_BACKGROUND_FALLBACK` are useful only
as historical diagnostics. They should not be production control flow. A region can have both strand
and density information. A region can also have only one of them. The fusion model should handle all
cases through likelihood multiplication and source uncertainty.

### Absence Is Different From Weak Evidence

There are two distinct situations:

1. A source is structurally inapplicable.
2. A source is applicable but weak or uncertain.

For example, a region with overlapping transcripts on opposite strands has no unique transcript
strand. Strand folding is not defined there, so strand evidence is absent. This is not a confidence
threshold; it is a structural missing channel.

By contrast, an unstranded library may still produce a fitted `StrandSummary`, but its contrast is
near zero or its confidence interval includes zero. In that case strand evidence is present as a
calculation but should be effectively flat over `D_r`.

### Zero Should Be a Posterior, Not a Gate

The PR 06 `force_zero_gdna_mass` flag should not survive as the long-term model. A no-gDNA dataset
should produce a density posterior concentrated near zero because the observed counts and
opportunities support that conclusion. A barely nonzero gDNA dataset should produce a low positive
posterior with finite variance. The transition from zero to tiny positive gDNA should be continuous.

An exact point mass at zero should be reserved for mathematical degeneracy, such as no possible gDNA
opportunity or an explicitly impossible support. Ordinary absence of observed gDNA-compatible anchor
fragments should become a low-density posterior with uncertainty determined by opportunity, not an
external override.

### Avoid Double Counting By Defining Conditional Evidence

Density evidence should model count/opportunity information:

```text
P(total compatible mass and boundary geometry | D_r)
```

Strand evidence should model orientation information conditional on the compatible mass:

```text
P(strand-folded orientation counts | D_r, total compatible mass)
```

This lets the same fragments contribute both geometry and orientation without double counting the
same likelihood term. The density channel should not pre-filter boundary counts through strand means
and then multiply by a second strand likelihood for the same signal.

## Current Code Audit

### Useful Pieces To Keep

- `strand_summary.py` already exposes global strand contrast and a confidence margin.
- `strand_deconv.py` already computes folded strand counts, beta-binomial strand likelihood pieces,
  per-region reliability, and exact/approximate paths.
- `density_model.py` already builds Gamma-rate posterior predictive density evidence and has
  `density_logpmf_grid()`.
- `integration.py` already contains a prototype `fuse_density_and_strand()` that multiplies density
  and strand evidence region-by-region.
- `density_observation.py`, `boundary_model.py`, and `boundary_sweep.py` already represent the
  geometry and boundary projection machinery needed by the density channel.
- `RegionUnsplicedMass` already provides the downstream mass-conserving handoff consumed by adaptive
  priors and exposure learning.

### Production Problems To Remove

- `_orchestrator.py` computes `_strand_summary_identifiable()` and hard-drops strand channels when
  the global summary is weak.
- `build_region_unspliced_mass()` implements Tier 1 -> Tier 2 -> Tier 3 ownership.
- `build_region_unspliced_mass()` marks strand as owning a region when reliability is only greater
  than `1e-6`, which creates cliffs.
- Near-unstranded strand deconvolution sets `mean_count = n_total` and `upper_count = n_total` while
  also setting precision to zero. That is dangerous because callers can mistake the mean for a real
  gDNA estimate.
- `build_boundary_local_posterior()` changes boundary alpha inputs depending on whether strand
  channels are supplied. PR 07 should separate geometry/count density evidence from orientation
  evidence instead of replacing density counts with strand-deconvolved counts.
- `calibration_m_step()` currently refits background from strand-contained means when strand exists.
  It should refit from fused regional gDNA mass and fused uncertainty.
- `estimate_region_exposure()` contains a temporary high-`p_unexpressed` pool gate from PR 06. PR 07
  should replace that with continuous weighting from expression probability and fused information.

## Target Data Contracts

### Strand Evidence Contract

Introduce a production evidence object that represents strand information as a likelihood or an
uncertain posterior summary, not as an owned mass estimate.

Recommended name:

```python
@dataclass(frozen=True, slots=True)
class StrandGdnaEvidence:
    n_total: np.ndarray                 # float64[R]
    applicable: np.ndarray              # bool[R]
    structural_absent: np.ndarray       # bool[R]
    information: np.ndarray             # float64[R], precision-like, >= 0
    mean_count: np.ndarray              # float64[R], diagnostic only when informative
    variance_count: np.ndarray          # float64[R], inf or very large when flat
    upper_count: np.ndarray             # float64[R]
    log_bayes_factor: np.ndarray        # float64[R]
    reliability: np.ndarray             # float64[R], diagnostic in [0, 1]
    flags: np.ndarray                   # uint16[R]
```

This object should also provide, directly or through helper functions:

```python
strand_logpmf_grid(evidence, region_idx, d_grid) -> np.ndarray
strand_normal_params(evidence, region_idx, n_observed) -> tuple[mean, variance]
```

#### Strand Applicability

`applicable` is true only when all structural requirements hold:

- the region has a unique transcript-strand class (`TS_POS` or `TS_NEG`),
- the library has orientation information that can in principle distinguish RNA from gDNA,
- the strand-folded counts have nonzero support.

`applicable` is false for:

- `TS_NONE`,
- `TS_AMBIG`, including overlapping transcripts on opposite strands,
- regions where no strand-folded observations are available.

When `applicable == false`, the strand likelihood is flat:

```text
log P_strand(D_r = d) = constant for all d
information = 0
variance = infinity
```

This makes density-only regions the easiest fusion case: one channel is absent, so the posterior is
just the density posterior bounded to `[0, T_r]`.

#### Weak Or Unstranded Libraries

For unstranded or weakly stranded data, do not hard-drop the strand object. Instead, make its
likelihood broad.

The global protocol uncertainty should scale strand information. A simple first implementation:

```text
c = abs(2 * p_r1_sense - 1)
se = signed_strand_contrast_se
global_info_scale = c^2 / (c^2 + (z * se)^2 + eps)
```

Then:

```text
effective_strand_information_r = local_orientation_information_r * global_info_scale
```

If `p_r1_sense` is near 0.5 or the contrast interval includes zero, `global_info_scale` approaches
zero. The strand likelihood becomes almost flat. If the library is strongly strand-specific and the
global estimate is precise, the strand likelihood remains sharp.

This scaling is a first implementation. A better later implementation can marginalize over the
global `p_r1_sense` posterior rather than using a scalar information scale.

#### Strand Output Semantics

Do not set near-unstranded `mean_count = n_total` as if all mass is gDNA. For flat strand evidence,
the mean is not a standalone estimate. Use one of these safer options:

- set `mean_count = nan` and require callers to use likelihood helpers, or
- set `mean_count = 0.5 * n_total` but mark `information = 0` and document it as diagnostic only.

The preferred production path is to avoid using `StrandGdnaEvidence.mean_count` for mass handoff at
all. Fusion should consume the likelihood or variance.

### Density Evidence Contract

Extend or clarify `DensityEvidence` so it exposes uncertainty in the same way:

```python
@dataclass(frozen=True, slots=True)
class DensityGdnaEvidence:
    mean_count: np.ndarray              # float64[R]
    variance_count: np.ndarray          # float64[R]
    upper_count: np.ndarray             # float64[R]
    information: np.ndarray             # float64[R], precision-like, >= 0
    alpha_post: np.ndarray              # float64[R]
    beta_post: np.ndarray               # float64[R]
    contained_leff: np.ndarray          # float64[R]
    boundary_count: np.ndarray          # float64[R]
    applicable: np.ndarray              # bool[R]
    flags: np.ndarray                   # uint16[R]
```

The existing `DensityEvidence` already has most of this:

- `mean_unbounded`,
- `upper_unbounded`,
- `variance_unbounded`,
- `alpha_post`,
- `beta_post`,
- `contained_leff`,
- `boundary_count`,
- `density_logpmf_grid()`.

PR 07 should make these fields part of the production contract and add an explicit information
array:

```text
information_r = 1 / max(variance_count_r, variance_floor)
```

For density-only absence of observed gDNA, use a finite low-density posterior whenever there is
finite opportunity. Avoid point-mass zero except for true degeneracy. This removes the `rho_ref ==
ZERO` behavior cliff while still allowing high-confidence near-zero density when opportunity is
large.

### Fused Evidence Contract

Introduce a fused regional posterior object:

```python
@dataclass(frozen=True, slots=True)
class FusedRegionGdnaEvidence:
    mean_count: np.ndarray              # float64[R]
    variance_count: np.ndarray          # float64[R]
    upper_count: np.ndarray             # float64[R]
    rna_lower_count: np.ndarray         # float64[R]
    p_nonzero: np.ndarray               # float64[R]
    observed_compatible_count: np.ndarray

    density_information: np.ndarray     # float64[R]
    strand_information: np.ndarray      # float64[R]
    density_weight: np.ndarray          # float64[R], diagnostic tau_density / tau_sum
    strand_weight: np.ndarray           # float64[R], diagnostic tau_strand / tau_sum

    density_applicable: np.ndarray      # bool[R]
    strand_applicable: np.ndarray       # bool[R]
    flags: np.ndarray                   # uint16[R]
```

`density_weight` and `strand_weight` are diagnostics, not gates. They explain which source supplied
information in the posterior. They should change smoothly as uncertainty changes.

### Regional Mass Contract

Keep `RegionUnsplicedMass` as the downstream handoff if that is the smallest code change, but change
its semantics:

```text
total_mass = observed compatible unspliced mass T_r
gdna_mass = fused posterior E[D_r]
rna_mass = T_r - E[D_r]
precision = fused information or effective sample size
```

Deprecate algorithmic use of `method`. Options:

1. Replace `method` with `source_mask` and source weights in a new dataclass.
2. Keep `method` temporarily for output compatibility but set it to a neutral `METHOD_FUSED` value.
3. Remove `method` from production and update tests/goldens.

Because the newcalib roadmap allows schema churn, option 1 is preferred:

```python
@dataclass(frozen=True, slots=True)
class RegionUnsplicedMass:
    total_mass: np.ndarray
    gdna_mass: np.ndarray
    rna_mass: np.ndarray
    region_size_bp: np.ndarray
    unspliced_counts: np.ndarray
    precision: np.ndarray
    density_information: np.ndarray
    strand_information: np.ndarray
    density_weight: np.ndarray
    strand_weight: np.ndarray
    flags: np.ndarray
```

If downstream churn is too large for one PR, keep `method` only as deprecated diagnostics and remove
all logic that branches on it.

## Fusion Algorithm

### Exact Fusion For Small Counts

For each region with observed compatible count `N_r`, define an integer grid:

```text
d in {0, 1, ..., round(N_r)}
```

Then compute:

```text
log_post(d) = log_density(d) + log_strand(d) + log_bounds(d)
```

where `log_bounds(d)` is zero inside `[0, N_r]` and `-inf` outside.

Summarize:

```text
mean = E[d]
variance = Var[d]
upper = posterior_quantile(confidence)
p_nonzero = P(d > 0)
rna_lower = max(N_r - upper, 0)
```

Use the exact path up to a configurable threshold, starting with the existing
`MAX_EXACT_POSTERIOR_N` or `ADAPTIVE_EXACT_MAX` values already used by `strand_deconv.py` and
`integration.py`.

### Approximate Fusion For Large Counts

For large counts, represent each channel as a normal or truncated normal over `D_r`:

```text
density: mean mu_d, variance var_d
strand:  mean mu_s, variance var_s
```

Convert variances to precisions:

```text
tau_d = 1 / var_d
tau_s = 1 / var_s
```

Absent or flat evidence has `tau = 0`. Degenerate evidence has very high precision and must be
handled before generic approximations to avoid optimizer fallback inventing mass.

Fuse:

```text
tau = tau_d + tau_s
mu = (tau_d * mu_d + tau_s * mu_s) / tau
var = 1 / tau
```

Then truncate to `[0, N_r]` and compute the requested posterior summaries. If both precisions are
zero, fall back to a broad bounded distribution and flag the region as low information.

### Source Weights

For diagnostics:

```text
density_weight = tau_d / (tau_d + tau_s)
strand_weight = tau_s / (tau_d + tau_s)
```

When only density is applicable, `density_weight = 1` and `strand_weight = 0`. When strand is
unstranded/flat, `tau_s = 0`, so the same result arises naturally. When both are strong, weights
reflect relative information.

### Deterministic And Near-Zero Density

Replace the PR 06 external `force_zero_gdna_mass` with density likelihood semantics:

- true mathematical zero: density likelihood is a point mass at `D_r = 0`, variance zero,
  information infinite or a large sentinel;
- observed-zero with finite opportunity: density likelihood is concentrated near zero but has finite
  variance;
- barely nonzero gDNA: density likelihood has low positive mean and finite variance.

This removes the cliff where strand suddenly starts owning the answer as soon as density is barely
nonzero.

## Calibration Flow After PR 07

Target production flow:

```text
payload + region table
    -> density observation
    -> density evidence with posterior uncertainty
    -> strand counts
    -> strand evidence with posterior uncertainty
    -> fused regional gDNA posterior
    -> RegionUnsplicedMass from fused mean and information
    -> expression-state update
    -> background density refit from fused mass, p_unexpressed, and fused information
    -> exposure refit from fused mass, p_unexpressed, and fused information
    -> adaptive prior and EM exposure handoff
```

Important changes from the current flow:

- `_strand_summary_identifiable()` no longer determines whether strand exists globally.
- `calibration_strand_channels = strand_channels if strand_usable else None` goes away.
- `build_boundary_local_posterior()` no longer substitutes strand-deconvolved boundary means as the
  density channel. Strand orientation remains a separate likelihood.
- `build_region_unspliced_mass()` is replaced or rewritten as `build_fused_region_unspliced_mass()`.
- `calibration_m_step()` refits background from fused regional mass, not from strand-only contained
  means.
- `estimate_background_density()` and `estimate_region_exposure()` use fused information/variance as
  continuous weights.
- The PR 06 `force_zero_gdna_mass` argument is deleted after fused density zero/near-zero behavior is
  in production.

## Ambiguous Opposite-Strand Regions

This is a first-class requirement.

Regions with overlapping transcripts on opposite strands cannot use strand orientation to separate
RNA from gDNA because there is no single region strand. In these regions:

```text
strand_applicable = false
strand_information = 0
strand log likelihood = constant
```

Fusion then reduces to density-only evidence:

```text
posterior(D_r) proportional to density_evidence(D_r)
```

This is not a fallback tier. It is the natural result of one missing likelihood channel.

Acceptance tests must include at least one opposite-strand overlap case where:

- strand evidence is structurally absent,
- density evidence remains active,
- fused output equals density-only output within tolerance,
- no strand flags or pseudo-means can move `D_r`.

## Exposure Refit After Fusion

The PR 06 high-`p_unexpressed` exposure gate should become a continuous weight.

Current temporary behavior:

```text
tau2 pool requires p_unexpressed >= 0.80
non-pool rows get omega = 1
```

PR 07 target behavior:

```text
row_weight_r = p_unexpressed_r * fused_information_r * support_weight_r
```

where `support_weight_r` can start as a bounded function of physical support count:

```text
support_weight_r = N_r / (N_r + N0)
```

or be folded into the fused posterior variance if that variance already captures support.

Rows with low `p_unexpressed` then contribute little, not nothing. Rows with ambiguous expression
contribute proportionally. Rows with strong fused gDNA information contribute more than noisy rows.

`omega` should still be neutralized for rows with no usable fused information, but that should be a
zero-information outcome rather than a hard expression-probability threshold.

## Background Density Refit After Fusion

Refit background density using fused mass and fused information:

```text
weight_r = p_unexpressed_r * fused_information_r * support_weight_r
rho0 = weighted robust estimate of gdna_mass_r / region_opportunity_r
```

This replaces the current source-method pool:

```text
method in {METHOD_STRAND, METHOD_BOUNDARY}
```

The source of information should not determine eligibility. The posterior uncertainty should.

Continue to preserve existing robust safeguards:

- exclude or downweight zero-opportunity regions,
- downweight spliced/expression evidence through `p_unexpressed`,
- trim or Huber-weight high-density outliers,
- carry prior/bootstrap when the effective pool is genuinely empty.

## Output And Diagnostics

Update calibration summaries to describe evidence fusion instead of tiers.

Recommended region-level diagnostics:

```text
gdna_mass
rna_mass
fused_variance
fused_precision
density_mean
density_variance
density_information
strand_mean
strand_variance
strand_information
density_weight
strand_weight
strand_applicable
strand_structural_absent
density_applicable
p_nonzero_gdna
```

Recommended summary JSON block:

```json
"strand_density_fusion": {
  "density_weight": "summary stats",
  "strand_weight": "summary stats",
  "fused_precision": "summary stats",
  "n_density_only": 0,
  "n_strand_informative": 0,
  "n_strand_structural_absent": 0,
  "n_near_unstranded": 0,
  "n_low_information": 0
}
```

Use explicit names like `density_weight` and `strand_weight`, not `tier` names.

## Implementation Phases

### Phase 1 - Evidence Contracts And Tests

- Add `StrandGdnaEvidence` and helpers in `strand_deconv.py` or a new `strand_evidence.py`.
- Add explicit density information/variance contract in `density_model.py`.
- Add tests for flat strand evidence, structural absence, and density uncertainty.
- Keep old production wiring during this phase.

Done when the evidence objects can be built and tested independently.

### Phase 2 - Fusion Engine

- Promote and revise `integration.py` into the production fusion module.
- Ensure exact fusion multiplies density and strand likelihoods.
- Ensure approximate fusion combines precisions and respects degenerate density cases.
- Add manual one-region tests where expected posterior can be computed directly.
- Add continuity tests over strand specificity and small positive gDNA density.

Done when fused outputs are smooth and source weights are diagnostic-only.

### Phase 3 - Production Mass Handoff

- Replace `build_region_unspliced_mass()` tier ownership with fused posterior mass.
- Delete or deprecate `METHOD_STRAND`, `METHOD_BOUNDARY`, and `METHOD_BACKGROUND_FALLBACK` from
  production logic.
- Remove `force_zero_gdna_mass` from `_orchestrator.py`, `calibration_e_step()`, and
  `run_calibration_iteration()`.
- Update `RegionUnsplicedMass` tests to assert fused behavior rather than tier promotion.

Done when calibration regional mass is always derived from `FusedRegionGdnaEvidence.mean_count`.

### Phase 4 - Iterative Calibration Refit

- Refit `BackgroundDensity` from fused mass and fused information.
- Refit `BackgroundModel` from fused mass instead of strand-only contained means.
- Refit exposure using continuous weights from `p_unexpressed` and fused information.
- Remove the PR 06 high-`p_unexpressed` tau2 gate once the continuous weighting tests pass.

Done when no scalar refit branches on the old source method.

### Phase 5 - Diagnostics, Goldens, And Cleanup

- Update output schemas and summaries to expose fusion diagnostics.
- Delete stale tier wording in docs and tests.
- Refresh goldens intentionally after targeted tests explain the behavior shift.
- Remove temporary compatibility aliases if they are no longer needed.

Done when outputs describe source uncertainty and fused posterior mass rather than tier ownership.

### Phase 6 - Scenario Validation And Benchmarks

- Run focused scenario tests for stranded, weak-stranded, unstranded, opposite-strand ambiguous,
  no-gDNA, tiny-gDNA, and boundary-only cases.
- Run full test suite.
- Run benchmark smoke where scratch benchmark data are available.
- Compare benchmark failures by source weights and fused uncertainty, not by tier method.

Done when the roadmap has evidence that remaining failures are biological identifiability or EM
competition issues, not calibration handoff cliffs.

## Test Plan

### Unit Tests

- `StrandSummary` weak contrast produces near-zero global information scale.
- Near-unstranded strand evidence returns flat likelihood and zero information.
- Opposite-strand ambiguous regions return `strand_applicable = false`.
- Strong strand-specific pure RNA region produces a sharp low-gDNA strand likelihood.
- Strong strand-specific mixed/gDNA region produces a sharp positive-gDNA strand likelihood.
- Density evidence exposes finite variance for observed-zero with finite opportunity.
- Sparse tiny-positive density evidence is positive but broad.
- Exact fusion equals manual grid multiplication.
- Density-only fusion equals bounded density posterior.
- Flat strand fusion equals density-only fusion.
- Strong strand can overcome a broad density prior when orientation evidence is decisive.
- Strong density can overcome weak/noisy strand when opportunity supports near-zero gDNA.
- Source weights change continuously as strand specificity moves from 0.50 to 1.00.
- Source weights change continuously as density moves from zero to tiny positive to low positive.

### Integration Tests

- No-gDNA RNA-only sentinel: calibration fused gDNA mass stays near zero without `force_zero`.
- Tiny-gDNA sentinel: fused gDNA mass is low positive, not snapped to zero or strand-owned.
- Low-strand-specificity nRNA case: strand weight is low and density controls mass.
- Strong-strand case: strand weight rises and improves RNA/gDNA separation.
- Unstranded boundary-only case: density imputes gDNA without strand.
- Opposite-strand ambiguous overlap: density-only fusion, no strand pseudo-signal.
- Hybrid capture skew: exposure learning receives fused mass and uncertainty.

### Continuity Tests

Continuity tests are critical for this PR. Add parameter sweeps that assert no large jumps in fused
mass or source weights unless the underlying posterior information changes sharply.

Recommended sweeps:

```text
p_r1_sense: 0.50, 0.55, 0.65, 0.80, 0.95, 1.00
gDNA fraction: 0, tiny, low, medium
boundary opportunity: none, sparse, strong
transcript strand class: unique, ambiguous opposite-strand
```

Assertions should focus on monotonicity, boundedness, and absence of cliffs rather than exact
golden numbers.

## Acceptance Criteria

- No production branch says "use strand if possible, else use density".
- Strand evidence always carries uncertainty or structural absence.
- Density evidence always carries uncertainty.
- Fused mass is the only regional gDNA mass source used by adaptive priors and exposure learning.
- Unstranded libraries produce flat or near-flat strand evidence, not confident gDNA masses.
- Opposite-strand ambiguous regions cannot use strand and reduce cleanly to density-only fusion.
- No-gDNA and tiny-gDNA cases transition smoothly.
- Former PR 06 `force_zero_gdna_mass` behavior is reproduced by density posterior concentration, not
  an external override.
- Tests cover source weights and continuity, not only final pass/fail scenario assertions.
- Full suite passes after intentional golden refresh.

## Risk Register

### Risk: Density And Strand Double Count The Same Fragments

Mitigation: define density as count/opportunity evidence and strand as orientation evidence
conditional on compatible mass. Avoid feeding strand-deconvolved boundary means into density and
then multiplying by strand likelihood again.

### Risk: Approximate Fusion Is Wrong Near Boundaries

Mitigation: use exact grid fusion for small/medium counts and adaptive exact fallback when the
approximate mode lands near 0 or `N_r`. Keep tests for large-count degenerate zero density.

### Risk: Zero-Density Cases Become Too Permissive

Mitigation: absence over large opportunity should yield a very low, high-information density
posterior. Validate no-gDNA sentinels and compare against PR 06 fixed behavior.

### Risk: Tiny Positive gDNA Gets Snapped To Zero

Mitigation: replace deterministic zero with finite low-density posteriors whenever opportunity is
finite. Add tiny-positive continuity tests.

### Risk: Exposure Refit Reintroduces Expression Leakage

Mitigation: use continuous `p_unexpressed * fused_information` weighting and retain diagnostics for
rows with high raw ratios but low expression confidence.

## Open Design Questions

1. Should `StrandGdnaEvidence.mean_count` be `nan` for flat evidence, or a bounded neutral diagnostic
   value? `nan` is safer but may require more schema handling.
2. Should global strand uncertainty be handled by a scalar information scale in PR 07, or should we
   immediately marginalize over the Beta posterior for `p_r1_sense`?
3. What finite prior should replace deterministic zero density when no positive anchors are observed
   but opportunity is large?
4. Should `RegionUnsplicedMass` be renamed now that it carries fused evidence diagnostics, or should
   the existing name survive for downstream stability?
5. How much fusion diagnostic output belongs in ordinary user-facing files versus summary JSON only?

## Done Means

The calibration pipeline no longer asks which estimator wins. It asks what each source says, how
uncertain each source is, and what posterior over regional gDNA mass follows from combining them.

Strand-specific data use strand powerfully. Unstranded data make strand uninformative. Ambiguous
opposite-strand regions omit strand structurally. Density remains active throughout. The regional
mass handed to priors and exposure is a fused posterior mean with inspectable uncertainty, not the
output of a tiered control-flow hierarchy.