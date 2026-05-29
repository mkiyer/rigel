# Adaptive Prior v2 - One-Knob Quantile RNA/gDNA Priors

Date: 2026-05-26
Status: implementation-ready design
Supersedes: `docs/prior/adaptive_prior_v1.md`
Incorporates review: `docs/prior/one_knob_parameterization.md`

## Executive Decision

Adopt the one-knob parameterization as the user-facing prior policy.

The grouped RNA/gDNA prior should be constructed from a posterior distribution
over the gDNA share of unspliced mass. The data determine both the posterior
location and its effective sample size. The user supplies only one scientific
preference:

```text
rna_confidence = q in (0, 1)
```

This is the posterior quantile level used to summarize the gDNA share when
building the grouped EM prior.

```text
share_for_prior_l = Quantile_q(phi_l | calibration)
```

where:

```text
phi_l = D_l / U_l
D_l   = latent gDNA-owned unspliced mass at locus l
U_l   = total prior-owned unspliced mass at locus l
```

The prior strength is not user-set. It is learned from posterior uncertainty:

```text
ess_l = posterior effective sample size of phi_l
alpha_gdna_add_l = ess_l * share_for_prior_l
alpha_rna_add_l  = ess_l * (1 - share_for_prior_l)
```

The v3 constants:

```text
aggregate_prior_strength
aggregate_prior_edge_count
aggregate_prior_max_count
gdna_prior_logit_bias
```

are demoted from user or conceptual policy to legacy compatibility and private
safety limits. They should not be exposed as ordinary configuration.

## Review Disposition

Accepted from `one_knob_parameterization.md`:

- The user should not tune prior strength, edge counts, max counts, or logit
  bias.
- The only legitimate scientific preference is asymmetric RNA-vs-gDNA risk.
- That risk is naturally represented as a posterior quantile of the gDNA share.
- ESS should be learned by moment matching from posterior variance.
- Sample-level empirical Bayes shrinkage should be learned from the sample.
- The old v3 constants should be deleted from the public surface and kept only
  as internal safety rails or legacy mode.

Clarified or modified:

- The quantile statement is exact for a continuous posterior share under
  asymmetric absolute loss. A binary RNA/gDNA classification loss would instead
  imply a posterior-probability threshold. Rigel is constructing a continuous
  prior share, so the quantile formulation is the correct one here.
- Existing `rna_lower_confidence` maps to the same upper-gDNA quantile level
  `q`, not to `1 - q`, when translated into the gDNA-share prior. If
  `R_lower(q)` is a lower credible bound for RNA mass, then
  `D_upper(q) = U - R_lower(q)` is the matching upper credible bound for gDNA
  mass.
- Internal numerical constants are not literally eliminated. They are removed
  from the user surface and constrained to safety caps, convergence guards, or
  debug-only overrides.
- Do not add a no-op `rna_confidence=0.5 -> v3 magic defaults` bridge as the
  final semantic contract. That would preserve behavior but lie about the
  meaning of the knob. Legacy v3 can remain as an explicit hidden compatibility
  mode during migration.

## User-Facing Contract

### Canonical Option

```text
--rna-confidence q
```

Range:

```text
0 < q < 1
```

Default:

```text
q = 0.5
```

Meaning:

```text
Use the q-th posterior quantile of the locus gDNA share as the grouped-prior
gDNA share.
```

Interpretation:

| `rna_confidence` | Summary | Behavior |
| --- | --- | --- |
| `0.5` | posterior median gDNA share | neutral under symmetric share-estimation loss |
| `0.9` | high gDNA-share quantile | conservative RNA; requires stronger evidence before prior protects RNA |
| `0.1` | low gDNA-share quantile | sensitive RNA; requires stronger evidence before prior protects gDNA |

The name is intentionally framed from the RNA reporting perspective: higher
`rna_confidence` means the tool uses a higher credible gDNA share when building
the prior, making RNA assignment more conservative in ambiguous unspliced
regions.

### Deprecated Aliases

For one release, keep these options as aliases:

```text
--rna-lower-confidence
--gdna-density-confidence
```

Alias rules:

- If `--rna-confidence` is provided, it is canonical.
- If only `--rna-lower-confidence` is provided, set `rna_confidence` to that
  value and emit a deprecation warning.
- If only `--gdna-density-confidence` is provided, set `rna_confidence` to that
  value and emit a deprecation warning.
- If multiple confidence flags are provided and differ by more than `1e-12`,
  fail fast with a clear error.

Do not expose `aggregate_prior_strength`, `aggregate_prior_edge_count`,
`aggregate_prior_max_count`, or `gdna_prior_logit_bias` on the CLI. Existing
YAML configs that contain them should be rejected or mapped only through an
explicit `legacy_v3` policy mode.

## Decision-Theoretic Basis

For a continuous share estimate `s` and true gDNA share `phi`, define the
asymmetric absolute loss:

```text
L_q(s, phi) = q       * max(phi - s, 0)
            + (1 - q) * max(s - phi, 0)
```

The Bayes action minimizing posterior expected loss is:

```text
s* = Quantile_q(phi | data)
```

Rigel's grouped prior needs exactly this kind of continuous share summary.
The posterior decides the distribution; `q` decides the user's risk preference.

This is cleaner than a strength-plus-bias system:

- the posterior mean and variance determine what the data know;
- the quantile determines the decision stance;
- the ESS determines how much the next EM stage should trust the summary.

## Core Posterior Object

For each locus, adaptive v2 constructs an approximate posterior:

```text
phi_l | calibration ~ Beta(a_l, b_l)
```

The final grouped prior is a decision-adjusted summary of this posterior:

```text
share_l = QBeta(q; a_l, b_l)
ess_l   = a_l + b_l

alpha_gdna_add_l = ess_l * share_l
alpha_rna_add_l  = ess_l * (1 - share_l)
```

Important distinction:

- `a_l` and `b_l` describe the calibration posterior over the gDNA share.
- `share_l` is a decision-theoretic quantile of that posterior.
- The EM receives a pseudo-count pair whose concentration is the posterior ESS
  but whose share is the risk-adjusted quantile.

This is intentional. The grouped EM prior is not a stored posterior sample; it
is a one-step summary passed from calibration into quantification under the
user's declared loss preference.

## Inputs Preserved From v3

Adaptive v2 keeps the v3 structural contract:

- `prior_unspliced_total_l`: only unspliced mass can feed prior alphas.
- `gdna_expected_count_l`: projected gDNA unspliced posterior mean.
- `rna_expected_count_l`: projected RNA unspliced posterior mean.
- `has_gdna_candidate_l`: native structural candidate availability.
- Spliced fragments remain likelihood evidence only.
- All transcript-like RNA components remain one aggregate prior group.

Adaptive v2 adds count-scale uncertainty:

```text
gdna_expected_var_l = Var(D_l | calibration)
gdna_share_var_l    = Var(phi_l | calibration)
```

## Required Data Structure Changes

### `RegionGdnaEstimate`

Add variance to the strand deconvolution summary:

```python
@dataclass(frozen=True, slots=True)
class RegionGdnaEstimate:
    n_total: np.ndarray
    mean_count: np.ndarray
    var_count: np.ndarray       # NEW: Var(D | strand evidence)
    upper_count: np.ndarray
    rna_lower_count: np.ndarray
    precision: np.ndarray
    flags: np.ndarray
```

The exact path already computes `sd_r`; store:

```text
var_count = sd_r^2
```

Because `D = N - R`, RNA and gDNA count variances are equal under fixed
observed `N`.

The normal path also computes `sd_r`; store the same quantity.

Degenerate cases:

- no observed unspliced mass: `mean_count = 0`, `var_count = 0`;
- strand-ineligible but preserved total: defer to density fallback if strand is
  globally unusable;
- near-unstranded exact strand path: mark uninformative and use a wide fallback
  variance if the value is ever consumed directly.

### `RegionGdnaChannelEstimate`

Add per-compartment variances:

```python
contained_var: np.ndarray
boundary_left_var: np.ndarray
boundary_right_var: np.ndarray
```

Regional variance under strand-informative calibration is:

```text
Var(D_r) = contained_var_r
         + boundary_left_var_r
         + boundary_right_var_r
```

This assumes independent compartment posterior approximations. Record this as
an approximation in diagnostics.

### `PriorMassDeconvolution`

Extend the mass-conserving prior payload:

```python
@dataclass(frozen=True, slots=True)
class PriorMassDeconvolution:
    unspliced_total: np.ndarray
    gdna_unspliced_mean: np.ndarray
    gdna_unspliced_var: np.ndarray       # NEW
    rna_unspliced_mean: np.ndarray
    method: np.ndarray
    precision: np.ndarray
    flags: np.ndarray
```

Optional convenience property:

```python
rna_unspliced_var = gdna_unspliced_var
```

The invariant remains:

```text
gdna_unspliced_mean + rna_unspliced_mean == unspliced_total
```

Variance must be finite and non-negative.

### `PriorTable`

Add adaptive diagnostics:

```python
gdna_expected_var: np.ndarray
prior_gdna_share_mean: np.ndarray
prior_gdna_share_var: np.ndarray
prior_gdna_share_quantile: np.ndarray
prior_ess_local: np.ndarray
prior_ess_global: np.ndarray
prior_ess_final: np.ndarray
prior_alpha_source: np.ndarray
prior_variance_inflation: np.ndarray
prior_conflict_score: np.ndarray | None
```

Keep v3 columns for one release:

```text
alpha_gdna_add
alpha_rna_add
prior_budget
prior_gdna_share_biased
```

but mark `prior_budget` as an alias for `prior_ess_final` and
`prior_gdna_share_biased` as a legacy name for `prior_gdna_share_quantile` in
adaptive mode.

## Density-Fallback Uncertainty

When strand information is unavailable, calibration currently builds:

```text
mu_gdna = p_captured * captured_mu
        + (1 - p_captured) * off_target_mu
```

Adaptive v2 must also build:

```text
Var(D) = p_captured * Var(D | captured)
       + (1 - p_captured) * Var(D | off_target)
       + p_captured * (1 - p_captured)
         * (captured_mu - off_target_mu)^2
```

For a Gamma-Poisson posterior predictive model:

```text
rho ~ Gamma(alpha, beta)
D | rho ~ Poisson(rho * L)

E[D]   = alpha * L / beta
Var[D] = E[D] + E[D]^2 / alpha
```

Off-target component:

```text
alpha_off = background.rho_off_alpha
beta_off  = background.rho_off_beta
L         = contained_leff
```

Captured component:

```text
alpha_cap = background.rho_off_alpha + sweep.alpha_excess
beta_cap  = background.rho_off_beta  + sweep.beta_excess
L         = contained_leff
```

If only a mean and upper bound are available during migration, use an
interval-derived variance only as a temporary fallback. Use the confidence
level that created the stored upper bound, not the user `rna_confidence`
quantile. This avoids the invalid `z = 0` case at `q = 0.5`.

```text
sd_fallback = max(upper - mean, 0) / z_interval
var_fallback = sd_fallback^2
```

and set a diagnostic flag:

```text
FLAG_PRIOR_VARIANCE_FROM_INTERVAL
```

## Projection To Locus Posterior Moments

Existing v3 allocation computes overlap shares `s_lr` from region `r` to locus
`l`. Extend it to project variance:

```text
U_l        = sum_r s_lr * U_r
mu_D_l     = sum_r s_lr * mu_D_r
sigma2_D_l = sum_r s_lr^2 * sigma2_D_r
```

The square on `s_lr` is required for independent variance propagation.

Inflate variance for projection ambiguity:

```text
allocation_gap_l = max(0, 1 - prior_allocated_fraction_l)
shared_mass_l    = multi_locus_region_mass_l + partial_coverage_region_mass_l
shared_frac_l    = shared_mass_l / max(U_l, eps)

sigma2_D_l += (allocation_gap_l^2 + shared_frac_l^2) * U_l^2
```

This lowers ESS; it does not gate gDNA.

## Moment Matching Helpers

Implement in a small pure-Python helper module, for example:

```text
src/rigel/calibration/adaptive_prior.py
```

Required functions:

```python
def beta_from_mean_var(
    mean: np.ndarray,
    var: np.ndarray,
    total: np.ndarray,
    *,
    max_ess: float,
) -> BetaMoments
def beta_quantile(alpha: np.ndarray, beta: np.ndarray, q: float) -> np.ndarray
def combine_beta(alpha1, beta1, alpha2, beta2) -> tuple[np.ndarray, np.ndarray]
```

Moment matching:

```text
p = clamp(mean, 0, 1)
v_max = p * (1 - p)
v_floor = p * (1 - p) / (total + 1)
v = clamp(var, v_floor, v_max)

tau = p * (1 - p) / v - 1
tau = clamp(tau, 0, max_ess)

alpha = tau * p
beta  = tau * (1 - p)
```

Boundary handling:

- If `total <= eps`, set `tau = 0`, `alpha = 0`, `beta = 0`.
- If `p <= eps` and `tau > 0`, treat the posterior as degenerate at zero for
  quantile evaluation unless a smoothing floor is explicitly requested for
  diagnostics.
- If `p >= 1 - eps` and `tau > 0`, treat it as degenerate at one.
- For nondegenerate `alpha, beta > 0`, compute quantiles with
  `scipy.special.betaincinv(alpha, beta, q)`.

Do not add hidden pseudo-count floors to `alpha_gdna_add` or `alpha_rna_add`.
Any floor used only for numerical quantile evaluation must not become modeled
prior mass.

## Sample-Level Empirical Bayes Shrinkage

Use regional or locus-level share moments to estimate a global gDNA-share
distribution:

```text
phi_l ~ Beta(psi * kappa, (1 - psi) * kappa)
```

Estimate `psi` from inverse-variance-weighted shares:

```text
p_i = mu_D_i / U_i
v_i = sigma2_D_i / U_i^2
w_i = 1 / max(v_i, p_i * (1 - p_i) / (U_i + 1), eps)
w_i = min(w_i, U_i + 1)

psi = (sum_i w_i * p_i + 0.5) / (sum_i w_i + 1.0)
```

Estimate heterogeneity after subtracting measurement variance:

```text
observed_var     = weighted_var(p_i, w_i)
measurement_var  = weighted_mean(v_i, w_i)
between_var      = max(observed_var - measurement_var, eps)
kappa            = psi * (1 - psi) / between_var - 1
```

Clamp only by private safety rails:

```text
kappa = clamp(kappa, 0, _MAX_GLOBAL_ESS)
```

If the method-of-moments estimate is degenerate, fall back to:

```text
kappa = min(sum_i w_i, _MAX_GLOBAL_ESS)
```

and record:

```text
global_ess_fallback_used = true
```

Do not use fixed method weights for strand vs density evidence. Method labels
are diagnostics. Their influence should come through posterior variance.

## Locus-Level Posterior Combination

For each locus, construct a local Beta approximation:

```text
a_local_l = tau_local_l * p_local_l
b_local_l = tau_local_l * (1 - p_local_l)
```

Construct the empirical-Bayes global contribution:

```text
tau_global_l = min(kappa, _MAX_GLOBAL_BORROW_ESS)
a_global_l = tau_global_l * psi
b_global_l = tau_global_l * (1 - psi)
```

Combine:

```text
a_post_l = a_local_l + a_global_l
b_post_l = b_local_l + b_global_l
```

Apply structural and reliability effects to ESS, not to share:

```text
structural_l = 1 if native has structural gDNA candidate else 0
prior_mass_l = 1 if U_l > eps else 0
allocation_reliability_l = clamp(prior_allocated_fraction_l, 0, 1)
convergence_reliability = 1 if calibration converged else 0.5
variance_completeness_l = 1 unless interval/precision fallback was required

reliability_l = structural_l
              * prior_mass_l
              * allocation_reliability_l
              * convergence_reliability
              * variance_completeness_l
```

Final ESS:

```text
ess_l = min(_MAX_ESS, reliability_l * (a_post_l + b_post_l))
```

Final share:

```text
share_l = QBeta(q; a_post_l, b_post_l)
```

Final grouped prior:

```text
alpha_gdna_add_l = ess_l * share_l
alpha_rna_add_l  = ess_l * (1 - share_l)
```

If `structural_l = 0` or `prior_mass_l = 0`, force:

```text
alpha_gdna_add_l = 0
alpha_rna_add_l  = 0
ess_l            = 0
```

This preserves v3 native candidate semantics.

## Why One-Sided Evidence Works Without `edge_count`

One-sided evidence is no longer a special case.

If calibration confidently says all unspliced mass is RNA:

```text
p_l ~= 0
Var(phi_l) small
tau_l large
QBeta(q) near 0 for any ordinary q < 1
alpha_rna_add large
```

If calibration weakly says all unspliced mass is RNA:

```text
p_l ~= 0
Var(phi_l) large
tau_l small
alpha_rna_add small unless global evidence is strong
```

The old `aggregate_prior_edge_count` was compensating for the lack of an
uncertainty-aware one-sided posterior. Quantile plus ESS makes it unnecessary.

## Internal Safety Rails

These constants may exist in code, but not on the ordinary CLI:

```python
_MAX_ESS = 3000.0
_MAX_GLOBAL_ESS = 3000.0
_MAX_GLOBAL_BORROW_ESS = 1000.0
_MIN_VARIANCE_EPS = 1.0e-12
```

They are not model-tuning parameters. They guard numerical pathologies and
should be audited by stress tests. Changing them should be rare and documented.

The only user-facing model knob is:

```text
rna_confidence
```

## Config And CLI Implementation

### `config.py`

Add canonical confidence to `CalibrationConfig`:

```python
rna_confidence: float = 0.5
```

Validation:

```python
0.0 < rna_confidence < 1.0
```

Keep these legacy fields for one release, but mark deprecated:

```python
rna_lower_confidence: float | None = None
gdna_density_confidence: float | None = None
```

Resolution helper:

```python
def resolved_rna_confidence(self) -> float:
    values = [rna_confidence if explicitly set]
    values += [rna_lower_confidence if not None]
    values += [gdna_density_confidence if not None]
    if values disagree: raise ValueError
    return canonical value or 0.5
```

For migration from the current non-optional fields, the resolver can initially
use sentinel defaults in CLI parsing rather than dataclass defaults, so old
YAML files continue to load with warnings.

### `cli.py`

Add:

```text
--rna-confidence
```

Keep existing confidence flags for one release but change help text to say:

```text
Deprecated alias for --rna-confidence.
```

Unknown or removed prior-tuning YAML keys should produce actionable errors:

```text
aggregate_prior_strength is no longer user-configurable; use --rna-confidence.
```

### Calibration Orchestrator

Use the resolved `q` for posterior summaries that are decision summaries:

```python
q = calibration_config.resolved_rna_confidence()
```

Do not use `q` as a generic numerical threshold when a step is model fitting
rather than posterior summarization. In particular, exon self-training and
calibration damping should eventually be converted to posterior-weighted or
ELBO-based procedures. During migration, keep their current internal defaults
but remove them from user surface.

This distinction prevents the one user knob from becoming a disguised tuning
constant for unrelated numerical procedures.

## Output Diagnostics

Add to `summary.json`:

```text
prior_policy: "adaptive_quantile"
rna_confidence: q
global_gdna_share: psi
global_gdna_ess: kappa
global_ess_fallback_used
prior_variance_source_histogram
prior_method_histogram
prior_ess_final summary
prior_gdna_share_quantile summary
prior_conflict_score summary
```

Add to `loci.feather`:

```text
prior_gdna_share_mean
prior_gdna_share_var
prior_gdna_share_quantile
prior_ess_local
prior_ess_global
prior_ess_final
alpha_gdna_add
alpha_rna_add
prior_variance_inflation
prior_reliability
prior_conflict_score
```

After EM, compute conflict:

```text
em_gdna_share_l = gdna_em_count_l / max(total_em_count_l, eps)
prior_conflict_score_l = abs(logit_clip(em_gdna_share_l)
                           - logit_clip(prior_gdna_share_quantile_l))
```

Large conflicts are diagnostics, not automatic gates.

## Implementation Phases

### Phase 0 - Golden-Free Unit Tests

Before implementation, add tests for pure helper math:

- `beta_from_mean_var()` returns monotone ESS as variance changes.
- ESS is finite and capped.
- Degenerate all-RNA and all-gDNA shares do not create opposite-side floors.
- `beta_quantile(q)` is monotone in `q`.
- `q=0.1`, `0.5`, `0.9` produce ordered gDNA prior shares for identical
  posterior inputs.
- Structural no-gDNA candidates force both alphas to zero.

### Phase 1 - Preserve Variance In Calibration

Implement `var_count` in strand deconvolution and density fallback.

Update tests:

- exact posterior variance matches direct discrete posterior variance;
- normal posterior variance is finite and decreases with stronger strand
  contrast;
- Gamma-Poisson variance is greater than or equal to mean;
- mixture variance includes the between-component term.

### Phase 2 - Extend Prior Projection

Project variance through `_allocate_counts_by_geometry()` using `s_lr^2`.

Add tests:

- disjoint region variance is conserved;
- split region variance is allocated by squared shares;
- allocation ambiguity inflates variance and reduces ESS.

### Phase 3 - Empirical-Bayes Global Profile

Add `PriorEvidenceProfile` with:

```python
global_gdna_share
global_gdna_ess
global_ess_fallback_used
n_regions_used
share_observed_variance
share_measurement_variance
share_between_region_variance
variance_source_histogram
```

Tests:

- no-gDNA high-confidence regions produce `psi` near zero;
- high-gDNA regions increase `psi`;
- heterogeneous regions lower `kappa`;
- missing variance triggers fallback diagnostics.

### Phase 4 - Adaptive Quantile Priors

Implement `compute_adaptive_quantile_grouped_prior_counts()` and dispatch from
`assemble_priors()`.

Keep a hidden `legacy_v3` path for comparison, but make adaptive quantile the
default once tests are green.

Tests:

- q monotonicity at fixed posterior;
- one-sided high-confidence RNA produces strong RNA alpha without edge count;
- one-sided low-confidence RNA produces weak alpha;
- global no-gDNA profile shrinks ambiguous loci toward RNA;
- strong local gDNA evidence can override global low-gDNA shrinkage.

### Phase 5 - CLI And Config Migration

- Add `--rna-confidence`.
- Deprecate `--rna-lower-confidence` and `--gdna-density-confidence`.
- Remove prior magic constants from public config/YAML.
- Write warnings for old confidence aliases.
- Fail on conflicting confidence flags.

### Phase 6 - Output And Golden Migration

- Add adaptive diagnostics to `loci.feather` and `summary.json`.
- Keep legacy columns for one release with documented alias semantics.
- Regenerate goldens intentionally.

### Phase 7 - Validation And Audit

Run:

```text
pytest tests/ -v
pytest tests/test_golden_output.py --update-golden -v
```

Run sentinel sweeps at:

```text
rna_confidence in {0.1, 0.5, 0.9}
```

Required monotonicity checks:

- `q=0.9` produces greater or equal gDNA prior share than `q=0.5` for every
  locus with nonzero ESS.
- `q=0.1` produces less or equal gDNA prior share than `q=0.5` for every locus
  with nonzero ESS.
- No retuning of hidden caps is allowed across q values.

Benchmark checks:

- zero-gDNA samples recover RNA without fixed logit bias;
- high-gDNA samples do not globally suppress gDNA;
- unstranded/no-global-evidence single-exon cases report low ESS and ambiguity;
- prior-vs-EM conflicts cluster in biologically or statistically ambiguous
  loci, not in clean sentinels.

## Failure Modes And Correct Behavior

### Miscalibrated Regional Posterior

If calibration confuses nRNA with gDNA, the quantile cannot repair the upstream
posterior. Correct behavior is to expose high conflict and method diagnostics,
then fix calibration.

### Truly Unidentifiable Loci

If local and global evidence are both weak, ESS should be low. The EM should be
allowed to decide from likelihood evidence, and outputs should report ambiguity.

### Strong Local Evidence Against Global Trend

Local high-ESS evidence must be able to override global shrinkage. This is a
required test case, especially for local gDNA hotspots in otherwise clean
libraries.

## Acceptance Criteria

The adaptive quantile plan is complete when:

- ordinary users see one prior-related model knob: `--rna-confidence`;
- v3 magic constants are unavailable in ordinary CLI/YAML config;
- prior ESS is derived from posterior variance;
- prior share is derived from the requested posterior quantile;
- confidence aliases resolve consistently and warn;
- all alphas remain structurally inactive when no gDNA candidate exists;
- q-response is monotone without retuning;
- diagnostics explain prior strength, prior share, shrinkage, and conflicts.

## Summary

Adaptive prior v2 turns the grouped RNA/gDNA prior into a decision-theoretic
posterior handoff:

```text
the data learn the posterior;
the posterior variance learns the ESS;
empirical Bayes learns the sample-level shrinkage;
the user supplies one quantile-level risk preference.
```

Everything else is legacy compatibility, diagnostics, or private numerical
safety.