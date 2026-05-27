# Adaptive Prior v1 - Uncertainty-Calibrated RNA/gDNA Group Priors

Date: 2026-05-26
Status: design proposal for the post-v3 prior operating policy
Builds on: `docs/prior/prior_redesign_v3.md`

## Executive Summary

The v3 prior implementation created the right native object: a grouped prior
over aggregate RNA vs aggregate gDNA. It also exposed the wrong operating
surface: four fixed constants that encode a global policy about how hard the
calibration prior should push every locus.

Adaptive Prior v1 replaces those constants as the primary policy with a
sample- and locus-specific posterior effective sample size (ESS). The prior
share and prior strength should be estimated from calibration uncertainty:

```text
calibration posterior for locus l
    -> posterior mean gDNA share p_l
    -> posterior variance Var(p_l)
    -> Beta-equivalent ESS tau_l
    -> alpha_gdna_add = tau_l * p_l
       alpha_rna_add  = tau_l * (1 - p_l)
```

The existing v3 constants remain only as internal safety limits for pathological
or uncertainty-missing cases. They are no longer user-facing conceptual knobs.

The user should not set `aggregate_prior_strength`,
`aggregate_prior_edge_count`, `aggregate_prior_max_count`, or
`gdna_prior_logit_bias`. The default policy should be:

```text
prior_policy = "adaptive"
```

and the output should explain what the sample taught Rigel: global gDNA rate,
regional calibration precision, per-locus prior ESS, shrinkage source, and
prior-vs-EM conflict diagnostics.

## Non-Negotiable Constraints

This design preserves the v3 modeling contract:

- Only unspliced mass contributes calibration-derived EM prior pseudocounts.
- Spliced fragments remain likelihood evidence and never become prior count
  mass.
- Internal inference remains transcript-first.
- All annotated mRNA and synthetic nRNA components remain one aggregate RNA
  group for the prior.
- Gene-level summaries remain output-only.
- Native gDNA availability remains structural: a gDNA component is active only
  where an unspliced fragment unit has finite gDNA likelihood.
- No posterior threshold gate such as `gdna_upper > eps` is reintroduced.

## Problem With v3 Operating Point

The current v3 budget is:

```text
m_g = projected gDNA unspliced mean
m_r = projected RNA unspliced mean
m   = m_g + m_r

balanced_budget = 2 * min(m_g, m_r)
edge_budget     = min(m, aggregate_prior_edge_count)
budget_raw      = max(balanced_budget, edge_budget)
budget_raw      = min(budget_raw, aggregate_prior_max_count)
budget          = min(aggregate_prior_strength * budget_raw,
                      aggregate_prior_max_count)

w_raw           = m_g / m
w_biased        = sigmoid(logit(w_raw) + gdna_prior_logit_bias)

alpha_gdna_add  = budget * w_biased
alpha_rna_add   = budget * (1 - w_biased)
```

This passed the sentinel suite after tuning:

```text
aggregate_prior_strength = 3.0
aggregate_prior_edge_count = 1000.0
aggregate_prior_max_count = 3000.0
gdna_prior_logit_bias = -6.0
```

But those values are not theoretically interpretable by a user. The four
constants combine three different concepts:

1. how much calibration evidence exists;
2. how uncertain that evidence is;
3. how conservative Rigel should be against false gDNA.

Adaptive v1 separates these concepts and estimates them from the sample.

## Statistical Target

For each `MultiLocus` `l`, define:

```text
U_l = prior-owned unspliced fragment mass allocated to locus l
D_l = latent gDNA-owned unspliced mass at locus l
R_l = U_l - D_l
phi_l = D_l / U_l
```

Calibration should provide an approximate posterior over `D_l`, or equivalently
over `phi_l`:

```text
E[D_l | calibration]      = mu_D_l
Var(D_l | calibration)    = sigma2_D_l
E[phi_l | calibration]    = p_l = mu_D_l / U_l
Var(phi_l | calibration)  = v_l = sigma2_D_l / U_l^2
```

Then approximate the posterior over `phi_l` with a Beta distribution:

```text
phi_l ~ Beta(a_l, b_l)
```

Moment matching gives:

```text
tau_l = a_l + b_l
p_l   = a_l / tau_l
v_l   = p_l * (1 - p_l) / (tau_l + 1)

tau_l = p_l * (1 - p_l) / v_l - 1
a_l   = tau_l * p_l
b_l   = tau_l * (1 - p_l)
```

The native grouped EM consumes additive pseudocounts, so:

```text
alpha_gdna_add[l] = a_l
alpha_rna_add[l]  = b_l
```

When variance is high, `tau_l` becomes small and the prior becomes weak. When
variance is low, `tau_l` becomes large, subject only to safety caps. No user
needs to guess a prior strength.

## Required Calibration Uncertainty Payload

The current `PriorMassDeconvolution` stores means and a coarse precision scalar:

```python
class PriorMassDeconvolution:
    unspliced_total: np.ndarray
    gdna_unspliced_mean: np.ndarray
    rna_unspliced_mean: np.ndarray
    method: np.ndarray
    precision: np.ndarray
    flags: np.ndarray
```

Adaptive v1 needs count-scale uncertainty. Extend it to carry at least:

```python
class PriorMassDeconvolution:
    unspliced_total: np.ndarray
    gdna_unspliced_mean: np.ndarray
    gdna_unspliced_var: np.ndarray
    rna_unspliced_mean: np.ndarray
    rna_unspliced_var: np.ndarray
    method: np.ndarray
    precision: np.ndarray
    flags: np.ndarray
```

Because `R = U - D`, `Var(R) = Var(D)` when `U` is treated as observed. Both
fields may be stored for clarity, but the implementation can store only
`gdna_unspliced_var` if the dataclass exposes an accessor for RNA variance.

### Strand-Informative Regions

The strand deconvolution already computes the needed posterior uncertainty.
For exact regions, `_summarize_exact()` computes `sd_r` from the discrete
posterior over RNA count. For large regions, `_summarize_normal()` computes a
normal-approximation `sd_r`.

Since `D = N - R`:

```text
Var(D) = Var(R) = sd_r^2
```

Implementation change:

- Add `sd_count` or `var_count` to `RegionGdnaEstimate`.
- Add compartment variance fields to `RegionGdnaChannelEstimate`:

```python
contained_var
boundary_left_var
boundary_right_var
```

- Build regional prior variance as:

```text
sigma2_D_r = contained_var[r]
           + boundary_left_var[r]
           + boundary_right_var[r]
```

This assumes independent compartment deconvolution. That is an approximation,
but it is substantially better than discarding uncertainty. If later needed,
the exact compartment covariance can be modeled at the ledger level.

### Density-Fallback Regions

When strand deconvolution is unavailable, `mu_gdna` is currently a mixture of
captured and off-target density expectations:

```text
mu_gdna = p_captured * captured_mu
        + (1 - p_captured) * off_target_mu
```

The variance should include both component count variance and state-mixture
uncertainty:

```text
Var(D) = p_captured * Var(D | captured)
       + (1 - p_captured) * Var(D | off_target)
       + p_captured * (1 - p_captured)
         * (captured_mu - off_target_mu)^2
```

The off-target component follows the background Gamma-Poisson posterior. If
`rho ~ Gamma(alpha, beta)` and `D | rho ~ Poisson(rho * L)`, then the posterior
predictive count variance is:

```text
mean = alpha * L / beta
var  = mean + mean^2 / alpha
```

The captured component uses the same Gamma-Poisson form with the sweep excess
posterior:

```text
alpha_captured = background.rho_off_alpha + sweep.alpha_excess
beta_captured  = background.rho_off_beta  + sweep.beta_excess
mean_captured  = alpha_captured * L / beta_captured
var_captured   = mean_captured + mean_captured^2 / alpha_captured
```

If a component only has a mean and one-sided upper bound during an intermediate
migration, approximate its standard deviation as:

```text
sd ~= max(upper - mean, 0) / z_confidence
var = sd^2
```

This interval-derived fallback should be marked with a flag. It is a migration
bridge, not the final uncertainty source.

### Minimum Variance Floor

Moment matching is unstable when `v_l` is zero. Use a variance floor that
reflects finite count resolution, not an arbitrary policy:

```text
v_floor_l = p_l * (1 - p_l) / (U_l + 1)
v_l       = max(v_l, v_floor_l, eps)
```

This prevents infinite ESS. It also ensures a locus with 10 unspliced fragments
cannot receive 3000 effective prior fragments merely because a posterior
approximation rounded its variance to zero.

## Projection From Regions To Loci

Existing geometry allocation projects means by region overlap share `s_lr`.
Adaptive v1 projects variance at the same time:

```text
mu_D_l      = sum_r s_lr * mu_D_r
U_l         = sum_r s_lr * U_r
sigma2_D_l = sum_r s_lr^2 * sigma2_D_r
```

The `s_lr^2` term is correct for independent region-level uncertainty.

Add allocation uncertainty when projection quality is imperfect:

```text
q_l = prior_allocated_fraction[l]
c_l = fraction of prior mass touching multiple loci or partial coverage
```

Then inflate variance:

```text
sigma2_D_l <- sigma2_D_l
              + allocation_penalty_l * U_l^2
```

where:

```text
allocation_penalty_l = max(0, 1 - q_l)^2 + c_l^2
```

This does not gate gDNA. It simply lowers prior ESS when region-to-locus
projection is ambiguous.

## Sample-Level Empirical Bayes Shrinkage

Local calibration can be weak for single-exon or unstranded loci. Those loci
should borrow strength from the sample's global gDNA profile, not from fixed
magic constants.

Estimate a sample-level gDNA share from high-quality prior-mass regions:

```text
p_r = mu_D_r / U_r
v_r = Var(D_r) / U_r^2
psi = sum_r w_r * p_r / sum_r w_r
```

Use inverse-variance weights capped by finite count information:

```text
w_r = 1 / max(v_r, p_r * (1 - p_r) / (U_r + 1), eps)
w_r = min(w_r, U_r + 1)
```

This is deliberately not a fixed method weight. Strand and density evidence
earn influence by carrying low posterior variance. Method labels remain
diagnostics and fallback indicators, not primary trust constants.

If a migration phase lacks variance for a method, use a conservative fallback
weight and mark the profile as variance-incomplete:

```text
fallback_weight_r = clamp(precision_r, 0, 1)^2 * min(U_r + 1, max_fallback_region_weight)
```

Do not use expression or gene labels to define the global prior. The estimate
comes from calibration states, strand deconvolution, density model evidence,
and unspliced regional mass.

Estimate global uncertainty with a Beta-equivalent ESS:

```text
psi_smooth = (sum w_r * p_r + 0.5) / (sum w_r + 1.0)
tau_global = effective information in the weighted regions
```

Primary implementation:

```text
tau_global = kappa_global
```

where `kappa_global` is the empirical-Bayes concentration of the sample-level
Beta distribution over regional gDNA shares. Estimate it by method of moments
after subtracting measurement variance:

```text
observed_var = weighted_var(p_r)
measurement_var = weighted_mean(v_r)
between_region_var = max(observed_var - measurement_var, eps)
kappa_global = psi * (1 - psi) / between_region_var - 1
```

This `kappa_global` describes real heterogeneity in gDNA share across regions.
It should not be forced large in heterogeneous libraries.

If the between-region variance estimate is unavailable or numerically
degenerate, fall back to a conservative global concentration and mark it as a
fallback:

```text
tau_global = min(sum_r w_r, max_global_ess)
```

and set `global_ess_fallback_used = true` in diagnostics.

## Locus Posterior Share

For each locus, compute local moment-matched values:

```text
p_local = clamp(mu_D_l / U_l, eps, 1 - eps)
v_local = clamp(sigma2_D_l / U_l^2, v_floor_l, p_local*(1-p_local))
tau_local = max(p_local * (1 - p_local) / v_local - 1, 0)
```

Shrink the local share toward the sample-level share on the logit scale:

```text
lambda_l = tau_local / (tau_local + tau_global_l)
eta_l    = lambda_l * logit(p_local)
         + (1 - lambda_l) * logit(psi)
p_l      = sigmoid(eta_l)
```

`tau_global_l` should be smaller for loci whose local calibration is already
precise and larger when local evidence is weak but the sample-level estimate is
strong. A practical first version:

```text
tau_global_l = min(tau_global, U_l, max_global_borrow_ess)
```

This replaces `gdna_prior_logit_bias`. In a no-gDNA sample, `psi` is near zero
with high confidence, so ambiguous loci get an RNA-conservative prior. In a
gDNA-heavy sample, `psi` is high, so Rigel does not globally suppress gDNA.

## Locus Prior ESS

The final ESS should represent posterior confidence, not raw expression scale.

Compute the posterior variance after combining local and global information.
For a first implementation, approximate precision additively:

```text
tau_combined = tau_local + tau_global_l
```

Then apply reliability multipliers:

```text
reliability_l = structural_l
              * allocation_reliability_l
              * convergence_reliability
              * variance_completeness_l
```

where:

```text
structural_l = 1 if native has a structural gDNA candidate else 0
allocation_reliability_l = clamp(prior_allocated_fraction_l, 0, 1)
convergence_reliability = 1 if calibration converged else e.g. 0.5
variance_completeness_l = 1 when count variance is available;
                          <1 only for interval/precision fallbacks
```

Do not multiply by the old `precision` scalar when count-scale variance is
available. The variance already controls `tau_local`. The precision scalar is
only a migration fallback and a human-readable diagnostic.

The final ESS is:

```text
tau_l = min(
    max_ess,
    reliability_l * tau_combined,
)
```

Then:

```text
alpha_gdna_add = tau_l * p_l
alpha_rna_add  = tau_l * (1 - p_l)
```

If no structural gDNA candidate exists, both alphas remain zero. This preserves
the v3 candidate contract.

## One-Sided Evidence Without Magic Edge Counts

The v3 `edge_count` exists because one-sided evidence can be important:
calibration may say nearly all unspliced mass is RNA, but `2 * min(m_g, m_r)`
would give no prior mass.

Adaptive v1 handles this through uncertainty, not a fixed edge count:

- If calibration confidently says `p_l ~= 0`, then `Var(p_l)` is small,
  `tau_l` is large, and `alpha_rna_add` is strong.
- If calibration weakly says `p_l ~= 0`, then `Var(p_l)` is large, `tau_l` is
  small, and the prior is weak unless the sample-level global gDNA estimate is
  itself precise.
- If both local and global evidence are weak, Rigel should report ambiguity
  rather than pretending a magic count can identify RNA or gDNA.

This is the correct replacement for `aggregate_prior_edge_count`.

## Safety Limits And Legacy Constants

The old constants remain internal fallback limits, not user-facing tuning
parameters:

| v3 constant | v4 role |
| --- | --- |
| `aggregate_prior_strength` | Legacy-mode multiplier only; in adaptive mode, at most a hidden cap on fallback interval-derived ESS. |
| `aggregate_prior_edge_count` | Maximum fallback one-sided ESS when variance is missing. Adaptive mode should not use it when posterior variance is available. |
| `aggregate_prior_max_count` | Absolute `max_ess` safety cap. This remains useful. |
| `gdna_prior_logit_bias` | Legacy-mode offset only. Adaptive mode replaces it with sample-level empirical Bayes shrinkage toward `psi`. |

Recommended config shape:

```python
class PriorPolicyConfig:
    mode: Literal["adaptive", "legacy_v3", "off"] = "adaptive"
    max_ess: float = 3000.0
    max_global_ess: float = 3000.0
    max_global_borrow_ess: float = 1000.0
    max_fallback_one_sided_ess: float = 1000.0
    max_fallback_region_weight: float = 100.0
    min_variance_eps: float = 1.0e-9
```

The CLI should expose only the mode, if anything:

```text
--prior-policy adaptive|legacy-v3|off
```

Do not expose the numeric limits to ordinary users. They are developer safety
rails and benchmark audit targets.

## Prior-vs-EM Conflict Diagnostics

Adaptive priors should produce their own audit trail. Add per-locus diagnostics:

```text
prior_policy
prior_gdna_share_local
prior_gdna_share_global
prior_gdna_share_final
prior_ess_local
prior_ess_global_borrowed
prior_ess_final
prior_variance_local
prior_variance_inflated
prior_reliability
prior_method_mix_strand
prior_method_mix_density
prior_conflict_score
```

After EM, compute:

```text
em_gdna_share = gdna_em_count_locus / max(total_em_count_locus, eps)
prior_conflict_score = abs(logit(em_gdna_share) - logit(prior_gdna_share_final))
```

Flag but do not auto-silence large conflicts:

```text
prior_conflict_high = prior_ess_final high and conflict_score high
```

These cases should appear in `summary.json` and optionally in
`locus_stats.feather`. They are exactly the loci users and developers should
inspect.

## Optional Two-Pass Self-Calibration

Adaptive v1 does not require rerunning EM, but v2 can add a safe second pass:

1. Build adaptive priors and run EM.
2. Identify high-confidence prior-vs-EM conflict classes.
3. Re-estimate global shrinkage dispersion or downweight low-reliability
   methods.
4. Re-run EM once.

This must not become a hidden threshold gate. It is empirical Bayes
hyperparameter refinement, and every adapted value must be recorded.

## Implementation Plan

### Phase 0 - Tests And Golden-Free Acceptance Criteria

Add tests before changing outputs:

- Moment matching converts `(mean, variance)` to finite Beta ESS.
- ESS decreases monotonically as variance increases.
- ESS is capped by `max_ess` and by finite-count variance floors.
- One-sided high-confidence RNA evidence yields large `alpha_rna_add` without
  fixed edge counts.
- One-sided low-confidence RNA evidence yields small ESS.
- High-confidence global no-gDNA profile shrinks ambiguous loci toward RNA.
- High-confidence global gDNA profile does not suppress local gDNA.
- No structural gDNA candidate still zeros both grouped alphas.

### Phase 1 - Preserve Regional Variance

Modify strand and density calibration outputs:

- Add `var_count` or `sd_count` to `RegionGdnaEstimate`.
- Add per-compartment variance to `RegionGdnaChannelEstimate`.
- Add `gdna_unspliced_var` to `PriorMassDeconvolution`.
- Compute density fallback variance from Gamma-Poisson component variance plus
  state-mixture variance.
- Preserve existing `precision` for summaries but stop using it as the primary
  uncertainty surface.

### Phase 2 - Project Means And Variances

Extend `assemble_priors()` allocation to produce:

```text
gdna_expected_count
rna_expected_count
prior_unspliced_total
gdna_expected_var
prior_gdna_share_local
prior_gdna_share_var
```

Use `s_lr^2` for variance projection and inflate variance for allocation
ambiguity.

### Phase 3 - Sample Evidence Profile

Add a `PriorEvidenceProfile` dataclass:

```python
class PriorEvidenceProfile:
    global_gdna_share: float
    global_gdna_ess: float
    strand_mass_fraction: float
    density_mass_fraction: float
    precision_summary: dict[str, float]
    allocation_summary: dict[str, float]
    fallback_used: bool
```

Store it in `PriorTable.to_summary_dict()`.

### Phase 4 - Adaptive Grouped Prior Counts

Replace `_compute_grouped_prior_counts()` with a policy dispatch:

```python
if policy.mode == "adaptive":
    grouped = compute_adaptive_grouped_prior_counts(...)
elif policy.mode == "legacy_v3":
    grouped = compute_legacy_v3_grouped_prior_counts(...)
else:
    grouped = zeros
```

The adaptive function consumes projected means, projected variances,
allocation diagnostics, structural candidate availability, and the sample
evidence profile.

### Phase 5 - Output Migration

Keep existing columns for one release:

```text
alpha_gdna_add
alpha_rna_add
prior_budget
prior_gdna_share_biased
```

but redefine or supplement them with adaptive names:

```text
prior_ess_final
prior_gdna_share_final
prior_gdna_share_local
prior_gdna_share_global
```

Mark `prior_budget` and `prior_gdna_share_biased` as legacy aliases in docs.

### Phase 6 - Validation

Validation should not tune on truth by hand. Use truth only to evaluate whether
the automatic policy generalizes.

Required checks:

- Existing full pytest suite.
- v3 sentinel suite with `prior_policy=adaptive`.
- Synthetic benchmark sweeps over:
  - true gDNA fraction: none, low, high;
  - strand specificity: 0.50, 0.65, 0.90, 1.00;
  - nRNA contamination: none, moderate, high;
  - capture profiles: off-target, targeted, hybrid-capture.
- Calibration diagnostics stratified by method and confidence.
- Prior-vs-EM conflict reports inspected for systematic method failure.

Success criteria:

- No ground-truth-tuned constants are needed to make zero-gDNA sentinels pass.
- In no-gDNA samples, global `psi` is near zero and high-confidence ambiguous
  loci receive RNA-conservative priors.
- In high-gDNA samples, global `psi` increases and gDNA is not globally
  suppressed.
- Low-confidence regions produce low ESS and visible ambiguity diagnostics.
- Output summary explains why each prior was strong or weak.

## Failure Modes And Expected Behavior

### Truly Unidentifiable Single-Locus Data

If the whole sample lacks strand information, has no high-confidence global
gDNA profile, and a locus is single-exon/unspliced-compatible, no method can
determine RNA vs gDNA from the data alone.

Correct behavior:

- low prior ESS;
- explicit ambiguity diagnostics;
- no hidden hard gate;
- no claim of confident RNA or gDNA recovery.

### Strong Global No-gDNA Evidence

If the sample has many high-confidence regions showing little or no gDNA, then
ambiguous loci should borrow that sample-level information.

Correct behavior:

- global `psi` near zero;
- substantial global ESS;
- RNA-conservative priors at ambiguous structural-gDNA loci;
- no need for fixed `gdna_prior_logit_bias=-6`.

### Strong Local gDNA Evidence In Low-gDNA Sample

If local evidence is precise and contradicts the global profile, local evidence
must win.

Correct behavior:

- `lambda_l` near one;
- final share near local `p_local`;
- large prior-vs-global conflict diagnostic;
- gDNA not suppressed merely because the sample is mostly RNA.

## Recommended Immediate Next Step

Implement the uncertainty payload first. Without real count variance, adaptive
v1 will be forced to reinterpret the current `precision` scalar, which risks
creating a new magic-number layer.

The first production milestone should therefore be:

```text
Preserve posterior variance from strand and density calibration,
project it to loci, and compute Beta-equivalent ESS.
```

Once that is in place, the old v3 constants can be demoted to safety caps and
the default prior policy can become adaptive.