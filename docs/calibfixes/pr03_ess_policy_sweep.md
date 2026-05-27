# PR 3: Adaptive Prior ESS Policy Sweep

## Goal

Evaluate whether the adaptive grouped prior should have a larger effective sample size in well-identified loci, after PR 1 and PR 2 have made source mass and exposure reliable. This PR is primarily diagnostics and benchmarking. Defaults should remain conservative until the benchmark gate proves a new policy is safe across all four RNA-seq strata.

## Dependency

Do not start this PR until both are merged or available on the branch:

1. PR 1 continuous strand reliability.
2. PR 2 empirical-Bayes capture exposure.

The ESS policy should use their diagnostics. Increasing ESS before source reliability and exposure identifiability are fixed can amplify prior noise.

## Non-goals

- Do not use ESS to compensate for bad source mass.
- Do not tune a linear constant until one synthetic case looks good.
- Do not add transcript-level RNA floors.
- Do not change native EM grouped-prior redistribution unless a focused numerical bug is found.
- Do not make the default policy more aggressive without all-strata benchmark evidence.

## Files to edit

| File | Purpose of edit |
|---|---|
| `src/rigel/config.py` | Add optional ESS policy configuration, defaulting to current behavior. |
| `src/rigel/cli.py` | Add config/CLI plumbing for experiments. |
| `src/rigel/calibration/adaptive_prior.py` | Factor ESS cap calculation and support policy variants. |
| `src/rigel/calibration/prior.py` | Compute locus-level reliability/identifiability inputs and pass them to the prior policy. |
| `src/rigel/calibration/_result.py` | Add summary diagnostics for policy, cap, and reliability inputs. |
| `scripts/debug/review_prior_noise_v4.py` or new script | Add counterfactual ESS policy report without requiring a quant rerun when possible. |
| `scripts/benchmarking/runner.py` | Allow named rigel configs to sweep ESS policy values. |
| `scripts/benchmarking/analysis.py` | Add ESS-policy comparison columns to reports. |
| `tests/test_adaptive_prior.py` | Unit tests for policy formulas and flags. |
| `tests/test_per_locus_gdna_mass.py` | Integration tests for prior table diagnostics. |

## Configuration contract

Add fields to `EMConfig` or a small nested prior config if the codebase already has a better home. Keep defaults equal to current behavior.

Suggested minimal fields:

```python
adaptive_prior_ess_policy: Literal[
    "current",
    "linear",
    "information_weighted",
] = "current"
adaptive_prior_max_ess: float = 3000.0
adaptive_prior_linear_fraction: float = 0.25
```

Validation:

- policy is one of the allowed strings;
- max ESS is finite and non-negative;
- linear fraction is finite and non-negative.

CLI flags:

```text
--adaptive-prior-ess-policy {current,linear,information_weighted}
--adaptive-prior-max-ess FLOAT
--adaptive-prior-linear-fraction FLOAT
```

The old module constant `MAX_ESS = 3000.0` can remain as the default value, but `compute_adaptive_prior()` should receive the chosen value explicitly.

## New diagnostics to compute before changing policy

For each MultiLocus, compute these inputs in `prior.py`:

```text
n_competing_fragments
source_reliability_locus
exposure_identifiability_locus
information_weight_locus
ess_cap
ess_policy
```

### `n_competing_fragments`

Count fragments or EM units where gDNA and RNA genuinely compete, not all locus fragments.

Recipe:

1. Start from `locus.unit_indices`.
2. Keep units with finite `em_data.gdna_log_liks`.
3. Keep unspliced units only, matching `enable_gdna_for_multilocus()` semantics.
4. If per-unit RNA likelihood availability is exposed, require at least one finite RNA likelihood too. If not exposed, use the locus having RNA components as the RNA side and document the approximation.
5. Count fragments if a unit count/multiplicity field exists; otherwise count units and name the diagnostic `n_competing_units` until fragment multiplicity is available.

Do not count spliced-only evidence as gDNA/RNA competing evidence.

### `source_reliability_locus`

Aggregate PR 1 region/compartment reliability over regions touched by the MultiLocus.

Recipe:

```text
source_reliability_locus = weighted_mean(region_source_reliability, weights=region_unspliced_total)
```

Use `PriorMassDeconvolution.unspliced_total` for weights. If all weights are zero, reliability is zero.

### `exposure_identifiability_locus`

Aggregate PR 2 exposure fit evidence over touched regions.

Recipe:

```text
exposure_identifiability_region = O_region / (O_region + global_beta)
exposure_identifiability_locus = weighted_mean(exposure_identifiability_region, weights=region_opportunity)
```

If PR 2 stores only unit-level `O_u`, map the unit value back to each region. If there is no panel/opportunity, identifiability is zero for capture exposure, but capture-off loci should not be penalized for needing no exposure. Use:

```text
if region_opportunity == 0: exposure_identifiability_region = 1 for capture-off/no-panel neutral exposure
```

This keeps the information-weighted policy from shrinking ordinary capture-off loci solely because no capture panel exists.

### `information_weight_locus`

Use a smooth product:

```text
information_weight = sqrt(source_reliability_locus * exposure_identifiability_locus)
```

Clip to `[0, 1]`. This is a candidate policy input, not a source mass value.

## Factor ESS cap calculation

In `adaptive_prior.py`, add a helper:

```python
def compute_ess_cap(
    *,
    locus_unspliced: np.ndarray,
    n_competing_fragments: np.ndarray,
    source_reliability: np.ndarray,
    exposure_identifiability: np.ndarray,
    policy: str,
    max_ess: float,
    linear_fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return cap and uint16 policy flags."""
```

Policies:

### Current

```text
cap = min(locus_unspliced, max_ess)
```

This must reproduce current output when `max_ess == MAX_ESS`.

### Linear stress-test

```text
linear_cap = linear_fraction * n_competing_fragments
cap = min(locus_unspliced, linear_cap)
```

This is for stress testing only. It can produce very large priors in mega-loci and must not become default without strong evidence.

### Information-weighted target

Use a smooth interpolation between current and linear caps:

```text
current_cap = min(locus_unspliced, max_ess)
linear_cap = min(locus_unspliced, linear_fraction * n_competing_fragments)
info = sqrt(source_reliability * exposure_identifiability)
cap = current_cap + info * max(linear_cap - current_cap, 0)
```

This policy never goes below current behavior and only grows toward the linear stress cap when both source and exposure evidence are reliable.

Add flags:

```python
PRIOR_ESS_POLICY_LINEAR = np.uint16(0x10)
PRIOR_ESS_POLICY_INFORMATION_WEIGHTED = np.uint16(0x20)
PRIOR_ESS_INFORMATION_LIMITED = np.uint16(0x40)
```

Preserve existing `PRIOR_ESS_CAPPED` semantics when total prior mass exceeds the final cap.

## Modify `compute_adaptive_prior()`

Add optional keyword-only inputs with safe defaults:

```python
n_competing_fragments: np.ndarray | None = None
source_reliability: np.ndarray | None = None
exposure_identifiability: np.ndarray | None = None
ess_policy: str = "current"
ess_linear_fraction: float = 0.25
```

Default behavior:

- if policy is `current`, output must be bit-for-bit or tolerance-equivalent to current tests;
- if optional arrays are missing for non-current policies, raise `ValueError` rather than silently using zeros;
- cap calculation should happen where `cap = np.minimum(locus_unspliced, cap_max)` currently lives.

Extend `AdaptivePriorResult`:

```python
ess_cap: np.ndarray
n_competing_fragments: np.ndarray
source_reliability: np.ndarray
exposure_identifiability: np.ndarray
information_weight: np.ndarray
```

Extend `PriorTable` with matching diagnostics.

## Counterfactual diagnostic script

Create a debug script, for example:

```text
scripts/debug/review_prior_ess_policies.py
```

The script should read existing quant outputs if possible and emit:

```text
diagnostics/prior_ess_policy_review/
  ess_policy_summary.tsv
  ess_policy_locus_detail.tsv
  ess_policy_gene_focus.tsv
  ess_policy_review.md
```

Required columns:

```text
condition
capture_enabled
strand_specificity
gdna_truth_label
policy
locus_id
locus_unspliced
n_competing_fragments
source_reliability
exposure_identifiability
information_weight
ess_cap
prior_ess_final
alpha_gdna_add
alpha_rna_add
truth_gdna_or_source_mass_if_available
```

This script is diagnostic. It should not be the only test of behavior because changing ESS changes EM output.

## Benchmark configs

Add named configs to `scripts/benchmarking/configs/default.yaml` or a separate calibration-fix config:

```yaml
rigel_configs:
  current:
    adaptive_prior_ess_policy: current
  ess_linear_025:
    adaptive_prior_ess_policy: linear
    adaptive_prior_linear_fraction: 0.25
  ess_info_weighted_025:
    adaptive_prior_ess_policy: information_weighted
    adaptive_prior_linear_fraction: 0.25
```

If current config naming uses CLI flag names, follow the existing runner convention rather than introducing a second schema.

## Tests

### `tests/test_adaptive_prior.py`

Add tests for:

1. Current policy reproduces the old cap exactly.
2. Linear policy uses `linear_fraction * n_competing_fragments` and respects `locus_unspliced`.
3. Information-weighted policy equals current when reliability or identifiability is zero.
4. Information-weighted policy equals linear when reliability and identifiability are one and linear is above current.
5. Information-weighted policy never exceeds `locus_unspliced`.
6. Missing diagnostic arrays raise `ValueError` for non-current policies.
7. Existing flags still report capped loci correctly.

### `tests/test_per_locus_gdna_mass.py`

Add a small integration test asserting `PriorTable.to_summary_dict()` includes:

```text
prior_policy_name
prior_ess_cap summary
n_competing_fragments summary
source_reliability summary
exposure_identifiability summary
information_weight summary
```

### CLI/config tests

Add tests that config files and CLI flags set the policy fields and preserve defaults.

## Validation commands

```bash
conda activate rigel && ruff check src/rigel/config.py src/rigel/cli.py src/rigel/calibration/adaptive_prior.py src/rigel/calibration/prior.py scripts/debug/review_prior_ess_policies.py scripts/benchmarking tests/test_adaptive_prior.py tests/test_per_locus_gdna_mass.py
conda activate rigel && pytest tests/test_adaptive_prior.py tests/test_per_locus_gdna_mass.py tests/test_cli.py -v
```

Run broader calibration tests:

```bash
conda activate rigel && pytest tests/test_calibrate.py tests/test_pipeline_wiring.py tests/test_golden_output.py -v
```

## Benchmark gate

Run the eight-condition hybrid-capture suite for at least:

- current;
- linear 0.25;
- information-weighted 0.25.

Acceptance requirements for changing the default away from current:

| Scenario group | Gate |
|---|---|
| no-gDNA capture-off | no new false gDNA; mRNA MARD does not regress beyond 5 percent relative |
| no-gDNA capture-on | no exposure-driven false gDNA amplification |
| high-gDNA stranded capture-off | maintain or improve mRNA and gDNA source error |
| high-gDNA stranded capture-on | improve or maintain PR 1 + PR 2 results |
| high-gDNA unstranded capture-off | no material regression |
| high-gDNA unstranded capture-on | no material regression; full repair remains deferred |

Additional transcript-collapse gate:

- identify transcripts that newly go to zero under a stronger policy;
- for each, report whether the isoform identifiability diagnostic is unavailable, ambiguous, or well-conditioned;
- do not merge a stronger default if well-conditioned transcripts collapse without an explanatory source/truth improvement.

## Merge decision

Three valid outcomes are allowed:

1. Keep default `current`, merge only diagnostics/config plumbing.
2. Add `information_weighted` as an opt-in experimental policy, keep default `current`.
3. Make `information_weighted` default only if the full benchmark gate passes.

The linear policy should remain a stress-test policy unless there is unexpectedly strong all-strata evidence.

## Review checklist

- [ ] Current policy is backward compatible by default.
- [ ] ESS policy uses source reliability and exposure identifiability from PR 1/PR 2.
- [ ] `n_competing_fragments` does not count spliced-only or non-competing evidence.
- [ ] No transcript floor or native EM redistribution change was added.
- [ ] Benchmark configs include current, linear, and information-weighted policies.
- [ ] Default policy decision is justified by all-strata benchmark output.
