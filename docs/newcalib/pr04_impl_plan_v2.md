# PR 04 v2 - EB Exposure Factor Model

## Status

Supersedes [pr04_impl_plan_v1.md](pr04_impl_plan_v1.md).

This v2 plan simplifies PR 04 after repo audit and design review. The implementation goal is now
narrower and cleaner:

1. Produce a per-region `RegionExposure` object with an empirical-Bayes-shrunk `omega_r`.
2. Use only the discrete physical unspliced support count `N_r` as the observation support.
3. Remove the stale `A_r` / `RegionExposure.from_density()` production path.
4. Do not let `omega_r` affect locus EM or `gdna_eff_len` in PR 04. PR 05 owns downstream
   denominator normalization.

The key correction from v1 is that PR 04 must not alias `RegionCalibration.A_r` to
`RegionExposure.omega`. The current downstream code multiplies gDNA effective length by exposure,
which is the wrong sign for enriched regions. PR 04 therefore exposes `omega_r` for inspection and
future PR 05 consumption, while keeping current EM behavior exposure-neutral.

## 1. Plain-Language Model

Each region has an estimated gDNA mass `M_r` from PR 03 and a region size `O_r`. PR 03 also learns a
library-wide gDNA density `rho0`. PR 04 asks one question:

```text
Is this region sampled more or less often than expected from the library-wide gDNA density?
```

The raw regional answer is:

```text
raw_ratio_r = regional density / global density
```

`raw_ratio_r = 1` means the region matches the global density. `raw_ratio_r = 10` means the region
looks ten times more visible to gDNA sampling than expected. `raw_ratio_r = 0.1` means it looks ten
times depleted.

Raw ratios are noisy when few physical fragments support the region. PR 04 therefore shrinks each
log ratio toward zero:

```text
log_raw_ratio_r = log(raw_ratio_r)
log_omega_r = shrink_weight_r * log_raw_ratio_r
omega_r = exp(log_omega_r)
```

`omega_r` is the production exposure factor. It is centered at 1.0, may be below or above 1.0, and
is not a probability.

## 2. What Is `tau2`?

`tau2` is the library-wide variance of real regional log-exposure.

More concretely:

- `log_raw_ratio_r` is the observed regional signal.
- Some spread in `log_raw_ratio_r` is just finite-count noise.
- Some spread is real sampling heterogeneity, such as capture enrichment, mappability, or
  accessibility.
- `tau2` estimates the real heterogeneity part after subtracting the expected finite-count noise.

Small `tau2` means the library behaves like ordinary RNA-seq: most regional differences are noise,
so `omega_r` shrinks strongly toward 1.0.

Large `tau2` means the library has real regional exposure variation: high-support regions with big
raw ratios should keep much of that signal.

PR 04 estimates `tau2` from dependable, non-imputed regions only. It does not learn `tau2` from
Tier 3 fallback rows because those rows are imputed from the global background and are not
independent observations.

### `tau2` Formula

For each dependable region in the pool:

```text
y_r = log_raw_ratio_r
v_obs_r = observation noise variance for y_r
```

PR 04 uses the method-of-moments estimator:

```text
empirical_var = weighted_mean(y_r^2)
mean_v_obs    = weighted_mean(v_obs_r)
tau2_hat      = max(empirical_var - mean_v_obs, 0.0)
```

The subtraction is the important part: if the observed spread is fully explained by finite-count
noise, `tau2_hat` becomes zero and all regions shrink to `omega_r = 1`.

## 3. Core Decisions In V2

### 3.1 No `A_r` Alias

`RegionCalibration.A_r` is removed, not aliased.

PR 04 must not pass `omega_r` into the existing `prior.py` exposure projection. Current `prior.py`
multiplies `gdna_eff_len` by an exposure weight, which penalizes high-exposure regions in native EM.
That sign is rejected by the roadmap.

The principled PR 04 fix is:

```text
region_exposure.omega exists for diagnostics and PR 05.
assemble_priors() remains exposure-neutral in PR 04.
gdna_em_exposure_weight = 1.0 for every locus in PR 04.
gdna_eff_len = gdna_eff_len_unweighted in PR 04.
```

PR 05 will introduce the correct downstream usage:

```text
gdna_eff_len_em = gdna_eff_len_unweighted / omega_locus
```

### 3.2 Discrete Support Count Is The ESS

The observation support for PR 04 is the discrete physical fragment count:

```text
N_r = RegionUnsplicedMass.unspliced_counts[r]
```

There is no fractional ESS and no mass-squared ESS. Fractional overlap affects `M_r`; it does not
change how many independent physical fragments were observed.

The observation variance is therefore intentionally simple:

```text
if N_r >= 1:
    v_obs_r = 1.0 / N_r
else:
    v_obs_r = no_support_v_obs   # large finite default, e.g. 1e6
```

This is not a heuristic fallback. It is the production definition of uncertainty for PR 04.

### 3.3 Tier 3 Is Exposure-Neutral

Tier 3 rows (`METHOD_BACKGROUND_FALLBACK`) are imputed from the global background. After clipping to
observed total mass, their diagnostic raw ratio may not equal exactly 1.0. That does not make them
independent exposure observations.

Therefore PR 04 handles Tier 3 as follows:

```text
raw_ratio_r      = computed and reported for diagnostics
log_raw_ratio_r  = log(raw_ratio_r), finite because of the pseudocount
pool_weight_r    = 0.0
tau2 contribution = none
shrink_weight_r  = 0.0
omega_r          = 1.0
```

The flag marks them as imputed fallback rows. This keeps the production exposure surface from
learning from its own background imputation.

### 3.4 Bootstrap Is Exposure-Neutral

The first calibration pass starts with `BackgroundDensity.from_bootstrap(...)`. That object has not
yet been fit from the PR 03 Tier 1/2 pool. When PR 04 sees this bootstrap density, it returns an
all-ones exposure surface:

```text
omega_r = 1.0
log_omega_r = 0.0
shrink_weight_r = 0.0
tau2_method = "bootstrap_neutral"
```

No bootstrap dispersion leaks into downstream state.

### 3.5 Damping Is Simple And Optional

When a previous real `RegionExposure` exists, PR 04 may damp `tau2` linearly:

```text
tau2 = (1.0 - tau2_damping) * previous.tau2 + tau2_damping * tau2_hat
```

Default `tau2_damping = 0.5`. If there is no previous real exposure, use `tau2_hat` directly.

There is no hidden floor in shrinkage. If `tau2 == 0`, then `shrink_weight == 0` and every non-clipped
`omega_r == 1`.

A tiny positive `tau2_floor_for_reporting` may be used only when computing reciprocal diagnostic
quantities. It must not affect `shrink_weight`.

## 4. Dataclass

Replace the existing `src/rigel/calibration/exposure.py` with the PR 04 production contract.
Delete `RegionExposure.uniform()`, `RegionExposure.from_density()`, `mode`, `A_r`, `rho_r`,
`rho_ref`, `reference_quantile`, and `eligible`.

```python
@dataclass(frozen=True, slots=True)
class RegionExposure:
    """Per-region multiplicative gDNA sampling exposure.

    Produced by PR 04. Consumed by PR 05 for locus-EM denominator normalization.
    PR 04 itself does not feed omega into EM.
    """

    # Primary output.
    omega: np.ndarray              # float64[R], positive exposure factor, center = 1.0
    log_omega: np.ndarray          # float64[R], log(omega)

    # Raw and shrunk signals.
    raw_ratio: np.ndarray          # float64[R], rho_hat_r / rho0, always > 0
    log_raw_ratio: np.ndarray      # float64[R], log(raw_ratio_r)
    shrink_weight: np.ndarray      # float64[R], tau2 / (tau2 + v_obs_r), in [0, 1]
    v_obs: np.ndarray              # float64[R], observation variance on log scale

    # Anchors and fit diagnostics.
    lambda_global: np.ndarray      # float64[R], rho0 * O_r
    rho0: float                    # rho0 used for this exposure fit
    tau2: float                    # fitted library-wide real log-exposure variance
    tau2_hat: float                # undamped method-of-moments estimate

    # Support carried through for diagnostics and downstream validation.
    support_count: np.ndarray      # uint64[R], copy of unspliced_counts

    # Pool diagnostics.
    tau2_pool_size: int            # number of regions with positive tau2 pool weight
    tau2_method: str               # see allowed values below

    # Region diagnostic bits.
    flags: np.ndarray              # uint16[R]
```

Allowed `tau2_method` values:

```text
"bootstrap_neutral"   BackgroundDensity was bootstrap; omega is all ones.
"no_pool_neutral"     No dependable positive-weight pool; omega is all ones.
"moment"              tau2_hat used directly.
"moment_damped"       tau2_hat blended with previous.tau2.
```

### Invariants

`RegionExposure.__post_init__` enforces:

```text
all float arrays are float64, C-contiguous, 1D, same shape (R,)
support_count is uint64, C-contiguous, shape (R,)
flags is uint16, C-contiguous, shape (R,)
omega is finite and > 0
log_omega is finite and consistent with log(omega)
raw_ratio is finite and > 0
log_raw_ratio is finite and consistent with log(raw_ratio)
0 <= shrink_weight <= 1
v_obs is finite and >= 0
lambda_global is finite and >= 0
rho0 is finite and > 0
tau2 and tau2_hat are finite and >= 0
tau2_pool_size is int and >= 0
tau2_method is one of the allowed strings
```

## 5. Flags

Use flags only when they have a clear production or diagnostic purpose.

```python
FLAG_EXPOSURE_NO_SUPPORT          = 1 << 0  # N_r == 0; v_obs uses no_support_v_obs
FLAG_EXPOSURE_NOT_TAU2_POOL       = 1 << 1  # region did not contribute to tau2 fit
FLAG_EXPOSURE_IMPUTED_TIER3       = 1 << 2  # method_r == METHOD_BACKGROUND_FALLBACK
FLAG_EXPOSURE_NUMERIC_FLOOR       = 1 << 3  # omega clipped at omega_floor
FLAG_EXPOSURE_NUMERIC_CEILING     = 1 << 4  # omega clipped at omega_ceiling
FLAG_EXPOSURE_BOOTSTRAP_NEUTRAL   = 1 << 5  # bootstrap density caused all-ones exposure
```

No separate `NOT_UNEXPRESSED` flag is needed. Low `p_unexpressed` affects only the tau2 pool weight;
if the resulting weight is zero, `FLAG_EXPOSURE_NOT_TAU2_POOL` is sufficient.

## 6. Per-Region Raw Signal

PR 04 reuses PR 03's Bayesian-shrunk density formula with the same pseudocount convention:

```text
pseudo_mass = alpha_floor
pseudo_size = alpha_floor / rho0
rho_hat_r   = (M_r + pseudo_mass) / (O_r + pseudo_size)
raw_ratio_r = rho_hat_r / rho0
log_raw_ratio_r = log(raw_ratio_r)
lambda_global_r = rho0 * O_r
```

Inputs:

```text
M_r = RegionUnsplicedMass.gdna_mass
O_r = RegionUnsplicedMass.region_size_bp
rho0 = BackgroundDensity.rho0_mean
```

The pseudocount makes `raw_ratio_r` strictly positive even when `M_r == 0`. There is no additive
`eps_mass` hack.

For Tier 3 rows, compute and report this raw diagnostic, but do not use it for shrinkage or `tau2`.
Final Tier 3 `omega_r` is exactly 1.0.

## 7. Observation Variance

The production observation variance is based only on the physical support count:

```text
support = float(N_r)
v_obs_r = 1.0 / support          if support >= 1
v_obs_r = no_support_v_obs       if support == 0
```

Default:

```text
no_support_v_obs = 1e6
```

This keeps all arrays finite and makes zero-support rows shrink essentially completely if they ever
reach the shrinkage formula. In practice, zero-support and Tier 3 rows are not part of the `tau2`
pool.

## 8. Learning `tau2`

### 8.1 Pool

Only dependable, independently estimated rows contribute to `tau2`:

```text
pool_weight_r = precision_r * float(N_r) * p_unexpressed_r
pool_mask_r = method_r in {METHOD_STRAND, METHOD_BOUNDARY}
            & O_r >= 1.0
            & N_r >= 1
            & pool_weight_r > 0
```

This is a soft gate. There is no hard `p_unexpressed` threshold. A region with
`p_unexpressed == 0` simply receives zero pool weight.

Tier 3 rows are excluded because they are imputed.

### 8.2 Winsorization For The Fit Only

Use PR 03 `BackgroundDensity.log_dispersion` as the robust scale for winsorization inside the tau2
fit only:

```text
clip_radius = winsorize_k * max(background_density.log_dispersion, tiny_positive)
y_fit_r = clip(log_raw_ratio_r, -clip_radius, +clip_radius)
```

Default:

```text
winsorize_k = 4.0
```

The final per-region `omega_r` uses the unwinsorized `log_raw_ratio_r`. Winsorization protects the
library-wide `tau2` estimate from one extreme row; it does not cap supported capture spikes.

### 8.3 Method Of Moments

For pool rows:

```text
empirical_var = weighted_mean(y_fit_r ** 2, pool_weight_r)
mean_v_obs    = weighted_mean(v_obs_r,      pool_weight_r)
tau2_hat      = max(empirical_var - mean_v_obs, 0.0)
```

If the pool is empty, return all-ones exposure with `tau2_method = "no_pool_neutral"`.

If a previous real exposure exists, damp linearly:

```text
tau2 = (1.0 - tau2_damping) * previous.tau2 + tau2_damping * tau2_hat
```

Otherwise:

```text
tau2 = tau2_hat
```

`tau2` is a variance, so it is allowed to be exactly zero.

## 9. Shrinkage

For non-bootstrap, non-Tier-3 rows:

```text
if tau2 == 0:
    shrink_weight_r = 0.0
else:
    shrink_weight_r = tau2 / (tau2 + v_obs_r)

log_omega_r = shrink_weight_r * log_raw_ratio_r
omega_r = exp(log_omega_r)
```

For bootstrap and no-pool neutral fits:

```text
omega_r = 1.0
log_omega_r = 0.0
shrink_weight_r = 0.0
```

For Tier 3 rows in otherwise real fits:

```text
omega_r = 1.0
log_omega_r = 0.0
shrink_weight_r = 0.0
```

### Numeric Guards

Use broad finite-float guards only:

```text
omega_floor = 1e-6
omega_ceiling = 1e6
```

After clipping, recompute `log_omega = log(omega)` so the two fields stay consistent.
Set floor/ceiling flags only for rows actually clipped.

## 10. Public API

```python
def estimate_region_exposure(
    region_unspliced_mass: RegionUnsplicedMass,
    background_density: BackgroundDensity,
    p_unexpressed: np.ndarray,
    *,
    previous: RegionExposure | None = None,
    alpha_floor: float = 1.0,
    tau2_damping: float = 0.5,
    winsorize_k: float = 4.0,
    no_support_v_obs: float = 1.0e6,
    omega_floor: float = 1.0e-6,
    omega_ceiling: float = 1.0e6,
) -> RegionExposure:
    ...
```

Validation rules:

```text
alpha_floor > 0
tau2_damping in [0, 1]
winsorize_k > 0
no_support_v_obs > 0 and finite
0 < omega_floor < 1 < omega_ceiling
p_unexpressed shape == (R,), finite, clipped or validated to [0, 1]
```

Use local helpers only where they clarify implementation:

```python
_weighted_mean(values, weights) -> float
_as_float64_vector(name, values, region_count) -> np.ndarray
```

No Gamma-Poisson alternate path. No capture-mode parameter. No labels.

## 11. Calibration Loop Wiring

PR 04 does not change `calibration_m_step()`.

The correct wiring point is `run_calibration_iteration()`:

```text
current_density = BackgroundDensity.from_bootstrap(background)
previous_region_exposure = None

for pass_index in ...:
    step = calibration_e_step(
        ...,
        background_density=current_density,
        previous_region_exposure=previous_region_exposure,
    )

    previous_region_exposure = step.region_exposure

    if not done:
        current_background, current_kappa = calibration_m_step(...)
        current_density = estimate_background_density(...)
```

`calibration_e_step()` builds `RegionUnsplicedMass`, then immediately calls
`estimate_region_exposure(...)` using the `BackgroundDensity` supplied for that pass.

Bootstrap pass behavior is deterministic because the supplied `BackgroundDensity` has
`fit_status == "fallback_bootstrap"` and `n_effective_regions == 0`.

Add one diagnostic to each pass:

```text
tau2
tau2_hat
relative_tau2_shift
tau2_method
tau2_pool_size
```

Do not add tau2 to the convergence stopping criterion in PR 04.

## 12. Files To Touch

### Replace

- `src/rigel/calibration/exposure.py`
  - Replace old `A_r` / `from_density` class with the PR 04 `RegionExposure` dataclass,
    flag constants, and `estimate_region_exposure()`.

### Update

- `src/rigel/calibration/calibration_iteration.py`
  - Import `RegionExposure` and `estimate_region_exposure`.
  - Add `region_exposure: RegionExposure` to `RegionCalibration`.
  - Add `region_exposure: RegionExposure` to `CalibrationStepResult`.
  - Remove the `A_r` fields from both dataclasses.
  - Remove the inline all-ones exposure vector in `calibration_e_step()`.
  - Add `previous_region_exposure` plumbing in `run_calibration_iteration()`.
  - Keep the `rho_off` property as the bridge to `background_density.rho0_mean`.

- `src/rigel/calibration/prior.py`
  - Remove `bp_weighted_mean_exposure_over_blocks` usage.
  - Set `gdna_em_exposure_weight` to ones in PR 04.
  - Set `gdna_eff_len` equal to `gdna_eff_len_unweighted` in PR 04.
  - Leave a short comment that PR 05 will apply `region_exposure.omega` with division.

- `src/rigel/calibration/_exposure.py`
  - Delete `bp_weighted_mean_exposure_over_blocks` if no longer used.
  - Keep geometric effective-length helpers: `l_eff_contained`,
    `fractional_boundary_side_exposure`, `gdna_eff_len_for_loci`, and related helpers.

- `src/rigel/calibration/_result.py`
  - Remove `A_r` summary.
  - Add `region_exposure` summary with `rho0`, `tau2`, `tau2_hat`, `tau2_method`,
    `tau2_pool_size`, and stats for `omega`, `raw_ratio`, `shrink_weight`, `v_obs`, and
    `support_count`.

- `src/rigel/pipeline.py`
  - Replace the `[CAL] A_r ...` log with an exposure log using `region_exposure` summary.

- `docs/newcalib/README.md`
  - Update the PR 04 link to this v2 plan when PR 04 starts.

### Delete Or Rewrite Tests

- `tests/test_exposure.py`
  - Rewrite for the PR 04 `RegionExposure` contract or delete if fully superseded.

- `tests/test_region_exposure_from_density.py`
  - Delete. `RegionExposure.from_density()` is removed.

### Update Existing Tests

- `tests/test_calibration_iteration.py`
- `tests/test_calibration_result.py`
- `tests/test_calibration_prior.py`
- `tests/test_per_locus_gdna_mass.py`
- `tests/test_bayesian_prior_acceptance.py`
- `tests/test_calibrate.py`
- `tests/test_pipeline_smoke.py`
- `tests/test_pipeline_wiring.py`

Fixture updates should construct `RegionExposure` where `RegionCalibration` is built directly.

## 13. New Test Plan

Create `tests/test_region_exposure.py`.

### Dataclass

1. Dtype and contiguity: all float arrays are `float64`, support is `uint64`, flags are `uint16`.
2. Shape validation: mismatched arrays raise.
3. Numeric validation: negative/NaN/inf `omega`, `raw_ratio`, `v_obs`, `tau2`, or `rho0` raises.
4. Log consistency: `log_omega` must match `log(omega)` and `log_raw_ratio` must match
   `log(raw_ratio)`.

### Raw Signal

5. Zero mass: `M_r = 0` gives positive finite `raw_ratio` and finite `log_raw_ratio`.
6. High mass: large `M_r / O_r` gives `raw_ratio > 1`.
7. Pseudocount consistency: tiny regions shrink raw density toward `rho0` before EB shrinkage.

### Support Variance

8. `N_r = 1` gives `v_obs = 1.0`.
9. `N_r = 1000` gives `v_obs = 0.001`.
10. `N_r = 0` gives `v_obs = no_support_v_obs` and sets `FLAG_EXPOSURE_NO_SUPPORT`.
11. No fractional support: changing `M_r` with fixed `N_r` does not change `v_obs`.

### Tau2

12. Bootstrap density returns all `omega == 1`, all `shrink_weight == 0`, method
    `"bootstrap_neutral"`.
13. Empty dependable pool returns all `omega == 1`, method `"no_pool_neutral"`.
14. Uniform high-support pool learns `tau2 == 0` or near zero and shrinks to `omega ~= 1`.
15. Heterogeneous high-support pool learns positive `tau2` and keeps supported enriched regions
    above 1.
16. Tier 3 rows do not contribute to `tau2_pool_size`.
17. Adding 1000 Tier 3 rows does not change `tau2` compared with the Tier 1/2-only fit.
18. Winsorization affects `tau2_hat` but not the final per-region `log_raw_ratio` used for omega.
19. Damping blends current `tau2_hat` with previous real `tau2` by the configured linear weight.
20. `tau2 == 0` gives `shrink_weight == 0` and `omega == 1` for all non-clipped rows.

### Shrinkage And Guards

21. With the same raw ratio and positive tau2, high support has larger `shrink_weight` than low
    support.
22. Extreme high omega clips to `omega_ceiling` and sets `FLAG_EXPOSURE_NUMERIC_CEILING`.
23. Extreme low omega clips to `omega_floor` and sets `FLAG_EXPOSURE_NUMERIC_FLOOR`.
24. Tier 3 rows have `omega == 1` and `shrink_weight == 0` even when diagnostic raw ratio is not
    exactly 1 because of clipping in PR 03 imputation.

### Plumbing

25. `RegionCalibration` and `CalibrationStepResult` expose `region_exposure` and do not expose an
    `A_r` field.
26. `assemble_priors()` is exposure-neutral in PR 04: non-one `region_exposure.omega` does not
    change `gdna_eff_len`, and `gdna_em_exposure_weight` is all ones.
27. Calibration summary contains `region_exposure` and no `A_r` key.
28. Pipeline smoke produces a populated `region_exposure` summary.
29. Grep-style guard: no production code references `RegionExposure.from_density`,
    `RegionExposure.uniform`, or `region_calibration.A_r`.

## 14. Validation Commands

Targeted:

```bash
conda activate rigel && pytest \
    tests/test_region_exposure.py \
    tests/test_region_unspliced_mass.py \
    tests/test_calibration_iteration.py \
    tests/test_calibration_result.py \
    tests/test_calibration_prior.py \
    tests/test_per_locus_gdna_mass.py \
    -v
```

Smoke and wiring:

```bash
conda activate rigel && pytest \
    tests/test_pipeline_smoke.py \
    tests/test_pipeline_wiring.py \
    tests/test_calibrate.py \
    -v
```

No native rebuild is required unless unrelated native files are edited.

## 15. Implementation Order

1. Replace `src/rigel/calibration/exposure.py` and land dataclass/raw-signal/support/tau2 unit
   tests.
2. Remove stale `from_density` tests and update direct `RegionExposure` tests.
3. Wire `RegionExposure` into `calibration_iteration.py` and update calibration iteration tests.
4. Remove `A_r` from production dataclasses and summaries.
5. Make `prior.py` exposure-neutral for PR 04 and update prior/gDNA effective-length tests.
6. Update pipeline logging and smoke tests.
7. Run targeted suite, then smoke/wiring suite.
8. Do a final grep for stale contracts:

```bash
rg "RegionExposure\.from_density|RegionExposure\.uniform|region_calibration\.A_r|\bA_r\b" \
    src/rigel tests
```

Expected result after PR 04: no production `A_r`; any remaining test/golden references are either
removed or intentionally updated to `region_exposure`.

## 16. Done Means

- Local repo is clean before implementation begins and remains free of the previously observed
  `PriorMassDeconvolution` clobber.
- `RegionExposure` is the only production exposure dataclass.
- `omega_r` is produced for every region in every calibration pass.
- Bootstrap/no-pool/Tier-3 rows are exposure-neutral (`omega == 1`).
- `v_obs` depends only on discrete physical support count `N_r`.
- `tau2` is fit only from dependable Tier 1/2 rows with positive soft pool weight.
- `A_r` is removed from production dataclasses, summaries, pipeline logs, and prior assembly.
- PR 04 does not change EM exposure behavior; `gdna_eff_len` stays unweighted until PR 05.
- Old `RegionExposure.from_density()` and `RegionExposure.uniform()` tests are gone or rewritten.
- Targeted and smoke validation commands pass.

## 17. Implementation Log

Implemented PR 04 on top of the restored `origin/main` calibration loop.

### Production Changes

- Replaced the stale `RegionExposure` API in `src/rigel/calibration/exposure.py` with the PR 04
  EB exposure dataclass, exposure flags, and `estimate_region_exposure()`.
- Removed `RegionCalibration.A_r` and `CalibrationStepResult.A_r`; both now carry
  `region_exposure: RegionExposure`.
- Wired `estimate_region_exposure()` immediately after `build_region_unspliced_mass()` in the
  E-step. `run_calibration_iteration()` now carries `previous_region_exposure` for damped `tau2`.
- Added pass diagnostics for `exposure_tau2`, `exposure_tau2_hat`,
  `exposure_relative_tau2_shift`, `exposure_tau2_method`, and `exposure_tau2_pool_size`.
- Made `assemble_priors()` exposure-neutral for PR 04: `gdna_em_exposure_weight` is all ones and
  `gdna_eff_len` equals the unweighted gDNA effective length.
- Deleted the stale `bp_weighted_mean_exposure_over_blocks()` helper from `_exposure.py`.
- Replaced the `A_r` calibration summary and pipeline log with nested `region_exposure`
  diagnostics.
- Updated `docs/newcalib/README.md` so the PR 04 roadmap link points to this v2 plan.

### Tests Updated

- Rewrote `tests/test_exposure.py` around the PR 04 estimator contract.
- Deleted `tests/test_region_exposure_from_density.py` because `RegionExposure.from_density()` no
  longer exists.
- Updated direct `RegionCalibration(...)` fixtures in calibration result, prior, per-locus prior,
  Bayesian prior acceptance, and calibrate wiring tests.
- Updated prior tests to assert PR 04 exposure neutrality even when fixture `omega` is non-one.

### Deviations From This Plan

- The new exposure tests live in the existing `tests/test_exposure.py` instead of a new
  `tests/test_region_exposure.py`. This keeps the previous exposure test location and avoids adding
  another near-duplicate test module.
- `estimate_region_exposure()` mirrors the PR 03 method integer values locally to avoid a runtime
  circular import with `calibration_iteration.py`. The public method constants remain owned by
  `calibration_iteration.py`.

### Validation Run

```bash
conda activate rigel && ruff check \
    src/rigel/calibration/exposure.py \
    src/rigel/calibration/calibration_iteration.py \
    src/rigel/calibration/prior.py \
    src/rigel/calibration/_exposure.py \
    src/rigel/calibration/_result.py \
    src/rigel/pipeline.py \
    tests/test_exposure.py \
    tests/test_calibration_iteration.py \
    tests/test_calibration_result.py \
    tests/test_calibration_prior.py \
    tests/test_per_locus_gdna_mass.py \
    tests/test_bayesian_prior_acceptance.py \
    tests/test_calibrate.py
```

Result: passed.

```bash
conda activate rigel && pytest \
    tests/test_exposure.py \
    tests/test_calibration_iteration.py \
    tests/test_calibration_result.py \
    tests/test_calibration_prior.py \
    tests/test_per_locus_gdna_mass.py \
    tests/test_bayesian_prior_acceptance.py \
    tests/test_calibrate.py \
    -v
```

Result: `31 passed`.

```bash
conda activate rigel && pytest tests/test_golden_output.py tests/test_pipeline_wiring.py -v
```

Result: `24 passed`.

```bash
conda activate rigel && pytest tests/ -q
```

Result: `1200 passed`.
