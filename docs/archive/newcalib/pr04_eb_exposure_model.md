# PR 04 - EB Exposure Factor Model

## Goal

Implement a per-region exposure factor `omega_r` that shrinks noisy regional density ratios toward
the global center while allowing strongly supported capture spikes to remain large.

This PR should introduce the production exposure model. It should not depend on capture labels.

## Inputs

From PR 03:

```text
M_r = regional gDNA mass
O_r = regional fractional-mass opportunity
N_r = physical unspliced fragment support
ESS_r = support effective sample size, preferably M_r^2 / sum_i w_ir^2
p_unexpressed_r = posterior probability the region is not expression-contaminated
rho0 = global gDNA density center
```

Define the expected mass at global exposure:

```text
lambda_r = rho0 * O_r
```

## Model

Use variance-aware EB shrinkage on the log exposure ratio. This is deliberately chosen over a pure
Gamma-Poisson update because `M_r` is fractional mass while `N_r` and `ESS_r` carry the missing
support information.

Raw regional signal:

```text
y_r = log((M_r + eps_mass) / (lambda_r + eps_mass))
```

Observation variance:

```text
v_obs_r = v_mass_r + v_support_r
v_mass_r = 1 / max(M_r, mass_floor)
v_support_r = c_support / max(ESS_r, 1)
```

`c_support` starts as a conservative constant and should be reported in summaries. If benchmarks show
over-shrinkage or under-shrinkage, learn it later from replicate simulations. Do not add a capture
mode switch.

Learn the library-wide exposure variance from clean regions:

```text
eligible_r = p_unexpressed_r high enough and O_r > 0 and ESS_r > 0
tau2 = robust_weighted_var(y_r over eligible regions) - median(v_obs_r over eligible regions)
tau2 = max(tau2, tau2_floor)
```

Use `p_unexpressed_r * min(ESS_r, ESS_cap)` as the fitting weight. Robustly winsorize or trim `y_r`
only while estimating `tau2`; do not winsorize the final regional `omega_r` except for numeric
finite guards.

Posterior shrinkage:

```text
shrink_weight_r = tau2 / (tau2 + v_obs_r)
log_omega_r = shrink_weight_r * y_r
omega_r = exp(log_omega_r)
```

The center is exactly `log(1.0) = 0.0`.

## Why This Uses Support

Two regions can have the same fractional mass and the same raw density ratio but very different
evidence:

- one fragment contributes 20 percent overlap,
- 200 fragments each contribute 0.1 percent overlap.

The raw fractional mass alone cannot distinguish those cases. `ESS_r` and `N_r` give the shrinkage
model the missing uncertainty scale.

## Relationship to Gamma-Poisson

The Gamma-Poisson posterior mean is a useful diagnostic:

```text
omega_gp_r = (nu + M_r) / (nu + lambda_r)
nu = 1 / tau2
```

Do not make this a second production path in the first implementation. It can be reported in tests or
debug output while validating the log-EB model. If the log model fails benchmark acceptance, replace
it cleanly rather than maintaining two exposure systems.

## New Dataclass

Add a production dataclass, likely in `src/rigel/calibration/exposure.py`:

```text
RegionExposure
    omega: float32[R]
    raw_ratio: float32[R]
    log_raw_ratio: float32[R]
    shrink_weight: float32[R]
    lambda_global: float32[R]
    support_count: uint32/uint64[R]
    support_ess: float64[R]
    tau2: float
    rho0: float
    flags: uint16[R]
```

Rename old `A_r` language to `omega` in new production code. If temporary aliases are needed during
the PR, keep them local and delete them before the PR is complete.

## Flags

Suggested flags:

```text
FLAG_EXPOSURE_LOW_OPPORTUNITY
FLAG_EXPOSURE_LOW_SUPPORT
FLAG_EXPOSURE_NOT_UNEXPRESSED_ELIGIBLE
FLAG_EXPOSURE_NUMERIC_FLOOR
FLAG_EXPOSURE_NUMERIC_CEILING
```

Numeric ceiling should be high enough not to clip realistic capture spikes. It is a finite-float
guard, not a modeling cap.

## Summaries

Summary JSON should include:

- `rho0`
- `tau2`
- effective `nu = 1 / tau2`
- `omega` min/p50/p95/p99/max/mean
- `raw_ratio` min/p50/p95/p99/max/mean
- `shrink_weight` min/p50/p95/p99/max/mean
- number of low-support regions
- number of exposure-training regions

## Tests

- Uniform synthetic ratios with low variance learn small `tau2` and keep `omega_r` near 1.0.
- Capture-like skew learns large `tau2` and allows high-support enriched regions to remain large.
- Low-support enriched region shrinks more strongly than high-support enriched region with the same
  raw ratio.
- Zero-mass region gives finite low `omega_r`, not NaN or infinity.
- High raw ratio with high support is not clipped to 1.0.

## Done Means

- `RegionCalibration` exposes `RegionExposure` or equivalent production field.
- Production code no longer uses `A_r` as a capture-derived value.
- The exposure model can be run on ordinary and capture libraries without choosing a capture mode.
