# PR 04 v1 - EB Exposure Factor Model (Implementation Ready)

## Status

Supersedes [pr04_eb_exposure_model.md](pr04_eb_exposure_model.md). Implementation-ready
blueprint after independent critique. Builds directly on PR 01 (two-state filter), PR 02a
(`region_unspliced_support`), and PR 03 v3 (`RegionUnsplicedMass` + `BackgroundDensity`).

The original PR 04 draft sketched the log-EB shrinkage idea correctly but left four
implementation traps unresolved:

1. **`eps_mass` was a hack for `log(M_r / lambda_r)` at zero mass.** PR 03 already solved
   this with the Bayesian-shrunk `rho_hat_r` (Section 6.3.1). PR 04 must reuse that same
   shrunk per-region density, not introduce a parallel additive epsilon.
2. **`v_obs_r = v_mass_r + v_support_r` mixed two heuristic constants** (a mass-floor and
   `c_support`) where a Gamma-posterior delta-method variance is the principled closed form
   and uses `ESS_r` cleanly.
3. **`tau2` and PR 03's `log_dispersion` are estimating the same multiplicative
   spread.** PR 04 must seed `tau2` from PR 03 (no double estimation) and refine; otherwise
   the two stages disagree on iteration 1.
4. **The Gamma-Poisson "diagnostic" path duplicated production state and invited drift.**
   Per the README, PR 04 ships **one** production exposure model. The GP form is referenced
   for future replacement only; it is not a maintained parallel path.

Additional critique fixes folded in:

- Eligibility is a soft weight (`p_unexpressed * support`), never a hard mask, matching
  the PR 01 / PR 03 contract.
- Numeric floor/ceiling are loose finite-float guards, not modeling caps. Capture spikes
  (omega ~ 1000x) survive.
- Iteration 1 has a deterministic story: PR 03's bootstrap `BackgroundDensity` carries a
  wide `log_dispersion` and `n_effective_regions == 0`, so PR 04 returns `omega_r == 1.0`
  identically until at least one Tier 1/2 refit has happened.
- The PR 04 output schema is what PR 05 consumes; this PR does **not** aggregate to loci.
- `RegionCalibration.A_r` is retired in this PR; the replacement is `RegionExposure.omega`.

## 1. Goals and Invariants

PR 04 introduces one production contract: a per-region multiplicative sampling exposure
factor with empirical-Bayes shrinkage toward a global center of 1.0.

```text
omega_r = regional gDNA sampling density / global gDNA sampling density
log_omega_r = E[log rho_r | M_r, O_r] - log(rho0)
log_omega_r := shrink(log_omega_r, weight=tau2 / (tau2 + v_obs_r))
omega_r := exp(log_omega_r)
```

Consumed by PR 05 (locus-EM gDNA denominator normalization). Not consumed by anything in
PR 04 itself.

### Invariants

1. **Center is exactly 1.0.** `omega_r == 1.0` iff `log_omega_r == 0.0`. The shrinkage
   target is `log(1.0) = 0.0` by construction; no separate "centering" pass.
2. **Finite and positive.** Every `omega_r` is finite, > 0, and within
   `[omega_floor, omega_ceiling]`. Floor/ceiling are wide finite guards
   (`1e-6`, `1e6`), not modeling caps.
3. **Capture spikes survive.** A high-support region with `raw_ratio_r = 1000` and
   `support >> tau2 / v_obs_r threshold` retains an `omega_r` near 1000, not 1.0.
4. **Low-support regions shrink hard.** A region with `support = 1` shrinks
   `log_omega_r` by `>= 0.9` of the raw signal in any realistic library (assumes
   `tau2 << v_obs_r` at that support).
5. **No capture latent state.** No `p_captured`, no capture-mode switch, no labels. The
   model is identical for ordinary and capture libraries; the data picks `tau2`.
6. **PR 03 reuse.** PR 04 consumes `RegionUnsplicedMass` and `BackgroundDensity` as
   produced by PR 03 and does not recompute `rho0` or `M_r`.
7. **No locus aggregation.** PR 04 emits a per-region `RegionExposure`. PR 05 owns the
   FL-window-aware projection to loci.

### Non-Goals

- Locus-EM denominator normalization (PR 05).
- A second per-unit `+log(omega)` likelihood term (explicitly rejected in the README
  exposure contract; would need to replace the denominator-only model, not stack on it).
- Capture-mode switches or labels.
- Feeding `omega_r` back into PR 03's `rho0` refit (avoid coupling; `rho0` is the anchor).
- Maintaining a Gamma-Poisson production path in parallel.

## 2. Notation

| Symbol | Source | Meaning |
|--------|--------|---------|
| `M_r` | `RegionUnsplicedMass.gdna_mass` | Per-region gDNA fractional mass (PR 03) |
| `O_r` | `RegionUnsplicedMass.region_size_bp` | Region width in bp (PR 03 canonical denominator) |
| `N_r` | `RegionUnsplicedMass.unspliced_counts` | Integer physical fragment ESS (PR 02a) |
| `precision_r` | `RegionUnsplicedMass.precision` | Reliability of the `M_r` estimator (PR 03) |
| `method_r` | `RegionUnsplicedMass.method` | Tier label (PR 03) |
| `p_unexpressed_r` | `RegionCalibration.p_unexpressed` | Soft unexpressed gate (PR 01) |
| `rho0` | `BackgroundDensity.rho0_mean` | Library-wide gDNA fragment mass per bp (PR 03) |
| `alpha0`, `beta0` | `BackgroundDensity` | Gamma conjugacy view of `rho0` (PR 03) |
| `log_dispersion_prev` | `BackgroundDensity.log_dispersion` | PR 03's robust log-MAD seed for `tau2` |
| `alpha_floor`, `beta_floor` | constants | Same values used by PR 03 (`1.0`, `1.0 / rho0_prev`) |

## 3. Dataclass

New file `src/rigel/calibration/exposure.py`:

```python
@dataclass(frozen=True, slots=True)
class RegionExposure:
    """Per-region multiplicative gDNA sampling exposure.

    Produced by PR 04. Consumed by PR 05 for locus-EM denominator normalization.
    """

    # Primary output.
    omega: np.ndarray              # float64[R]   exposure factor; center = 1.0
    log_omega: np.ndarray          # float64[R]   log(omega); for downstream log-domain math

    # Raw and shrunk signals (diagnostics).
    raw_ratio: np.ndarray          # float64[R]   rho_hat_r / rho0   (always > 0)
    log_raw_ratio: np.ndarray      # float64[R]   log(raw_ratio_r)
    shrink_weight: np.ndarray      # float64[R]   tau2 / (tau2 + v_obs_r)  in [0, 1]
    v_obs: np.ndarray              # float64[R]   per-region observation variance on log scale

    # Anchors (copies from the BackgroundDensity used).
    lambda_global: np.ndarray      # float64[R]   rho0 * O_r  (expected mass at global exposure)
    rho0: float                    # the rho0 used to compute raw_ratio
    tau2: float                    # library-wide exposure variance (nats^2)
    nu_effective: float            # 1.0 / tau2  (Gamma-Poisson-equivalent prior strength)

    # Support carried through for downstream diagnostics.
    support_count: np.ndarray      # uint64[R]    copy of unspliced_counts

    # Pool / fit diagnostics.
    tau2_pool_size: int            # number of regions that contributed to tau2 fit
    tau2_method: str               # "moment" | "seeded_from_prev" | "floor"

    # Region-level diagnostic bits.
    flags: np.ndarray              # uint16[R]
```

### Flag constants

```python
FLAG_EXPOSURE_LOW_SUPPORT          = 1 << 0  # N_r < tau2_pool_min_support
FLAG_EXPOSURE_NOT_UNEXPRESSED      = 1 << 1  # p_unexpressed_r below soft-pool threshold
FLAG_EXPOSURE_TIER3_FALLBACK       = 1 << 2  # method_r == METHOD_BACKGROUND_FALLBACK; raw_ratio == 1 by construction
FLAG_EXPOSURE_NUMERIC_FLOOR        = 1 << 3  # omega clipped at omega_floor (1e-6)
FLAG_EXPOSURE_NUMERIC_CEILING      = 1 << 4  # omega clipped at omega_ceiling (1e6)
FLAG_EXPOSURE_BOOTSTRAP_ITERATION  = 1 << 5  # tau2 fit deferred (BackgroundDensity is bootstrap); omega == 1
```

### Invariants enforced in `__post_init__`

```text
all primary float arrays are float64, C-contiguous, shape (R,)
support_count.dtype == np.uint64
omega > 0 and finite
log_omega is finite
omega_floor <= omega <= omega_ceiling
0.0 <= shrink_weight <= 1.0
v_obs >= 0 and finite
tau2 >= 0 and finite
nu_effective == 1.0 / max(tau2, tau2_floor)   # exact reciprocal
tau2_method in {"moment", "seeded_from_prev", "floor"}
```

## 4. The Per-Region Signal

### 4.1 Shrunk per-region density (reuse PR 03)

PR 03 already supplies a Bayesian-shrunk `rho_hat_r` formula. PR 04 reuses it verbatim
with the **same** `alpha_floor` so the two stages stay statistically consistent:

```text
pseudo_mass = alpha_floor                            # default 1.0
pseudo_size = alpha_floor / rho0                     # implied prior-equivalent bp

rho_hat_r   = (M_r + pseudo_mass) / (O_r + pseudo_size)
raw_ratio_r = rho_hat_r / rho0
log_raw_ratio_r = log(raw_ratio_r)
```

Properties (proved in PR 03 Section 6.3.1):

- Always finite and positive. No `eps_mass` epsilon hack.
- For tiny regions (`O_r << pseudo_size`), `rho_hat_r -> rho0`, so `raw_ratio_r -> 1` and
  the region contributes nothing to the shrinkage update.
- For large regions, `rho_hat_r -> M_r / O_r`, the data-only density.

For Tier 3 regions (`method_r == METHOD_BACKGROUND_FALLBACK`), `M_r = rho0 * O_r` by
construction, so `raw_ratio_r = 1.0` exactly (modulo the alpha_floor pseudocount, which
keeps it within `[1 - O(1/O_r), 1]`). Tier 3 regions therefore contribute negligibly to
any shrinkage update regardless of weighting; the `FLAG_EXPOSURE_TIER3_FALLBACK` is
informational only.

### 4.2 Observation variance (Gamma delta-method)

Replace the original draft's `v_mass + v_support` two-term heuristic with the closed-form
delta-method variance of `log(rho_hat_r)` under a Gamma posterior:

```text
alpha_post_r = alpha_floor + M_r
beta_post_r  = beta_floor  + O_r
# rho_hat_r = alpha_post_r / beta_post_r   (matches Section 4.1 with beta_floor = pseudo_size)

# Gamma delta-method: var(log X) ~= 1 / alpha for X ~ Gamma(alpha, beta).
# Quasi-Poisson scaling by support: when M_r aggregates fractional contributions
# from N_r physical fragments, multiply alpha_post_r by an effective scale before
# inverting. The natural scale is N_r itself (fragment count is the ESS).
v_obs_r = 1.0 / max(alpha_post_r * support_scale_r, alpha_floor)
support_scale_r = max(N_r, 1.0) / max(M_r + alpha_floor, alpha_floor)
                                                # mean per-fragment mass contribution
```

Why this works:

- Pure Gamma(alpha) delta-method `var(log) ~= 1/alpha` is correct when `M_r` is itself
  a fragment count. Here `M_r` is fractional mass; `support_scale_r` rescales alpha so
  that the variance reflects the integer number of independent fragments behind it, not
  the (smaller) fractional sum.
- For a region with `N_r = 1` and a tiny fractional overlap, `alpha_post_r * support_scale_r
  -> 1`, so `v_obs_r -> 1.0` (nat^2), reflecting one observation. Equivalent to a
  ~63% relative uncertainty on `rho_hat_r`.
- For a region with `N_r = 1000` heavily-contributing fragments, `v_obs_r -> 1 / 1000 =
  0.001`, reflecting tight uncertainty.
- The mass-floor disappears: the only floor is `alpha_floor` in the denominator, which is
  the same constant used everywhere in PR 03.

`v_obs_r` is reported in `RegionExposure.v_obs` for inspection.

**Rejected alternatives.** A pure `1.0 / N_r` was rejected: it ignores the relative
contribution of fractional mass per fragment, so a fragment contributing 1% overlap is
treated identically to one contributing 100%, biasing the variance low for marginal
overlaps. A `c_support` global constant was rejected: it is an unidentified hyperparameter
that the delta-method formula determines from first principles.

## 5. Learning `tau2`

### 5.1 Pool and weights

The same soft weight as PR 03:

```text
w_r = precision_r * float(N_r) * p_unexpressed_r
pool_mask_r = (method_r in {METHOD_STRAND, METHOD_BOUNDARY})
            & (O_r >= 1.0)
            & (N_r >= 1)
            & (w_r > 0)
```

Tier 3 regions are excluded from the `tau2` pool because `raw_ratio_r = 1` by
construction (Section 4.1) and would deflate the dispersion. This mirrors PR 03's
Tier 3 exclusion from `rho0` for the same reason.

**No hard threshold on `p_unexpressed_r`.** Soft weight only.

### 5.2 Method-of-moments `tau2`

The log-EB random-effects model assumes:

```text
log_raw_ratio_r ~ N(0, tau2 + v_obs_r)
```

(Centered at 0 because PR 03 already pinned `rho0` so that the pool-weighted
`mean(log_raw_ratio_r)` is approximately 0.)

Weighted method-of-moments estimator:

```text
empirical_var = sum_r w_r * (log_raw_ratio_r)^2 / sum_r w_r
mean_v_obs    = sum_r w_r * v_obs_r            / sum_r w_r
tau2_hat      = max(empirical_var - mean_v_obs, 0.0)
```

This is the standard random-effects MoM estimator (DerSimonian-Laird style). It can be
exactly zero when the observation variance already explains the empirical spread, in
which case shrinkage is total and every `omega_r -> 1`. That is the correct behavior for
an ordinary RNA-seq library with no real exposure heterogeneity.

### 5.3 Seeding from PR 03 and damping

PR 03 already produces `log_dispersion` (weighted MAD-derived sigma in nats). On
iteration 1 (when PR 03 has just done its first real fit), seed `tau2` from this:

```text
tau2_seed = log_dispersion_prev ** 2          # convert sigma -> variance
```

If the pool is empty or `BackgroundDensity.fit_status == "fallback_bootstrap"`, do not
fit `tau2_hat` at all; emit `omega_r = 1` everywhere and set
`FLAG_EXPOSURE_BOOTSTRAP_ITERATION` on every region. This guarantees iteration 1 is a
no-op for downstream EM and avoids feeding a spurious dispersion into PR 05 before any
data has spoken.

For subsequent iterations, damp `tau2` in log space to avoid oscillation:

```text
log_tau2_next = (1 - damping) * log(max(prev.tau2, tau2_floor))
              + damping       * log(max(tau2_hat,   tau2_floor))
tau2 = exp(log_tau2_next)
```

Default `damping = 0.5` (same as PR 03). `tau2_floor = (log(1.01))**2 ~= 9.9e-5` keeps
`nu_effective` finite without forcing measurable shrinkage at the lower end.

`tau2_method`:

- `"floor"`: pool empty, `tau2 = log_dispersion_prev**2` carried unchanged (or
  `tau2_floor` if even that is unavailable).
- `"seeded_from_prev"`: bootstrap iteration; `tau2 = log_dispersion_prev**2`, no MoM fit.
- `"moment"`: pool nonempty; MoM fit ran.

### 5.4 Robust trimming during the `tau2` fit only

To prevent a single mega-enriched region from inflating `tau2` past the bulk dispersion,
winsorize `log_raw_ratio_r` to `+/- k * log_dispersion_prev` (default `k = 4`) **only
inside the MoM sum**. The per-region `log_omega_r` (Section 6) uses the **unwinsorized**
value: capture spikes must survive in the final `omega_r`.

This is the only place trimming appears. There is no trimming on `omega_r`, no top-tail
cap, no per-region clip beyond the finite-float guards in Section 6.4.

## 6. Per-Region Shrinkage

### 6.1 Posterior

```text
shrink_weight_r = tau2 / (tau2 + v_obs_r)
log_omega_r     = shrink_weight_r * log_raw_ratio_r
omega_r         = exp(log_omega_r)
```

This is the standard normal-normal posterior mean when the prior is `N(0, tau2)` and the
observation is `N(log_raw_ratio_r, v_obs_r)`. Center is exactly 0 in log space, so
`omega_r = 1` is the no-shrinkage limit (`tau2 = 0`) and the no-data limit
(`v_obs_r -> infinity`).

### 6.2 Tier 3 shortcut

For Tier 3 regions, `log_raw_ratio_r ~ 0` by construction, so the formula already
returns `omega_r ~ 1`. No special case required. The flag is set for transparency.

### 6.3 Bootstrap iteration shortcut

If `FLAG_EXPOSURE_BOOTSTRAP_ITERATION` is set globally (Section 5.3), bypass the
formula and write `omega_r = 1.0`, `log_omega_r = 0.0`, `shrink_weight_r = 0.0`. This
prevents any bootstrap-derived `log_dispersion` from leaking into omega before a real
data-driven fit exists.

### 6.4 Finite-float guards

```text
omega_floor   = 1e-6
omega_ceiling = 1e6

omega_r := clip(omega_r, omega_floor, omega_ceiling)
log_omega_r := log(omega_r)                  # recompute to keep the two consistent

if omega_r == omega_floor:   flags_r |= FLAG_EXPOSURE_NUMERIC_FLOOR
if omega_r == omega_ceiling: flags_r |= FLAG_EXPOSURE_NUMERIC_CEILING
```

These bounds are 6 orders of magnitude on each side and exist only to keep float64 math
finite. A capture spike at omega = 5000 passes through untouched. If a real benchmark
ever needs more, raise the ceiling rather than introducing a modeling cap.

## 7. Coupling to Outer Loop

PR 04 does **not** feed back into PR 03. The data flow each outer iteration is:

```text
PR 03:  build_region_unspliced_mass   ->   RegionUnsplicedMass
        estimate_background_density    ->   BackgroundDensity (rho0, alpha0, beta0, log_dispersion)
                                           |
PR 04:  estimate_region_exposure       <---+
                                           |
                                           V
                                       RegionExposure (omega_r)
                                           |
PR 05:  project to loci, normalize gDNA denominator
        (no PR 04 change of state)
```

The outer `run_calibration_iteration()` already monitors `max_state_shift` and
`relative_rho_shift` from PR 03 for convergence. Add one summary diagnostic:
`relative_tau2_shift = |log(tau2_next) - log(tau2_prev)| / max(|log(tau2_prev)|, 1.0)`.
Do **not** add it to the convergence-stopping criterion in this PR; it is observational
only until PR 06 calibrates the tolerance.

## 8. Public API

New module `src/rigel/calibration/exposure.py`:

```python
def estimate_region_exposure(
    region_unspliced_mass: RegionUnsplicedMass,
    background_density: BackgroundDensity,
    p_unexpressed: np.ndarray,                # float32 or float64[R]
    *,
    previous: RegionExposure | None = None,
    alpha_floor: float = 1.0,
    beta_floor_strategy: str = "match_pr03",  # beta = alpha_floor / rho0
    tau2_damping: float = 0.5,
    tau2_floor: float = (np.log(1.01)) ** 2,
    winsorize_k: float = 4.0,
    omega_floor: float = 1e-6,
    omega_ceiling: float = 1e6,
) -> RegionExposure:
    ...
```

Single function, single output. No capture-mode parameter, no Gamma-Poisson alternate
return.

## 9. Files to Touch

### New
- `src/rigel/calibration/exposure.py` - `RegionExposure`, flag constants,
  `estimate_region_exposure()`, and `_weighted_mean()` helper local to the module.

### `src/rigel/calibration/calibration_iteration.py`
- Add `region_exposure: RegionExposure | None` to `RegionCalibration` and
  `CalibrationStepResult` (parallel to the PR 03 bridge fields).
- In `calibration_e_step()`: after `build_region_unspliced_mass()` runs and a real
  `BackgroundDensity` exists, call `estimate_region_exposure(...)` and stash the result
  on `CalibrationStepResult.region_exposure`.
- In `calibration_m_step()`: pass the previous iteration's `RegionExposure` (if any) as
  `previous=` so log-tau2 damping has a target.
- **Retire `A_r`.** Replace `RegionCalibration.A_r` and `CalibrationStepResult.A_r` with
  `region_exposure.omega.astype(np.float32)` for any consumer that has not yet migrated.
  A temporary `@property A_r -> region_exposure.omega.astype(np.float32)` is acceptable
  for the duration of this PR; the field itself is gone before merge. Updates to the
  bootstrap iteration ones-vector handoff: emit `omega_r = 1.0` from
  `estimate_region_exposure` instead of inlining the ones-vector.

### `src/rigel/calibration/_result.py`
- Add a `RegionExposureSummary` block to the calibration summary: `rho0`, `tau2`,
  `nu_effective`, `tau2_method`, `tau2_pool_size`, and percentiles of `omega`,
  `raw_ratio`, `shrink_weight`, `v_obs` (`p50, p95, p99, max, min`).
- Drop `A_r` summary fields.

### `src/rigel/calibration/prior.py`
- No change in PR 04. (PR 05 is the consumer of `omega`.) Confirm `prior.py` does not
  read `RegionCalibration.A_r` directly; if it does, route through the temporary
  property or update to read `region_exposure.omega`.

### `src/rigel/calibration/adaptive_prior.py`
- No change in PR 04. ESS plumbing landed in PR 03.

## 10. Test Plan

New `tests/test_region_exposure.py`.

### Sanity / dataclass
1. **dtype enforcement.** All float arrays in `RegionExposure` are `float64`,
   `support_count` is `uint64`, scalars are Python `float`/`int`/`str`.
2. **Invariant checks.** Construction with `omega < 0`, NaN, mismatched shapes, or
   `nu_effective != 1/tau2` raises.

### Per-region signal
3. **Tier 3 yields `raw_ratio == 1` modulo pseudocount.** A pool of pure-Tier 3 regions
   produces `raw_ratio_r` within `[1 - 1/O_r, 1]` and `omega_r == 1` exactly when the
   bootstrap shortcut fires; close to 1 otherwise.
4. **Zero-mass region.** `M_r = 0`, `O_r = 10000`: `raw_ratio_r` is small and positive
   (no `log(0)`); `v_obs_r` reflects support `N_r`.
5. **Variance scales inversely with support.** Two regions with identical `M_r`, `O_r`
   but `N_r = 1` vs `N_r = 1000`: `v_obs_r(N=1000) < v_obs_r(N=1) / 100`.

### `tau2` estimation
6. **Uniform pool learns small `tau2`.** Synthetic pool of 200 regions with
   `log_raw_ratio_r` drawn from `N(0, 0.01)` and `v_obs_r ~ 0.01`: MoM `tau2 -> 0` and
   most `omega_r in [0.95, 1.05]`. Locks the ordinary-RNA-seq behavior.
7. **Capture-like skew learns large `tau2`.** Pool with 90% near zero and 10% with
   `log_raw_ratio_r in [3, 7]` (omega 20-1000), all high support: `tau2 > 1.0`,
   high-support enriched regions retain `omega > 10`. Locks the capture behavior.
8. **Winsorization affects only `tau2`, not `omega_r`.** A single region with
   `log_raw_ratio_r = 20` (omega ~= 5e8) gets winsorized in the MoM sum but produces a
   final `omega_r` clipped only by `omega_ceiling = 1e6`, with
   `FLAG_EXPOSURE_NUMERIC_CEILING` set.
9. **Bootstrap iteration is a no-op.** When `background_density.fit_status ==
   "fallback_bootstrap"`, every `omega_r == 1.0`, every flag carries
   `FLAG_EXPOSURE_BOOTSTRAP_ITERATION`, `tau2_method == "seeded_from_prev"`.
10. **Damping caps tau2 jumps.** Inject a 10x `tau2_hat` shock vs `previous.tau2`;
    confirm `log(tau2_next) - log(tau2_prev) ~ damping * log(10)`.

### Shrinkage
11. **Low support shrinks harder than high support.** Two regions with identical
    `log_raw_ratio_r = 2.0`, one with `N_r = 1` and one with `N_r = 500`: the
    high-support `omega_r` is closer to `exp(2)`, the low-support one closer to `1`.
12. **Zero `tau2` shrinks everything to 1.** Force `tau2_hat = 0`: every
    `omega_r == 1.0` regardless of `log_raw_ratio_r`.
13. **No tier 3 leakage into tau2.** Pool with 1000 Tier 3 regions and 10 Tier 1
    regions: `tau2_pool_size == 10`, `tau2` matches the Tier-1-only fit within numerical
    tolerance.

### Plumbing
14. **No `A_r` left in production.** Grep-style assertion: `RegionCalibration.A_r`
    field is removed or aliased to `region_exposure.omega.astype(np.float32)`.
15. **`p_unexpressed` is a soft weight.** A region with `p_unexpressed_r = 0.0` has
    zero pool weight (does not contribute to `tau2`) but still receives a per-region
    `omega_r` computed from its `log_raw_ratio_r` and `v_obs_r`. The directional gate
    affects only the dispersion fit, not per-region shrinkage.
16. **Numeric guards fire on extreme inputs.** `omega_r` outside `[1e-6, 1e6]` is
    clipped and the corresponding flag is set.

### Update existing test files
- `tests/test_calibration_iteration.py`: add `region_exposure` field to fixtures.
- `tests/test_calibration_result.py`: add the exposure summary block.
- `tests/test_pipeline_smoke.py`: confirm pipeline still runs end-to-end with the new
  field populated.

## 11. Targeted Validation Commands

```bash
conda activate rigel && pytest \
    tests/test_region_exposure.py \
    tests/test_region_unspliced_mass.py \
    tests/test_calibration_iteration.py \
    tests/test_calibration_result.py \
    tests/test_calibration_prior.py \
    -v
```

After targeted tests pass:

```bash
conda activate rigel && pytest tests/test_pipeline_smoke.py tests/test_pipeline_wiring.py -v
```

No native rebuild required.

## 12. Recommended Implementation Order

1. Land `RegionExposure` dataclass with `__post_init__` validation and tests 1, 2.
2. Land Section 4 (`rho_hat_r`, `raw_ratio_r`, `v_obs_r`) as pure functions and tests
   3, 4, 5.
3. Land Section 5 (`tau2` MoM fit, seeding, damping, winsorization) and tests 6, 7, 8,
   9, 10, 13.
4. Land Section 6 shrinkage and `estimate_region_exposure()` and tests 11, 12, 16.
5. Wire into `calibration_e_step()` / `calibration_m_step()`; retire `A_r`; land test
   14, 15. Update mechanical fixtures in existing test files.
6. Update `_result.py` summaries. Run targeted suite then smoke suite.

## 13. Resolved Critiques (Final Checklist)

- [x] **Reused PR 03's Bayesian-shrunk density.** No `eps_mass` epsilon hack;
      `raw_ratio_r > 0` by construction (Section 4.1).
- [x] **Closed-form `v_obs_r`.** Gamma delta-method on `log(rho_hat_r)` with
      quasi-Poisson scaling by `N_r`. No `c_support` hyperparameter; no `mass_floor`
      heuristic (Section 4.2).
- [x] **`tau2` seeded from PR 03.** `log_dispersion_prev**2` is the warm start; MoM
      refines (Sections 5.2, 5.3). The two stages agree on iteration 1.
- [x] **One production exposure model.** Log-EB only. Gamma-Poisson is documented in
      the README as a future replacement, not a maintained parallel.
- [x] **Soft eligibility only.** `p_unexpressed_r` enters as a weight; no hard
      threshold (Section 5.1). Matches PR 01 / PR 03 contract.
- [x] **Numeric guards are wide.** Floor 1e-6, ceiling 1e6; capture spikes survive
      (Section 6.4; test 7 and 8).
- [x] **Deterministic iteration 1.** Bootstrap `BackgroundDensity` triggers
      `omega_r == 1` short-circuit; no leakage of bootstrap dispersion into downstream
      EM (Section 6.3, test 9).
- [x] **Tier 3 contributes nothing.** Excluded from `tau2` pool; `raw_ratio = 1` by
      construction (Sections 5.1, 6.2; test 13).
- [x] **PR 04 does not feed back into PR 03.** Strict topological data flow
      (Section 7).
- [x] **`A_r` retired.** Replaced everywhere by `RegionExposure.omega` (Section 9;
      test 14).
- [x] **No locus aggregation in PR 04.** Per-region output only; PR 05 owns the
      projection.

## 14. Done Means

- `RegionExposure` is produced for every region in every non-bootstrap calibration
  iteration; `omega_r = 1.0` everywhere in the bootstrap iteration with the flag set.
- `RegionCalibration.region_exposure` and `CalibrationStepResult.region_exposure` are
  populated; `A_r` is removed (or transitionally aliased to `omega.astype(np.float32)`).
- `tau2` is fit by weighted method-of-moments on the PR 03 Tier 1/2 pool, seeded from
  PR 03's `log_dispersion`, damped in log space.
- `v_obs_r` is the Gamma delta-method variance with `N_r`-based quasi-Poisson scaling.
- Per-region `omega_r` survives capture spikes (test 7) and shrinks low-support regions
  toward 1 (test 11).
- Ordinary RNA-seq simulation produces `omega_r` tightly concentrated near 1.0 (test 6).
- All tests in `tests/test_region_exposure.py` pass; existing calibration tests pass
  after mechanical fixture updates; pipeline smoke passes.
