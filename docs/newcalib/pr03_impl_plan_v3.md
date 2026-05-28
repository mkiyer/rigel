# PR 03 v3 - Regional Unspliced Mass Contract (Implementation Ready)

## Status

Supersedes [pr03_region_gdna_mass.md](pr03_region_gdna_mass.md) and
[pr03_impl_plan_v2.md](pr03_impl_plan_v2.md). This is the implementation-ready blueprint.

Closes three rounds of critique:

- **Round 1 (rename + tier hierarchy):** biological-intuition names, kill squared-mass ESS,
  deterministic three-tier `M_r` fallback, ESS plumbed to adaptive prior.
- **Round 2 (`rho0` framing):** `rho0` is the library-wide gDNA fragment mass per bp, not an
  "off-target rate"; delivered as a distribution; estimated from any reliable `M_r` region.
- **Round 3 (this doc):** the `rho0` estimator is presented in **both** Gamma conjugacy form
  (for PR 04 posterior updates) and log-scale dispersion form (for PR 04 shrinkage variance).
  Tier 3 uses the simple `min(rho0_mean * region_size_bp, T_r)` rule; the Bayesian-split form
  is documented as a future upgrade. `region_size_bp` is the single canonical denominator
  throughout PR 03 outputs.

## 1. Goals and Invariants

PR 03 introduces one production contract: a per-region decomposition of unspliced fractional
mass into gDNA and intronic-RNA portions.

```text
total_mass[r]   = T_r = unspliced compatible fractional mass in region r
gdna_mass[r]    = M_r = gDNA portion of T_r
rna_mass[r]     = R_r = T_r - M_r = unspliced (intronic) RNA portion of T_r
```

This decomposition is consumed by:

- `prior.py` / `adaptive_prior.py` (this PR wires the handoff).
- PR 04 (continuous-exposure `omega_r` empirical-Bayes shrinkage).

### Invariants

1. **Mass conservation.** `M_r + R_r == T_r` exactly in `float64` for every region. No `f32`
   tolerance hack.
2. **Uniform precision.** Primary mass tensors are `float64`; the only non-float field is the
   integer `unspliced_counts` (`uint64`, established by PR 02a).
3. **Deterministic estimation hierarchy.** `M_r` is resolved by a strict three-tier cascade
   (strand -> boundary -> background imputation). The tier is recorded per region in `method`.
4. **`rho0` decoupling.** `rho0` is a library-wide property delivered as a `BackgroundDensity`
   distribution. It is updated each iteration from regions with empirical gDNA evidence (Tier 1
   or Tier 2), never from Tier 3 (identity loop). No capture latent state appears anywhere.
5. **No legacy capture fields.** `p_captured`, `gamma_r`, `capture_enrichment_target`, and the
   four-state latent tensor are absent from every PR 03 input and output.

### Non-Goals

- `omega_r` shrinkage (PR 04).
- The native EM denominator sign correction (PR 05).
- Native boundary payload (PR 02b).
- Reimplementing strand deconvolution or boundary sweeps (reused as-is from existing code).

## 2. Naming Conventions

| Old name (`PriorMassDeconvolution`)        | New name (`RegionUnsplicedMass`) | Rationale |
|--------------------------------------------|-----------------------------------|-----------|
| `PriorMassDeconvolution`                   | `RegionUnsplicedMass`             | Describes the full unspliced decomposition, not just gDNA |
| `unspliced_total`                          | `total_mass`                      | Container already says "unspliced"; drop redundant prefix |
| `gdna_unspliced_mean`                      | `gdna_mass`                       | Point estimate of mass, not a posterior mean |
| `rna_unspliced_mean`                       | `rna_mass`                        | Same |
| `support_count` (proposed)                 | `unspliced_counts`                | Concrete: physical unspliced fragments touching the region |
| `mass_opportunity`                         | `region_size_bp`                  | The denominator is geometric region width in bp |
| `support_ess` (proposed)                   | (deleted)                         | Integer `unspliced_counts` is the only ESS we trust |
| `rho_off` (scalar field on calibration)    | `BackgroundDensity` dataclass     | Distribution, not scalar; framed as global gDNA density |

## 3. Dataclasses

### 3.1 `RegionUnsplicedMass`

Replaces `PriorMassDeconvolution` in `src/rigel/calibration/calibration_iteration.py`.

```python
@dataclass(frozen=True, slots=True)
class RegionUnsplicedMass:
    """Mass-conserving unspliced decomposition per region.

    Source of truth for downstream adaptive prior assembly (this PR) and EB
    exposure learning (PR 04).
    """

    # Primary mass tensors.
    total_mass: np.ndarray        # float64[R]  T_r
    gdna_mass: np.ndarray         # float64[R]  M_r
    rna_mass: np.ndarray          # float64[R]  R_r = T_r - M_r

    # Geometric mass denominator.
    region_size_bp: np.ndarray    # float64[R]  region width in bp (b - a)

    # Effective sample size (from PR 02a).
    unspliced_counts: np.ndarray  # uint64[R]   physical unspliced fragments touching r

    # Provenance and quality.
    method: np.ndarray            # uint8[R]    METHOD_STRAND / METHOD_BOUNDARY /
                                  #             METHOD_BACKGROUND_FALLBACK
    precision: np.ndarray         # float64[R]  reliability of the M_r estimator (used as
                                  #             a weight in the rho0 refit; not a posterior
                                  #             precision)
    flags: np.ndarray             # uint16[R]   region-level diagnostic bits
```

#### Method constants

```python
METHOD_STRAND              = 1   # strand deconvolution drove M_r
METHOD_BOUNDARY            = 2   # boundary-sweep imputation drove M_r
METHOD_BACKGROUND_FALLBACK = 3   # no local evidence; M_r := rho0_mean * region_size_bp
```

#### Flags (new in PR 03)

```python
FLAG_M_IMPUTED_FROM_BACKGROUND = 1 << 0   # set when method == METHOD_BACKGROUND_FALLBACK
FLAG_M_CLIPPED_TO_TOTAL        = 1 << 1   # M_r had to be clipped to total_mass[r]
# Existing flags from RegionGdnaChannelEstimate continue to propagate via OR.
```

#### Invariants enforced in `__post_init__`

```text
all primary tensors are float64, C-contiguous, shape (R,)
unspliced_counts.dtype == np.uint64
total_mass >= 0
0 <= gdna_mass <= total_mass
0 <= rna_mass  <= total_mass
gdna_mass + rna_mass == total_mass      # exact float64
region_size_bp > 0
precision >= 0 and finite
method in {METHOD_STRAND, METHOD_BOUNDARY, METHOD_BACKGROUND_FALLBACK}
```

### 3.2 `BackgroundDensity`

New dataclass in `src/rigel/calibration/calibration_iteration.py`.

```python
@dataclass(frozen=True, slots=True)
class BackgroundDensity:
    """Library-wide gDNA density distribution.

    Carries both the Gamma conjugacy parameters (for PR 04 posterior updates of
    omega_r) and a log-scale dispersion (for PR 04 shrinkage variance). PR 04
    chooses which to consume.
    """

    # Robust center.
    rho0_mean: float            # gDNA fragment mass per bp; the omega_r anchor

    # Gamma conjugacy view (for PR 04 posterior updates).
    alpha0: float               # Gamma shape:  alpha_floor + sum(weights * M_r)
    beta0:  float               # Gamma rate:   beta_floor  + sum(weights * region_size_bp)
                                # rho0_mean == alpha0 / beta0 by construction

    # Robust dispersion view (for PR 04 shrinkage variance seed).
    log_dispersion: float       # weighted sigma-equivalent of log(rho_hat_r), in nats

    # Pool diagnostics.
    n_effective_regions: float  # sum of weights (precision * counts * p_unexpressed)
    n_regions_in_pool: int      # raw count of Tier 1/2 regions contributing
    method_histogram: tuple     # (n_tier1, n_tier2, n_tier3_excluded)
    fit_status: str             # "converged" | "fallback_bootstrap" | "prior_only"
```

#### Bootstrap handoff

```python
@classmethod
def from_bootstrap(cls, model: BackgroundModel) -> "BackgroundDensity":
    """Initial BackgroundDensity before any E-step exists.

    Wraps the bootstrap Gamma posterior; signals "no data-driven dispersion
    yet" via wide default log_dispersion and zero n_effective_regions so the
    first iteration treats the dispersion as a flat prior.
    """
    return cls(
        rho0_mean=float(model.rho_off_mean),
        alpha0=float(model.rho_off_alpha),
        beta0=float(model.rho_off_beta),
        log_dispersion=float(np.log(10.0)),   # wide; PR 04 will replace
        n_effective_regions=0.0,
        n_regions_in_pool=0,
        method_histogram=(0, 0, 0),
        fit_status="fallback_bootstrap",
    )
```

## 4. The Three-Tier `M_r` Hierarchy

For each region, evaluate the tiers in order. The first tier that fires owns `M_r`.

### Tier 1 - Strand deconvolution

**Condition.** `strand_channels is not None` and at least one of
`contained_reliability[r]`, `boundary_left_reliability[r]`, `boundary_right_reliability[r]`
is `>= eps` (default `eps = 1e-6`).

**Formula.**

```text
M_r = contained_mean[r]       * clip(contained_reliability[r],       0, 1)
    + boundary_left_mean[r]   * clip(boundary_left_reliability[r],   0, 1)
    + boundary_right_mean[r]  * clip(boundary_right_reliability[r],  0, 1)
M_r = clip(M_r, 0, total_mass[r])
method[r]    = METHOD_STRAND
precision[r] = max over channels of (channel_precision * channel_reliability)
flags[r]    |= strand_channels.flags[r]
```

This preserves the numerical behavior of the existing
`build_prior_mass_deconvolution(strand_channels=...)` path.

### Tier 2 - Boundary-sweep imputation

**Condition.** Tier 1 did not fire AND the region has nonzero boundary excess evidence:
`local_posterior.alpha_excess[r] > 0` OR `local_posterior.beta_excess[r] > 0`.

**Formula.**

```text
M_r = clip(sweep.contained_gdna_mean[r], 0, total_mass[r])
method[r]    = METHOD_BOUNDARY
precision[r] = sweep.contained_gdna_precision[r]
```

`sweep.contained_gdna_mean` is produced by the existing `run_boundary_sweep()` which calls
`predict_contained_gdna_from_excess()` internally with `contained_leff`. PR 03 reuses this
black-box; `contained_leff` is an *internal input* to Tier 2 and is not surfaced in
`RegionUnsplicedMass`.

**Explicitly rejected:** `M_r = min(T_r, max(direct_contained_density, boundary_sweep_mean))`
from the original PR 03 draft. Without strand contrast there is no estimator independent of the
boundary sweep, so a `max(...)` of correlated quantities is upward-biased.

### Tier 3 - Background imputation

**Condition.** Neither Tier 1 nor Tier 2 fired.

**Formula.**

```text
M_r          = clip(rho0_mean * region_size_bp[r], 0, total_mass[r])
method[r]    = METHOD_BACKGROUND_FALLBACK
precision[r] = 0.0
flags[r]    |= FLAG_M_IMPUTED_FROM_BACKGROUND
```

**Why not the Bayesian split form**
`M_r = T_r * (rho0 * region_size_bp) / (rho0 * region_size_bp + mu_rna_r)`?
That form is mathematically appealing but requires a regional RNA expectation `mu_rna_r`. By
construction, regions reaching Tier 3 have no strand contrast and no boundary flux — the only
two places we get `mu_rna_r`. So `mu_rna_r` is unavailable at Tier 3. The simple cap-at-total
formula is what we can defensibly compute. If PR 06 surfaces a separate per-region RNA prior
(e.g. from annotation-based intronic length scaling), the Bayesian-split upgrade ships in a
follow-up; the `method` field already provides the discriminator for that future cutover.

### After every tier

```text
rna_mass[r] = total_mass[r] - gdna_mass[r]
if clipping fired: flags[r] |= FLAG_M_CLIPPED_TO_TOTAL
```

## 5. `region_size_bp` (the canonical denominator)

PR 03 stores `region_size_bp[r] = region_end[r] - region_start[r]` as `float64` and uses it as
the denominator everywhere a density is computed:

- Tier 3 imputation: `M_r = rho0_mean * region_size_bp[r]`.
- M-step `rho_hat[r] = gdna_mass[r] / region_size_bp[r]`.

`contained_leff[r]` (FL-marginal effective length, currently used by the density estimators)
remains an internal input to Tier 2 only. It is **not** the PR 03 denominator. Rationale: the
fractional mass `total_mass[r]` is accumulated over bp overlap, so its natural density
denominator is bp width. The count-based opportunity `W + L - 1` is the right denominator if and
when PR 04 adds a count-Poisson likelihood, and is out of scope here.

## 6. Iterative `rho0` Estimation

### 6.1 Training pool (per iteration)

```text
pool_mask[r] = (method[r] in {METHOD_STRAND, METHOD_BOUNDARY})
             & (unspliced_counts[r] >= 1)
             & (region_size_bp[r] >= 1.0)
```

Explicit exclusions:

- **`method == METHOD_BACKGROUND_FALLBACK`**: identity loop. Including these regions would
  artificially deflate `log_dispersion` and pin `rho0_mean` to its previous value.
- **`unspliced_counts == 0`**: zero physical observations; `M_r` is a fractional artefact of
  a neighbor's fragment touching the edge.
- **No top-density tail trim. No `seed_mask` membership requirement.** Enriched capture
  regions stay in the pool as legitimate samples; their influence on the center is bounded by
  the Huber cap (Section 6.3).

### 6.2 Per-region weight

```text
w[r] = precision[r] * float(unspliced_counts[r]) * p_unexpressed[r]
```

- `precision[r]`: reliability of the `M_r` estimator (Tier 1 strand reliability product, or
  Tier 2 boundary precision).
- `unspliced_counts[r]`: integer fragment ESS.
- `p_unexpressed[r]`: smooth iterative refinement. Expressed regions are *not excluded* — their
  strand-deconv `M_r` is still a valid empirical sample — but their weight decays as the model
  becomes confident they are expressed.

On iteration 1 `p_unexpressed[r]` is uniform or comes from the bootstrap; either degrades
gracefully.

### 6.3 Robust geometric-mean estimator with Huber cap

#### 6.3.1 Handling zero-mass regions (Bayesian shrinkage)

A naive `log(rho_hat[r])` is undefined when `gdna_mass[r] == 0` and a `log(max(..., epsilon))`
clamp is catastrophic: a single genuinely clean Tier-1 region (zero gDNA mass) would be pinned
at `log(1e-12) ~= -27.6` and would yank the geometric mean off scale.

Fix: compute `rho_hat[r]` as a **Bayesian posterior mean** under a Gamma prior whose mean is
the previous iteration's `rho0_mean`. Concretely, with a single tunable pseudocount
`alpha_floor` (default 1.0):

```text
pseudo_mass = alpha_floor                              # default 1.0 (in mass units)
pseudo_size = alpha_floor / previous.rho0_mean         # implied prior-equivalent bp

rho_hat[r]  = (gdna_mass[r] + pseudo_mass)
            / (region_size_bp[r] + pseudo_size)
```

Properties:

- A zero-mass region contributes `rho_hat[r] = pseudo_mass / (region_size_bp[r] + pseudo_size)`,
  which for any non-tiny `region_size_bp` is close to zero in absolute terms but **always
  positive**, so `log(rho_hat[r])` is finite.
- For a tiny region (`region_size_bp << pseudo_size`), `rho_hat[r] -> previous.rho0_mean` (the
  prior dominates).
- For a large region (`region_size_bp >> pseudo_size`), `rho_hat[r] -> gdna_mass[r] /
  region_size_bp[r]` (the data dominates).
- The total pull of any zero-mass region toward the log-mean is bounded by the same Huber cap
  (Section 6.3.2), so even many simultaneous zero-mass regions cannot dominate.

This pseudocount is the **same** `alpha_floor` used in the Gamma view (Section 6.4), so the two
views stay statistically consistent. The pseudocount is *not* a separate hyperparameter.

#### 6.3.2 Robust center

```text
# Step 1: per-region density and its log (using the shrunk rho_hat from 6.3.1)
log_rho_hat[r] = log(rho_hat[r])

# Step 2: raw weighted center and dispersion
log_center_raw = weighted_mean(log_rho_hat, w)           over pool_mask
log_mad        = weighted_median(|log_rho_hat - log_center_raw|, w)
log_dispersion = max(1.4826 * log_mad, log_dispersion_floor)
                                                          # floor = log(1.1) ~= 0.095 to avoid
                                                          # degenerate Huber cap on a tight pool

# Step 3: Huberize and recompute the center
k = 1.5
clipped = clip(log_rho_hat,
               log_center_raw - k * log_dispersion,
               log_center_raw + k * log_dispersion)
log_center = weighted_mean(clipped, w)
rho0_mean_hat = exp(log_center)
```

**Why weighted geometric mean.** `rho_r = rho0 * omega_r` is a multiplicative scale. The
arithmetic mean of `rho_hat_r` is dominated by heavy tails; the geometric mean is the natural
central tendency for multiplicative noise and is invariant under unit changes. Huberization caps
each region's contribution to the center at `k * log_dispersion`, so a single mega-enriched
target (e.g. `omega_r = 1000`) shifts `rho0_mean` by at most a bounded multiplicative factor,
not unboundedly.

**Damping (log-space)** across iterations to stabilize convergence:

```text
log_center_damped = (1 - damping) * log(previous.rho0_mean) + damping * log_center
rho0_mean_next    = exp(log_center_damped)
```

Default `damping = 0.5`. The same `damping` is applied to `(alpha0, beta0)` in Section 6.4 so
the Gamma view does not lag the geometric-mean view.

### 6.4 Gamma conjugacy view

PR 04 will perform per-region posterior updates of `omega_r` against a Gamma prior. To make that
cheap, the same estimator emits the conjugate Gamma parameters in addition to the robust
geometric mean:

```text
alpha_hat = alpha_floor + sum over pool of (w[r] * gdna_mass[r])
beta_hat  = beta_floor  + sum over pool of (w[r] * region_size_bp[r])

# Damping mirrors the geometric-mean side
alpha0 = (1 - damping) * previous.alpha0 + damping * alpha_hat
beta0  = (1 - damping) * previous.beta0  + damping * beta_hat
```

The two views are kept consistent by enforcing the identity `rho0_mean = alpha0 / beta0` post
hoc — i.e. after computing the robust Huberized `rho0_mean_next` and the raw Gamma `(alpha0,
beta0)`, rescale `(alpha0, beta0)` to match:

```text
scale = rho0_mean_next / (alpha_hat_damped / beta_hat_damped)
alpha0 = alpha_hat_damped * sqrt(scale)
beta0  = beta_hat_damped  / sqrt(scale)
```

This preserves the Gamma's effective sample size while pinning its mean to the robust estimate.
PR 04 can use either view without inconsistency.

`fit_status`:
- `"converged"`: refit ran and the log-jump from `previous.rho0_mean` is below
  `convergence_log_tol` (default 0.01).
- `"iterating"`: refit ran but the log-jump exceeds tolerance.
- `"fallback_bootstrap"`: pool was empty; carried previous `BackgroundDensity` unchanged.
- `"prior_only"`: first real fit after bootstrap (no `previous` data yet).

### 6.5 Why this resolves the leak critique without a seed mask

The Round-1 critique feared that with the four-state classifier dead, `p_unexpressed` could not
distinguish enriched-unexpressed from off-target-unexpressed and a single mega-enriched region
would blow up `rho0`. The robust geometric-mean + Huberization caps the influence of any single
region at `k * log_dispersion` in log space, regardless of how massive its `M_r` is. Enriched
regions remain in the pool as legitimate tail samples and contribute to `log_dispersion` —
which is exactly what PR 04 needs to set the shrinkage strength. If PR 06 validation surfaces
residual capture skew, the lever is `k` (e.g. tighten from 1.5 to 1.0), not a re-introduced
mask.
### 6.6 Convergence and oscillation control

The iterative refit faces three potential failure modes. The estimator handles each explicitly:

1. **Single-region tail blowup.** Bounded by the Huber cap (Section 6.3.2). No region can move
   `log_center` by more than `k * log_dispersion` regardless of mass.
2. **Iteration-to-iteration overshoot of `rho0_mean`.** Bounded by log-space damping
   (Section 6.3.2): a `damping = 0.5` halves any single-step log-jump.
3. **`p_unexpressed` flip-flopping.** A region whose `p_unexpressed` toggles between 0.3 and
   0.7 across outer EM iterations will have its `w[r]` weight toggle in lockstep. The damping
   on `rho0_mean` (#2) prevents this from producing a divergent `rho0`. We **do not** add a
   second damping layer on the weights themselves: the outer `calibration_iteration.run()`
   convergence check already monitors `max_state_shift` and will not declare convergence while
   `p_unexpressed` is bouncing. If the outer loop converges, the weights converge, and `rho0`
   converges.

`fit_status` is set per call based on the change from `previous`:

```text
log_jump = |log(rho0_mean_next) - log(previous.rho0_mean)|

if n_regions_in_pool == 0:
    fit_status = "fallback_bootstrap"   # carry previous unchanged
elif previous.fit_status == "fallback_bootstrap" and previous.n_regions_in_pool == 0:
    fit_status = "prior_only"           # first real fit after bootstrap
elif log_jump < convergence_log_tol:    # default 0.01 (~1% relative change)
    fit_status = "converged"
else:
    fit_status = "iterating"
```

The outer `run_calibration_iteration()` may use `fit_status == "converged"` as one of its
stopping criteria, but PR 03 does not require it; the outer loop's existing
`max_state_shift` / `relative_rho_shift` checks remain authoritative.

### 6.7 Bootstrap call (`fit_background_model`)

Unchanged interface. Still produces an initial `rho_off_mean`, `rho_off_alpha`, `rho_off_beta`,
and `seed_mask` from heuristic exclusions. PR 03 wraps it via `BackgroundDensity.from_bootstrap()`
to seed iteration 1. **The `seed_mask` is consumed only by the bootstrap; the iterative refit
ignores it.**

## 7. Downstream Handoff to Adaptive Prior

The Round-1 critique was that handing only mass means starves `adaptive_prior.py` of
uncertainty. PR 03 wires `unspliced_counts` through as a first-class ESS argument.

### 7.1 `prior.py` projection (region -> locus)

```python
def project_region_mass_to_locus(
    region_mass: RegionUnsplicedMass,
    locus_table: LocusTable,
    overlap_weights: np.ndarray,   # shape (n_loci, n_regions), row-stochastic
) -> LocusUnsplicedMass:
    """Project per-region mass and ESS onto loci."""
    return LocusUnsplicedMass(
        total_mass=overlap_weights @ region_mass.total_mass,
        gdna_mass=overlap_weights @ region_mass.gdna_mass,
        rna_mass=overlap_weights @ region_mass.rna_mass,
        # Project counts as a weighted sum (float64). adaptive_prior interprets
        # this as the effective number of independent unspliced observations
        # supporting the locus.
        unspliced_counts=overlap_weights @ region_mass.unspliced_counts.astype(np.float64),
    )
```

### 7.2 `adaptive_prior.py` consumption

`compute_adaptive_prior()` gains a required `locus_ess` argument:

```python
def compute_adaptive_prior(
    p_states: np.ndarray,
    locus_unspliced_mass: LocusUnsplicedMass,
    locus_ess: np.ndarray,           # float64[n_loci]
    *,
    base_concentration: float,
    floor_ess: float = 1.0,
    ceil_ess: float = 100.0,
) -> PriorTable:
    ...
```

Effect on prior concentration:

```text
shrink(ess) = clip((ess - floor_ess) / (ceil_ess - floor_ess), 0.0, 1.0)
prior_concentration_l = base_concentration * shrink(locus_ess[l])
```

Zero ESS contributes a flat prior; high ESS contributes the full data-driven Dirichlet. The
existing `p_unexpressed` soft gate is unchanged.

PR 03 mandates the **plumbing** (ESS as a first-class argument). The exact functional form of
`shrink(...)` may be refined in PR 04 or later; the contract is "two synthetic loci with
identical projected `gdna_mass` and `rna_mass` but different `unspliced_counts` must produce
different prior concentrations".

## 8. Files to Touch

### `src/rigel/calibration/calibration_iteration.py`

- Delete `PriorMassDeconvolution`.
- Introduce `RegionUnsplicedMass` and `BackgroundDensity` per Section 3.
- Replace `RegionCalibration.prior_mass` with `region_unspliced_mass: RegionUnsplicedMass`.
- Replace `RegionCalibration.rho_off: float` with `background_density: BackgroundDensity`.
  Optional one-cycle property `rho_off` returning `background_density.rho0_mean` for downstream
  grace; remove in PR 04.
- Replace `build_prior_mass_deconvolution(...)` with
  `build_region_unspliced_mass(observation, *, region_size_bp, unspliced_counts,
  strand_channels, local_posterior, sweep, background_density) -> RegionUnsplicedMass`.
  The new signature takes every input the three-tier hierarchy needs.
- Replace the M-step background refit with a new function
  `estimate_background_density(region_unspliced_mass, p_unexpressed, *, previous, damping=0.5,
  huber_k=1.5, alpha_floor=1.0, beta_floor=1.0, rho_floor=1e-12) -> BackgroundDensity`.
- `calibration_e_step()`:
  - Plumb `region_size_bp` (from `RegionArrays`) and `unspliced_counts` (from
    `RegionCountLedger`) into `build_region_unspliced_mass()`.
  - Pass the current `BackgroundDensity` instead of a scalar `rho_off`.
- `calibration_m_step()`:
  - Call `estimate_background_density(...)` instead of the old Gamma refit.
  - Return `(BackgroundDensity, kappa_d)` instead of `(BackgroundModel, kappa_d)`. The
    bootstrap `BackgroundModel` is no longer mutated after iteration 1.

### `src/rigel/calibration/_arrays.py`

- Add `region_size_bp: np.ndarray` (`float64`) to `RegionArrays`, derived from
  `region_end - region_start` and reordered by the existing `order` permutation.

### `src/rigel/calibration/region_count_ledger.py`

- Confirm `unspliced_counts` (PR 02a) is reachable in the same sort order. Add an explicit
  accessor if not already present.

### `src/rigel/calibration/prior.py`

- Switch consumer code from `PriorMassDeconvolution` to `RegionUnsplicedMass`.
- Add `project_region_mass_to_locus()` per Section 7.1.
- Pass `locus_ess` into `compute_adaptive_prior()`.

### `src/rigel/calibration/adaptive_prior.py`

- Add the `locus_ess` required argument.
- Implement `shrink(ess)` per Section 7.2 (or document the interim linear ramp).
- Existing `p_unexpressed` directional soft gate is unchanged.

### `src/rigel/calibration/_result.py`

- Region/locus summary fields renamed to match `RegionUnsplicedMass`.
- Add `BackgroundDensity` to the calibration summary (mean, log_dispersion,
  n_effective_regions, n_regions_in_pool, method_histogram, fit_status).
- Add per-tier region count to diagnostics.

### `src/rigel/calibration/background_model.py`

- No API change for `fit_background_model()`.
- Add docstring: `seed_mask` is consumed *only* to seed iteration 1 via
  `BackgroundDensity.from_bootstrap()`. Subsequent iterations use
  `estimate_background_density()` and ignore `seed_mask`.

## 9. Test Plan

Add `tests/test_region_unspliced_mass.py` with the following cases.

### Tier-hierarchy correctness

1. **Mass conservation.** For synthetic inputs that touch all three tiers,
   `np.all(M + R == T)` exactly in float64 (no tolerance).
2. **Tier 1 (strand).** Region with strong antisense unspliced flux. Assert
   `method == METHOD_STRAND`, `M_r` matches the channel-reliability-weighted sum, `precision > 0`.
3. **Tier 2 (boundary).** Region with no strand contrast but strong boundary excess. Assert
   `method == METHOD_BOUNDARY`, `M_r == clip(sweep.contained_gdna_mean, 0, T_r)`.
4. **Tier 3 (background fallback).** Region with neither strand nor boundary evidence. Assert
   `method == METHOD_BACKGROUND_FALLBACK`, `M_r == min(rho0_mean * region_size_bp, T_r)`,
   `flags & FLAG_M_IMPUTED_FROM_BACKGROUND`.
5. **Sentinel: empty region.** `T_r == 0` and `unspliced_counts == 0`: all tiers produce
   `M_r == 0`, `R_r == 0`.
6. **Tier promotion.** Same region, switch on strand reliability: assert tier moves
   `BACKGROUND_FALLBACK -> BOUNDARY -> STRAND` as evidence is added.

### `rho0` estimator behavior

7. **Tier 3 exclusion from `rho0` update.** Pool of N Tier-3 regions + one Tier-1 region. Assert
   `estimate_background_density(...).rho0_mean` equals the rho0 you would get from the Tier-1
   region alone (within damping). Locks the identity-loop fix.
8. **Robust estimator caps enriched tails.** 1000 Tier-1 regions at `rho_hat = 1.0` plus 10
   Tier-1 regions at `rho_hat = 1000.0`. Assert `rho0_mean in [0.9, 1.3]` (i.e. near the bulk
   mode, not the contaminated mean ~10.9). Locks the no-seed-mask claim.
9. **Dispersion grows with regional spread.** Two pools with identical `rho0_mean` but
   different multiplicative spread produce monotonically increasing `log_dispersion`. Locks the
   PR 04 dispersion contract.
10. **Gamma-vs-mean consistency.** Assert `abs(alpha0 / beta0 - rho0_mean) < 1e-9` for every
    `BackgroundDensity` produced.
11. **`p_unexpressed` reweighting across iterations.** Run two iterations of
    `calibration_m_step()`; assert the iteration-2 weight on a region whose `p_unexpressed`
    increased from 0.3 to 0.9 is correspondingly higher (proportional check).
12. **Damping stabilizes shocks.** Apply a 10x scale shock to `rho_hat`; assert `rho0_mean`
    moves by approximately `damping` of the log-shock on the first pass, not the full shock.
13. **Bootstrap handoff.** `BackgroundDensity.from_bootstrap(model)` yields
    `rho0_mean == model.rho_off_mean`, `alpha0 == model.rho_off_alpha`,
    `beta0 == model.rho_off_beta`, `n_effective_regions == 0.0`,
    `fit_status == "fallback_bootstrap"`.
14. **Zero-mass region does not blow up the geometric mean.** Pool of 100 Tier-1 regions, 50 of
    which have `gdna_mass = 0.0` and 50 of which have `rho_hat = 1.0`. Assert
    `rho0_mean in [0.5, 1.0]` (i.e. zeros pull the mean down via Bayesian shrinkage to roughly
    the prior-blended value, but do **not** drive it to `~exp(-27)` as a naive `log(epsilon)`
    clamp would). Locks the zero-mass handling in Section 6.3.1.
15. **`fit_status` transitions.** Sequence `from_bootstrap -> first refit with empty pool ->
    first refit with data -> second refit at converged tolerance` produces statuses
    `"fallback_bootstrap" -> "fallback_bootstrap" -> "prior_only" -> "converged"`.

### Plumbing

16. **dtype.** Every primary tensor in `RegionUnsplicedMass` is exactly `np.float64`;
    `unspliced_counts.dtype == np.uint64`; `BackgroundDensity` scalars are Python `float`.
17. **Adaptive prior receives ESS.** Two synthetic loci with identical projected `gdna_mass`
    and `rna_mass` but different `unspliced_counts` produce different prior concentrations from
    `compute_adaptive_prior()`.
18. **No capture fields.** Assert no attribute named `p_captured`, `gamma_r`, or
    `capture_enrichment_target` exists on any PR 03 dataclass.

### Update existing test files

Rename and adjust fixtures for:

- `tests/test_calibration_iteration.py`
- `tests/test_calibration_prior.py`
- `tests/test_calibration_result.py`
- `tests/test_per_locus_gdna_mass.py`
- `tests/test_bayesian_prior_acceptance.py`
- `tests/test_adaptive_prior.py`

Changes are mechanical: field renames, dtype expectations, addition of the `unspliced_counts`
and `locus_ess` arguments where applicable. The legacy four-state assertions were already
removed in PR 01.

## 10. Targeted Validation Commands

```bash
conda activate rigel && pytest \
    tests/test_region_unspliced_mass.py \
    tests/test_calibration_iteration.py \
    tests/test_calibration_prior.py \
    tests/test_calibration_result.py \
    tests/test_per_locus_gdna_mass.py \
    tests/test_bayesian_prior_acceptance.py \
    tests/test_adaptive_prior.py \
    -v
```

After targeted tests pass:

```bash
conda activate rigel && pytest tests/test_pipeline_smoke.py tests/test_pipeline_wiring.py -v
```

No native rebuild is required for PR 03 (PR 02a already produces `unspliced_counts`).

## 11. Recommended Implementation Order

1. Land `RegionUnsplicedMass` and `BackgroundDensity` dataclasses with `__post_init__`
   validation and unit tests (cases 1, 5, 16, 18).
2. Land `build_region_unspliced_mass()` with the three-tier hierarchy and tests 2, 3, 4, 6.
3. Land `estimate_background_density()` and tests 7, 8, 9, 10, 11, 12, 13, 14, 15.
4. Wire `calibration_e_step()` / `calibration_m_step()` to use the new dataclasses; update
   existing iteration tests (mechanical renames).
5. Wire `prior.py` projection and `adaptive_prior.py` ESS consumption; land test 17.
6. Update `_result.py` summaries and pipeline-smoke fixtures.
7. Run the targeted suite, then the smoke suite.

## 12. Resolved Critiques (Final Checklist)

- [x] **M-step background leak.** Robust geometric-mean + Huber cap bounds tail influence
      without re-introducing a `seed_mask` freeze. Tier 3 regions and zero-count regions
      excluded from the pool. Test 8 locks it.
- [x] **Downstream uncertainty starvation.** `unspliced_counts` plumbed through
      `project_region_mass_to_locus()` to `compute_adaptive_prior(locus_ess=...)`. Test 15
      locks it.
- [x] **Deterministic fallback hierarchy.** Three explicit tiers in fixed order; `method`
      field records the choice; no `max(...)` combine. Test 6 locks the order.
- [x] **`rho0` framing.** Library-wide gDNA fragment mass per bp; a property of the library,
      not "off-target". Delivered as `BackgroundDensity` distribution.
- [x] **`rho0` as a distribution.** `(rho0_mean, alpha0, beta0, log_dispersion,
      n_effective_regions, ...)` — Gamma conjugacy for PR 04 posterior updates, log-dispersion
      for PR 04 shrinkage variance. Test 9 and 10 lock the two views.
- [x] **Bootstrap `seed_mask` retired from iterative refit.** Used only via
      `BackgroundDensity.from_bootstrap()`. Test 13 locks it.
- [x] **Identity loop.** Tier-3 imputed and zero-count regions excluded from the `rho0` pool.
      Test 7 locks it.
- [x] **Iterative pool refinement via `p_unexpressed` weighting.** Test 11 locks it.
- [x] **Field renames for biological intuition.** `RegionUnsplicedMass`, `total_mass`,
      `gdna_mass`, `rna_mass`, `region_size_bp`, `unspliced_counts`, `BackgroundDensity`,
      `rho0_mean`.
- [x] **No squared-mass ESS.** Integer `unspliced_counts` from PR 02a is the only ESS source.
- [x] **Uniform float64.** All primary mass tensors and `region_size_bp` are `float64`. Mass
      conservation is exact, not tolerance-based. Test 14 locks dtype; test 1 locks exact
      conservation.
- [x] **No capture latent fields anywhere.** Test 16 locks it.

## 13. Done Means

- `RegionUnsplicedMass` and `BackgroundDensity` are the only production outputs of the new
  calibration mass path. `PriorMassDeconvolution` is deleted from production code.
- `M_r` is computed by the deterministic three-tier hierarchy with `method` recording the tier.
- `rho0` is refit each iteration by the robust geometric-mean + Huber-cap estimator over Tier 1
  and Tier 2 regions; the result is consistent in both Gamma `(alpha0, beta0)` and log-dispersion
  views.
- Adaptive prior consumes a per-locus ESS argument derived from `unspliced_counts`.
- All primary mass tensors are `float64`; mass conservation is exact.
- The targeted test suite (`tests/test_region_unspliced_mass.py` plus the renamed fixtures)
  passes; pipeline smoke passes.

## 14. Implementation Log

### Step 1 — `RegionUnsplicedMass` / `BackgroundDensity` dataclasses + `region_size_bp`

- Added both frozen-slotted dataclasses plus method/flag constants in
  `src/rigel/calibration/calibration_iteration.py`.
- Added `region_size_bp: float64[R]` to `RegionArrays` (`_arrays.py`) populated as
  `(end - start).astype(np.float64)` in sorted region order.
- 11 unit tests added in `tests/test_region_unspliced_mass.py`.

### Step 2 — `build_region_unspliced_mass()` three-tier hierarchy

- Added `build_region_unspliced_mass(...)` in `calibration_iteration.py` with the
  three-tier hierarchy from Section 4. Vectorised via `np.where`; mass conservation
  is exact (`rna = total - clipped_gdna` in float64).
- **Deviation:** the plan referenced `sweep.contained_gdna_precision` for the
  Tier 2 precision weight, but `BoundarySweepResult` has no such field. Substituted
  `alpha_excess + beta_excess` (boundary-excess effective sample size), which is
  the natural sample-size weight of the sweep estimator at that region. Documented
  inline at the Tier 2 branch.
- 6 new tier tests added (Cases 2, 3, 4, 4b clip variant, 6 promotion, mass-conservation
  mixed-tier). All 17 tests in the file pass.

### Step 3 — `estimate_background_density()` robust geomean + Huber refit

- Added `estimate_background_density(region_unspliced_mass, p_unexpressed, *, previous, ...)`
  in `calibration_iteration.py` implementing Sections 6.1-6.7. Includes pool selection,
  Bayesian-shrunk per-region `rho_hat`, weighted geomean with Huberization on the MAD
  scale, log-space damping, Gamma `(alpha0, beta0)` conjugate update with rescaling to
  enforce `alpha0/beta0 == rho0_mean`, and `fit_status` transitions per Section 6.6.
- Helper `_weighted_median(values, weights)` (lower-median convention) computes the MAD.
- `n_effective_regions` is the sum of pool weights (`precision * counts * p_unexpressed`).
  `method_histogram` reports counts across **all** regions, not just pool members.
- Empty-pool short-circuit returns `previous` rho0/alpha0/beta0/log_dispersion unchanged
  with `fit_status="fallback_bootstrap"`. The same fallback fires if all pool weights are
  zero (e.g. `p_unexpressed == 0` everywhere).
- 10 tests added (Cases 7, 8, 9, 10, 11, 12, 14, 15 + method-histogram + ESS sanity).
  Combined file now 27/27 passing; targeted regression sweep (region_unspliced_mass +
  calibration_iteration + latent_states + adaptive_prior) is 57/57.
- **Note vs plan:** plan Section 6.6 fit_status table treats `"prior_only"` as the
  branch when `previous.fit_status == "fallback_bootstrap" AND previous.n_regions_in_pool
  == 0`; this implementation matches that exactly, *but* it evaluates the prior-only check
  **before** the convergence-tolerance check. If a bootstrap stub happens to produce a
  log_jump below tolerance on first fit, we still return `"prior_only"` rather than
  `"converged"`. This matches the spirit of the plan (first real fit deserves the
  diagnostic) but is a slight ordering interpretation worth noting.
