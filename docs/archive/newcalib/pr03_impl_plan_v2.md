# PR 03 v2 - Regional Unspliced Mass Implementation Plan

## Status

Supersedes [pr03_region_gdna_mass.md](pr03_region_gdna_mass.md). Incorporates two rounds of
critique on 2026-05-28:

- **Round 1**: rename for biological intuition, kill the squared-mass ESS, define a deterministic
  fallback hierarchy for `M_r`, pass uncertainty (not just means) downstream, rewrite the
  unstranded fallback section.
- **Round 2**: redefine `rho0` as the **library-wide global gDNA density** (a property of the
  library, not of "unexpressed regions"). Estimate it from any region with a reliable `M_r`,
  not from the bootstrap `seed_mask`. Deliver it as a distribution `(mean, dispersion, ESS)`
  rather than a scalar, so PR 04 has the variance it needs for shrinkage. Use a **robust
  location estimator** to keep enriched capture tails from dominating the mean.

## Purpose

PR 03 introduces the single production contract for the regional unspliced mass decomposition:

```text
total_mass[r]  = T_r = total unspliced compatible fractional mass in region r
gdna_mass[r]   = M_r = portion of T_r attributable to gDNA
rna_mass[r]    = R_r = T_r - M_r = portion of T_r attributable to unspliced (intronic) RNA
```

This is the source of truth that

- PR 04 (`omega_r` empirical-Bayes exposure model) consumes to compute regional density ratios.
- `prior.py` / `adaptive_prior.py` consume to construct per-locus Dirichlet priors over
  `(mRNA, nRNA, gDNA)` states.

After PR 03 there must be exactly one production regional unspliced mass table, computed without
reference to any capture latent state.

## Non-Goals

- Do not implement `omega_r` shrinkage. That is PR 04.
- Do not touch the EM denominator sign. That is PR 05.
- Do not move the boundary payload to native. That is PR 02b.
- Do not rebuild strand deconvolution or boundary sweeps. They are reused as-is.

## Naming Conventions (Biological Intuition First)

The old `PriorMassDeconvolution` dataclass describes more than gDNA mass and uses opaque "support"
and "opportunity" language. PR 03 renames for clarity:

| Old name                       | New name              | Reason                                              |
|--------------------------------|-----------------------|-----------------------------------------------------|
| `PriorMassDeconvolution`       | `RegionUnsplicedMass` | Describes the full unspliced decomposition          |
| `unspliced_total`              | `total_mass`          | Already inside an `unspliced` table; redundant prefix |
| `gdna_unspliced_mean`          | `gdna_mass`           | These are point estimates of mass, not "means"      |
| `rna_unspliced_mean`           | `rna_mass`            | Same                                                |
| `support_count` (from PR 02a)  | `unspliced_counts`    | Concrete: physical unspliced fragments touching r   |
| `mass_opportunity`             | `region_size_bp`      | The denominator is region width in bp; say so       |
| `support_ess` (proposed)       | DELETED               | The integer `unspliced_counts` is the ESS           |

`method`, `precision`, and `flags` retain their names; they are not user-facing.

## What `rho0` Is, and Why It Stays

### Definition

`rho0` is the **library-wide expected gDNA fragment mass per bp**. It is a property of the
library (sequencing depth, gDNA contamination fraction, fragment-length distribution), not a
property of any particular class of regions. "Unspliced" is incidental: all gDNA fragments are
unspliced, and the unspliced fractional-mass channel is just our isolation pathway. The biology
that `rho0` summarizes is "how much sequenced gDNA mass per bp the library produces on average".

Under the new continuous-exposure framing, the regional gDNA density is

```text
rho_r = rho0 * omega_r
```

where `omega_r` is the per-region exposure factor with prior mean 1 (PR 04). So `rho0` is
exactly the scale anchor that makes `omega_r = 1` mean "average".

### Why we still need it

1. **PR 04 anchor.** The EB exposure model needs a center to shrink toward. Without `rho0` there
   is no `lambda_r = rho0 * region_size_bp[r]` to compute `raw_ratio_r = M_r / lambda_r`, so
   shrinkage has no scale.
2. **PR 04 dispersion seed.** The EB shrinkage strength `nu` in `Gamma(nu, nu)` is set from the
   variance of regional gDNA densities. PR 03 must therefore deliver not just a scalar `rho0`
   but a **dispersion** estimate computed over the same reliable-region pool. Without this,
   PR 04 has to re-estimate variance from scratch and cannot distinguish sampling noise from
   true exposure variance.
3. **Tier-3 imputation fallback for `M_r`.** When a region has no strand contrast and no
   boundary evidence, PR 03 must still emit a finite `M_r`. The only honest answer is
   `M_r = rho0 * region_size_bp[r]`. See the fallback hierarchy below.

`rho0` is **not** an off-target rate. It is **not** a classification boundary. It is **not** a
prior on capture status. In a non-capture library it equals the mean library gDNA density; in a
capture library it equals the median/robust-center of regional gDNA densities — enriched targets
are valid samples from the tail of the `omega_r` distribution, not contaminants to exclude.

### Why the bootstrap `seed_mask` is *not* the right gate any more

The old `seed_mask` (top-`t` density trim + spliced/strand-RNA exclusion + low-opportunity
exclusion) was built for the **old** goal: "find clean off-target seeds for `rho_off`". With the
new goal — estimate a global gDNA density distribution — the trim is the wrong tool:

- The trim excludes enriched capture regions, but those regions carry real gDNA mass and (when
  strand or boundary evidence is reliable) inform the *tail* of the gDNA density distribution
  that PR 04 wants to learn.
- The trim excludes high-density regions even in non-capture libraries, biasing `rho0` downward
  when natural gDNA density varies (e.g. mappability hotspots).
- The trim is an unprincipled threshold; PR 03 should use a principled robust estimator instead.

### The principled gate is *reliability of M_r*, not seed-mask membership

PR 03 drops the `seed_mask` freeze from the iterative `rho0` refit. The training pool becomes
"any region with a reliable `M_r`": Tier 1 (strand) or Tier 2 (boundary). Tier 3 imputed
regions are excluded (identity loop). Within that pool, the estimator is a **robust location +
dispersion** estimator weighted by reliability and ESS so that enriched tails do not dominate
the center. See "Iterative `rho0` Estimation" below for the formula.

The bootstrap `fit_background_model()` still produces an initial `rho0` from heuristic seeds
because we need a starting value before the first E-step has produced any `p_unexpressed` or
`M_r`. But after the first E-step the bootstrap mask is no longer used.

## Dataclass: `RegionUnsplicedMass`

Replaces `PriorMassDeconvolution` in `src/rigel/calibration/calibration_iteration.py`.

```python
@dataclass(frozen=True, slots=True)
class RegionUnsplicedMass:
    """Mass-conserving unspliced decomposition per region.

    Source of truth for downstream adaptive prior assembly (PR 03) and EB
    exposure learning (PR 04).
    """

    # Primary mass tensors (float64 throughout to avoid casting penalties
    # when projected into native EM, which consumes doubles).
    total_mass: np.ndarray          # float64[R]   T_r
    gdna_mass: np.ndarray           # float64[R]   M_r
    rna_mass: np.ndarray            # float64[R]   R_r = T_r - M_r
    region_size_bp: np.ndarray      # float64[R]   region width in bp (b - a)

    # Effective sample size for downstream uncertainty gating.
    # This is the integer count from PR 02a, exposed as the ESS field
    # the adaptive prior consumes. No squared-mass ESS is used.
    unspliced_counts: np.ndarray    # uint64[R]    N_r

    # Provenance and quality.
    method: np.ndarray              # uint8[R]     METHOD_STRAND / METHOD_BOUNDARY /
                                    #              METHOD_BACKGROUND_FALLBACK
    precision: np.ndarray           # float64[R]   reliability of the M_r estimator,
                                    #              not of the EB shrinkage
    flags: np.ndarray               # uint16[R]    region-level diagnostic bits
```

### Method constants

```python
METHOD_STRAND              = 1   # strand deconvolution drove M_r
METHOD_BOUNDARY            = 2   # boundary-sweep imputation drove M_r
METHOD_BACKGROUND_FALLBACK = 3   # no local evidence; M_r := rho0 * region_size_bp
```

### Invariants (enforced in `__post_init__`)

```text
shape(total_mass) == shape(gdna_mass) == shape(rna_mass)
                  == shape(region_size_bp) == shape(unspliced_counts)
                  == (R,)
total_mass >= 0
0 <= gdna_mass <= total_mass
0 <= rna_mass  <= total_mass
gdna_mass + rna_mass == total_mass   (float64 exact, not f32 tolerance)
region_size_bp > 0
unspliced_counts.dtype == np.uint64
precision >= 0 and finite
```

### dtype rationale

The previous design mixed `float32` for masses with `float64` for opportunity. Downstream native
EM consumes `double`. Keeping every primary tensor at `float64` avoids per-pass casting in
`prior.py` projection sweeps and removes the float32 conservation tolerance hack
(`total32 = total.astype(float32)` then re-deriving RNA) that the current code uses. Memory cost
is modest: at typical ~1M regions this is ~32 MB extra vs. f32; acceptable for the simplification.

## Computing `M_r`: Deterministic Fallback Hierarchy

The critique correctly demands an explicit, ordered cascade so the implementation is
deterministic. PR 03 commits to exactly three tiers, evaluated in order, per region.

### Tier 1 — Strand deconvolution (preferred)

Used when strand contrast is identifiable for region `r`. Identifiability is determined by the
existing `RegionGdnaChannelEstimate` reliability and precision fields, not by any new heuristic.

```text
M_r = contained_mean[r]       * clip(contained_reliability[r],       0, 1)
    + boundary_left_mean[r]   * clip(boundary_left_reliability[r],   0, 1)
    + boundary_right_mean[r]  * clip(boundary_right_reliability[r],  0, 1)
M_r = clip(M_r, 0, total_mass[r])
method[r]    = METHOD_STRAND
precision[r] = max over channels of (channel_precision * channel_reliability)
```

This is exactly what the current `build_prior_mass_deconvolution` does when `strand_channels` is
provided. PR 03 preserves the numerical behavior; only the surrounding contract changes.

### Tier 2 — Boundary-sweep imputation (when strand fails, evidence exists)

Used when strand deconvolution is not identifiable but the region has nonzero boundary evidence
(left or right flux), as reported by `BoundaryLocalPosterior` and `BoundarySweepResult`.

```text
boundary_contained = sweep.contained_gdna_mean[r]    # already conserves mass
M_r = clip(boundary_contained, 0, total_mass[r])
method[r]    = METHOD_BOUNDARY
precision[r] = sweep.contained_gdna_precision[r]
```

The previous design proposed `M_r = min(T_r, max(direct_contained_gdna_mean, boundary_sweep_mean))`.
PR 03 rejects the `max(...)` rule. Without strand contrast there is no "direct contained gDNA
mean" that is independent of the boundary sweep — the contained density estimator already folds
in the background. Taking a max of two correlated quantities biases `M_r` upward.

If strand contrast is partially available (e.g. one channel reliable, others not), Tier 1 still
fires; the unreliable channels contribute zero by construction (their reliability clips to ~0).

### Tier 3 — Background imputation (no local evidence)

Used when neither strand nor boundary evidence is informative for region `r`. Definition of "not
informative" (must hold for both):

```text
strand_channels is None for r   OR   all channel reliabilities < epsilon
local_posterior.alpha_excess[r] <= 0  AND  local_posterior.beta_excess[r] <= 0
```

Then:

```text
M_r          = clip(rho0_mean * region_size_bp[r], 0, total_mass[r])
method[r]    = METHOD_BACKGROUND_FALLBACK
precision[r] = 0.0
flags[r]    |= FLAG_M_IMPUTED_FROM_BACKGROUND
```

A new flag `FLAG_M_IMPUTED_FROM_BACKGROUND` is added. Both the `rho0` refit (this PR) and the
EB shrinkage training pool (PR 04) must exclude regions with this flag — they are an identity
loop with respect to `rho0`.

### After all tiers

```text
rna_mass[r] = total_mass[r] - gdna_mass[r]
```

Mass conservation holds exactly in `float64`.

## `region_size_bp` (the mass denominator)

The previous draft argued about whether `mass_opportunity` should be region width or
start-position count `W + L - 1`. Resolution:

- The fractional mass `total_mass[r]` is accumulated by the native fractional overlap pass. The
  natural denominator for turning that fractional mass into a density is **region width in bp**,
  `region_size_bp[r] = b - a`. A 1 bp region that accumulates fractional mass 0.05 has density
  0.05 per bp, regardless of how many fragments physically touched it.
- The start-position opportunity `W + L - 1` is the correct denominator for converting
  **integer fragment counts** into a count-based density. That denominator belongs to PR 04
  if/when it computes a count-Poisson likelihood. PR 03 does not need it.

So PR 03 stores only `region_size_bp[r] = b - a` (`float64`). No FL-marginal opportunity is
needed in this PR; if PR 06 validation shows edge effects at chromosome ends matter, that fix
ships separately and only touches a small bounded number of regions.

## Iterative `rho0` Estimation

### Output: `BackgroundDensity`

`rho0` is delivered as a small dataclass, not a scalar, so PR 04 has the dispersion it needs:

```python
@dataclass(frozen=True, slots=True)
class BackgroundDensity:
    """Library-wide gDNA density distribution.

    `rho0_mean` is the robust location of regional gDNA densities; PR 04 uses it
    as the anchor for `omega_r`. `log_dispersion` is a scale measure (standard
    deviation of log(rho_hat_r)) that PR 04 uses to initialize the EB shrinkage
    parameter `nu`. `n_effective_regions` is the sum of weights in the robust
    estimator and tells PR 04 how much to trust the dispersion.
    """

    rho0_mean: float           # robust center of regional gDNA densities (gDNA mass / bp)
    log_dispersion: float      # weighted MAD of log(rho_hat_r), in nats
    n_effective_regions: float # sum of weights (reliability * unspliced_counts)
    n_regions_in_pool: int     # raw count of regions contributing
    method_histogram: tuple    # counts of regions per M_r tier in the pool
```

The bootstrap `BackgroundModel` (`rho_off_alpha`, `rho_off_beta`, `rho_off_mean`,
`seed_mask`, ...) is preserved for the bootstrap call only. After the first E-step the
iteration carries `BackgroundDensity` and stops using the bootstrap Gamma posterior.

### Training pool (per iteration)

```text
pool_mask[r] = (method[r] in {METHOD_STRAND, METHOD_BOUNDARY})
             & (unspliced_counts[r] >= 1)
             & (region_size_bp[r] >= min_region_size)   # default: 1
```

Explicit exclusions:

- `method == METHOD_BACKGROUND_FALLBACK`: identity loop.
- `unspliced_counts == 0`: no physical observation in the region; `M_r` is a fractional artefact
  of a fragment touching a neighbor.
- No top-density tail trim. No bootstrap `seed_mask`. Enriched regions stay in the pool; the
  robust estimator (below) bounds their influence on the center.
- No spliced-evidence trim. Spliced fragments do not enter the unspliced mass channel anyway, so
  this is implicit.

### Per-region weight

```text
w[r] = precision[r] * unspliced_counts[r] * p_unexpressed[r]
```

Rationale per factor:

- `precision[r]`: how reliable the `M_r` estimator is for this region (strand reliability
  product, or boundary-sweep precision).
- `unspliced_counts[r]`: how many independent physical fragments support the region (the only
  ESS we trust).
- `p_unexpressed[r]`: as iterations refine which regions are unexpressed, the pool naturally
  shifts toward cleaner samples. Expressed regions are *not excluded* — their strand-deconv
  `M_r` is still a legitimate sample — but their weight decays smoothly as the model becomes
  more confident they are expressed.

On the first iteration `p_unexpressed[r]` may be uniform or may come from a heuristic; in either
case the formula degrades gracefully.

### Robust location and dispersion

Over `pool_mask` regions, compute per-region density and its log:

```text
rho_hat[r]      = gdna_mass[r] / region_size_bp[r]
log_rho_hat[r]  = log(max(rho_hat[r], rho_floor))     # rho_floor avoids -inf
```

The production estimator is a weighted geometric mean (equivalently, weighted mean of
`log_rho_hat`) with a Huber-style influence cap to bound tail influence:

```text
log_center_raw = weighted_mean(log_rho_hat, w)        over pool_mask
log_mad        = weighted_median(|log_rho_hat - log_center_raw|, w)
log_dispersion = 1.4826 * log_mad                     # convert MAD to sigma-equivalent

# Huberize: cap each point's deviation at k * log_dispersion, k = 1.5.
clipped = clip(log_rho_hat,
               log_center_raw - 1.5 * log_dispersion,
               log_center_raw + 1.5 * log_dispersion)
log_center = weighted_mean(clipped, w)
rho0_mean  = exp(log_center)
```

The weighted geometric mean is the right central tendency for a multiplicative scale parameter:
it is invariant under unit changes, behaves well on heavy-tailed positive data, and matches the
`Gamma`-with-log-link assumption underlying PR 04's `omega_r` prior. The Huberization step caps
the pull of enriched outliers without removing them from the pool.

Damping across iterations (carry the previous estimate to stabilize convergence):

```text
rho0_next = exp((1 - damping) * log(rho0_prev) + damping * log_center)
```

Default `damping = 0.5`.

### Why this resolves the leak critique without a seed mask

The original critique was: with the four-state classifier dead, `p_unexpressed` cannot
distinguish enriched-on-target-unexpressed from off-target-unexpressed, so a single mega-enriched
region could blow up `rho0`. The robust geometric-mean + Huberization estimator caps the
influence of any single region at a bounded multiplicative factor regardless of how massive its
`M_r` is. Enriched regions remain in the pool as legitimate samples from the tail of the gDNA
density distribution, and they contribute to `log_dispersion` (which PR 04 wants), but they
cannot single-handedly move `rho0_mean`.

If PR 06 validation surfaces a residual capture-library skew, the remedy is to tighten the Huber
cap `k` (e.g. from 1.5 to 1.0), not to re-introduce a fixed seed mask.

### Bootstrap call (`fit_background_model`)

Unchanged interface and behavior. It still runs once at the start of calibration to produce an
initial `rho0_mean` before any E-step exists. Its `seed_mask` output is no longer consumed by
the iterative refit. The bootstrap is allowed to keep its existing heuristics because it has no
better alternative on the first pass; subsequent passes use the principled estimator above.

## Downstream Handoff to Adaptive Prior

The critique correctly notes that handing only means starves `adaptive_prior.py` of variance
information. PR 03 binds both means and ESS.

### `prior.py` projection (region -> locus)

Currently `prior.py` projects per-region quantities onto loci via overlap weights and feeds the
result to `compute_adaptive_prior()`. PR 03 extends the projection signature to carry ESS:

```python
def project_region_mass_to_locus(
    region_mass: RegionUnsplicedMass,
    locus_table: LocusTable,
    overlap_weights: np.ndarray,
) -> LocusUnsplicedMass:
    """Project per-region mass and ESS onto loci."""
    return LocusUnsplicedMass(
        gdna_mass=overlap_weights @ region_mass.gdna_mass,
        rna_mass=overlap_weights @ region_mass.rna_mass,
        total_mass=overlap_weights @ region_mass.total_mass,
        # ESS projection: integer counts get a weighted-sum projection,
        # then floored at the dominant region's count to preserve discreteness.
        # Implementation detail: store as float64 and let adaptive_prior
        # interpret the value as the locus-level effective number of
        # independent unspliced observations.
        unspliced_counts=overlap_weights @ region_mass.unspliced_counts.astype(np.float64),
    )
```

### `adaptive_prior.py` consumption

`compute_adaptive_prior()` gains a new required argument `locus_ess` (the projected
`unspliced_counts`). The existing entropy/eligibility weight remains; the new ESS scales the
prior concentration:

```text
prior_concentration_l = base_concentration * shrink_toward_uniform(
    locus_ess_l,
    floor_ess,
    ceil_ess,
)
```

where `shrink_toward_uniform(ess, floor, ceil)` is a monotone increasing function in `[0, 1]`:
zero ESS contributes a flat prior (no signal), and high ESS contributes the full
data-driven Dirichlet. The exact functional form is left to the adaptive prior implementer, but
PR 03 mandates the **plumbing**: ESS must be a first-class argument, not derived inside
`adaptive_prior.py` from the means.

Acceptance test (PR 03 scope): two synthetic loci with identical projected `gdna_mass` and
`rna_mass` but different `unspliced_counts` (e.g. 5 vs. 500) must produce different prior
concentrations from `compute_adaptive_prior()`.

## Files to Touch

### `src/rigel/calibration/calibration_iteration.py`

- Delete `PriorMassDeconvolution`. Introduce `RegionUnsplicedMass` and `BackgroundDensity` per
  the dataclasses above.
- Update `RegionCalibration` and `CalibrationStepResult`: rename `prior_mass` field to
  `region_unspliced_mass`; replace scalar `rho_off` field with a `BackgroundDensity` field
  named `background_density`. The old scalar `rho_off` accessor may be retained as a property
  returning `background_density.rho0_mean` for one PR cycle if downstream churn is high; the
  preferred path is to migrate consumers to the dataclass immediately.
- Rewrite `build_prior_mass_deconvolution(...)` as
  `build_region_unspliced_mass(observation, *, region_size_bp, strand_channels, local_posterior,
  sweep, background_density, unspliced_counts) -> RegionUnsplicedMass`. The new signature
  takes the three estimator inputs and the background scalar so the fallback hierarchy is
  implementable in one place.
- Replace `calibration_m_step()` background refit with the robust geometric-mean estimator:
  - New function `estimate_background_density(region_unspliced_mass, p_unexpressed, *,
    previous: BackgroundDensity, damping=0.5, huber_k=1.5, rho_floor=1e-12) ->
    BackgroundDensity`.
  - Pool: `method in {STRAND, BOUNDARY}` AND `unspliced_counts >= 1`.
  - Weights: `precision * unspliced_counts * p_unexpressed`.
  - Compute weighted geometric mean of `rho_hat[r] = gdna_mass[r] / region_size_bp[r]` with
    Huber-capped influence; compute `log_dispersion` from weighted MAD.
  - Damp in log space against `previous.rho0_mean`.
- `calibration_e_step()` plumbs `region_size_bp` (from `RegionArrays`) and `unspliced_counts`
  (from PR 02a's `RegionCountLedger`) and the current `BackgroundDensity` into
  `build_region_unspliced_mass()`.

### `src/rigel/calibration/_arrays.py` / `region_count_ledger.py`

- Add `region_size_bp` as a derived field on `RegionArrays` if not already present. It is
  `region_end - region_start`, sorted by the same `order` permutation.
- Confirm `unspliced_counts` (from PR 02a) is reachable from the same sort order.

### `src/rigel/calibration/prior.py`

- Update to consume `RegionUnsplicedMass` instead of `PriorMassDeconvolution`.
- Add ESS projection per the handoff section.
- Pass `locus_ess` into `compute_adaptive_prior()`.

### `src/rigel/calibration/adaptive_prior.py`

- Accept and use `locus_ess`. Document the shrinkage function and its bounds.
- Existing `p_unexpressed` soft gate is unchanged.

### `src/rigel/calibration/_result.py`

- Region/locus summary fields renamed to match `RegionUnsplicedMass` field names.
- Add a `method` histogram (counts of regions per tier) to the diagnostics block.

### `src/rigel/calibration/background_model.py`

- No API change for the bootstrap `fit_background_model()` call; its `seed_mask` and Gamma
  posterior are still produced for the first pass.
- Add a docstring note that `seed_mask` is **only** used to initialize `rho0` before the first
  E-step. Subsequent iterations use `estimate_background_density()` from
  `calibration_iteration.py` and ignore `seed_mask`.
- Provide a small helper `BackgroundDensity.from_bootstrap(background_model) ->
  BackgroundDensity` that wraps the bootstrap Gamma posterior mean as `rho0_mean` and
  initializes `log_dispersion` to a wide default (e.g. `log(10.0)`) and
  `n_effective_regions = 0.0`, signaling "no data-driven dispersion yet" so the first E-step
  treats the dispersion as a flat prior.

## Tests

Add `tests/test_region_unspliced_mass.py` covering:

1. **Mass conservation.** `gdna_mass + rna_mass == total_mass` exactly in float64 for synthetic
   inputs that touch all three tiers.
2. **Tier 1 (strand).** Synthetic region with strong antisense unspliced flux: `M_r` matches the
   strand deconvolution prediction; `method == METHOD_STRAND`; `precision > 0`.
3. **Tier 2 (boundary).** Synthetic region with no strand contrast but strong boundary excess:
   `M_r` matches `sweep.contained_gdna_mean`; `method == METHOD_BOUNDARY`.
4. **Tier 3 (background fallback).** Synthetic region with neither strand contrast nor boundary
   evidence: `M_r == background_density.rho0_mean * region_size_bp`;
   `method == METHOD_BACKGROUND_FALLBACK`; `flags & FLAG_M_IMPUTED_FROM_BACKGROUND`.
5. **Sentinel: no-gDNA, no-RNA region.** Zero `total_mass`, zero `unspliced_counts`: all three
   tiers must produce `M_r == 0` and `R_r == 0`.
6. **Tier 3 exclusion from `rho0` update.** Construct a `RegionUnsplicedMass` with N tier-3
   regions and a single tier-1 region. Run `estimate_background_density()` once; assert that
   `rho0_next` equals what you would get from the tier-1 region alone. Locks the identity-loop
   fix.
7. **Robust estimator caps enriched-tail influence.** Construct a pool with 1000 tier-1
   regions at density `1.0` and 10 tier-1 regions at density `1000.0` (mimicking capture
   enrichment). Assert that `rho0_mean` lies within `[0.9, 1.3]` (i.e. close to the bulk mode,
   not the contaminated mean of ~10.9). Locks the no-seed-mask robustness claim.
8. **Dispersion grows with regional spread.** Two pools with identical `rho0_mean` but
   different multiplicative spread of `rho_hat_r` produce different `log_dispersion` values,
   monotonically increasing with the input spread. Locks the dispersion contract for PR 04.
9. **Iterative refinement uses `p_unexpressed`.** Run two iterations: after the second
   iteration, regions with `p_unexpressed -> 1` should contribute more weight than they did in
   iteration 1, and `rho0_mean` should move accordingly.
10. **Damping stabilizes.** Apply a synthetic shock (replace `rho_hat_r` by a 10x scaled
    version mid-iteration); confirm `rho0_mean` moves by approximately `damping` of the shock
    on the first pass, not the full shock.
11. **dtype.** Every primary tensor in `RegionUnsplicedMass` is exactly `np.float64` (or
    `np.uint64` for `unspliced_counts`); `BackgroundDensity` scalars are Python `float`.
12. **Adaptive prior receives ESS.** Two synthetic loci with identical projected `gdna_mass`
    and `rna_mass` but different `unspliced_counts` produce different prior concentrations.
13. **Bootstrap handoff.** `BackgroundDensity.from_bootstrap(...)` yields
    `rho0_mean == background_model.rho_off_mean` and `n_effective_regions == 0.0`.

Update existing tests that referenced `PriorMassDeconvolution` field names:

- `tests/test_calibration_iteration.py`
- `tests/test_calibration_prior.py`
- `tests/test_calibration_result.py`
- `tests/test_per_locus_gdna_mass.py`
- `tests/test_bayesian_prior_acceptance.py`

These changes are purely rename + dtype + the addition of an `unspliced_counts` argument; no
semantic test redesign is expected outside of the M-step freeze test (#7) and the ESS plumbing
test (#9).

## Targeted Validation Commands

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

No native rebuild is required for PR 03. PR 02a already exposes `unspliced_counts`.

## Resolved Critiques (Checklist)

- [x] **`rho0` redefined.** Library-wide gDNA fragment mass per bp; a property of the library,
      not of "unexpressed regions". Estimated from any reliable `M_r` region.
- [x] **`rho0` as a distribution.** Delivered as `BackgroundDensity(rho0_mean,
      log_dispersion, n_effective_regions, ...)` so PR 04 has the dispersion it needs.
- [x] **Bootstrap `seed_mask` retired from iterative refit.** Used only to seed the first
      `rho0`. Iterative refit uses reliability-weighted robust geometric mean over
      strand/boundary-tier regions.
- [x] **M-step background leak resolved by robust estimator.** Weighted geometric mean +
      Huberization caps tail influence so enriched capture targets cannot dominate the center;
      they still contribute to `log_dispersion`, which PR 04 wants.
- [x] **Identity loop closed.** Tier 3 imputed regions and zero-count regions are excluded from
      the `rho0` training pool.
- [x] **Iterative pool refinement.** Weighting by `p_unexpressed` lets the pool refine itself
      across calibration passes without hard-classifying regions.
- [x] **Downstream uncertainty starvation.** `RegionUnsplicedMass.unspliced_counts` is bound to
      the adaptive prior handoff via a `locus_ess` argument to `compute_adaptive_prior()`.
- [x] **Unstranded fallback hierarchy.** Three explicit tiers in fixed order: strand -> boundary
      sweep -> background imputation. No `max(strand, boundary)` combine.
- [x] **Field renames for biological intuition.** `RegionUnsplicedMass`, `total_mass`,
      `gdna_mass`, `rna_mass`, `region_size_bp`, `unspliced_counts`, `BackgroundDensity`.
- [x] **No squared-mass ESS.** PR 02a's integer `unspliced_counts` is the only ESS source.
- [x] **Uniform float64 across primary tensors.** Removes per-pass casting in `prior.py` and
      removes the float32 mass-conservation tolerance hack.

## Done Means

- One production regional unspliced mass table (`RegionUnsplicedMass`) consumed by both adaptive
  priors (PR 03 wiring) and EB exposure learning (PR 04, future work).
- `M_r` is computed by a deterministic three-tier fallback. The selected tier is recorded per
  region in `method`.
- `rho0` is delivered as `BackgroundDensity(rho0_mean, log_dispersion, n_effective_regions,
  ...)`. It is refit each iteration by a reliability-weighted robust geometric-mean estimator
  over strand/boundary-tier regions. The bootstrap `seed_mask` is used only to seed iteration 1.
- Adaptive prior accepts a per-locus ESS argument derived from `unspliced_counts`.
- All primary mass tensors are `float64`; mass conservation is exact (not tolerance-based).
- No capture-state field appears anywhere in PR 03 inputs or outputs.
- Enriched capture regions remain in the `rho0` estimation pool as legitimate tail samples;
  their influence on `rho0_mean` is bounded by the Huber cap, and they contribute to
  `log_dispersion` for PR 04.

## Implementation Log

(To be filled in as the PR lands.)
