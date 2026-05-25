# Rigel Calibration Fusion and EM Integration Implementation Plan (v2)

Date: 2026-05-25
Status: implementation plan
Supersedes: `docs/calibration/calibration_fusion_impl_plan_v1.md`
Uses as historical context only: older `docs/fineregions/calibration_roadmap_v*.md`

## 0. Purpose

v1 identified the right endpoint: Rigel needs to fuse regional density evidence with
local strand evidence, allocate the resulting gDNA prior to MultiLocus EM units, and
remove the final `quant_from_buffer` boundary so `rigel quant` is runnable end to end.

v2 keeps that endpoint, but rewrites the implementation contract against the code that
exists now. The central principle is simple:

```text
calibration owns regional evidence;
prior allocation owns MultiLocus projection;
pipeline owns orchestration only.
```

The implementation should feel boring in the best possible way. Existing modules keep
owning their domain logic. New modules are small, explicit adapters. No hidden state, no
parallel model stack, no policy encoded as duplicate flags.

## 1. Current Code Inventory

The following pieces are live and should be incorporated, not replaced wholesale.

| Area | Current owner | v2 action |
| --- | --- | --- |
| Sorted region/count views | `src/rigel/calibration/_arrays.py` | Keep. Use `RegionArrays` and `PayloadArrays` as the common regional coordinate system. |
| Region count channels | `src/rigel/calibration/region_count_ledger.py` | Keep. Use `RegionCountLedger`; do not create another ledger class unless this module is renamed in a separate cleanup. |
| Density observations | `src/rigel/calibration/density_observation.py` | Keep. Extend only if fusion needs additional explicit count/opportunity terms. |
| Density model | `src/rigel/calibration/density_model.py` | Extend with public predictive/log-PMF helpers. Do not reimplement NB math inside integration. |
| Strand deconvolution | `src/rigel/calibration/strand_deconv.py` | Keep. Refactor exact convolution into a reusable likelihood helper. Do not duplicate `_exact_posterior_R` logic. |
| Exposure surface | `src/rigel/calibration/exposure.py` | Keep. `RegionExposure.from_density` remains the first exposure surface. |
| FL-aware geometry | `src/rigel/calibration/_exposure.py` | Keep. Use `gdna_eff_len_for_loci` and `bp_weighted_mean_exposure_over_blocks` for EM denominators. |
| Calibration orchestration | `src/rigel/calibration/_orchestrator.py` | Extend. Fusion belongs here because this function already has the regional arrays, payload arrays, ledger, density evidence, and strand counts. |
| Calibration result | `src/rigel/calibration/_result.py` | Extend with compact optional fusion/prior products and summaries. |
| Locus graph | `src/rigel/locus.py` | Keep. `build_multi_loci` and `MultiLocus.loci` are the correct footprint primitives. |
| EM partitioning | `src/rigel/locus_partition.py` | Use `partition_and_free`. Do not invent `pr.build_partitions`. |
| Native EM | `src/rigel/native/em_solver.cpp` | Keep. It already accepts per-locus prior counts, gDNA effective lengths, and explicit gDNA eligibility. |
| Pipeline helpers | `src/rigel/pipeline.py` | Restore the existing helper flow in `quant_from_buffer`; do not bypass `_setup_geometry_and_estimator`, `_score_fragments`, or `_run_locus_em_partitioned`. |

## 2. What v2 Changes From v1

v1 remains the basis, but these fixes are required before coding:

1. There is no `calibration.calibration_context` field. v2 removes it. Fusion is run
   inside `calibrate()` while the concrete regional context is still in scope.
2. `DensityEvidence` does not currently expose the NB predictive parameters that v1
   assumes. v2 adds a public density predictive surface or helper methods to
   `density_model.py`.
3. `RegionGdnaEstimate` is a summary, not the full strand likelihood input. v2 fuses
   from `StrandRegionCounts` plus `kappa_d` and `p_r1_sense`, while still preserving
   `RegionGdnaEstimate` for summaries and backward-compatible diagnostics.
4. `TS_NONE` and `TS_AMBIG` do not get fake strand deconvolution. Their strand
   likelihood is neutral; density remains active.
5. The exact posterior threshold follows the current implementation: `N <= 200`.
6. `PriorTable` distinguishes expected counts, upper counts, and EM-facing prior
   pseudocounts. Native EM receives the EM-facing array only.
7. Hotspot-safe count allocation cannot be promised by base-pair overlap alone. v2
  starts with a simple conservation-aware geometry allocator plus mandatory
  diagnostics, and keeps midpoint/span-based allocation as an optional follow-up if
  diagnostics prove the smearing risk is material.
8. Pipeline wiring must first score fragments and build `em_data`; v1's snippet used
   `em_data` and `estimator` before constructing them.
9. The ignored tests in `tests/conftest.py` are part of the work. Passing by skipping
   the broken boundary is not a finished integration.

## 3. Design Locks

These rules should be treated as implementation invariants.

### 3.1 Model Ownership

- Density predictive math lives in `density_model.py`.
- Strand likelihood math lives in `strand_deconv.py`.
- Fusion lives in `integration.py` and only combines public density and strand APIs.
- MultiLocus allocation lives in `prior.py`.
- Pipeline code wires objects together but does not perform statistical work.

### 3.2 Density Is Predictive, Not Just A Point Estimate

Fusion needs a distribution over local gDNA count `D`, not only
`mean_unbounded` and `upper_unbounded`.

Add a public helper in `density_model.py` along these lines:

```python
def density_logpmf_grid(
    evidence: DensityEvidence,
    region_idx: int,
    d_grid: np.ndarray,
) -> np.ndarray:
    """Return log P_density(D=d) for one region and integer d_grid."""
```

The helper is the only code that knows how the current Gamma-rate posterior maps to an
NB predictive count distribution. `integration.py` must not reconstruct private
`alpha_post`, `beta_post`, or `p_nb` formulas ad hoc.

To make that helper possible, extend `DensityEvidence` with the minimal predictive
state it needs. Candidate fields:

```python
alpha_post: np.ndarray          # float64[R] or float32[R]
beta_post: np.ndarray           # float64[R] or float32[R]
contained_leff: np.ndarray      # float64[R] or float32[R]
boundary_count: np.ndarray      # float32[R]
variance_unbounded: np.ndarray  # float64[R] or float32[R]
tail_probability: np.ndarray    # float32[R]
expected_tail_count: np.ndarray # float32[R]
```

If memory pressure argues against storing all of these, keep `alpha_post`, `beta_post`,
`contained_leff`, and `boundary_count`, then derive the rest through methods. The
important requirement is that fusion uses a public density API, not duplicated private
math.

### 3.3 Strand Likelihood Is Reused, Not Rewritten

`strand_deconv.py` already implements the exact convolution in `_exact_posterior_R`.
Refactor that code into a reusable helper that returns likelihood over `D`:

```python
def strand_log_likelihood_d_grid(
    k_sense_obs: int,
    n_total: int,
    d_grid: np.ndarray,
    *,
    kappa_d: float,
    p_r1_sense: float,
) -> np.ndarray:
    """Return log L_strand(K | D=d, N) for integer d_grid."""
```

Then `_exact_posterior_R` can call the helper with a uniform prior for its existing
strand-only summary path. This gives one source of truth for the convolution.

### 3.4 Applicability Is Derived

Do not make region-type policy flags that duplicate `signature` or `ts_class`.

Strand is applicable only when all of the following hold:

```text
ts_class in {TS_POS, TS_NEG}
observed strand-folded count N > 0
abs(p_r1_sense - 0.5) >= STRAND_CONTRAST_NUMERICAL_FLOOR
```

`TS_NONE`, `TS_AMBIG`, zero-count regions, and near-unstranded libraries get a neutral
strand likelihood. Neutral means:

```text
log L_strand(D=d) = constant over d in [0, N]
```

Density remains active in all of those cases.

Existing strand flags such as `FLAG_INELIGIBLE`, `FLAG_NEAR_UNSTRANDED`, and
`FLAG_APPROX_NORMAL` remain diagnostics. They are not the primary fusion API.

### 3.5 Bounded Local Posterior, Explicit Tail Diagnostics

Density can predict more gDNA than was locally observed. EM cannot receive more local
prior count than the local observed compatible ledger supports without making the prior
look like observed data.

Fusion therefore produces two concepts:

```text
unbounded density prediction: model diagnostic and exposure-learning input
bounded fused posterior: EM-facing local count distribution over 0 <= D <= N
```

The bounded posterior is used for local prior allocation. Tail metrics are reported:

```text
tail_probability = P_density(D > N)
expected_tail_count = E_density[(D - N)+]
density_over_observed_ratio = mean_unbounded / max(N, eps)
```

### 3.6 EM Prior Pseudocount Semantics

Native EM treats `locus_gdna_prior_count` as an additive Dirichlet pseudocount for the
single gDNA component of a MultiLocus. A pseudocount should default to an expectation,
not a high-confidence upper bound.

v2 default:

```text
gdna_prior_count_em = allocated fused mean_count
```

Also record:

```text
gdna_upper_count = allocated fused upper_count
gdna_expected_count = allocated fused mean_count
```

If benchmarks later show systematic false RNA calls from underpowered priors, add a
configuration knob for `prior_count_mode = mean | upper | blend`. Do not silently bake
that policy into the first implementation.

### 3.7 gDNA Eligibility Is Decoupled From Prior Strength

A zero prior count must not disable the gDNA component. The native solver already has a
separate `enable_gdna` array. `assemble_priors` should compute it from EM data:

```text
enable_gdna[locus] = any unit in locus is unspliced and has finite gdna_log_lik
```

This preserves the Bayesian-prior redesign invariant and avoids resurrecting older
behavior where no prior meant no gDNA hypothesis.

### 3.8 Count Allocation Is Diagnostic-First

The first runnable implementation should not require new native span storage or a
complex interval engine before EM can run. The conservative v2 policy is:

1. Allocate regional count mass by normalized MultiLocus footprint overlap as the
   default MVP.
2. Record exactly how much fused mass came from regions touched by multiple
   MultiLoci, regions with partial footprint coverage, and loci where geometry-only
   allocation could smear a localized hotspot.
3. Use existing `ScoredFragments.genomic_midpoint` only for diagnostics or an optional
   allocation mode until benchmarks show the geometry allocator is insufficient.
4. Defer actual genomic fragment spans `(start, end)` unless midpoint diagnostics show
   that boundary-crossing fragments are common enough to justify a native payload
   change.

This intentionally favors a runnable, inspectable system. If the geometry allocator
over-allocates gDNA to a nearby locus, the likely first-order effect is conservative:
more local gDNA prior, less false RNA. That can still hurt sensitivity, so the
diagnostics are not optional. Base-pair overlap remains appropriate for opportunity
denominators regardless of which observed-count allocator is selected.

## 4. New And Extended Schemas

### 4.1 FusedRegionGdnaEvidence

New module: `src/rigel/calibration/integration.py`

```python
@dataclass(frozen=True, slots=True)
class FusedRegionGdnaEvidence:
    mean_count: np.ndarray                 # float32[R], E[D | density, strand, 0 <= D <= N]
    upper_count: np.ndarray                # float32[R], confidence quantile inside [0, N]
    variance_count: np.ndarray             # float32[R]
    rna_lower_count: np.ndarray            # float32[R], max(N - upper_count, 0)
    observed_compatible_count: np.ndarray  # float32[R]

    density_weight: np.ndarray             # float32[R], diagnostic information weight
    strand_weight: np.ndarray              # float32[R], diagnostic information weight
    density_applicable: np.ndarray         # bool[R]
    strand_applicable: np.ndarray          # bool[R]

    tail_probability: np.ndarray           # float32[R], density tail beyond observed bound
    expected_tail_count: np.ndarray        # float32[R]
    flags: np.ndarray                      # uint8[R]
```

Suggested flag bits:

```python
FUSED_DENSITY_ONLY = 1 << 0
FUSED_STRAND_USED = 1 << 1
FUSED_EXACT = 1 << 2
FUSED_APPROX = 1 << 3
FUSED_DENSITY_TAIL = 1 << 4
FUSED_DEGENERATE = 1 << 5
FUSED_FALLBACK = 1 << 6
FUSED_BOUNDARY_FALLBACK = 1 << 7
```

### 4.2 PriorTable

New module: `src/rigel/calibration/prior.py`

```python
@dataclass(frozen=True, slots=True)
class PriorTable:
    gdna_prior_count_em: np.ndarray        # float64[L], passed to native EM
    gdna_expected_count: np.ndarray        # float64[L], allocated fused mean
    gdna_upper_count: np.ndarray           # float64[L], allocated fused upper

    gdna_eff_len: np.ndarray               # float64[L], exposure-weighted EM denominator
    gdna_eff_len_unweighted: np.ndarray    # float64[L], raw FL-marginal denominator
    gdna_em_exposure_weight: np.ndarray    # float64[L], gdna_eff_len / unweighted
    enable_gdna: np.ndarray                # uint8[L] or bool[L]

    n_regions_touched: np.ndarray          # int32[L]
    n_units_used_for_diagnostics: np.ndarray # int64[L]
    count_allocation_mode: np.ndarray      # uint8[L], geometry / midpoint / span
    count_allocation_fallback: np.ndarray  # uint8[L], diagnostic fallback code
    multi_locus_region_mass: np.ndarray    # float64[L], mass from regions touching >1 locus
    partial_coverage_region_mass: np.ndarray # float64[L], mass from partially covered regions
```

Keep per-region exposure factors on `RegionExposure`, not on `PriorTable`. `PriorTable`
is a per-MultiLocus product.

### 4.3 CalibrationResult Additions

Extend `CalibrationResult` with optional products:

```python
fused_region_gdna: FusedRegionGdnaEvidence | None = None
prior_table: PriorTable | None = None
```

`to_summary_dict()` should include compact summaries only: histograms, totals,
quantiles, fallback counts, and `n_multi_loci`. Do not serialize large arrays into
`summary.json`.

If we later want regional calibration tables for debugging, write them as an explicit
optional artifact such as `calibration_regions.feather`. Do not hide multi-megabyte
arrays in JSON.

## 5. Fusion Algorithm

### 5.1 Inputs

Fusion is called from `_orchestrator.calibrate()` after these objects already exist:

```python
region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
ledger = build_region_count_ledger(payload_arrays)
observation = build_density_observation(region_arrays, ledger, fl_models.gdna)
density_evidence = fit_density_evidence(observation, ...)
strand_counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=...)
region_gdna = deconvolve_regions_by_strand(strand_counts, ...)
```

New call:

```python
fused_region_gdna = fuse_density_and_strand(
    region_arrays=region_arrays,
    ledger=ledger,
    density_observation=observation,
    density_evidence=density_evidence,
    strand_counts=strand_counts,
    strand_summary=strand_summary,
    kappa_d=kappa_d.kappa,
    confidence=rna_lower_confidence,
)
```

`region_gdna` remains in the result for strand-only diagnostics. Fusion does not depend
on `RegionGdnaEstimate.mean_count` as its sole input.

### 5.2 Exact Path For N <= 200

For each region:

1. Compute `N = round(observed_compatible_count)` from the unspliced compatible ledger.
2. If `N == 0`, emit zeros and `FUSED_DEGENERATE`.
3. Build `d_grid = np.arange(N + 1)`.
4. Get `log_density = density_logpmf_grid(density_evidence, region_idx, d_grid)`.
5. If strand is applicable, get `log_strand = strand_log_likelihood_d_grid(...)`.
6. If strand is neutral, initialize an explicit zero vector:

  ```python
  log_strand = np.zeros(N + 1, dtype=log_density.dtype)
  ```

  This is the additive log-likelihood for multiplying by `1.0`; it must not be
  represented as `None`, `-inf`, or a skipped normalization branch.
7. Normalize `log_density + log_strand` with log-sum-exp.
8. Emit mean, variance, upper quantile, `rna_lower_count`, and information weights.

Information weights are diagnostics, not model inputs. Start with curvature or variance
based weights where stable:

```text
w_density = tau_density / (tau_density + tau_strand)
w_strand = tau_strand / (tau_density + tau_strand)
```

If one channel is neutral, the other gets weight 1. If both are degenerate, both get 0.

### 5.3 Large-Count Path For N > 200

Use a bounded Laplace/truncated-normal approximation to the same fused log posterior:

```text
ell(d) = log_density(d) + log_strand(d), 0 <= d <= N
```

Implementation requirements:

- Optimize on `[0, N]` with a deterministic bounded scalar optimizer.
- Treat near-boundary modes as a separate regime. If the bounded mode satisfies
  `d_hat < 5` or `d_hat > N - 5`, do not rely on a symmetric Gaussian Hessian as the
  primary approximation.
- Estimate `-ell''(d_hat)` by finite differences on the count scale only for interior
  modes.
- Use `scipy.stats.truncnorm` or equivalent closed-form moments for mean, variance,
  and upper quantile.
- Set `FUSED_APPROX`.
- Validate approximation against exact fusion near the threshold and on downsampled
  fixtures.

Boundary and numerical fallbacks:

```text
if d_hat < 5 or d_hat > N - 5:
  if N <= ADAPTIVE_EXACT_MAX:
    evaluate the exact discrete posterior over 0..N
  else:
    evaluate a one-sided boundary window and approximate the remaining tail with
    an exponential/Gamma-like tail on the count scale
  set FUSED_BOUNDARY_FALLBACK

elif Hessian is nonpositive, near-zero, nonfinite, or yields invalid variance:
  expand the local evaluation window or fall back to exact when feasible
  otherwise use density-only bounded posterior for neutral/weak strand cases
  set FUSED_FALLBACK
```

`ADAPTIVE_EXACT_MAX` should be centralized in `integration.py` and tested. A first
reasonable value is `1000`, but the constant should be chosen from timing data rather
than scattered through call sites. The one-sided tail fallback is especially important
for strong strand-specific RNA regions where the posterior is pushed hard against
`D = 0`: a symmetric Hessian can become nonpositive, nearly flat on one side, or
numerically dominated by finite-difference noise. The fallback must be visible in
summary diagnostics and must never produce counts outside `[0, N]`.

## 6. Prior Allocation Algorithm

### 6.1 Function Contract

`assemble_priors` is called after scoring and MultiLocus construction:

```python
prior_table = assemble_priors(
    multi_loci=multi_loci,
    em_data=em_data,
    index=index,
    calibration=calibration,
)
```

It consumes:

- `calibration.fused_region_gdna`
- `calibration.region_exposure`
- `calibration.fl_models.gdna`
- `index.region_df` and `index.ref_name_to_id`
- `em_data.genomic_midpoint`, `em_data.is_spliced`, `em_data.gdna_log_liks`, and
  `MultiLocus.unit_indices`

### 6.2 gDNA Eligibility

For each MultiLocus:

```text
enable_gdna = any unit is unspliced and finite(gdna_log_lik)
```

This value is independent of `gdna_prior_count_em`.

### 6.3 Count Mass Allocation

The MVP allocation is geometry-based and conservation-aware. For each region `r`, find
MultiLocus footprints overlapping the region and compute:

```text
raw_share(r, L) = overlap_bp(region r, footprint L) / length(region r)
share(r, L) = raw_share(r, L) / sum_L raw_share(r, L)
```

Then:

```text
gdna_expected_count[L] += share(r, L) * fused.mean_count[r]
gdna_upper_count[L] += share(r, L) * fused.upper_count[r]
gdna_prior_count_em[L] = gdna_expected_count[L]
```

If `sum_L raw_share(r, L) == 0`, do not allocate that region's fused count mass to EM.
Record it as unallocated calibration mass.

The implementation must also compute diagnostics that tell us whether this simple rule
is good enough:

```text
n_regions_touching_multiple_loci
fused_mass_touching_multiple_loci
n_regions_partially_covered_by_loci
fused_mass_partially_covered_by_loci
n_loci_receiving_geometry_only_mass
top_loci_by_geometry_only_mass
```

This is the default because it is simple, deterministic, and uses no new native data.
It may overestimate local gDNA in rare hotspot-smearing cases; that is conservative for
false-RNA control, but the diagnostic mass must be reported so we can decide from real
data whether a more precise allocator is needed.

### 6.4 Optional Midpoint Allocator And Interval Indexing

If diagnostics show geometry-only mass is common or damaging, add an optional midpoint
allocator using the existing `em_data.genomic_midpoint` field. Do not use naive nested
loops over units and regions.

Per reference, region rows are sorted and non-overlapping. Midpoint-to-region mapping
should use binary search:

```python
local_idx = np.searchsorted(region_starts_for_ref, midpoint, side="right") - 1
valid = (local_idx >= 0) & (midpoint < region_ends_for_ref[local_idx])
```

This is `O(U log R_ref)` and can become `O(U + R_ref)` if unit midpoints are sorted per
reference before scanning. It is acceptable for diagnostics and likely acceptable for a
future allocation mode. It should be implemented in `prior.py` as vectorized per-ref
lookup, not Python nested loops.

Should we store midpoint or full genomic span?

- **Now:** keep the existing midpoint. It is already available, cheap, and sufficient to
  measure whether geometry-only allocation is risky.
- **Later:** store fragment `(start, end)` spans only if midpoint diagnostics show that
  boundary-crossing or within-region hotspot allocation changes important benchmark
  outcomes. Span storage requires native scanner/scorer changes and should not block
  the first runnable integration.

Should we use `cgranges`?

- **Not for midpoint lookup.** Sorted non-overlapping region bins plus
  `np.searchsorted` are simpler and faster to reason about.
- **Yes, possibly later for span overlap.** If we add full fragment spans or compact
  within-region bins, `cgranges` is the right tool for interval overlaps. That should be
  a focused allocator PR with diagnostics comparing geometry, midpoint, and span modes.

### 6.5 Effective Length Denominators

Use existing FL-aware helpers:

```python
gdna_eff_len_unweighted[L] = gdna_eff_len_for_loci(
    multi_locus.loci,
    index.ref_lengths,
    calibration.fl_models.gdna,
)

weight[L] = bp_weighted_mean_exposure_over_blocks(
    blocks=[(loc.ref_id, loc.start, loc.end) for loc in multi_locus.loci],
    region_arrays=region_arrays,
    exposure=calibration.region_exposure,
)

gdna_eff_len[L] = max(gdna_eff_len_unweighted[L] * weight[L], 1.0)
```

This matches the current production simplification recorded in the repository memory:
exact weighted denominators were evaluated previously and did not materially change the
known FLG2 denominator gap. If we later replace this with exact weighted denominators,
that should be a separate denominator PR with its own benchmarks.

### 6.6 Conservation Diagnostics

`assemble_priors` should report:

```text
sum_region_fused_mean
sum_locus_allocated_expected
sum_region_fused_upper
sum_locus_allocated_upper
n_regions_with_no_locus_units
n_loci_geometry_only_fallback
n_loci_enable_gdna_true
fused_mass_touching_multiple_loci
fused_mass_partially_covered_by_loci
```

The allocated total should be <= the regional fused total. It may be lower because some
regional calibration mass has no EM locus with eligible gDNA units.

## 7. Pipeline Wiring

`quant_from_buffer` should become a short, readable conduit.

Correct sequence:

```python
if calibration is None:
    raise ValueError(...)

em_config = em_config or EMConfig()
scoring_cfg = scoring or FragmentScoringConfig()

_geometry, estimator = _setup_geometry_and_estimator(
    index,
    calibration.fl_models.rna,
    em_config,
)

em_data = _score_fragments(
    buffer,
    index,
    strand_models,
    calibration.fl_models.rna,
    calibration.fl_models.gdna,
    stats,
    estimator,
    scoring_cfg,
    log_every,
    annotations,
)

multi_loci = build_multi_loci(em_data, index)
_assign_locus_ids(estimator, multi_loci)

prior_table = assemble_priors(
    multi_loci=multi_loci,
    em_data=em_data,
    index=index,
    calibration=calibration,
)

partitions = partition_and_free(em_data, multi_loci)

_run_locus_em_partitioned(
    estimator,
    partitions,
    multi_loci,
    index,
    gdna_prior_count_em=prior_table.gdna_prior_count_em,
    gdna_eff_len=prior_table.gdna_eff_len,
    enable_gdna=prior_table.enable_gdna,
    gdna_eff_len_unweighted=prior_table.gdna_eff_len_unweighted,
    gdna_prior_count_raw=prior_table.gdna_expected_count,
    gdna_em_exposure_weight=prior_table.gdna_em_exposure_weight,
    em_config=em_config,
    annotations=annotations,
    emit_locus_stats=emit_locus_stats,
)

calibration = dataclasses.replace(
    calibration,
    prior_table=prior_table,
    n_multi_loci=len(multi_loci),
)

return estimator, calibration
```

Edge cases:

- If `em_data.n_units == 0`, return the estimator with deterministic counts and an empty
  `PriorTable`; no native EM call is needed.
- If `calibration.fused_region_gdna is None`, raise a clear error. Do not silently fall
  back to strand-only priors in production wiring.
- If `index.region_df is None`, fail with the existing rebuild-index guidance.

## 8. Implementation Playbook

### PR 0: Contract Cleanup And Test Skeletons

Goal: make the work executable before changing model behavior.

Tasks:

1. Land this v2 plan.
2. Add empty or minimal tests for fusion and prior contracts.
3. Add TODO markers in `tests/conftest.py` for each ignored legacy test that must be
   rewritten before Phase 6 is considered done.
4. Confirm public signatures and names.

Acceptance:

```bash
conda activate rigel && pytest tests/test_density_model.py tests/test_strand_deconv.py -v
```

No behavior change expected.

### PR 1: Public Density And Strand Likelihood APIs

Goal: expose the model distributions without duplicating math.

Tasks:

1. Extend `DensityEvidence` or add a `DensityPredictiveSurface` so
   `density_logpmf_grid(...)` is possible.
2. Move density tail metrics that are currently transient into the summary surface.
3. Refactor `strand_deconv.py` to expose `strand_log_likelihood_d_grid(...)`.
4. Keep existing `deconvolve_regions_by_strand(...)` behavior and tests green.

Tests:

- Density log-PMF matches SciPy NB/Poisson sentinels.
- Density helper handles deterministic zero evidence.
- Strand log-likelihood reproduces `_exact_posterior_R` after applying a uniform prior.
- Perfect and near-unstranded strand cases remain stable.

### PR 2: Fusion Engine In Calibration

Goal: add `integration.py` and run fusion inside `calibrate()`.

Tasks:

1. Implement `FusedRegionGdnaEvidence`.
2. Implement exact fusion for `N <= 200`.
3. Implement bounded large-count approximation for interior modes with explicit
  boundary-mode fallback when `d_hat < 5` or `d_hat > N - 5`.
4. Add `fused_region_gdna` to `CalibrationResult` and compact summaries.
5. Update `_orchestrator.calibrate()` to build the fused surface immediately after
   strand deconvolution.

Tests:

- Unstranded input equals bounded density posterior.
- Strong strand-specific pure RNA drives fused gDNA near zero even when density prior is
  high.
- Flat density prior recovers strand-only behavior.
- `TS_NONE` and `TS_AMBIG` use neutral strand likelihood.
- Neutral strand likelihood is an explicit zero log-vector and leaves the density
  posterior unchanged.
- Fused outputs always stay in `[0, observed_compatible_count]`.
- Large-count approximation agrees with exact near the threshold.
- Boundary-mode large-count regions do not produce NaNs when the mode is near `0` or
  `N`.

### PR 3: Prior Allocation

Goal: create `calibration/prior.py` and produce native EM inputs.

Tasks:

1. Implement `PriorTable`.
2. Implement conservation-aware geometry allocation as the default count allocator.
3. Implement diagnostic midpoint lookup with `np.searchsorted` if inexpensive, but do
   not make midpoint allocation a correctness dependency for the MVP.
4. Implement denominator calculation through existing `_exposure.py` helpers.
5. Implement explicit `enable_gdna` from EM candidate data.
6. Add compact prior summaries to `CalibrationResult.to_summary_dict()`.

Tests:

- Partial footprint allocation conserves count mass.
- Regions touching multiple MultiLoci are counted and their fused mass is reported.
- A localized hotspot diagnostic fixture reports geometry-only smear risk; exact
  hotspot-safe allocation is not required for the MVP unless diagnostics show high mass.
- Midpoint-to-region lookup uses binary search/vectorized scanning, not nested loops.
- Uniform exposure yields `gdna_eff_len == gdna_eff_len_unweighted`.
- Zero prior count does not disable gDNA eligibility.

### PR 4: End-To-End Pipeline Wiring

Goal: remove the `NotImplementedError` and make `rigel quant` runnable.

Tasks:

1. Restore `quant_from_buffer` using the sequence in Section 7.
2. Use `partition_and_free`, not a nonexistent `build_partitions` helper.
3. Pass raw denominator and exposure diagnostics into `_run_locus_em_partitioned`.
4. Rewrite or unignore legacy tests in `tests/conftest.py` that block the new path.
5. Remove the NotImplementedError skip hook once Phase 6 tests pass.

Tests:

```bash
conda activate rigel && pytest tests/test_pipeline_wiring.py -v
conda activate rigel && pytest tests/test_pipeline_smoke.py -v
conda activate rigel && pytest tests/test_golden_output.py -v
```

Update goldens only after focused tests explain the expected output changes.

### PR 5: Benchmark Gate

Goal: prove the runnable pipeline is scientifically usable before declaring victory.

Commands:

```bash
conda activate rigel
python -m scripts.benchmarking status -c scripts/benchmarking/configs/default.yaml
python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml --dry-run
```

Then run selected conditions first:

```bash
python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml \
  --conditions gdna_none_ss_1.00_nrna_none gdna_high_ss_0.90_nrna_none

python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/default.yaml \
  -o results/benchmark_report_v2_smoke \
  --conditions gdna_none_ss_1.00_nrna_none gdna_high_ss_0.90_nrna_none
```

Report:

- total gDNA prior mass vs EM gDNA mass;
- fused density/strand weights;
- prior allocation fallback counts;
- gDNA false RNA and RNA false gDNA confusion;
- nRNA siphon changes;
- loci with largest density tail tension.

## 9. Test Cleanup Requirements

The following test posture is not acceptable for the final Phase 6 state:

- `tests/conftest.py` collection-ignores prior and pipeline tests.
- `tests/conftest.py` converts the Phase 6 `NotImplementedError` into skips.
- Legacy tests import deleted modules such as old `density_global` or `locus_prior` APIs.

Required cleanup:

1. Rewrite `tests/test_bayesian_prior_acceptance.py` around the new `PriorTable`.
2. Rewrite `tests/test_pipeline_wiring.py` around the current `BamScanner.set_regions`
   signature and the restored `quant_from_buffer` path.
3. Decide whether `tests/test_calibration_result.py` should be revived or replaced by
   focused summary tests for `CalibrationResult`.
4. Remove the Phase 6 skip hook after `quant_from_buffer` no longer raises.
5. Keep golden updates intentional and isolated.

## 10. Non-Goals For v2

These are important but should not block end-to-end runnable quantification:

- learned hybrid-capture exposure beyond `RegionExposure.from_density`;
- exact weighted gDNA effective-length denominators;
- expression-aware strand deconvolution for `TS_AMBIG`;
- BED-informed capture priors;
- full fragment span `(start, end)` storage for prior allocation;
- cgranges-based fragment-span allocation;
- C++ changes to the native EM solver.

If C++ changes become necessary, they must follow the native rebuild rule:

```bash
conda activate rigel && pip install --no-build-isolation -e .
```

## 11. Readiness Checklist

Implementation is ready to start when these are true:

- [ ] `density_model.py` has a public density log-PMF/predictive helper design.
- [ ] `strand_deconv.py` has a public strand likelihood helper design.
- [ ] `FusedRegionGdnaEvidence` fields are accepted.
- [ ] `PriorTable` fields match `_run_locus_em_partitioned` inputs.
- [ ] `calibrate()` is the chosen owner of regional fusion.
- [ ] `quant_from_buffer` wiring sequence is agreed.
- [ ] Tests for unstranded, strong stranded, `TS_AMBIG`, boundary-mode fallback,
      allocation diagnostics, and zero-prior eligibility exist.

## 12. Summary

v2 keeps the ambition of v1 but removes the fragile parts:

- no imaginary `calibration_context`;
- no duplicated density or strand math;
- no fake strand information for ambiguous regions;
- no uninstrumented base-pair smearing of localized count mass;
- no pipeline statistics hidden inside calibration code;
- no success condition based on skipped tests.

The result is a cleaner route to a runnable Rigel: calibration builds a fused regional
posterior, prior allocation projects it into the exact native EM contract, and pipeline
wiring becomes a short, obvious sequence of already-owned operations.
