# gDNA Regional Exposure Plan v1 Review

**Date**: 2026-05-18
**Reviewed plan**: `docs/calibration/gdna_exposure_plan_v1.md`

## Bottom line

The v1 plan has the right central idea: learn a regional gDNA exposure field from conservative calibration evidence, then use that field in the gDNA likelihood numerator and denominator. That is the correct modeling object for hybrid-capture gDNA.

It is not quite implementation-ready. A few details would either break the implementation, add unnecessary parameters, or make the math less internally consistent than it needs to be. The recommended v1.1 is smaller: build one exposure table, ship numerator and denominator together, vectorize lookups, and remove most hard thresholds by using exposure-weighted quantiles plus EB shrinkage that is already in the codebase.

## Must-fix issues before implementation

### 1. Stage 1 contradicts the mass-conservation decision

Section 0 says Stage 1 includes both the weighted gDNA denominator and `log A(midpoint)`. Section 4 says Stage 1 is production denominator-only and Stage 2 adds the numerator.

Do not ship an enabled denominator-only model. It changes `-log L_g` without the matching `+log A`, so it can create posterior drift that is not the target model. Either:

- implement denominator plus numerator in the same production stage, or
- land denominator plumbing behind a disabled flag and do not enable `regional_exposure` until the numerator is present.

The cleaner option is one production stage: exposure table, weighted `L_g`, and per-unit `log A` together.

### 2. A midpoint alone is insufficient for per-unit lookup

`RegionIndexPy.overlap()` needs `(ref_id, start, end)`. The proposed native ABI adds only `genomic_midpoint`, so `_apply_unit_gdna_weights()` cannot unambiguously look up `A(pos)`.

Minimal fix without adding a scanner ref column:

- emit `genomic_midpoint` from the scorer,
- infer `unit_ref_id` in Python from `em_data.locus_t_indices` and `index.t_df["ref"]`,
- document that multimappers use the best RNA candidate's reference as the representative gDNA reference in v1.

More exact but larger fix:

- propagate `genomic_ref_id` through `ResolvedFragment`, `FragmentAccumulator`, `_FinalizedChunk`, `StreamingScorer`, and `ScoredFragments`.

For v1, the inferred-ref approach is smaller and is acceptable for the no-multimap failure mode that motivated this work.

### 3. The weighted `L_g` algorithm is underspecified in a way that can be wrong

`gdna_eff_len_for_loci()` counts start positions whose length-`ell` interval overlaps the locus. For a locus `[a,b)`, the valid start window is `[a - ell + 1, b)`. If we weight by midpoint, the relevant region lookup is over midpoint positions `p + ell // 2`, not only over regions overlapping `[a,b)`.

So the weighted implementation should mirror the existing function exactly:

1. For each `ell`, build valid start windows for every `Locus`.
2. Clip to contig-valid starts.
3. Merge overlapping start windows per reference.
4. Convert each merged start window `[lo, hi)` to midpoint window `[lo + ell // 2, hi + ell // 2)`.
5. Intersect that midpoint window with `RegionArrays` and sum `overlap_bp * A_r`.
6. Weight by `pmf[ell]`.

When exposure mode is uniform, branch directly to `gdna_eff_len_for_loci()` so the no-op path is bit-exact instead of merely close.

### 4. The split-invariance proof is false under the proposed shrinkage

The plan states that splitting a region into two proportional subregions leaves the EB-shrunk density unchanged. With fixed `kappa`, this is not true:

```text
parent = (Y + kappa*rho) / (E + kappa)
half   = (Y/2 + kappa*rho) / (E/2 + kappa)
```

The half-region is pulled more strongly toward the prior unless `Y/E == rho`.

Do not make rebuilt-exposure split invariance a required test. The useful invariant is simpler: if two adjacent subregions have the same assigned `A`, weighted effective length must be unchanged after splitting. Test the geometry with fixed weights; test EB shrinkage separately.

### 5. `rho_global * A_r` is not an expected-count formula

If `A_r = rho_hat_r / rho_ref`, then `rho_global * A_r` generally equals neither `rho_hat_r` nor a correctly normalized local intensity. It double-uses the global average and can undercount heavily when `rho_global << rho_ref`.

For v1, prefer one of two cleaner choices:

- Minimal and safest: leave `gdna_prior_count` unchanged in production and add weighted-prior diagnostics only.
- If changing the prior now: compute expected gDNA count directly from the learned regional densities, `sum_r rho_hat_r * exposure_r`, with boundary terms using their boundary-channel shrunk densities.

Avoid calling `rho_global * A_r * exposure` an expected gDNA pseudocount unless the normalization is rederived.

### 6. Stage 0 cannot reconstruct per-region counts from `summary.json`

The saved summary has global density summaries, not `payload.per_region_counts`, `u_left`, `u_right`, or orientation-resolved counts. The diagnostic must rerun `scan_and_buffer` / `calibrate`, or a future production path must persist the calibration payload. The plan should remove the summary-only reconstruction path.

### 7. Per-unit Python lookups should be vectorized

The plan estimates `N_units * log(R)` scalar Python lookups as cheap. That is optimistic for tens of millions of units.

Add a vectorized method instead:

```python
weights_for_positions(ref_ids: np.ndarray, positions: np.ndarray) -> np.ndarray
```

Implementation can group by `ref_id`, use `np.searchsorted` against per-ref region starts/ends, and fill a contiguous float array. This is both faster and easier to test than a Python loop over units.

## Parameter reductions

The plan introduces too many knobs for a first patch. Most can be removed by leaning on existing EB shrinkage and exposure-weighted summaries.

### Remove `min_boundary_events_for_exon_refine`

Always compute the exon boundary-region density when boundary exposure is positive:

```text
rho_hat_exon = (Y_boundary + kappa_exon * rho_global_exon) / (E_boundary + kappa_exon)
```

If evidence is weak, EB shrinkage already pulls the estimate to the global baseline. A hard event threshold is redundant.

### Replace `min_exposure_bp` and `min_class_regions_for_quantile`

Use exposure-weighted quantiles of `log rho_hat` within each channel. Tiny regions then contribute tiny weight instead of requiring a threshold.

Fallback only when total exposure for a channel is zero. That is a structural condition, not a tunable parameter.

### Avoid `capture_index_threshold` as a model branch

The cleanest regression switch is explicit: `--regional-exposure off`.

For `auto`, prefer a continuous learned shrinkage of log-weights toward zero rather than a hard uniform/regional branch. A simple parameter-free option is:

```text
log_A_raw = min(log(rho_hat_r) - log(rho_ref), 0)
signal_k  = max(0, observed_spread_k - null_spread_k) / max(observed_spread_k, eps)
log_A_r   = signal_k * log_A_raw
```

where `observed_spread_k` is an exposure-weighted high-minus-median log-density spread and `null_spread_k` is derived from the EB posterior uncertainty already represented by `E_r + kappa_k`. Uniform libraries get `signal ~= 0`; capture libraries get `signal ~= 1`; no threshold is needed.

If that feels too much for v1, skip auto-uniform detection entirely and rely on the explicit off switch plus EB shrinkage.

### Remove `c_max` from the public surface

`lambda_floor = 1/c_max` is a cap parameter. It is not clearly biological. The EB formula already prevents arbitrary zero weights except in genuinely high-exposure zero-count regions.

For v1, use only a numerical log floor for safety, not a user-facing cap. If diagnostics show pathological overconfidence, learn a floor from the empirical lower tail of the exposure distribution and report it in `summary.json`.

### Keep CLI surface to one flag

Recommended v1 CLI:

```text
--regional-exposure {auto,off}
```

Do not expose `c_max`, threshold, minimum exposure, minimum class count, or boundary-event thresholds in the first production patch.

## Recommended v1.1 shape

### Config

Keep the config small:

```python
@dataclass(frozen=True)
class RegionalExposureConfig:
    enabled: bool = True
```

If diagnostics need fixed estimator constants such as a 0.95 reference quantile, keep them module-level constants in `_regional_exposure.py`, not user-facing configuration.

### Exposure table

`RegionalGdnaExposure` should own only the arrays and summaries needed downstream:

- `rho_hat: float64[R]`
- `log_weight: float64[R]`
- `weight: float64[R]` if needed for speed
- `mode: "uniform" | "regional"`
- `rho_ref` and per-class diagnostics
- `weights_for_positions(ref_ids, positions)` for vectorized unit weighting

The table should be aligned to `RegionArrays` sorted order and should validate array lengths at construction.

### Production pipeline order

Recommended order inside `quant_from_buffer`:

1. Score fragments.
2. Build multi-loci.
3. Assemble priors with `regional_exposure` available.
4. Apply per-unit `log A` while global `em_data` is still alive.
5. Partition and free.
6. Run EM.

This respects the existing gotcha that `assemble_priors()` and any global-unit mutation must happen before `partition_and_free()` nulls the arrays.

### First implementation slice

The smallest useful production slice is:

1. Build `RegionalGdnaExposure` in calibration.
2. Add weighted `L_g` with an exact uniform branch.
3. Emit `genomic_midpoint` from `StreamingScorer`.
4. Apply vectorized `log A` before partitioning.
5. Leave `gdna_prior_count` unchanged, but emit diagnostics comparing current and candidate weighted priors.

That slice directly tests the hybrid-capture denominator hypothesis with fewer moving parts. If it improves VCaP without uniform-library drift, the prior can be revisited with direct `rho_hat * exposure` rather than `rho_global * A`.

## Tests to adjust

Keep these tests:

- Uniform exposure returns the existing `gdna_eff_len_for_loci()` bit-exactly.
- Fixed-weight partition invariance: splitting adjacent equal-`A` regions does not change weighted `L_g`.
- Two-state geometry sanity: `A=1` over 10 kb and `A=0.01` over 90 kb gives the expected weighted opportunity away from FL boundary effects.
- Per-unit log-weight applies `log A` and leaves `-inf` gDNA likelihoods untouched.
- No-multimap VCaP hotspot diagnostic shows the predicted posterior movement before enabling by default.

Change these tests:

- Do not assert rebuilt EB exposure is invariant under arbitrary region splitting.
- Do not require a hard capture-index threshold to choose uniform mode.
- Do not require denominator-only Stage 1 benchmark gates.

## Final recommendation

Revise v1 into a smaller v1.1 before coding. The core model should be:

```text
rho_hat_r = EB_shrunk_conservative_gDNA_density(region r)
A_r       = learned relative exposure from rho_hat_r
log p_g(f) += log A(midpoint_f)
L_g       = sum over valid gDNA starts of A(midpoint)
```

Everything else should be diagnostic until the likelihood correction is validated on the five real-data hotspots and the uniform synthetic suite.