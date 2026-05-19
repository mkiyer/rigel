# gDNA Regional Exposure Plan v2 Review

**Date**: 2026-05-18
**Reviewed plan**: `docs/calibration/gdna_exposure_plan_v2.md`

## Verdict

v2 is close, and the core likelihood model is ready to prototype. It is not yet fully implementation-ready as written. The remaining issues are smaller than v1's, but a few would cause wrong behavior or unnecessary implementation churn if coding starts directly from the document.

The strongest parts of v2 are worth keeping:

- Denominator and numerator ship together.
- `gdna_prior_count` stays unchanged in v2.
- Weighted `L_g` is defined on midpoints of expanded start windows.
- The public CLI surface is a single on/off switch.
- Lookups are vectorized.

Before implementation, fix the items below.

## Must Fix Before Coding

### 1. Per-unit weighting is in the wrong pipeline order

The plan applies `_apply_unit_gdna_weights()` before `assemble_priors()`. That changes `em_data.gdna_log_liks` before `enable_gdna_for_multilocus()` runs inside `assemble_priors()`.

`enable_gdna_for_multilocus()` currently checks only finiteness, so the result probably stays the same, but this is a brittle hidden dependency. More importantly, the weight function needs `multi_loci` context for robust reference inference and validation.

Recommended order inside `quant_from_buffer`:

1. Score fragments.
2. Build `multi_loci`.
3. Assign locus IDs.
4. Assemble priors with the regional exposure, using unmodified `gdna_log_liks` for eligibility.
5. Apply per-unit `log A` while global `em_data` still exists.
6. Partition and free.
7. Run EM.

This preserves the existing `assemble_priors()` before `partition_and_free()` invariant and keeps eligibility independent of score perturbations.

### 2. `index.t_to_ref_arr` does not exist

The plan uses `index.t_to_ref_arr[first_t]`, but `TranscriptIndex.load()` currently materializes `t_to_g_arr`, `t_to_strand_arr`, and `g_to_strand_arr`, not a transcript-to-reference array.

Either add `t_to_ref_arr` during index load:

```python
t_ref_codes = index.t_df["ref"].cat.codes.values
cat_to_ref_id = np.array(
    [index.ref_name_to_id[str(name)] for name in index.t_df["ref"].cat.categories],
    dtype=np.int32,
)
index.t_to_ref_arr = cat_to_ref_id[t_ref_codes]
```

or compute this locally in `_apply_unit_gdna_weights()`. Adding the cached array is cleaner because `locus_prior.py` already reconstructs this mapping internally.

### 3. Inferred reference is not always a safe proxy for gDNA hit reference

For the no-multimap VCaP failure mode, inferring `ref_u` from the unit's best transcript is probably fine. For multimappers and multi-reference `MultiLocus` components, it is only an approximation.

The plan should make the v2 behavior explicit:

- v2 production lookup uses `locus_t_indices[u]` to infer a representative transcript reference.
- If a unit's candidate transcripts span multiple references, skip the per-unit numerator term for that unit or fall back to `A=1`.
- Count skipped cross-reference units in `summary.json`.

This is preferable to silently applying the wrong chromosome's exposure. Per-hit weighted gDNA LSE can remain deferred.

### 4. `genomic_midpoint` must be scattered only if partitions need it

The plan says `locus_partition.partition_and_free` scatters `genomic_midpoint`. But `_run_locus_em_partitioned()` does not pass midpoint to native EM, and after per-unit weighting it is no longer needed.

Do not scatter it in v2 unless a downstream diagnostic truly consumes it after partitioning. Keep the field on global `ScoredFragments`, apply weights before `partition_and_free()`, then let the global array be dropped. This avoids expanding `LocusPartition` and native scatter plumbing for no runtime benefit.

### 5. Regional diagnostics need a concrete data path

`loci.feather` is built from `estimator.locus_results`, not directly from `PriorTable`. The plan says to emit `gdna_eff_len_unweighted`, `gdna_eff_len_weight_ratio`, and `gdna_prior_count_regional`, but it does not specify how those values enter `estimator.locus_results`.

Recommended implementation:

- Extend `PriorTable` with optional arrays:
  - `gdna_eff_len_unweighted`
  - `gdna_prior_count_regional`
- In `_run_locus_em_partitioned()`, pass these arrays alongside `gdna_prior_count` and `gdna_eff_len` and store the per-locus values in `_build_locus_meta()`.
- Extend `AbundanceEstimator.get_loci_df()` columns to read those keys.

Do not try to modify `loci.feather` in `cli.py`; the estimator owns that output table.

### 6. `CalibrationResult` needs `regional_exposure` to survive `with_priors()`

Adding a field is straightforward, but the plan should explicitly require `dataclasses.replace()` in `with_priors()` to preserve it, and `build_calibration_result()` to require or construct a uniform exposure when disabled. This keeps `calibration.regional_exposure` always present, as v2 assumes.

### 7. Strand-corrected per-region counts are underspecified

The plan says intergenic and intron `Y_r` are strand-corrected via the existing `compute_global_densities` correction when active. The existing correction is global/channel-level, not a reusable per-region function. For introns, orientation arrays exist (`intron_counts_by_orient`), so a per-region corrected moment can be implemented. For intergenic, rows are strand-uninformative, so there is no strand correction to apply.

Make this explicit:

- Intergenic: use raw `PayloadArrays.intergenic_per_region`.
- Intron: use `_gdna_count_moment(same, opp, signed_strand_contrast=...)` on strand-identifiable rows; otherwise raw total.
- Exon boundary: use orientation-resolved boundary flux if strand active; otherwise raw eligible `u_left/u_right`.

This probably requires expanding `PayloadArrays` or passing the original payload into `_regional_exposure.build()`.

### 8. Auto-uniform attenuation is elegant but too much for first implementation

The null-spread closed form plus Poisson fallback is the least settled part of v2. It is parameter-light, but it is also a new statistical submodel with nontrivial test burden and edge cases around tiny `rho`, zero counts, and log transforms.

Recommendation: keep the formulas in the plan, but implement v2a with a simpler explicit mode:

- `--regional-exposure off`: uniform.
- `--regional-exposure auto`: build regional weights directly from EB densities and record observed/null spread diagnostics, but do not attenuate by `signal_k` yet.

Then gate default enablement on Stage-0 plus uniform synthetic benchmarks. If uniform drift is too high, add signal attenuation as v2b. This is a small retreat from maximal elegance, but it lowers implementation risk substantially.

If the team wants signal attenuation in the first patch, it needs a more precise numerical spec: epsilon value, zero-density log handling, whether the Poisson fallback is deterministic via fixed seed or analytic-only for tests, and exact exposure-weighted quantile semantics.

### 9. Weighted `L_g` tests need one correction

The W5 test says "degenerate exposure (all `A_r = 0`)" but v2 defines `A_r ∈ (0, 1]` with `LOG_A_FLOOR`, so all-zero exposure cannot occur. Rewrite W5 as either:

- no region covers the midpoint window, so the lookup contributes identity weight and the result follows the unweighted path for uncovered positions, or
- all regions have `A = exp(LOG_A_FLOOR)`, so the result is tiny but positive and then clamped to `min_value` only if below the floor.

The current wording would test an impossible state.

## Readiness Recommendation

Implement after one small v2.1 edit to the plan. The v2.1 should:

1. Move per-unit weighting after `assemble_priors()` and before `partition_and_free()`.
2. Add or compute `t_to_ref_arr` explicitly.
3. Define cross-reference unit behavior as skip / `A=1` with a counter.
4. Drop midpoint scattering into `LocusPartition` unless needed for a named diagnostic.
5. Add the `PriorTable -> _run_locus_em_partitioned -> estimator.get_loci_df()` diagnostic path.
6. Specify per-region strand correction exactly.
7. Either defer signal attenuation to v2b or fully specify deterministic null-spread numerics.

With those edits, the implementation path is clear and maintainable.