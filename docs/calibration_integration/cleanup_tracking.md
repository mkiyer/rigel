# Calibration Integration — Obsolete / Stale Code Tracking

**Date**: 2026-03-18 (updated 2026-03-19)
**Context**: Phase 1 integrated C++ region accumulation into the BAM scanner.
Phase 2 wired `calibrate_gdna()` into the pipeline and refactored `compute_eb_gdna_priors()`
with calibrated/uncalibrated code paths. Phase 3 replaced the crude intergenic density
anchor with calibrated per-ref density in the per-locus hybrid estimator. Phase 4 flows
the calibrated gDNA FL model into scoring with a discrimination-based quality guard.

---

## Code That Becomes Obsolete

### 1. `src/rigel/region_evidence.py` — `count_region_evidence()`

**Status**: Functionally replaced by C++ `RegionAccumulator` + `pipeline.py` extraction.

**Keep for now**: This module serves as the **oracle reference** for comparison tests
(`tests/test_region_evidence_cpp.py`). It should be retained until:
- Phase 2 is complete and the calibration pipeline is wired end-to-end.
- Confidence is established via large-scale BAM comparison (not just synthetic mini-BAMs).

**When to remove**: After Phase 2, remove the `count_region_evidence()` function and
the pysam-based BAM traversal. The unit tests in `tests/test_region_evidence.py`
(Sections 1–3: `TestAssembleFragment`, `TestCountFragment`, `TestFragLengthTabulation`)
can be retained since they test the internal counting logic independently.
Section 4 (`TestCountRegionEvidence`) can be removed when the oracle is no longer needed.

### 2. `tests/test_region_evidence.py::TestCountRegionEvidence` (Section 4)

**Status**: Integration tests for the Python pysam path. Overlap with
`tests/test_region_evidence_cpp.py::TestCppVsPythonRegionEvidence`.

**When to remove**: After `count_region_evidence()` is removed.

### 3. Index-side `region_cr` (Python cgranges)

**Status**: `TranscriptIndex.region_cr` (Python `_cgranges_impl` object built during
`load()`) is still used by the Python `count_region_evidence()`. The C++ path builds
its own `region_cr_` inside `FragmentResolver.build_region_index()`.

**When to remove**: After the Python `count_region_evidence()` is removed, the
`region_cr` attribute and the cgranges construction in `TranscriptIndex.load()`
(lines 988–992 of `index.py`) can be removed. The `region_df` must be retained
since it's passed to `calibrate_gdna()`.

---

## Code That Changed Signatures

### 1. `pipeline.py::scan_and_buffer()`

**Before**: Returns `tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]`

**After**: Returns `tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer, pd.DataFrame | None, pd.DataFrame | None]`

**Impact**: The only call site is `run_pipeline()` (updated). No external callers known.

---

## Phase 2: Newly Superseded Code Paths

### 4. `_compute_intergenic_density()` in `pipeline.py`

**Status**: Functionally superseded by `GDNACalibration.gdna_density_global` when calibration
is enabled. Still called in `_compute_priors()` as the intergenic density fallback anchor
for the uncalibrated path and as a supplementary signal in `compute_gdna_density_hybrid()`.

**When to consider removal**: When the uncalibrated path is no longer needed (all pipelines
use calibration). Until then, keep as fallback.

### 5. Strand-based global/ref density estimation in `compute_eb_gdna_priors()`

**Status**: The entire uncalibrated path (global `compute_gdna_density_hybrid()`, MoM κ_ref
estimation, `_compute_ref_gdna_densities()`) is bypassed when `calibration` is provided.
This code is preserved as the fallback path for when no region data is available.

**When to remove**: When all production pipelines use calibration. Retain for robustness
in edge cases (e.g., index built without region partition).

### 6. `EMConfig.gdna_kappa_ref` / `gdna_kappa_locus` parameters

**Status**: Still functional — they override the calibrated κ if explicitly set. However,
the calibrated κ from `calibrate_gdna()` is now the default source for both κ_ref and
κ_locus via the calibrated path. The EMConfig params serve as manual overrides.

**Action**: No removal needed. These remain useful for debugging / experimentation.

---

## Phase 3: Per-Locus Density Anchor Supersession

### 7. `intergenic_density` as per-locus density anchor in `_compute_per_locus_gdna_densities()`

**Status**: When calibration is enabled, the per-locus hybrid estimator now uses
`calibrated_ref_densities[primary_ref]` as the density anchor instead of the crude
`intergenic_density` scalar. The fallback to `intergenic_density` is preserved for
loci whose reference is missing from the calibrated dict, and for the uncalibrated path.

**When to remove**: The `intergenic_density` parameter remains needed as a fallback.
No removal planned.

---

## Phase 4: gDNA FL Model Replacement

### 8. `frag_length_models.gdna_model` (strand-deconvolved gDNA FL)

**Status**: Conditionally superseded by `calibration.gdna_fl_model` when calibration
is enabled AND the calibrated model passes the ESS gate (Σγ·w ≥ `min_fl_ess`).
The strand-deconvolved model is retained when: (a) calibration is disabled,
(b) too few eligible regions for calibration (`< min_gdna_regions`), or
(c) the ESS of γ-weighted FL observations is below `min_fl_ess`.

**When to remove**: The strand-deconvolved model remains needed as a fallback.
The `mix_models()` function in `frag_length_model.py` is unchanged. No removal planned.

### 9. `_FragmentLengthModels.mix_models()` strand-deconvolution logic

**Status**: Still active — produces the baseline `gdna_model` and `rna_model` before
calibration. The calibrated path conditionally replaces only `gdna_model`.
The `rna_model` from `mix_models()` is always used.

**When to remove**: Not planned. The strand-deconvolved FL model serves as the
fallback when calibration's ESS is insufficient.

---

## No Code Removed Yet

Phases 1–4 are additive — all new code, no deletions. Cleanup deferred until Phase 5+ validation.
