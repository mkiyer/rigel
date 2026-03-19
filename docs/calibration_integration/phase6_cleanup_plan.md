# Phase 6: Cleanup and Deprecation — Implementation Plan

**Date**: 2026-03-17  
**Status**: ✅ COMPLETE (995 tests passing)  
**Prerequisite**: Phases 1–4 complete (1085 tests passing)  
**Principle**: Calibration is mandatory. One code path. No fallbacks. Clear module contracts.

---

## 1. Design Philosophy

**Before Phase 6**: Two parallel code paths — a "calibrated" path and an "uncalibrated"
fallback — maintained in superposition. Every function that touches gDNA priors carries
`if calibration is not None` branches, `Optional` type annotations, and 8+ MoM kappa
parameters "just in case."

**After Phase 6**: Calibration always runs. `calibrate_gdna()` always returns a fully
populated `GDNACalibration` object with sensible internal fallbacks. Downstream code
receives the calibration result unconditionally — no None checks, no branching.

---

## 2. Module Contracts (Post-Cleanup)

### 2.1 `calibration.py` — The Authoritative Source for gDNA Parameters

**Contract**: `calibrate_gdna()` ALWAYS returns a valid `GDNACalibration`.
It handles its own edge cases internally:

| Scenario | Internal Behavior | Output |
|----------|-------------------|--------|
| Normal (hundreds+ regions) | Full EM deconvolution | Complete calibration |
| Too few regions (`< min_gdna_regions`) | Algebraic fallback (hard constraints) | Valid `GDNACalibration` with `gdna_fl_model=None`, hard-constraint posteriors |
| No eligible regions | Zero-density bailout | Valid `GDNACalibration` with `gdna_fl_model=None`, zeros |
| ESS too low for FL model | EM runs, FL gated | Valid `GDNACalibration` with `gdna_fl_model=None` |

The `gdna_fl_model` field is the ONLY nullable field. Everything else is always populated.
The bailout paths need enhancement (see §4.2) to provide useful `kappa`, `gdna_density_per_ref`,
and `gdna_density_global` even when the EM doesn't run — using algebraic strand-based
estimation as the fallback *inside calibration*, not outside it.

### 2.2 `locus.py` — Per-Locus gDNA Initialization

**Contract**: `compute_eb_gdna_priors()` receives a `GDNACalibration` (required, not Optional).
It uses `calibration.gdna_density_global`, `.gdna_density_per_ref`, `.kappa` directly.
It performs the per-locus hybrid density + strand estimation and hierarchical shrinkage.

Removed from this module:
- Global/ref-level strand deconvolution (moved into calibration bailout)
- MoM κ estimation (replaced by calibration's κ)
- `_compute_ref_gdna_densities()` (replaced by `calibration.gdna_density_per_ref`)
- `compute_gdna_density_from_strand()` (only caller was the uncalibrated path + priors.py)

Retained:
- `compute_gdna_density_hybrid()` — still used per-locus (strand + density anchor)
- `_compute_per_locus_gdna_densities()` — still used (simplified: `calibrated_ref_densities` always provided)
- Final locus→ref shrinkage loop

### 2.3 `pipeline.py` — Orchestration

**Contract**: `run_pipeline()` always calls `calibrate_gdna()` when region data is available.
`quant_from_buffer()` always receives a `GDNACalibration` (required).

- `CalibrationConfig.enabled` field: **deleted**
- `calibration: GDNACalibration | None` parameters: changed to `calibration: GDNACalibration`
- `_compute_priors()`: simplified — no MoM kappa threading

### 2.4 `config.py` — Simplified Configuration

**EMConfig**: 7 dead MoM kappa fields removed:
- `gdna_kappa_ref` → **DELETED** (calibration provides κ)
- `gdna_mom_min_evidence_ref` → DELETED
- `gdna_mom_min_evidence_locus` → DELETED
- `gdna_kappa_min` → DELETED
- `gdna_kappa_max` → DELETED
- `gdna_kappa_fallback` → DELETED
- `gdna_kappa_min_obs` → DELETED
- `gdna_kappa_locus` → **KEEP** as override (calibrated path already respects it)

**CalibrationConfig**: `enabled` field **DELETED**.

### 2.5 `priors.py` — **DELETE ENTIRE MODULE**

Both functions are dead:
- `compute_global_gdna_density()` — no production caller
- `estimate_kappa()` — only called by uncalibrated path (being deleted)

The `estimator.py` re-export (`from .priors import ...`) is also removed.

---

## 3. Change Inventory

### 3.1 Files to DELETE

| File | Lines | Reason |
|------|-------|--------|
| `src/rigel/priors.py` | ~110 | Both functions dead. Module fully superseded by calibration. |
| `src/rigel/region_evidence.py` | ~260 | Superseded by C++ `RegionAccumulator` (Phase 1). |

### 3.2 Files to EDIT (Production Code)

| File | Change | Est. Lines |
|------|--------|------------|
| `src/rigel/config.py` | Remove 7 MoM kappa fields from `EMConfig`; remove `CalibrationConfig.enabled` | −40 |
| `src/rigel/calibration.py` | Enhance bailout paths with algebraic density/κ fallbacks; strengthen `GDNACalibration` contract | +30, −5 |
| `src/rigel/locus.py` | Delete uncalibrated `else` block (~73 lines), `_compute_ref_gdna_densities()` (~35 lines), `compute_gdna_density_from_strand()` (~19 lines); make `calibration` required; remove MoM kappa params from `compute_eb_gdna_priors()` signature | −140, +5 |
| `src/rigel/pipeline.py` | Remove `cal_cfg.enabled` gate; make `calibration` required in `quant_from_buffer()`/`_compute_priors()`; remove MoM kappa threading; simplify CAL-FL block | −30, +5 |
| `src/rigel/cli.py` | Remove 7 `--gdna-kappa-*` / `--gdna-mom-*` CLI flags and their `_ParamSpec` entries | −50 |
| `src/rigel/estimator.py` | Remove `from .priors import ...` re-export line | −1 |

### 3.3 Files to EDIT (Tests)

| File | Change |
|------|--------|
| `tests/test_gdna.py` | Delete `TestComputeGdnaDensityFromStrand` class (~60 lines); delete any tests exercising `_compute_ref_gdna_densities`; update tests that set MoM kappa fields |
| `tests/test_calibration_integration.py` | Remove tests for `CalibrationConfig(enabled=False)` path; remove `enabled=False` from comparisons; simplify A/B tests |
| `tests/test_calibrated_density.py` | Remove uncalibrated comparison tests that exercise the deleted code path |
| `tests/test_gdna_fl_scoring.py` | Remove `calibration_enabled=False` test variants (or convert to `min_gdna_regions=1_000_000` to simulate no-calibration effect) |
| `tests/test_region_evidence.py` | Delete Section 4 (`TestCountRegionEvidence`) that tests the Python pysam path; keep Sections 1–3 for counting logic if useful, otherwise delete |
| `tests/test_estimator.py` | Remove references to `estimate_kappa`, `compute_global_gdna_density` |
| `docs/parameters.md` | Remove 7 MoM kappa parameter entries; update CalibrationConfig section (drop `enabled`) |

### 3.4 Files UNCHANGED

| File | Reason |
|------|--------|
| `src/rigel/frag_length_model.py` | `mix_models()` still produces baseline gdna/rna models. Unaffected. |
| `src/rigel/scoring.py` | Reads from already-replaced `frag_length_models.gdna_model`. No calibration awareness. |
| `src/rigel/native/*` (C++) | No changes. C++ has no calibration awareness. |

---

## 4. Detailed Changes

### 4.1 `config.py` — Simplify EMConfig and CalibrationConfig

**EMConfig** — Delete these 7 fields and their docstrings:
```python
# DELETE:
gdna_kappa_ref: float | None = None
gdna_mom_min_evidence_ref: float = 50.0
gdna_mom_min_evidence_locus: float = 30.0
gdna_kappa_min: float = 2.0
gdna_kappa_max: float = 200.0
gdna_kappa_fallback: float = 5.0
gdna_kappa_min_obs: int = 20
```

**KEEP** in EMConfig:
```python
gdna_prior_scale: float = 1.0       # Used in EM M-step
gdna_kappa_locus: float | None = None  # Manual override for per-locus shrinkage
```

**CalibrationConfig** — Delete `enabled: bool = True`. The remaining fields
(`max_iterations`, `convergence_tol`, `density_percentile`, `min_gdna_regions`, `min_fl_ess`)
stay as-is.

### 4.2 `calibration.py` — Self-Contained Fallbacks

The bailout paths (no eligible regions; too few regions) currently return degenerate values
(`gdna_density_global=0.0`, `gdna_density_per_ref={}`, `kappa=0.0`). These work but are
poor defaults for the per-locus shrinkage that runs downstream.

**Enhancement**: when bailing out, compute algebraic fallback densities and κ from the
available region stats. The algebraic fallback replaces what the uncalibrated path in
`locus.py` used to do, but scoped inside calibration:

```python
# In the "too few regions" bailout:
# 1. Fallback density: use unspliced counts / region length
fallback_density = _compute_fallback_density(stats, eligible)
fallback_per_ref = _compute_fallback_per_ref_density(stats, region_df, eligible)
# 2. Fallback kappa: uninformative default
fallback_kappa = 2.0  # Uniform Beta(1,1) — minimal shrinkage
```

This ensures `GDNACalibration` always provides non-degenerate values even on tiny datasets.
The pipeline never needs to know *how* calibration produced its numbers.

**Contract strengthening**:
- `gdna_density_global`: always ≥ 0 (can be 0.0 legitimately)
- `gdna_density_per_ref`: always a dict (can be empty `{}` on no-data bailout, but
  the too-few-regions bailout should provide per-ref values)
- `kappa`: always > 0 (default 2.0 on bailout)
- `gdna_fl_model`: `None` when ESS is insufficient (the ONLY nullable output)

### 4.3 `locus.py` — Delete Uncalibrated Path

Delete from `compute_eb_gdna_priors()`:
1. The entire `else:` block (lines 870–943): global strand decon, MoM κ_ref, ref shrinkage, MoM κ_locus
2. `_compute_ref_gdna_densities()` function (lines 656–690)
3. `compute_gdna_density_from_strand()` function (lines 565–588)
4. All MoM parameters from the function signature (`mom_min_evidence_ref`, `mom_min_evidence_locus`, `kappa_min`, `kappa_max`, `kappa_fallback`, `kappa_min_obs`)
5. The `kappa_params` dict construction
6. The `if calibration is not None` / `else` structure → just use calibration directly

**Simplified `compute_eb_gdna_priors()` signature**:
```python
def compute_eb_gdna_priors(
    loci: list[Locus],
    em_data: ScoredFragments,
    estimator: AbundanceEstimator,
    index: TranscriptIndex,
    strand_models: StrandModels,
    *,
    calibration: GDNACalibration,        # REQUIRED
    intergenic_density: float = 0.0,     # fallback density anchor
    kappa_locus: float | None = None,    # manual override
) -> list[float]:
```

Body becomes the former "calibrated path" only (~6 lines for setup) + the shared
per-locus section (unchanged).

### 4.4 `pipeline.py` — Make Calibration Mandatory

**`run_pipeline()`**: Remove the `cal_cfg.enabled` guard. The logic becomes:

```python
calibration = None
if region_counts is not None and fl_table is not None and index.region_df is not None:
    calibration = calibrate_gdna(...)
else:
    logger.warning("[CAL] No region data — cannot run calibration")
```

Decision point: should missing region data be a **warning** or an **error**?
- For maximum strictness: raise `RuntimeError` — require calibration-equipped indices.
- For pragmatism: warn and compute a minimal fallback `GDNACalibration` with  
  `gdna_density_global=intergenic_density, per_ref={}, kappa=2.0, gdna_fl_model=None`.

**Recommendation**: Error. If we're making calibration mandatory, require a region-equipped
index. Old indices without `regions.feather` should be rebuilt. This avoids a stealth
degradation path.

**`quant_from_buffer()`**: Change `calibration: GDNACalibration | None = None` to
`calibration: GDNACalibration`. Remove the `calibration is not None` guard on the
CAL-FL block — just check `calibration.gdna_fl_model is not None`.

**`_compute_priors()`**: Remove the 8 MoM kappa parameters. 
Simplified call to `compute_eb_gdna_priors()`:

```python
gdna_inits = compute_eb_gdna_priors(
    loci, em_data, estimator, index, strand_models,
    calibration=calibration,
    intergenic_density=intergenic_density,
    kappa_locus=em_config.gdna_kappa_locus,
)
```

### 4.5 `cli.py` — Remove Dead CLI Flags

Delete 7 `argparse` flag definitions and their `_ParamSpec` mappings:
- `--gdna-kappa-ref`
- `--gdna-mom-min-evidence-ref`
- `--gdna-mom-min-evidence-locus`
- `--gdna-kappa-min`
- `--gdna-kappa-max`
- `--gdna-kappa-fallback`
- `--gdna-kappa-min-obs`

**KEEP**: `--gdna-kappa-locus` (still functional as manual override).

### 4.6 `priors.py` — Delete Entire Module

Both functions are dead:
- `compute_global_gdna_density()` → no production caller
- `estimate_kappa()` → only called by the uncalibrated path being deleted

Also remove the re-export in `estimator.py`:
```python
# DELETE from estimator.py:
from .priors import compute_global_gdna_density, estimate_kappa  # noqa: F401
```

### 4.7 `region_evidence.py` — Delete Entire Module

The Python/pysam BAM traversal was superseded by C++ `RegionAccumulator` in Phase 1.
No production code imports it. Tests in `test_region_evidence_cpp.py` use the C++ path.

The internal counting helpers (`assemble_fragment`, `count_fragment`, `_tabulate_fragment_lengths`)
are tested in `test_region_evidence.py` Sections 1–3, but these are testing a dead code path.
Delete the module and the test file together.

---

## 5. Execution Order

The changes have dependencies — they must be executed in order:

| Step | What | Why |
|------|------|-----|
| **6a** | Enhance calibration bailouts with fallback density/κ | Ensures `GDNACalibration` always provides non-degenerate values before we remove downstream fallbacks |
| **6b** | Delete `CalibrationConfig.enabled`; make calibration mandatory in `pipeline.py` | Un-gates calibration; `calibration` is now always a `GDNACalibration` |
| **6c** | Simplify `locus.py`: delete uncalibrated path, `_compute_ref_gdna_densities`, `compute_gdna_density_from_strand`, MoM params | Remove the biggest dead code block |
| **6d** | Simplify `pipeline.py`: `_compute_priors()` drop MoM params, `quant_from_buffer()` require calibration | Clean up the orchestration |
| **6e** | Delete 7 MoM kappa fields from `EMConfig`, 7 CLI flags from `cli.py` | Config/CLI cleanup |
| **6f** | Delete `priors.py` module + `estimator.py` re-export | Module elimination |
| **6g** | Delete `region_evidence.py` module | Module elimination |
| **6h** | Update tests: delete dead-path tests, update signatures | Test cleanup |
| **6i** | Regenerate golden outputs, run full suite | Validation |
| **6j** | Update docs: `parameters.md`, `cleanup_tracking.md`, `implementation_plan.md` | Documentation |

---

## 6. Risk Assessment

| Risk | Mitigation |
|------|-----------|
| **Bailout produces degenerate calibration → bad gDNA priors on tiny datasets** | Step 6a: algebraic fallback density/κ inside calibration. Covered by `test_all_spliced_no_gdna` and new bailout tests. |
| **Old indices without `regions.feather` produce hard error** | Intentional. Users rebuild index with `rigel index`. Document in CHANGELOG. |
| **Removing CLI flags breaks existing YAML configs** | Accepted — these flags were never documented in the user manual and are only used for debugging. Document in CHANGELOG. |
| **Tests that compare calibrated vs uncalibrated paths break** | Step 6h: convert to single-path tests that validate calibration quality in isolation. |

---

## 7. Lines of Code Impact (Estimated)

| Category | Added | Deleted | Net |
|----------|-------|---------|-----|
| Production code | +40 | −350 | **−310** |
| Test code | +20 | −200 | **−180** |
| Docs | +20 | −80 | **−60** |
| **Total** | +80 | −630 | **−550** |

The codebase becomes ~550 lines lighter. One code path. Zero branching on calibration state.
