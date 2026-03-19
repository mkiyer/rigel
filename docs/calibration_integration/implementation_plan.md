# gDNA Calibration Integration — Implementation Plan

**Date**: 2026-03-17 (updated 2026-03-19)
**Status**: Phase 4 complete — Calibrated gDNA FL model used in scoring with discrimination guard
**Scope**: Integrate `calibration.py` region-based gDNA deconvolution into the main Rigel pipeline

---

## 1. Executive Summary

Rigel's main pipeline currently estimates gDNA contamination using a strand-based Empirical Bayes (EB) shrinkage hierarchy (global → per-reference → per-locus) computed from per-transcript unspliced sense/antisense counts. This system has known weaknesses:

1. **κ (strand symmetry concentration) is hardcoded** — the default `strand_symmetry_kappa=6` is catastrophically weak (~2 pseudo-counts vs thousands of fragments).
2. **Density estimation is crude** — relies on a single intergenic density scalar or strand deconvolution, both fragile.
3. **gDNA fragment length is not modeled** — the EM uses the global FL model for gDNA candidates rather than a gDNA-specific FL distribution.
4. **No pre-EM classification step** — gDNA/RNA separation is left entirely to the locus-level EM, which is vulnerable to mega-locus θ accumulation (the root cause of gDNA siphoning).

The calibration system (`calibration.py`) solves all four problems via a region-level Aggregate-First EM that classifies genomic regions as "not expressed" (gDNA-only) vs "expressed" using three convergent signals (density, strand, fragment length). Its outputs provide:

- **Calibrated gDNA density** (global + per-reference) — replaces crude intergenic density
- **Calibrated κ** — replaces hardcoded/auto-estimated κ
- **gDNA FL model** — provides gDNA-specific fragment length distribution
- **Per-region posteriors (γ)** — can seed EB priors with region-level certainty
- **Mixing proportion (π)** — global gDNA contamination rate

**However**, the calibration system currently exists as a **separate code path** requiring a second BAM traversal via pysam (`region_evidence.py`). This plan coalesces the two paths into a single unified pipeline.

---

## 2. Current Architecture: Two Parallel Paths

### Path A: Main Pipeline (`pipeline.py`)
```
BAM → C++ BamScanner → FragmentAccumulator → FragmentBuffer
    → FragmentScorer → FragmentRouter → ScoredFragments (CSR)
    → build_loci() → compute_eb_gdna_priors() → batch_locus_em()
```

**Data collected during C++ scan** (per fragment):
- Transcript indices + CSR structure
- splice_type, exon_strand, sj_strand, chimera_type
- genomic_footprint, genomic_start, read_length, NM
- Per-transcript: frag_lengths, exon_bp, intron_bp, unambig_intron_bp
- StrandObservations (spliced, exonic, intergenic) for model training
- FragLenObservations (spliced, unspliced, intergenic, same/opp strand)

**gDNA prior computation** (`locus.py: compute_eb_gdna_priors`):
- Uses per-transcript unspliced sense/antisense from scoring pass
- Strand-based deconvolution: `G = 2(A·SS − S·(1−SS)) / (2SS−1)`
- Hybrid density+strand weighting: `W = (2SS−1)²`
- Hierarchical shrinkage: global → ref → locus via MoM-estimated κ
- Final: `gdna_init = shrunk_density × exonic_bp_locus`

### Path B: Calibration (`region_evidence.py` → `calibration.py`)
```
BAM → pysam iteration → assemble_fragment() → count_fragment()
    → region_counts + fl_table
    → calibrate_gdna() → GDNACalibration
```

**Data collected during pysam scan** (per fragment → per region):
- Fractional counts: n_unspliced_pos, n_unspliced_neg, n_spliced_pos, n_spliced_neg
- Fragment-length observations (unspliced, single-region, unique mappers)

**Calibration EM** (`calibrate_gdna`):
- Region summary stats: strand_ratio, splice_rate, density, gene_strand
- Sense fraction (gene-strand-corrected)
- Log-density with Gaussian mixture
- Beta-Binomial strand model with shared κ MLE
- gDNA FL model (γ-weighted histogram)
- Three-channel E-step: density LLR + strand LLR + FL LLR → posterior γ

---

## 3. Redundancy Analysis: What Overlaps

### 3.1 BAM Traversal (Complete Redundancy)
Both paths parse the same BAM file. Path A uses C++ htslib (fast, multi-threaded). Path B uses pysam (slow, single-threaded). This is the most expensive redundancy.

### 3.2 Fragment Assembly (Substantial Redundancy)
Both paths:
- Group records by query name
- Walk CIGAR to extract aligned blocks
- Merge overlapping intervals
- Determine splice state from intron (N) operations
- Compute strand from read flags (R1/R2 orientation)
- Compute genomic footprint

**Path A** additionally: resolves fragments against the transcript interval tree, computes per-transcript overlap (exon_bp, intron_bp), determines chimeric status.

**Path B** additionally: overlaps blocks against the region partition (cgranges), fractionally distributes counts.

### 3.3 Strand Information (Partial Redundancy)
- Path A: collects per-fragment exon_strand + sj_strand + per-tx gene strand → trains StrandModels → computes strand_specificity (SS)
- Path B: needs per-region strand_ratio + gene_strand → computes sense_fraction → Beta-Binomial κ

Both need SS. Path A's StrandModels produce it. Path B's calibration *consumes* it.

### 3.4 Fragment Length (Partial Redundancy)
- Path A: collects FL from unique mappers (spliced, unspliced, intergenic, same/opp strand) → trains FragmentLengthModels
- Path B: collects FL from unspliced single-region unique mappers → builds gDNA FL model (γ-weighted)

Path B needs a subset of Path A's FL observations.

### 3.5 Interval-Overlap Queries (Different Targets, Same Mechanism)
- Path A: queries the **transcript/exon interval tree** (cgranges `cr_`) during `_resolve_core()` in C++
- Path B: queries the **region partition** (cgranges `region_cr`) during `count_fragment()` in Python

Both use cgranges. The interval trees are different but both are loaded from the same index.

---

## 4. Integration Strategy: Unified Single-Pass Architecture

### Core Insight
The C++ BAM scanner already sees every fragment and has access to per-fragment genomic coordinates. The region partition (`region_cr`) is already built during index loading and stored in `TranscriptIndex.region_cr`. **Region accumulation can be done inside the C++ scan loop alongside fragment resolution, at negligible incremental cost.**

### Target Architecture
```
BAM → C++ BamScanner (enhanced)
    ├── Fragment resolution → FragmentAccumulator (unchanged)
    ├── Model training → StrandObs + FragLenObs (unchanged)
    └── Region accumulation → RegionAccumulator (NEW)
          ├── per-region: [n_unspliced_pos, n_unspliced_neg, n_spliced_pos, n_spliced_neg]
          └── FL observations: [(region_id, frag_len), ...]

→ Finalize models (strand_specificity, FL models)
→ calibrate_gdna(region_counts, fl_table, region_df, SS) → GDNACalibration
→ Apply calibration outputs to pipeline:
    ├── gdna_density_global / per_ref → replace intergenic density in EB prior
    ├── kappa → replace auto-estimated / hardcoded κ
    ├── gdna_fl_model → replace global FL for gDNA candidate scoring
    └── region_posteriors → optional per-locus gDNA weighting
→ FragmentScorer → FragmentRouter → ScoredFragments
→ build_loci() → compute_eb_gdna_priors() (enhanced) → batch_locus_em()
```

---

## 5. Phased Implementation Plan

### Phase 1: C++ Region Accumulation in BAM Scanner ✅ COMPLETE

**Goal**: Eliminate the pysam-based `region_evidence.py` BAM traversal by accumulating region evidence inside the existing C++ scan loop.

**Status**: Implemented and validated (2026-03-18). All 1029 tests pass. 15 comparison tests confirm C++/Python numerical agreement.

**Implementation Summary**:

#### 1a. `RegionAccumulator` struct (resolve_context.h, lines 1137–1337) ✅

~200-line struct with:
- `counts`: `std::vector<double>` — `[n_regions × 4]` row-major (unspliced_pos, unspliced_neg, spliced_pos, spliced_neg)
- `fl_region_ids`, `fl_frag_lens`: `std::vector<int32_t>` — FL observations for unspliced, single-region, unique mappers
- `region_cr`: borrowed `cgranges_t*` pointer (read-only, thread-safe)
- `rgn_buf` / `rgn_buf_cap`: per-instance scratch buffer for cgranges queries
- `id_to_ref`: borrowed pointer to ref name lookup
- Move-only (non-copyable) with proper RAII for scratch buffer
- `accumulate()`: fractional overlap counting via cgranges, mirrors `region_evidence.py::count_fragment()`
- `merge_from()`: thread-safe merge via element-wise addition + FL concatenation
- `init()`: sets up counts array and borrows cgranges/ref pointers

**Bug found and fixed**: cgranges stores coordinates as packed `uint64_t x` after indexing:
`x = (start << 32) | end`. Must use `cr_st(iv)` / `cr_en(iv)` accessors, NOT `iv->x` / `iv->y` directly. The `y` field is repurposed for the augmented interval tree after `cr_index()`.

#### 1b. `build_region_index()` on FragmentResolver (resolve_context.h + resolve.cpp) ✅

Added `build_region_index(refs, starts, ends, ids)` method to `FragmentResolver`.
- Builds a separate `region_cr_` cgranges (owned by resolver, freed in destructor)
- Tracks `n_regions_` for array sizing
- Exposed via nanobind binding in `resolve.cpp`
- Also added: `has_region_index()`, `n_regions()` accessors

#### 1c. Region accumulation in `process_qname_group_threaded` (bam_scanner.cpp) ✅

Two accumulation sites, both guarded by `region_acc.enabled()`:
1. **Intergenic fragments** (resolve fails): uses `frag.exons` and `frag.has_introns()`
2. **Resolved non-chimeric fragments**: uses `frag.exons` and boolean `frag_spliced` from `result.splice_type`

Per-worker `WorkerState` contains its own `RegionAccumulator` (thread-local). Workers are initialized with `region_acc.init()` if the resolver has a region index. After scan, worker accumulators merge into `BamScanner::region_acc_` in the merge loop.

#### 1d. Region evidence in `build_result()` (bam_scanner.cpp) ✅

Result dict contains `result["region_evidence"]` when data was accumulated:
```python
result["region_evidence"] = {
    "counts": list[float],       # flat n_regions*4 doubles
    "n_regions": int,
    "fl_region_ids": list[int],
    "fl_frag_lens": list[int],
}
```

#### 1e. `scan_and_buffer()` updated (pipeline.py) ✅

- Before scanner creation: calls `resolve_ctx.build_region_index()` with `index.region_df` arrays
- After scan: extracts `result["region_evidence"]` → `region_counts` DataFrame and `fl_table` DataFrame
- Return signature extended: `tuple[..., pd.DataFrame | None, pd.DataFrame | None]`
- `run_pipeline()` call site updated to unpack `region_counts, fl_table`

**Files modified**:
- `src/rigel/native/resolve_context.h` — `RegionAccumulator` struct + `build_region_index()` + accessors
- `src/rigel/native/bam_scanner.cpp` — `WorkerState.region_acc`, merge, accumulation sites, `build_result()`
- `src/rigel/native/resolve.cpp` — nanobind binding for `build_region_index`
- `src/rigel/pipeline.py` — `scan_and_buffer()` region index setup + evidence extraction + return type
- `tests/test_region_evidence_cpp.py` — 15 comparison tests (11 C++ vs Python, 4 standalone)

---

### Phase 2: Wire Calibration into the Pipeline ✅ COMPLETE

**Goal**: Call `calibrate_gdna()` between model finalization and fragment scoring, feeding its outputs into subsequent pipeline stages.

**Status**: Implemented and validated (2026-03-19). All 1050 tests pass (1029 existing + 21 new integration tests).

**Implementation Summary**:

#### 2a. `CalibrationConfig` dataclass (config.py) ✅

New frozen `CalibrationConfig` with fields: `enabled` (default `True`), `max_iterations`, `convergence_tol`, `density_percentile`, `min_gdna_regions`. Added to `PipelineConfig` as `calibration: CalibrationConfig`.

#### 2b. Calibration call in `run_pipeline()` (pipeline.py) ✅

After model finalization and before `quant_from_buffer()`, `calibrate_gdna()` is called when:
- `cal_cfg.enabled` is True
- `region_counts`, `fl_table`, and `index.region_df` are all available

The `GDNACalibration` result is stored on `PipelineResult.calibration` and passed into `quant_from_buffer()`.

#### 2c. Threading calibration through the call chain ✅

- `quant_from_buffer()` accepts `calibration: GDNACalibration | None = None`
- `_compute_priors()` accepts and passes `calibration` to `compute_eb_gdna_priors()`
- `compute_eb_gdna_priors()` in `locus.py` accepts `calibration: GDNACalibration | None = None`

When calibration is provided (calibrated path):
- `global_density = cal.gdna_density_global`
- `ref_shrunk = cal.gdna_density_per_ref`
- `k_ref = cal.kappa`
- `k_locus = cal.kappa` (unless `kappa_locus` explicitly set)
- All strand-based deconvolution at global/ref level is bypassed
- Per-locus shrinkage toward calibrated refs is preserved

When calibration is `None` (uncalibrated path): original strand-based behavior preserved exactly.

#### 2d. Golden output regeneration ✅

Golden files regenerated with `--update-golden` to reflect calibration-enabled defaults. Verified uncalibrated path produces bit-exact match with pre-Phase-2 golden files.

**Files modified**:
- `src/rigel/config.py` — new `CalibrationConfig`, added to `PipelineConfig`
- `src/rigel/pipeline.py` — `PipelineResult.calibration`, calibration call in `run_pipeline()`, threading through `quant_from_buffer()` and `_compute_priors()`
- `src/rigel/locus.py` — `compute_eb_gdna_priors()` refactored with calibrated/uncalibrated paths
- `tests/test_calibration_integration.py` — 21 new tests (config, end-to-end, comparison, clean scenario)
- `tests/golden/` — regenerated golden output files

#### 2a. Add calibration call to `run_pipeline()`

After `strand_models.finalize()` / `frag_length_models.finalize()` and before `quant_from_buffer()`:

```python
# --- Calibration (post-model, pre-scoring) ---
if index.region_df is not None and region_counts is not None:
    from .calibration import calibrate_gdna
    cal = calibrate_gdna(
        region_counts,
        fl_table,
        index.region_df,
        strand_models.strand_specificity,
    )
    logger.info(
        f"[CAL] gDNA cal: π={cal.mixing_proportion:.3f}, "
        f"κ={cal.kappa:.1f}, λ_G={cal.gdna_density_global:.2e}"
    )
```

#### 2b. Add `GDNACalibration` to `PipelineConfig` or as separate config

Add optional calibration parameters to `PipelineConfig`:
```python
@dataclass(frozen=True)
class CalibrationConfig:
    enabled: bool = True
    max_iterations: int = 50
    convergence_tol: float = 1e-4
    density_percentile: float = 0.10
    min_gdna_regions: int = 100
```

#### 2c. Thread `GDNACalibration` outputs into `quant_from_buffer()`

The calibration outputs replace/augment data at specific points:

| Calibration Output | Replaces | Integration Point |
|---|---|---|
| `gdna_density_global` | `_compute_intergenic_density()` return | `_compute_priors()` → `intergenic_density` kwarg |
| `gdna_density_per_ref` | `ref_shrunk` dict in `compute_eb_gdna_priors` | `_compute_priors()` → new kwarg |
| `kappa` | Auto-estimated `k_ref` and `k_locus` | `_compute_priors()` → `kappa_ref`/`kappa_locus` kwargs or new calibrated path |
| `gdna_fl_model` | Global FL model used for gDNA scoring | `FragmentScorer` → new `gdna_fl_model` param |
| `mixing_proportion` | N/A (new information) | `gdna_prior_scale` tuning or diagnostic |
| `region_posteriors` | N/A (new information) | Optional: per-locus weighted gDNA initialization |

**Testing**: Integration tests comparing pipeline output with and without calibration. Golden output regeneration.

**Files modified**:
- `src/rigel/config.py` — new `CalibrationConfig` dataclass
- `src/rigel/pipeline.py` — calibration call + output threading
- `src/rigel/locus.py` — `compute_eb_gdna_priors()` accepts calibrated density + κ
- `tests/test_pipeline_smoke.py` — verify calibration runs end-to-end

---

### Phase 3: Replace EB gDNA Density with Calibrated Density ✅ COMPLETE

**Goal**: Use calibration's per-reference gDNA density as the authoritative density signal, replacing the crude intergenic density in the per-locus hybrid estimator.

**Status**: Implemented and validated (2026-03-19). All 1067 tests pass (1029 original + 21 Phase 2 + 17 Phase 3).

**Implementation Summary**:

#### 3a. Refactored `compute_eb_gdna_priors()` ✅ (done in Phase 2)

The function was refactored in Phase 2 to accept `calibration: GDNACalibration | None` and use calibrated global density, per-ref densities, and κ in the calibrated path.

#### 3b. Per-locus density anchor improved ✅

`_compute_per_locus_gdna_densities()` now accepts an optional `calibrated_ref_densities: dict[str, float] | None` parameter. When provided, each locus's `compute_gdna_density_hybrid()` call uses the calibrated per-ref density for its primary reference as the density anchor, instead of the crude intergenic density scalar.

This directly implements the Phase 3 design:
```
locus_density = W × strand_density + (1-W) × calibrated_ref_density
```

**Impact**: Most significant for weakly-stranded libraries (SS ≈ 0.5–0.7) where W = (2SS-1)² → 0 and the density component dominates. At high SS, the strand signal dominates and the anchor has minimal impact.

The fallback behavior is preserved: when the locus's reference isn't found in `calibrated_ref_densities`, `intergenic_density` is used.

#### 3c. Calibrated κ ✅ (done in Phase 2)

The calibrated κ from `calibrate_gdna()` replaces both `k_ref` and `k_locus` in the calibrated path.

**Testing**: 17 new tests in `tests/test_calibrated_density.py`:
- `TestHybridDensityAnchor` (4 tests): validates density component behavior at different SS levels
- `TestPerLocusCalibratedDensity` (7 tests): unit tests for `_compute_per_locus_gdna_densities()` with/without calibrated ref densities, including multi-ref, missing-ref fallback, and parent independence
- `TestCalibratedDensityEndToEnd` (6 tests): A/B pipeline comparison at SS=0.98 and SS=0.65

**Files modified**:
- `src/rigel/locus.py` — `_compute_per_locus_gdna_densities()` accepts `calibrated_ref_densities`, `compute_eb_gdna_priors()` passes calibrated refs when available
- `tests/test_calibrated_density.py` — 17 new tests
- `tests/golden/` — regenerated golden output files

---

### Phase 4: gDNA Fragment-Length Model in Scoring ✅ COMPLETE

**Goal**: Use the calibrated gDNA FL model for gDNA candidate scoring instead of the strand-deconvolved intergenic-based model.

**Status**: Implemented and validated (2026-03-19). All 1085 tests pass.

**Key Discovery**: The C++ scoring infrastructure already supports separate gDNA FL scoring (`gdna_frag_len_log_lik()`, `gdna_fl_log_prob_`, etc.) via `NativeFragmentScorer`. `FragmentScorer.from_models()` already reads `frag_length_models.gdna_model` for the gDNA FL lookup table. No C++ changes were needed — only the Python pipeline needed to flow the calibrated model.

**Implementation Summary**:

#### 4a. Model replacement in `quant_from_buffer()` ✅

Before calling `_score_fragments()`, the pipeline replaces `frag_length_models.gdna_model` with `calibration.gdna_fl_model` when calibration produced a valid model (non-None). `FragmentScorer.from_models()` automatically picks up the replacement.

#### 4b. ESS Gate in `build_gdna_fl_model()` ✅

The calibrated gDNA FL model is gated by Effective Sample Size (ESS):

- `build_gdna_fl_model()` accepts a `min_ess` parameter and returns `None` when `Σγ·w < min_ess` (the γ-weighted sum of FL observations is too low to fit a reliable distribution).
- `CalibrationConfig.min_fl_ess` (default 50) controls this threshold and is threaded through `calibrate_gdna()`.
- `GDNACalibration.gdna_fl_model` can be `None`, and the pipeline simply checks `is not None` before replacement.

This is a principled statistical approach: 50 observations clustered around a peak are sufficient to estimate the mean and variance of a unimodal distribution, but fewer are unreliable.

#### 4c. Region count bailout in `calibrate_gdna()` ✅

When the number of eligible regions is less than `min_gdna_regions` (default 100), the calibration EM lacks statistical power and bails out early with `gdna_fl_model=None`. This prevents the EM from running on tiny toy datasets where it would incorrectly assign RNA-origin fragments to gDNA.

#### 4d. C++ interface unchanged ✅

The existing C++ infrastructure handles gDNA FL scoring:
- `NativeFragmentScorer` stores `gdna_fl_log_prob_`, `gdna_fl_max_size_`, `gdna_fl_tail_base_`
- `gdna_frag_len_log_lik()` provides fast FL lookup for gDNA candidates
- `fused_score_buffer()` already uses the separate gDNA FL in gDNA log-likelihood computation

No C++ code changes or recompilation were required.

**Testing**: 18 new tests in `tests/test_gdna_fl_scoring.py`:
- `TestFragmentLengthModelReplacement` (4 tests): unit tests for model replacement mechanics
- `TestFragmentScorerGdnaFL` (3 tests): verify scorer picks up replaced model
- `TestGdnaFLEndToEnd` (5 tests): pipeline integration with calibration on/off
- `TestESSGate` (4 tests): ESS gate logic and pipeline behavior with low ESS
- `TestGdnaFLScoringImpact` (2 tests): end-to-end with distinctive gDNA FL

**Files modified**:
- `src/rigel/config.py` — `CalibrationConfig.min_fl_ess` parameter
- `src/rigel/calibration.py` — ESS gate in `build_gdna_fl_model()`, region count bailout, `min_fl_ess` threading
- `src/rigel/pipeline.py` — simple `None` check for calibrated gDNA FL replacement
- `tests/test_gdna_fl_scoring.py` — 18 tests (ESS gate replaces former discrimination guard tests)
- `tests/golden/` — regenerated golden output files

---

### Phase 5: Calibration-Informed EM Initialization

**Goal**: Use calibration's per-region posteriors to improve EM initialization beyond density-only priors.

**Changes**:

#### 5a. Region-to-locus posterior aggregation

Each locus spans multiple regions (via its transcripts' genomic coordinates). The calibration posteriors (γ_r) give P(not expressed) per region. Aggregate to per-locus:

```python
def compute_locus_gdna_weight(locus, region_posteriors, region_df, index):
    """Aggregate region posteriors into a locus-level gDNA weight."""
    # For each transcript in the locus, find overlapping regions
    # Weighted average of γ_r by overlap length
    ...
```

This gives a calibration-informed initial estimate of how much gDNA a locus should contain — independent of the strand deconvolution.

#### 5b. Warm-start EM θ_gDNA with calibration weight

Currently, `compute_ovr_prior_and_warm_start()` in `em_solver.cpp` warm-starts θ from coverage weights (OVR). The gDNA component starts from its prior alone.

With calibration, θ_gDNA can be warm-started:
```cpp
// In process_locus lambda:
double cal_gdna_weight = calibrated_gdna_weights[li];
theta_init[gdna_idx] = cal_gdna_weight * total_fragments;
```

This provides the EM with a much better initial point for gDNA, especially in mega-loci where the default initialization leads to θ accumulation runaway.

#### 5c. Adaptive prior scale from mixing proportion

The calibrated mixing proportion π gives the global gDNA contamination rate. This can inform the `gdna_prior_scale`:

```python
# If calibration says 10% gDNA, scale the prior accordingly
gdna_prior_scale = calibration.mixing_proportion / (1 - calibration.mixing_proportion)
```

**Testing**: EM convergence speed comparison. gDNA accuracy on synthetic data with known ground truth.

**Files modified**:
- `src/rigel/locus.py` — new `compute_locus_gdna_weight()`
- `src/rigel/native/em_solver.cpp` — accept calibrated warm-start
- `src/rigel/estimator.py` — pass calibration weights through

---

### Phase 6: Cleanup and Deprecation

**Goal**: Remove redundant code paths, make calibration the default.

**Changes**:

#### 6a. Deprecate `region_evidence.py`

Once C++ region accumulation (Phase 1) is verified to match, the Python pysam-based implementation becomes a test reference. Move to `tests/` or mark as deprecated.

#### 6b. Simplify `compute_eb_gdna_priors()` fallback path

The original strand-only EB path remains as a fallback when region_df is unavailable (e.g., old index without regions). But the calibrated path becomes the default.

#### 6c. Remove redundant config parameters

With calibration providing κ, the `gdna_kappa_ref` / `gdna_kappa_locus` / `gdna_kappa_min` / `gdna_kappa_max` / `gdna_kappa_fallback` parameters in `EMConfig` become less critical. They can be simplified or defaulted from calibration.

#### 6d. Update golden outputs

All golden output files regenerated with calibration enabled.

**Files modified**:
- `src/rigel/region_evidence.py` — deprecated / moved
- `src/rigel/config.py` — simplified EB config
- `tests/golden/` — regenerated

---

## 6. Synergies and Coalescing Opportunities

### 6.1 Single-Pass BAM Scan (Phase 1)
**Impact**: Major. Eliminates the most expensive redundancy (second BAM traversal). The C++ scanner processes fragments at ~5-10M/min; pysam at ~0.5-1M/min. Single-pass gives 2-10× faster calibration.

### 6.2 Shared Strand Specificity (Phase 2)
**Impact**: Medium. The calibration EM needs strand_specificity, which Path A already computes. No need to re-estimate from region data.

### 6.3 Shared FL Observations (Phase 1+4)
**Impact**: Medium. The C++ scanner already collects unspliced same/opp strand FLs. The calibration needs unspliced single-region FLs. Collecting the latter in C++ is a trivial extension of the existing FL collection.

### 6.4 Region cgranges in C++ (Phase 1)
**Impact**: Medium. The C++ resolver already uses cgranges for transcript overlap. Adding a second cgranges for regions reuses the same infrastructure. The cgranges query during region accumulation is lightweight (~2-5 lookups per fragment, each O(log n + k)).

### 6.5 Calibrated κ Unifies Two Estimators (Phase 3)
**Impact**: High. Currently there are THREE κ estimation methods:
1. `priors.estimate_kappa()` — MoM from density scatter (used by EB system)
2. `calibration.estimate_kappa_marginal()` — MLE from strand mixture (used by calibration)
3. Hardcoded defaults in config

With calibration integrated, method 2 (the most principled) replaces methods 1 and 3 for the strand symmetry parameter. Method 1 can be retained for shrinkage at the locus level.

### 6.6 Density Anchor Eliminates Circular Dependency (Phase 3)
**Impact**: High. The current EB system estimates gDNA density *during* the scoring/EM pipeline, which creates a subtle circularity (gDNA density depends on which fragments are "gDNA", which depends on EM output, which depends on gDNA density as prior). Calibration estimates density *before* EM from region-level evidence, breaking the loop.

### 6.7 gDNA FL Model Improves Scoring Discrimination (Phase 4)
**Impact**: High. Currently, gDNA and mRNA/nRNA candidates use the SAME global FL model for scoring. Since the gDNA FL distribution differs (it reflects insert-size of DNA, not mRNA transcript-constrained FL), using a gDNA-specific FL model improves the log-likelihood ratio between components, directly improving EM discrimination.

### 6.8 Warm-Start Mitigates Mega-Locus θ Accumulation (Phase 5)
**Impact**: Critical for problematic samples. The mega-locus gDNA siphoning occurs because θ_gDNA starts from prior alone and accumulates from ALL unspliced fragments. A calibration-informed warm start gives θ_gDNA a reasonable initial value, reducing the number of EM iterations needed and preventing runaway accumulation.

---

## 7. Data Flow Diagram: Before and After

### Before (Two Passes)
```
                    ┌─── C++ BamScanner ──────────────────────┐
                    │  Per-fragment:                           │
     BAM ──(1)────►│  - Transcript resolution                │──► FragmentBuffer
                    │  - StrandObs training                   │        │
                    │  - FragLenObs training                  │        │
                    └─────────────────────────────────────────┘        │
                                                                       ▼
                    ┌─── pysam iteration ─────────────────────┐    FragmentScorer
                    │  Per-fragment:                           │        │
     BAM ──(2)────►│  - Region overlap (cgranges)            │        ▼
                    │  - Fractional counting                  │    ScoredFragments
                    │  - FL observation tabulation            │        │
                    └──────────────┬──────────────────────────┘        ▼
                                   │                              build_loci()
                                   ▼                                   │
                            calibrate_gdna()                           ▼
                                   │                           compute_eb_gdna_priors()
                                   │    [NOT CONNECTED]        (strand-based only)
                                   ▼                                   │
                            GDNACalibration                            ▼
                            (unused by main)                   batch_locus_em()
```

### After (Single Pass, Unified)
```
                    ┌─── C++ BamScanner (enhanced) ───────────┐
                    │  Per-fragment:                           │
                    │  - Transcript resolution (unchanged)    │──► FragmentBuffer
     BAM ──(1)────►│  - StrandObs training (unchanged)       │
                    │  - FragLenObs training (unchanged)      │
                    │  - Region accumulation (NEW)            │──► region_counts + fl_table
                    └─────────────────────────────────────────┘
                                                                       │
                        ┌──────────────────────────────────────────────┘
                        │
                        ▼
        ┌── Finalize Models ──┐     ┌── Calibration ──────────────┐
        │ strand_specificity  │────►│ calibrate_gdna()            │
        │ frag_length_models  │     │ → density, κ, FL, γ, π     │
        └─────────────────────┘     └─────────┬───────────────────┘
                                              │
        ┌─── Scoring ────────────────────────┤
        │ FragmentScorer(gdna_fl_model=...)  │◄── gDNA FL model
        │ FragmentRouter → ScoredFragments   │
        └──────────┬─────────────────────────┘
                   │
                   ▼
        ┌─── Locus Construction + EB Priors ──────────────────┐
        │ build_loci()                                         │
        │ compute_eb_gdna_priors(calibration=cal)             │
        │   → calibrated density replaces intergenic estimate │
        │   → calibrated κ replaces MoM / hardcoded κ        │
        │   → shrinkage: locus → calibrated_ref → cal_global │
        └──────────┬───────────────────────────────────────────┘
                   │
                   ▼
        ┌─── batch_locus_em() ────────────────────────────────┐
        │ - Calibrated gDNA priors                            │
        │ - Calibrated θ_gDNA warm start (optional, Phase 5) │
        │ - gDNA FL-aware scoring                             │
        └─────────────────────────────────────────────────────┘
```

---

## 8. Risk Assessment and Mitigations

### Risk 1: C++ Region Accumulation Diverges from Python Reference
**Mitigation**: Add pairwise comparison tests in `test_region_evidence.py` that run both C++ and Python paths on the same BAM and assert agreement within tolerance. Keep Python implementation as test oracle until Phase 6.

### Risk 2: Region cgranges Construction Cost
**Mitigation**: `region_cr` is already built during `TranscriptIndex.load()`. Passing it to C++ adds zero runtime. The per-fragment cgranges query cost is O(log n + k) where k is typically 1-3 regions — negligible vs. the transcript interval query cost already paid.

### Risk 3: Calibration EM Non-Convergence on Edge Cases
**Mitigation**: `calibrate_gdna()` already handles edge cases: no-data bailout, pristine sample detection, minimum seed guarantees. The pipeline falls back to the original EB path when calibration fails or reports `converged=False`.

### Risk 4: Mega-Locus Regime Not Fully Addressed
**Mitigation**: Calibration addresses the *input* side (better density, κ, FL), but mega-locus θ accumulation is fundamentally an EM structural problem. Calibration significantly reduces the severity (correct density anchor prevents runaway), but may not fully eliminate it. A complementary fix (per-gene gDNA components or EM decomposition) may still be needed for the worst cases.

### Risk 5: Backward Compatibility
**Mitigation**: Calibration requires `regions.feather` in the index. Old indices without it fall back to the original EB path. The `CalibrationConfig.enabled` toggle allows users to disable calibration.

---

## 9. Testing Strategy

### Phase 1 Tests
- **Unit**: `RegionAccumulator` C++ tests (fractional counting correctness)
- **Integration**: Compare C++ vs Python region evidence on all scenario BAMs
- **Regression**: Existing golden outputs unchanged (region accumulation is additive)

### Phase 2 Tests
- **Smoke**: `test_pipeline_smoke.py` runs with calibration enabled
- **Integration**: Verify `GDNACalibration` object is populated in pipeline results

### Phase 3 Tests
- **A/B comparison**: Synthetic BAMs with known gDNA fractions — measure accuracy improvement
- **Regression**: Golden outputs regenerated with calibrated priors

### Phase 4 Tests
- **Unit**: gDNA FL log-likelihood lookup correctness
- **Integration**: Scoring output changes when gDNA FL differs from global FL

### Phase 5 Tests
- **EM convergence**: Verify fewer iterations needed with warm start
- **Accuracy**: gDNA count accuracy on stress-test scenarios

### Phase 6 Tests
- **Full suite**: All 1000+ tests pass with calibration as default
- **Backward compatibility**: Old indices (no regions.feather) still work

---

## 10. Implementation Priority and Dependencies

```
Phase 1 (C++ region accumulation)
    │
    ├──► Phase 2 (Wire calibration into pipeline) ─── depends on Phase 1
    │        │
    │        ├──► Phase 3 (Replace EB density) ──────── depends on Phase 2
    │        │
    │        ├──► Phase 4 (gDNA FL in scoring) ──────── depends on Phase 2
    │        │
    │        └──► Phase 5 (EM warm start) ───────────── depends on Phase 2+3
    │
    └──► Phase 6 (Cleanup) ─────────────────────────── depends on Phase 3+4+5
```

Phases 3, 4, and 5 are **independent of each other** and can be developed in parallel after Phase 2. Phase 6 should only happen after all prior phases are stable.

**Recommended implementation order**: 1 → 2 → 3 → 4 → 5 → 6

---

## 11. Estimated Scope per Phase

| Phase | New/Modified Files | C++ Changes | Python Changes | Risk |
|-------|--------------------|-------------|----------------|------|
| 1 | 5 | Major (new accumulator) | Minor (extract results) | Medium |
| 2 | 4 | None | Medium (wiring) | Low |
| 3 | 2 | None | Medium (refactor EB) | Low |
| 4 | 3 | Medium (FL lookup) | Medium (scorer) | Medium |
| 5 | 3 | Minor (warm start) | Medium (aggregation) | Medium |
| 6 | 5 | None | Minor (cleanup) | Low |

---

## 12. Open Questions for Review

1. **Should calibration be optional (toggle) or always-on?**
   Recommendation: Always-on when `regions.feather` exists in the index. This is the safest default — calibration strictly improves the information available to EM.

2. **Should calibration parameters be user-configurable?**
   Recommendation: Add a minimal `CalibrationConfig` (enabled, max_iterations). The internal parameters (density_percentile, min_gdna_regions) should remain defaults unless diagnostic investigation is needed.

3. **Should the gDNA FL model be used only for scoring (Phase 4) or also in the EM effective length correction?**
   Recommendation: Start with scoring only. The EM effective length system is already 1.0 for all components with per-fragment bias profiles handling the correction. Changing this is a deeper architectural shift.

4. **How to handle the `strand_symmetry_kappa` parameter in `BamScanConfig` vs calibrated κ?**
   Recommendation: The calibrated κ from `calibration.py` serves a different role (region-level strand classification) than the EM-side κ (per-locus gDNA coupling). They should be named differently. The EB shrinkage κ parameters (`gdna_kappa_ref`, `gdna_kappa_locus`) should be overridable by calibration output.

5. **Per-region posteriors for per-locus initialization (Phase 5) — use γ directly or only as density anchor?**
   Recommendation: Start with density anchor only (simpler, less risky). Per-region posteriors can later inform per-locus initialization if the density anchor alone is insufficient.
