# Calibration Integration — Code-Level Reference

**Companion to**: `implementation_plan.md`
**Purpose**: Maps each integration point to specific files, functions, and line ranges.

---

## Current Calibration Pipeline (Standalone)

### Entry Point
| File | Function | Lines | Role |
|------|----------|-------|------|
| `src/rigel/region_evidence.py` | `count_region_evidence()` | 232-290 | BAM → region counts + FL table (pysam) |
| `src/rigel/calibration.py` | `calibrate_gdna()` | 973-1211 | Main calibration EM orchestrator |

### Calibration Internal Functions
| Function | Lines | Role |
|----------|-------|------|
| `compute_region_stats()` | 165-237 | Per-region strand_ratio, density, splice_rate |
| `compute_sense_fraction()` | 244-278 | Gene-strand-corrected sense fraction |
| `_compute_strand_llr_binomial()` | 281-316 | Simple Binomial strand LLR (fallback) |
| `_compute_strand_llr_betabinom()` | 319-395 | Beta-Binomial strand LLR with shared κ |
| `estimate_kappa_marginal()` | 400-500 | MLE κ via marginal mixture likelihood |
| `compute_log_density()` | 594-618 | Log-density with ε pseudocount |
| `_compute_density_llr_gaussian()` | 621-641 | Gaussian density LLR per region |
| `_compute_fl_llr()` | 834-889 | Shape-normalized FL LLR per region |
| `_seed_initial_partition()` | 674-744 | Two-phase initialization (expressed + gDNA seeds) |
| `_e_step()` | 751-830 | Per-region posterior γ with hard constraints |
| `_m_step()` | 895-948 | Update π, λ_G, λ_E, Gaussian params |
| `build_gdna_fl_model()` | 548-589 | γ-weighted gDNA FL histogram |
| `_compute_per_ref_density()` | 521-542 | Per-reference γ-weighted density |
| `GDNACalibration` | 107-155 | Result container dataclass |

### Calibration Dependencies
| Dependency | Source | Used For |
|------------|--------|----------|
| `region_df` | `TranscriptIndex.region_df` (from `regions.feather`) | Region metadata (length, flags, ref) |
| `region_cr` | `TranscriptIndex.region_cr` (cgranges) | Spatial overlap queries |
| `strand_specificity` | `StrandModels.strand_specificity` | Strand channel SS parameter |
| `FragmentLengthModel` | `src/rigel/frag_length_model.py` | gDNA + RNA FL models |

---

## Current gDNA Prior Pipeline (Main Tool)

### Entry Points
| File | Function | Lines | Role |
|------|----------|-------|------|
| `src/rigel/pipeline.py` | `_compute_priors()` | 441-486 | Orchestrates nRNA init + EB gDNA priors |
| `src/rigel/locus.py` | `compute_eb_gdna_priors()` | 763-950 | Hierarchical EB gDNA initialization |

### EB Prior Internal Functions
| File | Function | Lines | Role |
|------|----------|-------|------|
| `locus.py` | `compute_gdna_density_from_strand()` | 588-610 | Strand deconvolution: G / exonic_bp |
| `locus.py` | `compute_gdna_density_hybrid()` | 613-660 | Hybrid density+strand with W=(2SS-1)² |
| `locus.py` | `_compute_ref_gdna_densities()` | 663-694 | Per-ref shrinkage toward global |
| `locus.py` | `_compute_per_locus_gdna_densities()` | 697-757 | Per-locus densities + parent refs |
| `priors.py` | `estimate_kappa()` | (MoM) | Shrinkage κ estimation from density scatter |
| `pipeline.py` | `_compute_intergenic_density()` | 100-135 | Intergenic density anchor |

### EB Prior Data Flow
```
estimator.transcript_unspliced_sense/antisense
    → compute_gdna_density_hybrid() [per transcript, per ref, per locus]
    → MoM estimate_kappa() for κ_ref, κ_locus
    → hierarchical shrinkage: locus → ref → global
    → gdna_init = shrunk_density × exonic_bp_locus
```

---

## Integration Touch Points (by Phase)

### Phase 1: C++ Region Accumulation ✅ COMPLETE

**C++ code added**:
- `src/rigel/native/resolve_context.h`:
  - `RegionAccumulator` struct (lines 1137–1337): counts, FL obs, accumulate, merge, init
  - `FragmentResolver::build_region_index()` (lines 530–553): builds `region_cr_` cgranges
  - `FragmentResolver::has_region_index()`, `n_regions()` accessors (lines 570–571)
  - Destructor updated to free `region_cr_` (line 524)
- `src/rigel/native/bam_scanner.cpp`:
  - `WorkerState` (line 338): `RegionAccumulator region_acc;` member
  - `merge_region_acc()` (line 438): static merge helper
  - `BamScanner::region_acc_` (line 947): class-level merged accumulator
  - Worker init (lines 1001–1007): `region_acc.init()` when region index available
  - Worker merge (lines 1118–1119): `merge_region_acc()` call
  - Intergenic accumulation (lines 1248–1249): `region_acc.accumulate()` for unresolved fragments
  - Resolved accumulation (lines 1351–1353): `region_acc.accumulate()` for non-chimeric resolved fragments
  - `build_result()` (lines 1436–1442): `result["region_evidence"]` dict
- `src/rigel/native/resolve.cpp` (lines 108–125): nanobind `.def("build_region_index", ...)`

**Python changes**:
- `src/rigel/pipeline.py::scan_and_buffer()`:
  - Return type extended to include `pd.DataFrame | None` × 2 (region_counts, fl_table)
  - Before scanner: `resolve_ctx.build_region_index()` with region_df arrays
  - After scan: extracts `result["region_evidence"]` → DataFrames (region_counts, fl_table)
- `src/rigel/pipeline.py::run_pipeline()`: updated unpack of `scan_and_buffer()` return

**Tests added**:
- `tests/test_region_evidence_cpp.py`: 15 tests
  - `TestCppVsPythonRegionEvidence` (11 tests): C++ vs Python comparison on identical BAMs
  - `TestCppRegionAccumulatorStandalone` (4 tests): shape, FL filtering, no-index guard

**Bug found during implementation**: cgranges `cr_intv_t.x` field stores packed
`(start << 32) | end` after `cr_index()`. Code initially accessed `iv->x` and `iv->y`
directly — must use `cr_st(iv)` / `cr_en(iv)` accessor macros.

### Phase 2: Wire Calibration

**Python changes**:
- `src/rigel/config.py` (after line ~155): add `CalibrationConfig` dataclass
- `src/rigel/config.py: PipelineConfig` (line 142): add `calibration: CalibrationConfig` field
- `src/rigel/pipeline.py: run_pipeline()` (after line ~730, post-model finalize):
  - Call `calibrate_gdna()` with region evidence
  - Store `GDNACalibration` result
- `src/rigel/pipeline.py: quant_from_buffer()` (line ~580): accept `calibration` kwarg
- `src/rigel/pipeline.py: _compute_priors()` (line 441): accept `calibration` kwarg

### Phase 3: Replace EB Density

**Python changes**:
- `src/rigel/locus.py: compute_eb_gdna_priors()` (line 763):
  - Add `calibration: GDNACalibration | None = None` parameter
  - When calibration is provided:
    - Skip global/ref strand deconvolution (lines 814-862)
    - Use `calibration.gdna_density_global` as `global_density`
    - Use `calibration.gdna_density_per_ref` as `ref_shrunk`
    - Use `calibration.kappa` for `k_ref` and potentially `k_locus`
  - Keep per-locus strand deconvolution (line 697+) but shrink toward calibrated ref

### Phase 4: gDNA FL Model

**C++ changes**:
- `src/rigel/native/scoring.cpp: fused_score_buffer()` (line ~1271):
  - Add `f64_1d gdna_fl_log_probs` parameter
  - In gDNA log-likelihood computation, use `gdna_fl_log_probs[fl_len]` instead of global FL

**Python changes**:
- `src/rigel/scoring.py: FragmentScorer` — add `gdna_fl_model` attribute
  - In `from_models()`: accept optional `gdna_fl_model` kwarg
  - Compute `gdna_fl_log_probs` array from model's count histogram
- `src/rigel/scan.py: FragmentRouter._scan_native()` — pass `gdna_fl_log_probs` to `fused_score_buffer()`

### Phase 5: EM Warm Start

**Python changes**:
- `src/rigel/locus.py` — new `compute_locus_gdna_weight(locus, region_posteriors, region_df, index)` function
- `src/rigel/pipeline.py: _compute_priors()` — compute locus-level calibration weights

**C++ changes**:
- `src/rigel/native/em_solver.cpp: batch_locus_em()` (line ~1606):
  - Add `f64_1d calibrated_gdna_weights` parameter (optional, can be zeros)
  - In `process_locus` lambda: use weight for θ_gDNA warm start

---

## Key Data Structures to Understand

### Region Partition (from `index.py: build_region_table`)
```
region_df columns: region_id, ref, start, end, length, exon_pos, exon_neg, tx_pos, tx_neg
```
- Non-overlapping bins cover entire reference set
- Flags from boundary-sweep algorithm
- Used by calibration to classify region types (intergenic, intronic, exonic)

### GDNACalibration (from `calibration.py`)
```python
region_posteriors: ndarray[n_regions]  # γ_r ∈ [0,1], P(not expressed)
gdna_density_global: float             # frags/bp from gDNA
gdna_density_per_ref: dict[str, float] # per-chromosome density
kappa: float                           # Beta-Binomial concentration
gdna_fl_model: FragmentLengthModel     # gDNA-specific FL distribution
mixing_proportion: float               # π = P(not expressed)
expressed_density: float               # frags/bp from RNA+gDNA
```

### Fragment Data Available in C++ Scanner (per fragment, at accumulation time)
```
From AssembledFragment:
  - exons: vector<ExonBlock> (ref, start, end, strand)
  - introns: vector<(ref, start, end)>
  - genomic_footprint(), has_introns()

From RawResolveResult:
  - t_inds, splice_type, exon_strand, sj_strand
  - frag_lengths[], exon_bp[], intron_bp[]

From scan context:
  - is_unique_mapper (from NH tag / secondary detection)
  - frag_id (sequential counter)
```

All of this is sufficient to replicate `count_fragment()` from `region_evidence.py`.
