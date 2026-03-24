# Phase 4 — Quantification-Stage Performance Results

**Date**: 2026-03-24
**Status**: Complete
**Prerequisite**: Phase 3c (streaming chunk architecture) complete

---

## Summary

Phase 4 optimized the Python quantification stages that dominated wall time
after Phase 3c's C++ scan improvements. Three sub-phases targeted specific
bottlenecks identified by cProfile on a real BAM (1.6M fragments, 457K
transcripts):

| Phase | Target | Before | After | Savings |
|-------|--------|--------|-------|---------|
| **4a** | Vectorized interval merge + cache | 4.60s | 0.73s | **3.87s (84%)** |
| **4b** | Bulk exon CSR construction | 2.44s | 0.97s | **1.47s (60%)** |
| **4c** | Dead code + gc.collect removal | 0.69s | 0.00s | **0.69s (100%)** |
| **Total** | | | | **6.03s** |

---

## Stage-by-Stage Comparison

| Stage | Phase 3c | Phase 4 | Change |
|-------|----------|---------|--------|
| **wall_time** | **22.70s** | **17.00s** | **-25.1%** |
| scan_and_buffer | 9.12s | 8.99s | -1.4% |
| calibration | 2.51s | 2.40s | -4.4% |
| **fragment_scorer** | **2.44s** | **0.97s** | **-60.2%** |
| **build_loci** | **2.21s** | **0.43s** | **-80.5%** |
| **eb_gdna_priors** | **2.39s** | **0.30s** | **-87.4%** |
| locus_em | 3.72s | 3.64s | -2.2% |
| fragment_router_scan | 0.29s | 0.27s | -6.9% |
| peak_rss_mb | 4039 MB | 3910 MB | -3.2% |
| **throughput** | **261K f/s** | **581K f/s** | **+2.23x** |

---

## Correctness Verification

- **Unit tests**: 1018 passed (0 failures)
- **Golden output tests**: 21 passed (0 regressions)

---

## Phase 4a — Vectorized Interval Merging & Locus Footprints

### What changed

- Deleted `_merged_intervals()` generator (69K calls, 5.3s cumulative)
- Pre-extracted Arrow-backed string column to numpy object array once
  (eliminated 2.2s `pyarrow.compute.take` overhead)
- Factorized ref strings to integer codes for fast sort/compare
- Merged intervals computed once per locus in `build_loci()`, cached
  on `Locus.merged_intervals`
- `compute_gdna_locus_gammas()` uses cached intervals (no recomputation)

### Files modified

- `src/rigel/scored_fragments.py` — Added `merged_intervals` field to `Locus`
- `src/rigel/locus.py` — Rewrote `build_loci()`, simplified
  `compute_gdna_locus_gammas()`, deleted `_merged_intervals()` and
  `defaultdict` import

### cProfile evidence

| Function | Before | After |
|----------|--------|-------|
| `_merged_intervals` | 5.32s (69K calls) | deleted |
| `pyarrow.compute.take` | 2.21s | eliminated |
| `build_loci` | 2.21s | 0.43s |
| `compute_gdna_locus_gammas` | 2.39s | 0.004s |

---

## Phase 4b — Bulk Exon CSR Construction

### What changed

- Added `TranscriptIndex.build_exon_csr()` — builds CSR arrays from the
  existing `_t_exon_intervals` dict in two vectorized passes
- Replaced 457K-iteration Python loop in `FragmentScorer.from_models()`
  with a single `build_exon_csr()` call
- Changed C++ `NativeFragmentScorer` constructor from dict-unpacking
  (~55 lines) to direct array copy (~20 lines)
- Removed `_t_exon_data` field from `FragmentScorer`

### Files modified

- `src/rigel/index.py` — Added `build_exon_csr()` method
- `src/rigel/scoring.py` — Removed 457K loop, removed `_t_exon_data` field
- `src/rigel/native/scoring.cpp` — Simplified constructor, updated bindings

### cProfile evidence

| Function | Before | After |
|----------|--------|-------|
| Per-transcript loop | 2.44s (457K iters) | deleted |
| `numpy.cumsum` | 0.47s (457K calls) | 0.16s (228K calls) |
| `build_exon_csr` | — | 0.012s |

---

## Phase 4c — Dead Code & Hygiene

### What changed

- Deleted `merge_accumulator_into()` — 53-line dead function in
  `bam_scanner.cpp`, never called since Phase 3b
- Removed `gc.collect()` after scoring in `pipeline.py` — saved 0.69s;
  CPython refcount GC handles the `del` + `release()` deterministically
- Kept `gc.collect()` after EM (pipeline.py L602) for final output stage

### Files modified

- `src/rigel/native/bam_scanner.cpp` — Deleted dead function
- `src/rigel/pipeline.py` — Removed `gc.collect()` after scoring

---

## What's NOT in Phase 4

- **C++ scan optimization** — requires py-spy profiling; separate investigation
- **C++ EM solver** — 3.6s, well-optimized with SQUAREM + OpenMP
- **Calibration optimizer** — 2.4s, dominated by scipy minimize calls
