# Memory Reduction Plan — Phase 5 Status Update

**Date:** 2026-03-06

## 1. Measured Results After Phase 5A + 5B

### Profiling Setup

- **Dataset:** 24M fragments (12M PE reads), simulated CCLE HeLa cervix
- **Transcripts:** 254,461 across 63,472 genes
- **Machine:** macOS arm64, 4 threads
- **Buffer spill:** 2 GiB budget; buffer was spilled to disk (0 MB in memory at scan time)

### Before vs After

| Metric | Baseline (pre-5A/5B) | After 5A+5B | Delta |
|---|---:|---:|---:|
| Peak RSS | 20,272 MB | 18,548 MB | **−1,724 MB (−8.5%)** |
| RSS before pipeline | 3,868 MB | 3,866 MB | ~0 |
| Wall time | 138.6s | 139.6s | ~0 |
| Throughput | 173K frags/s | 174K frags/s | ~0 |

### Stage Breakdown (--stages mode, 18,693 MB peak)

| Stage | Time | % |
|---|---:|---:|
| scan_and_buffer | 92.3s | 67% |
| locus_em | 31.4s | 23% |
| fragment_router_scan | 6.1s | 4% |
| build_loci | 3.9s | 3% |
| eb_gdna_priors | 2.4s | 2% |
| fragment_scorer | 1.5s | 1% |

### Why Only −1.7 GB Instead of Projected −5.7 GB?

The original plan projected −5.7 GB (−2 GB buffer + −3.7 GB ScoredFragments).
The actual savings were much smaller because:

1. **Buffer was already spilled to disk.** With only 2 GiB in-memory
   budget and 24M fragments × ~98 bytes/frag = ~2.3 GB, the buffer
   spills early.  Phase 5A freed metadata and file handles, not 2 GB
   of live data.

2. **Phase 5B frees ScoredFragments AFTER peak.** The peak RSS occurs
   during scanner accumulator → numpy conversion (FragmentRouter.scan),
   when both `array.array` accumulators and the new numpy arrays coexist.
   By the time `del em_data` runs (after EM), peak has already been set.

3. **builder (FragmentRouter) not freed.** The `builder` object holds
   ~1.3 GB of `array.array` accumulators that remain reachable as a
   local variable throughout the EM phase.  This was not addressed by
   5A or 5B.

**Key insight:** The dominant memory peak is during `FragmentRouter.scan()
→ _to_np()`, where both array.array and numpy copies coexist.  Phases
5A and 5B target post-peak cleanup, not peak reduction.  The high-
priority items are those that reduce the *peak itself*.


## 2. Critical Bug Fixed: C++ Type Mismatch

### Problem

`extract_locus_sub_problem` (em_solver.cpp:1245) declared `g_log_liks`,
`g_coverage_wts`, and `g_gdna_log_liks` as `const double*`, but
`batch_locus_em` (em_solver.cpp:1697) extracted `const float*` from
the `f32_1d` nanobind arrays and passed them directly.  This was a
**compile error** (not just UB) — the code could not build.

The mismatch was introduced when Phase 5E type narrowing converted
ScoredFragments storage from float64 to float32 in scan.py, but the
`extract_locus_sub_problem` signature was not updated to match.

### Fix Applied

Changed `extract_locus_sub_problem` parameter types from `const double*`
to `const float*` for the three affected arrays.  The function reads
individual float values into local `double` variables (implicit
widening), so EM numerical precision is unchanged.

**Build succeeds, all 766 tests pass.**


## 3. New Findings: Memory Waste Still Present

### 3a. builder (FragmentRouter) Not Freed — ~1.3 GB Leak

After `em_data = builder.scan(buffer, log_every)` at pipeline.py:490,
the `builder` variable remains in scope through the entire EM phase.
It holds 15 `array.array` accumulators that were the source data for
ScoredFragments.  These are now dead data but Python cannot collect
them because `builder` is a live local variable.

**Fix:** Add `del builder, ctx` immediately after scan, before
`gc.collect()`.  Same fix needed in profiler.py.

**Effort:** 2 lines.  **Savings:** ~1.3 GB at 24M frags.

### 3b. Triple-Copy Data Path in _scan_native — 3× Memory at Peak

The C++ → ScoredFragments path copies each data element three times:

```
C++ scan_chunk → std::vector<double>  [COPY 0: C++ native, float64]
  → .astype(np.float32).tobytes()     [COPY 1: numpy → bytes, float32]
  → array.array.frombytes()           [COPY 2: bytes → array.array, float32]
  → _to_np: frombuffer().copy()       [COPY 3: array.array → numpy, float32]
```

At the moment _to_np runs, both the array.array (copy 2) and the new
numpy array (copy 3) exist simultaneously = **2× ScoredFragments size
at peak**.  The C++ native copy (copy 0) is per-chunk and transient.

Additionally, C++ scoring.cpp produces `std::vector<double>` (float64)
for log_liks, coverage_weights, and gdna_log_liks.  These are then
downcast to float32 in Python via `.astype(np.float32)`.  Producing
float32 directly in C++ would save ~40% per-chunk memory and eliminate
one conversion step.

### 3c. ScoredFragments Docstring Drift

The ScoredFragments class docstring (estimator.py:92) says `float64`
for log_liks and coverage_weights, but actual arrays are `float32`
after scan.py narrowing.  Should be updated to match.

### 3d. TranscriptIndex Overhead — ~250-300 MB

| Structure | Type | Est. Size | CSR Alternative |
|---|---|---:|---|
| `_iv_t_set` | `list[frozenset[int]]` | ~50 MB | Flat int32 array + offsets |
| `sj_map` | `dict[tuple, frozenset[int]]` | ~150 MB | Sorted numpy arrays + binary search |
| `_t_exon_intervals` | `dict[int, np.ndarray]` | ~10 MB | Single flat int32 + offsets |
| Other (DataFrames, numpy) | mixed | ~50 MB | — |

Python frozenset has 56 bytes overhead per set.  For ~700K intervals
and ~1M splice junctions, this is substantial.  CSR-packing would
reduce `_iv_t_set` + `sj_map` from ~200 MB to ~30-40 MB.

**Priority:** Moderate.  The index is loaded once and persists for the
entire run.  Savings are fixed, not proportional to library size.  But
on a 16 GB machine, 200 MB matters.


## 4. Revised Priority Order

The original plan ordered phases as 5A → 5B → 5C → 5E → 5D → 5F.
Based on profiling, the revised order targets *peak reduction* first:

### Priority 1: Free builder After Scan (Trivial, ~1.3 GB)

Add `del builder, ctx` + `gc.collect()` in pipeline.py after scan.
This is the single easiest win remaining.

**Location:** pipeline.py:493 (after `buffer.cleanup()`), profiler.py:554

```python
em_data = builder.scan(buffer, log_every)
# -- Free scanner accumulators and buffer --
del builder, ctx          # NEW: free ~1.3 GB of array.array accumulators
buffer.cleanup()
buffer._chunks.clear()
buffer._memory_bytes = 0
gc.collect()
```

### Priority 2: Two-Pass Scanner — Phase 5D (High Impact, ~200 lines)

**This is the highest-impact optimization for scaling.**

The current scanner accumulates into `array.array` then copies to
numpy, doubling memory at peak.  A two-pass approach eliminates this:

1. **Count pass:** Iterate chunks, count per-unit candidates.
   Determine `n_units` and `n_candidates`.
2. **Pre-allocate:** Create final numpy arrays at exact size.
3. **Fill pass:** Re-iterate chunks, write directly into numpy
   at correct offsets.

**Peak memory drops from 2× ScoredFragments to 1× ScoredFragments.**

At 24M fragments with ~10M EM units × ~5 candidates:
- Current: ~2.6 GB (array.array) + ~2.6 GB (numpy) = **~5.2 GB peak**
- After: ~2.6 GB (numpy only) = **~2.6 GB peak**
- **Savings: ~2.6 GB at 24M frags, ~10+ GB at 100M frags**

The two-pass approach also eliminates the `array.array` intermediate
entirely, simplifying the data flow from:
```
C++ double → astype float32 → tobytes → frombytes → frombuffer → copy
```
to:
```
C++ float32 → direct write into pre-allocated numpy
```

**Sub-steps for Phase 5D:**

1. Add a lightweight counting scan pass to FragmentRouter that
   iterates buffer chunks and counts candidates per unit without
   storing any candidate data.
2. Pre-allocate all ScoredFragments numpy arrays at exact final size.
3. Modify the fill pass to write directly into the pre-allocated
   arrays at computed offsets.
4. Remove all `array.array` accumulators from FragmentRouter.
5. (Optional) Modify scoring.cpp `scan_chunk` to produce float32
   natively for log_liks, coverage_weights, gdna_log_liks — saves
   per-chunk memory and eliminates the `.astype(np.float32)` step.

### Priority 3: C++ Native float32 Output — Phase 5E-native (Moderate)

scoring.cpp `scan_chunk` produces `std::vector<double>` for log_liks,
coverage_weights, and gdna_log_liks, returning them as float64 numpy
arrays.  Python then immediately downcasts to float32 via
`.astype(np.float32)`.

Changing scoring.cpp to produce `std::vector<float>` directly would:
- Halve per-chunk C++ memory for these three arrays
- Eliminate the `.astype()` conversion allocation
- Naturally feed into Phase 5D's direct-write path

This pairs with Phase 5D and should be done concurrently.

### Priority 4: `--temp-dir` CLI — Phase 5C (Usability)

Wire the existing `BamScanConfig.spill_dir` parameter to the CLI.
Low effort, needed for production deployments on machines with
constrained `/tmp`.

### Priority 5: TranscriptIndex CSR Packing (Moderate, ~150 MB)

Replace `_iv_t_set` (list of frozenset) and `sj_map` (dict of
frozenset) with CSR-packed numpy arrays.  Fixed ~150-200 MB savings
regardless of library size.

### Priority 6: Memory-Mapped ScoredFragments — Phase 5F (Future)

Only needed if scaling to 200M+ fragments on 16 GB machines, or if
ScoredFragments exceeds physical RAM even after Phase 5D.


## 5. Revised Memory Projections

### Assumptions

- Index baseline: ~4 GB (TranscriptIndex + Python interpreter + htslib)
- EM units: 25% of fragments, avg 5 candidates/unit
- ScoredFragments: ~26 bytes/candidate + ~36 bytes/unit
- Buffer spills at 2 GiB (contributes 0 to peak RSS after spill)

### After Priority 1 (del builder) — Trivial

| Frags (M) | Current Peak | After | Savings |
|---:|---:|---:|---:|
| 24 | 18.5 GB | ~17.2 GB | ~1.3 GB |
| 50 | ~25 GB | ~23.5 GB | ~1.5 GB |
| 100 | ~36 GB | ~34 GB | ~2 GB |

### After Priority 1 + 2 (del builder + two-pass scanner)

| Frags (M) | Current Peak | After | Savings | Feasible On |
|---:|---:|---:|---:|---|
| 24 | 18.5 GB | ~13.5 GB | ~5 GB | **16 GB** (tight) |
| 50 | ~25 GB | ~16 GB | ~9 GB | **32 GB** |
| 100 | ~36 GB | ~22 GB | ~14 GB | 32 GB |
| 200 | ~53 GB | ~30 GB | ~23 GB | 64 GB |

### After Priority 1 + 2 + 3 + 5 (full plan excl. mmap)

| Frags (M) | Current Peak | After | Savings | Feasible On |
|---:|---:|---:|---:|---|
| 24 | 18.5 GB | ~12 GB | ~6.5 GB | **16 GB** |
| 50 | ~25 GB | ~14 GB | ~11 GB | **16 GB** |
| 100 | ~36 GB | ~20 GB | ~16 GB | **32 GB** |
| 200 | ~53 GB | ~28 GB | ~25 GB | 32 GB |


## 6. Implementation Checklist

| # | Change | Impact | Effort | Status |
|---|---|---|---|---|
| **P1** | `del builder, ctx` after scan | −1.3 GB | 2 lines | **Ready to implement** |
| **P2** | Two-pass scanner (Phase 5D) | Halve peak | ~200 lines Python | Not started |
| **P3** | C++ float32 output in scoring.cpp | −30% per-chunk | ~50 lines C++ | Not started |
| **P4** | `--temp-dir` CLI (Phase 5C) | Usability | ~20 lines | Not started |
| **P5** | TranscriptIndex CSR packing | −150 MB fixed | ~150 lines | Not started |
| **P6** | Memory-mapped ScoredFragments (Phase 5F) | Extreme scale | ~500 lines | Deferred |
| — | Phase 5A (free buffer after scan) | −metadata | 4 lines | **Done** ✓ |
| — | Phase 5B (free ScoredFragments after EM) | Post-peak | 2 lines | **Done** ✓ |
| — | Phase 5E (type narrowing in scan.py) | −25% CSR | — | **Done** ✓ (partially) |
| — | C++ extract_locus_sub_problem float* fix | Correctness | 3 lines | **Done** ✓ |


## 7. Correctness Issues to Track

1. **ScoredFragments docstring** (estimator.py:92-132): Says float64
   for log_liks and coverage_weights; actual dtype is float32.  Update
   when touching estimator.py.

2. **test_index_integrity.py**: 16 tests have pre-existing import
   errors (`from tests.conftest import MINI_GTF` fails with
   `ModuleNotFoundError: No module named 'tests'`).  Unrelated to
   memory work but should be fixed.

3. **pyfallback.py type alignment**: Fixed in this session — was using
   np.float64/np.int64 in `_to_np()` but array.array uses float32/int32
   typecodes.  Changed to np.float32/np.int32.  All tests pass.
