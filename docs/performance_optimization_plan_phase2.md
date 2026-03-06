# Performance Optimization Plan — Phase 2: Production-Ready Pipeline

## 1. Current Profile Summary (2026-03-06, Post Phase-2A)

**Dataset:** 12M buffered fragments, 254,461 transcripts, 63,472 genes, 29,266 loci
**Platform:** macOS ARM64 (Apple Silicon), Python 3.13, conda `hulkrna`

### Phase History

| Milestone | Total Wall | Key Change |
|---|---|---|
| Baseline (pre-Phase 1) | 416s | — |
| Phase 1 — C++ scan_chunk | 273s | fragment_router_scan 143.7s → 6.6s (22×) |
| Phase 2A — Batch C++ locus EM | **180.2s** | locus_em 199.8s → ~73s (2.7×) |
| Fallback path (for sub-stage reference) | 270.9s | Python per-locus baseline |

### Profile: Python Fallback Path (--stages mode)

This profile used the profiler's `--stages` mode, which **manually
calls the per-locus Python loop** to instrument sub-stage timings.
It therefore measures the *fallback* path, not the batch C++ path.
It remains valuable as a detailed decomposition of where time is spent.

| Stage | Wall Time | % Total | Notes |
|---|---|---|---|
| **locus_em** | **199.8s** | **73.8%** | Dominant bottleneck |
| scan_and_buffer | 54.8s | 20.2% | C++ BAM scan → buffer |
| fragment_router_scan | 6.4s | 2.4% | Phase 1 win (22× faster) |
| build_loci | 6.1s | 2.3% | scipy connected_components |
| eb_gdna_priors | 2.2s | 0.8% | Empirical Bayes gDNA |
| other | 1.6s | 0.6% | geometry, estimator, scorer, nrna init |
| **TOTAL** | **270.9s** | **100%** | |

**Memory:** 20,105 MB peak RSS (16,241 MB delta from 3,864 MB baseline)

### Locus EM Sub-stage Decomposition (29,266 loci)

| Sub-stage | Wall Time | cProfile Self | Per-call | Batch C++ Impact |
|---|---|---|---|---|
| `build_locus_em_data` | 15.4s | 9.1s | 0.53ms | **Eliminated** (reimplemented in C++) |
| `run_locus_em` (C++ solver) | 179.4s | 104.0s | 3.6ms | **Same solver; no Python overhead** |
| `assign_posteriors` | 5.0s | 2.5s | 0.17ms | **Eliminated** (reimplemented in C++) |

**Gap analysis: 179.4s wall vs 104.0s CPU for `run_locus_em`**

The 75.4s gap is attributable to:
- Python→C++ round-trip overhead × 29,266 calls (~30-40s)
- numpy array allocation/deallocation for inputs/outputs (~15-20s)
- GIL contention with the profiler's memory-sampling thread (~15-20s)

All of this overhead is eliminated by the batch C++ path.

### cProfile Hot Spots (top 10 by self-time)

| # | Function | Self Time | Calls | Notes |
|---|---|---|---|---|
| 1 | `pyfallback.run_locus_em` | 104.0s | 29,266 | C++ solver + marshaling |
| 2 | `_thread.lock.acquire` | 84.0s | 636 | Memory monitor thread |
| 3 | `locus.build_locus_em_data` | 9.1s | 29,267 | **Eliminated by batch C++** |
| 4 | `pyfallback.assign_locus_ambiguous` | 2.5s | 29,266 | **Eliminated by batch C++** |
| 5 | `numpy.full` | 1.3s | 87,776 | 3× per locus in build_locus_em_data |
| 6 | `threading.wait` | 1.1s | 159 | Memory monitor sleep |
| 7 | `numpy.ufunc.at` | 1.0s | 149,553 | Scatter ops in build/assign |
| 8 | `scipy.csr_sort_indices` | 0.9s | 1 | build_loci graph construction |
| 9 | `numpy.asarray` | 0.9s | 58,995 | Pervasive conversion overhead |
| 10 | `array.array.append` | 0.9s | 10.3M | Scanner accumulation |

### Memory Timeline Analysis

```
Time(s)   RSS(MB)   Event
───────   ───────   ─────────────────────────────────────
  0       3,864     Index loaded
 54       9,662     scan_and_buffer complete (+5.8 GB)
 59      14,154     ScoredFragments CSR allocated (+4.5 GB)
 63      19,546     Estimator + loci allocated (+5.4 GB)
 70      20,105     Peak — stable through EM (+0.6 GB)
271      20,105     End — nothing freed
```

Key observations:
- **5.8 GB** consumed by buffer chunks (BAM scan output)
- **10 GB** consumed by ScoredFragments CSR + estimator arrays
- **Memory never decreases** — buffer and CSR persist through EM
- Freeing buffer/CSR after use could recover **~8-10 GB**

### Profile: Batch C++ Path (simple mode, 2026-03-06)

Ran without `--stages` so `run_pipeline()` uses the batch C++ path
(`batch_locus_em()` — single call, no Python per-locus loop).

| Metric | Value |
|---|---|
| **Wall time** | **180.2s** |
| Peak RSS | 20,212 MB |
| Throughput | 133,162 frags/s |
| gDNA total | 1,967,791 |

**Estimated stage breakdown** (from log timestamps + cProfile):

| Stage | Time (est.) | % Total | Notes |
|---|---|---|---|
| scan_and_buffer | ~91s | 50.5% | C++ BAM scan → buffer (now **dominant**) |
| setup (geometry+scorer+scan+priors) | ~17s | 9.4% | build_loci, EB priors, nrna_frac |
| **batch C++ locus EM** | **~73s** | **40.5%** | Sequential, single-threaded |
| **TOTAL** | **~180s** | | |

**Phase 2A delivered: locus_em 199.8s → ~73s (2.7× for EM, 1.5× overall)**

The cProfile in simple mode is dominated by `_thread.lock.acquire`
(170.6s) — the memory-monitor thread blocking while the C++ batch EM
runs with GIL released.  Python-side work totals only ~9.6s, confirming
the pipeline is almost entirely C++ at this point.

**Remaining Python hot spots** (≥0.5s, from cProfile):

| Function | Self Time | Notes |
|---|---|---|
| `pd.ArrowArray.__getitem__` | 0.96s | Buffer chunk iteration |
| `array.array.append` | 0.87s | Scanner accumulation |
| `numpy.asarray` | 0.86s | Pervasive conversion |
| `list.append` | 0.76s | Scanner accumulation |
| `pipeline._build_locus_meta` | 0.68s | 29,266× Python dict construction |
| `locus.build_loci` | 0.62s | scipy connected_components |
| `buffer._load_chunk` | 0.62s | Feather deserialization |
| `ndarray.copy` | 0.51s | Buffer finalization copies |

**Memory timeline** (simple mode):

```
Time(s)   RSS(MB)   Event
───────   ───────   ─────────────────────────────────────
  0       3,851     Index loaded
 89       7,468     scan_and_buffer starting finalization
 91      11,867     Buffer chunks finalized (+8.0 GB)
 98      16,846     ScoredFragments + estimator (+5.0 GB)
105      20,212     Peak — stable through batch EM (+3.4 GB)
180      20,212     End — nothing freed
```

---

## 2. Phase 2A — Batch C++ Locus EM: COMPLETED ✓

**Status:** Implemented in em_solver.cpp (~500 lines C++), all 823 tests
passing on both batch C++ and Python fallback paths.

**What was done:**
- `batch_locus_em()` C++ function processes all loci in one call
- `extract_locus_sub_problem()` — C++ reimplementation of `build_locus_em_data`
- `assign_posteriors_native()` — C++ reimplementation of `assign_locus_ambiguous`
- Flat-array locus definitions (CSR offsets for transcripts/units per locus)
- Python fallback extracted to `pyfallback.py` for debugging

**Dispatch logic** (`pipeline.py`):
```python
if annotations is None and not os.environ.get("HULK_FORCE_PYTHON"):
    # Batch C++ path — single call for all loci
    total_gdna_em, mrna, nrna, gdna = estimator.run_batch_locus_em(...)
else:
    # Python per-locus fallback (for annotation support / debugging)
    from .pyfallback import run_per_locus_em_loop
    total_gdna_em = run_per_locus_em_loop(...)
```

**Measured Phase 2A results** (from batch C++ profile, 2026-03-06):

| Source | Fallback Path | Batch C++ (measured) | Savings |
|---|---|---|---|
| Python for-loop + round-trips | ~40s | 0 | −40s |
| `build_locus_em_data` numpy/pandas | 15.4s | <1s (C++) | −14.4s |
| `assign_locus_ambiguous` Python | 5.0s | <1s (C++) | −4.5s |
| C++ EM solver (sequential) | ~104s | ~71s (tighter C++ loop) | −33s |
| **locus_em total** | **199.8s** | **~73s** | **−127s** |

**Measured total: 180.2s** (was projected ~177s — close to estimate).
The EM solver itself ran faster than projected (~71s vs ~104s),
likely due to better cache locality in the flat-array C++ layout
vs the per-locus Python/numpy overhead in the fallback path.

**Phase 2A delivered: locus_em 199.8s → ~73s (2.7× EM speedup, 1.5× overall)**

---

## 3. Phase 2B — Thread-Parallel EM: NEXT PRIORITY ★

**Goal:** The C++ EM solver takes ~73s of CPU time processing 29,266
independent loci sequentially.  Add thread parallelism to cut this to
~10-15s on 8+ cores.

**This is the single highest-impact optimization remaining.**

### Architecture

```cpp
void batch_locus_em(...) {
    // ... setup ...

    nb::gil_scoped_release release;  // Already present

    std::vector<ThreadLocalAccum> tl_accums(n_threads);
    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < n_loci; ++i) {
        auto& acc = tl_accums[omp_get_thread_num()];
        auto sub = extract_locus_sub_problem(i, ...);
        auto [theta, alpha] = run_locus_em_core(sub, em_config);
        assign_posteriors(sub, theta, acc);
    }
    // Reduce: sum tl_accums → final output arrays
}
```

### Key Design Decisions

1. **`schedule(dynamic, 16)`** — Locus sizes vary enormously (1 transcript
   to 297 transcripts, 1 unit to 59,930 units).  Dynamic scheduling
   distributes large/small loci for load balance.

2. **Thread-local accumulator buffers** — Each thread gets its own
   `em_counts[N_t, n_cols]`, `nrna_em_counts[N_t]`, etc.  After
   completion, reduce (sum) across threads.  No atomic operations.
   Memory overhead: `N_t × n_cols × 8 bytes × N_threads` ≈ 130 MB
   for 8 threads — negligible vs 20 GB total.

3. **OpenMP on macOS** — Apple Clang lacks OpenMP.  Options:
   a. Use `libomp` from conda-forge (preferred — already in build env)
   b. Fall back to `std::thread` pool with work-stealing queue
   c. `CMakeLists.txt` detects OpenMP; gracefully degrades to sequential

4. **GIL release** — Already present (`nb::gil_scoped_release`) in
   `batch_locus_em()`.  No GIL contention since all work is pure C++.

5. **Determinism** — For a fixed thread count and `schedule(dynamic, 16)`,
   the work partition is deterministic.  Floating-point reduction order
   may differ; use Kahan summation in reducers for reproducibility.

### Expected Speedup

| Threads | Est. locus_em | Est. Total | Speedup (from baseline) |
|---|---|---|---|
| 1 (batch C++, no threads) | ~73s | ~180s | 2.3× |
| 2 | ~40s | ~148s | 2.8× |
| 4 | ~22s | ~130s | 3.2× |
| 8 | ~12s | ~120s | 3.5× |
| 12 (M-series 10-12 cores) | ~9s | ~117s | 3.6× |

### Implementation Plan

**Files:**
- `src/hulkrna/native/em_solver.cpp` — Add OpenMP parallel region
- `CMakeLists.txt` — Detect & link OpenMP (or libomp)
- `src/hulkrna/config.py` — Add `EMConfig.n_threads` parameter
- `src/hulkrna/estimator.py` — Pass thread count to `batch_locus_em()`

**Sub-tasks:**
1. Add `find_package(OpenMP)` to CMakeLists.txt with graceful fallback
2. Define `ThreadLocalAccum` struct in em_solver.cpp
3. Wrap locus loop with `#pragma omp parallel for schedule(dynamic, 16)`
4. Post-loop: reduce thread-local accumulators into final arrays
5. Add `n_threads` parameter to `batch_locus_em()` nanobind binding
6. Wire `EMConfig.n_threads` → `batch_locus_em(n_threads=...)`
7. Test: 823 tests must pass; numerical reproducibility within `rtol=1e-10`
8. Benchmark: 1, 2, 4, 8 thread configurations

---

## 4. Phase 3 — build_loci: C++ Union-Find (Target: −5s)

`build_loci` uses scipy COO→CSR→connected_components (6.1s).  Most of
the cost is scipy sparse matrix construction:

| scipy operation | Time |
|---|---|
| `coo_tocsr` | 0.39s |
| `csr_sort_indices` | 0.87s |
| `validate_graph` → `tocsr` | 1.78s |
| Python overhead | ~3.1s |

**Replace with C++ weighted union-find:**

```cpp
std::vector<int32_t> connected_components(
    const int64_t* offsets, const int32_t* t_indices,
    int32_t n_units, int32_t nrna_base, int32_t n_transcripts);
```

- O(N_candidates × α(N)) — near-linear, no sparse matrix needed
- Expected: 6.1s → ~0.3s (20× improvement)
- Low implementation risk — straightforward algorithm

---

## 5. Phase 4 — scan_and_buffer Optimization (Target: −15s)

After Phase 2B, `scan_and_buffer` (~91s) will be the **dominant
bottleneck** — potentially 65-75% of total runtime.

### 4A. htslib Multi-threaded BAM Decompression

Check if `hts_set_threads()` is already enabled.  If not, adding 2-4
decompression threads could save ~15-20s on I/O-bound workloads.
(Current throughput: ~264K reads/s = 24M reads / 91s.)

### 4B. Eliminate Redundant Copies in Buffer Finalization

`_FinalizedChunk` construction copies arrays that may already own their
memory via nanobind capsules.  Verify ownership semantics to eliminate
unnecessary `.copy()` calls.

### 4C. Pre-allocated Chunk Memory

Return pre-allocated numpy arrays from C++ accumulator instead of
byte-buffer → numpy conversion.

Expected total: 91s → ~75-80s

---

## 6. Phase 5 — Memory Reduction (Target: −8-10 GB)

Current peak: 20.2 GB.  Memory timeline shows no releases during the
pipeline.  Three independent opportunities:

### 5A. Free Buffer Chunks After Scan (−5-6 GB)

Buffer chunks are only needed during `fragment_router_scan` to build
the CSR.  After `scan()`, the chunk data is redundant — only the
ScoredFragments CSR is needed.  Explicit `buffer.cleanup()` or
`del buffer` after scan could free ~5.8 GB.

*The profiler already calls `buffer.cleanup()` — verify it happens
before the EM phase in the main pipeline path.*

### 5B. Free Global CSR After Batch EM (−3-4 GB)

After `batch_locus_em()`, ScoredFragments arrays (offsets, t_indices,
log_liks, coverage_weights, tx_starts, tx_ends, etc.) are no longer
needed.  Explicit `del em_data` post-EM.

### 5C. Compact Fragment Buffer Columns (−2-3 GB)

| Column | Current | Could be | Savings/12M frags |
|---|---|---|---|
| `frag_id` | int64 | int32 | 48 MB |
| `t_offsets` | int64 | uint32 | varies |
| `read_length` | uint32 | uint16 | 24 MB |

Lower priority — type narrowing requires careful validation.

### 5D. Stream-and-Discard Buffer Chunks

If the C++ scanner and CSR builder can be fused (scan directly into
ScoredFragments without materializing buffer chunks), the 5.8 GB
buffer stage can be eliminated entirely.  This is a larger refactor.

Expected total: **20.1 GB → ~12-14 GB** (40% reduction)

---

## 7. Summary: Projected Timeline

### Wall Time Projections

| Scenario | locus_em | scan_buf | other | Total | vs Baseline |
|---|---|---|---|---|---|
| **Baseline (pre-opt)** | 314s | 55s | 47s | **416s** | 1.0× |
| Phase 1 (scan_chunk) | 201s | 55s | 17s | **273s** | 1.5× |
| **Phase 2A (batch C++) ✓** | **~73s** | **91s** | **17s** | **180s** | **2.3×** |
| **Phase 2B (8 threads)** | **~12s** | **91s** | **17s** | **~120s** | **3.5×** |
| +Phase 3 (union-find) | ~12s | 91s | 12s | ~115s | 3.6× |
| +Phase 4 (scan opt) | ~12s | 75s | 12s | ~99s | 4.2× |

### Memory Projections

| Scenario | Peak RSS | Savings |
|---|---|---|
| Current | 20.2 GB | — |
| +Phase 5A (free buffer after scan) | ~14.4 GB | −5.8 GB |
| +Phase 5B (free CSR after EM) | ~10.6 GB | −3.8 GB |
| **Full plan (5A+5B)** | **~10.6 GB** | **−9.6 GB** |

### Recommended Execution Order

1. **★ Phase 2B** — Thread-parallel EM (highest ROI: ~61s savings → ~120s)
2. **Phase 5A+5B** — Memory freeing (easy wins: −9.6 GB, minimal code)
3. **Phase 3** — C++ union-find (clean −5s, low risk)
4. **Phase 4** — scan_and_buffer optimization (moderate complexity)

---

## 8. Profiler Notes

The `--stages` profiler mode manually runs the Python per-locus loop,
bypassing `batch_locus_em()`.  The Phase 2A measurement (180.2s) was
obtained by running the profiler in **simple mode** (no `--stages`),
which exercises the full `run_pipeline()` → `batch_locus_em()` path.

**Future profiling recommendations:**
1. Add batch-level timing to `quant_from_buffer()`:
   ```python
   import time
   t0 = time.perf_counter()
   total_gdna_em, ... = estimator.run_batch_locus_em(...)
   logger.info(f"batch_locus_em: {time.perf_counter() - t0:.3f}s")
   ```
2. Update the `--stages` profiler to support calling
   `estimator.run_batch_locus_em()` directly, with before/after
   timers, instead of the manual per-locus loop.
3. After Phase 2B (threading), profile with `OMP_NUM_THREADS=1,2,4,8`
   to measure scaling efficiency.

---

## 9. Risk Assessment

| Risk | Mitigation |
|---|---|
| OpenMP not available on macOS | Use `libomp` from conda-forge; detect at CMake time |
| Thread-safety of accumulation | Thread-local buffers; reduce (sum) at end |
| Numerical reproducibility under threading | Kahan summation in reducer; validate within rtol=1e-10 |
| Dynamic schedule variance | Use chunk size ≥16 to amortize scheduling overhead |
| Memory overhead from thread-local accumulators | ~130 MB for 8 threads (< 1% of peak 20 GB) |
| Large locus (297 tx, 59,930 units) as load-balance outlier | dynamic scheduling handles this; pre-sort loci by size (largest first) |

---

## 10. Testing Strategy

1. **Unit tests**: All 823 existing tests pass after each phase
2. **Parity tests**: 18 batch C++/Python parity tests (`test_batch_parity.py`)
3. **Thread determinism**: Same input → same output for fixed thread count
4. **Numerical tolerance**: `rtol=1e-10` between sequential and threaded
5. **Profiling**: Re-run profiler after each phase to confirm gains
6. **Edge cases**: Empty loci, single-transcript loci, max locus (297 tx,
   59,930 units), zero-gDNA loci, all-spliced loci
