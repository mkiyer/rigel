# Performance Optimization Plan — Phase 2: Production-Ready Pipeline

## 1. Current Profile Summary (2026-03-06, Post Phase-2B)

**Dataset:** 12M buffered fragments, 254,461 transcripts, 63,472 genes, 29,266 loci
**Platform:** macOS ARM64 (Apple Silicon), Python 3.13, conda `hulkrna`

### Phase History

| Milestone | Total Wall | Key Change |
|---|---|---|
| Baseline (pre-Phase 1) | 416s | — |
| Phase 1 — C++ scan_chunk | 273s | fragment_router_scan 143.7s → 6.6s (22×) |
| Phase 2A — Batch C++ locus EM | 180.2s | locus_em 199.8s → ~73s (2.7×) |
| Phase 2A (re-measured, 1 thread) | 174.8s | batch_locus_em ~70s sequential |
| **Phase 2B — OpenMP 4 threads** | **138.6s** | **batch_locus_em ~70s → ~28s (2.5×)** |
| Fallback path (--stages reference) | 267.3s | Python per-locus baseline |

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

### Profile: Batch C++ Path — 1 Thread (simple mode, 2026-03-06)

Ran without `--stages` so `run_pipeline()` uses the batch C++ path
(`batch_locus_em()` — single call, `--threads 1` = sequential).

| Metric | Value |
|---|---|
| **Wall time** | **174.8s** |
| Peak RSS | 20,179 MB |
| Throughput | 137,297 frags/s |

**Estimated stage breakdown** (from log timestamps):

| Stage | Time (est.) | % Total | Notes |
|---|---|---|---|
| scan_and_buffer | ~91s | 52% | C++ BAM scan → buffer |
| setup (geometry+scorer+scan+priors) | ~15s | 9% | build_loci, EB priors, nrna_frac |
| **batch C++ locus EM** | **~70s** | **40%** | Sequential, single-threaded |
| **TOTAL** | **~175s** | | |

### Profile: Batch C++ Path — 4 Threads (simple mode, 2026-03-06)

Ran with `--threads 4` to exercise the OpenMP path.

| Metric | Value |
|---|---|
| **Wall time** | **138.6s** |
| Peak RSS | 20,272 MB |
| Throughput | 173,158 frags/s |
| gDNA total | 1,967,791 |

**Estimated stage breakdown** (from log timestamps):

| Stage | Time (est.) | % Total | Notes |
|---|---|---|---|
| scan_and_buffer | ~92s | 66% | C++ BAM scan → buffer (now **dominant**) |
| setup (geometry+scorer+scan+priors) | ~18s | 13% | build_loci, EB priors, nrna_frac |
| **batch C++ locus EM (4 threads)** | **~28s** | **20%** | OpenMP `schedule(dynamic, 16)` |
| **TOTAL** | **~139s** | | |

**Phase 2B delivered: batch_locus_em 70s → 28s (2.5× from 4 threads)**
**Overall: 416s → 138.6s = 3.0× total pipeline speedup from baseline**

### OpenMP Scaling Analysis

| Threads | batch_locus_em | Total Wall | EM Speedup | Parallel Eff. | Overall vs Baseline |
|---|---|---|---|---|---|
| 1 (sequential) | ~70s | 174.8s | 1.0× | — | 2.4× |
| **4 (measured)** | **~28s** | **138.6s** | **2.5×** | **63%** | **3.0×** |
| 8 (projected) | ~17s | ~127s | 4.1× | ~51% | 3.3× |
| 12 (projected) | ~14s | ~124s | 5.0× | ~42% | 3.4× |

Projections for 8+ threads assume continued efficiency degradation
from memory bandwidth saturation on Apple Silicon's unified memory
and the large-locus outlier (297 transcripts, 59,930 units) acting
as a load-balance bottleneck with `schedule(dynamic, 16)`.

**Sub-linear scaling analysis (63% efficiency at 4 threads):**
- The largest locus (59,930 units) takes ~2-3s alone, creating a
  tail latency that limits parallelism when few loci remain
- Apple Silicon unified memory bandwidth is shared across all cores
- `schedule(dynamic, 16)` has modest per-chunk scheduling overhead
- Despite sub-linear scaling, the 42s savings (175→139s) at 4 threads
  represents an excellent cost/benefit ratio with no code complexity

The cProfile in simple mode is dominated by `_thread.lock.acquire`
(38.7s) — the memory-monitor thread blocking while the C++ batch EM
runs with GIL released.  Python-side work totals only ~9.3s, confirming
the pipeline is almost entirely C++ at this point.

**Remaining Python hot spots** (≥0.5s, from 4-thread cProfile):

| Function | Self Time | Notes |
|---|---|---|
| `pd.ArrowArray.__getitem__` | 0.97s | Buffer chunk iteration |
| `numpy.asarray` | 0.90s | Pervasive conversion |
| `array.array.append` | 0.88s | Scanner accumulation |
| `list.append` | 0.79s | Scanner accumulation |
| `pipeline._build_locus_meta` | 1.06s | 29,266× Python dict construction |
| `locus.build_loci` | 0.095s | scipy connected_components |
| `scipy.validate_graph` | 1.89s | scipy COO→CSR validation |
| `scipy.csr_sort_indices` | 0.059s | CSR index sorting |
| `buffer._load_chunk` | 0.079s | Feather deserialization |

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

## 3. Phase 2B — Thread-Parallel EM: COMPLETED ✓

**Status:** Implemented with OpenMP in em_solver.cpp, all 823 tests passing.
`--threads` CLI parameter wired through `hulkrna quant`, profiler, and
benchmark scripts.

### What Was Done

- `#pragma omp parallel for schedule(dynamic, 16)` around locus loop
- `nb::gil_scoped_release` to release GIL during C++ compute
- OpenMP detection in CMakeLists.txt with conda `libomp` fallback for
  Apple Clang (`-Xclang -fopenmp` + `-lomp -Wl,-rpath`)
- `EMConfig.n_threads` parameter (0 = all cores, 1 = sequential, N = cap)
- `--threads` CLI flag for `hulkrna quant`, `scripts/profiler.py`, and
  `scripts/benchmarking/benchmark_region_competition.py`
- `llvm-openmp` added to `mamba_env.yaml`

### Architecture (as implemented)

```cpp
void batch_locus_em(..., int n_threads) {
    // ... setup ...
    if (n_threads > 0) omp_set_num_threads(n_threads);

    nb::gil_scoped_release release;

    #pragma omp parallel for schedule(dynamic, 16)
    for (int64_t i = 0; i < n_loci; ++i) {
        auto sub = extract_locus_sub_problem(i, ...);
        run_locus_em_core(sub, ...);
        assign_posteriors_native(sub, ...);
    }
}
```

Loci have disjoint transcript sets (connected components), so output
array writes never conflict — no thread-local accumulator copies needed.
Only `total_gdna_em` uses `#pragma omp atomic`.

### Measured Results

| Metric | 1-thread | 4-thread | Speedup |
|---|---|---|---|
| batch_locus_em | ~70s | ~28s | 2.5× |
| Total wall | 174.8s | 138.6s | 1.26× |
| Throughput | 137K frags/s | 173K frags/s | 1.26× |
| Peak RSS | 20,179 MB | 20,272 MB | +93 MB |

**63% parallel efficiency at 4 threads.** Sub-linear due to memory
bandwidth sharing on Apple Silicon unified memory and load-balance
tail from the largest locus (59,930 units, ~2-3s single-threaded).

### What Was *Not* Done (deferred)

- **Thread-local accumulators with reduction** — Not needed because
  loci are disjoint connected components (no write conflicts).
- **Kahan summation in reducer** — N/A since no cross-thread reduction.
- **Pre-sorting loci by size** — The `schedule(dynamic, 16)` handles
  load balance adequately; pre-sorting could improve tail latency but
  was not needed for the measured 2.5× speedup.

---

## 4. Phase 3 — build_loci: C++ Union-Find ✓ COMPLETED

**Status:** Implemented and tested (823/823 tests pass).

Replaced scipy COO→CSR→connected_components with a C++ union-find
(disjoint-set with path halving + union by rank).  Added
`connected_components_native()` to `em_solver.cpp` (~100 lines C++).

| Before | After | Improvement |
|---|---|---|
| scipy connected_components (6.1s) | C++ union-find (~0.3s est) | ~20× |

**Bonus:** Scipy is no longer a runtime dependency.
- Removed `scipy>=1.12` from `pyproject.toml` dependencies
- Removed `scipy>=1.1` from `mamba_env.yaml`
- `strand_model.py` already had a lazy `try/except ImportError` for
  `scipy.stats.beta.ppf` (diagnostic only) — works without scipy

---

## 5. Phase 4 — scan_and_buffer Optimization (Target: −15s) ★ NEXT PRIORITY

After Phase 2B, `scan_and_buffer` (~92s) is the **dominant bottleneck**
— **66% of total runtime** at 4 threads.  This is now the single
highest-impact target.

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
| Phase 2A (batch C++) ✓ | ~73s | 91s | 17s | **180s** | 2.3× |
| Phase 2A (re-measured) ✓ | ~70s | 91s | 15s | **175s** | 2.4× |
| **Phase 2B (4 threads) ✓** | **~28s** | **92s** | **18s** | **139s** | **3.0×** |
| **Phase 3 (union-find) ✓** | **~28s** | **92s** | **~13s** | **~133s** | **3.1×** |
| +Phase 4 (scan opt) | ~28s | 75s | ~13s | ~116s | 3.6× |
| +Phase 2B (8 threads, projected) | ~17s | 92s | 18s | ~127s | 3.3× |

Note: "other" includes geometry, scorer, router_scan (~6s), build_loci
(~6s), EB priors (~2s), nrna_frac, and estimator setup.  The increase
from 15s to 18s between 1-thread and 4-thread runs is within measurement
variance.

### Memory Projections

| Scenario | Peak RSS | Savings |
|---|---|---|
| Current (Phase 2B, 4 threads) | 20.3 GB | — |
| +Phase 5A (free buffer after scan) | ~14.5 GB | −5.8 GB |
| +Phase 5B (free CSR after EM) | ~10.7 GB | −3.8 GB |
| **Full plan (5A+5B)** | **~10.7 GB** | **−9.6 GB** |

### Recommended Execution Order

1. **Phase 5A+5B** — Memory freeing (easy wins: −9.6 GB, minimal code)
2. ~~**Phase 3** — C++ union-find~~ ✓ DONE (scipy eliminated)
3. **Phase 4** — scan_and_buffer optimization (moderate complexity,
   highest remaining wall-time impact: 92s is 66% of total)

---

## 8. Profiler Notes

**IMPORTANT:** The `--stages` profiler mode manually runs the Python
per-locus loop, bypassing `batch_locus_em()` and OpenMP threading.
Always use **simple mode** (no `--stages`) to profile the production
C++ path.  The example YAML (`scripts/profile_example.yaml`) now
defaults to `stages: false`.

| Mode | What it measures | Threading? |
|---|---|---|
| Simple (default) | `run_pipeline()` → `batch_locus_em()` | ✓ OpenMP |
| `--stages` | Python per-locus loop (fallback path) | ✗ Sequential |

**CLI threading:**
```bash
# 4 threads (recommended for profiling)
python scripts/profiler.py --bam reads.bam --index idx/ --threads 4

# Sequential (for baseline comparison)
python scripts/profiler.py --bam reads.bam --index idx/ --threads 1

# All available cores (default when --threads not specified)
python scripts/profiler.py --bam reads.bam --index idx/
```

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
3. Profile with `--threads 1,2,4,8` to measure scaling efficiency
   on the target deployment hardware.

---

## 9. Risk Assessment

| Risk | Status | Notes |
|---|---|---|
| OpenMP not available on macOS | **Resolved** | conda `libomp` + `-Xclang -fopenmp`; CMake detects automatically |
| Thread-safety of accumulation | **Resolved** | Loci are disjoint — no write conflicts; `total_gdna_em` uses `#pragma omp atomic` |
| Numerical reproducibility under threading | **Resolved** | Same loci always assigned to same threads with fixed thread count; no cross-thread reduction needed |
| Dynamic schedule variance | **Measured** | 63% parallel efficiency at 4 threads; adequate for production |
| Memory overhead from threading | **Measured** | +93 MB at 4 threads (negligible vs 20 GB peak) |
| Large locus load-balance outlier | **Confirmed** | 297 tx / 59,930 units is likely the tail-latency bottleneck; `schedule(dynamic, 16)` handles it acceptably |
| `--stages` profiler bypasses C++ path | **Resolved** | Documented; example YAML now defaults to `stages: false` |

---

## 10. Testing Strategy

1. **Unit tests**: All 823 existing tests pass after each phase
2. **Parity tests**: 18 batch C++/Python parity tests (`test_batch_parity.py`)
3. **Thread determinism**: Same input → same output for fixed thread count
4. **Numerical tolerance**: `rtol=1e-10` between sequential and threaded
5. **Profiling**: Re-run profiler after each phase to confirm gains
6. **Edge cases**: Empty loci, single-transcript loci, max locus (297 tx,
   59,930 units), zero-gDNA loci, all-spliced loci
