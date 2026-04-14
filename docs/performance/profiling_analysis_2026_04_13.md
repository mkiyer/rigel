# Rigel Performance Analysis and Optimization Plan

**Date**: 2026-04-13
**Platform**: Linux 4.18.0 (RHEL 8), x86-64, 8 cores, 60 GB RAM
**Python**: 3.12.13 (conda-forge), NumPy 2.x
**BAM**: VCaP simulation, STAR-aligned, 19.4M fragments (1.8 GB, queryname-sorted)
**Index**: Human genome, 457,513 transcripts (254K coding + 203K synthetic nRNA), 63,472 genes
**Profiling tools**: Stage profiler (`profiler.py --stages`), py-spy 0.4.1 (99 Hz, `--native`, debug symbols)

---

## 1. Executive Summary

A performance regression was identified and fixed: an O(N×M) loop in `calibration.py:_build_gdna_fl_model()` caused a ~45-minute hang (4,600× slowdown). After the fix, the full pipeline completes in **~121 seconds** with an 8.5 GB peak RSS.

Profiling with py-spy (native C++ symbol resolution) reveals the top five optimization opportunities ordered by estimated impact:

| # | Optimization | Current Cost | Est. Savings | Effort |
|---|-------------|-------------|--------------|--------|
| 1 | Reduce GC pressure in partition | 20.0s (22%) | 15–18s | Medium |
| 2 | Pre-allocate StreamingScorer vectors | 5.6s (6%) | 3–5s | Low |
| 3 | Reduce nanobind list↔ndarray conversion in BAM scan | 6.4s (7%) | 3–5s | Medium |
| 4 | Eliminate scatter_candidates memmove | 2.2s (2.4%) | 1–2s | Low |
| 5 | Streaming CSR to cut peak memory | — | 3 GB RSS | High |

**Combined estimated impact**: 22–30s saved (18–25% speedup), 3 GB RSS reduction.

---

## 2. Stage Timing Baseline

Two independent profiling runs (with and without py-spy) agree closely:

| Stage | Profiler (s) | % | py-spy (s) | Notes |
|-------|-------------|---|------------|-------|
| `scan_and_buffer` | 81.0 | 67.0% | ~11* | C++ BAM scan; py-spy undercounts (multithreaded) |
| `locus_em` | 24.0 | 19.8% | 18.2 | C++ EM solver, 8 OpenMP threads |
| `partition` | — | (incl.) | 15.8 | Scatter + GC; profiler lumps into quant_from_buffer |
| `fragment_router_scan` | 11.2 | 9.2% | 13.8 | Scoring + CSR construction |
| `build_loci` | 1.6 | 1.3% | 2.5 | Connected components |
| `calibration` | 0.6 | 0.5% | — | Post-fix; was ~45 min before |
| **Total wall** | **121.0** | | **~90** | py-spy measures fewer wall seconds (thread merging) |

*py-spy collapses multithreaded samples; `scan_and_buffer` wall time is 81s but most is spent in C++ worker threads that py-spy samples less effectively.

## 3. Memory Profile

| Checkpoint | RSS (MB) | Δ |
|-----------|----------|---|
| After index load | 2,835 | baseline |
| After BAM scan | 5,474 | +2,639 (fragment buffer, 2.7 GB) |
| After scoring/routing | 8,622 | +3,148 (CSR arrays) |
| After buffer release | 5,563 | −3,059 freed |
| After EM | 5,282 | −281 |
| **Peak** | **8,622** | buffer + CSR coexist |

The peak RSS is dominated by the simultaneous presence of the fragment buffer (~2.7 GB) and the scored CSR arrays (~3.1 GB) before the buffer is freed.

---

## 4. py-spy Native Profile Findings

### 4.1 Self-Time (Leaf Frame) Distribution

| Category | Self Time | % of Total | Key Functions |
|---------|-----------|------------|---------------|
| **Python GC** | 20.0s | 22.2% | `gc_collect_main` |
| **NumPy operations** | 17.4s | 19.3% | Array creation, argsort, type discovery |
| **libc memory ops** | 7.5s | 8.3% | `memmove`, `memset`, `memcmp` |
| **Scoring C++** | ~5.6s | 6.2% | `compute_fragment_weight`, `push_back`, `genomic_to_tx_pos` |
| **EM C++** | ~5.5s | 6.1% | `em_step_kernel_range`, `KahanAccumulator`, `fast_exp` |
| **nanobind conversion** | ~2.8s | 3.1% | `list_caster::from_cpp` (vector→Python list) |
| **Index load** | ~2.0s | 2.2% | Feather I/O, string decode |
| **Other** | ~30s | 33% | Various (Python dispatch, htslib I/O, etc.) |

### 4.2 GC Attribution by Call Site

The garbage collector is triggered 13× during partition scatter plus additional times during EM cleanup:

| Call Site | GC Time | % |
|----------|---------|---|
| `partition.py:97` (unit scatter loop) | 5.5s | 6.1% |
| `partition.py:70` (candidate scatter loop) | 4.7s | 5.1% |
| `index.py:947` (index load) | 3.8s | 4.2% |
| `pipeline.py:616` (EM cleanup) | 1.0s | 1.1% |
| `index.py:1050` (index load) | 1.0s | 1.1% |
| `partition.py:74` (gc.collect after offsets free) | 0.9s | 1.0% |
| Other sites | ~3.1s | 3.5% |

**Root cause**: `partition_and_free()` explicitly calls `gc.collect()` after freeing each global array to bound peak memory. With ~5 GB of Python objects on the heap (buffer + CSR + per-locus results), each GC sweep takes ~1 second. With 13 arrays × 1s = ~13s of pure GC.

### 4.3 Scoring Vector Reallocation

The C++ `StreamingScorer::score_chunk_impl()` function uses 18 `std::vector`s that grow via `push_back()` with no `reserve()`. py-spy shows:

| Function | Self Time | Notes |
|---------|-----------|-------|
| `stl_vector.h:1081` (push_back/construct) | 1.75s | int and double vectors |
| `__memmove_avx` (within scoring) | 2.83s | Vector reallocation copies |
| `compute_fragment_weight` (scoring.cpp:125) | 1.24s | Fragment weight computation |

The memmove within scoring (2.83s) is almost entirely caused by vector reallocation during push_back. With 16M EM units × ~2.5 candidates = ~40M push_backs across 6 per-candidate vectors, the ~15 reallocation cycles per vector cause ~240M unnecessary element copies.

### 4.4 Nanobind Conversion in BAM Scan

The BAM scanner returns per-fragment results as Python dicts containing lists, which nanobind converts to numpy arrays. py-spy shows 6.4s spent in `PyArray_FromAny_int` / `DiscoverDTypeAndShape_Recursive` / `AssignFromCache_Recursive` within `scan_and_buffer`.

This is the cost of converting C++ `std::vector<int>` → Python `list` → `np.array()`:
- `nanobind::list_caster::from_cpp`: 2.76s creating Python lists from C++ vectors
- `PyArray_FromAny_int`: 6.39s converting those lists to numpy arrays

### 4.5 EM Solver Hot Path

The EM solver is well-optimized. The actual computation time is spent where expected:

| Function | Self Time | Notes |
|---------|-----------|-------|
| `_mm256_fmadd_pd` (FMA intrinsic) | 1.26s | Vectorized E-step computation |
| `em_step_kernel_range` (various lines) | 2.5s | E-step inner loop |
| `KahanAccumulator::add` | 0.70s | Compensated summation |
| `fast_exp_scalar` | 0.25s | exp() approximation |
| `__memset_avx2` | 1.04s | Zeroing accumulator buffers |

The memset (1.04s) is for zeroing the per-locus accumulator arrays before each E-step. The scatter functions (scatter_candidates, scatter_units) contribute 2.2s of memmove within the EM-attributed samples.

---

## 5. Optimization Proposals (Priority Order)

### P1: Reduce GC Pressure in Partition — Est. 15–18s savings

**Current**: `partition_and_free()` calls `gc.collect()` after freeing each of 13 global arrays to prevent the buffer + CSR peak from persisting. Each GC sweep costs ~1s on the 5+ GB heap.

**Proposal**: Replace per-array `gc.collect()` with a single call after all arrays are freed, or disable GC during the scatter loop entirely and call once at the end:

```python
gc.disable()
try:
    # ... scatter all arrays, set originals to None ...
    pass
finally:
    gc.enable()
    gc.collect()
```

**Risk**: Peak RSS may spike briefly (~1–2 GB) since freed arrays won't be reclaimed until the final gc.collect(). This is acceptable given the 8.5 GB current peak — the transient spike is bounded by the already-freed buffer.

**Alternative**: Batch the gc.collect() calls — do one after the candidate phase (6 arrays) and one after the unit phase (8 arrays), reducing from 13 to 2 GC sweeps (~2s instead of ~13s).

### P2: Pre-allocate StreamingScorer Vectors — Est. 3–5s savings

**Current**: All 18 `std::vector` members in `StreamingScorer` ([scoring.cpp](../../src/rigel/native/scoring.cpp)) are initialized with zero capacity and grow via `push_back()`. The chunk-scoring loop processes ~1M fragments per chunk, generating ~2.5M candidate-level entries and ~1M unit-level entries.

**Proposal**: Add `reserve()` calls in the `StreamingScorer` constructor or `reset()` method based on the known chunk size:

```cpp
// In StreamingScorer constructor or reset_for_chunk(n_frags):
size_t est_candidates = n_frags * 3;  // ~2.5× avg + headroom
size_t est_units = n_frags;
v_ti_->reserve(est_candidates);
v_ll_->reserve(est_candidates);
v_cw_->reserve(est_candidates);
// ... etc for all per-candidate vectors
v_offsets_->reserve(est_units);
v_fid_->reserve(est_units);
// ... etc for all per-unit vectors
```

**Impact**: Eliminates ~240M unnecessary element copies from reallocation. The 2.83s memmove within scoring and 1.75s push_back overhead should drop to near-zero reallocation cost.

### P3: Return numpy Arrays Directly from C++ BAM Scanner — Est. 3–5s savings

**Current**: The C++ BAM scanner accumulates results in `std::vector`s, which nanobind converts to Python `list` objects (2.76s), which are then converted to numpy arrays via `np.asarray()` (6.4s). This is a double conversion: C++ vector → Python list → numpy array.

**Proposal**: Return `nb::ndarray<numpy>` directly from the C++ scanner instead of `std::vector`. nanobind's ndarray support provides zero-copy return of C++ allocations as numpy arrays. This eliminates both the `list_caster` conversion and the `PyArray_FromAny` copy.

The BAM scanner's `Accumulator` class already maintains contiguous buffers. Wrapping them in `nb::ndarray` capsules on return avoids the Python list intermediary entirely.

**Risk**: Requires changes to the nanobind bindings in `bam_scanner.cpp`. The Accumulator already stores data contiguously, so the change is primarily in the return interface.

### P4: Batch scatter_candidates to Reduce memmove — Est. 1–2s savings

**Current**: The `scatter_candidates_impl` functions in `em_solver.cpp` are called once per array type (6 calls for candidate arrays). Each call iterates all loci and copies data from the global CSR to per-locus arrays. The 6 separate passes over the same index structure cause 2.2s of memmove.

**Proposal**: Fuse the 6 per-type scatter calls into a single pass that scatters all candidate arrays simultaneously. This halves the index iteration overhead and improves cache locality:

```cpp
// Instead of 6 calls to scatter_candidates_impl<T>(...):
scatter_all_candidates(offsets, locus_units, offsets_list,
                       t_indices, log_liks, cov_weights, tx_starts, tx_ends, count_cols,
                       n_loci);
```

### P5: Streaming CSR Construction to Reduce Peak RSS — Est. 3 GB savings

**Current**: The fragment buffer (~2.7 GB) must coexist with the scored CSR (~3.1 GB) during the `fragment_router_scan` phase, creating a peak of 8.6 GB. The buffer is freed immediately after scoring completes.

**Proposal**: Process and free buffer chunks one at a time during scoring, spilling scored CSR results to temporary files (Arrow IPC with LZ4). Reload for the partition phase. This caps peak at ~max(buffer_chunk, CSR_chunk) + baseline ≈ 5.5 GB.

**Risk**: I/O overhead from spill/reload. With LZ4, the CSR (~3.1 GB) compresses well, and sequential write/read is fast on local storage.

---

## 6. Additional Optimization Opportunities (Lower Priority)

### P6: EM Solver — Increase MAX_K_STACK from 512 to 1024

The EM solver uses a stack-allocated buffer (`MAX_K_STACK = 512`) for the row-local accumulator. With nRNA transcripts doubling component counts, loci with >512 components hit the heap fallback. Increasing to 1024 covers the p99 case for human genome data.

### P7: EM Solver — memset Overhead in E-step

The E-step zeros the accumulator buffer (1.04s for `__memset_avx2`). For mega-loci (98K components), this is unavoidable. For small loci (<100 components), the zero-fill could be skipped by tracking dirty ranges.

### P8: Deterministic-Unambiguous Fast Path Recovery

Only 3.5% of scoring units (680K / 19.4M) take the deterministic fast path. The expected rate for mRNA-only indexes is 60–80%. The 203K synthetic nRNA transcripts double the candidate set per fragment, collapsing the fast path. Coalescing overlapping nRNA spans or increasing the merge tolerance (currently 20bp) in index construction could recover the fast path.

### P9: int32→int64 Upcast in buffer.py

`to_scoring_arrays()` upcasts `t_offsets` from int32→int64 (4 MB → 8 MB per chunk × ~16 chunks = ~64 MB extra copies). Accept int32 offsets in the C++ scoring kernel.

---

## 7. Locus EM Hotspot Analysis

| Metric | Value |
|--------|-------|
| Total loci | 29,458 |
| Mega-loci (>50K transcripts) | 1 |
| Total EM wall time (8 threads) | 24.0s |
| Top 10 loci | 17.7s (24% of EM) |

### Mega-Locus Detail (98,608 transcripts, 4.06M units, 140K ECs)

| Sub-stage | Time | % |
|-----------|------|---|
| SQUAREM iterations (59 iters) | 9.69s | 69.6% |
| Build equivalence classes | 1.43s | 10.3% |
| Assign posteriors | 1.06s | 7.6% |
| Extract fragment data | 0.95s | 6.8% |
| Bias correction | 0.69s | 4.9% |

The mega-locus is dominated by the SQUAREM solver (69.6%), which is well-optimized with AVX2 FMA intrinsics and Kahan summation. The main opportunities are reducing the locus size through nRNA index optimization (P8) or improving EC construction cache efficiency.

---

## 8. Completed Fix: Calibration O(N×M) Regression

**Commit**: `26001ab` introduced a per-region Python for-loop in `_build_gdna_fl_model()` that iterated 542K regions, each scanning 13.3M fragment-length observations. The O(N×M) complexity caused ~45 minutes of wall time.

**Fix**: Replaced with three vectorized NumPy broadcasts using fancy indexing. Calibration dropped from ~45 min to 0.58s (4,600× speedup). Validated by 1,013/1,013 tests passing.

---

## 9. Profiling Setup Reference

### Stage Profiler

```bash
conda activate rigel
python scripts/profiling/profiler.py \
  --bam <input.bam> \
  --index <rigel_index/> \
  --outdir <output_dir> \
  --threads 8 --stages
```

### py-spy with Native C++ Symbols

To get C++ function names in py-spy output, the .so files must have `-g` debug info and must *not* be linked with `-Wl,-s` (nanobind default). Temporarily modify CMakeLists.txt:

1. Add `-g` to all `target_compile_options` calls
2. Disable LTO (`set(ipo_supported FALSE)`)
3. Rebuild: `pip install --no-build-isolation -e .`
4. Manually re-link without strip:
   ```bash
   cd build/cp312-abi3-linux_x86_64
   for mod in _em_impl _scoring_impl _bam_impl; do
       cmd=$(cat CMakeFiles/${mod}.dir/link.txt | sed 's/-Wl,-s //g')
       eval "$cmd" -o ${mod}_debug.abi3.so
       cp ${mod}_debug.abi3.so $(python -c "import rigel.${mod#_}; print(rigel.${mod#_}.__file__)" 2>/dev/null || echo /path/to/site-packages/rigel/${mod}.abi3.so)
   done
   ```
5. Run py-spy:
   ```bash
   py-spy record --native --rate 99 --format raw \
     -o pyspy_raw.txt -- python driver.py ...
   ```

### Analyzing py-spy Raw Output

The raw collapsed stack format (`func1;func2;...;leafN count`) can be analyzed with custom Python scripts or converted to flamegraphs with [flamegraph.pl](https://github.com/brendangregg/FlameGraph).

---

## 10. Summary

| Metric | Before Fix | After Fix | After All P1–P5 (est.) |
|--------|-----------|-----------|------------------------|
| Wall time | >45 min | 121s | ~90–95s |
| Peak RSS | ~8.5 GB | 8.6 GB | ~5.5 GB |
| Throughput | <6K frags/s | 345K frags/s | ~430K frags/s |
| GC overhead | ~20s (22%) | 20s (22%) | ~2s (2%) |
| Calibration | ~45 min | 0.6s | 0.6s |
