# Rigel Performance Analysis & Optimization Plan

**Date:** March 27, 2026
**System:** 2× Intel Xeon Gold 6154 (36 cores total), 186 GB RAM, RHEL 8.10
**Python:** 3.12.13 (conda-forge), GCC 14.3.0
**Index:** Human GENCODE, 457,513 transcripts, 63,472 genes

---

## 1. Profiling Summary

### Test Datasets

| Metric | Small BAM (LBX0338) | Large BAM (MI_1190) | Ratio |
|--------|---------------------|---------------------|-------|
| BAM size | 203 MB | 1.2 GB | 5.9× |
| Total fragments | 4,145,225 | 24,149,525 | 5.8× |
| Buffered fragments | 728,785 | 9,250,691 | 12.7× |
| Unique alignments | 562,699 | 8,691,392 | 15.4× |
| Multimappers | 166,086 | 559,299 | 3.4× |
| **Total wall time** | **50.7s** | **947.7s** | **18.7×** |
| Throughput | 81,792 frags/s | 25,483 frags/s | 0.31× |
| Peak RSS | 2,860 MB | 6,142 MB | 2.1× |

### Stage Breakdown

| Stage | Small BAM | % | Large BAM | % | Scaling |
|-------|-----------|---|-----------|---|---------|
| **index_load** | 15.1s | 23.0% | 14.4s | 1.5% | 1.0× (fixed) |
| **scan_and_buffer** | 27.3s | 41.6% | 209.0s | 22.0% | 7.7× |
| **calibration** | 6.8s | 10.4% | 40.1s | 4.2% | 5.9× |
| **fragment_scorer** | 2.1s | 3.2% | 1.6s | 0.2% | 0.8× |
| **fragment_router_scan** | 0.9s | 1.4% | 9.4s | 1.0% | 10.4× |
| **build_loci** | 0.7s | 1.1% | 1.0s | 0.1% | 1.4× |
| **eb_gdna_priors** | 0.6s | 0.9% | 0.5s | 0.0% | 0.8× |
| **locus_em** | **12.1s** | **18.4%** | **686.1s** | **72.4%** | **56.6×** |

### Key Observation: Super-Linear EM Scaling

The EM stage exhibits catastrophically super-linear scaling:
- 12.7× more buffered fragments → **56.6× more EM time**
- Root cause: **mega-loci**. The largest locus in the large BAM has:
  - **315,030 transcripts** (51,656 in small → 6.1× more)
  - **5,777,023 units** (179,704 in small → 32× more)
- EM complexity scales as O(units × components × iterations), so a 32× increase in units with more components leads to the 56.6× slowdown

### Memory Profile (Large BAM)

| Phase | RSS (MB) | Growth |
|-------|----------|--------|
| After index load (baseline) | 2,860 | — |
| After scan_and_buffer | 3,882 | +1,022 MB (116 bytes/fragment) |
| After router_scan | 5,494 | +1,612 MB (CSR data structures) |
| After buffer release | 3,774 | −1,720 MB (buffer freed) |
| After locus_em | 4,814 | +1,039 MB (per-locus working memory) |
| After cleanup | 3,202 | −1,612 MB (CSR freed) |
| **Peak** | **6,142** | **+3,283 MB above baseline** |

**Peak occurs during router_scan** when both the buffer (1.5 GB) and CSR arrays coexist.

---

## 2. Bottleneck Analysis

### Bottleneck #1: Locus EM — 72.4% of total (large BAM)

**The problem:** One or two mega-loci with hundreds of thousands of transcripts and millions of fragment-units dominate EM runtime. The batch EM has two phases:
1. **Phase 1 (mega-loci):** Process sequentially with all 8 threads collaborating on intra-locus E-step
2. **Phase 2 (normal loci):** Dynamic work-stealing across threads

For the large BAM, almost all 686s is Phase 1 — iterating SQUAREM on a single mega-locus with 5.7M units × 315K transcripts.

**Root causes:**
- **No AVX2/AVX-512 exp() kernel on x86_64.** The E-step inner loop calls `fast_exp_scalar()` — a NEON 2-wide path exists for ARM64 but x86 falls back to scalar. This is the single hottest instruction sequence in Rigel.
- **Equivalence class column-sum accumulation is column-major on row-major data.** The E-step kernel iterates columns in the outer loop and rows inner, causing cache-line thrashing for wide equivalence classes.
- **Per-locus EC builder allocates fresh hash maps.** For Phase 2's thousands of small loci, this causes millions of small heap allocations.
- **SQUAREM convergence for mega-loci may be slow.** With 315K components (2×t + 1 = 630K+1), convergence requires many iterations.

**py-spy evidence:** `_em_impl.abi3.so` accounted for **1,033 out of 1,151 samples** (89.7%) in the py-spy profile. The top 5 hot addresses within the module are concentrated in 2-3 functions (likely `em_step_kernel_range` and `build_equiv_classes`).

### Bottleneck #2: BAM Scan — 22.0% of total (large BAM)

**The problem:** 209s for 9.25M fragments = ~44K fragments/s effective throughput (far below I/O limits).

**Architecture:** 3-tier producer-consumer pipeline:
- Reader thread (htslib + BGZF decompression with 2 sub-threads)
- N worker threads (resolve + accumulate)
- Main thread (Python callback with GIL)

**Likely sub-bottlenecks:**
- **htslib BGZF decompression** is a known bottleneck for BAM reading; only 2 decompression threads by default
- **Fragment resolution (cgranges queries)** for each fragment against 457K transcripts
- **Per-fragment `build_fragment()`** uses `std::unordered_map` + `std::set` for common 2-4 exon block cases
- **Queue backpressure:** Bounded queues (default sizes) may cause producer stalls

### Bottleneck #3: Calibration — 4.2% but 40s absolute

**The problem:** Python-level EM loop runs 17 iterations (large BAM), each calling:
- `estimate_kappa_marginal()` → Brent's optimization with ~19 evaluations of `_marginal_loglik()` over all regions
- `build_gdna_fl_model()` × 2 → re-bins millions of FL observations per iteration

**cProfile evidence (small BAM):**
- `_marginal_loglik`: 3.46s / 765 calls — the hot inner function
- `build_gdna_fl_model`: 0.42s / 83 calls
- `_compute_strand_llr_betabinom`: 0.49s / 41 calls
- `np.cumsum`: 227,982 calls / 0.35s — from locus build, not calibration

### Bottleneck #4: Memory — 6.1 GB peak for 1.2 GB BAM

**The problem:** Peak RSS is 5× the input BAM size. Breakdown:
- **Index:** ~2.8 GB fixed baseline (457K transcripts, splice junctions, interval trees, region data)
- **Fragment buffer:** ~1.0 GB (116 bytes/fragment × 9.25M fragments, CSR layout)
- **Scored CSR arrays:** ~1.6 GB (fragment_router_scan output — parallel arrays for unit scores, component indices, offsets)
- **EM working memory:** ~1.0 GB (per-locus EC flat arrays, scratch buffers)

---

## 3. Optimization Recommendations

### Priority 1 — CRITICAL: AVX2 Vectorized exp() for x86_64

**Impact:** Up to 4× speedup on E-step kernel (the single hottest loop)
**Effort:** Medium
**Risk:** Low (isolated change in `fast_exp.h`)

The current `fast_exp.h` has:
- `fast_exp_scalar()` — polynomial approximation (used on x86)
- NEON 2-wide path (ARM64 only)
- AVX2/AVX-512 macros defined but **no vectorized exp() implementations**

**Recommendation:** Implement `fast_exp_avx2()` processing 4 doubles/cycle and `fast_exp_avx512()` processing 8 doubles/cycle. The Cephes-style polynomial approximation already used in `fast_exp_scalar()` vectorizes trivially with `_mm256_*` intrinsics.

**Expected gains:**
- E-step kernel: 3-4× faster (4-wide AVX2 vs scalar)
- Locus EM overall: ~2.5-3× faster (accounting for non-exp work)
- Large BAM total: 686s → ~230-275s for EM alone

### Priority 2 — HIGH: Mega-Locus Decomposition / Pruning

**Impact:** Dramatic reduction in worst-case EM time
**Effort:** High
**Risk:** Medium (affects quantification accuracy for complex loci)

The single mega-locus with 315K transcripts dominates runtime. Options:

**A. Component pruning before EM:**
- After scoring, remove transcript components with negligible total likelihood (e.g., all fragments score below threshold for that transcript)
- This eliminates "phantom" transcripts from the EM that receive near-zero posterior anyway
- Could reduce effective components from 315K to <50K in many cases

**B. Locus decomposition via sparse connectivity:**
- Within a mega-locus, find weakly-connected sub-components in the fragment-transcript bipartite graph
- Fragments that share no transcript candidates don't need to be in the same EM problem
- Use BFS/DFS on the equivalence class structure to find disconnected sub-problems

**C. Iterative refinement:**
- Run initial EM with reduced iterations → identify near-zero components → prune → re-run
- Similar to the "online EM" approach used by Salmon

**D. Early termination for converged components:**
- Track per-component convergence; freeze components that stop changing
- Reduces effective EM width over iterations

### Priority 3 — HIGH: Optimize E-step Memory Access Pattern

**Impact:** 20-50% improvement on E-step cache efficiency
**Effort:** Low
**Risk:** Low

The column-sum accumulation in `em_step_kernel_range()` iterates `j` (columns/components) in the outer loop and rows in the inner loop. For EC data stored row-major in `ll_flat[n × k]`, this causes:
- Each inner row iteration strides by `k` doubles → poor spatial locality for large `k`

**Fix:** Transpose to row-major accumulation — iterate rows in outer loop, accumulate all columns per row. This reads each cache line once per row instead of `k` times.

For mega-loci with `k=630K`, this change alone could be transformative since L1 cache is typically 32KB and each row is `630K × 8B = 5MB`.

### Priority 4 — MEDIUM: Calibration C++ Port

**Impact:** Calibration 5-10× faster (40s → 4-8s)
**Effort:** Medium
**Risk:** Low

Move the calibration EM loop to C++:
1. **Pre-bin FL observations per region** (one-time O(n_fragments) → per-region histogram matrix)
2. **Fuse E-step channels** (density, strand, FL log-likelihood ratios) into single pass over regions
3. **Inline Brent's optimization** for kappa with fast `betaln` approximation (avoid 765 Python→C round trips)
4. **Warm-start Brent's** with previous iteration's κ value

Alternatively, a cheaper Python-level fix:
- Pre-compute per-region FL histograms once (O(n_regions × max_fl) matrix)
- Rebuild FL model via matrix-vector multiply instead of re-binning millions of observations
- This alone would save ~80% of `build_gdna_fl_model` time

### Priority 5 — MEDIUM: BAM Scanner Optimizations

**Impact:** scan_and_buffer 1.5-2× faster
**Effort:** Medium-High
**Risk:** Low

**A. Increase BGZF decompression threads:**
Current default is 2 sub-threads for BAM decompression. With 36 cores available, using 4-6 decompression threads would better saturate I/O.

**B. Small-buffer-optimized fragment assembly:**
Replace `std::unordered_map` + `std::set` in `build_fragment()` with stack-allocated small buffers for the common case (2-4 exon blocks per read pair). The scoring path already uses `MrnaScored[64]` stack allocation — apply same pattern here.

**C. Configurable worker thread count:**
Ensure the number of resolution workers scales with available cores. Profile the queue fill levels to identify if the producer or consumers are the bottleneck.

**D. Consider memory-mapped BAM reading:**
For BAMs on fast SSDs or in filesystem cache, mmap-based reading could reduce syscall overhead.

### Priority 6 — MEDIUM: Memory Reduction

**Impact:** 20-40% peak RSS reduction
**Effort:** Medium
**Risk:** Low

**A. Fragment buffer compaction:**
Current: 116 bytes/fragment. Potential savings:
- Use 32-bit indices instead of 64-bit where feasible (transcript indices fit in 32-bit for <4B transcripts)
- Pack boolean flags into bitfields
- Use delta encoding for sorted arrays (e.g., reference positions)
- Target: 64-80 bytes/fragment (30-45% savings)

**B. Streaming CSR construction:**
Currently, the full buffer is kept in memory until scoring is complete, then CSR arrays are built. Peak occurs when both coexist (+1.6 GB in the large BAM). Options:
- Process buffer chunks incrementally, building CSR on the fly and freeing processed buffer chunks
- Use memory-mapped temporary files for the buffer (already supported via `spill_dir` but not used by default)

**C. Index memory reduction:**
The index baseline is ~2.8 GB for human. This is dominated by:
- Interval tree data (1.2M collapsed intervals)
- Splice junction hash map (404K unique → 1.4M total)
- Transcript DataFrames (457K rows)
Possible savings via more compact representations (e.g., packed interval arrays instead of DataFrame, FlatHashMap instead of std::unordered_map for SJ lookup).

### Priority 7 — LOW: Index Load Time

**Impact:** 15s → 3-5s
**Effort:** Low-Medium
**Risk:** Low

Index load is 15s of fixed overhead. Current format uses feather files read via pyarrow/pandas. Options:
- Use numpy `.npy` for pure array data (skip pyarrow/pandas overhead)
- Memory-map the interval arrays instead of loading them fully
- Parallelize loading of independent feather files
- Cache the cgranges interval tree in serialized form

### Priority 8 — LOW: gc.collect() Overhead

**Impact:** 0.9s savings
**Effort:** Trivial
**Risk:** Low

`gc.collect()` takes 0.885s in the cProfile output. This is called explicitly after buffer release. Consider:
- Disabling GC during the compute-heavy phases
- Using `gc.freeze()` / `gc.unfreeze()` around the EM phase
- Or simply accept this as acceptable cleanup cost

---

## 4. Implementation Roadmap

### Phase 1: Quick Wins (1-2 days each)

| Item | Expected Impact | Effort |
|------|----------------|--------|
| AVX2 `fast_exp()` in `fast_exp.h` | EM 3× faster | 1 day |
| E-step row-major accumulation | EM 20-50% cache improvement | 0.5 day |
| Pre-bin FL observations for calibration | Calibration 2× faster | 0.5 day |
| Warm-start Brent's κ optimization | Calibration 30% fewer evals | 0.5 day |
| Increase BGZF decompression threads | Scan 10-20% faster | 0.5 day |

### Phase 2: Medium Effort (3-5 days each)

| Item | Expected Impact | Effort |
|------|----------------|--------|
| Mega-locus pruning (remove near-zero components) | EM 5-10× for worst-case loci | 3 days |
| Calibration C++ port | Calibration 5-10× faster | 3 days |
| Small-buffer fragment assembly | Scan 10-15% faster | 2 days |
| Fragment buffer compaction (32-bit indices) | 30% buffer memory savings | 2 days |

### Phase 3: Major Architecture (1-2 weeks)

| Item | Expected Impact | Effort |
|------|----------------|--------|
| Locus decomposition (sparse sub-problems) | EM 10-50× for mega-loci | 1 week |
| Streaming CSR construction | 30-40% peak RSS reduction | 1 week |
| AVX-512 E-step kernel | Additional 2× on AVX-512 servers | 3 days |

### Projected Impact

With Phase 1 + Phase 2 optimizations on the large BAM:

| Stage | Current | Projected | Improvement |
|-------|---------|-----------|-------------|
| scan_and_buffer | 209s | 170s | 1.2× |
| calibration | 40s | 8s | 5× |
| locus_em | 686s | 70-140s | 5-10× |
| **Total** | **948s** | **~260-330s** | **~3-3.5×** |

With all phases:

| Stage | Current | Projected | Improvement |
|-------|---------|-----------|-------------|
| scan_and_buffer | 209s | 120s | 1.7× |
| calibration | 40s | 4s | 10× |
| locus_em | 686s | 30-70s | 10-20× |
| **Total** | **948s** | **~170-210s** | **~4.5-5.5×** |

---

## 5. Profiling Infrastructure Notes

### What Worked
- `profiler.py --stages` gives excellent per-stage visibility
- cProfile captures Python-level hotspots well (calibration functions)
- py-spy with `--native` captures C++ module contribution percentages
- Memory timeline RSS sampling identifies peak memory phases

### What's Missing
- **C++ function-level profiling:** py-spy shows unresolved addresses because LTO strips debug symbols from installed `.so` files. Fix: either preserve debug info in a separate `.debug` file, or use `perf record` with `--call-graph dwarf` before LTO stripping.
- **Per-locus timing:** Add timing instrumentation to the batch EM to log per-locus wall time and identify which specific loci are the mega-locus bottlenecks.
- **Thread utilization:** No visibility into thread idle time during EM Phase 1 (how well does the E-step work partitioning utilize all 8 threads?). Add thread utilization counters.
- **I/O vs compute split in BAM scan:** Unclear if the scan bottleneck is I/O (htslib/BGZF) or computation (resolve/accumulate). Add per-stage counters within C++ scanner.
- **Equivalence class statistics:** Log EC count, median/max EC width, and compression ratio per locus to understand EM data characteristics.

### Recommended Profiling Enhancements
1. Add `--emit-locus-stats` flag to output per-locus timing and size data
2. Build a debug `.so` variant with preserved symbols for py-spy/perf symbol resolution
3. Add internal C++ timers (e.g., `std::chrono`) for sub-phases within scan and EM
4. Instrument thread pool utilization (idle time, task imbalance)
