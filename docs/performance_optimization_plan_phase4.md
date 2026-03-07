# Performance Optimization Plan — Phase 4: Native Profiling Analysis & Roadmap

**Date:** 2026-03-07
**Dataset:** Oracle BAM (1.4 GB, 24M reads / 12M fragments, 254K transcripts, 29K loci)
**Platform:** macOS 26.3 ARM64 (Apple Silicon), AppleClang 17.0.0, Python 3.12.13

---

## 1. Profiling Summary

### 1.1 Stage Timings (profiler.py --stages)

| Stage                  | Time (s) |   %   | Notes                          |
|------------------------|----------|-------|--------------------------------|
| **scan_and_buffer**    |   80.71  | 58.4% | BAM I/O + fragment resolution  |
| **locus_em**           |   47.96  | 34.7% | Batch C++ EM (4 threads)       |
| fragment_router_scan   |    4.31  |  3.1% | Python routing + scoring       |
| eb_gdna_priors         |    3.16  |  2.3% | Empirical Bayes gDNA           |
| fragment_scorer        |    1.49  |  1.1% | Score matrix construction      |
| Other (6 stages)       |    0.46  |  0.3% | build_loci, geometry, etc.     |
| **TOTAL**              | **138.1**|       | 174K frags/s, 18 GB peak RSS   |

### 1.2 cProfile (Python-side, GIL-held only)

Python sees 59.7 s cumulative, but 54.0 s is `_thread.lock.acquire` — the Python
main thread is blocked on the GIL while C++ extensions run. Meaningful Python
hotspots are small: `pandas ArrowArray.__getitem__` (1.9 s), `gc.collect` (0.35 s),
`scoring.py genexprs` (0.56 s), `numpy.asarray` (0.78 s).

### 1.3 Native Sampling (macOS `sample`, 180 s @ 1 ms)

66K-line native call-stack capture covering the full BAM scan + EM phases with
per-instruction sample counts. Key threads analyzed below.

---

## 2. Native Profiling Deep Dive

### 2.1 BAM Scanner Main Thread (80.7 s phase, 70,680 samples in `_bam_impl`)

The main reader thread pushes `QnameGroup` items one-at-a-time through a
`BoundedQueue<QnameGroup>` (capacity = `2 × n_workers`). The native sample
breaks down as follows:

| Operation                              | Samples | % of BAM time |
|----------------------------------------|---------|---------------|
| `std::mutex::lock()`                   | 22,640  | 32.0%         |
| `std::condition_variable::notify_one()`| 22,129  | 31.3%         |
| `std::mutex::unlock()`                 | 20,128  | 28.5%         |
| BAM parsing + qname logic              |  5,783  |  8.2%         |
| **Total in `_bam_impl`**              | 70,680  | 100%          |

**Conclusion:** 91.8% of the BAM scanner's time is spent on mutex
lock/unlock/notify — not on actual I/O or parsing. The `BoundedQueue` with
single-item push is the dominant bottleneck in the entire pipeline.

Worker threads (Arrow ThreadPool, 16 threads) show 121,544 / 121,720 samples
in `condition_variable::wait` — they are **starved**, spending >99.8% of their
time idle waiting for work.

### 2.2 htslib Thread Pool (BAM Decompression)

htslib `tpool_worker` threads decompress BGZF blocks via libdeflate:

| Operation                        | Samples | % of thread |
|----------------------------------|---------|-------------|
| `_pthread_cond_wait` (idle)      | 64,538  | 92.4%       |
| `bgzf_decode_func` → libdeflate | 4,962   |  7.1%       |
| `pthread_cond_broadcast`         |   311   |  0.4%       |
| Memory freeing                   |    29   |  0.0%       |

BAM decompression is not the bottleneck — it completes quickly between reads.
The serial `sam_read1()` in the main reader thread is what gates throughput.

### 2.3 EM Solver Worker Threads (48.0 s phase)

Each EM worker thread (e.g., Thread_7431265 with 28,982 samples):

| Operation                        | Samples | % of thread |
|----------------------------------|---------|-------------|
| `em_step_kernel` inner loop      |  9,585  | 33.5%       |
| — of which `std::exp()` calls    |  5,441  | 19.0%       |
| — of which EM arithmetic         |  4,144  | 14.5%       |
| Work-stealing + locus setup      | 19,032  | 66.5%       |
| `std::log()` (digamma/VBEM)      |    ~12  |  <0.1%      |

Within the EM kernel, `std::exp()` (from libsystem_m) is the single most
expensive call — it appears at **9 distinct call sites** within the
log-sum-exp normalize loop, each consuming 400–950 samples. The tight loop
structure (`for j in 0..k: scratch[j] = exp(scratch[j] - max_val)`) would
benefit enormously from vectorized exp.

The 66.5% in "work-stealing + locus setup" includes:
- `extract_locus_sub_problem()` data gathering (indirect memory access through CSR)
- Atomic `fetch_add` for work stealing (CHUNK_SIZE=16)
- SQUAREM outer iterations and convergence checks
- Per-locus vector allocations inside `LocusSubProblem`

---

## 3. Optimization Phases (Priority Order)

### Phase 1 — BAM Scanner Queue Redesign [CRITICAL]

**Impact estimate:** 40–60 s savings (50–75% of scan_and_buffer)
**Effort:** Medium

The native profiling reveals that the `BoundedQueue` single-item push pattern
causes 92% overhead from mutex contention. The queue capacity is only
`2 × n_workers` (e.g., 8 items for 4 workers), meaning the producer must
acquire the lock, push one `QnameGroup`, notify, and unlock — **for every
single fragment**.

#### 1a. Batch-Push with Configurable Chunk Size

Replace single-item `queue.push(group)` with batched pushes:

```cpp
// Producer side: accumulate a batch before pushing
constexpr size_t BATCH_SIZE = 128;
std::vector<QnameGroup> batch;
batch.reserve(BATCH_SIZE);

while (sam_read1(fp, hdr, b) >= 0) {
    // ... parse and build current_group ...

    if (/* qname boundary */) {
        batch.push_back(std::move(current_group));
        current_group = QnameGroup{};

        if (batch.size() >= BATCH_SIZE) {
            queue.push_batch(std::move(batch));
            batch.clear();
            batch.reserve(BATCH_SIZE);
        }
    }
}
// Flush remaining
if (!batch.empty()) queue.push_batch(std::move(batch));
```

Add `push_batch()` to `BoundedQueue`:

```cpp
void push_batch(std::vector<T> items) {
    std::unique_lock<std::mutex> lock(mutex_);
    not_full_.wait(lock, [this, &items] {
        return items_.size() + items.size() <= capacity_;
    });
    for (auto& item : items) {
        items_.push_back(std::move(item));
    }
    not_empty_.notify_all();  // Wake all workers at once
}
```

This reduces lock acquisitions by ~128×, and `notify_all` replaces per-item
`notify_one`.

#### 1b. Increase Queue Capacity

Change capacity from `2 * n_workers` to `1024` or more. The current tiny queue
forces constant synchronization. Larger capacity allows the producer to run
ahead of workers:

```cpp
BoundedQueue<QnameGroup> queue(1024);
```

Combined with batching, this would reduce mutex contention from ~65K samples
to ~500 samples (a 130× improvement).

#### 1c. Lock-Free Alternative (Optional)

If batching alone is insufficient, replace `BoundedQueue` with a lock-free
SPMC ring buffer (e.g., `moodycamel::ConcurrentQueue` or a custom SPSC ring
if we switch to one consumer that re-distributes). This eliminates all mutex
overhead but adds implementation complexity.

#### 1d. Worker-Side Batch Pop

Workers should also pop in batches to reduce lock contention on the consumer
side:

```cpp
bool pop_batch(std::vector<T>& out, size_t max_items) {
    std::unique_lock<std::mutex> lock(mutex_);
    not_empty_.wait(lock, [this] {
        return !items_.empty() || closed_;
    });
    if (items_.empty()) return false;
    size_t n = std::min(items_.size(), max_items);
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        out.push_back(std::move(items_.front()));
        items_.pop_front();
    }
    not_full_.notify_one();
    return true;
}
```

---

### Phase 2 — CMake Compiler Flags [HIGH PRIORITY]

**Impact estimate:** 10–20% overall (15–25 s), primarily in EM
**Effort:** Low

The current build uses only `-O2` (CMake Release default). No `-O3`,
`-march=native`, `-ffast-math`, or LTO. For a numerical workload dominated by
`std::exp()` and floating-point arithmetic, this is a significant missed
opportunity.

#### 2a. Add Optimization Flags

In `CMakeLists.txt`, add after the existing target definitions:

```cmake
# Aggressive optimization for numerical C++ extensions
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(RIGEL_OPT_FLAGS -O3 -march=native -DNDEBUG)

    # -ffast-math enables:
    #   - Reciprocal approximation instead of division
    #   - FMA (fused multiply-add) contraction
    #   - Reordering of floating-point operations
    #   - Assumes no NaN/Inf (we check explicitly in EM kernel)
    # Only apply to performance-critical extensions:
    target_compile_options(_em_impl PRIVATE ${RIGEL_OPT_FLAGS} -ffast-math)
    target_compile_options(_scoring_impl PRIVATE ${RIGEL_OPT_FLAGS} -ffast-math)
    target_compile_options(_bam_impl PRIVATE ${RIGEL_OPT_FLAGS})
    target_compile_options(_resolve_impl PRIVATE ${RIGEL_OPT_FLAGS})
    target_compile_options(_cgranges_impl PRIVATE ${RIGEL_OPT_FLAGS})
endif()
```

**Why `-ffast-math` matters for EM:**
- Allows the compiler to use `fmla` (fused multiply-add) ARM instructions
- Enables auto-vectorization of the inner `exp` loop with NEON
- May replace `std::exp` with a faster approximate implementation
- The EM kernel already guards against NaN/Inf with explicit checks

#### 2b. Enable Link-Time Optimization (LTO)

```cmake
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
```

LTO allows cross-TU inlining (e.g., inlining `em_step_kernel` at call sites)
and eliminates dead code across compilation units.

#### 2c. Verify OpenMP Linkage

The Phase 2 plan mentions OpenMP but the current code uses `std::thread` with
manual work-stealing, not OpenMP pragmas. If OpenMP is intended (e.g., for the
EM loop), verify it's actually linked and used. If the current `std::thread`
approach is intentional, remove unused OpenMP references from CMakeLists.txt.

---

### Phase 3 — EM Solver: Vectorized exp() [HIGH PRIORITY]

**Impact estimate:** 15–25 s savings (30–50% of locus_em)
**Effort:** Medium

The EM `em_step_kernel` inner loop (em_solver.cpp lines 237–243) is:

```cpp
for (int j = 0; j < k; ++j) {
    row[j] = std::exp(row[j] - max_val);
    row_sum += row[j];
}
```

This calls scalar `std::exp()` per element. On ARM64 this invokes Apple's
libsystem_m `exp()` which is accurate to full precision but slow for bulk use.

#### 3a. SIMD Vectorized exp with ARM NEON

Replace the scalar loop with a NEON-vectorized version:

```cpp
#include <arm_neon.h>

// Fast approximate exp using 6th-order polynomial (< 1 ULP error)
// Process 2 doubles at a time with NEON
static inline float64x2_t fast_exp_neon(float64x2_t x) {
    // Clamping and range reduction: x = n*ln(2) + r
    // Use Horner's method for polynomial evaluation
    // ...actual implementation from Cephes or Sleef library
}

// In em_step_kernel:
int j = 0;
double row_sum = 0.0;
for (; j + 2 <= k; j += 2) {
    float64x2_t v = vld1q_f64(row + j);
    v = vsubq_f64(v, vdupq_n_f64(max_val));
    v = fast_exp_neon(v);
    vst1q_f64(row + j, v);
    row_sum += vgetq_lane_f64(v, 0) + vgetq_lane_f64(v, 1);
}
// Scalar tail
for (; j < k; ++j) {
    row[j] = std::exp(row[j] - max_val);
    row_sum += row[j];
}
```

#### 3b. Use Sleef or xsimd Library

Rather than hand-rolling NEON intrinsics, consider using the
[Sleef](https://sleef.org/) or [xsimd](https://github.com/xtensor-stack/xsimd)
library which provides portable vectorized math functions:

```cpp
#include <sleef.h>
// Sleef_expd2_u10 — 2-wide double exp, 1.0 ULP accuracy
for (int j = 0; j + 2 <= k; j += 2) {
    double2 v = Sleef_expd2_u10(vld1q_f64(row + j));
    vst1q_f64(row + j, v);
}
```

This is especially effective because the EM kernel has **9 unrolled exp call
sites** (visible in the native sample as 9 distinct instruction addresses
hitting `exp()`).

#### 3c. Loop Restructuring for Better Vectorization

The current column-sum accumulation (lines 256–262) iterates rows-then-cols
with stride `k`. Transposing to column-major storage or tiling the accumulation
would improve cache locality and enable wider SIMD:

```cpp
// Current: column-major accumulation with indirect cidx
for (int j = 0; j < k; ++j) {
    double col_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        col_sum += scratch[i * k + j];    // stride-k access
    }
    em_totals[cidx[j]] += col_sum;
}

// Better: row-major accumulation into a local buffer
std::array<double, 512> local_sums{};  // stack-allocated, fits L1
for (int i = 0; i < n; ++i) {
    const double* row = scratch + i * k;
    for (int j = 0; j < k; ++j) {
        local_sums[j] += row[j];         // sequential access
    }
}
for (int j = 0; j < k; ++j) {
    em_totals[cidx[j]] += local_sums[j];
}
```

---

### Phase 4 — EM Solver: Memory Access & Work Distribution [MEDIUM]

**Impact estimate:** 5–15 s savings
**Effort:** Medium-High

The native sample shows 66.5% of EM thread time outside the inner kernel,
split between `extract_locus_sub_problem()` (CSR data gathering) and
work-stealing overhead. This phase targets those costs.

#### 4a. Pre-Extract Locus Data at Build Time

Currently, every EM iteration re-extracts the sub-problem data via indirect
CSR indexing in `extract_locus_sub_problem()`. For iterative EM (multiple
SQUAREM outer iterations), this redundant extraction is wasteful.

Move sub-problem extraction to a one-time setup phase before the EM loop:

```cpp
// Pre-extract all locus sub-problems once
std::vector<LocusSubProblem> sub_problems(n_loci);
parallel_for(n_loci, [&](int li) {
    extract_locus_sub_problem(sub_problems[li], ...);
});

// EM iterations reuse pre-extracted data
for (int iter = 0; iter < max_iter; ++iter) {
    parallel_for(n_loci, [&](int li) {
        em_step_kernel(sub_problems[li].equiv_classes, ...);
    });
}
```

**Caveat:** This increases peak memory by storing all sub-problems
simultaneously. For 29K loci, this should be manageable.

#### 4b. Sort Loci by Size (Descending)

The current work-stealing uses `CHUNK_SIZE=16` with no ordering. The largest
locus has 300 transcripts × 60K units = massive computation. If it's assigned
near the end, other threads idle waiting.

Sort loci largest-first so the big loci are scheduled early:

```cpp
// Sort indices by unit count descending
std::vector<int> order(n_loci);
std::iota(order.begin(), order.end(), 0);
std::sort(order.begin(), order.end(), [&](int a, int b) {
    return (lu_off[a+1] - lu_off[a]) > (lu_off[b+1] - lu_off[b]);
});
```

#### 4c. Reduce Per-Locus Allocation

`LocusSubProblem` contains `std::vector` members that are allocated/freed for
every locus. Thread-local scratch with maximum-size pre-allocation would
eliminate these allocations entirely. The worker thread already has a
`LocusSubProblem sub` but the vectors inside it may be resized per locus.

Ensure all vectors use `reserve()` once at thread start based on max locus
size, and use `clear()` + index tracking instead of resize:

```cpp
// At thread start:
sub.ll_flat.reserve(max_units * max_components);
sub.comp_idx.reserve(max_components);
sub.scratch.reserve(max_units * max_components);
```

---

### Phase 5 — BAM Scanner: Parallel Reading [MEDIUM]

**Impact estimate:** 10–20 s savings (requires architectural change)
**Effort:** High

Even after fixing the queue bottleneck, `sam_read1()` is serial. Parallel BAM
reading requires one of:

#### 5a. Multi-Region Parallel Scan

Split the genome into N regions, open N independent BAM file handles with
`hts_set_threads()`, and scan in parallel. Merge fragment buffers afterward.
Requires coordinate-sorted BAM (already the case for the oracle BAM).

```cpp
// Pseudocode
std::vector<Region> regions = split_genome(chromosomes, n_workers);
std::vector<std::future<ScanResult>> futures;
for (auto& region : regions) {
    futures.push_back(std::async([&] {
        return scan_region(bam_path, region);
    }));
}
// Merge results
for (auto& f : futures) merge(global_result, f.get());
```

#### 5b. htslib Worker Thread Tuning

The native profiling shows htslib's thread pool is 92.4% idle. This is because
`sam_read1()` is called sequentially and decompression is not the bottleneck.
Adding more htslib threads via `hts_set_threads()` won't help unless the
reading pattern is changed to request multiple records in parallel.

However, increasing `hts_set_threads()` may provide marginal benefit by
pre-decompressing BGZF blocks ahead of the reader.

---

### Phase 6 — Python Layer Optimization [LOW]

**Impact estimate:** 2–5 s savings
**Effort:** Low-Medium

The cProfile shows several Python-side costs that are small individually but
add up:

#### 6a. Eliminate `gc.collect()` Call (0.35 s)

If `gc.collect()` is called explicitly between stages:

```python
# Remove or defer:
# gc.collect()  # 0.35s wasted
```

#### 6b. Vectorize scoring.py Generator Expressions (0.56 s)

Lines 180–181 of `scoring.py` have generator expressions consuming 0.56 s
total across 1.9M calls. Replace with bulk numpy operations.

#### 6c. Optimize Pandas Arrow Access (1.9 s)

`ArrowArray.__getitem__` accounts for 1.0 s over 745K calls. If this is
per-element access, switch to bulk `.to_numpy()` conversion:

```python
# Instead of: [arr[i] for i in range(len(arr))]
# Use:        arr.to_numpy()
```

---

## 4. Expected Impact Summary

| Phase | Target                          | Est. Savings | New Total | Effort |
|-------|---------------------------------|--------------|-----------|--------|
| 1     | BAM queue batching              | 40–60 s      | 78–98 s   | Medium |
| 2     | CMake -O3/-ffast-math/LTO       | 15–25 s      | 60–80 s   | Low    |
| 3     | Vectorized exp() in EM          | 15–25 s      | 50–65 s   | Medium |
| 4     | EM memory access + scheduling   | 5–15 s       | 40–55 s   | Med-Hi |
| 5     | Parallel BAM reading            | 10–20 s      | 30–45 s   | High   |
| 6     | Python layer cleanup            | 2–5 s        | 28–40 s   | Low    |

**Optimistic floor: ~30 s (4.6× speedup from 138 s)**
**Conservative target: ~55 s (2.5× speedup)**

Phases 1–3 alone are expected to deliver a 2× speedup with moderate effort.

---

## 5. Measurement Plan

Each phase should be validated with the same profiling methodology:

```bash
# Stage timing (always)
python scripts/profiler.py \
    --bam .../sim_oracle.bam \
    --index .../rigel_index \
    --outdir .../profile_phaseN \
    --stages

# Native sampling (for C++ changes)
python scripts/profiler.py ... &
BGPID=$!; sleep 2
sample $BGPID 180 1 -file sample_phaseN.txt
wait $BGPID
```

Key metrics to track:
1. `scan_and_buffer` time (target: <30 s)
2. `locus_em` time (target: <20 s)
3. Peak RSS (target: <16 GB)
4. Throughput (target: >400K frags/s)

---

## 6. Relationship to Prior Plans

This document builds on and supersedes the analysis portions of:
- `performance_optimization_plan.md` — Original P1–P7 roadmap
- `performance_optimization_plan_phase2.md` — Batch EM + OpenMP (completed)
- `performance_optimization_plan_phase3.md` — Fused score buffer (completed)

The key new insight from native profiling is that **BAM scanner mutex
contention — not I/O, not parsing, not decompression — is the #1 bottleneck**,
accounting for 92% of the scan phase. This re-prioritizes targets: queue
redesign (Phase 1 here) should come before all other optimizations.

For the EM solver, the new data confirms that `std::exp()` is the dominant
arithmetic cost (~19% of thread time) but that memory access and work
distribution (66.5%) are equally important — suggesting that compiler flags
(Phase 2) and data structure optimization (Phase 4) may yield more benefit
than SIMD alone (Phase 3).
