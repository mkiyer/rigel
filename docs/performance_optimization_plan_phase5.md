# Performance Optimization Plan — Phase 5

## 1  Current Performance Baseline

All measurements below were collected on Apple M4 Max (ARM64), macOS 26.3,
`rigel` v0.1.0-dev built with `-O3 -march=native -ffp-contract=fast -flto`,
4 EM threads, Python 3.12.13 (conda `rigel` env).

### 1.1  Stage-Level Profile (3 Datasets)

| Dataset | Total | scan\_buffer | locus\_em | frag\_router | Peak RSS | Loci | Max t×u |
|---------|-------|-------------|-----------|-------------|----------|------|---------|
| Pristine sim (20 M reads) | 92.2 s | 39.9 s (43%) | 46.1 s (50%) | 3.3 s (4%) | 12,596 MB | 18,441 | 260 × 59,706 |
| Contaminated sim (24 M reads) | 179.9 s | 42.3 s (24%) | 129.9 s (72%) | 4.3 s (2%) | 18,569 MB | 28,992 | 300 × 60,211 |
| Real clinical (21.6 M reads) | 15.3 s | 3.7 s (24%) | 5.7 s (38%) | 3.8 s (25%) | 3,805 MB | 14,635 | 42,617 × 292,595 |

Key observations:
- **EM dominates contaminated data** (72% of wall time) because gDNA / nRNA
  contamination creates many more ambiguous fragments requiring iterative EM.
- **Scan dominates pristine data** (43%) because the EM has less work.
- **Real clinical data** has a massive mega-locus (42,617 transcripts × 292,595
  units) but most reads are filtered out by the scanner (only 1.1 M fragments
  buffered from 21.6 M reads), keeping overall runtime low.
- **Peak RSS is very high**: 18.6 GB for contaminated, 12.6 GB for pristine.

### 1.2  Native Self-Time Profile (Contaminated Dataset)

Collected via macOS `sample` (1 ms interval, 300 s) attached to a 4-thread
contaminated run (191.2 s total with sampling overhead).

#### EM Worker Threads (4 threads, 461,754 self-time samples)

| Category | Samples | % of EM | Projected wall-time |
|----------|---------|---------|---------------------|
| EM arithmetic (`_em_impl.abi3.so`) | 270,703 | 58.6% | ~81 s |
| `std::exp()` (`libsystem_m.dylib`) | 160,177 | 34.7% | ~48 s |
| DYLD-STUB$$exp (dynamic linker) | 27,575 | 6.0% | ~8 s |
| `std::log()` | 555 | 0.1% | ~0.2 s |
| `_platform_memset` | 1,042 | 0.2% | — |
| malloc/free | ~800 | 0.2% | — |

**exp() total (including DYLD stub): 40.7% of EM worker CPU → ~56 s of 138.8 s
EM wall time.**

#### Main Thread (147,842 self-time samples)

| Category | Samples | % |
|----------|---------|---|
| `__ulock_wait` (joining EM workers) | 115,440 | 78.1% |
| `__psynch_cvwait` (condvar) | 8,835 | 6.0% |
| `__psynch_mutexwait` | 7,500 | 5.1% |
| `_scoring_impl` | 2,591 | 1.8% |

Main thread is almost entirely blocked on EM thread join — good utilization.

#### BAM Worker Threads (OTHER category, 5 threads)

| Category | Samples | % |
|----------|---------|---|
| `__workq_kernreturn` (idle) | 147,842 | 57.7% |
| `__psynch_mutexwait` (lock contention) | 18,790 | 7.3% |
| `__psynch_cvwait` | 13,240 | 5.2% |
| `_bam_impl` code | 11,163 | 4.4% |
| malloc/free | 37,000 | ~14% |

BAM workers have significant lock contention and malloc pressure, but since
scan is not the bottleneck on contaminated data, the impact is modest.

#### Idle Threads (37 threads, 4.6 M samples — 99.8% `__psynch_cvwait`)

15 of these are **numpy OpenMP threads** (`libomp.dylib`) that never do useful
work during the EM phase.  The remaining are HTS decompression pool threads,
Arrow threadpool workers, and Python GC threads — all idle during EM.

---

## 2  Optimization Opportunities (Ranked by Expected Impact)

### P1  Vectorized exp() — NEON/AVX2 Intrinsics for the E-Step Kernel

**Impact: ~56 s → ~14–20 s for EM exp() (save 36–42 s on contaminated)**

`std::exp()` consumes 40.7% of EM worker CPU time.  The system `libsystem_m`
`exp()` is a high-accuracy scalar implementation.  For posterior normalization,
4–5 ULP accuracy is more than sufficient.

**Approach:** Implement a `fast_exp_neon()` (ARM64) and `fast_exp_avx2()`
(x86-64) using the standard polynomial approximation for `2^(x/ln2)`:

1. Split `x` into integer part `n` and fractional part `f = x - n·ln2`.
2. Compute `2^n` via bit manipulation (set exponent field in IEEE 754).
3. Evaluate `exp(f)` via degree-5 or degree-6 minimax polynomial.
4. Multiply: `2^n × exp(f)`.

This processes 2 doubles/cycle on NEON (2-wide `float64x2_t`), or 4 doubles
per AVX2 `__m256d` register.  Expected 3–4× speedup over scalar `std::exp()`.

**Changes required:**
- Add `fast_exp.h` with NEON and AVX2 implementations, guarded by
  `simd_detect.h` macros.
- Modify `em_step_kernel_range()` Pass 1 to use vectorized exp:
  - Process `k` columns in SIMD-width batches (2 or 4 at a time).
  - Handle tail elements with scalar fallback.
- Include unit tests comparing fast_exp vs std::exp for the relevant input
  range (typically `[-700, 0]` for log-probability differences).
- Benchmark accuracy: verify max relative error < 1e-12 across the range.

**DYLD-STUB elimination:** The dynamic linker stub for `exp()` adds 6% on top.
With inline SIMD intrinsics, the stub is eliminated entirely — both the 34.7%
exp() and 6.0% DYLD-STUB overhead are replaced.

**Risk:** Low.  The E-step only needs posteriors that sum to 1 per row; ULP
differences in exp() are absorbed by normalization.

### P2  Thread Pool for parallel_estep — Eliminate Spawn/Join Overhead

**Impact: Moderate — eliminates per-E-step thread creation overhead**

Currently `parallel_estep()` spawns fresh `std::thread` objects every call.
For a mega-locus with 300 SQUAREM iterations × 3 EM steps = 900 E-step calls,
that's 900 × (n_threads−1) = 2,700 thread create/join cycles.  Each
`pthread_create` + `pthread_join` costs ~20–50 µs on macOS.

**Approach:** Create a lightweight thread pool (persistent worker threads +
task queue) that lives for the duration of `batch_locus_em()`:

1. Spawn `n_threads - 1` worker threads at the start of `batch_locus_em()`.
2. Workers wait on a condition variable for task assignments.
3. `parallel_estep()` posts row-range tasks and signals workers.
4. Main thread participates as worker 0.
5. After all Phase 1 mega-loci are done, signal workers to exit.

This also reduces the main thread's `__ulock_wait` time (78% of main thread
samples) since join overhead is eliminated.

**Changes required:**
- Add a simple `ThreadPool` class with `submit_batch()` / `wait()` interface.
- Modify `parallel_estep()` to accept a `ThreadPool&` instead of spawning.
- Modify `batch_locus_em()` to own a `ThreadPool` for the batch scope.

**Risk:** Low.  Thread pool pattern is well-understood.

### P3  Suppress numpy OpenMP Idle Threads

**Impact: Eliminates 15+ idle threads wasting resources**

The native profile shows 15 threads from `libomp.dylib` (numpy's OpenMP
runtime) sitting idle in `__kmp_fork_barrier` for the entire EM phase.  These
threads consume kernel scheduling resources and ~4.6 M idle sample events.

**Approach:** Set `OMP_NUM_THREADS=1` before importing numpy, or call
`omp_set_num_threads(1)` from C++ at module init.  This prevents numpy from
spawning its threadpool since rigel's BAM scan and EM already manage their
own parallelism.

**Changes required:**
- One line in rigel's C++ module init or Python `__init__.py`:
  `os.environ.setdefault("OMP_NUM_THREADS", "1")`.

**Risk:** Minimal.  Should be gated to avoid affecting users who import rigel
alongside numpy-heavy workloads.  Use `setdefault` so users can override.

### P4  Memory Reduction — Release Builder & Context After Scan

**Impact: ~1.3 GB peak RSS reduction (already designed, not yet implemented)**

Per `docs/memory_reduction_plan_phase5_update.md` P1:
- `del builder` after `scan_and_buffer()` returns — frees the fragment builder.
- `del ctx` after finalize — frees the scan context.
- Expected savings: ~1.3 GB on 24 M fragment contaminated dataset.

**Changes required:**
- 2 lines in `pipeline.py` (already designed).

### P5  Two-Pass Fused Scanner — Halve Peak Memory

**Impact: ~2.6 GB peak RSS reduction (designed in `fused_score_buffer_plan.md`)**

Currently the scan + buffer pipeline involves 3 memory copies:
Python list → numpy concat → typed arrays.  A two-pass C++ scanner would:

1. **Count pass:** Determine exact array sizes.
2. **Fill pass:** Write directly into pre-allocated numpy arrays.

This eliminates the intermediate Python list of tuples and the concat step,
roughly halving the scan's memory footprint from ~5.2 GB to ~2.6 GB.

**Changes required:**
- Significant C++ refactor of `_bam_impl` scan interface.
- See `docs/fused_score_buffer_plan.md` for full design.

### P6  Column-Major Layout for E-Step Pass 2

**Impact: Moderate — improves cache locality for Kahan column sums**

The E-step kernel's Pass 2 computes column sums over `scratch[i*k + j]` for
each column `j`.  This is stride-`k` access pattern.  For large equivalence
classes (e.g., k=300, n=200), each column access skips 2,400 bytes — poor
cache line utilization.

**Approach:** After Pass 1 normalization, transpose the `n×k` scratch matrix
to `k×n` layout, then do row-wise sums for Pass 2.  Alternatively, fuse the
column accumulation into Pass 1 using thread-local column accumulators.

The fused approach is preferable:
```cpp
// Pass 1 with fused column accumulation
for (int i = row_start; i < row_end; ++i) {
    // ... exp + normalize as before ...
    for (int j = 0; j < k; ++j) {
        col_acc[j].add(row[j]);  // Kahan accumulate during normalize
    }
}
```

This eliminates Pass 2 entirely and halves the memory traffic (scratch is
only written, never re-read).

**Changes required:**
- Modify `em_step_kernel_range()` to fuse column sums into Pass 1.
- Remove separate Pass 2 loop.
- Verify numerical equivalence (Kahan order changes slightly).

**Risk:** Low-medium.  The Kahan sum order changes from column-major to
row-major, which may produce ULP-level differences.  Need to verify test
suite still passes.

### P7  Reduce Hash Map Overhead in build_equiv_classes

**Impact: Small — reduces per-locus setup overhead**

`build_equiv_classes()` uses `unordered_map<vector<int32_t>, ...>` which
does many small heap allocations for the key vectors.  For loci with 10K+
units, this creates significant allocator pressure.

**Approach options:**
1. **Cache EC structure across SQUAREM iterations:** The equivalence class
   grouping depends only on the CSR structure, not on `theta`, so it can be
   computed once per locus and reused across all EM iterations.  Currently
   this is already the case (built once in `process_locus`), so no change
   needed there.
2. **Replace `vector<int32_t>` keys with hash of the key:** Use a 128-bit
   hash as the map key instead of storing the full key vector.  This
   eliminates per-key allocations.
3. **Pre-size the hash map:** Use `reserve(n_units)` (already done).

**Risk:** Low.  Option 2 has a theoretical collision risk but with 128-bit
hashes it's negligible.

### P8  BAM Worker Lock Contention

**Impact: Small on contaminated data, moderate on pristine data**

BAM workers spend 7.3% of their time on `__psynch_mutexwait`.  On pristine
data where scan accounts for 43% of wall time, reducing lock contention
could save 1–2 seconds.

**Approach:** Profile the specific mutexes involved.  Likely candidates:
- `BoundedQueue` SPMC queue contention between producers.
- Fragment builder lock during buffer append.

Potential fixes:
- Increase queue capacity to reduce producer blocking.
- Use per-worker staging buffers with periodic bulk commits.

---

## 3  Implementation Priority

| Priority | Item | Expected Impact | Effort | Risk |
|----------|------|----------------|--------|------|
| **P1** | Vectorized exp() | **-36–42 s** (contaminated) | Medium | Low |
| **P2** | Thread pool for parallel_estep | -5–10 s, better scaling | Medium | Low |
| **P3** | Suppress OpenMP idle threads | Cleaner resource usage | Trivial | Minimal |
| **P4** | `del builder`/`ctx` | **-1.3 GB** RSS | Trivial | Minimal |
| **P5** | Two-pass fused scanner | **-2.6 GB** RSS | High | Medium |
| **P6** | Fused column sums | -5–10 s (cache improvement) | Low | Low-medium |
| **P7** | Hash map overhead | -1–2 s on large loci | Low | Low |
| **P8** | BAM lock contention | -1–2 s (pristine only) | Low | Low |

### Recommended Implementation Order

1. **P1 + P3** — Highest impact: vectorized exp saves ~40 s on contaminated,
   OMP_NUM_THREADS=1 is a one-liner.
2. **P6** — Fused column sums: moderate impact, low effort.
3. **P2** — Thread pool: good engineering improvement, measurable on mega-loci.
4. **P4** — Memory: two-line fix for 1.3 GB reduction.
5. **P5** — Memory: larger project, 2.6 GB reduction.
6. **P7 + P8** — Diminishing returns.

### Projected Total Impact

With P1 + P2 + P6 implemented, contaminated dataset projection:
- EM: 138.8 s → ~80–90 s (exp speedup + cache improvements + thread pool)
- Total: 179.9 s → ~120–130 s (~30% improvement)

With all items:
- Total: 179.9 s → ~110–120 s
- Peak RSS: 18.6 GB → ~14.7 GB (P4 + P5)

---

## 4  Profiling Methodology Notes

### Data Locations

| Dataset | Path | Size | Reads |
|---------|------|------|-------|
| Pristine | `sim_ccle_hela_cervix_short/gdna_none_ss_1.00_nrna_none/sim_oracle.bam` | 967 MB | 20 M |
| Contaminated | `sim_ccle_hela_cervix_short/gdna_high_ss_0.90_nrna_low/sim_oracle.bam` | 1.4 GB | 24 M |
| Real clinical | `mctp_LBX0069_SI_42153_HFFFMDRX7/bam/star.collate.markdup.bam` | 1.7 GB | 21.6 M |
| Index | `rigel_index/` | — | 254,461 tx, 63,472 genes |

### Tools Used

- **Stage profiler:** `scripts/profiler.py --stages --threads 4`
- **Native sampling:** macOS `sample <pid> 300 1` (1 ms interval, 300 s)
- **Sample analysis:** `scripts/_analyze_em_threads.py` (self-time extraction
  from hierarchical macOS sample output)
- **Key insight:** macOS `sample` reports inclusive counts in a call tree;
  accurate self-time requires computing `self = inclusive - sum(children)`.

### Reproducing

```bash
# Stage-level profile
conda run -n rigel python scripts/profiler.py \
  --bam <bam_path> --index <index_path> \
  --outdir <output_dir> --stages --threads 4

# Native sample (attach to running profiler PID)
sudo sample <pid> 300 1 -file <output.txt>
```
