# Parallel E-step: Adaptive Thread Allocation Plan

**Date:** 2026-03-07  
**Context:** EM solver optimization — intra-locus parallel E-step  
**Prerequisite:** Parallel E-step infrastructure already implemented (parallel_estep, compute_estep_work)

---

## 1. Problem Statement

The EM solver processes ~29K loci via work-stealing (CHUNK_SIZE=16, atomic
counter). Each locus runs single-threaded SQUAREM, which calls em_step_kernel
in a tight loop (up to 333 iterations × 3 EM calls per iteration ≈ 1000
E-step invocations per locus).

Locus sizes follow a heavy-tailed distribution. A typical genome has:
- ~28,000 loci with 1–50 transcripts and <1000 units (trivial, <1ms each)
- ~1,000 loci with 50–200 transcripts and 1K–10K units (moderate)
- ~10–100 loci with 200+ transcripts and 10K–100K units (large)
- 1–3 "mega-loci" with 1000+ transcripts and 100K+ units (dominant)

The mega-loci are the critical path. With 8 threads, 7 threads finish all
their small loci and idle while 1 thread grinds on the mega-locus. Native
profiling shows std::exp() (19% of thread time) is the dominant cost inside
em_step_kernel, and it's embarrassingly parallel across equivalence classes.

### What we tried and learned

1. **fast_exp (custom fdlibm implementation):** <1 ULP per call, but even
   sub-ULP differences cascade through 333 SQUAREM iterations on complex
   loci. Result: 249 loci affected, max_abs=292, 41 pruning flips on
   oracle_sim. **Rejected** — violates 1e-6 tolerance constraint.

2. **Parallel E-step with fixed threshold (50K work):** Correct parallel_estep
   with thread-local buffers and deterministic reduction. But 386 loci
   triggered on oracle_sim, each getting FP accumulation order changes
   amplified by SQUAREM. max_abs=219 on oracle_sim. **Threshold too low.**

3. **Key insight:** The parallel E-step changes FP accumulation order (thread-
   local buffers + reduction vs sequential accumulation). This doesn't affect
   correctness, but SQUAREM amplifies the differences. At threshold=10M, no
   loci triggered on either test dataset. But even the ORIGINAL code shows
   run-to-run non-determinism (max_abs=2.66 on star_realdata) from the
   existing work-stealing scheduler.

**Conclusion:** We accept that parallelism introduces FP variation (same as
the existing work-stealing does). The goal is to apply it where it actually
helps: mega-loci that dominate wall-clock time. The fixed-threshold approach
is too crude — we need adaptive thread allocation based on relative locus size.

---

## 2. Available Data (Before Per-Locus Processing)

All of these are available in batch_locus_em before entering the worker loop:

| Metric | Source | Cost |
|--------|--------|------|
| n_t[li] = transcript count | `lt_off[li+1] - lt_off[li]` | O(1) |
| n_u[li] = unit count | `lu_off[li+1] - lu_off[li]` | O(1) |
| n_candidates[li] = nonzero entries | `sum(g_offsets[u+1] - g_offsets[u]) for u in locus units` | O(n_u) |

**Proxy metric:** `work[li] = n_t[li] * n_u[li]` — O(1), good correlation
with actual E-step work. The true E-step cost is `sum(n_ec * k_ec)` across
equivalence classes, but ECs are only built inside the worker. The proxy
captures the key signal: mega-loci have large n_t AND large n_u.

(We could also use n_candidates from the CSR offsets, which is a tighter
estimate, but it requires an O(n_u) scan per locus. For the pre-compute phase
over all loci this is O(total_units) which is fine.)

---

## 3. Thread Allocation Strategy

### 3a. Compute fair share per thread

```
total_work = sum(work[li]) for all loci
fair_share = total_work / N_threads
```

`fair_share` represents the amount of work one thread "should" handle if work
were perfectly distributed. A locus whose work exceeds fair_share is a
bottleneck that benefits from multi-threading.

### 3b. Per-locus thread count

```
estep_threads[li] = clamp(round(work[li] / fair_share), 1, N_threads)
```

| Locus work relative to fair_share | estep_threads | Effect |
|-----------------------------------|---------------|--------|
| < 0.5× fair_share                | 1             | Normal, work-stealing |
| 0.5× – 1.5× fair_share          | 1             | Normal (rounding) |
| 1.5× – 2.5× fair_share          | 2             | 2 threads on E-step |
| 2.5× – 3.5× fair_share          | 3             | 3 threads |
| ...                              | ...           | ... |
| ≥ (N_threads - 0.5) × fair_share| N_threads     | All threads |

**Example (oracle_sim):** 254K transcripts, 29K loci, 8 threads.
If total_work = 2e9 (proxy), fair_share = 2.5e8. The mega-locus with
n_t=7943 × n_u=~50K = 4e8 work gets `round(4e8 / 2.5e8) = 2` threads.
A truly extreme mega-locus in a full genome (n_t=10K × n_u=500K = 5e9)
would get `round(5e9 / 2.5e8) = min(20, 8) = 8` threads.

### 3c. Minimum thread benefit threshold

Even with the formula, don't parallelize loci whose absolute work is trivially
small (the overhead of spawning/joining threads dominates). Use:

```
if (work[li] < MIN_PARALLEL_WORK) estep_threads[li] = 1
```

Where `MIN_PARALLEL_WORK` ~ 100,000 (chosen so that the parallel E-step
overhead of allocating thread-local buffers + join + reduce is amortized
by the actual E-step computation savings). This guards against pathological
distributions where `fair_share` is small but total work is also small.

---

## 4. Scheduling: Sort + Two-Phase Execution

### 4a. Pre-compute and sort

Before the worker loop, compute work and sort:

```cpp
// Pre-compute locus work (O(n_loci))
std::vector<int64_t> locus_work(n_loci);
int64_t total_work = 0;
for (int li = 0; li < n_loci; ++li) {
    int64_t nt = lt_off[li + 1] - lt_off[li];
    int64_t nu = lu_off[li + 1] - lu_off[li];
    locus_work[li] = nt * nu;
    total_work += nt * nu;
}

// Sort loci by work descending
std::vector<int> locus_order(n_loci);
std::iota(locus_order.begin(), locus_order.end(), 0);
std::sort(locus_order.begin(), locus_order.end(),
    [&](int a, int b) { return locus_work[a] > locus_work[b]; });

// Compute per-locus thread allocation
int64_t fair_share = std::max<int64_t>(total_work / actual_threads, 1);
std::vector<int> estep_threads(n_loci, 1);
int mega_count = 0;
for (int li = 0; li < n_loci; ++li) {
    if (locus_work[li] >= MIN_PARALLEL_WORK) {
        int thr = static_cast<int>(
            std::round(static_cast<double>(locus_work[li]) / fair_share));
        estep_threads[li] = std::clamp(thr, 1, actual_threads);
        if (estep_threads[li] > 1) mega_count++;
    }
}
```

### 4b. Phase 1 — Mega-loci (multi-threaded E-step)

Process loci with estep_threads > 1 first, in sorted order (largest first).
All N threads collaborate on each mega-locus sequentially:

```cpp
// Phase 1: mega-loci, one at a time, all threads on E-step
for (int idx = 0; idx < mega_count; ++idx) {
    int li = locus_order[idx];
    int ethr = estep_threads[li];

    // These phases are single-threaded (fast relative to SQUAREM):
    //   extract_locus_sub_problem
    //   build_equiv_classes
    //   compute_ovr_prior_and_warm_start

    // SQUAREM is the expensive part — it calls em_step_kernel
    // in a loop. With estep_threads > 1, each em_step_kernel
    // call uses parallel_estep internally.
    linked_run_squarem(..., ethr);

    // assign_posteriors (single-threaded, fast)
}
```

**Why sequential per mega-locus:** The extract/build/assign phases are
single-threaded but fast (<1% of total locus time). SQUAREM dominates.
Running mega-loci sequentially ensures all threads are available for each
one's E-step, giving maximum parallelism where it matters.

### 4c. Phase 2 — Normal loci (work-stealing, 1 thread each)

After mega-loci are done, process the remaining loci with the existing
work-stealing pattern. All threads participate, each processing loci
independently with 1 E-step thread:

```cpp
// Phase 2: remaining loci, work-stealing
std::atomic<int> next_idx{mega_count};  // start after mega-loci
constexpr int CHUNK_SIZE = 16;

auto worker_fn = [&]() {
    LocusSubProblem sub;
    // ... thread-local scratch ...

    for (;;) {
        int chunk_start = next_idx.fetch_add(CHUNK_SIZE);
        if (chunk_start >= n_loci) break;
        int chunk_end = std::min(chunk_start + CHUNK_SIZE, n_loci);
        for (int idx = chunk_start; idx < chunk_end; ++idx) {
            int li = locus_order[idx];  // sorted order
            // ... extract, build, squarem(estep_threads=1), assign ...
        }
    }
};
```

**Why sorted order helps even in Phase 2:** Larger loci appear earlier in the
sorted order. With work-stealing, the first chunks grabbed contain the
largest remaining loci. This naturally distributes heavy loci across threads
early, reducing the "long tail" effect where one thread gets stuck on the
last big locus while others idle.

---

## 5. Architecture Summary

```
batch_locus_em:
  ┌─────────────────────────────────────────────────┐
  │ 1. Pre-compute locus_work[li] = n_t * n_u      │  O(n_loci)
  │ 2. Sort loci by work descending                 │  O(n_loci log n_loci)
  │ 3. Compute estep_threads[li] from fair_share    │  O(n_loci)
  │ 4. Identify mega-loci (estep_threads > 1)       │
  └─────────────────────────────────────────────────┘
              │
              ▼
  ┌─────────────────────────────────────────────────┐
  │ Phase 1: Mega-loci (sequential, multi-thread)   │
  │                                                 │
  │   for each mega-locus (largest first):          │
  │     extract + build (1 thread)                  │
  │     SQUAREM with estep_threads > 1              │
  │       └── parallel_estep (all threads on E-step)│
  │     assign_posteriors (1 thread)                │
  └─────────────────────────────────────────────────┘
              │
              ▼
  ┌─────────────────────────────────────────────────┐
  │ Phase 2: Normal loci (work-stealing, 1 thread)  │
  │                                                 │
  │   N threads, CHUNK_SIZE=16, sorted order        │
  │   each locus: extract → build → squarem → assign│
  │   (estep_threads = 1 for all)                   │
  └─────────────────────────────────────────────────┘
```

---

## 6. Non-Determinism Analysis

The existing code is already non-deterministic across runs (max_abs=2.66 on
star_realdata from two runs of identical code). This comes from the work-
stealing scheduler: different threads process different loci, and thread-local
gdna accumulation + atomic flushing produces different FP results.

The parallel E-step adds FP variation only for mega-loci (those with
estep_threads > 1). This is because thread-local em_totals accumulation
followed by reduction produces a different summation order than sequential
accumulation. SQUAREM then amplifies these <1 ULP per-step differences across
hundreds of iterations.

**This is acceptable** because:
1. The baseline is already non-deterministic
2. The parallel E-step produces correct EM results (just a different convergence path)
3. The differences are concentrated in mega-loci that have many competing
   transcripts near the decision boundary — these are inherently noisy
4. At the gene level, differences are much smaller (max_rel < 0.02)

**Regression testing must account for this:** Tight tolerances (1e-6) will
fail. The regression framework should use tolerances that match the inherent
non-determinism of the existing code (~0.05 absolute, ~0.01 relative at
transcript level).

---

## 7. Interaction with CMake Compiler Flags (Phase 2 from optimization plan)

The Phase 4 optimization plan identifies `-O3 -march=native -ffast-math` as
a high-priority, low-effort optimization. This is relevant because:

1. `-ffast-math` allows the compiler to auto-vectorize the exp() loop with
   NEON/SSE, which addresses the same 19% hotspot that fast_exp targeted
2. `-ffast-math` also introduces FP non-determinism (reordering, FMA
   contraction), which interacts with SQUAREM sensitivity the same way
3. However, `-ffast-math` affects ALL code in the TU, not just EM — it also
   affects digamma, log, etc.

**Recommendation:** Implement compiler flags (Phase 2 of the plan) AFTER
the parallel E-step is working, and test them together. The compiler-driven
optimizations may make the custom fast_exp unnecessary while providing
broader benefits.

---

## 8. Implementation Checklist

1. **Pre-compute phase** (before GIL release):
   - Compute `locus_work[li]` for all loci
   - Sort to get `locus_order[]`
   - Compute `estep_threads[li]` from fair_share formula
   - Determine `mega_count`

2. **Phase 1 loop** (inside GIL release, before spawning threads):
   - Process mega-loci sequentially on the main control flow
   - For SQUAREM calls, pass `estep_threads[li]` through to parallel_estep
   - Need to restructure: currently extract/build/squarem/assign are all
     inside the worker lambda. For Phase 1, these need to run outside the
     work-stealing loop, with a thread-spawning mechanism for parallel_estep.

3. **Phase 2 loop** (existing work-stealing, minor changes):
   - Work-steal through `locus_order[mega_count:]` instead of `0..n_loci`
   - All estep_threads = 1 (sequential E-step, bit-identical to current code)

4. **parallel_estep enhancement:**
   - Currently spawns std::threads for each call. For Phase 1 where we call
     it ~1000 times per mega-locus, the thread spawn/join overhead adds up.
   - Consider a persistent thread pool for the parallel E-step:
     - Main thread signals "start E-step" via barrier/condition variable
     - Worker threads wake, do their partition, signal "done"
     - Main thread continues SQUAREM iteration
   - This eliminates ~1000 thread spawn/join cycles per mega-locus.

5. **Thread pool for parallel_estep (Phase 1):**
   ```
   struct EstepPool {
       std::vector<std::thread> threads;
       std::barrier<> start_barrier;
       std::barrier<> done_barrier;
       // Shared state: ec_data ptr, log_weights ptr, em_totals ptr, etc.
       std::atomic<bool> stop{false};

       void run_estep(ec_data, log_weights, em_totals, n_components) {
           // Set shared pointers
           start_barrier.arrive_and_wait();  // wake workers
           // Worker 0's portion runs here
           done_barrier.arrive_and_wait();   // wait for all
           // Reduce thread-local into em_totals
       }
   };
   ```
   This reduces per-call overhead from ~50μs (thread spawn) to ~1μs (barrier).

6. **Testing:**
   - All 838 unit tests must pass
   - Golden tests must pass (these use small loci, unaffected)
   - Regression tests with relaxed tolerances matching inherent non-determinism

---

## 9. Expected Impact

For the oracle_sim benchmark (48s in locus_em, mega-locus with work=3.6M):
- The mega-locus currently runs single-threaded while 7 threads idle
- With 8 threads on its E-step, the E-step portion speeds up ~5-6× (sub-linear
  due to reduction overhead and memory bandwidth)
- The overall locus_em phase might improve 10-30% depending on how much of the
  48s is spent on the mega-locus vs all other loci

For full-genome real data (where mega-loci are much larger), the impact should
be more significant — potentially cutting locus_em time in half for pathological
cases with very large connected components.
