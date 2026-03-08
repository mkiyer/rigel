# E-Step Kernel Rewrite Plan

## Motivation

Two performance/correctness issues in the current E-step implementation:

1. **Thread load imbalance in `parallel_estep`**: Partitions ECs by count, not
   by work. A single EC with n=50K, k=20 dominates one thread while others
   idle. Profiling shows exponential slowdowns on mega-loci.

2. **Precision swamping in `em_step_kernel` column sums**: Sequential
   accumulation of 500K+ normalized probabilities into a growing double drops
   low-order bits. Over hundreds of SQUAREM iterations, these errors amplify
   into non-deterministic trajectory divergence.

## Changes

### A. `KahanAccumulator` struct

New helper struct (~10 lines):

```cpp
struct KahanAccumulator {
    double sum = 0.0;
    double c = 0.0;
    inline void add(double val) {
        double y = val - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
};
```

Location: insert before `em_step_kernel`, after `EmEquivClass` and the EC builder.

### B. Replace `em_step_kernel` with `em_step_kernel_range`

Replace the current `em_step_kernel(ec, log_weights, em_totals)` with a
range-based variant that processes rows `[row_start, row_end)`:

```cpp
static inline void em_step_kernel_range(
    const EmEquivClass& ec,
    int row_start, int row_end,
    const double* log_weights,
    double* em_totals);
```

Key differences from current kernel:
- Takes `[row_start, row_end)` row range instead of processing all n rows.
- Fused pass: combines log-weight addition + log-sum-exp normalize in one
  loop over rows (better cache locality, same as current structure).
- **Kahan summation** on column-sum accumulation in Pass 2.
- Old `em_step_kernel` is deleted entirely.

Thread safety: multiple threads processing non-overlapping row ranges of the
same EC write to disjoint regions of `ec.scratch` (which is `mutable`).
Read-only access to `ec.ll_flat` and `ec.comp_idx`. Safe, no data races.

### C. Update all single-threaded E-step call sites

Four functions currently call `em_step_kernel` in a sequential loop:

| Function | Line | Current call |
|---|---|---|
| `map_em_step` | ~432 | `em_step_kernel(ec, ...)` |
| `vbem_step` | ~482 | `em_step_kernel(ec, ...)` |
| `linked_map_em_step` | ~547 (else branch) | `em_step_kernel(ec, ...)` |
| `parallel_estep` | ~374 (worker) | `em_step_kernel(ec_data[i], ...)` |

All sequential call sites change to:
```cpp
em_step_kernel_range(ec, 0, ec.n, log_weights.data(), em_totals);
```

The `parallel_estep` call site is replaced by the task-based worker (see D).

### D. Rewrite `parallel_estep` with task-based load balancing

Replace the current EC-count-based partitioning with granular task chunking:

1. **Task creation**: Break each EC into row-range tasks of ~4096/k rows,
   capped so each task fits comfortably in L1 cache. Each task records
   `{ec_idx, row_start, row_end, cost}` where cost = (row_end - row_start) * k.

2. **Cost-based partitioning**: Single greedy pass assigns contiguous task
   ranges to threads such that each thread gets approximately
   `total_cost / n_threads` work. Because tasks are bounded in size by the
   chunk cap, the greedy assignment produces near-optimal balance.

3. **Thread-local em_totals**: Each thread accumulates into its own
   `double[n_components]` buffer (unchanged from current design).

4. **Workers call `em_step_kernel_range`**: Each thread loops over its
   assigned tasks, calling `em_step_kernel_range(ec, row_start, row_end, ...)`.

5. **Deterministic Kahan reduction**: Cross-thread reduction sums
   thread-local buffers in thread-index order using `KahanAccumulator`
   (unchanged determinism guarantee, improved precision).

### E. Wire `estep_threads` through the non-linked path

Currently only `linked_map_em_step` → `linked_run_squarem` has `estep_threads`.
The classic (non-linked) path needs it too:

1. Add `int estep_threads = 1` parameter to `map_em_step`.
2. Add `int estep_threads = 1` parameter to `vbem_step`.
3. Add `int estep_threads = 1` parameter to `run_squarem`.
4. In `map_em_step` and `vbem_step`: branch on `estep_threads > 1` to call
   `parallel_estep` vs the sequential loop (same pattern as
   `linked_map_em_step`).
5. In `run_squarem`: thread `estep_threads` through to its `map_em_step` /
   `vbem_step` calls.
6. In `process_locus` lambda: pass `estep_thr` to `run_squarem` in the
   non-linked branch.

### F. Remove dead code

- Delete `em_step_kernel` (replaced by `em_step_kernel_range`).
- Delete `compute_estep_work` (no longer used — work proxy is computed
  from CSR offsets in the locus scheduler, not from EC data).

## Files modified

- `src/rigel/native/em_solver.cpp` — all changes above

## Validation

1. Build: `conda run -n rigel pip install --no-build-isolation -ve .`
2. Unit tests: `conda run -n rigel python -m pytest tests/ -x -q` (838 tests)
3. Benchmark against simulated truth to confirm accuracy is preserved.

## Interactions with existing two-phase locus scheduler

The two-phase sorted locus scheduling (Phase 1: mega-loci with all-thread
E-step, Phase 2: work-stealing with single-thread E-step) remains unchanged.
The improvements here make Phase 1 faster (better intra-locus load balance)
and more numerically stable (Kahan summation). Phase 2 benefits from Kahan
column sums in the sequential path.
