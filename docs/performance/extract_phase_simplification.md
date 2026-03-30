# Extract-Phase Simplification & Dead Code Removal

## Summary

Three-phase change to simplify and optimize the locus EM extract phase:

1. **Rewrite test infrastructure** to exercise the live partitioned EM path
   instead of the dead non-partitioned path.
2. **Delete ~870 lines of dead code** across C++ and Python.
3. **Simplify the live extract function** by removing the epoch-based dedup
   machinery, pre-allocating output arrays, and using a reusable sort buffer.

## Motivation

### Dead code tested, live code under-tested

The production pipeline exclusively uses `batch_locus_em_partitioned()` →
`extract_locus_sub_problem_from_partition()`. However, 70 unit tests call
`batch_locus_em()` → `extract_locus_sub_problem()` through the conftest
helper `_run_and_assign()`. The non-partitioned path is dead production
code kept alive only by tests.

### Over-engineered extract logic

The current `extract_locus_sub_problem_from_partition()` uses an
epoch-based sparse-array approach that is both complex and cache-hostile:

- **`seen_epoch[nc]` + `best_buf[nc]`**: Two arrays sized to `n_components`
  (up to 160K entries ≈ 5 MB combined). Randomly accessed per unit,
  destroying L1/L2 cache for every fragment.
- **`dirty_comps` + per-unit `std::sort`**: A dynamically-sized component
  list that is sorted after each unit. Sorting overhead applied to 18.9M
  units in the mega-locus.
- **Six `push_back()` calls per candidate**: Hidden heap reallocations in
  the innermost loop.

### Dedup solves a non-existent problem

The epoch machinery deduplicates candidates with the same local component
index within a unit. However:

- The scorer (`scoring.cpp`) already merges per-transcript duplicates via
  `st.mm_merged` (a hashmap keyed by transcript index).
- The global→local mapping is injective (one global transcript maps to
  exactly one local component).
- Therefore, duplicates in the input CSR are impossible by construction.

The dedup machinery adds complexity and cache pressure for zero benefit.

## Phase 1: Rewrite test infrastructure

### Change: `tests/conftest.py` — `_run_and_assign()`

Replace the call to `rc.run_batch_locus_em()` with:

```python
from rigel.partition import partition_and_free

def _run_and_assign(rc, em_data, ...):
    # Unpack tuple from _make_locus_em_data
    if isinstance(em_data, tuple):
        em_data, loci, locus_gammas, index = em_data

    _ensure_estimator_geometry(rc)

    # Partition ScoredFragments into per-locus LocusPartition objects
    partitions = partition_and_free(em_data, loci)

    # Build the 12-tuples and transcript index lists
    partition_tuples = [
        (p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
         p.tx_starts, p.tx_ends, p.count_cols,
         p.is_spliced, p.gdna_log_liks, p.genomic_footprints,
         p.locus_t_indices, p.locus_count_cols)
        for p in [partitions[i] for i in range(len(loci))]
    ]
    locus_t_lists = [l.transcript_indices for l in loci]
    gdna_spans = np.array([l.gdna_span for l in loci], dtype=np.int64)

    total_gdna, _, _ = rc.run_batch_locus_em_partitioned(
        partition_tuples, locus_t_lists,
        locus_gammas, gdna_spans, index,
        em_iterations=em_iterations,
    )
    rc._gdna_em_total += total_gdna
    return {"mrna": ..., "nrna": ..., "gdna": ...}
```

### Verification

All 977 tests pass. The 70 tests in `test_estimator.py` and `test_gdna.py`
now exercise the live `batch_locus_em_partitioned()` C++ code path.

## Phase 2: Delete dead code

| File | Function / Block | Lines (approx) |
|------|------------------|-----------------|
| `em_solver.cpp` | `extract_locus_sub_problem()` | ~200 |
| `em_solver.cpp` | `batch_locus_em()` | ~470 |
| `em_solver.cpp` | nanobind binding for `batch_locus_em` | ~85 |
| `estimator.py` | `run_batch_locus_em()` | ~140 |
| `pipeline.py` | `_run_locus_em()` | ~60 |
| `native.py` | `batch_locus_em` import/export | ~3 |
| **Total** | | **~958** |

### Verification

All 977 tests pass after dead code removal.

## Phase 3: Optimize the live extract function

### 3a. Add no-duplicate-candidates test

Add a test in `test_em_impl.py` that constructs a multi-transcript locus
with ambiguous fragments, partitions the data, and verifies that no unit
in the partition CSR contains duplicate global transcript indices.

This proves the scorer contract and justifies removing the dedup machinery.

### 3b. Pre-allocate output arrays

Before the unit loop, compute the worst-case output size:

```
max_out = pv.n_candidates + pv.n_units   // every unit could get a gDNA candidate
```

Call `resize(max_out)` on all six output vectors (`t_indices`, `log_liks`,
`coverage_wts`, `tx_starts`, `tx_ends`, `count_cols`). Replace all
`push_back()` calls with direct indexed writes via a `write_cursor`
integer. After the loop, call `resize(write_cursor)` to trim.

This eliminates all dynamic allocation in the innermost loop.

### 3c. Replace epoch arrays with a reusable sort buffer

Declare a single `std::vector<LocalCandidate>` before the unit loop.
Reuse it across all units — `resize()` is a no-op when the existing
capacity is sufficient, so after processing the first unit there are
effectively zero allocations. The buffer stays cache-hot because it is
small (typically 2-3 elements) and reused every iteration.

Per-unit logic:

```
width_in  = offsets[ui+1] - offsets[ui]     // input RNA candidates
has_gdna  = (!is_spliced && isfinite(gdna_log_lik)) ? 1 : 0
width_out = width_in + has_gdna

buf.resize(width_out)

// Fill: remap global→local and store in buf
for j in [offsets[ui] .. offsets[ui+1]):
    local = local_map[global_t]
    buf[k++] = {local, log_lik, cov_wt, tx_start, tx_end, count_col}

// Append gDNA if eligible (component = n_t, always the largest index)
if has_gdna:
    buf[k++] = {n_t, gdna_ll, 1.0, 0, footprint, 0}

// Sort by local component index (required by equiv-class builder)
if width_out > 1:
    std::sort(buf.begin(), buf.end(), by local_comp)

// Write directly to pre-allocated output at write_cursor
for entry in buf:
    out_t_indices[cursor] = entry.local_comp
    out_log_liks[cursor]  = entry.log_lik
    ...
    cursor++

offsets[ui+1] = cursor
```

### What is deleted

- `BestCandidate` struct — replaced by the existing `LocalCandidate`
- `seen_epoch` vector (sized `n_components`, up to 160K entries)
- `best_buf` vector (sized `n_components`)
- `dirty_comps` vector + per-unit `std::sort` on it
- `current_epoch` counter
- All `push_back()` calls on the six output vectors

### Why this is fast

1. **No cache thrashing**: The reusable buffer is typically 2-3 entries
   (~80 bytes), fitting entirely in L1 cache. The epoch arrays were
   160K entries (~5 MB), blowing L1/L2 every iteration.

2. **Zero dynamic allocation**: Output arrays pre-allocated; sort buffer
   reused. The old code did 6× `push_back()` per candidate.

3. **Minimal sorting**: `std::sort` on 2-3 elements compiles to branchless
   compare-and-swap. The old code sorted a `dirty_comps` vector plus
   performed indirect reads from `best_buf` — two levels of indirection.

4. **Linear memory access**: Writing to pre-allocated output at a cursor
   is purely sequential. The old scatter pattern (epoch check → write to
   `best_buf[local]` → later read back sorted) had poor locality.

### Verification

- All 977 tests pass.
- Golden outputs (`tests/golden/`) are bit-for-bit identical.
- New test proving no duplicate candidates per (unit, transcript).

## Files modified

| File | Change |
|------|--------|
| `tests/conftest.py` | Rewrite `_run_and_assign()` to use partitioned path |
| `tests/test_em_impl.py` | Add no-duplicate-candidates test |
| `src/rigel/native/em_solver.cpp` | Delete dead functions; optimize extract |
| `src/rigel/estimator.py` | Delete `run_batch_locus_em()` |
| `src/rigel/pipeline.py` | Delete `_run_locus_em()` |
| `src/rigel/native.py` | Remove `batch_locus_em` import/export |
