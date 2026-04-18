# P4: Fused Partition Scatter

**Date**: 2026-04-14
**Status**: Plan
**Estimated savings**: 2–3s wall time + significant code simplification
**Current partition cost**: 4.93s (4.5% of 108.8s total)

---

## 1. Problem Statement

The `partition_and_free()` function scatters 14 global CSR arrays into
29,458 per-locus partitions via **15 separate C++ calls** — each one
acquiring the GIL, extracting `n_loci` raw pointers from Python lists,
performing its scatter loop, wrapping results as numpy arrays, and
releasing the GIL.  The 6 candidate-scatter calls each iterate the same
CSR offset structure, and the 8 unit-scatter calls each iterate the same
unit-index arrays, producing 14 redundant passes over the locus index.

py-spy profiling attributes 2.2s of `memmove` to these scatter calls,
plus overhead from 15 GIL cycles, 15 × 29K pointer-extraction loops,
and poor cache reuse across passes.

### Current architecture: 15 calls, 14 passes

```
Python partition_and_free()
│
├─ build_partition_offsets()           1 C++ call, 1 pass
│
├─ scatter_candidates_f64(log_liks)    ─┐
├─ scatter_candidates_f64(cov_wts)      │  6 C++ calls
├─ scatter_candidates_i32(t_indices)    │  6 passes over same
├─ scatter_candidates_i32(tx_starts)    │  (offsets, locus_units,
├─ scatter_candidates_i32(tx_ends)      │   partition_offsets)
├─ scatter_candidates_u8(count_cols)   ─┘
│
├─ scatter_units_f64(gdna_log_liks)    ─┐
├─ scatter_units_i32(genomic_fp)        │  8 C++ calls
├─ scatter_units_i32(locus_t_idx)       │  8 passes over same
├─ scatter_units_u8(locus_count_cols)   │  locus_units
├─ scatter_units_u8(is_spliced)         │
├─ scatter_units_i64(frag_ids)          │
├─ scatter_units_u8(frag_class)         │
├─ scatter_units_u8(splice_type)       ─┘
│
└─ Python loop: assemble LocusPartition objects (29K iterations)
   └─ Pipeline: re-pack into 12-tuples for C++ EM (29K iterations)
```

**Total**: 15 C++ calls, 14 redundant passes, 2 Python reshuffling loops.

---

## 2. Key Insight: The Double Transpose

The scatter functions are column-oriented: each call produces one
`list[ndarray]` keyed by array name.  The EM solver is row-oriented: it
expects one `tuple[12 arrays]` per locus.  The current code performs an
explicit transpose through two Python loops:

1. **Assemble LocusPartition**: for each locus, pick column `li` from
   each of 14 result lists → dataclass (29K × 16 getattr/setattr ops).
2. **Pack 12-tuples**: for each locus, extract 12 fields from LocusPartition
   → tuple (29K × 12 attribute accesses).

A fused C++ function that scatters row-major — producing one complete
per-locus tuple in each iteration — eliminates both transpose loops and
the `LocusPartition` intermediate entirely from the EM hot path.

---

## 3. Design: Single `partition_all` Function

Replace all 15 C++ scatter functions and the Python orchestration with
**one** C++ function that performs the entire partition in two internal
passes (candidates + units) and returns results in exactly the format
the EM solver expects.

### Signature

```cpp
static nb::tuple partition_all(
    // Global CSR offsets
    i64_1d g_offsets,
    // Per-candidate globals (6 arrays)
    f64_1d log_liks,
    f64_1d coverage_weights,
    i32_1d t_indices,
    i32_1d tx_starts,
    i32_1d tx_ends,
    u8_1d  count_cols,
    // Per-unit globals (5 EM arrays)
    u8_1d  is_spliced,
    f64_1d gdna_log_liks,
    i32_1d genomic_footprints,
    i32_1d locus_t_indices,
    u8_1d  locus_count_cols,
    // Per-unit globals (3 metadata arrays)
    i64_1d frag_ids,
    u8_1d  frag_class,
    u8_1d  splice_type,
    // Locus structure
    nb::list locus_units,
    int n_loci
) -> nb::tuple  // (em_partition_tuples, meta_partition_tuples)
```

### Return value

```python
(
    em_partition_tuples,    # list[tuple[12]]  — ready for batch_locus_em_partitioned
    meta_partition_tuples,  # list[tuple[3]]   — (frag_ids, frag_class, splice_type)
)
```

Each `em_partition_tuples[li]` is the exact 12-tuple that `batch_locus_em_partitioned` 
expects:

```
(offsets, t_indices, log_liks, coverage_weights,
 tx_starts, tx_ends, count_cols,
 is_spliced, gdna_log_liks, genomic_footprints,
 locus_t_indices, locus_count_cols)
```

Each `meta_partition_tuples[li]` holds the annotation-only arrays:

```
(frag_ids, frag_class, splice_type)
```

### Internal structure: 2 passes, not 14

```cpp
// Phase 1: Extract all locus_units pointers (1 loop over n_loci)
struct LocusInfo {
    const int32_t* units;
    int n_units;
    int64_t n_candidates;
    int64_t* offsets;      // allocated once
    // Candidate output pointers (allocated once per locus)
    double*  log_liks;
    double*  cov_wts;
    int32_t* t_indices;
    int32_t* tx_starts;
    int32_t* tx_ends;
    uint8_t* count_cols;
    // Unit output pointers
    uint8_t* is_spliced;
    double*  gdna_log_liks;
    int32_t* genomic_fp;
    int32_t* locus_t_idx;
    uint8_t* locus_count_cols;
    int64_t* frag_ids;
    uint8_t* frag_class;
    uint8_t* splice_type;
};
// Allocate all arrays for all loci in Phase 1.

// Phase 2: Single candidate scatter pass
for (int li = 0; li < n_loci; ++li) {
    auto& info = infos[li];
    // Compute offsets (same as build_partition_offsets)
    // memcpy 6 candidate arrays in one unit loop
    for (int k = 0; k < info.n_units; ++k) {
        int u = info.units[k];
        int64_t g_start = goff[u];
        int64_t seg_len = goff[u + 1] - g_start;
        if (seg_len > 0) {
            memcpy(info.log_liks   + info.offsets[k], src_ll  + g_start, seg_len * 8);
            memcpy(info.cov_wts    + info.offsets[k], src_cw  + g_start, seg_len * 8);
            memcpy(info.t_indices  + info.offsets[k], src_ti  + g_start, seg_len * 4);
            memcpy(info.tx_starts  + info.offsets[k], src_ts  + g_start, seg_len * 4);
            memcpy(info.tx_ends    + info.offsets[k], src_te  + g_start, seg_len * 4);
            memcpy(info.count_cols + info.offsets[k], src_cc  + g_start, seg_len * 1);
        }
    }
}

// Phase 3: Single unit gather pass
for (int li = 0; li < n_loci; ++li) {
    auto& info = infos[li];
    for (int k = 0; k < info.n_units; ++k) {
        int u = info.units[k];
        info.is_spliced[k]       = src_is[u];
        info.gdna_log_liks[k]    = src_gl[u];
        info.genomic_fp[k]       = src_gf[u];
        info.locus_t_idx[k]      = src_lt[u];
        info.locus_count_cols[k] = src_lc[u];
        info.frag_ids[k]         = src_fi[u];
        info.frag_class[k]       = src_fc[u];
        info.splice_type[k]      = src_st[u];
    }
}

// Phase 4: Wrap all allocated arrays as numpy ndarrays,
//          assemble into per-locus tuples
```

### Why this is faster

| Factor | Before (14 calls) | After (1 call) |
|--------|-------------------|-----------------|
| GIL acquire/release | 15× | 1× |
| Pointer extraction loops | 15 × 29K = 435K | 1 × 29K = 29K |
| Candidate index iterations | 6 × Σ(n_units) | 1 × Σ(n_units) |
| Unit index iterations | 8 × Σ(n_units) | 1 × Σ(n_units) |
| Cache: `g_offsets[u]` lookups | 6 × 16M | 1 × 16M |
| Cache: `units[k]` lookups | 14 × 16M | 2 × 16M |
| Python LocusPartition assembly | 29K × 16 attrs | 0 |
| Python 12-tuple packing | 29K × 12 attrs | 0 |

The candidate scatter fuses 6 `memcpy` calls per unit-segment into a
single iteration, so the offset lookup (`g_offsets[u]`, `g_offsets[u+1]`)
and destination offset (`partition_offsets[k]`) are computed once and
reused for all 6 arrays.  At 16M units with ~2.5 candidates each, this
eliminates ~80M redundant offset lookups.

The unit gather fuses 8 indexed reads per unit into one loop, so the
indirection `units[k]` is resolved once per unit instead of 8×.

---

## 4. Python Side: Radical Simplification

### Before: partition.py (~120 lines)

```python
from .native import (
    build_partition_offsets,
    scatter_candidates_f64, scatter_candidates_i32, scatter_candidates_u8,
    scatter_units_f64, scatter_units_i32, scatter_units_i64, scatter_units_u8,
)

def partition_and_free(em_data, loci):
    # ... 14 scatter calls in loops, incremental free, assemble LocusPartition ...
```

### After: partition.py (~30 lines)

```python
from .native import partition_all

def partition_and_free(em_data, loci):
    locus_units = [locus.unit_indices for locus in loci]

    em_tuples, meta_tuples = partition_all(
        em_data.offsets,
        em_data.log_liks, em_data.coverage_weights,
        em_data.t_indices, em_data.tx_starts, em_data.tx_ends,
        em_data.count_cols,
        em_data.is_spliced.view(np.uint8) if em_data.is_spliced.dtype == np.bool_ else em_data.is_spliced,
        em_data.gdna_log_liks, em_data.genomic_footprints,
        em_data.locus_t_indices, em_data.locus_count_cols,
        em_data.frag_ids,
        em_data.frag_class.view(np.uint8) if em_data.frag_class.dtype == np.int8 else em_data.frag_class,
        em_data.splice_type,
        locus_units, len(loci),
    )

    # Free all global arrays at once
    for attr in ("offsets", "t_indices", "log_liks", "count_cols",
                 "coverage_weights", "tx_starts", "tx_ends",
                 "locus_t_indices", "locus_count_cols", "is_spliced",
                 "gdna_log_liks", "genomic_footprints", "frag_ids",
                 "frag_class", "splice_type"):
        setattr(em_data, attr, None)

    return em_tuples, meta_tuples
```

### Pipeline changes

`_run_locus_em_partitioned` currently receives `dict[int, LocusPartition]`
and re-packs into tuples.  With the fused scatter, it receives
`em_tuples` directly.  The change is:

```python
# Before:
partitions = partition_and_free(em_data, loci)
_run_locus_em_partitioned(estimator, partitions, loci, ...)

# After:
em_tuples, meta_tuples = partition_and_free(em_data, loci)
_run_locus_em_partitioned(estimator, em_tuples, meta_tuples, loci, ...)
```

Inside `_run_locus_em_partitioned`:

```python
# Before: _call_batch_em builds tuples from LocusPartition objects
partition_tuples = [
    (p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
     p.tx_starts, p.tx_ends, p.count_cols, p.is_spliced,
     p.gdna_log_liks, p.genomic_footprints, p.locus_t_indices,
     p.locus_count_cols)
    for p in parts
]

# After: em_tuples already in the right format
partition_tuples = [em_tuples[loc.locus_id] for loc in batch_loci]
```

`_populate_em_annotations` changes from accessing `p.frag_ids` etc. to
using `meta_tuples`:

```python
# Before:
frag_ids = np.concatenate([p.frag_ids for p in batch_parts])
frag_class = np.concatenate([p.frag_class for p in batch_parts])
splice_type = np.concatenate([p.splice_type for p in batch_parts])

# After:
frag_ids = np.concatenate([meta_tuples[loc.locus_id][0] for loc in batch_loci])
frag_class = np.concatenate([meta_tuples[loc.locus_id][1] for loc in batch_loci])
splice_type = np.concatenate([meta_tuples[loc.locus_id][2] for loc in batch_loci])
```

---

## 5. Memory Analysis

### Transient peak during scatter

The fused function holds all 14 global arrays alive for the duration of
the single call.  The current incremental approach frees each global
after its scatter.  The difference:

| | Current (incremental) | Fused (single call) |
|---|---|---|
| Globals alive during last scatter | 1 array (~16 MB) | All 14 (~930 MB) |
| Scattered results | All 14 lists (~930 MB) | All 14 lists (~930 MB) |
| Transient overhead | ~0 | ~900 MB |

The transient overhead is bounded at ~900 MB for ~2 seconds (the scatter
call duration) and drops to zero when the call returns and Python frees
the globals.  Given the current peak RSS of 7,044 MB and index baseline
of 2,836 MB, this transient spike is well within the already-allocated
memory headroom — the peak occurs later during EM, not during partition.

After the fused call returns and globals are freed, RSS matches the
current post-partition state exactly.

### Long-term: LocusPartition elimination saves object overhead

Eliminating 29,458 LocusPartition dataclass instances (each with 16
numpy array references) removes ~4.7M Python object-header bytes and
~470K reference-counting operations during the partition→EM handoff.

---

## 6. Code Cleanup Summary

### Removed

| Item | Lines | Location |
|------|-------|----------|
| `scatter_candidates_impl<T>` template + 3 bindings | ~65 | em_solver.cpp |
| `scatter_units_impl<T>` template + 4 bindings | ~55 | em_solver.cpp |
| `build_partition_offsets` function + binding | ~55 | em_solver.cpp |
| 7 scatter imports in native.py | 7 | native.py |
| 7 scatter imports in partition.py | 7 | partition.py |
| CAND_ARRAYS/UNIT_ARRAYS loops + assembly | ~70 | partition.py |
| LocusPartition from EM hot path | ~20 | pipeline.py |
| 7 scatter imports + tests for individual functions | ~100 | test_partition.py |

**Total removed**: ~380 lines across 5 files.

### Added

| Item | Lines | Location |
|------|-------|----------|
| `partition_all` function | ~130 | em_solver.cpp |
| 1 binding for `partition_all` | ~20 | em_solver.cpp |
| 1 import in native.py | 1 | native.py |
| Simplified `partition_and_free` | ~25 | partition.py |
| Updated `_run_locus_em_partitioned` | ~10 | pipeline.py |
| Tests for `partition_all` | ~80 | test_partition.py |

**Total added**: ~265 lines.

**Net**: −115 lines, fewer abstractions, one function replaces fifteen.

---

## 7. Implementation Steps

### Step 1: Add `partition_all` to em_solver.cpp

Write the fused C++ function with the signature from §3.  Internal
structure:

1. Extract all raw pointers under GIL (one loop over `locus_units`)
2. Release GIL
3. Candidate pass: for each locus, compute local offsets + `memcpy` all
   6 candidate arrays in one inner loop
4. Unit pass: for each locus, gather all 8 unit arrays in one inner loop
5. Re-acquire GIL
6. Wrap all allocated arrays as numpy ndarrays with capsule deleters
7. Assemble into per-locus tuples; return `(em_tuples, meta_tuples)`

Add the nanobind `m.def("partition_all", ...)` binding.

### Step 2: Update native.py

Replace:
```python
from ._em_impl import build_partition_offsets
from ._em_impl import scatter_candidates_f64
from ._em_impl import scatter_candidates_i32
from ._em_impl import scatter_candidates_u8
from ._em_impl import scatter_units_f64
from ._em_impl import scatter_units_i32
from ._em_impl import scatter_units_i64
from ._em_impl import scatter_units_u8
```

With:
```python
from ._em_impl import partition_all
```

Keep the old imports available temporarily for backward compatibility
with tests (removed in Step 5).

### Step 3: Rewrite partition.py

Replace the 14-call orchestration loop with a single call to
`partition_all`.  Change the return type from `dict[int, LocusPartition]`
to `tuple[list, list]` (em_tuples, meta_tuples).

### Step 4: Update pipeline.py

Update `_run_locus_em_partitioned` to accept `(em_tuples, meta_tuples)`
instead of `dict[int, LocusPartition]`.  Remove the tuple-packing loop
in `_call_batch_em`.  Update `_populate_em_annotations` to use
`meta_tuples`.

Update `quant_from_buffer` call site to pass the new return format.

Update the profiler similarly (same call pattern).

### Step 5: Update tests

- Rewrite `test_partition.py` to test `partition_all` directly and
  via the `partition_and_free` wrapper.
- Remove tests for individual `scatter_candidates_*` / `scatter_units_*`
  functions.
- Update `conftest.py` fixtures that call `partition_and_free` to handle
  the new return type.
- Run golden output tests to verify EM results are bit-identical.

### Step 6: Remove old scatter functions from em_solver.cpp

Delete `scatter_candidates_impl`, `scatter_units_impl`,
`build_partition_offsets`, and their 8 nanobind bindings.  Remove old
imports from native.py.

### Step 7: Compile, test, profile

```bash
pip install --no-build-isolation -e .
pytest tests/ -v
python scripts/profiling/profiler.py \
  --bam <BAM> --index <INDEX> --outdir /tmp/rigel_p4_profile \
  --threads 8 --stages
```

Verify:
- All 1028 tests pass
- Golden outputs unchanged (bit-identical EM results)
- Partition stage drops from ~5s to ~2–3s
- No RSS regression

---

## 8. Risk Assessment

| Risk | Mitigation |
|------|------------|
| Transient RSS spike (~900 MB) during fused call | Acceptable: well within current peak headroom, transient for ~2s |
| Changing `partition_and_free` return type breaks callers | Only 3 callers: pipeline.py, profiler.py, conftest.py — all updated in same PR |
| Removing LocusPartition from EM path loses debugging info | LocusPartition dataclass remains in scored_fragments.py for optional use; meta_tuples still carry annotation data |
| Bit-identical EM results | Same data, same layout, same C++ EM path — only the scatter transport changes |

---

## 9. Estimated Impact

| Metric | Before | After (est.) |
|--------|--------|--------------|
| Partition wall time | 4.93s | 2–3s |
| GIL acquire/release cycles | 15 | 1 |
| Passes over locus index | 14 | 2 |
| Python tuple-packing loops | 2 × 29K | 0 |
| Lines of scatter code | ~380 | ~265 |
| C++ scatter functions | 15 | 1 |
| Python scatter imports | 8 | 1 |
