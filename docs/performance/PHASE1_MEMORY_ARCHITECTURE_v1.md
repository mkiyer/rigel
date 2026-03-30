# Phase 1: Memory Architecture — Implementation Plan

## Executive Summary

Replace the monolithic `ScoredFragments` global CSR with a **locus-partitioned CSR** that allows incremental processing and memory release. Peak RSS drops from **28.8 GB → ~17 GB** for CAPAN-1, eliminating OOM on 32 GB machines.

Two key design decisions:

1. **Array-by-array scatter** — Python drives a loop that partitions one global array at a time, calling a C++ scatter function for each, then immediately frees the global array. Peak swap-over memory is bounded to `base + global_remaining + sizeof(largest_array)`.

2. **Unified partition-native EM** — A single new C++ function `batch_locus_em_partitioned()` replaces `batch_locus_em()`. It accepts per-locus partition data (as Python lists of numpy arrays), extracts data pointers under GIL, then uses the same OpenMP mega/normal scheduling internally. No re-concatenation of partitions.

---

## Current Architecture

```
FragmentBuffer (54.5M frags, 3.2 GB)
    ↓ fused_score_buffer (76s)
ScoredFragments — MONOLITHIC GLOBAL CSR (50.9M units, 15.3 GB)
    ↓ connected_components (4s, reads only offsets + t_indices = 1.5 GB)
Loci[14,298]  (references into global CSR via unit_indices)
    ↓ batch_locus_em (378s)
    │   For each locus (C++, GIL released):
    │     extract_locus_sub_problem → scattered reads from global CSR
    │     build_equiv_classes → SQUAREM → assign_posteriors
    │     LocusSubProblem is reused (vectors cleared, not freed)
    ↓
del ScoredFragments
```

### Why Peak RSS Is 28.8 GB

The monolithic ScoredFragments (15.3 GB) must persist *throughout* `batch_locus_em` because each locus references arbitrary, non-contiguous ranges within the global CSR. During mega-locus processing:

| Component | Memory |
|-----------|--------|
| Base (Python, index, models) | 2.9 GB |
| ScoredFragments (global CSR) | 15.3 GB |
| Mega-locus LocusSubProblem extraction | ~4.5 GB |
| Mega-locus EC data + SQUAREM scratch | ~2.5 GB |
| Allocator fragmentation / overhead | ~3.6 GB |
| **Peak** | **28.8 GB** |

### Data Flow Audit

**Connected components** touches only 2 of 16 arrays:
- `offsets` (int64, 0.41 GB) + `t_indices` (int32, 1.14 GB) = **1.5 GB**

**`extract_locus_sub_problem`** (C++) reads 12 global arrays:
- Per-candidate (×284M): `t_indices`, `log_liks`, `coverage_wts`, `tx_starts`, `tx_ends`, `count_cols`
- Per-unit (×50.9M): `offsets`, `is_spliced`, `gdna_log_liks`, `genomic_footprints`, `locus_t_indices`, `locus_count_cols`

**`compute_gdna_locus_gammas`** does NOT access ScoredFragments — uses only `Locus.merged_intervals` + calibration data.

**Post-EM arrays** (`frag_ids`, `frag_class`, `splice_type`) are used only for annotation output, NOT by the EM. They consume 1.8 GB and can be handled separately.

---

## Proposed Architecture

```
FragmentBuffer (54.5M frags, 3.2 GB)
    ↓ fused_score_buffer (76s) — UNCHANGED
ScoredFragments — MONOLITHIC GLOBAL CSR (15.3 GB) — same as today
    ↓ connected_components (4s) — CHANGED: also returns comp_u_offsets/comp_u_flat
Loci[14,298] + comp_u_offsets + comp_u_flat
    ↓ compute_gdna_locus_gammas — UNCHANGED (no ScoredFragments access)
    ↓ array-by-array scatter (NEW, ~8s)
    │   Python loop: for each array in ScoredFragments:
    │     C++ scatter → per-locus partition arrays
    │     delete global array, gc.collect()
    ↓
dict[locus_id → LocusPartition] (~10 GB total)
    ↓ batch_locus_em_partitioned (NEW unified C++ function)
    │   Under GIL: extract all partition pointers → PartitionView[]
    │   Release GIL:
    │     Phase 1: mega-loci (all threads on intra-locus E-step)
    │     Phase 2: normal loci (OpenMP work-stealing)
    │   Same extract → build_ec → SQUAREM → assign pipeline
    ↓
Python: free partition dict; gc.collect()
```

For mega-loci that benefit from per-locus memory recovery, the Python driver
calls `batch_locus_em_partitioned` once per mega-locus, frees the partition
via `dict.pop()`, then calls again for all remaining normal loci:

```
for mega_locus in mega_loci:
    batch_locus_em_partitioned([mega_partition], ...)
    del mega_partition; gc.collect()       # ← incremental memory recovery

batch_locus_em_partitioned(normal_partitions, ...)
del partitions; gc.collect()
```

This naturally collapses to a single code path — the same C++ function
handles 1 partition (mega) or 14,000 partitions (normal). The scheduling
logic inside the function adapts automatically based on `n_loci` and work
distribution.

### Memory Timeline (CAPAN-1)

```
Time →     Scoring      CC    Scatter            Mega-EM        Normal-EM

Global    ║ 15.3 GB ║──────╗
CSR       ║ (all    ║      │ shrinks array-by-array
          ║ arrays) ║      │ (← freed after each scatter)
                           ▼
                      ╔══════════╗
Per-locus             ║ grows to ║─────────╗
partitions            ║ ~10 GB   ║         │
                      ╚══════════╝         │
                                           ▼
                                   ╔════════════╗  ╔══════════╗
                                   ║ 1 mega     ║  ║ all      ║
                                   ║ partition   ║  ║ normal   ║
                                   ║ + EM work   ║  ║ partns   ║
                                   ╚════════════╝  ╚══════════╝
                                      pop+free        free all
                                     (~−4.2 GB)

RSS:     18.2 GB   ≤18.2 ≤17.5 GB     ~17 GB       ~10 GB     → 2.9 GB
```

### Peak Memory Accounting

The scatter processes one global array at a time, largest first. At each step,
the per-locus copy for that array is allocated, then the global array is freed:

| Step | Array | Size | Global Remaining | Partitions | Total (+ base 2.9 GB) |
|------|-------|------|------------------|------------|------------------------|
| 0 (plan) | build offsets | — | 15.3 GB | +0.41 GB | 18.6 GB |
| 1 | `log_liks` | 2.27 GB | 13.0 | 2.68 | **17.6 GB** ← scatter peak |
| 2 | `coverage_weights` | 2.27 GB | 10.7 | 4.95 | 16.4 |
| 3 | `t_indices` | 1.14 GB | 9.6 | 6.09 | 15.4 |
| 4 | `tx_starts` | 1.14 GB | 8.5 | 7.23 | 15.5 |
| 5 | `tx_ends` | 1.14 GB | 7.3 | 8.37 | 15.5 |
| 6 | `count_cols` | 0.28 GB | 7.0 | 8.65 | 15.5 |
| — | free `g_offsets` | — | 6.6 | 8.65 | 15.3 |
| 7 | `gdna_log_liks` | 0.41 GB | 6.2 | 9.06 | 15.0 |
| 8 | `genomic_footprints` | 0.20 GB | 6.0 | 9.26 | 15.0 |
| 9 | `is_spliced` | 0.05 GB | 5.9 | 9.31 | 15.0 |
| 10 | `locus_t_indices` | 0.20 GB | 5.7 | 9.51 | 15.0 |
| 11 | `locus_count_cols` | 0.05 GB | 5.7 | 9.56 | 15.0 |
| — | free remaining | — | 0 | 9.56 | 12.5 |

Scatter peak = **17.6 GB** (step 1, during `log_liks` swap-over).

**During mega-locus EM** (after scatter complete, 12.5 GB baseline):
- Mega-partition data: 4.2 GB (in the partition dict)
- LocusSubProblem extraction: ~4.5 GB (dedup, reindex, gDNA — smaller than raw partition)
- EC data + SQUAREM scratch: ~2.5 GB
- **Peak ≈ 17 GB** (extraction + EC transient, overlapping with partition data)

**After mega-locus freed** (partition popped from dict):
- Remaining partitions: 9.56 − 4.2 = 5.36 GB
- Normal-locus EM working set: ~1–2 GB
- **Peak ≈ 10 GB**

**Overall peak: ~17.6 GB** (during first scatter or mega-locus EM) — **39% reduction from 28.8 GB**.

---

## Phase 1A: Array-by-Array Scatter

### Design

The scatter is driven by a Python loop that processes one ScoredFragments
array at a time. For each array, a C++ function copies the relevant
elements into per-locus numpy arrays, and Python immediately frees the
global array. This ensures the peak swap-over memory never exceeds
`base + global_remaining + sizeof(largest_array)`.

The scatter requires two kinds of copy:

1. **Per-candidate arrays** (6 arrays: `t_indices`, `log_liks`, `count_cols`,
   `coverage_weights`, `tx_starts`, `tx_ends`) — elements are indexed via
   the CSR `offsets` array. For global unit `u`, candidates occupy
   `global_arr[g_offsets[u] : g_offsets[u+1]]`. The global `offsets` array
   must stay alive until all per-candidate arrays are scattered.

2. **Per-unit arrays** (5 arrays: `is_spliced`, `gdna_log_liks`,
   `genomic_footprints`, `locus_t_indices`, `locus_count_cols`) —
   one element per unit, simple gather by unit index.

Before scattering any data arrays, we must build the **per-locus offsets**
from `g_offsets` and the connected components assignment. These local offsets
define the CSR structure within each partition and are needed to allocate
the correct sizes for per-candidate arrays.

### Scatter Plan: `build_partition_offsets()`

A lightweight C++ function that reads `g_offsets` and the CC assignment to
produce per-locus offset arrays:

```cpp
nb::list build_partition_offsets(
    i64_1d g_offsets,          // [n_global_units + 1]
    i64_1d comp_u_offsets,     // [n_loci + 1]
    i32_1d comp_u_flat,        // unit indices per locus (from CC)
    int n_loci
);
// Returns: Python list of int64 numpy arrays, one per locus
//          partition_offsets[li] has shape [n_local_units + 1]
//          partition_offsets[li][0] = 0
//          partition_offsets[li][k+1] = partition_offsets[li][k]
//                                    + (g_offsets[u+1] - g_offsets[u])
//          where u = comp_u_flat[comp_u_offsets[li] + k]
```

This is O(total_units), costs ~0.4 GB for the output (same as global offsets),
and runs in <1s.

### C++ Scatter Functions (Templated)

Two scatter functions, templated on element type, handle per-candidate and
per-unit arrays respectively:

```cpp
template <typename T, typename NB_T>
nb::list scatter_candidates(
    NB_T global_arr,            // 1D array of type T
    i64_1d g_offsets,           // global CSR offsets (still alive)
    i64_1d comp_u_offsets,      // [n_loci + 1]
    i32_1d comp_u_flat,         // unit indices per locus
    int n_loci
);
// Returns: list of T[] numpy arrays, one per locus
// Each output array has length = partition_offsets[li][-1]

template <typename T, typename NB_T>
nb::list scatter_units(
    NB_T global_arr,            // 1D array of type T
    i64_1d comp_u_offsets,      // [n_loci + 1]
    i32_1d comp_u_flat,         // unit indices per locus
    int n_loci
);
// Returns: list of T[] numpy arrays, one per locus
// Each output array has length = n_local_units
```

**Implementation pattern** (same for both):
```
for each locus li:
    allocate output array (new T[n_elements])
    copy elements (memcpy for per-candidate, gather for per-unit)
    wrap in numpy ndarray with capsule ownership
    append to Python list
```

Each template is instantiated for 4 types (float64, int32, int64, uint8)
and exposed to Python as typed functions:
```cpp
m.def("scatter_candidates_f64", &scatter_candidates<double, f64_1d>, ...);
m.def("scatter_candidates_i32", &scatter_candidates<int32_t, i32_1d>, ...);
m.def("scatter_candidates_u8",  &scatter_candidates<uint8_t, u8_1d>, ...);
m.def("scatter_units_f64", &scatter_units<double, f64_1d>, ...);
m.def("scatter_units_i32", &scatter_units<int32_t, i32_1d>, ...);
m.def("scatter_units_u8",  &scatter_units<uint8_t, u8_1d>, ...);
```

GIL is released during the copy loops (all data is extracted as raw pointers
under GIL first). The capsule-based ownership pattern matches the existing
`connected_components` return style.

### Python Driver: `partition_and_free()`

```python
def partition_and_free(
    em_data: ScoredFragments,
    comp_u_offsets: np.ndarray,
    comp_u_flat: np.ndarray,
    n_loci: int,
) -> dict[int, LocusPartition]:
    """Scatter global CSR into per-locus partitions, freeing each global
    array immediately after scatter. Returns dict[locus_id → LocusPartition]."""

    # Step 0: Build per-locus offsets (needs g_offsets)
    offsets_list = build_partition_offsets(
        em_data.offsets, comp_u_offsets, comp_u_flat, n_loci
    )

    # Step 1–6: Scatter per-candidate arrays (largest first)
    # g_offsets must stay alive during this phase.
    CAND_ARRAYS = [
        ("log_liks",          scatter_candidates_f64),
        ("coverage_weights",  scatter_candidates_f64),
        ("t_indices",         scatter_candidates_i32),
        ("tx_starts",         scatter_candidates_i32),
        ("tx_ends",           scatter_candidates_i32),
        ("count_cols",        scatter_candidates_u8),
    ]
    cand_results = {}
    for attr, scatter_fn in CAND_ARRAYS:
        global_arr = getattr(em_data, attr)
        cand_results[attr] = scatter_fn(
            global_arr, em_data.offsets, comp_u_offsets, comp_u_flat, n_loci
        )
        setattr(em_data, attr, None)
        del global_arr
        gc.collect()

    # Free g_offsets (no longer needed for per-candidate scatter)
    em_data.offsets = None
    gc.collect()

    # Step 7–11: Scatter per-unit arrays
    UNIT_ARRAYS = [
        ("gdna_log_liks",     scatter_units_f64),
        ("genomic_footprints", scatter_units_i32),
        ("is_spliced",        scatter_units_u8),   # bool in Python, view as uint8
        ("locus_t_indices",   scatter_units_i32),
        ("locus_count_cols",  scatter_units_u8),
    ]
    unit_results = {}
    for attr, scatter_fn in UNIT_ARRAYS:
        global_arr = getattr(em_data, attr)
        # is_spliced is stored as bool; view as uint8 for typed scatter
        if global_arr.dtype == np.bool_:
            global_arr = global_arr.view(np.uint8)
        unit_results[attr] = scatter_fn(
            global_arr, comp_u_offsets, comp_u_flat, n_loci
        )
        setattr(em_data, attr, None)
        del global_arr
        gc.collect()

    # Assemble per-locus LocusPartition objects
    partitions = {}
    for li in range(n_loci):
        partitions[li] = LocusPartition(
            locus_id=li,
            n_units=len(offsets_list[li]) - 1,
            n_candidates=int(offsets_list[li][-1]),
            offsets=offsets_list[li],
            t_indices=cand_results["t_indices"][li],
            log_liks=cand_results["log_liks"][li],
            count_cols=cand_results["count_cols"][li],
            coverage_weights=cand_results["coverage_weights"][li],
            tx_starts=cand_results["tx_starts"][li],
            tx_ends=cand_results["tx_ends"][li],
            is_spliced=cand_results["is_spliced"][li],
            gdna_log_liks=unit_results["gdna_log_liks"][li],
            genomic_footprints=unit_results["genomic_footprints"][li],
            locus_t_indices=unit_results["locus_t_indices"][li],
            locus_count_cols=unit_results["locus_count_cols"][li],
        )

    return partitions
```

**Why `t_indices` stay in global transcript space:** Mapping to local
transcript indices during scatter would require O(n_transcripts) reverse-lookup
construction per locus, destroying the O(N) performance of the scatter step.
The global→local remapping is handled inside
`extract_locus_sub_problem_from_partition()`, exactly as the current
`extract_locus_sub_problem()` does.

### New Data Structure: `LocusPartition`

```python
@dataclass(slots=True)
class LocusPartition:
    """Per-locus CSR subset. All arrays are contiguous, 0-indexed."""
    locus_id: int
    n_units: int
    n_candidates: int

    # CSR structure
    offsets: np.ndarray          # int64[n_units + 1]

    # Per-candidate arrays (indexed by offsets)
    t_indices: np.ndarray        # int32[n_candidates] — GLOBAL transcript indices
    log_liks: np.ndarray         # float64[n_candidates]
    count_cols: np.ndarray       # uint8[n_candidates]
    coverage_weights: np.ndarray # float64[n_candidates]
    tx_starts: np.ndarray        # int32[n_candidates]
    tx_ends: np.ndarray          # int32[n_candidates]

    # Per-unit arrays
    is_spliced: np.ndarray       # uint8[n_units] (stored as bool, viewed as uint8 for C++)
    gdna_log_liks: np.ndarray    # float64[n_units]
    genomic_footprints: np.ndarray # int32[n_units]
    locus_t_indices: np.ndarray  # int32[n_units]
    locus_count_cols: np.ndarray # uint8[n_units]
```

No `PartitionedCSR` wrapper class — a plain `dict[int, LocusPartition]`
suffices. Popping entries for mega-locus freeing is native dict behavior.

---

## Phase 1B: Unified Partition-Native EM

### Why No Re-concatenation

The original plan (Option A) proposed concatenating normal-locus partitions
back into a mini-global CSR to reuse the existing `batch_locus_em()`. This
is a trap:

- It creates a **temporary memory spike**: ~5.7 GB for the concatenated CSR
  co-existing with the partition dict before freeing.
- It wastes CPU time on a throwaway data transformation.
- It introduces `concatenate_partitions()` — a function that exists only to
  bridge an API mismatch and would be deleted once the real fix ships.

Instead, we write **one new C++ function** that replaces `batch_locus_em()`,
accepting per-locus partition data directly. The same function handles both
mega-loci (called with `n_loci=1`) and normal loci (called with
`n_loci=14K`). The scheduling logic (mega vs normal, thread allocation)
adapts automatically based on work distribution.

### `extract_locus_sub_problem_from_partition()`

A new C++ function that replaces `extract_locus_sub_problem()`. The key
simplification: with partitioned data, units are numbered 0..n_units−1
and the partition's `offsets` array IS the local CSR — no `g_offsets[global_u]`
indirection needed.

```cpp
static void extract_locus_sub_problem_from_partition(
    LocusSubProblem& sub,
    const PartitionView& pv,       // ← replaces 12+ global CSR pointers
    double locus_gamma,
    int64_t gdna_span,
    const double*  all_unambig_row_sums,
    const int64_t* all_t_starts,
    const int64_t* all_t_ends,
    const int64_t* all_t_lengths,
    int32_t* local_map,
    int local_map_size);
```

**Changes from `extract_locus_sub_problem()`:**

| Aspect | Old (global CSR) | New (partition) |
|--------|-----------------|-----------------|
| Unit iteration | `u = u_arr[ui]` (scattered global indices) | `ui = 0..n_units-1` (contiguous) |
| Candidate range | `g_offsets[u]..g_offsets[u+1]` | `pv.offsets[ui]..pv.offsets[ui+1]` |
| Candidate access | `g_t_indices[j]`, `g_log_liks[j]` | `pv.t_indices[j]`, `pv.log_liks[j]` |
| Per-unit access | `g_is_spliced[u]` | `pv.is_spliced[ui]` |
| Cache behavior | Scattered reads across global arrays | Sequential reads within partition |

The dedup, global→local remapping, gDNA candidate addition, prior
construction, and bias profile setup are **identical** to the current code.
This is a mechanical refactor: replace `src[global_index]` with
`pv.field[local_index]`.

### C++ Data Bridge: `PartitionView`

```cpp
struct PartitionView {
    const int64_t* offsets;
    const int32_t* t_indices;
    const double*  log_liks;
    const double*  coverage_wts;
    const int32_t* tx_starts;
    const int32_t* tx_ends;
    const uint8_t* count_cols;
    const uint8_t* is_spliced;
    const double*  gdna_log_liks;
    const int32_t* genomic_footprints;
    const int32_t* locus_t_indices;
    const uint8_t* locus_count_cols;
    int n_units;
    int64_t n_candidates;

    const int32_t* transcript_indices;  // global transcript indices for this locus
    int n_transcripts;
};
```

Under GIL, the function iterates over the Python partition list and
extracts raw pointers into a `std::vector<PartitionView>`. The numpy
arrays are kept alive by the Python list (a function argument), so the
raw pointers remain valid throughout the GIL-released compute phase.

```cpp
// Under GIL (pseudocode):
std::vector<PartitionView> views(n_loci);
for (int i = 0; i < n_loci; ++i) {
    nb::object part = partitions[i];    // Python LocusPartition object
    auto off = nb::cast<i64_1d>(part.attr("offsets"));
    views[i].offsets = off.data();
    views[i].n_units = static_cast<int>(off.shape(0)) - 1;
    views[i].n_candidates = off.data()[views[i].n_units];
    // ... same pattern for all 12 fields
    auto ti = nb::cast<i32_1d>(locus_transcript_indices[i]);
    views[i].transcript_indices = ti.data();
    views[i].n_transcripts = static_cast<int>(ti.shape(0));
}
// 14,298 loci × 13 casts ≈ 186K Python attribute accesses (~50ms)
```

### `batch_locus_em_partitioned()` — Unified C++ Function

**Signature:**

```cpp
nb::tuple batch_locus_em_partitioned(
    // Per-locus partitions
    nb::list partitions,                    // list[LocusPartition]
    nb::list locus_transcript_indices,      // list[np.ndarray(int32)]
    f64_1d   locus_gammas,                  // [n_loci]
    i64_1d   gdna_spans,                    // [n_loci]
    // Per-transcript data (GLOBAL, shared across all loci)
    f64_2d   unambig_counts,                // [N_T, n_cols]
    i64_1d   t_starts,
    i64_1d   t_ends,
    i64_1d   t_lengths,
    // Mutable output accumulators (GLOBAL, disjoint writes per locus)
    f64_2d_mut em_counts_out,               // [N_T, n_cols]
    f64_2d_mut gdna_locus_counts_out,       // [N_T, n_cols]
    f64_1d_mut posterior_sum_out,            // [N_T]
    f64_1d_mut n_assigned_out,              // [N_T]
    // EM config (same params as batch_locus_em)
    int max_iterations, double convergence_delta,
    double total_pseudocount, bool use_vbem,
    int assignment_mode, double assignment_min_posterior,
    uint64_t rng_seed, int n_transcripts_total,
    int n_splice_strand_cols, int n_threads,
    bool emit_locus_stats
);
```

**Internal structure** (mirrors current `batch_locus_em` exactly):

```cpp
{
    // ---- Under GIL ----
    // 1. Extract PartitionView[] from Python objects
    // 2. Pre-compute per-transcript unambig row sums
    // 3. Allocate scratch, output vectors, locus profiling

    nb::gil_scoped_release release;

    // ---- GIL released ----
    auto process_locus = [&](int li, int estep_thr,
                             LocusSubProblem& sub,
                             std::vector<int32_t>& local_map,
                             bool mega,
                             EStepThreadPool* pool) -> double {
        // 1. extract_locus_sub_problem_from_partition(sub, views[li], ...)
        // 2. apply_bias_correction_uniform(...)
        // 3. build_equiv_classes(...)
        // 4. compute_ovr_prior_and_warm_start(...)
        // 5. run_squarem(...)
        // 6. assign_posteriors(...)
        // Identical logic to current process_locus lambda
    };

    // Phase 1: mega-loci (all threads on intra-locus E-step)
    // Phase 2: normal loci (OpenMP work-stealing, chunks of 16)
    // Same scheduling as current batch_locus_em
}
```

The only difference from the current `batch_locus_em` is step 1 of
`process_locus`: it calls `extract_locus_sub_problem_from_partition()`
reading from `views[li]` instead of `extract_locus_sub_problem()` reading
from global CSR pointers. **Everything downstream is identical.**

### Python EM Driver

```python
def _run_locus_em_partitioned(
    estimator, partitions, loci, index, locus_gammas, em_config,
    *, emit_locus_stats=False,
):
    """Run locus EM from partitioned data with incremental memory freeing."""
    n_threads = em_config.n_threads or os.cpu_count() or 1

    # Classify loci into mega vs normal (same heuristic as C++ scheduling)
    total_work = sum(
        len(l.transcript_indices) * len(l.unit_indices) for l in loci
    )
    fair_share = total_work // n_threads if n_threads > 1 else total_work + 1

    mega_loci = sorted(
        [l for l in loci
         if len(l.transcript_indices) * len(l.unit_indices) >= fair_share],
        key=lambda l: len(l.transcript_indices) * len(l.unit_indices),
        reverse=True,
    )
    mega_ids = {l.locus_id for l in mega_loci}

    # --- Phase A: Mega-loci one at a time (incremental memory freeing) ---
    for locus in mega_loci:
        partition = partitions.pop(locus.locus_id)
        _call_batch_em_partitioned(
            estimator, [partition], [locus], index,
            locus_gammas, em_config,
            n_threads=n_threads,
            emit_locus_stats=emit_locus_stats,
        )
        del partition
        gc.collect()

    # --- Phase B: All normal loci in one batched call ---
    normal_loci = [l for l in loci if l.locus_id not in mega_ids]
    if normal_loci:
        normal_partitions = [partitions[l.locus_id] for l in normal_loci]
        _call_batch_em_partitioned(
            estimator, normal_partitions, normal_loci, index,
            locus_gammas, em_config,
            n_threads=n_threads,
            emit_locus_stats=emit_locus_stats,
        )
    del partitions
    gc.collect()
```

The `_call_batch_em_partitioned()` helper prepares the arguments and calls
the C++ `batch_locus_em_partitioned()` function — the same code path
regardless of whether it's processing 1 mega-locus or 14,000 normal loci.

---

## Pipeline Integration

### Modified `quant_from_buffer()` Flow

```python
def quant_from_buffer(buffer, index, ...):
    # Phase 1–2: geometry + scoring — UNCHANGED
    geometry, estimator = _setup_geometry_and_estimator(...)
    em_data = _score_fragments(buffer, ...)  # → ScoredFragments (15.3 GB)

    # Phase 3: Locus construction — UNCHANGED
    loci = build_loci(em_data, index)        # CC reads offsets + t_indices only
    locus_gammas = _compute_priors(estimator, loci, index, calibration)

    # Phase 4 (NEW): Array-by-array scatter + incremental free
    # Preserve CC output arrays for scatter (comp_u_offsets, comp_u_flat
    # are already in the Locus objects as unit_indices ranges)
    partitions = partition_and_free(em_data, loci)
    # At this point em_data arrays are all None, global CSR is freed

    # Extract annotation arrays before destroying em_data shell
    annotation_arrays = AnnotationArrays(
        frag_ids=em_data.frag_ids,
        frag_class=em_data.frag_class,
        splice_type=em_data.splice_type,
    )
    del em_data; gc.collect()

    # Phase 5 (NEW): Streaming locus EM with incremental partition freeing
    _run_locus_em_partitioned(
        estimator, partitions, loci, index,
        locus_gammas, em_config,
        emit_locus_stats=emit_locus_stats,
    )

    # Phase 6: cleanup
    gc.collect()
```

### Backward Compatibility

The old `batch_locus_em()` C++ function and `_run_locus_em()` Python function
remain in the codebase during the transition. They can be removed once the
new path is fully validated and golden tests pass. This enables A/B testing
between old and new paths.

---

## Implementation Steps

### Step 1: `LocusPartition` Dataclass
- **File:** `src/rigel/scored_fragments.py`
- **What:** Add `LocusPartition` dataclass (as shown above)
- **Risk:** None (additive, no existing code changes)

### Step 1.5: Modify `build_loci()` to Return CC Arrays
- **File:** `src/rigel/locus.py`
- **What:** Change `build_loci()` to return
  `(loci, comp_u_offsets, comp_u_flat)` in addition to the loci list.
  Currently `comp_u_offsets`/`comp_u_flat` from
  `connected_components_native()` are consumed internally — each locus
  copies its slice into `Locus.unit_indices`. The scatter functions need
  the CSR-form assignment directly for efficient array-by-array scatter.
- **Callers:** Update `pipeline.py` to capture the new return values.
- **Risk:** None — additive, no logic changes

### Step 2: C++ Scatter Functions
- **File:** `src/rigel/native/em_solver.cpp`
- **What:**
  - `build_partition_offsets()` — from `g_offsets` + CC assignment → per-locus offset arrays
  - `scatter_candidates<T>()` — scatter per-candidate array by locus
  - `scatter_units<T>()` — scatter per-unit array by locus
  - Typed instantiations for float64, int32, int64, uint8
- **Bindings:** Expose via `_em_impl` module
- **Testing:** Round-trip test: scatter each array, then for each locus verify that gathering back using the original `Locus.unit_indices` produces the same data as indexing the original global array
- **Risk:** None — pure data reorganization, no algorithmic changes

### Step 3: `partition_and_free()` Python Function
- **File:** `src/rigel/locus.py` (or new `src/rigel/partition.py`)
- **What:** Python driver that calls `build_partition_offsets`, then
  iterates scatter functions array-by-array, freeing each global array
- **Testing:** Verify all global CSR arrays are freed (set to None) and
  per-locus data matches original global indexing
- **Risk:** None

### Step 4: `PartitionView` + `extract_locus_sub_problem_from_partition()` in C++
- **File:** `src/rigel/native/em_solver.cpp`
- **What:** `PartitionView` struct + extraction function that reads from
  contiguous per-locus arrays instead of scattered global CSR. Mechanical
  refactor of `extract_locus_sub_problem()`.
- **Testing:** For each locus, verify the `LocusSubProblem` produced by
  the new function is bit-for-bit identical to the old function. This is
  the critical correctness gate.
- **Risk:** Low — same algorithm, different data layout

### Step 5: `batch_locus_em_partitioned()` C++ Function
- **File:** `src/rigel/native/em_solver.cpp`
- **What:** Fork of `batch_locus_em()` that accepts Python list of
  `LocusPartition` objects, extracts `PartitionView[]` under GIL, then
  runs the same scheduling + `process_locus` pipeline with the new
  extraction function.
- **Bindings:** Expose via `_em_impl` module
- **Testing:** Verify bit-for-bit identical EM output counts vs old
  `batch_locus_em()` path for same input data
- **Risk:** Low — same EM logic, different data access pattern

### Step 6: Python Pipeline Integration
- **Files:** `src/rigel/pipeline.py`, `src/rigel/estimator.py`
- **What:**
  - Replace `_run_locus_em()` with `_run_locus_em_partitioned()`
  - Modify `quant_from_buffer()` flow (Phase 4–5 changes)
  - Update `run_batch_locus_em()` in estimator to call the new C++ function
- **Testing:** Golden output tests must produce identical results
- **Risk:** Low — same data, same algorithms, different memory lifecycle

### Step 7: Validation
- All existing tests pass (967 tests, bit-for-bit where possible)
- Golden output regression tests pass without `--update-golden`
- Memory profiling on CAPAN-1 to verify peak RSS reduction
- Wall-time benchmark to verify overhead is < 10s

---

## Invariants & Correctness Properties

1. **Data equivalence:** For every locus, the `LocusSubProblem` produced by
   `extract_locus_sub_problem_from_partition()` must be identical to the one
   produced by `extract_locus_sub_problem()` given the same input data. This
   is the single most important invariant and must be verified per-locus in
   unit tests.

2. **Locus isolation:** Each transcript and each unit belongs to exactly one
   locus (guaranteed by connected components). Partitioning is lossless —
   no data is shared between loci, and the sum of all per-locus arrays
   equals the global array.

3. **Annotation preservation:** `frag_ids`, `frag_class`, `splice_type` are
   preserved separately (not partitioned) and remain accessible for post-EM
   annotation output.

4. **Thread safety:** Identical to current `batch_locus_em`. Different loci
   write to disjoint transcript indices in the output accumulators. Shared
   `total_gdna_em` uses `std::atomic`. LocusSubProblem and local_map are
   thread-private.

5. **Determinism:** Results are identical regardless of processing order.
   Per-locus RNG seed = `rng_seed ^ (locus_id * golden_ratio_constant)`,
   same as existing.

6. **Memory monotonicity during scatter:** After each scatter+free step,
   total memory (global remaining + partition accumulated) is constant
   (≈10 GB for CAPAN-1). The swap-over peak during any single step is
   bounded by `total + sizeof(current_array)`.

---

## Complexity & Performance

| Operation | Time Complexity | Expected Wall Time (CAPAN-1) |
|-----------|----------------|------------------------------|
| `build_partition_offsets` | O(total_units) | <1s |
| 6 × `scatter_candidates` + free | O(total_candidates) each | ~5s total |
| 5 × `scatter_units` + free | O(total_units) each | ~2s total |
| LocusPartition assembly | O(n_loci) | <1s |
| `batch_locus_em_partitioned` (mega) | Same as today | ~310s |
| `batch_locus_em_partitioned` (normal) | Same as today | ~68s |
| **Total overhead** | | **~8s** |

The 8s overhead is 1% of total wall time (757s). The scatter step may even
*improve* mega-locus EM performance due to better cache locality — the
partition data is contiguous vs scattered reads across a 15 GB global array.

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Scatter peak exceeds 32 GB RAM | **Impossible** — peak is 17.6 GB by construction (array-by-array guarantees monotonic swap) | — | Architecture prevents this |
| GIL-held pointer extraction too slow | Very low (~50ms for 186K attr reads) | Low | Profile; if needed, pass parallel lists instead of objects |
| `extract_from_partition` produces different LocusSubProblem | Low | **Critical** | Per-locus equivalence test comparing old vs new extraction |
| Wall-time regression | Very low (~8s on 757s total) | Low | Benchmark before/after |
| Python GC doesn't free arrays immediately | Low | Medium | Explicit `gc.collect()` after each `del`; verify with RSS tracking |

---

## Summary

This plan eliminates the monolithic ScoredFragments bottleneck through two
clean, composable mechanisms:

1. **Array-by-array scatter**: A Python-driven loop calls typed C++ scatter
   functions one global array at a time, immediately freeing each. The
   swap-over peak is bounded by construction to `base + 10 GB + 2.27 GB`
   = **17.6 GB**, with no reliance on "brief peaks" or GC timing.

2. **Unified partition-native EM**: A single `batch_locus_em_partitioned()`
   C++ function replaces the old `batch_locus_em()`. It accepts per-locus
   data via Python lists, extracts raw pointers under GIL, and runs the
   same scheduling pipeline with a simplified extraction function that
   benefits from contiguous, cache-friendly partition data.

There is no re-concatenation, no throwaway bridge functions, and no API
mismatch to work around. The same C++ function handles both mega-loci
(called with n=1 for incremental freeing) and normal loci (called with
n=14K for OpenMP batching).

Expected outcome: **peak RSS drops from 28.8 GB to ~17.6 GB** (39%
reduction), with **~8s** wall-time overhead and **zero** risk of
intermittent OOM on 32 GB machines.
