# Phase 1: Memory Architecture — Implementation Plan

## Executive Summary

Replace the monolithic `ScoredFragments` global CSR with a **locus-partitioned
CSR** that allows incremental processing and memory release.  Peak RSS drops
from **28.8 GB → ~20 GB** for CAPAN-1, eliminating OOM on 32 GB machines.

Three key design decisions:

1. **Array-by-array scatter** — Python drives a loop that partitions one
   global array at a time via a typed C++ scatter function, then immediately
   frees the global array.  Peak swap-over memory is bounded by construction.

2. **Unified partition-native EM** — A single new C++ function
   `batch_locus_em_partitioned()` replaces `batch_locus_em()`.  It accepts
   per-locus partition data as a list of tuples (matching the existing
   `fused_score_buffer` convention), extracts raw pointers under GIL, then
   runs the same OpenMP mega/normal scheduling internally.  No
   re-concatenation.

3. **Zero scaffold debt** — No `build_loci()` modification, no
   `AnnotationArrays` class, no `concatenate_partitions()` function.  Every
   new artifact exists permanently in the final architecture.

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

The monolithic ScoredFragments (15.3 GB) must persist *throughout*
`batch_locus_em` because each locus references arbitrary, non-contiguous
ranges within the global CSR.  During mega-locus processing:

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
- Per-candidate (×284M): `t_indices`, `log_liks`, `coverage_wts`,
  `tx_starts`, `tx_ends`, `count_cols`
- Per-unit (×50.9M): `offsets`, `is_spliced`, `gdna_log_liks`,
  `genomic_footprints`, `locus_t_indices`, `locus_count_cols`

**`compute_gdna_locus_gammas`** does NOT access ScoredFragments — uses
only `Locus.merged_intervals` + calibration data.

**Dead-weight arrays** (`frag_ids` int64, `frag_class` int8,
`splice_type` uint8 — total 0.51 GB) are never read after
`ScoredFragments` construction.  The `AnnotationTable` populated during
scanning has its own independent copies.  These arrays are pure waste and
can be freed immediately.

**Dead parameters in `extract_locus_sub_problem`:** The C++ function
accepts `all_t_starts` and `all_t_ends` pointers but never dereferences
them.  Only `all_t_lengths` is used (for `bias_profiles`).  The new
extraction function drops these dead parameters.

---

## Proposed Architecture

```
FragmentBuffer (54.5M frags, 3.2 GB)
    ↓ fused_score_buffer (76s) — UNCHANGED
ScoredFragments — MONOLITHIC GLOBAL CSR (15.3 GB) — same as today
    ↓ connected_components (4s) — UNCHANGED
Loci[14,298]  (each holds unit_indices, transcript_indices, etc.)
    ↓ compute_gdna_locus_gammas — UNCHANGED (no ScoredFragments access)
    ↓ free dead-weight arrays (−0.51 GB)
    ↓ array-by-array scatter (NEW, ~8s)
    │   Python loop: for each EM array in ScoredFragments:
    │     C++ scatter_T(global_arr, g_offsets, locus_units) → per-locus arrays
    │     setattr(em_data, attr, None); gc.collect()
    ↓
dict[locus_id → LocusPartition] (~10 GB total)
    ↓ batch_locus_em_partitioned (NEW unified C++ function)
    │   Under GIL: extract pointers from list of tuples → PartitionView[]
    │   Release GIL:
    │     Phase 1: mega-loci (all threads collaborate on intra-locus E-step)
    │     Phase 2: normal loci (work-stealing, chunks of 16)
    │   Same extract → build_ec → SQUAREM → assign pipeline
    ↓
Python: free partition dict; gc.collect()
```

For mega-loci, Python calls `batch_locus_em_partitioned` once per
mega-locus, `pop()`s the partition, then calls once more for all remaining
normal loci:

```python
for mega_locus in mega_loci:
    part = partitions.pop(mega_locus.locus_id)
    batch_em_partitioned([part_tuple], [mega_locus], ...)
    del part; gc.collect()    # ← incremental recovery

batch_em_partitioned(normal_part_tuples, normal_loci, ...)
del partitions; gc.collect()
```

The same C++ function handles 1 partition (mega) or 14,000 partitions
(normal).  Scheduling adapts internally.

### Memory Timeline (CAPAN-1)

```
Time →   Scoring    CC   Dead-free  Scatter       Mega-EM     Normal-EM

Global  ║ 15.3 GB ║────╗
CSR     ║         ║    │ −0.51
        ║         ║    │  (dead wt)
                       │  shrinks array-by-array
                       ▼
                  ╔══════════╗
Per-locus         ║ grows to ║────────╗
partitions        ║ ~10 GB   ║        │
                  ╚══════════╝        ▼
                                ╔══════════╗  ╔══════════╗
                                ║ 1 mega   ║  ║ all      ║
                                ║ + EM     ║  ║ normal   ║
                                ╚══════════╝  ╚══════════╝
                                  pop+free       free all
                                 (−4.2 GB)

RSS:   18.2 GB  18.2  17.7     ≤20 GB      ≤20 GB    ~10 GB  → 2.9 GB
```

### Peak Memory Accounting

**Step 0: Free dead-weight arrays**

Before any scatter, free `frag_ids` (0.41 GB), `frag_class` (0.05 GB),
`splice_type` (0.05 GB).  These are never read and their AnnotationTable
copies live independently.

| Array | Type | Count | Size |
|-------|------|-------|------|
| `frag_ids` | int64 | 50.9M | 0.41 GB |
| `frag_class` | int8 | 50.9M | 0.05 GB |
| `splice_type` | uint8 | 50.9M | 0.05 GB |
| **Dead-weight total** | | | **0.51 GB** |

Post-free global CSR: 15.3 − 0.51 = **14.8 GB** (RSS ≈ 17.7 GB with base).

**Array-by-array scatter peak**

Each scatter step: C++ allocates per-locus arrays (total = array size),
returns to Python, Python frees the global array.  During the C++ call,
both the global and per-locus copies coexist:

```
peak_during_step = RSS_before_step + sizeof(array_being_scattered)
```

Since each swap is size-neutral (free as much as allocate), RSS before
each step is ≈ constant at 17.7 GB.  The peak occurs during the largest
array swap:

```
Scatter peak = 17.7 + 2.27 (log_liks) ≈ 20.0 GB
```

| Step | Array | Size | Notes |
|------|-------|------|-------|
| 0 | `build_partition_offsets` | 0.41 GB | from `g_offsets` |
| 1 | `log_liks` | 2.27 GB | **scatter peak ≈ 20.0 GB** |
| 2 | `coverage_weights` | 2.27 GB | peak ≈ 20.0 GB |
| 3 | `t_indices` | 1.14 GB | |
| 4 | `tx_starts` | 1.14 GB | |
| 5 | `tx_ends` | 1.14 GB | |
| 6 | `count_cols` | 0.28 GB | |
| — | free `g_offsets` | −0.41 GB | no longer needed |
| 7–11 | per-unit arrays | 0.91 GB total | |

After all scatter + global freed: base (2.9) + partitions (10.1)
= **13.0 GB**.

**During mega-locus EM** (after scatter complete):

| Component | Memory |
|-----------|--------|
| Partitions (all) | 10.1 GB |
| EM extract + EC + SQUAREM for mega-locus | ~7.0 GB |
| Base | 2.9 GB |
| **Peak** | **~20 GB** |

**After mega-locus freed** (partition popped):

| Component | Memory |
|-----------|--------|
| Remaining partitions | ~5.9 GB |
| EM working set (normal loci) | ~2 GB |
| Base | 2.9 GB |
| **Peak** | **~11 GB** |

**Overall peak: ~20 GB** (during scatter of largest array or during
mega-locus EM) — **31% reduction from 28.8 GB**.

---

## Phase 1A: Array-by-Array Scatter

### Design

The scatter is driven by a Python loop that processes one ScoredFragments
array at a time.  For each array, a typed C++ scatter function copies the
relevant elements into per-locus numpy arrays.  Python immediately frees
the global array.

Two kinds of scatter:

1. **Per-candidate arrays** (6 arrays: `t_indices`, `log_liks`, `count_cols`,
   `coverage_weights`, `tx_starts`, `tx_ends`) — elements are indexed via
   the CSR `offsets` array.  For global unit `u`, candidates occupy
   `global_arr[g_offsets[u] : g_offsets[u+1]]`.  The global `offsets` array
   must stay alive until all per-candidate arrays are scattered.

2. **Per-unit arrays** (5 arrays: `is_spliced`, `gdna_log_liks`,
   `genomic_footprints`, `locus_t_indices`, `locus_count_cols`) —
   one element per unit, simple gather by unit index.

Before scattering data arrays, we compute **per-locus offsets** from
`g_offsets` and the locus unit assignments.  These define the CSR structure
within each partition and are needed to allocate per-candidate arrays.

### Key Design Choice: Use `Locus.unit_indices` Directly

Each `Locus` object already holds `unit_indices` — a contiguous int32 array
of global unit indices belonging to that locus.  These are copies of slices
from the connected-components output, already allocated.  We pass them
directly to the scatter functions as a `nb::list` of `i32_1d` arrays:

```python
locus_units = [locus.unit_indices for locus in loci]
```

This eliminates any need to modify `build_loci()` or reconstruct CSR arrays
from the connected-components output.  The list is tiny (~0.2 GB, already
allocated in the Locus objects).

### `build_partition_offsets()` C++ Function

```cpp
nb::list build_partition_offsets(
    i64_1d g_offsets,       // [n_global_units + 1]
    nb::list locus_units,   // list of int32[] arrays (one per locus)
    int n_loci
);
// Returns: Python list of int64 numpy arrays, one per locus
//   partition_offsets[li][0] = 0
//   partition_offsets[li][k+1] = partition_offsets[li][k]
//       + (g_offsets[u+1] - g_offsets[u])
//   where u = locus_units[li][k]
```

O(total_units).  Produces ~0.4 GB of per-locus offset arrays.  < 1s.

### Typed C++ Scatter Functions

Two function templates, instantiated for each element type:

```cpp
template <typename T>
nb::list scatter_candidates(
    nb::ndarray<const T, nb::ndim<1>, nb::c_contig> global_arr,
    i64_1d g_offsets,
    nb::list locus_units,       // list of int32[]
    nb::list partition_offsets,  // from build_partition_offsets
    int n_loci
);
// For each locus, copies candidate data into a contiguous array.
// Inner loop per unit: memcpy(dst + p_off[ui], src + g_off[u], n * sizeof(T))

template <typename T>
nb::list scatter_units(
    nb::ndarray<const T, nb::ndim<1>, nb::c_contig> global_arr,
    nb::list locus_units,
    int n_loci
);
// For each locus, gathers per-unit elements: dst[k] = src[u_arr[k]]
```

Implementation pattern (same for both):
```
Under GIL: extract all raw pointers (g_offsets, each locus_units[li], etc.)
Release GIL:
  for each locus li (sequential — I/O bound, parallelism not worthwhile):
    allocate output array (new T[n_elements])
    copy elements:
      scatter_candidates: memcpy per-unit candidate segments
      scatter_units: element-wise gather
Under GIL: wrap each in numpy ndarray with capsule ownership, append to list
```

Exposed to Python as typed functions:
```cpp
m.def("scatter_candidates_f64", &scatter_candidates<double>, ...);
m.def("scatter_candidates_i32", &scatter_candidates<int32_t>, ...);
m.def("scatter_candidates_u8",  &scatter_candidates<uint8_t>, ...);
m.def("scatter_units_f64",      &scatter_units<double>, ...);
m.def("scatter_units_i32",      &scatter_units<int32_t>, ...);
m.def("scatter_units_u8",       &scatter_units<uint8_t>, ...);
```

The `scatter_candidates` functions take `partition_offsets` (output of
`build_partition_offsets`) so that each call knows allocation sizes and
write positions without recomputing from `g_offsets`.

### Python Driver: `partition_and_free()`

```python
def partition_and_free(
    em_data: ScoredFragments,
    loci: list[Locus],
) -> dict[int, LocusPartition]:
    """Scatter global CSR into per-locus partitions, freeing each global
    array immediately after scatter."""
    n_loci = len(loci)
    locus_units = [locus.unit_indices for locus in loci]

    # ---- Free dead-weight arrays (never used after construction) ----
    em_data.frag_ids = None
    em_data.frag_class = None
    em_data.splice_type = None
    gc.collect()

    # ---- Build per-locus CSR offsets ----
    offsets_list = build_partition_offsets(
        em_data.offsets, locus_units, n_loci
    )

    # ---- Scatter per-candidate arrays (largest first) ----
    # g_offsets must stay alive during this phase.
    CAND_ARRAYS = [
        ("log_liks",         scatter_candidates_f64),
        ("coverage_weights", scatter_candidates_f64),
        ("t_indices",        scatter_candidates_i32),
        ("tx_starts",        scatter_candidates_i32),
        ("tx_ends",          scatter_candidates_i32),
        ("count_cols",       scatter_candidates_u8),
    ]
    cand_results = {}
    for attr, scatter_fn in CAND_ARRAYS:
        global_arr = getattr(em_data, attr)
        cand_results[attr] = scatter_fn(
            global_arr, em_data.offsets,
            locus_units, offsets_list, n_loci,
        )
        setattr(em_data, attr, None)
        del global_arr
        gc.collect()

    # g_offsets no longer needed
    em_data.offsets = None
    gc.collect()

    # ---- Scatter per-unit arrays ----
    UNIT_ARRAYS = [
        ("gdna_log_liks",     scatter_units_f64),
        ("genomic_footprints", scatter_units_i32),
        ("locus_t_indices",   scatter_units_i32),
        ("locus_count_cols",  scatter_units_u8),
        ("is_spliced",        scatter_units_u8),
    ]
    unit_results = {}
    for attr, scatter_fn in UNIT_ARRAYS:
        global_arr = getattr(em_data, attr)
        if global_arr.dtype == np.bool_:
            global_arr = global_arr.view(np.uint8)
        unit_results[attr] = scatter_fn(global_arr, locus_units, n_loci)
        setattr(em_data, attr, None)
        del global_arr
        gc.collect()

    # ---- Assemble LocusPartition objects ----
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
            is_spliced=unit_results["is_spliced"][li],
            gdna_log_liks=unit_results["gdna_log_liks"][li],
            genomic_footprints=unit_results["genomic_footprints"][li],
            locus_t_indices=unit_results["locus_t_indices"][li],
            locus_count_cols=unit_results["locus_count_cols"][li],
        )

    return partitions
```

**Why `t_indices` stay in global transcript space:** Mapping to local
transcript indices during scatter would require building a reverse-lookup
table per locus, destroying the O(N) scatter performance.  The
global→local remapping is deferred to `extract_locus_sub_problem_from_partition()`,
exactly as the current code does.

### `LocusPartition` Dataclass

```python
@dataclass(slots=True)
class LocusPartition:
    """Per-locus CSR subset.  All arrays contiguous, 0-indexed."""
    locus_id: int
    n_units: int
    n_candidates: int

    # CSR structure
    offsets: np.ndarray          # int64[n_units + 1]

    # Per-candidate arrays (indexed by offsets)
    t_indices: np.ndarray        # int32  — GLOBAL transcript indices
    log_liks: np.ndarray         # float64
    count_cols: np.ndarray       # uint8
    coverage_weights: np.ndarray # float64
    tx_starts: np.ndarray        # int32
    tx_ends: np.ndarray          # int32

    # Per-unit arrays
    is_spliced: np.ndarray       # uint8 (bool viewed as uint8 for C++)
    gdna_log_liks: np.ndarray    # float64
    genomic_footprints: np.ndarray # int32
    locus_t_indices: np.ndarray  # int32
    locus_count_cols: np.ndarray # uint8
```

No wrapper class — a plain `dict[int, LocusPartition]` suffices.  Popping
entries for mega-locus freeing is native dict behavior.

---

## Phase 1B: Unified Partition-Native EM

### Why No Re-concatenation

Re-concatenating normal-locus partitions back into a mini-global CSR to
reuse the existing `batch_locus_em()`:

- Creates a pointless memory spike (~6 GB during re-concat).
- Wastes CPU on a throwaway transformation.
- Introduces `concatenate_partitions()` — a function that exists only to
  bridge an API mismatch and would be deleted later.

Instead: **one new C++ function** that accepts per-locus partition data
directly.  Same function handles both mega-loci (n=1) and normal loci
(n=14K).

### Tuple-of-Arrays Convention for C++ Interface

Following the established `fused_score_buffer` pattern in this codebase
(which already passes chunk data as a `nb::list` of tuples), the new
function receives partition data as a **list of tuples**, one per locus:

```python
# Python helper builds the tuple list once:
partition_tuples = []
for partition in partitions:
    partition_tuples.append((
        partition.offsets,          # [0] int64
        partition.t_indices,        # [1] int32
        partition.log_liks,         # [2] float64
        partition.coverage_weights, # [3] float64
        partition.tx_starts,        # [4] int32
        partition.tx_ends,          # [5] int32
        partition.count_cols,       # [6] uint8
        partition.is_spliced,       # [7] uint8
        partition.gdna_log_liks,    # [8] float64
        partition.genomic_footprints, # [9] int32
        partition.locus_t_indices,  # [10] int32
        partition.locus_count_cols, # [11] uint8
    ))
```

C++ unpacks under GIL:
```cpp
for (int i = 0; i < n_loci; ++i) {
    nb::tuple tup = nb::borrow<nb::tuple>(partition_tuples[i]);
    views[i].offsets          = nb::cast<i64_1d>(tup[0]).data();
    views[i].t_indices        = nb::cast<i32_1d>(tup[1]).data();
    views[i].log_liks         = nb::cast<f64_1d>(tup[2]).data();
    views[i].coverage_wts     = nb::cast<f64_1d>(tup[3]).data();
    views[i].tx_starts        = nb::cast<i32_1d>(tup[4]).data();
    views[i].tx_ends          = nb::cast<i32_1d>(tup[5]).data();
    views[i].count_cols       = nb::cast<u8_1d>(tup[6]).data();
    views[i].is_spliced       = nb::cast<u8_1d>(tup[7]).data();
    views[i].gdna_log_liks    = nb::cast<f64_1d>(tup[8]).data();
    views[i].genomic_footprints = nb::cast<i32_1d>(tup[9]).data();
    views[i].locus_t_indices  = nb::cast<i32_1d>(tup[10]).data();
    views[i].locus_count_cols = nb::cast<u8_1d>(tup[11]).data();
    views[i].n_units = static_cast<int>(nb::cast<i64_1d>(tup[0]).shape(0)) - 1;
    views[i].n_candidates = views[i].offsets[views[i].n_units];
}
// 14,298 tuples × 12 casts ≈ 172K nanobind casts → ~40ms
```

This is cleaner than passing 12 parallel lists (verbose), cleaner than
Python attribute access (no string lookups in C++), and matches a proven
pattern already in the codebase.  The tuple positional order is
documented by the Python helper function.

### `PartitionView` (C++ struct)

```cpp
struct PartitionView {
    // Per-locus CSR data (contiguous, 0-indexed)
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
    int     n_units;
    int64_t n_candidates;

    // Locus transcript membership (from Locus.transcript_indices)
    const int32_t* transcript_indices;
    int n_transcripts;
};
```

### `extract_locus_sub_problem_from_partition()`

Replaces `extract_locus_sub_problem()`.  Key simplifications:

| Aspect | Old (global CSR) | New (partition) |
|--------|-----------------|-----------------|
| **Signature** | 20+ parameters | **8 parameters** (PartitionView replaces 12 CSR ptrs) |
| Unit iteration | `u = u_arr[ui]` (scattered) | `ui = 0..n_units−1` (sequential) |
| Candidate range | `g_offsets[u]..g_offsets[u+1]` | `pv.offsets[ui]..pv.offsets[ui+1]` |
| Per-unit access | `g_is_spliced[u]` | `pv.is_spliced[ui]` |
| Dead params | `all_t_starts`, `all_t_ends` (never used) | **Dropped** |
| **Cache behavior** | Scattered reads across 15 GB global | **Sequential reads in ~4 GB partition** |

```cpp
static void extract_locus_sub_problem_from_partition(
    LocusSubProblem& sub,
    const PartitionView& pv,       // replaces 12 global CSR pointers
    double locus_gamma,
    int64_t gdna_span,
    const double*  all_unambig_row_sums,
    const int64_t* all_t_lengths,  // only t_lengths needed (not t_starts/t_ends)
    int32_t* local_map,
    int local_map_size)
{
    int n_t = pv.n_transcripts;
    int n_u = pv.n_units;
    const int32_t* t_arr = pv.transcript_indices;

    sub.n_t = n_t;
    sub.n_local_units = n_u;
    sub.n_components = n_t + 1;
    sub.gdna_idx = n_t;
    int nc = sub.n_components;

    // --- Build global→local mapping ---
    // (identical logic to existing function)

    // --- Gather + dedup candidates from partition CSR ---
    for (int ui = 0; ui < n_u; ++ui) {
        // Partition offsets are directly usable — no indirection
        auto p_start = pv.offsets[ui];
        auto p_end   = pv.offsets[ui + 1];

        sub.locus_t_arr[ui] = pv.locus_t_indices[ui];
        sub.locus_ct_arr[ui] = pv.locus_count_cols[ui];

        // Gather and dedup RNA candidates
        for (auto j = p_start; j < p_end; ++j) {
            int32_t global_t = pv.t_indices[j];
            // ... same dedup logic as existing code ...
            // but reads pv.log_liks[j], pv.coverage_wts[j], etc.
            // instead of g_log_liks[j], g_coverage_wts[j]
        }

        // Add gDNA candidate for unspliced units
        bool is_spliced = (pv.is_spliced[ui] != 0);
        double gdna_ll = pv.gdna_log_liks[ui];
        // ... identical gDNA logic ...
    }

    // --- Per-component arrays (identical) ---
    // unambig_totals, bias_profiles, prior, eligible — same logic
    // Clean up local_map scratch — same
}
```

The dedup, global→local remapping, gDNA candidate addition, prior
construction, and bias profile setup are **mechanically identical** to
the existing code.  Only the data *source* changes:
`pv.field[local_index]` instead of `g_field[global_index]`.

**Cache locality bonus:** With the old code, accessing
`g_log_liks[g_offsets[u]]` for scattered global unit `u` causes cache
misses when the mega-locus's 21.7M units are distributed across a
15 GB array.  With partition data, the mega-locus's candidates are
contiguous — sequential reads with excellent prefetcher behavior.
This may produce a measurable speedup in the extraction phase beyond
the memory savings.

### `batch_locus_em_partitioned()` — C++ Function

```cpp
static std::tuple<double, nb::ndarray<...>, nb::ndarray<...>, nb::list>
batch_locus_em_partitioned(
    // Per-locus partition data (list of 12-tuples)
    nb::list partition_tuples,
    // Per-locus transcript membership (list of int32[])
    nb::list locus_transcript_indices,
    // Per-locus scalars
    f64_1d   locus_gammas,           // [n_loci]
    i64_1d   gdna_spans,             // [n_loci]
    // Per-transcript globals (shared across all loci)
    f64_2d   unambig_counts,         // [N_T, n_cols]
    i64_1d   t_lengths,              // [N_T]
    // Mutable output accumulators (disjoint writes per locus)
    f64_2d_mut em_counts_out,        // [N_T, n_cols]
    f64_2d_mut gdna_locus_counts_out,// [N_T, n_cols]
    f64_1d_mut posterior_sum_out,     // [N_T]
    f64_1d_mut n_assigned_out,       // [N_T]
    // EM config
    int max_iterations, double convergence_delta,
    double total_pseudocount, bool use_vbem,
    int assignment_mode, double assignment_min_posterior,
    uint64_t rng_seed, int n_transcripts_total,
    int n_splice_strand_cols, int n_threads,
    bool emit_locus_stats
)
```

**vs. current `batch_locus_em` signature:**

| Parameter group | Old | New |
|----------------|-----|-----|
| CSR data | 12 global arrays | `nb::list partition_tuples` |
| Locus definitions | 4 CSR arrays | `nb::list locus_transcript_indices` |
| Locus scalars | `locus_gammas`, `gdna_spans` | Same |
| Per-transcript | `unambig_counts`, `t_starts`, `t_ends`, `t_lengths` | `unambig_counts`, `t_lengths` (drops dead `t_starts`/`t_ends`) |
| Output | 4 mutable accumulators | Same |
| Config | 10 params | Same |

**Internal structure** (mirrors current `batch_locus_em`):

```cpp
{
    // ---- Under GIL ----
    // 1. Extract PartitionView[] from partition_tuples + locus_transcript_indices
    // 2. Pre-compute per-transcript unambig row sums
    // 3. Allocate scratch, output vectors

    nb::gil_scoped_release release;

    // ---- GIL released ----
    // Identical scheduling to current batch_locus_em:
    //   - Compute approximate work per locus (n_t × n_u)
    //   - Sort by descending work
    //   - Phase 1: mega-loci (all threads collaborate)
    //   - Phase 2: normal loci (work-stealing, chunks of 16)
    //
    // The only difference is process_locus step 1:
    //   extract_locus_sub_problem_from_partition(sub, views[li], ...)
    // instead of:
    //   extract_locus_sub_problem(sub, t_arr, n_t, u_arr, n_u, ..., goff, gti, ...)
}
```

### Python EM Driver

```python
def _run_locus_em_partitioned(
    estimator, partitions, loci, index,
    locus_gammas, em_config, *, emit_locus_stats=False,
):
    """Run locus EM from partitioned data with incremental memory freeing."""
    n_threads = em_config.n_threads or os.cpu_count() or 1

    # Classify mega vs normal (same heuristic as C++ scheduler)
    total_work = sum(
        len(l.transcript_indices) * p.n_units
        for l, p in ((l, partitions[l.locus_id]) for l in loci)
    )
    fair_share = total_work // n_threads if n_threads > 1 else total_work + 1

    mega_loci = sorted(
        [l for l in loci
         if len(l.transcript_indices) * partitions[l.locus_id].n_units >= fair_share],
        key=lambda l: len(l.transcript_indices) * partitions[l.locus_id].n_units,
        reverse=True,
    )
    mega_ids = {l.locus_id for l in mega_loci}

    # ---- Phase A: Mega-loci (one at a time, free after each) ----
    for locus in mega_loci:
        part = partitions.pop(locus.locus_id)
        _call_batch_em_partitioned(
            estimator, [part], [locus], index,
            locus_gammas, em_config,
            n_threads=n_threads,
            emit_locus_stats=emit_locus_stats,
        )
        del part
        gc.collect()

    # ---- Phase B: Normal loci (one batched call) ----
    normal_loci = [l for l in loci if l.locus_id not in mega_ids]
    if normal_loci:
        normal_parts = [partitions[l.locus_id] for l in normal_loci]
        _call_batch_em_partitioned(
            estimator, normal_parts, normal_loci, index,
            locus_gammas, em_config,
            n_threads=n_threads,
            emit_locus_stats=emit_locus_stats,
        )

    del partitions
    gc.collect()


def _call_batch_em_partitioned(
    estimator, partitions, loci, index,
    locus_gammas, em_config, *, n_threads, emit_locus_stats,
):
    """Thin wrapper: pack tuples, call C++, record results."""
    # Pack partition data as tuples (matching fused_score_buffer convention)
    partition_tuples = [
        (p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
         p.tx_starts, p.tx_ends, p.count_cols,
         p.is_spliced, p.gdna_log_liks, p.genomic_footprints,
         p.locus_t_indices, p.locus_count_cols)
        for p in partitions
    ]
    locus_t_lists = [l.transcript_indices for l in loci]
    gammas = np.array([locus_gammas[l.locus_id] for l in loci], dtype=np.float64)
    spans = np.array([l.gdna_span for l in loci], dtype=np.int64)

    total_gdna, locus_mrna, locus_gdna, stats = _batch_locus_em_partitioned(
        partition_tuples,
        locus_t_lists,
        gammas,
        spans,
        estimator.unambig_counts,
        index.t_df["length"].values.astype(np.int64),
        estimator.em_counts,
        estimator.gdna_locus_counts,
        estimator._em_posterior_sum,
        estimator._em_n_assigned,
        em_config.iterations,
        em_config.convergence_delta,
        em_config.prior_pseudocount,
        em_config.mode == "vbem",
        _ASSIGNMENT_MODE_MAP[em_config.assignment_mode],
        em_config.assignment_min_posterior,
        int(estimator._rng.integers(0, 2**63)),
        estimator.num_transcripts,
        NUM_SPLICE_STRAND_COLS,
        n_threads,
        emit_locus_stats,
    )

    # Record results (same as current _run_locus_em)
    estimator._gdna_em_total = (
        getattr(estimator, '_gdna_em_total', 0.0) + total_gdna
    )
    # ... locus_results append, locus_stats recording ...
```

---

## Pipeline Integration

### Modified `quant_from_buffer()` Flow

```python
def quant_from_buffer(buffer, index, ...):
    # Phase 1–2: geometry + scoring — UNCHANGED
    geometry, estimator = _setup_geometry_and_estimator(...)
    em_data = _score_fragments(buffer, ...)  # → ScoredFragments (15.3 GB)

    # Phase 3: Locus construction + priors — UNCHANGED
    loci = build_loci(em_data, index)  # CC reads offsets + t_indices only
    locus_gammas = _compute_priors(estimator, loci, index, calibration)

    # Phase 4 (NEW): Array-by-array scatter + incremental free
    partitions = partition_and_free(em_data, loci)
    del em_data; gc.collect()  # shell is now empty (all arrays are None)

    # Phase 5 (NEW): Streaming locus EM with incremental partition freeing
    _run_locus_em_partitioned(
        estimator, partitions, loci, index,
        locus_gammas, em_config,
        emit_locus_stats=emit_locus_stats,
    )

    # Phase 6: cleanup
    gc.collect()
```

**No `AnnotationArrays` needed.** `frag_ids`, `frag_class`, `splice_type`
in ScoredFragments are dead weight — never read after construction.  The
AnnotationTable (populated during scanning) has its own independent copies.
These three arrays are freed as the first step of `partition_and_free()`.

### Backward Compatibility

The old `batch_locus_em()` C++ function and `_run_locus_em()` Python
function remain in the codebase during the transition for A/B testing.
They can be removed once golden tests pass with the new path.

---

## Implementation Steps

### Step 1: `LocusPartition` Dataclass
- **File:** `src/rigel/scored_fragments.py`
- **What:** Add `LocusPartition` dataclass (as shown above)
- **Risk:** None (additive)

### Step 2: C++ Scatter Functions
- **File:** `src/rigel/native/em_solver.cpp`
- **What:**
  - `build_partition_offsets()` — from `g_offsets` + `locus_units` list
  - `scatter_candidates<T>()` — per-candidate scatter using partition offsets
  - `scatter_units<T>()` — per-unit gather by unit index
  - Typed instantiations exposed to Python
- **Bindings:** `_em_impl` module
- **Testing:** Round-trip: for each locus, verify that the scattered
  per-locus array matches `global_arr[unit_global_indices]` (per-unit)
  or the concatenation of `global_arr[g_offsets[u]:g_offsets[u+1]]`
  segments (per-candidate)
- **Risk:** None — pure data reorganization

### Step 3: `partition_and_free()` Python Function
- **File:** `src/rigel/partition.py` (new file)
- **What:** Python driver that frees dead-weight arrays, builds offsets,
  iterates scatter functions array-by-array with gc.collect() between
  each, assembles LocusPartition objects
- **Testing:** Verify all global CSR arrays freed; per-locus data
  matches original
- **Risk:** None

### Step 4: `PartitionView` + `extract_locus_sub_problem_from_partition()`
- **File:** `src/rigel/native/em_solver.cpp`
- **What:** `PartitionView` struct (8 params instead of 20+) + extraction
  function that reads from contiguous partition data.  Drops dead params
  `all_t_starts`/`all_t_ends`.
- **Testing:** **Critical gate:** for each locus, verify the
  `LocusSubProblem` is bit-for-bit identical to the old function's
  output given the same input data
- **Risk:** Low — same algorithm, different data layout

### Step 5: `batch_locus_em_partitioned()` C++ Function
- **File:** `src/rigel/native/em_solver.cpp`
- **What:** New function that unpacks list-of-tuples under GIL (matching
  `fused_score_buffer` convention), then runs identical scheduling +
  `process_locus` pipeline with the new extraction function
- **Bindings:** `_em_impl` module
- **Testing:** Bit-for-bit identical EM output vs old `batch_locus_em()`
- **Risk:** Low

### Step 6: Python Pipeline Integration
- **Files:** `src/rigel/pipeline.py`, `src/rigel/estimator.py`
- **What:** Replace `_run_locus_em()` with `_run_locus_em_partitioned()`,
  modify `quant_from_buffer()` flow, add `_call_batch_em_partitioned()`
  wrapper
- **Testing:** Golden output tests must pass without `--update-golden`
- **Risk:** Low

### Step 7: Validation
- All existing tests pass (bit-for-bit where possible)
- New unit tests for scatter functions and extraction equivalence
- Memory profiling on CAPAN-1: verify peak RSS ≤ 20 GB
- Wall-time benchmark: verify overhead < 10s
- RSS tracking at each gc.collect() to verify arrays are actually freed

---

## Invariants & Correctness Properties

1. **Data equivalence:** For every locus, the `LocusSubProblem` from
   `extract_locus_sub_problem_from_partition()` must be identical to
   `extract_locus_sub_problem()`.  This is the critical invariant.

2. **Locus isolation:** Each transcript and each unit belongs to exactly
   one locus (connected components guarantee).  Partitioning is lossless.

3. **No dead weight:** `frag_ids`, `frag_class`, `splice_type` are freed
   immediately.  No `AnnotationArrays` class exists.

4. **Thread safety:** Identical to `batch_locus_em`.  Different loci write
   to disjoint output indices.  Shared `total_gdna_em` uses
   `std::atomic<double>`.

5. **Determinism:** Per-locus RNG seed = `rng_seed ^ (locus_id × 0x9e3779b97f4a7c15ULL)`.

6. **Memory monotonicity during scatter:** After each scatter+free step,
   total memory = constant (≈17.7 GB for CAPAN-1).  Swap-over peak =
   constant + sizeof(current_array).

---

## Complexity & Performance

| Operation | Complexity | Wall Time (CAPAN-1) |
|-----------|-----------|---------------------|
| Free dead-weight arrays | O(1) | < 1s |
| `build_partition_offsets` | O(total_units) | < 1s |
| 6 × `scatter_candidates` + free | O(total_candidates) each | ~5s total |
| 5 × `scatter_units` + free | O(total_units) each | ~2s total |
| LocusPartition assembly | O(n_loci) | < 1s |
| Tuple packing + GIL pointer extraction | O(n_loci × 12) | ~40ms |
| `batch_locus_em_partitioned` (mega) | Same as today | ~310s |
| `batch_locus_em_partitioned` (normal) | Same as today | ~68s |
| **Total overhead** | | **~8s** |

The 8s overhead is 1% of total wall time (757s).

**Bonus: cache locality improvement.**  The mega-locus extract phase
currently does scattered reads across 15 GB.  With partition data, the
mega-locus's 4.2 GB of data is contiguous — the L3 cache and hardware
prefetcher can serve sequential reads efficiently.  Profiling data from
CAPAN-1 shows the extract phase takes 11.2s for the mega-locus; cache
locality may reduce this.

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Scatter peak exceeds 32 GB | **Impossible** — bounded at ~20 GB by construction | — | Architecture prevents |
| GIL pointer extraction too slow | Very low (~40ms) | Low | Profile; switch to parallel lists if needed |
| Extraction produces different LocusSubProblem | Low | **Critical** | Per-locus bit-for-bit test |
| Wall-time regression from scatter | Very low (8s / 757s) | Low | Benchmark before/after |
| Python GC doesn't free arrays | Low | Medium | `gc.collect()` + RSS tracking |

---

## What This Plan Does NOT Do

- **Modify `build_loci()`** — uses `Locus.unit_indices` directly
- **Create `AnnotationArrays`** — dead-weight arrays freed, not preserved
- **Create `concatenate_partitions()`** — no re-concatenation needed
- **Create `run_single_locus_em()`** — unified function handles n=1 and n=14K
- **Change any EM algorithm** — identical dedup, SQUAREM, assignment logic
- **Pass `t_starts`/`t_ends` to extraction** — dead parameters dropped

Every new artifact (scatter functions, PartitionView, partition_and_free,
batch_locus_em_partitioned) is permanent in the final architecture.

---

## Summary

This plan eliminates the monolithic ScoredFragments bottleneck through two
composable mechanisms:

1. **Array-by-array scatter**: Python drives a loop calling typed C++ scatter
   functions, freeing each global array immediately.  Dead-weight arrays
   (`frag_ids`, `frag_class`, `splice_type`) are freed before scattering
   begins.  Peak memory is bounded by construction.

2. **Unified partition-native EM**: A single `batch_locus_em_partitioned()`
   C++ function replaces `batch_locus_em()`, accepting per-locus data via
   the established list-of-tuples convention.  The new extraction function
   has a clean 8-parameter signature (vs 20+), drops dead parameters, and
   benefits from contiguous partition data for better cache behavior.

There is no re-concatenation, no throwaway bridge code, no `build_loci()`
modification, and no `AnnotationArrays` wrapper.  Every new artifact
exists permanently.

Expected outcome: **peak RSS drops from 28.8 GB to ~20 GB** (31%
reduction), with **~8s** wall-time overhead and **zero** risk of OOM
on 32 GB machines.
