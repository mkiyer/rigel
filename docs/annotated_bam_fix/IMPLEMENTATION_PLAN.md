# Fix Annotated BAM Writer for EM-Routed Fragments

## Problem Statement

The annotated BAM writer only annotates ~5% of fragments (deterministic-unambiguous + chimeric). The remaining ~95% of fragments—those routed through the locus-level EM solver—receive default `ZP=intergenic, ZT=".", ZG="."` tags because their per-fragment assignments are never written to the `AnnotationTable`.

### Root Cause

During memory optimization refactoring, the `frag_id` primary key was discarded too early in the pipeline, severing the link between the EM solver's mathematical output and the physical BAM records.

Specifically:
1. **`partition.py:47`** — `em_data.frag_ids = None` destroys the frag_id array before EM runs.
2. **`assign_posteriors()` in `em_solver.cpp:1283`** — The C++ EM solver computes a per-unit `winner` (the transcript or gDNA component with highest posterior) but only scatters results to aggregate `em_counts_2d` count arrays. No per-fragment assignment output exists.
3. **`_run_locus_em_partitioned()` in `pipeline.py:423`** — Returns no per-fragment data to the caller; only aggregate `(gdna_em, locus_mrna, locus_gdna)`.

As a result, the `AnnotationTable` is only populated by the scoring pass (`scan.py:_scan_native`) for two fragment categories:
- Det-unambig fragments (SPLICED_ANNOT + FRAG_UNAMBIG): ~689K fragments
- Chimeric fragments: ~few hundred

The ~12.2M EM-routed fragments are never annotated.

### Why BAM Synchronization Works Without Order Preservation

The `frag_id` is the absolute binding key between Pass 1 (scan) and Pass 2 (annotation):

1. **Pass 1**: The reader thread assigns `frag_id = 0, 1, 2, ...` strictly sequentially as it pulls read-name groups from the BAM file (`bam_scanner.cpp:1006`). Even though worker threads resolve fragments out of order, the frag_id is assigned by the reader *before* dispatch.
2. **EM Phase**: Processes fragments grouped by locus, completely out of original BAM order.
3. **Pass 2**: `BamAnnotationWriter` (`bam_scanner.cpp:1508`) re-reads the same BAM file sequentially and recreates the same sequential frag_id counter. Because the BAM file hasn't changed, the Nth fragment in Pass 2 is guaranteed to be the same fragment that received `frag_id = N` in Pass 1.

Therefore, as long as the `AnnotationTable` contains the correct annotation keyed by `frag_id`, the order of EM processing is entirely irrelevant. The `AnnotationTable.frag_id_to_row` dict acts as a perfect hash map for O(1) lookup.

## Architectural Approach

Two coordinated changes:

1. **Wire per-unit metadata through the partition system.** Stop discarding `frag_ids`, `frag_class`, and `splice_type` in `partition_and_free()`. Instead, scatter them into per-locus `LocusPartition` objects alongside existing per-unit arrays. This keeps data flowing forward through the pipeline naturally—no `.copy()`, no save-and-reconstruct, and per-locus freeing for mega-loci. The only new C++ code is a `scatter_units_i64` template instantiation for the int64 `frag_ids` array.

2. **Add conditional per-unit assignment output to the C++ EM solver.** Behind an `emit_assignments` flag (zero overhead when annotation is not requested), `assign_posteriors()` writes the winning transcript index and posterior for each unit into position-indexed output arrays. After C++ returns, Python reads the EM output alongside the `frag_ids`/`frag_class`/`splice_type` from the `LocusPartition` objects (which are still alive) and populates the `AnnotationTable` in batch.

### Key Design Property

`frag_ids` do **not** need to enter the C++ EM code path. The C++ side only knows about position-indexed units. The Python side holds the `frag_id -> physical BAM record` mapping via `LocusPartition` objects. After C++ returns position-indexed assignment arrays, Python joins them with partition metadata to build annotation entries. This keeps the C++ changes minimal and localized.

### What Changes vs. What Doesn't

**Modified files (7):**
- `src/rigel/native/em_solver.cpp` — `scatter_units_i64` binding, `assign_posteriors()`, `batch_locus_em_partitioned()`, nanobind binding
- `src/rigel/native.py` — import `scatter_units_i64`
- `src/rigel/scored_fragments.py` — add 3 fields to `LocusPartition`
- `src/rigel/partition.py` — stop discarding 3 arrays; scatter them
- `src/rigel/estimator.py` — `run_batch_locus_em_partitioned()` expanded return
- `src/rigel/pipeline.py` — `_run_locus_em_partitioned()`, `_call_batch_em()`, annotation population
- `src/rigel/annotate.py` — `AnnotationTable.add_batch()`, `_grow_to()`

**NOT modified:**
- C++ `PartitionView`, `LocusSubProblem` structs — frag_ids stay in Python land
- `scan.py` — scoring/routing logic unchanged
- `bam_scanner.cpp` — `BamAnnotationWriter` is already correct; it just needs data in the table

---

## Phase 1: Scatter Per-Unit Metadata Through Partitions

### 1A. Add `scatter_units_i64` Template Instantiation

**File:** `src/rigel/native/em_solver.cpp`, near line 2494

The `scatter_units_impl<T>` template (line 1591) already handles all types. Add a new instantiation for `int64_t`:

```cpp
    m.def("scatter_units_i64", &scatter_units_impl<int64_t>,
          nb::arg("global_arr"), nb::arg("locus_units"), nb::arg("n_loci"),
          "Scatter per-unit int64 array into per-locus arrays.");
```

### 1B. Export from `native.py`

**File:** `src/rigel/native.py`

Add `scatter_units_i64` to the import from `_em_impl` and to `__all__`.

### 1C. Add Fields to `LocusPartition`

**File:** `src/rigel/scored_fragments.py`, `LocusPartition` dataclass (line ~174)

Add three new per-unit fields:

```python
@dataclass(slots=True)
class LocusPartition:
    # ... existing fields ...

    # Per-unit annotation metadata (not passed to C++ EM)
    frag_ids: np.ndarray       # int64[n_units] — buffer frag_id per unit
    frag_class: np.ndarray     # uint8[n_units] — fragment class code (stored as uint8)
    splice_type: np.ndarray    # uint8[n_units] — SpliceType enum value
```

### 1D. Stop Discarding Arrays; Scatter Them

**File:** `src/rigel/partition.py`

**Remove** the three discard lines (currently line 47-49):
```python
    # REMOVE these lines:
    em_data.frag_ids = None
    em_data.frag_class = None
    em_data.splice_type = None
    gc.collect()
```

**Add** `scatter_units_i64` to imports:
```python
from .native import (
    build_partition_offsets,
    scatter_candidates_f64,
    scatter_candidates_i32,
    scatter_candidates_u8,
    scatter_units_f64,
    scatter_units_i32,
    scatter_units_i64,      # <-- NEW
    scatter_units_u8,
)
```

**Add** 3 entries to `UNIT_ARRAYS`:
```python
    UNIT_ARRAYS = [
        ("gdna_log_liks", scatter_units_f64),
        ("genomic_footprints", scatter_units_i32),
        ("locus_t_indices", scatter_units_i32),
        ("locus_count_cols", scatter_units_u8),
        ("is_spliced", scatter_units_u8),
        # --- NEW: annotation metadata ---
        ("frag_ids", scatter_units_i64),
        ("frag_class", scatter_units_u8),
        ("splice_type", scatter_units_u8),
    ]
```

The existing dtype conversion block already handles `bool_ -> uint8` via `.view()`. Extend it to also handle `int8 -> uint8` for `frag_class`:

```python
    for attr, scatter_fn in UNIT_ARRAYS:
        global_arr = getattr(em_data, attr)
        if global_arr.dtype == np.bool_:
            global_arr = global_arr.view(np.uint8)
        elif global_arr.dtype == np.int8:           # <-- NEW
            global_arr = global_arr.view(np.uint8)  # <-- bit-exact reinterpretation
        unit_results[attr] = scatter_fn(global_arr, locus_units, n_loci)
        setattr(em_data, attr, None)
        del global_arr
        gc.collect()
```

**Add** 3 fields to the `LocusPartition` constructor in the assembly loop:

```python
    for li in range(n_loci):
        partitions[li] = LocusPartition(
            # ... existing fields ...
            locus_t_indices=unit_results["locus_t_indices"][li],
            locus_count_cols=unit_results["locus_count_cols"][li],
            # --- NEW ---
            frag_ids=unit_results["frag_ids"][li],
            frag_class=unit_results["frag_class"][li],
            splice_type=unit_results["splice_type"][li],
        )
```

---

## Phase 2: C++ Per-Unit Assignment Output (`em_solver.cpp`)

### 2A. Modify `assign_posteriors()` Signature

**File:** `src/rigel/native/em_solver.cpp`, line 1288

Add 4 parameters at the end for optional per-unit output:

```cpp
static void assign_posteriors(
    const LocusSubProblem& sub,
    const double* theta,
    int assignment_mode,
    double min_posterior,
    SplitMix64& rng,
    // Output accumulators (existing)
    double* em_counts_2d,
    double* gdna_locus_counts_2d,
    double* posterior_sum,
    double* n_assigned,
    double& mrna_total,
    double& gdna_total,
    int N_T_TOTAL,
    int n_cols,
    // --- NEW: per-unit assignment output (nullable) ---
    int32_t* out_winner_tid,   // nullptr = skip output
    float*   out_winner_post,
    int16_t* out_n_candidates,
    int      out_offset)
```

**Semantics of output arrays** (when `out_winner_tid != nullptr`):
- `out_winner_tid[out_offset + ui]`: **global** transcript index of the winning component, or `-2` if gDNA wins, or `-1` if no assignment
- `out_winner_post[out_offset + ui]`: posterior probability of the winner (float32)
- `out_n_candidates[out_offset + ui]`: number of candidates for this unit (int16)

### 2B. Add Per-Unit Output Logic

Insert between the `// Track max mRNA posterior` comment and the `// Scatter assignment weights` loop. The key insight: for MAP/sample modes, the `winner` variable is already computed. For fractional mode, find the MAP winner solely for annotation display:

```cpp
        // --- Per-unit annotation output ---
        if (out_winner_tid != nullptr) {
            // Find annotation winner:
            // - fractional mode: MAP of raw posteriors (for display only)
            // - MAP/sample modes: reuse already-computed winner
            int ann_winner = -1;
            if (assignment_mode == ASSIGN_FRACTIONAL) {
                double best_post = -1.0;
                for (int j = 0; j < seg_len; ++j) {
                    if (posteriors[j] > best_post) {
                        best_post = posteriors[j];
                        ann_winner = j;
                    }
                }
            } else {
                ann_winner = winner;
            }

            int write_idx = out_offset + ui;
            out_n_candidates[write_idx] = static_cast<int16_t>(
                std::min(seg_len, static_cast<int>(INT16_MAX)));

            if (ann_winner >= 0) {
                int32_t comp = sub.t_indices[s + ann_winner];
                if (comp < n_t) {
                    out_winner_tid[write_idx] = local_to_global[comp];
                } else {
                    out_winner_tid[write_idx] = -2;  // gDNA
                }
                out_winner_post[write_idx] = static_cast<float>(
                    posteriors[ann_winner]);
            } else {
                out_winner_tid[write_idx] = -1;
                out_winner_post[write_idx] = 0.0f;
            }
        }
```

### 2C. Modify `batch_locus_em_partitioned()` — Allocation and Return

**File:** `src/rigel/native/em_solver.cpp`, line 1840

#### Return type — expand from 4-tuple to 7-tuple:

```cpp
static std::tuple<
    double,                                          // total_gdna_em
    nb::ndarray<nb::numpy, double, nb::ndim<1>>,     // locus_mrna[n_loci]
    nb::ndarray<nb::numpy, double, nb::ndim<1>>,     // locus_gdna[n_loci]
    nb::list,                                        // locus_stats
    nb::object,                                      // out_winner_tid (ndarray or None)
    nb::object,                                      // out_winner_post (ndarray or None)
    nb::object                                       // out_n_candidates (ndarray or None)
>
batch_locus_em_partitioned(
    // ... existing params ...
    bool   emit_locus_stats,
    bool   emit_assignments)           // <-- NEW parameter
```

#### Allocate output arrays — after extracting PartitionViews, before `gil_scoped_release`:

```cpp
    // --- Per-unit assignment output arrays ---
    int total_units = 0;
    std::vector<int> locus_write_offsets(n_loci, 0);
    for (int i = 0; i < n_loci; ++i) {
        locus_write_offsets[i] = total_units;
        total_units += views[i].n_units;
    }

    int32_t* out_tid_ptr   = nullptr;
    float*   out_post_ptr  = nullptr;
    int16_t* out_ncand_ptr = nullptr;

    std::vector<int32_t> out_tid_vec;
    std::vector<float>   out_post_vec;
    std::vector<int16_t> out_ncand_vec;

    if (emit_assignments && total_units > 0) {
        out_tid_vec.assign(total_units, -1);
        out_post_vec.assign(total_units, 0.0f);
        out_ncand_vec.assign(total_units, 0);
        out_tid_ptr   = out_tid_vec.data();
        out_post_ptr  = out_post_vec.data();
        out_ncand_ptr = out_ncand_vec.data();
    }
```

**Thread safety:** Each locus writes to `[locus_write_offsets[li], locus_write_offsets[li] + n_units)` — non-overlapping sections, no synchronization needed.

#### Update `assign_posteriors` call in `process_locus` lambda:

```cpp
            assign_posteriors(
                sub, result.theta.data(),
                assignment_mode, assignment_min_posterior, locus_rng,
                em_out, gdna_out,
                psum_out, nass_out,
                locus_mrna, locus_gdna,
                N_T, N_COLS,
                out_tid_ptr, out_post_ptr, out_ncand_ptr,
                locus_write_offsets[li]);
```

#### Construct return arrays — after GIL reacquired, alongside existing `stats_list`:

```cpp
    nb::object py_out_tid   = nb::none();
    nb::object py_out_post  = nb::none();
    nb::object py_out_ncand = nb::none();

    if (emit_assignments && total_units > 0) {
        size_t u_shape[1] = {static_cast<size_t>(total_units)};

        auto* tid_copy = new int32_t[total_units];
        std::memcpy(tid_copy, out_tid_vec.data(), total_units * sizeof(int32_t));
        nb::capsule tid_owner(tid_copy, [](void* p) noexcept {
            delete[] static_cast<int32_t*>(p); });
        py_out_tid = nb::cast(nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(
            tid_copy, 1, u_shape, std::move(tid_owner)));

        auto* post_copy = new float[total_units];
        std::memcpy(post_copy, out_post_vec.data(), total_units * sizeof(float));
        nb::capsule post_owner(post_copy, [](void* p) noexcept {
            delete[] static_cast<float*>(p); });
        py_out_post = nb::cast(nb::ndarray<nb::numpy, float, nb::ndim<1>>(
            post_copy, 1, u_shape, std::move(post_owner)));

        auto* ncand_copy = new int16_t[total_units];
        std::memcpy(ncand_copy, out_ncand_vec.data(), total_units * sizeof(int16_t));
        nb::capsule ncand_owner(ncand_copy, [](void* p) noexcept {
            delete[] static_cast<int16_t*>(p); });
        py_out_ncand = nb::cast(nb::ndarray<nb::numpy, int16_t, nb::ndim<1>>(
            ncand_copy, 1, u_shape, std::move(ncand_owner)));
    }

    return std::make_tuple(
        total_gdna_em_val,
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(
            mrna_copy, 1, shape, std::move(mrna_owner)),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(
            gdna_copy, 1, shape, std::move(gdna_owner)),
        stats_list,
        py_out_tid,
        py_out_post,
        py_out_ncand
    );
```

### 2D. Update Nanobind Module Binding

**File:** `src/rigel/native/em_solver.cpp`, line 2499

Add the new parameter:

```cpp
    m.def("batch_locus_em_partitioned", &batch_locus_em_partitioned,
          // ... existing args ...
          nb::arg("n_threads") = 0,
          nb::arg("emit_locus_stats") = false,
          nb::arg("emit_assignments") = false,   // <-- NEW
          "Run locus EM from per-locus partition data.\n\n"
          "Returns (total_gdna_em, locus_mrna, locus_gdna, locus_stats,\n"
          " out_winner_tid, out_winner_post, out_n_candidates).");
```

---

## Phase 3: Python Estimator (`estimator.py`)

### 3A. Update `run_batch_locus_em_partitioned()`

**File:** `src/rigel/estimator.py`, line 204

Add `emit_assignments` parameter, pass through to C++, unpack and return expanded tuple:

```python
    def run_batch_locus_em_partitioned(
        self,
        # ... existing params ...
        emit_locus_stats: bool = False,
        emit_assignments: bool = False,       # <-- NEW
    ) -> tuple[float, np.ndarray, np.ndarray,
               np.ndarray | None, np.ndarray | None, np.ndarray | None]:

        # ... existing setup code ...

        total_gdna_em, locus_mrna, locus_gdna, locus_stats_raw, \
            out_winner_tid, out_winner_post, out_n_candidates = (
            _batch_locus_em_partitioned(
                # ... existing args ...
                emit_locus_stats,
                emit_assignments,
            )
        )

        # ... existing locus_stats handling ...

        return (
            total_gdna_em,
            np.asarray(locus_mrna),
            np.asarray(locus_gdna),
            out_winner_tid,
            out_winner_post,
            out_n_candidates,
        )
```

---

## Phase 4: Pipeline Integration (`pipeline.py`)

### 4A. Update `_call_batch_em()` Inner Function

**File:** `src/rigel/pipeline.py`, inside `_run_locus_em_partitioned()`

```python
    emit_assignments = annotations is not None

    def _call_batch_em(parts, batch_loci, batch_gammas, batch_spans):
        partition_tuples = [
            (
                p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
                p.tx_starts, p.tx_ends, p.count_cols,
                p.is_spliced, p.gdna_log_liks, p.genomic_footprints,
                p.locus_t_indices, p.locus_count_cols,
            )
            for p in parts
        ]
        locus_t_lists = [l.transcript_indices for l in batch_loci]

        return estimator.run_batch_locus_em_partitioned(
            partition_tuples,
            locus_t_lists,
            batch_gammas,
            batch_spans,
            index,
            em_iterations=em_config.iterations,
            em_convergence_delta=em_config.convergence_delta,
            emit_locus_stats=emit_locus_stats,
            emit_assignments=emit_assignments,
        )
```

### 4B. Add Annotation Population Helper

New module-level function in `pipeline.py`:

```python
def _populate_em_annotations(
    batch_parts,
    out_winner_tid,
    out_winner_post,
    out_n_candidates,
    annotations,
    index,
):
    if out_winner_tid is None:
        return

    t_to_g = index.t_to_g_arr
    is_syn = index.t_df["is_synthetic_nrna"].values

    # Concatenate per-locus metadata in same order as C++ output
    frag_ids = np.concatenate([p.frag_ids for p in batch_parts])
    frag_class = np.concatenate([p.frag_class for p in batch_parts])
    splice_type_arr = np.concatenate([p.splice_type for p in batch_parts])

    n = len(frag_ids)
    best_tid = np.asarray(out_winner_tid, dtype=np.int32)
    posteriors = np.asarray(out_winner_post, dtype=np.float32)
    n_cand = np.asarray(out_n_candidates, dtype=np.int16)

    # Gene index
    valid_t = best_tid >= 0
    best_gid = np.full(n, -1, dtype=np.int32)
    best_gid[valid_t] = t_to_g[best_tid[valid_t]]

    # Pool code
    pool = np.full(n, POOL_CODE_INTERGENIC, dtype=np.uint8)
    pool[best_tid == -2] = POOL_CODE_GDNA
    if valid_t.any():
        is_nrna = is_syn[best_tid[valid_t]].astype(bool)
        pool[valid_t] = np.where(
            is_nrna, POOL_CODE_NRNA, POOL_CODE_MRNA
        ).astype(np.uint8)

    # Clean sentinel: -2 (gDNA) -> -1
    best_tid_clean = best_tid.copy()
    best_tid_clean[best_tid == -2] = -1

    annotations.add_batch(
        frag_ids=frag_ids,
        best_tids=best_tid_clean,
        best_gids=best_gid,
        pools=pool,
        posteriors=posteriors,
        frag_classes=frag_class.view(np.int8),
        n_candidates=n_cand,
        splice_types=splice_type_arr,
    )
```

### 4C. Call After Each Batch EM

**Phase A (mega-loci):**
```python
    for locus in mega_loci:
        part = partitions.pop(locus.locus_id)
        gdna_em, mrna_arr, gdna_arr, out_tid, out_post, out_ncand = _call_batch_em(
            [part], [locus], ...)
        total_gdna_em += gdna_em
        if annotations is not None:
            _populate_em_annotations(
                [part], out_tid, out_post, out_ncand, annotations, index)
        ...
        del part
        gc.collect()
```

**Phase B (normal loci):**
```python
    if normal_loci:
        normal_parts = [partitions[l.locus_id] for l in normal_loci]
        gdna_em, mrna_arr, gdna_arr, out_tid, out_post, out_ncand = _call_batch_em(
            normal_parts, normal_loci, normal_gammas, normal_spans)
        total_gdna_em += gdna_em
        if annotations is not None:
            _populate_em_annotations(
                normal_parts, out_tid, out_post, out_ncand, annotations, index)
        ...
```

### 4D. Update Signatures

`_run_locus_em_partitioned()` — add `annotations` param.
`quant_from_buffer()` — pass `annotations=annotations` through.

---

## Phase 5: Batch Insert Method (`annotate.py`)

### 5A. Add `add_batch()` to `AnnotationTable`

```python
    def add_batch(self, frag_ids, best_tids, best_gids, pools,
                  posteriors, frag_classes, n_candidates, splice_types):
        n = len(frag_ids)
        if n == 0:
            return
        needed = self._size + n
        if needed > self.capacity:
            self._grow_to(max(self.capacity * 2, needed))
        start = self._size
        end = start + n
        self.frag_ids[start:end] = frag_ids
        self.best_tid[start:end] = best_tids
        self.best_gid[start:end] = best_gids
        self.pool[start:end] = pools
        self.posterior[start:end] = posteriors
        self.frag_class[start:end] = frag_classes
        self.n_candidates[start:end] = n_candidates
        self.splice_type[start:end] = splice_types
        self.frag_id_to_row.update(zip(frag_ids.tolist(), range(start, end)))
        self._size = end
```

### 5B. Add `_grow_to()` Method

```python
    def _grow_to(self, new_cap):
        if new_cap <= self.capacity:
            return
        for attr in ("frag_ids", "best_tid", "best_gid", "pool",
                     "posterior", "frag_class", "n_candidates", "splice_type", "locus_id"):
            old = getattr(self, attr)
            new = np.empty(new_cap, dtype=old.dtype)
            new[:self.capacity] = old
            if attr in ("frag_ids", "best_tid", "best_gid", "frag_class", "locus_id"):
                new[self.capacity:] = -1
            else:
                new[self.capacity:] = 0
            setattr(self, attr, new)
        self.capacity = new_cap

    def _grow(self):
        self._grow_to(max(self.capacity * 2, 1024))
```

---

## Data Flow Summary

```
Pass 1: BAM Scan
  |-- Reader thread: frag_id = 0, 1, 2, ... (sequential)
  |-- Workers: resolve fragments
  |-- Scoring (scan.py):
  |   |-- Det-unambig -> annotations.add(frag_id, tid, ..., pool=mRNA/nRNA)
  |   |-- Chimeric -> annotations.add(frag_id, ..., pool=chimeric)
  |   +-- EM-routed -> ScoredFragments (with frag_ids per unit)
  |
  |-- partition_and_free():
  |   |-- Scatters frag_ids/frag_class/splice_type into per-locus LocusPartition
  |   +-- Frees global arrays
  |
  |-- EM batch (C++ with emit_assignments=true):
  |   +-- assign_posteriors -> writes out_winner_tid[offset+ui], etc.
  |
  |-- Post-EM annotation (Python):
  |   |-- Reads frag_ids from partition.frag_ids (still alive)
  |   |-- Reads frag_class/splice_type from partition (still alive)
  |   |-- Joins with C++ position-indexed output
  |   +-- annotations.add_batch(frag_ids, best_tid, ...)
  |
  +-- AnnotationTable now has ~100% coverage

Pass 2: BAM Annotation (BamAnnotationWriter)
  |-- Re-reads BAM sequentially, frag_id = 0, 1, 2, ...
  |-- annotations.frag_id_to_row[frag_id] -> row -> annotation
  +-- Stamps ZT, ZG, ZP, ZW, ZC, ZH, ZN, ZS tags
```

---

## Verification Plan

### 1. Compile

```bash
conda activate rigel && pip install --no-build-isolation -e .
```

### 2. Unit Tests

```bash
pytest tests/test_annotate.py -v     # annotation-specific tests
pytest tests/ -v                     # full suite, no regressions
```

### 3. Integration Validation

```bash
rigel quant \
  --bam /scratch/.../sim_oracle.bam \
  --index /scratch/.../rigel_index \
  -o /tmp/test_annot \
  --annotated-bam /tmp/test_annot/annotated.bam \
  --threads 8 --seed 42 2>&1
```

Check `write_annotated_bam` summary:
- **Before fix:** `n_annotated ~ 689K, n_intergenic ~ 12.2M`
- **After fix:** `n_annotated ~ 12.9M, n_intergenic ~ few thousand` (true intergenic only)

### 4. Count Invariance

`test_annotated_bam_counts_match` verifies annotation does not change quantification results.

---

## Design Decisions

| Decision | Rationale |
|----------|-----------|
| Wire metadata through partition, not save+reconstruct | Eliminates `.copy()`, enables per-locus freeing for mega-loci, data flows forward naturally |
| frag_ids stay in Python land (not in C++ PartitionView) | C++ EM has no use for frag_ids; minimizes C++ surface area |
| `scatter_units_i64` template instantiation | Trivial — reuses existing template; 2 lines of C++ |
| `frag_class` (int8) scattered as uint8 | Same bit pattern; `.view(np.int8)` restores signed type |
| `nb::object` return for optional arrays | `nb::none()` when not emitting; Python `None` |
| `emit_assignments` flag | Zero overhead when not requested; default=false |
| MAP winner for fractional mode annotation | Annotation needs single winner for ZT/ZG tags |
| `add_batch()` with numpy slicing | Avoids 12M individual `add()` calls |
| Sentinel -2 for gDNA in C++ | Distinguishes gDNA from unassigned (-1) and transcript (>=0) |
