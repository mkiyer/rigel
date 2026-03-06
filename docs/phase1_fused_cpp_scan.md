# Phase 1: Fused C++ Scan — Detailed Implementation Plan

## 1. Executive Summary

**Goal:** Eliminate ~100–120s of Python overhead in `FragmentRouter.scan()` by
moving the per-fragment inner loop into a single C++ function that processes an
entire `_FinalizedChunk` at once, operating directly on the columnar NumPy
arrays without creating 11.2M Python `BufferedFragment` objects.

**Estimated savings:**

| Source of overhead | Current cost | Eliminated? |
|-|-|-|
| `chunk[i]` → `BufferedFragment` creation | 27.4s | ✅ Yes |
| `array.append` (92.8M calls) | 10.0s | ✅ Yes |
| `array.frombytes` (61.9M calls) | 6.0s | ✅ Yes |
| `_add_single_fragment` Python dispatch | 35.7s (self) | ✅ Yes |
| `_finalize_unit_metadata` | 5.0s | ✅ Yes |
| Pre-EM accumulation (is_antisense + loops) | ~20s | ✅ Yes |
| `assign_unambig` Python calls | 2.0s | ✅ Yes |
| Stats counting (Python attribute access) | ~2s | ✅ Yes |
| **Total estimated savings** | **~108s** | |

**What remains in Python:**
- Chunk iteration (`buffer.iter_chunks()`) — ~50 chunks, negligible
- Multimapper group handling (`_flush_mm_group`) — called ~50K times, complex
  dict-merge logic, only ~5s total; not worth initial C++ port
- Chimeric annotation recording (rare, annotation-table only)
- Final `_to_np()` conversion — already fast

## 2. Current Architecture

```
FragmentBuffer.iter_chunks()
    │
    ▼ (per chunk, ~250K frags each)
┌──────────────────────────────────────────────────────────────────────┐
│  Python for loop: for i in range(chunk.size)                        │
│                                                                      │
│  1. Classify: fc = frag_classes[i]  (vectorized, fast)              │
│  2. bf = chunk[i]          ← Creates BufferedFragment (27.4s total) │
│  3. Pre-EM accumulation    ← is_antisense + 6 array updates        │
│  4. If SPLICED_ANNOT+UNAMBIG: estimator.assign_unambig(bf)          │
│  5. If MULTIMAPPER: group by frag_id, flush via _flush_mm_group()   │
│  6. Else: _add_single_fragment(bf, chunk, i, fc)                    │
│     └── nc.score_emit_fragment(bf.*)  ← C++ per-fragment            │
│         └── 9× array.frombytes() to append results                  │
│  7. Append per-unit metadata (8 array.append calls)                 │
└──────────────────────────────────────────────────────────────────────┘
    │
    ▼
ScoredFragments (CSR arrays)
```

The per-fragment C++ call `score_emit_fragment` is already fast (~2µs) but is
called 10.3M times from Python, with ~5µs of Python overhead per call
(creating `BufferedFragment`, unpacking, `frombytes`, metadata appends).

## 3. Target Architecture

```
FragmentBuffer.iter_chunks()
    │
    ▼ (per chunk, ~250K frags each)
┌──────────────────────────────────────────────────────────────────────┐
│  Python: ONE call per chunk                                         │
│                                                                      │
│  1. Classify: frag_classes = chunk.fragment_classes  (vectorized)    │
│  2. Identify MM fragments (frag_classes == FRAG_MULTIMAPPER)        │
│  3. Call C++: scan_chunk(chunk_arrays..., frag_classes, ...)        │
│     └── C++ processes ALL non-MM, non-chimeric fragments            │
│         └── Scoring + CSR emit + pre-EM accum + assign_unambig      │
│         └── Returns: CSR arrays + metadata arrays + updated accums  │
│  4. Python: flush MM groups (unchanged)                             │
│  5. Python: record chimeric annotations (unchanged, rare)           │
└──────────────────────────────────────────────────────────────────────┘
    │
    ▼
ScoredFragments (CSR arrays)
```

All non-multimapper, non-chimeric fragments (~95% of total) are processed
entirely in C++ with zero per-fragment Python calls.

## 4. C++ API Design

### 4.1. New Method: `NativeFragmentScorer::scan_chunk`

Added to the existing `NativeFragmentScorer` class in `scoring.cpp`.

```cpp
nb::tuple scan_chunk(
    // ---- Chunk columnar arrays (from _FinalizedChunk) ----
    i32_1d  t_offsets,         // int64[N+1] -> CSR offsets into t_indices etc
    i32_1d  t_indices,         // int32[M]   -> flat transcript indices
    i32_1d  frag_lengths,      // int32[M]   -> per-candidate fragment lengths
    i32_1d  exon_bp,           // int32[M]   -> per-candidate exon base pairs
    i32_1d  intron_bp,         // int32[M]   -> per-candidate intron base pairs
    i32_1d  unambig_intron_bp, // int32[M]   -> per-cand unambig intron bp
    u8_1d   splice_type,       // uint8[N]   -> per-fragment splice type
    u8_1d   exon_strand,       // uint8[N]   -> per-fragment exon strand
    u8_1d   frag_classes,      // uint8[N]   -> classification (0-4)
    i64_1d  frag_id,           // int64[N]   -> fragment IDs
    u32_1d  read_length,       // uint32[N]  -> per-fragment read length
    i32_1d  genomic_footprint, // int32[N]   -> per-fragment genomic footprint
    i32_1d  genomic_start,     // int32[N]   -> per-fragment genomic start
    u16_1d  nm,                // uint16[N]  -> NM tag values

    // ---- Index arrays for pre-EM strand classification ----
    i8_1d   t_to_strand_arr,   // int8[T]   -> transcript strand

    // ---- Pre-EM accumulator arrays (modified in-place) ----
    f64_mut exonic_sense,           // float64[T] writable
    f64_mut exonic_antisense,       // float64[T] writable
    f64_mut unspliced_sense,        // float64[T] writable
    f64_mut unspliced_antisense,    // float64[T] writable
    f64_mut intronic_sense,         // float64[T] writable
    f64_mut intronic_antisense,     // float64[T] writable
    f64_2d_mut unambig_counts,     // float64[T, 6] writable

    // ---- gDNA scoring parameters ----
    double  gdna_splice_penalty_unspliced,  // typically 1.0
    double  gdna_splice_penalty_unannot     // typically 0.01
);
```

**Returns** a tuple of numpy arrays + scalar stats:
```
(
    csr_offsets,         // int64[K+1]  — CSR offsets for K new EM units
    csr_t_indices,       // int32[C]    — flattened candidate transcript indices
    csr_log_liks,        // float64[C]  — flattened log-likelihoods
    csr_count_cols,      // uint8[C]    — flattened count columns
    csr_coverage_weights,// float64[C]  — flattened coverage weights
    csr_tx_starts,       // int32[C]    — transcript start positions
    csr_tx_ends,         // int32[C]    — transcript end positions
    locus_t,             // int32[K]    — best transcript per unit
    locus_ct,            // uint8[K]    — best count_col per unit
    is_spliced,          // int8[K]     — is_spliced flag per unit
    gdna_ll,             // float64[K]  — gDNA log-likelihood per unit
    genomic_footprints,  // int32[K]    — genomic footprint per unit
    frag_ids,            // int64[K]    — fragment ID per unit
    frag_class,          // int8[K]     — fragment class per unit
    splice_types,        // uint8[K]    — splice type per unit
    // Stats
    n_det_unambig,       // int   — deterministic unambig count
    n_em_unambig,        // int   — EM-routed unambig count
    n_em_ambig_same,     // int   — EM-routed ambig-same-strand count
    n_em_ambig_opp,      // int   — EM-routed ambig-opp-strand count
    n_gated_out,         // int   — fragments gated out (no winners)
    n_chimeric,          // int   — chimeric fragments skipped
)
```

### 4.2. Input Array Type Aliases

New aliases needed for mutable and wider-type arrays:

```cpp
using i64_1d   = nb::ndarray<const int64_t,  nb::ndim<1>, nb::c_contig>;
using u8_1d    = nb::ndarray<const uint8_t,  nb::ndim<1>, nb::c_contig>;
using u16_1d   = nb::ndarray<const uint16_t, nb::ndim<1>, nb::c_contig>;
using u32_1d   = nb::ndarray<const uint32_t, nb::ndim<1>, nb::c_contig>;
using f64_mut  = nb::ndarray<double, nb::ndim<1>, nb::c_contig>;
using f64_2d_mut = nb::ndarray<double, nb::ndim<2>, nb::c_contig>;
```

### 4.3. Internal Data Flow (per fragment, inside C++)

For each fragment index `i` in `[0, N)`:

1. **Skip** if `frag_classes[i] == FRAG_CHIMERIC` (4) or `frag_classes[i] == FRAG_MULTIMAPPER` (3).

2. **Slice** the per-candidate arrays: `t_inds[start:end]`, `exon_bp[start:end]`,
   etc. using `t_offsets[i]` and `t_offsets[i+1]`.

3. **Pre-EM accumulation** (for `FRAG_UNAMBIG` or `FRAG_AMBIG_SAME_STRAND`):
   - Get `first_t = t_indices[start]`, `t_strand = t_to_strand_arr[first_t]`
   - Compute `is_anti` using the same strand-likelihood logic already in
     `score_emit_fragment` (compare `exon_strand[i]` to `t_strand`).
   - For each candidate `k` in `[start, end)`:
     - `t_idx = t_indices[k]`
     - `weight = 1.0 / n_cand`
     - If `unambig_intron_bp[k] <= 0`: accumulate `exonic_sense/antisense[t_idx] += weight`
     - If `splice_type[i] == SPLICE_UNSPLICED`: accumulate `unspliced_sense/antisense[t_idx] += weight`
     - If `unambig_intron_bp[k] > 0`: accumulate `intronic_sense/antisense[t_idx] += weight`

4. **Deterministic unambig**: If `frag_classes[i] == FRAG_UNAMBIG` AND
   `splice_type[i] == SPLICE_ANNOT`:
   - Compute `is_anti` (same as above)
   - `col = splice_type * 2 + is_anti`
   - `unambig_counts[t_idx, col] += 1.0`
   - Increment `n_det_unambig` counter
   - **Do NOT emit** to CSR (deterministic units bypass EM)

5. **EM-routed scoring** (all other non-MM, non-chimeric fragments):
   - Call the existing `score_emit_fragment` logic **inline** (not as a method
     call — just reuse the same code path already in the class).
   - This produces: `out_ti`, `out_ll`, `out_ct`, `out_cw`, `out_ts`, `out_te`,
     `best_ll`, `best_t`, `best_ct`.
   - If `out_ti` is non-empty (fragment not gated out):
     - Append all candidate data to the output vectors.
     - Record offset, metadata (`is_spliced`, `gdna_ll`, `genomic_footprint`,
       `frag_id`, `frag_class`, `splice_type`, `locus_t`, `locus_ct`).
   - If empty: increment `n_gated_out`.
   - Increment appropriate stats counter based on `frag_classes[i]`.

6. **NOTE on deterministic-unambig + EM-routed for FRAG_UNAMBIG**:
   In the current Python code, `FRAG_UNAMBIG + SPLICE_ANNOT` fragments go to
   `assign_unambig` AND skip `_add_single_fragment`. All other `FRAG_UNAMBIG`
   (non-`SPLICE_ANNOT`) go to `_add_single_fragment`. The C++ must reproduce
   this exact routing.

### 4.4. gDNA Log-Likelihood Computation (inline)

For unspliced fragments (not `SPLICE_ANNOT` and not `SPLICE_UNANNOT`):

```cpp
double gdna_ll = LOG_HALF + frag_len_log_lik(genomic_footprint[i])
               + log(max(gdna_splice_pen, 1e-10));
```

For spliced fragments: `gdna_ll = -inf`.

This replaces `_gdna_log_lik()` and `_finalize_unit_metadata()`.

### 4.5. `is_antisense` Logic (duplicated from scoring)

The strand classification is already computed inside `score_emit_fragment` for
the mRNA scoring branch. For the pre-EM accumulation, we need the same logic
applied to the **first** candidate only:

```cpp
bool compute_is_antisense(int exon_strand, int t_strand) const {
    if (exon_strand == 1 || exon_strand == 2) {
        bool same_strand = (exon_strand == t_strand);
        double p_sense = same_strand ? exp(log_p_sense_) : exp(log_p_antisense_);
        // But we only need: p_sense < 0.5
        // Since log_p_sense > log(0.5) iff the model has strand specificity,
        // same_strand → p_sense = exp(log_p_sense) which is > 0.5 for SS libs
        // → NOT antisense.  !same_strand → p_sense = exp(log_p_antisense)
        // which is < 0.5 → IS antisense.
        //
        // Actually, the Python code uses:
        //   p = strand_model.strand_likelihood_int(exon_strand, ref_strand)
        //   return p < 0.5
        // where strand_likelihood_int returns the sense probability.
        //
        // For the C++ scorer, the relevant parameters are already available:
        //   same_strand → p = exp(log_p_sense_)
        //   !same_strand → p = exp(log_p_antisense_)
        //   is_anti = p < 0.5
        //
        // Since exp(log_p_sense_) ≥ 0.5 always (SS model is calibrated),
        // and exp(log_p_antisense_) < 0.5 always,
        // this simplifies to: is_anti = !same_strand (for SS libraries)
        // and for unstranded: both are exp(-0.693) = 0.5, is_anti = false.
        //
        // But to be safe, use the actual comparison:
        bool same = (exon_strand == static_cast<int>(t_strand));
        double log_p = same ? log_p_sense_ : log_p_antisense_;
        return std::exp(log_p) < 0.5;
    }
    return false;  // no strand info → not antisense
}
```

**Key nuance**: The Python `is_antisense` uses the `exonic_spliced` strand
model. The C++ scorer already holds `log_p_sense_` and `log_p_antisense_` which
are constructed from the same strand model during `FragmentScorer.from_models()`.
So these values are identical.

### 4.6. `assign_unambig` Logic

The Python `assign_unambig` does:
```python
t_idx = int(next(iter(bf.t_inds)))
t_strand = int(index.t_to_strand_arr[t_idx])
anti = self.is_antisense(bf.exon_strand, t_strand, sm)
col = splice_type * 2 + int(anti)
self.unambig_counts[t_idx, col] += 1.0
```

In C++, this becomes a simple array write inside `scan_chunk`. The
`unambig_counts` array is passed as a mutable `float64[T, 6]` and modified
in-place.

## 5. Output Strategy

### 5.1. Pre-allocated Output Vectors

The C++ function uses `std::vector` for all output arrays, growing dynamically.
After processing all fragments, the vectors are converted to NumPy arrays via
nanobind's `nb::ndarray` (zero-copy when possible, or via `nb::bytes` +
`frombytes` as currently done).

**Better approach**: Return raw numpy arrays directly. Nanobind supports
constructing `nb::ndarray` from `std::vector` data with a capsule-based
deleter. This avoids the intermediate `bytes` + `frombytes` overhead entirely.

```cpp
// Helper: wrap std::vector<T> as a numpy array, transferring ownership
template<typename T>
nb::ndarray<nb::numpy, T, nb::ndim<1>> vec_to_numpy(std::vector<T>&& vec) {
    size_t n = vec.size();
    T* data = vec.data();
    // Create a capsule that owns the vector
    auto* owner = new std::vector<T>(std::move(vec));
    nb::capsule deleter(owner, [](void* p) noexcept {
        delete static_cast<std::vector<T>*>(p);
    });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(data, {n}, deleter);
}
```

### 5.2. Python-Side Concatenation

After each chunk's `scan_chunk` call returns CSR arrays, the Python scan loop
concatenates them:

```python
all_offsets = []
all_t_indices = []
# ... etc

for chunk in buffer.iter_chunks():
    result = nc.scan_chunk(chunk.t_offsets, chunk.t_indices, ...)
    # Unpack result arrays
    all_offsets.append(result.csr_offsets)
    all_t_indices.append(result.csr_t_indices)
    # ... handle MMs ...

# Final concatenation
offsets = np.concatenate(all_offsets)
t_indices = np.concatenate(all_t_indices)
```

The offset arrays from each chunk need rebasing (adding the cumulative
candidate count). This is a simple vectorized operation.

## 6. Multimapper Handling

Multimappers (`FRAG_MULTIMAPPER`, ~4% of fragments) are **excluded** from
`scan_chunk` and continue to be handled by the existing Python
`_flush_mm_group()` method. Reasons:

1. MM handling requires grouping by `frag_id` across multiple alignment records,
   dict-merge of WTA scores, and complex splice-type priority logic.
2. Total MM cost is only ~5s (~3.5% of scan time).
3. The existing `_score_wta_mrna`/`_score_wta_nrna` C++ kernels already
   accelerate the scoring part.

The `scan_chunk` C++ function simply skips `frag_classes[i] == FRAG_MULTIMAPPER`,
and the Python loop collects MM fragments for group-flush.

**Chimeric fragments** are also skipped in C++ (`frag_classes[i] == FRAG_CHIMERIC`).
The Python loop records chimeric annotations if the annotation table is active.

## 7. Implementation Steps

### Step 1: Add `scan_chunk` to `scoring.cpp`

1. Add new type aliases for `u8_1d`, `u16_1d`, `u32_1d`, `i64_1d`, `f64_mut`,
   `f64_2d_mut`.

2. Add a private helper `compute_is_antisense(int exon_strand, int t_strand)`.

3. Add the `scan_chunk` method to `NativeFragmentScorer`:
   - Iterate over all fragments in the chunk
   - Skip chimeric and multimapper
   - Pre-EM accumulation for FRAG_UNAMBIG / FRAG_AMBIG_SAME_STRAND
   - Deterministic unambig for FRAG_UNAMBIG + SPLICE_ANNOT
   - Score+emit for everything else (reusing the existing inline scoring code)
   - Track stats counters

4. Register the method in the nanobind module definition.

### Step 2: Refactor `FragmentRouter.scan()` in `scan.py`

1. If `self._native_ctx` is available, use the new chunk-level C++ path:
   ```python
   for chunk in buffer.iter_chunks():
       frag_classes = chunk.fragment_classes

       # --- C++ processes all non-MM, non-chimeric ---
       result = self._native_ctx.scan_chunk(
           chunk.t_offsets, chunk.t_indices,
           chunk.frag_lengths, chunk.exon_bp,
           chunk.intron_bp, chunk.unambig_intron_bp,
           chunk.splice_type, chunk.exon_strand,
           frag_classes, chunk.frag_id,
           chunk.read_length, chunk.genomic_footprint,
           chunk.genomic_start, chunk.nm,
           index.t_to_strand_arr,
           estimator.transcript_exonic_sense,
           estimator.transcript_exonic_antisense,
           estimator.transcript_unspliced_sense,
           estimator.transcript_unspliced_antisense,
           estimator.transcript_intronic_sense,
           estimator.transcript_intronic_antisense,
           estimator.unambig_counts,
           gdna_pen_unspliced, gdna_pen_unannot,
       )
       # Unpack CSR arrays, append to accumulation lists
       # Handle MMs (unchanged)
       # Handle chimeric annotations (unchanged)
   ```

2. The existing per-fragment Python path remains as a fallback when
   `self._native_ctx` is None.

3. After all chunks:
   - Concatenate per-chunk CSR arrays with offset rebasing
   - Concatenate per-chunk metadata arrays
   - Merge MM CSR data (from `_flush_mm_group`, which still uses `array.array`)
   - Build final `ScoredFragments`

### Step 3: Update `FragmentScorer.from_models()` in `scoring.py`

Pass additional parameters to `NativeFragmentScorer` constructor if needed.
Currently the constructor already receives all necessary index arrays. No
changes expected unless we need gDNA splice penalty values persisted in the
C++ object (see Step 1).

**Alternative**: Pass gDNA splice penalties as constructor args so they're
available in `scan_chunk` without extra per-call overhead. Add two new
constructor params: `gdna_splice_penalty_unspliced` and
`gdna_splice_penalty_unannot`.

### Step 4: Test & Validate

1. Run all 805 existing tests — must all pass.
2. Verify numerical equivalence: add a test that runs the same buffer through
   both the old per-fragment path (with `_native_ctx = None`) and the new
   chunk-level path, comparing all output arrays element-by-element.
3. Run the profiler on the same 10M-read dataset to measure speedup.

## 8. Key Design Decisions

### 8.1. Modify Estimator Arrays In-Place

The 6 pre-EM accumulator arrays (`exonic_sense/antisense`,
`unspliced_sense/antisense`, `intronic_sense/antisense`) + `unambig_counts`
are passed as **mutable** numpy arrays to C++. The C++ function writes to them
directly. This avoids allocating per-chunk delta arrays and summing them later.

Nanobind supports this via `nb::ndarray<double, nb::ndim<1>>` (without `const`).

### 8.2. Return CSR Arrays as NumPy (Not Bytes)

Unlike `score_emit_fragment` which returns `bytes` for `array.frombytes()`, the
new `scan_chunk` returns proper NumPy arrays. The Python side concatenates these
with `np.concatenate`, which is efficient for ~50 chunks.

### 8.3. Stats as Return Values

Pipeline stats counters (deterministic_unambig_units, em_routed_*, n_gated_out)
are returned as integers in the tuple, not modified in-place. The Python side
adds them to the `PipelineStats` object.

### 8.4. CSR Offset Handling

Each chunk produces its own local CSR offsets starting from 0. The Python merger
rebases them:

```python
cumulative = 0
for chunk_offsets, chunk_t_indices in zip(all_offsets, all_t_indices):
    chunk_offsets += cumulative
    cumulative += len(chunk_t_indices)
```

### 8.5. Annotation Table

Chimeric annotations require the Python `AnnotationTable.add()` method. The
C++ function does not handle annotations — it simply counts chimeric fragments.
The Python loop handles chimeric annotation recording in a small separate pass:

```python
if annotations is not None:
    chimeric_mask = (frag_classes == FRAG_CHIMERIC)
    for i in np.where(chimeric_mask)[0]:
        annotations.add(
            frag_id=int(chunk.frag_id[i]),
            best_tid=-1, best_gid=-1,
            pool=POOL_CODE_CHIMERIC, posterior=0.0,
            frag_class=FRAG_CHIMERIC, n_candidates=0,
            splice_type=int(chunk.splice_type[i]),
        )
```

Similarly, deterministic-unambig annotations (if active) need a separate pass:

```python
if annotations is not None:
    # Collect det-unambig fragment info from C++ result
    det_t_indices = result.det_unambig_t_indices  # int32[D]
    det_frag_ids = result.det_unambig_frag_ids    # int64[D]
    for j in range(len(det_t_indices)):
        tid = int(det_t_indices[j])
        gid = int(t_to_g[tid])
        annotations.add(
            frag_id=int(det_frag_ids[j]),
            best_tid=tid, best_gid=gid,
            pool=POOL_CODE_MRNA, posterior=1.0,
            frag_class=FRAG_UNAMBIG, n_candidates=1,
            splice_type=int(SPLICE_ANNOT),
        )
```

For this to work, `scan_chunk` also returns:
- `det_unambig_t_indices`: int32 array of transcript indices for det-unambig
- `det_unambig_frag_ids`: int64 array of fragment IDs for det-unambig

These are small (~800K entries) and only needed when annotations are active.
Always return them regardless (negligible overhead).

## 9. Memory Considerations

The current approach uses `array.array` with dynamic growth. The C++ approach
uses `std::vector` with the same growth pattern. For a typical 250K-fragment
chunk with ~5 candidates/fragment, the output vectors are ~1.25M entries
(~10 MB per candidate array). This is negligible.

Total additional memory: ~50 MB peak (one chunk's output arrays at a time).

## 10. Risk Mitigation

| Risk | Mitigation |
|-|-|
| Numerical mismatch (C++ vs Python) | Comparison test: run both paths, assert allclose |
| Off-by-one in CSR offsets | Test with known small inputs where manual verification is feasible |
| Strand model parameter mismatch | Assert `log_p_sense`/`log_p_antisense` match in Python and C++ |
| Build failure (new nanobind types) | Test incrementally — add type aliases first, then method |
| MM group correctness | MM path is unchanged — existing tests cover it |

## 11. Files Modified

| File | Change |
|-|-|
| `src/hulkrna/native/scoring.cpp` | Add `scan_chunk` method + helpers (~200 lines) |
| `src/hulkrna/scan.py` | Refactor `scan()` to call `scan_chunk` per chunk (~50 lines changed) |
| `src/hulkrna/scoring.py` | Pass gDNA penalties to NativeFragmentScorer constructor (~5 lines) |
| `CMakeLists.txt` | No change (scoring.cpp already built) |
| `tests/test_*.py` | Add numerical equivalence test |

## 12. Success Criteria

1. All 805 existing tests pass.
2. `fragment_router_scan` time drops from ~144s to ~40s (multimappers + Python overhead for chunk iteration).
3. Total pipeline time drops from ~416s to ~310s (conservatively).
4. Zero numerical differences in output (exact float64 match for scoring, allclose for accumulations due to summation order).

## 13. Implementation Order

1. **scoring.cpp**: Add type aliases → `compute_is_antisense` helper → `scan_chunk` method → nanobind `.def()` registration.
2. **scoring.py**: Add gDNA penalty constructor args if needed.
3. **scan.py**: Add chunk-level C++ branch in `scan()`, keeping Python fallback.
4. **Build & test**: `pip install --no-build-isolation -e .` → `pytest tests/ -x -q`.
5. **Profile**: Run profiler to verify speedup.
