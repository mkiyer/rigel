# Performance Phase 2 — Optimization Plan

**Date**: 2026-04-14  
**Baseline**: 109.75s wall | 7930 MB peak RSS | 380K frags/s  
**Dataset**: VCaP simulated BAM, 19.4M buffered fragments, 457K transcripts  

## Current Profile

| Hotspot | Time | % | RSS Impact |
|---------|------|---|------------|
| `scan_and_buffer` | 71.0s | 64.7% | +2730 MB (fragment buffer) |
| `locus_em` | 22.2s | 20.3% | +2373 MB jump at partition |
| `fragment_router_scan` (scoring) | 11.7s | 10.6% | ~155 MB temp/chunk |
| `partition` | 4.8s | 4.4% | Peak RSS spike point |
| Other | 4.8s | 4.4% | — |

RSS snapshots (MB): before=2836, after_scan=5566, after_router_scan=5442, after_locus_em=5251, peak=7930 (during partition).

---

## R0: Revert P4v2 (Fused Partition Scatter)

### Explanation

The P4v2 optimization replaced 15 individual scatter functions with two fused C++ calls (`partition_candidates` + `partition_units`). While this simplified the code, profiling reveals it was a net negative:

- **RSS regression**: Peak RSS increased from 7044 MB (P5 baseline) to 7930 MB (+886 MB). The two-phase design holds all candidate globals alive during phase 1, causing a larger simultaneous memory footprint than the original incremental-free approach.
- **Negligible runtime gain**: Partition time went from ~5.0s to ~4.8s — within noise. Total wall time actually increased slightly (108.8s → 109.75s).
- **Root cause**: The original profiling analysis overestimated the memmove cost in scatter. Most memmove overhead was actually in the scoring phase (P2), not partition. The 15 individual scatter calls were already efficient because each freed its global array immediately after scatter, keeping peak RSS bounded.

### Implementation Plan

1. **Restore original scatter functions** in `em_solver.cpp`:
   - Re-add the 15 individual `scatter_*` functions that were removed in P4.
   - Remove `partition_candidates()` and `partition_units()` C++ functions.

2. **Restore `partition.py`** to pre-P4 state:
   - Replace the two-phase `partition_and_free()` with the original implementation that calls individual scatter functions sequentially, freeing each global array immediately after scatter.

3. **Restore `native.py` imports**:
   - Remove `partition_candidates`, `partition_units` exports.
   - Re-add the 8 individual scatter function imports.

4. **Restore `test_partition.py`**:
   - Re-add individual scatter function tests that were removed in P4.
   - Remove `TestPartitionAll` tests specific to fused API.

5. **Update `pipeline.py`** if the partition return type changed.

6. **Verify**: Run full test suite (`pytest tests/ -v`), confirm 1028 tests pass. Re-run profiler to confirm RSS returns to ~7044 MB baseline.

**Expected gain**: −886 MB peak RSS. No runtime change.

---

## O1: StreamingScorer Vector `reserve()`

### Explanation

The `StreamingScorer` in `scoring.cpp` manages 19 `std::vector` pointers that accumulate the scored CSR output. These vectors are allocated with **zero `reserve()` calls** and grow entirely via `push_back()`. For a typical run with ~48.5M candidate push_backs (19.4M fragments × ~2.5 candidates), this causes:

- **Repeated reallocations**: `std::vector` grows by 1.5–2× each time capacity is exceeded. For 48.5M elements, this means ~25–30 reallocation events per vector, each copying the entire contents.
- **Memory fragmentation**: Freed old allocations leave holes in the heap. With 19 vectors reallocating independently, this creates severe allocator fragmentation — a significant contributor to the gap between actual data size and observed RSS.
- **Cache pollution**: Reallocation touches cold memory, evicting hot working-set data from L1/L2 cache.

The fix is straightforward: after the first `score_chunk` call provides a chunk size, estimate the total output size and pre-allocate all vectors. This eliminates reallocation overhead entirely.

### Implementation Plan

1. **Add a `reserve_vectors()` helper** to `StreamingScorer` in `scoring.cpp`:
   ```cpp
   void reserve_vectors(int64_t est_units, int64_t est_cands) {
       v_offsets_->reserve(est_units + 1);
       v_fid_->reserve(est_units);
       v_fc_->reserve(est_units);
       v_st_->reserve(est_units);
       // Per-candidate vectors:
       v_ti_->reserve(est_cands);
       v_ll_->reserve(est_cands);
       v_ct_->reserve(est_cands);
       v_cw_->reserve(est_cands);
       v_ts_->reserve(est_cands);
       v_te_->reserve(est_cands);
       v_lt_->reserve(est_cands);
       v_lct_->reserve(est_cands);
       v_isp_->reserve(est_cands);
       v_gll_->reserve(est_cands);
       v_gfp_->reserve(est_cands);
       // Deterministic vectors — harder to estimate, use units as proxy:
       v_dti_->reserve(est_units / 2);
       v_dfid_->reserve(est_units / 2);
       v_chim_fid_->reserve(est_units / 10);
       v_chim_stype_->reserve(est_units / 10);
   }
   ```

2. **Call `reserve_vectors()` from `score_chunk()`** on the first invocation:
   - After the first chunk is processed, we know `chunk.N` (fragments per chunk) and can observe the candidate-to-fragment ratio.
   - Use a `bool reserved_ = false` flag. On first call: estimate total fragments from the hint passed at construction (or from the first chunk's size × expected number of chunks), compute `est_cands = est_units * avg_cands_per_unit`, call `reserve_vectors()`.
   - Alternatively, accept an `estimated_fragments` parameter in the `StreamingScorer` constructor and derive estimates upfront.

3. **Pass fragment count estimate from Python**:
   - In `scan.py`, the `FragmentRouter` knows `buffer.n_total` before scoring begins. Pass this to the `StreamingScorer` constructor or via a new `set_size_hint(n_fragments)` method.
   - Estimate candidates as `n_fragments * 2.5` (empirical average from profiling data).

4. **Verify**: Run tests, confirm identical output (reserve does not change behavior). Profile to measure fragmentation reduction and scoring speedup.

**Expected gain**: 2–5% scoring speedup (~0.3–0.6s), reduced heap fragmentation (indirect RSS benefit).

---

## O2: Eliminate int32→int64 Upcasts in Scoring Path

### Explanation

The fragment buffer stores `t_offsets` and `frag_id` as int32 (sufficient for per-chunk sizes up to 2B), but the C++ `score_chunk()` expects int64 arrays. The `to_scoring_arrays()` method in `buffer.py` performs:

```python
np.ascontiguousarray(self.t_offsets, dtype=np.int64)  # int32 → int64 copy
np.ascontiguousarray(self.frag_id, dtype=np.int64)    # int32 → int64 copy
```

For each chunk of ~1M fragments:
- `t_offsets`: 1M × 8 bytes = 8 MB copy (was 4 MB as int32)
- `frag_id`: 1M × 8 bytes = 8 MB copy (was 4 MB as int32)
- **Per chunk**: ~16 MB of unnecessary temporary allocation
- **Across ~19 chunks**: ~300 MB total temp allocations, plus the copy time

More importantly, the C++ scoring kernel only uses `t_offsets` for CSR indexing (values fit in int32 within a chunk) and `frag_id` as an opaque identifier (also fits int32 within a chunk). There is no arithmetic reason to require int64.

### Implementation Plan

1. **Change `ChunkPtrs` in `scoring.cpp`**:
   - Change `t_off` from `const int64_t*` to `const int32_t*`.
   - Change `f_id` from `const int64_t*` to `const int32_t*`.

2. **Update `score_chunk()` casts**:
   ```cpp
   cp.t_off = nb::cast<i32_1d>(chunk_arrays[0]).data();  // was i64_1d
   cp.f_id  = nb::cast<i32_1d>(chunk_arrays[8]).data();  // was i64_1d
   ```

3. **Update `score_chunk_impl()`** internal usage:
   - Anywhere `cp.t_off[i]` is used as an index, it's already int32-safe (chunk offsets < 2^31).
   - Anywhere `cp.f_id[i]` is used, it's stored into `st_.v_fid_` which is `vector<int64_t>` — the widening happens at push_back time (implicit `int32→int64` promotion), which is free.

4. **Update `to_scoring_arrays()` in `buffer.py`**:
   ```python
   # Remove dtype=np.int64 upcasts:
   np.ascontiguousarray(self.t_offsets, dtype=np.int32),  # no copy needed
   np.ascontiguousarray(self.frag_id, dtype=np.int32),    # no copy needed
   ```

5. **Verify**: Run tests — output arrays from `finish()` should still be int64 (the widening happens on push_back into `v_offsets_` and `v_fid_`). Check that `v_offsets_` accumulates correctly using global offsets (it uses `push_back(v_ti_->size())` which is int64).

6. **Handle the global offset accumulation in `v_offsets_`**: The StreamingScorer's `v_offsets_` is int64 because it tracks the **global** CSR offset across all chunks (can exceed int32 for >2B candidates). The chunk-local `t_off` is int32, but `v_offsets_` must remain int64. Verify that the offset accumulation adds the chunk-local offset to a running int64 base — this is already the case since it uses `v_ti_->size()`.

**Expected gain**: ~155 MB peak RSS reduction per chunk, ~50ms total copy overhead eliminated.

---

## O3: Index String Column Integer Coding

### Explanation

The rigel index loads transcript metadata into Pandas DataFrames (`t_df`, etc.) with string columns for reference names (`ref`), gene IDs (`g_id`), gene names (`g_name`), and gene types (`g_type`). With 457K transcripts:

- **String columns**: Each unique string is a Python `str` object (56+ bytes overhead). With ~25 unique `ref` values, ~63K unique `g_id` values, and ~63K unique `g_name` values, the 457K-row DataFrame stores 457K Python string object references per column — even though most values are duplicates.
  - Estimated waste: 100–150 MB across 4 string columns.

- **frozenset objects in `_iv_t_set`**: The interval-to-transcript-set mapping uses `frozenset` objects. Each frozenset has 232+ bytes overhead. With ~200K–400K intervals, this consumes 50–100 MB.

- **`sj_map` string keys**: Splice junction map uses `(ref_name, pos1, pos2)` tuples with string ref names as keys. Hash table overhead with string hashing costs 30–50 MB.

The index baseline RSS is 2836 MB. While much of this is allocator overhead/fragmentation, ~200–350 MB is recoverable through integer coding.

### Implementation Plan

1. **Convert string columns to categorical/integer codes in `t_df`**:
   - In `index.py` (or wherever `t_df` is constructed), convert `ref`, `g_id`, `g_name`, `g_type` to `pd.Categorical`:
     ```python
     for col in ('ref', 'g_id', 'g_name', 'g_type'):
         t_df[col] = t_df[col].astype('category')
     ```
   - Categorical columns store one integer code per row (int8/int16) plus a small dictionary of unique values. For 457K rows × 4 columns: ~2 MB vs ~150 MB.

2. **Update downstream consumers**:
   - Audit all code that reads `t_df['ref']`, `t_df['g_id']`, etc. Categorical columns support `.values`, `.iloc[]`, equality checks, and groupby natively — most code should work unchanged.
   - Watch for code that calls `.str` accessor or does `isinstance(val, str)` checks on individual elements — these may need adjustment since `.iloc[i]` on a categorical returns the original string.

3. **Convert `_iv_t_set` from dict-of-frozensets to CSR**:
   - Replace `Dict[interval_key, frozenset[int]]` with two arrays:
     ```python
     iv_t_offsets: np.ndarray  # int32, length = n_intervals + 1
     iv_t_indices: np.ndarray  # int32, sorted transcript indices
     ```
   - Build CSR during index construction. Lookups become `t_indices[offsets[i]:offsets[i+1]]`.
   - Update all callers (primarily in `resolve.cpp` / `_resolve_impl`) to use the CSR interface.

4. **Convert `sj_map` to integer-keyed lookup**:
   - Replace string ref name with integer ref index (from the categorical encoding above).
   - Key becomes `(ref_idx, pos1, pos2)` — all integers, faster hashing.

5. **Persist categorical encoding in feather files**:
   - Feather (Arrow IPC) natively supports dictionary-encoded columns. When saving/loading `t_df`, categorical columns will round-trip efficiently.

6. **Verify**: Run full test suite. Compare index load RSS (expect ~200–350 MB reduction from 2836 MB baseline).

**Expected gain**: 200–350 MB RSS reduction in index baseline.

---

## O4: Mega-Locus Convergence Investigation

### Explanation

The EM solver exhibits two convergence anomalies visible in the profiling data:

1. **Mega-locus** (locus 0): 98,608 transcripts, 4.06M units, 140,407 equivalence classes. Takes **14.2s** with 59 SQUAREM iterations (9.9s in SQUAREM alone). This single locus consumes 64% of total EM time. Each iteration costs ~240ms due to the problem size (400+ components × 5K–10K units per EC pass).

2. **Slow small loci**: Locus 24525 has only 15 transcripts but requires **260 SQUAREM iterations** — far more than the mega-locus. This suggests a convergence pathology: poorly conditioned initialization, near-degenerate likelihood landscape, or adversarial prior configuration causing the solver to oscillate.

The SQUAREM accelerator is supposed to converge in fewer iterations than plain EM, but 260 iterations for a 15-transcript problem indicates the acceleration is failing (likely due to negative step-length clipping or repeated restarts).

### Implementation Plan

1. **Add convergence diagnostics to EM solver**:
   - In `em_solver.cpp`, add optional per-iteration logging controlled by a debug flag:
     ```cpp
     struct ConvergenceDiag {
         int iteration;
         double delta;        // max |α_new - α_old|
         double step_length;  // SQUAREM step length (r, v, alpha)
         bool was_clipped;    // step length was clipped to [-1, ∞)
         double log_lik;      // observed data log-likelihood
     };
     ```
   - Store diagnostics in a thread-local vector, return alongside EM results.
   - Add a Python-accessible flag in `EMConfig` to enable/disable diagnostic collection.

2. **Run diagnostic profiling**:
   - Enable diagnostics for the VCaP benchmark.
   - Extract convergence curves for: (a) the mega-locus, (b) locus 24525 (260 iterations), (c) a sample of "normal" loci for comparison.
   - Plot delta vs iteration, step-length vs iteration, log-likelihood vs iteration.

3. **Analyze failure modes**:
   - **SQUAREM step-length collapse**: If `alpha` is repeatedly clipped to small values, the accelerator degenerates to plain EM. Root cause: non-monotone extrapolation.
   - **Oscillation**: If delta oscillates rather than monotonically decreasing, the prior or initialization is adversarial.
   - **Near-zero components**: If many components have near-zero abundance, the solver wastes iterations on numerically insignificant updates.

4. **Potential fixes** (dependent on diagnostic findings):
   - **Component pruning**: After N iterations, drop components with abundance < ε threshold. Reduces per-iteration cost for mega-loci.
   - **Adaptive convergence threshold**: Use relative delta (`max(|Δα| / α)`) instead of absolute delta for faster convergence on small-abundance components.
   - **Initialization improvement**: Use coverage-weighted initialization instead of uniform priors for faster initial convergence.
   - **SQUAREM damping**: Adjust step-length bounds to prevent oscillation in pathological cases.

5. **Verify**: Any convergence changes must not alter quantification accuracy. Validate against golden outputs and benchmark truth abundances.

**Expected gain**: If mega-locus iterations reduce from 59→40: ~3–5s saved. If pathological small loci are fixed: ~1–2s additional. Total potential: 4–7s (18–31% of EM time).

---

## O5: Equivalence Class Sorting Pre-Check

### Explanation

The `build_equiv_classes()` function in `em_solver.cpp` groups scored candidates into equivalence classes (units sharing the same transcript set). After grouping, it sorts the log-likelihood and coverage-weight matrices to match the deterministic EC ordering. This sorting step allocates full temporary copies of the `ll_flat` and `wt_flat` matrices:

```cpp
// Allocates new_ll and new_wt: n_units × n_components × sizeof(double)
// For mega-locus: ~4M × 400 × 8 bytes ≈ 12.8 GB (theoretical max, actual is sparse)
```

In practice, the data is often **already in correct order** (when the input is naturally sorted by transcript set). The allocation and copy are wasted work in these cases.

### Implementation Plan

1. **Add identity-permutation check** before allocating sort temporaries:
   ```cpp
   // Compute sort permutation
   std::vector<int> perm = compute_ec_sort_order(ec_data);
   
   // Check if already sorted (identity permutation)
   bool already_sorted = true;
   for (int i = 0; i < (int)perm.size(); ++i) {
       if (perm[i] != i) { already_sorted = false; break; }
   }
   
   if (!already_sorted) {
       // Existing reorder logic: allocate new_ll, new_wt, apply permutation
       apply_permutation(perm, ll_flat, wt_flat);
   }
   ```

2. **Locate the exact reorder code** in `em_solver.cpp`:
   - Find `build_equiv_classes()` or the EC construction path.
   - Identify where the sort permutation is computed and where temporaries are allocated.
   - Insert the identity check between permutation computation and allocation.

3. **Verify**: Run tests to confirm identical EC construction. Profile to measure skip rate (what fraction of loci have already-sorted ECs).

**Expected gain**: 1–2% EM speedup (~0.2–0.4s). More impactful for mega-loci where the temporary allocation is large.

---

## Execution Order

| Priority | Item | Type | Expected Gain | Effort |
|----------|------|------|---------------|--------|
| 1 | R0 — Revert P4v2 | Revert | −886 MB RSS | Medium (restore old code) |
| 2 | O1 — Vector reserve() | Quick win | 2–5% scoring, RSS fragmentation | ~15 lines C++ |
| 3 | O2 — int32→int64 upcast | Quick win | −155 MB/chunk, −50ms | ~30 lines C++ + Python |
| 4 | O4 — Convergence investigation | Diagnostic | 4–7s potential (18–31% EM) | Investigation first |
| 5 | O3 — Index string coding | Refactor | −200–350 MB RSS baseline | Medium refactor |
| 6 | O5 — EC sort pre-check | Quick win | 1–2% EM | ~20 lines C++ |

**Cumulative target**: ~100–105s wall time, ~6200–6500 MB peak RSS (vs current 109.75s / 7930 MB).
