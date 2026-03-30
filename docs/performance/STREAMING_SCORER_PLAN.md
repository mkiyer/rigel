# Streaming Chunk Scoring: Implementation Plan

**Goal:** Eliminate the ~10 GB chunk-array peak during fragment scoring by
streaming chunks one at a time through a stateful C++ `StreamingScorer`.

**Expected savings:** Router phase peak RSS drops from ~23.6 GB to ~16.4 GB
(CAPAN-1 dataset, 14 GB BAM, 57.4M fragments, 54 spilled chunks).

## Background

### Current Architecture (all chunks in memory)

`FragmentRouter._scan_native()` loads ALL buffer chunks into a Python list,
then passes the entire list to the C++ `fused_score_buffer()` method:

```python
# scan.py — current code
chunk_arrays = []
for chunk in buffer.iter_chunks():
    chunk_arrays.append(chunk.to_scoring_arrays())  # ~10 GB total

result = native_scorer.fused_score_buffer(chunk_arrays, ...)
```

The C++ side extracts raw pointers for every chunk into a `vector<ChunkPtrs>`
and stores them in `FillState::all_cps`. Multimapper (MM) fragments accumulate
`(chunk_idx, row_idx)` pairs in `FillState::mm_members`. When the MM group is
flushed, `flush_mm_group()` reads raw chunk data via `cps[chunk_idx]` — this is
why ALL chunk arrays must remain alive simultaneously.

### Why This Is Wasteful

The only reason chunks are retained after processing is the MM back-reference
pattern. Multimapper alignments are ~8% of fragments (4.6M of 57.4M for
CAPAN-1), and MM groups are tiny (median ~3 alignments). The ~10 GB chunk
array retention is driven entirely by the `(chunk_idx, row_idx)` storage
pattern in `mm_members`.

### Key Invariant: Name-Sorted Contiguity

Because Rigel requires name-sorted/collated BAM input, all alignments for a
given `frag_id` are strictly contiguous. Chunk boundaries may split an MM group
(one group spans at most 2 chunks), but **the moment a new `frag_id` appears,
the previous group is guaranteed complete**. The existing code already flushes
on `frag_id` transition — it just reads chunk data during the flush, which is
the only thing keeping chunks alive.

## Design: Eager MM Scoring

The core insight: instead of storing `(chunk_idx, row_idx)` back-references
and re-reading chunk data at flush time, we **score each MM alignment eagerly**
as it arrives and merge into running accumulators. At flush time, no chunk data
is needed — only the accumulators.

### What `flush_mm_group()` Currently Does (3 phases)

1. **Score + merge** — iterates `mm_members`, reads chunk data via `cps[ci]`,
   scores each candidate, merges into `unordered_map<int32_t, MergedMrna>`
2. **Prune + emit** — pool-separated likelihood pruning, push to CSR vectors
3. **Finalize unit metadata** — aggregate splice type, compute best gDNA LL

### What Changes: All 3 Phases Become Incremental

Phase 1 moves to a new `score_mm_alignment()` method called during
`score_chunk_impl()` — it reads the *current* chunk's live arrays and merges
into persistent accumulators in `FillState`. Phase 3 metadata is also
accumulated eagerly. The flush path does only phase 2 (prune + emit).

### FillState MM Accumulators (replace `mm_members` + `all_cps`)

| Field | Type | Updated | Purpose |
|-------|------|---------|---------|
| `mm_merged` | `unordered_map<int32_t, MergedMrna>` | Each MM alignment | Scored candidates merged across alignments |
| `mm_is_any_spliced` | `bool` | Each MM alignment | OR'd splice check |
| `mm_best_stype` | `int` | Each MM alignment | max(ANNOT > UNANNOT > UNSPLICED) |
| `mm_best_gdna_ll` | `double` | Each unspliced MM alignment | Best gDNA log-likelihood |
| `mm_best_gdna_fp` | `int32_t` | Each unspliced MM alignment | Footprint for best gDNA |
| `mm_first_gfp` | `int32_t` | First MM alignment | Fallback footprint for spliced groups |
| `mm_n_members` | `int` | Each MM alignment | Count (replaces members.size()) |

**Memory re-use:** `mm_merged.clear()` retains heap allocation for bucket
reuse across groups — virtually zero-allocation per group.

**Memory footprint:** Active MM group holds ~100 alignments max (~4 KB),
vs ~10 GB for all chunk arrays. Six orders of magnitude improvement.

## Stateful API: `StreamingScorer`

```cpp
class StreamingScorer {
    const NativeFragmentScorer& scorer_;  // borrowed reference
    FillState st_;                        // growing CSR + MM accumulators
    // Growing output vectors (heap-allocated, capsule-transferred to Python)
    std::vector<int64_t>*  v_offsets_;
    // ... (17 output vectors total)
    // Chimeric fragment accumulator
    std::vector<int64_t>*  v_chim_fid_;
    std::vector<uint8_t>*  v_chim_stype_;
    // Statistics
    int64_t stat_det_, stat_em_u_, ...;
    double  gdna_log_sp_;

public:
    StreamingScorer(NativeFragmentScorer& scorer,
                    i8_1d t_strand_arr,
                    f64_2d_mut unambig_counts,
                    double gdna_log_splice_pen_unspliced);

    void score_chunk(nb::tuple chunk_arrays);  // score one chunk, free after
    nb::tuple finish();                         // flush final MM, return CSR
};
```

### Python-Side Usage

```python
# scan.py — new streaming code
from .native import StreamingScorer

scorer = StreamingScorer(native_scorer, t_strand_arr,
                         estimator.unambig_counts, gdna_log_penalty)
n_processed = 0
for chunk in buffer.iter_chunks():
    arrays = chunk.to_scoring_arrays()
    scorer.score_chunk(arrays)
    n_processed += chunk.size
    # 'arrays' goes out of scope → GC immediately

result = scorer.finish()
# result has same 24-element tuple as fused_score_buffer + chimeric arrays
```

### Chimeric Fragment Handling

Chimeric fragments (`FRAG_CHIMERIC`) are currently `continue`'d past in the
C++ scoring loop, then their annotations are extracted from chunk arrays in
Python after scoring. With streaming, chunk arrays are freed before annotation
extraction.

Solution: during `score_chunk_impl()`, push `f_id[i]` and `s_type[i]` into
small vectors in FillState (`v_chim_fid`, `v_chim_stype`). Return these from
`finish()` alongside the CSR. Chimeric fragments are a tiny fraction (~0.1%),
so these vectors consume < 1 MB.

## Implementation Steps

### Phase A: C++ Refactor (`src/rigel/native/scoring.cpp`)

**Step 1: FillState MM accumulators**
Replace `mm_members` (vector<pair<int,int>>) and `all_cps` (pointer) with
the 7 eager accumulators listed above. Add chimeric accumulators
(`v_chim_fid`, `v_chim_stype`). Update FillState constructor.

**Step 2: `score_mm_alignment()` private method**
Extract the per-alignment scoring loop (currently inside `flush_mm_group`
lines 459-540) into a new method on `NativeFragmentScorer`:
```cpp
void score_mm_alignment(const ChunkPtrs& cp, int row,
                        FillState& st, double gdna_log_sp) const;
```
Reads from the chunk's live arrays. Merges candidates into `st.mm_merged`.
Updates `mm_is_any_spliced`, `mm_best_stype`, `mm_best_gdna_ll/fp`,
`mm_first_gfp`, `mm_n_members`.

**Step 3: Simplify `flush_mm_group()`**
Remove the scoring loop and all chunk data reads. Keep only:
- Prune `mm_merged` → emit to CSR vectors
- Write unit metadata from accumulators (`mm_best_stype`, etc.)
- Reset all accumulators: `mm_merged.clear()` (preserves buckets),
  reset scalars to defaults

**Step 4: Update MM branch in `score_chunk_impl()`**
Replace `st.mm_members.emplace_back(chunk_idx, i)` with
`score_mm_alignment(cp, i, st, gdna_log_sp)`. Add chimeric accumulation:
push `f_id[i]` and `s_type[i]` when `fclass == FRAG_CHIMERIC`.
Remove `chunk_idx` parameter (no longer meaningful).

**Step 5: `StreamingScorer` class**
New class that owns the growing CSR vectors and FillState. Three-method API:
constructor, `score_chunk()`, `finish()`. The scorer reference is borrowed
(must outlive the StreamingScorer).

**Step 6: Nanobind bindings**
Add `nb::class_<StreamingScorer>` with `.def(nb::init<...>)`,
`.def("score_chunk", ...)`, `.def("finish", ...)`.
Keep existing `fused_score_buffer` for backward compatibility.

### Phase B: Python Integration

**Step 7: Update `scan.py`**
Rewrite `_scan_native()` to use `StreamingScorer`. The streaming loop replaces
bulk chunk loading. Handle chimeric and det-unambig annotations from the
`finish()` result. Progress logging inline using `chunk.size`.

**Step 8: Update `native.py`**
Export `StreamingScorer` from `_scoring_impl`.

### Phase C: Cleanup

**Step 9: Remove `fused_score_buffer`**
After the streaming path is validated, remove the old method and all
associated dead code (`ChunkPtrs` vector, `all_cps`, `mm_members`).

## Verification

1. `pytest tests/ -v` — all tests pass with bit-for-bit identical results
2. `pytest tests/test_golden_output.py -v` — golden outputs identical
3. Unit test for MM groups spanning chunk boundaries → identical scoring
4. CAPAN-1 profiling: peak RSS drops from ~23.6 GB to ~16.4 GB
5. `python scripts/benchmark/analyze_golden.py` — no benchmark regression

## Memory Budget (CAPAN-1 Estimates)

| Component | Current | After Streaming |
|-----------|---------|-----------------|
| Index | 2.1 GB | 2.1 GB |
| Buffer (in-memory chunks + spill) | 2.0 GB | 2.0 GB |
| All chunks loaded for scoring | **~10.0 GB** | **0 GB** |
| Current chunk being scored | — | ~0.2 GB |
| Active MM group accumulators | — | ~4 KB |
| Growing CSR output vectors | ~14.2 GB | ~14.2 GB |
| **Router phase peak** | **~23.6 GB** | **~16.4 GB** |

## Files Modified

| File | Change |
|------|--------|
| `src/rigel/native/scoring.cpp` | FillState refactor, score_mm_alignment, StreamingScorer class, nanobind bindings |
| `src/rigel/scan.py` | `_scan_native()` rewrite to streaming loop |
| `src/rigel/native.py` | Export `StreamingScorer` |
