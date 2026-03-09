# P1 Detailed Implementation Plan: Multimapper Scoring in C++

## Overview

Move the entire multimapper (MM) scoring and CSR accumulation path from
Python (`_flush_mm_group` and its callees in `scan.py`) into the existing
C++ `fused_score_buffer` two-pass framework in `scoring.cpp`.

After implementation, the Python hot-path contribution drops from ~57s
to ~2s, eliminating ~200M Python function calls and ~55s of wall time.

---

## 1. Current Architecture (Before)

### Data Flow

```
_scan_native()
  │
  ├── Prepare chunk_arrays (14 columnar numpy arrays per chunk)
  │
  ├── nc.fused_score_buffer(chunk_arrays, ...)     ─── C++ TWO-PASS ───┐
  │     • Pass 1: count EM units + candidates, accumulate strand counts │
  │     • Pass 2: fill pre-allocated numpy arrays                       │
  │     • Skips FRAG_MULTIMAPPER fragments (continues past them)        │
  │     • Returns 17-element tuple of numpy arrays + 6 stat counters    │
  │                                                                     │
  │   Result: (cpp_offsets, cpp_ti, cpp_ll, ..., det_tids, det_fids,   │
  │            n_det, n_em_u, n_em_as, n_em_ao, n_gated, n_chim)      │
  │                                                                     │
  ├── Python MM loop ─── SLOW PATH ────────────────────────────────────┘
  │     for each chunk:
  │       mm_indices = np.where(frag_classes == FRAG_MULTIMAPPER)
  │       for idx in mm_indices:
  │         bf = chunk[idx]            ← creates BufferedFragment dataclass
  │         group by frag_id
  │         _flush_mm_group()          ← per-group Python scoring loop
  │           _score_wta_mrna(bf)      ← nc.score_wta_mrna() returns dict
  │           _score_wta_nrna(bf)      ← nc.score_wta_nrna() returns dict
  │           merge dicts across hits (best per transcript, global WTA)
  │           _emit_mrna(final_mrna)   ← appends to array.array CSR lists
  │           _emit_nrna(final_nrna)   ← appends to array.array CSR lists
  │           _gdna_log_lik(bf)        ← per-hit gDNA computation
  │           append per-unit metadata to array.array lists
  │
  ├── Merge: np.concatenate(cpp_arrays, mm_arrays)
  │
  └── Return ScoredFragments
```

### Key Functions Called on the MM Python Path

| Function | Location | Calls | Self (s) | What it does |
|---|---|---:|---:|---|
| `_flush_mm_group` | scan.py:510 | 885K | 9.9 | Orchestrates per-group scoring |
| `__getitem__` | buffer.py:186 | 3.39M | 9.3 | Creates BufferedFragment dataclass |
| `_score_wta_mrna` | scan.py:116 | 3.39M | 2.2 | mRNA WTA → returns `nb::dict` from C++ |
| `_score_wta_nrna` | scan.py:243 | 2.88M | 3.0 | nRNA WTA → returns `nb::dict` from C++ |
| `_emit_mrna` | scan.py:367 | 885K | 1.4 | Dict → array.array appends |
| `_emit_nrna` | scan.py:399 | 731K | 6.2 | Dict → array.array appends |
| `_gdna_log_lik` | scan.py:422 | 2.36M | 1.0 | Per-fragment gDNA score |
| `array.array.append` | (builtin) | 108M | 9.6 | CSR accumulation |
| `frag_len_log_lik` | scoring.py:254 | 2.36M | 0.4 | LUT lookup for gDNA path |

### Current C++ Methods on NativeFragmentScorer

- `score_wta_mrna(...)` → `nb::dict` of per-transcript scored results
- `score_wta_nrna(...)` → `nb::dict` of per-transcript scored results
- `score_emit_fragment(...)` → fused single-fragment scoring (dead code)
- `fused_score_buffer(...)` → two-pass bulk scoring (production path)

### Data Structures

**Chunk (14-element tuple, all numpy):**
```
[0]  t_offsets          int64[N+1]    CSR offsets into candidate arrays
[1]  t_indices          int32[M]      transcript indices (flat)
[2]  frag_lengths       int32[M]      per-candidate fragment lengths
[3]  exon_bp            int32[M]      exonic basepairs
[4]  intron_bp          int32[M]      intronic basepairs
[5]  unambig_intron_bp  int32[M]      unambiguous intron bp
[6]  splice_type        uint8[N]      per-fragment splice type
[7]  exon_strand        uint8[N]      per-fragment exon strand
[8]  frag_class         uint8[N]      computed fragment classification
[9]  frag_id            int64[N]      multimapper group ID
[10] read_length        uint32[N]     read length
[11] genomic_footprint  int32[N]      genomic span
[12] genomic_start      int32[N]      genomic start position
[13] nm                 uint16[N]     mismatch count
```

**Multimapper grouping:** Fragments with `frag_class == FRAG_MULTIMAPPER`
(3) and same `frag_id` form one MM group. Within a chunk, MM fragments
with the same frag_id are contiguous (name-sorted BAM → same molecule's
alignments are adjacent).

---

## 2. Target Architecture (After)

### Data Flow

```
_scan_native()
  │
  ├── Prepare chunk_arrays (same 14 arrays per chunk, no changes)
  │
  ├── nc.fused_score_buffer(chunk_arrays, ...)     ─── C++ TWO-PASS ───┐
  │     • Pass 1: count EM units + candidates for ALL fragment types    │
  │       INCLUDING FRAG_MULTIMAPPER (currently skipped)                │
  │       - Identify MM groups (contiguous same-frag_id runs)           │
  │       - Count their output units and candidates                     │
  │     • Pass 2: fill pre-allocated numpy arrays for ALL types         │
  │       - Score MM groups inline using the same mRNA/nRNA WTA logic   │
  │       - Merge across alignments (best per transcript, global WTA)   │
  │       - Compute gDNA log-lik (best unspliced hit)                   │
  │       - Write unit metadata (splice_type, is_spliced, gDNA, etc.)   │
  │     • Returns same 17-element tuple + 7 stat counters (add n_mm)    │
  │                                                                     │
  ├── NO PYTHON MM LOOP                                                 │
  │   (all MM scoring now handled inside fused_score_buffer)            │
  │                                                                     │
  ├── NO MERGE STEP                                                     │
  │   (fused_score_buffer returns unified arrays for all fragment types) │
  │                                                                     │
  └── Return ScoredFragments directly from C++ arrays
```

### Key Design Decisions

1. **Handle MM within the existing `score_chunk_impl<FillMode>` template**
   rather than adding a separate pass. The two-pass framework (count, then
   fill) already iterates all fragments — we just need to stop skipping
   `FRAG_MULTIMAPPER`.

2. **MM group detection in C++:** Scan linearly — when we encounter a
   `FRAG_MULTIMAPPER`, start accumulating. When the `frag_id` changes
   (or we hit a non-MM fragment, or end of chunk), flush the accumulated
   group. This mirrors the Python logic exactly.

3. **Cross-chunk MM groups:** The current Python code handles MM groups
   that span chunk boundaries (same frag_id at end of chunk N and start
   of chunk N+1). The C++ implementation must handle this too. Approach:
   keep a `pending_mm` state in `FillState` that persists across chunks.

4. **MM merge logic (critical correctness invariants):**
   - Per-alignment: score mRNA WTA and nRNA WTA independently (reuse
     existing `score_chunk_impl` mRNA/nRNA scoring logic)
   - Cross-alignment merge: for each transcript, keep the best score
     (lower overhang wins; ties broken by higher log-lik)
   - Global WTA gate: after merge, keep only minimum-overhang transcripts
   - nRNA gate: skip nRNA if any hit has `SPLICE_ANNOT`
   - gDNA: best gDNA log-lik across unspliced hits only
   - Unit splice_type: ANNOT > UNANNOT > UNSPLICED (most informative)
   - Unit footprint: best unspliced hit's footprint (if not spliced),
     else first hit's footprint

---

## 3. Implementation Steps

### Step 1: Add MM Group Scoring to `score_chunk_impl`

**File:** `src/rigel/native/scoring.cpp`

Currently at line ~887:
```cpp
if (fclass == FRAG_CHIMERIC) { ++stat_chim; continue; }
if (fclass == FRAG_MULTIMAPPER) continue;  // ← REMOVE THIS SKIP
```

**New logic:** Replace the MM skip with MM group accumulation.

```
Algorithm (inside score_chunk_impl loop):

  if (fclass == FRAG_MULTIMAPPER) {
      if (f_id[i] != mm_current_fid) {
          // New group — flush previous group if any
          flush_mm_group(...)
          mm_current_fid = f_id[i]
          mm_group_start = i     // first index in this group
          mm_group_size = 1
      } else {
          mm_group_size++
      }
      continue  // don't score individually, wait for flush
  }

  // If previous fragment was MM and this one isn't, flush
  if (mm_group_size > 0) {
      flush_mm_group(...)
      mm_group_size = 0
  }

  // ... existing non-MM scoring code ...
```

At end of `score_chunk_impl` (after the loop), flush any pending MM group.

**Cross-chunk handling:** `mm_current_fid`, `mm_group_start`, and
`mm_group_size` are stored in `FillState` and persist across chunk calls.
When a new chunk starts and the first MM fragment has the same frag_id as
the pending group from the previous chunk, we must handle the continuation.

**Approach for cross-chunk:** Before processing a new chunk, record the
previous chunk's pending state. After the loop, if there's a pending group,
don't flush it — let it carry over to the next chunk. After all chunks are
processed, flush any final pending group. This requires storing the chunk
pointer index alongside the group state.

**Simpler approach:** Since MM fragments with the same frag_id are always
contiguous within a chunk (from the name-sorted BAM), and chunks are
processed in order, we only need to track a single `(frag_id, start_idx,
count)` triple across chunks. If the first MM of the next chunk has the
same frag_id, extend the group; otherwise flush and start new.

**Implementation note:** For cross-chunk groups, we need access to the
previous chunk's data during the flush. This means we must store the
group members' indices (chunk_idx, row_idx pairs) rather than just a
start/count within one chunk. Use a small `std::vector<std::pair<int,int>>`
for pending MM group members.

### Step 2: Implement `flush_mm_group` in C++

**New private method on `NativeFragmentScorer`:**

```cpp
template<bool FillMode>
void flush_mm_group(
    const std::vector<ChunkPtrs>& all_chunks,
    const std::vector<std::pair<int,int>>& mm_members, // (chunk_idx, row_idx)
    FillState& st,
    double gdna_log_sp,
    int64_t& stat_mm) const
```

**Logic (mirrors `_flush_mm_group` in scan.py exactly):**

```
1. Determine splice status across all hits:
   - is_any_spliced = any hit has ANNOT or UNANNOT
   - is_annot_spliced = any hit has ANNOT

2. For each hit in mm_members:
   a. Extract per-fragment data from chunk arrays
      (t_inds, exon_bp, frag_lengths, etc.)
   b. Score mRNA WTA (same logic as existing inline scoring)
      → produces per-candidate (t_idx, overhang, log_lik, count_col,
        coverage_wt, tx_start, tx_end)
   c. Merge into merged_mrna map: for each transcript, keep best
      (lower overhang wins, ties broken by higher log_lik)

3. Global mRNA WTA: find min overhang across merged_mrna, keep only
   candidates with that overhang

4. If !is_annot_spliced:
   For each hit (skipping ANNOT ones):
     a. Score nRNA WTA
     b. Merge into merged_nrna map (same merge logic)
   Global nRNA WTA: min overhang filter

5. Emit winners:
   CountMode:  increment st.cand_cur for each winner
   FillMode:   write candidates into output arrays at st.cand_cur

6. If any candidates emitted:
   a. Record unit metadata:
      - offsets[unit_cur+1] = cand_cur
      - locus_t = best_t (highest mRNA log_lik among winners)
      - locus_ct = best_count_col
      - fid = frag_id of the group
      - fclass = FRAG_MULTIMAPPER
      - stype = most informative (ANNOT > UNANNOT > UNSPLICED)
      - is_spliced = is_any_spliced
      - gDNA: if !is_any_spliced, compute best gDNA log-lik across
        unspliced hits; else -inf
      - genomic_footprint: best unspliced hit's footprint, or first
        hit's footprint if spliced
   b. Increment unit_cur and stat_mm
   else:
     Increment stat_gated
```

**Data structure for merged candidates:** Use
`std::unordered_map<int32_t, MergedCandidate>` where:
```cpp
struct MergedCandidate {
    int32_t overhang;
    double log_lik;
    int32_t count_col;  // mRNA only
    double coverage_wt;
    int32_t tx_start, tx_end;
};
```

For typical MM groups (3-4 hits × 5-20 candidates = 15-80 entries), an
unordered_map is fine. For very large groups, consider a flat sorted
vector, but profiling suggests this isn't needed.

### Step 3: Update `fused_score_buffer` Return Value

Add `stat_mm` (multimapper unit count) to the returned tuple. Currently
returns 23 elements; add one more for `n_mm`. Update the Python unpacking.

### Step 4: Simplify `_scan_native` in Python

Remove:
- The `mm_indices` / `__getitem__` / `_flush_mm_group` loop
- The array.array CSR accumulators on `FragmentRouter.__init__`
- The merge step (`np.concatenate(cpp_arrays, mm_arrays)`)
- The `_to_np` helper

The method becomes simply:
```python
def _scan_native(self, chunks, buffer, log_every):
    result = nc.fused_score_buffer(chunk_arrays, ...)
    (cpp_offsets, cpp_ti, ..., n_det, n_em_u, ..., n_mm) = result
    # Update stats
    # Handle annotations (det-unambig + chimeric)
    # Return ScoredFragments directly from C++ arrays
```

### Step 5: Update Tests

The existing tests in `test_pipeline_routing.py` exercise MM scoring
through the full `FragmentRouter.scan()` → `_scan_native()` path.
They use mock chunks with `__getitem__` support. After P1:

- These tests still work (same input, same output)
- The mock `_Chunk.__getitem__` is no longer called by production code
  but may still be used by tests for setup — verify and clean up
- Add new tests that specifically verify C++ MM scoring correctness:
  - MM group with mixed splice types (ANNOT + UNSPLICED)
  - MM group where one hit has better overhang than others
  - MM group where two hits map to the same transcript (merge)
  - Cross-chunk MM group (same frag_id at chunk boundary)
  - Large MM group (> SCORED_STACK_CAPACITY candidates)
  - MM group where all candidates are gated out (zero WTA survivors)

---

## 4. Concerns and Design Issues

### Concern 1: Cross-Chunk MM Groups

**Issue:** MM fragments with the same frag_id may span chunk boundaries.
The current Python code handles this naturally (it accumulates `mm_pending`
across chunks). In C++, the two-pass `score_chunk_impl` processes one
chunk at a time, so we need state that persists across chunk calls.

**Mitigation:** Add `mm_pending_fid`, `mm_pending_members` (vector of
`(chunk_idx, row_idx)` pairs) to `FillState`. At each chunk boundary,
check if the pending group continues. After all chunks processed, flush
the final pending group. Both count and fill passes must handle this
identically.

**Risk:** Medium. The two-pass design means we must track the same
cross-chunk state in both passes. The count pass doesn't have output
arrays yet, so it only increments cursors. The fill pass writes into
pre-allocated arrays. Both must flush at exactly the same points to
maintain count/fill agreement.

### Concern 2: Count/Fill Pass Agreement

**Issue:** The two-pass design requires that pass 1 (count) and pass 2
(fill) produce identical unit/candidate counts. Any divergence causes
buffer overflows or incorrect results.

**Mitigation:** The existing `score_chunk_impl` template already handles
this — pass 1 has `FillMode=false` and pass 2 has `FillMode=true`. The
same template code runs both passes. The MM scoring logic must be
similarly templated so both passes execute identical control flow.

**Risk:** Low if we keep the MM logic inside the same template function.

### Concern 3: MM Merge Uses Unordered Map

**Issue:** `std::unordered_map<int32_t, MergedCandidate>` has overhead
per insert/lookup compared to a flat array. For MM groups with ~3.8 hits
× ~10 candidates = ~38 entries, this is negligible. But pathological cases
(e.g., a 100-hit MM group mapping to 5000 transcripts) could be slow.

**Mitigation:** Profile after implementation. If needed, switch to a flat
`std::vector` sorted by t_idx with binary search, or a stack-allocated
flat hash table. For now, `std::unordered_map` is fine for correctness.

### Concern 4: Strand Accumulation for MM Fragments

**Issue:** The current `fused_score_buffer` accumulates strand counts
(`acc_es`, `acc_ea`, etc.) for `FRAG_UNAMBIG` and `FRAG_AMBIG_SAME_STRAND`
only. MM fragments are **not** included in strand accumulation (the Python
path also doesn't accumulate them). This is correct and should not change.

**Mitigation:** Just make sure the MM flush does NOT accumulate strand
counts. The existing `if (fclass == FRAG_UNAMBIG || fclass == FRAG_AMBIG_SAME_STRAND)`
guard already excludes MM.

### Concern 5: `__getitem__` and BufferedFragment Still Needed for Tests

**Issue:** The test mock `_Chunk.__getitem__` returns `_BF` instances
(not `BufferedFragment`). The production `_FinalizedChunk.__getitem__`
returns `BufferedFragment` instances. After P1, the production code no
longer calls `__getitem__` on chunks — but tests may still construct
`_BF` / `BufferedFragment` for other purposes.

**Mitigation:** Keep `BufferedFragment` and `__getitem__` as a debug/test
utility. Mark with a comment noting they are not on the production path.
Remove from production hot path only.

### Concern 6: `_gdna_log_lik` Uses Full Splice Penalty Map

**Issue:** The Python `_gdna_log_lik` looks up `ctx.gdna_splice_penalties`
by `bf.splice_type`. The C++ fused scorer only passes in the UNSPLICED
penalty. For the MM path, the gDNA computation only runs on unspliced hits
(the `_flush_mm_group` code explicitly skips ANNOT/UNANNOT before calling
`_gdna_log_lik`). So using only the UNSPLICED penalty is correct.

**Mitigation:** In the C++ MM flush, only compute gDNA for unspliced MM
hits, using the same `gdna_log_sp` parameter already passed to
`fused_score_buffer`. No additional penalty map needed.

### Concern 7: Python `_score_wta_mrna`/`_score_wta_nrna` Already Dispatch to C++

**Issue:** The Python methods `_score_wta_mrna` and `_score_wta_nrna`
already dispatch to `nc.score_wta_mrna()` and `nc.score_wta_nrna()` on
the `NativeFragmentScorer`. So the MM path currently goes:
`Python(_flush_mm_group) → C++(score_wta_mrna) → Python(dict merge) →
Python(dict iteration for emit)`. The overhead is in the Python↔C++
round-trips (creating `nb::dict` objects 3.39M times) and the Python
dict merge/emit loop, not in the scoring math.

**Implication for P1:** The C++ MM scoring can reuse the **same inline
scoring logic** already in `score_chunk_impl` (the mRNA/nRNA scoring
blocks). No need to call `score_wta_mrna()`/`score_wta_nrna()` — just
copy the inline loop. This avoids `nb::dict` construction entirely.

---

## 5. Code Deletion Inventory

After P1 is implemented and verified, the following code becomes dead
and should be removed:

### scan.py — Methods to Delete

| Lines | Method | Reason |
|---|---|---|
| 116–241 | `_score_wta_mrna` | Only called by `_add_single_fragment` (dead) and `_flush_mm_group` (replaced by C++) |
| 243–365 | `_score_wta_nrna` | Same — only used by dead/replaced paths |
| 367–397 | `_emit_mrna` | Only called by `_add_single_fragment` (dead) and `_flush_mm_group` (replaced by C++) |
| 399–419 | `_emit_nrna` | Same |
| 422–432 | `_gdna_log_lik` | Only called by `_finalize_unit_metadata` (dead) and `_flush_mm_group` (replaced by C++) |
| 435–445 | `_finalize_unit_metadata` | Only called by `_add_single_fragment` (dead) |
| 451–498 | `_add_single_fragment` | **Already dead code** — never called anywhere |
| 510–642 | `_flush_mm_group` | Replaced by C++ MM scoring in `fused_score_buffer` |

**Total: ~530 lines of Python deleted from scan.py**

### scan.py — `__init__` Attributes to Delete

The following `array.array` accumulators on `FragmentRouter` are used
exclusively by the MM Python path and the dead `_add_single_fragment`:

```python
# All of these (lines 91-114) can be deleted:
self.offsets = _array('q', [0])
self.t_indices_list = _array('i')
self.log_liks_list = _array('d')
self.count_cols_list = _array('B')
self.coverage_weights_list = _array('d')
self.tx_starts_list = _array('i')
self.tx_ends_list = _array('i')
self.locus_t_list = _array('i')
self.locus_ct_list = _array('B')
self.is_spliced_list = _array('b')
self.gdna_ll_list = _array('d')
self.genomic_footprints_list = _array('i')
self.frag_id_list = _array('i')
self.frag_class_list = _array('b')
self.splice_type_list = _array('B')

self.mm_pending: list = []
self.mm_fid: int = -1
```

### scan.py — `_scan_native` Code to Delete

| Lines | Code | Reason |
|---|---|---|
| 788–801 | MM index extraction + `__getitem__` + `_flush_mm_group` loop | Replaced by C++ |
| 817–820 | Final `_flush_mm_group` call | Replaced by C++ |
| 824 | `del chunk_arrays` | No longer needed at this point |
| 827–828 | `mm_n_units`, `mm_n_cands` computation | No MM arrays to merge |
| 830–833 | `_to_np` helper function | No array.array conversion needed |
| 835–908 | Entire `if mm_n_units > 0:` merge block + `else:` block | No merge — use C++ arrays directly |

### scan.py — Imports to Remove

```python
# These imports are only used by deleted methods:
from .scoring import frag_len_log_lik               # used by _gdna_log_lik
from .scoring import genomic_to_transcript_pos_bisect  # used by _score_wta_mrna
from .scoring import compute_fragment_weight         # used by _score_wta_mrna/_score_wta_nrna
```

### scoring.cpp — Methods to Delete

| Lines | Method | Reason |
|---|---|---|
| 305–428 | `score_wta_mrna` | Returns `nb::dict` — only called from Python MM path. The inline scoring logic in `score_chunk_impl` replaces it. |
| 430–540 | `score_wta_nrna` | Same — returns `nb::dict`, only for Python MM path |
| 542–780 | `score_emit_fragment` | **Already dead code** — only called from dead `_add_single_fragment` |

**Total: ~475 lines of C++ deleted from scoring.cpp**

### scoring.cpp — Nanobind Bindings to Delete

```cpp
// These binding declarations can be removed:
.def("score_wta_mrna", ...)     // ~8 lines
.def("score_wta_nrna", ...)     // ~6 lines
.def("score_emit_fragment", ...) // ~6 lines
```

### buffer.py — Code to Demote (Not Delete)

`BufferedFragment` dataclass and `_FinalizedChunk.__getitem__` are no
longer called on the production hot path but are still useful for:
- Test construction (`test_pipeline_routing.py` uses `_Chunk.__getitem__`)
- Debugging / interactive exploration
- The `_BF` test helper mirrors `BufferedFragment` fields

**Action:** Keep `BufferedFragment` and `__getitem__` but add a comment:
```python
# NOTE: Not on the production hot path after P1 (C++ MM scoring).
# Retained for testing and debugging.
```

### scoring.py — Code to Demote

| Function | Lines | Status |
|---|---|---|
| `frag_len_log_lik` | 254–263 | No longer called from production path (C++ has its own). Keep for tests. |
| `genomic_to_transcript_pos_bisect` | 331–378 | No longer called from production path. Keep for tests. |
| `compute_fragment_weight` (Python) | 268–329 | Already has C++ version. Keep for tests. |
| `_t_exon_data` construction | 182–195 | Still needed by NativeFragmentScorer constructor (passes dict to C++). |

### Summary of Deletions

| File | Lines Deleted | What |
|---|---:|---|
| scan.py methods | ~530 | 8 methods (scoring, emit, flush, metadata) |
| scan.py __init__ | ~25 | 17 array.array accumulators + MM state |
| scan.py _scan_native | ~90 | MM loop, merge, helper |
| scan.py imports | ~3 | Unused imports |
| scoring.cpp methods | ~475 | 3 methods (wta_mrna, wta_nrna, score_emit) |
| scoring.cpp bindings | ~20 | 3 binding declarations |
| **Total** | **~1,143** | |

---

## 6. Verification Strategy

### Phase A: C++ Implementation with Python Cross-Check

Before deleting the Python path:

1. **Add the C++ MM scoring** to `fused_score_buffer` and return MM
   results alongside the existing non-MM results.

2. **In Python `_scan_native`**, run BOTH paths:
   - New: extract MM results from the C++ return value
   - Old: run the existing Python MM loop

3. **Assert equality** between C++ and Python MM outputs:
   ```python
   # For each MM unit, verify:
   assert np.array_equal(cpp_mm_offsets, py_mm_offsets)
   assert np.array_equal(cpp_mm_t_indices, py_mm_t_indices)
   assert np.allclose(cpp_mm_log_liks, py_mm_log_liks, atol=1e-12)
   assert np.array_equal(cpp_mm_count_cols, py_mm_count_cols)
   # etc.
   ```

4. **Run full test suite** (852 tests) + profiler benchmark with
   cross-check enabled.

### Phase B: Switch to C++ Only

Once Phase A passes:

1. Remove the Python cross-check code
2. Delete all dead Python and C++ methods (per Section 5)
3. Run full test suite again
4. Run benchmark to verify performance improvement

### Phase C: Golden Output Regression

- Run `test_golden_output.py` to verify end-to-end quantification
  results are numerically identical
- Run benchmark suite to verify accuracy metrics (correlation, MARD)
  are unchanged

---

## 7. Testing Additions

### New Unit Tests for C++ MM Scoring

Add to `test_pipeline_routing.py`:

```python
def test_mm_group_mixed_splice_types():
    """MM group with ANNOT + UNSPLICED hits: nRNA gated, gDNA = -inf."""

def test_mm_group_best_overhang_wins():
    """MM group where hit 2 has lower mRNA overhang than hit 1."""

def test_mm_group_same_transcript_merge():
    """Two MM hits map to the same transcript — best score kept."""

def test_mm_group_all_gated():
    """MM group where exon_bp=0 for all candidates → gated out."""

def test_mm_group_gdna_best_unspliced():
    """gDNA log-lik picks the best among unspliced hits."""

def test_mm_group_splice_type_priority():
    """Unit splice_type = ANNOT when any hit is ANNOT."""

def test_mm_group_cross_chunk():
    """MM group spanning two chunks (same frag_id at boundary)."""
```

---

## 8. Implementation Order

```
Step 1: Add MM group tracking to score_chunk_impl       [scoring.cpp]
Step 2: Implement flush_mm_group in C++                  [scoring.cpp]
Step 3: Update fused_score_buffer return value            [scoring.cpp]
Step 4: Update Python _scan_native unpacking              [scan.py]
Step 5: Add Phase A cross-check (temporary)               [scan.py]
Step 6: Run tests + benchmark with cross-check
Step 7: Add new MM-specific unit tests                    [tests/]
Step 8: Remove Phase A cross-check                        [scan.py]
Step 9: Delete dead Python code                           [scan.py]
Step 10: Delete dead C++ code                              [scoring.cpp]
Step 11: Clean up imports, demote utility functions         [scan.py, scoring.py, buffer.py]
Step 12: Run full test suite + golden output regression
Step 13: Run profiler benchmark to measure speedup
```

---

## 9. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Cross-chunk MM group handling bug | Medium | High (silent wrong results) | Phase A cross-check; dedicated test |
| Count/fill pass disagreement | Low | High (buffer overflow / crash) | Same template handles both passes |
| Floating-point divergence (C++ vs Python) | Low | Low (cosmetic) | Use `atol=1e-12` in cross-check |
| Large MM groups (>64 candidates) | Low | Low (correctness, not perf) | Heap fallback already in pattern |
| Strand accumulation incorrectly including MM | Low | Medium | Explicit guard; test verification |
