# HulkIndex Interval Query Optimization Research

**Date:** 2026-02-28  
**Benchmark:** PVT1 locus, 100k reads, 379 transcripts, 35 genes  
**Current bottleneck:** `resolve_fragment` 12.1s self / `query_exon_with_coords` 5.9s  

---

## 1. Current Architecture

### Index structure

HulkIndex stores a single **cgranges** interval tree containing every EXON,
INTRON, and INTERGENIC interval across all transcripts. Each interval is
inserted with a unique integer label that indexes into three parallel NumPy
lookup arrays:

```
cr.add(ref, start, end, label)

_iv_t_index[label] → transcript index (int, -1 for intergenic)
_iv_g_index[label] → gene index      (int, -1 for intergenic)
_iv_type[label]    → IntervalType     (EXON=0, INTRON=1, INTERGENIC=2)
```

### Query pattern

For each aligned exon block of a fragment, `query_exon_with_coords` calls
`cr.overlap(ref, start, end)` and constructs a Python list of 5-tuples:

```python
for h_start, h_end, label in self.cr.overlap(ref, start, end):
    hits.append((
        int(self._iv_t_index[label]),   # numpy scalar → python int
        int(self._iv_g_index[label]),   # numpy scalar → python int
        int(self._iv_type[label]),      # numpy scalar → python int
        h_start,
        h_end,
    ))
```

### Downstream consumer

`resolve_fragment` iterates over these hits to build:
- `block_exon_t` / `block_intron_t` — frozensets of transcript indices per block
- `_t_exon_bp` / `_t_intron_bp` / `_t_unambig_intron_bp` — per-transcript clipped
  overlap BP accumulators (3 defaultdicts)

The `g_idx` field is **never used** in the hot loop — it's unpacked and
discarded. Gene indices are derived later from `index.t_to_g_arr[t]`.

---

## 2. Profiling Findings

### Interval table statistics (PVT1 benchmark)

| Metric | Value |
|---|---|
| Total intervals in cgranges | **3,238** |
| EXON intervals | 1,802 (55.6%) |
| INTRON intervals | 1,423 (43.9%) |
| INTERGENIC intervals | 13 (0.4%) |
| Unique (ref, start, end, type) boundaries | **993** |
| **Redundancy ratio** | **3.25×** |
| Max transcripts sharing one boundary | 97 |

**69.3% of all interval entries** are duplicates of the same genomic
boundary with a different `t_index`. They return identical `(h_start,
h_end)` coordinates from cgranges — only the label differs.

### Actual query hit distribution (from 100k BAM fragments)

| Metric | Value |
|---|---|
| Total exon block queries | 187,021 |
| **Total hits returned** | **13,683,327** |
| Mean hits per query | **73.2** |
| Median hits per query | 80 |
| P99 hits per query | 159 |
| Max hits per query | 243 |

After collapsing identical `(h_start, h_end)` boundaries:

| Metric | Value |
|---|---|
| **Collapsed hits** | **2,581,391** |
| Collapsed mean per query | **13.8** |
| **Hit reduction** | **81.1%** |

### Timing breakdown (extrapolated to full 100k fragments)

| Component | Time |
|---|---|
| Raw cgranges C library overlap | **0.65s** |
| Python tuple construction (int() + tuple + list.append) | **2.94s** |
| **Total `query_exon_with_coords`** | **5.9s** (measured by cProfile) |

The cgranges C library accounts for only **11%** of `query_exon_with_coords`
time. The remaining **89% is pure Python overhead**: 13.7M iterations of
NumPy scalar indexing, `int()` conversion, tuple construction, and
`list.append`.

### Equivalence class structure

Exon boundaries:
- 628 unique boundaries, mapping to 445 unique transcript sets
- 480 boundaries (76%) are transcript-unique (set size = 1)
- 23 boundaries have ≥10 transcripts sharing the same exon

Intron boundaries:
- 365 unique boundaries, mapping to 278 unique transcript sets
- 217 boundaries (59%) are transcript-unique

---

## 3. Candidate Optimizations

### A. Collapsed Interval Index ★★★ RECOMMENDED — HIGH IMPACT

**Concept:** Replace per-transcript intervals with per-boundary intervals.
Each unique `(ref, start, end)` is stored **once** in cgranges with a single
label. The label maps to a pre-computed `(itype, frozenset_of_t_indices)`
instead of a single `(t_index, g_index, itype)`.

**Implementation sketch:**

```python
# During load():
# Group intervals by (ref, start, end) → merge into collapsed entries
collapsed = {}  # (ref, start, end) → (itype, set_of_t_indices)
for label, row in enumerate(iv_df.itertuples(index=False)):
    key = (row.ref, row.start, row.end)
    if key not in collapsed:
        collapsed[key] = (row.interval_type, set())
    collapsed[key][1].add(row.t_index)

# Insert collapsed into cgranges (993 entries instead of 3,238)
for i, ((ref, start, end), (itype, t_set)) in enumerate(collapsed.items()):
    cr.add(ref, start, end, i)
    _collapsed_type[i] = itype
    _collapsed_t_set[i] = frozenset(t_set)  # or tuple for speed
```

```python
# New query method:
def query_collapsed(self, exon):
    hits = []
    for h_start, h_end, label in self.cr.overlap(exon.ref, exon.start, exon.end):
        hits.append((
            h_start,
            h_end,
            self._collapsed_type[label],    # int, not numpy
            self._collapsed_t_set[label],   # frozenset[int]
        ))
    return hits
```

```python
# In resolve_fragment:
for h_start, h_end, itype, t_set in index.query_collapsed(exon_block):
    if itype == IntervalType.EXON:
        block_exon_t |= t_set  # frozenset union
        bp = min(block_end, h_end) - max(block_start, h_start)
        if bp > 0:
            for t_idx in t_set:
                _t_exon_bp[t_idx] += bp
            global_exon_intervals.append(
                (max(block_start, h_start), min(block_end, h_end))
            )
    elif itype == IntervalType.INTRON:
        block_intron_t |= t_set
        bp = min(block_end, h_end) - max(block_start, h_start)
        if bp > 0:
            for t_idx in t_set:
                _t_intron_bp[t_idx] += bp
            per_t_intron_intervals_collapsed.append(...)
```

**Expected savings:**

| Component | Current | Collapsed | Savings |
|---|---|---|---|
| cgranges entries | 3,238 | 993 | 69% fewer |
| cgranges query time | 0.65s | ~0.35s | ~0.30s |
| Python iterations in query | 13.7M | 2.6M | 81% fewer |
| Tuple construction overhead | 2.94s | ~0.55s | ~2.4s |
| resolve_fragment type checks/clipping | (part of 12.1s) | 5.3× fewer | ~2s |
| **Estimated total savings** | | | **4–6s** |

**Compatibility:** The downstream API changes (returns transcript sets
instead of individual transcript tuples), requiring updates to
`resolve_fragment` and any callers of `query_exon_with_coords`.
The standalone `compute_overlap_profile` (retained for tests) would
also need updating, or we keep both query methods.

**Risks:** Low. The collapsed representation is a strict superset of the
current per-transcript representation. Correctness is easily verified — the
frozenset union in `resolve_fragment` produces identical `block_exon_t` and
`block_intron_t` sets. The per-transcript BP accumulators receive identical
values because clipped overlap is the same for all transcripts sharing a
boundary.

**Note on intron/exon overlap:** A given (ref, start, end) boundary is
always either EXON or INTRON within a single transcript. However, the
same boundary can be EXON for transcript A and INTRON for transcript B
in rare cases of overlapping genes with different exon structures. The
profiling shows that the current data has 993 unique `(ref, start, end,
type)` boundaries — the same as 993 unique `(ref, start, end)` — meaning
**no boundary is both EXON and INTRON** in the PVT1 dataset. For safety,
the implementation should group by `(ref, start, end, type)` to handle
pathological cases where exon/intron classifications differ.

---

### B. Breakpoint Partition (Sorted-Array Binary Search) ★★☆ MEDIUM IMPACT

**Concept:** Pre-compute a set of "elementary intervals" — maximal genomic
regions where the set of overlapping reference intervals is constant. These
are defined by the sorted breakpoints (starts and ends) of all intervals.

Between consecutive breakpoints, every query touching that region returns
the same result. A query `[qstart, qend)` is answered by:

1. Binary search for the first breakpoint ≥ qstart: O(log B)
2. Walk breakpoints until > qend: O(k) where k is the number of elementary
   intervals overlapping the query
3. Emit the pre-computed (type, t_set, elem_start, elem_end) for each

Where B is the number of breakpoints (at most 2 × 993 = 1,986 for the
PVT1 locus).

**Advantages:**
- Eliminates cgranges entirely — pure Python + NumPy with binary search
- O(log B + k) per query with very small constants
- Pre-computed results mean zero per-hit computation

**Disadvantages:**
- Complex implementation (breakpoint sweep, stable sorting of interval events)
- The elementary intervals may be more numerous than collapsed intervals in
  dense gene regions (every distinct start/end creates a new partition)
- Memory: must store a transcript set per elementary interval
- Requires per-chromosome arrays, adding bookkeeping

**Estimated savings:** Similar to Option A for the query phase (3–5s), but
potentially faster for the resolve_fragment loop since results are
pre-formatted. However, the implementation complexity is significantly
higher.

**Recommendation:** Defer unless Option A proves insufficient. The cgranges
C library is already highly optimized for interval overlap; the bottleneck
is the Python iteration, which Option A addresses directly.

---

### C. Eliminate int() Conversion Overhead ★★☆ COMPLEMENTARY

**Concept:** The current `query_exon_with_coords` does three `int(numpy_scalar)`
conversions per hit. Storing the lookup arrays as Python lists of native ints
(instead of NumPy arrays) eliminates the numpy-to-Python conversion:

```python
# During load():
self._iv_t_index = iv_df["t_index"].values.tolist()  # list[int]
self._iv_g_index = iv_df["g_index"].values.tolist()  # list[int]
self._iv_type = iv_df["interval_type"].values.tolist()
```

This makes `self._iv_t_index[label]` return a native Python int directly.

**Savings:** ~1.5s (estimated from 13.7M × 3 × ~35ns per int() conversion).
If combined with Option A, this becomes less impactful (only 2.6M
iterations), saving ~0.3s.

**Recommendation:** Apply as a quick win if Option A is not implemented.
With Option A, the lookup arrays become `_collapsed_type` (Python list of
ints) and `_collapsed_t_set` (Python list of frozensets), which are
natively Python — the issue disappears.

---

### D. Transcript Span Model ★☆☆ LOW IMPACT

**Concept:** Instead of storing separate INTRON intervals per transcript,
store one "transcript span" `[tx_start, tx_end)` per transcript. To
determine whether a query hit is exonic or intronic, check the exon
intervals (already cached in `_t_exon_intervals`).

This reduces genic intervals from 3,225 to 2,181 (32% fewer = 1,044 fewer
intervals). However, it introduces a per-transcript exon-membership lookup
for every span hit.

**Analysis:** The 32% interval reduction is modest compared to Option A's 69%.
Moreover, the per-hit exon-membership check adds O(log n_exons) overhead.
The combination does not clearly win.

**Recommendation:** Not recommended as a standalone change. Option A
subsumes most of this benefit through boundary collapsing, which already
eliminates the intron interval redundancy.

---

### E. Pre-computed Per-Query Caching ★☆☆ LOW APPLICABILITY

**Concept:** Cache query results so that fragments with identical exon
blocks reuse previous results.

**Analysis:** In the PVT1 benchmark with 100k oracle-aligned fragments
landing on 379 transcripts, exon block coordinates are highly diverse
(different fragment positions within each transcript). Cache hit rates
would be very low. This approach works better for k-mer or pseudo-alignment
pipelines where equivalence classes are defined by a small alphabet.

**Recommendation:** Not applicable to the current alignment-based pipeline.

---

### F. Restructured resolve_fragment Loop ★★☆ COMPLEMENTARY TO A

**Concept:** Once Option A provides collapsed hits with frozensets, the
resolve_fragment inner loop can be further optimized:

1. **Drop g_idx from query:** The gene index is unused in the hot loop.
   Already handled by Option A's new return format.

2. **Batch BP accumulation:** Instead of per-transcript `_t_exon_bp[t_idx]
   += bp` in a Python loop, accumulate into a NumPy array:

   ```python
   # Pre-allocate: _bp_accum = np.zeros(n_transcripts, dtype=np.int32)
   # Per collapsed hit:
   t_arr = np.array(list(t_set), dtype=np.int32)
   np.add.at(_bp_accum, t_arr, bp)
   ```

   This avoids Python-level dict operations for the per-transcript
   accumulation.

3. **Pre-compute clipped intervals as NumPy operations:** Collect all
   `(h_start, h_end)` from the query, then vectorize `max/min` clipping:

   ```python
   h_starts = np.array([h[0] for h in hits])
   h_ends   = np.array([h[1] for h in hits])
   clipped_lo = np.maximum(block_start, h_starts)
   clipped_hi = np.minimum(block_end, h_ends)
   bp = np.maximum(clipped_hi - clipped_lo, 0)
   ```

**Estimated savings:** 1–2s on top of Option A.

---

## 4. Unambiguous Intron BP Handling

The `_t_unambig_intron_bp` computation requires special attention.
Currently, for each block:

1. Collect all EXON hit intervals (across ALL transcripts) → `global_exon_intervals`
2. Merge into non-overlapping union
3. For each transcript's INTRON hit, subtract the exon union to get
   "unambiguously intronic" bases

With collapsed intervals, step 1 becomes simpler (fewer iterations), but
step 3 still requires per-transcript processing because different
transcripts may have different intron boundaries even if some share the
same boundary.

**Key insight:** With collapsing, the `per_t_intron_intervals` dict can be
built directly from the collapsed t_set. All transcripts in a collapsed
intron hit share the *same* `(h_start, h_end)`, so `_subtract_intervals`
returns the same value for all of them:

```python
# Before (per-transcript):
for t_idx, intron_ivs in per_t_intron_intervals.items():
    for iv in intron_ivs:
        _t_unambig_intron_bp[t_idx] += _subtract_intervals(iv, merged_exon)

# After (collapsed):
for (h_start, h_end), t_set in collapsed_intron_hits:
    unambig_bp = _subtract_intervals((h_start, h_end), merged_exon)
    if unambig_bp > 0:
        for t_idx in t_set:
            _t_unambig_intron_bp[t_idx] += unambig_bp
```

The `_subtract_intervals` call is done **once per collapsed boundary**
instead of once per transcript — saving up to 97× for the most redundant
boundaries.

---

## 5. Recommendations & Implementation Roadmap

### Priority 1: Collapsed Interval Index (Option A)

**Impact: HIGH (estimated 4–6s savings, 10–14% of pipeline)**

This is the single highest-impact optimization available for the interval
query system. It addresses the root cause — 81% of cgranges hit iterations
are redundant — and reduces Python overhead proportionally.

Implementation steps:
1. Add `_build_collapsed_index()` method to HulkIndex.load()
2. Create `query_collapsed()` returning `(h_start, h_end, itype, t_set)` tuples
3. Update `resolve_fragment` to consume collapsed hits
4. Update `compute_overlap_profile` (test-only path)
5. Verify all 654 tests pass and benchmark MAE unchanged

### Priority 2: Restructured resolve_fragment Loop (Option F)

**Impact: MEDIUM (estimated 1–2s additional)**

Once collapsed hits are available, vectorize the BP accumulation and
clipping in resolve_fragment with NumPy operations.

### Priority 3: Native Python lookup arrays (Option C)

**Impact: LOW (subsumed by Option A)**

If Option A is implemented, this comes for free. If not, it's a quick
standalone ~1.5s win.

### Not recommended:

- **Option B (breakpoint partition):** Higher complexity for similar gains.
  Worth revisiting only if Python cgranges becomes a bottleneck after
  Option A.
- **Option D (transcript spans):** Marginal benefit, adds complexity.
- **Option E (query caching):** Low hit rate for alignment-based pipelines.

---

## 6. Projected Impact

Starting from current baseline: **43.4s** (after P1 + P2)

| Optimization | Expected savings | Running total |
|---|---|---|
| P3: Uniform bias fast-path | ~10s | ~33s |
| Collapsed interval index (A) | 4–6s | ~28s |
| Restructured resolve loop (F) | 1–2s | ~27s |
| **Combined** | **15–18s** | **~26–28s** |

With all three: **1.6–1.7× additional speedup** from the current 43.4s,
or **~2.0× total speedup** from the original 52.7s baseline.

This brings the pipeline into the ~26s range for 100k fragments with pure
Python/NumPy — no C extensions or Cython required.

---

## Appendix: Diagnostic Data

### Interval type breakdown

```
EXON           :    1,802 intervals  median=162bp  mean=265bp
INTRON         :    1,423 intervals  median=11,457bp  mean=26,677bp
INTERGENIC     :       13 intervals  median=28,748bp
```

### Per-transcript interval counts

```
Transcripts:           379
Intervals/transcript:  median=9, mean=8.5, max=25
Exons/transcript:      median=5, mean=4.8, max=13
Introns/transcript:    median=4, mean=4.0, max=12
Single-exon tx:        23  (6.1%)
Multi-exon tx:         356 (93.9%)
```

### Boundary redundancy

```
Total genic intervals:        3,225
Unique (ref,start,end):         993
Redundancy:                   3.25×
Max transcripts at one boundary: 97
```

### Actual query patterns (20k BAM fragment sample, extrapolated to 100k)

```
Queries:                187,021
Total hits:          13,683,327
Mean hits/query:           73.2
Median hits/query:           80
P99:                        159
Max:                        243

After collapsing identical boundaries:
Collapsed hits:       2,581,391
Mean collapsed/query:      13.8
Hit reduction:            81.1%
```

### Timing (extrapolated to 100k fragments)

```
cgranges C overlaps:     0.65s  (11% of query_exon_with_coords)
Python tuple overhead:   2.94s  (50%)
Other Python overhead:   2.31s  (39%)
Total:                   5.90s
```
