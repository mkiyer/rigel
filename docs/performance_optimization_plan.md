# Performance Optimization Plan for hulkrna

## Profiling Summary

**Dataset:** Real BAM — 2,437,978 fragments, 254,461 transcripts, 63,472 genes  
**Wall time:** 16.37s pipeline + 51.5s index load = **67.9s total**  
**Throughput:** 148,904 fragments/sec  

### Phase Breakdown (Pipeline Only — 16.37s)

| Phase | Time (s) | % | Description |
|-------|----------|---|-------------|
| `quant_from_buffer` | 12.33 | 75.3% | Scoring + EM quantification |
| `scan_and_buffer` | 4.00 | 24.4% | C++ BAM scan + model replay |
| Overhead | 0.05 | 0.3% | Pipeline glue code |

### Bottleneck Functions (by self time)

| Rank | Function | Self (s) | Cum (s) | Calls | Description |
|------|----------|----------|---------|-------|-------------|
| 1 | `scan_and_buffer` | 3.25 | 4.00 | 1 | C++ BAM scan + Python model replay |
| 2 | `run_locus_em` | 1.35 | 1.35 | 6,670 | C++ EM solver dispatch |
| 3 | `buffer.__getitem__` | 1.03 | 1.08 | 401,556 | Python object creation per fragment |
| 4 | `scan` (FragmentRouter) | 0.92 | 5.42 | 1 | Python fragment iteration loop |
| 5 | `from_models` | 0.85 | 1.55 | 1 | FragmentScorer initialization |
| 6 | `_add_single_fragment` | 0.80 | 1.47 | 220,734 | Per-fragment scoring dispatch |
| 7 | `build_locus_em_data` | 0.62 | 1.90 | 6,670 | Per-locus numpy array construction |
| 8 | `ArrowArray.__getitem__` | 0.60 | 0.89 | 414,241 | Pandas Arrow indexing overhead |
| 9 | `array.append` | 0.55 | 0.55 | 5,399,703 | CSR list building |
| 10 | `_flush_mm_group` | 0.52 | 1.75 | 72,237 | Multimapper Python scoring |

---

## Optimization Opportunities

### Priority 1: Eliminate Python Per-Fragment Loop in `quant_from_buffer` (est. 4-5s savings)

**Current cost: ~7.5s** (scan + __getitem__ + _add_single_fragment + accumulation + _flush_mm_group)

The `FragmentRouter.scan()` iterates every fragment in Python:
```
for chunk in buffer.iter_chunks():
    for i in range(chunk.size):  # 401K iterations
        bf = chunk[i]    # creates BufferedFragment dataclass: 1.03s
        ...dispatch...
```

Each `chunk[i]` constructs a `BufferedFragment` dataclass with 17 fields, slicing
multiple numpy arrays into Python objects. Then `_add_single_fragment` calls into C++
via `score_emit_fragment` — but the per-call overhead of creating the Python object
and crossing the Python→C++ boundary 220K times is enormous.

**Fix:** Move the entire `scan()` loop into C++/nanobind. The C++ function would:
1. Accept the entire chunk's columnar arrays at once (already numpy).
2. Iterate fragments within C++, calling the scoring+emit logic.
3. Return the assembled CSR arrays and accumulator counters.

This eliminates:
- 401K `BufferedFragment` dataclass constructions (1.03s)
- 220K Python→C++ crossing overhead in `_add_single_fragment` (0.80s)
- 72K `_flush_mm_group` Python iterations with per-fragment scoring (0.52s)
- Per-fragment Python branching in `scan()` (0.92s)
- Millions of `array.append` calls (0.55s)
- Pre-EM accumulation (exonic/intronic/unspliced sense/antisense) Python loops (embedded in scan)

The fragment classification (`fragment_classes`) is already vectorized. The entire
loop should be a single C++ call per chunk.

### Priority 2: Vectorize `from_models` Exon Data Pre-computation (est. 1.0s savings)

**Current cost: 1.55s** (0.85s self, 254K `get_exon_intervals` calls)

```python
for t_idx in range(index.num_transcripts):   # 254K iterations
    exon_ivs = index.get_exon_intervals(t_idx)
    starts = tuple(int(x) for x in exon_ivs[:, 0])
    ends = tuple(int(x) for x in exon_ivs[:, 1])
    ...cumsum...
    t_exon_data[t_idx] = (starts, ends, tuple(cb))
```

This Python loop creates 254K dict entries with Python tuples, calling into the
index's dict-based exon storage 254K times.

**Fix A (quick):** Store exon interval data in a flat CSR (offsets + starts + ends) in
the index itself (computed once at index build time or at load). Pass the flat arrays
to C++ NativeFragmentScorer instead of a Python dict.

**Fix B (deeper):** Move this pre-computation into the C++ NativeFragmentScorer
constructor, accepting the CSR arrays directly.

### Priority 3: Accelerate `build_locus_em_data` (est. 0.5-1.0s savings)

**Current cost: 1.90s cumulative** across 6,670 loci.

Each locus call does:
- `np.full` for local_map array (up to 508K elements)
- Vectorized gather + dedup + gDNA concatenation
- `np.argsort`, `np.lexsort`, `np.concatenate` per locus

**Fix A:** Reuse pre-allocated scratch buffers across loci instead of allocating
fresh numpy arrays per locus (the `local_map` array is re-created every call at
`max_global` size which can be ~500K elements).

**Fix B:** Move the locus EM data extraction into C++ — given that `run_locus_em`
is already C++, the entire locus loop (build_locus_em_data + run_locus_em +
assign_locus_ambiguous) could be a single C++ call over all loci.

### Priority 4: Accelerate `scan_and_buffer` Model Replay (est. 0.7s savings)

**Current cost:** 0.75s in `_replay_strand_observations` (0.32s) +
`_replay_fraglen_observations` (0.42s)

These functions iterate Python lists element-by-element:
```python
for o, t in zip(obs, truth):
    strand_models.exonic_spliced.observe(Strand(int(o)), Strand(int(t)))
```

The `Strand(int(o))` call triggers enum construction (0.20s for 744K `Enum.__call__`).

**Fix A (quick):** Replace `Strand(int(o))` with direct int passing — the `observe()`
methods only need integers internally. Add `observe_batch()` methods that accept
numpy arrays and update counts directly.

**Fix B (deeper):** Move strand/fraglen model training entirely into C++ scanner,
returning only the finalized model parameters (counts, mean, var) rather than
raw observations.

### Priority 5: Optimize Index Loading (est. 30-40s savings)

**Current cost: 51.5s** to load TranscriptIndex.

This is by far the largest wall-time cost. While it's a one-time cost, for
interactive use it dominates the experience.

**Fix A:** Profile `TranscriptIndex.load()` separately to identify whether the
bottleneck is Feather I/O, interval tree construction, or DataFrame operations.

**Fix B:** Cache the constructed interval tree / ContainmentIndex in a binary
format (e.g., pickle or a custom binary format) that can be memory-mapped.

**Fix C:** Lazy-load components that aren't needed for all runs (e.g., SJ data
only needed when SJ tag detection is required).

### Priority 6: Reduce Pandas Overhead in Locus Loop (est. 0.5s savings)

**Current cost:** 0.67s in `pandas DataFrame.__getitem__` (33,646 calls) +
associated Arrow array indexing and Series operations.

The locus EM loop accesses `index.t_df["start"]`, `index.t_df["end"]`,
`index.t_df["ref"]` repeatedly per locus.

**Fix:** Pre-extract these as numpy arrays once before the locus loop
(some are already extracted, but `t_refs` uses `index.t_df["ref"].values`
inside the loop). The `build_locus_em_data` function accesses
`index.t_df["start"].values[t_arr]` and `index.t_df["end"].values[t_arr]`
per locus — these `.values` lookups go through Pandas each time.

---

## Implementation Roadmap

### Phase A: Quick Wins (est. 2-3s savings, ~1 day)

1. **Vectorize model replay** — Replace element-by-element `observe()` with
   `observe_batch(array)` methods that update internal counts via numpy
   vectorized operations. Eliminate `Strand(int(o))` enum construction.

2. **Pre-extract DataFrame columns** — Extract `index.t_df["start"].values`,
   `["end"].values`, `["ref"].values`, `["length"].values` as numpy arrays
   before entering the locus loop, and pass them as arguments to
   `build_locus_em_data`. Same for `quant_from_buffer`.

3. **Reuse `local_map` buffer** — Allocate once at max size and zero/reset
   per locus instead of `np.full(max_global, -1)` per locus.

### Phase B: Fuse Scoring Loop into C++ (est. 4-5s savings, ~3-5 days)

1. **Extend `NativeFragmentScorer`** to accept entire chunk columnar arrays
   and perform the full scan loop in C++, returning packed CSR arrays +
   accumulator counters.

2. **Move multimapper flush logic** into C++ (currently 72K Python-side
   `_flush_mm_group` calls with nested scoring).

3. **Move pre-EM accumulation** (exonic/intronic/unspliced sense/antisense
   counting) into the C++ scan loop.

### Phase C: Fuse Locus Loop into C++ (est. 1-2s savings, ~2-3 days)

1. **Move `build_locus_em_data` + `run_locus_em` + `assign_locus_ambiguous`**
   into a single C++ call that iterates all loci, eliminating per-locus
   numpy allocation overhead and Python dispatch.

### Phase D: Index Load Optimization (est. 30-40s savings, ~2-3 days)

1. **Profile index loading** to identify specific bottleneck.
2. **Serialize interval tree** to binary format for fast reload.
3. **Lazy-load** SJ data and exon intervals.

---

## Expected Impact

| Phase | Savings (s) | New Pipeline (s) | New Total (s) | Effort |
|-------|-------------|-------------------|----------------|--------|
| Baseline | — | 16.4 | 67.9 | — |
| Phase A | 2-3 | 13-14 | 65-66 | 1 day |
| Phase B | 4-5 | 8-10 | 60-62 | 3-5 days |
| Phase C | 1-2 | 7-8 | 59-60 | 2-3 days |
| Phase D | 30-40 | 7-8 | 19-28 | 2-3 days |
| **All** | **~40-50** | **~7-8** | **~19-28** | **~9-12 days** |

Pipeline throughput improvement: **148K → ~300-350K frags/sec** (Phases A-C)  
Total user-facing improvement: **68s → ~20-28s** (with Phase D)
