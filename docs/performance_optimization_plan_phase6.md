# Phase 6: Python → C++ Migration & Memory Optimization Plan

## Profiling Baseline

**Dataset:** `gdna_high_ss_0.90_nrna_low` minimap2 alignment (worst-case: high
gDNA, low nRNA, 29.4M fragments, 254K transcripts, 63K genes).

**Platform:** Apple M4 Max, macOS 26.3, Python 3.12.13, 8 threads.

### Stage Breakdown (wall clock)

| Stage | Time (s) | % | Notes |
|---|---:|---:|---|
| `scan_and_buffer` (C++ BAM scan) | 76.1 | 30.6% | Already fully in C++ |
| `fragment_router_scan` (scoring + routing) | 57.5 | 23.1% | C++ fused scorer + Python MM loop |
| `locus_em` (C++ batch EM) | 110.0 | 44.3% | Already fully in C++ |
| `eb_gdna_priors` | 2.9 | 1.2% | Python |
| `fragment_scorer` (setup) | 1.5 | 0.6% | Python one-shot setup |
| `build_loci` | 0.4 | 0.1% | C++ connected_components |
| Other | 0.1 | <0.1% | |
| **Total** | **248.5** | | |

### Memory

| Metric | Value |
|---|---|
| RSS before pipeline | 3,851 MB |
| RSS after pipeline | 22,136 MB |
| RSS delta | 18,284 MB |
| Peak RSS | 22,136 MB |

### cProfile Top Python Hotspots

All times are **cumulative** (include callees) unless marked **(self)**.

| Rank | Function | Self (s) | Cum (s) | Calls | Issue |
|---|---|---:|---:|---:|---|
| 1 | `scan.py:_flush_mm_group` | 9.9 | 61.1 | 885K | Pure Python MM loop, drives all below |
| 2 | `buffer.py:__getitem__` | 9.3 | 10.6 | 3.39M | Creates `BufferedFragment` dataclass per row |
| 3 | `{array.array.append}` | 9.6 | 9.6 | 108M | CSR accumulation for MM path |
| 4 | `scan.py:_emit_nrna` | 6.2 | 13.5 | 731K | Per-dict-item CSR append loop |
| 5 | `scan.py:_score_wta_nrna` | 3.0 | 3.0 | 2.88M | Python nRNA scoring (MM hits) |
| 6 | `scan.py:_score_wta_mrna` | 2.2 | 2.2 | 3.39M | Python mRNA scoring (MM hits) |
| 7 | `{dict.get}` | 1.9 | 1.9 | 26.7M | gdna_splice_penalties + merge lookups |
| 8 | `scan.py:_gdna_log_lik` | 1.0 | 2.5 | 2.36M | Per-fragment gDNA log-lik (MM path) |
| 9 | `scan.py:_emit_mrna` | 1.4 | 3.0 | 885K | Per-dict-item CSR append loop |
| 10 | `scoring.py:frag_len_log_lik` | 0.4 | 0.4 | 2.36M | Python LUT lookup (MM path) |

### Key Observations

1. **The multimapper (MM) Python loop dominates.** `_flush_mm_group` is
   called 885K times and collectively accounts for ~57.5s of the 57.5s
   `fragment_router_scan` stage. The C++ `fused_score_buffer` handles
   non-MM fragments; MM fragments fall back to per-row Python.

2. **`buffer.py:__getitem__` creates 3.4M `BufferedFragment` dataclasses**
   solely to serve the MM Python path — each call constructs 17 fields
   from NumPy column indexing. This is pure overhead that vanishes when
   MM scoring moves to C++.

3. **CSR accumulation for MMs** uses `array.array.append` (108M calls, 9.6s).
   The C++ fused scorer already pre-allocates and fills arrays in two passes;
   the MM path should do the same.

4. **`locus_em` at 110s (44.3%)** is entirely in C++ (`batch_locus_em`).
   The Python overhead is only 0.083s for the wrapper call. Further gains
   require C++ algorithmic improvements (e.g., mega-locus handling, SQUAREM
   acceleration tuning, memory layout), not Python changes.

5. **`scan_and_buffer` at 76s (30.6%)** is entirely in C++. No Python work.

6. **Memory: 22 GB peak RSS** is very high. The CSR arrays for 15M+ scored
   fragments × 16 parallel arrays dominate. Spilling helped (1 chunk spilled
   to Arrow), but the in-memory EM data is still enormous. Opportunities:
   - Compact CSR column types (e.g., float32 log-liks, int16 tx positions)
   - Release buffer immediately after fused_score_buffer completes
   - Streaming EM that processes loci without holding all CSR data

---

## P1: Move Multimapper Scoring to C++ [HIGH — ~55s savings]

**Current state:** `fused_score_buffer` handles `FRAG_UNAMBIG`,
`FRAG_AMBIG_SAME_STRAND`, and `FRAG_AMBIG_OPP_STRAND` fragments in C++.
`FRAG_MULTIMAPPER` fragments are extracted per-fragment via
`chunk[idx]` (creating a `BufferedFragment`) and scored through the Python
`_flush_mm_group → _score_wta_mrna → _score_wta_nrna → _emit_mrna →
_emit_nrna → _gdna_log_lik` call chain.

**Problem:** 885K multimapper groups × ~3.8 alignments/group = 3.39M
individual fragment scores, all in pure Python. This accounts for
~55s of wall time (the entire `fragment_router_scan` stage minus C++
`fused_score_buffer` overhead).

**Plan:**

### P1a: C++ Multimapper Scoring Kernel

Add a `score_multimapper_group()` method to `NativeFragmentScorer` in
`scoring.cpp` that:

1. Accepts CSR-style multimapper group data: arrays of `(t_offsets, t_indices,
   exon_bp, intron_bp, frag_lengths, exon_strand, splice_type, nm,
   read_length, genomic_footprint, genomic_start)` for all alignments
   in one MM group.
2. Internally performs:
   - Per-alignment mRNA WTA scoring → merge across alignments (best per
     transcript, global WTA gate)
   - Per-alignment nRNA WTA scoring → merge (same)
   - gDNA log-lik computation (best unspliced hit)
   - Output: same `(offsets, t_indices, log_liks, count_cols, ...)` arrays
     as the fused scorer produces
3. Returns pre-sized numpy arrays (two-pass or single-pass with
   stack-allocated scratch).

### P1b: Integrate MM Groups into Fused Buffer Pass

Currently the fused scorer skips MM fragments and the Python loop
picks them up afterward. Better approach:

1. **Pre-group** MM fragments in C++ during the fused pass 1 (count phase).
   Track `(frag_id, chunk_offset)` pairs in a side buffer.
2. In pass 2, score MM groups inline using the kernel from P1a.
3. Eliminate the Python MM loop entirely — _flush_mm_group, __getitem__,
   _emit_mrna, _emit_nrna, _gdna_log_lik all become dead code on the
   hot path.

**Complexity:** Medium-high. The MM merge logic (best-per-transcript,
global WTA, splice priority) has several invariants that must be
carefully replicated.

**Expected savings:** ~55s → ~1-2s (C++ loop over 885K groups × 3.8
hits/group on M4 Max should be negligible). Eliminates ~200M Python
function calls.

**Memory savings:** Eliminates `array.array` CSR accumulators for MM
path (~500 MB for this dataset). Pre-allocated numpy arrays from fused
pass replace them.

---

## P2: Eliminate `BufferedFragment` / `__getitem__` Overhead [HIGH — folded into P1]

**Current state:** `__getitem__` is called 3.39M times (10.6s cumulative),
constructing a 17-field `@dataclass(slots=True)` from NumPy column indexing.
Every call does 17 `int()` conversions + 6 slice operations + 1 dataclass
`__init__`.

**This is eliminated automatically by P1** — once MM scoring moves to C++,
there is no per-fragment Python access. The `__getitem__` path becomes
test-only / debug-only.

If P1 is delayed, a cheaper interim fix:
- Replace `BufferedFragment` with a `namedtuple` (faster __init__)
- Batch `__getitem__` calls: add `__getitems__(indices)` that returns
  a columnar sub-chunk for all MM hits at once, avoiding per-row overhead

---

## P3: Compact CSR Column Types [MEDIUM — ~3-5 GB RSS savings]

**Current state:** The CSR arrays from `fused_score_buffer` use:
- `log_liks`: float64 (8 bytes/candidate)
- `coverage_weights`: float64 (8 bytes/candidate)
- `tx_starts`, `tx_ends`: int32 (4 bytes/candidate)
- `count_cols`: uint8 (1 byte/candidate)
- `t_indices`: int32 (4 bytes/candidate)
- `offsets`: int64 (8 bytes/unit)

For 15M units and ~30M+ candidates, this is ~750 MB+ for the EM data alone.

**Opportunities:**

### P3a: float32 log-likelihoods

Log-likelihoods are relative scores — float32 has ~7 decimal digits,
far more than needed for WTA gating and EM posterior computation.
Savings: 4 bytes/candidate × 30M = ~120 MB.

Requires: Change C++ output type, modify EM solver to accept float32
input (or upcast at locus level — the per-locus float arrays are small).

### P3b: float32 coverage weights

Coverage weights are in [0, 1]. float32 is sufficient.
Savings: ~120 MB.

### P3c: int16 transcript-relative positions

`tx_starts` and `tx_ends` encode transcript-relative positions. Most
transcripts are <32,767 bp (int16 range). Use int16 with an overflow
sentinel for the rare long transcripts.
Savings: ~120 MB.

### P3d: Packed count_cols

`count_cols` is 0–5 (6 values, 3 bits). Currently uint8. Could pack
4 values per byte, but complexity may not justify the savings.
Skip unless memory is critical.

**Total potential savings:** ~360 MB on EM data for this dataset.
Larger datasets (100M+ fragments) would see proportionally more.

---

## P4: Early Buffer Release [MEDIUM — ~3-4 GB RSS savings]

**Current state:** The `FragmentBuffer` (spilled chunks) is loaded into
memory by `iter_chunks()` and held until the `scan()` method returns,
even though `fused_score_buffer` has already consumed all chunk data.

**Plan:** After `fused_score_buffer` returns and MM groups are extracted,
immediately:
1. Delete `chunk_arrays` (already done).
2. Call `buffer.release()` to free spilled chunks.
3. Delete chunk objects.

The profiler shows `buffer.release()` is called after `scan()` returns
in `profile_stages()`, but the main `quant_from_buffer` path may
hold the buffer longer than necessary.

**Expected savings:** 3-4 GB (the spilled ~15M fragment chunk is ~3.3 GB).

---

## P5: `FragmentScorer.from_models` Build Cost [LOW — 1.5s]

**Current state:** `from_models` takes 1.5s, mostly in building
`_t_exon_data` (a Python dict of tuples for 254K transcripts). This calls
`index.get_exon_intervals` 254K times (0.053s for the calls, but tuple
construction and dict insertion add up).

**Plan:** Move `_t_exon_data` construction to C++ as part of the
`NativeFragmentScorer` constructor. The C++ side already receives
the transcript arrays — it can build an internal flat exon lookup
table without Python dict/tuple overhead.

If the native scorer already has the exon data from index load, the
Python `_t_exon_data` dict for the Python fallback scoring path can be
eliminated entirely.

**Expected savings:** ~1.0-1.3s.

---

## P6: `_compute_intergenic_density` [LOW — 0.15s]

**Current state:** Reads `ref_lengths.feather`, then iterates
`index.t_df.groupby("ref")` to compute genic union span via a
Python merger loop. Takes 0.15s.

**Plan:** Pre-compute genic union span at index build time and store it
as a scalar in the index metadata. At quant time, just read the
pre-computed value. Alternatively, add a C++ helper that computes the
genic union in a single pass over pre-sorted transcript coordinates.

**Expected savings:** ~0.14s (minor but clean).

---

## P7: `compute_eb_gdna_priors` + `compute_hybrid_nrna_frac_priors` [LOW — 2.9s]

**Current state:** 2.9s in Python, mostly in NumPy vectorized operations
over 254K transcripts × 8.7K loci. Functions include hierarchical MoM
estimation, per-locus rate computation, and Beta prior assignment.

**Plan:** Not worth porting to C++ unless the transcript/locus count
grows significantly. The code is clean NumPy vectorized — the overhead
is startup cost (groupby, slicing), not inner-loop Python.

**If needed:** Bundle into the C++ `batch_locus_em` entrypoint so that
prior computation and EM run in one C++ call, eliminating the Python
locus-iteration overhead.

---

## P8: Locus EM C++ Algorithmic Improvements [OPTIONAL — 110s]

The locus EM is already entirely in C++ but consumes 44% of total time.
This is algorithmic complexity, not Python overhead. Potential C++
improvements (separate from Python cleanup):

### P8a: Mega-locus SQUAREM acceleration

The single 217K-transcript, 10.6M-unit mega-locus dominates. SQUAREM
is already used but the budget_divisor tuning and convergence detection
could be improved for this extreme case.

### P8b: Memory layout for E-step cache efficiency

The parallel E-step accesses `log_weights[t_indices[j]]` with
scattered indices. Reordering equivalence classes by connected component
and renumbering transcript indices to be dense within each locus would
improve L2/L3 cache hit rates.

### P8c: Sparse prior optimization

Many transcripts in the mega-locus have negligible posterior mass.
A sparse-theta representation that prunes near-zero components during
EM iterations could dramatically reduce E-step cost.

---

## P9: Streaming EM / Reduce Peak Memory [OPTIONAL — target 50% RSS reduction]

**Current state:** All 15M+ scored fragments are held in memory as
dense CSR arrays throughout the EM. Peak RSS is 22 GB.

**Plan (long-term):** Process loci in streaming fashion:
1. Sort scored fragments by locus during `fused_score_buffer` (group
   candidates by connected component).
2. Process one locus at a time: extract local CSR → run EM → write
   posteriors → free local arrays.
3. Only global accumulators (per-transcript counts) stay resident.

This would reduce peak memory from O(total_candidates) to
O(max_locus_candidates) + O(transcripts).

**Complexity:** High. Requires significant refactoring of the CSR
build + locus EM pipeline.

---

## Priority Summary

| ID | Description | Savings | Complexity | Priority |
|---|---|---|---|---|
| **P1** | MM scoring to C++ | ~55s (22%) | Medium-high | **HIGH** |
| **P2** | Eliminate BufferedFragment | Folded into P1 | — | **HIGH** |
| **P3** | Compact CSR types | ~360 MB | Low | **MEDIUM** |
| **P4** | Early buffer release | ~3-4 GB | Low | **MEDIUM** |
| **P5** | from_models to C++ | ~1.3s | Low | LOW |
| **P6** | Intergenic density precompute | ~0.14s | Low | LOW |
| **P7** | EB priors to C++ | ~2.5s | Medium | LOW |
| **P8** | Locus EM C++ algorithmic | Variable | High | OPTIONAL |
| **P9** | Streaming EM | ~11 GB RSS | High | OPTIONAL |

## Implementation Order

1. **P1** (MM to C++) — Single biggest win. Eliminates nearly all
   remaining Python hot-path code, ~200M function calls, and ~55s wall time.
2. **P4** (early buffer release) — Quick RSS win, low risk.
3. **P3** (compact types) — Straightforward type narrowing.
4. **P5** (from_models) — Clean up Python scoring setup.
5. **P6–P7** — Marginal gains, nice-to-have.

## Post-P1 Projected Profile

After P1 (MM to C++), the Python contribution to the hot path
drops to effectively zero:

| Stage | Current (s) | Projected (s) | Notes |
|---|---:|---:|---|
| scan_and_buffer | 76.1 | 76.1 | Unchanged (C++) |
| fragment_router_scan | 57.5 | ~2-3 | fused_score_buffer only (C++) |
| locus_em | 110.0 | 110.0 | Unchanged (C++) |
| Other | 5.0 | 3.5 | P5 saves ~1.5s |
| **Total** | **248.5** | **~192** | **~23% faster** |

The remaining optimization frontier is then purely C++ algorithmic
(P8: mega-locus EM acceleration).
