# Phase 3 Implementation Plan: Early Pipeline (Scan & Buffer) Performance

**Date:** 2026-03-24  
**Status:** Approved Plan  
**Prerequisite:** Phases 1–2 complete (calibration + medium-impact optimizations)  
**Goal:** Reduce scan+buffer wall time by 30–40%, reduce peak RSS by 5–10 GB  

---

## Executive Summary

The scan-and-buffer stage is the dominant bottleneck in the rigel pipeline, consuming 43% of wall time (110s) and driving peak RSS to ~21 GB. This plan details a phased implementation to:

1. **Phase 3a** — Quick wins: expose BGZF decompression threads, pre-allocate accumulators, lazy-load index components. (**~10–20s wall time savings**)
2. **Phase 3b** — Zero-copy finalization and merge elimination. (**~8–15s wall time, ~1.3 GB memory savings**)
3. **Phase 3c** — Streaming chunk architecture: workers emit bounded chunks, Python spills while scan continues. (**~10–20s wall time, ~5–10 GB peak RSS savings**)
4. **Phase 3d** — Router single-pass scoring and pipeline overlap. (**~5–9s wall time savings**)

**Projected total savings: 35–65s wall time (30–50% of scan+router), 6–11 GB peak RSS.**

---

## Current Architecture (Baseline)

### Data Flow

```
BAM File (htslib + 2 BGZF threads)
    ↓ sam_read1()
[Reader: main thread]
    parse_bam_record() → ParsedAlignment
    Group by qname
    ↓
[BoundedQueue<QnameGroup>]  capacity = 2 × n_workers
    ↓
[N Worker Threads]
    Each has WorkerState:
      - ResolverScratch
      - FragmentAccumulator   ← grows unbounded
      - BamScanStats
      - StrandObservations
      - FragLenObservations
      - RegionAccumulator
    ↓
[Main thread: merge N accumulators → 1]  ← O(total_data) copy
    ↓
[scanner.finalize_accumulator()]
    to_bytes() copies all vectors → nb::bytes  ← copy #1
    ↓
[Python: _FinalizedChunk.from_raw()]
    np.frombuffer().copy() on each array  ← copy #2
    ↓
[FragmentBuffer.inject_chunk()]
    Single chunk in memory (no spill for typical workloads)
```

### Key Numbers (simulated 14.7M fragment workload)

| Metric | Value |
|--------|-------|
| BAM records | 38.7M |
| Fragment groups | 14.7M |
| Buffered units (CSR expansion) | 17.7M |
| Per-fragment storage | ~87 bytes (37 fixed + 50 CSR) |
| Accumulator data | ~1.3 GB |
| Peak RSS (scan phase) | ~17.3 GB |
| Peak RSS (total) | ~20.8 GB |
| Scan wall time (default threads) | 109.8s |
| Router wall time | 18.1s |

### Bottleneck Analysis

| Sub-stage | Estimated Time | Notes |
|-----------|---------------|-------|
| BAM I/O + BGZF decompress | 25–35s | Only 2 decompression threads |
| Record parsing + qname grouping | 15–20s | `parse_bam_record()` extracts tags, cigar, coords |
| Fragment resolution (cgranges) | 30–40s | Per-fragment interval queries |
| Accumulator append + merge | 15–20s | Vector growth + O(N) merge copy |
| Finalize + from_raw copies | 5–10s | Two full data copies |

### Memory Lifecycle (Current)

```
Phase                              C++ Heap    nb::bytes    NumPy     Total Δ
──────────────────────────────────────────────────────────────────────────────
Workers accumulating (scan)        N × ~1.3/N GB                     ~1.3 GB
Merge: main + workers              ~2.6 GB                           ~2.6 GB
Workers freed                      ~1.3 GB                           ~1.3 GB
finalize → nb::bytes               ~1.3 GB    ~1.3 GB               ~2.6 GB
from_raw → np.copy()               ~1.3 GB    ~1.3 GB    ~1.3 GB    ~3.9 GB
GC frees acc + bytes                                      ~1.3 GB    ~1.3 GB
──────────────────────────────────────────────────────────────────────────────
Observed RSS delta: ~11 GB (includes macOS allocator retention, vector
over-allocation up to 2×, and temporary duplication boundaries)
```

---

## Phase 3a: Quick Wins

**Effort:** Low  
**Risk:** Minimal  
**Expected impact:** ~10–20s wall time reduction, ~1–2s startup reduction

### Step 3a.1 — Expose `n_decomp_threads` in BamScanConfig

**Problem:** BGZF decompression is hardcoded to 2 threads in the nanobind binding (`nb::arg("n_decomp_threads") = 2`). On Apple M-series (10 cores), this leaves decompression bandwidth on the table.

**Files to modify:**
- [src/rigel/config.py](../../src/rigel/config.py) — Add `n_decomp_threads: int = 4` to `BamScanConfig`
- [src/rigel/pipeline.py](../../src/rigel/pipeline.py) — Pass `n_decomp_threads=scan.n_decomp_threads` to `scanner.scan()`
- [src/rigel/cli.py](../../src/rigel/cli.py) — Expose `--decomp-threads` CLI arg (optional)

**Implementation:**

```python
# config.py — BamScanConfig
n_decomp_threads: int = 4
"""Number of htslib BGZF decompression threads (default 4).
Set to 0 to disable multi-threaded decompression."""
```

```python
# pipeline.py — scan_and_buffer()
result = scanner.scan(
    bam_path,
    n_workers=n_scan,
    n_decomp_threads=scan.n_decomp_threads,
)
```

**Validation:**
- Run profiler with `n_decomp_threads` = 1, 2, 4, 8 on the benchmark BAM
- Measure scan wall time at each setting
- Expected: 5–15s improvement at 4 threads vs current 2

### Step 3a.2 — Pre-allocate Accumulator Vectors with `reserve()`

**Problem:** `FragmentAccumulator` vectors grow via amortized doubling. For 14.7M fragments with ~1.2× CSR expansion: ~17 vectors × ~24 doublings = significant reallocation overhead and 2× over-allocation at peak.

**Files to modify:**
- [src/rigel/native/resolve_context.h](../../src/rigel/native/resolve_context.h) — Add `reserve()` method to `FragmentAccumulator`
- [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp) — Call `reserve()` on worker accumulators after estimating fragment count

**Implementation:**

```cpp
// resolve_context.h — FragmentAccumulator
void reserve(int64_t n_fragments, int64_t n_candidates) {
    splice_type_.reserve(n_fragments);
    exon_strand_.reserve(n_fragments);
    sj_strand_.reserve(n_fragments);
    num_hits_.reserve(n_fragments);
    merge_criteria_.reserve(n_fragments);
    chimera_type_.reserve(n_fragments);
    frag_id_.reserve(n_fragments);
    read_length_.reserve(n_fragments);
    genomic_footprint_.reserve(n_fragments);
    genomic_start_.reserve(n_fragments);
    nm_.reserve(n_fragments);
    t_indices_.reserve(n_candidates);
    t_offsets_.reserve(n_fragments + 1);
    frag_lengths_.reserve(n_candidates);
    exon_bp_.reserve(n_candidates);
    intron_bp_.reserve(n_candidates);
}
```

**Estimation strategy (in `scan()`):**
Use a two-phase approach:
1. Read first 100K BAM records to estimate fragments/record ratio and candidates/fragment ratio
2. Check total BAM record count from header (or file size heuristic)
3. `reserve()` per-worker accumulator = estimated_total / n_workers

Alternatively, simpler: use BAM file size as crude estimate — typical BAM yields ~20 bytes/record, so `file_size / 20` ≈ record count. With ~0.38 fragments/record (14.7M/38.7M), this gives a reasonable estimate.

```cpp
// bam_scanner.cpp — scan() method, after creating worker states
int64_t est_records = /* file size heuristic or first-batch sampling */;
int64_t est_frags_per_worker = (est_records * 4 / 10) / n_workers;  // 0.4 frag/rec
int64_t est_cands_per_worker = est_frags_per_worker * 12 / 10;      // 1.2 cand/frag
for (auto& ws : worker_states) {
    ws.accumulator.reserve(est_frags_per_worker, est_cands_per_worker);
}
```

**Validation:**
- Compare RSS and wall time before/after on benchmark BAM
- Verify exact output match (golden test)
- Expected: 2–5s wall time improvement, ~30% RSS reduction from accumulator over-allocation

### Step 3a.3 — Lazy-load `region_df` and `nrna_df` in TranscriptIndex

**Problem:** `TranscriptIndex.load()` eagerly reads `regions.feather` and processes nRNA data, adding ~0.5–1s to startup even when not needed for the scan phase.

**Files to modify:**
- [src/rigel/index.py](../../src/rigel/index.py) — Convert region_df, nrna_df, and derived structures to `@cached_property` or similar lazy pattern

**Implementation:**

```python
# index.py — TranscriptIndex
# Instead of loading in __init__:
#   self._region_df = pf.read_feather(regions_path)
# Use:
@cached_property
def region_df(self):
    regions_path = self._index_dir / "regions.feather"
    if not regions_path.exists():
        return None
    return pf.read_feather(regions_path)
```

**Constraint:** `region_df` is only needed during `scan_and_buffer()` for `build_region_index()` and later in calibration. The lazy load defers I/O until first access.

**Validation:**
- Index load time drops by ~0.5–1s
- All tests pass (access patterns unchanged)
- No regression in pipeline output

### Step 3a.4 — Memory-map Feather Files

**Problem:** `pf.read_feather()` copies all data into the process heap. With `memory_map=True`, Arrow can zero-copy map columns directly from disk.

**Files to modify:**
- [src/rigel/index.py](../../src/rigel/index.py) — Pass `memory_map=True` to all `read_feather()` / `read_table()` calls

**Implementation:**

```python
# index.py — TranscriptIndex.load()
t_table = pf.read_table(t_path, memory_map=True)
iv_table = pf.read_table(iv_path, memory_map=True)
sj_table = pf.read_table(sj_path, memory_map=True)
```

**Validation:**
- Measure index load time reduction
- Ensure subsequent operations that materialize columns still work
- Expected: 0.5–1s savings

### Phase 3a Completion Criteria

- [ ] All four changes implemented
- [ ] `pytest tests/ -v` passes (no regressions)
- [ ] Profiler run shows measurable improvement
- [ ] Golden output match confirmed

---

## Phase 3b: Zero-Copy Finalization & Merge Elimination

**Effort:** Medium  
**Risk:** Moderate (changes C++/Python data ownership boundary)  
**Expected impact:** ~8–15s wall time, ~1.3 GB peak RSS reduction

### Step 3b.1 — Zero-Copy Finalization via nanobind Capsules

**Problem:** The current finalization path copies data twice:
1. `to_bytes()`: C++ vector → `nb::bytes` (memcpy of all vector data)
2. `np.frombuffer().copy()`: bytes → independent NumPy array

**Solution:** Transfer ownership of C++ vectors directly to NumPy arrays via nanobind's capsule/ndarray mechanism, eliminating both copies.

**Files to modify:**
- [src/rigel/native/resolve_context.h](../../src/rigel/native/resolve_context.h) — New `finalize_zero_copy()` method returning `nb::dict` of `nb::ndarray`
- [src/rigel/buffer.py](../../src/rigel/buffer.py) — Update `_FinalizedChunk.from_raw()` to skip `.copy()` when arrays are already owned

**Implementation approach:**

Use nanobind's `nb::ndarray` with a capsule-backed deleter that owns the underlying `std::vector`. When the Python ndarray is garbage-collected, the capsule destructor frees the C++ vector.

```cpp
// resolve_context.h — helper template
template <typename T>
nb::ndarray<nb::numpy, T> vec_to_ndarray(std::vector<T>&& vec) {
    size_t n = vec.size();
    auto* data = new std::vector<T>(std::move(vec));
    nb::capsule owner(data, [](void* p) noexcept {
        delete static_cast<std::vector<T>*>(p);
    });
    return nb::ndarray<nb::numpy, T>(
        data->data(), {n}, owner);
}

// FragmentAccumulator::finalize_zero_copy()
nb::dict finalize_zero_copy(const std::vector<int32_t>& t_strand_arr) {
    // Compute ambig_strand (same as current finalize())
    // ...
    
    nb::dict result;
    result["splice_type"]   = vec_to_ndarray(std::move(splice_type_));
    result["exon_strand"]   = vec_to_ndarray(std::move(exon_strand_));
    result["sj_strand"]     = vec_to_ndarray(std::move(sj_strand_));
    result["num_hits"]      = vec_to_ndarray(std::move(num_hits_));
    result["merge_criteria"]= vec_to_ndarray(std::move(merge_criteria_));
    result["chimera_type"]  = vec_to_ndarray(std::move(chimera_type_));
    result["t_indices"]     = vec_to_ndarray(std::move(t_indices_));
    result["t_offsets"]     = vec_to_ndarray(std::move(t_offsets_));
    result["frag_lengths"]  = vec_to_ndarray(std::move(frag_lengths_));
    result["exon_bp"]       = vec_to_ndarray(std::move(exon_bp_));
    result["intron_bp"]     = vec_to_ndarray(std::move(intron_bp_));
    result["ambig_strand"]  = vec_to_ndarray(std::move(ambig_strand));
    result["frag_id"]       = vec_to_ndarray(std::move(frag_id_));
    result["read_length"]   = vec_to_ndarray(std::move(read_length_));
    result["genomic_footprint"] = vec_to_ndarray(std::move(genomic_footprint_));
    result["genomic_start"] = vec_to_ndarray(std::move(genomic_start_));
    result["nm"]            = vec_to_ndarray(std::move(nm_));
    result["size"]          = nb::cast(size_);
    return result;
}
```

```python
# buffer.py — _FinalizedChunk.from_raw()
@classmethod
def from_raw(cls, raw: dict) -> "_FinalizedChunk":
    # If values are already numpy arrays (zero-copy path), use directly
    # If values are bytes (legacy path), use np.frombuffer().copy()
    def _to_array(val, dtype):
        if isinstance(val, np.ndarray):
            # Zero-copy: already a numpy array via capsule
            return val if val.dtype == dtype else val.astype(dtype)
        return np.frombuffer(val, dtype=dtype).copy()
    
    return cls(
        splice_type=_to_array(raw["splice_type"], np.uint8),
        exon_strand=_to_array(raw["exon_strand"], np.uint8),
        # ... etc
    )
```

**Key design decisions:**
- The `std::move` transfers vector ownership to the heap-allocated copy inside the capsule
- The accumulator is invalidated after `finalize_zero_copy()` — it cannot be reused
- The `from_raw()` path must still handle both bytes (old API) and ndarray (new API) for backward compatibility during the transition

**Validation:**
- Golden output test must match exactly
- Verify no double-free or use-after-free via ASAN build: `CMAKE_BUILD_TYPE=Debug CXXFLAGS="-fsanitize=address" pip install --no-build-isolation -e .`
- Measure RSS reduction: should see ~1.3 GB drop at finalize boundary

### Step 3b.2 — Eliminate Worker Accumulator Merge

**Problem:** After scan completes, `merge_accumulator_into()` copies all per-worker accumulator data into the main scanner accumulator. For N workers, this is an O(total_data) copy that doubles peak C++ memory momentarily.

**Current flow (bam_scanner.cpp lines 1635–1643):**
```cpp
for (auto& ws : worker_states) {
    stats_.merge_from(ws.stats);
    merge_accumulator_into(accumulator_, ws.accumulator);  // O(N) copy
    merge_strand_obs(strand_obs_, ws.strand_obs);
    merge_fraglen_obs(fraglen_obs_, ws.fraglen_obs);
    merge_region_acc(region_acc_, ws.region_acc);
}
```

**Solution:** Instead of merging, have each worker finalize its own accumulator independently and return N chunks to Python. Python receives a list of chunks, each of which is an independent `_FinalizedChunk`.

**Files to modify:**
- [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp) — New `finalize_worker_accumulators()` method returning list of dicts
- [src/rigel/native/resolve_context.h](../../src/rigel/native/resolve_context.h) — No changes (finalize per-accumulator already works)
- [src/rigel/pipeline.py](../../src/rigel/pipeline.py) — Loop over returned chunks
- [src/rigel/buffer.py](../../src/rigel/buffer.py) — `inject_chunk()` already handles multiple chunks

**Implementation:**

```cpp
// bam_scanner.cpp — BamScanner
nb::list finalize_worker_accumulators(
    const std::vector<int32_t>& t_strand_arr
) {
    nb::list chunks;
    // Finalize main scanner accumulator (for single-thread case)
    if (accumulator_.get_size() > 0) {
        chunks.append(accumulator_.finalize_zero_copy(t_strand_arr));
    }
    // Finalize each worker accumulator independently
    for (auto& ws : worker_states_) {
        if (ws.accumulator.get_size() > 0) {
            chunks.append(ws.accumulator.finalize_zero_copy(t_strand_arr));
        }
    }
    return chunks;
}
```

**Design changes required:**
- `worker_states_` must survive past `scan()` return (currently local). Move to class member or return alongside result dict.
- The `merge_accumulator_into()` call is removed from the scan loop.
- Stats, strand_obs, fraglen_obs, region_acc still need merging (lightweight, not accumulators).
- The main accumulator is only non-empty in single-thread mode.

```python
# pipeline.py — scan_and_buffer()
# Old:
#   raw = scanner.finalize_accumulator(t_strand_arr)
#   chunk = _FinalizedChunk.from_raw(raw)
#   buffer.inject_chunk(chunk)
#
# New:
chunks = scanner.finalize_worker_accumulators(t_strand_arr)
for raw in chunks:
    chunk = _FinalizedChunk.from_raw(raw)
    buffer.inject_chunk(chunk)
```

**Validation:**
- Golden output test (fragment ordering may change → need order-independent comparison or canonical sort)
- **Critical:** Fragment ordering will differ from single-accumulator path. Verify that downstream EM is order-invariant (it should be — EM operates on equivalence classes).
- RSS profile: peak during merge phase should drop by ~1.3 GB

### Step 3b.3 — Single-Pass `fused_score_buffer()`

**Problem:** The C++ `fused_score_buffer()` makes two passes over all fragment data:
- **Pass 1 (count):** Count total output units and candidates to pre-allocate exact-size arrays
- **Pass 2 (fill):** Score and write to pre-allocated arrays

This doubles the scoring work. For 17.7M buffered units, this is ~18s total (two 9s passes).

**Solution:** Switch to a single-pass design with growing vectors, then transfer to NumPy.

**Files to modify:**
- [src/rigel/native/scoring.cpp](../../src/rigel/native/scoring.cpp) — Rewrite `fused_score_buffer()` to single-pass

**Implementation approach:**

```cpp
// scoring.cpp — fused_score_buffer() single-pass variant
//
// Instead of:
//   Pass 1: count → exact sizes
//   Pass 2: fill pre-allocated
//
// Do:
//   Single pass: append to growing vectors (with reserve heuristic)
//   Transfer ownership to NumPy via capsule

// Pre-allocate with estimate: total_fragments * 1.2 (typical expansion)
int64_t est_units = 0;
for (auto& chunk : chunks) est_units += chunk.size;
int64_t est_cands = est_units * 12 / 10;

std::vector<int32_t> unit_offsets;
std::vector<int32_t> t_indices;
// ... other output vectors
unit_offsets.reserve(est_units);
t_indices.reserve(est_cands);
// ... reserve all output vectors

// Single scoring pass
for (auto& chunk : chunks) {
    for (int i = 0; i < chunk.size; i++) {
        // Score fragment, append to vectors
    }
}

// Transfer ownership via capsules (same pattern as 3b.1)
```

**Trade-offs:**
- ~10% vector over-allocation (vs exact pre-allocation)
- Eliminates redundant pass over data
- Net: 40–50% wall time reduction for router stage

**Validation:**
- Exact output match for scoring results (order and values)
- Router wall time should drop from ~18s to ~9–13s
- Memory may slightly increase due to over-allocation, but capsule transfer avoids the copy

### Phase 3b Completion Criteria

- [ ] Zero-copy finalization implemented and tested
- [ ] Worker merge eliminated, multi-chunk finalization working
- [ ] Single-pass scoring implemented
- [ ] `pytest tests/ -v` passes
- [ ] Golden output match confirmed
- [ ] ASAN clean build with no memory errors
- [ ] Profiler shows measurable improvement

---

## Phase 3c: Streaming Chunk Architecture

**Effort:** Medium–High  
**Risk:** Moderate (well-constrained threading model with clean GIL boundaries)  
**Expected impact:** ~10–20s wall time, ~5–10 GB peak RSS reduction

### Architecture

The streaming design replaces the two-call `scan()` + `finalize_worker_accumulators()` API with a single `scan()` that delivers chunks to Python via callback as they are produced. Three thread roles with clean separation of concerns:

```
┌─────────────────────┐
│    Main Thread       │  Python's thread. Owns the GIL.
│                      │  Drains output queue, converts C++ → Python,
│    (GIL consumer)    │  calls chunk_callback.
└──────────┬──────────┘
           │ output_queue.pop(acc)      ← blocks (GIL released)
           │ finalize_zero_copy(acc)    ← GIL held
           │ chunk_callback(dict)       ← GIL held
           │
┌──────────┴──────────┐
│   Output Queue       │  BoundedQueue<FragmentAccumulator>
│   (MPSC, capacity N) │  Workers push, main thread pops.
└──────────┬──────────┘
           │
     ┌─────┴─────┬─────────┬─── ... ───┐
     │           │         │            │
┌────┴───┐ ┌────┴───┐ ┌───┴────┐ ┌────┴───┐
│Worker 0│ │Worker 1│ │Worker 2│ │Worker N│  Pure C++ threads.
│        │ │        │ │        │ │        │  Pop input queue, resolve
│        │ │        │ │        │ │        │  fragments, push full
│        │ │        │ │        │ │        │  accumulators to output.
└────┬───┘ └────┬───┘ └────┬───┘ └────┬───┘
     │          │          │           │
     └─────┬────┴──────────┴─── ... ───┘
           │
┌──────────┴──────────┐
│   Input Queue        │  BoundedQueue<QnameGroup>
│   (SPMC, capacity 2N)│  Reader pushes, workers pop.
└──────────┬──────────┘
           │
┌──────────┴──────────┐
│   Reader Thread      │  Background C++ thread.
│                      │  BAM I/O → qname grouping → push to input queue.
│   (BAM I/O)          │  After EOF: close input → join workers →
│                      │  merge stats → close output queue.
└─────────────────────┘
```

**Key design principle:** No thread ever needs the GIL except the main thread. Workers and reader are pure C++. The GIL is released while the main thread blocks on `output_queue.pop()` and re-acquired before creating Python objects.

### Thread Lifecycle

```
Main Thread                    Reader Thread                  Workers
═══════════                    ═════════════                  ═══════
scan() called
 ├─ Open BAM, create queues
 ├─ Launch workers             ─────────────────────────────→ start
 ├─ Launch reader              ─→ BAM read loop:
 ├─ Release GIL                   sam_read1 → parse →
 │                                input_queue.push(group) ──→ pop(group)
 │                                                            resolve fragment
 │                                                            accumulator.append()
 │                                                            if size ≥ chunk_size:
 │                                                              output_queue.push(acc)
 │                                                              acc = fresh
 ├─ output_queue.pop(acc)  ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←  ↑
 ├─ Acquire GIL                   ...
 ├─ finalize_zero_copy(acc)       ...
 ├─ chunk_callback(dict)          ...
 ├─ Release GIL                   ...
 │   ... repeat ...               BAM EOF
 │                                input_queue.close() ──────→ drain remaining
 │                                                            push final partial acc
 │                                join workers ←←←←←←←←←←←←  exit
 │                                merge stats/strand/fraglen
 │                                output_queue.close()
 ├─ pop() returns false ←←←←←←←  ↑
 ├─ Acquire GIL
 ├─ Join reader thread ←───────  exit
 ├─ Build result dict
 └─ Return
```

### GIL Safety Proof

| Thread | Touches Python objects? | GIL needed? | How? |
|--------|------------------------|-------------|------|
| Main | Yes (finalize_zero_copy, callback) | Yes, cycled | Release before pop(), acquire before finalize |
| Reader | No (htslib, C++ queues, stats) | Never | Pure C++ |
| Workers | No (resolve, accumulate, push C++ objects) | Never | Pure C++ |

The GIL is released while the main thread blocks on `output_queue.pop()`. This means:
- Python signal handlers can fire (Ctrl-C works)
- Workers and reader run freely (no GIL contention)
- Backpressure: if Python callback is slow, `output_queue` fills → workers block on push → natural throttling

### Deadlock Freedom Proof

The system has two queues forming a pipeline: `input_queue → workers → output_queue → main`.

- **Input queue full, workers blocked on pop:** Cannot happen simultaneously — if workers are blocked on pop, input queue is empty, not full.
- **Output queue full, workers blocked on push:** Main thread is either (a) blocked on pop (queue is full → pop succeeds immediately) or (b) executing callback → will return to pop soon. Workers unblock.
- **Reader blocked on input push, output queue full:** Reader blocks on input push only when input queue is full. Workers are blocked on output push (queue full). Main thread pops output → frees worker → worker pops input → frees reader. Chain unblocks.
- **All blocked:** Impossible — main thread's pop on a full output queue always succeeds immediately.

### Step 3c.1 — `scan()` Becomes Streaming

Replace the current two-call API with a unified `scan()` that accepts a callback and chunk size. One call does everything.

**New `scan()` signature (C++):**

```cpp
nb::dict scan(
    const std::string& bam_path,
    nb::callable chunk_callback,
    const std::vector<int32_t>& t_strand_arr,
    int n_workers = 1,
    int n_decomp_threads = 4,
    int64_t chunk_size = 1000000
)
```

**New Python usage:**

```python
result = scanner.scan(
    bam_path,
    chunk_callback=on_chunk,
    t_strand_arr=index.t_to_strand_arr.tolist(),
    n_workers=8,
    n_decomp_threads=4,
    chunk_size=1_000_000,
)
# result contains stats, strand_observations, frag_length_observations
# No separate finalize call — chunks already delivered via callback.
```

**Files to modify:**
- [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp) — Rewrite `scan()`, add reader thread, remove `finalize_worker_accumulators()` and `finalize_accumulator()`

**Implementation:**

```cpp
nb::dict scan(
    const std::string& bam_path,
    nb::callable chunk_callback,
    const std::vector<int32_t>& t_strand_arr,
    int n_workers,
    int n_decomp_threads,
    int64_t chunk_size)
{
    if (n_workers < 1) n_workers = 1;

    // ── Setup (main thread, GIL held) ─────────────────────
    // Open BAM, build tid mapping, create worker states, reserve accumulators.
    // Create input_queue (capacity: 2 * n_workers).
    // Create output_queue (capacity: n_workers).
    BoundedQueue<QnameGroup> input_queue(n_workers * 2);
    BoundedQueue<FragmentAccumulator> output_queue(n_workers);

    // ── Launch workers ────────────────────────────────────
    for (int i = 0; i < n_workers; i++) {
        workers.emplace_back([&input_queue, &output_queue, ctx, chunk_size, i, ...]() {
            WorkerState& ws = *worker_states_[i];
            QnameGroup group;
            while (input_queue.pop(group)) {
                process_qname_group_threaded(group, *ctx, ws, include_multimap);

                // Emit chunk when accumulator reaches threshold
                if (ws.accumulator.get_size() >= chunk_size) {
                    output_queue.push(std::move(ws.accumulator));
                    ws.accumulator = FragmentAccumulator();
                    ws.accumulator.reserve(chunk_size, chunk_size * 3 / 2);
                }
            }
            // Flush remaining data
            if (ws.accumulator.get_size() > 0) {
                output_queue.push(std::move(ws.accumulator));
            }
        });
    }

    // ── Launch reader thread ──────────────────────────────
    BamScanStats reader_stats;
    std::thread reader_thread([&]() {
        // BAM read loop (identical to current scan() reader loop)
        bam1_t* b = bam_init1();
        QnameGroup current_group;
        std::string current_qname;
        int64_t frag_id = 0;

        while (sam_read1(fp, hdr, b) >= 0) {
            reader_stats.total++;
            // ... filter, group by qname, push to input_queue ...
        }
        // Flush last group
        if (!current_group.records.empty()) {
            current_group.frag_id = frag_id++;
            input_queue.push(std::move(current_group));
        }
        input_queue.close();

        // Wait for all workers to finish
        for (auto& w : workers) w.join();

        // Merge lightweight observations (NOT accumulators — already pushed)
        for (auto& ws_ptr : worker_states_) {
            stats_.merge_from(ws_ptr->stats);
            merge_strand_obs(strand_obs_, ws_ptr->strand_obs);
            merge_fraglen_obs(fraglen_obs_, ws_ptr->fraglen_obs);
            if (ws_ptr->region_acc.enabled())
                merge_region_acc(region_acc_, ws_ptr->region_acc);
        }
        stats_.merge_from(reader_stats);
        worker_states_.clear();

        // Signal main thread: all chunks delivered
        output_queue.close();

        bam_destroy1(b);
        bam_hdr_destroy(hdr);
        hts_close(fp);
    });

    // ── Main thread: drain output queue ───────────────────
    {
        nb::gil_scoped_release release;    // free the GIL
        FragmentAccumulator acc;
        while (output_queue.pop(acc)) {    // blocks until chunk or closed
            nb::gil_scoped_acquire acquire;  // re-acquire for Python work
            nb::dict chunk = acc.finalize_zero_copy(t_strand_arr);
            chunk_callback(chunk);
        }
    }

    reader_thread.join();
    return build_result();
}
```

**Removed methods:**
- `finalize_accumulator()` — folded into streaming scan
- `finalize_worker_accumulators()` — folded into streaming scan

**Removed members:**
- `accumulator_` — vestigial (unused since Phase 3b)

**Key subtlety — reader stats:** The reader thread maintains its own `BamScanStats reader_stats` for record-level counters (total, qc_fail, unmapped, duplicate, secondary, supplementary). These are disjoint from worker stats (which track fragment-level counters). The reader merges both into `stats_` after joining workers.

### Step 3c.2 — Worker Chunk Emission

Workers push full accumulators to the output queue and create fresh ones. This is the inner loop of each worker thread:

```cpp
// Inside worker lambda, after process_qname_group_threaded():
if (ws.accumulator.get_size() >= chunk_size) {
    output_queue.push(std::move(ws.accumulator));
    ws.accumulator = FragmentAccumulator();
    ws.accumulator.reserve(chunk_size, chunk_size * 3 / 2);
}
```

**Why this works:**
- `FragmentAccumulator` is implicitly moveable (all members are `std::vector`).
- After `std::move`, the source accumulator is empty. The `= FragmentAccumulator()` creates a fresh one with the leading `t_offsets_[0] = 0` sentinel.
- `reserve(chunk_size, chunk_size * 3/2)` pre-allocates for the next chunk. This is exact (not a heuristic) since we know the target.
- The `BoundedQueue::push()` blocks if the output queue is full — this is the backpressure mechanism that bounds C++ memory.

**Worker exit:** After the input queue is drained, each worker pushes its remaining partial accumulator (if non-empty) and exits. The reader thread joins all workers, then closes the output queue.

### Step 3c.3 — Pipeline Integration

**Files to modify:**
- [src/rigel/pipeline.py](../../src/rigel/pipeline.py) — Simplify `scan_and_buffer()` to use callback

**New `scan_and_buffer()` (complete rewrite of the scan section):**

```python
def scan_and_buffer(bam_path, index, scan):
    stats = PipelineStats()
    strand_models = StrandModels()
    frag_length_models = FragmentLengthModels(max_size=scan.max_frag_length)
    buffer = FragmentBuffer(
        t_strand_arr=index.t_to_strand_arr,
        chunk_size=scan.chunk_size,
        max_memory_bytes=scan.max_memory_bytes,
        spill_dir=scan.spill_dir,
    )

    # ... resolver setup (unchanged) ...

    scanner = _NativeBamScanner(resolve_ctx, sj_spec, ...)

    def _on_chunk(raw_dict):
        chunk = _FinalizedChunk.from_raw(raw_dict)
        buffer.inject_chunk(chunk)

    n_scan = scan.n_scan_threads or os.cpu_count() or 1
    result = scanner.scan(
        bam_path,
        chunk_callback=_on_chunk,
        t_strand_arr=index.t_to_strand_arr.tolist(),
        n_workers=n_scan,
        n_decomp_threads=scan.n_decomp_threads,
        chunk_size=scan.chunk_size,
    )

    # Replay stats / strand / fraglen observations (unchanged)
    _apply_scan_stats(stats, result["stats"])
    _replay_strand_observations(result["strand_observations"], strand_models)
    _replay_fraglen_observations(result["frag_length_observations"], frag_length_models)

    # No finalize step — chunks already delivered via callback.

    # ... region evidence extraction (unchanged) ...

    return stats, strand_models, frag_length_models, buffer, region_counts, fl_table
```

What changed:
- The `scanner.scan()` call gains `chunk_callback`, `t_strand_arr`, and `chunk_size` parameters
- The `finalize_worker_accumulators()` loop is deleted entirely
- The `accumulator_size` check is deleted
- One call does everything: scan, resolve, train, deliver chunks

### Step 3c.4 — Chunk Size Configuration

The existing `BamScanConfig.chunk_size` field (default 1,000,000) is already the right knob. It is now passed through to the C++ `scan()` as the streaming threshold.

**Tuning guidance:**

| chunk_size | Per-chunk memory | Chunks for 15M frags | Notes |
|------------|-----------------|---------------------|-------|
| 500K | ~44 MB | ~30 | Lower peak, more GIL cycles |
| **1M** | **~87 MB** | **~15** | **Default. Good balance.** |
| 2M | ~174 MB | ~8 | Higher peak, fewer GIL cycles |

The GIL acquisition cost per chunk is negligible (~1 µs) compared to the ~87 MB of data being finalized (~5 ms). Even with 30 chunks (500K), the overhead is < 0.1% of scan time.

### Error Handling

If the Python callback raises an exception:

1. The main thread catches the exception (nanobind propagates Python exceptions as C++ exceptions from `chunk_callback()`).
2. The main thread sets `std::atomic<bool> abort_{false}` on the scanner.
3. The main thread continues draining the output queue **discarding chunks** until it closes.
4. Meanwhile, the reader thread checks `abort_` periodically (every N records) and breaks out of the BAM read loop early, closes the input queue.
5. Workers drain naturally, push remaining data (which is discarded by main).
6. Reader joins workers, closes output queue.
7. Main thread joins reader, re-throws the original Python exception.

The abort path adds one atomic load per BAM record on the reader thread (~0.5 ns, uncontended). No performance impact on the happy path.

### Memory Profile

```
Steady-state during scan (N=8 workers, chunk_size=1M, ~87 bytes/frag):

  C++ worker accumulators:  8 × ≤87 MB  =  ≤ 696 MB  (each capped at chunk_size)
  C++ output queue:         8 × ≤87 MB  =  ≤ 696 MB  (capacity N, worst case full)
  C++ total:                             ≤  1.4 GB

  Python buffer:            bounded by max_memory_bytes (default 2 GB)
  Python spill:             excess chunks written to Arrow IPC on disk

  Total peak:               ≤ 3.4 GB    (vs current ~11 GB observed)
  Expected peak:            ~ 2.5 GB    (mix of full/partial accumulators)
```

**Why this is dramatically better:** The current architecture lets accumulator memory grow unbounded (all fragments in one accumulator per worker, ~1.3 GB total), then doubles it during finalize-to-numpy. The streaming architecture caps each accumulator at `chunk_size` fragments and pipelines delivery to Python, which can spill excess to disk.

### Test Migration

Tests that call `scan()` and `finalize_worker_accumulators()` separately need updating. The migration is mechanical:

```python
# OLD (two calls):
result = scanner.scan(bam_path, n_workers=4, n_decomp_threads=2)
chunks = scanner.finalize_worker_accumulators(t_strand_arr)
for raw in chunks:
    buffer.inject_chunk(_FinalizedChunk.from_raw(raw))

# NEW (one call, callback):
chunks = []
result = scanner.scan(
    bam_path,
    chunk_callback=chunks.append,
    t_strand_arr=t_strand_arr,
    n_workers=4,
    n_decomp_threads=2,
    chunk_size=999_999_999,  # large → one chunk per worker (batch behavior)
)
for raw in chunks:
    buffer.inject_chunk(_FinalizedChunk.from_raw(raw))
```

Setting `chunk_size` to a very large value replicates the old "accumulate everything" behavior — workers never hit the threshold during scanning and only push their final partial accumulator at exit.

### Phase 3c Completion Criteria

- [ ] `scan()` rewritten with reader thread + output queue + callback
- [ ] Workers emit bounded chunks at configurable threshold
- [ ] Main thread drains with proper GIL acquire/release cycling
- [ ] `finalize_worker_accumulators()` and `finalize_accumulator()` removed
- [ ] `pipeline.py` uses the new single-call API
- [ ] `pytest tests/ -v` passes (all tests migrated)
- [ ] Golden output match (order-independent — EM is order-invariant)
- [ ] Peak RSS reduced by ≥ 5 GB on benchmark workload
- [ ] Profiler confirms wall time improvement
- [ ] No deadlocks (verified with chunk_size=1000 stress test)
- [ ] Ctrl-C works during scan (GIL released while blocking)

---

## Phase 3d: Router Parallelization & Pipeline Overlap

**Effort:** Medium–High  
**Risk:** Moderate  
**Expected impact:** ~5–9s wall time reduction  
**Prerequisite:** Phase 3b (single-pass scoring) and Phase 3c (streaming chunks)

### Step 3d.1 — Parallel `fused_score_buffer()` Across Chunks

**Problem:** `fused_score_buffer()` processes all chunks sequentially in a single thread. With streaming chunks, multiple chunks can be scored in parallel.

**Files to modify:**
- [src/rigel/native/scoring.cpp](../../src/rigel/native/scoring.cpp) — Add OpenMP parallel-for over chunks
- [src/rigel/scan.py](../../src/rigel/scan.py) — Pass chunks individually or in batches

**Implementation approach:**

```cpp
// scoring.cpp — parallel chunk scoring
// Each chunk is independently scorable (no shared mutable state)
#pragma omp parallel for schedule(dynamic)
for (int c = 0; c < n_chunks; c++) {
    score_single_chunk(chunks[c], thread_local_output[c]);
}
// Merge thread-local outputs into final CSR
```

**Constraint:** The output CSR must be assembled in deterministic order. Use per-chunk local outputs, then concatenate in chunk order.

### Step 3d.2 — Score-While-Scanning Pipeline (Deferred to Phase 4)

Scoring requires finalized strand and fragment-length models, which are only available after the scan completes. True scan-score overlap would need preliminary models — high complexity, marginal benefit. Defer.

---

## Implementation Sequence & Dependencies

```
Phase 3a.1 (decomp threads)  ──── independent ────┐
Phase 3a.2 (reserve)         ──── independent ────┤
Phase 3a.3 (lazy index)      ──── independent ────┤
Phase 3a.4 (mmap feather)    ──── independent ────┤
                                                   ├─→ Profile & validate 3a
Phase 3b.1 (zero-copy)       ──── depends on nothing extra ─┐
Phase 3b.2 (merge elim)      ──── depends on 3b.1 ──────────┤
Phase 3b.3 (single-pass)     ──── independent ───────────────┤
                                                              ├─→ Profile & validate 3b
Phase 3c.1 (streaming scan)  ──── depends on 3b.1 ─────────────┐
Phase 3c.2 (worker emission) ──── part of 3c.1 ────────────────┤
Phase 3c.3 (pipeline integ)  ──── depends on 3c.1 ─────────────┤
Phase 3c.4 (chunk config)    ──── depends on 3c.1 ─────────────┤
                                                                ├─→ Profile & validate 3c
Phase 3d.1 (parallel score)  ──── depends on 3b.3 ──────────────→ Profile & validate 3d
Phase 3d.2 (pipelining)      ──── defer to Phase 4
```

### Risk Matrix

| Step | Risk | Mitigation |
|------|------|------------|
| 3a.1 | None | Simple config pass-through |
| 3a.2 | Low — estimates may be inaccurate | Over-estimate by 50%; worse case = unused reserved memory |
| 3a.3 | Low — may break callers that assume eager load | Search all `region_df` access sites; test thoroughly |
| 3a.4 | Low — some Arrow operations may copy anyway | Measure; if no improvement, revert |
| 3b.1 | Medium — capsule lifecycle errors | ASAN build, stress test, fuzz with small chunks |
| 3b.2 | Medium — fragment ordering changes | Verify EM is order-invariant; add canonical sort to golden test comparisons |
| 3b.3 | Medium — vector over-allocation waste | Limit to 10% over-reserve; shrink_to_fit before transfer |
| 3c.1 | Medium — reader thread + GIL cycling | Deadlock-free by construction (pipeline topology); stress test with chunk_size=1000 |
| 3c.2 | Low — worker flush is a simple threshold check | Accumulator move semantics are trivially correct (all vector members) |
| 3c.3 | Low — mechanical pipeline.py rewrite | One-call API is simpler than the two-call API it replaces |
| 3d.1 | Medium — determinism of parallel scoring | Per-chunk local buffers, deterministic merge |

### Regression Testing Strategy

**For every step:**

1. Run `pytest tests/ -v` — all existing tests must pass
2. Run `pytest tests/test_golden_output.py -v` — golden output must match
3. For steps that change fragment ordering (3b.2, 3c.*):
   - Add order-independent comparison mode to golden test
   - Or sort output by (gene, transcript) before comparison
4. Profile with `scripts/profiler.py` on the benchmark BAM
5. Verify correctness with `scripts/synthetic_sim_sweep.py` on at least one config

### Profiling Gates

After each phase, run the profiler and record:

| Metric | Phase 3a Target | Phase 3b Target | Phase 3c Target |
|--------|----------------|----------------|----------------|
| scan_and_buffer (s) | ≤ 95s | ≤ 85s | ≤ 75s |
| fragment_router_scan (s) | 18s (unchanged) | ≤ 12s | ≤ 12s |
| Peak RSS (MB) | ≤ 20,500 | ≤ 19,000 | ≤ 14,000 |
| Total wall time (s) | ≤ 240s | ≤ 225s | ≤ 200s |

---

## Appendix A: File Change Matrix

| File | 3a.1 | 3a.2 | 3a.3 | 3a.4 | 3b.1 | 3b.2 | 3b.3 | 3c.1 | 3c.3 | 3d.1 |
|------|------|------|------|------|------|------|------|------|------|------|
| config.py | ✓ | | | | | | | | | |
| pipeline.py | ✓ | | | | | ✓ | | ✓ | ✓ | |
| buffer.py | | | | | ✓ | | | | | |
| index.py | | | ✓ | ✓ | | | | | | |
| scan.py | | | | | | | | | | ✓ |
| resolve_context.h | | ✓ | | | ✓ | | | | | |
| bam_scanner.cpp | | ✓ | | | | ✓ | | ✓ | | |
| scoring.cpp | | | | | | | ✓ | | | ✓ |
| cli.py | ✓ | | | | | | | | | |

## Appendix B: Projected Performance Summary

| Metric | Current | After 3a | After 3b | After 3c | After 3d |
|--------|---------|----------|----------|----------|----------|
| Scan wall time | 109.8s | ~95s | ~85s | ~75s | ~75s |
| Router wall time | 18.1s | 18s | ~12s | ~12s | ~8s |
| Index load time | 8.2s | ~6.5s | ~6.5s | ~6.5s | ~6.5s |
| Peak RSS | 20.8 GB | ~20.5 GB | ~19 GB | ~14 GB | ~14 GB |
| Total pipeline | 254.3s | ~240s | ~225s | ~200s | ~195s |
| Speedup vs baseline (392.8s) | 1.54× | 1.64× | 1.75× | 1.96× | 2.01× |

## Appendix C: Beyond Phase 3

Remaining wall time (~195s) dominated by:
- **Locus EM: ~109s** — Requires algorithmic changes (locus decomposition, transcript pruning) → Phase 4
- **BAM I/O: ~25–35s** — htslib sequential name-sorted read; BGZF decompression bound → limited optimization potential
- **Fragment resolution: ~30–40s** — cgranges queries, batch optimization potential → Phase 4
