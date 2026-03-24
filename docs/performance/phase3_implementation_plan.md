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

**Effort:** High  
**Risk:** High (fundamental architecture change to C++/Python boundary)  
**Expected impact:** ~10–20s wall time, ~5–10 GB peak RSS reduction

### Overview

The streaming architecture changes the fundamental contract between C++ scanner and Python buffer. Instead of accumulating ALL fragments in C++ and returning one giant result, workers periodically emit bounded chunks to Python while scanning continues.

```
Current:  scan(ALL) → merge → finalize → Python
Proposed: scan(chunk₁) → Python(chunk₁)
          scan(chunk₂) → Python(chunk₂)   (concurrent with spill of chunk₁)
          scan(chunk₃) → Python(chunk₃)
          ...
          scan_done → finalize remaining
```

### Step 3c.1 — Worker-Side Chunk Emission

**Problem:** Each worker's `FragmentAccumulator` grows without bound during the entire scan. For 14.7M fragments / N workers, each accumulator holds ~1.3GB / N.

**Solution:** When a worker's accumulator reaches a configurable chunk size (e.g., 1M fragments), it finalizes and emits the chunk, then creates a fresh accumulator.

**Files to modify:**
- [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp) — Add chunk emission logic to worker loop
- [src/rigel/native/resolve_context.h](../../src/rigel/native/resolve_context.h) — No changes
- [src/rigel/native/thread_queue.h](../../src/rigel/native/thread_queue.h) — Add output queue for finalized chunks

**Implementation approach:**

Two architectural options:

#### Option A: Callback-Based (Recommended)

Workers push finalized chunks to a shared output queue. A dedicated Python-side consumer thread (or the main thread after scan) drains the queue.

```cpp
// bam_scanner.cpp — WorkerState additions
struct WorkerState {
    // ... existing fields ...
    int64_t chunk_size_limit;
    BoundedQueue<nb::dict>* output_queue;  // shared output queue
    std::vector<int32_t>* t_strand_arr;    // for finalize
    
    void maybe_emit_chunk() {
        if (accumulator.get_size() >= chunk_size_limit) {
            auto chunk = accumulator.finalize_zero_copy(*t_strand_arr);
            output_queue->push(std::move(chunk));
            accumulator = FragmentAccumulator();
            accumulator.reserve(chunk_size_limit, chunk_size_limit * 12 / 10);
        }
    }
};

// In process_qname_group_threaded(), after accumulator.append():
ws.maybe_emit_chunk();
```

**Issue with Option A:** `nb::dict` and nanobind Python objects cannot be created/manipulated from non-Python threads without holding the GIL. This makes a callback/queue approach with Python objects complex.

#### Option B: C++ Output Queue with Deferred Finalization (Recommended)

Workers push raw C++ accumulators (moved) to a shared output queue. The main thread (holding GIL) drains the queue after scan and finalizes each to Python objects.

```cpp
// bam_scanner.cpp
BoundedQueue<FragmentAccumulator> chunk_queue_;  // capacity = 2 × n_workers

// Worker loop:
void process_in_worker(QnameGroup& group, WorkerState& ws) {
    // ... resolve + append ...
    if (ws.accumulator.get_size() >= ws.chunk_size_limit) {
        chunk_queue_.push(std::move(ws.accumulator));
        ws.accumulator = FragmentAccumulator();
        ws.accumulator.reserve(ws.chunk_size_limit, ...);
    }
}

// After scan:
nb::list finalize_all_chunks(const std::vector<int32_t>& t_strand_arr) {
    nb::list result;
    // Drain output queue
    FragmentAccumulator acc;
    while (chunk_queue_.pop(acc)) {
        result.append(acc.finalize_zero_copy(t_strand_arr));
    }
    // Finalize remaining per-worker accumulators
    for (auto& ws : worker_states_) {
        if (ws.accumulator.get_size() > 0) {
            result.append(ws.accumulator.finalize_zero_copy(t_strand_arr));
        }
    }
    return result;
}
```

**Memory bound:** At any point, there are at most `n_workers` in-flight accumulators (each ≤ chunk_size_limit) plus `queue_capacity` drained accumulators. With chunk_size_limit = 1M fragments and 87 bytes/fragment:

```
Peak accumulator memory = (n_workers + queue_capacity) × 1M × 87 bytes
                        = (8 + 16) × 87 MB
                        = ~2.1 GB  (vs current ~1.3 GB + merge copy overhead)
```

For the actual data, peak is bounded. The drained accumulators are finalized immediately to Python arrays and freed.

#### Option C: True Streaming with Python Callback (Most Complex, Highest Reward)

A Python callback is invoked from C++ (on the main thread) every time the output queue has a chunk. This allows Python to spill chunks to disk while scanning continues, overlapping I/O.

```python
# pipeline.py
def _on_chunk(raw_dict):
    chunk = _FinalizedChunk.from_raw(raw_dict)
    buffer.inject_chunk(chunk)  # may spill to disk

result = scanner.scan_streaming(
    bam_path,
    n_workers=n_scan,
    n_decomp_threads=scan.n_decomp_threads,
    chunk_callback=_on_chunk,
    chunk_size=scan.chunk_size,
)
```

```cpp
// bam_scanner.cpp — scan_streaming()
// Reader thread also acts as chunk drainer:
// After pushing each qname_group and periodically:
while (chunk_queue_.try_pop(acc)) {
    nb::dict chunk = acc.finalize_zero_copy(t_strand_arr);
    callback(chunk);  // Call Python callback (holding GIL)
}
```

**GIL management:** The reader thread holds the GIL when calling the Python callback. Workers release the GIL during resolution work. The reader must briefly acquire the GIL for callbacks. This is feasible because the reader thread alternates between BAM I/O (no GIL needed) and queue pushing.

### Step 3c.2 — Python Buffer Integration with Streaming

**Files to modify:**
- [src/rigel/pipeline.py](../../src/rigel/pipeline.py) — New `_scan_streaming()` path
- [src/rigel/buffer.py](../../src/rigel/buffer.py) — Chunk injection already works; add incremental memory tracking

**Implementation:**

```python
# pipeline.py — scan_and_buffer() modified for streaming
def scan_and_buffer(bam_path, index, scan_config):
    buffer = FragmentBuffer(...)
    
    def _on_chunk(raw_dict):
        chunk = _FinalizedChunk.from_raw(raw_dict)
        buffer.inject_chunk(chunk)
        logger.debug(f"Chunk: {chunk.size:,} fragments, "
                     f"buffer: {buffer.memory_bytes / 1024**2:.0f} MB")
    
    result = scanner.scan_streaming(
        bam_path,
        n_workers=n_scan,
        n_decomp_threads=scan_config.n_decomp_threads,
        chunk_callback=_on_chunk,
        chunk_size=scan_config.chunk_size,
    )
    
    # Stats/observations replay (unchanged)
    _apply_scan_stats(stats, result["stats"])
    # ...
```

### Step 3c.3 — Configurable Chunk Size

**Files to modify:**
- [src/rigel/config.py](../../src/rigel/config.py) — Already has `chunk_size: int = 1_000_000`
- Wire through to C++ `scan_streaming()` as `chunk_size` parameter

**Guidance for tuning:**
- **1M fragments**: ~87 MB per chunk, good balance of memory and overhead
- **500K fragments**: ~44 MB per chunk, lower peak memory
- **2M fragments**: ~174 MB per chunk, less overhead but higher peak

### Step 3c.4 — Backward Compatibility

Maintain the existing `scan()` + `finalize_accumulator()` API alongside the new `scan_streaming()` API. The old path is useful for:
- Small BAMs where streaming overhead isn't justified
- Testing and debugging
- Gradual migration

**Selection logic:**
```python
# pipeline.py
if scan_config.streaming:  # new config flag, default True
    # streaming path
else:
    # legacy path
```

### Phase 3c Memory Profile (Expected)

```
Phase                              C++ Heap    Python      Total Δ
──────────────────────────────────────────────────────────────────
Workers accumulating               N × ~87 MB              ~700 MB (8 workers)
Chunk emitted → queue              +87 MB                  ~787 MB
Python callback → buffer           -87 MB C++  +87 MB Py   ~787 MB
Buffer spill (if threshold)                    -87 MB      ~700 MB
... repeats for each chunk ...
──────────────────────────────────────────────────────────────────
Peak C++: ~700 MB  (vs current ~2.6 GB during merge)
Peak Python buffer: bounded by max_memory_bytes config (default 2 GB)
Total peak: ~2.7 GB  (vs current ~3.9 GB, observed ~11+ GB with overhead)
```

### Phase 3c Completion Criteria

- [ ] Streaming scan implemented (Option B or C)
- [ ] Chunk emission from workers at configurable threshold
- [ ] Python buffer receives chunks incrementally
- [ ] Peak RSS reduced by ≥ 5 GB on benchmark workload
- [ ] `pytest tests/ -v` passes
- [ ] Golden output match (order-independent)
- [ ] Profiler confirms wall time improvement
- [ ] No deadlocks under heavy workloads (stress test with small chunk sizes)

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

### Step 3d.2 — Score-While-Scanning Pipeline

**Problem:** Currently, scoring waits for the entire scan to finish. With streaming chunks, scoring can start as chunks arrive.

**Design:**

This is an extension of Step 3c.3 (streaming callback). The callback not only buffers the chunk but also queues it for scoring:

```python
# pipeline.py — conceptual
def _on_chunk(raw_dict):
    chunk = _FinalizedChunk.from_raw(raw_dict)
    buffer.inject_chunk(chunk)
    scoring_queue.put(chunk)  # score in parallel
```

However, scoring requires the strand and fragment-length models, which are only finalized after the scan completes. So true overlap is only possible for chunks scored with a *preliminary* model, then re-scored with the final model only for affected chunks.

**Decision:** Defer this to Phase 4. The complexity of preliminary model scoring is high and the benefit is marginal given that the scan dominates wall time.

### Phase 3d Completion Criteria

- [ ] Parallel chunk scoring implemented
- [ ] Deterministic output ordering maintained
- [ ] Router wall time drops from ~18s to ~9–13s with single-pass, or ~5–8s with parallelism

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
Phase 3c.1 (streaming)       ──── depends on 3b.1 + 3b.2 ──┐
Phase 3c.2 (buffer integ)    ──── depends on 3c.1 ──────────┤
Phase 3c.3 (chunk config)    ──── depends on 3c.1 ──────────┤
Phase 3c.4 (compat)          ──── depends on 3c.1 ──────────┤
                                                              ├─→ Profile & validate 3c
Phase 3d.1 (parallel score)  ──── depends on 3b.3 ────────────→ Profile & validate 3d
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
| 3c.1 | High — threading + GIL + queue design | Start with Option B (deferred finalization); Option C later |
| 3c.2 | Medium — incremental buffer bookkeeping | Existing inject_chunk handles this |
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

| File | 3a.1 | 3a.2 | 3a.3 | 3a.4 | 3b.1 | 3b.2 | 3b.3 | 3c.1 | 3c.2 | 3d.1 |
|------|------|------|------|------|------|------|------|------|------|------|
| config.py | ✓ | | | | | | | | | |
| pipeline.py | ✓ | | | | | ✓ | | ✓ | ✓ | |
| buffer.py | | | | | ✓ | | | | ✓ | |
| index.py | | | ✓ | ✓ | | | | | | |
| scan.py | | | | | | | | | | ✓ |
| resolve_context.h | | ✓ | | | ✓ | | | | | |
| bam_scanner.cpp | | ✓ | | | | ✓ | | ✓ | | |
| thread_queue.h | | | | | | | | ✓ | | |
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
