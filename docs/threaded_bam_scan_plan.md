# Threaded BAM Scan — Implementation Plan

## Executive Summary

Parallelize the BAM scan stage (`scan_and_buffer`) using a producer-consumer
architecture. The reader thread handles serial BAM I/O + qname grouping
while a pool of worker threads performs all compute-heavy parsing, overlap
resolution, and buffer accumulation. Combined with htslib's built-in BGZF
decompression threading, this targets a **2–3× speedup** of the scan stage
(91.5 s → 30–45 s).

## Architecture

```
┌──────────────────────────┐
│  htslib BGZF thread pool │   (hts_set_threads, transparent)
│  Decompresses blocks     │
│  ahead of reader thread  │
└──────────┬───────────────┘
           │ pre-decompressed blocks
           ▼
┌──────────────────────────┐
│  Reader Thread           │
│  sam_read1() loop        │
│  Flag filter             │
│  bam_dup1() deep copy    │
│  Group by qname          │
│  Assign frag_id          │
│  Enqueue to BoundedQueue │
└──────────┬───────────────┘
           │ QnameGroup (vector<bam1_t*>, frag_id)
           ▼
┌──────────────────────────┐
│  Bounded SPMC Queue      │
│  capacity ~2048 groups   │
│  ~2 MB in-flight data    │
└──────────┬───────────────┘
           │ dequeue
     ┌─────┼─────┬─────┐
     ▼     ▼     ▼     ▼
┌────────┐┌────────┐┌────────┐┌────────┐
│Worker 0││Worker 1││Worker 2││Worker 3│
│        ││        ││        ││        │
│ Parse  ││ Parse  ││ Parse  ││ Parse  │
│ tags,  ││ tags,  ││ tags,  ││ tags,  │
│ CIGAR  ││ CIGAR  ││ CIGAR  ││ CIGAR  │
│        ││        ││        ││        │
│ Group  ││ Group  ││ Group  ││ Group  │
│ by hit ││ by hit ││ by hit ││ by hit │
│        ││        ││        ││        │
│ Build  ││ Build  ││ Build  ││ Build  │
│ frag   ││ frag   ││ frag   ││ frag   │
│        ││        ││        ││        │
│ Resolve││ Resolve││ Resolve││ Resolve│
│ core() ││ core() ││ core() ││ core() │
│        ││        ││        ││        │
│ Accum  ││ Accum  ││ Accum  ││ Accum  │
│ ulate  ││ ulate  ││ ulate  ││ ulate  │
└────────┘└────────┘└────────┘└────────┘
  Each worker has:
  - Own ResolverScratch (cgranges query bufs + bp counters)
  - Own FragmentAccumulator
  - Own BamScanStats
  - Own StrandObservations
  - Own FragLenObservations
```

### After all workers join:
1. Merge per-worker stats (sum counters)
2. Concatenate per-worker observation vectors
3. Concatenate per-worker accumulators (append all vectors)
4. Call `finalize()` on merged accumulator → bytes for Python

## Design Decisions

### Decision 1: Producer does NO BAM parsing

**The producer copies raw `bam1_t` records and groups by qname. Workers
do ALL parsing (CIGAR, tags, resolution).**

Rationale — cost comparison for 24M records:

| Operation | Per-record cost | Total (24M) |
|-----------|----------------|-------------|
| `bam_dup1()` (deep copy ~400 bytes) | ~20 ns | ~0.5 s |
| Tag extraction (NM, NH, HI, XS) | ~150 ns | ~3.6 s |
| CIGAR parsing → exon blocks | ~100 ns | ~2.4 s |
| `read_sj_strand()` with fallback | ~50 ns | ~1.2 s |

The deep copy is **~15× cheaper** than parsing. Every nanosecond saved in
the producer is a nanosecond off the serial critical path. Parsing in the
producer would add ~7 s to the bottleneck; the copy adds ~0.5 s. The
trade-off is clear: copy the raw bytes, push parsing to N workers where
it divides by N.

The `bam1_t` record is self-contained: `core` struct (fixed fields) +
`data` pointer (qname + CIGAR + seq + qual + aux, `l_data` bytes). The
`bam_dup1()` function allocates a new `bam1_t` and copies everything.

Memory in flight: with a 2048-slot queue and ~2 records/group × 400
bytes/record, the queue holds ~1.6 MB. Negligible.

### Decision 2: htslib BGZF decompression threads

Call `hts_set_threads(fp, n_decomp)` (where `n_decomp` = 2 by default)
to let htslib decompress BGZF blocks ahead of the reader in a background
thread pool. This is built-in, well-tested, and transparently turns
`sam_read1()` from a decompression-bound call into a fast buffer read.

Without `hts_set_threads`: reader pays ~1.0 µs/record (decompress + parse).
With `hts_set_threads(2)`: reader pays ~0.3–0.5 µs/record (read from
pre-decompressed buffer + minimal bookkeeping).

### Decision 3: What the reader extracts

The reader needs just enough information to:
1. **Filter** records (flag check for QC fail, unmapped, duplicate)
2. **Detect qname boundaries** (string compare on `bam_get_qname`)
3. **Assign `frag_id`** (simple counter increment)

All of this operates on the `bam1_t` core fields and qname, which are
immediately accessible without parsing CIGAR or aux data. The reader does
NOT call `bam_aux_get`, `parse_cigar`, or any tag/CIGAR functions.

The only operation that requires accessing `bam1_t::data` is `bam_get_qname`,
which is a pointer offset (zero-copy). The `bam_dup1` deep copy preserves
the data blob for workers.

### Decision 4: Thread allocation

Default thread counts (configurable via `n_scan_threads` parameter):

| Role | Threads | Notes |
|------|---------|-------|
| Reader | 1 | Serial BAM I/O |
| BGZF decompressors | 2 | Via `hts_set_threads()` |
| Workers | N (default 4) | Parse + resolve + accumulate |
| **Total** | **3 + N = 7** | |

The BGZF threads are managed internally by htslib and don't need our
coordination. The reader and workers communicate via the bounded queue.

## Expected Speedup

### Time budget (current serial, 24M records, 91.5 s)

| Component | Est. time | Thread |
|-----------|-----------|--------|
| BGZF decompression | ~15 s | htslib pool (hidden) |
| BAM record parsing from buffer | ~8 s | Reader |
| Flag filter + qname compare | ~1 s | Reader |
| `bam_dup1` copies | ~0.5 s | Reader |
| **Reader total** | **~10 s** | |
| Tag extraction (NM, NH, HI, XS/ts) | ~7 s | Workers |
| CIGAR parse → exon blocks | ~2.5 s | Workers |
| `group_records_by_hit` | ~2.5 s | Workers |
| `build_fragment` (merge exons+introns) | ~7 s | Workers |
| `_resolve_core` (cgranges overlap) | ~29 s | Workers |
| `compute_frag_lengths` (SJ gap lookup) | ~12 s | Workers |
| Model training obs + accum write | ~7 s | Workers |
| **Worker total** | **~67 s** | **/N** |
| Queue overhead, merge, cache effects | ~3 s | |

### Projected wall time

$$T = \max(T_{reader},\; T_{worker}/N) + T_{overhead}$$

| N workers | Reader | Workers/N | Overhead | **Total scan** | **vs. 91.5 s** |
|-----------|--------|-----------|----------|----------------|-----------------|
| 2 | 10 s | 33.5 s | 3 s | **~37 s** | **2.5×** |
| 4 | 10 s | 16.8 s | 4 s | **~21 s** | **4.4×** |
| 6 | 10 s | 11.2 s | 5 s | **~15 s** | **6.1×** |
| 8 | 10 s | 8.4 s | 5 s | **~15 s** | **6.1×** (reader-bound) |

At N=4, total pipeline drops from 131.5 s → **~62 s** (scan 21 + EM 31 + other 10).
At N=6, diminishing returns — reader becomes the bottleneck at ~10 s.

## Implementation Plan

### Phase 1: Refactor `FragmentResolver` scratch state (prerequisite)

**Goal**: Separate read-only index data from per-call mutable scratch
buffers, enabling multiple threads to call `_resolve_core()` concurrently.

**File**: `resolve_context.h`

#### 1a. Create `ResolverScratch` struct

```cpp
struct ResolverScratch {
    // Per-transcript bp counters (sparse-cleaned)
    std::vector<int32_t> t_exon_bp;
    std::vector<int32_t> t_transcript_bp;
    std::vector<int32_t> t_unambig_intron_bp;
    std::vector<uint8_t> t_dirty;
    std::vector<int32_t> dirty_indices;

    // cgranges query buffers (reusable per-call)
    int64_t* buf = nullptr;
    int64_t  buf_cap = 0;
    int64_t* sj_buf = nullptr;
    int64_t  sj_buf_cap = 0;

    explicit ResolverScratch(int32_t n_transcripts)
        : t_exon_bp(n_transcripts, 0),
          t_transcript_bp(n_transcripts, 0),
          t_unambig_intron_bp(n_transcripts, 0),
          t_dirty(n_transcripts, 0)
    {
        dirty_indices.reserve(512);
    }

    ~ResolverScratch() {
        free(buf);
        free(sj_buf);
    }

    // Non-copyable, movable
    ResolverScratch(const ResolverScratch&) = delete;
    ResolverScratch& operator=(const ResolverScratch&) = delete;
    ResolverScratch(ResolverScratch&& o) noexcept;
    ResolverScratch& operator=(ResolverScratch&& o) noexcept;

    void mark_dirty(int32_t t_idx) {
        if (!t_dirty[t_idx]) {
            t_dirty[t_idx] = 1;
            dirty_indices.push_back(t_idx);
        }
    }

    void clean() {
        for (int32_t t : dirty_indices) {
            t_exon_bp[t] = 0;
            t_transcript_bp[t] = 0;
            t_unambig_intron_bp[t] = 0;
            t_dirty[t] = 0;
        }
        dirty_indices.clear();
    }
};
```

#### 1b. Add overloads to `FragmentResolver` that accept external scratch

Add new method signatures that take `ResolverScratch&`:

```cpp
bool _resolve_core(
    const std::vector<ExonBlock>& exons,
    const std::vector<IntronBlock>& introns,
    int32_t genomic_footprint,
    RawResolveResult& cr,
    ResolverScratch& scratch);   // <-- NEW parameter

std::unordered_map<int32_t, int32_t> compute_frag_lengths(
    const std::vector<ExonBlock>& exons,
    const std::vector<IntronBlock>& introns,
    const std::vector<int32_t>& t_inds,
    ResolverScratch& scratch) const;   // <-- NEW parameter
```

The old signatures (no scratch param) remain as wrappers that forward to
the internal scratch — preserving backward compatibility for
`BamAnnotationWriter` and Python-exposed `resolve()` / `resolve_fragment()`.

```cpp
// Backward-compatible wrapper (uses internal scratch)
bool _resolve_core(
    const std::vector<ExonBlock>& exons,
    const std::vector<IntronBlock>& introns,
    int32_t genomic_footprint,
    RawResolveResult& cr)
{
    return _resolve_core(exons, introns, genomic_footprint, cr, scratch_);
}
```

The existing member scratch buffers (`t_exon_bp_`, etc.) become a private
`ResolverScratch scratch_` member — zero behavioral change for existing
callers.

**Memory per scratch**: 254K transcripts × (4+4+4+1) bytes = ~3.3 MB.
For 4 workers = ~13 MB total. Negligible.

### Phase 2: Threading infrastructure

**File**: New header `src/rigel/native/thread_queue.h`

#### 2a. Bounded SPMC queue

```cpp
template <typename T>
class BoundedQueue {
    std::deque<T> items_;
    std::mutex mutex_;
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    size_t capacity_;
    bool done_ = false;  // producer signals completion

public:
    explicit BoundedQueue(size_t capacity);

    // Producer calls: blocks if queue is full
    void push(T item);

    // Producer calls: signals no more items
    void close();

    // Consumer calls: returns false when queue closed AND empty
    bool pop(T& item);
};
```

#### 2b. QnameGroup work unit

```cpp
struct QnameGroup {
    std::vector<bam1_t*> records;  // deep-copied via bam_dup1
    int64_t frag_id;               // assigned by reader

    ~QnameGroup() {
        for (auto* b : records) {
            if (b) bam_destroy1(b);
        }
    }

    // Move-only
    QnameGroup(QnameGroup&& o) noexcept
        : records(std::move(o.records)), frag_id(o.frag_id)
    {
        o.records.clear();
    }
    QnameGroup& operator=(QnameGroup&& o) noexcept;
};
```

#### 2c. WorkerState — per-thread mutable state

```cpp
struct WorkerState {
    ResolverScratch scratch;
    FragmentAccumulator accumulator;
    BamScanStats stats;
    StrandObservations strand_obs;
    FragLenObservations fraglen_obs;

    explicit WorkerState(int32_t n_transcripts)
        : scratch(n_transcripts) {}
};
```

### Phase 3: Refactor `BamScanner::scan()` into threaded version

**File**: `bam_scanner.cpp`

#### 3a. Add `scan_threaded()` method alongside existing `scan()`

The existing `scan()` method is preserved for single-threaded operation
and backward compatibility. The new method:

```cpp
nb::dict scan_threaded(const std::string& bam_path,
                       int n_workers = 4,
                       int n_decomp_threads = 2);
```

#### 3b. Reader thread function

```cpp
void reader_loop(
    const std::string& bam_path,
    BoundedQueue<QnameGroup>& queue,
    const std::vector<int32_t>& tid_to_ref_id,
    SJTagMode sj_tag_mode,
    bool skip_duplicates,
    bool include_multimap,
    BamScanStats& reader_stats,  // minimal: total, qc_fail, unmapped, dup
    int n_decomp_threads)
{
    htsFile* fp = hts_open(bam_path.c_str(), "rb");
    hts_set_threads(fp, n_decomp_threads);
    bam_hdr_t* hdr = sam_hdr_read(fp);
    bam1_t* b = bam_init1();

    QnameGroup current_group;
    std::string current_qname;
    int64_t frag_id = 0;

    while (sam_read1(fp, hdr, b) >= 0) {
        reader_stats.total++;

        uint16_t flag = b->core.flag;
        if (flag & BAM_FQCFAIL) { reader_stats.qc_fail++; continue; }
        if (flag & BAM_FUNMAP)  { reader_stats.unmapped++;  continue; }
        if (flag & BAM_FDUP) {
            reader_stats.duplicate++;
            if (skip_duplicates) continue;
        }

        const char* qname = bam_get_qname(b);

        if (!current_group.records.empty() && current_qname != qname) {
            current_group.frag_id = frag_id++;
            queue.push(std::move(current_group));
            current_group = QnameGroup{};  // reset
        }

        current_qname = qname;
        current_group.records.push_back(bam_dup1(b));
    }

    // Flush last group
    if (!current_group.records.empty()) {
        current_group.frag_id = frag_id++;
        queue.push(std::move(current_group));
    }

    queue.close();

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(fp);
}
```

**Key points:**
- `bam_dup1(b)` deep-copies each record (~20 ns, ~400 bytes)
- Reader does NO tag extraction, NO CIGAR parsing
- `hts_set_threads` handles BGZF decompression in background
- `frag_id` is assigned sequentially by the reader
- `QnameGroup` owns its `bam1_t*` pointers — destroyed
  after workers process them

#### 3c. Worker thread function

```cpp
void worker_loop(
    BoundedQueue<QnameGroup>& queue,
    FragmentResolver& resolver,      // read-only shared index
    WorkerState& state,              // per-thread mutable state
    const std::vector<int32_t>& tid_to_ref_id,
    SJTagMode sj_tag_mode,
    bool include_multimap)
{
    QnameGroup group;
    while (queue.pop(group)) {
        process_qname_group_threaded(
            group, resolver, state, tid_to_ref_id,
            sj_tag_mode, include_multimap);
        // group destructor frees bam1_t* records
    }
}
```

#### 3d. `process_qname_group_threaded` — worker-side processing

This is essentially the existing `process_qname_group` but:

1. **Accepts raw `bam1_t*`** instead of pre-parsed `ParsedAlignment`
2. **Parses tags and CIGAR** from `bam1_t*` data (moved from reader)
3. **Uses `ResolverScratch&`** for thread-safe `_resolve_core` calls
4. **Writes to `WorkerState`** (accumulator, stats, observations)

```cpp
void process_qname_group_threaded(
    QnameGroup& group,
    FragmentResolver& resolver,
    WorkerState& state,
    const std::vector<int32_t>& tid_to_ref_id,
    SJTagMode sj_tag_mode,
    bool include_multimap)
{
    if (group.records.empty()) return;

    // 1. Parse bam1_t* records into ParsedAlignments
    std::vector<ParsedAlignment> parsed;
    parsed.reserve(group.records.size());
    for (bam1_t* b : group.records) {
        // Same parsing logic as current reader loop:
        // extract flag, pos, tags (NM, NH, HI, XS/ts), parse CIGAR
        ParsedAlignment rec;
        rec.flag = b->core.flag;
        rec.ref_id = /* mapped via tid_to_ref_id */;
        rec.ref_start = b->core.pos;
        rec.mate_ref_id = /* mapped */;
        rec.mate_ref_start = b->core.mpos;
        rec.nm = /* bam_aux_get(b, "NM") */;
        rec.nh = /* bam_aux_get(b, "NH") */;
        rec.hi = /* bam_aux_get(b, "HI") */;
        rec.sj_strand = read_sj_strand(b, sj_tag_mode);
        parse_cigar(b, rec.ref_id, rec.sj_strand, rec.exons, rec.sjs);
        parsed.push_back(std::move(rec));
    }

    // 2. Same logic as existing process_qname_group():
    //    group_records_by_hit, build_fragment, resolve, train, accumulate
    //    But using state.scratch for _resolve_core
    //    and writing to state.accumulator, state.stats, etc.
    // ...
}
```

This function is **identical in logic** to the existing `process_qname_group`
after the parsing step. The only differences are:
- Parsing from `bam1_t*` instead of reading from the current htslib buffer
- Using `state.scratch` for `_resolve_core` and `compute_frag_lengths`
- Writing to `state.stats` / `state.strand_obs` / `state.fraglen_obs` /
  `state.accumulator` instead of `this->stats_` etc.

#### 3e. Merge after join

```cpp
// After all workers join:
BamScanStats merged_stats = reader_stats;  // start with reader's counters
FragmentAccumulator merged_accum;
StrandObservations merged_strand;
FragLenObservations merged_fraglen;

for (auto& ws : worker_states) {
    // Sum stats
    merged_stats.n_read_names     += ws.stats.n_read_names;
    merged_stats.unique           += ws.stats.unique;
    merged_stats.multimapping     += ws.stats.multimapping;
    merged_stats.n_fragments      += ws.stats.n_fragments;
    // ... all integer counters ...

    // Concatenate accumulator vectors (append)
    merged_accum.splice_type_.insert(
        merged_accum.splice_type_.end(),
        ws.accumulator.splice_type_.begin(),
        ws.accumulator.splice_type_.end());
    // ... same for all accumulator vectors ...
    // Fix t_offsets: shift by current merged total
    // ...

    // Concatenate observation vectors
    merged_strand.exonic_spliced_obs.insert(...);
    merged_fraglen.lengths.insert(...);
    // ... etc ...
}
```

**Accumulator merge detail**: The tricky part is `t_offsets_`. Each worker's
`t_offsets_` starts at 0 and indexes into that worker's `t_indices_`. When
merging, the second worker's offsets must be shifted by the first worker's
`t_indices_.size()`. Implementation:

```cpp
void merge_accumulator_into(FragmentAccumulator& dst,
                            FragmentAccumulator& src)
{
    if (src.size_ == 0) return;

    int64_t t_base = static_cast<int64_t>(dst.t_indices_.size());

    // Append per-fragment columns
    dst.splice_type_.insert(dst.splice_type_.end(),
                            src.splice_type_.begin(),
                            src.splice_type_.end());
    // ... repeat for all per-fragment vectors ...

    // Append CSR data (t_indices, frag_lengths, exon_bp, etc.)
    dst.t_indices_.insert(dst.t_indices_.end(),
                          src.t_indices_.begin(),
                          src.t_indices_.end());
    dst.frag_lengths_.insert(dst.frag_lengths_.end(),
                             src.frag_lengths_.begin(),
                             src.frag_lengths_.end());
    // ... etc for exon_bp, intron_bp, unambig_intron_bp ...

    // Shift and append t_offsets (skip src's leading 0)
    for (int32_t i = 1; i <= src.size_; i++) {
        dst.t_offsets_.push_back(src.t_offsets_[i] + t_base);
    }

    dst.size_ += src.size_;
}
```

### Phase 4: Python-side integration

**File**: `pipeline.py`

#### 4a. Add `n_scan_threads` to `BamScanConfig`

```python
@dataclass
class BamScanConfig:
    # ... existing fields ...
    n_scan_threads: int = 4  # Worker threads for BAM scan (0 = auto)
```

#### 4b. Update `scan_and_buffer()` to pass thread count

```python
# In scan_and_buffer():
if scan.n_scan_threads > 1:
    result = scanner.scan_threaded(
        bam_path,
        n_workers=scan.n_scan_threads,
        n_decomp_threads=2,
    )
else:
    result = scanner.scan(bam_path)  # existing single-threaded path
```

#### 4c. CLI exposure

Add `--scan-threads` to the CLI argument parser (default 4).

### Phase 5: Testing

1. **Golden test parity**: All 21 golden tests must pass bit-exact.
   Fragment ordering in the buffer may change (workers process in
   non-deterministic order), but **the EM and final quantification
   are order-independent**. The golden tests compare final abundances,
   not intermediate buffer contents.

2. **Correctness invariant**: The merged stats, strand observations,
   and fragment-length observations must be identical to the
   single-threaded path (up to observation ordering, which doesn't
   affect model fitting).

3. **Stress test**: Run with 1, 2, 4, 8 workers and verify identical
   final results.

4. **Edge cases**: Empty BAM, BAM with only intergenic reads, BAM
   with multimappers, BAM with chimeric reads.

## Thread-Safety Audit

### Read-only shared state (safe for concurrent access):

| Data | Location | Accessed by |
|------|----------|-------------|
| `cgranges_t* cr_` | `FragmentResolver` | `_resolve_core` (overlap query) |
| `cgranges_t* sj_cr_` | `FragmentResolver` | `compute_frag_lengths` |
| `SJMap sj_map_` | `FragmentResolver` | `sj_lookup` |
| `t_set_data_`, `t_set_offsets_` | `FragmentResolver` | `_resolve_core` |
| `sj_map_data_` | `FragmentResolver` | `sj_lookup` |
| `sj_t_index_`, `sj_strand_` | `FragmentResolver` | `compute_frag_lengths` |
| `t_to_g_arr_`, `t_strand_arr_` | `FragmentResolver` | `_resolve_core` |
| `iv_type_` | `FragmentResolver` | `_resolve_core` |
| `ref_to_id_`, `id_to_ref_` | `FragmentResolver` | `_resolve_core` |
| `tid_to_ref_id_` | `BamScanner` | reader + workers |

**Note on cgranges**: `cr_overlap()` writes ONLY to the caller-supplied
`buf` pointer. The tree structure itself (`cr_->root`, `cr_->a`, etc.)
is read-only after `cr_index()`. Each worker thread supplies its own
`scratch.buf` / `scratch.sj_buf`, so there is no data race.

### Per-thread mutable state (isolated, no sharing):

| Data | Where | Thread |
|------|-------|--------|
| `ResolverScratch` (bp counters, query bufs) | `WorkerState` | Each worker |
| `FragmentAccumulator` | `WorkerState` | Each worker |
| `BamScanStats` | `WorkerState` + reader | Each has own |
| `StrandObservations` | `WorkerState` | Each worker |
| `FragLenObservations` | `WorkerState` | Each worker |
| `bam1_t* b` (htslib read buffer) | Reader loop | Reader only |
| `htsFile* fp` | Reader loop | Reader only |

### Synchronization points:

1. **BoundedQueue**: mutex + 2 condition variables (push/pop)
2. **Thread join**: `std::thread::join()` before merge
3. **No other locks needed**

## Files Modified

| File | Change |
|------|--------|
| `resolve_context.h` | Add `ResolverScratch` struct; add overloaded `_resolve_core` and `compute_frag_lengths` accepting scratch; move existing scratch into internal `scratch_` member |
| `bam_scanner.cpp` | Add `scan_threaded()` method; add `process_qname_group_threaded()`; add `reader_loop()` and `worker_loop()` functions; add `merge_*` functions; add nanobind binding |
| `thread_queue.h` (new) | `BoundedQueue<T>` template, `QnameGroup` struct, `WorkerState` struct |
| `pipeline.py` | Add `n_scan_threads` to `BamScanConfig`; call `scan_threaded` when threads > 1 |
| `config.py` | Add `n_scan_threads` field to `BamScanConfig` |
| `CMakeLists.txt` | Link pthreads (macOS: implicit; Linux: `-lpthread`) |
| `profiler.py` | Pass `n_scan_threads` through to scanner |

## Implementation Order

1. **Phase 1** (resolve_context.h): `ResolverScratch` extraction.
   Build + run all tests. This is a pure refactor with zero behavior change.

2. **Phase 2** (thread_queue.h): Queue + WorkerState + QnameGroup.
   Header-only, no integration yet.

3. **Phase 3** (bam_scanner.cpp): Implement `scan_threaded()`.
   Wire up reader → queue → workers → merge. Expose via nanobind.

4. **Phase 4** (Python side): Config + pipeline + CLI integration.

5. **Phase 5**: Test with golden suite, profile, tune thread counts.

## Risk Mitigations

| Risk | Mitigation |
|------|------------|
| cgranges has hidden global state | Audit cgranges.c — confirmed: `cr_overlap` uses only caller's buffer. Tree is read-only after index. |
| `bam_dup1` allocation pressure | Records are freed immediately after worker processes them. In-flight memory is bounded by queue capacity × group size. Pool allocator can be added later if needed. |
| Fragment order changes break golden tests | Golden tests compare final transcript abundances and EM outputs, not buffer order. Verified: EM is order-independent. |
| Reader becomes I/O-bound (SSD vs HDD) | With `hts_set_threads`, decompression is pipelined. For HDD, this could be the limit. NVMe/tmpfs should be fine. Profile to confirm. |
| GIL interference | All threading is in C++ (no GIL held). nanobind releases GIL for the `scan_threaded` call. |
| htslib not built with thread support | Already verified: conda htslib links libSystem (pthreads on macOS). `hts_set_threads` is available. |
