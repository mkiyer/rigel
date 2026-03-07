# P1 — BAM Scan: Eliminate `bam_dup1` Allocation Overhead

## Current Architecture

```
┌──────────────────────┐     BoundedQueue<QnameGroup>     ┌────────────────────┐
│    Producer (main)    │ ──push(move(group))────────────▸ │   Worker threads   │
│                       │                                  │                    │
│ sam_read1(fp, hdr, b) │     capacity = 2 × n_workers     │ pop(group)         │
│ bam_dup1(b) ──── ✗    │                                  │ parse bam1_t*      │
│ group.records.push()  │                                  │ resolve + buffer   │
│                       │                                  │ ~group destroys    │
│ on qname boundary:    │                                  │   bam_destroy1()   │
│   queue.push(group)   │                                  │                    │
└──────────────────────┘                                  └────────────────────┘
```

### What Happens Per Record (24M reads)

1. **Producer** calls `sam_read1()` into a reused `bam1_t* b`.
2. On each record, producer calls **`bam_dup1(b)`** — this `malloc`s a new
   `bam1_t` (header + variable-length data: sequence, quality, CIGAR, aux tags)
   and `memcpy`s everything. ~200–500 bytes per record.
3. The pointer is pushed into `QnameGroup.records` (a `std::vector<bam1_t*>`).
4. On qname boundary (~2 records/group for PE), the group is **moved** into
   the bounded queue.
5. **Worker** pops the group, parses each `bam1_t*` into a `ParsedAlignment`
   (extracting flag, pos, CIGAR, tags — all scalar fields + small vectors).
6. After parsing, the worker is **done** with the `bam1_t*` — it works
   exclusively with `ParsedAlignment` objects for the rest of processing.
7. `~QnameGroup()` calls `bam_destroy1()` on each pointer → `free()`.

### Cost of `bam_dup1`

For 24M records:
- 24M `malloc` calls (~200–500 bytes each)
- 24M `memcpy` calls
- 24M `free` calls (in worker threads)

Total transient allocation: ~6–12 GB. The malloc/free churn creates
allocator pressure and TLB/cache pollution. On macOS with the system
allocator, this is a significant cost.

### Key Observation

**Workers only need `bam1_t*` for the initial parse step** — to extract
scalars + CIGAR into `ParsedAlignment`. After that, the raw BAM record is
dead weight. The parse is pure extraction (no mutation of `bam1_t`).

---

## Approach A: Pool-Based `bam1_t` Recycling

### Design

Pre-allocate a pool of `bam1_t` objects. Producer borrows from pool, reads
into borrowed record, pushes group. Worker finishes parsing, returns records
to pool.

```
 ┌──────────────┐    BoundedQueue<QnameGroup>    ┌──────────────────┐
 │   Producer    │────push(group)───────────────▸ │     Worker       │
 │               │                                │                  │
 │ pool.borrow() │◂───────return(bam1_t*)────────│ parse bam1_t* →  │
 │ sam_read1(b)  │        RecordPool              │  ParsedAlignment │
 │ group.push(b) │      (lock-free ring)          │ pool.return(b)   │
 └──────────────┘                                └──────────────────┘
```

**Pool implementation**: Lock-free SPMC ring buffer, or simple mutex-guarded
`std::vector<bam1_t*>` (contention is low — workers return records at high
rate, producer borrows at the same rate).

htslib's `sam_read1()` calls `bam_set1()` internally which **reuses** the
`bam1_t`'s `data` buffer if it's large enough (doubling-growth via
`realloc`). So a pooled `bam1_t` will quickly grow to the max record size
and then get reused with zero allocation.

### Changes Required

1. **New class `Bam1Pool`**: ~40 lines.
   - `borrow() → bam1_t*` (pop from free list, or `bam_init1()` if empty)
   - `return_to_pool(bam1_t*)` (push to free list)
   - Destructor calls `bam_destroy1()` on all pooled records.

2. **Modify `QnameGroup`**:
   - Add `Bam1Pool* pool` pointer.
   - Destructor returns records to pool instead of `bam_destroy1()`.
   - Or: workers return records explicitly after parsing, before processing.

3. **Modify `scan_threaded()`**:
   - Create `Bam1Pool` before the loop.
   - Replace `bam_dup1(b)` with `pool.borrow()` + `sam_read1(fp, hdr, pooled)`.
   - Pass pool reference to workers (or embed in QnameGroup).

### Pros
- Eliminates 24M malloc + 24M free.
- htslib's `bam_set1` reuses internal buffer → zero-copy after warmup.
- Records returned to pool after parse → bounded pool size ≈ queue_depth × avg_group_size.
- Simple, well-understood pattern.

### Cons
- Lifetime management: must ensure workers return records before pool is destroyed.
- `QnameGroup` now has a non-owning relationship with records (or needs pool ptr).
- Pool must be thread-safe (producer borrows, workers return).
- Workers must explicitly return records — adds a step after parsing.
- If a worker is slow, pool can starve (producer blocks). Need enough pool capacity.

### Complexity
- ~60 lines new code (pool class).
- ~20 lines changed (scan_threaded reader loop + worker lambda).
- Thread safety: one mutex for pool, or lock-free ring.

---

## Approach B: Parse-in-Producer (Eliminate `bam1_t*` Transfer)

### Design

Move the `bam1_t → ParsedAlignment` extraction into the **producer thread**.
The producer reads `sam_read1()` into a single reused `bam1_t`, immediately
parses it into a `ParsedAlignment`, and puts `ParsedAlignment` objects into
the queue. Workers never see raw `bam1_t*` at all.

```
 ┌────────────────────────┐  BoundedQueue<ParsedGroup>  ┌──────────────────┐
 │     Producer (main)     │────push(group)────────────▸ │     Worker       │
 │                         │                             │                  │
 │ sam_read1(fp, hdr, b)   │                             │ pop(group)       │
 │ ParsedAlignment = parse │                             │ resolve + buffer │
 │ group.records.push(pa)  │                             │                  │
 │ on boundary: push(group)│                             │                  │
 └────────────────────────┘                             └──────────────────┘
```

**No `bam_dup1`, no pool, no record lifecycle management.**

The single-threaded `scan()` path already does exactly this pattern: it reads
into a reused `b`, immediately builds a `ParsedAlignment`, and pushes it into
`current_group`. The threaded path was unnecessarily deferring this parse to
workers.

### Changes Required

1. **New struct `ParsedGroup`** (replaces `QnameGroup`):
   ```cpp
   struct ParsedGroup {
       std::vector<ParsedAlignment> records;
       int64_t frag_id = 0;
       // Move-only, no bam1_t* lifecycle, trivially destructible
   };
   ```
   ~15 lines, replacing ~35 lines of `QnameGroup` (eliminates destructor,
   move ops, bam_destroy1 calls).

2. **Modify `scan_threaded()` reader loop**:
   - Reuse single `bam1_t* b` (already exists for stats).
   - After early filters, build `ParsedAlignment` inline (copy the 15 lines
     from the single-threaded `scan()` path — flag, pos, tags, parse_cigar).
   - Push `ParsedAlignment` into `current_parsed_group.records`.
   - On qname boundary: `queue.push(move(current_parsed_group))`.

3. **Modify `process_qname_group_from_bam1()`**:
   - Rename to `process_parsed_group()`.
   - Remove the 25-line `bam1_t → ParsedAlignment` parsing block at the top.
   - Receives `ParsedGroup&` instead of `QnameGroup&`.
   - Rest of function (resolve, buffer, stats) is **unchanged**.

4. **Update `BoundedQueue` type**:
   - `BoundedQueue<ParsedGroup>` instead of `BoundedQueue<QnameGroup>`.
   - No changes to `BoundedQueue` itself (it's a template).

### Pros
- **Zero raw BAM records cross the thread boundary** — no pool, no lifecycle.
- No `bam_dup1`, no `bam_destroy1`, no pool synchronization.
- `ParsedAlignment` is a value type with `std::vector` members — moves are trivial.
- The parse (CIGAR + tags) is cheap (~100–200ns/record) vs. resolve (~1–5µs).
- Matches the single-threaded code path exactly — easier to reason about.
- `ParsedGroup` destructor is trivial (no raw pointers to manage).
- Removes the `bam1_t*` concept from the threading layer entirely.
- **Thread safety is free** — no shared mutable state beyond the queue.

### Cons
- Producer does slightly more work (CIGAR parse + tag extraction per record).
  However, this was always sequential anyway — the question is whether the
  producer can keep up with htslib I/O + parse vs. worker consumption.
- If CIGAR parsing is slow, the producer becomes the bottleneck. But:
  - `parse_cigar()` is ~100ns (cache-hot small-CIGAR PE reads).
  - `sam_read1()` + BGZF decompression dominates the producer.
  - Workers do resolve + accumulate (~1–5µs/fragment) — much heavier.
  - So the producer has ~10–50× headroom.

### Complexity
- ~15 lines new code (ParsedGroup struct).
- ~35 lines deleted (QnameGroup with its bam1_t* lifecycle).
- ~15 lines moved (parse block from worker → producer).
- ~10 lines changed (worker function signature).
- Net: **~25 fewer lines of code**.

---

## Comparison

| Criterion                    | A: Pool Recycling             | B: Parse-in-Producer           |
|------------------------------|-------------------------------|--------------------------------|
| **malloc/free eliminated**   | ✅ Yes (pool reuse)           | ✅ Yes (no bam1_t* at all)     |
| **memcpy eliminated**        | ❌ No (sam_read1 fills pool)  | ✅ Yes (tag/cigar → scalars)   |
| **New synchronization**      | ⚠️ Pool mutex/lock-free       | ✅ None (queue-only)           |
| **Lifetime management**      | ⚠️ Pool + worker return       | ✅ Trivial (value types)       |
| **Code complexity**          | +60 lines (pool class)       | −25 lines (simpler)            |
| **Producer work increase**   | None (still does dup)         | +100–200ns/record (parse)      |
| **Producer bottleneck risk** | Low                           | Low (10–50× headroom)          |
| **Matches single-thread**    | ❌ Different from `scan()`    | ✅ Identical pattern           |
| **Queue item size**          | ~16 bytes (2 ptrs)           | ~200 bytes (ParsedAlignment)   |
| **Memory pattern**           | Pool sits in cache            | Value-type moves               |
| **Risk of bugs**             | Medium (pool lifetime)        | Low (no raw pointers)          |

---

## Recommendation: **Approach B — Parse-in-Producer**

Approach B is clearly superior for this codebase:

1. **Simplicity**: It eliminates the `bam1_t*` concept from the threading
   layer entirely. No pool, no lifecycle, no synchronization beyond the
   existing bounded queue. The `QnameGroup` struct (with its destructor,
   move ops, and `bam_destroy1` calls) is replaced by a simpler value-type
   struct.

2. **Correctness**: Value-type `ParsedAlignment` objects in
   `std::vector<ParsedAlignment>` have trivial move semantics. No risk of
   use-after-free, double-free, or pool exhaustion.

3. **Already proven**: The single-threaded `scan()` path does exactly this —
   parse in the reader loop, process `ParsedAlignment` objects. The threaded
   path just needs the same structure.

4. **Smaller diff**: Net deletion of ~25 lines. The parse block is copied
   from the single-threaded path (or factored into a shared helper).

5. **No producer bottleneck**: `parse_cigar()` + tag extraction is ~100–200ns
   per record. Workers do resolve + accumulate at ~1–5µs per fragment. With
   4 workers, the consumer side processes at ~4M fragments/s. The producer
   (sam_read1 + parse) runs at ~20–30M records/s (CIGAR-parsing benchmark)
   but is I/O-bound by BGZF decompression at ~10–15M records/s. So the
   producer easily keeps workers fed.

### Implementation Plan

**Files to modify**: `thread_queue.h`, `bam_scanner.cpp`

**Step 1**: Add `ParsedGroup` struct to `thread_queue.h` (or `bam_scanner.cpp`):
```cpp
struct ParsedGroup {
    std::vector<ParsedAlignment> records;
    int64_t frag_id = 0;
};
```

**Step 2**: Factor out `bam1_t → ParsedAlignment` into a static helper:
```cpp
static ParsedAlignment parse_bam_record(
    bam1_t* b,
    const std::vector<int32_t>& tid_to_ref_id,
    SJTagMode sj_tag_mode);
```
This extracts the 15-line block that's duplicated between `scan()` and
`process_qname_group_from_bam1()`. Both call sites become one-liners.

**Step 3**: Modify `scan_threaded()` reader loop:
- Change queue type: `BoundedQueue<ParsedGroup>`
- Replace `bam_dup1(b)` + group.records.push_back with
  `parsed_group.records.push_back(parse_bam_record(b, tid_to_ref_id_, sj_tag_mode_))`
- On qname boundary: `queue.push(std::move(current_parsed_group))`

**Step 4**: Modify worker function:
- Rename `process_qname_group_from_bam1` → `process_parsed_group`
- Accept `ParsedGroup&` instead of `QnameGroup&`
- Remove the 25-line `bam1_t → ParsedAlignment` parsing loop at the top
- Rest of function unchanged

**Step 5**: Remove `QnameGroup` from `thread_queue.h` (now unused).

**Step 6**: Update single-threaded `scan()` to also use `parse_bam_record()`
helper for consistency (optional but nice).

### Expected Outcomes

- **Runtime**: Eliminate ~5–15s from scan_and_buffer (from 48s baseline).
  Conservative estimate: 10–20% reduction in scan time.
- **Memory**: Eliminates transient 6–12 GB of bam_dup1 allocations.
  Pool-equivalent memory is ~0 (ParsedAlignment lives briefly in queue).
- **Code**: Net reduction of ~25 lines. Simpler threading model.
- **Test impact**: No behavioral change — same ParsedAlignment data reaches
  workers. Downstream resolve + accumulate is identical.
