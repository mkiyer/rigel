# P5: Streaming Buffer Release — Eliminating the Peak Memory Cliff

**Date**: 2026-04-14
**Status**: Proposed
**Estimated savings**: 2.6 GB peak RSS (8.8 GB → 6.2 GB), ~30% reduction

---

## 1. The Problem

Rigel's memory profile has a cliff:

```
                                  ┌─── 8.8 GB peak (buffer + CSR coexist)
                                  │
   8 GB ─                    ┌────┤
                             │CCCC│    C = CSR arrays (3.2 GB)
   6 GB ─  ┌────────────────┤CCCC│
           │BBBBBBBBBBBBBBBB│CCCC├────┐
   4 GB ─  │BBBBBBBBBBBBBBBB│CCCC│    │    B = fragment buffer (2.7 GB)
           │BBBBBBBBBBBBBBBB│CCCC│    │
   2 GB ─  │BBBBBBBBBBBBBBBB│CCCC│    │
       ────┘BBBBBBBBBBBBBBBB└────┘    └────
           ↑                ↑    ↑
         scan          score+route  buffer.release()
```

The fragment buffer (~2.7 GB) and scored CSR arrays (~3.2 GB) coexist during the `fragment_router_scan` phase. The buffer persists through the entire scoring loop even though each chunk is consumed exactly once and never read again. Only after `scorer.finish()` returns the complete CSR does `buffer.release()` free the buffer — at which point peak has already been set.

This is wasteful. The scoring loop already processes chunks one at a time in sequence. Each chunk is converted to scoring arrays, passed to C++, then `del arrays` frees the temporary copy. But the original chunk data sits in `buffer._chunks` untouched, waiting for a bulk `release()` call that comes too late.

For the VCaP benchmark (19.4M fragments):

| Phase | RSS | Components |
|-------|-----|-----------|
| After BAM scan | 5,573 MB | index (2,835) + buffer (2,737) |
| After scoring/routing | 8,792 MB | index + buffer + CSR (3,220) |
| After buffer release | 5,467 MB | index + CSR |
| After EM | 5,249 MB | index + results |

The 8.8 GB peak is set by ~10 seconds of coexistence during scoring. All 2.7 GB of buffer would be freeable during this window.

---

## 2. The Insight

The scoring loop in `FragmentRouter._scan_native()` already streams:

```python
for chunk in buffer.iter_chunks():
    arrays = chunk.to_scoring_arrays()
    scorer.score_chunk(arrays)
    del arrays
```

Each chunk is visited exactly once, in order, never revisited. The `del arrays` line frees the *upcast copies* created by `to_scoring_arrays()`, but the **original chunk** (the `_FinalizedChunk` object with its 17 numpy arrays) stays alive in `buffer._chunks`.

The fix is to release each chunk **after it has been scored**, not after all chunks are scored. This turns the memory cliff into a gentle handoff — as the CSR grows, the buffer shrinks.

```
   After (streaming release):

   6 GB ─            ┌────────────
                     │CCCCCCCCCCCC    C = CSR (grows)
   4 GB ─  ┌────┬───┤CCCCCCCCCCCC
           │BBBB│B  C│CCCCCCCCCCCC
   2 GB ─  │BBBB│B C │CCCCCCCCCCCC
       ────┘BBBB└BC──┘CCCCCCCCCCCC────
           ↑          ↑
         scan     score (streaming)
```

Peak memory: index + one chunk + growing CSR. At the end of scoring, the CSR is full (~3.2 GB) but the buffer is already gone. New peak ≈ 6.2 GB.

---

## 3. Design

### 3.1 The core change: `iter_chunks_consuming()`

Add a single new method to `FragmentBuffer` that yields chunks and releases them:

```python
def iter_chunks_consuming(self) -> Iterator[_FinalizedChunk]:
    """Yield chunks one at a time, releasing each after the caller advances.

    After this method returns, the buffer is empty. Spilled chunk files
    are deleted as they are consumed.

    This is the streaming counterpart to iter_chunks(). Use it when
    chunks will not be revisited — the typical case during scoring.
    """
    while self._chunks:
        chunk_ref = self._chunks[0]
        if isinstance(chunk_ref, Path):
            chunk = _load_chunk(chunk_ref)
            chunk_ref.unlink(missing_ok=True)
        else:
            chunk = chunk_ref
        self._chunks[0] = None       # drop reference from list
        self._chunks.pop(0)           # shrink list
        self._memory_bytes -= chunk.memory_bytes
        yield chunk
        # Caller has advanced past this chunk.
        # The chunk's numpy arrays are freed when 'chunk'
        # goes out of scope at the next iteration.
```

**Why a new method, not modifying `iter_chunks()`?**

`iter_chunks()` is a pure, non-destructive view — it yields without side effects. This is the correct default behavior. A consuming iterator is a fundamentally different contract: it destroys the buffer as a side effect of reading. Callers opt in explicitly by choosing `iter_chunks_consuming()`. The two methods coexist cleanly.

**Refinement**: Using `pop(0)` on a list is O(n). With 16 chunks this is negligible, but for elegance we can use `collections.deque`:

```python
# In __init__:
self._chunks: deque[_FinalizedChunk | Path] = deque()

# In iter_chunks_consuming:
chunk_ref = self._chunks.popleft()  # O(1)
```

This also makes `_accept_chunk` (which appends) and `_spill_oldest` (which iterates from left) more natural.

### 3.2 Wire into the scoring loop

In `FragmentRouter._scan_native()`, replace the buffer iteration:

```python
# Before
for chunk in buffer.iter_chunks():
    arrays = chunk.to_scoring_arrays()
    scorer.score_chunk(arrays)
    ...
    del arrays

# After
for chunk in buffer.iter_chunks_consuming():
    arrays = chunk.to_scoring_arrays()
    scorer.score_chunk(arrays)
    ...
    del arrays, chunk     # explicit: chunk freed here
```

The `del chunk` is technically redundant (the loop variable is overwritten at the next iteration), but it documents intent: this chunk will not be revisited.

### 3.3 Remove the separate buffer.release() call

In `pipeline.py::_score_fragments()`:

```python
# Before
em_data = builder.scan(buffer, log_every)
del builder, ctx
buffer.release()           # ← was the only place buffer got freed
return em_data

# After
em_data = builder.scan(buffer, log_every)
del builder, ctx
# buffer is already consumed — chunks released during scan
return em_data
```

The buffer object itself (the shell with an empty deque) can be garbage-collected when `_score_fragments` returns. No explicit release needed.

### 3.4 Memory accounting

`FragmentBuffer.memory_bytes` already tracks the sum of in-memory chunk sizes. `iter_chunks_consuming()` decrements this counter as each chunk is released. By the time the iterator is exhausted, `memory_bytes == 0`. External monitoring (profiler memory snapshots) will see a gradual decline instead of a cliff.

---

## 4. What Changes In Each File

### buffer.py — 3 changes

1. **Change `_chunks` from `list` to `deque`** — O(1) popleft for consuming iteration. The existing `_spill_oldest()` method iterates the deque looking for the first in-memory chunk, which is still O(n) in the worst case but unchanged in practice.

2. **Add `iter_chunks_consuming()` method** — ~15 lines. Pops from the left, loads spilled chunks, yields, decrements memory counter.

3. **Adjust `release()` to be idempotent** — Already is, but verify it handles an empty deque gracefully. Add a guard for `_memory_bytes` already being zero.

### scan.py — 2 changes

1. **`_scan_native()`: Replace `buffer.iter_chunks()` with `buffer.iter_chunks_consuming()`**. Add `del chunk` after `del arrays`. No other changes.

2. **Pass buffer by reference** (already the case) so the consuming iterator modifies the caller's buffer object.

### pipeline.py — 1 change

1. **`_score_fragments()`: Remove `buffer.release()`** — the buffer is already consumed. Keep `del builder, ctx`.

### profiler.py — 1 change

1. **Adjust RSS snapshot placement** — The `after_buffer_release` snapshot should be taken after `builder.scan()` returns, since the buffer is now consumed during scanning. The `after_router_scan` snapshot already captures this. Consider adding intermediate snapshot labels for clarity.

---

## 5. The Bigger Picture: Architectural Coherence

This change aligns with a broader principle already emerging in the codebase: **data flows forward, never backward**. Each stage consumes its input and produces output for the next stage. Nothing is revisited.

The pipeline already follows this in several places:
- The BAM scanner streams fragments, never re-reads
- `partition_and_free()` destroys ScoredFragments as it scatters
- Mega-locus partitions are popped and freed after EM

The buffer → CSR transition was the remaining exception. After this change, the full pipeline becomes:

```
BAM file  →  BAM scanner  →  [buffer chunks]  →  scorer  →  [CSR]  →  loci  →  [partitions]  →  EM
              consume           consume             consume           consume           consume
```

Every arrow is a consuming handoff. No stage's output outlives the next stage's start. This is the minimum-memory pipeline.

### The `list` → `deque` migration

Changing `_chunks` from `list` to `deque` is a small structural improvement with broad benefits:

| Operation | list | deque |
|-----------|------|-------|
| `append()` (right) | O(1) amortized | O(1) |
| `pop(0)` / `popleft()` | **O(n)** | O(1) |
| `[i]` random access | O(1) | O(1)* |
| iteration | O(n) | O(n) |

\* `deque` supports O(1) indexing for reasonable sizes via `__getitem__`.

The only place that does indexed access is `_spill_oldest()` (`self._chunks[i] = path`), which replaces a chunk with its spill path. With a deque this remains valid since deque supports `__setitem__`.

---

## 6. Estimated Impact

### Memory

| Metric | Before | After | Δ |
|--------|--------|-------|---|
| Peak RSS | 8,792 MB | ~6,226 MB | **−2,566 MB (−29%)** |
| Buffer at peak | 2,737 MB | ~171 MB (1 chunk) | −2,566 MB |
| CSR at peak | 3,220 MB | 3,220 MB | unchanged |

The new peak occurs at the end of scoring: index (2,835 MB) + last chunk (~171 MB) + full CSR (3,220 MB) = ~6,226 MB.

### Performance

**Expected: neutral or slight improvement.**

The scoring loop does the same work in the same order. The only difference is that chunk objects are freed earlier, which reduces the allocator's live-set and may improve cache behavior. The `deque.popleft()` is cheaper than not freeing at all.

Spilled chunks are loaded identically — `_load_chunk()` reads from Arrow IPC. The only new cost is `unlink()` on the spill file after load, which is negligible (one syscall per spilled chunk, ~0 chunks for the VCaP benchmark which fits in the 2 GB budget).

### Scaling

The savings scale linearly with dataset size:

| BAM size | Buffer | CSR | Peak (before) | Peak (after) | Savings |
|----------|--------|-----|---------------|--------------|---------|
| 19M frags (VCaP) | 2.7 GB | 3.2 GB | 8.8 GB | 6.2 GB | 2.6 GB |
| 50M frags | ~7 GB | ~8 GB | ~18 GB | ~11 GB | ~7 GB |
| 100M frags | ~14 GB | ~16 GB | ~33 GB | ~19 GB | ~14 GB |

For large clinical samples (100M+ fragments), this is the difference between fitting in 32 GB and requiring 64 GB.

---

## 7. What About P5-Original: Streaming CSR to Disk?

The original P5 proposal in the profiling report suggested a more aggressive approach: spill the *scored CSR* to disk (Arrow IPC with LZ4) during scoring, then reload for partition. This would cap peak at ~max(chunk, CSR\_chunk) ≈ 5.5 GB.

We deliberately chose not to pursue that approach for the following reasons:

1. **Complexity**: Spilling the CSR requires designing a new serialization format for the 15 scored arrays, managing temp files, and handling the reload path. The consuming iterator adds ~15 lines.

2. **I/O cost**: The CSR for VCaP is 3.2 GB. Even with LZ4, writing and reading this from disk adds several seconds (sequential I/O at ~2 GB/s = ~3s write + 3s read). The consuming iterator adds 0 I/O.

3. **Diminishing returns**: The consuming iterator already recovers 2.6 GB of the 2.7 GB buffer. The remaining chunk (~171 MB) is noise. The additional ~170 MB savings from CSR spill doesn't justify the complexity.

4. **Simplicity**: The consuming iterator is an obvious refinement — it makes the existing streaming pattern complete. CSR spill is a fundamentally new mechanism.

If even larger datasets require further memory reduction, CSR spill can be added later as an orthogonal layer. The consuming iterator is a prerequisite for it regardless.

---

## 8. Edge Cases and Safety

**Empty buffer**: `iter_chunks_consuming()` on an empty deque yields nothing. The scorer loop body never executes. `scorer.finish()` returns empty arrays. Everything works.

**Single chunk**: One chunk is popped, scored, freed. CSR contains all units. No issue.

**All chunks spilled**: Each spilled chunk is loaded from disk, yielded, and its file is deleted. The `_temp_dir` may become empty. `cleanup()` handles this (rmtree of an empty dir is fine).

**Mixed spilled + in-memory**: The deque contains a mix of `_FinalizedChunk` objects and `Path` references. `popleft()` returns whichever is next. Type dispatch in the iterator handles both.

**Caller doesn't exhaust the iterator**: If someone calls `iter_chunks_consuming()` and breaks early, remaining chunks stay in the deque. The buffer remains usable (but partially consumed). `release()` still works as a fallback. No leak.

**Thread safety**: Not required — the pipeline is single-threaded in Python. The C++ scoring happens under the GIL.

---

## 9. Implementation Checklist

1. **buffer.py**: Change `_chunks` from `list` to `deque`. Add `iter_chunks_consuming()`. Verify `release()`, `_spill_oldest()`, `_accept_chunk()`, `summary()`, `inject_chunk()` work with deque.

2. **scan.py**: In `_scan_native()`, replace `buffer.iter_chunks()` with `buffer.iter_chunks_consuming()`. Add `del chunk` after `del arrays`.

3. **pipeline.py**: In `_score_fragments()`, remove `buffer.release()`.

4. **profiler.py**: Adjust `after_buffer_release` RSS snapshot to be taken after `builder.scan()` (since buffer is consumed during scan, not after).

5. **Tests**: Existing `test_buffer.py` tests cover `iter_chunks()`, spill, and cleanup. Add tests for `iter_chunks_consuming()`:
   - Buffer is empty after consuming iteration
   - `memory_bytes == 0` after consuming iteration
   - Spilled files are deleted during consuming iteration
   - Partial consumption leaves remaining chunks available

6. **Profile**: Rerun stage profiler on VCaP BAM. Verify peak RSS drops from ~8.8 GB to ~6.2 GB. Verify wall time is unchanged.
