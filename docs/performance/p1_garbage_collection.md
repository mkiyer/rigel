# P1: Garbage Collection — Analysis and Implementation Plan

**Date**: 2026-04-14
**Status**: Proposed
**Estimated savings**: 8–12s (7–10% of total wall time)

---

## 1. Should We Be Managing Garbage Collection Explicitly?

**No.** The explicit `gc.collect()` calls are unnecessary and actively harmful. Here's why:

### The root misunderstanding

The current code assumes that setting an attribute to `None` and calling `del` is not sufficient to free a large numpy array — that an explicit `gc.collect()` is needed to "force" deallocation. This is incorrect.

CPython uses **reference counting** as its primary memory management mechanism. When the last reference to an object drops (via `del`, `= None`, or scope exit), the object is deallocated **immediately** — synchronously, deterministically, at that exact line of code. No garbage collector intervention is needed.

Python's cyclic garbage collector (`gc` module) exists for exactly one purpose: reclaiming objects trapped in **reference cycles** — situations where A→B→A and neither can reach refcount zero despite being unreachable. The collector discovers these cycles by scanning all tracked objects and identifying unreachable clusters.

### Why the current gc.collect() calls are pure waste

The objects being freed in `partition_and_free()` are:

| Object | Type | GC-tracked? | Can form cycles? |
|--------|------|-------------|------------------|
| `ScoredFragments` | `@dataclass(slots=True)` | Yes | **No** — slots hold only `ndarray` and `int` |
| numpy arrays (the actual memory) | `np.ndarray` | **No** | No — C buffer, no PyObject refs |
| `int`, `float` scalars | built-in | No | No |

A numpy array holds a C-level `void*` data buffer. It does not contain pointers to other Python objects. It is not tracked by the cyclic GC. When its refcount reaches zero — which happens the instant `setattr(em_data, attr, None)` executes — the buffer is freed by `PyMem_Free`. No `gc.collect()` delay, no "pending" state.

The `ScoredFragments` dataclass uses `slots=True`, which prevents `__dict__` creation. Its slots hold only `ndarray` and `int` — neither of which is GC-tracked. Without `__dict__` and without any tracked child objects, `ScoredFragments` itself **cannot participate in a reference cycle**. It, too, is freed by refcount alone.

### What gc.collect() actually does when called

`gc.collect()` performs a **full generational sweep** across all three generations (gen0 + gen1 + gen2). This means it must visit every GC-tracked object on the heap. In rigel's pipeline, the heap contains:

| Source | Tracked objects | Notes |
|--------|----------------|-------|
| Collapsed interval index | 1,220,804 frozensets | Transcript sets per collapsed interval |
| Splice junction map | 404,168 frozensets | Transcript sets per junction |
| Module/class infrastructure | ~50,000 | Functions, types, descriptors, etc. |
| **Total** | **~1,675,000** | |

Measured cost: **567ms per `gc.collect()`** with this heap.

The current code calls `gc.collect()` **15 times** inside `partition_and_free()` plus 4 more in `pipeline.py`:

```
partition_and_free:  6 (candidate arrays) + 1 (offsets) + 8 (unit arrays) = 15
pipeline.py:         1 (after partition) + 1 (after mega-locus) + 1 (after normal loci) + 1 (cleanup) = 4
Total:               19 × 567ms ≈ 10.8s
```

The profiler measured the partition phase at **13.2s**, of which ~4–5s is the actual scatter work (memcpy) and **~8–9s is pure gc.collect() overhead** — scanning 1.6M frozensets that have nothing to do with the arrays being freed.

### What about automatic GC?

Python's automatic GC uses generational thresholds (default: 700/10/10). Gen0 collections trigger after 700 net new allocations and scan only *young* objects — not the 1.6M frozensets in gen2. Gen0 sweeps cost microseconds, not hundreds of milliseconds.

The 1.6M frozensets migrated to gen2 during index load (they survived many collections). They will only be revisited during automatic gen2 sweeps, which trigger every ~70,000 net allocations — well outside the partition hot path.

**Conclusion**: The explicit `gc.collect()` calls should be removed entirely. Refcounting already provides immediate, deterministic deallocation. The GC sweeps serve no purpose and cost ~9s per pipeline run.

---

## 2. The Elegant Solution

The correct approach is not to "manage GC better" — it is to **stop fighting the runtime**. CPython's reference counting already does the right thing. The solution removes code rather than adding it.

### Principle: Trust refcounting for acyclic objects

If an object graph is acyclic — and ours is (numpy arrays can't form cycles, slots dataclasses holding only arrays and ints can't form cycles) — then refcounting provides:

1. **Immediate deallocation** at the exact point where the last reference drops
2. **Deterministic timing** — no GC pause jitter
3. **Zero overhead** — no heap scan, no generational bookkeeping

This is strictly superior to explicit `gc.collect()` in every dimension: faster, more predictable, and simpler.

### The change

Remove all `gc.collect()` calls from `partition.py`. Remove the `gc.collect()` calls from `pipeline.py` that exist solely to reclaim these same objects. Keep the refcount-based cleanup patterns (`setattr(em_data, attr, None)`, `del`, `= None`) — these are correct and sufficient.

Specifically:

**partition.py** — Remove `import gc` and all 15 `gc.collect()` calls. The `setattr(em_data, attr, None)` + `del global_arr` pattern already frees by refcount:

```python
# Before (current)
for attr, scatter_fn in CAND_ARRAYS:
    global_arr = getattr(em_data, attr)
    cand_results[attr] = scatter_fn(...)
    setattr(em_data, attr, None)
    del global_arr
    gc.collect()          # ← 567ms of wasted work, 15 times

# After
for attr, scatter_fn in CAND_ARRAYS:
    global_arr = getattr(em_data, attr)
    cand_results[attr] = scatter_fn(...)
    setattr(em_data, attr, None)  # refcount drops → immediate free
    del global_arr                # belt-and-suspenders (optional)
```

**pipeline.py** — Remove the 4 `gc.collect()` calls in `quant_from_buffer()` and `_run_locus_em_partitioned()`:

| Location | Current | After |
|----------|---------|-------|
| After `partition_and_free` | `del em_data; gc.collect()` | `del em_data` |
| After mega-locus `del part` | `del part; gc.collect()` | `del part` |
| After `del partitions` | `del partitions; gc.collect()` | `del partitions` |
| Phase 6 cleanup | `gc.collect()` | *(remove)* |

### What about edge cases?

**Q: What if a future code change introduces a cycle involving these objects?**
The automatic generational GC will catch it — that's what it's for. Cycles involving young objects are collected by gen0 within ~700 allocations, at negligible cost. No explicit calls needed.

**Q: Could RSS spike without the intermediate gc.collect() calls?**
No. The arrays are freed immediately when their refcount drops. The only case where RSS wouldn't decrease is if the OS doesn't actually reclaim the pages (common with `mmap`-backed allocations). But these are heap-allocated numpy arrays using `PyMem_Malloc` / `malloc` — they are returned to the allocator immediately and available for reuse. The RSS measurement reflects this.

**Q: What about the `LocusPartition` objects created during assembly?**
These are also `@dataclass(slots=True)` holding numpy arrays + ints. When the partitions dict is popped or deleted, refcounting handles it. The auto-GC handles any theoretical cycles.

**Q: Why were the gc.collect() calls added in the first place?**
Likely a precautionary measure during development, when the memory lifecycle was less clear. It's a common Python pattern — "when in doubt, gc.collect()." In most Python programs this costs 1–5ms and is harmless. Here, with 1.6M GC-tracked frozensets from the index, it costs 567ms per call.

---

## 3. Estimated Impact

| Metric | Before | After | Δ |
|--------|--------|-------|---|
| `partition` stage | 13.2s | ~4.5s | **−8.7s (−66%)** |
| Pipeline `gc.collect()` | ~2.3s | 0s | **−2.3s** |
| **Total wall** | 112.1s | ~101s | **−11s (−10%)** |
| Peak RSS | 8,774 MB | 8,774 MB | unchanged |
| Code lines | 19 gc.collect() + import | 0 | −20 lines |
| Throughput | 372K frags/s | ~413K frags/s | **+11%** |

The partition stage estimate (4.5s) is the observed 13.2s minus the measured gc overhead (~8.5s from 15 × 567ms) rounded conservatively. The pipeline gc.collect() calls contribute an additional ~2.3s (4 × 567ms).

### Risk assessment

**Risk: None.** This change removes dead code — calls that scan objects uninvolved in the freed memory. The numpy arrays are already freed by refcount at the `setattr(em_data, attr, None)` line. The gc.collect() runs *after* the memory is already freed and accomplishes nothing.

---

## 4. A Note On The Index Frozensets

The 1.6M frozensets (collapsed intervals + splice junctions) are not a bug — they are a valid data structure for transcript set lookups. However, they are the indirect cause of the gc.collect() overhead: each full-generation sweep must visit all of them.

A separate future optimization could convert these to a more GC-friendly representation (e.g., numpy arrays of transcript indices, or a C++ data structure). This would reduce the base cost of any accidental gc.collect() call and speed up the rare automatic gen2 sweeps. But this is orthogonal to P1 — the correct fix is to stop calling gc.collect() where it isn't needed.

---

## 5. Implementation Checklist

1. **partition.py**: Remove `import gc`. Remove all `gc.collect()` calls from both scatter loops and the offsets-free line. Keep `setattr(..., None)` and `del` statements.

2. **pipeline.py**: Remove the 4 `gc.collect()` calls in:
   - `quant_from_buffer()` after `del em_data` (line ~786)
   - `_run_locus_em_partitioned()` after `del part` (line ~613)
   - `_run_locus_em_partitioned()` after `del partitions` (line ~657)
   - `quant_from_buffer()` Phase 6 cleanup (line ~811)
   - Remove `import gc` if no other uses remain.

3. **Tests**: Run full suite (`pytest tests/ -v`). No test should depend on gc.collect() timing.

4. **Profile**: Rerun stage profiler on VCaP BAM. Verify partition stage drops from ~13s to ~4–5s.

5. **Validate RSS**: Confirm peak RSS is unchanged via memory timeline CSV.
