# Fused Score Buffer — Two-Pass C++ Architecture

**Date:** 2026-03-06  
**Status:** Approved, implementing

## Problem

The current `FragmentRouter._scan_native()` path copies every scored
candidate **three times** between C++ and Python:

```
C++ scan_chunk → numpy (f64)                      [C++ allocation]
  → .astype(f32).tobytes() → array.array.frombytes [COPY 1+2: into Python accumulator]
  → np.frombuffer(...).copy()                      [COPY 3: into final numpy]
```

At 24M fragments (~50M candidates), the `array.array` accumulators and
the final numpy arrays **coexist** during finalization, doubling peak
memory to ~5.2 GB for the scanner output alone.

## Solution: Two-Pass Fused Scoring in C++

One C++ function replaces the entire Python accumulation loop.

### Pass 1 — Count

Iterate all buffer chunks through the existing scoring + WTA gating
logic.  Only **count** survivors (n_units, n_candidates).  Unambiguous
fragment counts (`unambig_counts`, `exonic_sense`, etc.) are accumulated
as a side effect — this is the only mutation.

### Allocate

Pre-allocate all 15 output arrays at their exact final sizes using
C++ `std::vector`.

### Pass 2 — Fill

Re-iterate the same chunks, re-run scoring + WTA gating, and write
results directly into the pre-allocated arrays at write cursors.
Skip unambig count accumulation (already done in pass 1).

### Return

Hand the `std::vector` storage to Python as numpy arrays via
nanobind capsules (zero-copy).  Python wraps them in `ScoredFragments`.

## Memory Budget (24M fragments)

| Stage | Current | After |
|---|---:|---:|
| array.array accumulators | ~2,600 MB | **0 MB** |
| _to_np transient copy | ~2,600 MB | **0 MB** |
| ScoredFragments (final) | ~2,600 MB | ~2,600 MB |
| **Peak during scan** | **~5,200 MB** | **~2,600 MB** |

## Float32 Throughout

All floating-point output arrays (`log_liks`, `coverage_weights`,
`gdna_log_liks`) are produced as `float` (32-bit) natively in C++.
The EM solver already accepts `const float*`.

## What Changes

### New in C++ (`scoring.cpp`)

- `fused_score_buffer()` — top-level function exposed to Python.
  Accepts list of chunks + scoring context + estimator accumulators.
  Returns tuple of 15+ numpy arrays.
- Refactors `scan_chunk` internals into reusable per-fragment helpers
  so both count and fill passes share the same scoring logic without
  code duplication.

### Modified in Python (`scan.py`)

- `FragmentRouter._scan_native()` calls `fused_score_buffer()` once
  instead of the per-chunk accumulation loop.
- Handles multimapper groups in Python (appends to result arrays).
- Constructs `ScoredFragments` from returned numpy views.

### Deleted from Python (`scan.py`)

- All 15 `array.array` accumulators in `FragmentRouter.__init__`.
- `_to_np()` helper.
- Per-chunk `.astype().tobytes()` → `.frombytes()` accumulation loop.

### Unchanged

- `FragmentBuffer` and its spill/reload infrastructure.
- `ScoredFragments` dataclass.
- `FragmentRouter._scan_python()` fallback for tests/mocks.
- `FragmentRouter._flush_mm_group()` (Python multimapper path).

## Multimapper Handling — Future Work

Multimapper fragments (NH > 1) are currently handled entirely in
Python via `FragmentRouter._flush_mm_group()`.  These fragments
accumulate across chunks (same `frag_id`, different alignments) and
require cross-chunk state.  They represent < 1% of fragments in
typical RNA-seq libraries.

For this iteration, multimappers remain in Python.  The C++ fill pass
leaves headroom for MM units, and Python appends them after the C++
pass.

**Future iteration:** Port `_flush_mm_group()` and the cross-chunk MM
state machine to C++, eliminating the last Python-side accumulation.
This will require the C++ function to maintain a pending MM buffer
across chunk boundaries and flush on frag_id transitions.

## Implementation Checklist

1. [ ] Add `fused_score_buffer()` to `scoring.cpp` (count + alloc + fill)
2. [ ] Expose via nanobind module definition
3. [ ] Update `scan.py` `_scan_native()` to call new function
4. [ ] Handle MM append in Python
5. [ ] Delete dead Python accumulation code
6. [ ] Golden test validation (21/21)
7. [ ] Full test suite validation
8. [ ] Profile peak RSS reduction
