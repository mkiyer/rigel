# PR 01: Remove Synchronous Scan Spill From The Critical Path

## Summary

Make temporary `FragmentBuffer` spill non-blocking or avoid it entirely when the
next stage can consume chunks once. The goal is to keep the C++ scanner running
while older chunks are written, loaded, scored, or released.

## Motivation

The profiling runs show that spill is not a minor accounting detail. The same
VCAP workload moved from 200.9s with 14 spilled chunks to 135.7s when the memory
budget was large enough to avoid spill. The current spill path runs inside the
Python chunk callback:

```text
C++ scanner -> finalize_zero_copy -> Python _on_chunk -> buffer.inject_chunk
             -> _accept_chunk -> _spill_oldest -> pyarrow.feather.write_feather
```

While this executes, the scan pipeline cannot accept output chunks normally.
Live samples showed main-thread stacks in Arrow IPC serialization and LZ4
compression while scanner threads were waiting.

## Current Code

- Chunk callback wiring: [src/rigel/pipeline.py](../../../src/rigel/pipeline.py)
- Buffer spill logic: [src/rigel/buffer.py](../../../src/rigel/buffer.py)
- Native output queue and finalization:
  [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)

## Proposed Scope

This PR should implement the smallest design that removes synchronous spill from
the scan callback while preserving the current public `FragmentBuffer` API.

Recommended design:

1. Add a background spill writer to `FragmentBuffer`.
2. Replace spilled chunk entries with a small pending-spill handle containing:
   path, completion event, exception slot, and memory accounting metadata.
3. In `_spill_oldest`, enqueue the chunk for the writer and immediately subtract
   in-memory bytes from the buffer budget.
4. In `iter_chunks`, `iter_chunks_consuming`, and `cleanup`, wait on pending
   spill handles only when their data is needed.
5. Propagate writer exceptions on the next buffer operation.
6. Join the writer during `cleanup`, `release`, and context-manager exit.

Secondary option for the same PR or a follow-up:

- Add a `spill_compression` option with default `lz4`, but allow `None` for
  temporary local spill. If disk bandwidth is sufficient, uncompressed IPC may be
  faster and use less CPU than LZ4.

Larger follow-up, not required in this PR:

- Stream scan chunks directly into scoring and release them, avoiding temporary
  disk serialization for the normal quantification path.

## Implementation Steps

1. Introduce a private `_PendingSpill` dataclass in `buffer.py`.
2. Add a writer thread and queue initialized lazily when the first spill occurs.
3. Update `_chunks` to allow `_FinalizedChunk`, `Path`, and `_PendingSpill`.
4. Update `_spill_oldest` to enqueue rather than write synchronously.
5. Update `_load_chunk` call sites to wait if the chunk reference is pending.
6. Make `cleanup` robust when the writer has pending chunks.
7. Add logging that distinguishes `Queued spill` from `Finished spill`.

## Tests

Add focused tests in [tests/test_buffer.py](../../../tests/test_buffer.py):

- Spilling with a tiny memory budget produces pending or completed spill entries.
- Iteration returns the same chunks and fragment rows as synchronous spill.
- Writer exceptions are surfaced to the caller.
- `cleanup()` removes files and joins the writer.
- `iter_chunks_consuming()` deletes spilled files after use.

Run:

```bash
conda activate rigel
pytest tests/test_buffer.py tests/test_pipeline_smoke.py tests/test_golden_output.py -v
```

## Benchmark Plan

Run the scan profiler with a memory budget that forces spill:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr01 \
  --name-prefix pr01 \
  --n-scan-threads 4 8 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 4
```

Compare against:

- `s8 d2, 4GiB spill`: 200.9s, 14 spills
- `s4 d2, 12GiB no-spill`: 135.7s, 0 spills

## Acceptance Criteria

- Golden outputs are unchanged.
- Forced-spill scan wall time improves materially versus 200.9s.
- The scan main thread no longer shows sustained Arrow/LZ4 stacks in a live
  sample during spill bursts.
- Peak RSS remains bounded by the configured memory budget plus writer queue
  headroom.

## Risks

- Async writer errors can otherwise be delayed or hidden. Treat this as a core
  correctness requirement.
- A writer queue with too much buffering can silently increase RSS. Keep the
  queue bounded.
- If scoring immediately needs a chunk whose spill is still running, the wait is
  only moved, not removed. The main win is preserving scan throughput.

## Non-Goals

- Do not introduce Polars here. The issue is synchronous serialization, not the
  dataframe engine.
- Do not redesign the scoring pipeline in this PR unless the async spill design
  proves insufficient.
