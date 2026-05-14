# PR 04: Async Chunk Store For Spill

**Position in roadmap:** Fourth. Land after PR 01–03 establish the new
scan baseline, or pull earlier if forced-spill workloads are the current
release blocker.

## Summary

Keep Rigel's single-pass scan/training/calibration architecture, but
remove Arrow/LZ4 spill writes from the scanner callback path. The C++
scanner already emits finalized chunks through a Python callback; this
PR makes `FragmentBuffer.inject_chunk(...)` cheap by queuing spill work
to a bounded background writer and surfacing any writer error before the
buffer is consumed.

This replaces the direct scan-to-score streaming plan. That plan is not
valid yet because scoring depends on finalized strand models and v6 FL /
gDNA calibration, all of which are produced after the scan completes.

## Motivation

Current path:

```text
C++ scanner output queue
  -> finalize_zero_copy
  -> Python chunk_callback
  -> _FinalizedChunk.from_raw
  -> FragmentBuffer.inject_chunk
  -> _spill_oldest
  -> pyarrow.feather.write_feather(..., compression="lz4")
```

When the memory budget is exceeded, the callback blocks inside Arrow IPC
serialization and LZ4 compression. The forced-spill VCAP run paid about
65 s versus the no-spill run.

We cannot simply score chunks as they arrive because the scorer needs
models and calibration derived from the full scan. The right local fix is
therefore an asynchronous chunk store: preserve the current semantics,
but move temporary chunk serialization off the scanner's critical path.

## Current Code

* Scanner callback already exists: [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)
* Python callback and buffer injection: [src/rigel/pipeline.py](../../../src/rigel/pipeline.py)
* Spill / load implementation: [src/rigel/buffer.py](../../../src/rigel/buffer.py)
* Router consumer: [src/rigel/scan.py](../../../src/rigel/scan.py)

## Proposed Change

Add a private pending-spill representation to `FragmentBuffer`:

```python
@dataclass(slots=True)
class _PendingSpill:
    path: Path
    done: threading.Event
    error: BaseException | None = None
    memory_bytes: int = 0
```

`FragmentBuffer._spill_oldest()` should enqueue the chunk to a single
background writer and immediately replace the in-memory chunk with a
`_PendingSpill`. The writer serializes the chunk to Feather/IPC, sets
the completion event, and records any exception.

Consumers (`iter_chunks`, `iter_chunks_consuming`, `cleanup`, and
`release`) must wait on pending spills only when they need the data or
when resources are being cleaned up.

### Backpressure

The writer queue must be bounded. If scan output outruns the writer for
long enough to fill the queue, the callback may block. That is acceptable
and explicit: it means disk serialization is truly the bottleneck. The
goal is to remove normal spill bursts from the hot path, not to allow
unbounded RAM growth.

Initial defaults:

* one writer thread
* writer queue capacity: 2 chunks
* spill compression: keep `lz4` by default; optionally add a private
  config escape hatch for uncompressed spill benchmarking

## Implementation Steps

1. Add `_PendingSpill` and a private `_SpillWriter` helper in
   `buffer.py`. Use `threading.Thread`, `queue.Queue(maxsize=2)`, and a
   sentinel for shutdown.
2. Change `FragmentBuffer._chunks` from `deque[_FinalizedChunk | Path]`
   to `deque[_FinalizedChunk | Path | _PendingSpill]`.
3. Start the writer lazily on first spill.
4. In `_spill_oldest`, enqueue the chunk, replace it with
   `_PendingSpill`, subtract `chunk.memory_bytes` from `_memory_bytes`,
   and increment `_n_spilled` immediately.
5. In `_load_chunk` call sites, first resolve pending handles with a
   `_wait_pending_spill(...)` method that joins the event and re-raises
   writer exceptions.
6. In `cleanup()` / `release()`, stop the writer and wait for all pending
   work before deleting the temporary directory.
7. Make `summary()` report pending/in-progress spill count separately
   from completed on-disk chunks.
8. Add logging for both enqueue and finish events so profile traces can
   distinguish callback time from writer time.

## Tests

```bash
conda activate rigel
pytest tests/test_buffer.py tests/test_pipeline_smoke.py \
       tests/test_pipeline_routing.py -v
pytest tests/test_golden_output.py -v
```

Add focused tests in [tests/test_buffer.py](../../../tests/test_buffer.py):

* forced-spill iteration returns byte-for-byte identical chunks versus
  the synchronous path
* `iter_chunks_consuming()` waits for pending spill and deletes the file
* writer exception is re-raised on the next buffer operation
* `cleanup()` waits for the writer before removing the temp directory
* `release()` is idempotent with pending and completed spills

## Benchmark Plan

Run the forced-spill case, then compare to no-spill:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr04 \
  --name-prefix pr04 \
  --threads 4 8 \
  --scan-bgzf-threads 2 \
  --scan-fragments-per-chunk 1000000 \
  --scan-buffer-size 4
```

Compare against:

* `spill_s8_d2`: 200.9 s, 14 spills
* `nospill_s4_d2`: 135.7 s, 0 spills

Capture a short native/Python sample during spill bursts. The scanner
callback should no longer show sustained Arrow/LZ4 frames; those should
move to the writer thread.

## Acceptance Criteria

* Golden outputs unchanged.
* Forced-spill wall time improves materially versus 200.9 s and moves
  toward the no-spill baseline.
* Scanner callback samples are dominated by `_FinalizedChunk.from_raw` /
  bookkeeping, not Arrow IPC or LZ4.
* Peak RSS remains bounded by the configured memory budget plus at most
  `writer_queue_capacity × chunk_size` worth of chunk data.
* Writer exceptions are never swallowed.

## Risks

* Async error propagation is easy to get wrong. Treat it as a correctness
  requirement, not polish.
* A queue that is too deep hides disk bottlenecks by increasing RSS.
  Keep it small and explicit.
* Deleting the temp directory before pending writes finish can corrupt
  cleanup. `cleanup()` must always drain/stop the writer first.

## Deferred: Direct Scan-To-Score Streaming

True scan-to-score streaming would require one of these larger design
changes:

* online calibration and scoring models that are valid before the full
  scan finishes;
* provisional scoring plus later correction/rescoring;
* a two-pass BAM design where pass 1 trains/calibrates and pass 2 streams
  chunks directly to the scorer.

All three are broader than a spill-removal PR. Revisit only after the
single-pass async chunk store is measured and the remaining wall-time
profile justifies the architectural cost.
