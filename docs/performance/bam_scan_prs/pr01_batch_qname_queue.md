# PR 01: Batch Qname Queue Work

**Position in roadmap:** First. Every other PR's wall-time numbers
depend on this landing.

## Summary

Change the scanner input queue from one `QnameGroup` per push/pop to
small batches of qname groups. Removes per-read-name mutex /
condition-variable traffic and restores positive scaling with worker
threads.

## Motivation

The no-spill thread sweep shows a 50% regression from 4 → 8 workers and
another 10% from 8 → 12:

| Run | Wall time | Read names/s |
|---|---:|---:|
| s4 d2 nospill | 135.7s | 235.6k |
| s8 d2 nospill | 204.2s | 156.7k |
| s12 d2 nospill | 226.9s | 141.0k |

Live samples show worker threads spending major fractions of time in
`BoundedQueue<QnameGroup>::pop` (mutex lock, condition variable wait,
`pop_front`) and the reader spending matching time in `push`. With ~32M
pushes and ~32M pops on the VCAP workload, queue overhead is a fundamental
ceiling.

## Current Code

* Queue: [src/rigel/native/thread_queue.h](../../../src/rigel/native/thread_queue.h)
* Scanner threading: [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)

## Proposed Change

Introduce a `QnameBatch` work unit:

```cpp
struct QnameBatch {
    std::vector<QnameGroup> groups;
};
```

* Reader thread accumulates groups and pushes one batch.
* Workers pop one batch and process every group locally.
* Queue type becomes `BoundedQueue<QnameBatch>`; capacity counted in
  batches.

Defaults: `qname_batch_size = 512`, queue capacity expressed in batches
(initial guess: 4 × n_workers, tune in benchmark). Keep the batch size
private; do not expose it to users until benchmarks prove they need it.

## Implementation Steps

1. Add `QnameBatch` next to `QnameGroup` in `bam_scanner.cpp`.
2. Replace `BoundedQueue<QnameGroup>` with `BoundedQueue<QnameBatch>` in
   `BamScanner::scan`.
3. Add a small reader helper that appends the current group to an
   in-progress batch and pushes it when full.
4. Flush the final partial batch before closing the input queue.
5. Update workers to drain the batch in a tight inner loop, calling the
   existing per-group processing without changes.
6. Reuse the batch's internal `groups` vector (`groups.clear()` after
   processing) to avoid reallocation.
7. Keep `frag_id` assignment exactly where it is today.

## Tests

The existing scanner tests cover the per-group invariants. Force-batch
edge cases via two new short tests:

* `qname_batch_size = 1` produces output identical to the default size on
  a small fixture.
* A workload smaller than one batch finishes correctly (final-flush path).

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_pipeline_smoke.py tests/test_pipeline_routing.py \
       tests/test_cross_chunk.py -v
pytest tests/test_golden_output.py -v
```

## Benchmark Plan

Run the no-spill sweep so spill noise doesn't contaminate queue
measurements:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr01 \
  --name-prefix pr01 \
  --n-scan-threads 4 8 12 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

## Acceptance Criteria

* Golden outputs unchanged.
* **`s8 d2 nospill` ≤ 1.10 × `s4 d2 nospill`.** This is the headline
  criterion: positive scaling restored within 10% slack.
* `s12 d2 nospill` ≤ `s8 d2 nospill`. (Additional workers don't make the
  scanner slower.)
* Live samples show `BoundedQueue::pop` / `push` frames absent from the
  top of the worker and reader stacks.
* Peak RSS does not grow by more than `qname_batch_size × n_workers ×
  sizeof(QnameGroup)` (ballpark a few MB).

## Risks

* Large batches plus pathological qname distributions (one batch full of
  100-way multimappers) can cause load imbalance. Mitigate with a
  moderate batch size; revisit only if benchmarks show standing
  imbalance.
* Output chunk latency grows by at most one batch worth of work. With 1M
  fragments per output chunk and 512 groups per batch, this is
  negligible.

## Non-Goals

* Do not replace the queue with a lock-free implementation. Batching
  removes the lock pressure.
* Do not change resolution semantics or multimapper pairing.
