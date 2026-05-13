# PR 02: Batch Qname Queue Work

## Summary

Change the scanner input queue from one `QnameGroup` per push/pop to batches of
qname groups. This reduces mutex, condition-variable, move, and destructor
overhead in the reader-to-worker pipeline.

## Motivation

The no-spill sweep showed severe negative scaling:

| Run | Wall time | Read names/s |
|---|---:|---:|
| s4 d2 no-spill | 135.7s | 235.6k |
| s8 d2 no-spill | 204.2s | 156.7k |
| s12 d2 no-spill | 226.9s | 141.0k |

Live samples showed worker threads spending large fractions of time in
`BoundedQueue<QnameGroup>::pop`, including `std::mutex::lock`,
`std::condition_variable::wait`, `notify_one`, and `pop_front`. The reader also
spent heavy time in `BoundedQueue<QnameGroup>::push` when worker counts changed.

The current design performs one queue operation per read name, roughly 32M queue
items for the VCAP workload. That is too fine-grained.

## Current Code

- Queue implementation: [src/rigel/native/thread_queue.h](../../../src/rigel/native/thread_queue.h)
- Scanner threading: [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)

## Proposed Scope

Introduce a `QnameBatch` work unit:

```cpp
struct QnameBatch {
    std::vector<QnameGroup> groups;
};
```

The reader thread accumulates groups into a batch and pushes the batch to the
input queue. Workers pop one batch and process its groups in a local loop.

Initial defaults:

- `qname_batch_size`: 512 or 1024 groups
- input queue capacity: measured in batches, not groups
- output queue unchanged for this PR

Expose the batch size as an internal scanner argument only if needed for
benchmarking. Otherwise keep it private until there is evidence users need to
tune it.

## Implementation Steps

1. Add `QnameBatch` next to `QnameGroup` in `bam_scanner.cpp`.
2. Change `BoundedQueue<QnameGroup>` to `BoundedQueue<QnameBatch>` in
   `BamScanner::scan`.
3. Add a reader helper that appends the current group to the current batch and
   pushes the batch when full.
4. Ensure the final partial batch is flushed before closing the input queue.
5. Update workers to process every group in the popped batch.
6. Keep `frag_id` assignment exactly where it is today so output identity is
   unchanged.
7. Consider `groups.clear()` reuse to reduce batch vector reallocation.

## Tests

Existing scanner behavior tests should catch most regressions. Add a focused
test if needed with a small name-sorted BAM where batch size is forced to 1 and
then to a larger value; outputs should match.

Run:

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_pipeline_smoke.py tests/test_pipeline_routing.py tests/test_cross_chunk.py -v
pytest tests/test_golden_output.py -v
```

## Benchmark Plan

Use the no-spill sweep to isolate queue behavior:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr02 \
  --name-prefix pr02 \
  --n-scan-threads 4 8 12 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

## Acceptance Criteria

- Golden outputs are unchanged.
- `s8 d2 no-spill` is no longer dramatically slower than `s4 d2 no-spill`.
- Live samples show far fewer `BoundedQueue<QnameGroup>::pop` and
  `BoundedQueue<QnameGroup>::push` frames.
- Queue batching does not increase peak RSS by more than the expected batch
  buffer size.

## Risks

- Larger batches can create load imbalance if a batch contains many expensive
  multimapper groups. Start with moderate batch sizes and benchmark.
- Larger batches can increase latency before output chunks fill. This is likely
  negligible compared with 1M-fragment output chunks.

## Non-Goals

- Do not replace the queue implementation with a lock-free queue in this PR.
  Batching should remove most of the lock pressure with far less risk.
- Do not alter fragment resolution semantics or multimapper pairing.
