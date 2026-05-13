# PR 09: Parallel Streaming Scorer

**Position in roadmap:** Fifth. Land after PR07/PR08 reduce payload size
and scatter cost.

## Summary

Parallelize `fragment_router_scan` by scoring buffer chunks concurrently
into per-worker CSR builders, then deterministically merging the chunk
results into one `ScoredFragments` object.

## Motivation

`fragment_router_scan` is 14.72 s, the largest non-scan stage. The
current `StreamingScorer.score_chunk` path is native but single-threaded
across chunks.

Chunks are mostly independent after scan has finalized strand and
fragment-length models. This makes chunk-level parallelism the clearest
wall-time opportunity outside scan.

## Current code

- Router driver: [src/rigel/scan.py](../../../src/rigel/scan.py)
- Native scorer state: [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp)
- Buffer consuming iterator: [src/rigel/buffer.py](../../../src/rigel/buffer.py)
- Estimator unambiguous counts: [src/rigel/estimator.py](../../../src/rigel/estimator.py)

## Design constraints

- Output must be deterministic.
- Multimapper groups must not be split incorrectly. The scan buffer must
  guarantee all alignments for one molecule share a chunk, or the scorer
  must carry boundary state between adjacent chunks.
- Unambiguous count accumulation must not race on shared estimator arrays.
- Annotation output order must remain stable or be explicitly sorted.
- Peak RSS must stay bounded; parallel workers cannot all retain huge
  completed chunk outputs indefinitely.

## Proposed design

Use a two-level builder:

1. Each worker owns a local native `ChunkScorer` with local output
   vectors and local unambiguous count arrays.
2. Python submits chunks from `iter_chunks_consuming()` to a bounded
   worker pool.
3. Results are returned with a monotonically increasing `chunk_id`.
4. A deterministic merge step concatenates chunk outputs by `chunk_id`,
   fixes CSR offsets, and sums local unambiguous count arrays into the
   estimator.

Keep `NativeFragmentScorer` immutable and shared. Keep per-worker mutable
state local.

## Implementation steps

1. Split `StreamingScorer` internals into:
   - immutable `NativeFragmentScorer`,
   - per-output mutable scorer/builder state,
   - finish method for one chunk or chunk sequence.
2. Add a native function/class that scores exactly one finalized chunk
   into local arrays without touching shared estimator counts.
3. Replace direct writes to `estimator.unambig_counts` with per-worker or
   per-chunk unambiguous count buffers.
4. In [src/rigel/scan.py](../../../src/rigel/scan.py), add a parallel path
   controlled by `n_threads` or a new scoring config field. Default can
   remain serial until validated.
5. Use a bounded queue/future pool. Keep at most `2 * n_workers` chunk
   results in memory.
6. Merge chunk results in chunk order:
   - concatenate `t_indices`, payload arrays, and per-unit arrays,
   - rebuild global `offsets` by adding cumulative candidate counts,
   - concatenate deterministic and chimeric annotation arrays,
   - sum unambiguous counts.
7. Preserve the existing serial path for comparison while the new path is
   being validated. Remove it later only after confidence is high.
8. Add profile counters: scoring workers, chunks scored, merge time,
   max in-flight chunk results.

## Multimapper boundary check

Before implementing parallel scoring, add an assertion or test proving
multimapper groups are not split across chunks. If they can be split,
there are two options:

- make the buffer chunker flush only at qname group boundaries, or
- add a serial boundary reconciliation path for the first/last frag_id of
  each chunk.

Do not proceed without this invariant.

## Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_buffer.py tests/test_pipeline_routing.py tests/test_pipeline_smoke.py -v
pytest tests/test_golden_output.py -v
pytest tests/test_em_impl.py tests/test_estimator.py -v
```

Add focused tests:

- Serial and parallel scoring produce identical `ScoredFragments` arrays
  on a multi-chunk fixture.
- Parallel unambiguous counts match serial counts exactly.
- Multimapper groups at chunk boundaries are handled or rejected with a
  clear invariant failure.
- Annotation rows are deterministic.

## Benchmark plan

Run full profiles at scoring worker counts 1, 2, 4, and 8 if exposed:

```bash
conda activate rigel
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr09_parallel_scorer \
  --stages --threads 8 --memory-interval 250
```

Collect clean wall times and peak RSS. Use cProfile only to attribute
merge overhead.

## Acceptance criteria

- Golden outputs unchanged.
- Parallel path produces array-identical `ScoredFragments` to serial on
  focused tests.
- `fragment_router_scan` improves by at least 2x on VCAP at 8 threads,
  with a stretch target of 3-5x.
- Peak RSS does not grow by more than the configured in-flight chunk
  result bound.
- Merge time is reported separately and is not more than 25% of the old
  router stage.

## Risks

- Shared estimator count races. Avoid by using local count buffers.
- Memory growth from too many completed chunk outputs. Use bounded
  in-flight work and ordered merge.
- Multimapper chunk-boundary semantics. Prove or fix before parallelism.

## Non-goals

- Do not change scoring formulas.
- Do not combine with float32 payload conversion.
- Do not redesign calibration or scan-to-score streaming in this PR.