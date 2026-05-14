# PR 09: Parallel Streaming Scorer

**Position in roadmap:** Fifth. Land after PR07 (smaller payloads) and
PR08 (cheaper scatter) so the parallel scorer doesn't immediately
saturate memory bandwidth.

## Summary

Score buffer chunks concurrently into per-worker CSR builders, then
merge in chunk order into one `ScoredFragments`. Eliminate the
estimator-shared mutable `unambig_counts` write inside the scoring
kernel.

## Motivation

`fragment_router_scan` is 14.72 s, the largest non-scan stage.
Today the loop in
[src/rigel/scan.py](../../../src/rigel/scan.py) (`FragmentRouter`) is
single-threaded across chunks: it drains
`buffer.iter_chunks_consuming()` and calls a single
`StreamingScorer.score_chunk` instance whose vectors grow monotonically.

Chunks are independent in two specific senses that make
chunk-parallelism legal:

* Strand and FL models are frozen before scoring begins.
* The scanner already ends chunks at qname-group boundaries (an
  invariant noted in `scoring.cpp` ~line 1003 — "cross-chunk references
  not possible"). Every alignment for one molecule lives in one chunk.

What is **not** independent today is `unambig_counts`. The native
scorer writes deterministic-fragment counts directly into a single
estimator-owned `f64_2d_mut` array passed through
[src/rigel/scan.py](../../../src/rigel/scan.py) ~line 129. Any parallel
design must privatize that array per worker and merge.

## Current code

* Driver: [src/rigel/scan.py](../../../src/rigel/scan.py)
  (`FragmentRouter`, lines ~30–150).
* Native scorer: [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp)
  (`StreamingScorer`, `score_chunk` ~line 1128, `finish` ~line 1157).
* Buffer iterator: [src/rigel/buffer.py](../../../src/rigel/buffer.py)
  (`iter_chunks_consuming`).
* Estimator counts: [src/rigel/estimator.py](../../../src/rigel/estimator.py)
  (`unambig_counts` ~line 192).

## Design

Three separable pieces:

### A. Privatize unambig counts

Change the native scoring entrypoint so it writes into a
**locally-owned** `unambig_counts` accumulator (one per chunk or per
worker). The driver sums these into the estimator after the scoring
phase finishes. This is a prerequisite for any concurrent design and
also makes the existing serial path easier to reason about.

### B. Per-chunk scorer

Refactor `StreamingScorer` into:

* `NativeFragmentScorer` (immutable, shared) — index, models,
  scoring config. Holds nothing mutable.
* `ChunkScorer` (per-call) — owns the local output vectors and local
  `unambig_counts`. Constructed cheaply, destroyed when the chunk's
  output is consumed.
* `score_chunk(chunk_arrays) -> ChunkResult` returns a Python tuple
  of numpy arrays + a chunk_id. No global state mutation.

### C. Concurrent driver and ordered merge

Python driver:

```python
in_q  = queue.Queue(maxsize=2 * n_workers)   # chunks awaiting work
out_q = queue.Queue(maxsize=2 * n_workers)   # completed ChunkResults

# producer: drain buffer.iter_chunks_consuming(), assign monotone chunk_id, push
# workers: pop chunk, call native scorer, push (chunk_id, result)
# merger: pop in chunk_id order (heap or expected_id counter), append to global builders
```

Use `concurrent.futures.ThreadPoolExecutor`. The native call releases
the GIL, so threads scale.

The global builders are append-only Python lists of per-chunk arrays;
final assembly concatenates and rebuilds CSR offsets in a single pass
once all chunks are merged. Concatenation is a reduction step, not a
per-chunk cost.

## Multimapper boundary check

Add an explicit invariant test before parallelization is enabled:

```python
def test_chunks_dont_split_qname_groups():
    # Construct a buffer where chunk size forces a multimapper group near
    # the boundary, finalize, iterate, assert every frag_id appears in
    # exactly one chunk.
```

If it ever fails, the buffer chunker (not the scorer) is wrong; fix
there.

## Implementation steps

1. **Privatize `unambig_counts` (independent commit).** Change the
   native scorer to fill a local accumulator and return it from
   `finish()`. Driver sums into the estimator. Keep the rest of the
   scorer single-threaded for this commit. Validate against goldens.
2. Split `StreamingScorer` per the A/B refactor above. The per-call
   `ChunkScorer` should be cheap to construct (no model rebuild); all
   immutable state stays in `NativeFragmentScorer`.
3. Bind `score_chunk` as a free function on `NativeFragmentScorer`
   that returns its own ChunkResult tuple.
4. Add the parallel driver in `scan.py` behind a config flag
   (`scoring.n_threads`, default 1 for the first PR landing).
5. Implement ordered merge using a min-heap keyed by chunk_id.
6. Concatenate per-chunk arrays in chunk order. Rebuild global CSR
   offsets by scanning concatenated `count_cols`.
7. Sum per-chunk `unambig_counts` into the estimator at end.
8. Bound in-flight memory: cap `out_q` and total queued
   ChunkResult bytes via PR05's `scoring_csr.candidate_bytes`
   accounting (or a simpler cap on number of in-flight chunks).
9. Add profile counters: scoring workers, chunks processed, merge
   wall, peak in-flight chunk-result bytes.

## Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_buffer.py tests/test_pipeline_routing.py \
       tests/test_pipeline_smoke.py tests/test_cross_chunk.py -v
pytest tests/test_em_impl.py tests/test_estimator.py -v
pytest tests/test_golden_output.py -v
```

Add focused tests:

* Privatized `unambig_counts` (commit 1) gives bit-identical estimator
  totals to the previous shared-write design.
* `scoring.n_threads = 1` and `n_threads = 4` produce array-identical
  `ScoredFragments` and identical `unambig_counts`.
* Multimapper-boundary invariant test (above).
* Annotation arrays preserve the expected order (sorted by chunk_id,
  then by within-chunk index).

## Benchmark plan

```bash
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr09_parallel_scorer \
  --stages --threads 8 --memory-interval 250
```

Sweep `scoring.n_threads` ∈ {1, 2, 4, 8}. Compare wall, peak RSS, and
the merge wall reported by the new profile counter.

## Acceptance criteria

* Golden outputs unchanged.
* Parallel and serial paths produce array-identical `ScoredFragments`
  and `unambig_counts` on focused tests.
* `fragment_router_scan` improves by ≥ 2× at `scoring.n_threads = 8`
  on VCAP. Stretch target: 3×.
* Peak RSS does not exceed
  `pre_PR09_peak + max_inflight_chunks × per_chunk_bytes`.
* Merge wall reported separately and ≤ 25% of pre-PR09
  `fragment_router_scan` time.
* Multimapper-boundary invariant test added and passing.

## Risks

* The `unambig_counts` privatization is the largest correctness change
  in the PR. Land it in commit 1 so any regression is bisectable.
* Memory growth from queued ChunkResults. Bound the queue and instrument
  the bound.
* Determinism. Always merge in chunk_id order. Never let workers append
  directly to global builders.

## Non-goals

* No scoring-formula changes.
* No FragmentBuffer-layer changes (PR10 owns those).
* No fusion of scoring with partitioning. Scoring still feeds an
  intermediate `ScoredFragments`; fusing further is a separate
  conversation.
