# PR 02: Resolver Scratch Vectors

## Summary

Replace per-fragment `unordered_set<int32_t>` and `unordered_map<int32_t,
…>` allocations inside `_resolve_core` with thread-local scratch
`std::vector<int32_t>` buffers, and use STL set algorithms (`std::sort`,
`std::unique`, `std::set_intersection`, `std::set_union`,
`std::binary_search`) instead of hash containers.

## Motivation

Live samples on the VCAP workload show malloc/free traffic from
`unordered_set` and `unordered_map` construction inside the resolver
dominating CPU in the per-fragment hot path. Each ambiguous read cluster
constructs and destroys multiple hash containers; for 32M fragments this
is the biggest single source of allocator pressure in the scan stage.

Sorted-vector set operations are uniformly faster than `unordered_set`
for the small-K, integer-key case typical of resolver work, *and* they
free us from amortised-rehash spikes.

## Current Code

* Resolver core: [src/rigel/native/resolve.cpp](../../../src/rigel/native/resolve.cpp)
* Resolve context types: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)

Hot constructions to replace:

* `std::unordered_set<int32_t>` for per-exon transcript-ID sets, then
  merged across exons.
* `std::unordered_map<int32_t, int32_t>` for per-fragment FL accumulation
  (covered separately in PR 03 for the FL case; this PR covers the
  resolver-side maps).
* Repeated allocation of small temporary vectors for transcript-ID
  intersection / union.

## Proposed Change

Add a `ResolverScratch` struct, owned by the worker thread (one per
thread), reset per fragment:

```cpp
struct ResolverScratch {
    std::vector<int32_t> exon_t_set;       // sorted-unique
    std::vector<int32_t> transcript_t_set; // sorted-unique
    std::vector<int32_t> tmp_a;
    std::vector<int32_t> tmp_b;
    std::vector<int32_t> tmp_out;

    void reset_per_fragment() noexcept {
        exon_t_set.clear();
        transcript_t_set.clear();
    }
};
```

Pass `ResolverScratch&` through `_resolve_core` and helpers.

Replace hash-set constructions:

```cpp
// Before: std::unordered_set<int32_t> ts(begin, end);
// After:
scratch.tmp_a.assign(begin, end);
sort_unique(scratch.tmp_a);
```

Replace merge / intersect operations:

```cpp
// merge two sorted-unique vectors into out
std::set_union(a.begin(), a.end(), b.begin(), b.end(),
               std::back_inserter(out));
// intersect
std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                      std::back_inserter(out));
// membership
bool present = std::binary_search(a.begin(), a.end(), x);
```

Only one helper is justified, because `std::sort` + `std::unique` +
`erase` is a 3-line idiom we'll write a dozen times:

```cpp
inline void sort_unique(std::vector<int32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}
```

Do **not** add `intersect_sorted_into`, `union_sorted_into`,
`append_unique_sorted`, `contains_sorted`, `merge_sets_into`,
`detect_chimera` scratch-aware variants, or any other rigel-specific
wrappers. STL is sufficient and more idiomatic.

## Implementation Steps

1. Add `ResolverScratch` and `sort_unique` to `resolve.cpp` (private,
   anonymous namespace).
2. Hold one `ResolverScratch` per worker thread (the existing per-thread
   resolver state already exists; add this as a member).
3. Walk `_resolve_core` and helpers; replace each `unordered_set<int32_t>`
   construction with `scratch.tmp_X.assign(...) + sort_unique`.
4. Replace each set-merge / intersect with `std::set_union` /
   `std::set_intersection` writing to a fresh scratch vector.
5. Replace `unordered_map<int32_t, T>` lookups inside the resolver hot
   path with parallel sorted-vector + lower_bound lookups, *only* where
   profiling justifies it. Cold paths can stay.
6. Call `reset_per_fragment()` at the top of `_resolve_core` (the per-
   fragment `clear()` is implicit in `assign`, but doing it explicitly
   keeps lifetime obvious).

## Tests

This is a pure representation change; existing tests are the contract:

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_resolution.py tests/test_pipeline_routing.py \
       tests/test_pipeline_smoke.py tests/test_orient_routing.py \
       tests/test_layout_iter.py -v
pytest tests/test_golden_output.py -v
```

Add one new test that constructs an ambiguous-multimapper fixture (≥ 8
candidate transcripts overlapping ≥ 4 exons) and checks the resolver
output equals a snapshot. This is regression insurance for the merge
logic.

## Benchmark Plan

Run after PR 01 has landed (otherwise queue noise dominates):

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr02 \
  --name-prefix pr02 \
  --n-scan-threads 8 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

Compare against `pr01` baseline at the same `s8 d2` setting.

## Acceptance Criteria

* All resolver / routing tests pass; golden outputs unchanged.
* `_resolve_core` self-time drops by ≥ 30% on the s8 nospill run.
* Live sampling no longer shows `unordered_set::insert` /
  `_M_rehash_aux` frames in resolver stacks.
* No new heap allocations in steady-state worker loops (verified by
  spot-check with `MallocStackLogging` or instrumented allocator if
  desired; not gating).

## Risks

* Sorted-vector set operations require sorted-unique invariants. Audit
  every assignment to the scratch vectors. Add `assert(std::is_sorted
  (...))` debug-only checks at key boundaries.
* `std::set_intersection` requires both inputs sorted. The `sort_unique`
  helper enforces this; document the precondition next to each call.

## Non-Goals

* Don't touch the resolver's actual algorithm.
* Don't move scratch ownership to the BoundedQueue work item; keep it on
  the worker thread.
