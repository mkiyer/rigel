# PR 02: Resolver Scratch Vectors

## Summary

Replace per-fragment `unordered_set<int32_t>` and allocation-returning
set helpers inside `_resolve_core` with thread-local scratch
`std::vector<int32_t>` buffers. Use STL set algorithms (`std::sort`,
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

* Resolver core and `ResolverScratch`: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)
* Existing allocation-returning set helpers: [src/rigel/native/constants.h](../../../src/rigel/native/constants.h)

Hot constructions to replace:

* `std::unordered_set<int32_t>` for per-exon transcript-ID sets, then
  merged across exons.
* `std::unordered_set<int32_t>` for nRNA parent collection and
  `all_overlap_t`.
* `merge_sets(...)` currently returns vectors by value and uses an
  `unordered_set` for the union branch.
* `detect_chimera(...)` allocates union-find bookkeeping plus
  `unordered_map` / `unordered_set` temporaries.
* Repeated allocation of small temporary vectors for transcript-ID
  intersection / union.

## Proposed Change

Extend the existing `ResolverScratch` struct in `resolve_context.h`. It
is already owned by the worker thread; this PR adds set-operation
buffers to it and resets those buffers per fragment:

```cpp
struct ResolverScratch {
    std::vector<std::vector<int32_t>> exon_t_sets;       // per block, sorted-unique
    std::vector<std::vector<int32_t>> transcript_t_sets; // per block, sorted-unique
    std::vector<std::vector<int32_t>> sj_t_sets;         // per intron, sorted-unique
    std::vector<int32_t> tmp_a;
    std::vector<int32_t> tmp_b;
    std::vector<int32_t> tmp_out;
    std::vector<int32_t> tmp_union;
    std::vector<int32_t> all_overlap_t;
    std::vector<int32_t> nrna_t;

    void reset_per_fragment() noexcept {
        for (auto& v : exon_t_sets) v.clear();
        for (auto& v : transcript_t_sets) v.clear();
        for (auto& v : sj_t_sets) v.clear();
        tmp_a.clear(); tmp_b.clear(); tmp_out.clear();
        tmp_union.clear(); all_overlap_t.clear(); nrna_t.clear();
    }
};
```

    `ResolverScratch&` is already passed through `_resolve_core`; keep that
    ownership model.

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

Do **not** add a broad rigel-specific wrapper layer. STL is sufficient.
However, do not leave the existing hot helpers untouched: either rewrite
`merge_sets` / `detect_chimera` in place to avoid hash containers, or
bypass them from `_resolve_core` with local STL operations. The end
state should have one implementation of each algorithm, not old and new
variants that drift.

## Implementation Steps

1. Extend the existing `ResolverScratch` in `resolve_context.h` with the
  reusable set-operation buffers above.
2. Add `sort_unique` near the resolver hot path as the only new helper.
3. Walk `_resolve_core` and helpers; replace each `unordered_set<int32_t>`
   construction with `scratch.tmp_X.assign(...) + sort_unique`.
4. Replace each set-merge / intersect with `std::set_union` /
   `std::set_intersection` writing to a fresh scratch vector.
5. Replace `ref_set`, `nrna_set`, and `all_overlap_t` with sorted-vector
   scratch buffers.
6. Rewrite or bypass the hot `merge_sets` union branch so it no longer
   constructs an `unordered_set`.
7. Rewrite `detect_chimera` in place with vector-backed components and
   sorted unique strand collection; do not create a second chimera
   algorithm.
8. Leave index-load maps alone (`SJMap`, blacklist maps, ref lookup).
   They are not per-fragment allocations.
9. Call `reset_per_fragment()` at the top of `_resolve_core` (the per-
   fragment `clear()` is implicit in `assign`, but doing it explicitly
   keeps lifetime obvious).

## Tests

This is a pure representation change; existing tests are the contract:

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_resolution.py tests/test_pipeline_routing.py \
       tests/test_pipeline_smoke.py tests/test_orient_routing.py \
  tests/test_layout_iter.py tests/test_splice_blacklist.py -v
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
* `git grep` shows no local `std::unordered_set` construction inside
  `_resolve_core`, `merge_sets`, or `detect_chimera`. Index-level maps
  and other cold-path maps may remain.
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
