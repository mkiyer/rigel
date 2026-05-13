# PR (deferred): Fragment Assembly Allocator Cleanup

**Status:** Deferred. Revisit only if `build_fragment` and friends are
still hot in profiles after PR 02 (Resolver Scratch Vectors) lands.

## Why deferred

The original PR-07 forked `build_fragment` into a "fast path" (typical
unique mapper, hand-rolled allocation) and a "slow path" (the existing
multimapper logic). That doubles the test surface for one function and
creates drift risk: in two years, when someone adds a feature, they fix
it on the fast path, miss the slow path, and we ship a subtle
inconsistency.

The underlying problem (small `unordered_map` and `std::set` allocations
inside `build_fragment` and adjacent helpers) is the same problem PR 02
solves in `_resolve_core`. The same fix — thread-local scratch vectors,
sorted-unique invariants, STL set algorithms — applies here too. There
is no need for two paths.

The right sequence is:

1. Land PR 02 (resolver scratch vectors), validate it.
2. Re-profile. If `build_fragment` is no longer in the top frames,
   there is no work to do.
3. If it is, apply the same scratch-vector pattern *in place* in
   `build_fragment`. One path, one set of tests, no fast/slow guard.

## Preserved analysis from the original PR-07

The functions of interest are in
[src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)
(`build_fragment`, `group_records_by_hit`, `pair_multimapper_reads`) and the data structures
in [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h).

Hot operations observed in the original profile:

* `build_fragment`: `std::unordered_map<int64_t, vector<pair<int32_t,
  int32_t>>>` for `(ref_id, strand)` exon groups.
* `build_fragment`: `std::set<tuple<int32_t, int32_t, int32_t,
  int32_t>>` for intron de-duplication.
* `group_records_by_hit`: `std::unordered_map<int32_t, ...>` for HI-tag
  grouping.
* `pair_multimapper_reads`: `std::unordered_set<int>` for paired-location
  bookkeeping plus repeated `build_fragment` calls for secondary
  locations.
* Transient `std::vector` allocations for per-fragment/per-hit metadata.

If these persist after PR 02, the in-place rewrite is:

* Replace tiny maps with sorted parallel vectors or a small vector of
  pairs plus `std::sort` / `std::lower_bound`.
* Replace `std::set<tuple<...>>` with `std::vector<IntronBlock> + sort +
  unique`.
* Move per-thread scratch into a `BuildFragmentScratch` struct owned
  by the worker, mirroring `ResolverScratch`.

## Acceptance Criteria (if revived)

* All scanner / pipeline tests pass.
* Golden outputs unchanged.
* `build_fragment` self-time drops by ≥ 30% on the post-PR-02 baseline.
* No new code paths; the change is representational only.

## Non-Goals

* No fast/slow path split. One implementation, always.
* No public API changes.
