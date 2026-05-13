# PR 03: Rewrite Resolver Set Logic With Reusable Scratch Vectors

## Summary

Replace hot per-fragment `unordered_set`, `unordered_map`, and nested vector
allocations inside `FragmentResolver::_resolve_core` with reusable sorted-vector
scratch buffers.

## Motivation

Native samples repeatedly showed allocator-heavy stacks in `_resolve_core`:

- `std::__hash_table<int>::__emplace_unique_key_args`
- `operator new`
- `free_tiny`
- `_platform_memset`
- `std::vector<int>::insert`
- `merge_sets`
- `detect_chimera`

This work occurs for millions of fragments. The current implementation creates
new containers for each call:

- `std::unordered_set<int32_t> block_exon_t`
- `std::unordered_set<int32_t> block_transcript_t`
- `std::vector<std::vector<int32_t>> exon_t_sets`
- `std::vector<std::vector<int32_t>> transcript_t_sets`
- `std::unordered_set<int32_t> nrna_set`
- `std::unordered_set<int32_t> all_overlap_t`
- additional sets/maps in `merge_sets` and `detect_chimera`

The transcript candidate sets are already sorted before most downstream use.
Using append/sort/unique on reusable vectors should be much cheaper than hashing
small integers into fresh heap nodes.

## Current Code

- Resolver core: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)
- Set helpers: [src/rigel/native/constants.h](../../../src/rigel/native/constants.h)

## Proposed Scope

This PR should preserve the exact candidate sets and merge criteria while
changing only the internal representation.

Add reusable buffers to `ResolverScratch`, for example:

```cpp
std::vector<std::vector<int32_t>> exon_t_sets;
std::vector<std::vector<int32_t>> transcript_t_sets;
std::vector<std::vector<int32_t>> sj_t_sets;
std::vector<int32_t> tmp_candidates;
std::vector<int32_t> tmp_union;
std::vector<int32_t> tmp_intersection;
std::vector<int32_t> tmp_overlap;
std::vector<int32_t> tmp_nrna;
```

Introduce helper operations:

- `sort_unique(vector<int32_t>&)`
- `append_unique_sorted(dst, value)` only where incremental insertion is truly
  cheaper than append/sort/unique
- `intersect_sorted_into(a, b, out)`
- `union_sorted_into(parts, out)`
- `contains_sorted(vec, value)`

## Implementation Steps

1. Extend `ResolverScratch` with reusable set buffers and a `reset_sets(n_exons)`
   method that clears vectors without freeing capacity.
2. Replace per-block `unordered_set` construction with appending to per-block
   vectors followed by `sort_unique`.
3. Replace `nrna_set` with a scratch vector plus `sort_unique`.
4. Replace `all_overlap_t` with a sorted union scratch vector and binary search.
5. Add a scratch-aware `merge_sets_into(...)` that writes into caller-provided
   output vectors and avoids constructing temporary unions via `unordered_set`.
6. Add a scratch-aware `detect_chimera` or a low-allocation variant for the
   common small-block case.
7. Keep the old helpers temporarily if tests or Python-visible APIs still use
   them, but move scanner hot paths to the scratch-aware versions.

## Tests

Add resolver equivalence tests that compare old and new behavior on synthetic
block patterns before deleting old code, or keep a test-only reference helper.
Cases should include:

- single-exon unique mapper
- multi-exon annotated splice
- unannotated splice
- synthetic nRNA parent injection
- intergenic empty `t_inds`
- cis and trans chimera detection
- opposite-strand overlaps

Run:

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_resolution.py tests/test_splice.py tests/test_splice_blacklist.py -v
pytest tests/test_native_gdna_eligibility.py tests/test_pipeline_integration_v6.py -v
```

## Benchmark Plan

Use both scan profiler and a future resolver microbenchmark if one is added.
The scan profiler should be enough for PR acceptance:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr03 \
  --name-prefix pr03 \
  --n-scan-threads 4 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

## Acceptance Criteria

- Golden outputs are unchanged.
- Native samples show less allocator time under `_resolve_core`.
- `nospill_s4_d2` wall time improves versus 135.7s.
- Candidate counts, splice type counts, intergenic counts, and chimera counts are
  identical to baseline on the VCAP scan.

## Risks

- Sorted-vector semantics must exactly match set semantics, including duplicate
  removal and deterministic order.
- nRNA parent injection currently inserts into every block set while preserving
  sorted uniqueness. The replacement must preserve that behavior.
- Chimera detection is biologically meaningful; do not simplify it without
  targeted tests.

## Non-Goals

- Do not alter cgranges interval indexing in this PR.
- Do not change transcript-centric resolution semantics.
