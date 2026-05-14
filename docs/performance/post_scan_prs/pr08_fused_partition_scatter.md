# PR 08: Fused Partition Scatter

**Position in roadmap:** Fourth. Pairs cleanly with PR07; if PR07 lands
first the fused scatter writes float32 payloads.

## Summary

Replace the 11 sequential native scatter calls in `partition_and_free`
with one fused per-locus pass implemented on top of the existing
templated `scatter_*_impl<T>` helpers. Walk each locus's unit list
**once**, allocating and filling all destination arrays for that locus
together.

## Motivation

`partition_and_free` is 6.15 s, ~8.4% of the post-scan pipeline. The
current implementation in
[src/rigel/locus_partition.py](../../../src/rigel/locus_partition.py)
calls `build_partition_offsets` once and then issues 11 scatter calls,
each of which re-walks the same locus / unit structure. The work is
memory-bandwidth bound; one fused traversal is the obvious fix.

## Current code

* Python orchestration:
  [src/rigel/locus_partition.py](../../../src/rigel/locus_partition.py)
  (`partition_and_free`, lines ~27–100).
* Templated scatter:
  [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp)
  (`scatter_candidates_impl<T>` ~line 1371,
  `scatter_units_impl<T>` ~line 1436;
  bindings `scatter_candidates_{f64,i32,u8}`,
  `scatter_units_{f64,i32,u8,i64}` ~lines 2342–2373).
* Per-locus container: [src/rigel/scored_fragments.py](../../../src/rigel/scored_fragments.py)
  (`LocusPartition`).

## Proposed change

Add **one** native entrypoint, `partition_scatter_fused`, that consumes
all source arrays plus the per-locus `units` index lists and returns
ready-to-construct `LocusPartition` arrays. Reuse the existing
templates internally — do not duplicate the per-array copy logic.

### Native signature (shape-only)

```cpp
nb::list partition_scatter_fused(
    // global CSR (offsets-aware payloads)
    int64_1d offsets,
    int32_1d t_indices,
    /* candidate payload */    nb::handle log_liks,         // f32 or f64
    u8_1d count_cols,
    /* candidate payload */    nb::handle coverage_weights, // f32 or f64
    // global per-unit fields
    /* unit payload */         nb::handle gdna_log_liks,    // f32 or f64
    int32_1d locus_t_indices,
    u8_1d locus_count_cols,
    i8_1d is_spliced,
    int64_1d frag_ids,
    i8_1d frag_class,
    u8_1d splice_type,
    // partition spec
    nb::list locus_units   // list[int64_1d]; length == n_loci
);
```

* Returns a Python `list` of length `n_loci`; each entry is a tuple of
  numpy arrays in the order that `LocusPartition.__init__` expects.
* Payload dtypes are detected per-array (float32 vs float64) so this
  PR does not depend on PR07 having landed but benefits when it does.

### Implementation pattern

```cpp
// Under GIL: cache pointer + dtype + length for every input.
// Release GIL.
//
// For locus k in 0 .. n_loci - 1:
//   const int64_t* u = locus_units[k];
//   size_t n_units = sizes[k];
//   // candidate offsets within this locus, computed in-place
//   compute_local_offsets(u, n_units, offsets, /*out*/ local_off);
//   // allocate destination arrays
//   alloc_locus_arrays(...);
//   // single pass, copying both candidate ranges and unit fields
//   for (size_t i = 0; i < n_units; ++i) {
//     int64_t unit = u[i];
//     copy_candidate_range(unit, dst_cand_arrays);
//     copy_unit_fields    (unit, dst_unit_arrays);
//   }
//
// Re-acquire GIL only to wrap arrays in capsules and append to result.
```

The per-array `copy_candidate_range` / `copy_unit_fields` helpers are
the bodies of the existing `scatter_*_impl<T>` templates, lifted to
take an explicit destination index. Do not introduce new templates;
share code.

### Caller change

`partition_and_free` collapses to one call:

```python
parts = native.partition_scatter_fused(
    offsets, t_indices, log_liks, count_cols, coverage_weights,
    gdna_log_liks, locus_t_indices, locus_count_cols,
    is_spliced, frag_ids, frag_class, splice_type,
    locus_units,
)
em_data.release_global_arrays()  # existing free path
return {locus_id: LocusPartition(*parts[k]) for k, locus_id in enumerate(locus_ids)}
```

Keep `release_global_arrays` exactly as it is today.

## Implementation steps

1. Add `partition_scatter_fused` in `em_solver.cpp` next to the
   existing scatter helpers; refactor `scatter_*_impl<T>` so its inner
   copy body is callable from both the legacy single-array binding
   *and* the fused entrypoint.
2. Bind it via nanobind. Detect payload dtype with `dtype_kind` /
   `dtype_size` checks; dispatch to f32 vs f64 branches.
3. Release the GIL for all allocation and copy work.
4. Wire `partition_and_free` to call the fused entrypoint.
5. Keep the legacy single-array bindings (`scatter_candidates_f64`
   etc.) in place for now; mark them `// retained for tests` in the
   nanobind module so they can be deleted in a follow-up after every
   call site is on the fused path.
6. Optional second commit (same PR, behind a config flag, off by
   default): parallelize the per-locus loop with OpenMP. Skip if
   single-thread fused already meets the acceptance criterion.

## Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_locus_partition.py tests/test_partition_units.py -v
pytest tests/test_em_impl.py tests/test_estimator.py -v
pytest tests/test_pipeline_smoke.py tests/test_golden_output.py -v
```

Add focused tests:

* Construct a small `ScoredFragments` with several loci of varying
  shape (empty, one-unit, variable candidate widths). Compare the
  fused output array-by-array against the legacy 11-call scatter.
* Verify dtype detection: pass float32 and float64 payloads; assert
  output dtypes match input dtypes.
* Repeat under a tight GC loop to catch capsule ownership bugs.

## Benchmark plan

```bash
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr08_fused_scatter \
  --stages --threads 8 --memory-interval 250
```

## Acceptance criteria

* Per-locus arrays are byte-for-byte identical to the legacy scatter on
  focused tests; golden outputs unchanged.
* `partition_and_free` improves by ≥ 25% on the VCAP profile.
* No peak-RSS regression > 250 MB.
* No EM regression > 5%.
* All `LocusPartition` arrays remain contiguous and have expected
  dtypes (verified via PR05's profile JSON).

## Risks

* Many-array native bindings invite capsule / lifetime bugs. Lift the
  copy bodies as functions; do not embed them in lambdas.
* Parallelizing per locus before single-thread is correct invites two
  bug categories at once. Land single-thread first, parallel as a
  follow-on commit.
* Removing the legacy bindings without an audit may break tests.
  Retain them in this PR; delete in a separate PR.

## Non-goals

* No EM extraction-semantics changes.
* No locus-construction changes.
* No payload-dtype changes (PR07 owns that).
