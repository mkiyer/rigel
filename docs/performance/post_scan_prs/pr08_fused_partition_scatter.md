# PR 08: Fused Partition Scatter

**Position in roadmap:** Fourth. Best after PR07 so the fused scatter can
use final payload dtypes.

## Summary

Replace the current array-by-array partition scatter with a fused native
scatter that walks each locus's unit list once and fills all destination
arrays for that locus together. Optionally parallelize across loci after
the fused single-thread implementation is correct.

## Motivation

`partition_and_free` is 6.15 s, about 8.4% of the current full profile.
The current implementation performs one offsets pass, four candidate
array scatters, and seven unit array scatters. Each pass repeatedly walks
the same locus/unit structure.

This is memory-bandwidth work. Fusing the pass makes the code simpler at
the pipeline level and reduces repeated traversal.

## Current code

- Python orchestration: [src/rigel/locus_partition.py](../../../src/rigel/locus_partition.py)
- Native scatter helpers: [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp)
- Per-locus container: [src/rigel/scored_fragments.py](../../../src/rigel/scored_fragments.py)

## Proposed native API

Add one native function, name subject to taste:

```python
partition_scatter_fused(
    offsets,
    t_indices,
    log_liks,
    count_cols,
    coverage_weights,
    is_spliced,
    gdna_log_liks,
    locus_t_indices,
    locus_count_cols,
    frag_ids,
    frag_class,
    splice_type,
    locus_units,
    n_loci,
) -> list[tuple]
```

Each returned tuple should contain the arrays needed to construct one
`LocusPartition`.

The implementation should allocate all destination arrays for one locus,
then copy unit and candidate data in one loop over that locus's units.

## Implementation steps

1. Keep the existing scatter functions in place initially.
2. Add a fused native helper in
   [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp).
3. Under the GIL, extract raw pointers and lengths for each `locus_units`
   array into a small metadata vector.
4. Release the GIL for allocation and copying.
5. For each locus:
   - compute local partition offsets,
   - allocate candidate arrays,
   - allocate unit arrays,
   - walk `units[k]` once,
   - copy candidate segments with `memcpy`,
   - gather unit fields by direct indexing.
6. Wrap the result arrays with capsules and return a Python list.
7. Replace `partition_and_free(...)` internals with one fused call while
   preserving the public return type: `dict[int, LocusPartition]`.
8. Null out global `em_data` arrays after the fused call, as today.
9. Once correct, add OpenMP or `std::thread` parallelism across loci if
   the single-thread fused pass leaves enough wall time. Keep this as a
   second commit inside the same PR only if the first commit is already
   measured.
10. Remove old scatter functions only after all call sites and tests are
    updated. If external tests import them, deprecate first rather than
    delete.

## Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_locus_partition.py tests/test_partition_units.py -v
pytest tests/test_em_impl.py tests/test_estimator.py -v
pytest tests/test_pipeline_smoke.py tests/test_golden_output.py -v
```

Add focused tests:

- Construct a small `ScoredFragments` with several loci and compare the
  old scatter output to fused output field-by-field.
- Include empty loci, one-unit loci, and a locus with variable candidate
  widths.
- Include bool/int8 `is_spliced` input to verify dtype handling.

## Benchmark plan

Use the full profile and PR05 metrics:

```bash
conda activate rigel
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr08_fused_scatter \
  --stages --threads 8 --memory-interval 250
```

Compare `partition_and_free`, peak RSS, and `locus_em`. The EM stage
should not regress; if it does, inspect whether candidate ordering or
array contiguity changed.

## Acceptance criteria

- Partition output is byte-for-byte identical to the old implementation
  on focused tests.
- Golden outputs unchanged.
- `partition_and_free` improves by at least 25% on the VCAP profile, or
  the PR documents why the remaining cost is allocation rather than
  traversal.
- No peak RSS regression above 250 MB.
- Candidate and unit arrays remain contiguous and have expected dtypes.

## Risks

- This function returns many arrays, so ownership/capsule bugs are the
  main correctness risk. Keep wrappers simple and test with repeated
  garbage collection.
- Parallel scatter can increase peak allocation if every thread handles a
  large locus at once. Start with fused single-thread; add bounded
  parallelism only after measuring.

## Non-goals

- Do not change EM extraction semantics.
- Do not alter locus construction.
- Do not combine this with float32 or buffer dtype changes unless those
  PRs have already landed.