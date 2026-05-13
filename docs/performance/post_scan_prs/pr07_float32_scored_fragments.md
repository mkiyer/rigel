# PR 07: Float32 Scored-Fragment Payloads

**Position in roadmap:** Third. Land before parallel scoring and fused
scatter so later PRs operate on the smaller payload.

## Summary

Store per-candidate scoring payloads as float32 in the global and
partitioned CSR representations while preserving double-precision
accumulation in numerically sensitive EM reductions.

Target arrays:

- `log_liks`
- `coverage_weights`
- `gdna_log_liks` if validation shows no numerical issue

## Motivation

The scoring/router stage creates the second major RSS jump: 9.81 GB to
16.20 GB. The biggest individual arrays are per-candidate float64
payloads. At roughly 260 M candidates in the initial profile,
`log_liks` and `coverage_weights` are about 2.08 GB each. Converting
both to float32 saves about 2.1 GB before considering partition copies
and bandwidth effects.

This also makes downstream scatter and extraction lighter. The EM solver
can still promote values to double in row-local scratch and accumulation.

## Current code

- Python container dtype docs: [src/rigel/scored_fragments.py](../../../src/rigel/scored_fragments.py)
- Native scoring output vectors: [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp)
- Partition scatter bindings: [src/rigel/locus_partition.py](../../../src/rigel/locus_partition.py), [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp)
- EM partition view and extraction: [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp)

## Numerical position

`log_liks` are row-normalized by log-sum-exp style pivots before they
affect posterior probabilities. Float32 precision is sufficient for the
stored likelihood payload if the EM inner loop promotes to double for
accumulation. `coverage_weights` are warm-start/prior weights, not used
in the EM E-step hot loop. Float32 is sufficient for those values.

This PR must prove that with tests, not assume it.

## Implementation steps

1. Add nanobind aliases for float32 one-dimensional arrays in the native
   modules, for example `f32_1d`.
2. In [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp),
   change `StreamingScorer` vectors for `v_ll_` and `v_cw_` from
   `std::vector<double>` to `std::vector<float>`.
3. Push values with explicit `static_cast<float>(...)`. Keep local
   scoring arithmetic in double.
4. Decide whether `v_gll_` also becomes `std::vector<float>` in the same
   PR. Preferred: include it if tests pass, because it is a pure stored
   likelihood payload. If not, document why it stays double.
5. Add `scatter_candidates_f32` and, if needed, `scatter_units_f32` native
   bindings in [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp).
6. Update [src/rigel/locus_partition.py](../../../src/rigel/locus_partition.py)
   to scatter float32 payloads with the new bindings.
7. Update `PartitionView` input pointer types for float32 payloads while
   keeping `LocusSubProblem.log_liks` and `coverage_wts` as double if
   those arrays are reused heavily by SQUAREM.
8. In extraction, cast float32 input to double when writing row-local
   `LocusSubProblem` arrays.
9. Update dtype annotations and docstrings in
   [src/rigel/scored_fragments.py](../../../src/rigel/scored_fragments.py).
10. Update tests and fixtures that explicitly construct float64
    `log_liks`, `coverage_weights`, or `gdna_log_liks` arrays.
11. Add a small native/Python test that asserts the production scorer
    returns float32 arrays for the selected payloads.

## Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_em_impl.py tests/test_estimator.py tests/test_locus_partition.py -v
pytest tests/test_pipeline_smoke.py tests/test_pipeline_routing.py -v
pytest tests/test_golden_output.py -v
pytest tests/ -q
```

Add or update tests to compare float32 and float64 fixture paths on small
synthetic loci. Use tight but realistic tolerances and inspect any golden
diff before accepting it.

## Benchmark plan

Run clean and cProfile staged profiles:

```bash
conda activate rigel
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr07_float32 \
  --stages --threads 8 --memory-interval 250
```

Compare:

- peak RSS
- `fragment_router_scan`
- `partition_and_free`
- `locus_em`
- output abundance diffs against the pre-PR profile run

## Acceptance criteria

- Golden outputs pass without regeneration, or numeric changes are
  explained and below accepted output tolerances.
- Peak RSS drops by at least 1.5 GB on the VCAP full profile.
- `fragment_router_scan`, `partition_and_free`, and `locus_em` do not
  regress by more than 5%.
- Profile JSON from PR05 shows candidate payload bytes are roughly halved.
- EM uses double accumulation for reductions and final abundance counts.

## Risks

- ABI changes ripple through tests and fixtures.
- Some small probabilities may be sensitive if float32 values are used
  directly inside repeated SQUAREM updates. Mitigate by promoting to
  double inside per-locus subproblem arrays or the E-step.
- Mixed precision can become confusing. Keep clear type names and
  docstrings.

## Non-goals

- Do not parallelize scoring in this PR.
- Do not remove `coverage_weights`; only shrink representation.
- Do not change EM convergence thresholds.