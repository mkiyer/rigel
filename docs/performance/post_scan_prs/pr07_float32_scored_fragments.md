# PR 07: Float32 Scored-Fragment Payloads

**Position in roadmap:** Third. Land before PR08 (fused scatter) and
PR09 (parallel scoring) so they operate on the smaller payload.

## Summary

Convert the per-candidate / per-unit *log-likelihood* and
*coverage-weight* payloads from `double` to `float32` in the global and
partitioned CSR. Promote to `double` only where reductions need it
(EM E-step accumulation and SQUAREM updates).

Target arrays:

* `log_liks` (per candidate)
* `coverage_weights` (per candidate)
* `gdna_log_liks` (per unit) — included if numerical validation passes

## Motivation

Scoring/router materializes the second-largest resident plateau (9.81
GB → 16.20 GB). At ~263 M candidates,
`log_liks` and `coverage_weights` are ~2.1 GB each as `double`. Halving
the dtype saves ~2.1 GB before considering partition copies and memory
bandwidth wins.

These payloads are normalized (log-sum-exp pivots) before they affect
posteriors, and `coverage_weights` is now used only for the EM
warm-start `θ_init` (see
[docs/cleanup/rna_prior_cleanup.md](../../cleanup/rna_prior_cleanup.md)).
Float32 is sufficient for the stored values; the EM solver already
performs Kahan-compensated double accumulation.

## Current code

* Python container: [src/rigel/scored_fragments.py](../../../src/rigel/scored_fragments.py)
  — `log_liks`, `coverage_weights`, `gdna_log_liks` are all `np.float64`.
* Native scorer: [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp)
  — `StreamingScorer` holds `std::vector<double>* v_ll_`, `v_cw_`,
  `v_gll_`. `finish()` returns an 18-element tuple including these as
  `double` arrays.
* Scatter: [src/rigel/locus_partition.py](../../../src/rigel/locus_partition.py)
  + the templated `scatter_candidates_impl<T>` /
  `scatter_units_impl<T>` in
  [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp).
  Today only `f64`, `i32`, `u8`, `i64` instantiations are bound.
* EM consumers of `coverage_weights`: warm-start `θ_init` in
  `compute_ovr_prior_and_warm_start` (em_solver.cpp ~lines 765–787) and
  the secondary tie-break in `build_equiv_classes` (~lines 281–296).
  No SQUAREM E-step path reads `coverage_weights` directly.

## Numerical position

* `log_liks`: stored values are pre-normalized log-likelihoods. Their
  absolute value is bounded by the scoring kernel and they are
  immediately exponentiated against a row pivot. Float32 dynamic range
  (~38 decades) is plenty; precision (~7 decimal digits) is sufficient
  given that posteriors are normalized.
* `coverage_weights`: only used to seed warm-start θ. Float32 is
  obviously sufficient for a starting point.
* `gdna_log_liks`: same family as `log_liks` but per-unit. Almost
  certainly safe; the PR must demonstrate it on golden outputs.

## Implementation steps

1. Add `f32_1d` and `f32_1d_mut` nanobind aliases in `scoring.cpp` and
   `em_solver.cpp` (and any shared header).
2. In `StreamingScorer`, change `v_ll_`, `v_cw_` (and conditionally
   `v_gll_`) from `std::vector<double>*` to `std::vector<float>*`.
   Push values with explicit `static_cast<float>(...)`. Keep all
   intermediate scoring arithmetic in `double`.
3. Update `finish()` to wrap these vectors as `float32` numpy arrays.
   Document the shape change in the function docstring.
4. Add the `float` instantiation of `scatter_candidates_impl<T>` and,
   if `gdna_log_liks` becomes float32, `scatter_units_impl<T>`. Bind
   them as `scatter_candidates_f32` / `scatter_units_f32`.
5. Update [src/rigel/scored_fragments.py](../../../src/rigel/scored_fragments.py)
   dtype declarations and any docstrings.
6. Update `partition_and_free` to dispatch the f32 scatter for the
   payload columns.
7. In `em_solver.cpp`'s subproblem extraction (`PartitionView` →
   `LocusSubProblem`), accept float32 inputs and cast to `double` when
   filling `LocusSubProblem.log_liks` / `coverage_wts`. Preserve the
   `LocusSubProblem` representation as `double` to avoid any change to
   the SQUAREM hot loop.
8. Update fixtures and tests that explicitly construct `np.float64`
   payloads.
9. Add a regression test that asserts the production scorer returns
   `float32` for the converted columns.

## Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_em_impl.py tests/test_estimator.py \
       tests/test_locus_partition.py -v
pytest tests/test_pipeline_smoke.py tests/test_pipeline_routing.py -v
pytest tests/test_golden_output.py -v
pytest tests/ -q
```

If golden outputs shift, do **not** regenerate without a numeric
analysis. Acceptable noise floor: per-transcript abundance relative
error < 1e-5 on the synthetic scenarios; if a real-data diff is
larger, the PR is wrong.

## Benchmark plan

```bash
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr07_float32 \
  --stages --threads 8 --memory-interval 250
```

Compare against PR05/PR06 baseline at the same scan settings:

* peak RSS
* `fragment_router_scan`, `partition_and_free`, `locus_em`
* PR05's `scoring_csr.candidate_bytes.{log_liks,coverage_weights}`

## Acceptance criteria

* Golden outputs unchanged within the 1e-5 relative-error noise floor.
* Peak RSS drops by ≥ 1.5 GB on VCAP.
* `scoring_csr.candidate_bytes.log_liks` is exactly half its previous
  value.
* `fragment_router_scan`, `partition_and_free`, and `locus_em` do not
  regress by more than 5%.
* All EM reductions still use `double` accumulation. Document this in
  the relevant function comments.

## Risks

* Mixed-precision pipelines drift. Use clear aliases (`f32_1d`,
  `f64_1d`) and put a one-line type comment at every conversion site.
* `gdna_log_liks` is the riskiest of the three because the
  per-fragment gDNA likelihood enters the EM E-step directly. If
  golden outputs shift outside the noise floor, **leave it as
  `double`** in this PR and follow up.
* Tests hard-code `np.float64`; expect ~10–20 fixture updates.

## Non-goals

* No parallel scoring (PR09).
* No removal of `coverage_weights` (the warm-start consumer is alive).
* No EM tolerance changes.
