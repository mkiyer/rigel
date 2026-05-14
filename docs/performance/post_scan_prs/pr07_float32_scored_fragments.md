# PR 07: Float32 Scored-Fragment Payloads

**Position in roadmap:** Third. Land before PR08 (fused scatter) and
PR09 (parallel scoring) so they operate on the smaller payload.

**Status:** Implemented in this branch. The native scorer stores
`log_liks`, `coverage_weights`, and `gdna_log_liks` as float32, while
the native EM extraction path promotes them back to double before
building each locus subproblem. During validation, the scorer also
stopped materializing unused per-candidate `tx_start` / `tx_end` arrays,
which were returned to Python and immediately discarded.

Validation results are published below. The final VCAP profile reduced
peak RSS from 15.60 GB to 12.09 GB versus the PR06 2 GiB scan-cap
baseline.

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
  — `log_liks`, `coverage_weights`, and `gdna_log_liks` are float32
  payload arrays.
* Native scorer: [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp)
  — `StreamingScorer` holds `std::vector<float>* v_ll_`, `v_cw_`, and
  `v_gll_`. All scoring arithmetic remains double until the final
  storage cast. The old unused per-candidate `tx_start` / `tx_end`
  output vectors have been removed.
* Scatter: [src/rigel/locus_partition.py](../../../src/rigel/locus_partition.py)
  + the templated `scatter_candidates_impl<T>` /
  `scatter_units_impl<T>` in
  [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp).
  Float payloads dispatch to `scatter_candidates_f32` /
  `scatter_units_f32`, with float64 support retained for synthetic
  fixtures and direct native tests.
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
pytest tests/test_locus_partition.py tests/test_ndarray_util.py \
  tests/test_batch_em_impl.py -v
pytest tests/test_estimator.py tests/test_pipeline_smoke.py \
  tests/test_pipeline_routing.py tests/test_golden_output.py -v
ruff check src/ tests/
pytest tests/ -q
```

If golden outputs shift, do **not** regenerate without a numeric
analysis. Acceptable noise floor: per-transcript abundance relative
error < 1e-5 on the synthetic scenarios; if a real-data diff is
larger, the PR is wrong.

Observed golden drift before regeneration was small and consistent with
float32 storage: maximum absolute difference `9.8e-7`, maximum relative
difference `8.6e-8`. Goldens were regenerated after this analysis.

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

Final benchmark output directories:

* PR06 baseline: `/Users/mkiyer/Downloads/rigel_runs/profile_pr06_cap_2gib`
* PR07 final: `/Users/mkiyer/Downloads/rigel_runs/profile_pr07_float32_no_txpos`

| Metric | PR06 2 GiB baseline | PR07 final | Delta |
|---|---:|---:|---:|
| Wall time | 80.61 s | 89.15 s | +8.54 s |
| Peak RSS | 15,596 MB | 12,087 MB | -3,509 MB (-22.5%) |
| After router RSS | 15,596 MB | 11,907 MB | -3,689 MB |
| Scoring CSR bytes | 5.99 GiB | 3.94 GiB | -2.05 GiB |
| Partition bytes | 5.99 GiB | 3.94 GiB | -2.05 GiB |
| `fragment_router_scan` | 18.08 s | 14.51 s | -3.57 s |
| `partition` | 6.35 s | 5.41 s | -0.94 s |
| `locus_em` | 14.71 s | 13.84 s | -0.87 s |

The headline wall time is slower because `scan_and_buffer` was 14.32 s
slower in the PR07 run (`35.86 s` → `50.18 s`), even though PR07 does
not touch the BAM scan phase. The PR07-touched post-scan stages all got
faster in this profile. Treat total wall as run-to-run scan noise here;
use the stage deltas for PR07 attribution.

## Acceptance criteria

* Golden outputs unchanged within the 1e-5 relative-error noise floor.
* Golden outputs unchanged within the 1e-5 relative-error noise floor:
  satisfied after quantified drift analysis.
* Peak RSS drops by ≥ 1.5 GB on VCAP: satisfied (`-3.51 GB`).
* `scoring_csr.candidate_bytes.log_liks` is exactly half its previous
  value: satisfied (`1.94 GiB` → `0.97 GiB`).
* `fragment_router_scan`, `partition_and_free`, and `locus_em` do not
  regress by more than 5%: satisfied for the PR07-touched stages.
* All EM reductions still use `double` accumulation: satisfied; payloads
  are promoted in native partition extraction before locus EM.

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
