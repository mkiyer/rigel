# PR 12: EM High-Iteration Workset

**Position in roadmap:** Eighth. Investigative. Run after the
memory / scatter / scoring PRs unless EM becomes the dominant stage.

## Summary

Build a focused diagnostic and microbenchmark harness for high-iteration
loci, then evaluate **model-preserving** EM improvements (warm-start
quality, active-set pruning of components with persistently negligible
posterior mass, deterministic-EC shortcuts). This PR produces evidence
and at most one small, low-risk solver tweak.

## Motivation

`locus_em` is 14.30 s. The largest single locus is only 3.80 s of that;
the rest is many medium loci with high SQUAREM iteration counts (some
reported at 333 iterations).

| Substage (summed across loci) | Time |
|---|---:|
| `squarem_us` | 72.76 s |
| `assign_us` | 4.66 s |
| `build_ec_us` | 4.37 s |
| `extract_us` | 1.87 s |
| `warm_start_us` | 0.29 s |

(Wall = 14.30 s because EM is parallel across loci.)

Loosening `convergence_delta` would reduce iterations and ship the wrong
answer. The right next step is to characterize the hard loci and search
for model-preserving wins.

## Current code

* Batch driver: [src/rigel/estimator.py](../../../src/rigel/estimator.py).
* Native EM: [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp).
* Locus stats emission: [scripts/profiling/profiler.py](../../../scripts/profiling/profiler.py).

## Workstream

### 1. Workset

Add a script under `scripts/profiling/em_workset.py` that reads
`locus_stats_default.json` and emits a compact ranked table:

* top-N loci by total EM time
* top-N by SQUAREM iterations
* top-N by EC count
* top-N where `assign_us + build_ec_us` exceeds `squarem_us`

Keep N = 20 by default; CLI override for full dump.

### 2. Replay harness

Add a script that can re-run EM for a single locus_id from a saved
`PartitionBundle`. Saving the bundle is heavyweight; instead implement
a lighter path that:

* runs the full pipeline up to `partition_and_free`,
* selects the locus IDs of interest,
* invokes the native single-locus EM with each,
* reports per-iteration log-likelihood, posterior mass per component,
  and final assignment.

Output goes to `scripts/profiling/em_replay/<locus_id>/`.

### 3. Hypotheses to test (in this order)

#### a. Warm-start quality

Compare the current warm start (coverage-weighted θ_init) against:

* uniform θ
* prior-weighted θ (scaled by `gdna_prior_count`)
* posterior from one EM iteration of a coarsened model

Report iterations to convergence and final log-likelihood. Choose the
warm start that minimizes iterations on the workset *without changing
the converged answer beyond `convergence_delta`*.

#### b. Active-set pruning

For loci where many components hold persistently small posterior mass
(say < 1e-12 for ≥ 10 iterations), test removing them from the
remaining iterations. Re-introduce only if any of their candidate
likelihoods exceeds a re-entry threshold. Verify converged θ matches
the un-pruned solution within `convergence_delta`.

#### c. Deterministic-EC shortcut

Identify ECs whose rows are effectively deterministic after likelihood
pruning (single eligible component per row). Move their contribution
to a closed-form update instead of the SQUAREM loop.

### 4. Solver change (optional)

Implement *at most one* of (a)/(b)/(c) in this PR, gated behind a
config flag, default off. The PR can be diagnostics-only and still
land — write up the negative results.

## Implementation steps

1. Build the workset script. Land it.
2. Build the replay harness. Land it.
3. Run experiments. Document results in
   [docs/performance/em_workset_findings.md](../em_workset_findings.md).
4. If a tweak is justified, implement it behind a config flag and
   add an explicit numerical test on a synthetic hard locus.

## Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_em_impl.py tests/test_estimator.py \
       tests/test_em_prior_weight.py -v
pytest tests/test_golden_output.py -v
```

Any solver change requires a new test: a synthetic hard locus where
the new path produces an answer numerically equivalent to the old path
within `convergence_delta`.

## Benchmark plan

Two levels:

1. **Microbench** on the workset, repeated 5×, take median.
2. **Full VCAP staged profile** to confirm macro impact.

```bash
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr12_em_workset \
  --stages --threads 8 --memory-interval 250
```

## Acceptance criteria

* Workset report lands and is reproducible from any profile run.
* Replay harness lands and produces per-iteration diagnostics.
* If a solver change ships:
  * Goldens unchanged.
  * `locus_em` improves by ≥ 5% on VCAP, *or* the change documents a
    specific hard-locus class that benefits.
  * `convergence_delta` is **not** loosened.
* If no solver change ships, the PR delivers the documented
  experimental record and the harness for future work.

## Risks

* Solver changes can silently bias abundance estimates. Required
  guard: per-component posterior mass must match the un-changed path
  within `convergence_delta` on the synthetic hard-locus test.
* Microbenchmark wins can vanish under full-pipeline scheduling. The
  full profile is the binding test.

## Non-goals

* No SQUAREM rewrite.
* No model or prior changes.
* No tolerance loosening.
