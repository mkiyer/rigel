# PR 12: EM High-Iteration Workset

**Position in roadmap:** Eighth and explicitly investigative. Do this
after memory/scatter/scoring work unless EM becomes the dominant stage.

## Summary

Build a focused diagnostic and benchmark harness for high-iteration loci,
then evaluate low-risk EM improvements such as better warm starts,
active-set pruning, and deterministic-equivalence-class shortcuts.

This PR should produce evidence and small safe improvements, not broad
solver rewrites.

## Motivation

`locus_em` is 14.30 s. It is already threaded, and the largest mega-locus
accounts for only 3.80 s of summed locus time. Several smaller loci take
0.3-0.7 s because they require many SQUAREM iterations; some report 333
iterations.

Changing convergence tolerances would be easy and wrong. The right next
step is to characterize the hard loci and find model-preserving ways to
reduce work.

## Current code

- Batch EM driver: [src/rigel/estimator.py](../../../src/rigel/estimator.py)
- Native EM solver: [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp)
- Locus stats emission: [scripts/profiling/profiler.py](../../../scripts/profiling/profiler.py)

## Proposed work

Create a repeatable high-iteration workset:

- top loci by total EM time,
- top loci by SQUAREM iterations,
- top loci by equivalence-class elements,
- top loci where assignment/build-EC dominates rather than SQUAREM.

Then test targeted solver improvements on that workset before touching
the production path.

## Implementation steps

1. Add a script under `scripts/profiling/` or `scripts/debug/` that reads
   `locus_stats_default.json` and emits a compact workset table.
2. Add optional profile output that includes enough locus identity to
   reproduce a specific locus: multi-locus id, transcript count, unit
   count, EC count, iteration count, and timing breakdown.
3. Add a debug entry point that can run EM for selected locus ids from a
   saved partition bundle or a freshly built profile run. If saving
   partitions is too heavy, document the limitation and use full-run
   filtering first.
4. Investigate warm-start quality for high-iteration loci:
   - compare current warm start to uniform,
   - compare to prior-weighted starts,
   - report iterations and final log likelihood.
5. Investigate active-set pruning only after warm-start data is collected:
   components with persistently tiny posterior mass may be removable from
   subsequent iterations if final assignment remains unchanged within
   tolerance.
6. Add deterministic-EC shortcut candidates: loci where every EC has a
   single eligible component, or where rows are effectively deterministic
   after likelihood pruning.
7. Implement only the smallest proven improvement. Keep diagnostics from
   unsuccessful ideas in the PR notes.

## Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_em_impl.py tests/test_estimator.py tests/test_em_prior_weight.py -v
pytest tests/test_golden_output.py -v
```

For any solver change, add a focused test that demonstrates identical or
numerically equivalent results on a synthetic hard locus.

## Benchmark plan

Use two levels:

1. Workset microbenchmarks for selected loci, repeated several times.
2. Full VCAP staged profile to confirm macro impact.

Full profile command:

```bash
conda activate rigel
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr12_em_workset \
  --stages --threads 8 --memory-interval 250
```

## Acceptance criteria

- The workset report identifies and ranks hard loci reproducibly.
- Any production solver change preserves golden outputs or has a clearly
  bounded numerical diff.
- No convergence tolerance is loosened.
- If a solver optimization lands, `locus_em` improves by at least 5% on
  VCAP or the PR remains diagnostics-only.

## Risks

- Solver changes can silently bias abundance estimates. Keep changes
  small, measured, and backed by golden outputs.
- Microbenchmark wins may not survive full-run scheduling. Require a full
  profile before declaring victory.

## Non-goals

- Do not rewrite SQUAREM.
- Do not change the model or priors.
- Do not tune tolerances as a performance fix.