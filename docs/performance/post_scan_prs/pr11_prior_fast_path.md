# PR 11: Prior Fast Path

**Position in roadmap:** Seventh. Lower priority than scoring / scatter
work because `compute_eb_gdna_priors` is now 4.23 s. Still a clean
Python-bound target with a clear separation between EM-required outputs
and diagnostic outputs.

## Summary

Add a config-gated fast path that computes only the EM-consumed prior
outputs (`gdna_prior_count` per multi-locus, `enable_gdna` per
multi-locus). Skip the per-locus locoregional diagnostic estimates that
EM does not need. Move `enable_gdna` computation earlier in the
pipeline so it does not require a second Python gather over the global
CSR.

## Motivation

cProfile breakdown of the remaining 4.23 s:

| Function | cProfile cumulative |
|---|---:|
| `_compute_locus_scratch` | 2.57 s |
| `estimate_locus_gdna` | 1.30 s |
| `contained_exposure_clipped` | 1.25 s |
| `partition_units_to_loci` | 1.09 s |
| `compute_all_transcript_eff_lens` | 0.99 s |
| `enable_gdna_for_multilocus` | 0.66 s |

EM consumes the global `eta_g` prior count and the per-multi-locus
`enable_gdna` flag. The locoregional `estimate_locus_gdna` output is a
diagnostic; it is *not* on the EM path.

## Sequencing constraint (do not skip)

`enable_gdna_for_multilocus` reads `is_spliced` and `gdna_log_liks`
from the global CSR. `partition_and_free` consumes those arrays.
**Eligibility must be computed before scatter**, and PR08's fused
scatter must not change that ordering.

The fast path proposed here keeps the existing call ordering. The only
change is that eligibility can be computed by a vectorized native
helper instead of the current Python gather.

## Current code

* Prior assembly: [src/rigel/calibration/locus_prior.py](../../../src/rigel/calibration/locus_prior.py)
  (`assemble_priors`, lines ~778–942; `enable_gdna_for_multilocus`
  ~line 725).
* Unit→locus binning: [src/rigel/calibration/_locus_n_obs.py](../../../src/rigel/calibration/_locus_n_obs.py).
* Pipeline call site: [src/rigel/pipeline.py](../../../src/rigel/pipeline.py).
* Result schema: [src/rigel/calibration/_result.py](../../../src/rigel/calibration/_result.py)
  (`PriorTable`).

## Proposed change

### Config

Add to `CalibrationConfig` in [src/rigel/config.py](../../../src/rigel/config.py):

```python
emit_prior_diagnostics: bool = True   # default preserves current output
```

### Refactor `assemble_priors`

Split into three functions:

* `_prior_setup(...)` — region arrays, payload arrays, region index,
  boundary-crossing exposure, transcript refs/starts. Shared by both
  paths.
* `_assemble_priors_full(setup) -> PriorTable` — current behaviour.
* `_assemble_priors_fast(setup) -> PriorTable` — global `eta_g` only;
  diagnostic fields filled with documented sentinels (`np.nan` for
  floats; explicit empty arrays for per-locus diagnostics).

`assemble_priors(config, ...)` dispatches based on
`config.emit_prior_diagnostics`.

### Vectorize `enable_gdna_for_multilocus`

Today this is a Python loop over multi-loci that gathers per-unit
`is_spliced` and `gdna_log_liks` slices. Move it to a single vectorized
pass that:

* slices once with `partition_units_to_loci`'s output (already cached
  if the caller is `assemble_priors`),
* applies the existing eligibility predicate column-wise.

If the predicate is non-trivial, add a small native helper. Prefer
NumPy first; reach for native only if NumPy stays > 1 s after PR05's
profiler shows it.

### PriorTable representation

* Always has `gdna_prior_count` and `enable_gdna` populated.
* Diagnostic fields use `np.nan` (float) or zero-length arrays (lists)
  in fast mode. **Do not** fill diagnostic fields with misleading
  zeros.
* Output writers must check for the sentinels and either omit the
  diagnostic columns or write them as nulls. Document the schema
  change in [docs/MANUAL.md](../../MANUAL.md).

## Implementation steps

1. Add `emit_prior_diagnostics` to `CalibrationConfig`. Default
   `True`.
2. Extract `_prior_setup`, `_assemble_priors_full`,
   `_assemble_priors_fast`. Land this refactor as commit A; assert
   `_assemble_priors_full` is bit-identical to the old function.
3. Vectorize `enable_gdna_for_multilocus`. Land as commit B.
4. Update output writers to handle sentinel diagnostic fields.
5. Add CLI flag only if the team wants this user-facing. Otherwise
   keep it internal for benchmarking until output policy is settled.

## Tests

```bash
conda activate rigel
pytest tests/test_calibration_result.py tests/test_calibrate_orchestrator.py -v
pytest tests/test_locus_partition.py tests/test_per_locus_gdna_mass.py -v
pytest tests/test_pipeline_smoke.py tests/test_golden_output.py -v
```

Add focused tests:

* `_assemble_priors_full` produces byte-identical output to the
  pre-refactor function on synthetic fixtures.
* Fast and full modes produce identical `gdna_prior_count` and
  `enable_gdna` arrays on synthetic fixtures.
* Output writers handle sentinel diagnostic values without crashing
  and without writing misleading zeros.

## Benchmark plan

```bash
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr11_prior_fast \
  --stages --threads 8 --memory-interval 250
```

Compare `compute_eb_gdna_priors` between full and fast modes.

## Acceptance criteria

* Default mode (`emit_prior_diagnostics=True`) preserves goldens.
* Fast mode produces identical EM inputs (`gdna_prior_count`,
  `enable_gdna`).
* Fast mode reduces `compute_eb_gdna_priors` by ≥ 50% on VCAP. (4.23 s
  → ≤ 2.0 s.)
* Output schema in fast mode is documented and tested.

## Risks

* Diagnostics are scientifically useful. Do not silently drop them in
  the default user path.
* PriorTable schema changes can ripple to downstream tools. Sentinel
  values are explicit; prefer them over the misleading-zero
  alternative.
* Vectorizing `enable_gdna_for_multilocus` must respect the
  pre-`partition_and_free` ordering constraint. Add an assertion in
  the pipeline that gDNA likelihood arrays are still live when
  eligibility is computed.

## Non-goals

* No prior-formula changes.
* No EM gDNA-eligibility-semantics changes.
* No removal of diagnostic functionality, only of unconditional
  computation.
