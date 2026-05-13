# PR 11: Prior Fast Path And Eligibility Fusion

**Position in roadmap:** Seventh. Lower priority than scoring/scatter
because prior assembly is now 4.23 s, but still a clean Python-heavy
target.

## Summary

Add a production fast path for gDNA prior assembly that computes only the
EM-consumed prior count and gDNA eligibility when detailed locoregional
diagnostics are not requested. Move or fuse gDNA eligibility computation
so it does not perform a separate Python gather over every multi-locus.

## Motivation

`compute_eb_gdna_priors` improved from 9.5 s to 4.23 s after cleanup, but
the remaining cost is still visible. cProfile shows time in
`_compute_locus_scratch`, `estimate_locus_gdna`,
`contained_exposure_clipped`, `partition_units_to_loci`, and
`enable_gdna_for_multilocus`.

The key modeling fact: EM consumes the global `eta_g` prior count and an
eligibility flag. The locoregional `estimate_locus_gdna` output is useful
diagnostic material, but it is not needed to run EM.

## Current code

- Prior assembly: [src/rigel/calibration/locus_prior.py](../../../src/rigel/calibration/locus_prior.py)
- Unit-to-locus partitioning: [src/rigel/calibration/_locus_n_obs.py](../../../src/rigel/calibration/_locus_n_obs.py)
- Pipeline call site: [src/rigel/pipeline.py](../../../src/rigel/pipeline.py)
- Calibration result schema: [src/rigel/calibration/_result.py](../../../src/rigel/calibration/_result.py)

## Proposed change

Introduce a config-controlled prior diagnostic mode:

```python
class CalibrationConfig:
    emit_prior_diagnostics: bool = True
```

Default can stay `True` initially to preserve output. The fast path can
be enabled in benchmarks and later considered as default if output files
remain satisfactory.

Fast path behavior:

- compute `gdna_prior_count_arr` from global density only,
- compute `enable_gdna_arr`, preferably during partition/scatter or with
  a vectorized native helper,
- populate minimal prior summary fields,
- skip per-locus locoregional diagnostic estimates unless requested.

## Implementation steps

1. Add `emit_prior_diagnostics` or similarly named config field. Keep the
   default preserving current output.
2. Split `assemble_priors(...)` into two internal paths:
   - full diagnostic path: current behavior,
   - fast EM path: global `eta_g` plus eligibility only.
3. Factor shared setup out of both paths: `RegionArrays`,
   `PayloadArrays`, `RegionIndexPy`, `b_cross`, transcript refs/starts.
4. Implement a vectorized or native `compute_enable_gdna_by_locus(...)`
   helper that consumes `multi_loci`, `em_data.is_spliced`, and
   `em_data.gdna_log_liks` once. If PR08 has landed, compute this during
   fused scatter and return it as partition metadata.
5. Ensure `CalibrationResult.with_priors(...)` can accept a minimal prior
   table or a prior table with diagnostic fields marked unavailable.
6. Update output writers to handle missing diagnostics gracefully when the
   fast path is enabled.
7. Add CLI/config visibility only if this mode is user-facing. Otherwise
   keep it internal for profiling until behavior is settled.

## Tests

```bash
conda activate rigel
pytest tests/test_calibration_result.py tests/test_calibrate_orchestrator.py -v
pytest tests/test_locus_partition.py tests/test_per_locus_gdna_mass.py -v
pytest tests/test_pipeline_smoke.py tests/test_golden_output.py -v
```

Add focused tests:

- Full diagnostic mode produces byte-compatible prior diagnostics with
  current outputs.
- Fast mode produces identical `gdna_prior_count` and `enable_gdna` to
  full mode on synthetic fixtures.
- Output writing either omits diagnostic-only columns in fast mode or
  fills them with documented nulls.

## Benchmark plan

Run full staged profiles with diagnostics on and off:

```bash
conda activate rigel
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr11_prior_fast \
  --stages --threads 8 --memory-interval 250
```

Compare `compute_eb_gdna_priors`, output schema, and final quant outputs.

## Acceptance criteria

- Full diagnostic mode preserves existing golden outputs.
- Fast mode produces identical EM inputs: `gdna_prior_count` and
  `enable_gdna` arrays.
- Fast mode reduces `compute_eb_gdna_priors` by at least 25% on VCAP.
- The output schema behavior is documented and tested.

## Risks

- Diagnostics are scientifically useful. Do not silently drop them in the
  default user path without an explicit decision.
- Prior table schema may assume diagnostic fields exist. Keep the minimal
  representation explicit rather than using misleading zeros.

## Non-goals

- Do not change the prior formula.
- Do not change EM gDNA eligibility semantics.
- Do not remove diagnostic functionality.