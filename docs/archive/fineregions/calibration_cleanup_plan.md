# Calibration Cleanup Plan

Date: 2026-05-24
Status: implementation-ready
Goal: strip current calibration down to a small, boring substrate so the next
density model and a working EM path can land cleanly.

## 0. Why

Calibration today carries scaffolding from several abandoned designs:
strand-deconvolution-as-prior, regional exposure, kappa-d self-training,
`RegionExposure`, `RegionGdnaEstimate`, multi-flag region taxonomies, and a
`NotImplementedError` boundary in `pipeline.py`. The roadmaps (v1-v3) keep
proposing more layers on top.

We are stuck because the substrate is heavy. The cleanup below removes the
non-essential parts in one focused pass, then leaves exactly the pieces a
density-first calibration needs:

- payload arrays (already good);
- region geometry (already good);
- FL models (already good);
- one global gDNA density scalar per class (already good, will be replaced);
- an empty seat where the density model lands next.

Strand deconvolution is preserved as a **diagnostic-only** output until the
density model is working end-to-end. It does not feed priors. It does not
gate behavior.

## 1. Outcome After Cleanup

```text
src/rigel/calibration/
  __init__.py
  _arrays.py            kept
  _diagnostics.py       kept
  _exposure.py          trimmed
  _fl_sources.py        kept
  _orchestrator.py      simplified
  _result.py            simplified
  density_global.py     trimmed (pooled scalars only)
  fl.py                 kept
  fractional_evidence.py kept
  regions.py            kept
  scan_payload.py       kept
  signature.py          kept
  strand_summary.py     kept
  strand_deconv.py      kept, demoted to diagnostic
```

Removed concepts:

- `RegionExposure`, `exposure.py`;
- per-region exposure weights and "regional exposure" anything;
- kappa_d exon self-training;
- `region_gdna` as a required field on `CalibrationResult`;
- the Phase 6 boundary `NotImplementedError` in `pipeline.py`;
- all flag taxonomies past `FLAG_INELIGIBLE` and `FLAG_NEAR_UNSTRANDED`.

After cleanup the calibration result carries: FL models, pooled global gDNA
densities per class, optional strand diagnostic, and `n_multi_loci`. Nothing
more.

## 2. Non-Goals

This cleanup does not:

- design or implement the density model (separate doc);
- change C++ scanner code;
- change `RegionArrays`, `PayloadArrays`, `scan_payload.py`,
  `signature.py`, or FL training;
- alter pipeline EM internals beyond removing the
  `NotImplementedError` and wiring a trivial uniform prior.

## 3. Concrete Targets

### 3.1 Files to delete

```text
src/rigel/calibration/exposure.py
tests/test_exposure.py                       (if it exists)
tests/test_region_exposure*.py               (any leftovers)
tests/test_regional_exposure*.py             (any leftovers)
```

### 3.2 Files to trim

`src/rigel/calibration/_result.py`:

- drop `region_exposure` field and all references;
- make `region_gdna: RegionGdnaEstimate | None = None` (was required);
- drop `_strand_deconv_summary` requirement; instead conditionally include
  the strand block only when `region_gdna is not None`;
- summary key for strand becomes `strand_deconv_diagnostic` to make the
  demotion explicit;
- remove `region_exposure` key from `to_summary_dict`.

`src/rigel/calibration/_orchestrator.py`:

- drop `RegionExposure.uniform(...)` and the `region_exposure` argument
  threaded into `build_calibration_result`;
- keep strand deconvolution call, but mark its output as diagnostic;
- add a single config knob slot `enable_strand_diagnostic: bool = True` so
  strand work can be turned off cheaply when iterating on density;
- ensure the orchestrator returns even when strand cannot be computed (no
  hard failure on near-unstranded libraries).

`src/rigel/calibration/strand_deconv.py`:

- remove `screen_no_rna_exons` and `kappa_d` exon self-training; replace
  `estimate_kappa_d` with the pooled MoM-only version that uses
  `estimate_strand_balance` and falls back to `kappa_d = 2.0` when MoM
  fails;
- delete unused flags past `FLAG_INELIGIBLE`, `FLAG_NEAR_UNSTRANDED`,
  `FLAG_KAPPA_FALLBACK`, `FLAG_APPROX_NORMAL`;
- delete the normal-approximation branch if it is significantly more code
  than the exact path, and keep only exact discrete over
  `n_total <= 50`. For larger `n_total` set `rna_lower_count = 0`,
  `mean_count = n_total`, `upper_count = n_total`, `precision = 0`, and
  raise `FLAG_INELIGIBLE`. This is a deliberate simplification; the
  density model owns the prior anyway;
- target: file under 250 lines after trim.

`src/rigel/calibration/density_global.py`:

- keep `l_eff_contained`, `compute_global_densities`, `GlobalGdnaDensity`,
  `GlobalDensityTable`, `estimate_strand_balance`, `StrandBalanceEstimate`;
- remove anything related to exon composite densities, kappa fields that are
  unused, and any helper that was only called by `screen_no_rna_exons`;
- target: file under 250 lines after trim.

`src/rigel/calibration/_exposure.py`:

- delete `bp_weighted_mean_exposure_over_blocks`. With `RegionExposure`
  gone there is no exposure to weight by;
- keep `_merged_blocks`, `gdna_eff_len_for_loci`, `contained_exposure_clipped`,
  `fractional_boundary_side_exposure`. These are pure FL/geometry helpers
  the density model will reuse.

`src/rigel/pipeline.py`:

- remove the `NotImplementedError("rigel quant: locus EM lands in Phase 6")`
  raise;
- add a tiny helper `build_uniform_prior_table(loci, fl_models)` that
  returns the four arrays EM consumes:

```text
gdna_prior_count_em = zeros(L,  float64)
gdna_eff_len        = gdna_eff_len_for_loci(...)
enable_gdna         = ones(L,   uint8)
gdna_prior_count    = zeros(L,  float64)   # same as _em until density lands
```

- wire `_run_locus_em_partitioned` with these arrays;
- do not introduce any density logic in `pipeline.py`; the density module
  will provide a replacement function with the same signature
  (`build_prior_table(loci, fl_models, region_arrays, payload_arrays,
  fl_models, density_model) -> PriorTable`).

### 3.3 CLI / config knob trims

`src/rigel/config.py`:

- delete every field referencing exposure, regional exposure, capture, or
  kappa self-training that is no longer read;
- keep `rna_lower_confidence` only if the demoted strand diagnostic still
  uses it. Otherwise delete it too.

`src/rigel/cli.py`:

- delete any flag whose config field was removed;
- add no new flags in cleanup. The density model will introduce one or two
  later if needed.

### 3.4 Tests to delete

```text
tests/test_exposure.py
tests/test_region_exposure*.py
tests/test_regional_exposure*.py
tests/test_prior_table*.py             (if exposure-driven)
tests/test_assemble_priors.py          (already deleted in Phase 1, verify)
tests/test_calibrate_orchestrator.py   (already deleted, verify)
```

### 3.5 Tests to keep

```text
tests/test_strand_deconv.py            keep, but skip any cases that
                                       depended on exon self-training
tests/test_calibration_*.py            keep, drop exposure-dependent cases
tests/test_buffer.py, test_scan.py     unaffected
tests/golden/*                         regenerate if pipeline output text
                                       changed
```

A targeted re-run of `pytest tests/ -q` after each step is required.

## 4. Step-By-Step Order

Each step is a single PR-sized change.

### Step 1: Remove `RegionExposure`

- Delete `src/rigel/calibration/exposure.py` and any test importing it.
- Drop `region_exposure` from `_result.py`, `_orchestrator.py`,
  `_diagnostics.py`, and summary JSON.
- Delete `bp_weighted_mean_exposure_over_blocks` and its tests.
- Replace any `region_exposure.uniform(...)` call sites with nothing.
- Run `pytest tests/`. Triage and delete tests that only existed to
  validate `RegionExposure` plumbing.

### Step 2: Demote strand deconvolution

- Remove `screen_no_rna_exons` from `strand_deconv.py`.
- Replace `estimate_kappa_d` with a MoM-only function from
  `estimate_strand_balance`, with a hard fallback to `kappa_d = 2.0`.
- Remove `FLAG_EXON_SELF_TRAIN` and any related logic.
- Make `region_gdna` optional everywhere it appears.
- Update `to_summary_dict` to write `strand_deconv_diagnostic` only when
  the field is populated.
- Update or skip strand tests that exercised exon self-training; do not add
  new tests in this step.

### Step 3: Trim density_global

- Delete `ExonCompositeDensity` if still present and any related summary
  fields.
- Reduce `density_global.py` to the public surface listed in 3.2.
- Ensure `compute_global_densities` still returns INTERGENIC and INTRON-ONLY
  rows with `rho` and `n_regions_used` populated. These are the seeds the
  density model will consume.

### Step 4: Unblock pipeline.py

- Delete the `NotImplementedError` raise.
- Add `build_uniform_prior_table` in `pipeline.py` as described.
- Call it from `quant_from_buffer` so EM runs end-to-end with zero gDNA
  prior counts. `enable_gdna` is 1; `gdna_prior_count_em` is 0.
- This is the temporary path. The density model will replace
  `build_uniform_prior_table` with a real prior table without changing the
  EM call site.

### Step 5: Add a thin pipeline end-to-end test

- New test `tests/test_pipeline_uniform_prior.py`: tiny synthetic BAM,
  run `quant_from_buffer`, assert no exception, quant output has expected
  shape and `enable_gdna` is honored.
- This becomes the regression anchor that lets us iterate on density.

### Step 6: Sweep dead code

- Run `grep -R "RegionExposure" src/ tests/` and remove stragglers.
- Run `grep -R "FractionalCutover" src/ tests/` and remove stragglers.
- Run `grep -R "regional_exposure" src/ tests/` and remove stragglers.
- Run `grep -R "screen_no_rna_exons" src/ tests/` and remove stragglers.
- `ruff check src/ tests/` clean.
- `pytest tests/ -q` green.

## 5. Acceptance Gates

After all six steps:

```text
grep -R "RegionExposure"           src/ tests/   -> no matches
grep -R "regional_exposure"        src/ tests/   -> no matches
grep -R "screen_no_rna_exons"      src/ tests/   -> no matches
grep -R "NotImplementedError.*locus EM" src/    -> no matches
wc -l src/rigel/calibration/*.py                 -> total < 3000 lines
pytest tests/ -q                                 -> green
rigel quant on a tiny synthetic BAM              -> produces quant.feather
```

`rigel quant` works end-to-end with **zero gDNA prior**. RNA may be slightly
overestimated, gDNA may be underestimated, and that is acceptable for the
duration of the cleanup. The density model will fix the numbers next.

## 6. What Cleanup Explicitly Does NOT Touch

- The C++ scanner and payload schema.
- `RegionArrays.signature`, `ts_class`, payload channel layout.
- FL model training and pool quality logic.
- The EM kernel.

If a "small" cleanup edit looks like it needs to change any of these, stop
and call it out. It does not belong in this cleanup pass.

## 7. Risks and Counter-Risks

- **Risk:** removing `region_gdna` from `CalibrationResult` breaks a
  downstream consumer.
  **Mitigation:** keep the field optional, default None. Downstream reports
  must guard against None.
- **Risk:** removing kappa_d self-training degrades strand diagnostics.
  **Mitigation:** acceptable. The strand path is demoted to diagnostic-only
  and will be revisited once the density model is in place.
- **Risk:** end-to-end EM with zero gDNA prior over-assigns to RNA on real
  data.
  **Mitigation:** acceptable for one PR cycle. Density model lands next.

## 8. Definition of Done

When all of the following are true:

```text
- the calibration/ directory matches the layout in Section 1
- rigel quant runs end-to-end on a tiny synthetic BAM
- the pipeline accepts a PriorTable-shaped object whose four arrays are
  the only contract
- the rest of the calibration module surface area is at least one third
  smaller than before
- there is no exposure code, no kappa exon self-training, no fractional
  cutover sentinel, and no NotImplementedError boundary
```

Cleanup is finished. The density model can land against a small, stable
substrate.
