# Phase IV - Clean Cutover To RegionCalibration

Date: 2026-05-26
Status: implementation plan
Scope: Replace the live calibration handoff with the v6 `RegionCalibration` path. No `use_v6_calibration` flag, no old density/fusion state, and no compatibility-gated behavior.

---

## Overview

Phase IV wires the Phase I/II/III primitives into the production calibration and EM pipeline.

This is a clean cutover. The current density/fusion path is not preserved behind a flag because it is not a useful fallback. The new path becomes the only calibration path:

```text
BAM calibration payload
-> RegionArrays / PayloadArrays / RegionCountLedger
-> DensityObservation
-> compartment strand deconvolution
-> BackgroundModel
-> BoundaryTable
-> BoundaryLocalPosterior + BoundarySweepResult
-> four-state RegionCalibration
-> PriorTable from RegionCalibration
-> locus EM
```

The old path:

```text
fit_density_evidence
-> fuse_density_and_strand
-> RegionExposure.from_density
-> calibration.fused_region_gdna
```

is removed from the live pipeline.

---

## 1. Cutover Rules

### 1.1 No feature flag

Do not add `use_v6_calibration` or any equivalent flag. Phase IV makes `RegionCalibration` the production contract.

### 1.2 No old calibration state in `CalibrationResult`

Remove these live result fields:

```python
density_evidence
region_gdna
region_exposure
fused_region_gdna
```

Replace them with:

```python
region_calibration: RegionCalibration
strand_channels: RegionGdnaChannelEstimate | None
background_model: BackgroundModel
boundary_local: BoundaryLocalPosterior
boundary_sweep: BoundarySweepResult
```

`strand_channels`, `background_model`, `boundary_local`, and `boundary_sweep` are retained for diagnostics/debugging, not for downstream decision logic. Downstream EM consumes `region_calibration`.

### 1.3 No calibration-side gDNA on/off switching

The v6 concept keeps gDNA smooth. Calibration should not disable gDNA based on a posterior threshold.

The native EM still has a technical `enable_gdna` input. In Phase IV this field means only:

```text
this locus has at least one unspliced unit with finite gDNA likelihood
```

It must not depend on `RegionCalibration.upper_gdna`, `p_captured`, or any calibration evidence threshold. A zero prior with native gDNA eligibility still leaves the gDNA component available to absorb likelihood-supported mass.

### 1.4 No old golden preservation requirement

Existing golden outputs can change. Phase IV should update tests and, when appropriate, golden data to reflect the new calibration path. There is no requirement that v6-off outputs remain unchanged because v6-off no longer exists.

---

## 2. Calibration Orchestrator Cutover

File: `src/rigel/calibration/_orchestrator.py`

### 2.1 Remove old imports/path

Remove live use of:

```python
fit_density_evidence
fuse_density_and_strand
RegionExposure.from_density
build_strand_region_counts/deconvolve_regions_by_strand as final handoff
```

`build_strand_region_counts()` and `deconvolve_regions_by_strand()` may remain used only if needed for kappa estimation or diagnostics, but the final output must be compartment-aware `RegionGdnaChannelEstimate` plus `RegionCalibration`.

### 2.2 New orchestrator flow

```python
region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
ledger = build_region_count_ledger(payload_arrays)
observation = build_density_observation(region_arrays, ledger, fl_models.gdna)

strand_summary = strand_summary or StrandSummary.uninformative()
strand_counts = build_strand_region_counts(...)
kappa_d = estimate_kappa_d(...)
compartment_counts = build_compartment_strand_counts(...)
strand_channels = deconvolve_compartments_by_strand(
    compartment_counts,
    kappa_d=kappa_d.kappa,
    rna_lower_confidence=rna_lower_confidence,
)

background = fit_background_model(
    observation,
    strand_channels=strand_channels if stranded_identifiable else None,
    top_t_fraction=background_trim_fraction,
    min_eff_length=background_min_eff_length,
)

boundaries = build_boundary_table(
    region_arrays,
    ledger,
    observation.boundary_left_leff,
)
local = build_boundary_local_posterior(
    observation,
    boundaries,
    background,
    strand_channels=strand_channels if stranded_identifiable else None,
    confidence=gdna_density_confidence,
)

region_calibration = run_calibration_iteration(
    region_arrays,
    observation,
    boundaries,
    background,
    strand_channels=strand_channels if stranded_identifiable else None,
    local_posterior=local,
    max_calibration_passes=max_calibration_passes,
    confidence=gdna_density_confidence,
)
```

### 2.3 Stranded vs unstranded input

If the library strand model is not identifiable, pass `strand_channels=None` to background, boundary, and tensor components that require strand evidence. Still compute compartment counts where useful for diagnostics, but do not let unidentifiable strand estimates drive calibration.

Concrete rule:

```text
strand_usable = _strand_summary_identifiable(strand_summary)
calibration_strand_channels = strand_channels if strand_usable else None
```

This keeps unstranded data broad and uncertainty-aware.

### 2.4 Build result

`build_calibration_result()` should receive the new objects and reject missing `region_calibration`.

---

## 3. CalibrationResult Contract

File: `src/rigel/calibration/_result.py`

### 3.1 New dataclass

```python
@dataclass(frozen=True, slots=True)
class CalibrationResult:
    fl_models: FLModels
    diagnostics: Diagnostics
    region_calibration: RegionCalibration
    strand_channels: RegionGdnaChannelEstimate | None
    background_model: BackgroundModel
    boundary_local: BoundaryLocalPosterior
    boundary_sweep: BoundarySweepResult
    prior_table: PriorTable | None = None
    n_multi_loci: int = 0
    rna_lower_confidence: float = 0.95
```

Remove old fields from the dataclass, summary, tests, and constructors:

```python
density_evidence
region_gdna
region_exposure
fused_region_gdna
```

### 3.2 Summary output

`to_summary_dict()` should emit:

```text
calibration_config
fl_models
diagnostics
region_calibration
background_model
boundary_local
boundary_sweep
strand_channels
prior_table
n_multi_loci
```

`region_calibration` summary should include:

```text
n_regions
rho_off
kappa_d
capture_enrichment_target
n_passes
converged
state_mass: sums/means of each state probability
p_expressed/p_captured summaries
mu_gdna/upper_gdna/rna_lower summaries
A_r/gamma_r summaries
flag histogram
pass_diagnostics
```

`boundary_local` and `boundary_sweep` summaries can be compact: count, mean/max of predicted gDNA, number of sparse/no-evidence/swept regions, and transfer-weight summaries.

### 3.3 No compatibility properties

Do not add properties like `.fused_region_gdna` or `.region_exposure` that recreate old state. Tests should be updated to the new contract.

---

## 4. Prior Assembly Cutover

File: `src/rigel/calibration/prior.py`

### 4.1 Inputs

`assemble_priors()` should consume:

```python
rc = calibration.region_calibration
region_mean = rc.mu_gdna
region_upper = rc.upper_gdna
exposure = rc.A_r
```

Remove the requirement for `calibration.fused_region_gdna`.

### 4.2 Count allocation

Keep the existing geometry allocation machinery, but pass `rc.mu_gdna` and `rc.upper_gdna` instead of fused evidence:

```python
unallocated_expected, unallocated_upper = _allocate_counts_by_geometry(
    region_arrays=region_arrays,
    multi_loci=multi_loci,
    fused_mean=np.asarray(rc.mu_gdna, dtype=np.float64),
    fused_upper=np.asarray(rc.upper_gdna, dtype=np.float64),
    ...
)
```

Rename internal variable names from `fused_*` to `region_*` during cleanup to avoid preserving old mental model.

### 4.3 Denominator/exposure

Use `RegionCalibration.A_r` directly in `bp_weighted_mean_exposure_over_blocks()`. That helper already accepts any object with an `A_r` attribute, so either:

```python
exposure=calibration.region_calibration
```

or add a clearer wrapper if needed. Do not rebuild `RegionExposure`.

### 4.4 Native gDNA eligibility

Change:

```python
evidence_enable = gdna_upper > _GDNA_ENABLE_UPPER_EPS
enable_gdna = candidate_enable & evidence_enable
```

to:

```python
enable_gdna = candidate_enable.astype(np.uint8)
```

Reason: calibration posterior mass should not act as an on/off switch. Native eligibility remains a technical unit-likelihood requirement.

### 4.5 PriorTable naming

Keep `PriorTable.enable_gdna` if native EM requires that name. Document it as native technical eligibility, not calibration state. If renaming is too invasive for Phase IV, defer the rename but update comments/tests.

---

## 5. Exposure Module Cutover

File: `src/rigel/calibration/exposure.py`

`RegionExposure` is old state. Options:

1. Delete `exposure.py` if no tests/imports need it after prior assembly moves to `RegionCalibration.A_r`.
2. If deletion is too disruptive in one pass, stop importing it from the live calibration path and mark it removed/deprecated in tests. Do not keep it on `CalibrationResult`.

Preferred Phase IV action: remove live imports and update tests first. Delete the file only after grep confirms no production imports remain.

---

## 6. Pipeline Cutover

File: `src/rigel/pipeline.py`

### 6.1 quant_from_buffer

Replace the old precondition:

```python
if calibration.fused_region_gdna is None: ...
```

with:

```python
if calibration.region_calibration is None: ...
```

The result should say `RegionCalibration`, not fused evidence.

### 6.2 Prior table replacement

`_replace(calibration, prior_table=prior_table, n_multi_loci=len(multi_loci))` remains valid if `CalibrationResult` is still frozen/dataclass-compatible.

### 6.3 EM call

Continue passing `prior_table.enable_gdna` to native EM because the native solver needs it. Its semantics change to technical eligibility only.

### 6.4 Logging

Change calibration logs from old `v6 quality` plus density/fusion summaries to:

```text
[CAL] region states: background=..., captured=..., expressed=...
[CAL] rho_off=... kappa_d=... passes=... converged=...
[CAL] A_r mean/p99/max=...
```

---

## 7. Config And CLI Cutover

Files: `src/rigel/config.py`, `src/rigel/cli.py`

### 7.1 Remove old density-exposure knob

Remove or repurpose:

```python
density_max_exposure
```

Do not clip `A_r` in Phase IV. If later benchmarks need clipping, it should be a new explicit v6 knob with justification.

### 7.2 Keep confidence knob but update semantics

Keep `gdna_density_confidence` for now to avoid unnecessary CLI churn, but update documentation:

```text
confidence for regional gDNA upper bounds in boundary/local/sweep predictions and RegionCalibration
```

### 7.3 Add calibration pass controls

Add:

```python
background_trim_fraction: float = 0.01
max_calibration_passes: int = 5
```

Potentially keep `density_min_eff_length` but update text to describe background seed effective-length eligibility, not old density model fitting.

### 7.4 No `use_v6_calibration`

Do not add this field. Remove it from docs if present.

---

## 8. Test Migration Plan

### 8.1 Update existing tests

Expected updates:

- `tests/test_calibrate.py`: expect `region_calibration`, `background_model`, boundary summaries, and no `fused_region_gdna` / `region_exposure` assertions.
- `tests/test_calibration_prior.py`: construct `CalibrationResult` with `RegionCalibration`; validate count allocation from `mu_gdna`/`upper_gdna`; validate `enable_gdna` no longer depends on positive upper evidence.
- `tests/test_pipeline_wiring.py`: mock calibration with `region_calibration`, not old fused/exposure fields.
- `tests/test_calibration_result.py`: update dataclass and summary expectations.
- Golden tests: update intentionally after smoke tests demonstrate deterministic output.

### 8.2 Add new tests

Add focused Phase IV tests:

1. `tests/test_calibration_cutover.py`
   - `calibrate()` returns `CalibrationResult.region_calibration`.
   - Old fields are absent.
   - Summary includes state mass and pass diagnostics.

2. `tests/test_prior_region_calibration.py`
   - `assemble_priors()` consumes `RegionCalibration.mu_gdna`, `upper_gdna`, and `A_r`.
   - gDNA effective length uses bp-weighted `A_r`.
   - `enable_gdna` is true for finite unspliced gDNA likelihood even when `upper_gdna == 0`.

3. `tests/test_pipeline_cutover.py`
   - `quant_from_buffer()` requires `region_calibration`.
   - Pipeline passes `RegionCalibration`-derived priors into EM.

### 8.3 Delete or rewrite old tests

Tests whose only purpose is old density/fusion behavior should be removed or rewritten. Do not preserve failing old-path tests by adding compatibility shims.

Likely affected:

```text
tests/test_density_model.py
tests/test_density_global.py
tests/test_density_observation.py
tests/test_exposure.py
tests/test_region_exposure_from_density.py
tests/test_calibration_integration.py
tests/test_calibration_prior.py
tests/test_calibrate.py
```

Review individually; keep low-level math tests only if the module remains used.

---

## 9. Implementation Order

### Step 1 - Result contract

- Rewrite `CalibrationResult` and `build_calibration_result()` around `RegionCalibration`.
- Add compact summary helpers for `RegionCalibration`, `BackgroundModel`, boundary local/sweep, and strand channels.
- Update direct result tests.

### Step 2 - Orchestrator cutover

- Replace old density/fusion/exposure path with Phase I/II/III flow.
- Ensure stranded/unstranded handling is explicit.
- Remove old imports from `_orchestrator.py`.
- Run `tests/test_calibrate.py` after updating expectations.

### Step 3 - Prior assembly cutover

- Update `assemble_priors()` to consume `calibration.region_calibration`.
- Change `enable_gdna` to technical native eligibility only.
- Update prior tests and Bayesian prior acceptance tests.

### Step 4 - Pipeline cutover

- Replace fused evidence preconditions.
- Update pipeline wiring tests.
- Run pipeline smoke tests.

### Step 5 - Config/CLI/docs cleanup

- Remove `density_max_exposure` from config/CLI.
- Add `background_trim_fraction` and `max_calibration_passes`.
- Update CLI help text and docs.

### Step 6 - Old-state cleanup

- Grep for old fields/modules.
- Remove live imports and stale tests.
- Delete old modules only after production references are gone, or leave them temporarily only if tests still exercise reusable low-level math. They must not be part of live calibration state.

---

## 10. Validation Plan

Run in escalating order:

```bash
conda activate rigel

ruff check src/ tests/

pytest tests/test_latent_states.py tests/test_calibration_iteration.py -v
pytest tests/test_calibration_cutover.py tests/test_prior_region_calibration.py -v
pytest tests/test_calibrate.py tests/test_calibration_prior.py tests/test_pipeline_wiring.py -v
pytest tests/test_pipeline_smoke.py tests/test_golden_output.py -v
```

Then broader validation:

```bash
pytest tests/ -v
```

Expected outcome: old golden tests may fail until intentional golden updates. Treat those as acceptance-review checkpoints, not as reasons to preserve the old path.

---

## 11. Open Issues To Resolve During Implementation

### R1. Exact strand tensor remains deferred

The clean cutover can ship with the conservative strand summary term, but Phase IV should not make it harder to pass raw folded strand counts later. Keep function signatures extensible.

### R2. `kappa_d` live refit needs raw counts

Phase III carries `kappa_d` rather than refitting it. Phase IV should decide whether to pass the needed folded counts to `run_calibration_iteration()` now or keep kappa fixed until benchmark evidence justifies the additional loop coupling.

Initial recommendation: keep fixed for cutover, log it, and benchmark.

### R3. Config naming

`gdna_density_confidence` is still an old-ish name. Keep it for one cutover to avoid CLI churn, but update semantics. Rename later only if docs/users are confused.

### R4. Old module deletion blast radius

`density_model.py`, `integration.py`, and `exposure.py` have many tests. Deleting them in the same PR may create noisy mechanical churn. The hard requirement is no live production dependency and no old state on `CalibrationResult`. Physical deletion can be a cleanup checkpoint after the cutover tests are green.

---

## 12. Acceptance Criteria

Phase IV is complete when:

- `calibrate()` returns `CalibrationResult.region_calibration` and no old density/fusion/exposure result fields.
- `quant_from_buffer()` only requires `region_calibration` for calibration priors.
- `assemble_priors()` uses `mu_gdna`, `upper_gdna`, and `A_r` from `RegionCalibration`.
- Calibration no longer gates the gDNA component off using posterior evidence thresholds.
- No `use_v6_calibration` flag or equivalent exists.
- Focused calibration, prior, and pipeline tests pass under the new contract.
- Running log documents all deviations, deleted old assumptions, and any remaining residual issues.
