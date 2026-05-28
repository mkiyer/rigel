# PR 01 - Two-State Calibration Teardown

## Goal

Replace the current four-state `(expressed yes/no) x (captured yes/no)` tensor with a two-state
expression model. Capture is no longer a latent class. Regional exposure will be learned later as a
continuous factor from calibrated gDNA mass.

This PR should remove production capture-state behavior rather than hiding it behind compatibility
switches.

## Current State

`src/rigel/calibration/latent_states.py` defines:

```text
STATE_UNEXPRESSED_OFFTARGET
STATE_UNEXPRESSED_CAPTURE
STATE_EXPRESSED_CAPTURE
STATE_EXPRESSED_OFFTARGET
STATE_IS_CAPTURED
STATE_IS_EXPRESSED
```

`src/rigel/calibration/calibration_iteration.py` consumes those states to compute:

- `p_states` with four columns
- `p_captured`
- `mu_gdna = p_captured * captured_mu + (1 - p_captured) * off_target_mu`
- `gamma_r`
- `capture_enrichment_target`
- pass diagnostics such as `n_regions_captured`

This is the classification task we want to remove.

## New State Model

The new state constants are:

```text
STATE_UNEXPRESSED = 0
STATE_EXPRESSED = 1
N_STATES = 2
STATE_NAMES = ("unexpressed", "expressed")
```

The E-step should answer one question:

```text
Is this region expression-contaminated enough that it should not train the gDNA background/exposure
profile as cleanly as an unexpressed region?
```

The posterior output is:

```text
p_unexpressed = p_states[:, STATE_UNEXPRESSED]
p_expressed = p_states[:, STATE_EXPRESSED]
```

## What Changes

### `latent_states.py`

- Delete capture constants and `STATE_IS_CAPTURED`.
- Delete `build_logbf_capture`.
- Delete the capture-vs-off-target use of `build_logbf_gdna_density` as a state factor.
- Keep an expression factor that uses spliced evidence and strand-derived RNA lower bounds.
- Keep an annotation prior that protects intergenic and intron-only regions from spurious expression.
- Rename functions around expression if that makes the module clearer. This file should no longer
  mention capture.

Expected functions after the PR:

```text
build_expression_log_prior(...)
build_logbf_expression(...)
build_logbf_strand_expression(...)
build_state_log_tensor(...)
normalize_state_log_tensor(...)
build_state_tensor(...)
```

### `calibration_iteration.py`

- Rename the top-level dataclass from `RegionCalibration` only if it helps readability; otherwise
  keep the name and change its fields.
- Remove `gamma_r`, `p_captured`, and `capture_enrichment_target`.
- Keep `PriorMassDeconvolution` temporarily if PR 03 has not landed, but stop deriving `mu_gdna`
  from captured/off-target mixing.
- `calibration_m_step` should refit `rho0` from regions weighted by `p_unexpressed`, plus background
  seed masks and RNA-exclusion masks.
- Pass diagnostics should report expression quantities only:

```text
n_regions_expressed
n_regions_unexpressed
max_state_shift
rho0
relative_rho_shift
sum_log_evidence
```

### `adaptive_prior.py`

- Update entropy weighting for `N_STATES = 2`.
- Verify prior mass projection still uses `gdna_unspliced_mean` and `rna_unspliced_mean`, not capture
  states.

### `_result.py`

- Update summaries to remove capture fields.
- State mass summaries should enumerate `unexpressed` and `expressed` only.

## What Stays

The following modules should remain production code and should not be deleted in this PR:

- `strand_deconv.py`
- `boundary_model.py`
- `boundary_sweep.py`
- `density_observation.py`
- `background_model.py`, after its seed weighting is adjusted to use `p_unexpressed`

They are gDNA mass-estimation machinery, not capture-classification machinery.

## Acceptance Tests

- Rewrite `tests/test_calibration_iteration.py` around two-state rows and row-sum invariants.
- Update latent-state tests so entropy, normalization, and expression evidence use `N_STATES = 2`.
- Add a regression where strong spliced or strand-RNA evidence raises `p_expressed`.
- Add a regression where intergenic unspliced-only evidence stays high `p_unexpressed`.
- Remove tests that assert capture-state labels or `p_captured` fields.

## Done Means

- Production code contains no capture latent constants.
- Calibration still produces mass-conserving gDNA/RNA prior evidence.
- Summary JSON no longer reports capture state mass.
- The smallest relevant test command passes:

```bash
conda activate rigel && pytest tests/test_calibration_iteration.py tests/test_calibration_prior.py -v
```
