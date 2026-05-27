# PR 2: Empirical-Bayes Local Capture Exposure

## Goal

Replace the current `A_r = mu_gdna / off_target_expectation` ratio with capture exposure built from physical capture opportunity and learned local enrichment. The model should allow strongly captured probe neighborhoods to have high gDNA exposure, weakly captured neighborhoods to shrink toward the panel mean, and no-capture data to return `A_r = 1` exactly.

This PR depends on PR 1. Source-reliable gDNA mass and reliability weights must be available before local exposure is fit.

## Non-goals

- Do not use latent capture states as RNA/gDNA source labels.
- Do not learn exposure from raw local observed counts alone.
- Do not add fixed `A_r` caps, fixed source pseudocounts, or ratio guards that become hidden model assumptions.
- Do not solve unstranded capture-on source splitting. This PR can estimate exposure there, but source mass remains limited by available evidence.
- Do not make production calibration depend on `rigel.sim` internals without an explicit extraction or adapter.

## Files to edit

| File | Purpose of edit |
|---|---|
| `src/rigel/config.py` | Add optional capture-panel calibration configuration. |
| `src/rigel/cli.py` | Add CLI/config-file plumbing for the optional capture-panel fields. |
| `src/rigel/calibration/capture_exposure.py` | New EB exposure model and diagnostics. |
| `src/rigel/calibration/_exposure.py` | FL-aware capture-opportunity helpers. |
| `src/rigel/calibration/calibration_iteration.py` | Call the exposure fit and stop deriving `A_r` from `mu_gdna`. |
| `src/rigel/calibration/prior.py` | Continue consuming `region_calibration.A_r`, now backed by EB exposure, for gDNA EM effective length. |
| `src/rigel/calibration/_result.py` | Include compact capture-exposure diagnostics in `summary.json`. |
| `src/rigel/sim/capture.py` or new shared panel module | Reuse or extract probe parsing/overlap semantics for production opportunity calculation. |
| `scripts/benchmarking/runner.py` | Pass simulation capture probe files to `rigel quant` when a benchmark condition has a panel. |
| `tests/test_capture_exposure.py` | Unit tests for opportunity and EB shrinkage. |
| `tests/test_per_locus_gdna_mass.py` | Integration test that gDNA effective length uses EB `A_r`. |

## Configuration contract

Add optional fields to `CalibrationConfig`:

```python
capture_probes: str | Path | None = None
capture_probe_format: Literal["auto", "transcript", "transcript_tsv", "bed12"] = "auto"
capture_gdna_split_penalty: float = 0.2
capture_min_overlap: int = 1
```

Do not add `binding_per_base` as a model parameter for calibration exposure. The EB model learns enrichment strength. The panel geometry should contribute a normalized opportunity `c_r` in `[0, 1]`, not a fixed capture strength.

Validation:

- `capture_probes is None` disables panel opportunity and returns neutral exposure.
- path fields are converted to strings in `PipelineConfig.to_dict()`.
- split penalty is finite and non-negative.
- min overlap is an integer >= 1.

Add CLI flags under `rigel quant`:

```text
--capture-probes PATH
--capture-probe-format {auto,transcript,transcript_tsv,bed12}
--capture-gdna-split-penalty FLOAT
--capture-min-overlap INT
```

Add `_ParamSpec` mappings so YAML config files can set these under `calibration`.

## New module: `capture_exposure.py`

Add flags:

```python
CAPTURE_EXPOSURE_NO_PANEL = np.uint16(0x1)
CAPTURE_EXPOSURE_NO_OPPORTUNITY = np.uint16(0x2)
CAPTURE_EXPOSURE_WEAK_GLOBAL = np.uint16(0x4)
CAPTURE_EXPOSURE_WEAK_LOCAL = np.uint16(0x8)
CAPTURE_EXPOSURE_SOURCE_UNRELIABLE = np.uint16(0x10)
```

Add dataclasses:

```python
@dataclass(frozen=True, slots=True)
class CaptureOpportunity:
    contained: np.ndarray              # float32[R], in [0, 1]
    boundary_left: np.ndarray          # float32[R], in [0, 1]
    boundary_right: np.ndarray         # float32[R], in [0, 1]
    unit_id: np.ndarray                # int32[R], -1 when no unit
    unit_names: tuple[str, ...]
    unit_opportunity: np.ndarray       # float64[U]

@dataclass(frozen=True, slots=True)
class CaptureExposureFit:
    A_r: np.ndarray                    # float32[R]
    gamma_r: np.ndarray                # float32[R], captured-only opportunity multiplier for diagnostics
    region_theta: np.ndarray           # float32[R]
    region_opportunity: np.ndarray     # float32[R]
    unit_theta_post: np.ndarray        # float64[U]
    unit_lambda_post: np.ndarray       # float64[U]
    unit_X: np.ndarray                 # float64[U]
    unit_O: np.ndarray                 # float64[U]
    unit_reliability_mass: np.ndarray  # float64[U]
    global_theta_mean: float
    global_alpha: float
    global_beta: float
    flags: np.ndarray                  # uint16[R]
    unit_flags: np.ndarray             # uint16[U]
```

Add `to_summary_dict()` or a private summary helper reporting:

```text
global_theta_mean
global_alpha
global_beta
n_units
n_units_weak_local
n_regions_no_panel
n_regions_no_opportunity
A_r summary stats
region_opportunity summary stats
unit_lambda_post summary stats
unit_X and unit_O summary stats
```

## Opportunity calculation contract

Capture opportunity is geometric and FL-aware. For each region/compartment, compute:

```text
c_{r,c} = E_{fragment length} E_{legal start in compartment window} best_probe_overlap_fraction(start, length)
```

where `best_probe_overlap_fraction` is in `[0, 1]` and uses non-stacking semantics: overlapping probes do not add; take the best single probe group contribution. For BED12 split probes, use `capture_gdna_split_penalty` for genomic/pre-mRNA blocks exactly as simulation does.

No counts enter this calculation.

Implementation recipe:

1. Extract the reusable panel parser/landscape logic from `src/rigel/sim/capture.py` into a production-safe helper, for example `src/rigel/capture_panel.py`, or add a narrow adapter that imports only pure parsing/interval utilities.
2. Keep `CaptureSampler` simulation behavior unchanged; update it to use the shared parser if extraction is chosen.
3. In `_exposure.py`, add helpers:

```python
def build_capture_opportunity(
    *,
    region_arrays: RegionArrays,
    boundaries: BoundaryTable,
    fl_model: FragmentLengthModel,
    panel: CapturePanel | None,
) -> CaptureOpportunity: ...
```

4. Use the gDNA FL model for gDNA exposure opportunity. Do not use RNA FL for `A_r`.
5. For no panel, return arrays of zeros, `unit_id=-1`, `unit_names=()`, and set no-panel flags later.
6. For regions overlapping multiple probe groups, choose the dominant unit by total opportunity for the first implementation. Store total opportunity as the non-stacking best-overlap value. A later PR can support fractional multi-unit attribution if needed.

## EB model contract

For each region and compartment:

```text
e_{r,c} = rho_off * leff_{r,c}
c_{r,c} = capture opportunity in [0, 1]
D_{r,c} = source-reliable gDNA mean from PR 1
w_{r,c} = source reliability from PR 1
```

Aggregate to exposure unit `u`:

```text
O_u = sum_{r,c in u} w_{r,c} * e_{r,c} * c_{r,c}
X_u = sum_{r,c in u} w_{r,c} * max(D_{r,c} - e_{r,c}, 0)
```

Then:

```text
theta_u = lambda_u - 1
theta_u >= 0
X_u | theta_u ~ Poisson(O_u * theta_u)
theta_u ~ Gamma(a_global, b_global)
theta_u_post = (a_global + X_u) / (b_global + O_u)
A_r = 1 + sum_c c_{r,c} * theta_{unit(r)}_post
```

For the region-level `A_r`, use contained opportunity as the primary multiplier for EM gDNA effective length. Keep boundary opportunities in the fit, not in the region multiplier, unless `_exposure.py` already has a clear region-level boundary-aware exposure convention. Record the chosen convention in the docstring.

## Estimating the global prior

Use source-reliable units only:

```text
eligible_global = unit_O > 0 and unit_reliability_mass > 0
```

Compute panel mean:

```text
theta_panel = max(sum(unit_X) / sum(unit_O), 0)
```

When at least three eligible units exist, use method of moments on unit rates:

```text
theta_hat_u = unit_X / unit_O
poisson_var_u = theta_panel / unit_O
between_var = weighted_var(theta_hat_u, weights=unit_O)
tau2 = max(between_var - weighted_mean(poisson_var_u, weights=unit_O), 0)
```

If `theta_panel > 0` and `tau2 > 0`:

```text
a_global = theta_panel^2 / tau2
b_global = theta_panel / tau2
```

Otherwise use a broad exponential prior centered on the panel mean:

```text
a_global = 1
b_global = 1 / max(theta_panel, eps)
```

Set `CAPTURE_EXPOSURE_WEAK_GLOBAL` in this fallback path. `eps` is only a numerical guard such as `1e-12`; it is not an exposure pseudocount.

If there is no panel or no opportunity, return neutral exposure:

```text
A_r = 1
gamma_r = 1
region_theta = 0
```

and set `CAPTURE_EXPOSURE_NO_PANEL` or `CAPTURE_EXPOSURE_NO_OPPORTUNITY`.

## Integration into `calibration_iteration.py`

Current code derives:

```python
denominator = rho_off * contained_leff
exposure = _safe_exposure(mu_gdna, denominator)
captured_exposure = _safe_exposure(captured_mu, denominator)
```

Replace this with:

1. Build or accept `CaptureOpportunity` before/inside `calibration_e_step()`.
2. Call `fit_capture_exposure()` after `prior_mass` is built, because the fit needs source-reliable mass.
3. Set `A_r = capture_fit.A_r` and `gamma_r = capture_fit.gamma_r`.
4. Store the fit on `CalibrationStepResult` and final `RegionCalibration` as `capture_exposure`.
5. Keep `_safe_exposure()` only if another call site needs it; otherwise delete it in a cleanup commit.

Suggested signature changes:

```python
def calibration_e_step(..., capture_opportunity: CaptureOpportunity | None = None, ...) -> CalibrationStepResult:
```

```python
def run_calibration_iteration(..., capture_opportunity: CaptureOpportunity | None = None, ...) -> RegionCalibration:
```

`calibrate()` in `_orchestrator.py` should build `capture_opportunity` once after `observation` and `boundaries` are available, then pass it into the iteration loop.

## Integration into `prior.py`

`assemble_priors()` already uses:

```python
bp_weighted_mean_exposure_over_blocks(..., exposure=region_calibration)
```

Keep that API stable if possible. It should now read `region_calibration.A_r`, whose source is EB exposure. Add diagnostics to `PriorTable.to_summary_dict()` so reviewers can see:

```text
gdna_em_exposure_weight summary
capture_exposure_active boolean
capture_exposure_global_weak boolean
```

Do not multiply gDNA effective length by any latent-state source probability.

## Benchmark runner plumbing

The synthetic suite manifest already records capture configuration for each condition. Update `scripts/benchmarking/runner.py` so a condition with `capture_config.probes` passes:

```text
--capture-probes <probe path>
--capture-probe-format <format>
--capture-gdna-split-penalty <value>
--capture-min-overlap <value>
```

If the condition is capture-off, pass no capture flags. This is important: PR 2 should be disabled by absence of panel/opportunity, not by condition-name parsing.

## Tests

### `tests/test_capture_exposure.py`

Add tests for:

1. No panel: all opportunity arrays are zero, `A_r == 1`, `gamma_r == 1`, and no-panel flags are set.
2. Panel with no overlapping regions: `A_r == 1`, no-opportunity flags are set.
3. No source-reliable excess: positive opportunity but `D <= e` gives local `X_u == 0`; posterior shrinks to the global prior and does not inflate no-gDNA regions.
4. Strong local enrichment: one unit with large `X_u` and `O_u` has `lambda_u_post` near `1 + X_u/O_u`.
5. Weak local enrichment: tiny `O_u` stays close to `1 + a_global/b_global`.
6. Boundary evidence contributes to `unit_X` and `unit_O` when boundary opportunity and reliability are positive.
7. Overlapping probes do not stack: opportunity is based on the best single probe group.
8. Split BED12 probe uses `capture_gdna_split_penalty` for genomic opportunity.
9. Very large supported enrichment returns a large finite `A_r`; there is no arbitrary maximum cap.

### Existing tests to update

- Update any `RegionCalibration` or `CalibrationStepResult` constructors in tests to include `capture_exposure` if the dataclass gains a required field. Prefer defaulting it to `None` for compatibility where possible.
- Update summary tests to assert a compact capture-exposure block exists.
- Update CLI/config tests for the new quant flags.

## Validation commands

```bash
conda activate rigel && ruff check src/rigel/config.py src/rigel/cli.py src/rigel/calibration src/rigel/sim/capture.py scripts/benchmarking tests/test_capture_exposure.py tests/test_per_locus_gdna_mass.py
conda activate rigel && pytest tests/test_capture_exposure.py tests/test_per_locus_gdna_mass.py tests/test_calibration_iteration.py tests/test_cli.py -v
```

Run broader calibration and pipeline tests before merge:

```bash
conda activate rigel && pytest tests/test_calibrate.py tests/test_calibration_result.py tests/test_pipeline_wiring.py tests/test_golden_output.py -v
```

Golden outputs should be updated only if the changed diagnostics or quantification are intentional and explained in the PR.

## Benchmark gate

After PR 1 + PR 2, rerun the eight-condition hybrid-capture suite with capture probes passed to `rigel quant` for capture-on conditions.

Required outcomes:

- capture-off conditions: `A_r` summary should be exactly or effectively `1`; mRNA MARD does not regress beyond 5 percent relative;
- no-gDNA capture-on: low source-reliable excess should not produce large local exposure;
- high-gDNA stranded capture-on: learned local exposure rises in probe-overlapping units and gDNA fraction moves toward truth;
- high-gDNA unstranded capture-on: no material regression, but full repair remains deferred to FL-plus-boundary source splitting;
- summary diagnostics report local/global weak flags instead of hiding fallbacks.

## Review checklist

- [ ] `A_r` no longer comes from `mu_gdna / denominator`.
- [ ] No fixed exposure cap or count pseudocount was introduced.
- [ ] No panel/no opportunity returns neutral `A_r = 1`.
- [ ] Local enrichment shrinks to a learned global prior.
- [ ] Boundary evidence is included in the EB fit.
- [ ] `prior.py` consumes EB `A_r` only as exposure, not source mass.
- [ ] Benchmark runner passes capture probes for capture-on synthetic conditions.
