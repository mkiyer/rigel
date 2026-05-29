# PR 05 v2 - EM Exposure Normalization and Prior Wiring Cleanup

> Superseded by [pr05_impl_plan_v3.md](pr05_impl_plan_v3.md). Do not implement this v2 plan:
> it applies exposure only to the gDNA component and is therefore asymmetric for capture libraries.

## Status

Builds on [pr04_impl_plan_v2.md](pr04_impl_plan_v2.md) and replaces the implementation details in
[pr05_downstream_em_exposure.md](pr05_downstream_em_exposure.md). The roadmap goal is unchanged:
make high regional gDNA exposure increase gDNA competitiveness in locus EM.

PR 04 intentionally produced `region_exposure.omega` but kept EM exposure-neutral. PR 05 is the
first PR that consumes `omega` downstream.

## 1. Goal

Implement a clean, inspectable path from calibration output to native EM inputs:

```text
RegionUnsplicedMass.gdna_mass/rna_mass  ->  additive Dirichlet mass split
RegionExposure.omega                    ->  gDNA denominator normalization
```

The important cleanup is to stop treating those as one mixed wiring concern. Regional gDNA/RNA mass
answers the Bayesian prior split question. Regional exposure answers the gDNA visibility denominator
question.

## 2. Core Decisions

### 2.1 Exposure Does Not Move Dirichlet Alpha Mass In PR 05

The additive EM prior split remains:

```text
alpha_gdna_add comes from M_r = RegionUnsplicedMass.gdna_mass
alpha_rna_add  comes from R_r = RegionUnsplicedMass.rna_mass
```

with the existing `p_unexpressed` soft gate, locus overlap projection, ESS cap, structural gDNA gate,
RNA call bias, and locus ESS shrinkage.

`omega` does not convert RNA mass into gDNA mass and does not change `alpha_gdna_add` or
`alpha_rna_add` in PR 05. It changes only the gDNA component denominator passed to native EM.

### 2.2 Exposure Is Denominator-Only

Native EM already uses the gDNA effective length as a denominator:

```text
log P(read | gDNA component) includes -log(gdna_eff_len)
```

Therefore PR 05 uses:

```text
gdna_eff_len_em = max(raw_gdna_eff_len_unweighted / gdna_exposure_factor, 1.0)
```

High exposure (`gdna_exposure_factor > 1`) decreases the denominator. Low exposure
(`gdna_exposure_factor < 1`) increases it.

Do not add a per-unit `+log(omega)` term in PR 05.

### 2.3 Locus Exposure Must Be FL-Window-Aware

A locus can contain many fine regions, and a small high-exposure subregion inside a long candidate
window must not dominate by simple region count or transcript span count. PR 05 should aggregate
`omega` over the same fragment-start opportunity windows used by `gdna_eff_len_for_loci()`.

The locus factor is:

```text
gdna_exposure_factor = exposed_start_opportunity / unweighted_start_opportunity
```

where exposure defaults to `1.0` for start positions not covered by the region table.

## 3. Current Post-PR04 State

Relevant production state after PR 04:

- `RegionCalibration.region_exposure` is a `RegionExposure` object with per-region `omega`.
- `RegionCalibration.A_r` is gone.
- `RegionExposure.from_density()` and `RegionExposure.uniform()` are gone.
- `assemble_priors()` computes the additive mass split from `RegionUnsplicedMass` and then sets:

```text
gdna_eff_len = gdna_eff_len_unweighted
gdna_em_exposure_weight = 1.0
```

- Output metadata still uses PR04 placeholder names:

```text
gdna_eff_len
gdna_eff_len_unweighted
gdna_eff_len_weight_ratio
gdna_em_exposure_weight
```

PR 05 replaces those names with sign-explicit names.

## 4. Target Contracts

### 4.1 PriorTable Fields

Update `PriorTable` in `src/rigel/calibration/prior.py`:

```python
@dataclass(frozen=True, slots=True)
class PriorTable:
    alpha_gdna_add: np.ndarray
    alpha_rna_add: np.ndarray

    # Native EM gDNA denominator.
    gdna_eff_len_em: np.ndarray

    # Denominator diagnostics.
    gdna_eff_len_unweighted: np.ndarray
    gdna_exposure_factor: np.ndarray
    gdna_eff_len_adjustment_ratio: np.ndarray

    enable_gdna: np.ndarray
    ... existing adaptive-prior diagnostics ...
```

Delete from the final PR05 contract:

```text
gdna_eff_len
gdna_em_exposure_weight
gdna_eff_len_weight_ratio
```

`gdna_eff_len_adjustment_ratio` is always computed as:

```text
gdna_eff_len_adjustment_ratio = gdna_eff_len_em / gdna_eff_len_unweighted
```

It is usually close to `1 / gdna_exposure_factor`, but it can differ when the post-adjustment
minimum floor of `1.0` is active.

### 4.2 Locus Metadata / Output Columns

Update locus metadata dictionaries in `src/rigel/pipeline.py` and user-facing loci output in
`src/rigel/estimator.py`.

New columns:

```text
gdna_eff_len_em
gdna_eff_len_em_per_bp
gdna_eff_len_unweighted
gdna_exposure_factor
gdna_eff_len_adjustment_ratio
```

Remove these output columns:

```text
gdna_eff_len
gdna_eff_len_per_bp
gdna_eff_len_weight_ratio
gdna_em_exposure_weight
```

`gdna_eff_len_em` is the value passed to native EM. `gdna_eff_len_unweighted` is the baseline
opportunity before exposure normalization. `gdna_exposure_factor` is the FL-window aggregation of
regional `omega`.

### 4.3 Summary JSON

Update `PriorTable.to_summary_dict()`:

```text
gdna_eff_len_em
gdna_eff_len_unweighted
gdna_exposure_factor
gdna_eff_len_adjustment_ratio
```

Update CLI summary in `src/rigel/cli.py` from the current `gdna_eff_len.value/per_bp` shape to a
sign-explicit block:

```python
"gdna_eff_len": {
    "em": _locus_series_summary("gdna_eff_len_em"),
    "em_per_bp": _locus_series_summary("gdna_eff_len_em_per_bp"),
    "unweighted": _locus_series_summary("gdna_eff_len_unweighted"),
    "exposure_factor": _locus_series_summary("gdna_exposure_factor"),
    "adjustment_ratio": _locus_series_summary("gdna_eff_len_adjustment_ratio"),
}
```

## 5. Geometry Algorithm

### 5.1 Low-Level Helper

Add FL-window-aware helpers in `src/rigel/calibration/_exposure.py`.

Recommended public helper:

```python
@dataclass(frozen=True, slots=True)
class GdnaExposureOpportunity:
    raw_unweighted: float
    raw_exposed: float
    exposure_factor: float


def gdna_exposure_opportunity_for_loci(
    loci: tuple | list,
    ref_lengths: Mapping[str | int, int] | Sequence[int],
    fl: FragmentLengthModel,
    region_arrays: RegionArrays,
    exposure: RegionExposure,
    *,
    default_omega: float = 1.0,
    min_factor: float = 1.0e-6,
    max_factor: float = 1.0e6,
) -> GdnaExposureOpportunity:
    ...
```

A convenience wrapper may also be added:

```python
def gdna_exposure_factor_for_loci(...) -> float:
    return gdna_exposure_opportunity_for_loci(...).exposure_factor
```

`RegionArrays` should be imported back into `_exposure.py` for this helper.

### 5.2 Exact Start-Window Integration

For each positive fragment length `ell` with probability `pmf[ell]`:

1. For each locus block `(ref_id, start, end)`, compute the same valid start window used by
   `gdna_eff_len_for_loci()`:

```text
valid_hi = ref_len - ell + 1
lo = max(start - ell + 1, 0)
hi = min(end, valid_hi)
```

Skip windows with `hi <= lo`.

2. Merge overlapping start windows per reference so multimember loci do not double count start
   positions.

3. Add unweighted opportunity:

```text
raw_unweighted += pmf[ell] * sum(width(merged_windows))
```

4. Add exposed opportunity:

```text
raw_exposed += pmf[ell] * sum(integral_omega_over_start_window(ref_id, lo, hi))
```

5. Compute:

```text
if raw_unweighted > 0:
    exposure_factor = raw_exposed / raw_unweighted
else:
    exposure_factor = 1.0
```

Then clip only as a numeric guard:

```text
exposure_factor = clip(exposure_factor, min_factor, max_factor)
```

### 5.3 Interval Exposure Integral

For an interval `[lo, hi)` on a reference, integrate exposure over start-coordinate positions:

```text
integral = uncovered_bp * default_omega + sum(overlap_bp_with_region_r * omega_r)
```

Implementation detail:

- Use `region_arrays.ref_offsets` to select regions on the reference.
- Use `np.searchsorted(region_arrays.end[lo_ref:hi_ref], lo, side="right")` and
  `np.searchsorted(region_arrays.start[lo_ref:hi_ref], hi, side="left")` to find overlapping
  regions.
- Track `covered_bp` from overlaps and add `(interval_width - covered_bp) * default_omega` for gaps.
- Assume fine regions are non-overlapping within a reference, matching the index contract.
- Validate `exposure.omega.shape == region_arrays.start.shape`.
- Validate all `omega` values are finite and positive.

This default-omega behavior is important because FL start windows can extend outside annotated fine
regions near locus boundaries.

### 5.4 Refactor Existing Unweighted Helper Carefully

`gdna_eff_len_for_loci()` must keep its existing public behavior, including the default
`min_value=1.0`.

For implementation, either:

1. Keep `gdna_eff_len_for_loci()` as-is and implement the exact exposed helper separately, or
2. Extract a shared internal opportunity-window iterator and make both helpers call it.

Preferred cleanup is option 2 if it stays readable:

```python
def _iter_gdna_start_windows(loci, ref_lengths, fl):
    yield ell_i, pmf_weight, ref_key, merged_windows
```

Do not reintroduce the old bp-weighted exposure helper.

## 6. Denominator Assembly

Add a focused builder in `src/rigel/calibration/prior.py`:

```python
@dataclass(frozen=True, slots=True)
class LocusGdnaDenominator:
    gdna_eff_len_em: np.ndarray
    gdna_eff_len_unweighted: np.ndarray
    gdna_exposure_factor: np.ndarray
    gdna_eff_len_adjustment_ratio: np.ndarray


def compute_locus_gdna_denominators(
    *,
    multi_loci: list[MultiLocus],
    index: TranscriptIndex,
    gdna_fl: FragmentLengthModel,
    region_arrays: RegionArrays,
    region_exposure: RegionExposure,
) -> LocusGdnaDenominator:
    ...
```

For each multi-locus:

```text
opp = gdna_exposure_opportunity_for_loci(...)
raw_unweighted = opp.raw_unweighted
omega_locus = opp.exposure_factor

unweighted = max(raw_unweighted, 1.0)
em = max(raw_unweighted / omega_locus, 1.0)
ratio = em / unweighted
```

Fill arrays by `multi_locus_id`. Validate contiguous locus IDs as `compute_adaptive_prior()` does.

Important: use the raw unweighted opportunity for the exposure division, then apply the native EM
floor. Do not divide the already-floored unweighted denominator. This preserves the intended order:
physical opportunity first, exposure normalization second, numeric floor last.

## 7. Prior Wiring Cleanup

Refactor `assemble_priors()` into three readable phases:

```text
1. Build common region arrays and native gDNA eligibility.
2. Build additive prior mass split from RegionUnsplicedMass.
3. Build exposure-adjusted gDNA denominators from RegionExposure.
```

Suggested shape:

```python
def assemble_priors(...):
    region_arrays = RegionArrays.from_region_df(...)
    has_gdna_candidate = ...

    adaptive = compute_adaptive_prior(...)
    denominators = compute_locus_gdna_denominators(...)

    return PriorTable(
        alpha_gdna_add=adaptive.alpha_gdna_add,
        alpha_rna_add=adaptive.alpha_rna_add,
        gdna_eff_len_em=denominators.gdna_eff_len_em,
        gdna_eff_len_unweighted=denominators.gdna_eff_len_unweighted,
        gdna_exposure_factor=denominators.gdna_exposure_factor,
        gdna_eff_len_adjustment_ratio=denominators.gdna_eff_len_adjustment_ratio,
        ...
    )
```

No exposure factor should be passed into `compute_adaptive_prior()` in PR 05.

## 8. Pipeline And Estimator Wiring

### 8.1 Pipeline

Update `_run_locus_em_partitioned()` signature:

```python
def _run_locus_em_partitioned(
    ...,
    gdna_eff_len_em: np.ndarray,
    ...,
    gdna_eff_len_unweighted: np.ndarray | None = None,
    gdna_exposure_factor: np.ndarray | None = None,
    gdna_eff_len_adjustment_ratio: np.ndarray | None = None,
) -> None:
```

Inside `_call_batch_em()`, pass the EM denominator to the native wrapper:

```python
estimator.run_batch_locus_em_partitioned(
    ...,
    gdna_eff_len=batch_gdna_eff_len_em,
    ...,
)
```

The native wrapper can keep the parameter name `gdna_eff_len` because that API represents the actual
native denominator. The pipeline-level name should be explicit.

Update `_build_locus_meta()` to write:

```text
gdna_eff_len_em
gdna_eff_len_em_per_bp
gdna_eff_len_unweighted
gdna_exposure_factor
gdna_eff_len_adjustment_ratio
```

### 8.2 Estimator Output

Update `AbundanceEstimator.get_loci_df()`:

- Update `cols`.
- Read new metadata keys.
- Remove old fallback aliases unless needed briefly within a single patch.
- Update docstring so `gdna_eff_len_em` is described as the denominator passed to EM.

### 8.3 CLI Summary

Update the summary JSON block in `src/rigel/cli.py` to use the new loci columns.

## 9. Native EM

No native code change is required.

Native EM already floors and logs the gDNA effective length it receives. PR 05 changes the Python
value passed into that parameter.

Because no native files are touched, no native rebuild is required for PR 05 unless later
implementation discovers a native test fixture that needs C++ changes.

## 10. Tests

### 10.1 New Geometry Tests

Add or extend tests for `_exposure.py`, probably in `tests/test_exposure.py` or a new focused
`tests/test_gdna_exposure_denominator.py`.

Required cases:

1. Uniform `omega == 1` gives `gdna_exposure_factor == 1`.
2. Whole-window `omega == 50` gives factor `50` and EM denominator `unweighted / 50`.
3. Whole-window `omega == 0.1` gives factor `0.1` and EM denominator larger than unweighted.
4. Mixed-exposure single-locus case with delta FL length 1:

```text
locus = [0, 1000)
regions = [0, 100) omega=10, [100, 1000) omega=1
expected factor = (100*10 + 900*1) / 1000 = 1.9
```

5. Gap handling:

```text
locus = [0, 1000)
only region [0, 100) has omega=10
uncovered [100, 1000) defaults to omega=1
expected factor = 1.9
```

6. Boundary expansion with delta FL length > 1 uses start windows, not just locus bp span.
7. Multi-block loci merge overlapping start windows and do not double count.
8. Empty/no-valid-window loci return factor `1.0` and denominator `1.0`.

### 10.2 Prior Assembly Tests

Update `tests/test_calibration_prior.py` and `tests/test_per_locus_gdna_mass.py`:

1. Existing mass-split tests still pass and assert alpha values are unchanged by exposure.
2. High exposure fixture now asserts:

```text
gdna_exposure_factor > 1
gdna_eff_len_em < gdna_eff_len_unweighted
gdna_eff_len_adjustment_ratio < 1
```

3. Low exposure fixture asserts:

```text
gdna_exposure_factor < 1
gdna_eff_len_em > gdna_eff_len_unweighted
gdna_eff_len_adjustment_ratio > 1
```

4. Unit test from roadmap:

```text
raw_unweighted = 100000
omega_locus = 50
gdna_eff_len_em = 2000, not 5000000
```

This can be tested through `compute_locus_gdna_denominators()` with a delta-FL/mock region setup.

5. Minimum floor test:

```text
raw_unweighted / omega_locus < 1
gdna_eff_len_em == 1
ratio may differ from 1 / omega_locus
```

### 10.3 Pipeline / Estimator Tests

Update `tests/test_pipeline_wiring.py`:

- `prior_table` fixture uses `gdna_eff_len_em`, `gdna_exposure_factor`, and
  `gdna_eff_len_adjustment_ratio`.
- `_run_locus_em_partitioned()` receives and records the EM denominator field.

Update estimator/locus output tests:

- New columns exist.
- Removed old columns do not exist.
- `gdna_eff_len_adjustment_ratio == gdna_eff_len_em / gdna_eff_len_unweighted`.

### 10.4 Native-Level Behavioral Smoke

Add a focused EM behavior test, likely in `tests/test_batch_em_impl.py` or a small prior/pipeline
smoke:

- Build the same ambiguous unspliced locus twice with identical log likelihoods and priors.
- Run once with uniform exposure-equivalent denominator.
- Run once with high exposure-equivalent denominator (`gdna_eff_len_em` smaller).
- Assert gDNA posterior/count increases in the high-exposure run.

This test verifies the sign without needing a full synthetic capture simulation.

### 10.5 Golden Outputs

Golden loci TSV headers will change. Update golden outputs only after the sign tests pass and the
new columns are inspected.

Expected header replacement:

```text
remove: gdna_eff_len, gdna_eff_len_per_bp, gdna_eff_len_weight_ratio, gdna_em_exposure_weight
add:    gdna_eff_len_em, gdna_eff_len_em_per_bp, gdna_exposure_factor,
        gdna_eff_len_adjustment_ratio
keep:   gdna_eff_len_unweighted
```

## 11. Validation Commands

Targeted:

```bash
conda activate rigel && pytest \
    tests/test_exposure.py \
    tests/test_calibration_prior.py \
    tests/test_per_locus_gdna_mass.py \
    tests/test_pipeline_wiring.py \
    tests/test_batch_em_impl.py \
    -v
```

Output contract:

```bash
conda activate rigel && pytest tests/test_golden_output.py -v
```

Full suite:

```bash
conda activate rigel && pytest tests/ -q
```

Lint touched Python files:

```bash
conda activate rigel && ruff check \
    src/rigel/calibration/_exposure.py \
    src/rigel/calibration/prior.py \
    src/rigel/pipeline.py \
    src/rigel/estimator.py \
    src/rigel/cli.py \
    tests/test_exposure.py \
    tests/test_calibration_prior.py \
    tests/test_per_locus_gdna_mass.py \
    tests/test_pipeline_wiring.py \
    tests/test_batch_em_impl.py
```

## 12. Implementation Order

1. Add exact FL-window exposure opportunity helper and tests.
2. Add `LocusGdnaDenominator` and `compute_locus_gdna_denominators()` in `prior.py`.
3. Update `PriorTable` fields and `to_summary_dict()`.
4. Refactor `assemble_priors()` into mass-split and denominator phases.
5. Update prior tests for high/low exposure sign and mass-split invariance.
6. Update `_run_locus_em_partitioned()` and `quant_from_buffer()` wiring to pass
   `gdna_eff_len_em`.
7. Update `AbundanceEstimator.get_loci_df()` columns and docstring.
8. Update CLI summary block.
9. Add native-level sign smoke test.
10. Update golden outputs intentionally.
11. Run targeted tests, golden tests, then full suite.
12. Final grep:

```bash
grep -R --exclude-dir='__pycache__' -n -E \
  'gdna_em_exposure_weight|gdna_eff_len_weight_ratio|unweighted \*|\* .*exposure' \
  src/rigel tests
```

Expected: no production path multiplies gDNA effective length by exposure.

## 13. Open Issues To Confirm

### Issue 1 - Should `omega` Ever Affect Additive Alpha Mass?

Recommendation: no for PR 05. `M_r/R_r` is the mass split evidence. `omega` is sampling visibility
and should affect the gDNA likelihood denominator only. This matches the roadmap and avoids double
counting exposure.

Clarification needed only if you want PR 05 to change `alpha_gdna_add`/`alpha_rna_add` directly.

### Issue 2 - Start-Coordinate Exposure Versus Fragment-Interval Exposure

Recommendation: aggregate `omega` over fragment start windows because that exactly mirrors
`gdna_eff_len_for_loci()` and the PR05 roadmap says to use the same opportunity windows.

Alternative: weight each candidate fragment by average exposure across its full interval. That is a
different likelihood model and starts to resemble the future local-exposure likelihood. It should not
be mixed into this denominator-only PR unless explicitly chosen.

### Issue 3 - Exposure For Start Windows Outside Region Table

Recommendation: uncovered start positions use `default_omega = 1.0`.

Reason: the region table may not cover every expanded gDNA start window near locus boundaries. Treat
missing exposure as neutral rather than dropping opportunity or treating it as zero exposure.

### Issue 4 - Output Backward Compatibility

Recommendation: no backward compatibility. The roadmap says no compatibility is required, and the
new output names are much clearer. This will require intentional golden updates.

### Issue 5 - Performance On Mega-Loci

Recommendation: implement the exact helper first, with shared window iteration if practical. The
expected cost is roughly `n_loci * n_positive_fragment_lengths * n_blocks_per_locus`, similar to the
existing exact unweighted helper, plus interval-over-region lookups. If profiling shows this is a
hotspot, add caching or per-reference prefix integrals of `omega` in a later optimization PR.

Do not fall back to bp-weighted exposure in production PR 05 unless this becomes a measured blocker.

## 14. Done Means

- `alpha_gdna_add` and `alpha_rna_add` remain a mass-conserving projection of
  `RegionUnsplicedMass`.
- `gdna_exposure_factor` is computed from exact FL start-window opportunity, not a bp-weighted span
  average.
- High exposure decreases `gdna_eff_len_em`.
- Low exposure increases `gdna_eff_len_em`.
- Native EM receives `gdna_eff_len_em`.
- Locus output names make the sign inspectable.
- No production path multiplies gDNA effective length by exposure.
- Golden outputs are updated only after sign tests pass.
