# gDNA Regional Exposure Class-Normalization Fix Plan

**Status**: implementation-ready  
**Date**: 2026-05-18  
**Applies to**: implemented regional exposure plan v3  
**Primary module**: `src/rigel/calibration/_regional_exposure.py`

## 1. Summary

The v3 regional exposure implementation is theoretically correct in shape: use
the same exposure field `A(x)` in the gDNA likelihood numerator and denominator.
The real VCaP hybrid-capture regression comes from the way the implementation
normalizes that exposure field.

Current code computes a per-class high-density reference, then collapses the
references into one global maximum:

```python
rho_ref = max(rho_ref_per_class.values())
log_rho_ref = np.log(rho_ref)
```

That compares intron and intergenic density estimates to the EXON-INTRON
boundary channel. In the VCaP sample, the class references differ by orders of
magnitude:

```text
EXON-INTRON rho_ref_class  = 1.928126e-01
INTRON rho_ref_class       = 3.792302e-04
INTERGENIC rho_ref_class   = 1.065998e-09
```

As a result, an intron at its own 95th percentile receives `A ~= 0.0061`
instead of `A = 1.0`, and many ordinary introns are pushed to the `1e-4` floor.
The EM then sees true unspliced gDNA fragments as locally implausible gDNA and
routes them to RNA, mostly nRNA.

The fix is to make exposure weights class-relative: each region class must use
its own `rho_ref_class` when converting `rho_hat` to `A_r`. This preserves the
existing model and adds no new user-facing parameter.

## 2. Goals

1. Normalize each region class against its own exposure reference.
2. Preserve the v3 numerator/denominator contract:
   - per-unit gDNA score gets `log A(midpoint)`;
   - per-locus gDNA effective length integrates the same `A(x)`.
3. Keep `gdna_prior_count` unchanged; `gdna_prior_count_regional` remains a
   diagnostic only.
4. Add a regression test that fails under global cross-class normalization.
5. Validate on the real VCaP 50/50 RNA/gDNA mix using the full truth-derived
   confusion matrix.

## 3. Non-Goals

1. Do not add a new CLI flag or config parameter.
2. Do not change the native EM solver or scoring ABI.
3. Do not change the regional effective-length denominator formula.
4. Do not tune the `LOG_A_FLOOR`, `REFERENCE_QUANTILE`, or signal attenuation
   constants in this patch.
5. Do not introduce BED/capture-target support in this patch.

This is a surgical implementation fix, not a redesign of the exposure model.

## 4. Root Cause In Code

The problematic section is in `RegionalGdnaExposure.build()`:

```python
rho_ref = max(rho_ref_per_class.values()) if rho_ref_per_class else 0.0
...
log_rho_ref = float(np.log(rho_ref))
for cname in _CLASS_ORDER:
    ...
    raw = np.where(rho_c > 0.0, np.log(rho_c) - log_rho_ref, LOG_A_FLOOR)
```

`rho_hat` values are estimated within three different evidence channels:

- `INTERGENIC`: contained intergenic fragments per contained exposure;
- `INTRON`: intron-contained evidence per contained exposure;
- `EXON-INTRON`: exon-boundary crossing flux per boundary exposure.

These densities are not on an exchangeable numerical scale. Taking the maximum
across classes lets the boundary-flux channel define the normalization for all
other classes.

## 5. Implementation Plan

### Step 1: Use Class-Specific References For Weights

Change only the weight-construction block in
`src/rigel/calibration/_regional_exposure.py`.

Keep the existing global `rho_ref` field as a run-level summary diagnostic, but
do not use it for per-region weights. Instead, move reference lookup inside the
per-class loop:

```python
# Run-level diagnostic only.
rho_ref = max(rho_ref_per_class.values()) if rho_ref_per_class else 0.0

log_weight = np.zeros(R, dtype=np.float64)
for cname in _CLASS_ORDER:
    mask = class_mask[cname]
    if not mask.any():
        continue
    signal = signal_per_class[cname]
    if signal <= 0.0:
        continue

    rho_ref_c = float(rho_ref_per_class[cname])
    if rho_ref_c <= 0.0:
        continue
    log_rho_ref_c = float(np.log(rho_ref_c))

    rho_c = rho_hat[mask]
    with np.errstate(invalid="ignore", divide="ignore"):
        raw = np.where(rho_c > 0.0, np.log(rho_c) - log_rho_ref_c, LOG_A_FLOOR)
    raw = np.minimum(raw, 0.0)
    log_weight[mask] = np.maximum(signal * raw, LOG_A_FLOOR)
```

The behavior by class becomes:

- a region at or above its class `REFERENCE_QUANTILE` has `A_r = 1`;
- lower-density regions are attenuated within that class only;
- classes with `signal <= 0` remain identity exposure;
- classes with invalid or zero references remain identity exposure.

### Step 2: Improve Summary Diagnostics

The existing `per_class` summary already records `rho_ref_class`. Add the
weight distribution after `log_weight` is constructed:

```python
"weight_q05": ...,
"weight_q50": ...,
"weight_q95": ...,
"n_at_floor": ...,
```

These fields should be exposure-weighted quantiles using `E[mask]`, matching the
existing density quantiles. They are diagnostic only and should not affect EM.

Expected VCaP sanity check after the fix:

```text
INTRON rho_ref_class remains 3.792302e-04
INTRON weight_q95 should be 1.0 or very close to 1.0
INTRON n_at_floor should drop substantially from the global-reference run
```

### Step 3: Leave Existing Pipeline Order Untouched

Do not move any pipeline calls. The current order is correct:

1. `assemble_priors(...)` computes weighted `gdna_eff_len` from the exposure;
2. `_apply_unit_gdna_weights(...)` adds `log A(midpoint)` to finite gDNA units;
3. `partition_and_free(...)` scatters the already-weighted data;
4. native EM uses the supplied `gdna_eff_len`.

This patch changes the exposure table only. The numerator and denominator will
automatically use the corrected class-relative field.

### Step 4: No Native Rebuild Required

This fix touches Python only. A native rebuild is not required unless a later
implementation changes `src/rigel/native/`.

## 6. Unit Tests

### Test 1: Cross-Class Reference Does Not Suppress Intron Q95

Add to `tests/test_regional_exposure.py`.

Construct a small mixed-class region set with:

- several INTRON regions with a clear high-density intron at the class q95;
- several EXON regions with a much larger EXON-INTRON q95;
- global densities and payload counts chosen so both classes have nonzero
  signal.

Assert:

```python
exp.mode == "regional"
exp.per_class["EXON-INTRON"]["rho_ref_class"] > 100 * exp.per_class["INTRON"]["rho_ref_class"]
intron_high_weight == pytest.approx(1.0)
```

This test should fail under the current global-reference implementation because
the high intron would be normalized by the exon reference and strongly
attenuated.

### Test 2: Class With Zero Signal Remains Identity

Extend or add a test where EXON-INTRON has strong signal but INTERGENIC has
`signal == 0`. Assert all intergenic weights remain `1.0`. This protects the
existing auto-uniform behavior.

### Test 3: Summary Contains Per-Class Weight Diagnostics

Assert `to_summary_dict()["per_class"][cname]` contains:

- `rho_ref_class`;
- `weight_q05`;
- `weight_q50`;
- `weight_q95`;
- `n_at_floor`.

### Test 4: Uniform Mode Parity

Keep the existing uniform/no-op tests. The fix must not perturb:

- `RegionalGdnaExposure.uniform()`;
- `RegionalGdnaExposure.build(..., enabled=False)`;
- `weighted_gdna_eff_len_for_loci(..., uniform_exposure)`.

## 7. Focused Test Commands

Run these after implementation:

```bash
conda activate rigel
pytest tests/test_regional_exposure.py tests/test_weighted_eff_len.py -v
pytest tests/test_pipeline_smoke.py tests/test_pipeline_integration_v6.py -v
```

Then run the full suite before calling the fix complete:

```bash
conda activate rigel && pytest tests/ -v
```

## 8. Synthetic Validation

The synthetic 24-condition run previously passed because it did not stress the
real-data class-scale separation strongly enough. Still, the fix must preserve
synthetic behavior.

Run the existing synthetic regional rerun/analyzer scripts:

```bash
conda activate rigel && python scripts/debug/run_synthetic24_regional_v3.py
conda activate rigel && python scripts/debug/analyze_synthetic24_regional_v3.py
```

Acceptance criteria:

1. All 24 conditions complete.
2. No material regression versus `regional_off` in strand-specific uniform
   scenarios.
3. Pool confusion matrix remains comparable to the prior v3 synthetic report.
4. Any drift is explained by explicit changes in regional weights, not by
   output instability or missing annotations.

## 9. VCaP Validation

This is the decisive gate for the fix.

### Inputs

```text
input BAM: /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam
index:     /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index
truth:     C6EL5ANXX = RNA, H7MFFDSXY = gDNA
```

### Runs

Keep the same settings as the failed comparison. Write to a new output directory
so the current evidence is preserved:

```bash
conda activate rigel && rigel quant \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  -o /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_auto_classnorm_fix \
  --regional-exposure auto \
  --include-multimap \
  --seed 20260518 \
  --threads 8 \
  --annotated-bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_auto_classnorm_fix/annotated.bam \
  --emit-locus-stats
```

Compare against the preserved same-version control:

```text
/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_off_v3_with_mm
```

Use the existing diagnostic scripts:

```bash
conda activate rigel && python scripts/debug/analyze_vcap_regional_v3_confusion.py \
  --before /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_off_v3_with_mm/annotated.bam \
  --after /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_auto_classnorm_fix/annotated.bam \
  --out-dir /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_classnorm_fix_confusion

conda activate rigel && python scripts/debug/diagnose_vcap_regional_transitions.py \
  --before /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_off_v3_with_mm/annotated.bam \
  --after /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_auto_classnorm_fix/annotated.bam \
  --after-loci /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_auto_classnorm_fix/loci.feather \
  --out-dir /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/regional_classnorm_fix_confusion
```

### Acceptance Criteria

The failed v3 regional-auto run produced:

```text
gDNA -> RNA before/control: 2,913,875 / 18,228,677 = 15.99%
gDNA -> RNA regional-auto:  7,265,279 / 18,228,677 = 39.86%
```

The fixed run should satisfy:

1. `gDNA -> RNA` must fall far below the failed regional-auto rate.
2. The dominant new transition `gDNA before -> nRNA after` must collapse from
   the failed value of `4,253,739`.
3. `RNA -> gDNA` should remain no worse than the same-version uniform control;
   if it improves, record that as a benefit but do not trade it for massive
   gDNA false negatives.
4. `summary.json["calibration"]["regional_exposure"]["per_class"]["INTRON"]`
   should report class-relative high-intron weights near identity.
5. `loci.feather.gdna_eff_len_weight_ratio` may remain small for genuinely
   depleted loci, but top false gDNA-to-RNA loci should no longer be dominated
   by ratios caused by cross-class EXON-INTRON scaling.

## 10. Failure Escalation If Class-Normalization Is Insufficient

If the VCaP confusion matrix remains poor after class-relative normalization,
do not tune constants immediately. First run a second diagnostic focused on
local numerator/denominator mismatch:

1. For newly false gDNA-to-RNA fragments, record `log A(midpoint)` by class.
2. For their loci, record `log(gdna_eff_len / gdna_eff_len_unweighted)`.
3. Compute the net exposure term:

   ```text
   delta = log A(midpoint) - log(gdna_eff_len / gdna_eff_len_unweighted)
   ```

4. If `delta` remains strongly negative in true gDNA fragments, the residual
   issue is local support mismatch, not cross-class scaling.

Only then consider a second patch, such as smoother regional weights or an
exposure floor learned from the locus opportunity distribution. That would be a
separate plan.

## 11. Rollback

Rollback is trivial: restore the old global `rho_ref` use inside the weight
loop. No persisted index format or output schema depends on the new behavior.

The added summary diagnostics are backward-compatible extra JSON fields.

## 12. Implementation Checklist

1. Patch `_regional_exposure.py` to use `rho_ref_per_class[cname]` inside the
   weight loop.
2. Add per-class weight distribution diagnostics.
3. Add the cross-class regression test.
4. Run focused tests.
5. Run full tests.
6. Rerun synthetic 24-condition validation.
7. Rerun VCaP class-normalized regional-auto quantification.
8. Publish the before/after confusion matrix and transition diagnostics.
