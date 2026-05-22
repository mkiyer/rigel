# Exact Exposure-Weighted gDNA Effective Length Plan

Date: 2026-05-21

Status: implementation plan, not yet implemented

Related diagnostics:

- `docs/benchmarks/flg2_hotspot_diagnostics_2026-05-21.md`
- `docs/benchmarks/locus3_exposure_audit_2026-05-21.md`

## 1. Executive Summary

The FLG2 hotspot failure is not caused by EM being unable to solve the mixture.
The EM is being given a gDNA denominator that is still far too large for the
local capture opportunity.

For MultiLocus 3 in the VCAP run:

| Quantity | Value |
| --- | ---: |
| merged MultiLocus footprint | 353,317,609 bp |
| unweighted FL-marginal gDNA effective length | 353,762,179 bp |
| current scalar exposure multiplier | 0.261212 |
| gDNA effective length passed to EM | 92,407,095 bp |
| reduction versus unweighted | 3.83x |

The `92.4M` denominator is already exposure-weighted. The problem is that the
current weighting scheme only reduces the huge locus by a scalar average over the
whole mega-locus footprint. That average is not close to the local/capture
opportunity relevant to the FLG2 pileup.

There are two separable issues:

1. **Projection issue:** production uses a scalar footprint approximation even
   though the repository contains an exact FL-aware weighted gDNA effective
   length helper.
2. **Exposure-field issue:** the current `rho_ref` / cap-at-1 normalization
   causes too much exonic mega-locus territory to look fully exposed.

The first issue is a theory-consistency fix and should be implemented first.
However, for a 353 Mb mega-locus, exact FL-aware projection alone is unlikely to
shrink the denominator by orders of magnitude, because the fragment-length
expansion contributes only about 0.13% extra length. If exact projection leaves
MultiLocus 3 near `92M`, the next fix is exposure-field normalization and
saturation, not another EM heuristic.

## 2. What Is The Current Footprint Exposure Weight?

Current production code in `assemble_priors()` computes:

```text
L_g_raw = gdna_eff_len_for_loci(ml.loci, ref_lengths, gdna_fl)
W_g     = footprint_exposure_weight(ml.loci, regional_exposure)
L_g_em  = max(L_g_raw * W_g, 1.0)
eta_em  = eta_raw * L_g_em / L_g_raw
```

`footprint_exposure_weight(...)` is a bp-weighted mean of the per-region exposure
field over the unexpanded genomic footprint:

```text
W_g = integral_{x in MultiLocus footprint} A(x) dx
      / length(MultiLocus footprint)
```

Implementation properties:

- input blocks are half-open genomic intervals `(ref_id, start, end)`;
- overlapping blocks on the same reference are merged before integration;
- `RegionalGdnaExposure.weighted_length_on_ref(...)` computes
  `integral A(x) dx` over each merged block;
- the result is clipped to `[1e-4, 1]`;
- uniform exposure returns exactly `1`.

This rule is simple and consistent across component classes, but it is not the
same operation as the gDNA effective-length denominator used by EM.

## 3. Why Are We Computing It This Way Today?

This was an intentional simplification in the v4.3 denominator-only plan, not a
low-level typo.

History from the repository:

- `3aada0e` introduced the v3 regional exposure plan and the exact
  `weighted_gdna_eff_len_for_loci(...)` helper.
- `7a9caf9` implemented the v4.3 plan and changed production to the scalar
  footprint rule.
- `docs/calibration/gdna_exposure_plan_v4.3.md` explicitly says the footprint
  average is "not exact fragment-midpoint integration" and that the exact
  refinement was deferred until the simple model was tested.

So the current behavior is best classified as a **model simplification that is
now outside its valid regime**, not a coding accident. The simplification is
acceptable for small loci or slowly varying exposure fields. It is fragile for
mega-loci created by multimapper connectivity, especially in hybrid capture data
where local exposure can vary sharply across the component.

## 4. What Does `weighted_gdna_eff_len_for_loci(...)` Do Differently?

The unweighted gDNA effective length is not just the locus span. For a fragment
of length `ell`, a genomic start position contributes if the length-`ell`
fragment overlaps any locus interval. For a locus interval `[a, b)`, the valid
start window is approximately:

```text
[a - ell + 1, b)
```

after clipping to valid contig starts. The unweighted denominator averages the
number of valid starts over the gDNA fragment-length PMF.

The exact weighted helper preserves that geometry and weights the actual
fragment opportunity:

1. Iterate over fragment lengths with nonzero `h_gDNA(ell)`.
2. Expand every locus interval into its valid start window for that `ell`.
3. Merge overlapping start windows per reference to avoid double-counting.
4. Shift start windows to midpoint windows using `ell // 2`.
5. Integrate the exposure field over those midpoint windows:

   ```text
   L_g_weighted = sum_ell h(ell) * integral_{m in midpoint windows(ell)} A(m) dm
   ```

6. Return the weighted sum, floored at `min_value`.

The scalar approach instead computes:

```text
L_g_scalar = L_g_unweighted * mean_A(unexpanded footprint)
```

That is equivalent only when `A(x)` is effectively constant over the locus and
its fragment-length shoulders. The exact helper is the mathematically correct
projection of a spatial exposure field onto the same denominator geometry used
by gDNA likelihood normalization.

## 5. What Difference Should We Expect For MultiLocus 3?

For MultiLocus 3, the saved values imply:

```text
unweighted FL effective length - footprint span
= 353,762,179 - 353,317,609
= 444,570 bp
= 0.126% of footprint span
```

Because this locus is huge, exact FL-aware integration can differ from scalar
footprint integration mainly in the fragment-length shoulders and in midpoint
alignment near exposure discontinuities. Those differences are small relative to
353 Mb unless the exposure field itself is extremely sparse inside the footprint.

If we keep the current exposure field and only replace scalar projection with
exact projection, a rough bound is:

| Quantity | Value |
| --- | ---: |
| current scalar denominator | 92,407,095 bp |
| scalar-weighted unexpanded footprint | 92,290,968 bp |
| maximum possible added shoulder mass if shoulders have `A=1` | 444,570 bp |
| minimum possible added shoulder mass if shoulders are at floor `A=1e-4` | 44 bp |

So exact projection alone is not expected to turn `92M` into `kb` for this
particular mega-locus. It is still the correct first fix because the current
implementation violates the denominator geometry, but it may not be sufficient.

The likely larger failure is that the current exposure field makes a broad
fraction of the mega-locus look exposed. The VCAP run has:

```text
rho_ref = 3.609e-4
EXON-COMPOSITE rho_q50 = 6.098e-4
EXON-COMPOSITE rho_q95 = 6.281e-3
```

Since weights are capped as `A_r = min(rho_hat_r / rho_ref, 1)`, many exonic
regions saturate at `A=1`. This is why every annotated RNA component in locus 3
has `em_exposure_weight=1.0`, and why the gDNA denominator remains enormous.

## 6. Is This A Bug Or A Misconception?

It is a misconception exposed by a real failure mode.

The scalar footprint average was adopted to make a consistent denominator-only
model quickly after numerator weighting caused regressions. That was a reasonable
engineering move for the v4.3 checkpoint. But it conflates two different
objects:

- a component-level average exposure over a genomic footprint;
- the FL-marginal gDNA opportunity that EM actually normalizes against.

For small loci the difference is negligible. For mega-loci, the scalar average
also hides the deeper problem: the exposure field is not sparse enough. The
current implementation answers "what is the mean exposure over this entire
connected component?" The EM needs "what is the effective exposed gDNA
opportunity for fragments that can explain this component under the assay?"

## 7. Implementation Plan

### Phase 0: Add An Exact-Versus-Scalar Audit Hook

Before changing EM behavior, add a diagnostic path that computes both values
for every MultiLocus:

```text
gdna_eff_len_unweighted
gdna_eff_len_scalar
gdna_eff_len_exact
gdna_eff_len_scalar_weight
gdna_eff_len_exact_weight_ratio = gdna_eff_len_exact / gdna_eff_len_unweighted
gdna_eff_len_exact_over_scalar = gdna_eff_len_exact / gdna_eff_len_scalar
```

Persist these fields to `loci.feather` and, when `emit_locus_stats` is enabled,
to `locus_stats.feather` or a dedicated diagnostic table.

Expected result for uniform exposure:

```text
gdna_eff_len_exact == gdna_eff_len_unweighted
gdna_eff_len_scalar == gdna_eff_len_unweighted
```

Expected result for simple constant subunit exposure:

```text
gdna_eff_len_exact ~= gdna_eff_len_scalar
```

Expected result for sharp exposure transitions:

```text
gdna_eff_len_exact differs from scalar according to FL midpoint geometry
```

### Phase 1: Use Exact Weighted gDNA Effective Length In `assemble_priors()`

Change imports in `src/rigel/calibration/locus_prior.py`:

```python
from ._exposure import weighted_gdna_eff_len_for_loci
```

For each `MultiLocus`:

```python
unweighted_eff_len = gdna_eff_len_for_loci(...)
scalar_weight = footprint_exposure_weight(...)
scalar_eff_len = max(unweighted_eff_len * scalar_weight, 1.0)
exact_eff_len = weighted_gdna_eff_len_for_loci(
    ml.loci,
    ref_lengths_arr,
    gdna_fl,
    regional_exposure,
    min_value=1.0,
)
gdna_eff_len = exact_eff_len
gdna_prior_count_em = eta_g * gdna_eff_len / unweighted_eff_len
```

Keep `gdna_prior_count` as the canonical unweighted/global diagnostic. Keep the
current prior-density invariant:

```text
gdna_prior_count_em / gdna_eff_len
== gdna_prior_count / gdna_eff_len_unweighted
```

This preserves the current denominator-only contract: regional exposure changes
the gDNA opportunity, not the global expected density.

### Phase 2: Tests

Update or add tests in `tests/test_weighted_eff_len.py` and
`tests/test_assemble_priors.py`:

1. Uniform exposure remains bit-exact with unweighted.
2. Constant exposure weight gives exact ratio.
3. Sharp transition synthetic case proves exact projection differs from scalar
   when FL midpoint windows extend into a different exposure region.
4. `assemble_priors()` uses `weighted_gdna_eff_len_for_loci(...)`, not
   `unweighted * footprint_weight`, for the EM denominator.
5. `gdna_prior_count_em / gdna_eff_len` invariant holds.
6. Diagnostics persist scalar and exact values.

### Phase 3: Focused FLG2/VCAP Validation

Run a focused VCAP quant with exact weighted denominator:

```bash
conda activate rigel
rigel quant \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  -o /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exact_weighted_leff_v1 \
  --include-multimap \
  --emit-locus-stats \
  --annotated-bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exact_weighted_leff_v1/annotated.bam
```

Then regenerate:

- FLG2 hotspot confusion;
- full gDNA/RNA confusion matrix;
- top false-RNA windows;
- locus 3 denominator audit.

Primary success criteria:

- `gdna_eff_len_exact` is persisted and consumed by EM.
- Global gDNA-to-RNA false assignments improve versus
  `exon_strand_deconv_v1`.
- FLG2 window gDNA-to-RNA assignment rate drops materially from 71.8%.

Important interpretation criterion:

- If `gdna_eff_len_exact` remains close to `92M`, exact projection is correct
  but not the root cause of the remaining failure. Move immediately to Phase 4.

### Phase 4: Exposure-Field Normalization And Capping Sweep

If exact projection is not enough, test reference-scale alternatives for
converting `rho_hat` to `A_r`.

Current rule:

```text
rho_ref = global E-weighted Q95(rho_hat)
A_r = clip(rho_hat_r / rho_ref, floor, 1)
```

This saturates many exonic regions in VCAP. Explore:

1. Global Q99 / Q99.5 / Q99.9.
2. Class-aware reference scales, especially EXON-COMPOSITE-specific reference
   values.
3. Max normalizer as a diagnostic upper bound, not the first production choice:

   ```text
   rho_ref = max(rho_hat)
   A_r = max(rho_hat_r / rho_ref, floor)
   ```

4. Robust high-tail normalizers such as trimmed max or weighted Q99.9.
5. Optional soft cap instead of hard `min(..., 1)`, if the high tail is unstable.

The maximum normalizer is mathematically attractive because nothing saturates,
but it is likely too strict in real capture data: a few probe outliers, CNV
regions, or duplicate-heavy regions can make almost the entire genome near zero.
Use it as a bounding experiment. If max-normalized VCAP fixes FLG2 while
destroying normal loci, the production target is a robust high-tail normalizer,
not raw max.

### Phase 5: Decide Based On Denominator Decomposition

For MultiLocus 3, record the exact denominator by exposure bin:

```text
sum contribution where A in [1e-4, 1e-3)
sum contribution where A in [1e-3, 1e-2)
sum contribution where A in [1e-2, 1e-1)
sum contribution where A in [1e-1, 1]
sum contribution where A == 1
```

This answers the central question: is the giant denominator coming from many
moderately exposed regions or from saturation at `A=1`? The next implementation
should be chosen from that decomposition, not from aggregate confusion alone.

## 8. Verification Commands

Focused unit tests after Phase 1:

```bash
conda activate rigel && pytest \
  tests/test_weighted_eff_len.py \
  tests/test_assemble_priors.py \
  tests/test_bayesian_prior_acceptance.py \
  -v
```

Lint changed files:

```bash
conda activate rigel && ruff check \
  src/rigel/calibration/_exposure.py \
  src/rigel/calibration/locus_prior.py \
  tests/test_weighted_eff_len.py \
  tests/test_assemble_priors.py
```

Full benchmark validation should compare at least:

- current `exon_strand_deconv_v1`;
- exact weighted denominator;
- exact weighted denominator plus alternate `rho_ref` candidates.

## 9. Open Questions

1. Does exact FL-aware projection materially differ from scalar projection on
   VCAP mega-loci? Geometry suggests probably not for MultiLocus 3, but we need
   persisted exact diagnostics to confirm.
2. How much of MultiLocus 3's denominator comes from saturated `A=1` regions?
3. Should `rho_ref` be global, class-aware, or a robust high-tail value for
   hybrid capture?
4. Should exposure normalization use a hard cap at 1, or should the model carry
   relative exposure above 1 and normalize elsewhere? Current EM denominator
   expects `A <= 1`, so this would be a larger design change.

## 10. Recommendation

Implement exact weighted gDNA effective length first because it is the correct
projection and restores consistency with the original v3 theory. Do not expect
it alone to make MultiLocus 3 kb-scale. Treat its diagnostics as the decision
point:

- If exact `L_gDNA` drops sharply, validate and keep the fix.
- If exact `L_gDNA` remains near `92M`, proceed directly to exposure-field
  normalization and capping experiments, with max normalization included only as
  a bounding diagnostic.
