# EXON Strand Deconvolution Implementation Plan v1

**Date**: 2026-05-20
**Status**: planned
**Companion audit log**: `docs/calibration/exon_strand_deconv_audit_log_v1.html`

## 1. Goal

Use strand-specific deconvolution on fully contained EXON fragments when the library has identifiable strand signal, and reconcile that evidence with the existing EXON-INTRON boundary-flux estimator for exonic gDNA density.

The change should improve large or boundary-sparse EXON regions with unambiguous region strand without weakening the current boundary-flux behavior in unstranded libraries, ambiguous-strand regions, or samples where strand training is underpowered.

Important v1 limitation: `RegionStrand.AMBIG` rows cannot use the existing SAME/OPP orientation model because both Python and native orientation helpers deliberately classify AMBIG rows as uninformative. Dense antisense loci are therefore a known partial miss for v1. The implementation must report how much EXON opportunity gains contained-strand power, split by POS/NEG versus AMBIG, and a follow-up design should evaluate per-strand sub-regions or per-strand contained payloads for AMBIG rows.

## 2. Current State

The calibration path currently has three gDNA density channels:

| Channel | Evidence | Opportunity | Strand correction today |
| --- | --- | --- | --- |
| `INTERGENIC` | `per_region_counts[:, MASK_INTERGENIC]` | FL-weighted containment | No, unstranded by definition |
| `INTRON` | `intron_counts_by_orient` | FL-weighted containment | Yes |
| `EXON-INTRON` | `u_left_by_orient`, `u_right_by_orient` on eligible EXON sides | boundary-crossing exposure | Yes |

EXON-contained fragments are currently not used for gDNA density estimation. Exonic gDNA mass is projected from the `EXON-INTRON` density in:

- `estimate_global_gdna_fragments()` for whole-genome gDNA projection;
- `expected_gdna_count_global()` for the canonical global-only EM prior;
- `_boundary_term_prorated()` for locoregional diagnostic `n_gdna_exon_only`.

The strand math already exists and should be reused:

- `StrandSummary` in `src/rigel/calibration/_orient.py`;
- `strand_correction_usable()` in `src/rigel/calibration/density_global.py`;
- `_gdna_count_moment()` in `src/rigel/calibration/density_global.py`;
- orientation routing through `orient::classify()` in native code.

## 3. Design Principles

1. Boundary flux remains the conservative baseline. It is the only usable EXON evidence when strand contrast is unavailable.
2. EXON-contained evidence is a strand-conditioned gDNA proxy, not a raw count proxy. Never use uncorrected EXON-only counts as gDNA evidence.
3. Fuse evidence by statistical power in count/opportunity space. Do not replace boundary density with contained density, and do not average rates without denominators or overdispersion.
4. Preserve the canonical global-only prior contract. Local payload counts may affect diagnostics and the regional exposure field, but `expected_gdna_count_global()` must remain payload-free.
5. Do not add EXON-only counts to the gDNA fragment-length pool in this PR. `extract_gdna_counts()` should continue to exclude `MASK_EXON`.
6. Keep gene-level outputs downstream-only. All inference changes stay at transcript/region/locus level.

## 4. Evidence Model

For an unambiguous-strand EXON region and an identifiable strand model:

```text
same = fragments whose observed read-1 orientation matches the EXON region strand
opp  = fragments whose observed read-1 orientation opposes the EXON region strand
c    = signed_strand_contrast = 2 * p_r1_sense - 1

gDNA_count_moment = same + opp - (same - opp) / c
```

This is already implemented by `_gdna_count_moment()`.

Define a scalar strand reliability weight:

```text
strand_power = c * c, if strand_correction_usable(strand_summary)
strand_power = 0, otherwise
```

This is the first implementation of strand reliability for EXON-contained evidence, not a full inverse-variance weight:

- SS=0.5 gives `c=0`, so `strand_power=0` and the model falls back to boundary flux.
- Perfect stranded libraries give `abs(c)=1`, so `strand_power=1` and contained evidence can get full precision weight.
- Intermediate stranded libraries contribute smoothly and conservatively.

The critique of the first draft correctly noted that `c * c` alone is a heuristic. It damps deconvolution noise but does not account for Poisson/NB variance. v1 should therefore use `strand_power` as only one factor in a precision-like opportunity weight.

Use the existing NB kappa estimate as a precision cap. `estimate_kappa()` returns a Gamma-Poisson/NB alpha in count units. Whenever it is used as an opportunity length, convert it to a bp-equivalent beta:

```text
beta = alpha / max(rho_ref, rho_floor)
```

Then convert raw opportunity into a precision opportunity:

```text
precision_opportunity(E, beta) = E * beta / (E + beta)
```

This behaves like `E` when Poisson noise dominates and saturates at `beta` when overdispersion dominates. It prevents a huge EXON contained opportunity from automatically overwhelming a tighter boundary estimate.

The v1 fusion rule for exonic contained density is therefore:

```text
P_boundary  = precision_opportunity(E_boundary, beta_boundary)
P_contained = strand_power * precision_opportunity(E_contained, beta_contained)

scale_boundary  = P_boundary / E_boundary, if E_boundary > 0 else 0
scale_contained = P_contained / E_contained, if E_contained > 0 else 0

Y_exon = scale_boundary * Y_boundary + scale_contained * Y_contained
E_exon = P_boundary + P_contained
rho_exon = Y_exon / E_exon, if E_exon > 0 else fallback
```

where `Y` is a strand-corrected gDNA count moment and `E` is the matching opportunity. Scaling `Y` and `E` by the same channel factor preserves each channel's rate while giving the fusion rule a precision-like denominator.

For global densities, aggregate strand-corrected moments first and then clip the aggregate to nonnegative, matching the current `_compute_density()` behavior. For regional exposure and locoregional diagnostics, do not introduce new clip-before-shrink behavior for EXON-contained evidence. Use signed corrected moments through the EB shrinkage step and apply the nonnegative floor only to final rates/weights. This avoids the upward bias caused by `np.maximum(corrected, 0)` on many small rows.

## 5. Native Scanner Changes

### 5.1 Add Payload Field

Modify `src/rigel/native/calibration/accumulator.h`:

```cpp
std::vector<int64_t> exon_contained_counts_by_orient;  // shape (n_regions, orient::N)
```

Initialize it in the constructor with `n_regions * orient::N` zeros, merge it in `merge_from()`, and expose it in `src/rigel/native/bam_scanner.cpp` as a 2-D ndarray named `exon_contained_counts_by_orient`.

The v1 field should stay `int64_t` for consistency with the existing counter schema and to avoid overflow surprises in high-depth targeted libraries. The extra memory is real: `(n_regions, 3) * int64` is about 16 MB for a 700k-region index per worker. If this becomes a bottleneck, evaluate an `int32_t` payload migration as a separate ABI/schema change with an explicit maximum-count guarantee.

### 5.2 Counting Gate

In `CalibrationAccumulator::observe()`, add a subset counter in the per-region loop. The recommended gate is:

```text
flux_eligible == true                       # splice_type == SPLICE_UNSPLICED
obs_mask == mask::EXON                      # no intron/intergenic aligned-block hit
regions.type_mask(rid) has mask::EXON
frag_start >= regions.start(rid)
frag_end <= regions.end(rid)
```

Then increment:

```cpp
payload_.exon_contained_counts_by_orient[orient_idx]++;
```

Use exact half-open containment, not merely "not boundary crossing", because the denominator is `l_eff_contained()`. The intended convention is `[frag_start, frag_end)` contained in `[region_start, region_end)`, so `frag_end == regions.end(rid)` counts as contained. At `q=1`, that same fragment is not a right-boundary crossing under the existing predicate `frag_end >= region_end + q`, which preserves disjointness between contained and boundary evidence.

### 5.3 What Not To Change

Do not add these fragments to `fl_hist[MASK_EXON]` for gDNA FL training. The existing `fl_hist` already records global mask bins, and `extract_gdna_counts()` must continue to exclude `MASK_EXON` because EXON-only is dominated by mature RNA and can contain spliced genomic spans.

## 6. Python Payload and Array Changes

### 6.1 `CalibrationScanPayload`

Modify `src/rigel/calibration/scan_payload.py`:

- add dataclass field `exon_contained_counts_by_orient: np.ndarray`;
- validate dtype `int64` and shape `(n_regions, ORIENT_N)`;
- validate row sums are less than or equal to `per_region_counts[:, MASK_EXON]`;
- document that strict inequality is expected and correct, because `per_region_counts[:, MASK_EXON]` includes spliced and non-contained EXON-only fragments;
- reject only row sums greater than the EXON mask count, which would imply a native counter bug.

### 6.2 `PayloadArrays`

Modify `src/rigel/calibration/_arrays.py`:

- add sorted view `exon_contained_by_orient`;
- populate it using `payload.exon_contained_counts_by_orient[region_arrays.order]`;
- add a regression test where `region_df` order differs from `(ref_id, start)` order so the new payload view cannot silently drift from `RegionArrays`.

All tests and synthetic payload builders that instantiate `CalibrationScanPayload` need a zero-filled field.

## 7. Global Density Changes

### 7.1 Add EXON-CONTAINED Channel

Modify `src/rigel/calibration/density_global.py`:

- extend `DensityType` with `"EXON-CONTAINED"`;
- add `_channel_exon_contained(region_mask, counts_by_orient, leff, strands)`;
- call `_compute_density()` for EXON-contained rows;
- use `l_eff_contained(exon_spans, gdna_fl)` as the denominator;
- compute the full `leff` vector once per call and reuse it for INTERGENIC, INTRON, EXON-INTRON diagnostics, and EXON-CONTAINED so the new channel does not duplicate fragment-length cache work.

The channel must only use EXON rows. If strand correction is not usable, this raw channel may still be reported diagnostically, but it must not contribute to the composite exonic gDNA density.

### 7.2 Add EXON Composite Diagnostics

Add a small frozen dataclass, for example:

```python
@dataclass(frozen=True, slots=True)
class ExonCompositeDensity:
      rho: float
      boundary_rho: float
      contained_rho: float
      boundary_precision: float
      contained_precision: float
      strand_power: float
      contained_active: bool
      precision_model: str
```

Extend `GlobalDensityTable` with:

```python
exon_contained: GlobalGdnaDensity
exon_composite: ExonCompositeDensity
```

Build the composite as:

```text
beta_boundary = opportunity_beta(exon_intron.kappa, exon_intron.rho)
beta_contained = opportunity_beta(exon_contained.kappa, exon_contained.rho)

boundary_precision = precision_opportunity(exon_intron.eff_length_bp, beta_boundary)
contained_precision = (
    strand_power * precision_opportunity(exon_contained.eff_length_bp, beta_contained)
    if exon_contained.strand_active
    else 0
)

rho = (exon_intron.rho * boundary_precision + exon_contained.rho * contained_precision)
   / (boundary_precision + contained_precision)
```

Fallbacks:

- if `contained_precision == 0`, `rho = exon_intron.rho`;
- if both precision terms are zero, `rho = 0`;
- if boundary precision is zero but contained precision is positive, use contained density.

Add `EXON-CONTAINED` and `EXON-COMPOSITE` sections to `GlobalDensityTable.to_summary_dict()` while keeping the existing `EXON-INTRON` section for compatibility and diagnostics.

The composite density is for EXON-contained projection. It should not erase the boundary density: downstream boundary-crossing expected counts should continue to use `global_densities.exon_intron.rho`.

## 8. Regional Exposure Changes

Modify `src/rigel/calibration/_regional_exposure.py`.

For EXON rows, compute two evidence streams:

```text
E_boundary = eligible_sides * B_cross
Y_boundary = current strand-corrected boundary count moment

E_contained = l_eff_contained(exon_span, gdna_fl)
Y_contained = strand-corrected exon_contained_by_orient count moment
```

Use contained evidence only for rows with POS or NEG region strand and a usable `StrandSummary`. AMBIG rows remain boundary-only in v1 because the current orientation payload is UNINF for those rows.

Then set:

```text
P_boundary  = precision_opportunity(E_boundary, beta_boundary)
P_contained = strand_power * precision_opportunity(E_contained, beta_contained)

Y[is_exon] = (P_boundary / E_boundary) * Y_boundary
           + (P_contained / E_contained) * Y_contained
E[is_exon] = P_boundary + P_contained
```

Use zero scale factors when a raw opportunity is zero. Keep the current global shrinkage and `rho_ref` logic after `Y` and `E` are built, but do not clip EXON-contained row moments before shrinkage. If signed moments produce negative post-shrink rates, floor the final `rho_hat` before log-weight calculation and report the number of rows floored.

Add per-class summaries for:

- `EXON-INTRON` boundary-only diagnostics;
- `EXON-CONTAINED` strand-only diagnostics;
- `EXON-COMPOSITE` actual evidence used in the exposure field.

This is the main place where giant EXON regions benefit, but only when precision supports it. A giant contained opportunity should not dominate if the contained channel is overdispersed or strand contrast is weak.

## 9. Locus Prior Changes

### 9.1 Diagnostic Locoregional Estimate

Modify `src/rigel/calibration/locus_prior.py`.

Add `_exon_contained_term_prorated()` analogous to `_density_term_prorated()`, but using `payload_arrays.exon_contained_by_orient` and strand correction instead of a single full-count column.

For each locus:

1. Compute the existing boundary branch:
   - local boundary observed count;
   - local boundary opportunity;
   - boundary-only `rho_b_loco` shrunk to `global_densities.exon_intron.rho`.
2. Compute contained branch when strand is usable:
   - strand-corrected contained gDNA count moment inside the locus;
   - contained opportunity from `eff_clip_core` on EXON rows;
   - contained `rho_c_loco` shrunk to `global_densities.exon_contained.rho`.
3. Fuse the two local rates:

```text
w_b = precision_opportunity(local_boundary_opportunity, beta_boundary)
w_c = strand_power * precision_opportunity(local_contained_opportunity, beta_contained)
rho_ex_loco = weighted_mean(rho_b_loco, rho_c_loco, w_b, w_c)
```

Fallbacks:

- if `w_c == 0`, preserve current boundary behavior;
- if `w_b == 0` and `w_c > 0`, use contained strand evidence;
- if both are zero, use `global_densities.exon_composite.rho`.

Keep `n_gdna_boundary_observed` as the raw local observed boundary count. Compute `n_gdna_exon_only` from the fused exonic density:

```text
n_gdna_exon_only = rho_ex_loco * l_core_exon
```

Add diagnostics as a nested frozen dataclass, for example `ExonGdnaDiagnostics`, instead of widening the EM-facing fields one by one. Flatten selected fields into `per_locus_gdna_df` only if they are needed for QC. Candidate fields:

- `rho_exon_boundary`;
- `rho_exon_contained`;
- `rho_exon_composite`;
- `exon_boundary_precision`;
- `exon_contained_precision`;
- `n_exon_contained_observed`;
- `n_exon_contained_estimated`.

When using kappa in `shrink_to_loco()`, convert `GlobalGdnaDensity.kappa.value` from NB alpha to bp-equivalent opportunity beta before passing it as the shrinkage opportunity. This is the same units issue already handled in `_regional_exposure.py`.

### 9.2 Canonical Global-Only Prior

Keep `expected_gdna_count_global()` payload-free. Keep the boundary-crossing term tied to the boundary density and use the composite only for EXON-contained projection:

```text
boundary_crossing_expected = rho_exon_boundary * s_locus * B_cross
exon_contained_expected    = rho_exon_composite * L_exon_locus
```

This lets global EXON-contained strand evidence improve the canonical prior without leaking local payload counts into the helper.

This still intentionally changes global-only priors in stranded libraries through `exon_contained_expected`. Report old boundary-only versus new composite exonic projections in tests/diagnostics so baseline shifts are explainable.

### 9.3 Whole-Genome Projection

Update `estimate_global_gdna_fragments()` so EXON-contained projection uses `global_densities.exon_composite.rho` instead of `global_densities.exon_intron.rho`. Include a diagnostic delta from the previous boundary-only projection in summary or debug output so regression shifts are not silent.

## 10. Output and Diagnostics

Update summaries and diagnostics in a minimal, explicit way:

- `summary.json["calibration"]["global_densities"]` includes `EXON-CONTAINED` and `EXON-COMPOSITE`.
- `summary.json["calibration"]["regional_exposure"]["per_class"]` includes composite EXON evidence diagnostics, and `_CLASS_ORDER` is updated without dropping existing `INTERGENIC`, `INTRON`, or `EXON-INTRON` keys.
- `summary.json` includes the fraction of EXON opportunity with active contained power, split by unambiguous-strand rows and AMBIG rows.
- `per_locus_gdna_df` includes flattened local fusion diagnostics only if needed; otherwise keep the nested diagnostics inside `LocusGdnaEstimate`/`CalibrationResult` internals.
- Existing `loci.feather` columns can remain unchanged for v1 unless we decide to surface `rho_exon_composite` there. The EM-facing fields `gdna_prior_count`, `gdna_prior_count_em`, and `gdna_eff_len` already capture the changed prior and denominator behavior.
- Update debug/analysis scripts that assume the old three-key `global_densities` or `regional_exposure.per_class` layout, especially scripts under `scripts/debug/analyze_*` and benchmark report readers.

## 11. Tests

### 11.1 Native Accumulator and Payload

Update `tests/test_calibration_accumulator.py`:

- payload shape includes `exon_contained_counts_by_orient`;
- worker-merge equality compares the new array;
- payload validation rejects wrong shape/dtype;
- a fully contained EXON-only unspliced fragment increments the proper SAME/OPP orientation;
- a boundary-crossing EXON fragment increments boundary flux but not contained counts;
- a fragment with `frag_end == region.end` is counted as contained and not as a right-boundary crossing under the half-open convention;
- a spliced EXON fragment does not increment contained counts;
- row-sum less than or equal to `per_region_counts[:, MASK_EXON]` is validated, including a case where a spliced EXON-only fragment makes the inequality strict;
- sorted `PayloadArrays` keep `exon_contained_by_orient` aligned when `region_df` original order differs from sorted `(ref_id, start)` order.

After native edits, rebuild before running tests:

```bash
conda activate rigel && pip install --no-build-isolation -e .
```

### 11.2 Global Density

Update `tests/test_density_global.py`:

- EXON-contained pure RNA in perfect stranded data gives contained rho zero with a tight tolerance;
- EXON-contained pure gDNA in perfect stranded data matches uncorrected density;
- near-unstranded data reports contained diagnostics but composite falls back to boundary;
- boundary and contained fusion weights use kappa-capped precision opportunity plus strand power;
- ambiguous EXON rows do not contribute contained power;
- read-1 antisense mirror remains invariant.

### 11.3 Regional Exposure

Update `tests/test_regional_exposure.py`:

- large EXON with contained strand evidence gets higher evidence weight than sparse boundary sides only when precision supports it;
- boundary with strong evidence on few sides is not overwhelmed by sparse contained counts on a huge span;
- unstranded EXON-contained counts do not alter exposure relative to boundary-only;
- composite `Y/E` uses precision scaling and remains finite;
- signed contained moments avoid clip-before-shrink bias; include a comparison test bounding per-row clipping versus sum-then-clip behavior;
- existing uniform/no-spread behavior is preserved.

### 11.4 Locus Prior

Update `tests/test_per_locus_gdna_mass.py`, `tests/test_assemble_priors.py`, and `tests/test_bayesian_prior_acceptance.py`:

- `expected_gdna_count_global()` still has no payload argument;
- perturbing local payload counts still does not change `gdna_prior_count` when `global_densities` is held fixed;
- global EXON composite rho changes only the EXON-contained prior projection, while boundary-crossing expected counts remain tied to `exon_intron.rho`;
- diagnostic `estimate_locus_gdna()` uses contained local evidence only when strand-active;
- no-eligible-boundary single-exon locus can use contained evidence in stranded mode;
- unstranded single-exon locus preserves boundary/global fallback;
- kappa alpha is converted to opportunity beta before use in precision fusion or `shrink_to_loco()`.

### 11.5 Golden and Benchmark Gates

After unit tests pass, run focused scenario/golden checks that include stranded and unstranded cases. If golden output shifts intentionally, update with:

```bash
conda activate rigel && pytest tests/ --update-golden
```

For benchmark interpretation, compare at least:

- synthetic stranded gDNA/RNA mixtures with large exons;
- unstranded synthetic cases, which should remain close to boundary-only behavior;
- VCaP capture runs, focusing on gDNA false RNA, EXON exposure weights, and single-exon/large-exon loci.
- the fraction of EXON bp/opportunity that remains AMBIG and therefore cannot receive contained power in v1.

## 12. Risk Register

| Risk | Why it matters | Mitigation |
| --- | --- | --- |
| Double-counting boundary and contained fragments | Would inflate exonic gDNA density | Use exact containment gate and keep boundary flux separate |
| Using raw EXON-only counts in unstranded data | Would interpret mature RNA as gDNA | Set contained power to zero unless strand correction is usable |
| Adding EXON-only to gDNA FL training | Would contaminate gDNA FL with RNA and spliced genomic spans | Leave `extract_gdna_counts()` unchanged |
| Ambiguous-strand exons | Dense antisense loci may be exactly where extra EXON evidence is needed | Exclude AMBIG rows in v1, report lost opportunity, and plan per-strand AMBIG sub-regions/payloads |
| Local payload leakage into canonical prior | Breaks Bayesian-prior redesign invariants | Keep `expected_gdna_count_global()` payload-free |
| Overtrusting moderate strand contrast | Deconvolution variance rises as contrast falls | Use strand power plus kappa-capped precision opportunity, not raw opportunity alone |
| Clip-before-sum bias | Many small EXON rows can turn noisy negative moments into an upward rate bias | Keep signed moments through aggregation/shrinkage and floor final rates |
| Silent prior baseline shifts | Composite EXON density changes stranded-library projections | Keep boundary-crossing priors on boundary density and report old-vs-new exonic projection deltas |
| Payload memory growth | New `(n_regions, 3)` int64 array adds per-worker memory | Keep int64 for v1 consistency; measure memory and consider int32 only as a separate schema change |
| Payload sort drift | New array must stay aligned with sorted `RegionArrays` | Add explicit shuffled-order regression test for `PayloadArrays.from_payload()` |
| Schema churn across many tests | Many fixtures instantiate `CalibrationScanPayload` and `GlobalDensityTable` | Add focused helper updates first, then broaden fixture migration |

## 13. Implementation Order

1. Add native payload field and Python payload validation with zero-filled test fixtures.
2. Add shared strand/precision helpers: strand power, alpha-to-beta conversion, and kappa-capped precision opportunity.
3. Add EXON-contained global density and summary diagnostics.
4. Add EXON composite density and update whole-genome plus canonical EXON-contained prior projection.
5. Add regional exposure fusion in `_regional_exposure.py`, including signed-moment handling and final-rate flooring.
6. Add locoregional diagnostic fusion in `locus_prior.py`.
7. Update output diagnostics, debug-script consumers, and tests.
8. Rebuild native extension and run focused tests:

```bash
conda activate rigel && pip install --no-build-isolation -e .
conda activate rigel && pytest tests/test_calibration_accumulator.py -v
conda activate rigel && pytest tests/test_density_global.py -v
conda activate rigel && pytest tests/test_regional_exposure.py -v
conda activate rigel && pytest tests/test_per_locus_gdna_mass.py tests/test_assemble_priors.py tests/test_bayesian_prior_acceptance.py -v
```

Then run the broader suite or benchmark subset requested for the PR.

## 14. Acceptance Criteria

- In SS=0.5 or strand-unidentifiable libraries, EXON-contained counts have zero fusion power and behavior is boundary-only apart from new diagnostics.
- In SS=1.0 libraries with unambiguous EXON rows, pure RNA contained fragments deconvolve to approximately zero gDNA and pure gDNA contained fragments deconvolve to the uncorrected rate.
- Large EXON regions with many contained fragments can outweigh sparse boundary sides only when strand power and kappa-capped precision support it.
- The summary reports the fraction of EXON bp/opportunity that receives contained power, and the fraction blocked by AMBIG strand rows.
- Existing boundary flux diagnostics remain available and interpretable.
- Canonical boundary-crossing prior terms remain tied to `global_densities.exon_intron.rho`; only EXON-contained projection uses `global_densities.exon_composite.rho`.
- `extract_gdna_counts()` continues to exclude `MASK_EXON`.
- `expected_gdna_count_global()` remains payload-free.
- Native worker merge remains deterministic and byte-identical across one versus many scan workers.
- Half-open containment is pinned by tests, including the `frag_end == region.end` edge.
