# Strand Model Implementation Phase 5: Minimal Prior Assembly

Status: design document only. No Phase 5 code has been implemented here.

This document elaborates Phase 5 of `strand_model_impl_plan_v5.md`: converting the
per-region calibration artifacts produced in Phase 4 into per-MultiLocus gDNA prior
inputs for the locus EM. It also corrects three weaknesses in the original Phase 5
sketch:

1. Regions that were not identifiable by strand deconvolution must not be dropped.
2. The EM prior should use the conservative gDNA upper-count surface controlled by
   `rna_lower_confidence`, not the posterior mean count by default.
3. Hybrid-capture effective-length normalization cannot be treated as a solved detail;
   Phase 5 must preserve enough structure for Phase 7 to replace the uniform exposure
   scaffold with a real capture-aware denominator.

## Phase 5 Goal

Phase 4 gives calibration at the fine-region level:

- `CalibrationResult.region_gdna`: per-region strand deconvolution result.
- `CalibrationResult.region_exposure`: per-region exposure weights. This is uniform
  for now and becomes real in Phase 7.
- `index.region_df`: the genome-wide fine-region tiling.

Phase 5 should produce a per-MultiLocus table with the arrays needed by Phase 6 EM:

- `gdna_prior_count`: physical gDNA prior pseudocount for the locus EM.
- `gdna_eff_len`: gDNA effective length/opportunity denominator for the locus EM.
- `enable_gdna`: geometry/data eligibility for allowing the gDNA component.
- Diagnostics that explain how much prior mass came from deconvolved, conservative,
  ineligible, near-unstranded, and exposure-weighted sources.

The key modeling stance is conservative but explicit: if a region contains observed
mass that cannot be confidently assigned to RNA, Phase 5 should carry that mass into
the gDNA prior upper-bound surface and label it. It should not hide the mass by
filtering the region away.

## Current Plan Correction

The original Phase 5 sketch said to filter out rows flagged `FLAG_INELIGIBLE` or
`FLAG_NEAR_UNSTRANDED`, then sum `mean_count` into each locus. That is not adequate.

Filtering flagged rows has the wrong semantics. A flagged row does not mean "no gDNA";
it means "strand information was insufficient to make a direct RNA/gDNA split." Those
regions still contain observed fragments and often are exactly the rows where a
conservative gDNA prior matters.

Using `mean_count` also does not match the intended high-confidence RNA survivor model.
Phase 3/4 already compute both the posterior mean and the confidence-controlled upper
bound:

```text
N_r        = region_gdna.n_total[r]
D_mean_r   = region_gdna.mean_count[r]
R_lower_r  = region_gdna.rna_lower_count[r]
D_upper_r  = region_gdna.upper_count[r] = N_r - R_lower_r
```

For Phase 5, the EM prior count source should be:

```text
D_prior_r = D_upper_r
```

`D_mean_r` should remain in diagnostics and may be used later to estimate exposure,
but it should not be the default EM prior count if the goal is to preserve only
high-confidence RNA.

## Count Semantics

`rna_lower_confidence` controls how conservative `D_prior_r` is.

- At `rna_lower_confidence = 0.95`, `R_lower_r` is the amount of RNA that survives
  with high posterior support, and `D_upper_r = N_r - R_lower_r` becomes a conservative
  gDNA upper count.
- At `rna_lower_confidence = 0.99`, the gDNA prior is even more conservative.
- At `rna_lower_confidence = 0.50`, the prior becomes much closer to a central
  posterior bound, but it is still a quantile-derived bound, not exactly the posterior
  mean.

Phase 5 should therefore keep two separate ideas visible:

```text
gdna_prior_count   = sum allocated D_upper_r       # used by EM
gdna_expected_count = sum allocated D_mean_r       # diagnostic/statistical mean
rna_lower_count    = sum allocated R_lower_r       # high-confidence RNA survivor count
observed_count     = sum allocated N_r             # total allocated observed mass
```

This means the current implementation state after Phase 4 is halfway there: the
needed arrays exist, but the Phase 5 plan must consume `upper_count`, not `mean_count`,
for the EM prior. A literal "use posterior mean as EM prior" mode would require a
separate policy knob. It should not be smuggled into `rna_lower_confidence`, because
that knob already has a clear confidence-bound meaning.

## Handling Excluded Or Unidentified Regions

Phase 5 should never treat a strand-deconvolution exclusion flag as permission to
discard observed mass.

### Exclusion Taxonomy

The current `FLAG_INELIGIBLE` name is too abstract for Phase 5. In the Phase 3
implementation it is emitted when `StrandRegionCounts.eligible` is false or when the
strand-folded total is zero. That conflates several biologically different cases:
zero observed mass, regions with no transcript strand, and true opposite-strand
annotation ambiguity. Phase 5 should not branch on `FLAG_INELIGIBLE` alone. It should
derive an explicit category from `region_arrays.signature`, `region_arrays.ts_class`,
the payload observed counts, and the global strand model.

Recommended categories:

| Category | Detection | Why strand deconvolution is unavailable | Phase 5 strategy |
| --- | --- | --- | --- |
| `ZERO_OBSERVED` | observed region mass is zero | There are no fragments to decompose | Direct EM prior count is zero. Keep the row as exposure/zero-count information for later density fitting, but do not invent prior fragments. |
| `NO_TRANSCRIPT_STRAND` | `ts_class == TS_NONE`, usually intergenic signature `0x0` | There is no transcript-relative sense/antisense axis | Treat observed unspliced mass as high-purity gDNA evidence and as density-training data. This is not a failure case; it is an anchor class. |
| `STRAND_AMBIGUOUS` | `ts_class == TS_AMBIG`, i.e. both POS and NEG transcript bits are active | POS/NEG channels cannot be mapped to a single RNA-major/RNA-minor axis | Use density prediction, not strand deconvolution, for the gDNA prior. Subclass by signature: bistrand intron-only, bistrand exon-only, and mixed exon/intron opposite-strand regions. |
| `NEAR_UNSTRANDED_LIBRARY` | global strand contrast is effectively zero or too small to be useful | RNA and gDNA have nearly the same strand distribution, so the deconvolution denominator is unstable | Disable strand-deconvolution-derived RNA lower bounds globally. Use density/direct-anchor logic instead. |
| `LOW_PRECISION_DECONV` | eligible single-strand region, but posterior precision is near zero | Strand deconvolution ran, but the posterior is broad | Keep the conservative upper bound, but do not use the row as a strong density-training anchor. |
| `KAPPA_FALLBACK_DECONV` | `FLAG_KAPPA_FALLBACK` is set | The gDNA strand overdispersion was not learned from enough seed data | Keep the conservative output, but label it as fallback and downweight for density training. |

The first three categories are the real exclusions from per-region strand
deconvolution. The last two are not exclusions; they are quality labels on rows where
the strand model technically ran but should not be treated as high-confidence training
evidence.

This also exposes a current implementation hazard: `RegionGdnaEstimate.n_total` is the
strand-folded deconvolution total, not necessarily the full observed region total. For
`TS_NONE` and `TS_AMBIG` rows, Phase 3 currently carries zeros through the
strand-folded fields because those rows cannot be expressed as sense/antisense. Phase 5
therefore must either add an explicit `observed_total`/`observed_unspliced_total` field
to the region estimate or accept `PayloadArrays` alongside `RegionGdnaEstimate` and
compute observed mass directly from the payload. It must not infer "zero observed" from
`region_gdna.n_total == 0` unless the payload total is also zero.

### What Near-Unstranded Means

Near-unstranded is a property of the library, not a local region class. The strand
deconvolution model works because RNA and gDNA have different strand expectations:

```text
gDNA: P(POS) ~= 0.5
RNA:  P(RNA-major) = q, where q = max(p_r1_sense, 1 - p_r1_sense)
```

The usable signal is the contrast:

```text
contrast = q - 0.5 = abs(p_r1_sense - 0.5)
```

When that contrast approaches zero, the moment equation divides by a tiny number.
Balanced POS/NEG counts are consistent with both pure gDNA and unstranded RNA, and an
800/200 split no longer has a protocol-defined interpretation if the protocol itself is
not directional. In that setting, strand deconvolution should not claim confident RNA
lower bounds.

The current code uses a numerical floor:

```text
abs(p_r1_sense - 0.5) < 1e-3
```

That catches exactly unstranded or nearly exactly unstranded libraries. For production
modeling, Phase 5/6 should treat this as a practical identifiability question, not only
a floating-point guard. A library with `p_r1_sense = 0.55` has a nonzero contrast, but
the posterior intervals will be very wide and many regions will have little usable
deconvolution precision. Recommended diagnostics:

```text
signed_contrast = 2 * p_r1_sense - 1
contrast_margin_99 = strand_summary.signed_strand_contrast_margin(confidence=0.99)
strand_identifiable = abs(signed_contrast) > max(contrast_margin_99, min_practical_contrast)
```

where `min_practical_contrast` should be chosen empirically. Until that threshold is
benchmarked, the numerical floor can remain a hard guard, but the summary should not
describe all weak-strand libraries as equally identifiable.

### Strand-Ambiguous Regions Need Density

Transcript-strand-ambiguous regions are the central unsolved case. Suppose two
transcripts have overlapping exons on opposite strands and the region has 1000 total
fragments: 800 POS and 200 NEG. With a strongly stranded protocol, this implies a real
strand imbalance, but it does not tell us how much is gDNA, POS-strand RNA, or
NEG-strand RNA.

In simplified form, the two observed channels are trying to explain three unknowns:

```text
POS = 0.5 * D + q       * R_pos + (1 - q) * R_neg
NEG = 0.5 * D + (1 - q) * R_pos + q       * R_neg
N   = D + R_pos + R_neg
```

The total equation is redundant with POS + NEG. Without another source of information,
`D`, `R_pos`, and `R_neg` are underdetermined. This is exactly where an orthogonal gDNA
density model is needed. The density model should predict plausible gDNA count from
region size, FL-aware opportunity, region class/signature, local/capture exposure, and
training regions where gDNA is identifiable. The EM can then use transcript likelihoods
to decide how the remaining RNA mass splits across transcripts.

The temporary Phase 5 answer should not be "all ambiguous mass is gDNA" except as a
last-resort fallback when no density model is available. That fallback is safe against
false RNA, but it will suppress real RNA in exactly the difficult opposite-strand exon
cases. A lightweight density predictor is therefore part of the minimal useful Phase 5,
even if the full capture-aware model remains a later phase.

The required minimal policy is:

1. Include every region that overlaps a MultiLocus footprint and has nonzero allocated
   observed mass.
2. Use `upper_count` for directly strand-deconvolved, single-strand rows.
3. Use direct observed gDNA-compatible mass for intergenic/no-transcript-strand anchor
   rows.
4. Use a density-predicted upper count for strand-ambiguous rows and for globally
   near-unstranded libraries.
5. If density prediction is unavailable, fall back to `D_upper_r = N_r` only with an
   explicit `DENSITY_MISSING_ALL_GDNA_FALLBACK` diagnostic flag.
6. Track those rows separately in diagnostics so the output tells us whether a locus
   prior was informed by identifiable strand imbalance or by conservative fallback.

Recommended per-locus diagnostics:

```text
observed_count
gdna_prior_count
gdna_expected_count
rna_lower_count
deconvolved_prior_count
conservative_prior_count
ineligible_prior_count
near_unstranded_prior_count
ambiguous_signature_prior_count
density_predicted_prior_count
density_fallback_prior_count
zero_observed_region_bp
allocated_region_bp
n_regions_total
n_regions_deconvolved
n_regions_conservative
n_regions_density_predicted
n_regions_strand_ambiguous
```

If a future Phase 3 change ever emits `NaN` or missing values for unavailable rows,
Phase 5 must coerce the prior count for those rows to the conservative bound, not to
zero. Zero would mean "known RNA or no contamination," which is not what an
unidentified row means.

### Alternatives Considered

Dropping flagged rows is rejected because it loses observed mass and biases the prior
toward RNA in exactly the regions where strand information is weak.

Treating every excluded row as all-gDNA is also rejected as the default for
strand-ambiguous exon regions. It is a useful emergency bound, but it throws away the
main biological distinction between "unknown because no reads" and "unknown because
two opposite-strand RNAs could explain the imbalance."

The recommended compromise is a minimal density predictor before full Phase 5 EM
wiring. It can be much simpler than the eventual Phase 7 capture-exposure model, but it
must provide bounded predictions for strand-ambiguous regions. The all-gDNA bound then
remains a flagged last resort, not the ordinary path.

## Region To MultiLocus Allocation

The EM unit is `MultiLocus`, not a gene and not a single transcript. A `MultiLocus`
contains one or more genomic `Locus` intervals:

```text
MultiLocus(
    multi_locus_id,
    transcript_indices,
    unit_indices,
    gdna_span,
    loci: tuple[Locus, ...],
)
```

Phase 5 should aggregate fine-region calibration rows over the union of the `loci`
intervals for each `MultiLocus`.

For each region `r` and MultiLocus `m`, define:

```text
B_m             = union of genomic intervals in multi_locus.loci
overlap_bp_mr   = bp overlap between B_m and region r on the same reference
raw_share_mr    = overlap_bp_mr / region_length_r
```

Most region boundaries should align with transcript-span boundaries because the
fine-region builder uses exon and intron events. The overlap code should still handle
partial overlaps, because intergenic rows can be large and future index changes should
not require rewriting the prior assembler.

### Avoiding Double Counting Across Overlapping MultiLoci

Different MultiLoci can overlap the same genomic bases, especially in opposite-strand
or otherwise ambiguous annotation contexts. Blindly adding the full region count to
every overlapping MultiLocus duplicates calibration mass.

Phase 5 should use a conservation-aware share:

```text
sum_raw_r = sum_m raw_share_mr
scale_r   = 1                        if sum_raw_r <= 1
          = 1 / sum_raw_r             if sum_raw_r > 1
share_mr  = raw_share_mr * scale_r
```

Then for every contribution:

```text
prior_count_m    += share_mr * D_upper_r
expected_count_m += share_mr * D_mean_r
rna_lower_m      += share_mr * R_lower_r
observed_m       += share_mr * N_r
```

This has three useful properties:

- If one locus covers half of a large region, it receives half the region's count;
  the rest remains outside the EM footprint.
- If two disjoint loci cover different parts of the same region, each receives its
  bp-proportional share.
- If two overlapping MultiLoci both cover the same bases, the row's count is split
  rather than duplicated.

The per-region invariant should be:

```text
sum_m share_mr <= 1
```

and therefore:

```text
sum_m contribution_mr <= region_count_r
```

That invariant should be tested directly.

## PriorTable Schema

Phase 5 should introduce a small immutable table, probably in a new module such as
`src/rigel/calibration/prior.py`.

Proposed schema:

```python
@dataclass(frozen=True, slots=True)
class PriorTable:
    n_multi_loci: int

    # EM-facing count surface.
    gdna_prior_count: np.ndarray        # float64[M], sum of allocated upper_count
    gdna_prior_count_em: np.ndarray     # float64[M], same as above in Phase 5

    # Diagnostics and future exposure fitting.
    gdna_expected_count: np.ndarray     # float64[M], sum of allocated mean_count
    rna_lower_count: np.ndarray         # float64[M], sum of allocated rna_lower_count
    observed_count: np.ndarray          # float64[M], sum of allocated n_total
    density_predicted_count: np.ndarray # float64[M], density-model contribution
    density_upper_count: np.ndarray     # float64[M], density-model upper contribution

    # Effective-length denominator.
    gdna_eff_len_raw: np.ndarray        # float64[M], unweighted FL-marginal denominator
    gdna_eff_len: np.ndarray            # float64[M], exposure-weighted denominator
    gdna_em_exposure_weight: np.ndarray # float64[M], gdna_eff_len / gdna_eff_len_raw

    # Eligibility and diagnostics.
    enable_gdna: np.ndarray             # uint8[M]
    flags: np.ndarray                   # uint32[M]
    n_regions_total: np.ndarray         # int32[M]
    n_regions_conservative: np.ndarray  # int32[M]
    n_regions_deconvolved: np.ndarray   # int32[M]
    conservative_prior_count: np.ndarray
    deconvolved_prior_count: np.ndarray
    ineligible_prior_count: np.ndarray
    near_unstranded_prior_count: np.ndarray
    no_transcript_strand_prior_count: np.ndarray
    strand_ambiguous_prior_count: np.ndarray
    density_predicted_prior_count: np.ndarray
    density_fallback_prior_count: np.ndarray
```

`gdna_prior_count_em` is included because the current pipeline metadata already has
both raw and EM-facing prior-count fields. In Phase 5 they should be identical. Later
phases may choose to regularize or cap the EM-facing count, but any difference must be
explicitly visible in diagnostics. Exposure should generally affect the denominator,
not silently rescale the observed-count numerator.

### gDNA Eligibility

Do not define `enable_gdna` as `gdna_prior_count > 0`. That would disable the gDNA
component in a locus with zero prior but valid unspliced gDNA likelihood evidence.

The safer contract is:

```text
enable_gdna_geometry_m = gdna_eff_len_m > 0 and locus has nonempty genomic footprint
```

Phase 6 can combine this with the partition-level condition already used by
`AbundanceEstimator.run_batch_locus_em_partitioned`: at least one unspliced EM unit
has finite `gdna_log_lik`. In other words, the prior table should not use a zero
prior count as a hard biological impossibility statement.

## Effective Length And Hybrid Capture

There are two separate quantities:

```text
count numerator       = how much observed mass is available as gDNA prior
opportunity denominator = how much gDNA-accessible sequence could have generated it
```

For directly strand-deconvolved rows, the numerator comes from
`region_gdna.upper_count`. For strand-ambiguous rows, it comes from the density-v0
upper count capped by observed local mass. Hybrid-capture correction belongs primarily
in the denominator and exposure model, not by dropping observed counts.

### Phase 5 Minimal Denominator

For Phase 5, use the existing FL-aware geometry helper:

```text
gdna_eff_len_raw_m = gdna_eff_len_for_loci(multi_locus.loci, index.ref_lengths, gdna_fl_model)
```

Then apply the current exposure scaffold:

```text
A_m                 = bp_weighted_mean_exposure_over_blocks(multi_locus.loci, region_arrays,
                                                            calibration.region_exposure)
gdna_eff_len_m       = max(gdna_eff_len_raw_m * A_m, min_value)
gdna_em_exposure_weight_m = gdna_eff_len_m / gdna_eff_len_raw_m
```

Because Phase 4 exposure is uniform, Phase 5 should effectively produce:

```text
A_m = 1
gdna_eff_len_m = gdna_eff_len_raw_m
```

This is acceptable only as scaffolding. It is not a complete hybrid-capture solution.
The key Phase 5 requirement is to preserve `gdna_eff_len_raw`, `A_m`, and
`gdna_eff_len` separately so Phase 7 can replace `A_m` without changing the EM
contract.

### Why Uniform Exposure Is Not Enough

In hybrid-capture data, a large genomic locus may span 1 Mb while only 10 kb of that
span is efficiently captured. A raw locus-span denominator treats the full 1 Mb as
equally available to gDNA. That makes the gDNA component's effective length too large
relative to the actual capture-accessible footprint and can distort the abundance
scale. Exact FL geometry alone does not fix this if the exposure field is still flat.

Therefore Phase 7 should not be framed as a cosmetic improvement. It is the stage that
must learn the accessible gDNA opportunity surface.

## Phase 7 Exposure Vision Needed By Phase 5

Phase 5 should be designed so that Phase 7 can replace the scalar uniform exposure
with a capture-aware exposure surface. The desired future model is:

```text
D_mean_r ~ count model with mean eta_g * A_r * L_eff_r
```

where:

- `D_mean_r` is the deconvolved posterior mean count, not the conservative upper bound.
- `eta_g` is a global or library-wide gDNA density scale.
- `A_r` is a relative capture/accessibility weight for region `r`.
- `L_eff_r` is the fragment-length-aware effective opportunity for region `r`.

Important design guidance for Phase 7:

1. Use `mean_count` to estimate exposure, because exposure is an opportunity model.
   Use `upper_count` for conservative EM prior counts.
2. Downweight or label low-precision and fallback rows when learning `A_r`, but do not
   erase them from the count ledger.
3. Avoid a plain high-quantile cap as the whole normalization strategy. Previous
   experiments showed that exact denominators or high-tail caps alone do not solve the
   leakage problem. The hard part is estimating the exposure field and its scale.
4. Enforce count closure as a diagnostic: the exposure model should predict the total
   deconvolved gDNA mean count within tolerance when summed over regions.
5. For large captured loci, `A_m` should reflect the captured subfootprint. A 1 Mb locus
   with 10 kb of exposed target should get a denominator on the 10 kb-accessible scale,
   not the raw 1 Mb span scale.

The first Phase 7 implementation can use `gdna_eff_len_raw * bp_weighted_mean(A_r)` as
the scalar approximation. The interface should still leave room for a more exact
fragment-length marginal integral:

```text
L_gdna_m(A) = sum_l P(length=l) * integral 1(fragment_start, l overlaps B_m)
                                      * A(fragment_start, l) d fragment_start
```

The scalar approximation is probably fine for Phase 5/7 bring-up; the critical part is
that `A_r` must be biologically meaningful and normalized, not that the final integral
is implemented first.

## Recommended Roadmap

The best path forward is not to implement Phase 5 as originally scoped and hope the
ambiguous regions are rare enough. The opposite-strand ambiguous regions are exactly
where strand deconvolution has no mathematical answer, and they need an orthogonal
density signal. However, we also do not need the full final density/exposure model
before any prior assembly work can proceed.

Recommended organization:

```text
Phase 5A: classify region inference mode
Phase 5B: implement density-v0 for non-strand-deconvolvable rows
Phase 5C: merge strand and density surfaces into PriorTable
Phase 6:  wire PriorTable into EM
Phase 7:  replace density-v0/uniform exposure with capture-aware exposure modeling
```

### Phase 5A: Region Inference Mode Table

Add an explicit per-region inference-mode table before aggregating anything to loci.
This can be a lightweight dataclass or extra arrays on a region-prior helper:

```text
ZERO_OBSERVED
DIRECT_INTERGENIC_GDNA
STRAND_DECONVOLVED_SINGLE_STRAND
DENSITY_PREDICTED_STRAND_AMBIGUOUS
DENSITY_PREDICTED_NEAR_UNSTRANDED
CONSERVATIVE_ALL_GDNA_FALLBACK
```

Inputs:

- `region_arrays.signature` and `region_arrays.ts_class`.
- payload observed unspliced totals and, where useful, spliced totals.
- `region_gdna.flags`, `mean_count`, `upper_count`, and `precision`.
- global `StrandSummary` identifiability diagnostics.

Outputs:

- a count source for each region: strand upper bound, direct gDNA anchor, density
  prediction, or zero;
- a quality/source label for diagnostics and downstream training exclusion;
- an explicit distinction between "zero observed" and "not strand-identifiable."

### Phase 5B: Density-v0 For Ambiguous Rows

Implement a deliberately small density model before full prior assembly. Its job is
not to be the final hybrid-capture solution; its job is to avoid the false choice
between dropping ambiguous regions and treating all ambiguous mass as gDNA.

Training rows should come from identifiable gDNA evidence:

- intergenic/no-transcript-strand rows with observed unspliced mass;
- single-strand intron-only rows with low spliced evidence;
- no-RNA exon self-training rows accepted by the Phase 3 screen;
- high-precision strand-deconvolved single-strand rows, using `mean_count` for the
  density mean and `upper_count` only for conservative bounds.

Rows to exclude or downweight for density fitting:

- `STRAND_AMBIGUOUS` rows, because those are the prediction target;
- low-precision strand-deconvolved rows;
- hard `kappa_d` fallback rows;
- rows with strong spliced RNA evidence unless specifically selected as no-RNA exons.

A useful density-v0 model can be class-stratified and simple:

```text
rho_class = sum_r D_mean_r / sum_r L_eff_r       over training rows in class
D_pred_r  = rho_class(signature_class_r) * L_eff_r
D_upper_r = q_alpha_count_model(D_pred_r, dispersion_class)
D_prior_r = min(observed_unspliced_r, D_upper_r)
```

where `signature_class` can start as intergenic, intron-only, exon-only, and mixed;
`L_eff_r` is the FL-aware regional opportunity; and the upper bound can initially be a
Poisson or negative-binomial quantile with conservative dispersion. The cap by observed
unspliced mass is important: the density model predicts how much gDNA is plausible
inside the observed local fragment ledger, not permission to inject more local prior
fragments than were observed in that region.

If density-v0 has no valid training data for a class, fall back in this order:

1. Same-reference density from all identifiable anchors.
2. Genome-wide density from all identifiable anchors.
3. All-gDNA upper bound, flagged as `DENSITY_MISSING_ALL_GDNA_FALLBACK`.

For the concrete 1000-fragment, 800/200 opposite-strand exon example, density-v0 would
not try to decide whether the 800-channel excess is POS-strand RNA or NEG-strand RNA.
It would estimate a plausible `D_prior` from size/type/exposure. If it predicts 120
gDNA fragments, Phase 5 passes roughly 120 gDNA prior fragments to EM and leaves the
remaining mass to transcript likelihoods. If it predicts 500, EM sees a much stronger
gDNA prior. The strand counts alone cannot choose between those outcomes.

### Phase 5C: Merge Surface For Prior Assembly

After 5A and 5B, Phase 5 prior assembly is straightforward:

```text
single-strand identifiable rows -> region_gdna.upper_count
intergenic/no-transcript rows   -> direct observed gDNA-compatible count
strand-ambiguous rows           -> density_upper_count, capped by observed mass
near-unstranded rows            -> density_upper_count/direct anchors only
zero-observed rows              -> zero count
```

Then aggregate those region-level counts into MultiLoci using the region-to-locus
allocation algorithm described below.

### Should Phase 5 Pause?

Do not wire Phase 5 directly into EM until at least density-v0 exists. Otherwise the
only honest options for strand-ambiguous rows are all-gDNA or zero, and both are bad
defaults. But also do not pause for a full capture-aware exposure model. The organized
middle path is:

1. Implement region inference-mode classification and diagnostics.
2. Implement density-v0 with simple class/reference/global fallback.
3. Merge strand and density counts into `PriorTable`.
4. Run targeted simulations where ambiguous opposite-strand exons are the only hard
   case.
5. Then wire into EM and iterate on Phase 7 exposure if hybrid-capture loci still show
   denominator-driven leakage.

## Implementation Outline

1. Add region inference-mode classification and diagnostic constants. Do not use
   `FLAG_INELIGIBLE` as a control-flow category.
2. Add `observed_unspliced_total` or pass `PayloadArrays` into prior assembly so
   `TS_NONE` and `TS_AMBIG` rows do not look falsely empty.
3. Add density-v0 for strand-ambiguous and near-unstranded rows.
4. Add `PriorTable` and `assemble_locus_priors(...)` in `rigel.calibration.prior`.
5. Build a compact region-array representation from `index.region_df` if the existing
   `RegionArrays` object is not sufficient for overlap iteration.
6. For each reference, sweep sorted MultiLocus intervals against sorted region rows.
7. Emit raw `(multi_locus_id, region_id, overlap_bp, raw_share)` contributions.
8. Normalize per-region shares only when the raw shares sum above 1.
9. Aggregate region counts by source: strand upper bounds, direct anchors,
   density-predicted upper counts, conservative fallbacks, and zero-observed rows.
10. Compute raw and exposure-weighted gDNA effective lengths.
11. Fill diagnostics and return `PriorTable`.
12. In Phase 6, pass `prior_table.gdna_prior_count_em` and `prior_table.gdna_eff_len`
   to `_run_locus_em_partitioned`.

## Tests And Acceptance Criteria

Unit tests should cover:

1. A directly deconvolved region contributes `upper_count`, not `mean_count`, to
   `gdna_prior_count`.
2. `mean_count` is still aggregated into `gdna_expected_count`.
3. Increasing `rna_lower_confidence` increases or preserves `gdna_prior_count` for the
   same region evidence.
4. `FLAG_INELIGIBLE` rows are included conservatively, not filtered out.
5. `FLAG_NEAR_UNSTRANDED` rows are included conservatively, not filtered out.
6. `ZERO_OBSERVED` rows produce zero prior count even when they have nonzero exposure.
7. `TS_NONE` intergenic rows with observed unspliced mass become direct gDNA anchors,
   not `FLAG_INELIGIBLE` losses.
8. `TS_AMBIG` rows use density-predicted prior count, not strand-deconvolved count.
9. The 1000-fragment 800/200 opposite-strand example is underdetermined by strand
   equations alone and follows the density path.
10. Density-v0 caps predicted local prior count by observed local unspliced mass.
11. Partial region overlap prorates counts by `overlap_bp / region_length`.
12. Overlapping MultiLoci split duplicated region share so `sum_m share_mr <= 1`.
13. MultiLocus objects with multiple `Locus` intervals aggregate across all intervals
   without double-counting their own footprint.
14. Uniform exposure gives `gdna_eff_len == gdna_eff_len_raw` and exposure weight 1.
15. Zero prior count does not automatically force `enable_gdna` to 0 when the geometry
    is valid.

End-to-end acceptance for Phase 5:

- No observed region mass is lost because of strand-deconvolution flags.
- The default EM prior count is the confidence-controlled gDNA upper bound.
- Strand-ambiguous regions are routed through density-v0 or explicitly flagged
   all-gDNA fallback; they are not silently treated as strand-deconvolved rows.
- Per-region contribution conservation is enforced when MultiLoci overlap.
- The resulting table contains enough diagnostics to explain why a locus has a large
  gDNA prior.
- Hybrid-capture exposure remains explicitly marked as uniform scaffold, with raw and
  exposure-weighted denominators both retained for Phase 7.

## Open Questions

1. Should Rigel eventually expose an explicit `gdna_prior_count_source = upper|mean`
   knob? The current recommendation is no for Phase 5. Keep the single confidence
   knob and use `upper_count` for the EM prior.
2. Should overlapping MultiLoci split region mass or receive independent copies of the
   same prior evidence? This design recommends conservation-aware splitting to avoid
   prior inflation. If that underpowers specific opposite-strand cases, we should
   change it with targeted diagnostics rather than duplicating mass by default.
3. Should Phase 7 use scalar `A_m` or exact FL-marginal exposure integration first?
   Start with scalar `A_m` for simplicity, but keep the schema ready for exact
   integration. The exposure estimate itself is the bigger modeling risk.
