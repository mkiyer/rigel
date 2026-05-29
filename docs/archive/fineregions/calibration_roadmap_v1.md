# Integrated Calibration Roadmap v1

Date: 2026-05-24
Status: roadmap document only
Supersedes/organizes: `strand_model_impl_plan_v5.md`, `strand_model_impl_phase5.md`, `density_model_plan_v1.md`

Note: phase numbers in this document define a new integrated calibration roadmap. They
are intentionally reset around the density/strand/exposure integration problem.

## 0. Executive Summary

Rigel now has the beginning of two complementary calibration evidence channels:

1. **Strand deconvolution**: local evidence. It asks whether the observed POS/NEG
   balance in a fine region is better explained by RNA, gDNA, or a mixture.
2. **Density modeling**: opportunity evidence. It asks how much gDNA is plausible in a
   region given fragment-length-aware opportunity, region class, boundary flux, and
   eventually capture exposure.

These channels are orthogonal. The strand model learns from directional imbalance; the
density model learns from where gDNA could have arisen and from how much gDNA-like
flux is observed in anchor regions. Neither dominates all cases.

The immediate implementation priority is **Phase 0: integration cleanup and synthesis**.
Before adding more modeling, calibration should expose one coherent region-evidence
system that owns:

- scan payload arrays,
- region geometry,
- library strand diagnostics,
- current strand deconvolution output,
- current uniform exposure scaffold,
- placeholders for density predictions,
- a single integration contract consumed by future prior assembly.

After Phase 0, implement density as a first-class evidence channel, then fuse density
and strand into a region-level gDNA prior surface, then aggregate that surface to
MultiLoci and wire EM.

## 1. First Principles

The calibration problem is not "estimate one global gDNA fraction." It is:

```text
For each genomic fine region r:
  observed fragments_r = RNA_r + gDNA_r + noise/unsupported classes

We need enough evidence to decide how much local observed mass should enter the
per-locus gDNA prior, and how much gDNA opportunity should enter the EM denominator.
```

There are two different questions:

```text
Count numerator:
  How much observed mass in this region is plausibly gDNA?

Opportunity denominator:
  How much accessible sequence could have generated gDNA fragments here?
```

Strand deconvolution helps the numerator by using local POS/NEG imbalance. Density
modeling helps both numerator and denominator by predicting plausible gDNA counts from
opportunity. Exposure modeling modifies the denominator by learning which genomic
opportunities were actually accessible in a capture library.

A correct calibration system must keep these concepts separate:

```text
observed mass        -> count ledger
strand prediction    -> local RNA/gDNA split evidence
density prediction   -> opportunity-based gDNA plausibility
exposure             -> accessible-opportunity correction
prior table          -> MultiLocus EM contract
```

## 2. The Four RNA-seq Use Cases

Rigel does not know the use case in advance. It can learn strand specificity early, and
it must learn or infer capture/exposure later.

### Case 1: unstranded, no hybrid capture

```text
strand model: OFF
density model: ON
exposure: uniform is usually acceptable
```

Strand deconvolution has no usable contrast because RNA and gDNA have nearly the same
strand distribution. Density owns the gDNA prior. Intergenic and intronic contained
regions are representative anchors, and boundary flux adds useful support.

Implementation stance:

- Do not emit confident RNA lower bounds from strand.
- Fit density from gDNA-compatible anchors.
- Use density upper counts for exonic and ambiguous regions.
- Use uniform exposure until Phase 7 unless diagnostics indicate otherwise.

### Case 2: unstranded, hybrid capture

```text
strand model: OFF
density model: ON, boundary-driven
exposure: learned in Phase 7
```

Intergenic and intronic contained regions are depleted by capture and are not
representative of on-target gDNA density. Boundary flux at exon-intron and
exon-intergenic edges becomes essential because it samples gDNA-like fragments near
captured regions.

Implementation stance:

- Use boundary flux as a first-class density observation, not a diagnostic afterthought.
- Fit density-v0 from boundary and target-proximal anchors when background contained
  anchors are inconsistent with boundary anchors.
- Before Phase 7, preserve raw and exposure-weighted denominators separately and mark
  exposure as uniform scaffold.
- Phase 7 learns `RegionExposure.A_r` so a large raw genomic span does not imply a large
  accessible gDNA denominator.

### Case 3: strand-specific, no hybrid capture

```text
strand model: ON
density model: ON
exposure: uniform is usually acceptable
```

This is the cleanest synergy case. Strand gives a local RNA/gDNA split in single-strand
regions. Density gives an independent prediction from region class and opportunity,
handles regions where strand is weak, and acts as a stabilizer when strand posteriors
are broad.

Implementation stance:

- Use strand deconvolution directly for single-transcript-strand regions.
- Use density for `TS_AMBIG`, `TS_NONE`, near-unstranded, and low-precision rows.
- Compare strand and density predictions in regions where both are available.
- Eventually perform Bayesian fusion: density becomes the prior over `D`, strand counts
  become the likelihood.

### Case 4: strand-specific, hybrid capture

```text
strand model: ON
density model: ON, boundary-aware
exposure: learned in Phase 7
```

The beauty of strand-specific capture data is that entire exons can be deconvolved from
antisense/sense information without relying on boundary flux. But density remains
necessary for strand-ambiguous regions and for denominator/exposure normalization.

Implementation stance:

- Use strand deconvolution for identifiable single-strand exons and introns.
- Use density for opposite-strand ambiguous regions.
- Train initial density from exon-intron and exon-intergenic boundary flux because
  intergenic/intronic contained background is depleted.
- Use strand-derived `mean_count` as an additional exposure-learning signal in Phase 7.

## 3. Current Implementation Inventory

### Implemented and useful

- `scan_payload.py`: owns the fractional 12-channel scan payload, including contained
  and boundary left/right unspliced POS/NEG counts.
- `_arrays.py`: builds `RegionArrays` and `PayloadArrays`, sorted by reference. This is
  the right shared geometry/count view for both density and strand.
- `strand_summary.py`: owns library strand summary, `p_r1_sense`, signed contrast, and
  the numerical floor for near-unstranded behavior.
- `fl.py`: builds RNA/gDNA/global fragment-length models used by both scoring and
  calibration geometry.
- `density_global.py`: contains the current contained-only global density estimates,
  `l_eff_contained`, and `estimate_strand_balance` for gDNA strand overdispersion.
- `_exposure.py`: contains FL-aware geometry primitives:
  `contained_exposure_clipped`, `fractional_boundary_side_exposure`,
  `gdna_eff_len_for_loci`, and `bp_weighted_mean_exposure_over_blocks`.
- `strand_deconv.py`: implements region-level strand folding, kappa_d estimation,
  exon self-training, and `RegionGdnaEstimate` with `mean_count`, `upper_count`,
  `rna_lower_count`, precision, and flags.
- `exposure.py`: defines `RegionExposure`; currently only `uniform` is implemented.
- `_orchestrator.py`: wires FL models, current global densities, strand counts,
  kappa_d, region strand deconvolution, and uniform exposure into `CalibrationResult`.
- `_result.py`: summarizes calibration config, FL models, diagnostics, strand
  deconvolution, and region exposure.
- `pipeline.py`: already has the lower-level EM function signature ready for
  `gdna_prior_count_em`, `gdna_eff_len`, `enable_gdna`, raw prior counts, unweighted
  denominators, and exposure weights. `quant_from_buffer` still stops before prior
  assembly and EM.

### Missing pieces

- A unified per-region evidence table that joins payload totals, region class,
  strand prediction, density prediction, and exposure state.
- A density observation builder that includes boundary flux and FL-aware boundary
  opportunity.
- A density model module with per-region predictions, upper bounds, diagnostics, and
  fallback behavior.
- A formal integration layer that decides how strand and density cooperate per region.
- A `PriorTable` and conservation-aware region-to-MultiLocus allocation.
- End-to-end pipeline wiring past the Phase 6 boundary.
- Unsupervised exposure learning for hybrid capture.
- Validation across all four use cases.

## 4. Calibration Architecture Target

The target architecture should look like this:

```text
scan payload + index region_df
  -> RegionArrays + PayloadArrays
  -> FLModels
  -> StrandSummary
  -> RegionEvidenceTable
       observed ledgers
       region classes
       strand prediction
       density observation
       density prediction
       exposure state
       integration decision
  -> IntegratedRegionGdnaSurface
       D_prior_upper_r
       D_expected_r
       R_lower_r
       count source flags
       density/source diagnostics
  -> PriorTable over MultiLoci
       gdna_prior_count_em
       gdna_prior_count_raw
       gdna_expected_count
       gdna_eff_len_raw
       gdna_eff_len
       gdna_em_exposure_weight
       enable_gdna
       diagnostics
  -> locus EM
```

The important shift is that neither `strand_deconv.py` nor `density_model.py` should
own prior assembly. They each produce evidence. A separate integration layer owns the
policy for combining evidence and producing the region surface consumed by prior
assembly.

## 5. Region Evidence Types

Every fine region should be assigned an evidence mode. This should not be inferred from
`FLAG_INELIGIBLE` alone.

Recommended modes:

```text
ZERO_OBSERVED
DIRECT_GDNA_ANCHOR_NO_TRANSCRIPT_STRAND
STRAND_DECONVOLVED_SINGLE_STRAND
DENSITY_PREDICTED_STRAND_AMBIGUOUS
DENSITY_PREDICTED_NEAR_UNSTRANDED
DENSITY_PREDICTED_LOW_PRECISION
CONSERVATIVE_ALL_GDNA_FALLBACK
```

Key rule:

```text
A strand-deconvolution flag never means "drop observed mass."
```

It means either no strand axis existed, strand was not identifiable, or the posterior
was broad. The count ledger remains real and must be carried to density, prior, or a
flagged conservative fallback.

## 6. How Strand And Density Should Integrate

### Independent evidence surfaces

For each region, construct both possible predictions when possible:

```text
StrandPrediction:
  D_mean_strand
  D_upper_strand
  R_lower_strand
  precision_strand
  strand_flags

DensityPrediction:
  D_mean_density
  D_upper_density
  density_leff
  density_class
  density_quality
  density_flags
```

`D_mean_strand` and `D_mean_density` are not interchangeable:

- strand mean is local evidence from POS/NEG imbalance;
- density mean is opportunity evidence from similar gDNA-compatible regions.

`D_upper_strand` and `D_upper_density` are conservative count surfaces suitable for
RNA-protective prior construction.

### Integration v0: source-gated fusion

Use simple source gating before implementing full Bayesian fusion:

```text
ZERO_OBSERVED:
  D_prior = 0
  D_expected = 0

DIRECT_GDNA_ANCHOR_NO_TRANSCRIPT_STRAND:
  D_prior = observed_unspliced_total
  D_expected = observed_unspliced_total
  training_role = density anchor

STRAND_DECONVOLVED_SINGLE_STRAND:
  D_prior = D_upper_strand
  D_expected = D_mean_strand
  density_prediction retained as diagnostic and exposure-learning input

DENSITY_PREDICTED_STRAND_AMBIGUOUS:
  D_prior = min(observed_compatible_count, D_upper_density)
  D_expected = min(observed_compatible_count, D_mean_density)

DENSITY_PREDICTED_NEAR_UNSTRANDED:
  D_prior = min(observed_compatible_count, D_upper_density)
  D_expected = min(observed_compatible_count, D_mean_density)

DENSITY_PREDICTED_LOW_PRECISION:
  D_prior = max(D_upper_strand, min(observed_compatible_count, D_upper_density))
  D_expected = blend_or_report_only(D_mean_strand, D_mean_density)

CONSERVATIVE_ALL_GDNA_FALLBACK:
  D_prior = observed_compatible_count
  D_expected = observed_compatible_count
  flag loudly
```

The low-precision row is the only intentionally conservative `max` path. If strand ran
but was uninformative, we should not let a weak density prediction erase observed mass.
This can be tightened after validation.

### Integration v1: Bayesian fusion

The cleaner long-term model is:

```text
P(D | strand counts, density model)
  proportional to P(strand counts | D, N, p_r1_sense, kappa_d)
               * P_density(D | L_eff, class, exposure)
```

This keeps density and strand in their natural roles:

- density supplies the prior over `D`;
- strand supplies the likelihood from POS/NEG imbalance;
- the fused posterior yields `D_mean_fused`, `D_upper_fused`, and `R_lower_fused`.

Do not block the roadmap on this. Build integration v0 first, but design schemas so v1
can replace source-gated fusion without changing the EM contract.

## 7. Ordered Implementation Roadmap

### Phase 0: Integration cleanup and synthesis

**Goal:** Consolidate the current implemented pieces into one integrative calibration
system without adding a new statistical model yet.

**Why first:** We already have FL, strand summary, region arrays, global contained
rho, strand deconvolution, and uniform exposure. They are wired, but not organized as
complementary evidence channels. Phase 0 creates the place where density will plug in.

**Scope:**

- Add a small module, likely `calibration/evidence.py`.
- Define a `RegionEvidenceTable` or equivalent dataclass containing:
  - `RegionArrays`, or references to its structural arrays;
  - observed contained/boundary unspliced totals;
  - observed all-compartment totals for strand accounting;
  - region class / transcript strand class;
  - `StrandRegionCounts`;
  - `RegionGdnaEstimate`;
  - current `RegionExposure`;
  - placeholder density fields filled with `NaN` or zeros and clear mode flags.
- Move evidence construction out of `_orchestrator.py` into named helpers.
- Add a `CalibrationModeSummary`:
  - strand identifiable / not identifiable;
  - exposure mode (`uniform`, future `unsupervised`);
  - density mode (`not_built` in Phase 0);
  - counts of `TS_NONE`, `TS_POS`, `TS_NEG`, `TS_AMBIG`, zero-observed rows.
- Keep `CalibrationResult.region_gdna` and `region_exposure` stable.

**Tests:**

- Evidence table row count matches `RegionArrays` and `RegionGdnaEstimate`.
- `TS_NONE` and `TS_AMBIG` rows preserve observed unspliced mass even when
  strand-folded `region_gdna.n_total` is zero.
- Strand identifiable summary agrees with `StrandSummary` contrast diagnostics.
- Existing calibration tests still pass.

**Acceptance:**

- Calibration summary explicitly says strand model active/inactive, density model
  absent, and exposure uniform.
- No prior assembly yet.
- No behavior change to `quant_from_buffer` boundary.

### Phase 1: Unified density observations

**Goal:** Build per-region density observations `(count, L_eff)` that include contained
and boundary gDNA-compatible evidence.

**Scope:**

- Add `density_observation.py` or extend `density_global.py` with:

```text
compute_density_observations(region_arrays, payload_arrays, gdna_fl)
```

- Output arrays:
  - contained unspliced count;
  - left/right boundary unspliced count;
  - total density count;
  - contained effective length;
  - left/right boundary effective opportunity;
  - total density effective opportunity;
  - class labels: intergenic, intron-only, exon-only, exon-intron, exon-intergenic,
    strand-ambiguous.
- Use `_exposure.fractional_boundary_side_exposure` for boundary opportunity.
- Keep current contained-only `compute_global_densities` diagnostics, but add the
  unified observation block alongside it.

**Tests:**

- With zero boundary counts, unified pooled rho matches contained-only rho within
  rounding.
- Boundary counts increase both numerator and denominator.
- Boundary side opportunity is nonnegative and FL-sensitive.
- Unstranded and stranded payloads collapse POS/NEG correctly for density counts.

**Acceptance:**

- Density observation arrays exist for all regions.
- Boundary flux is visible in summary diagnostics.
- No density model fit yet.

### Phase 2: Density-v0 model

**Goal:** Implement the first orthogonal gDNA density predictor active in all four use
cases.

**Scope:**

- New module: `calibration/density_model.py`.
- Start with a robust density-v0:
  - class-stratified pooled mean using unified `(count, L_eff)`;
  - Poisson or negative-binomial upper quantile;
  - reference-level fallback;
  - genome-wide fallback;
  - all-gDNA fallback only when no density anchor exists.
- Keep the NB-GLM from `density_model_plan_v1.md` as the Phase 2B upgrade:

```text
log mean = beta0 + beta1 * log L_eff
log kappa = gamma0 + gamma1 * log L_eff
```

- Train from anchors by family:
  - no-transcript-strand/intergenic contained anchors;
  - intron-only anchors;
  - exon-intron and exon-intergenic boundary anchors;
  - high-precision strand-deconvolved rows when strand is active;
  - no-RNA exon self-training rows when available.
- Exclude or downweight:
  - `TS_AMBIG` prediction targets;
  - low-precision strand rows;
  - kappa fallback rows;
  - strongly spliced RNA rows unless accepted as no-RNA exons.

**Tests:**

- Density predicts nonzero upper counts for `TS_AMBIG` rows when anchors exist.
- Prediction is capped by observed compatible mass when used as a local prior count.
- Boundary-only anchor data can train an exonic/boundary density when contained
  intergenic anchors are depleted.
- Missing class anchors fall back to reference/global before all-gDNA fallback.

**Acceptance:**

- `DensityPrediction` exists per region with mean, upper, quality, class, and flags.
- The density model is active even when strand is off.
- The summary compares contained-anchor density to boundary-anchor density; large
  disagreement is reported as capture-like evidence, not yet as a final exposure model.

### Phase 3: Region-level integration surface

**Goal:** Fuse strand and density predictions into one region-level gDNA prior surface.

**Scope:**

- New module: `calibration/integration.py`.
- Implement evidence modes from Section 5.
- Produce `IntegratedRegionGdnaSurface`:
  - `gdna_prior_count_upper`;
  - `gdna_expected_count`;
  - `rna_lower_count`;
  - `observed_compatible_count`;
  - `strand_count_used`;
  - `density_count_used`;
  - source flags;
  - training/exposure eligibility flags.
- Use source-gated integration v0.
- Preserve independent strand and density predictions for diagnostics.

**Tests:**

- Unstranded/no-capture fixture routes all non-anchor rows through density.
- Unstranded/capture fixture can use boundary-trained density when contained anchors
  are depleted.
- Strand-specific/no-capture fixture uses strand for single-strand rows and density for
  `TS_AMBIG` rows.
- Strand-specific/capture fixture uses strand for identifiable exons and density for
  ambiguous exons.
- No observed mass is lost because a strand flag is set.

**Acceptance:**

- Every region has exactly one integration mode.
- Sum of integrated prior counts is bounded by observed compatible mass plus explicitly
  flagged model-injected fallback mass. Prefer no model-injected mass in v0.
- Summary reports how much mass came from strand, density, direct anchors, and fallback.

### Phase 4: PriorTable and conservation-aware MultiLocus allocation

**Goal:** Aggregate the integrated region surface to the EM unit, `MultiLocus`.

**Scope:**

- New module: `calibration/prior.py`.
- Define `PriorTable` with:
  - `gdna_prior_count`;
  - `gdna_prior_count_em`;
  - `gdna_expected_count`;
  - `rna_lower_count`;
  - `observed_count`;
  - `gdna_eff_len_raw`;
  - `gdna_eff_len`;
  - `gdna_em_exposure_weight`;
  - `enable_gdna`;
  - rich diagnostics by source.
- Implement region-to-MultiLocus overlap allocation:

```text
raw_share_mr = overlap_bp_mr / region_length_r
scale_r = 1 if sum_m raw_share_mr <= 1 else 1 / sum_m raw_share_mr
share_mr = raw_share_mr * scale_r
```

- Compute `gdna_eff_len_raw` using `_exposure.gdna_eff_len_for_loci`.
- Compute exposure-weighted denominator using
  `_exposure.bp_weighted_mean_exposure_over_blocks` and current `RegionExposure`.
- Do not define `enable_gdna` as `gdna_prior_count > 0`; use geometry and later EM
  finite-likelihood evidence.

**Tests:**

- Region contributions are conserved across overlapping MultiLoci.
- Partial overlaps prorate counts.
- Multiple `Locus` intervals in one `MultiLocus` do not double-count their own
  footprint.
- Uniform exposure gives `gdna_eff_len == gdna_eff_len_raw`.
- Zero prior does not force `enable_gdna = 0` when geometry is valid.

**Acceptance:**

- Prior table exists and is independent of pipeline EM wiring.
- Diagnostics explain large gDNA priors by source.

### Phase 5: Pipeline wiring to EM

**Goal:** Remove the `quant_from_buffer` boundary and run locus EM end-to-end.

**Scope:**

- Build loci and partitions as the pipeline already does.
- Call prior assembly after locus construction.
- Pass to `_run_locus_em_partitioned`:
  - `prior_table.gdna_prior_count_em`;
  - `prior_table.gdna_eff_len`;
  - `prior_table.enable_gdna`;
  - `prior_table.gdna_eff_len_raw`;
  - `prior_table.gdna_prior_count`;
  - `prior_table.gdna_em_exposure_weight`.
- Add prior diagnostics to `loci.feather` or summary/locus metadata as appropriate.

**Tests:**

- Tiny unstranded/no-capture pipeline run completes.
- Tiny strand-specific/no-capture pipeline run completes.
- Ambiguous opposite-strand fixture routes ambiguous rows through density, then EM.
- No `NotImplementedError("rigel quant: locus EM lands in Phase 6")` remains in
  production flow.

**Acceptance:**

- `rigel quant` runs end-to-end with uniform exposure.
- Density and strand source summaries are written.
- Golden output updates are intentional and explained.

### Phase 6: Bayesian strand-density fusion

**Goal:** Replace source-gated strand/density combination with a real fused posterior
where both evidence channels are available.

**Scope:**

- Extend `deconvolve_regions_by_strand` or add a sibling function that accepts a
  density prior over `D`.
- For exact regions:

```text
posterior(D) proportional to strand_likelihood(D) * density_prior(D)
```

- For approximate regions, combine moments or use a stable numerical approximation.
- Keep source-gated v0 as fallback when density prior is unavailable or invalid.

**Tests:**

- With a flat density prior, fused output equals current strand-only output.
- With strong density prior and weak strand evidence, fused output moves toward density.
- With strong strand evidence and broad density prior, fused output remains strand-led.
- `TS_AMBIG` remains density-led because the strand likelihood is underdetermined.

**Acceptance:**

- Case 3 synergy is real, not just side-by-side diagnostics.
- Case 4 identifiable exons remain strand-powered while ambiguous regions remain
  density-powered.

### Phase 7: Unsupervised exposure learning

**Goal:** Learn the hybrid-capture accessibility surface `A_r`.

**Scope:**

- Replace uniform `RegionExposure` with an unsupervised model.
- Use density and strand outputs correctly:
  - use `D_mean` for exposure learning;
  - use `D_upper` for conservative EM prior counts;
  - downweight low-precision/fallback rows;
  - retain all rows in diagnostics.
- Learn relative exposure from residuals against density opportunity:

```text
D_mean_r approximately eta_g * A_r * L_eff_r
```

- Boundary-trained density is essential for capture libraries because background
  intergenic/intronic contained anchors can be depleted.
- Preserve count closure diagnostics:

```text
sum predicted D_mean_r approximately sum deconvolved/training D_mean_r
```

**Tests:**

- Synthetic no-capture data learns near-uniform exposure.
- Synthetic capture data learns higher exposure on target/boundary-proximal regions.
- Strand-specific capture exons use strand counts while exposure corrects denominator.
- Ambiguous captured regions use density plus exposure.

**Acceptance:**

- Exposure summary reports target-like concentration without requiring a BED input.
- Captured loci get accessible-scale denominators, not raw genomic-span denominators.
- No regression in no-capture datasets.

### Phase 8: Validation matrix

**Goal:** Validate the integrated calibration system across the four use cases.

**Matrix:**

```text
1. unstranded, no capture
2. unstranded, capture
3. strand-specific, no capture
4. strand-specific, capture
```

For each case report:

- strand identifiable diagnostics;
- density anchor family contributions;
- boundary vs contained density agreement;
- exposure mode and learned exposure summaries;
- prior mass by source;
- gDNA/RNA assignment error on truth data;
- ambiguous-region behavior;
- failure modes and fallback rates.

**Acceptance:**

- Strand-off cases still estimate gDNA using density.
- Capture cases do not rely on depleted intergenic/intronic contained regions alone.
- Strand-specific cases use strand where identifiable.
- Ambiguous regions are never silently dropped and never silently called all-gDNA
  without a diagnostic fallback flag.

## 8. First Priority: What To Do Next

Start with **Phase 0**.

The next PR should not implement the NB-GLM yet and should not wire EM. It should make
calibration structurally honest:

1. Build a `RegionEvidenceTable` from `RegionArrays`, `PayloadArrays`, current strand
   predictions, and current uniform exposure.
2. Preserve observed totals for rows where strand folding is impossible.
3. Add mode/capability diagnostics to `summary.json`.
4. Keep density prediction fields empty but typed/owned.
5. Add tests proving no observed mass disappears for `TS_NONE` and `TS_AMBIG` rows.

Once this exists, density model implementation becomes a plug-in to a known surface,
not another parallel path.

## 9. Design Rules To Preserve

- Genes are output conveniences only. Calibration, density, prior construction, loci,
  and EM must remain transcript/region/locus based.
- Strand specificity is learned early and can disable strand deconvolution.
- Density is active in all four use cases.
- Boundary flux is first-class density evidence, especially under capture.
- Exposure affects opportunity denominators, not observed count ledgers.
- `upper_count` is the conservative EM count surface; `mean_count` is the statistical
  expectation and exposure-learning surface.
- Unidentifiable does not mean zero.
- Ambiguous opposite-strand regions require density or an explicit fallback.
- Prior assembly must conserve observed region mass across overlapping MultiLoci.
- The EM contract should not care whether a prior count came from strand, density, or
  direct anchors; diagnostics should care deeply.

## 10. Open Questions

1. In source-gated integration v0, should low-precision rows use
   `max(strand_upper, density_upper)` or density upper alone with a separate fallback
   threshold? The conservative recommendation is `max` until validation says it is too
   aggressive.
2. Should density-v0 start pooled/class-stratified or jump directly to NB-GLM? The
   roadmap recommends pooled/class-stratified first, then NB-GLM as Phase 2B, because
   integration and diagnostics are the larger risk.
3. How should boundary anchor families be selected before exposure learning? The
   roadmap recommends fitting/reporting both contained-anchor and boundary-anchor
   density families, then letting validation decide the default predictor per region
   class.
4. When should density become a prior inside strand deconvolution rather than a
   post-hoc integration surface? Phase 6, after end-to-end uniform-exposure EM runs.
