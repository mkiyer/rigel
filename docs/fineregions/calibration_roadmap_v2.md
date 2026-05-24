# Integrated Calibration Roadmap v2

Date: 2026-05-24
Status: roadmap document only
Supersedes/updates: `calibration_roadmap_v1.md`
Inputs synthesized: `strand_model_impl_plan_v5.md`, `strand_model_impl_phase5.md`, `density_model_plan_v1.md`

## 0. What Changed From v1

v1 correctly identified that strand and density are orthogonal evidence channels, but it
pushed too much policy into per-region evidence types. That was inelegant and partly
redundant:

- zero-observed status is derived from count arrays, and there are multiple relevant
  count surfaces: contained, boundary-left, boundary-right, spliced, unspliced, and
  all-compartment totals;
- transcript-strand applicability is derived from the fine-region `signature` and
  `ts_class`, so it should not be duplicated as a flag;
- near-unstranded behavior is a library-level context, not a region category;
- low precision and fallback are model-quality diagnostics, better represented as
  reliability/variance/quality scores than as a growing enum.

v2 replaces the region-type taxonomy with a cleaner decomposition:

```text
Region facts        = coordinates, signature, ts_class, region class
Count ledgers       = contained/boundary/spliced/unspliced POS/NEG observed mass
Library context     = strand specificity, strand strength, exposure mode, density mode
Model evidence      = strand predictive distribution + density predictive distribution
Integration layer   = fused posterior over regional gDNA count D
Prior assembly      = conservation-aware aggregation to MultiLoci
```

The key design rule is now:

```text
Do not encode library-level or signature-derived facts as per-region policy flags.
Use them to compute model applicability and model reliability.
```

## 1. First Principles

For each fine region `r`, calibration sees a ledger of observed fragments. It needs to
infer a distribution over the local gDNA count `D_r` and eventually aggregate those
regional distributions into per-MultiLocus EM priors.

The two evidence channels answer different questions:

```text
Strand channel:
  Given local POS/NEG counts and library strand specificity, which values of D_r are
  compatible with the observed strand imbalance?

Density channel:
  Given region opportunity, region class, boundary flux, and exposure state, which
  values of D_r are plausible before looking at local strand imbalance?
```

The natural fusion target is therefore:

```text
P(D_r | density, strand, observed ledger)
  proportional to P_density(D_r | L_eff_r, class_r, exposure_r)
               * P_strand(strand_counts_r | D_r, ledger_r, SS, kappa_d)
```

When strand specificity is weak, `P_strand` becomes broad. It should naturally receive
low weight because its likelihood is flat, not because every region carries an
`UNSTRANDED` flag.

When strand is inapplicable, such as no transcript-relative strand axis or opposite-
strand ambiguity, `P_strand` is neutral or not constructed. Density then owns the
prediction.

When density is weak, such as poor training support or uncertain capture exposure,
`P_density` is broad. Strand can then dominate in identifiable regions.

## 2. The Four RNA-seq Use Cases

Rigel learns strand specificity early. Exposure/capture state is harder and is learned
later. The density channel remains active in all cases.

### Case 1: unstranded, no hybrid capture

```text
strand likelihood: neutral or extremely broad
density prior: active
exposure: usually uniform
```

The strand channel is OFF because `SS` is close to 0.5. Density drives gDNA estimates
using contained intergenic/intronic anchors plus boundary flux.

Implementation consequence:

- no per-region unstranded flag is needed;
- library context records low strand strength;
- integration treats `strand_precision_r = 0` or near zero;
- fused posterior is density-led.

### Case 2: unstranded, hybrid capture

```text
strand likelihood: neutral or extremely broad
density prior: active, boundary-heavy
exposure: learned in Phase 7
```

Background intergenic/intronic contained regions may be depleted and unrepresentative.
Boundary flux at exon-intron and exon-intergenic edges becomes the essential density
signal because it samples gDNA-like fragments near captured sequence.

Implementation consequence:

- density observation construction must include boundary counts and boundary
  opportunity from the beginning;
- density summaries must compare contained-anchor and boundary-anchor estimates;
- exposure remains uniform until Phase 7, but raw and exposure-weighted denominators
  are kept separate.

### Case 3: strand-specific, no hybrid capture

```text
strand likelihood: active in applicable regions
density prior: active
exposure: usually uniform
```

This is the main synergy case. Density gives a prior over `D_r`; strand counts give a
local likelihood. Weak or moderate strand specificity, for example `SS = 0.6`, should
produce a flatter strand likelihood and therefore a smaller effective fusion weight.
Strong strand specificity produces a sharper likelihood.

Implementation consequence:

- use exact posterior fusion where possible;
- report diagnostic weights such as `w_strand_r` and `w_density_r` from posterior
  precision or information, not from hand-written region categories;
- keep density predictions even when strand is strong, because they are useful for
  diagnostics, ambiguous regions, and exposure learning.

### Case 4: strand-specific, hybrid capture

```text
strand likelihood: active in applicable regions
density prior: active, boundary-aware
exposure: learned in Phase 7
```

Strand-specific capture data lets entire single-strand exons be deconvolved from local
sense/antisense imbalance without relying on boundary flux. But density is still needed
for strand-ambiguous regions, weak-strand regions, and denominator/exposure
normalization.

Implementation consequence:

- identifiable exons can be strand-powered;
- opposite-strand ambiguous regions are density-powered unless a future transcript-level
  expression prior can make strand informative;
- initial density must train from exon-adjacent boundary flux because depleted
  background anchors can understate on-target gDNA density;
- Phase 7 uses `D_mean` surfaces to learn exposure and uses conservative upper surfaces
  for EM priors.

## 3. Current Implementation Inventory

Implemented and useful today:

- `scan_payload.py`: fractional 12-channel payload, including contained and boundary
  left/right unspliced POS/NEG counts.
- `_arrays.py`: `RegionArrays` and `PayloadArrays`, the shared sorted region/count
  view for density and strand.
- `strand_summary.py`: library-level strand summary, `p_r1_sense`, signed contrast,
  contrast margin, and numerical floor.
- `fl.py`: RNA/gDNA/global fragment-length models.
- `density_global.py`: contained-only global density diagnostics,
  `l_eff_contained`, and gDNA strand-balance kappa estimation.
- `_exposure.py`: FL-aware geometry helpers, including boundary-side opportunity and
  locus gDNA effective length.
- `strand_deconv.py`: current strand-only regional deconvolution producing
  `RegionGdnaEstimate.mean_count`, `upper_count`, `rna_lower_count`, precision, and
  flags.
- `exposure.py`: `RegionExposure.uniform` scaffold.
- `_orchestrator.py`: wires FL, contained densities, strand counts, kappa_d, strand
  deconvolution, and uniform exposure into `CalibrationResult`.
- `pipeline.py`: lower-level EM path already accepts prior counts, effective lengths,
  raw denominators, exposure weights, and `enable_gdna`, but `quant_from_buffer` still
  stops before prior assembly.

Missing:

- a unified evidence ledger that preserves observed mass independent of strand folding;
- density observations with boundary flux and FL-aware boundary opportunity;
- density predictive distributions;
- a fusion layer that combines density prior and strand likelihood;
- MultiLocus prior aggregation;
- end-to-end EM wiring;
- unsupervised exposure learning;
- validation across all four use cases.

## 4. Evidence Schema, Without Redundant Flags

### 4.1 Region facts

These are intrinsic and already derivable from the index:

```python
@dataclass(frozen=True, slots=True)
class RegionFacts:
    ref_id: np.ndarray
    start: np.ndarray
    end: np.ndarray
    signature: np.ndarray
    ts_class: np.ndarray
```

Do not add flags like `NO_TRANSCRIPT_STRAND` or `STRAND_AMBIGUOUS` as policy state.
Those are derived by helper functions from `signature` and `ts_class` when needed.

### 4.2 Count ledger

This is the most important Phase 0 object. It prevents false zeros when strand folding
is not possible.

```python
@dataclass(frozen=True, slots=True)
class RegionCountLedger:
    contained_unspliced_pos: np.ndarray
    contained_unspliced_neg: np.ndarray
    boundary_left_unspliced_pos: np.ndarray
    boundary_left_unspliced_neg: np.ndarray
    boundary_right_unspliced_pos: np.ndarray
    boundary_right_unspliced_neg: np.ndarray
    spliced_pos: np.ndarray
    spliced_neg: np.ndarray

    contained_unspliced_total: np.ndarray
    boundary_unspliced_total: np.ndarray
    total_unspliced: np.ndarray
    total_spliced: np.ndarray
    total_observed: np.ndarray
```

Zero-observed status is not a mode. It is a query:

```text
total_observed_r == 0
contained_unspliced_total_r == 0
boundary_unspliced_total_r == 0
total_unspliced_r == 0
```

Different callers care about different zeros. Density cares about gDNA-compatible
unspliced/boundary counts. Strand deconvolution cares about the strand-foldable ledger.
Diagnostics can report all of them without making them evidence types.

### 4.3 Library context

Library-level decisions live once, not once per region.

```python
@dataclass(frozen=True, slots=True)
class CalibrationContext:
    p_r1_sense: float
    strand_specificity: float
    signed_strand_contrast: float
    signed_strand_contrast_margin_99: float
    strand_strength: float        # [0, 1], model reliability scale
    kappa_d: float
    kappa_d_fallback_used: bool
    density_mode: str             # none, pooled_v0, nb_glm, ...
    exposure_mode: str            # uniform, unsupervised
```

`strand_strength` should be derived from contrast and uncertainty, for example:

```text
raw_contrast = abs(2 * p_r1_sense - 1)
noise_floor = max(contrast_margin_99, min_practical_contrast)
strand_strength = clip((raw_contrast - noise_floor) / (1 - noise_floor), 0, 1)
```

The exact formula can be calibrated later. The key point is that weak strand specificity
is a library context that changes the strand likelihood or its precision. It is not a
per-region flag.

### 4.4 Model evidence surfaces

Each model returns a distributional statement and a reliability statement.

```python
@dataclass(frozen=True, slots=True)
class StrandEvidence:
    applicable: np.ndarray        # derived from ts_class and context
    mean: np.ndarray              # strand-only E[D]
    upper: np.ndarray             # strand-only conservative D upper
    variance: np.ndarray          # strand-only uncertainty over D
    precision: np.ndarray         # 1 / variance, clipped and finite
    log_likelihood_shape: str     # exact, normal_approx, neutral
    quality: np.ndarray           # diagnostic [0, 1]
```

```python
@dataclass(frozen=True, slots=True)
class DensityEvidence:
    applicable: np.ndarray
    mean: np.ndarray
    upper: np.ndarray
    variance: np.ndarray
    precision: np.ndarray
    leff: np.ndarray
    anchor_support: np.ndarray
    quality: np.ndarray
    family: np.ndarray            # contained, boundary, mixed, fallback
```

Flags can still exist at the low-level model boundary for debugging, but integration
should consume applicability, moments/distributions, and quality. The primary API
should be continuous reliability, not a taxonomy of cases.

## 5. Fusion Design

### 5.1 Exact fusion target

For a region with a single transcript-relative strand axis:

```text
P(D | data) proportional to P_density(D | L_eff, class, exposure)
                       * P_strand(k_sense, k_antisense | D, N, p_r1_sense, kappa_d)
```

Where:

- density supplies the prior over `D`;
- strand supplies the local likelihood;
- `D` is bounded to `[0, observed_compatible_count]`;
- `SS` enters through `p_r1_sense` and therefore changes likelihood curvature.

If `SS = 0.6`, the likelihood is broad. If `SS = 0.99`, it is sharp. No special
per-region unstranded flag is needed.

### 5.2 Approximate fusion

When exact discrete fusion is too expensive, use moment/precision fusion as an
approximation:

```text
tau_density = 1 / var_density
tau_strand  = 1 / var_strand
mean_fused  = (tau_density * mean_density + tau_strand * mean_strand)
              / (tau_density + tau_strand)
var_fused   = 1 / (tau_density + tau_strand)
```

Diagnostic channel weights are then:

```text
w_density = tau_density / (tau_density + tau_strand)
w_strand  = tau_strand  / (tau_density + tau_strand)
```

This is the principled version of the earlier inverse-variance weighting idea. The
important improvement is that variances come from the actual density predictive model
and the actual strand likelihood curvature, not directly from a single library-level
`SS` scalar. `SS` influences `var_strand`; it does not directly choose a policy branch.

### 5.3 Applicability handling

Applicability is derived, not flagged:

```text
strand applicable:
  ts_class in {TS_POS, TS_NEG}
  and strand_strength > 0
  and local strand-folded count > 0

density applicable:
  density model has a valid prediction for the region class/reference/fallback path
```

For `TS_NONE`, there is no transcript-relative RNA strand axis. These regions can be
excellent density anchors, but strand likelihood is neutral.

For `TS_AMBIG`, the two observed strand channels are insufficient to identify total
`D` without additional assumptions about POS-strand RNA and NEG-strand RNA. Density is
therefore the main channel. A future transcript-expression prior could make strand
information useful here, but Phase 0-5 should not pretend it is identified.

For zero observed mass, both likelihoods are trivial and the fused count is zero by
ledger bounds.

### 5.4 Fusion output

```python
@dataclass(frozen=True, slots=True)
class FusedRegionGdnaEvidence:
    mean: np.ndarray              # E[D | density, strand]
    upper: np.ndarray             # confidence-controlled conservative count
    rna_lower: np.ndarray         # observed_compatible_count - upper
    variance: np.ndarray
    w_density: np.ndarray
    w_strand: np.ndarray
    observed_compatible_count: np.ndarray
    density_applicable: np.ndarray
    strand_applicable: np.ndarray
    quality: np.ndarray
```

`upper` is the EM-facing conservative count surface. `mean` is the statistical
expectation and the main exposure-learning surface.

## 6. Density Model Role

Density is active in all four cases. It should be implemented as a predictive
distribution, not just a point estimate.

### 6.1 Density observations

Per region:

```text
N_density = contained_unspliced_total
          + boundary_left_unspliced_total
          + boundary_right_unspliced_total

L_density = contained_effective_length
          + boundary_left_effective_opportunity
          + boundary_right_effective_opportunity
```

Boundary opportunity uses `_exposure.fractional_boundary_side_exposure` and is
especially important for hybrid capture.

### 6.2 Density-v0

Start with class/reference stratified pooled predictors:

```text
rho_class = sum N_density / sum L_density
mean_r    = rho_class * L_density_r
upper_r   = count_quantile(mean_r, dispersion_class, alpha)
```

Use fallback order:

```text
class + reference
class genome-wide
all anchors on reference
all anchors genome-wide
broad all-gDNA fallback, diagnostic only
```

### 6.3 Density quality

Density evidence should report quality continuously:

```text
anchor_support_r
fallback_depth_r
prediction_variance_r
exposure_uncertainty_r
quality_r in [0, 1]
```

The integration layer uses variance/precision. Human-facing diagnostics can summarize
quality and fallback depth.

## 7. Strand Model Role

The current strand model is already useful, but it should expose likelihood/reliability
more explicitly for fusion.

Needed additions:

- preserve current `RegionGdnaEstimate` for backward-compatible summaries;
- expose or reconstruct `var_strand_r` for each region;
- expose `strand_applicable_r` as a derived mask, not a new region type;
- include library `strand_strength` in the calibration summary;
- for weak `SS`, produce broad strand evidence rather than a hard ON/OFF branch.

A hard numerical floor is still useful to avoid unstable arithmetic, but production
behavior should be governed by likelihood width and model quality.

## 8. Revised Implementation Roadmap

### Phase 0: Evidence ledger and context synthesis

**Goal:** Consolidate current calibration into a clean evidence system before adding
new model behavior.

**Implement:**

- New module: `calibration/evidence.py`.
- Dataclasses:
  - `RegionFacts` or a thin wrapper around `RegionArrays`;
  - `RegionCountLedger`;
  - `CalibrationContext`;
  - `RegionEvidenceBundle` holding facts, ledger, context, current strand output, and
    uniform exposure.
- Helper functions:
  - `build_region_count_ledger(payload_arrays)`;
  - `build_calibration_context(strand_summary, kappa_d, region_exposure, density_mode)`;
  - `strand_applicability(region_facts, context, ledger)`;
  - `density_anchor_masks(region_facts, ledger)`.
- Orchestrator change: build this bundle and include compact summaries in
  `CalibrationResult`.

**Do not implement yet:** density fitting, prior assembly, EM wiring.

**Tests:**

- `TS_NONE` and `TS_AMBIG` rows retain nonzero observed ledgers even when
  `RegionGdnaEstimate.n_total` is zero.
- Zero status is derivable separately for contained, boundary, unspliced, and total
  ledgers.
- Strand applicability is derived from `ts_class` and context.
- Weak strand specificity reduces `strand_strength` without writing a per-region
  near-unstranded flag.

**Acceptance:**

- Summary reports strand strength, density mode `not_built`, exposure mode `uniform`,
  count-ledger totals, and applicability counts.
- No behavior change after calibration; `quant_from_buffer` still stops before EM.

### Phase 1: Unified density observations

**Goal:** Build density observations from count ledgers and FL-aware opportunity.

**Implement:**

- New module: `calibration/density_observation.py`.
- `DensityObservation` dataclass with `N_density`, `L_density`, contained and boundary
  components, class labels, and anchor masks.
- Boundary opportunity from `_exposure.fractional_boundary_side_exposure`.
- Summary comparing contained-only and boundary-inclusive density totals.

**Tests:**

- Boundary-free data reproduces current contained-only pooled density.
- Boundary counts and boundary opportunity are both included.
- Capture-like fixtures can show depleted contained anchors but nonzero boundary
  anchors.

**Acceptance:** density observations exist for every region and are not tied to strand
model output.

### Phase 2: Density predictive distributions

**Goal:** Implement density as an orthogonal predictive distribution over `D`.

**Implement:**

- New module: `calibration/density_model.py`.
- `DensityEvidence` dataclass with mean, upper, variance, precision, quality, family,
  and fallback depth.
- Start with class/reference pooled density-v0 and Poisson or NB upper quantiles.
- Keep NB-GLM with empirical quantile calibration as Phase 2B, after v0 is validated.

**Tests:**

- Density predicts for `TS_AMBIG` rows using anchors.
- Predictions are capped to local observed compatible mass when converted to local
  prior counts.
- Boundary anchors can drive exonic/capture-proximal density when contained anchors are
  depleted.
- Fallback order is deterministic and diagnosed.

**Acceptance:** density is active when strand is off and produces per-region mean,
upper, variance, and quality.

### Phase 3: Probabilistic density-strand fusion

**Goal:** Integrate the two channels by uncertainty rather than by region-type flags.

**Implement:**

- New module: `calibration/integration.py`.
- `StrandEvidence` wrapper around current strand output, including variance/precision
  and applicability.
- `FusedRegionGdnaEvidence` dataclass.
- Exact discrete fusion for small counts where strand likelihood is applicable.
- Precision/moment fusion fallback for large counts.
- Neutral strand likelihood for inapplicable regions.

**Tests:**

- With flat density prior, fused output matches strand-only output.
- With neutral strand likelihood, fused output matches density output.
- With weak `SS`, `w_strand` is lower than with strong `SS` for the same local counts.
- `TS_AMBIG` does not use a fake strand split; density dominates.
- No observed compatible mass is lost.

**Acceptance:** every region gets a fused mean, upper, variance, and diagnostic
`w_density`/`w_strand`.

### Phase 4: PriorTable and MultiLocus allocation

**Goal:** Aggregate fused regional evidence to MultiLocus EM inputs.

**Implement:**

- New module: `calibration/prior.py`.
- `PriorTable` with raw and EM-facing prior counts, expected counts, RNA lower counts,
  raw and exposure-weighted effective lengths, exposure weights, `enable_gdna`, and
  source/weight diagnostics.
- Conservation-aware region-to-MultiLocus allocation:

```text
raw_share_mr = overlap_bp_mr / region_length_r
share_mr = raw_share_mr / max(1, sum_m raw_share_mr)
```

- Uniform exposure path should yield `gdna_eff_len == gdna_eff_len_raw`.

**Tests:**

- Partial overlaps prorate counts.
- Overlapping MultiLoci conserve region mass.
- Multiple loci inside one MultiLocus do not double-count the same footprint.
- Zero prior does not disable gDNA when geometry and likelihood support exist.

**Acceptance:** prior table exists independently of EM wiring and explains its inputs.

### Phase 5: Pipeline wiring to EM

**Goal:** Remove the current pipeline boundary and run end-to-end with uniform exposure.

**Implement:**

- Build fused region evidence during calibration.
- Build `PriorTable` after MultiLocus construction.
- Pass `gdna_prior_count_em`, `gdna_eff_len`, `enable_gdna`, raw prior count,
  unweighted denominator, and exposure weight into `_run_locus_em_partitioned`.

**Tests:**

- End-to-end unstranded/no-capture fixture.
- End-to-end strand-specific/no-capture fixture.
- Opposite-strand ambiguous fixture shows density-dominated ambiguous regions.

**Acceptance:** `rigel quant` runs end-to-end before unsupervised exposure lands.

### Phase 6: Density refinement

**Goal:** Improve density-v0 after the integrated system is runnable.

**Implement:**

- NB-GLM density model:

```text
log mean = beta0 + beta1 * log L_eff
log kappa = gamma0 + gamma1 * log L_eff
```

- Empirical quantile calibration by log-length bins.
- Optional covariates later: GC, mappability, signature class.

**Tests:**

- NB-GLM recovers pooled density when slope is constrained to 1 and dispersion is
  constant.
- Quantile calibration improves empirical coverage.
- Predictions remain monotone and nonnegative.

**Acceptance:** density quality improves without changing downstream interfaces.

### Phase 7: Unsupervised exposure learning

**Goal:** Learn the hybrid-capture exposure surface `A_r`.

**Implement:**

- Replace uniform `RegionExposure` with an unsupervised exposure estimate.
- Use fused or strand/density `mean` surfaces for exposure learning.
- Use conservative `upper` surfaces for EM priors.
- Keep count closure diagnostics:

```text
sum predicted D_mean approximately sum observed/fused D_mean over training rows
```

**Tests:**

- No-capture synthetic data learns near-uniform exposure.
- Capture synthetic data learns higher exposure over target/boundary-proximal regions.
- Strand-specific capture still deconvolves identifiable exons locally.
- Ambiguous capture regions use density plus exposure.

**Acceptance:** captured loci get accessible-scale denominators without requiring BED
input.

### Phase 8: Validation matrix

Validate across:

```text
1. unstranded, no capture
2. unstranded, capture
3. strand-specific, no capture
4. strand-specific, capture
```

Report:

- strand strength and applicability;
- density anchor families and quality;
- boundary vs contained density disagreement;
- exposure mode and learned exposure distribution;
- fused `w_density` and `w_strand` distributions;
- prior mass and EM gDNA/RNA assignment error;
- fallback rates.

## 9. First Priority

Start with **Phase 0**.

The next implementation step is not density fitting and not EM wiring. It is to create
the evidence ledger and library context so every later model reads the same facts:

1. Preserve observed ledgers for every region, including rows where strand folding is
   impossible.
2. Derive applicability from `signature`, `ts_class`, count ledgers, and library
   context.
3. Represent library strand state once as `CalibrationContext`, including continuous
   `strand_strength`.
4. Keep model uncertainty as variance/precision/quality, not as redundant region-type
   flags.
5. Add summary diagnostics that show what the calibration system can and cannot infer
   before prior assembly.

## 10. Design Rules

- Do not duplicate `signature` and `ts_class` as policy flags.
- Do not encode near-unstranded as a per-region category.
- Do not make zero-observed a region type; derive it from the relevant count ledger.
- Use flags for low-level debugging, not as the primary integration API.
- Integrate density and strand through distributions, variance, and likelihood
  curvature.
- Let weak strand specificity reduce strand precision naturally.
- Density remains active in all four use cases.
- Boundary flux is first-class density evidence.
- Exposure changes opportunity denominators, not observed count ledgers.
- `upper` is the conservative EM count surface; `mean` is the expectation and exposure
  learning surface.
- Ambiguous opposite-strand regions need density or a future transcript-level prior;
  they should not receive a fake strand deconvolution.

## 11. Open Questions

1. What formula should map library strand contrast and uncertainty to
   `strand_strength`? Start simple and calibrate empirically.
2. Should Phase 3 exact fusion use a density prior over `D` directly, or a truncated
   density prior conditioned on local observed compatible mass? The conservative choice
   is to truncate to the local observed ledger before EM prior construction.
3. For weak but nonzero strand specificity, when is exact strand likelihood still worth
   computing? The answer may be count-dependent, because high-depth regions can carry
   information even at moderate `SS`.
4. Can `TS_AMBIG` regions eventually use strand information by adding a transcript-side
   RNA prior? Not in the near-term roadmap, but the fusion framework leaves room for it.
