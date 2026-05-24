# Integrated Calibration Roadmap v3

Date: 2026-05-24
Status: roadmap document only
Supersedes/updates: `calibration_roadmap_v2.md`
Inputs synthesized: `strand_model_impl_plan_v5.md`, `strand_model_impl_phase5.md`,
`density_model_plan_v1.md`, and review feedback on `calibration_roadmap_v2.md`

## 0. Executive Summary

v2 made the right architectural move: calibration should not be a taxonomy of region
flags. Strand and density are continuous, overlapping evidence channels, and the
integration layer should combine their uncertainty instead of choosing hard modes.

v3 keeps that architecture but tightens three things before implementation:

1. The strand model likelihood must own the effect of weak strand specificity. A
   `strand_strength` scalar may be reported as a diagnostic, but it must not multiply
   or double-discount the strand likelihood.
2. Fusion must be mathematically explicit. Exact discrete fusion is the default for
   small observed counts. Large-count approximation starts simple with a bounded
   Laplace/truncated-normal approximation to the fused log posterior, not a free-form
   Gaussian average that can escape the count bounds.
3. Prior allocation to MultiLoci must not smear regional count mass by base-pair
  overlap when observed reads are spatially clustered. Phase 5 now requires a
   footprint-aware allocation path before EM wiring.

The implementation order changes accordingly:

```text
Phase 0: design locks before plumbing
Phase 1: memory-light evidence ledger and calibration context
Phase 2: unified density observations
Phase 3: density predictive distributions
Phase 4: probabilistic density-strand fusion
Phase 5: footprint-aware PriorTable and MultiLocus allocation
Phase 6: pipeline wiring to EM
Phase 7: density refinement
Phase 8: unsupervised exposure learning
Phase 9: validation and deferred expression-aware TS_AMBIG work
```

## 1. Core Model

For each fine region `r`, define:

```text
N_r = observed gDNA-compatible count in the relevant local ledger
D_r = latent observed gDNA count assigned within the region
0 <= D_r <= N_r
```

Density and strand answer different questions:

```text
Density channel:
  What distribution over D_r is plausible from effective opportunity, region class,
  boundary flux, and exposure state?

Strand channel:
  Given N_r and the local POS/NEG count split, which values of D_r are compatible
  with the library strand model?
```

The fused target is:

```text
p(D_r | data_r, 0 <= D_r <= N_r)
  proportional to p_density(D_r | L_eff_r, class_r, exposure_r)
             * L_strand(counts_r | D_r, N_r, p_r1_sense, kappa_d)
             * 1[0 <= D_r <= N_r]
```

This formulation deliberately separates three things:

```text
unbounded density prediction       -> model expectation before local count bounds
bounded local posterior over D_r   -> EM-facing regional prior count surface
tail/tension diagnostics           -> evidence that density/exposure is mismatched
```

A local EM prior count cannot exceed observed compatible mass. But the density model's
unbounded prediction should not be silently overwritten. If `p_density` places large
mass above `N_r`, the bounded posterior is still used for EM, while the lost tail is
reported as model tension and used to improve density or exposure modeling.

## 2. Four RNA-seq Use Cases

Rigel learns library strand state early. Exposure/capture state is learned later. The
density channel remains active in all cases.

### Case 1: unstranded, no hybrid capture

```text
strand likelihood: neutral or very broad
density prior: active
exposure: usually uniform
```

No per-region unstranded flag is needed. The strand likelihood is naturally flat
because `p_r1_sense` is close to the gDNA strand balance. Density drives the posterior.

### Case 2: unstranded, hybrid capture

```text
strand likelihood: neutral or very broad
density prior: active, boundary-heavy
exposure: learned later
```

Contained background anchors may be depleted by capture. Boundary flux near targeted
sequence becomes first-class density evidence. Uniform exposure remains the initial
scaffold, but raw and exposure-weighted denominators must remain separate.

### Case 3: strand-specific, no hybrid capture

```text
strand likelihood: active where transcript-relative strand is identifiable
density prior: active
exposure: usually uniform
```

This is the clean synergy case. Density supplies the prior over `D_r`; strand counts
supply local likelihood curvature. Weak but nonzero strand specificity still carries
information at sufficient depth because the likelihood includes both `SS` and `N_r`.

### Case 4: strand-specific, hybrid capture

```text
strand likelihood: active in identifiable regions
density prior: active, boundary-aware
exposure: learned later
```

Strand-specific capture can deconvolve identifiable single-strand exons. Density is
still required for opposite-strand ambiguity, weak strand evidence, denominator
normalization, and exposure learning.

## 3. Implemented Inventory

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

- a memory-light observed ledger that preserves mass independent of strand folding;
- formal exact and approximate fusion math;
- density observations with boundary flux and FL-aware boundary opportunity;
- density predictive distributions;
- footprint-aware MultiLocus prior allocation;
- end-to-end EM wiring;
- unsupervised exposure learning with bias/exposure separation;
- validation across all four use cases.

## 4. Design Locks Before Implementation

The following decisions must be locked before writing the Phase 1 `evidence.py`
plumbing.

### 4.1 `strand_strength` is diagnostic only

Do not use a scalar `strand_strength` to scale the strand likelihood or strand
precision. That would double-discount weak strand libraries because weak `SS` already
flattens the likelihood.

The model uses:

```text
p_r1_sense
kappa_d or the gDNA strand-balance model
local strand counts
N_r
```

The likelihood curvature then determines strand information naturally.

A summary metric may still be useful:

```text
strand_contrast = abs(p_r1_sense - p_gdna_sense)
strand_information_diagnostic = normalized Fisher information under typical N
```

But this is for reporting and QA only. It is not an input multiplier to fusion.

### 4.2 Exact fusion threshold

Start simple and explicit:

```text
If N_r <= 50:
  compute exact discrete posterior over d = 0, ..., N_r

If N_r > 50:
  use the large-count approximation, unless diagnostics force exact/adaptive fallback
```

The threshold can be configurable later, but the first implementation should hard-code
or centralize a single default so tests are deterministic.

### 4.3 Large-count approximation

The first approximation should be a bounded Laplace approximation to the fused log
posterior, not a direct Gaussian average of two independent means.

Define:

```text
ell(d) = log p_density(d) + log L_strand(counts | d)
```

Find the bounded mode:

```text
d_hat = argmax ell(d), d in [0, N_r]
```

Approximate locally:

```text
D_r | data approx Normal(d_hat, sigma2_hat) truncated to [0, N_r]
sigma2_hat = 1 / max(epsilon, -ell_second_derivative(d_hat))
```

Then compute:

```text
mean_fused  = E[D_r | 0 <= D_r <= N_r]
upper_fused = conditional upper quantile over [0, N_r]
var_fused   = Var[D_r | 0 <= D_r <= N_r]
```

This preserves non-negativity and the observed-count bound. It also keeps the first
implementation concise. Do not add Gamma or log-normal fusion in Phase 4 unless this
bounded Laplace approximation fails validation.

Fallback conditions that should trigger exact/adaptive treatment:

```text
mode at 0 or N_r
sigma_hat is large relative to N_r
posterior quality flags indicate strong skew or multimodality
density tail above N_r is high
```

For v3, this is a design requirement, not a fully specified numerical optimizer. The
implementation phase should write the exact formulas against the concrete density and
strand likelihoods before coding.

### 4.4 Truncation and tail accounting

Do not truncate the density prior before fusion as though the tail never existed.

Use the unbounded density model in the fused posterior expression, then condition the
posterior on the physically possible local event `0 <= D_r <= N_r` for EM-facing
counts:

```text
Z_in   = integral_or_sum p_density(d) * L_strand(d) over [0, N_r]
Z_tail = integral_or_sum p_density(d) * L_strand(d) over (N_r, infinity)
posterior_for_EM(d) = p_density(d) * L_strand(d) / Z_in, d in [0, N_r]
```

If `L_strand(d)` is not defined for `d > N_r`, use density tail diagnostics:

```text
p_density_tail = P_density(D_r > N_r)
expected_tail  = E_density[(D_r - N_r)+]
```

Record tail/tension metrics rather than silently losing them:

```text
tail_probability_r
expected_tail_count_r
density_over_observed_ratio_r
```

The bounded posterior supplies local EM priors. Tail diagnostics inform density,
exposure, and model-fit reporting. They should not be fed back as observed count mass.

### 4.5 Memory layout for ledgers

Do not implement `RegionCountLedger` as many persistent float64 derived arrays.

Phase 1 should use a columnar structure-of-arrays layout:

```python
@dataclass(frozen=True, slots=True)
class RegionCountLedger:
    contained_pos: np.ndarray       # float32
    contained_neg: np.ndarray       # float32
    boundary_left_pos: np.ndarray   # float32
    boundary_left_neg: np.ndarray   # float32
    boundary_right_pos: np.ndarray  # float32
    boundary_right_neg: np.ndarray  # float32
    spliced_pos: np.ndarray         # float32
    spliced_neg: np.ndarray         # float32
```

Derived totals are methods, not stored columns:

```python
def contained_total(out: np.ndarray | None = None) -> np.ndarray: ...
def boundary_total(out: np.ndarray | None = None) -> np.ndarray: ...
def unspliced_total(out: np.ndarray | None = None) -> np.ndarray: ...
def observed_total(out: np.ndarray | None = None) -> np.ndarray: ...
```

Rules:

- store raw fractional counts as `float32` unless tests prove precision loss;
- perform global reductions in `float64`;
- allocate derived totals into caller-provided scratch buffers when possible;
- do not create per-region Python objects;
- keep masks as derived `np.bool_` arrays or packed diagnostics, not persistent flag
  matrices;
- add a memory-footprint test or diagnostic for million-region scale.

This keeps the initial ledger around 32 bytes per region for the eight raw channels,
excluding existing region metadata and transient scratch arrays.

### 4.6 Footprint-aware prior allocation

Do not allocate observed count mass to MultiLoci by base-pair overlap alone unless the
full fine region is covered or the ledger is known to be uniform at that resolution.

Base-pair prorating is dangerous when a large fine region contains a small read-depth
hotspot. The roadmap therefore requires one of these allocation substrates before
Phase 5 EM wiring:

```text
Preferred:
  FragmentFootprintLedger with fragment span, strand bucket, fractional weight, and
  region id, intersected against MultiLoci using cgranges or equivalent interval joins.

Acceptable first fallback:
  A compact binned ledger inside large fine regions, with bin size chosen small enough
  that count smearing cannot dominate locus priors.

Restricted geometry-only fallback:
  Base-pair prorating may allocate opportunity denominators, but not observed prior
  count mass, unless tests prove the region is fully covered or count-uniform.
```

Phase 5 must include tests with a clustered read hotspot inside a long fine region.
The MultiLocus that misses the hotspot must not inherit hotspot prior mass.

## 5. Evidence Schema

The evidence schema stays close to v2, with memory and model-boundary changes.

### 5.1 Region facts

Region facts are intrinsic and derived from the index:

```python
@dataclass(frozen=True, slots=True)
class RegionFacts:
    ref_id: np.ndarray
    start: np.ndarray
    end: np.ndarray
    signature: np.ndarray
    ts_class: np.ndarray
```

Do not duplicate `signature` or `ts_class` as policy flags. Helper functions derive
applicability masks when needed.

### 5.2 Count ledger

The count ledger stores only primitive observed channels. Totals are derived.

```python
@dataclass(frozen=True, slots=True)
class RegionCountLedger:
    contained_pos: np.ndarray
    contained_neg: np.ndarray
    boundary_left_pos: np.ndarray
    boundary_left_neg: np.ndarray
    boundary_right_pos: np.ndarray
    boundary_right_neg: np.ndarray
    spliced_pos: np.ndarray
    spliced_neg: np.ndarray
```

Zero-observed status is a query against the relevant derived total:

```text
contained_total_r == 0
boundary_total_r == 0
unspliced_total_r == 0
observed_total_r == 0
```

Different consumers use different ledgers. Density cares about gDNA-compatible
unspliced and boundary mass. Strand cares about strand-foldable mass. Diagnostics can
report all relevant zeros without introducing region types.

### 5.3 Calibration context

Library-level state lives once:

```python
@dataclass(frozen=True, slots=True)
class CalibrationContext:
    p_r1_sense: float
    p_gdna_sense: float
    signed_strand_contrast: float
    signed_strand_contrast_margin_99: float
    strand_information_diagnostic: float
    kappa_d: float
    kappa_d_fallback_used: bool
    density_mode: str
    exposure_mode: str
```

`strand_information_diagnostic` is report-only. The model consumes `p_r1_sense`,
`p_gdna_sense`, and the strand likelihood parameters directly.

### 5.4 Model evidence surfaces

```python
@dataclass(frozen=True, slots=True)
class StrandEvidence:
    applicable: np.ndarray
    mean: np.ndarray
    upper: np.ndarray
    variance: np.ndarray
    precision: np.ndarray
    quality: np.ndarray
    likelihood_family: str
```

```python
@dataclass(frozen=True, slots=True)
class DensityEvidence:
    applicable: np.ndarray
    mean_unbounded: np.ndarray
    upper_unbounded: np.ndarray
    variance_unbounded: np.ndarray
    precision_unbounded: np.ndarray
    leff: np.ndarray
    anchor_support: np.ndarray
    quality: np.ndarray
    family: np.ndarray
    fallback_depth: np.ndarray
```

Density reports unbounded predictions. The integration layer creates bounded fused
outputs for EM.

### 5.5 Fused output

```python
@dataclass(frozen=True, slots=True)
class FusedRegionGdnaEvidence:
    mean_bounded: np.ndarray
    upper_bounded: np.ndarray
    variance_bounded: np.ndarray
    rna_lower: np.ndarray
    observed_compatible_count: np.ndarray
    w_density_info: np.ndarray
    w_strand_info: np.ndarray
    density_applicable: np.ndarray
    strand_applicable: np.ndarray
    tail_probability: np.ndarray
    expected_tail_count: np.ndarray
    quality: np.ndarray
```

`upper_bounded` is the conservative EM count surface. `mean_bounded` is the local
expected observed gDNA count. Density `mean_unbounded` remains available for density
model diagnostics and exposure learning, but exposure learning must respect the bias
controls in Phase 8.

## 6. Strand Likelihood

For a strand-identifiable region, let:

```text
N = observed strand-foldable compatible count
K = observed sense-like count after transcript-relative folding
D = latent gDNA count
R = N - D
p_r = RNA sense probability from p_r1_sense and region orientation
p_d = gDNA sense probability from kappa_d or gDNA strand-balance model
```

A simple starting likelihood is binomial or beta-binomial with mixture probability:

```text
q(D) = (R * p_r + D * p_d) / N
K | D ~ Binomial(N, q(D))
```

or an overdispersed equivalent if current `kappa_d` semantics require it.

Important consequences:

- if `p_r` is close to `p_d`, `q(D)` changes slowly with `D`, so the likelihood is
  broad;
- if `N` is high, even moderate contrast can still carry information;
- no external `strand_strength` multiplier is needed.

`TS_NONE` uses a neutral strand likelihood. `TS_AMBIG` also uses a neutral likelihood
until a post-1.0 expression-aware model is designed.

## 7. Density Model Role

Density is active in all four use cases. It must produce a predictive distribution,
not only a point estimate.

### 7.1 Density observations

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

### 7.2 Density-v0

Start with class/reference pooled predictors:

```text
rho_class = sum N_density / sum L_density
mean_r    = rho_class * L_density_r
upper_r   = count_quantile(mean_r, dispersion_class, alpha)
```

Use deterministic fallback:

```text
class + reference
class genome-wide
all anchors on reference
all anchors genome-wide
broad all-gDNA fallback, diagnostic only
```

### 7.3 Density quality

Density evidence reports:

```text
anchor_support_r
fallback_depth_r
prediction_variance_r
exposure_uncertainty_r
quality_r in [0, 1]
tail_probability_r after local bounding
```

The integration layer uses the distribution. Human-facing summaries report quality and
fallback depth.

## 8. Unsupervised Exposure Learning Guardrails

Exposure learning is deferred until after uniform-exposure EM is runnable. When it
lands, it must not confuse capture accessibility with baseline sequencing biases.

The model should eventually distinguish:

```text
observed gDNA density = baseline_bias_r * capture_exposure_r * opportunity_r
```

where baseline bias may include GC, mappability, local sequence composition, and
alignment artifacts. Capture exposure should represent broad accessibility/enrichment,
not isolated bias hotspots.

A safe bootstrapping strategy:

1. Initialize `A_r = 1`.
2. Fit density-v0 or density-v1 on robust anchors with uniform exposure.
3. Add baseline-bias covariates before or during exposure fitting, starting with GC and
   mappability if available.
4. Estimate exposure from residual broad regional enrichment using high-quality fused
   `D_mean` rows only.
5. Shrink exposure strongly toward 1 and constrain global mean exposure.
6. Smooth or regularize exposure so isolated hotspots are not interpreted as capture
   targets.
7. Refit density with exposure as an offset and iterate only if diagnostics improve.

Rules:

- do not train exposure from low-quality ambiguous regions without strong density or
  strand support;
- do not let the same unexplained hotspot simultaneously define density, exposure, and
  EM prior mass;
- report capture/exposure tension separately from GC/mappability bias;
- BED support remains deferred, but the design should allow an optional BED prior later.

## 9. Revised Implementation Roadmap

### Phase 0: Design locks

**Goal:** Freeze the minimum statistical and engineering contracts before adding
`evidence.py`.

**Decisions to write down in code-adjacent docs/tests:**

- `strand_information_diagnostic` is report-only and never scales likelihoods.
- Exact fusion uses `N_r <= 50` as the first threshold.
- Large-count approximation is bounded Laplace/truncated normal.
- Density predictions are unbounded model outputs; fused outputs are bounded local
  posteriors.
- Tail probability and expected tail count are recorded.
- Ledger stores eight primitive float32 count channels; totals are derived lazily.
- MultiLocus count allocation requires footprint-aware or binned observed mass, not
  naive base-pair smearing.

**Acceptance:** no model code required, but implementation tests can be written against
these decisions as expected contracts.

### Phase 1: Memory-light evidence ledger and context synthesis

**Goal:** Consolidate current calibration facts without changing inference behavior.

**Implement:**

- New module: `calibration/evidence.py`.
- Dataclasses:
  - `RegionFacts` or a thin wrapper around `RegionArrays`;
  - memory-light `RegionCountLedger` with eight primitive channels;
  - `CalibrationContext`;
  - `RegionEvidenceBundle` holding facts, ledger, context, current strand output, and
    uniform exposure.
- Helper functions:
  - `build_region_count_ledger(payload_arrays)`;
  - `build_calibration_context(strand_summary, kappa_d, region_exposure, density_mode)`;
  - `strand_applicability(region_facts, context, ledger)`;
  - `density_anchor_masks(region_facts, ledger)`.
- Summary diagnostics:
  - count-ledger totals from derived methods;
  - strand applicability counts;
  - memory footprint estimate;
  - report-only strand information diagnostic.

**Do not implement yet:** density fitting, fusion, prior assembly, or EM wiring.

**Tests:**

- `TS_NONE` and `TS_AMBIG` rows retain nonzero observed ledgers even when strand-folded
  totals are zero.
- Contained, boundary, unspliced, and total zeros are derivable separately.
- Weak strand specificity changes diagnostics but does not alter ledger construction.
- Ledger dtype and derived totals are as expected.

**Acceptance:** no behavior change after calibration; `quant_from_buffer` still stops
before EM.

### Phase 2: Unified density observations

**Goal:** Build density observations from count ledgers and FL-aware opportunity.

**Implement:**

- New module: `calibration/density_observation.py`.
- `DensityObservation` dataclass with `N_density`, `L_density`, contained components,
  boundary components, class labels, and anchor masks.
- Boundary opportunity from `_exposure.fractional_boundary_side_exposure`.
- Summary comparing contained-only and boundary-inclusive density totals.

**Tests:**

- Boundary-free data reproduces current contained-only pooled density.
- Boundary counts and boundary opportunity are both included.
- Capture-like fixtures show depleted contained anchors but nonzero boundary anchors.

**Acceptance:** density observations exist for every region and are independent of
strand deconvolution output.

### Phase 3: Density predictive distributions

**Goal:** Implement density as an orthogonal predictive distribution over unbounded
regional gDNA count.

**Implement:**

- New module: `calibration/density_model.py`.
- `DensityEvidence` with unbounded mean, upper, variance, precision, quality, family,
  and fallback depth.
- Class/reference pooled density-v0 with Poisson or NB upper quantiles.
- Tail diagnostics after comparing unbounded predictions to local observed compatible
  mass.

**Tests:**

- Density predicts for `TS_AMBIG` rows using anchors.
- Unbounded predictions can exceed local observations without being silently clipped.
- Tail probability is reported when predictions exceed observed compatible mass.
- Boundary anchors can drive exonic/capture-proximal density when contained anchors are
  depleted.
- Fallback order is deterministic and diagnosed.

**Acceptance:** density is active when strand is off and reports a predictive
distribution plus quality diagnostics.

### Phase 4: Probabilistic density-strand fusion

**Goal:** Integrate density and strand through the fused posterior over bounded `D_r`.

**Implement:**

- New module: `calibration/integration.py`.
- `StrandEvidence` wrapper around current strand output, including variance/precision
  derived from the likelihood, not from an external multiplier.
- `FusedRegionGdnaEvidence` with bounded mean, upper, variance, information weights,
  and tail diagnostics.
- Exact discrete fusion for `N_r <= 50`.
- Bounded Laplace/truncated-normal approximation for large counts.
- Neutral strand likelihood for `TS_NONE`, `TS_AMBIG`, and inapplicable regions.

**Tests:**

- With neutral strand likelihood, fused output matches bounded density output.
- With flat density prior, fused output matches strand-only likelihood behavior.
- With weak `SS`, strand information is lower because likelihood curvature is lower,
  not because a scalar multiplier was applied.
- High-depth moderate-SS cases can still produce nonzero strand information.
- Fused mean and upper are always in `[0, N_r]`.
- Bounded large-count approximation matches exact fusion on sampled fixtures.

**Acceptance:** every region gets bounded fused evidence and diagnostics, with no lost
observed compatible mass.

### Phase 5: Footprint-aware PriorTable and MultiLocus allocation

**Goal:** Aggregate fused regional evidence to MultiLocus EM inputs without smearing
localized observed mass.

**Implement:**

- New module: `calibration/prior.py`.
- `PriorTable` with raw and EM-facing prior counts, expected counts, RNA lower counts,
  raw and exposure-weighted effective lengths, exposure weights, `enable_gdna`, and
  source/weight diagnostics.
- Observed count allocation from a footprint-aware substrate:
  - preferred: fragment footprint interval joins against MultiLoci;
  - acceptable first fallback: compact within-region bins;
  - base-pair overlap only for opportunity denominators or fully covered/count-uniform
    regions.
- Conservation diagnostics by region and MultiLocus.

**Tests:**

- Partial overlaps conserve count mass.
- Overlapping MultiLoci do not double-count the same observed footprint.
- A long fine region with a small read hotspot does not allocate hotspot prior mass to
  a MultiLocus that misses the hotspot.
- Uniform exposure path yields `gdna_eff_len == gdna_eff_len_raw`.
- Zero prior does not disable gDNA when geometry and likelihood support exist.

**Acceptance:** prior table exists independently of EM wiring and explains exactly how
count mass and denominators were allocated.

### Phase 6: Pipeline wiring to EM

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
- Local prior count bounds hold after MultiLocus aggregation.

**Acceptance:** `rigel quant` runs end-to-end before unsupervised exposure lands.

### Phase 7: Density refinement

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
- Tail/tension diagnostics decrease or remain explainable.

**Acceptance:** density quality improves without changing downstream interfaces.

### Phase 8: Unsupervised exposure learning

**Goal:** Learn hybrid-capture exposure `A_r` without confusing it with baseline bias.

**Implement:**

- `RegionExposure` learned mode.
- Bias-aware bootstrapping strategy from Section 8.
- Strong shrinkage and smoothing for exposure.
- Count closure diagnostics:

```text
sum predicted bounded D_mean approximately sum observed/fused D_mean over training rows
```

**Tests:**

- No-capture synthetic data learns near-uniform exposure.
- Capture synthetic data learns higher broad exposure over target-proximal regions.
- GC or mappability hotspots are not mistaken for capture exposure.
- Strand-specific capture still deconvolves identifiable exons locally.
- Ambiguous capture regions use density plus exposure.

**Acceptance:** captured loci get accessible-scale denominators without requiring BED
input.

### Phase 9: Validation and deferred expression-aware ambiguity

Validate across:

```text
1. unstranded, no capture
2. unstranded, capture
3. strand-specific, no capture
4. strand-specific, capture
```

Report:

- strand likelihood information and diagnostic strand contrast;
- density anchor families, fallback depth, and quality;
- density tail/tension against observed ledgers;
- boundary vs contained density disagreement;
- exposure mode and learned exposure distribution;
- fused density/strand information weights;
- footprint allocation conservation;
- prior mass and EM gDNA/RNA assignment error;
- fallback rates.

Expression-aware use of `TS_AMBIG` strand information is explicitly deferred to
post-1.0. Phase 0 through Phase 8 should rely on density for `TS_AMBIG` rather than
introducing a transcript-side expression prior that can destabilize EM.

## 10. Answered Open Questions

### 10.1 How should strand contrast map to `strand_strength`?

It should not map to an active model weight. Replace active `strand_strength` with
report-only diagnostics. The strand likelihood itself carries the effect of weak `SS`
through `p_r1_sense`, `p_gdna_sense`, and local count depth.

### 10.2 Should exact fusion use a density prior directly or a truncated one?

Use the unbounded density model in the fused posterior, then condition the posterior on
`0 <= D_r <= N_r` for EM-facing counts. Record density tail probability and expected
tail count rather than silently clipping the density model.

### 10.3 When is exact strand likelihood worth computing for weak `SS`?

Always for small counts, starting with `N_r <= 50`. For larger counts, use the bounded
large-count approximation. Weak `SS` does not automatically mean zero information;
high-depth regions can still be informative.

### 10.4 Can `TS_AMBIG` regions use strand information with a transcript-side prior?

Yes in principle, but defer it to post-1.0. Near-term calibration should not ask EM to
solve expression ambiguity and gDNA deconvolution at the same time in `TS_AMBIG`
regions.

## 11. Immediate Next Steps

Before Phase 1 implementation:

1. Write the concrete likelihood/fusion math spec for exact and bounded approximate
   fusion against the current `strand_deconv.py` semantics.
2. Decide whether Phase 5 will use fragment footprints or compact bins for observed
   count allocation.
3. Add memory-footprint targets for million-region ledgers.
4. Add contract tests for diagnostic-only strand information, bounded fused counts,
   density tail accounting, and hotspot-safe prior allocation.

Only after those locks should implementation start with `calibration/evidence.py`.
