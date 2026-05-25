# gDNA Density Model Implementation Plan v1

Date: 2026-05-24
Status: implementation plan
Companion: `calibration_roadmap_v3.md`
Supersedes/refines: `initial_density_model_plan_v2.md`,
`rnaseq_mode_aware_gdna_density_plan.md`

## 0. Executive Summary

The proposed Poisson-Gamma shrinkage model is the right foundation for Rigel's
independent gDNA density channel. It has the key property we need for both
hybrid-capture and non-capture RNA-seq: local boundary flux can update a
region-specific latent density without forcing a hard target/off-target
classification.

The implementation should keep one important refinement: do not collapse
contained and boundary counts into one global training count too early. In
hybrid capture, boundary counts are exactly where target enrichment appears.
If they are pooled into the global background prior, the off-target prior can be
pulled upward and the model loses the contrast it was meant to learn.

The v1 implementation should therefore use this separation:

```text
contained DNA-dominant anchors  -> fit off-target Gamma prior
boundary flux in each region    -> update local lambda_r posterior
contained opportunity           -> predict held-out contained gDNA/background
observed boundary mass          -> direct local unspliced-background evidence
lambda_r / lambda_ref           -> relative exposure denominator for EM
```

This keeps the model simple, conjugate, and continuous. It also keeps the
density channel independent from strand deconvolution while giving the later
fusion layer a real predictive distribution, not just a point estimate.

## 1. Evaluation of the Proposed Formulation

### 1.1 What is strong

The core latent-density model is well matched to the data Rigel now collects:

```text
lambda_r ~ Gamma(alpha_0, beta_0)
N_b,r | lambda_r ~ Poisson(lambda_r * L_b,r)
lambda_r | N_b,r ~ Gamma(alpha_0 + N_b,r, beta_0 + L_b,r)
N_c,r | N_b,r ~ NegativeBinomial(
    shape = alpha_0 + N_b,r,
    p = (beta_0 + L_b,r) / (beta_0 + L_b,r + L_c,r),
)
```

This naturally handles the two important regimes:

- Off-target or non-captured regions have sparse boundary evidence, so they
  shrink to the background prior.
- Captured regions have elevated boundary evidence, so the local posterior
  density increases without needing a hard on-target label.

The Negative Binomial predictive distribution is especially useful because the
upper count bound propagates both sampling noise and uncertainty in the local
exposure estimate.

### 1.2 Required refinements

First, the global background prior should be fit from contained intergenic and
intron-only anchors by default. Boundary counts should be retained separately
and used as local exposure evidence. A later robust mixture can learn explicit
off-target/on-target boundary components, but v1 should not let enriched
boundary counts redefine the off-target prior.

Second, the density model should call its target carefully. If boundary counts
are used to update `lambda_r`, they should not also be treated as an independent
future draw from the same predictive distribution. For a region-level density
surface, v1 should emit:

```text
predicted_contained_mean      = E[N_c,r | boundary]
predicted_contained_upper     = Q_conf[N_c,r | boundary]
observed_boundary_count       = N_b,r
mean_unbounded                = observed_boundary_count + predicted_contained_mean
upper_unbounded               = observed_boundary_count + predicted_contained_upper
```

The later fusion layer will bound this against observed compatible mass. The
density model itself should also report tail diagnostics such as
`P(predicted_contained > observed_contained)`.

Third, the nascent-RNA interpretation should be explicit. Boundary-driven
lambda inflation is conservative for mature RNA protection because it can
absorb unspliced background. But it is not always pure gDNA, especially in
unstranded transcribed introns. The density output should remain a gDNA-density
evidence channel with quality/tension flags; strand fusion should be allowed to
move mass away from gDNA when strand evidence indicates RNA or nRNA.

Fourth, exposure weights should be relative opportunities, not capped
probabilities. The useful EM denominator is:

```text
A_r = E[lambda_r | boundary] / lambda_ref
gdna_eff_len_weighted = sum_r A_r * raw_gdna_opportunity_r
```

`A_r` can be greater than 1 in captured regions. The current
`RegionExposure` docstring says `(0, 1]`; that contract should be revised when
the density-derived exposure surface becomes an EM input.

## 2. Current Code Surface

Already implemented and useful:

- `calibration._arrays.RegionArrays`: sorted `(ref_id, start)` geometry,
  `signature`, `ts_class`, and `ref_offsets`.
- `calibration._arrays.PayloadArrays`: sorted 12-channel count view with hot
  unspliced contained and boundary POS/NEG totals.
- `calibration.density_global.l_eff_contained`: FL-PMF-weighted contained
  opportunity.
- `calibration._exposure.fractional_boundary_side_exposure`: FL-PMF-weighted
  one-sided fractional boundary opportunity. The input is the receiving
  region's length, not the neighbor gap.
- `calibration.strand_deconv`: independent strand-only per-region gDNA/RNA
  posterior.
- `calibration._result.CalibrationResult`: current result object with
  `region_gdna` from strand deconvolution and uniform `region_exposure`.
- `pipeline._run_locus_em_partitioned`: already accepts raw/EM gDNA prior
  counts, weighted/unweighted gDNA effective lengths, exposure weight, and
  explicit `enable_gdna`.

Important current boundary:

- `quant_from_buffer` intentionally raises `NotImplementedError` until
  footprint-aware prior assembly and EM wiring land. Density model v1 should
  not unblock EM by itself.

## 3. Statistical Contract

### 3.1 Primitive observations

For each region `r`:

```text
C_r   = contained_unspliced_pos + contained_unspliced_neg
B_l   = boundary_left_unspliced_pos + boundary_left_unspliced_neg
B_r   = boundary_right_unspliced_pos + boundary_right_unspliced_neg
B_tot = B_l + B_r

L_c  = contained effective length under gDNA FL
L_l  = left boundary-side opportunity under gDNA FL
L_r  = right boundary-side opportunity under gDNA FL
L_b  = L_l + L_r
```

Use only unspliced channels for density. Spliced channels remain RNA evidence
and diagnostics.

### 3.2 Boundary opportunity

For each side of a region with length `S`, the denominator must match the
fractional side mass emitted by the accumulator:

```text
L_side(S) = sum_ell h(ell) * min((ell - 1) / 2, S / 2)
```

This is already implemented by `fractional_boundary_side_exposure(lengths_bp,
gdna_fl)`. Use the current region length for each side. Do not use adjacent
region gap lengths, and do not resurrect `splicing_anchor_tolerance` for
calibration.

### 3.3 Background prior fit

Fit Gamma priors to DNA-dominant contained anchors:

```text
lambda_class ~ Gamma(alpha_class, beta_class)
C_i | lambda_i ~ Poisson(lambda_i * L_c,i)
```

Initial anchor classes:

```text
INTERGENIC: signature == 0x0
INTRON:     intron bits present and no exon bits (0x4, 0x8, 0xC)
```

Eligibility:

```text
L_c,i > min_effective_length
C_i is finite and nonnegative
```

The first implementation should fit class-level genome-wide priors and include
the deterministic fallback hooks needed for later reference-specific pooling:

```text
class + reference
class genome-wide
all anchors on reference
all anchors genome-wide
weak broad fallback
```

The first version may only populate class genome-wide and all-anchors
genome-wide, but the dataclasses should carry `family` and `fallback_depth` so
the interface does not have to change.

### 3.4 Gamma prior parameterization

Use shape/rate form:

```text
E[lambda] = alpha / beta
Var(lambda) = alpha / beta^2
```

Fit by marginal Gamma-Poisson likelihood over anchor regions:

```text
log p(C_i | alpha, beta, L_i)
  = gammaln(C_i + alpha) - gammaln(alpha) - gammaln(C_i + 1)
    + alpha * log(beta / (beta + L_i))
    + C_i * log(L_i / (beta + L_i))
```

Counts are fractional in Rigel, but the same expression is valid as the Gamma
extension of the count likelihood and works with `gammaln(C_i + 1)`. Optimize
`log(alpha)` and `log(beta)` with bounded deterministic starting values.

Fallback when fitting is ill-conditioned:

```text
rho = sum C_i / sum L_i
alpha = max(min_prior_count, rho * fallback_prior_opportunity)
beta = max(min_prior_opportunity, fallback_prior_opportunity)
```

The fallback is intentionally a broad prior, not a hard point estimate.

### 3.5 Local boundary update

For every region, choose the best available prior family and update it with
local boundary evidence:

```text
alpha_post,r = alpha_prior,r + B_tot,r
beta_post,r  = beta_prior,r  + L_b,r
rho_post,r   = alpha_post,r / beta_post,r
```

If `L_b,r == 0` and `B_tot,r == 0`, the posterior equals the background prior.
If `B_tot,r` is high relative to `L_b,r`, the posterior learns capture enrichment
locally.

### 3.6 Contained predictive distribution

The held-out contained count has a Negative Binomial predictive distribution:

```text
N_c,r^gdna | boundary ~ NB(
    shape = alpha_post,r,
    p = beta_post,r / (beta_post,r + L_c,r),
)
```

Mean and variance:

```text
mean_c = (alpha_post / beta_post) * L_c
var_c  = mean_c * (1 + L_c / beta_post)
```

Upper confidence:

```text
upper_c = scipy.stats.nbinom.ppf(confidence, shape, p)
```

`scipy.stats.nbinom` supports non-integer shape. If performance becomes a
problem, add a normal or saddlepoint approximation for large means after the
exact path is validated.

### 3.7 Region-level density evidence

The density model should emit unbounded predictions and local tension
diagnostics:

```text
mean_unbounded  = B_tot + mean_c
upper_unbounded = B_tot + upper_c
var_unbounded   = var_c

observed_compatible_count = B_tot + C_r
tail_probability = P(N_c^gdna > C_r | boundary)
expected_tail_count = E[(N_c^gdna - C_r)+ | boundary]
```

Do not clip `mean_unbounded` or `upper_unbounded` in the density model. Phase 4
fusion will condition on the physically possible range `0 <= D_r <= N_r` for
EM-facing counts. Tail metrics should be carried forward as model tension.

### 3.8 Exposure weight

Let `lambda_ref` be the off-target reference density for the selected prior
family, usually `alpha_prior / beta_prior` or a robust all-anchor background.

```text
A_r = rho_post,r / max(lambda_ref, eps)
```

This `A_r` is the density-derived relative exposure candidate. Store it in
`DensityEvidence`. Do not force it into `RegionExposure` or EM until Phase 5
prior allocation defines the denominator aggregation path.

## 4. New Modules and Dataclasses

### 4.1 `calibration/evidence.py`

Roadmap v3 Phase 1 should land before density modeling so count channels have a
single memory-light contract.

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

Derived totals should be methods with optional `out` scratch buffers:

```python
contained_total(out=None)
boundary_total(out=None)
unspliced_total(out=None)
observed_total(out=None)
```

Keep raw count arrays as `float32`; use `float64` for reductions and model
fitting.

### 4.2 `calibration/density_observation.py`

```python
@dataclass(frozen=True, slots=True)
class DensityObservation:
    contained_count: np.ndarray          # float32[R]
    boundary_left_count: np.ndarray      # float32[R]
    boundary_right_count: np.ndarray     # float32[R]
    boundary_count: np.ndarray           # float32[R]
    observed_compatible_count: np.ndarray # float32[R]

    contained_leff: np.ndarray           # float64[R]
    boundary_left_leff: np.ndarray       # float64[R]
    boundary_right_leff: np.ndarray      # float64[R]
    boundary_leff: np.ndarray            # float64[R]
    total_leff: np.ndarray               # float64[R]

    class_code: np.ndarray               # uint8[R]
    anchor_intergenic: np.ndarray        # bool[R]
    anchor_intron: np.ndarray            # bool[R]
```

Builder:

```python
build_density_observation(region_arrays, ledger_or_payload_arrays, gdna_fl)
```

Use `l_eff_contained(region_arrays.end - region_arrays.start, gdna_fl)` and
`fractional_boundary_side_exposure(region_lengths, gdna_fl)` for each side.

### 4.3 `calibration/density_model.py`

```python
@dataclass(frozen=True, slots=True)
class GammaRatePrior:
    family: np.ndarray | str
    alpha: float
    beta: float
    mean_density: float
    n_regions: int
    n_fragments: float
    eff_length: float
    fallback_depth: int
    fit_status: str
```

```python
@dataclass(frozen=True, slots=True)
class DensityEvidence:
    mean_unbounded: np.ndarray
    upper_unbounded: np.ndarray
    variance_unbounded: np.ndarray
    precision_unbounded: np.ndarray
    leff_contained: np.ndarray
    leff_boundary: np.ndarray
    anchor_support: np.ndarray
    relative_exposure: np.ndarray
    tail_probability: np.ndarray
    expected_tail_count: np.ndarray
    family: np.ndarray
    fallback_depth: np.ndarray
    quality: np.ndarray
    confidence: float
```

Public builder:

```python
fit_density_evidence(
    observation,
    region_arrays,
    *,
    confidence: float,
    min_effective_length: float,
) -> DensityEvidence
```

## 5. Calibration Integration

### 5.1 Orchestrator sequence

The v1 calibration flow should become:

```text
build FL models
build RegionArrays and PayloadArrays
build RegionCountLedger / evidence bundle
compute current contained global densities for summary compatibility
build strand counts and strand-only RegionGdnaEstimate
build density observations
fit density evidence
return CalibrationResult with strand and density surfaces
```

The density model should not consume `region_gdna` from strand deconvolution.
The two evidence channels meet later in `calibration/integration.py`.

### 5.2 CalibrationResult changes

Add optional fields without renaming the current strand output:

```python
density_observation: DensityObservation | None = None
density_evidence: DensityEvidence | None = None
```

`to_summary_dict()` should include compact density summaries:

```text
confidence
n_regions
n_anchor_intergenic
n_anchor_intron
background_prior families and fallback counts
mean/upper/tail summaries
relative_exposure min/mean/p95/max
tail_probability p50/p95/max
```

Avoid serializing full per-region arrays to `summary.json`.

### 5.3 Config changes

Add to `CalibrationConfig`:

```python
gdna_density_confidence: float = 0.95
density_min_effective_length: float = 1.0
density_fit_min_regions: int = 20
density_fit_min_fragments: float = 20.0
```

Expose `--gdna-density-confidence` in the CLI because the user-facing upper
gDNA confidence level is part of the model contract. The other knobs can remain
internal until validation proves they need CLI exposure.

## 6. Relationship to Strand Fusion

Density evidence is independent of strand deconvolution. It should not use
`strand_strength`, `p_r1_sense`, or `kappa_d`.

Phase 4 fusion will consume:

```text
DensityEvidence.mean_unbounded / variance_unbounded / upper_unbounded
RegionGdnaEstimate or StrandEvidence likelihood summaries
observed compatible count from DensityObservation
```

Rules for fusion:

- Neutral strand likelihood for `TS_NONE`, `TS_AMBIG`, and near-unstranded
  regions.
- Exact bounded fusion for small `N_r`.
- Bounded Laplace/truncated-normal approximation for large `N_r`.
- Density tail is reported, not silently clipped.

This separation is important: density can solve unstranded and capture exposure
problems; strand can deconvolve gDNA/RNA when transcript-relative orientation is
identifiable.

## 7. Relationship to Prior Allocation and EM

Density model v1 should not assemble locus priors directly. Phase 5 will build
`calibration/prior.py` and decide how region-level evidence intersects
MultiLoci without smearing localized read hotspots.

When Phase 5 lands, the denominator should use the density-derived relative
exposure:

```text
gdna_eff_len_raw      = sum raw gDNA opportunity over locus footprint
gdna_eff_len_weighted = sum A_r * raw gDNA opportunity over locus footprint
gdna_prior_count_raw  = footprint-allocated fused gDNA count expectation
gdna_prior_count_em   = uncertainty-adjusted EM-facing prior count
```

Observed prior count mass must come from footprint-aware or binned allocation,
not naive base-pair overlap. Geometry-only prorating is acceptable for
denominators, not for observed count mass, unless the region is fully covered or
known count-uniform.

## 8. Tests

### 8.1 Evidence ledger tests

- Raw ledger arrays are `float32`; reductions are `float64`.
- Contained, boundary, unspliced, spliced, and observed totals are derived
  correctly.
- `TS_NONE` and `TS_AMBIG` rows retain observed counts even when strand-folded
  totals are zero.
- Million-region memory footprint stays within the roadmap target.

### 8.2 Density observation tests

- Boundary-free fixtures reproduce contained-only density observations.
- Boundary counts and side opportunities are both included.
- Boundary side opportunity uses region length and matches hand-computed
  `sum h(ell) * min((ell - 1) / 2, S / 2)`.
- Intergenic and intron-only anchor masks are deterministic.
- Exon and mixed signatures receive observations but are not background anchors.

### 8.3 Gamma fit tests

- Simulated Gamma-Poisson anchors recover mean density within tolerance.
- Poisson-like anchors fall back to a high-concentration or broad stable prior
  rather than failing.
- Sparse anchors trigger deterministic fallback depth.
- Fractional counts are accepted and produce finite parameters.

### 8.4 Local posterior tests

- With zero boundary evidence, posterior density equals the background prior.
- With high boundary density, posterior density and relative exposure increase.
- Non-capture-like sparse boundary data does not inflate off-target regions.
- `mean_unbounded` and `upper_unbounded` are not clipped to observed counts.
- Tail probability is positive when predicted contained gDNA exceeds observed
  contained mass.

### 8.5 Calibration integration tests

- `calibrate()` returns density observation and density evidence blocks.
- Summary includes prior family, fallback, tail, and relative exposure
  diagnostics.
- Strand-specific weak libraries alter strand diagnostics only; density evidence
  is unchanged for the same region payload.
- Existing strand deconvolution tests continue to pass without consuming density.

## 9. Implementation Order

1. Add `calibration/evidence.py` and tests for the memory-light ledger.
2. Add `calibration/density_observation.py` and tests for count/opportunity
   construction.
3. Add `calibration/density_model.py` with Gamma prior fitting, local boundary
   updates, NB predictive moments/quantiles, exposure candidates, and tail
   diagnostics.
4. Add `CalibrationConfig.gdna_density_confidence` and internal defaults.
5. Extend `CalibrationResult` summary with density blocks.
6. Wire density observation/evidence into `_orchestrator.calibrate()`.
7. Validate on synthetic no-capture and hybrid-capture calibration fixtures.
8. Only after this is stable, implement `integration.py` for probabilistic
   density-strand fusion.

## 10. Known Limitations of v1

- Boundary evidence may be weak for long captured exons whose probes sit far
  from exon-intron or exon-intergenic edges. Later exposure learning should
  smooth/borrow across neighboring regions or consume optional capture BED/probe
  priors.
- Intronic contained anchors can contain nascent RNA. Strand fusion and future
  robust anchor filtering should reduce this contamination; v1 should report
  anchor quality rather than pretending it is absent.
- Mappability-corrected denominators are still deferred. The density interface
  should leave room for `contained_leff_mappable` and `boundary_leff_mappable`
  without changing downstream consumers.
- The first exposure surface is local and evidence-driven, not a full
  unsupervised capture-accessibility model. Phase 8 can replace or smooth it
  while preserving the `relative_exposure` contract.
