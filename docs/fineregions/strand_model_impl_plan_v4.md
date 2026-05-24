# Strand-Based gDNA / RNA Deconvolution - Implementation Plan v4

Date: 2026-05-23
Status: implementation-ready after top-down critique
Supersedes: `strand_model_impl_plan_v3.md`
Design basis: `rnaseq_mode_aware_gdna_density_plan.md`, `gdna_exposure_plan_v4.3.md`
Primary use case: strand-specific hybrid-capture RNA-seq

## 0. Top-Down Critique Of v3

v3 was much closer to the right implementation shape, but it still had four
serious problems.

1. **It deleted too much of the capture/exposure story.** The current
   `_regional_exposure.py` file is a stub and should go, but the exposure
   concept is not stale. Hybrid capture needs either known targets or an
   inferred per-region exposure field. The useful historical idea is a
   denominator-side exposure weight `A_r`, not the current cutover scaffolding.
2. **It replaced old complexity with new magic numbers.** `STRAND_ALPHA`,
   `EXON_SCREEN_N_MIN`, `EXON_SCREEN_R_FRAC`, and
   `EXON_SCREEN_P_REJECT` were statistical policy choices disguised as
   constants. The user must control the RNA lower-bound confidence. Exon
   self-training should use that same posterior decision, not its own private
   thresholds.
3. **It made the minimal prior too small.** A per-locus prior must still know
   the region's observed strand count and any capture/exposure weight. v3's
   `RegionGdnaEstimate` did not carry `n_total` or exposure, then the prior
   section referenced `region_gdna.n_total` anyway.
4. **It did not distinguish deletion from archival.** We can delete current
   modules from production code, but the implementation should record where the
   exposure work came from so we can draw from it later without keeping dead
   imports alive.

v4 fixes these issues while keeping the core implementation small.

## 1. System Shape

```text
BAM
  -> C++ scanner
  -> FragmentBuffer + StrandModel + fractional region accumulator
  -> calibration.calibrate(...)
       1. FL models
       2. CaptureProfile: uniform now, BED soon, inferred later
       3. strand kappa_d seed from high-purity regions
       4. exon self-training screen using the same RNA confidence rule
       5. strand kappa_d refit
       6. RegionGdnaEstimate with exposure_weight
  -> prior assembly
  -> EM
```

The calibration package must produce one downstream contract:

```python
@dataclass(frozen=True, slots=True)
class RegionGdnaEstimate:
    n_total: np.ndarray          # observed unspliced strand-model count per region
    mean_count: np.ndarray       # posterior mean gDNA count
    upper_count: np.ndarray      # upper bound on gDNA count
    rna_lower_count: np.ndarray  # lower bound on RNA count
    exposure_weight: np.ndarray  # A_r in [0, 1], uniform == 1
    precision: np.ndarray        # prior strength / confidence proxy in [0, 1]
    flags: np.ndarray
    kappa_d: float
    kappa_d_n_seed_regions: int
    kappa_d_n_exon_self_training: int
    p_r1_sense: float
    rna_lower_confidence: float
```

This restores the exposure seam without reintroducing regional-exposure
spaghetti. If capture is unknown, `exposure_weight` is all ones. If a BED file
is provided later, `exposure_weight` is a region-level target-overlap/capture
weight. If an unsupervised model lands later, it writes the same array.

## 2. User-Facing Policy Knobs

There should be exactly one new statistical user knob in the first
implementation:

```python
@dataclass(frozen=True)
class CalibrationConfig:
    rna_lower_confidence: float = 0.95
    capture_targets_bed: str | None = None  # accepted by CLI later; None => uniform
```

CLI:

```bash
rigel quant ... --rna-lower-confidence 0.95
rigel quant ... --capture-targets targets.bed
```

Meaning of `rna_lower_confidence`:

- Higher values are more conservative about calling gDNA-derived fragments RNA.
- `rna_lower_count` is the lower confidence bound on RNA count.
- `upper_count = n_total - rna_lower_count` is the corresponding conservative
  upper bound on gDNA count.
- Exon self-training uses the same confidence: an exon is accepted for gDNA
  training only if `P(R > 0 | observed counts) <= 1 - rna_lower_confidence`.

No separate `STRAND_ALPHA`. No separate exon p-value. No separate RNA fraction
threshold. The user's false-positive tolerance has one place to live.

## 3. What Is Not A User Parameter

The first implementation should not expose these as configuration:

- exact-vs-normal numerical switch;
- kappa fallback prior family;
- low-count exon threshold;
- spliced-count threshold;
- posterior-grid work budget.

If a value is only about runtime or numerical hygiene, keep it private and name
it as such. If a value changes biological/statistical behavior, either derive it
from `rna_lower_confidence` or make it a real config field. Do not hide policy
inside all-caps constants.

## 4. Exon Self-Training Without Magic Thresholds

The exon self-training rule is now simple and explainable.

Candidate filters:

1. Region is strand-informative: `ts_class` is `TS_POS` or `TS_NEG`.
2. Region has exon annotation and is strand-pure: `is_exon_any(signature)` and
   not `has_both_strands(signature)`.
3. Both transcript-relative channels are observed:
   `k_sense > 0` and `k_antisense > 0`.
4. There is no spliced evidence in the region. Use
   `np.isclose(spliced_sense + spliced_antisense, 0.0)` to avoid a new
   tolerance parameter.

Posterior rule:

```text
Accept exon for self-training iff
P(R > 0 | k_sense, k_antisense, p_r1_sense, kappa_d_seed)
    <= 1 - rna_lower_confidence
```

Interpretation: the same confidence level that governs RNA lower-bound calls
also governs which exons are safe enough to treat as gDNA-only training data.
At 95% confidence, the screen accepts an exon only when the posterior
probability of any RNA contribution is at most 5%.

No `EXON_SCREEN_N_MIN`. No `EXON_SCREEN_R_FRAC`. No separate
`EXON_SCREEN_P_REJECT`.

## 5. Kappa Training

### 5.1 Seed

Use high-purity regions:

```text
intergenic OR intron-only
contained unspliced POS/NEG counts only
```

Call the existing `estimate_strand_balance(...)` machinery. That function
already computes the symmetric beta-binomial method-of-moments fit for DNA
strand imbalance.

### 5.2 Fallback

If the method-of-moments estimate is not identified, use a named prior rather
than a magic concentration:

```text
DNA strand probability pi_d ~ Beta(1, 1)
```

That implies beta-binomial concentration `kappa_d = 2.0`. The code should name
this as `DNA_STRAND_UNIFORM_PRIOR`, not as an arbitrary fallback number.

This fallback is intentionally conservative: it allows broad DNA strand
imbalance when the data do not identify `kappa_d`.

### 5.3 Refit With Exons

After the exon screen:

1. Pool seed regions and accepted exons.
2. Recompute the same method-of-moments sufficient statistics.
3. Use the refit `kappa_d` for all region deconvolution.
4. Record accepted-exon count in `RegionGdnaEstimate` and `summary.json`.

Do not iterate. One seed pass, one exon screen, one refit.

## 6. Capture / Exposure Support

### 6.1 Keep The Concept, Delete The Stub

Delete current production `RegionalGdnaExposure` scaffolding from
`_regional_exposure.py`; it is a fail-fast stub. Do not delete the idea.
Archive references in the plan and commit message:

- `docs/calibration/gdna_exposure_plan_v4.3.md` - denominator-only exposure
  design;
- git commits `7a9caf9` and `d8f7a38` - historical implementation with
  `RegionalGdnaExposure.build`, `log_weights_for_positions`,
  `weighted_length_on_ref`, and `_weighted_quantile`;
- memory note `regional-gdna-exposure-v3-complete-2026-05` - prior behavior
  and gotchas.

The current production code should not import those old names. Git history is
the code archive; this plan records the handles.

### 6.2 New Minimal Capture Contract

Add a small future-proof object when capture support is wired:

```python
@dataclass(frozen=True, slots=True)
class CaptureProfile:
    mode: Literal["uniform", "bed", "inferred"]
    region_weight: np.ndarray  # float32[R], clipped to [0, 1]
    source: str                # "none", BED path, or inferred model name
```

M2 can start with `CaptureProfile.uniform(region_arrays)` only. The important
thing is that `RegionGdnaEstimate.exposure_weight` exists from day one, so the
BED and inferred paths do not change the downstream schema.

### 6.3 BED Path

The simplest known-capture implementation is:

```text
region_weight[r] = fraction of region r overlapped by any capture BED interval
```

Later refine this to FL-aware footprint exposure, but do not start there. The
first BED mode should be transparent and testable.

### 6.4 Inferred Path

The unsupervised path should not be reimplemented from the old shrinkage code.
The future model is:

```text
region features -> predicted capture/exposure weight A_r
```

Possible features:

- strand-deconvolved gDNA density;
- region class and length;
- local neighborhood density;
- mappability / GC / probe-like sequence features if available;
- uncertainty from `precision` and count depth.

This replaces empirical-Bayes shrink-to-global. The old EB shrinkage can be
deleted now because the future object is a predictive exposure model, not a
closed-form global-mean shrinker.

## 7. Cleanup Plan

### 7.1 Delete Or Replace

Delete:

```text
src/rigel/calibration/_regional_exposure.py
src/rigel/calibration/_locus_n_obs.py
src/rigel/calibration/_region_index_py.py
src/rigel/calibration/density_loco.py
src/rigel/calibration/_kappa.py
src/rigel/calibration/locus_prior.py
src/rigel/calibration/errors.py
```

Replace:

```text
src/rigel/calibration/locus_prior.py -> src/rigel/calibration/prior.py
```

Trim, do not blindly delete:

```text
src/rigel/calibration/_exposure.py
```

Keep the geometry helpers that remain conceptually useful:

- `gdna_eff_len_for_loci(...)`;
- `contained_exposure_clipped(...)` if prior assembly still needs clipped
  contained opportunity;
- `fractional_boundary_side_exposure(...)` if boundary evidence remains;
- `_merged_blocks(...)` as private geometry.

Delete or defer the old regional-weight application functions tied to
`RegionalGdnaExposure` until the new `CaptureProfile` is introduced.

### 7.2 CLI / Config Cleanup

Remove the stale regional-exposure switches:

```text
--regional-exposure
--no-regional-exposure
--regional-exposure-reference-quantile
CalibrationConfig.regional_exposure_enabled
CalibrationConfig.regional_exposure_reference_quantile
```

Add:

```text
--rna-lower-confidence FLOAT
--capture-targets BED       # optional; can land after uniform mode
CalibrationConfig.rna_lower_confidence
CalibrationConfig.capture_targets_bed
```

Validation:

```text
grep -R FractionalCutoverPending src/ tests/  # no matches
```

## 8. Strand Deconvolution API

Create one module:

```text
src/rigel/calibration/strand_deconv.py
```

Public functions:

```python
def build_strand_region_counts(...) -> StrandRegionCounts: ...
def estimate_kappa_d(...) -> StrandBalanceEstimate: ...
def screen_no_rna_exons(...) -> np.ndarray: ...
def deconvolve_regions_by_strand(...) -> RegionGdnaEstimate: ...
```

The deconvolution function accepts:

```python
def deconvolve_regions_by_strand(
    counts: StrandRegionCounts,
    *,
    kappa_d: float,
    rna_lower_confidence: float,
    exposure_weight: np.ndarray | None = None,
) -> RegionGdnaEstimate:
    ...
```

If `exposure_weight is None`, fill ones.

## 9. Prior Assembly

The prior should use region estimates and exposure, without reviving the old
multi-module prior machinery.

For each locus footprint:

1. Find overlapping regions.
2. Ignore rows with `FLAG_INELIGIBLE` or `FLAG_NEAR_UNSTRANDED` for the
   strand-derived prior contribution.
3. Weight each row by `precision * exposure_weight`.
4. Sum weighted `mean_count` into `gdna_prior_count`.
5. Sum weighted `(n_total - mean_count)` into `rna_prior_count`.
6. Carry a raw denominator/effective-length field for EM if the existing EM
   interface requires it.

Do not implement empirical-Bayes shrinkage here. If priors are too strong,
we weaken them through `precision` and the posterior uncertainty, not through a
separate shrink-to-global module.

## 10. Tests

### 10.1 Strand deconvolution

- `test_observations_ts_pos_and_ts_neg_basic`.
- `test_protocol_symmetry`.
- `test_near_unstranded_is_conservative`.
- `test_rna_lower_confidence_controls_bound` - higher confidence never
  increases `rna_lower_count`.
- `test_exon_screen_accepts_balanced_no_spliced_exon`.
- `test_exon_screen_rejects_clear_rna_exon`.
- `test_exon_screen_rejects_spliced_exon`.
- `test_exon_screen_requires_both_strand_channels`.
- `test_kappa_refit_records_exon_training_count`.

### 10.2 Capture seam

- `test_uniform_capture_profile_is_all_ones`.
- `test_region_gdna_estimate_defaults_exposure_to_one`.
- later with BED: `test_capture_targets_overlap_fraction`.

### 10.3 Prior

- `test_prior_weights_by_precision_and_exposure`.
- `test_prior_ignores_unstranded_rows`.
- `test_prior_preserves_required_em_denominator_fields`.

## 11. Ordered Implementation

1. **Cleanup and archive references.** Delete current fail-fast modules, trim
   `_exposure.py`, remove stale CLI/config flags, and mention the exposure
   archive handles from §6.1 in the PR description.
2. **Config knob.** Add `CalibrationConfig.rna_lower_confidence` and CLI
   `--rna-lower-confidence`; remove hard-coded alpha from the plan and code.
3. **Strand module.** Add `strand_deconv.py`, including two-pass exon
   self-training driven by `P(R > 0) <= 1 - rna_lower_confidence`.
4. **Result schema.** Add required `region_gdna` and `strand_region_counts`;
   include `exposure_weight` and `rna_lower_confidence` in the summary.
5. **Minimal prior.** Add `prior.py` consuming `RegionGdnaEstimate`.
6. **Pipeline.** Remove `FractionalCutoverPending`; run EM end-to-end.
7. **Capture BED.** Add `--capture-targets` once the uniform strand path is
   green. This should touch only `CaptureProfile` construction, not the prior
   or EM interface.
8. **Unsupervised exposure research.** Build the predictive `A_r` model after
   BED mode and validation, drawing on `gdna_exposure_plan_v4.3.md` but not
   restoring the old shrinkage implementation.

## 12. Done When

- There are no `FractionalCutoverPending` references in production code.
- `rigel quant` runs on a strand-specific fixture.
- The user can set `--rna-lower-confidence` and see the expected monotone
  behavior in RNA lower bounds.
- Exon self-training is active and reported in `summary.json`.
- `RegionGdnaEstimate.exposure_weight` exists and is all ones in uniform mode.
- Old exposure code is deleted from imports, but archived by doc/commit handles
  for later capture work.
