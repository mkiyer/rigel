# Calibration Redesign: Continuous Exposure Factors

This directory is the implementation plan for rebuilding calibration around a continuous
per-region exposure factor. It builds on the theory note in
[new_calib_plan_v1.md](new_calib_plan_v1.md). The intent is to tear down the brittle
capture-classification path, keep the useful gDNA mass estimators, and rebuild the production path
with a smaller set of clear contracts.

No backward compatibility is required. Legacy behavior should be deleted from production code once
the replacement is in place. If old code is useful for reference, archive it outside the active
module path or in a clearly marked design note.

## Core Decision

Calibration should keep a latent expression state and remove the latent capture state.

The old v6 model has four states: unexpressed/off-target, unexpressed/capture,
expressed/capture, and expressed/off-target. That design makes capture a classification problem.
The new model treats capture, mappability, accessibility, and other sampling distortions as a
continuous regional exposure factor:

```text
omega_r = regional gDNA sampling density / global gDNA sampling density
```

`omega_r` is shrunk toward 1.0, but it is not a probability and it is not capped at 1.0. In capture
libraries it may be far below 1.0 in depleted regions and hundreds or thousands in targeted
regions. In ordinary RNA-seq it should remain tightly concentrated around 1.0 because the data will
learn a small exposure variance.

## Current Code Audit

The current production path already has useful pieces, but the contracts are mixed together:

- `src/rigel/calibration/latent_states.py` defines the four-state tensor and capture log Bayes
  factors.
- `src/rigel/calibration/calibration_iteration.py` iterates the four-state model, computes
  `RegionCalibration.A_r`, and carries `p_captured`, `gamma_r`, and `capture_enrichment_target`.
- `src/rigel/calibration/strand_deconv.py` provides the strongest existing gDNA/RNA separation for
  strand-informative libraries.
- `src/rigel/calibration/boundary_model.py` and `src/rigel/calibration/boundary_sweep.py` provide
  the boundary-to-contained gDNA projection used when strand information is absent or weak.
- `src/rigel/native/calibration/accumulator.cpp` accumulates fractional per-region mass, but it
  does not currently expose the number of physical fragment observations supporting each region.
- `src/rigel/calibration/prior.py` now assembles adaptive prior mass and component exposure inputs
  separately. Exposure is projected onto transcript and gDNA EM effective lengths after the prior
  mass split.

## New Production Contracts

### 1. Regional gDNA Mass

For each fine region `r`, calibration must produce:

```text
M_r = estimated unspliced gDNA mass in the region
T_r = total unspliced compatible fractional mass in the region
R_r = T_r - M_r
```

`M_r` is the source of truth for both adaptive prior construction and exposure learning. It should
come from strand deconvolution when strand contrast is identifiable, otherwise from the boundary and
contained density estimator. It must conserve unspliced mass after clipping:

```text
0 <= M_r <= T_r
0 <= R_r <= T_r
M_r + R_r = T_r
```

### 2. Regional Opportunity and Support

Exposure learning needs two different denominators:

- `O_r`: the mass opportunity used to turn `M_r` into a regional density. This must match the
  fractional mass contract. For fractional overlap mass, the first implementation should treat this
  as the exact region mass opportunity, not as the number of possible fragment starts.
- `N_r`: the physical fragment support count used to estimate uncertainty. A tiny region can have
  small fractional mass but many physical observations. The native accumulator must expose this.

The overlap start-site count for a region of width `W = b - a` and fragment length `L` is:

```text
start positions from a - L + 1 through b - 1
opportunity = W + L - 1
```

That start-site opportunity explains why tiny regions can have many physical observations. It is not
the same thing as the fractional mass opportunity. The implementation should keep these concepts
separate.

### 3. Exposure Factor

For each region:

```text
lambda_r = rho0 * O_r
raw_ratio_r = M_r / lambda_r
omega_r = EB_shrink(raw_ratio_r, support=N_r, center=1.0)
```

`rho0` is the global gDNA mass density learned from high-confidence unexpressed regions. The EB
model learns how much regional exposure variance exists in the library. Non-capture libraries should
learn low variance and shrink hard toward 1.0. Capture libraries should learn high variance and let
supported regions move far away from 1.0.

The planned production shrinkage family is a log-normal random-effects model: estimate regional
`log(raw_ratio_r)`, learn the library-wide exposure variance, and shrink the regional log ratio
toward `log(1.0) = 0` according to support-derived observation variance. A Gamma-Poisson posterior
mean can remain a diagnostic during validation, but PR 04 should ship one production exposure model,
not parallel capture-specific modes.

### 4. Downstream EM Usage

Exposure is sequence visibility, not a gDNA-only property. Hybrid-capture baits and other regional
sampling distortions do not know whether the molecule is genomic or transcriptomic, so downstream EM
must apply the same regional exposure surface to every overlapping component.

PR 05 represents exposure as component-level opportunity scaling:

```text
component_exposure_factor = bp_weighted_average(omega_r over component blocks)
eff_len_em = max(unweighted_eff_len * component_exposure_factor, 1.0)
```

This applies to annotated transcripts, synthetic nRNA transcripts, and locus gDNA components. The
adaptive prior mass split remains exposure-free and continues to use the calibrated gDNA/RNA mass
contract. If a future PR moves exposure into fragment-local likelihoods, it must replace this
component-denominator contract with a tested equivalent formulation rather than stacking both models.

## PR Sequence

1. [PR 01 - Two-State Calibration Teardown](pr01_v2.md)
2. [PR 02 - Native Observation Support and Boundary Payload](pr02_native_support_and_boundaries.md)
3. [PR 03 - Regional gDNA Mass Contract](pr03_region_gdna_mass.md)
4. [PR 04 - EB Exposure Factor Model](pr04_impl_plan_v2.md)
5. [PR 05 - Fair Component Exposure for EM](pr05_impl_plan_v3.md)
6. [PR 06 - Validation and Benchmarks](pr06_validation_and_benchmarks.md)

The sequence is intentionally staged so every PR has a local acceptance test. The native support
payload can land before the calibration math changes. The downstream sign fix can be tested with a
mock exposure surface before the full EB model is complete, but it should not be released as the
final behavior until `omega_r` is produced by the new calibration path.

## What We Keep

- Fractional region mass accumulation.
- Strand deconvolution and its reliability estimates.
- Boundary-local posterior and left/right boundary sweep, after detaching them from capture states.
- Mass-conserving `PriorMassDeconvolution`, renamed or reshaped around the new `RegionGdnaMass`
  contract.
- Grouped adaptive priors in native EM, after replacing entropy with a directional
  `p_unexpressed` soft gate.

## What We Delete or Replace

- Four-state latent constants and all capture-state posterior fields.
- `build_logbf_capture` and capture-specific density log Bayes factors.
- `p_captured`, `gamma_r`, and `capture_enrichment_target` from production outputs.
- gDNA-only exposure weighting or `gdna_eff_len / omega` denominator adjustments.
- Python-side transitional boundary construction once the native boundary payload is available.

## Acceptance Shape

A successful rebuild should show:

- Non-capture libraries: `omega_r` concentrated near 1.0, no widespread false gDNA rescue.
- Capture libraries: on-target regions can reach large exposure factors without clipping; off-target
  regions can fall below 1.0.
- Low-support regions: raw ratios are visibly shrunk toward 1.0.
- High-support regions: raw ratios largely survive shrinkage.
- Locus EM: regional exposure is applied fairly to every overlapping RNA and gDNA component.
- No reliance on gene-level constructs inside calibration, priors, locus construction, scoring, or EM.
