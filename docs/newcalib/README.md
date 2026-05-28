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
- `src/rigel/calibration/prior.py` currently multiplies gDNA effective length by the regional
  exposure average. Because native EM subtracts `log(gdna_eff_len)`, this has the wrong sign for
  enriched regions.

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

### 4. Downstream EM Usage

Exposure increases sampling visibility. In the initial production model, exposure is represented as
a denominator adjustment only:

```text
gdna_eff_len_em = gdna_eff_len_unweighted / omega_locus
```

Do not multiply the effective length by exposure. Native EM applies `-log(gdna_eff_len)` to the gDNA
component, so multiplying by `omega_locus` penalizes exactly the high-exposure regions where gDNA
should be most competitive.

Do not stack a denominator adjustment and a per-fragment `+log(omega)` adjustment in the same PR. If
we later want a fully normalized local-exposure likelihood, that must replace the denominator-only
model with a tested equivalent contract.

## PR Sequence

1. [PR 01 - Two-State Calibration Teardown](pr01_two_state_calibration.md)
2. [PR 02 - Native Observation Support and Boundary Payload](pr02_native_support_and_boundaries.md)
3. [PR 03 - Regional gDNA Mass Contract](pr03_region_gdna_mass.md)
4. [PR 04 - EB Exposure Factor Model](pr04_eb_exposure_model.md)
5. [PR 05 - Downstream EM Exposure Normalization](pr05_downstream_em_exposure.md)
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
- Grouped adaptive priors in native EM, after updating their entropy calculation for two states.

## What We Delete or Replace

- Four-state latent constants and all capture-state posterior fields.
- `build_logbf_capture` and capture-specific density log Bayes factors.
- `p_captured`, `gamma_r`, and `capture_enrichment_target` from production outputs.
- Effective-length multiplication by exposure.
- Python-side transitional boundary construction once the native boundary payload is available.

## Acceptance Shape

A successful rebuild should show:

- Non-capture libraries: `omega_r` concentrated near 1.0, no widespread false gDNA rescue.
- Capture libraries: on-target regions can reach large exposure factors without clipping; off-target
  regions can fall below 1.0.
- Low-support regions: raw ratios are visibly shrunk toward 1.0.
- High-support regions: raw ratios largely survive shrinkage.
- Locus EM: enriched gDNA regions are no longer suppressed by inflated effective length.
- No reliance on gene-level constructs inside calibration, priors, locus construction, scoring, or EM.
