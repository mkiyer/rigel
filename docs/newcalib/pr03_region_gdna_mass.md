# PR 03 - Regional gDNA Mass Contract

## Goal

Make regional gDNA mass a first-class calibration output independent of capture states.

This PR defines the handoff that exposure learning and adaptive priors both consume:

```text
T_r = total unspliced compatible fractional mass
M_r = estimated gDNA portion of T_r
R_r = estimated RNA portion of T_r
```

The implementation should preserve the good estimators already in the codebase and remove the
capture-state mixing around them.

## Current State

`PriorMassDeconvolution` in `src/rigel/calibration/calibration_iteration.py` already enforces:

```text
gdna_unspliced_mean + rna_unspliced_mean = unspliced_total
```

The inputs are useful:

- `DensityObservation.observed_compatible_count` contains contained plus boundary unspliced mass.
- `RegionGdnaChannelEstimate` provides compartment-wise strand deconvolution when strand contrast is
  identifiable.
- `BoundaryLocalPosterior` and `BoundarySweepResult` project boundary excess into contained gDNA
  predictions.

The problem is that `mu_gdna` is currently mixed using `p_captured`, which will no longer exist.

## New Dataclass

Introduce a dataclass with an explicit name, either replacing or renaming `PriorMassDeconvolution`:

```text
RegionGdnaMass
    unspliced_total: float32[R]
    gdna_mass: float32[R]
    rna_mass: float32[R]
    mass_opportunity: float64[R]
    support_count: uint32/uint64[R]
    support_ess: float64[R]
    method: uint8[R]
    precision: float32[R]
    flags: uint16[R]
```

`support_count` comes from PR 02. `support_ess` should prefer the squared-mass ESS when available;
otherwise use `support_count` as the fallback.

## Mass Opportunity

The mass opportunity must match the fractional mass being normalized. Do not use possible start-site
count as the density denominator for fractional mass.

First implementation:

```text
mass_opportunity_r = region width in bp, with chromosome-edge handling if needed
```

If tests show edge effects matter, replace this with an exact FL-marginal fractional opportunity.
The key is that the expected fractional mass of a tiny region is small even though many fragments can
touch it.

Boundary and contained effective lengths remain useful for the boundary projection model, but they
should not be confused with the mass denominator for exposure learning.

## Estimator Logic

### Strand-Informative Libraries

When strand contrast is identifiable:

```text
M_r = contained_mean * contained_reliability
    + boundary_left_mean * boundary_left_reliability
    + boundary_right_mean * boundary_right_reliability
```

Then clip to `[0, T_r]` and set `R_r = T_r - M_r`.

The precision field should reflect the best reliable strand channel, as the current code does, but
it should be interpreted as mass-estimator precision, not capture precision.

### Weak-Strand or Unstranded Libraries

When strand contrast is not identifiable:

1. Fit `rho0` from high-confidence unexpressed seed regions.
2. Build boundary-local posterior from boundary-compatible unspliced evidence.
3. Run boundary sweeps to impute contained gDNA mass.
4. Combine direct contained evidence and boundary-imputed evidence into `M_r`.
5. Clip to `[0, T_r]` and conserve mass.

The exact combine rule should be conservative in the first implementation. A good starting rule is:

```text
M_r = min(T_r, max(direct_contained_gdna_mean, boundary_sweep_mean))
```

If this proves too eager in no-gDNA libraries, use posterior precision to blend rather than take a
max. The important contract is that `M_r` is the final regional estimate, not an intermediate
capture-state quantity.

## Background Refit

`rho0` should be fit from unexpressed regions, weighted by `p_unexpressed` from PR 01:

```text
rho0 = sum_r p_unexpressed_r * M_r / sum_r p_unexpressed_r * O_r
```

Retain the existing safeguards:

- exclude low-opportunity regions from training when they would dominate noise,
- exclude regions with spliced evidence,
- exclude strand-RNA lower-bound regions,
- trim the highest-density seed tail during the initial background fit.

For capture data, trimming and expression exclusion should prevent on-target outliers from defining
the global center. For ordinary data, the center should be stable and close to the observed global
gDNA density.

## Adaptive Prior Handoff

`src/rigel/calibration/prior.py` and `src/rigel/calibration/adaptive_prior.py` should consume:

```text
gdna_unspliced_mean = region_gdna_mass.gdna_mass
rna_unspliced_mean = region_gdna_mass.rna_mass
unspliced_total = region_gdna_mass.unspliced_total
```

No capture posterior should enter adaptive prior construction.

## Tests

- Mass conservation after float32 conversion.
- Strand-informative synthetic region where antisense-compatible mass is assigned to gDNA.
- Unstranded boundary-only synthetic region where boundary evidence imputes contained gDNA.
- No-gDNA/no-RNA sentinel where `M_r` stays zero or near zero.
- Low-opportunity region where support exists but mass opportunity prevents density blow-up from
  training `rho0`.

## Done Means

- There is one production regional gDNA mass table.
- Exposure learning and adaptive priors use the same `M_r` source.
- Capture state fields are not needed to compute regional gDNA mass.
