# Current gDNA Prior Construction

Date: 2026-05-26
Status: current-behavior audit after Phase IV clean cutover

This note describes how calibration output is currently translated into locus-level EM inputs for the gDNA component. It is descriptive, not a proposal. The goal is to make the existing numerator, denominator, pseudocounts, floors, warm start, and eligibility behavior explicit before designing a principled RNA/gDNA prior.

## Executive Summary

The current EM gDNA input is not just an arbitrary pseudocount of `1.0`.

For each `MultiLocus`, Python constructs three live inputs:

```text
gdna_prior_count_em[locus] = allocated expected gDNA count from RegionCalibration.mu_gdna
gdna_eff_len[locus]        = FL-marginal gDNA locus effective length * mean regional A_r
enable_gdna[locus]         = technical likelihood eligibility bit
```

The allocated count is the numerator-like calibration prediction. The effective length is the denominator-like exposure term used by native EM through `-log(gdna_eff_len)`. The eligibility bit only asks whether the locus has at least one unspliced unit with finite gDNA likelihood.

There are, however, several independent floors/baselines:

- The background density model uses a Gamma prior with `alpha_floor = 1.0` and `beta_floor = 1.0`.
- The gDNA locus effective length is floored at `1.0` after exposure weighting.
- The exposure weight has a minimum of `1.0e-4`.
- Native VBEM adds a `0.5` Jeffreys-style baseline to every eligible component, including RNA components and the gDNA component.
- Native gDNA prior count is guarded by `max(gdna_prior_count, 1.0e-300)`, which is numerical protection, not a meaningful count.
- Native `EM_PRIOR_EPSILON = 1.0e-10` is used as an eligibility/bootstrap sentinel before the final prior vector is built; it is not the modeled per-locus gDNA prior.

The most important current asymmetry is that calibration contributes an explicit extra prior count only to gDNA. RNA components only receive the native VBEM baseline in VBEM mode, or no modeled prior count in MAP mode. There is not yet a calibrated RNA-vs-gDNA prior pair.

## End-To-End Data Flow

The live path is:

```text
calibration payload
-> DensityObservation
-> BackgroundModel
-> BoundaryTable / BoundarySweepResult
-> RegionCalibration
-> PriorTable
-> AbundanceEstimator.run_batch_locus_em_partitioned()
-> native batch_locus_em_partitioned()
-> compute_gdna_prior_and_warm_start()
-> SQUAREM MAP/VBEM locus EM
```

The old live density/fusion/exposure handoff is no longer used. `CalibrationResult` now carries `region_calibration`, `background_model`, `boundary_local`, `boundary_sweep`, and diagnostics. Downstream prior assembly consumes `region_calibration.mu_gdna`, `region_calibration.upper_gdna`, and `region_calibration.A_r` directly.

## Regional Calibration Quantities

The four-state calibration produces `RegionCalibration` arrays aligned to `index.region_df`:

```text
p_states[r, :]   posterior state probabilities
mu_gdna[r]       expected regional gDNA count
upper_gdna[r]    conservative upper gDNA count bound
rna_lower[r]     conservative RNA lower bound; currently diagnostic for prior assembly
A_r[r]           unitless regional gDNA exposure multiplier
gamma_r[r]       captured-only exposure diagnostic
rho_off          off-target gDNA density mean
```

For each region, the E-step computes:

```text
p_captured[r] = P(gdna_only_capture) + P(expressed_capture)

off_target_mu[r] = rho_off_mean * contained_leff[r]
captured_mu[r]   = boundary_sweep.mu_sweep[r]

mu_gdna[r] = p_captured[r] * captured_mu[r]
           + (1 - p_captured[r]) * off_target_mu[r]
```

`upper_gdna` is mixed the same way, using `boundary_sweep.upper_sweep` for captured regions and a background predictive upper bound for off-target regions.

The exposure multiplier is then:

```text
A_r[r] = mu_gdna[r] / max(rho_off_mean * contained_leff[r], eps)
```

with numerical clipping to a finite non-negative range.

Interpretation:

- `mu_gdna` is the regional expected-count prediction that later becomes the per-locus gDNA prior numerator after geometric allocation.
- `A_r` is a density multiplier relative to the off-target background. It later scales the locus-level gDNA effective length.
- `upper_gdna` is currently preserved for diagnostics but no longer gates gDNA eligibility.
- `rna_lower` is currently not used in `assemble_priors()` after the cleanup.

## Where The `1.0` Background Pseudocount Enters

The background model is a Gamma posterior for the off-target gDNA density `rho_off`:

```text
rho_off_alpha = alpha_floor + background_seed_gdna_count
rho_off_beta  = beta_floor  + background_seed_effective_length
rho_off_mean  = rho_off_alpha / rho_off_beta
```

The default floors are:

```text
alpha_floor = 1.0
beta_floor  = 1.0
```

The same default floors are used again when the calibration M-step refits `rho_off` from posterior-weighted background evidence:

```text
alpha_hat = 1.0 + weighted_background_gdna_count
beta_hat  = 1.0 + weighted_background_effective_length
```

Then the refit is damped against the previous background posterior.

This `1.0` is upstream of the per-locus EM prior. It affects `rho_off_mean`, which affects `off_target_mu`, `mu_gdna`, and `A_r`. It is not itself directly passed as `gdna_prior_count_em`.

## Projection From Regions To Loci

`assemble_priors()` builds one `PriorTable` row per `MultiLocus`.

### Count Numerator

The expected regional count array is:

```text
region_mean = region_calibration.mu_gdna
```

Each region's `mu_gdna` is allocated to overlapping `MultiLocus` genomic blocks by base-pair overlap. If multiple loci overlap a region, the region count is split by normalized overlap share:

```text
raw_share[locus] = bp_overlap(region, locus_blocks) / region_length
share[locus]     = raw_share[locus] / sum(raw_share over overlapping loci)

gdna_expected[locus] += share[locus] * mu_gdna[region]
gdna_upper[locus]    += share[locus] * upper_gdna[region]
```

Regions that touch no locus are accumulated into `unallocated_expected_count` and `unallocated_upper_count` diagnostics.

The EM prior count is currently exactly:

```text
gdna_prior_count_em[locus] = gdna_expected[locus]
```

There is no additional Python-side per-locus pseudocount added here. `gdna_upper_count` remains diagnostic under the clean baseline; it no longer gates or suppresses gDNA.

### Effective-Length Denominator

Python computes an unweighted gDNA effective length for each `MultiLocus`:

```text
gdna_eff_len_unweighted[locus]
    = FL-PMF-marginal number of valid gDNA fragment starts overlapping the locus blocks
```

For a single interior interval, this reduces to the closed-form span-based effective length used by `FragmentLengthModel.compute_all_transcript_eff_lens()`. For multiple intervals or edge cases, it sums valid fragment-start windows over each positive fragment length in the gDNA FL PMF and merges overlapping windows by reference.

Then Python computes an exposure weight:

```text
gdna_em_exposure_weight[locus]
    = bp-weighted mean of region_calibration.A_r over the locus blocks
```

The weight is floored at `1.0e-4` and is not clipped above `1.0`.

The final native denominator is:

```text
gdna_eff_len[locus]
    = max(gdna_eff_len_unweighted[locus] * gdna_em_exposure_weight[locus], 1.0)
```

This is the exposure-adjusted gDNA effective length passed to native EM. Native EM stores `log(gdna_eff_len)` and subtracts it from gDNA component log weights, matching the RNA component effective-length correction.

### Technical Eligibility

Python computes:

```text
enable_gdna[locus] = any unit in the locus is unspliced and has finite gDNA log-likelihood
```

This is intentionally independent of `mu_gdna`, `upper_gdna`, `A_r`, or `p_captured`. The old behavior equivalent to `gdna_upper > eps` is not present.

Consequences:

- A locus with `gdna_prior_count_em == 0` can still have `enable_gdna == 1`.
- A technically eligible gDNA component can absorb likelihood-supported mass even when calibration predicted zero or near-zero gDNA.
- All-spliced loci, or loci with no finite gDNA likelihoods, are technically ineligible.

## Python-To-Native Handoff

`pipeline.quant_from_buffer()` calls `assemble_priors()` after fragment scoring and `build_multi_loci()`. It passes these arrays into `_run_locus_em_partitioned()`:

```text
prior_table.gdna_prior_count_em
prior_table.gdna_eff_len
prior_table.enable_gdna
prior_table.gdna_eff_len_unweighted       diagnostics
prior_table.gdna_expected_count           diagnostics/raw count mirror
prior_table.gdna_em_exposure_weight       diagnostics
```

`AbundanceEstimator.run_batch_locus_em_partitioned()` then sends to native:

```text
locus_gdna_prior_count = gdna_prior_count_em
locus_enable_gdna      = enable_gdna
locus_gdna_eff_lens    = gdna_eff_len
t_eff_lens_arr         = per-transcript RNA effective lengths
use_vbem               = em_config.mode == "vbem"
```

If no explicit `enable_gdna` is supplied, the wrapper can recompute the same technical eligibility from partition data. In the live pipeline, `PriorTable.enable_gdna` is supplied explicitly.

## Native Locus Construction

For each locus, native EM builds `n_transcripts + 1` components:

```text
components 0..n_t-1 = RNA transcript components
component n_t       = single locus-level gDNA component
```

For each fragment/unit:

- RNA candidates are remapped to local transcript component IDs.
- A gDNA candidate is appended when the unit is unspliced and has finite gDNA log-likelihood.
- The gDNA candidate uses the scorer-supplied `gdna_log_lik`; native EM applies the locus-level `-log(gdna_eff_len)` correction separately.

Per-component effective lengths are:

```text
RNA:  log_eff_len[t]    = log(max(transcript_eff_len[t], 1.0))
gDNA: log_eff_len[gdna] = log(max(locus_gdna_eff_len, 1.0))
```

The initial `sub.prior` array is filled with `EM_PRIOR_EPSILON = 1.0e-10` so all RNA components, and an enabled gDNA component, are marked eligible. If `enable_gdna == 0`, the gDNA slot is set to `0.0` and therefore excluded from the eligibility mask. In the live pipeline, `enable_gdna == 0` also means no unit should have a finite gDNA candidate.

## Native Prior Vector And Warm Start

Native code then calls `compute_gdna_prior_and_warm_start()`.

### Warm Start

The warm start is not the modeled prior. It initializes EM state as:

```text
theta_init = unambiguous RNA totals
```

Then each ambiguous equivalence-class row is distributed across eligible components using the current coverage weights normalized within that row:

```text
theta_init[component] += row_coverage_weight_share(component)
```

This means technically eligible gDNA can receive initial mass from ambiguous unspliced rows even before iterations begin. That is initialization, not `prior_out`.

### Final Prior Vector

The final native prior vector is mode-dependent:

```text
baseline = 0.5 if VBEM else 0.0

prior[disabled component] = 0.0
prior[RNA component]      = baseline
prior[gDNA component]     = baseline + max(gdna_prior_count_em[locus], 1.0e-300)
```

In default `EMConfig`, mode is `"vbem"`, so every eligible RNA component receives `0.5`, and the eligible gDNA component receives:

```text
0.5 + gdna_prior_count_em[locus]
```

up to the numerical `1.0e-300` guard when the count is exactly zero.

In MAP mode, eligible RNA components receive `0.0`, and gDNA receives approximately:

```text
gdna_prior_count_em[locus]
```

again with only a numerical `1.0e-300` guard at zero.

## How The Prior Enters EM Updates

For MAP EM, the E-step log weight is:

```text
log_weight[c] = log(theta[c] + 1.0e-300) - log_eff_len[c]
```

and the M-step is:

```text
theta_new[c] = normalize(unambiguous_count[c] + ambiguous_posterior_count[c] + prior[c])
```

For VBEM, the E-step log weight is:

```text
log_weight[c] = digamma(alpha[c]) - digamma(sum(alpha)) - log_eff_len[c]
```

and the update is:

```text
alpha_new[c] = unambiguous_count[c] + ambiguous_posterior_count[c] + prior[c]
```

Thus the gDNA calibration count acts as a Dirichlet-style extra count for the gDNA component, while `gdna_eff_len` acts as the component-specific normalization/exposure denominator in the likelihood term.

## What Currently Does Not Happen

The current clean baseline does not do the following:

- It does not disable gDNA based on `upper_gdna`, `mu_gdna`, `A_r`, `p_captured`, or any posterior threshold.
- It does not use `rna_lower` to protect RNA components from gDNA siphoning.
- It does not construct a matched calibrated RNA prior count.
- It does not expose a gDNA aggressiveness knob.
- It does not use gene-level summaries or gene IDs in prior construction.
- It does not use `gdna_upper_count` in native EM; that field is diagnostic under the current baseline.

## Current Interpretation And Open Problem

The current prior has the following structure:

```text
gDNA numerator evidence:
    allocated RegionCalibration.mu_gdna count

gDNA denominator/exposure evidence:
    FL-marginal locus gDNA effective length weighted by mean A_r

technical availability:
    any unspliced finite-gDNA-likelihood unit

native generic baseline:
    0.5 for every eligible component in VBEM, 0.0 in MAP
```

This explains why the model can still leak mass into gDNA in zero-gDNA scenarios after removing the old `gdna_upper > eps` gate. The gDNA component is technically available, receives a warm start from likelihood-compatible unspliced ambiguity, receives the generic VBEM baseline, and may receive a positive calibration-derived count through `mu_gdna`. RNA does not receive an analogous calibrated prior count at the same projection layer.

The next design should therefore treat RNA and gDNA as a paired prior problem rather than only suppressing or enabling gDNA. A principled replacement likely needs:

- a calibrated RNA prior numerator and denominator,
- a calibrated gDNA prior numerator and denominator,
- an explicit operating-point/aggressiveness parameter,
- benchmark sweeps over true-zero gDNA and true-positive gDNA conditions.

## Code Map

- `src/rigel/calibration/background_model.py`: initial off-target Gamma background posterior with `alpha_floor = 1.0`, `beta_floor = 1.0`.
- `src/rigel/calibration/calibration_iteration.py`: four-state `RegionCalibration`, `mu_gdna`, `upper_gdna`, `A_r`, and iterative background refit.
- `src/rigel/calibration/prior.py`: projection from regional calibration arrays to `PriorTable` per `MultiLocus`.
- `src/rigel/calibration/_exposure.py`: gDNA effective length and bp-weighted `A_r` exposure weighting.
- `src/rigel/pipeline.py`: handoff from `PriorTable` to locus EM.
- `src/rigel/estimator.py`: Python wrapper for native `batch_locus_em_partitioned()`.
- `src/rigel/native/em_solver.cpp`: native gDNA component construction, warm start, prior vector, and MAP/VBEM updates.