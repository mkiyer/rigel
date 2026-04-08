# nRNA-gDNA Misclassification: Root Cause Analysis

**Date**: 2026-03-31
**Status**: Fix implemented and verified (1013/1013 tests pass)

## Problem Statement

Stress testing revealed 29 simulation runs where nRNA fragments were entirely misclassified as gDNA. The worst cases (SS=0.5–0.7, NTA=100–500) showed:
- nRNA_observed = 0 (expected ~4900)
- gDNA_observed = ~5000 (expected 0)
- Calibration γ = 0.91–0.94 (should be ~0)

## Root Cause

### Primary: Zero-Count Region Exclusion in Density Pathway

The density pathway in `calibrate_gdna()` computes λ_G (global gDNA background density) from the 10th percentile of unspliced fragment density across genomic regions. The eligibility filter:

```python
eligible = (stats["n_total"] > 0) & (stats["region_length"] > 0)
```

**excluded intergenic regions with zero reads** from the percentile calculation.

When there is no actual gDNA, intergenic regions have exactly 0 reads, so they are removed from the density distribution. This eliminates the "zero baseline" that the density pathway relies on to distinguish background from transcription.

With only gene-body regions remaining (all saturated with nRNA unspliced fragments), the 10th percentile returns a HIGH density value, causing massive gDNA overestimation.

#### Trace through for worst case (NO_GDNA_SS05):

| Step | What happens |
|------|-------------|
| **Regions** | 9 total: 2 intergenic (0,8) with 0 reads, 7 gene-body with density 0.51–0.58 |
| **Eligibility** | Regions 0,8 excluded (n_total=0) → only 7 high-density regions remain |
| **Percentile** | 10th percentile of [0.51, 0.58, 0.55, 0.56, 0.56, 0.55, 0.53] → **λ_G = 0.521** |
| **Blending** | w=(2×0.58−1)²=0.025, so 97.5% density → blended **λ_G = 0.508** |
| **Calibration** | γ = 0.937, α_gDNA = 4.685, α_RNA = 0.315 |
| **EM** | Strong gDNA prior + indistinguishable nRNA/gDNA likelihoods → all nRNA → gDNA |

#### With the fix (zero-count regions included):

| Step | What happens |
|------|-------------|
| **Density eligibility** | All 9 regions with length > 0 included |
| **Values** | [0, 0.51, 0.58, 0.55, 0.56, 0.56, 0.55, 0.53, 0] |
| **Weights** | [1k, 1k, 2.5k, 1k, 1.5k, 0.5k, 1.5k, 1k, **40k**] |
| **Zero-density** | Regions 0+8 = 41,000 bp = 83.7% of total length |
| **Percentile** | 10th percentile → **λ_G = 0.0** |
| **Calibration** | γ = 0.025, α_gDNA = 0.125, α_RNA = 4.875 |

### Secondary: Noisy SS Estimate from Low Spliced Observation Count

With NTA=500 and low mRNA (abundance=37), 98% of fragments are unspliced nRNA. Only ~19 fragments are spliced (mRNA with exon-exon junctions). SS is trained exclusively from spliced reads, giving:
- True SS=0.5 → Trained SS=0.5789 (±0.22 at 19 obs)
- True SS=0.7 → Trained SS=0.5263

The strand pathway formula `E[gDNA] = (N_anti − N_unspliced·(1−SS)) / (SS−0.5)` has denominator (SS−0.5) ≈ 0.08 at trained SS=0.58, amplifying any noise by ~12×. However, with the density fix, the strand pathway only contributes ~2.5% (w=0.025 at SS=0.58), so its overestimate is harmless.

### Tertiary: Fundamental Identifiability at SS=0.5

At true SS=0.5, nRNA and gDNA are **theoretically indistinguishable**: both are unspliced, both have 50/50 strand distribution, and both cover gene body regions. Even with perfect calibration (γ=0.025), the EM cannot discriminate at the fragment level. This explains why the SS=0.5 case still fails after the fix (the calibration is correct but the EM has no discriminating signal).

## Fix

**One-line conceptual change**: Use `density_eligible = (region_length > 0)` instead of the existing `eligible = (n_total > 0) & (region_length > 0)` for the density percentile calculation.

File: `src/rigel/calibration.py`, density pathway section (~line 546)

Before:
```python
d_unspliced = np.zeros(n_regions, dtype=np.float64)
pos = (region_length > 0) & eligible
d_unspliced[pos] = n_unspliced[pos] / region_length[pos]

if n_eligible >= _MIN_DENSITY_REGIONS:
    lambda_g_density = _length_weighted_percentile(
        d_unspliced[eligible], region_length[eligible], density_percentile)
```

After:
```python
density_eligible = region_length > 0
d_unspliced = np.zeros(n_regions, dtype=np.float64)
pos = density_eligible & (n_unspliced > 0)
d_unspliced[pos] = n_unspliced[pos] / region_length[pos]

if int(density_eligible.sum()) >= _MIN_DENSITY_REGIONS:
    lambda_g_density = _length_weighted_percentile(
        d_unspliced[density_eligible], region_length[density_eligible], density_percentile)
```

Zero-count regions naturally have d_unspliced=0, providing the "zero baseline" the percentile needs.

## Verification Results

### Diagnostic Cases (6-case comparison matrix)

| Case | Before γ | Before nRNA/gDNA | After γ | After nRNA/gDNA | Status |
|------|----------|-------------------|---------|------------------|--------|
| NO_GDNA_SS05 | 0.937 | 0 / 4959 | **0.025** | 0 / 4965 | Cal fixed, EM limited |
| NO_GDNA_SS07 | 0.933 | 0 / 4947 | **0.000** | **4929 / 0** | **FIXED** |
| NO_GDNA_SS09 | 0.274 | 4204 / 721 | **0.088** | 4213 / 716 | Improved |
| NO_GDNA_SS10 | 0.000 | 4944 / 0 | 0.000 | 4933 / 0 | No change |
| CONTROL_SS05 | 0.334 | 5000 mRNA | **0.000** | 5000 mRNA | **FIXED** |
| GDNA_ONLY | 0.224 | 0/1498 | 0.224 | 18/1525 | No change |

### Golden Test Regression Analysis

6 golden output tests changed, all improvements:
- **antisense_contained_ss90**: gDNA dropped from 514→0 in a zero-gDNA scenario (was the old density bug)
- 5 others: tiny numerical shifts (<0.35% relative error)

All 1013 tests pass after golden update.

## Remaining Limitations

### SS=0.5 Identifiability Wall
At true SS=0.5, nRNA and gDNA produce identical observable signatures (unspliced + unstranded). No calibration fix can overcome this — it requires discriminating signals the EM doesn't currently use, such as:
1. **Fragment length distribution**: nRNA follows RNA FL; gDNA may differ
2. **Coverage profile**: nRNA has 5'→3' bias; gDNA is uniform
3. **Proximity to promoters**: nRNA concentrates near TSS

### SS=0.9 Residual from Strand Noise
With only 19 spliced observations, the strand pathway SS estimate is noisy. This creates small residual gDNA overestimation (~700/5000 fragments at SS=0.9). More spliced reads → better SS estimate → smaller residual.

## Impact on Real Data

- **Stranded libraries (SS>0.7)**: Major improvement — density pathway no longer overestimates for nRNA-rich loci
- **Unstranded libraries (SS≈0.5)**: Improved calibration but EM identifiability remains a fundamental limit
- **Real gDNA present**: No change — gDNA fills intergenic regions → all regions eligible → no exclusion
