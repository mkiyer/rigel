# gDNA-Everywhere Calibration Design 2026-05-27

## Core Reframe

The calibration latent states should not be interpreted as source labels. gDNA is present in every
region and every state. The four states should stratify regions by two orthogonal properties:

- RNA expression: unexpressed vs expressed.
- Capture exposure: off-target vs capture-enriched.

Under this interpretation, the current state names are misleading:

| Current label | Better semantic label | Meaning |
| --- | --- | --- |
| `background` | `unexpressed_offtarget` | no RNA expression signal, off-target density |
| `gdna_only_capture` | `unexpressed_capture` | no RNA expression signal, capture-enriched density |
| `expressed_offtarget` | `expressed_offtarget` | RNA expression signal, off-target density |
| `expressed_capture` | `expressed_capture` | RNA expression signal, capture-enriched density |

All four states can have gDNA. The state chooses a region's RNA/capture stratum; the pool split comes
from a separate deconvolution of unspliced mass into gDNA and RNA.

## Model Invariant

For each region `r`, calibration should estimate two nonnegative masses:

```text
U_r = G_r + R_r
```

where:

- `U_r` is observed compatible unspliced mass.
- `G_r` is gDNA mass and is allowed in every latent state.
- `R_r` is unspliced RNA mass and is concentrated in expressed states.

The adaptive EM prior consumes only:

```text
E[G_r | data], E[R_r | data], confidence_r
```

Latent state probabilities can inform `confidence_r`, exposure, and diagnostics, but they must not be
collapsed into `G_r` versus `R_r`.

## Observation Channels

The deconvolution should combine independent evidence channels rather than forcing a state label to
stand in for source identity.

### Strand-Specific Evidence

In strand-specific data, orientation channels are the strongest local source split. This is why the
high-gDNA SS0.99 capture-on case has accurate regional prior mass even though the dominant latent
label is `expressed_capture`.

Target behavior:

```text
G_r = strand-deconvolved gDNA mean
R_r = U_r - G_r
confidence_r = strand posterior precision, discounted by state entropy and sparse-boundary flags
```

### Density And Capture Evidence

Density evidence estimates expected gDNA mass from background/off-target rate, capture enrichment,
effective length, and probe overlap. It should estimate `G_r` in all four expression/capture strata.

Target behavior:

```text
E[G_r | density] = function(rho_off, capture_exposure_r, contained_leff_r, local boundary evidence)
```

Expression should not zero out density-derived gDNA. It should only explain the residual unspliced
mass as RNA after gDNA is accounted for.

### Spliced Evidence

Spliced fragments are direct evidence for RNA expression, not direct evidence against gDNA. They
increase confidence that `R_r > 0`, but they do not imply `G_r = 0`.

Target behavior:

```text
spliced evidence -> expression probability and RNA lower bound
spliced evidence !-> gDNA mass suppression
```

### Fragment-Length Evidence

In unstranded capture-on mixtures, strand evidence is unavailable and density can be confounded by
captured RNA. Fragment-length evidence should become a local source-split channel whenever RNA and
gDNA FL distributions are separable enough.

Target behavior:

```text
E[G_r | FL] = local mixture posterior using global RNA/gDNA FL models
confidence_r = separation-aware posterior precision
```

## Proposed Architecture

### 1. Rename Or Alias State Semantics

Keep the four-state tensor shape for compatibility, but expose semantic aliases in code and reports:

```text
unexpressed_offtarget
unexpressed_capture
expressed_offtarget
expressed_capture
```

The old `gdna_only_capture` name should be retired from diagnostics because it implies source purity
that is false.

### 2. Make `PriorMassDeconvolution` The Source-Split Contract

The calibration-to-EM handoff should be:

```text
prior_mass.unspliced_total
prior_mass.gdna_unspliced_mean
prior_mass.rna_unspliced_mean
prior_mass.precision
prior_mass.method
prior_mass.flags
```

Adaptive priors should use this contract everywhere. There should be no production fallback that
reconstructs source mass from state labels.

### 3. Separate Mean From Strength

The posterior mean split and prior strength should be independent concepts.

Mean:

```text
n_gdna_region = gdna_unspliced_mean
n_rna_region = rna_unspliced_mean
```

Strength:

```text
w_r = state_entropy_confidence_r * source_split_precision_r * structural_confidence_r
```

Then:

```text
weighted_gdna_region = w_r * n_gdna_region
weighted_rna_region = w_r * n_rna_region
```

This directly addresses the high-gDNA SS0.99 capture-off regression: the all-region mean is more
truthful, but prior strength can over-influence isoform-ambiguous loci. The fix should discount the
strength when source precision, structural confidence, or EM identifiability is weak; it should not
change the mass semantics.

### 4. Keep Pool Priors Pool-Level

The adaptive prior should constrain the gDNA-vs-RNA pool split without distorting isoform allocation
inside the RNA pool. If a locus has many nearly indistinguishable isoforms, a stronger gDNA prior
should not arbitrarily collapse one RNA isoform into another.

Design target:

- apply `alpha_gdna_add` to the gDNA component;
- apply `alpha_rna_add` as RNA-pool support in a way that is neutral across compatible expressed RNA
  components;
- audit whether current EM prior injection couples pool prior strength to isoform choice.

GENE0008/locus 7 in high-gDNA SS0.99 capture-off is the sentinel case for this behavior.

### 5. Add A Mixed-Mass Calibration Diagnostic

Every region report should expose:

```text
p_expressed
p_captured
prior_mass_gdna_fraction
state_implied_gdna_fraction
state_mass_disagreement
source_split_method
source_split_precision
```

`state_mass_disagreement` should be a diagnostic only. It should flag where state labels and source
mass strongly disagree, which is expected in expressed regions containing gDNA.

## High-gDNA Strand-Specific Focus

### Capture On: Remaining Failure Mode

Observed behavior with all-region prior mass:

- Regional `prior_mass` is accurate: true gDNA 100,000; prior gDNA 98,698.
- Adaptive prior rescue is large: estimated gDNA rises from 39.01% to 41.52%.
- Final EM still undercalls gDNA relative to true 50%.

Interpretation:

The handoff mean is no longer the main failure. The remaining gap is likely downstream:

- gDNA likelihood/exposure may be too weak in captured expressed exons;
- RNA likelihood may dominate unspliced fragments despite accurate source prior;
- prior ESS may still be too low relative to ambiguous fragment likelihood mass;
- gDNA effective length or capture exposure may be mismatched to the probe-enriched geometry.

Next diagnostics:

- per-fragment RNA-vs-gDNA log-likelihood margins in probe-overlap exons;
- per-locus prior-vs-likelihood posterior movement;
- gDNA effective length and exposure audit for captured exons;
- compare EM posterior gDNA fraction against regional `prior_mass` fraction by locus.

### Capture Off: Sentinel Regression

Observed behavior with all-region prior mass:

- Pool gDNA improves: gDNA delta moves from -1,739 to -1,184.
- Transcript allocation regresses in GENE0008/locus 7.

Interpretation:

The source mean is better, but the current prior strength or RNA-pool injection can perturb isoform
allocation in an ambiguous locus. This is not a reason to restore latent-label mass splitting. It is a
reason to decouple pool prior strength from isoform choice.

Next diagnostics:

- run prior ESS sweeps under all-region `prior_mass` only;
- inspect EM equivalence classes for GENE0008/locus 7;
- compare transcript posterior changes with and without gDNA prior while holding RNA-pool support
  neutral;
- audit whether `alpha_rna_add` is distributed in a way that favors one isoform.

## Implementation Phases

### Phase A: Clean Handoff

- Remove the capture-gated split weight.
- Remove any production state-derived source-mass fallback.
- Require `PriorMassDeconvolution` arrays in adaptive prior computation.
- Keep state entropy as confidence only.

### Phase B: Confidence-Weighted Source Priors

- Define density precision for density-derived `prior_mass`, not only strand-derived mass.
- Combine source precision with state entropy and sparse-boundary flags.
- Benchmark ESS/precision sweeps on high-gDNA SS0.99 capture on/off.

### Phase C: Pool-Neutral RNA Prior Injection

- Audit how `alpha_rna_add` affects RNA isoform allocation.
- Redesign prior injection if needed so RNA support is pool-level and isoform-neutral.
- Use GENE0008/locus 7 as the sentinel regression test.

### Phase D: Unstranded Capture-On Deconvolution

- Add FL-informed local source splitting for unstranded mixtures.
- Consider an explicit mixed expressed-plus-gDNA observation model.
- Add a calibration acceptance test where captured expressed probe exons contain known ~50% gDNA.

## Acceptance Criteria

For high-gDNA SS0.99 capture-on:

- regional `prior_mass` remains near truth;
- final gDNA fraction improves beyond the current all-region result;
- gDNA-to-RNA assignment decreases without increasing mRNA-to-gDNA materially.

For high-gDNA SS0.99 capture-off:

- pool gDNA improvement from all-region `prior_mass` is preserved;
- GENE0008/locus 7 no longer shows major isoform collapse;
- transcript Spearman and MARD stay near baseline or improve.

For the full eight-condition suite:

- no production mode switches by capture status;
- no state-label-derived source mass;
- no major regression masked by special-case logic.