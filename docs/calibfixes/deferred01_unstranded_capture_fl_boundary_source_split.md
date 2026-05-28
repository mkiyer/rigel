# Deferred PR: Unstranded Capture-On FL Plus Boundary Source Split

## Status

Deferred from the v5 calibration-prior repair. Do not fold this into PR 1, PR 2, or PR 3. Unstranded capture-on data lacks strand evidence, and expressed captured exons can overwhelm density evidence. It needs a separate source-identification model.

## Goal

Produce strand-free regional gDNA/RNA source mass for capture-on data by combining fragment length, boundary/intronic compatibility, and capture geometry. The output should feed the same `PriorMassDeconvolution` contract used by the stranded repair:

```text
gdna_unspliced_mean + rna_unspliced_mean = unspliced_total
```

## Non-goals

- Do not use latent states as source labels.
- Do not use capture exposure alone as source evidence.
- Do not classify expressed captured exon mass as gDNA solely because exposure is high.
- Do not replace transcript-level EM likelihoods with this regional model.

## Files likely to edit

| File | Purpose |
|---|---|
| `src/rigel/calibration/fl.py` and `_fl_sources.py` | Build or expose captured-gDNA FL training pools. |
| `src/rigel/calibration/density_observation.py` | Provide per-compartment counts/effective lengths needed for a strand-free mixture. |
| `src/rigel/calibration/capture_exposure.py` | Reuse target weights and source-reliable exposure diagnostics from PR 2. |
| `src/rigel/calibration/unstranded_source_split.py` | New module for FL plus boundary regional source split. |
| `src/rigel/calibration/calibration_iteration.py` | Use strand-free source split when strand reliability is inactive and capture opportunity is present. |
| `tests/test_unstranded_capture_source_split.py` | Unit tests for model behavior. |

## Model sketch

For a fragment or regional evidence unit `i`:

```text
P(i | source = RNA)  uses RNA FL + transcript/exonic compatibility
P(i | source = gDNA) uses captured-gDNA FL + genomic/boundary compatibility + capture opportunity
```

At the region level:

```text
logit P(gDNA | i) = log prior_gdna_exposure_i
                   + log P(FL_i | gDNA_capture)
                   + log P(boundary_i | gDNA)
                   - log P(FL_i | RNA)
                   - log P(exonic_i | RNA)
```

Aggregate fragment posteriors to region/compartment source mass. If only regional count summaries are available, use a binned FL mixture rather than a per-fragment model.

## Training captured-gDNA FL

Candidate training evidence:

1. source-reliable stranded capture-on boundary fragments from PR 1 + PR 2;
2. off-target intronic/intergenic capture-on fragments with high gDNA posterior;
3. synthetic oracle labels for diagnostics only, never production fitting.

Use empirical-Bayes shrinkage toward global gDNA FL when the captured-gDNA pool is sparse. Emit pool quality diagnostics just like existing FL models.

## Boundary evidence

Boundary-crossing and intronic compatibility are required. A fragment spanning exon-intron boundary or lying in intronic sequence is direct evidence against mature RNA. The model should separate at least:

```text
contained exonic-compatible
boundary left
boundary right
intronic/intergenic
```

Do not collapse these into a single density count before source splitting.

## Acceptance tests

1. Unstranded capture-off remains unchanged or improves only through density/FL evidence.
2. Unstranded capture-on no-gDNA does not produce large false gDNA in highly expressed captured exons.
3. Unstranded capture-on high-gDNA recovers gDNA mass better than PR 1 + PR 2 alone.
4. Boundary-only gDNA evidence is recognized even when contained exon expression is high.
5. Sparse captured-gDNA FL training falls back smoothly to global gDNA FL and emits weak-pool diagnostics.

## Benchmark gate

This deferred PR should be judged on the unstranded capture-on scenarios that PR 1 + PR 2 explicitly do not promise to solve. It must also preserve the seven other strata/condition groups.
