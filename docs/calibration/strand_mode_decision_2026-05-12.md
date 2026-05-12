# Strand Mode Decision Rule

Status: implemented, 2026-05-12.

This note is the authoritative description of when calibration uses strand
orientation data and when it runs the unstranded count/exposure estimator.
Unstranded RNA-seq is a first-class supported mode. Do not describe it as a
legacy mode.

## User-facing parameters

There is no longer a public `--min-strand-contrast` CLI parameter, and
`min_strand_contrast` is no longer part of `CalibrationConfig`.

Before this cleanup, that parameter existed with default `1e-3`. The benchmark
regression showed that this was the wrong user-facing abstraction: deciding
whether strand data are usable is a statistical identifiability question, not a
knob users should tune per run. Rigel now makes the decision automatically.

The value `1e-3` remains only as the internal constant
`STRAND_CONTRAST_NUMERICAL_FLOOR`. It prevents division by a nearly zero signed
contrast after the statistical check. It is not a CLI or YAML setting.

## Strand summary

Calibration uses a compact `StrandSummary`:

```text
p = p_r1_sense
n = n_observations
signed_strand_contrast = c = 2p - 1
SE(c) = 2 * sqrt(p * (1 - p) / n)
margin_99 = 2.575829 * SE(c)
```

If `n == 0`, `SE(c)` is infinite and the summary is not identifiable.

A summary is identifiable only if:

```text
abs(c) >= max(STRAND_CONTRAST_NUMERICAL_FLOOR, margin_99)
```

This means the 99% normal-approximation interval for signed strand contrast
does not include zero, and the contrast is not numerically tiny.

## Which strand model calibration uses

Rigel builds two relevant strand summaries during the BAM scan:

1. `exonic_spliced`: the primary model, trained from annotated spliced RNA
   evidence. This remains the authoritative model for scoring.
2. `exonic`: a diagnostic model from all exonic strand evidence. This can be
  useful for QC, but it is not RNA-pure and is not used for calibration.

Calibration chooses the summary as follows:

1. Always use `exonic_spliced` as the calibration strand summary.
2. If `exonic_spliced` is identifiable by the rule above, strand-corrected
   channels may use it.
3. If `exonic_spliced` is not identifiable, calibration does not substitute any
   unspliced or exonic diagnostic model. Channels run unstranded unless future
   user-declared protocol validation is added.

The diagnostic `exonic` model is never used for calibration or scoring. This is
intentional. Annotated spliced fragments are the only RNA-pure strand-training
source available to Rigel. Unspliced exonic fragments can include gDNA and
nascent RNA, so using them to rescue calibration would make high-contamination
or atypical libraries look artificially well-stranded. If `exonic_spliced` is
not identifiable but the diagnostic `exonic` model is strongly strand-specific,
Rigel warns: this pattern can indicate sparse splice evidence with high gDNA or
nascent RNA signal, and the run should be inspected rather than silently
rescued.

## Per-channel density decision

Each global gDNA density channel makes its own final decision:

- `INTERGENIC` has no transcript strand, so it always uses the unstranded
  count/exposure estimator.
- `INTRON` and `EXON-INTRON` can use strand correction only on exposure whose
  region strand is `POS` or `NEG`.

For a channel, define:

```text
T = n_same + n_opp
D = n_same - n_opp
c = signed_strand_contrast
```

The channel uses strand correction only when both are true:

1. The selected strand summary is identifiable.
2. The channel has positive strand-identifiable exposure.

When active, the gDNA count moment is:

```text
G = sum(T - D / c) over strand-identifiable rows or boundary sides
E = sum(exposure) over the same strand-identifiable rows or boundary sides
rho = max(0, G / E)
```

Only the final aggregate ratio is clipped at zero. Individual signed moments
are not clipped per row.

When inactive, the channel uses the unstranded count/exposure estimator:

```text
G = sum(raw observed counts) over all eligible rows or boundary sides
E = sum(exposure) over all eligible rows or boundary sides
rho = G / E
```

Unknown-strand observations on strand-identifiable rows are diagnostic data.
They are reported in `n_uninf_observed` and reflected in `rho_uncorrected`, but
they are not added as raw gDNA counts to a strand-corrected numerator.

## Summary diagnostics

`summary.json` reports the channel-level decision and audit fields:

- `strand_active`: whether that density channel used the signed correction.
- `rho_uncorrected`: the unstranded count/exposure density for comparison.
- `strand_correction_fragments`: projected corrected count minus raw observed
  count.
- `n_strand_informative_regions`: rows or boundary sides with usable region
  strand.
- `strand_informative_exposure_fraction`: fraction of eligible exposure that
  was strand-identifiable.
- `n_uninf_observed`: observations whose fragment or region strand was not
  informative.

For unstranded libraries, `strand_active` should be false for `INTRON` and
`EXON-INTRON`, and the reported density is the unstranded count/exposure
estimate. This behavior is intentional support for unstranded data.