# Prior-Mass Handoff Root-Cause Plan 2026-05-27

## Position

The capture-gated prior-mass handoff is not the right design. It creates two semantic modes:
captured regions use `PriorMassDeconvolution`, while off-target regions can still infer gDNA/RNA
mass from four-state latent labels. That is not theoretically clean.

The clean invariant should be:

1. `PriorMassDeconvolution.gdna_unspliced_mean` and `rna_unspliced_mean` are the only source of
   regional gDNA/RNA mass split used by adaptive EM priors.
2. Four-state latent probabilities may contribute confidence/entropy weighting, capture/expression
   diagnostics, and state-level calibration, but they must not be reinterpreted as pool-mass labels.
3. Any regression from this single approach must be explained at the calibration, prior-strength, or
   EM/scoring layer; it should not be hidden by switching mass semantics by capture status.

## Diagnostics Run

Generated report:

```
/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb/diagnostics/prior_mass_handoff_root_causes/
```

Key files:

- `snapshot_metrics.tsv`: baseline state split vs all-region prior mass vs capture-gated prior mass.
- `regional_truth_vs_signals.tsv`: oracle regional truth vs `prior_mass` vs state-implied
  unexpressed mass.
- `locus_prior_summary.tsv`: locus-level prior strength and EM outcomes for all snapshots.
- `capture_off_transcript_regression_focus.tsv`: transcript-level source of the high-SS capture-off
  all-region regression.

## Root Causes

### 1. Latent state labels are not gDNA/RNA mass labels

The four-state model now has `unexpressed_offtarget`, `unexpressed_capture`,
`expressed_capture`, and `expressed_offtarget`. It has no mixed expressed-plus-gDNA state.
Therefore, a region can be labeled `expressed_capture` while containing large gDNA mass.

Observed in high-gDNA SS0.99 capture-on:

- True regional gDNA: 100,000
- `prior_mass` gDNA: 98,698
- State-implied unexpressed mass: 5,809
- gDNA in expressed states: 93,748
- Probe-overlap exon true gDNA: 93,146
- Probe-overlap exon `prior_mass` gDNA: 92,315
- Probe-overlap exon state-implied unexpressed mass: 1.45

This proves the state-derived mass split is the wrong abstraction. The adaptive prior must consume
the deconvolved prior mass directly.

Solved by the single approach: high-gDNA SS0.99 capture-on changes from 39.01% estimated gDNA to
41.52% under all-region prior mass; gDNA-to-RNA assignment drops from 24.73% to 19.98%; mRNA MARD
drops from 674.9% to 319.3%.

### 2. Unstranded capture-on is an upstream calibration deconvolution failure

For high-gDNA SS0.50 capture-on, both state labels and `prior_mass` are wrong in captured expressed
probe exons.

- True regional gDNA: 100,000
- `prior_mass` gDNA: 2,834
- State-implied unexpressed mass: 6,400
- gDNA in expressed states: 93,244
- Probe-overlap exon true gDNA: 93,126
- Probe-overlap exon `prior_mass` gDNA: 130
- Probe-overlap exon state-implied unexpressed mass: 7.77

This is not fixable by choosing a different handoff. The calibration layer is failing to deconvolve
captured expressed unstranded RNA/gDNA mixtures. The downstream prior can only propagate the wrong
mass split.

Solved by future calibration work, not by adaptive-prior switching.

### 3. The all-region clean handoff regression is a prior-strength/EM interaction, not evidence that
latent mass labels are correct

The all-region prior-mass handoff improves pool-level gDNA calls in high-gDNA capture-off scenarios:

- SS0.99 capture-off gDNA delta improves from -1,739 to -1,184.
- SS0.50 capture-off gDNA delta improves from -2,145 to -1,645.

The concerning regression is transcript-level in SS0.99 capture-off:

- mRNA MARD worsens from 46.7% to 87.0%.
- Spearman worsens from 0.935 to 0.754.
- Absolute transcript error rises from 3,008 to 5,268.

The regression is localized, especially in GENE0008/locus 7:

- `GENE0008.2` true count 340; baseline 824; all-region 1,972.
- `GENE0008.3` true count 4,119; baseline 1,047; all-region 0.
- Locus 7 local prior gDNA/RNA changes from 245/12,746 to 2,396/10,595.
- Locus 7 additive alpha changes from 78 gDNA / 2,922 RNA to 566 gDNA / 2,434 RNA.
- The locus-level pool call moves only slightly, but isoform allocation collapses for one transcript.

This means the all-region mass split exposes a separate weakness: the current adaptive prior strength
can be too influential for transcript-level EM in an isoform-ambiguous locus. The right response is to
diagnose and fix prior strength, uncertainty, and/or EM likelihood balance. The wrong response is to
infer mass from latent labels in off-target regions.

## Proposed Clean Implementation Direction

### Phase 1: Clean the handoff semantics

- Remove `prior_mass_split_weight`.
- Remove the state-derived gDNA/RNA fallback from production adaptive-prior assembly.
- Keep latent-state entropy as a confidence weight only.
- Feed `prior_mass.gdna_unspliced_mean` and `prior_mass.rna_unspliced_mean` everywhere.
- Update tests so they assert the invariant: mass split comes from `PriorMassDeconvolution`, not
  `STATE_IS_EXPRESSED`.

Expected result: one simple, intuitive implementation path for all eight scenarios.

### Phase 2: Diagnose prior strength independently of mass semantics

Run controlled prior-strength sweeps using all-region prior mass only:

- `MAX_ESS`: current 3000 vs smaller values.
- entropy-only weighting vs entropy times `prior_mass.precision` where available.
- local-only prior vs local-plus-global shrinkage for high-ambiguity loci.
- per-locus impact diagnostics for GENE0008/locus 7 and neighboring high-error loci.

Acceptance criterion: no major transcript-level regression in SS0.99 capture-off while preserving the
SS0.99 capture-on gDNA rescue.

### Phase 3: Root-cause unstranded capture-on calibration

For high-gDNA SS0.50 capture-on, prior mass itself is wrong. The diagnostic target is the captured
expressed probe-exon group where true gDNA is 93,126 but `prior_mass` gDNA is only 130.

The next calibration work should test:

- whether density `mu_gdna` is suppressed by captured expressed RNA density;
- whether the capture/background model needs an explicit mixed expressed-plus-gDNA state;
- whether fragment-length evidence can supply a pool split when strand evidence is unavailable;
- whether probe-overlap/exon groups need a bounded gDNA floor derived from off-target or intergenic
  calibration rather than expression-state labels.

Acceptance criterion: prior-mass gDNA in unstranded high-gDNA capture-on probe exons moves toward
oracle truth before adaptive priors are involved.

## Benchmark Interpretation

The all-region prior-mass handoff is the theoretically correct semantic direction. It is already
better for pool-level deconvolution in the high-gDNA strand-specific scenarios, including capture-off.
The observed capture-off transcript regression is real, but it is not evidence for preserving latent
state mass splitting. It is evidence that prior ESS/shrinkage/likelihood balance needs its own root
cause analysis under the clean mass model.

The plan is therefore to replace the gated workaround with one prior-mass path, then solve the two
remaining issues at their actual layers: prior strength for isoform-ambiguous loci, and unstranded
capture-on calibration deconvolution.