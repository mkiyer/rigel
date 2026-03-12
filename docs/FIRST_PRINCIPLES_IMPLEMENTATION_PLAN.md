# Rigel First-Principles Implementation Plan

This document is the implementation-planning companion to:

- `docs/THEORETICAL_MODEL.md`
- `docs/THEORETICAL_EM_OBJECTIVE.md`

Detailed workstream expansion now lives in:

- `docs/rigel_revision_plan/README.md`

It is intentionally redesign-oriented rather than patch-oriented. The purpose is
to sequence implementation work from the theory outward, keeping conceptual
purity and minimizing accidental complexity.

This is an initial plan. It is expected to be revised repeatedly.

## 1. Redesign Goal

Implement a cleaner abundance model with the following core properties:

1. one collapsed gDNA abundance component per locus
2. one free locus-specific gDNA strand-balance parameter $\phi_\ell$
3. one global empirical-Bayes symmetry hyperparameter
   $\kappa_{\mathrm{sym}}$
4. calibration of gDNA nuisance parameters from soft-weighted, gDNA-dominant
   genomic regions
5. a clean separation among:
   - structural likelihoods
   - abundance parameters
   - biological priors
   - initialization or calibration procedures

## 2. Guiding Principles

### 2.1 Preserve the public problem definition

The output target remains:

- transcript-level mature RNA abundance
- nRNA-span-level or transcript-fanned-out nascent RNA abundance
- locus-level total gDNA abundance

No strand-resolved gDNA output is required.

### 2.2 Keep the generative story simple

The main locus model should stay close to the compact formulation in
`THEORETICAL_EM_OBJECTIVE.md`.

### 2.3 Move library-level nuisance learning upstream

The gDNA fragment-length model and gDNA symmetry hyperparameter should be
learned in a calibration stage, not improvised inside the locus EM.

### 2.4 Treat heuristics as approximations to latent statistical objects

If practical approximations are needed, they should correspond to clearly stated
latent targets such as posterior regional gDNA purity.

## 3. Scope of the First Redesign

The first redesign should focus on the following deliverables.

### 3.1 In scope

- redesign the gDNA symmetry model around free $\phi_\ell$ and global
  $\kappa_{\mathrm{sym}}$
- define and implement a region-based gDNA calibration stage
- define and implement a soft regional gDNA-purity weighting framework
- refactor gDNA fragment-length estimation to use the same calibration logic
- simplify the role of initialization versus prior in the main EM

### 3.2 Explicitly out of scope for the first pass

- redesign of all RNA scoring terms
- changes to the public file formats
- full rethinking of transcript or nRNA indexing
- aggressive optimization before the statistical structure is validated

## 4. Major Workstreams

The redesign naturally decomposes into six workstreams.

### Workstream A: Annotation-context partition

Goal:

Build or expose the genomic interval partition needed for calibration.

Required capabilities:

- partition each reference into intervals with constant annotation context
- preserve exon, intron, intergenic, and strand-specific annotation labels
- expose strand ambiguity explicitly
- expose effective interval lengths for density calculations

Primary output:

- a region table keyed by `(ref, start, end, context)`

### Workstream B: Calibration fragment extraction

Goal:

Map parsed fragment objects onto calibration regions and accumulate regional
evidence summaries.

Required capabilities:

- operate on fragment objects, not raw BAM records
- restrict to unspliced fragments for gDNA-specific calibration signals
- optionally use nearby spliced fragments as RNA evidence for the same region or
  neighborhood
- retain only context-unambiguous fragment-to-region assignments in the first
  implementation

Primary output:

- a per-region evidence table with counts and densities

### Workstream C: Regional evidence model

Goal:

Compute the three fragment-level RNA evidence channels at region level.

Evidence families:

1. splice evidence
2. density-contrast evidence
3. strand-asymmetry evidence

Primary output:

- region-level evidence summaries suitable for purity modeling

### Workstream D: Soft gDNA-purity model

Goal:

Estimate continuous regional gDNA purity or a practical approximation to it.

Preferred target:

- latent continuous purity variable $\pi_r \in [0,1]$
- posterior weight $w_r = \mathbb{E}[\pi_r \mid y_r]$

Fallbacks:

- two-state mixture model
- heuristic score explicitly framed as posterior approximation

Primary output:

- one weight per region for gDNA-dominance calibration

### Workstream E: Empirical-Bayes nuisance learning

Goal:

Learn library-level gDNA nuisance parameters from the weighted calibration set.

Primary targets:

- gDNA fragment-length distribution
- gDNA symmetry hyperparameter $\kappa_{\mathrm{sym}}$

Primary output:

- calibrated gDNA nuisance parameter bundle for the run

### Workstream F: Main locus EM integration

Goal:

Refit the locus model around the cleaner parameterization.

Primary changes:

- replace current gDNA symmetry treatment with free locus-specific $\phi_\ell$
- feed calibrated $\kappa_{\mathrm{sym}}$ into the gDNA strand prior
- ensure total gDNA abundance remains one collapsed locus component
- simplify or remove older gDNA-specific initialization machinery that is no
  longer needed once calibration is explicit

## 5. Proposed Build Order

The work should proceed in the following order.

### Phase 0: Freeze theory and terminology

Deliverables:

- finalized theory docs
- stable vocabulary for:
  - calibration region
  - context ambiguity
  - regional purity
  - gDNA-dominance weight
  - locus-specific strand-balance parameter

Exit criteria:

- no open ambiguity about the statistical target

### Phase 1: Region partition infrastructure

Deliverables:

- region partition schema
- mapping from existing index structures to calibration intervals
- clear policy for ambiguous annotation contexts

Exit criteria:

- genome can be deterministically partitioned into calibration-ready intervals

### Phase 2: Regional evidence extraction

Deliverables:

- per-region fragment counts
- spliced evidence tallies
- strand-resolved unspliced counts
- exonic and nonexonic effective-length-normalized densities

Exit criteria:

- a run can emit a region evidence table independent of the main EM

### Phase 3: Purity-model prototype

Deliverables:

- first working soft-weight model
- region weight diagnostics
- sanity plots or summaries showing that evidently RNA-rich regions receive low
  weights and evidently RNA-poor regions receive high weights

Exit criteria:

- region weights behave qualitatively correctly on synthetic and real data

### Phase 4: Empirical-Bayes calibration of gDNA nuisance parameters

Deliverables:

- weighted gDNA fragment-length estimation
- weighted empirical-Bayes estimate of $\kappa_{\mathrm{sym}}$
- self-consistent alternation if the purity weights depend on the current
  symmetry model

Exit criteria:

- calibrated nuisance parameters are stable across repeated runs and reasonable
  under known controls

### Phase 5: Main EM refactor

Deliverables:

- clean integration of calibrated nuisance parameters into the locus EM
- free locus-specific $\phi_\ell$ under symmetric Beta prior
- simplified gDNA initialization semantics

Exit criteria:

- main EM runs end-to-end with the new calibration stack

### Phase 6: Validation and ablation

Deliverables:

- synthetic validation suites targeted at:
  - gDNA symmetry dispersion
  - gDNA fragment-length recovery
  - robustness under low intergenic coverage
  - capture-enriched or exon-biased settings
- real-data diagnostics comparing old and new calibration behavior

Exit criteria:

- redesign improves interpretability and does not regress pristine-case accuracy

## 6. Detailed Design Targets by Workstream

### 6.1 Workstream A: Region partition design

The region table should contain at minimum:

- `region_id`
- `ref`
- `start`
- `end`
- `length`
- `is_genic`
- `has_exon_plus`
- `has_exon_minus`
- `has_intron_plus`
- `has_intron_minus`
- `strand_context`
- `annotation_context`
- `is_context_ambiguous`

Open design questions:

- whether to keep context as a compact enum or bitmask
- whether to merge tiny adjacent intervals with identical context for stability
- whether to precompute neighborhood relations for splice-evidence borrowing

### 6.2 Workstream B: Fragment-to-region assignment design

For each fragment, the calibration path needs:

- whether the fragment is spliced or unspliced
- its observed strand
- its genomic footprint
- the set of overlapped calibration regions

First-pass admissibility rule:

- retain only fragments with one interpretable calibration context
- hold out context-ambiguous fragments for the calibration stage

This should be a conscious first-order simplification, not an accidental loss of
information.

### 6.3 Workstream C: Splice evidence design

Splice evidence has two roles.

Direct role:

- any spliced fragment is pure RNA evidence

Indirect role:

- spliced signal in a transcriptional neighborhood implies expected unspliced
  RNA signal in that same neighborhood

The first implementation may use only the direct role if needed. A later
refinement can convert spliced evidence into an expected local RNA burden using
transcript geometry and fragment-length information.

### 6.4 Workstream C: Density evidence design

The density signal should be based on effective-length-normalized regional
fragment density.

Preferred first contrast:

- exonic density versus nonexonic density

where nonexonic is intronic plus intergenic or the subset available and reliable
in the given library design.

This choice explicitly acknowledges that modern annotations may leave little
truly intergenic space, whereas intronic space often remains the more informative
RNA-poor background.

### 6.5 Workstream C: Strand evidence design

For strand-resolved non-ambiguous regions, accumulate:

- positive-strand count
- negative-strand count
- total count

The strand-evidence term should not be a crude imbalance ratio alone. It should
be a divergence relative to the current symmetric Beta-Binomial reference model.

### 6.6 Workstream D: Continuous purity model candidates

Candidate family A: latent continuous purity regression

- latent $\pi_r \in [0,1]$
- evidence likelihoods conditioned on $\pi_r$
- posterior mean used as weight

Candidate family B: discretized latent mixture approximation

- region is gDNA-dominant versus contaminated
- posterior state probability used as weight

Candidate family C: heuristic posterior surrogate

- monotone transformation of evidence channels
- calibrated to approximate posterior purity ordering

Recommended starting target:

- design the system so that family C can be implemented first if needed, while
  preserving interfaces compatible with family A

### 6.7 Workstream E: Self-consistent empirical Bayes

The nuisance-learning step should support alternating refinement:

1. initialize provisional purity weights
2. estimate gDNA fragment-length distribution
3. estimate $\kappa_{\mathrm{sym}}$
4. recompute strand-based purity evidence
5. update weights
6. iterate until stable

This is likely the cleanest practical architecture even if the first pass uses
simple submodels internally.

### 6.8 Workstream F: Main EM changes

The main EM should receive from calibration:

- a gDNA fragment-length model
- a global $\kappa_{\mathrm{sym}}$

The main EM should estimate within each locus:

- transcript and nRNA abundance weights
- one total gDNA abundance weight
- one locus-specific strand-balance parameter $\phi_\ell$

This preserves the collapsed gDNA abundance representation while making the
symmetry model explicit and statistically clean.

## 7. Validation Strategy

Validation should be treated as part of the redesign, not as a final add-on.

### 7.1 Synthetic requirements

Need simulations that vary independently:

- total gDNA fraction
- global strand symmetry dispersion
- regional heterogeneity in gDNA strand balance
- RNA expression sparsity
- capture or exonic enrichment bias
- fragment-length separation between RNA and gDNA

### 7.2 Real-data diagnostic requirements

Need summary outputs showing:

- distribution of region weights
- distribution of estimated $\phi_\ell$
- fitted global $\kappa_{\mathrm{sym}}$
- calibration-set composition by annotation context
- comparison of gDNA fragment-length estimates before and after redesign

### 7.3 Ablation requirements

At minimum, compare:

- old model versus new model
- hard region filtering versus soft weighting
- fixed $\phi = 1/2$ versus free $\phi_\ell$
- fixed hand-tuned $\kappa_{\mathrm{sym}}$ versus empirical-Bayes estimate
- heuristic weights versus probabilistic purity approximation

## 8. Initial Implementation Recommendation

The most conservative path that still respects the theory is:

1. expose the region partition from existing index structures
2. build a standalone region evidence table from fragments
3. implement a simple soft purity approximation first
4. use it to calibrate weighted gDNA fragment lengths and
   $\kappa_{\mathrm{sym}}$
5. integrate free $\phi_\ell$ into the main locus EM
6. only then decide whether the purity model itself needs to be made more fully
   probabilistic

This sequencing keeps the first implementation tractable while preserving the
first-principles structure.

## 9. Open Questions for Future Iteration

This plan intentionally leaves several questions open for future turns.

1. What is the cleanest first implemented approximation to the latent continuous
   purity model?
2. How should local splice evidence be propagated to nearby unspliced intervals?
3. Should calibration-region evidence be accumulated per base interval, merged
   neighborhood, or annotation cluster?
4. How should context-ambiguous fragments be reincorporated, if at all, after
   the first pass?
5. What is the cleanest interface between calibration outputs and the main native
   EM code?
