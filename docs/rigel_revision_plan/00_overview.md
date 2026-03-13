# Revision Overview

This document converts the high-level redesign into a concrete dependency-aware
implementation roadmap.

## 1. Core Observation from the Current Code

Several ingredients needed by the redesign already exist.

### 1.1 Region substrate is partially present already

`src/rigel/index.py` already builds a tiled interval representation of the
genome in `intervals.feather`, including exon, transcript-span, unambiguous
intron, splice-junction, and intergenic structure.

This means Workstream A can build a dedicated flattened calibration-region
index from existing transcript and exon coordinates rather than inventing a
separate annotation source from scratch.

### 1.2 Fragment-level calibration signals are already present during scan

The scan and buffer path already exposes per-fragment information sufficient for
the first calibration pass:

- splice state
- exon strand
- genomic start
- genomic footprint
- exonic overlap summary
- intronic overlap summary
- unambiguous intronic overlap summary

These appear in the chunk/buffer representation and are already consumed by the
native scan path.

### 1.3 Current gDNA calibration logic is localized enough to replace cleanly

The existing gDNA prior path is largely concentrated around:

- `src/rigel/locus.py`
- `src/rigel/priors.py`
- `src/rigel/pipeline.py`

That is favorable. It implies the redesign can first add a new upstream
calibration product without immediately disturbing the rest of the EM stack.

## 2. Recommended Execution Strategy

The codebase suggests a stronger recommendation than the earlier theory-only
plan.

The first implementation should not start with the purity model itself. It
should start by creating the observable regional evidence layer and a highly
conservative calibration seed set.

That order is preferable because:

1. the evidence substrate is already close to implementable
2. the purity approximation should be informed by real observed region-level
   distributions rather than assumed up front
3. the same evidence table will support both gDNA symmetry and gDNA
   fragment-length calibration
4. it avoids prematurely hard-coding a purity model that later turns out to be
   misaligned with the empirical data

The updated recommendation is more specific:

- build and persist a flattened calibration-region index as `regions.feather`
- derive that index from transcript-span and exon boundaries at index-build
   time
- store exactly four region flags: `exon_pos`, `exon_neg`, `tx_pos`, `tx_neg`
- use a two-state purity model as the first practical implementation target
- estimate gDNA fragment length and $\kappa_{\mathrm{sym}}$ jointly from the
  same purity-weighted regions
- preserve the one-sided robustness insight from the current targeted-excess
  penalty when designing the new gDNA-pair solver constraint

## 3. Public Model Versus Internal Parameterization

Publicly, Rigel should still report one collapsed gDNA abundance per locus.

Internally, the chosen implementation target is:

1. a `T + N + 2` model with `g_pos` and `g_neg` as ordinary mixture
   components
2. public reporting that collapses them to total locus gDNA abundance

The March 2026 methodology review makes a strong engineering case for this
choice because it keeps the EM closer to a textbook simplex update, avoids a
custom collapsed-gDNA nuisance-parameter M-step, and makes the strand-symmetry
prior a coupling prior on an ordinary pair of components.

The current recommendation is therefore:

- keep the public collapsed gDNA output unchanged
- implement the solver around `T + N + 2`
- reserve the collapsed one-gDNA model as a marginal reference model rather
   than the preferred implementation target

## 4. Workstream Dependencies

The practical dependency graph is:

1. Workstream A produces a flattened calibration-region table
2. Workstream B/C produces a per-region evidence table
3. Workstream D maps region evidence to gDNA-dominance weights
4. Workstream E consumes those weights to estimate gDNA nuisance parameters
5. Workstream F injects calibrated nuisance parameters into the locus EM

In code terms, the upstream interface that matters most is:

- input: fragment stream plus region table
- output: one per-region evidence record with counts, densities, and optional
  neighborhood summaries

Once that exists, the purity model can be iterated independently.

## 5. Definitive Versus Deferred Workstreams

The current repository supports a fairly definitive plan for these workstreams.

### 4.1 Definitive now

- Workstream A: flattened region-index build and load path
- Workstream B/C: fragment-to-region evidence extraction
- Workstream E: weighted calibration interface and estimators
- much of Workstream F: integration points and deprecation targets
- conservative calibration seed-set definition
- bridge analysis from current targeted-excess penalty to target prior model

### 4.2 Intentionally deferred

- the exact feature family and likelihood parameterization inside the first
   two-state purity model
- whether splice evidence should borrow from a strict overlap window or a wider
  neighborhood model

## 6. Proposed Milestones

### Milestone 1: Standalone regional evidence export

Deliver a run mode or internal function that emits region evidence without
changing locus EM behavior.

Success criteria:

- deterministic `regions.feather` exists and loads through `TranscriptIndex`
- fragments can be summarized by region
- evidence summaries are inspectable on synthetic and real data

### Milestone 2: Seeded gDNA calibration bundle

Deliver a calibration object that contains at least:

- `kappa_sym`
- gDNA fragment-length model or weighted fragment-length histogram
- diagnostics describing the effective calibration mass

Success criteria:

- stable across repeated runs
- interpretable on low-RNA controls and simulated libraries
- initially supported by a deliberately conservative calibration seed set

### Milestone 3: Main EM migration

Replace the current embedded gDNA heuristics with calibrated nuisance inputs.

Success criteria:

- one collapsed gDNA abundance per locus remains unchanged publicly
- the internal `T + N + 2` solver parameterization is implemented explicitly
   and justified
- older gDNA init semantics are reduced to calibration-informed gating or
  removed entirely