# Transcript-Centric Cleanup Plan

**Date**: 2026-03-31
**Status**: Draft for review

---

## Why this cleanup is needed

The codebase has drifted into an inconsistent middle state:

1. The intended model is transcript-centric, with gene outputs as a convenience layer.
2. The current implementation still encodes nascent RNA partly as synthetic transcripts appended to the transcript table.
3. That representation leaks gene ownership into places where the annotation is not trustworthy and where many-to-many relationships are the real biological situation.
4. Annotated single-exon nascent-equivalent transcripts are treated as metadata in some places and as model objects in others.
5. Lower-level API contracts and scheduling heuristics no longer match the real pipeline requirements.

The goal of this cleanup is not cosmetic refactoring. The goal is to make the core data model, solver inputs, and reporting layers agree on one coherent design.

The previous version of this plan got one central point wrong: it tried to separate annotated single-exon transcripts into distinct mature and nascent objects. That would be a modeling regression, not an improvement. If two states have the same genomic span and no likelihood term can distinguish them, splitting them creates a non-identifiable EM problem, increases state count, and forces arbitrary mass-splitting. That is worse for both accuracy and runtime.

---

## Non-negotiable design decisions

These are the architectural rules for the cleanup.

### 0. Annotation status and biochemical state are orthogonal

Rigel needs to model two separate dichotomies:

1. annotated versus unannotated
2. mature RNA versus nascent RNA

Those axes must not be collapsed into one label.

That gives four biologically meaningful cases:

- annotated mature RNA
- unannotated mature RNA
- annotated nascent RNA
- unannotated nascent RNA

The cleanup in this document focuses on getting cases 1, 3, and 4 represented coherently. Case 2, especially unannotated mature RNA with unannotated splice junctions, should remain explicitly deferred for future work rather than being smuggled into the current data model.

Orthogonal biological labels do not imply separate EM states. Separate components should only exist when they represent genuinely distinct modeled intervals.

### 1. Core quantification is transcript-centric

Rigel should continue to quantify mature RNA at the transcript level.
Gene groupings are downstream conveniences only and must not shape the core model.

### 2. Keep one canonical transcript state space

Rigel should continue to use one unified transcript dataframe as the canonical EM state space.

Nascent RNA models are still transcript models. The important distinctions are:

- annotated versus Rigel-added
- distinct genomic span versus already-represented genomic span

The cleanup should preserve that single canonical state space rather than introducing a second nRNA entity space that duplicates single-exon annotated transcripts.

### 3. Annotated nascent-equivalent single-exon transcripts must collapse duplicate states

If an annotated single-exon transcript already covers the implied nascent span, Rigel must not create a second row or EM component for that same interval.

That single annotated transcript row is sufficient because the mature and nascent interpretations of that exact single-exon span are not identifiable from the data. Duplicating the row would only:

- enlarge the state space
- destabilize the EM by introducing a symmetric duplicate
- produce arbitrary attribution with no accuracy gain

### 4. Gene ownership must be removed from implied nascent transcript rows

If a user wants gene-oriented summaries, they must be derived after quantification. The core index and solver must not force an implied nascent transcript to belong to one gene just because one representative transcript happened to be seen first.

### 5. Calibration is mandatory

`quant_from_buffer()` should reflect reality: calibration is required before scoring and EM.

### 6. Scheduler behavior must be robust for both tiny and massive loci

The locus execution strategy should avoid both extremes:

- treating every small locus as a mega-locus
- starving threads behind one or two large loci

---

## Target end state

After cleanup, the core data model should look like this.

### Core reference tables

1. `t_df`
   One canonical transcript table used by the EM.

   It includes:

   - annotated transcripts from the input annotation
   - Rigel-added implied nascent transcripts when no equivalent annotated transcript already exists

   It must not contain duplicate rows for the same single-exon genomic span just to represent mature and nascent interpretations separately.

   Suggested metadata cleanup:

   - `is_annotated`
   - `is_synthetic_nrna`
   - `is_nascent_equiv`
   - optional explicit role flags only if they do not increase state cardinality

2. Optional helper metadata tables

   Helper tables may still be useful for reporting or debug, but they must not define a second canonical EM state space.

3. Gene convenience mapping

   Synthetic implied nascent transcripts should no longer inherit a single `g_id` or `g_name` from one representative transcript.

### Core solver layout

For a locus with `T` transcript rows in the canonical table:

- each transcript row is one component
- a single-exon annotated nascent-equivalent transcript still contributes only one component
- a synthetic implied nascent transcript contributes a separate component only when it represents a genuinely distinct genomic span
- one additional locus gDNA component is included when appropriate

The critical invariant is no duplicate mature/nascent components for the same single-exon interval.

### Output contract

1. `quant.feather` / transcript output
   Annotated transcript output should contain annotated transcripts only.

   Annotated single-exon nascent-equivalent transcripts remain here once, not twice.

2. Nascent reporting

   Any nascent-oriented report must be derived from the same canonical transcript rows, not from duplicate EM components.

   If a separate `nrna_quant` file is retained, it should be understood as a reporting view over canonical transcript rows and synthetic implied rows, not as evidence that annotated single-exon transcripts were modeled twice.

3. `gene_quant.feather`
   Mature transcript convenience aggregation only. Synthetic implied nascent rows must not inflate transcript counts or gene transcript counts via fake single-gene ownership.

4. Optional convenience link outputs
   If needed, export helper metadata tables rather than forcing nascent counts into a one-gene ownership model.

---

## Cleanup phases

## Phase 0 - Lock the intended behavior with tests

This phase should happen before any structural refactor.

### Add tests for transcript-centric semantics

1. Same-strand overlapping or interleaved annotations where one nRNA entity is supported by transcripts carrying different gene labels.
   Expected result: the implied nascent transcript exists without inheriting a single gene owner.

2. Annotated single-exon nascent-equivalent transcript.
   Expected result: no duplicate transcript/component is created, the annotated single-exon transcript remains the sole row for that span, and EM never has to split mass across indistinguishable mature/nascent duplicates.

3. Gene convenience outputs.
   Expected result: synthetic or explicit nRNA entities do not inflate `num_transcripts`, `n_transcripts`, or transcript-level mRNA outputs.

4. Direct `quant_from_buffer()` call without calibration.
   Expected result: immediate, explicit failure at the Python API boundary.

5. Scheduler edge cases.
   Expected result: small workloads are batched, mega-loci are isolated only when warranted, and `total_work < n_threads` does not force all loci down the mega path.

### Files likely touched

- `tests/test_index_integrity.py`
- `tests/test_pipeline_routing.py`
- `tests/test_estimator.py`
- `tests/test_pipeline_smoke.py`
- new scenario tests under `tests/scenarios/`

### Exit criteria

The new tests fail under current behavior for the right reasons and define the target semantics.

---

## Phase 1 - Clean up the unified transcript index contract

This is the most important cleanup phase.

### 1A. Keep `t_df` as the canonical EM state space

Do not replace the unified transcript dataframe with a separate nRNA entity table.

### 1B. Keep implied nascent transcripts only when they add a real state

Rigel-added implied nascent transcripts should remain transcript rows only when they create a distinct modeled span that is not already present as an annotated single-exon transcript.

### 1C. Handle annotated equivalents by collapse, not duplication

When a single-exon annotated transcript fully covers a nascent span:

- do not append a synthetic transcript row
- point any `nrna_t_index` or equivalent links at that existing transcript row
- optionally mark the transcript row as `is_nascent_equiv` or with equivalent metadata

This preserves the state while avoiding non-identifiable duplication.

### 1D. Remove single-gene inheritance

Synthetic implied nascent transcripts must not copy `g_id`, `g_name`, or `g_type` from one representative transcript.

Any gene-oriented relationship must be derived later from transcript membership or helper metadata, never stored as fake ownership on the implied transcript row.

### 1E. Make index format changes explicit

This is an index schema change. There is NO need for backward compatibility. We can move forward and disregard the prior behavior.

### Files likely touched

- `src/rigel/index.py`
- `src/rigel/transcript.py`
- index writing and loading code
- integrity tests and golden index expectations

### Exit criteria

The loaded index has one canonical transcript row per modeled state. Annotated-equivalent single-exon spans are collapsed to one row, and implied nascent transcripts no longer carry fake single-gene ownership.

---

## Phase 2 - Refactor scoring and routing around one canonical transcript state space

Once the index is clean, the hot path must stop inferring component type from transcript flags.

### 2A. Remove `is_synthetic_nrna` from scoring semantics

The native scorer should stop relying on flags whose meaning mixes origin, annotation status, and biochemical interpretation.

### 2B. Prevent duplicate mature/nascent candidates for identical single-exon spans

For annotated single-exon nascent-equivalent transcripts, there should be one candidate/component, not separate mature and nascent duplicates.

For distinct implied nascent transcript rows, the scorer should continue treating them as separate transcript candidates when they truly represent distinct genomic spans.

### 2C. Preserve output annotations by component type

The annotation layer should classify winners from transcript-row metadata that distinguishes annotated versus implied rows without inventing duplicate single-exon states.

### Files likely touched

- `src/rigel/scoring.py`
- `src/rigel/scan.py`
- `src/rigel/scored_fragments.py`
- `src/rigel/native/scoring.cpp`
- related native bindings

### Exit criteria

Candidate generation never duplicates identical single-exon spans, and transcript-row metadata is sufficient to interpret annotated versus implied models.

---

## Phase 3 - Stabilize the locus/EM contract around identifiable transcript states

The current EM path needs semantic cleanup, but the right target is not a separate transcript space plus nRNA entity space for identical single-exon regions.

### 3A. Update per-locus data structures

`Locus`, `ScoredFragments`, `LocusPartition`, and related structures should represent one component per canonical transcript row, plus gDNA where appropriate.

### 3B. Remove non-identifiable duplicate states

Any path that would create separate mature/nascent components for the same single-exon annotated span must be removed.

### 3C. Update native EM interfaces

The C++ EM solver should operate on a cleaner transcript-row contract where each component is an identifiable modeled state.

This includes:

- priors by transcript role where needed
- bias profiles for transcript rows versus gDNA
- annotation winner decoding

### 3D. Keep reporting semantics separate from EM state

If the code needs separate mature or nascent reporting views, they must be derived from transcript-row metadata after quantification rather than by introducing duplicate EM states.

### Files likely touched

- `src/rigel/locus.py`
- `src/rigel/scored_fragments.py`
- `src/rigel/partition.py`
- `src/rigel/estimator.py`
- `src/rigel/native/em_solver.cpp`

### Exit criteria

The solver operates on identifiable transcript states only, with no duplicate mature/nascent components for the same single-exon interval.

---

## Phase 4 - Rebuild reporting around one canonical transcript table

This phase removes the downstream bookkeeping hacks while preserving one canonical transcript table.

### 4A. Transcript output

`get_counts_df()` should report annotated transcript rows only.

Annotated single-exon nascent-equivalent transcripts remain here once because they are annotated transcript rows.

### 4B. Nascent reporting

If `nrna_quant` is retained, it should be a reporting view over canonical transcript rows and synthetic implied rows.

It must not imply that annotated single-exon transcripts were modeled as two separate EM components.

### 4C. Gene convenience output

`get_gene_counts_df()` should aggregate mature transcript outputs from annotated transcripts only.

It should not:

- count implied nascent transcript rows as annotated transcripts
- derive gene-level state from fake single-gene ownership on implied nascent rows
- assume nRNA belongs to a single gene

If a gene-oriented nRNA convenience export is still wanted later, it should be a separate derived table with explicit many-to-many semantics.

### 4D. Derived gene table cleanup

The base gene table should count annotated transcripts only. It should never include synthetic implied nascent rows.

### Files likely touched

- `src/rigel/estimator.py`
- `src/rigel/index.py`
- `src/rigel/cli.py`
- golden output tests and output docs

### Exit criteria

Transcript, nascent reporting, and gene outputs can each be explained directly from one canonical transcript table without duplicate single-exon EM states.

---

## Phase 5 - Make calibration mandatory in the public API

This is a smaller but important contract fix.

### 5A. Update function signatures

`quant_from_buffer()` should require a calibration object instead of defaulting to `None`.

### 5B. Fail fast and clearly

If an internal caller reaches quantification without calibration, raise an immediate `ValueError` or `TypeError` with a message that says calibration is mandatory.

### 5C. Align docs and tests

Remove stale optional wording from docstrings and docs.

### Files likely touched

- `src/rigel/pipeline.py`
- tests covering direct quantification entry points
- docs that describe pipeline stages

### Exit criteria

The public and internal API contracts match actual runtime requirements.

---

## Phase 6 - Replace the mega-locus heuristic with balanced work scheduling

The current `fair_share` rule is too brittle.

### 6A. Estimate locus cost explicitly

Use a work estimate tied to real solver cost, for example:

- `n_units * mean_candidates_per_unit`
- or `n_candidates`
- or equivalence-class work if that data is already available cheaply

The exact metric can be chosen after a quick profile comparison, but it must correlate with EM wall time better than `n_transcripts * n_units` alone.

### 6B. Batch with greedy LPT scheduling

Use a Largest Processing Time style packing strategy:

1. sort loci by descending estimated work
2. keep a small set of work bins or batches
3. assign each locus to the currently lightest batch

This behaves well for both small and skewed workloads and avoids the `fair_share = 0` failure mode.

### 6C. Preserve single-locus isolation for true outliers

If one locus exceeds a configurable fraction of total work, run it alone. That should be based on explicit work ratios, not an incidental integer-division threshold.

### 6D. Measure scheduler quality

Profile with:

- many tiny loci
- mixed real-world loci
- one or two mega-loci plus a long tail

The scheduler should reduce idle time without adding complicated execution modes.

### Files likely touched

- `src/rigel/pipeline.py`
- potentially `src/rigel/native/em_solver.cpp` only if deeper task-queue changes become necessary
- scheduler-focused tests

### Exit criteria

Small jobs stay batched, large jobs stay balanced, and scheduler behavior is explainable from the code.

---

## Phase 7 - Remove stale documentation and compatibility hacks

After the model and outputs are stable:

1. update `docs/METHODS.md`
2. update `docs/MANUAL.md`
3. fix stale docstrings in `pipeline.py`, `locus.py`, `estimator.py`, and `scored_fragments.py`
4. remove temporary compatibility branches that only existed to support the old synthetic-transcript design

This phase should be last. Documentation should describe the code that actually exists.

---

## Recommended execution order

This cleanup is large enough that it should be split into reviewable chunks.

### PR 1

Phase 0 tests plus Phase 5 calibration contract fix.

### PR 2

Phase 1 index redesign and index-format update.

### PR 3

Phase 2 scoring/routing refactor.

### PR 4

Phase 3 locus/EM refactor.

### PR 5

Phase 4 reporting cleanup and golden output refresh.

### PR 6

Phase 6 scheduler redesign and performance validation.

### PR 7

Phase 7 documentation cleanup and dead-code removal.

---

## Risks to manage explicitly

1. Index schema churn.
   The cleanup changes on-disk reference structures. Keep rebuild requirements explicit.

2. Golden output churn.
   Output files will change once transcript-row semantics are cleaned up. Expect intentional golden updates.

3. Performance regressions during the data-model transition.
   The right design should reduce duplicate non-identifiable states, but routing and EM extraction must still be kept efficient.

4. Partial migrations.
   The largest risk is leaving the code half-converted again. Each phase needs a crisp ownership model and clear exit criteria.

---

## Success criteria

The cleanup is complete when all of the following are true:

1. The canonical EM state space remains one transcript table.
2. Rigel-added implied nascent transcripts exist only when they add a distinct modeled state.
3. Annotated-equivalent single-exon transcripts never create duplicate mature/nascent EM components.
4. Synthetic implied nascent transcripts do not carry fake single-gene ownership.
5. Gene outputs are convenience aggregations, not core state.
6. `quant_from_buffer()` requires calibration.
7. The scheduler behaves sensibly across both small and pathological loci.
8. The docs and docstrings describe the implementation that actually runs.
