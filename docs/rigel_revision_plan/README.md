# Rigel Revision Plan

This folder expands the high-level redesign in
`docs/FIRST_PRINCIPLES_IMPLEMENTATION_PLAN.md` into concrete workstream-level
plans tied to the current codebase.

The purpose is to separate:

- theory and model definition
- implementation sequencing
- per-workstream design notes

The intended workflow is:

1. agree on the workstream-level implementation target
2. execute one workstream at a time in small, testable steps
3. revise downstream workstreams as the upstream interfaces stabilize

## Document Map

- `00_overview.md`
  - repository-level roadmap and dependency graph among workstreams
- `01_workstream_A_region_partition.md`
  - annotation-context partition and region-table design
- `02_workstream_BC_region_evidence.md`
  - fragment-to-region assignment and regional evidence extraction
- `03_workstream_E_calibration.md`
  - empirical-Bayes calibration of gDNA nuisance parameters
- `04_workstream_F_em_integration.md`
  - later integration of calibrated nuisance parameters into the locus EM
- `05_purity_model_options.md`
  - practical purity approximations and recommended first implementation path
- `06_methodology_review.md`
  - assessment of the March 2026 methodology critique and resulting revisions

## Current Recommendation

The most definitive implementation path is:

1. Workstream A: expose a calibration-ready region partition from the existing
   index interval machinery
2. Workstream B/C: accumulate a standalone region evidence table from fragment
   objects already available during scan
3. Workstream E: calibrate gDNA nuisance parameters from weighted regional
   evidence
4. Workstream F: migrate the locus solver to the chosen internal `T + N + 2`
  gDNA-pair model while keeping collapsed public gDNA output

Workstream D, the purity model, should remain flexible for now. The revised
recommendation is to bootstrap calibration from a highly conservative seed set
of gDNA-dominant regions first, then optionally expand to soft weighting once
the empirical region-evidence distributions are visible.

## Near-Term Coding Order

1. add a calibration-region schema and export path on top of `TranscriptIndex`
2. add a calibration evidence accumulator independent of the main EM
3. emit diagnostic tables or summaries for region evidence
4. implement a first weighted calibration pass for gDNA symmetry and fragment
   length
5. refactor the locus EM only after upstream calibration outputs are stable