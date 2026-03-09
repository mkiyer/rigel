# Rigel Production Implementation Plan

Date: 2026-03-08

Source review:
- `docs/production_code_review_2026-03-08.md`

## Goal

Turn the production review into an execution plan that can be implemented in small, low-risk stages. The plan is designed to improve production safety first, then simplify boundaries, then reduce code size and long-term maintenance cost.

This plan assumes the scientific model remains intact. The focus is engineering quality: reliability, packaging, lifecycle safety, maintainability, and code clarity.

## Execution principles

1. Fix correctness and release-risk issues before beautification work.
2. Keep each stage small enough to validate with targeted tests.
3. Do not mix architectural refactors with algorithm changes unless required.
4. Preserve output semantics unless a stage explicitly redefines them.
5. Add tests before or alongside every boundary change.

## Current status update

The plan below has been updated after the first implementation pass.

Already landed:

- `_bam_impl` is now treated as required at build time.
- Publish validation imports the native extensions instead of only `import rigel`.
- Multimapper scoring has been moved into the fused C++ scan path, and the Python multimapper hot path has been removed from production scanning.

Partially landed:

- The scan path no longer does `list(buffer.iter_chunks())`, but it still accumulates all converted chunk arrays before calling the fused C++ scorer. This means Stage 2 is not complete.

Still unimplemented:

- A true end-to-end spill-safe scan interface.
- Cross-chunk multimapper regression coverage.
- A wheel smoke test that runs a minimal quant path.
- Output semantics, lifecycle cleanup, resolver boundary cleanup, and decode centralization.

Validation snapshot for the current implementation pass:

- `pytest tests/test_pipeline_routing.py tests/test_buffer.py -q`
- Result: 52 tests passed in the `rigel` conda environment.

## Stage 0: Baseline and guardrails

Status: partially complete

### Objective

Establish a stable baseline so the cleanup work does not regress correctness, determinism, or packaging.

### Scope

- Freeze a representative smoke-test dataset for end-to-end validation.
- Define the canonical verification matrix for each stage.
- Add one reusable smoke-test entry point for CI and release workflows.

### Deliverables

- A tiny end-to-end fixture or scenario that exercises:
  - index build/load
  - BAM scan
  - EM quantification
  - optional annotated BAM path if feasible
- A small script or test target used by both CI and publish workflows.
- A short contributor note describing required verification before merge.

### Validation

- Existing targeted tests pass.
- New smoke test passes locally and in CI.
- Smoke test is runnable from built wheels, not just editable installs.

### Status notes

- Targeted local validation exists for the changed scan/buffer paths.
- The canonical shared smoke test for CI and publish is still missing.
- Publish remains import-level validation only; it does not yet execute a tiny quant path.

### Exit criteria

- There is a single canonical smoke test used for release validation.
- The project can validate more than `import rigel` before publishing.

## Stage 1: Packaging and runtime contract for native modules

Status: largely complete

### Objective

Resolve the `_bam_impl` contract so builds, imports, and runtime behavior are consistent and predictable.

### Problem being solved

`_bam_impl` is optional at build time in `CMakeLists.txt`, but the Python pipeline imports it unconditionally. That creates a broken release contract.

### Decision required

Choose one of two supported models:

#### Option A: `_bam_impl` is required

- Build fails if `htslib` is unavailable.
- Documentation and installation guidance state that native BAM support is mandatory.
- Release artifacts are guaranteed to include the extension.

#### Option B: `_bam_impl` is optional

- Python imports become lazy or guarded.
- Runtime raises a clear capability error when BAM-backed features are invoked without native support.
- A fallback implementation exists, or functionality is explicitly unavailable.

### Recommended direction

Option A unless there is a strong distribution reason to support degraded installs. The current codebase is built around native BAM scanning, so pretending it is optional adds complexity without much product value.

### Scope

- Make build-time behavior explicit in `CMakeLists.txt`.
- Align Python imports with the chosen contract.
- Update release smoke tests to import quant-related modules.
- Update installation documentation.

### Deliverables

- Explicit native dependency policy.
- Consistent build/import/runtime behavior.
- Publish workflow smoke tests that import at least:
  - `rigel.pipeline`
  - `rigel.annotate`
  - native-backed modules used by quantification

### Status notes

- The required-native policy is now implemented in `CMakeLists.txt`.
- The publish workflow now imports multiple native extensions directly.
- Remaining work is documentation cleanup and upgrading publish smoke from import-only to minimal end-to-end quant validation.

### Validation

- Build succeeds or fails for the right reason.
- Wheel smoke tests exercise actual quant-related imports.
- Documentation matches real behavior.

### Exit criteria

- No release artifact can pass validation while shipping a broken quant import path.

Practical interpretation after current changes:

- The broken import-path problem is mostly resolved.
- The remaining gap is quant-path execution coverage on built artifacts.

## Stage 2: Fix spill-path memory behavior

Status: partially complete

### Objective

Make `FragmentBuffer` spilling actually protect memory during quantification.

### Problem being solved

`FragmentRouter.scan()` re-materializes all chunks with `list(buffer.iter_chunks())`, which defeats on-demand spill loading.

Updated status:

- The `list(buffer.iter_chunks())` materialization is gone.
- However, `_scan_native()` still appends every converted chunk into `chunk_arrays` before invoking `fused_score_buffer()`, so all converted chunk data still accumulates in memory.
- This stage should remain open until the fused scan path can operate without retaining the full converted input at once.

### Scope

- Refactor the chunk-consumption path to stream instead of eagerly materialize.
- If the fused native scorer truly requires all chunks together, redesign the interface so the requirement is explicit and spill-aware.
- Add regression coverage around forced spilling.

### Preferred implementation order

1. Document the intended memory contract in `buffer.py` and `scan.py`.
2. Measure current peak memory with forced spilling.
3. Refactor the fused scorer boundary so chunk conversion and native consumption can proceed incrementally.
4. Re-measure memory behavior.

### Deliverables

- Streaming-safe scan path.
- Forced-spill test case.
- Clear invariant: spill loading must not rehydrate all chunks at once.

Additional deliverable now required:

- Remove the all-chunk `chunk_arrays` accumulation from `scan.py` or make that buffering requirement explicit and justified in the design.

### Validation

- Buffer tests continue to pass.
- New forced-spill regression test passes.
- Peak memory is materially lower or at least bounded as intended.

### Exit criteria

- The chunk iterator contract is true in practice, not just in docstrings.

Not sufficient for closure:

- Merely removing `list(buffer.iter_chunks())` is not enough if the converted arrays are still accumulated before scoring.

## Stage 2A: Regression hardening for the new fused scan path

### Objective

Lock down the correctness of the new C++ multimapper path and the revised scan contract before more refactors land.

### Problem being solved

The multimapper path is now in C++, which is the right architecture, but the current focused tests do not cover some of the highest-risk boundaries.

### Scope

- Add a cross-chunk multimapper regression where the same `frag_id` spans a chunk boundary.
- Add a mixed-splice multimapper regression across that boundary.
- Add a forced-spill scan test that validates the intended memory contract at the scan layer.

### Deliverables

- Dedicated multimapper chunk-boundary tests.
- Dedicated forced-spill scan-path regression test.
- Clear documentation of what the fused scan path promises about chunk retention.

### Validation

- `tests/test_pipeline_routing.py` covers cross-chunk multimapper behavior.
- Spill-focused tests cover scan behavior, not only buffer round-tripping.

### Exit criteria

- The new fused scan path has explicit regression coverage for its most failure-prone boundaries.

## Stage 3: Clean up lifecycle ownership and native boundary decoding

### Objective

Reduce fragility around cleanup and Python/C++ transport boundaries.

### Problems being solved

- `FragmentBuffer` cleanup relies on `__del__`.
- Native raw-byte decoding into `_FinalizedChunk` is duplicated.
- Integer narrowing is implicit rather than enforced.

### Scope

- Remove or demote destructor-based cleanup.
- Centralize native raw-bytes to `_FinalizedChunk` decoding.
- Audit and fix unsafe narrowing conversions.

### Deliverables

- A single decode entry point, for example:
  - `_FinalizedChunk.from_native_raw(...)`
  - or a dedicated decoder helper in `buffer.py`
- Explicit cleanup ownership model using context managers and/or `weakref.finalize`.
- Checked conversions or preserved `int64` where needed.

### Validation

- Buffer tests pass.
- Spill cleanup tests pass under normal and exception paths.
- Large-ID assumptions are explicit and tested.

### Exit criteria

- There is one authoritative decode path.
- Cleanup behavior does not depend on CPython destructor timing.

## Stage 4: Make resolver ownership explicit

### Objective

Stop leaking private resolver state across the codebase.

### Problem being solved

Multiple modules rely on `TranscriptIndex._resolver`, which is private implementation detail used as a shared global capability.

### Scope

- Introduce a public API boundary for resolution-related operations.
- Move call sites away from direct `_resolver` access.

### Possible designs

#### Option A: Public property

- `index.native_resolver`
- Minimal change, but still exposes the low-level object.

#### Option B: Public methods on `TranscriptIndex`

- `index.resolve_fragment(...)`
- `index.set_gene_strands(...)`
- `index.set_transcript_strands(...)`
- Better encapsulation.

#### Option C: Dedicated service wrapper

- `ResolutionContext` or similar object owned by the index.
- Best separation if more refactoring is planned.

### Recommended direction

Option B now, Option C only if the resolution subsystem is going to grow significantly.

### Deliverables

- No production modules reach into `index._resolver` directly.
- The resolver lifecycle stays owned by the index/runtime index layer.

### Validation

- Resolution and pipeline tests pass.
- No direct `_resolver` references remain outside the index internals.

### Exit criteria

- The runtime index exposes a clean public boundary for native resolution capabilities.

## Stage 5: Resolve unfinished output semantics

### Objective

Eliminate placeholder logic from production outputs and define stable semantics.

### Problem being solved

Gene-level `locus_id` attribution is unfinished when a gene spans multiple transcript loci.

### Scope

- Decide the semantics for gene-level locus representation.
- Implement the policy.
- Add tests covering ambiguous multi-locus cases.

### Decision options

#### Option A: Dominant transcript policy

- Gene gets the `locus_id` of the transcript contributing the most mRNA.

#### Option B: Multi-locus sentinel

- Gene gets `-1` or a dedicated sentinel when transcripts span multiple loci.

#### Option C: Stop forcing a single scalar

- Keep gene table simple and emit a separate gene-to-loci relation table.

### Recommended direction

Option C is the cleanest model if downstream users need truth. Option A is the simplest if a single scalar is required for compatibility.

### Deliverables

- Defined semantics in code and docs.
- Tests for multi-locus genes.
- No `pass` placeholder in the output path.

### Validation

- Output tests pass.
- Golden-output or fixture-based tests reflect the chosen semantics.

### Exit criteria

- Every published output field has explicit, tested semantics.

## Stage 6: Split the orchestration layer

### Objective

Reduce cognitive load by turning `pipeline.py` into a true composition root.

### Problem being solved

`pipeline.py` is doing orchestration, transport decoding, model setup, annotation handling, and output shaping in one large module.

### Scope

Refactor without changing scientific behavior.

### Suggested module split

- `pipeline_scan.py`
  - BAM scan orchestration
  - model replay
  - buffer construction
- `pipeline_quant.py`
  - transcript geometry
  - estimator construction
  - EM handoff
- `pipeline_annotation.py`
  - annotation table setup
  - second-pass annotated BAM logic
- `pipeline_native.py`
  - native bridge helpers and decode shims

### Deliverables

- Smaller orchestration modules.
- `run_pipeline()` remains short and readable.
- Internal helpers grouped by responsibility.

### Validation

- Pipeline and scenario tests pass.
- Public API of `run_pipeline()` remains stable unless intentionally changed.

### Exit criteria

- `pipeline.py` reads like composition, not implementation detail accumulation.

## Stage 7: Split the estimator/output subsystem

### Objective

Break `estimator.py` into smaller units that match the real architecture.

### Problem being solved

`estimator.py` currently mixes:

- EM data structures
- prior estimation helpers
- batch EM integration
- output-frame construction
- posterior/confidence metrics

### Suggested split

- `estimator_core.py`
  - `AbundanceEstimator`
  - core state and accumulation
- `estimator_priors.py`
  - `estimate_kappa`
  - `compute_hybrid_nrna_frac_priors`
  - related prior logic
- `estimator_em.py`
  - batch EM bridge
  - locus EM integration helpers
- `estimator_outputs.py`
  - transcript, gene, locus, and detail DataFrame builders

### Deliverables

- Reduced file size and clearer API boundaries.
- Output generation separated from inference state transitions.

### Validation

- Estimator tests pass.
- Output schemas remain stable or are intentionally versioned.

### Exit criteria

- Inference code and reporting code are no longer interleaved in one giant module.

## Stage 8: Split build-time and runtime index responsibilities

### Objective

Untangle index construction from runtime query/index loading behavior.

### Problem being solved

`index.py` handles both index generation and runtime loading/query setup. Those are different responsibilities with different maintenance pressures.

### Suggested split

- `index_build.py`
  - reference length loading
  - GTF parsing and transcript extraction
  - interval and splice-junction table generation
  - on-disk writing
- `index_runtime.py`
  - `TranscriptIndex`
  - load-time validation
  - runtime lookup structures
  - native resolver setup

### Deliverables

- Cleaner build/runtime separation.
- Easier future support for alternate index formats or versioning.

### Validation

- Index integrity tests pass.
- Build/load round-trip tests pass.

### Exit criteria

- Runtime index loading can evolve without carrying build-pipeline clutter.

## Stage 9: Native code modularization

### Objective

Reduce the risk of large-scale native files becoming unreviewable and hard to modify safely.

### Problem being solved

Several native files are very large, especially `em_solver.cpp`, `bam_scanner.cpp`, `scoring.cpp`, and `resolve_context.h`.

### Scope

Do not change algorithms in this stage. Only split by responsibility.

### Suggested split areas

#### `em_solver.cpp`

- equivalence-class construction
- bias correction
- E-step/M-step kernels
- SQUAREM / convergence logic
- nanobind module interface

#### `bam_scanner.cpp`

- BAM iteration and htslib interaction
- fragment grouping
- stats accumulation
- model observation collection
- accumulator export

#### `scoring.cpp`

- fragment weighting and mapping helpers
- scorer object definition
- fused buffer scoring
- binding layer

#### `resolve_context.h`

- data types/constants
- overlap index building
- splice-junction helpers
- core fragment resolution logic

### Deliverables

- Smaller native translation units and headers.
- Better locality for profiling and code review.

### Validation

- Native extension tests pass.
- Performance benchmark does not regress materially.
- Profiling symbols remain usable.

### Exit criteria

- No critical native subsystem depends on one monolithic source file for maintainability.

## Stage 10: CI, linting, typing, and contributor ergonomics

### Objective

Raise the baseline engineering bar for future changes.

### Scope

- Run `ruff` in CI.
- Add type checking for Python where practical.
- Add contributor guidance for Python/native boundary work.
- Consider a small style guide for new modules.

### Deliverables

- Lint gate in CI.
- Initial type-check coverage or scoped enforcement.
- `docs/contributing_engineering.md` or similar guidance.

### Validation

- CI fails on lint regressions.
- Type checking is actionable and not pure noise.

### Exit criteria

- New cleanup debt is harder to introduce accidentally.

## Recommended PR sequence

Use this as the implementation order.

1. Finish Stage 0 with one canonical smoke test shared by CI and publish.
2. Close the remaining work in Stage 1 by upgrading publish smoke from import-only to minimal quant execution.
3. Finish Stage 2 by removing all-chunk accumulation from the fused scan entry point.
4. Complete Stage 2A with cross-chunk multimapper and forced-spill regressions.
5. Stage 3: cleanup ownership, decode centralization, integer-width audit.
6. Stage 4: public resolver boundary.
7. Stage 5: output semantics for gene-level locus attribution.
8. Stage 6: split `pipeline.py`.
9. Stage 7: split `estimator.py`.
10. Stage 8: split `index.py`.
11. Stage 9: native modularization.
12. Stage 10: CI/lint/type/contributor hardening.

## Suggested milestone structure

### Milestone A: Production safety

- Stages 0 through 3, including Stage 2A

Outcome:
- shipping and runtime behavior are reliable
- memory behavior is defensible
- boundary/lifecycle bugs are reduced

### Milestone B: Public architecture cleanup

- Stages 4 through 8

Outcome:
- module boundaries become explicit
- large Python files become easier to maintain
- output semantics are fully defined

### Milestone C: Long-term maintainability

- Stages 9 through 10

Outcome:
- native code becomes more reviewable
- future contributors have better guardrails

## Do not combine in one PR

Avoid combining these categories in a single change:

- packaging/runtime policy changes with module-splitting refactors
- output-semantic changes with native performance work
- boundary cleanup with unrelated algorithm adjustments
- large native file splits with feature additions

## Final recommendation

Given the current implementation state, the next two things to do should be:

1. Finish the spill-path fix by removing all-chunk converted-array accumulation before fused scoring.
2. Add regression coverage for cross-chunk multimappers and one minimal wheel-level quant smoke test.

The packaging/runtime contract has already improved substantially. The highest remaining production risk is the gap between the intended streaming/spill behavior and what the fused scan path actually retains in memory. Once that gap and the missing regression coverage are closed, the codebase is in a much better position for the structural cleanup work that makes it beautiful, understandable, and safe to extend.