# Hulkrna Cleanup Backlog (2026-02-18)

## How to use this backlog
- Priority scale: P0 critical correctness, P1 high-value stabilization, P2 maintainability and developer velocity.
- Estimate scale: S (0.5-1 day), M (2-4 days), L (1-2 weeks).
- Definition of done for all tickets:
  - Code changes merged with tests.
  - Relevant docs updated.
  - No regression in targeted benchmark scenarios.

## Epic A: Correctness and Contract Safety (P0/P1)

### HULK-001 — Remove dead condition in EM shadow routing
- Priority: P0
- Estimate: S
- Files:
  - [src/hulkrna/pipeline.py](src/hulkrna/pipeline.py)
  - [tests/test_gdna.py](tests/test_gdna.py)
- Problem:
  - Dead always-true branch makes logic misleading and risks future bugs.
- Scope:
  - Replace the always-true condition in multimapper shadow generation with explicit intended logic.
  - Add focused tests for expected shadow creation behavior.
- Acceptance criteria:
  - No always-true or equivalent dead branch remains in this code path.
  - Tests explicitly validate shadow behavior for spliced_annot and non-spliced_annot multimappers.
  - Existing gDNA tests continue to pass.

### HULK-002 — Fix unique-count stats semantics
- Priority: P0
- Estimate: M
- Files:
  - [src/hulkrna/pipeline.py](src/hulkrna/pipeline.py)
  - [src/hulkrna/stats.py](src/hulkrna/stats.py)
  - [tests/test_stats.py](tests/test_stats.py)
- Problem:
  - n_counted_truly_unique is incremented in multiple routes, including non-deterministic paths.
- Scope:
  - Define explicit counters:
    - deterministic_unique_units
    - em_routed_unique_units
    - isoform_ambig_units
    - gene_ambig_units
    - multimapper_units
  - Update output schema and internal increment points.
- Acceptance criteria:
  - Each buffered unit increments exactly one route counter.
  - Counter definitions are documented in stats schema.
  - Unit tests verify route exclusivity and conservation.

### HULK-003 — Replace assert-based index integrity checks with runtime exceptions
- Priority: P1
- Estimate: S
- Files:
  - [src/hulkrna/index.py](src/hulkrna/index.py)
  - [tests/test_sim.py](tests/test_sim.py)
- Problem:
  - assert can be skipped in optimized runs.
- Scope:
  - Replace assert checks in load path with explicit validation and descriptive exceptions.
- Acceptance criteria:
  - No assert-based runtime data validation remains in index load path.
  - Failing index integrity emits actionable exception text.

### HULK-004 — Guard overlap and fragment-length dtype bounds
- Priority: P1
- Estimate: M
- Files:
  - [src/hulkrna/buffer.py](src/hulkrna/buffer.py)
  - [tests/test_buffer.py](tests/test_buffer.py)
- Problem:
  - int16 and uint16 may overflow for long fragments or high overlap values.
- Scope:
  - Audit and update dtypes for overlap_bp, intron_bp, frag_length.
  - Add boundary tests near type limits.
- Acceptance criteria:
  - No overflow in synthetic high-length stress tests.
  - Memory increase is measured and documented.

### HULK-005 — Add strict GTF parse mode
- Priority: P1
- Estimate: M
- Files:
  - [src/hulkrna/gtf.py](src/hulkrna/gtf.py)
  - [src/hulkrna/transcript.py](src/hulkrna/transcript.py)
  - [src/hulkrna/cli.py](src/hulkrna/cli.py)
  - [tests/test_sim.py](tests/test_sim.py)
- Problem:
  - Current parser warns and skips malformed lines, potentially hiding data corruption.
- Scope:
  - Introduce parse mode flag: strict or warn-skip.
  - Expose mode in CLI and default to strict for production.
- Acceptance criteria:
  - Strict mode fails fast with line-numbered error.
  - warn-skip mode preserves current behavior and reports skipped line counts.

## Epic B: Simplification and Refactoring (P1/P2)

### HULK-006 — Split pipeline into focused components
- Priority: P1
- Estimate: L
- Files:
  - [src/hulkrna/pipeline.py](src/hulkrna/pipeline.py)
  - [src/hulkrna/counter.py](src/hulkrna/counter.py)
- Problem:
  - One large orchestration file makes reasoning and changes risky.
- Scope:
  - Extract into modules:
    - candidate_builder
    - scoring
    - em_solver
    - assignment
  - Keep public run_pipeline API stable.
- Acceptance criteria:
  - Each component has its own unit tests.
  - run_pipeline behavior remains backward-compatible.
  - File-level complexity and function length are reduced.

### HULK-007 — Introduce typed EM candidate schema
- Priority: P1
- Estimate: L
- Files:
  - [src/hulkrna/pipeline.py](src/hulkrna/pipeline.py)
  - [src/hulkrna/counter.py](src/hulkrna/counter.py)
- Problem:
  - Parallel list assembly is error-prone and hard to validate.
- Scope:
  - Add dataclass or typed builder object for candidate batches before CSR conversion.
  - Enforce validation on lengths, offsets, and component ranges.
- Acceptance criteria:
  - Candidate construction invariants are validated in one place.
  - EMData creation has explicit schema checks.

### HULK-008 — Remove or fully implement gdna_intronic channel
- Priority: P1
- Estimate: M
- Files:
  - [src/hulkrna/counter.py](src/hulkrna/counter.py)
  - [src/hulkrna/pipeline.py](src/hulkrna/pipeline.py)
  - [tests/test_gdna.py](tests/test_gdna.py)
- Problem:
  - gdna_intronic fields are reported but not actively assigned in current route.
- Scope:
  - Choose one:
    - remove from schema and outputs, or
    - implement full assignment path and tests.
- Acceptance criteria:
  - No partially-implemented state remains exposed in outputs.
  - Decision and rationale documented in benchmark docs.

### HULK-009 — Prune stale stats fields and add schema version
- Priority: P1
- Estimate: M
- Files:
  - [src/hulkrna/stats.py](src/hulkrna/stats.py)
  - [src/hulkrna/cli.py](src/hulkrna/cli.py)
  - [tests/test_stats.py](tests/test_stats.py)
- Problem:
  - Legacy fields linger without active updates; output semantics are unclear.
- Scope:
  - Add stats_schema_version.
  - Deprecate or remove stale fields.
  - Add migration note for downstream users.
- Acceptance criteria:
  - stats.json includes schema version.
  - All fields are either populated or explicitly deprecated.

### HULK-010 — Clean strand model surface
- Priority: P2
- Estimate: M
- Files:
  - [src/hulkrna/strand_model.py](src/hulkrna/strand_model.py)
  - [src/hulkrna/pipeline.py](src/hulkrna/pipeline.py)
  - [tests/test_strand_model.py](tests/test_strand_model.py)
- Problem:
  - Multiple stored models but one active scoring path increases cognitive load.
- Scope:
  - Clarify which models are active vs diagnostic only.
  - Rename fields or split active and diagnostic models.
- Acceptance criteria:
  - API makes active scoring model obvious.
  - Diagnostic-only models are clearly tagged and tested.

## Epic C: Testing and Quality Gates (P1/P2)

### HULK-011 — Add CLI contract tests
- Priority: P1
- Estimate: M
- Files:
  - [src/hulkrna/cli.py](src/hulkrna/cli.py)
  - [tests/test_cli.py](tests/test_cli.py)
- Problem:
  - No direct tests for CLI argument behavior and output contract.
- Scope:
  - Add parser and command tests for key flags and defaults.
  - Validate generated config and stats schema structure.
- Acceptance criteria:
  - Core count/index command options are covered.
  - Invalid parameter combinations fail with clear messages.

### HULK-012 — Add targeted regression set from diagnostic failures
- Priority: P1
- Estimate: M
- Files:
  - [tests/test_diagnostic.py](tests/test_diagnostic.py)
  - [tests/test_gdna.py](tests/test_gdna.py)
  - [tests/test_counter.py](tests/test_counter.py)
- Problem:
  - Known failure modes are diagnosed but not all are codified as stable regression gates.
- Scope:
  - Convert selected diagnostic scenarios into deterministic regression tests with thresholds.
- Acceptance criteria:
  - At least 3 historically unstable regimes become gating tests.
  - Thresholds justified by benchmark evidence.

### HULK-013 — Add overflow and conservation invariants suite
- Priority: P2
- Estimate: M
- Files:
  - [tests/test_buffer.py](tests/test_buffer.py)
  - [tests/test_counter.py](tests/test_counter.py)
  - [tests/test_pipeline.py](tests/test_pipeline.py)
- Scope:
  - Add invariants:
    - Unit conservation per EM route.
    - Pool conservation mRNA+nRNA+gDNA per unit.
    - Buffer numeric bounds.
- Acceptance criteria:
  - Invariant failures produce clear diagnostics.
  - Tests run in CI-friendly time.

## Epic D: Documentation and Developer UX (P2)

### HULK-014 — Publish architecture and decision record for new model
- Priority: P2
- Estimate: M
- Files:
  - [README.md](README.md)
  - [docs](docs)
- Problem:
  - The main architecture and model rationale are not discoverable from README.
- Scope:
  - Add concise architecture doc for mRNA/nRNA/gDNA model.
  - Add glossary and route definitions.
- Acceptance criteria:
  - New contributor can trace fragment flow from BAM to final counts.
  - All count outputs are defined in one reference page.

### HULK-015 — Add benchmark regression playbook
- Priority: P2
- Estimate: S
- Files:
  - [docs/regional_benchmark.md](docs/regional_benchmark.md)
  - [docs/benchmark_hulkrna_vs_salmon.md](docs/benchmark_hulkrna_vs_salmon.md)
- Scope:
  - Standardize benchmark commands, expected artifacts, and pass criteria.
- Acceptance criteria:
  - One reproducible command set for pre-merge performance and accuracy checks.
  - Explicit fail conditions for contamination and isoform bias metrics.

## Suggested execution order
1. HULK-001, HULK-002, HULK-003
2. HULK-004, HULK-005
3. HULK-009, HULK-011, HULK-012
4. HULK-006, HULK-007, HULK-008
5. HULK-010, HULK-013, HULK-014, HULK-015

## Milestone proposal
- Milestone 1 (Stabilize): HULK-001 to HULK-005
- Milestone 2 (Make it simpler): HULK-006 to HULK-010
- Milestone 3 (Harden and document): HULK-011 to HULK-015
