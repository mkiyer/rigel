# hulkrna Code Review and Cleanup Plan

Date: 2026-03-06

## Executive Summary

`hulkrna` already contains strong scientific ideas:

- transcript abundance estimation from RNA-seq,
- explicit deconvolution of genomic DNA contamination,
- separation of nascent RNA from mature RNA,
- strand specificity modeled as a continuum rather than a binary flag.

Those innovations are worth protecting. The cleanup plan below is therefore **not** a rewrite. It is a staged refactor focused on:

1. reducing duplicated and redundant code paths,
2. removing hard-coded parameters and magic numbers,
3. establishing one clear, well-documented home for user-facing settings,
4. promoting appropriate settings into advanced CLI/YAML options,
5. renaming vague symbols to make the code readable,
6. improving maintainability without materially regressing performance.

The main theme is:

> preserve the scientific model and the current performance shape, while making the codebase understandable, configurable, and production-ready.

---

## What I Observed

This plan is grounded in the current structure of the repository, especially these files:

- `src/hulkrna/cli.py`
- `src/hulkrna/config.py`
- `src/hulkrna/pipeline.py`
- `src/hulkrna/estimator.py`
- `src/hulkrna/locus.py`
- `src/hulkrna/scan.py`
- `src/hulkrna/scoring.py`
- `src/hulkrna/strand_model.py`
- `src/hulkrna/buffer.py`
- `src/hulkrna/annotate.py`
- `tests/test_cli.py`
- `docs/TODO.md`

### 1. Configuration is only partially centralized

`src/hulkrna/config.py` says it is the “single source of truth”, but in practice defaults and semantics are still split across multiple places:

- `src/hulkrna/config.py` dataclass defaults
- `src/hulkrna/cli.py` `_DEFAULTS` in `quant_command()`
- `src/hulkrna/cli.py` argparse help text
- `src/hulkrna/tests/test_cli.py` `_SMALL_DEFAULTS`
- duplicated fallback defaults in `src/hulkrna/estimator.py`, `src/hulkrna/locus.py`, `src/hulkrna/index.py`, and `src/hulkrna/pyfallback.py`

This is the single biggest maintainability problem because it guarantees drift.

### 2. The top-level quant command is doing too much

`quant_command()` currently:

- resolves CLI and YAML values,
- generates defaults,
- creates output paths,
- loads the index,
- normalizes special cases like `sj_strand_tag`,
- writes a config snapshot,
- builds `PipelineConfig`,
- invokes the pipeline,
- writes output artifacts.

That orchestration is functional, but it is crowded and hard to reason about.

### 3. Scientific defaults and engineering constants are mixed together

There are several categories of numbers in the codebase, but they are not clearly separated:

#### User-tunable scientific/model parameters
Examples:

- `em_prior_alpha = 0.01`
- `em_prior_gamma = 1.0`
- `prune_threshold = 0.1`
- `confidence_threshold = 0.95`
- `tss_window = 200`
- `nrna_frac_*`
- `gdna_kappa_*`
- `strand_prior_kappa = 2.0`
- scoring penalties like overhang and mismatch alpha

#### Internal numerical safeguards
Examples:

- `LOG_SAFE_FLOOR = 1e-10` in `src/hulkrna/scoring.py`
- `EM_PRIOR_EPSILON` from the native EM layer
- `_NEG_INF`
- low-level clamp and fallback values

#### Operational/performance constants
Examples:

- `_ANNOTATION_TABLE_PADDING = 1024` in `src/hulkrna/pipeline.py`
- `_ANNOTATION_TABLE_MIN_CAPACITY = 4096` in `src/hulkrna/pipeline.py`
- annotation table growth floor `1024` in `src/hulkrna/annotate.py`
- `chunk_size = 1_000_000`
- `max_memory_bytes = 2 * 1024**3`

These categories should not live in the same conceptual bucket.

### 4. Some modules are cleanly designed, but the design is not consistently applied

Good examples already exist:

- `src/hulkrna/config.py` uses typed dataclasses.
- `src/hulkrna/scoring.py` explicitly moves closure state into `FragmentScorer`.
- `src/hulkrna/pipeline.py` is trying to act as an orchestrator.

But other parts still retain legacy shape:

- duplicated default values,
- helper functions with overlapping responsibility,
- long functions with mixed orchestration and domain logic,
- fallback/native code paths that appear conceptually aligned but not always structurally unified.

### 5. Naming quality is uneven

At the module/API level, naming is often good: `StrandModel`, `FragmentScorer`, `PipelineConfig`, `AbundanceEstimator`.

At the implementation level, readability drops quickly because of short, overloaded names such as:

- `ctx`, `bf`, `fc`, `ll`, `ct`, `cw`, `rl`, `ebp`, `ibp`, `tx_s`, `tx_e`, `t_arr`, `u_arr`, `gt`, `lid`, `lbl`
- `mm_fid`, `gdna_lls`, `v_ll`, `v_cc`, `v_cw`

These are understandable in hot loops, but they make nontrivial scientific code much harder to audit.

### 6. Tests exist, but the test structure should become a refactor tool

The project already has meaningful test coverage, especially for CLI/config behavior and many scientific behaviors. That is a major asset.

However, the cleanup work will be much safer if tests are reorganized around the new architecture:

- parameter registry behavior,
- CLI/YAML/config round-tripping,
- backend parity,
- output schema contracts,
- benchmark guardrails for speed and memory.

### 7. Documentation of parameters is not yet user-grade

A power user can infer a lot from docstrings and CLI help, but there is not yet one obvious document answering:

- what each parameter means biologically,
- what unit it uses,
- what default means,
- when a user should touch it,
- whether it is safe/common vs advanced/expert-only.

That is critical for a tool with as much modeling depth as `hulkrna`.

---

## Refactor Principles

These principles should govern the whole cleanup effort.

### Principle A — Do not rewrite the science when the real problem is structure

The gDNA/nRNA/mRNA model is the product. Cleanup should mostly change code organization, naming, configuration plumbing, and documentation.

### Principle B — One source of truth per concern

- one source of truth for user-facing parameters,
- one source of truth for scientific defaults,
- one source of truth for low-level implementation constants,
- one source of truth for CLI/YAML/config serialization.

### Principle C — Separate user-facing knobs from internal constants

Not every number should become a CLI option. Some values are:

- model defaults worth exposing,
- advanced tuning knobs worth exposing carefully,
- internal implementation details that should remain private.

### Principle D — Refactor under test and benchmark protection

Every phase should preserve:

- numerical parity where intended,
- similar runtime,
- similar memory use,
- stable output schema unless intentionally versioned.

### Principle E — Prefer thin orchestration and pure domain modules

The healthiest architecture here is:

- CLI resolves input,
- config layer validates and normalizes,
- pipeline orchestrates,
- domain modules compute,
- reporting/output modules serialize.

---

## Target End State

By the end of this cleanup program, the codebase should have:

1. a **single parameter registry** that drives CLI, YAML, defaults, validation, saved config, and docs,
2. a **small top-level `quant` entrypoint** that mostly wires together components,
3. a **clear separation** between public model parameters and private implementation constants,
4. a **documented advanced-options surface** for expert tuning,
5. **fewer duplicated code paths** between Python/native/fallback layers,
6. **clear naming and module boundaries**,
7. **benchmark-backed confidence** that cleanup did not damage performance.

---

## Multi-Phase Implementation Plan

## Phase 0 — Establish Guardrails Before Moving Code

### Goal
Create a safe envelope for refactoring.

### Deliverables

- A pinned “refactor baseline” test command for fast local validation.
- A small benchmark suite for representative workloads:
  - small synthetic scenario,
  - moderate locus-heavy scenario,
  - one real-data smoke benchmark if available.
- Baseline metrics recorded in docs:
  - wall clock time,
  - peak memory if measurable,
  - key output parity summaries.
- A short architectural glossary document for core concepts:
  - mRNA,
  - nRNA,
  - gDNA,
  - locus,
  - TSS group,
  - strand specificity,
  - unambig vs multimapping.

### Why first
A spaghetti codebase becomes much riskier when touched without baseline measurements.

### Likely files

- `tests/`
- `scripts/benchmark.py`
- `docs/`

### Acceptance criteria

- Fast smoke suite exists.
- Benchmark entrypoints are documented.
- At least one parity fixture is available for before/after comparison.

---

## Phase 1 — Create a Real Single Source of Truth for Parameters

### Goal
Replace scattered defaults and help text with a unified parameter registry.

### Core problem being solved
Right now the same settings live in multiple places. That makes drift inevitable.

### Proposed design
Introduce a parameter schema module, for example:

- `src/hulkrna/parameter_registry.py`

Each parameter record should include:

- internal name,
- CLI flag name,
- YAML key,
- default value,
- type,
- category,
- whether it is basic or advanced,
- valid range/choices,
- short help text,
- long explanation,
- biological/algorithmic meaning,
- whether it is user-facing or internal.

Example categories:

- `io`
- `em`
- `strand_model`
- `fragment_scoring`
- `nrna_prior`
- `gdna_prior`
- `memory_and_scaling`
- `annotation_output`

### Concrete refactor steps

1. Move all user-facing defaults out of `src/hulkrna/cli.py` `_DEFAULTS`.
2. Make `src/hulkrna/config.py` consume the registry rather than hardcoding duplicated values.
3. Make CLI argument construction derive help/default metadata from the registry.
4. Make YAML resolution and saved `config.json` derive key names from the same registry.
5. Replace `tests/test_cli.py` `_SMALL_DEFAULTS` with tests that assert against the registry.
6. Identify duplicated default values still living in:
   - `src/hulkrna/estimator.py`
   - `src/hulkrna/locus.py`
   - `src/hulkrna/index.py`
   - `src/hulkrna/pyfallback.py`
   and either remove them or explicitly mark them as internal-only.

### Important distinction
Not every dataclass field in `PipelineConfig` needs to be a user-facing registry entry. Derived/internal values can stay internal.

### Likely files

- `src/hulkrna/cli.py`
- `src/hulkrna/config.py`
- new `src/hulkrna/parameter_registry.py`
- `tests/test_cli.py`
- possibly `README.md` and generated docs

### Acceptance criteria

- There is exactly one authoritative default for every user-facing quant parameter.
- CLI help, YAML loading, config snapshot output, and tests all agree.
- Adding a new parameter requires editing one place, not four.

---

## Phase 2 — Split `quant_command()` Into Clear Layers

### Goal
Turn `quant_command()` into a thin entrypoint.

### Current issue
`quant_command()` in `src/hulkrna/cli.py` is functioning as parser, config merger, validator, serializer, and orchestrator.

### Proposed target structure

Create dedicated helpers/modules such as:

- `resolve_quant_config()`
- `build_output_paths()`
- `write_run_manifest()`
- `run_quant_workflow()`
- `write_quant_outputs()`

or equivalent modules like:

- `src/hulkrna/cli_quant.py`
- `src/hulkrna/output_writer.py`
- `src/hulkrna/run_manifest.py`

### Concrete refactor steps

1. Extract argument/YAML normalization from `quant_command()`.
2. Extract output path creation and manifest writing.
3. Extract `PipelineConfig` construction into a dedicated conversion function.
4. Keep `quant_command()` focused on:
   - validate paths,
   - resolve config,
   - load index,
   - run workflow,
   - write outputs,
   - return exit code.

### Benefits

- easier unit testing,
- easier future CLI extension,
- much better readability,
- less risk of subtle config drift.

### Likely files

- `src/hulkrna/cli.py`
- new helper modules in `src/hulkrna/`
- `tests/test_cli.py`
- `tests/test_pipeline_routing.py`

### Acceptance criteria

- `quant_command()` becomes short and readable.
- config resolution can be unit-tested without running the pipeline.
- output manifest writing can be unit-tested independently.

---

## Phase 3 — Separate Scientific Parameters From Internal Implementation Constants

### Goal
Remove magic-number ambiguity.

### Current issue
The codebase mixes:

- scientifically meaningful defaults,
- internal numerical floors,
- performance heuristics,
- allocation sizing constants.

That makes it hard to know what should be configurable.

### Proposed structure

#### A. User-facing scientific defaults
Keep in the parameter registry and config dataclasses.

#### B. Internal numerical constants
Move into explicit constants modules, for example:

- `src/hulkrna/constants_numeric.py`
- `src/hulkrna/constants_performance.py`

#### C. Backend implementation constants
If shared with C++/nanobind code, keep a documented backend contract describing which values must stay mirrored.

### Immediate candidates to classify

#### User-facing or advanced-user-facing

- EM prior and convergence settings
- pruning threshold
- scoring penalties
- `tss_window`
- `nrna_frac_*`
- `gdna_*`
- `strand_prior_kappa`
- thread count
- memory spill settings such as `chunk_size` and `max_memory_bytes` if you want power users to control them

#### Internal-only

- `LOG_SAFE_FLOOR`
- `_NEG_INF`
- annotation growth floor constants
- C++ epsilon-style implementation safeguards

#### Needs a decision

- `max_frag_length`
- `chunk_size`
- annotation-table sizing parameters
- growth heuristics in `AnnotationTable`

These are operational knobs, not biological ones. They may belong in advanced configuration, but only if there is a real use case.

### Concrete refactor steps

1. Inventory all repeated constants.
2. Label each as:
   - public basic,
   - public advanced,
   - private internal,
   - backend-contract.
3. Remove duplicated declarations where possible.
4. Add comments explaining why each remaining private constant exists.

### Likely files

- `src/hulkrna/config.py`
- `src/hulkrna/scoring.py`
- `src/hulkrna/locus.py`
- `src/hulkrna/estimator.py`
- `src/hulkrna/pipeline.py`
- `src/hulkrna/annotate.py`
- `src/hulkrna/frag_length_model.py`
- `src/hulkrna/pyfallback.py`

### Acceptance criteria

- Every nontrivial number is either documented or removed.
- User-facing settings are easy to discover.
- Private numerical constants are obviously private.

---

## Phase 4 — Reduce Redundant Code Paths and Tighten Module Boundaries

### Goal
Reduce structural duplication while preserving performance.

### Highest-value targets

#### 4.1 Config/default duplication
This is the easiest and highest ROI redundancy to remove.

#### 4.2 Native/Python fallback coordination
There appears to be a good-faith effort to centralize some behavior, but the project still has multiple layers carrying similar defaults and assumptions.

The refactor goal is not to eliminate the fallback path. It is to ensure both backends consume the same contracts:

- same parameter objects,
- same named constants where appropriate,
- same result schema,
- same tests.

#### 4.3 Output/reporting duplication
Config snapshots, summaries, and output-file naming should be centralized instead of assembled inline.

#### 4.4 Buffer/annotation sizing logic
Sizing and growth heuristics are currently spread across `pipeline.py` and `annotate.py`. They should be defined in one place or clearly separated by concern.

### Concrete refactor steps

1. Define clean interfaces between:
   - scan stage,
   - EM/locus stage,
   - output/annotation stage.
2. Move “workflow glue” out of domain modules.
3. Audit for repeated transformations of the same data (e.g. repeated normalization of tags, repeated default conversion, repeated config serialization logic).
4. Add contract tests for backend parity.

### Likely files

- `src/hulkrna/pipeline.py`
- `src/hulkrna/scan.py`
- `src/hulkrna/annotate.py`
- `src/hulkrna/estimator.py`
- `src/hulkrna/pyfallback.py`
- native bindings under `src/hulkrna/native/` or related extension code

### Acceptance criteria

- fewer entrypoints with overlapping responsibilities,
- clearer ownership of each transformation,
- parity tests cover shared behavior.

---

## Phase 5 — Naming Cleanup and Readability Pass

### Goal
Make the code explain itself.

### Strategy
This should be an incremental pass, not a giant rename bomb.

### Renaming priorities

#### A. High-level public API names
These should change only if genuinely misleading.

#### B. Internal function names with vague purpose
Candidates should be renamed when their true role is clearer than the current name.

#### C. Local variable cleanup in dense scientific code
Examples of current readability pain points:

- `ctx` → `scoring_context`
- `bf` → `buffered_fragment`
- `fc` → `fragment_class`
- `ll` → `log_likelihood`
- `ct` → `count_column`
- `cw` → `coverage_weight`
- `rl` → `read_length_bp`
- `ebp` → `exonic_bases`
- `ibp` → `intronic_bases`
- `tx_s` / `tx_e` → `tx_start` / `tx_end`
- `t_arr` → `transcript_indices`
- `u_arr` → `unit_indices`
- `lid` → `locus_id`
- `lbl` → `component_label`

### Important nuance
In true hot loops, some short names may still be acceptable if performance and vectorized patterns matter. But scientific logic should not read like code golf.

### Supporting work
Create and enforce a naming glossary for recurring concepts:

- `unambig` vs `multimapping`
- `nrna_frac`
- `strand_specificity`
- `coverage_weight`
- `effective_length`
- `locus`
- `component`
- `candidate`

### Likely files

- `src/hulkrna/scan.py`
- `src/hulkrna/locus.py`
- `src/hulkrna/estimator.py`
- `src/hulkrna/pipeline.py`
- `src/hulkrna/scoring.py`

### Acceptance criteria

- New contributors can follow the main data flow without reverse-engineering abbreviations.
- Variable names distinguish biology, math, and implementation details.

---

## Phase 6 — Publish First-Class Parameter Documentation

### Goal
Make settings understandable to users.

### Proposed outputs

Create generated or semi-generated docs such as:

- `docs/parameters.md`
- `docs/advanced_parameters.md`
- `docs/config_schema.yaml` or `docs/example_config.yaml`

Each parameter should document:

- name,
- CLI flag,
- YAML key,
- type,
- default,
- allowed range,
- what it affects,
- when to change it,
- whether it is advanced,
- whether changing it may alter runtime or memory.

### Special recommendation
For this tool, documentation should explicitly separate:

- **scientific behavior knobs**,
- **performance/resource knobs**,
- **debugging or developer knobs**.

That distinction will massively reduce user confusion.

### Likely files

- new docs under `docs/`
- `README.md`
- generated examples referenced from CLI help

### Acceptance criteria

- A user can discover all important quant parameters without reading source.
- Advanced options are clearly marked and justified.
- CLI help and docs remain synchronized.

---

## Phase 7 — Production Hardening and Cleanup Sweep

### Goal
Finish the conversion from research code to production-ready code.

### Candidate work

- remove stale helper code and compatibility shims,
- remove dead comments that reference old architecture,
- standardize module docstrings,
- ensure consistent logging style,
- reduce giant functions where still present,
- improve README with real usage/config examples,
- add lint/type-check targets if desired,
- optionally add a small architectural overview diagram.

### Acceptance criteria

- no obvious dead paths,
- no duplicate parameter declarations,
- no orphan docs that describe outdated behavior,
- new contributors can identify the main workflow quickly.

---

## Priority-Ordered Implementation Sequence

If this work is done over multiple turns/PRs, this is the order I recommend.

### PR 1 — Parameter registry foundation

Scope:

- create parameter registry,
- migrate quant defaults into it,
- update CLI resolution,
- update tests.

Why first:

- highest leverage,
- removes drift immediately,
- unlocks better docs and advanced-option handling.

### PR 2 — Split `quant_command()` and output manifest handling

Scope:

- extract config resolution/build helpers,
- shorten CLI orchestration,
- centralize config snapshot writing.

Why second:

- makes later refactors easier,
- reduces cognitive load immediately.

### PR 3 — Constant classification and magic-number cleanup

Scope:

- classify constants,
- centralize internal constants,
- expose only justified advanced knobs.

Why third:

- now that parameter ownership is clear, you can safely separate public vs private numbers.

### PR 4 — Naming cleanup in `scan.py`, `locus.py`, `estimator.py`

Scope:

- high-impact readability renames,
- no behavioral changes beyond naming and helper extraction.

Why fourth:

- by this point interfaces are stable enough for safe renaming.

### PR 5 — Backend parity and duplication reduction

Scope:

- tighten Python/native contracts,
- remove repeated assumptions/defaults,
- add parity tests.

### PR 6 — Documentation generation and advanced config UX

Scope:

- docs pages,
- example YAML,
- README upgrades,
- perhaps `hulkrna quant --write-example-config` in the future.

---

## Which Parameters Should Become “Advanced” Options?

My recommendation is to explicitly divide user controls into three layers.

## Basic user options
These are reasonable for everyday use.

- input/output paths
- `seed`
- `include_multimap`
- `keep_duplicates`
- `annotated_bam`
- `threads`

Potentially also:

- `em_prior_alpha`
- `em_prior_gamma`
- `em_iterations`

if you expect scientifically sophisticated users.

## Advanced scientific options
These are valid user controls, but should be clearly marked expert-only.

- `prune_threshold`
- `em_convergence_delta`
- `overhang_alpha`
- `mismatch_alpha`
- `gdna_splice_penalty_unannot`
- `tss_window`
- all `nrna_frac_*`
- all `gdna_kappa_*` and `gdna_mom_*`
- `strand_prior_kappa`

## Advanced operational options
These affect runtime/memory more than model interpretation.

- `tmpdir`
- `chunk_size` if exposed
- `max_memory_bytes` if exposed
- `max_frag_length` if exposed

## Internal-only
These should remain code-level constants unless a real need emerges.

- log safety floors
- raw epsilon constants
- array growth floors
- sentinel values and internal clamps used only for implementation safety

---

## Specific Cleanup Opportunities Worth Targeting Early

### Opportunity 1 — Make `config.py` actually authoritative

This is the single most important cleanup step.

### Opportunity 2 — Turn parameter docs into generated artifacts

Once a registry exists, docs, CLI help, YAML examples, and saved manifests should all be derivable from it.

### Opportunity 3 — Unify run manifest / config snapshot behavior

The current `config.json` writing in `quant_command()` is useful and should stay, but it should become a deliberate, reusable “run manifest” layer.

### Opportunity 4 — Separate domain math from plumbing

Modules like `scan.py`, `locus.py`, and `estimator.py` contain real model logic. That logic becomes much easier to trust when orchestration and serialization are moved elsewhere.

### Opportunity 5 — Use the test suite as architectural scaffolding

Refactors should add tests for:

- parameter registry completeness,
- CLI/YAML precedence,
- serialization of config snapshots,
- backend parity,
- no-regression runtime envelopes.

---

## Risks and How To Control Them

## Risk 1 — Refactor changes numerical behavior accidentally

Mitigation:

- snapshot parity tests,
- benchmark comparisons,
- small-scenario golden outputs,
- avoid combining naming cleanup with math changes.

## Risk 2 — Exposing too many knobs confuses users

Mitigation:

- basic vs advanced split,
- strong docs,
- stable defaults.

## Risk 3 — Performance regresses while making code cleaner

Mitigation:

- preserve hot-loop structure where needed,
- benchmark after each major phase,
- refactor orchestration first, hot loops second.

## Risk 4 — Native and Python code paths drift further apart

Mitigation:

- contract tests,
- shared parameter objects,
- explicit backend interface ownership.

---

## Definition of Done

This cleanup initiative should be considered successful when all of the following are true:

- there is one obvious home for every user-facing parameter,
- defaults are not duplicated across CLI/tests/core modules,
- advanced options are clearly identified and documented,
- top-level workflow code is short and understandable,
- major modules have clearer names and smaller responsibilities,
- magic numbers are classified and documented,
- tests and benchmarks protect refactors,
- README/docs explain how to configure the model without reading source.

---

## Recommended Next Step

Start with **Phase 1: parameter registry and config centralization**.

That work has the best ratio of impact to risk because it immediately improves:

- duplication,
- discoverability,
- CLI/YAML behavior,
- test maintainability,
- future documentation generation.

It also creates the foundation for the rest of the plan.
