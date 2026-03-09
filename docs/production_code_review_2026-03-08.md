# Rigel Production Code Review

Date: 2026-03-08

Scope reviewed:
- Python production code under `src/rigel` (excluding most of `src/rigel/sim`)
- Native code under `src/rigel/native`
- Packaging and release configuration (`pyproject.toml`, `CMakeLists.txt`, CI workflows)

## Executive summary

Rigel has a strong technical core. The project already has several production-positive traits: a clear domain model, extensive tests, a documented algorithm, a mostly coherent Python/C++ split, and attention to determinism in the EM implementation. The main risks are no longer algorithmic correctness at a glance; they are maintainability, packaging robustness, resource lifetime, and a few concrete implementation defects that will become expensive once this code is widely deployed.

Update after the first implementation pass:

- The `_bam_impl` packaging/runtime contract has been materially improved. The build now fails if `htslib` is unavailable, and the publish workflow imports the native extensions instead of only `import rigel`.
- The multimapper scoring path has been moved deeper into C++, which is a real reduction in Python hot-path complexity and a meaningful maintainability/performance improvement.
- The spill-path memory issue is only partially addressed. The code no longer does `list(buffer.iter_chunks())`, but it still accumulates all converted chunk arrays in memory before calling the fused C++ scorer.

The biggest remaining production blockers are:

1. `FragmentBuffer`'s spill-to-disk design is still partially defeated in the scan stage.
2. At least one production reporting path contains unfinished placeholder logic.
3. Several boundary and lifecycle issues remain unresolved (`_resolver` leakage, destructor cleanup, duplicated native decode logic, integer narrowing).

The codebase is closer to production than it was at the start of review, but it still needs a deliberate cleanup pass focused on true spill-safe scanning, output semantics, lifecycle ownership, and boundary design.

## Validation snapshot

Validated in the `rigel` conda environment against the current workspace copy:

- `pytest tests/test_pipeline_routing.py tests/test_buffer.py -q`
- Result: 52 tests passed

What this validation proves:

- The new C++ multimapper routing path is consistent with the current focused test suite.
- Buffer spill/frag-id round-trip behavior covered by `tests/test_buffer.py` still passes.

What it does not yet prove:

- True bounded-memory behavior under forced spill during end-to-end scan/quant.
- Cross-chunk multimapper-group correctness.
- End-to-end wheel validation of a real quantification path.

## Findings

### 1. High: spill-path memory behavior is improved but not actually fixed end to end

Files:
- `src/rigel/scan.py:124`
- `src/rigel/scan.py:126`
- `src/rigel/scan.py:127`
- `src/rigel/buffer.py:600`
- `src/rigel/buffer.py:607`

The earlier `chunks = list(buffer.iter_chunks())` materialization is gone, which is an improvement. However, `_scan_native()` still iterates through every chunk, converts each chunk into contiguous numpy arrays, and appends those arrays to `chunk_arrays` before invoking the fused C++ scorer once at the end.

Why this matters:
- The raw spilled chunks are no longer all retained at once, but the converted arrays still are.
- Large production runs can still end up rehydrating the full scored input in memory before C++ processing starts.
- The code and docstrings now imply stronger streaming guarantees than the implementation actually provides.

Recommendation:
- Treat this finding as partially resolved, not closed.
- Either refactor the fused scorer to consume chunks incrementally, or introduce an explicit spill-aware staging interface on the native side.
- Add a regression test that forces spill and verifies bounded memory growth through the scan path, not just correct round-trip data reload.

### 2. Resolved at the primary contract level: `_bam_impl` is now required at build time

Files:
- `src/rigel/pipeline.py:50`
- `src/rigel/pipeline.py:51`
- `src/rigel/annotate.py:280`
- `CMakeLists.txt:163`
- `.github/workflows/publish.yml:89`

This was one of the most important issues in the original review, and the main contract problem has now been fixed. `CMakeLists.txt` fails fast when `htslib` is missing, which matches the unconditional runtime imports in `pipeline.py` and `annotate.py`.

Why this matters:
- The build/import/runtime behavior is now internally consistent.
- This removes a major class of post-release import failures.

Recommendation:
- Keep the required-native policy.
- Update installation and contributor docs to say this explicitly.
- Do not reopen the optional-native path unless there is a real distribution requirement.

### 3. Medium: release smoke testing is stronger, but still import-only

Files:
- `.github/workflows/publish.yml:89`

The publish workflow no longer only imports `rigel`. It now imports multiple native extensions directly, which is a meaningful improvement. However, it still does not run a minimal end-to-end quantification path on the built wheel.

Why this matters:
- Native import validation is better than before, but it still does not prove the packaged quant path works.
- Packaging problems can still hide in module wiring, runtime data loading, or small integration seams that import tests miss.

Recommendation:
- Keep the native import smoke test.
- Add one canonical wheel smoke test that builds or loads a tiny index and runs a minimal quant invocation.
- Reuse the same smoke asset in CI and publish.

### 4. Medium: multimapper C++ migration is a strong improvement, but coverage is still missing at one critical seam

Files:
- `src/rigel/native/scoring.cpp`
- `src/rigel/scan.py`
- `tests/test_pipeline_routing.py`

Moving multimapper scoring into `fused_score_buffer()` is the right direction. It removes a large amount of Python hot-path logic, eliminates several Python/C++ round-trips, and simplifies `scan.py` substantially.

Why this matters:
- This is one of the best changes in the current implementation set.
- It reduces hot-path complexity and narrows the production path to one native scoring implementation.

Remaining risk:
- I did not find a dedicated regression test for multimapper groups that span chunk boundaries.
- The C++ implementation explicitly carries pending multimapper state across chunks, which is correct in design, but that boundary is exactly where a silent counting bug would hide.

Recommendation:
- Add a targeted test where the same `frag_id` appears at the end of one chunk and the start of the next.
- Add at least one test that exercises mixed splice-type multimapper groups across that boundary.

### 5. Medium: gene-level `locus_id` output contains unfinished placeholder logic

File:
- `src/rigel/estimator.py:1363`

`AbundanceEstimator.get_gene_counts_df()` contains a `pass  # keep current assignment` in the branch that should reconcile multiple transcript loci for one gene. This is not just a style issue; it means a published output field is backed by incomplete conflict-resolution logic.

Why this matters:
- Output semantics are unclear when a gene spans multiple transcript loci.
- The current behavior is effectively “first assignment wins”, but the code comments describe something more deliberate.
- Production reporting code should not contain placeholders in ambiguous-data paths.

Recommendation:
- Define a real policy for gene-level locus attribution.
- Examples: dominant transcript by mRNA mass, multi-locus sentinel, or emit a separate gene-to-loci relation table instead of forcing a single `locus_id`.
- Add tests that cover multi-locus gene cases and pin the expected semantics.

### 6. Medium: private native resolver state leaks across module boundaries

Files:
- `src/rigel/index.py:905`
- `src/rigel/pipeline.py:291`
- `src/rigel/resolution.py:177`
- `src/rigel/annotate.py:290`

`TranscriptIndex` constructs native resolver state as a private attribute (`self._resolver`), but the rest of the codebase treats that private field as a shared service locator. Multiple modules reach inside the index to use it directly.

Why this matters:
- `TranscriptIndex` leaks internal construction details into unrelated modules.
- The object graph is harder to reason about and harder to refactor safely.
- This is exactly how “works now, impossible to untangle later” code evolves.

Recommendation:
- Expose a public API for the capabilities being used: for example `index.resolve_fragment(...)`, `index.native_resolver`, or a dedicated `ResolutionContext` wrapper.
- Keep ownership and lifecycle in one place.
- Treat private attributes as truly private again.

### 7. Medium: resource cleanup depends on `__del__`, which is not a reliable lifecycle boundary

File:
- `src/rigel/buffer.py:633`

`FragmentBuffer` uses `__del__` to call `cleanup()`. Destructor-based cleanup is non-deterministic and can behave poorly during interpreter shutdown, exception unwinding, reference cycles, and partial object initialization.

Why this matters:
- Temporary spill directories are real filesystem resources, not trivial Python objects.
- Production cleanup should be explicit and deterministic.
- The current code is mostly safe because `cleanup()` is idempotent and uses `ignore_errors=True`, but the ownership model is still fragile.

Recommendation:
- Remove `__del__` and rely on explicit cleanup plus context-manager usage.
- If a last-resort fallback is still desired, use `weakref.finalize` rather than destructor behavior.
- Add tests for interrupted runs and exception paths to verify spill directory cleanup.

### 8. Medium: native-buffer decoding logic is duplicated in multiple places

Files:
- `src/rigel/buffer.py:511`
- `src/rigel/buffer.py:515`
- `src/rigel/pipeline.py:325`
- `src/rigel/pipeline.py:329`

The byte-dict returned from native accumulation is decoded into `_FinalizedChunk` in both `FragmentBuffer._finalize_native()` and `pipeline.scan_and_buffer()`. The field-by-field conversions are nearly identical.

Why this matters:
- Changes to the native schema can drift between the two decode paths.
- Duplicated marshaling code is high-friction and easy to break.
- This is boundary code between Python and C++; it should exist exactly once.

Recommendation:
- Introduce a single decoder function or `_FinalizedChunk.from_native_raw(...)` constructor.
- Keep dtype validation there.
- Make the native output schema explicit and documented.

### 9. Medium: large-run identifiers are narrowed to `int32` in buffer serialization

File:
- `src/rigel/buffer.py:528`

`frag_id` is read from native output as `int64` and immediately cast to `int32`. `t_offsets` is also narrowed to `int32` elsewhere in the buffer serialization path.

Why this matters:
- It silently bakes in data-size ceilings that are not documented.
- For production batch systems or pooled runs, fragment counts can become very large.
- Silent narrowing is much worse than an explicit checked limit.

Recommendation:
- Keep `frag_id` and offsets in `int64` unless there is a measured reason not to.
- If a hard upper bound is intentional, enforce it with a checked conversion and a clear error.
- Document the expected maximum chunk and dataset sizes.

## Architecture and maintainability observations

These are not all defects, but they are where future spaghetti risk is coming from.

### Oversized modules

The following files are too large for long-term maintainability:

- `src/rigel/estimator.py` (1485 lines)
- `src/rigel/index.py` (994 lines)
- `src/rigel/cli.py` (936 lines)
- `src/rigel/locus.py` (786 lines)
- `src/rigel/pipeline.py` (763 lines)
- `src/rigel/buffer.py` (649 lines)
- `src/rigel/native/em_solver.cpp` (2683 lines)
- `src/rigel/native/bam_scanner.cpp` (1929 lines)
- `src/rigel/native/scoring.cpp` (1510 lines)
- `src/rigel/native/resolve_context.h` (1091 lines)

This is the single biggest maintainability smell in the codebase. Even where the code is technically sound, the cognitive load is too high. New contributors will need too much context to make safe changes.

Recommended module split directions:

- `cli.py`: split parser construction, YAML/config merging, subcommand handlers, and output writing.
- `pipeline.py`: split scan stage, quant stage, output/annotation stage, and native bridge helpers.
- `estimator.py`: split priors, EM batch bridge, output tables, and statistics/reporting.
- `index.py`: split build-time index generation from load-time runtime index construction.
- `em_solver.cpp`: split equivalence-class construction, bias correction, EM kernels, SQUAREM, and nanobind module glue.
- `bam_scanner.cpp`: split BAM reading, fragment grouping, statistics, and accumulator export.

### Boundary objects are under-specified

The Python/C++ seam relies heavily on raw dicts, NumPy arrays, and informal field contracts. That is efficient, but it is not self-describing.

Recommended cleanup:

- Define explicit boundary schemas for native outputs.
- Centralize all decode and validation logic.
- Keep Python-side ownership of only domain objects, not transport details.

### “Thin orchestrator” modules are carrying too much domain logic

`pipeline.py` describes itself as a thin orchestrator, but it currently knows about:

- native scanner wiring
- strand and fragment-length model replay
- buffer decoding
- transcript geometry construction
- TSS grouping
- empirical Bayes prior setup
- locus metadata formatting
- annotation-table lifecycle

That is a service layer plus data marshaling layer plus partial domain layer all in one file.

Recommended cleanup:

- Introduce explicit stage objects or pure helper modules with narrow responsibilities.
- Keep `run_pipeline()` as a readable composition root.

## Detailed cleanup and improvement opportunities

### Priority 1: correctness and release safety

1. Finish the spill-path fix by removing all-chunk accumulation before the fused C++ call.
2. Replace the `pass` placeholder in gene-level locus attribution with real semantics.
3. Add a wheel smoke test that runs one minimal quant path, not just native imports.
4. Add a regression test for cross-chunk multimapper groups.

### Priority 2: simplify boundaries and ownership

1. Replace raw native result dict decoding with one canonical adapter.
2. Remove direct `_resolver` access from other modules.
3. Remove destructor-based cleanup and use explicit ownership/finalizers.
4. Audit all narrowing casts across Python/native boundaries and document invariants.

### Priority 3: reduce module size and improve navigability

1. Split `cli.py` into parser, config resolution, and command execution modules.
2. Split `pipeline.py` into scan, quant, annotations, and report-writing modules.
3. Split `estimator.py` into priors, EM bridge, outputs, and metrics.
4. Split `index.py` into build pipeline and runtime index loader.
5. Split native files by algorithmic stage rather than “one massive extension file per subsystem”.

### Priority 4: strengthen production engineering quality

1. Add a formatter/linter gate to CI if not already enforced elsewhere. `ruff` is configured in `pyproject.toml`, but I did not see it run in CI.
2. Add a type-checking pass for the Python layer, especially around pandas/NumPy/native interfaces.
3. Add a dedicated smoke-test dataset small enough to run in CI and publish workflows.
4. Add explicit compatibility notes for required native dependencies and failure modes.
5. Make line between public and private APIs explicit in module design.

### Priority 5: documentation consistency

1. Align installation/runtime documentation with the real `_bam_impl` and `htslib` requirements.
2. Document supported failure modes when native modules are absent.
3. Document size assumptions for chunking, offsets, and fragment IDs.
4. Add a contributor guide describing module boundaries and native/Python interface rules.

## Positive notes

The codebase already has several foundations worth preserving:

- The algorithmic intent is well documented in code and docs.
- There is substantial test coverage, including native-path tests and determinism checks.
- The project has clearly invested in performance engineering rather than layering premature abstractions everywhere.
- Many of the current issues are fixable with structural cleanup rather than a redesign of the scientific model.
- The current implementation pass made real progress: the `_bam_impl` contract is tighter, publish validation is better, and the Python multimapper hot path has been removed from production scanning.

## Suggested productionization sequence

1. Finish the spill behavior work for real, including a forced-spill regression.
2. Add missing regression coverage for cross-chunk multimappers and wheel-level quant smoke.
3. Clean up output semantics and lifecycle ownership next.
4. Then do one focused refactor pass on module boundaries, starting with `pipeline.py`, `estimator.py`, and `index.py`.

The code does not need to be rewritten. It needs its boundaries tightened, its contracts made explicit, and its large modules carved into smaller units that match the actual architecture already present in the system.
