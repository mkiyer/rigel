# Rigel: Cleanup Implementation Plan

> Goal: A single, non-redundant code path where most implementation is in C++, and the Python→C++ boundary is obvious to an outside reviewer.

---

## Table of Contents

1. [Current State Assessment](#1-current-state-assessment)
2. [Priority 1 — Eliminate Redundancy & Clarify Boundaries](#2-priority-1--eliminate-redundancy--clarify-boundaries)
3. [Priority 2 — Consolidate the Python↔C++ Interface](#3-priority-2--consolidate-the-pythonc-interface)
4. [Priority 3 — Move Remaining Hot Python to C++](#4-priority-3--move-remaining-hot-python-to-c)
5. [Priority 4 — Structural Cleanup](#5-priority-4--structural-cleanup)
6. [Priority 5 — Documentation & Testing](#6-priority-5--documentation--testing)
7. [Files Untouched (No Changes Needed)](#7-files-untouched-no-changes-needed)
8. [Implementation Order](#8-implementation-order)

---

## 1. Current State Assessment

### What's Good

The legacy Python implementations (bam.py, em.py, fragment.py) have already been removed. All hot paths are in C++. The codebase is functional and correct.

### What's Messy

| Issue | Where | Impact |
|-------|-------|--------|
| **Duplicated chunk construction** | pipeline.py and buffer.py both build `_FinalizedChunk` from raw C++ output with subtly different dtype conversions | Maintenance risk, confusing to reviewers |
| **`_resolver` private attribute leaks across modules** | index.py creates it; pipeline.py, resolution.py, annotate.py, scoring.py all access `index._resolver` directly | Hidden dependency, unclear ownership |
| **scoring.py has dead imports and test-only functions mixed in** | `compute_fragment_weight` imported but unused; `genomic_to_transcript_pos` and `score_gdna_standalone` are test utilities living in production code | Confusing — reviewer can't tell what's production vs test |
| **`_FinalizedChunk` is ~300 lines of boilerplate** | Duplicated numpy bookkeeping for 18 arrays plus classification, iteration, sorting | Verbose, fragile if schema changes |
| **buffer.py uses `__del__` for cleanup** | Destructor-based cleanup is unreliable in Python (no guaranteed ordering) | Resource leak risk |
| **scan.py constructs contiguous numpy arrays from chunk** | ~30 lines of `.ascontiguousarray()` calls before passing to C++ | Boilerplate that could be eliminated |
| **pipeline.py is a "god function"** | `quant_from_buffer()` is ~250 lines doing scoring setup + memory management + EB priors + locus EM + cleanup | Hard to follow, hard to test in isolation |
| **estimator.py mixes data structures, EM dispatch, prior computation, and output formatting** | 1,470 lines with 5 different responsibilities | Largest file, most entangled |
| **strand_model.py is 563 lines for a Beta-Binomial** | Verbose property implementations, duplicated likelihood paths (_int vs enum) | Overengineered for what it does |
| **`build_locus_em_data()` is 260 lines of Python numpy** | Per-locus CSR extraction + deduplication + renumbering | Performance-sensitive, pure Python, could be C++ |

### Code Path Clarity Score

An outside reviewer would need to trace through **4–5 files** to understand the main pipeline, and would encounter:

- Private attribute access patterns (`index._resolver`)
- Two different chunk construction paths (pipeline.py vs buffer.py)
- Test utilities in production scoring module
- 1,470-line estimator file mixing concerns
- No "map" of what calls what (until now — see CODE_PATH.md)

---

## 2. Priority 1 — Eliminate Redundancy & Clarify Boundaries

### 2.1 Unify `_FinalizedChunk` construction

**Problem:** pipeline.py lines 328–361 and buffer.py `_finalize_native()` lines 505–532 both construct `_FinalizedChunk` from raw C++ dict output. They use different dtype conversions:
- pipeline.py: `np.frombuffer(..., dtype=np.int64).copy()` for t_offsets/frag_id
- buffer.py: `np.frombuffer(..., dtype=np.int64).astype(np.int32)` for t_offsets/frag_id

**Fix:** Add a `_FinalizedChunk.from_raw(raw_dict)` class method in buffer.py that is the **single** way to construct a chunk from C++ output. Use it in both pipeline.py and buffer.py.

**Scope:** ~50 lines changed across 2 files.

### 2.2 Make `_resolver` a public attribute

**Problem:** `index._resolver` is accessed directly by 4 modules (pipeline.py, resolution.py, annotate.py, and indirectly scoring.py via the C++ context). The underscore-prefix convention signals "don't touch this," but everyone touches it.

**Fix:** Rename to `index.resolver` (public attribute). Add a `@property` if lazy initialization is needed. Document in `TranscriptIndex` docstring that this is the C++ `FragmentResolver` instance.

**Scope:** ~15 lines changed across 5 files.

### 2.3 Move test-only functions out of scoring.py

**Problem:** `genomic_to_transcript_pos()` and `score_gdna_standalone()` are only used by tests but live in the production scoring module. `compute_fragment_weight` is imported but never used.

**Fix:**
- Move `genomic_to_transcript_pos()` and `score_gdna_standalone()` to `tests/conftest.py` or a `tests/helpers.py` module.
- Remove the dead `compute_fragment_weight` import from scoring.py.

**Scope:** ~80 lines moved, ~2 lines deleted.

### 2.4 Resolve the int32 frag_id narrowing

**Problem:** buffer.py narrows frag_id from int64 to int32 (max ~2.1 billion fragments). pipeline.py keeps it as int64. This inconsistency is undocumented and would silently corrupt data for very large BAM files.

**Fix:** Pick one dtype and use it consistently. If int32 is intentional (it is — BAM files rarely exceed 2B fragments), document it in the `_FinalizedChunk` docstring and add a runtime check/warning.

**Scope:** ~10 lines.

---

## 3. Priority 2 — Consolidate the Python↔C++ Interface

### 3.1 Create an explicit `native.py` interface module

**Problem:** C++ module imports are scattered across pipeline.py, scoring.py, buffer.py, estimator.py, locus.py, index.py, and annotate.py. Each imports different pieces from `_bam_impl`, `_scoring_impl`, `_resolve_impl`, `_em_impl`, `_cgranges_impl`. An outside reviewer has no single place to see "what does the C++ layer provide?"

**Fix:** Create `src/rigel/native.py` (or `_native_api.py`) that re-exports all C++ symbols used by Python code:

```python
"""Public interface to Rigel's C++ native extensions.

All C++ functionality used by Python code is imported through this module.
"""
# BAM scanning
from ._bam_impl import BamScanner, BamAnnotationWriter, detect_sj_strand_tag

# Fragment resolution
from ._resolve_impl import FragmentResolver, FragmentAccumulator, ResolvedFragment

# Scoring
from ._scoring_impl import NativeFragmentScorer

# EM solver
from ._em_impl import (
    batch_locus_em,
    connected_components,
    EM_PRIOR_EPSILON,
)

# Interval overlap
from ._cgranges_impl import cgranges
```

Other modules then import from `rigel.native` instead of scattered `_*_impl` modules. This makes the C++ boundary a single, visible seam.

**Scope:** New file (~30 lines) + ~20 import lines changed across 7 files.

### 3.2 Eliminate the numpy contiguous-array boilerplate in scan.py

**Problem:** scan.py lines ~100–140 convert every chunk array to contiguous numpy before passing to C++. This is pure boilerplate that exists because the chunk arrays *might* not be contiguous.

**Fix:** Either (a) guarantee contiguity at chunk construction time (in `_FinalizedChunk.from_raw()`), or (b) move the contiguity enforcement into the C++ `fused_score_buffer` (nanobind can handle this automatically with `nb::ndarray<>` flags).

**Scope:** ~40 lines removed from scan.py.

---

## 4. Priority 3 — Move Remaining Hot Python to C++

These are Python functions that run in the quantification critical path and would benefit from C++ implementation. They are candidates for future migration, not immediate cleanup.

### 4.1 Move `build_locus_em_data()` into C++ (locus.py → em_solver.cpp)

**Problem:** This 260-line function runs once per locus, extracting a local CSR sub-problem from the global CSR. It involves deduplication, renumbering, gDNA candidate injection, and array sorting — all in numpy. For datasets with many loci, this is a significant fraction of total Python time.

**Current state:** The C++ `batch_locus_em()` already receives flattened locus definitions. The per-locus data extraction could be folded into the C++ side.

**Benefit:** Eliminates the last significant Python loop over loci. The entire path from `ScoredFragments` → EM results would be a single C++ call.

**Scope:** ~300 lines moved from Python to C++. Medium complexity.

### 4.2 Move EB prior computation into C++ (locus.py + estimator.py → em_solver.cpp)

**Problem:** `compute_eb_gdna_priors()` (~160 lines) and `compute_nrna_frac_priors()` (~365 lines) are pure numpy. They're not hot-path (run once), but they break the "C++ does all computation" narrative and force data to round-trip between Python and C++.

**Benefit:** The entire quantification phase (after model training) would be a single C++ call: `run_quantification(scored_fragments, config) → results`. This is the cleanest possible architecture.

**Risk:** These functions encode complex statistical logic (3-tier EB shrinkage, Method-of-Moments κ estimation). Porting them to C++ increases maintenance burden and makes the math harder to iterate on.

**Recommendation:** Defer this unless profiling shows these functions are a bottleneck. They're clean, correct, and well-tested in Python.

**Scope:** ~525 lines moved. High complexity. Low urgency.

### 4.3 Move fragment-length effective-length computation into C++

**Problem:** `frag_length_model.py:compute_all_transcript_eff_lens()` is called during Phase 3 setup. It's vectorized numpy (Salmon-style eCDF convolution) and runs once.

**Recommendation:** Defer. Single call, already vectorized.

---

## 5. Priority 4 — Structural Cleanup

### 5.1 Split estimator.py into focused modules

**Problem:** estimator.py is 1,470 lines with 5 responsibilities:
1. `ScoredFragments` data structure (CSR container)
2. `Locus` and `LocusEMInput` data structures
3. `AbundanceEstimator` class (state + EM dispatch)
4. Prior computation functions (`compute_nrna_frac_priors`, `estimate_kappa`, etc.)
5. Output formatting (`get_counts_df`, `get_gene_counts_df`, etc.)

**Fix:** Split into:
- `scored_fragments.py` — `ScoredFragments` dataclass (~85 lines)
- `locus_types.py` — `Locus`, `LocusEMInput` dataclasses (~125 lines), currently in estimator — move to locus.py alongside `build_loci` and `build_locus_em_data`
- `priors.py` — EB prior functions: `compute_nrna_frac_priors`, `estimate_kappa`, `compute_global_gdna_density`, etc. (~400 lines)
- `estimator.py` — `AbundanceEstimator` class (EM dispatch + output) (~860 lines, down from 1,470)

**Alternative (minimal):** Keep as one file but reorganize with clear section headers and move `Locus`/`LocusEMInput` to locus.py where they logically belong.

**Scope:** 4 new files, ~100 import updates. Medium effort.

### 5.2 Decompose `quant_from_buffer()` in pipeline.py

**Problem:** This function is ~250 lines doing 5 distinct phases (scoring setup, memory management, EB priors, locus EM, cleanup). Hard to test in isolation.

**Fix:** Extract into named sub-functions that pipeline.py calls in sequence:

```python
def quant_from_buffer(buffer, index, ...):
    scorer, router, estimator = _setup_scoring(buffer, index, ...)
    em_data = _score_fragments(router, buffer)
    _compute_priors(estimator, em_data, index, ...)
    _run_em(estimator, em_data, ...)
    return estimator
```

**Scope:** ~50 lines of function extraction, no logic changes.

### 5.3 Replace `__del__` with `weakref.finalize` in buffer.py

**Problem:** `FragmentBuffer.__del__` calls `self.cleanup()` but `__del__` is unreliable — no guaranteed call order, may not run if interpreter is shutting down.

**Fix:** Use `weakref.finalize(self, shutil.rmtree, self._spill_dir)` in `__init__` instead.

**Scope:** ~5 lines.

### 5.4 Simplify strand_model.py

**Problem:** 563 lines for a Beta-Binomial model with duplicated likelihood paths (`strand_likelihood` vs `strand_likelihood_int`).

**Fix:**
- Remove `strand_likelihood()` (Strand enum version) — only `strand_likelihood_int()` is used in production (via the C++ scorer).
- Consolidate `observe()` and `observe_batch()` — if `observe()` is never called with single values in production (check: it's replayed from C++ batch arrays), remove the single-observation path.
- Remove `posterior_95ci()` if unused in production (check: only used in `log_summary()`).

**Scope:** ~50 lines removed.

---

## 6. Priority 5 — Documentation & Testing

### 6.1 Add module-level docstrings

Every `.py` and `.cpp` file should have a 2–3 line module docstring stating:
- What this module does
- What phase of the pipeline it belongs to
- Whether it's Python-only or a C++ interface

### 6.2 Add architecture diagram to README

Include a simplified version of the data flow diagram from CODE_PATH.md.

### 6.3 Add a smoke test for the full quant pipeline

Currently tests cover individual components. A single integration test that runs `rigel quant` on a small BAM + index and validates output schema would catch import/interface breakage.

### 6.4 Add cross-chunk multimapper regression test

Flagged in prior review — multimapper scoring happens in C++ `fused_score_buffer()`, but there's no test covering fragments that span multiple buffer chunks.

---

## 7. Files Untouched (No Changes Needed)

These files are clean, focused, and don't need modification:

| File | Reason |
|------|--------|
| `cli.py` | Well-organized, declarative parameter registry, single responsibility |
| `config.py` | Frozen dataclasses, clean defaults, no logic |
| `types.py` | Pure type definitions, used everywhere, stable |
| `splice.py` | Pure enum + column indexing, minimal |
| `stats.py` | Clean dataclass counter, no logic complexity |
| `gtf.py` | Streaming parser, well-tested |
| `transcript.py` | Simple dataclass + GTF reader |
| `bias.py` | Tight prefix-sum implementation, correct |
| `sim/*` | Test/benchmarking infrastructure, separate concern |
| `cgranges/*` | Vendored C library + nanobind wrapper |
| `native/constants.h` | Enum mirrors, must stay in sync with types.py |
| `native/fast_exp.h` | Optimized math, self-contained |
| `native/thread_pool.h` | Self-contained thread pool |
| `native/thread_queue.h` | Self-contained bounded queue |
| `native/simd_detect.h` | Build-time detection, no logic |

---

## 8. Implementation Order

Recommended sequencing to minimize disruption:

### Phase A — Quick wins (1–2 days)

1. **2.1** Unify `_FinalizedChunk` construction → `from_raw()` class method
2. **2.3** Move test-only functions out of scoring.py
3. **2.4** Document/standardize int32 frag_id
4. **5.3** Replace `__del__` with `weakref.finalize`
5. **2.2** Make `_resolver` public

### Phase B — Interface clarity (2–3 days)

6. **3.1** Create `native.py` C++ interface module
7. **3.2** Eliminate contiguous-array boilerplate in scan.py
8. **5.2** Decompose `quant_from_buffer()`
9. **5.4** Simplify strand_model.py

### Phase C — Structural reorganization (3–5 days)

10. **5.1** Split estimator.py (move `Locus`/`LocusEMInput` to locus.py, extract priors.py)
11. **6.1–6.4** Documentation and testing

### Phase D — Future C++ migration (optional, larger effort)

12. **4.1** Move `build_locus_em_data()` into C++
13. **4.2** Move EB prior computation into C++ (only if profiling justifies)

---

## Summary

The codebase is **functionally correct** and the legacy Python→C++ migration is **mostly complete**. The remaining cleanup is about **clarity**, not correctness:

- **5 changes** eliminate redundancy and leaky abstractions (Priority 1)
- **2 changes** make the C++ boundary visible and clean (Priority 2)
- **4 changes** improve module structure and reduce file sizes (Priority 4)
- **2 optional changes** would complete the C++ migration (Priority 3)

After this cleanup, an outside reviewer would see:
1. `cli.py` → `pipeline.py` → clear phase sequence
2. `native.py` → single seam showing all C++ exports
3. Each `.py` file < 900 lines with a single responsibility
4. CODE_PATH.md explaining the full algorithm
