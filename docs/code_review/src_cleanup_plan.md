# Source Code Cleanup Plan

**Date:** 2026-04-01  
**Scope:** `src/rigel/` — Python and C++ modules  
**Goal:** Remove dead code, consolidate redundancies, and reduce complexity

---

## Overall Assessment

The codebase is well-maintained with no egregious technical debt. No TODO/FIXME/HACK markers exist in production code. All imports are used. The findings below are categorized by priority.

---

## Priority 1 — Dead Code Removal

These items are confirmed dead and can be removed immediately with no behavioral change.

### 1.1 Remove unused `_splice_type_label()` in `annotate.py`

**File:** `src/rigel/annotate.py` line 282  
**Issue:** Function defined but never called anywhere in `src/` or `tests/`.  
**Action:** Delete the function (6 lines).

### 1.2 Remove legacy `resolve()` method from `resolve_context.h`

**File:** `src/rigel/native/resolve_context.h` lines 1134–1195  
**Issue:** The 13-element-tuple legacy entry point is exposed via nanobind in `resolve.cpp` line 92 but never called from Python. All callers use `resolve_fragment()` instead.  
**Action:** Remove the `resolve()` method (~60 lines) from `resolve_context.h` and its `.def("resolve", ...)` binding in `resolve.cpp`. Update any docstrings that reference it.

### 1.3 Remove unused `#include <sys/stat.h>` in `bam_scanner.cpp`

**File:** `src/rigel/native/bam_scanner.cpp` line 44  
**Issue:** No `stat()`, `fstat()`, or `S_IS*` macros are used anywhere in the file.  
**Action:** Delete the include line.

### 1.4 Delete the `deprecated/` directory

**Location:** Top-level `deprecated/` directory  
**Issue:** Contains only `benchmark_old/` — an obsolete version of the benchmarking framework. Not imported or referenced by any current code. The current benchmarking system lives in `scripts/benchmarking/`.  
**Action:** Remove the entire `deprecated/` directory.

---

## Priority 2 — Consolidate Redundant Code

### 2.1 Extract shared TPM calculation helper in `estimator.py`

**File:** `src/rigel/estimator.py` lines 370, 490, 584  
**Issue:** The same TPM calculation pattern is repeated 3 times:
```python
rpk = count / eff
rpk_sum = rpk.sum()
tpm = (rpk / rpk_sum * 1e6) if rpk_sum > 0 else np.zeros_like(rpk)
```
Appears in `get_counts_df()`, `get_gene_counts_df()`, and `get_nrna_counts_df()`.  
**Action:** Extract a `_compute_tpm(counts, effective_lengths)` helper and call it from all three methods.

### 2.2 Consolidate ground-truth counting methods in `sim/scenario.py`

**File:** `src/rigel/sim/scenario.py` lines 87–230  
**Issue:** Multiple near-duplicate method pairs for counting fragments from FASTQ vs BAM:
- `ground_truth_gdna_count()` (FASTQ) / `ground_truth_gdna_count_from_bam()` (BAM)  
- `ground_truth_nrna_count()` (FASTQ) / `ground_truth_nrna_count_from_bam()` (BAM)  
- `ground_truth_from_fastq()` / `ground_truth_from_bam()` (mRNA counts)

The `ground_truth_auto()` method already demonstrates the dispatch pattern. The FASTQ and BAM variants share logic (prefix matching, counting) but differ only in source.  
**Action:** Refactor to use a shared `_count_fragments(prefix_filter, source)` internal method, or at minimum use `_parse_bam_qnames()` as the unified qname iterator and add a parallel `_parse_fastq_qnames()`, then implement the counters as thin wrappers.

### 2.3 Consolidate duplicate probability computation in `sim/reads.py`

**File:** `src/rigel/sim/reads.py` lines 256–310  
**Issue:** `_compute_probs()` and `_compute_nrna_probs()` are structurally identical — both compute `max(0, length - frag_len + 1) × abundance` and normalize. They differ only in input arrays (`_t_lengths`/`abundance` vs `_premrna_lengths`/`nrna_abundance`).  
**Action:** Merge into a single parameterized `_compute_probs(lengths, abundances, frag_len)` method.

### 2.4 Centralize read-name format documentation

**Files:** `src/rigel/sim/oracle_bam.py` lines 14–22, `src/rigel/sim/reads.py` lines 16–27  
**Issue:** The simulation read-name encoding format (`{t_id}:{start}-{end}:{strand}:{idx}`) is documented independently in both files.  
**Action:** Move the canonical format documentation to `src/rigel/sim/__init__.py` or a shared docstring, and reference it from both files.

---

## Priority 3 — Simplify Complex Functions

These are optional improvements for long-term maintainability. Each is self-contained and can be done incrementally.

### 3.1 Break down `build_locus_em_data()` in `locus.py`

**File:** `src/rigel/locus.py` lines 126–381 (~255 lines)  
**Issue:** Single function with 8+ logical sections: nRNA deduplication, candidate array building, gDNA candidate construction, CSR layout, etc.  
**Action:** Extract 2–3 focused helpers along existing section boundaries (the function already has section-style comments). Candidates: nRNA span deduplication, gDNA candidate scoring, CSR assembly.

### 3.2 Break down `get_nrna_counts_df()` in `estimator.py`

**File:** `src/rigel/estimator.py` lines 509–631 (~122 lines)  
**Issue:** Dense DataFrame assembly with children-lookup loop and multiple aggregation steps.  
**Action:** Extract the per-entity children aggregation into a helper.

### 3.3 Break down `create_nrna_transcripts()` in `index.py`

**File:** `src/rigel/index.py` lines 156–256 (~100 lines)  
**Issue:** Dense vectorized NumPy operations for nRNA merging with nested loops and multiple intermediate data structures.  
**Action:** Split into 3 focused helpers: (1) per-transcript span computation, (2) span merging/deduplication, (3) synthetic transcript creation.

---

## Priority 4 — C++ Cleanup (Lower Risk)

These require recompilation and testing. Each change is small and isolated.

### 4.1 Consolidate `compute_frag_lengths()` overloads in `resolve_context.h`

**File:** `src/rigel/native/resolve_context.h` lines 284–307  
**Issue:** Two overloads with near-identical logic, differing only in whether a `ResolverScratch&` is passed. The backward-compatible wrapper (line 289) uses internal scratch.  
**Action:** Collapse into a single implementation, defaulting to internal scratch if none is provided.

### 4.2 Evaluate `SplitMix64` PRNG in `em_solver.cpp`

**File:** `src/rigel/native/em_solver.cpp` around line 101  
**Issue:** `SplitMix64` struct is defined and used only in `ASSIGN_SAMPLE` assignment mode (~line 2177), which is typically disabled in production (`ASSIGN_FRACTIONAL` is the default).  
**Action:** No immediate removal needed — this is intentional for optional sampling. Add a brief comment noting it is exercised only in sampling mode.

---

## Out of Scope (Not Recommended)

The following were evaluated and intentionally excluded:

- **`gc.collect()` calls in `partition.py`**: These are deliberate memory-management optimizations for mega-loci. Profiling confirms they help.
- **`_cache` dict pattern in `locus.py`**: Functional and tested. Refactoring to a TypedDict would be cosmetic and risks breaking the pipeline.
- **C++ function complexity** (`score_chunk_impl` at 380 lines, `batch_locus_em_partitioned` at 300 lines): These are performance-critical hot paths where inlining decisions matter. Splitting would risk performance regressions without measurable maintainability gain.
- **`finalize()` method in `buffer.py`**: Used in tests. Removing it would require test rewrites with no production benefit.
- **`AnnotationTable.locus_id` deferred initialization**: By design — locus IDs aren't known until after EM completes. The current pattern is correct.

---

## Execution Order

1. **P1 items** (dead code removal) — zero risk, immediate cleanup
2. **P2 items** (consolidation) — low risk, reduces duplication
3. **P3 items** (simplification) — moderate effort, improves readability
4. **P4 items** (C++ cleanup) — requires recompile + test cycle

Run `pytest tests/ -v` after each change to verify no regressions.
