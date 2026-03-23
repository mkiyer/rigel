# Phase C: Remove nRNA-Specific Machinery — Detailed Implementation Plan

## Overview

Phase B (completed) integrated synthetic nRNA transcripts into the index as regular transcripts. Phase C removes the legacy nRNA shadow system: the separate nRNA index space, nRNA scoring path, nRNA fan-out in union-find, and nRNA-specific EM components. After Phase C, the EM layout simplifies from `[mRNA(0..n_t), nRNA(n_t..n_t+n_nrna), gDNA]` to `[transcript(0..n_t), gDNA(n_t)]`.

**All tests must pass after each sub-phase.** Recompile (`pip install --no-build-isolation -e .`) after every C++ change.

---

## Dependencies and Ordering

The sub-phases are ordered to maintain a compilable, testable state at each step. Independence is noted where sub-phases can be parallelized.

```
C.1 (scoring.cpp) ─────┐
                        ├── C.5 (estimator.py + pipeline.py) ── C.7 (output/annotate)
C.2 (em_solver.cpp) ────┘                                              │
C.3 (scoring.py) ──── depends on C.1                                   │
C.4 (locus.py) ────── depends on C.2                                   │
C.6 (scan.py + scored_fragments.py) ── depends on C.1, C.5             │
C.8 (index.py cleanup) ── depends on C.2, C.5                          │
C.9 (tests + golden) ──── depends on all above ────────────────────────┘
```

**Recommended execution order** (minimize recompiles, catch issues early):
1. C.1 + C.2 together (both C++, one recompile)
2. C.3 + C.4 + C.5 + C.6 together (all Python, no recompile)
3. C.7 (output/annotate adjustments)
4. C.8 (index cleanup)
5. C.9 (tests + golden regeneration)

---

## C.1 — scoring.cpp: Remove nRNA Scoring Path

**Goal**: Eliminate the ~270 lines of nRNA-specific scoring code. Synthetic nRNA transcripts score as regular mRNA candidates (they ARE single-exon transcripts in the index).

### What to Remove

1. **Member variables** (constructor + class body):
   - `nrna_base_` (int32_t)
   - `t_to_nrna_` (vector<int32_t>)
   - `n_nrna_` (int32_t)
   - `nrna_span_` (vector<int32_t>)
   - `nrna_start_` (vector<int32_t>)

2. **Constructor parameters** (~5 params):
   - `nrna_base`, `t_to_nrna_arr`, `nrna_span_arr`, `nrna_start_arr`
   - Constructor body: remove array copy loops for these arrays

3. **MergedNrna struct** (lines ~446-451): Data structure for nRNA multimapper dedup.

4. **Multimapper nRNA block** in `flush_mm_group()` (lines ~616-755):
   - nRNA dedup hash map construction
   - nRNA scoring loop (overhang, fragment length, coverage weight, strand)
   - nRNA hard overhang gate selection
   - nRNA candidate emission with `nrna_base_ + ni` indexing
   - ~140 lines

5. **Singlemapper nRNA block** in `score_chunk_impl()` (lines ~1055-1184):
   - `NrnaMerged` struct and `nrna_map` construction
   - nRNA scoring loop (same computations as multimapper path)
   - nRNA hard overhang gate selection
   - nRNA candidate emission with `nrna_base_ + ni` indexing
   - ~130 lines

6. **Nanobind bindings**: Remove nRNA-related constructor args from the `NativeFragmentScorer` class binding.

### What to Keep

- `POOL_CODE_NRNA` constant and its `m.attr()` export — needed for annotate.py pool labeling.
- `case POOL_CODE_NRNA: return "nRNA";` in bam_scanner.cpp — output label.

### Verification

After removing, synthetic nRNA transcripts will ONLY score through the existing mRNA path. Since they are single-exon transcripts in the index, `exon_bp = span` for all fragments overlapping them. This produces the same scoring behavior as the old nRNA path (verified in bench analysis).

---

## C.2 — em_solver.cpp: Remove nRNA Section from EM

**Goal**: Simplify the EM from `2*n_t + 1` to `n_t + 1` components. Remove all nRNA-specific indexing, fan-out, and accumulation.

### What to Remove

1. **`extract_locus_sub_problem()`** — nRNA section (~80 lines):
   - Remove `nrna_base` parameter
   - Remove all nRNA local mapping: `unique_nrna`, `gnrna_to_local`, `local_t_to_local_nrna`, `nrna_to_t CSR` build
   - Remove nRNA entries from `local_map` construction
   - Simplify: `gdna_idx = n_t`, `n_components = n_t + 1`
   - Remove `LocusSub` fields: `n_nrna`, `local_t_to_local_nrna`, `nrna_to_t_off`, `nrna_to_t_idx`, `local_to_global_nrna`

2. **nRNA bias profiles** (lines ~1386-1398):
   - Remove `bias_profiles[n_t + n] = genomic_span` for nRNA entities
   - Synthetic transcripts already have their own bias profiles in the mRNA section (indices 0..n_t)

3. **All-single-exon prior zeroing** (lines ~1410-1424):
   - Remove entirely — not needed, no dummy entities created

4. **`assign_locus_posteriors()`** — nRNA scatter (~15 lines):
   - Remove `nrna_em_counts` parameter
   - Remove `nrna_total` parameter
   - Remove `N_NRNA_TOTAL` parameter
   - Remove the branch that scatters counts from nRNA components to `nrna_em_counts[]`
   - All RNA counts go to `em_counts_2d[global_t, col]`

5. **`batch_locus_em()`** — nRNA parameters (~20 lines):
   - Remove `nrna_base_index` parameter
   - Remove `t_to_nrna_arr` parameter
   - Remove `num_nrna_total` parameter
   - Remove `nrna_em_counts_out` parameter
   - Remove `locus_nrna_vec`/`locus_nrna_data` accumulators
   - Simplify return type (remove nRNA counts from tuple)

6. **`connected_components_native()`** — nRNA fan-out (~30 lines):
   - Remove `nrna_base` parameter
   - Remove `nrna_to_t_offsets_arr`, `nrna_to_t_indices_arr` parameters
   - **Remove the fan-out block** (lines ~2166-2185): `if (t >= nrna_base)` → fan out to all sharing transcripts. This is the primary driver of mega-loci.
   - Remove nRNA→first transcript mapping in unit-to-component assignment (lines ~2237-2242)

7. **Nanobind bindings**: Remove all nRNA-related args from `run_locus_em_native`, `batch_locus_em`, `connected_components_native`.

### Key Design Decision

After Phase C, the `gdna_idx` is simply `n_t` (the last component). All transcript-level counts (including synthetic nRNA transcripts) go directly into `em_counts_2d[global_t, col]`. The Python layer extracts nRNA counts by filtering on `is_synthetic_nrna`.

### Verification

The EM behavior changes subtly:
- Synthetic transcripts compete as regular EM components (correct)
- No forced fan-out in union-find → smaller loci (beneficial)
- gDNA component is at index `n_t` instead of `n_t + n_nrna` (internal change, transparent to output)

---

## C.3 — scoring.py: Remove nRNA Parameters

**Goal**: Clean up the Python scorer interface.

### What to Remove

1. `FragmentScorer` dataclass: remove `nrna_base: int` field
2. `from_models()`: remove `nrna_base=estimator.nrna_base_index`
3. Native constructor call: remove `nrna_base`, `t_to_nrna_arr`, `nrna_span_arr`, `nrna_start_arr` arguments

**Depends on**: C.1 (scoring.cpp nRNA params removed first)

---

## C.4 — locus.py: Remove nRNA CSR and Fan-out

**Goal**: Simplify locus building — remove nRNA→transcript mapping and simplify EM data construction.

### What to Remove

1. **`build_loci()`** (~22 lines):
   - Remove `nrna_base = em_data.nrna_base_index`
   - Remove nRNA→transcript CSR construction: `t_to_nrna`, `nrna_counts`, `nrna_to_t_offsets`, `nrna_to_t_indices`
   - Remove nRNA args to `_cc_native()`(or `connected_components`)

2. **`build_locus_em_data()`** (~40 lines):
   - Remove `nrna_base = em_data.nrna_base_index`
   - Remove `global_t_to_nrna`, `unique_global_nrna`, `local_t_to_local_nrna`, `nrna_to_t CSR` construction
   - Simplify: `gdna_idx = n_t`, `n_components = n_t + 1`
   - Remove nRNA portion of `local_map` construction
   - Remove nRNA bias profile computation
   - Remove all-single-exon prior zeroing
   - Simplify `LocusEMInput` construction (no nRNA fields)

3. **Component layout**: Update comments from `[mRNA, nRNA, gDNA]` to `[transcript, gDNA]`

**Depends on**: C.2 (em_solver.cpp signature changes)

---

## C.5 — estimator.py + pipeline.py: Remove nRNA Fan-out

**Goal**: Eliminate nRNA-specific accumulation and post-EM fan-out.

### estimator.py — What to Remove

1. Constructor: remove `num_nrna`, `t_to_nrna` parameters and stored attributes
2. Remove `self.nrna_base_index = num_transcripts`
3. Remove `self.nrna_em_counts = np.zeros(...)` array
4. Remove `_nrna_per_transcript()` method entirely

### estimator.py — What to Modify

1. `nrna_em_count` property: Instead of `self.nrna_em_counts.sum()`, compute nRNA total from `em_counts` by summing rows for synthetic transcripts:
   ```python
   @property
   def nrna_em_count(self) -> float:
       if self._synthetic_mask is None:
           return 0.0
       return float(self.em_counts[self._synthetic_mask].sum())
   ```
   Where `_synthetic_mask` is a boolean array set during construction from `is_synthetic_nrna` flags.

2. `get_counts_df()`: For the "nrna" column, use synthetic transcript mRNA counts directly (they ARE nRNA):
   - For synthetic transcripts: their mRNA count IS their nRNA count
   - For annotated transcripts: nRNA = 0 (no separate nRNA channel)
   - Or: rethink the output schema — synthetic transcripts could have a "pool" column indicating nRNA

3. `run_batch_locus_em()`: Remove nRNA-specific parameters from `_batch_locus_em` call. Update return tuple unpacking (no `locus_nrna`).

### pipeline.py — What to Modify

1. `AbundanceEstimator()` construction: remove `num_nrna`, `t_to_nrna` args
2. `_build_locus_meta()`: extract nRNA from em_counts for synthetic transcripts
3. Return tuple unpacking: adapt to new `batch_locus_em` return (no `locus_nrna_arr`)

---

## C.6 — scan.py + scored_fragments.py: Remove nRNA Base Index

**Goal**: Clean up data classes.

### scan.py
- Remove `nrna_base_index=self.ctx.nrna_base` from `ScoredFragments` construction

### scored_fragments.py
- Remove `ScoredFragments.nrna_base_index` field
- Remove `LocusEMInput.n_nrna` field
- Remove `LocusEMInput.local_to_global_nrna` field
- Remove `LocusEMInput.local_t_to_local_nrna` field
- Remove `LocusEMInput.nrna_to_t_offsets` field
- Remove `LocusEMInput.nrna_to_t_indices` field
- Update component layout comments

---

## C.7 — Output and Annotate Adjustments

**Goal**: Ensure output formats correctly label synthetic nRNA transcripts.

### annotate.py
- **Keep** `POOL_CODE_NRNA` and pool labeling. Fragments assigned to synthetic nRNA transcripts should be labeled as nRNA pool in BAM annotations.
- **Modify** the pool assignment logic: instead of checking `component >= nrna_base`, check whether the assigned transcript has `is_synthetic_nrna = True`.

### Output DataFrames
- In `get_counts_df()`: the "nrna" column for annotated transcripts should show the count assigned to the corresponding synthetic nRNA transcript(s). This requires a mapping from annotated transcripts to their associated synthetic.
- Design decision: In Phase C, keep it simple — the "nrna" column shows 0 for annotated transcripts and the mRNA count for synthetics. Gene-level aggregation sums synthetics' mRNA counts as gene nRNA. This matches the conceptual model: synthetic transcripts represent nascent RNA.

---

## C.8 — index.py Cleanup

**Goal**: Remove the nRNA table infrastructure.

### What to Remove

1. `NRNA_FEATHER`, `NRNA_TSV` constants
2. `compute_nrna_table()` function entirely (~45 lines)
3. `TranscriptIndex.__init__`: remove `nrna_df`, `t_to_nrna_arr`, `num_nrna` attributes
4. `build()`: remove `compute_nrna_table()` call and nrna.feather/tsv writes
5. `load()`: remove nrna.feather loading and attribute population

### What to Keep

- `create_nrna_transcripts()` — Phase B infrastructure, still needed
- `_cluster_coordinates()` — helper for TSS/TES clustering
- `NRNA_MERGE_TOLERANCE` — configurable tolerance
- `is_synthetic_nrna` / `is_nascent_equiv` columns in transcripts.feather

### transcript.py

- Remove `nrna_idx: int = -1` field
- Remove `nrna_idx` from `to_dict()`
- Keep `is_synthetic_nrna`, `is_nascent_equiv`, `nrna_abundance`

---

## C.9 — Tests and Golden Outputs

**Goal**: Update all tests that reference nRNA-specific machinery.

### Test files requiring updates

| File | Changes |
|------|---------|
| `conftest.py` | Remove mock `nrna_base`, `num_nrna`, `t_to_nrna_arr` from estimator/index mocks |
| `test_em_impl.py` | Remove `t_to_nrna` params from direct C++ EM calls; simplify EM component expectations |
| `test_estimator.py` | Remove `nrna_base`, `nrna_em_counts` assertions; update `nrna_em_count` property tests |
| `test_cross_chunk.py` | Update `nrna_em_count` assertions |
| `test_gdna.py` | Update `nrna_em_counts.sum()` in total assertions |
| `test_pipeline_routing.py` | Remove mock `nrna_df`, `nrna_base_index` checks |
| `test_gdna_fl_scoring.py` | Remove `num_nrna`, `t_to_nrna` from estimator construction |
| `test_golden_output.py` | Regenerate goldens with `--update-golden` |
| `test_opp_strand_order.py` | Update `nrna_em_counts[ti]` references |
| `test_equiv_class_sort_ab.py` | Update `nrna_em_counts[ti]` references |
| `tests/golden/*_scalars.json` | Update `nrna_em_count` field values |

### Strategy

1. Fix test infrastructure (conftest.py) first
2. Fix C++ binding tests (test_em_impl.py) — these call C++ directly
3. Fix Python unit tests (estimator, cross_chunk, gdna, etc.)
4. Regenerate golden outputs last
5. Full suite validation

---

## Estimated Scope

| Category | Lines removed | Lines modified |
|----------|--------------|----------------|
| C++ (scoring.cpp) | ~270 | ~10 |
| C++ (em_solver.cpp) | ~250 | ~30 |
| Python (locus.py) | ~60 | ~15 |
| Python (estimator.py) | ~30 | ~40 |
| Python (scoring.py, scan.py, scored_fragments.py) | ~20 | ~5 |
| Python (pipeline.py, index.py) | ~60 | ~15 |
| Tests | ~50 | ~100 |
| **Total** | **~740** | **~215** |

---

## Risk Assessment

### Low Risk
- C.1 (scoring.cpp): Clean removal, well-isolated code blocks
- C.3, C.6: Parameter removal, straightforward
- C.8 (index cleanup): Removing dead infrastructure

### Medium Risk
- C.2 (em_solver.cpp): Core EM changes, multiple interacting pieces. Must maintain correctness of:
  - Component indexing (gdna_idx shift)
  - Prior computation
  - Posterior scatter
  - Connected components (no fan-out)
- C.5 (estimator.py): Output format changes. Must ensure:
  - TPM normalization still correct
  - Gene-level aggregation handles synthetics
  - nRNA counts correctly extracted from em_counts

### Mitigations
- Recompile + test after every sub-phase
- Keep benchmark configs from Phase B as regression tests
- Golden output comparison catches any count drift
- Scenario tests validate end-to-end accountability

---

## Post-Phase C: Phase D (Prior for Synthetic nRNAs)

After Phase C, synthetic nRNA transcripts compete as regular transcripts in the EM. In the absence of real nascent RNA signal, they may absorb fragments due to their large single-exon geometry (higher coverage weight for intronic fragments). Options:

1. **No special treatment**: Let the OVR prior handle it (coverage-weighted prior already penalizes large-span components)
2. **Flag-based prior penalty**: Apply a configurable sparsity penalty to `is_synthetic_nrna` transcripts
3. **Data-driven**: Use benchmark results to determine if a prior is needed

Decision deferred to Phase D benchmarking.
