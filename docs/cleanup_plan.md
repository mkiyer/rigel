# hulkrna Cleanup Plan

Date: 2026-03-06

## Goal

Make the codebase clean, readable, and production-ready without changing
scientific behavior or performance characteristics.

## Principles

1. The frozen dataclasses in `config.py` are already the right abstraction —
   just remove the duplicate copies everywhere else.
2. Every phase must pass all 766+ tests and 21 golden scenarios bit-exact.
3. No new abstraction layers unless they eliminate concrete duplication.
4. Internal numerical constants are fine as module-level names — just name
   and document them. No need for separate "constants modules."

---

## Phase 1 — Config as Single Source of Truth ✅ DONE

**Problem:** The same default values live in three places that drift independently:

- `config.py` dataclass field defaults (authoritative)
- `cli.py` `_DEFAULTS` dict (~35 keys, manually synchronized)
- `test_cli.py` `_SMALL_DEFAULTS` dict (~25 keys, manually synchronized)

Additionally, `scoring.py` defines `DEFAULT_OVERHANG_ALPHA`, `DEFAULT_MISMATCH_ALPHA`,
and `DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT` as module constants, and
`FragmentScoringConfig` stores `None` to mean "use those module defaults."
Meanwhile `estimator.py` defines `_KAPPA_MIN`, `_KAPPA_MAX`, `_KAPPA_FALLBACK`,
`_KAPPA_MIN_OBS` as module constants that duplicate `EMConfig` fields.

**Fix:**

1. **Delete `_DEFAULTS`** from `cli.py`. Rewrite `_resolve_quant_args` to
   pull defaults from config dataclass fields via `dataclasses.fields()`.
2. **Delete `_SMALL_DEFAULTS`** from `test_cli.py`. Tests derive defaults
   from config.py.
3. **Move scoring defaults** into `FragmentScoringConfig` as real field
   defaults instead of `None` + module-level constants. Keep the alpha
   defaults in scoring.py only as the conversion functions.
4. **Delete `_KAPPA_*` module constants** from `estimator.py`. The
   `estimate_kappa()` function should receive these from the config that
   already carries them.
5. **Audit `EB_K_LOCUS` / `EB_K_CHROM`** in `locus.py` — these are
   vestigial from before MoM auto-estimation was added. Verify whether
   they're still used; if not, remove them.

**Files:** config.py, cli.py, test_cli.py, scoring.py, estimator.py, locus.py
**Risk:** Low — pure data plumbing, no algorithm changes

---

## Phase 2 — Slim quant_command() + Classify Constants ✅ DONE

**Problem:** `quant_command()` is ~200 lines mixing path validation, config
construction, index loading, pipeline invocation, and output writing.
Several magic numbers across the codebase are unnamed or undocumented.

**Fix:**

1. Extract from `quant_command()`:
   - `_build_pipeline_config(args, sj_strand_tag) → PipelineConfig`
   - `_write_quant_outputs(result, index, output_dir, args)`
   - `_write_run_config(args, output_dir)`
2. Document or name remaining magic numbers:
   - `_ANNOTATION_TABLE_PADDING=1024` (pipeline.py) — internal allocation
   - `_TAIL_DECAY_LP = log(0.99)` (frag_length_model.py) — internal
   - `0x7FFFFFFF` sentinel (scan.py) → `_OH_SENTINEL`
   - `_MIN_CI_OBSERVATIONS=10` (strand_model.py) — internal
   - Clarify `STRAND_DENOM_MIN=0.2` vs `_STRAND_DENOM_EPS=0.01` (locus.py)

**Files:** cli.py, pipeline.py, scan.py, locus.py, strand_model.py
**Risk:** Low

---

## Phase 3 — Dead Code Removal + Function Decomposition ✅ DONE

**Completed:**

1. **Deleted `pyfallback.py`** (630 lines) and `test_batch_parity.py` (18 tests).
   - Inlined `scan_python` per-fragment iteration into `FragmentRouter._scan_python`
     (used by test mocks with simple `_Chunk` objects).
   - Inlined `run_locus_em` and `assign_locus_ambiguous` into `AbundanceEstimator`
     methods (used by unit tests).
   - Removed `HULK_FORCE_PYTHON` env var dispatch from `pipeline.py`.
   - Batch C++ path is now the sole production EM path.
2. **Extracted `_aggregate_nrna_frac_by_group()`** in estimator.py — replaces
   duplicate 15-line aggregation blocks at TSS-group and locus-strand levels.
3. **Extracted `_compute_chrom_gdna_rates()` and `_compute_per_locus_gdna_rates()`**
   in locus.py — isolate chrome-level and per-locus gDNA rate computation.
4. 748/748 tests passing.

---

## Phase 4 — Naming Cleanup + Parameter Documentation ✅ DONE

**Completed:**

1. **Renamed scan.py abbreviations:** `oh`→`overhang`, `ct`→`count_col`,
   `cov_wt`→`coverage_wt`, `rl`→`read_len`, `ebp`→`exon_bp_val`,
   `tx_s`/`tx_e`→`tx_start`/`tx_end`, `best_ct`→`best_count_col`,
   `min_oh`→`min_overhang` across 6 functions.
2. **Renamed locus.py abbreviations:** `ss`→`strand_spec`,
   `r`/`n`→`rate`/`evidence`, `W`→`inv_var_weight`,
   `w`→`shrink_weight` across 6 functions.
3. **Renamed estimator.py abbreviations:** `s`→`strand_spec`,
   deleted `_z = np.zeros` alias (replaced with direct `np.zeros` calls),
   `W`→`inv_var_weight` across 3 functions.
4. **Documented `SpliceStrandCol`:** Added explicit column formula
   docstring: `col = splice_type * 2 + int(is_antisense)`.
5. **Wrote `docs/parameters.md`:** Comprehensive reference for all 34
   CLI parameters organized by category, plus config dataclass reference.
6. Also renamed `ss`→`strand_spec` in `pipeline.py`.
7. 748/748 tests passing.

---

## What We're NOT Doing (and Why)

| Suggestion from outside review | Why skip |
|---|---|
| New `parameter_registry.py` module | The frozen dataclasses ARE the registry. Deleting the duplicate dicts is the fix. |
| Phase 0 guardrails | 766 tests + 21 golden scenarios already exist. |
| Backend parity tests | C++ paths are tested by the same suite. |
| Generated docs from registry | Manual `parameters.md` is fine at this scale. |
| Separate constants modules | Module-level names with comments work. Extra indirection hurts readability. |
| Architectural glossary | The existing docstrings cover this. |
