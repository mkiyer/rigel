# Fragment Length Distribution Simplification Plan

**Date:** 2026-03-17  
**Status:** Draft  
**Prerequisite:** Phase 6 complete (calibration mandatory, uncalibrated code paths deleted)

## Motivation

The fragment length (FL) modeling code carries two coexisting approaches for
building the gDNA and RNA FL distributions:

1. **OLD — Strand deconvolution** (`mix_models()` in `frag_length_model.py`):
   Splits unspliced genic fragments into sense/antisense pools and uses strand
   specificity to estimate gDNA vs. RNA fractions.  Builds gDNA FL from
   `intergenic + certainty-weighted gDNA fraction of unspliced genic` and RNA FL
   from `annotated-spliced + certainty-weighted RNA fraction of unspliced genic`.

2. **NEW — Calibration EM** (`build_gdna_fl_model()` in `calibration.py`):
   Weights per-region FL observations by γ_r (posterior probability of "not
   expressed") from the calibration EM.  Uses ALL regions (genic + intergenic),
   not just intergenic.  This is the superior approach.

With calibration now mandatory, `mix_models()` serves only as a fallback when
the calibration EM's ESS gate rejects the gDNA FL model.  Meanwhile:

- The RNA FL model is contaminated with certainty-weighted unspliced genic
  fragments (via `mix_models()`), instead of being the pure annotated-spliced
  "gold standard."
- Five sub-models exist solely to support `mix_models()` strand deconvolution.
- C++ BAM scanner collects `unspliced_same_strand_lengths` and
  `unspliced_opp_strand_lengths` exclusively for this purpose.
- `observe_intergenic_batch()` writes to `gdna_model` which is immediately
  overwritten by `mix_models()`.
- `global_model` (which includes intergenic fragments) is used for effective
  length computation, where only RNA FL should be used.

## Design Principles

1. **RNA FL = annotated-spliced only.** This is the gold standard.  No mixing.
2. **gDNA FL owned by calibration.** The calibration module provides the gDNA FL
   model using γ-weighted region posteriors.  If the EM-based model fails (low
   ESS), calibration itself provides the fallback — a pure intergenic FL model.
3. **Diagnostic sub-models retained.** Category histograms (unspliced,
   spliced-annotated, spliced-unannotated) and the intergenic histogram are
   retained for logging and JSON diagnostics.  They are NOT used in scoring.
4. **Effective length uses RNA FL.**  The `global_model` is replaced by
   `rna_model` for transcript effective length computation.

## Detailed Changes

### Step 1: Simplify `FragmentLengthModels` — Remove strand mixing

**Files:** `src/rigel/frag_length_model.py`

- **Delete `mix_models()` method** (~100 lines).  This is the core deletion.
- **Delete `unspliced_same_strand` and `unspliced_opp_strand` sub-models.**
  Remove from `__init__`, `finalize()`, `to_dict()`, `write_json()`,
  `log_summary()`.
- **Stop `observe_intergenic_batch()` from writing to `gdna_model`.**  Route to
  `global_model` + `intergenic` only (diagnostic).
- **Set `rna_model = category_models[SPLICED_ANNOT]`.**  The RNA model is now
  simply the SPLICED_ANNOT histogram — no mixing.  In `__init__`, make
  `rna_model` a property alias or set it equal to the SPLICED_ANNOT model after
  training.  Since the model is a mutable reference, we can either:
  - (a) Make `rna_model` a `@property` returning `category_models[SPLICED_ANNOT]`,
    or
  - (b) In a new `build_scoring_models()` method (replacing `mix_models()`),
    simply assign `self.rna_model = self.category_models[SpliceType.SPLICED_ANNOT]`.
  Option (b) is preferred — it avoids hidden aliasing and keeps the explicit
  sequencing of observe → build → finalize.
- **Remove `gdna_model` from the init.** The gDNA model will be injected by the
  calibration module, not built inside `FragmentLengthModels`. Add a method or
  setter to accept the calibrated gDNA model.
- **Keep `global_model`, `intergenic`, and all `category_models` for
  diagnostics.**

The resulting `FragmentLengthModels` holds:

| Sub-model | Purpose |
|-----------|---------|
| `global_model` | Diagnostic summary stats (CLI, JSON, logging) |
| `category_models[UNSPLICED]` | Diagnostic only |
| `category_models[SPLICED_ANNOT]` | **Source of `rna_model`** — used in scoring |
| `category_models[SPLICED_UNANNOT]` | Diagnostic only |
| `intergenic` | Diagnostic; also **fallback source for gDNA FL** in calibration |
| `rna_model` | **Scoring** — alias/copy of SPLICED_ANNOT after build |
| `gdna_model` | **Scoring** — injected from calibration |

### Step 2: Move gDNA FL fallback into calibration

**Files:** `src/rigel/calibration.py`

When the calibration EM produces `gdna_fl_model = None` (ESS below threshold),
the calibration module should provide its own fallback rather than relying on
the pipeline to keep the strand-deconvolved model.  The fallback is a pure
**intergenic FL model** — the simplest, assumption-free gDNA proxy.

**Implementation:**

- Add `intergenic_fl_model` as a parameter to `calibrate_gdna()` (or pass the
  intergenic histogram counts through `fl_table`).
- In the `GDNACalibration` return path, when `final_fl is None` (ESS gate
  rejects the EM-based model), substitute the intergenic FL model:
  ```python
  if final_fl is None and intergenic_fl_model is not None:
      final_fl = intergenic_fl_model
      logger.info("[CAL-FL] Using intergenic FL as fallback (low ESS)")
  ```
- Similarly, in the early-bailout path (not enough eligible regions), use the
  intergenic FL model as fallback.
- **Result:** `GDNACalibration.gdna_fl_model` is ALWAYS non-None when any FL
  observations exist.  The pipeline no longer needs a fallback branch.

### Step 3: Simplify pipeline wiring

**Files:** `src/rigel/pipeline.py`

- **Remove the `mix_models()` call** at L798.  Replace with a new call to set
  up the RNA model:
  ```python
  frag_length_models.build_scoring_models()  # sets rna_model = SPLICED_ANNOT
  ```
- **Remove the CAL-FL replacement block** at L660-674.  Instead, unconditionally
  inject the calibrated gDNA FL model:
  ```python
  frag_length_models.gdna_model = calibration.gdna_fl_model
  ```
  (This is safe because Step 2 guarantees `gdna_fl_model` is always non-None.)
- **Pass `frag_length_models.intergenic` to `calibrate_gdna()`** so the
  calibration module has the intergenic FL model available for its fallback.
- **Remove `_replay_fraglen_observations()` entries for
  `unspliced_same_strand_lengths` and `unspliced_opp_strand_lengths`.**
- **Use `rna_model` for effective length** instead of `global_model`:
  ```python
  effective_lengths = frag_length_models.rna_model.compute_all_transcript_eff_lens(...)
  ```
  Keep `global_model.mean` for the `mean_frag` fallback, or use `rna_model.mean`.

### Step 4: Remove C++ strand-direction FL collection

**Files:** `src/rigel/native/bam_scanner.cpp`

- **Delete the unspliced same/opp strand fragment length collection block**
  (~L1323-1346).  This includes the `result.get_is_same_strand()` check, the
  `ctx.t_strand_arr_` lookup, and the push to
  `unspliced_same_strand_lengths` / `unspliced_opp_strand_lengths`.
- **Delete the vector declarations** (`unspliced_same_strand_lengths`,
  `unspliced_opp_strand_lengths` at ~L215-216).
- **Delete the merge logic** that concatenates these vectors across chunks
  (~L428-435).
- **Delete the nanobind dict entries** that export these vectors (~L1427-1432).

### Step 5: Update tests

**Files:** `tests/test_gdna_frag_length.py`, others

- **Delete or rewrite `tests/test_gdna_frag_length.py`.**  This file is
  entirely dedicated to testing `mix_models()` strand deconvolution.  With
  `mix_models()` deleted, all tests in this file become obsolete.  Replace with
  a small set of tests verifying:
  - `build_scoring_models()` correctly sets `rna_model` = SPLICED_ANNOT.
  - The intergenic fallback in calibration works.
  - `build_gdna_fl_model()` produces correct γ-weighted histograms.
- **Update any tests in `test_calibration.py`** that reference the old gDNA FL
  replacement logic.
- **Update `test_gdna.py`** helpers that call `mix_models()` on
  `FragmentLengthModels`.
- **Update `test_pipeline_routing.py`** if it references `mix_models`.
- **Update golden outputs** (`pytest tests/ --update-golden`).

### Step 6: Update CLI / diagnostics / JSON

**Files:** `src/rigel/cli.py`, `src/rigel/frag_length_model.py`

- **CLI summary stats:** Switch from `global_model` to `rna_model` for
  `frag_length_mean`, `frag_length_median`, etc. — or keep both and label
  clearly.  Keep `global_model` stats too if they're useful for diagnostics.
- **`to_dict()` / `write_json()`:** Remove `unspliced_same_strand` and
  `unspliced_opp_strand` entries.  Keep all other sub-models for diagnostics.
- **`log_summary()`:** Remove strand sub-model logging.

## Execution Order

| Step | Description | Risk | Reversible |
|------|-------------|------|------------|
| 1 | Simplify `FragmentLengthModels` | Medium — scoring model inputs change | Yes |
| 2 | Move gDNA FL fallback into calibration | Low — additive | Yes |
| 3 | Simplify pipeline wiring | Medium — integration point | Yes |
| 4 | Remove C++ strand-direction FL collection | Low — data collection only | Yes (recompile) |
| 5 | Update tests + golden outputs | Low | Yes |
| 6 | Update CLI / diagnostics / JSON | Low — cosmetic | Yes |

Steps 1-3 are tightly coupled and should be done together.  Step 4 (C++) is
independent and can be done in any order relative to 1-3, but requires
recompilation.  Steps 5-6 are cleanup after the core changes.

## Lines of Code Impact (Estimates)

| Action | Lines deleted | Lines added |
|--------|--------------|-------------|
| Delete `mix_models()` | ~100 | 0 |
| Delete strand sub-models | ~30 | 0 |
| `build_scoring_models()` | 0 | ~10 |
| Calibration fallback | 0 | ~15 |
| Pipeline simplification | ~40 | ~10 |
| C++ deletion | ~40 | 0 |
| Test rewrite | ~600 | ~80 |
| **Total** | **~810** | **~115** |

Net reduction: ~700 lines across Python + C++.

## Risks and Mitigations

1. **RNA FL from spliced-only may have fewer observations** in low-depth
   samples.  Mitigation: `FragmentLengthModel` uses Laplace smoothing, so it
   gracefully handles sparse histograms.  The SPLICED_ANNOT category is already
   the dominant category for most RNA-seq libraries.

2. **Intergenic FL fallback may differ from the strand-deconvolved model.**
   Mitigation: For any library with reasonable coverage, the calibration EM will
   produce a γ-weighted model (ESS > 50).  The intergenic fallback only
   activates for extremely low-coverage or unusual libraries, and intergenic
   fragments are a reasonable gDNA proxy.

3. **Effective length change from global → RNA FL.**  This shifts effective
   lengths slightly because intergenic (typically longer) fragments no longer
   inflate the distribution.  The change is correct — effective length should
   reflect the fragment length distribution of fragments that could have come
   from a transcript (i.e., RNA).  This may cause minor golden output changes.

4. **C++ changes require recompilation.**  This is standard per project
   instructions.

## Validation

After implementation:

```bash
conda activate rigel && pip install --no-build-isolation -e .
pytest tests/ -v                 # all tests must pass
pytest tests/ --update-golden    # regenerate golden outputs
ruff check src/ tests/
ruff format src/ tests/
```

Expect golden output changes due to:
- RNA FL model changing (pure spliced vs. mixed)
- Effective length computation changing (RNA FL vs. global FL)
- gDNA FL model potentially changing in low-ESS edge cases
