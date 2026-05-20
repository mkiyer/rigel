# RNA Bayesian Prior Cleanup Plan v3

**Date:** 2026-05-13
**Status:** Production-path-first rewrite of [rna_prior_cleanup_v2.md](rna_prior_cleanup_v2.md).
**Principle:** If a code path exists only because old tests call it, delete or rework both. Tests should protect the product Rigel actually ships, not preserve archaeological APIs.

This plan supersedes v2 in one important way: it stops treating legacy-native compatibility as something to preserve. The production quantification path is:

```text
quant_from_buffer -> _run_locus_em_partitioned -> AbundanceEstimator.run_batch_locus_em_partitioned -> _em_impl.batch_locus_em_partitioned
```

The Python-visible native helper `run_locus_em_native` is not on that path. It is used by tests and historical docs only. It should be removed, and the tests that currently exercise it should be rewritten to call the production batch interface with small synthetic partitions.

---

## 0. What v3 Changes Relative To v2

1. **Delete legacy-only native API.**
   Remove `_em_impl.run_locus_em_native` from C++ and nanobind. Do not rename it, preserve a compatibility wrapper, or keep tests just to exercise it.

2. **Make one coherent cleanup PR for the dead RNA prior system.**
   `prior_weight_rna`, `nrna_weight`, `alpha_rna`, and `rna_prior_count` are one conceptual mistake now: a dead RNA-prior allocation surface. Remove them together so no intermediate codebase still tells readers that an RNA prior budget exists.

3. **Do not preserve MAP behavior solely for legacy tests.**
   If MAP has no modeled RNA prior, say so in code. Numerical floors belong inside log/normalization guards, not in the prior model as fake RNA pseudocounts. Default VBEM remains the production recommendation and should remain numerically unchanged.

4. **Tests must target production contracts.**
   Replace single-locus native tests with production-batch tests: small `partition_tuples`, real `AbundanceEstimator.run_batch_locus_em_partitioned`, real gDNA eligibility, real effective lengths, real assignment behavior.

5. **Favor fewer names and fewer layers.**
   Use one Python-side name for the calibration prior mass: `gdna_prior_count`. Use local `alpha` only for actual Dirichlet vectors inside the native solver.

---

## 1. Target Shape

After cleanup, the production prior model is:

1. **gDNA:** eligible gDNA component receives calibration mass: `prior[gdna] = baseline + max(gdna_prior_count, EM_LOG_EPSILON)`.
2. **RNA:** RNA components receive only the mode baseline: `0.5` in VBEM, `0.0` in MAP.
3. **Eligibility:** eligibility is a separate `0/1` component mask; ineligible components get prior `0.0` and no unambiguous total.

The solver can still use numerical guards such as `log(theta + EM_LOG_EPSILON)` and `digamma(max(alpha, EM_LOG_EPSILON))`. Those are stability guards, not prior mass. Keep that distinction visible in names and comments.

The native helper should become a private implementation detail, not a Python API:

```cpp
static void compute_gdna_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,
    double        gdna_prior_count,
    int           gdna_idx,
    double*       prior_out,
    double*       theta_init_out,
    int           n_components,
    bool          use_vbem)
{
    std::copy(unambig_totals, unambig_totals + n_components, theta_init_out);

    for (const auto& ec : ec_data) {
        const int n = ec.n;
        const int k = ec.k;
        const int32_t* cidx = ec.comp_idx.data();
        const double* wt = ec.wt_flat.data();
        for (int i = 0; i < n; ++i) {
            double row_sum = 0.0;
            for (int j = 0; j < k; ++j) {
                row_sum += wt[i * k + j] * eligible[cidx[j]];
            }
            if (row_sum == 0.0) row_sum = 1.0;
            const double inv_row_sum = 1.0 / row_sum;
            for (int j = 0; j < k; ++j) {
                theta_init_out[cidx[j]] +=
                    (wt[i * k + j] * eligible[cidx[j]]) * inv_row_sum;
            }
        }
    }

    const double baseline = use_vbem ? 0.5 : 0.0;
    for (int i = 0; i < n_components; ++i) {
        if (eligible[i] <= 0.0) {
            prior_out[i] = 0.0;
        } else if (i == gdna_idx) {
            prior_out[i] = baseline + std::max(gdna_prior_count, EM_LOG_EPSILON);
        } else {
            prior_out[i] = baseline;
        }
    }
}
```

This is cleaner than v2 because it does not smuggle a legacy MAP RNA floor into the prior vector. If any production MAP tests change, update them as an intentional MAP cleanup. The default VBEM path is the path that must remain stable.

---

## 2. PR-A: Excise The Dead RNA Prior System

**Goal:** Remove `prior_weight_rna`, `nrna_weight`, `alpha_rna`, `rna_prior_count`, and the Python-visible `run_locus_em_native` API in one coherent change.

**Expected numeric impact:** no change for default VBEM production quantification. MAP may change where it previously received `EM_LOG_EPSILON` as fake RNA prior mass; that is intended cleanup, not a regression.

### Native EM

| File | Change |
|---|---|
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Remove `run_locus_em_native` entirely: function, nanobind binding, module docstring references, and comments. |
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Rename `compute_ovr_prior_and_warm_start` to `compute_gdna_prior_and_warm_start` while changing the signature. It should accept only `gdna_prior_count`, `gdna_idx`, eligibility, and mode. |
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Remove `alpha_rna`, `prior_weight_rna`, `coverage_totals`, `total_rna_coverage`, `n_rna_eligible`, and the gDNA warm-start ratio override. |
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Remove `PartitionView::prior_weight_rna`, `pwr_keepalive`, and the `locus_prior_weight_rna` nanobind argument. |
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Remove `locus_alpha_rna` from `batch_locus_em_partitioned`; use only `locus_gdna_prior_count` plus `locus_enable_gdna`. |

Native design notes:

- `build_equiv_classes`, `run_squarem`, `assign_posteriors`, and the batch scheduling logic stay.
- The batch interface is the only Python-visible EM entry point.
- If future tests need to inspect internal priors, expose a small production-adjacent diagnostic in `locus_stats` rather than reintroducing a standalone native helper.

### Python Prior Assembly And Pipeline

| File | Change |
|---|---|
| [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py) | Delete `build_prior_weight_rna`. Remove `alpha_rna`, `rna_prior_count`, and `prior_weight_rna` from `PriorTable`. Remove `rna_prior_count` from `MultiLocusPrior`. Keep `gdna_prior_count` as the single calibration prior mass. |
| [src/rigel/calibration/_result.py](../../src/rigel/calibration/_result.py) | Delete `CalibrationResult.alpha_rna` and `CalibrationResult.prior_weight_rna`. Remove `rna_prior_count` from `_MULTI_LOCUS_COLUMNS`. |
| [src/rigel/calibration/__init__.py](../../src/rigel/calibration/__init__.py) | Remove `build_prior_weight_rna` from imports and `__all__`. |
| [src/rigel/config.py](../../src/rigel/config.py) | Remove `CalibrationConfig.nrna_weight`. |
| [src/rigel/cli.py](../../src/rigel/cli.py) | Remove `--cal-nrna-weight`, `_PARAM_SPECS` mapping, summary/config YAML output for that argument, and help text. |
| [src/rigel/pipeline.py](../../src/rigel/pipeline.py) | Remove `nrna_weight` from `quant_from_buffer`, remove `cal_cfg.nrna_weight` from `run_pipeline`, and pass only `gdna_prior_count` plus `enable_gdna` into `_run_locus_em_partitioned`. |
| [src/rigel/estimator.py](../../src/rigel/estimator.py) | Remove `alpha_rna` and `locus_prior_weight_rna` from `run_batch_locus_em_partitioned`. Remove `alpha_rna` from `locus_results` and `get_loci_df`. Consider renaming user-facing `alpha_gdna` to `gdna_prior_count` in the same PR for clarity. |
| [scripts/profiling/profiler.py](../../scripts/profiling/profiler.py) | Mirror the production pipeline change. Profiling scripts are tooling, but they should follow production signatures. |
| [scripts/debug/](../../scripts/debug) | Delete stale RNA-prior debug scripts or update them to the new gDNA-only model. Do not keep dead fields to satisfy debug output. |

Simulator note: `src/rigel/sim/reads.py` uses local `nrna_weight` variables for simulated abundance weighting. Rename them to `nrna_abundance_weight` if grep clarity matters, but do not confuse them with calibration config.

### Tests To Delete

Delete tests whose target is the removed API or removed knob:

- [tests/test_em_prior_weight.py](../../tests/test_em_prior_weight.py)
- [tests/test_prior_weight_rna_policy.py](../../tests/test_prior_weight_rna_policy.py)
- The `build_prior_weight_rna` / `nrna_weight` block in [tests/test_code_review_fixes_2026_05.py](../../tests/test_code_review_fixes_2026_05.py)
- The `nrna_weight` flow test in [tests/test_pipeline_integration_v6.py](../../tests/test_pipeline_integration_v6.py)
- Any `alpha_rna is zero` acceptance test whose only purpose is preserving a dead field

### Tests To Rebuild Around Production Paths

Create or rewrite tests so they exercise `AbundanceEstimator.run_batch_locus_em_partitioned` or full `run_pipeline`, not `run_locus_em_native`.

Recommended replacements:

| Test focus | Production-path replacement |
|---|---|
| Native EM symmetric split | Small one-locus `partition_tuples` through `run_batch_locus_em_partitioned`; assert RNA/gDNA totals and assignment conservation. |
| VBEM Jeffreys baseline | Batch EM with equal RNA components under `mode="vbem"`; assert no component collapses in a symmetric locus. |
| MAP mode behavior | Batch EM with `mode="map"`; assert the warning is emitted and no fake RNA prior columns exist. Do not assert legacy epsilon floors. |
| gDNA eligibility | Keep [tests/test_native_gdna_eligibility.py](../../tests/test_native_gdna_eligibility.py), but remove any case that only exists for positive `alpha_rna` warm-start override. |
| Prior table schema | [tests/test_assemble_priors.py](../../tests/test_assemble_priors.py), [tests/test_calibration_result.py](../../tests/test_calibration_result.py), and [tests/test_calibrate_orchestrator.py](../../tests/test_calibrate_orchestrator.py) should assert the new compact schema: `gdna_prior_count`, `enable_gdna`, diagnostics. |
| End-to-end stability | [tests/test_golden_output.py](../../tests/test_golden_output.py) remains the main numeric guard for production defaults. |

`tests/test_em_impl.py` should either be deleted outright or rewritten as `tests/test_batch_em_impl.py` with production-batch fixtures. The latter is preferred if it still covers equivalence-class determinism, SQUAREM convergence, effective length correction, and VBEM behavior.

### PR-A Acceptance

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_golden_output.py -v
pytest tests/ -v
ruff check src/ tests/
```

String checks:

```bash
rg "run_locus_em_native|prior_weight_rna|build_prior_weight_rna|alpha_rna|rna_prior_count|cal_nrna_weight" src/rigel tests scripts
```

Expected result: no live-code/test matches. Historical docs may still match only after PR-B marks them historical.

---

## 3. PR-B: Rename Outputs And Publish The Prior Model

**Goal:** Make names match the simpler model.

### Naming

Preferred names after PR-B:

| Current | Preferred | Reason |
|---|---|---|
| `alpha_gdna` in Python schemas | `gdna_prior_count` | This is a calibration-derived count, not one member of a symmetric alpha pair anymore. |
| `gdna_prior` in loci output | `gdna_prior_rate` | It is `gdna_prior_count / n_em_fragments`, not the count itself. |
| `compute_ovr_prior_and_warm_start` | `compute_gdna_prior_and_warm_start` | The OVR allocation is gone. |

If renaming public dataframe columns is too disruptive for one PR, keep `alpha_gdna` as a one-release output alias while making internal Python code use `gdna_prior_count`. Do not keep `alpha_rna` as a matching alias.

### Documentation

- Create `docs/calibration/prior_model.md` as the canonical source of truth.
- Replace live docs in [docs/parameters.md](../parameters.md), [docs/MANUAL.md](../MANUAL.md), and [docs/METHODS.md](../METHODS.md).
- Add a banner to older implementation plans in `docs/calibration/` and `docs/bayesian_prior/`: `Historical design note. Current implementation: docs/calibration/prior_model.md.`
- Update [.github/copilot-instructions.md](../../.github/copilot-instructions.md) and [CLAUDE.md](../../CLAUDE.md) only if they mention removed knobs.

Canonical wording:

> Rigel uses a per-multi-locus asymmetric Dirichlet prior where only the eligible gDNA component receives calibration-derived mass. RNA components receive the mode baseline only: Jeffreys `0.5` in VBEM, zero in MAP. Eligibility is independent of prior strength and is computed from whether the locus contains unspliced units with finite gDNA likelihoods.

### PR-B Acceptance

```bash
rg "OVR prior|prior_weight_rna|alpha_rna|rna_prior_count|nrna_weight|run_locus_em_native" README.md CLAUDE.md .github docs src/rigel
```

Expected result: no live-doc/source matches, except historical docs with the banner and cleanup plans.

---

## 4. PR-C: MAP/VBEM Mode Discipline

**Goal:** Make MAP explicitly research/debug and keep VBEM as the default production mode.

Implementation:

- In [src/rigel/config.py](../../src/rigel/config.py), annotate `EMConfig.mode` as `Literal["vbem", "map"]`; default remains `"vbem"`.
- In [src/rigel/cli.py](../../src/rigel/cli.py), make `--em-mode` use argparse choices `vbem` and `map`.
- Emit one warning per run when `mode == "map"`: `EM mode 'map' selected: RNA components receive no modeled prior mass; 'vbem' is recommended for production quantification.`
- Test the warning through an existing pipeline or CLI path, not by calling a removed native helper.

Acceptance:

- Golden outputs unchanged under defaults.
- Config and CLI tests prove invalid modes are rejected.
- MAP warning appears once per run.

---

## 5. PR-D: `coverage_weights` Audit

**Goal:** Decide whether the remaining `coverage_weights` array is worth carrying. This is independent of the dead RNA-prior cleanup once PR-A has landed.

After PR-A, `coverage_weights` has two reasons to exist:

1. Row-normalized warm-start shares for `theta_init`.
2. Deterministic row ordering when equal log-likelihood rows need a secondary key.

Experiment:

| Config | Behavior | Decision use |
|---|---|---|
| A | Current trapezoid weights | Baseline. |
| B | Trapezoid weights stored as float32 | Keep signal, reduce memory. |
| C | Uniform warm start, no stored array | Remove signal and memory. |

Rules:

- Use a temporary benchmark branch or short-lived compile switch. Do not merge a permanent mode selector unless the final architecture needs it.
- If Config C ships, replace the sort tie-break with a deterministic key independent of weights.
- Benchmark VCAP plus at least three synthetic conditions: no/low/high gDNA at `ss=0.90`, including an `nrna_rand` condition.
- Report Spearman correlation, max relative drift on non-negligible counts, total mRNA/nRNA/gDNA fractions, per-locus SQUAREM iterations, wall time, and peak RSS.

Decision:

- Ship C if drift is scientifically negligible and deterministic tests pass.
- Else ship B if B is indistinguishable from A.
- Else keep A and document why the array is structurally important.

Publish the decision in `docs/performance/`.

---

## 6. Final Definition Of Done

The cleanup is complete when:

1. There is exactly one Python-visible native EM entry point: `batch_locus_em_partitioned`.
2. Live code and tests contain no `run_locus_em_native`, `prior_weight_rna`, `build_prior_weight_rna`, calibration `nrna_weight`, `alpha_rna`, or `rna_prior_count`.
3. Tests that previously called legacy native helpers now exercise production batch EM or full pipeline paths.
4. The production prior model has one calibration mass: `gdna_prior_count`, plus independent gDNA eligibility.
5. Default VBEM golden outputs are unchanged after PR-A through PR-C.
6. MAP behavior is documented as debug/research and tested through production configuration.
7. `coverage_weights` has a benchmark-backed keep/narrow/remove decision.
8. The docs contain one canonical current prior model, and older design docs are clearly historical.

---

## 7. Why This Is Healthier

The elegant codebase is the one where every public knob has a real effect, every public native entry point is used by the product, and every test explains a production guarantee. This cleanup removes a false nRNA-prior control surface, a legacy RNA-only native EM binding, and tests that were keeping both alive. What remains is smaller, easier to read, and more honest about Rigel's actual model: transcript-centric EM with a gDNA calibration prior, not a hidden RNA allocation heuristic.