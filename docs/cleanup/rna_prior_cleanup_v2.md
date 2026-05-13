# RNA Bayesian Prior Cleanup Plan v2

**Date:** 2026-05-13
**Status:** Reviewed against the live code before implementation.
**Scope:** Remove the dead RNA-prior allocation path feeding the native EM, make the surviving gDNA-only prior model explicit, and defer any `coverage_weights` math change to a benchmarked follow-up.

This v2 plan keeps the good core of [rna_prior_cleanup.md](rna_prior_cleanup.md) but tightens the sequencing and fixes several risky details. The most important change is a sharper distinction between **numeric EM identity** and **schema/API churn**. Removing dead priors should leave transcript, gene, and locus numeric quantities unchanged; removing dead columns, CLI flags, and compatibility properties is still an intentional public surface change and should be tested as such.

---

## 0. Critique of v1

The original plan is directionally right: `prior_weight_rna` and `alpha_rna` are dead in production because calibration pins `alpha_rna` to zero, while `coverage_weights` remains live only for warm start and deterministic ordering. The improvements below are about making the cleanup easier to merge and harder to misunderstand later.

1. **The bit-identity claim is too broad.**
   Quantification math should be unchanged for the production batch path, but deleting `alpha_rna`, `prior_weight_rna`, `nrna_weight`, or `rna_prior_count` changes Python object schemas, config YAML, summary command arguments, and some diagnostic dataframes. Verification should say "numeric outputs unchanged" unless a PR explicitly preserves every output column.

2. **The empty-locus MAP snippet in v1 is not bit-safe.**
   Today the legacy single-locus native path gives each eligible component an `EM_LOG_EPSILON` floor when `alpha_rna == 0`, so an otherwise empty eligible locus normalizes to a uniform `theta`. Replacing the floor with exact zero returns all-zero `theta`. That path is test/legacy rather than production batch EM, but the plan should not call it bit-safe.

3. **`run_locus_em_native` is a legacy/test-only API and needs separate handling.**
   In production, the batch native path consumes per-locus gDNA priors. The single-locus native test helper has no gDNA component and uses `alpha_rna` as a generic RNA prior budget. Removing `alpha_rna` there is a test/API cleanup, not the same operation as removing dead production calibration plumbing.

4. **The test inventory is incomplete.**
   Besides `test_em_prior_weight.py`, `test_prior_weight_rna_policy.py`, and `test_pipeline_integration_v6.py`, the cleanup must update `test_code_review_fixes_2026_05.py`, `test_bayesian_prior_acceptance.py`, `test_assemble_priors.py`, `test_calibration_result.py`, `test_native_gdna_eligibility.py`, `tests/conftest.py`, and scenario/debug tests that print `alpha_rna`.

5. **The script inventory is incomplete.**
   `scripts/profiling/profiler.py` directly mirrors `quant_from_buffer` prior plumbing and must be updated with the production code. Debug scripts under `scripts/debug/` may be retired or made tolerant of missing `alpha_rna` diagnostics.

6. **`scan.py` does not forward `alpha_rna`.**
   It forwards `coverage_weights`. The PR-A table in v1 should not list `scan.py` for `alpha_rna` removal.

7. **Naming should optimize reader clarity, not only internal symmetry.**
   `alpha_gdna` is concise in the native EM, but `gdna_prior_count` is clearer in Python diagnostics and user-facing tables. v2 recommends one semantic name at the Python boundary, then a thin native argument named for the same concept.

8. **The `coverage_weights` experiment should not leave permanent build knobs unless they ship.**
   A temporary benchmark branch can use a compile-time or local runtime toggle. The merged code should contain only the selected design: keep current weights, narrow them, or remove them.

---

## 1. Target Model

After PR-A and PR-B, Rigel's production EM prior model should be explainable in three lines:

1. **gDNA component:** if the component is eligible, `prior[gdna] = vbem_baseline + max(gdna_prior_count, EM_LOG_EPSILON)`, where `vbem_baseline` is `0.5` in VBEM and `0.0` in MAP.
2. **RNA components:** each eligible RNA component gets only the mode baseline: `0.5` in VBEM, and `EM_LOG_EPSILON` in MAP as a numerical floor, not a modeled RNA pseudocount.
3. **Eligibility:** component eligibility is independent of prior strength. Ineligible components get prior `0.0` and have their unambiguous totals zeroed before EM.

In Python-facing names, prefer `gdna_prior_count` for the calibration-derived mass. Use `alpha_gdna` only where the native EM vocabulary already requires it, and do not preserve an `alpha_rna` sibling just for symmetry.

The native helper should end up close to this shape:

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

    const double rna_floor = use_vbem ? 0.5 : EM_LOG_EPSILON;
    const double gdna_baseline = use_vbem ? 0.5 : 0.0;
    for (int i = 0; i < n_components; ++i) {
        if (eligible[i] <= 0.0) {
            prior_out[i] = 0.0;
        } else if (i == gdna_idx) {
            prior_out[i] = gdna_baseline + std::max(gdna_prior_count, EM_LOG_EPSILON);
        } else {
            prior_out[i] = rna_floor;
        }
    }
}
```

This intentionally preserves the MAP epsilon floor for eligible RNA components. It is negligible for real EM updates, but it avoids changing legacy empty-locus behavior while still removing the dead RNA allocation model.

---

## 2. Revised PR Sequence

### PR-A: Remove `prior_weight_rna` and `nrna_weight`

**Goal:** Delete the per-component RNA-prior allocation weight end to end. Numeric EM output should be unchanged because every production consumer multiplied it by `alpha_rna == 0`.

This is the highest-confidence cleanup and captures the per-locus float32 allocation savings without touching the larger `alpha_rna` schema yet.

| Area | Required change |
|---|---|
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Drop `prior_weight_rna` from `compute_ovr_prior_and_warm_start`, `PartitionView`, `batch_locus_em_partitioned`, and `run_locus_em_native`. Keep `alpha_rna` for this PR only, and distribute any positive RNA prior by coverage alone. |
| [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py) | Delete `build_prior_weight_rna`, remove `PriorTable.prior_weight_rna`, remove `assemble_priors(..., nrna_weight=...)`, and stop allocating one float32 vector per locus. |
| [src/rigel/calibration/__init__.py](../../src/rigel/calibration/__init__.py) | Remove the `build_prior_weight_rna` export. |
| [src/rigel/config.py](../../src/rigel/config.py) | Remove `CalibrationConfig.nrna_weight`. |
| [src/rigel/cli.py](../../src/rigel/cli.py) | Remove `--cal-nrna-weight` and its `_PARAM_SPECS` mapping. Config YAML written after this PR should no longer contain `cal_nrna_weight`. |
| [src/rigel/pipeline.py](../../src/rigel/pipeline.py) | Remove the `nrna_weight` argument to `quant_from_buffer`, the `cal_cfg.nrna_weight` call site in `run_pipeline`, and all `prior_weight_rna_per_locus` plumbing. |
| [src/rigel/estimator.py](../../src/rigel/estimator.py) | Remove `locus_prior_weight_rna` from `run_batch_locus_em_partitioned` and the native call. |
| [scripts/profiling/profiler.py](../../scripts/profiling/profiler.py) | Mirror the production pipeline signature changes so profiling remains usable. |
| [src/rigel/sim/reads.py](../../src/rigel/sim/reads.py) | Leave the local `nrna_weight` variables alone or rename them to `nrna_abundance_weight`; they are simulation abundance weights, not calibration prior knobs. |

Tests to delete or rewrite in PR-A:

- Delete [tests/test_em_prior_weight.py](../../tests/test_em_prior_weight.py).
- Delete [tests/test_prior_weight_rna_policy.py](../../tests/test_prior_weight_rna_policy.py).
- Remove the `nrna_weight` flow test from [tests/test_pipeline_integration_v6.py](../../tests/test_pipeline_integration_v6.py).
- Rewrite `test_code_review_fixes_2026_05.py` by deleting the `build_prior_weight_rna` and `nrna_weight` regression block; keep unrelated fixes.
- Update `test_calibrate_orchestrator.py`, `test_calibration_result.py`, and `test_assemble_priors.py` fixtures that still construct `PriorTable(prior_weight_rna=...)`.

PR-A verification:

```bash
conda activate rigel && pip install --no-build-isolation -e .
conda activate rigel && pytest tests/test_golden_output.py -v
conda activate rigel && pytest tests/ -v
conda activate rigel && ruff check src/ tests/
```

Acceptance checks:

- `rg "prior_weight_rna|build_prior_weight_rna|cal_nrna_weight|CalibrationConfig\(.*nrna_weight" src tests scripts` returns no production/test references.
- `nrna_weight` appears only in simulator abundance contexts, or those locals have been renamed.
- Transcript, gene, locus numeric golden comparisons pass without `--update-golden`.

### PR-B: Remove `alpha_rna` and `rna_prior_count`

**Goal:** Delete the dead RNA prior budget from the production native ABI and Python prior schema. Numeric EM output should remain unchanged for production because the deleted value is always zero.

This PR is where v1's "bit-safe" argument applies, but only with the MAP epsilon caveat from section 1.

| Area | Required change |
|---|---|
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Remove `alpha_rna` from `compute_ovr_prior_and_warm_start`, `batch_locus_em_partitioned`, and the nanobind signature. Delete `coverage_totals` and the gDNA warm-start ratio override. Preserve the `EM_LOG_EPSILON` RNA floor in MAP. |
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Decide explicitly what to do with `run_locus_em_native`: either remove the `alpha_rna` argument and make it a baseline-only test helper, or keep a renamed private `rna_prior_budget` only if tests truly need a non-production RNA-only prior. Prefer removal for clarity. |
| [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py) | Remove `PriorTable.alpha_rna`, `PriorTable.rna_prior_count`, `MultiLocusPrior.rna_prior_count`, and the `eta_r = 0.0` allocation. |
| [src/rigel/calibration/_result.py](../../src/rigel/calibration/_result.py) | Remove `CalibrationResult.alpha_rna`; remove `rna_prior_count` from `_MULTI_LOCUS_COLUMNS` unless a one-release compatibility column is deliberately retained. |
| [src/rigel/pipeline.py](../../src/rigel/pipeline.py) | Pass only `gdna_prior_count`/`alpha_gdna` and `enable_gdna` into batch EM. Remove `alpha_rna` from `locus_results`. |
| [src/rigel/estimator.py](../../src/rigel/estimator.py) | Remove `alpha_rna` from `run_batch_locus_em_partitioned`, `locus_results`, and `get_loci_df`. Consider renaming the user-facing `alpha_gdna` column to `gdna_prior_count` in this PR or PR-C. |
| [tests/conftest.py](../../tests/conftest.py) | Update batch EM fixtures to return `(em_data, loci, gdna_prior_count, index)` instead of carrying `alpha_rna`. |
| [scripts/debug/](../../scripts/debug) | Retire or update scripts that report `alpha_rna_sum` or `locus_alpha_rna`. Debug scripts are not a reason to keep dead production fields. |

Tests to update in PR-B:

- [tests/test_bayesian_prior_acceptance.py](../../tests/test_bayesian_prior_acceptance.py): replace "`alpha_rna` is always zero" with "no RNA prior budget exists" and assert the gDNA prior count plus eligibility fields instead.
- [tests/test_assemble_priors.py](../../tests/test_assemble_priors.py): remove `rna_prior_count` assertions and fixtures.
- [tests/test_calibration_result.py](../../tests/test_calibration_result.py): update locked dataframe columns and alias properties.
- [tests/test_calibrate_orchestrator.py](../../tests/test_calibrate_orchestrator.py): update `PriorTable.empty()` and `with_priors` expectations.
- [tests/test_native_gdna_eligibility.py](../../tests/test_native_gdna_eligibility.py): delete cases whose only purpose is the positive-`alpha_rna` warm-start override; keep coverage for gDNA eligibility independent of `gdna_prior_count`.
- [tests/test_em_impl.py](../../tests/test_em_impl.py): remove OVR/RNA-prior-budget tests from the single-locus helper. Preserve tests for equivalence classes, effective length, SQUAREM convergence, VBEM baseline, and deterministic behavior.
- [tests/scenarios/test_gdna_diagnosis.py](../../tests/scenarios/test_gdna_diagnosis.py): stop printing `alpha_rna`.

PR-B verification:

```bash
conda activate rigel && pip install --no-build-isolation -e .
conda activate rigel && pytest tests/test_golden_output.py -v
conda activate rigel && pytest tests/ -v
conda activate rigel && ruff check src/ tests/
```

Acceptance checks:

- `rg "alpha_rna|rna_prior_count" src/rigel tests scripts` returns no references, except explicitly historical docs if those are marked as historical.
- Numeric golden comparisons pass without `--update-golden`.
- If any output/dataframe column is removed or renamed, schema-specific tests are updated intentionally. Do not describe that part as byte-identical.

### PR-C: Rename and Document the Surviving Prior

**Goal:** Make the codebase say what it now does: gDNA-only asymmetric prior plus mode baseline, not OVR prior allocation.

Recommended renames:

| Old | New | Notes |
|---|---|---|
| `compute_ovr_prior_and_warm_start` | `compute_gdna_prior_and_warm_start` | Rename after PR-B so the name is honest. |
| Python `alpha_gdna` fields where practical | `gdna_prior_count` | Prefer semantic names in Python and docs. Keep native/local `alpha` names only where they refer to Dirichlet vectors. |
| `gdna_prior` locus column | Keep as rate, but document as `gdna_prior_count / n_em_fragments` | Avoid overloading the count and rate. |

Documentation work:

- Create `docs/calibration/prior_model.md` as the source of truth. Keep it short and paste the three-line model from section 1.
- Update the module docstring in [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py).
- Update the native helper comment in [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp).
- Update user-facing parameter docs in [docs/parameters.md](../parameters.md), [docs/MANUAL.md](../MANUAL.md), and the relevant prior section of [docs/METHODS.md](../METHODS.md).
- Mark older implementation plans in `docs/calibration/` and `docs/bayesian_prior/` as historical rather than rewriting every old milestone document. Add a top banner such as: `Historical design note. The current prior model is documented in docs/calibration/prior_model.md.`
- Update [.github/copilot-instructions.md](../../.github/copilot-instructions.md) and [CLAUDE.md](../../CLAUDE.md) only if they mention removed knobs.

PR-C verification:

- `rg "OVR prior|prior_weight_rna|alpha_rna|rna_prior_count|nrna_weight" README.md CLAUDE.md .github docs src/rigel` returns either no matches or matches only inside clearly marked historical docs and this cleanup plan.
- A reader can derive the native prior vector from `docs/calibration/prior_model.md` alone.

### PR-D: MAP/VBEM Mode Discipline

**Goal:** Make the regularization mode explicit at the configuration boundary. The default remains VBEM, so production default math is unchanged.

Implementation:

- In [src/rigel/config.py](../../src/rigel/config.py), annotate `EMConfig.mode` as `Literal["vbem", "map"]` and keep the default `"vbem"`.
- In [src/rigel/cli.py](../../src/rigel/cli.py), make `--em-mode` an argparse `choices=("vbem", "map")` option with concise help text: `vbem` is the default and recommended mode; `map` is for research/debug runs and has no modeled RNA prior beyond numerical floors.
- Emit one warning per run when `mode == "map"`: `EM mode 'map' selected: RNA components have no modeled prior mass; results may differ materially from the recommended 'vbem' mode.`
- Document the tradeoff in [docs/parameters.md](../parameters.md) in one paragraph.

PR-D verification:

- Unit test the config validation and CLI choices.
- Run a small existing pipeline/CLI fixture with `--em-mode map` and assert the warning appears once.
- Golden tests remain unchanged because the default mode is unchanged.

### PR-E: `coverage_weights` Audit and Decision

**Goal:** Decide whether the remaining `coverage_weights` array is worth its memory cost. This is the only planned stage that may change EM numerics.

After PR-B, `coverage_weights` has two native consumers:

1. Row-normalized warm-start shares in `compute_gdna_prior_and_warm_start`.
2. The secondary sort key in `build_equiv_classes` for deterministic row order when log-likelihood rows tie.

Experiment design:

| Config | Behavior | Purpose |
|---|---|---|
| A | Current trapezoid `coverage_weights` as float64 | Baseline. |
| B | Current trapezoid weights narrowed to float32 end to end | Tests dtype savings without removing the signal. |
| C | Uniform warm start, no stored coverage-weight array | Tests whether the signal can be removed. |

Implementation guidance:

- Use a temporary benchmark branch or short-lived build flag. Do not merge a permanent `RIGEL_COVERAGE_WEIGHT_MODE` switch unless the final shipped design genuinely needs a user/developer knob.
- If Config C is tested, replace the secondary sort key with a deterministic key that does not depend on removed weights. The equivalence-class comparator should remain stable across repeated runs.
- Coordinate with [docs/performance/tier1_plan.md](../performance/tier1_plan.md), because dtype narrowing of `log_liks` and `coverage_weights` overlaps this work.

Metrics and decision rule:

- Compare transcript-level counts/TPM against Config A on VCAP plus at least three synthetic benchmark conditions: no gDNA, low gDNA, high gDNA at `ss=0.90`, including at least one `nrna_rand` condition.
- Report Spearman correlation, max relative drift on non-negligible counts, total mRNA/nRNA/gDNA mass fractions, per-locus SQUAREM iteration distributions, wall time, and peak RSS.
- Ship Config C only if drift is scientifically negligible and deterministic tests pass. Otherwise ship Config B if it is indistinguishable from A. If both diverge, keep A and document why the array is structurally important.
- Publish the decision in `docs/performance/` with command lines, conditions, metrics, and the selected code path.

---

## 3. Verification Matrix

Run all commands with the project environment active:

```bash
conda activate rigel
```

For PR-A through PR-D:

```bash
pip install --no-build-isolation -e .
pytest tests/test_golden_output.py -v
pytest tests/ -v
ruff check src/ tests/
```

Known pre-existing exception: `tests/test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna` may fail for unrelated reasons. Do not hide new failures behind that known issue.

Additional numeric identity check for PR-A through PR-D:

- Run one existing benchmark condition before and after each PR, preferably `gdna_high_ss_0.90_nrna_none`.
- Compare `quant.feather`, `gene_quant.feather`, `nrna_quant.feather`, and numeric columns of `loci.feather` row by row.
- If schemas changed, compare the intersection of intended numeric columns and record the schema change in the PR summary.

For PR-E:

- Run the benchmark matrix from section 2.
- If chosen output changes exceed golden tolerances, regenerate goldens with `--update-golden` and add a short `tests/golden/CHANGELOG.md` entry explaining the intentional coverage-weight decision.

---

## 4. Final Definition of Done

The cleanup is complete when:

1. `prior_weight_rna`, `build_prior_weight_rna`, calibration `nrna_weight`, `alpha_rna`, and `rna_prior_count` are absent from live `src/rigel/` code.
2. Any remaining `nrna_weight` name is either gone or clearly local to simulation abundance weighting, not calibration prior weighting.
3. The production native EM ABI accepts one calibration prior quantity per locus: the gDNA prior count, plus independent gDNA eligibility.
4. The current prior model is documented once in `docs/calibration/prior_model.md`; older design docs are marked historical instead of competing with it.
5. PR-A through PR-D show unchanged quantification numerics on goldens and at least one benchmark smoke condition.
6. PR-E has a published benchmark decision: keep, narrow, or remove `coverage_weights`.
7. The full suite and ruff pass, aside from the documented unrelated strand-LLR failure if it is still present.

---

## 5. Expected Payoff

| Stage | Main simplification | Numeric change? | Expected benefit |
|---|---|---|---|
| PR-A | Remove `prior_weight_rna` / `nrna_weight` | No | Deletes dead config, CLI, tests, and one float32 vector per locus. |
| PR-B | Remove `alpha_rna` / `rna_prior_count` | No for production batch EM | Deletes dead native ABI path, dead warm-start override, and misleading diagnostics. |
| PR-C | Rename and document | No | Makes the code readable without knowing the historical OVR design. |
| PR-D | MAP/VBEM discipline | No for default | Prevents accidental use of poorly regularized MAP mode. |
| PR-E | Decide `coverage_weights` | Maybe | Up to ~2 GB RSS reduction if uniform warm start ships; otherwise at least informs dtype narrowing. |

The strategic win is not just memory. The cleanup removes a false control surface: users and future contributors should no longer see an nRNA prior weight knob that cannot affect production EM under the current prior model.