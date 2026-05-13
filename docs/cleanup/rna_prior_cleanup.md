# RNA Prior + `coverage_weights` Cleanup Plan

**Date:** 2026-05-13
**Scope:** Strip the dead RNA-prior allocation machinery from rigel and align the code, tests, and documentation with what calibration v6 + the Bayesian-prior redesign actually do today.

This plan is the implementation half of the deep-dive interrogation that precedes it. It is organised as four PRs (A → D). PRs A, B, C are designed to produce **bit-identical EM output**; PR D is the only one that may change numbers and is gated on a benchmark.

---

## 0. Findings recap (the ground truth in production)

What the EM actually sees today, traced from [`assemble_priors`](../../src/rigel/calibration/locus_prior.py) → [`compute_ovr_prior_and_warm_start`](../../src/rigel/native/em_solver.cpp):

| Quantity | Calibration writes | C++ EM applies (MAP) | C++ EM applies (VBEM, default) |
|---|---|---|---|
| `alpha_gdna[i]` | `η_g(ℓ)` (calibration-derived) | `prior[gdna] = η_g` | `prior[gdna] = 0.5 + η_g` |
| `alpha_rna[i]` | **`0.0`** (hard-coded, [locus_prior.py:1039](../../src/rigel/calibration/locus_prior.py#L1039)) | `prior[t] = EM_LOG_EPSILON ≈ 1e-300` | `prior[t] = 0.5` (Jeffreys) |
| `prior_weight_rna[i]` | `1.0` mRNA / `nrna_weight` synth-nRNA / `1.0` gDNA | only used in `α_rna · w · cov / total` ⇒ **× 0** | same ⇒ **× 0** |

Cross-references in C++ for verification:
- Coverage-share warm start (alive): [em_solver.cpp:765-787](../../src/rigel/native/em_solver.cpp#L765-L787).
- RNA-allocation prior (dead): [em_solver.cpp:807-822](../../src/rigel/native/em_solver.cpp#L807-L822).
- gDNA warm-start override (dead, guarded by `α_rna > 0`): [em_solver.cpp:827-848](../../src/rigel/native/em_solver.cpp#L827-L848).
- Empty-locus uniform branch (dead, multiplies α_rna): [em_solver.cpp:1131](../../src/rigel/native/em_solver.cpp#L1131).

**The only live consumers of `coverage_weights` in production are:** (1) the warm-start `θ_init` in `compute_ovr_prior_and_warm_start`, and (2) the secondary-key tie-break inside `build_equiv_classes` ([em_solver.cpp:281-296](../../src/rigel/native/em_solver.cpp#L281-L296)).

---

## 1. Why this cleanup is bit-safe

Every RNA-prior term in the C++ EM is multiplied by `alpha_rna`. Production calibration sets `alpha_rna ≡ 0`. Therefore:

- **Removing the RNA-allocation arithmetic is mathematically equivalent to keeping it.** Every removed expression evaluates to `0` (modulo `EM_LOG_EPSILON ≈ 1e-300`, which is below double-precision relevance after the M-step normalisation).
- **Removing `prior_weight_rna`** is bit-safe for the same reason: its only consumer is multiplied by `α_rna`.
- **Removing the `(α_gdna / α_rna) * others` warm-start override** is bit-safe because its guard requires `α_rna > 0`.

The only place where this argument fails is if a downstream caller passes a non-zero `α_rna`. Rigel has exactly one such caller in production: it's `assemble_priors`, which sets it to `0`. The unit tests in [tests/test_em_prior_weight.py](../../tests/test_em_prior_weight.py) and [tests/test_em_impl.py](../../tests/test_em_impl.py) deliberately pass non-zero `α_rna` to exercise the dead path; these tests must be retired or rewritten.

PR-D (coverage_weights) is the only non-bit-safe step and is explicitly flagged as such.

---

## 2. PR-A — Cleanup A + B: strip dead `α_rna` / `prior_weight_rna` machinery

### Goal
Remove the entire RNA-allocation prior pathway from the ABI. **Bit-identical EM output.**

### Files touched

| File | Change |
|---|---|
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Drop `alpha_rna` + `prior_weight_rna` args from `compute_ovr_prior_and_warm_start`, `run_locus_em_native`, `run_batch_locus_em_partitioned`. Replace RNA-allocation block with a 3-branch `prior_out[i]` assignment. Drop the dead warm-start override. Drop `coverage_totals` accumulation. |
| [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py) | Drop `alpha_rna`, `rna_prior_count`, `prior_weight_rna` from `PriorTable`. Delete `build_prior_weight_rna`. Drop `nrna_weight` from `assemble_priors` signature. Drop `eta_r = 0.0` and the `rna_prior_count` field on `MultiLocusPrior` (or keep it for diagnostics only; see §2.4). |
| [src/rigel/calibration/__init__.py](../../src/rigel/calibration/__init__.py) | Remove `build_prior_weight_rna` from re-exports. |
| [src/rigel/config.py](../../src/rigel/config.py) | Delete `CalibrationConfig.nrna_weight`. |
| [src/rigel/cli.py](../../src/rigel/cli.py) | Delete the `--nrna-weight` flag if present. |
| [src/rigel/pipeline.py](../../src/rigel/pipeline.py) | Drop `prior_weight_rna_per_locus` plumbing (lines ~592, 595, 647, 677, 707-709, 762-764, 974, 991). Drop `nrna_weight` argument (line ~827). Drop `alpha_rna` from `_call_batch_em` and downstream batch tuples. |
| [src/rigel/estimator.py](../../src/rigel/estimator.py) | Drop `alpha_rna` and `locus_prior_weight_rna` from `run_batch_locus_em_partitioned` signature; drop the `alpha_rna` cast at line ~343 and the `locus_prior_weight_rna` cast at line ~362; drop the `"alpha_rna"` summary column. |
| [src/rigel/scan.py](../../src/rigel/scan.py) | Anywhere `alpha_rna` is forwarded — drop. |
| [src/rigel/sim/reads.py](../../src/rigel/sim/reads.py) | The two `nrna_weight` locals (lines 608, 645) are the simulator's own per-transcript abundance weights, **unrelated to the prior**. Verify by reading context and rename if confusion is likely; otherwise leave alone. |
| [tests/test_prior_weight_rna_policy.py](../../tests/test_prior_weight_rna_policy.py) | Delete. The policy it pinned is enforced structurally by removing the parameter. |
| [tests/test_em_prior_weight.py](../../tests/test_em_prior_weight.py) | Either delete or refactor to test the gDNA-only prior. Recommend delete; it adds no coverage that gDNA-prior tests don't already provide. |
| [tests/test_pipeline_integration_v6.py](../../tests/test_pipeline_integration_v6.py) | Delete `test_nrna_weight_does_not_break` (lines ~223-242). |
| [tests/test_calibrate_orchestrator.py](../../tests/test_calibrate_orchestrator.py) | Delete the `prior_weight_rna=[...]` and `t.prior_weight_rna == []` assertions (lines ~232, 256). |
| [tests/test_em_impl.py](../../tests/test_em_impl.py) | Drop the `alpha_rna` and `prior_weight_rna` arguments from the `run_locus_em_native` shim and any tests that pass non-zero `α_rna`. |
| [tests/conftest.py](../../tests/conftest.py) | Drop any `prior_weight_rna` plumbing if present. |
| [docs/](../../docs/) | Search-and-update: rg `prior_weight_rna|alpha_rna|nrna_weight|OVR prior` and rewrite. Defer to PR-B for the bulk of the docs work. |

### Concrete C++ replacement for `compute_ovr_prior_and_warm_start`

Replace the function body lines ~743-849 with this. Note: signature loses both `alpha_rna` and `prior_weight_rna`; rename to `compute_gdna_prior_and_warm_start` happens in PR-B.

```cpp
static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,    // [n_components] 1.0 if eligible, 0.0 otherwise
    double        alpha_gdna,        // calibration gDNA prior (physical count)
    int           gdna_idx,          // index of gDNA component (-1 if none)
    double*       prior_out,         // [n_components] output
    double*       theta_init_out,    // [n_components] output
    int           n_components,
    bool          use_vbem)
{
    // --- Coverage-weighted warm start (unchanged math; coverage_totals removed) ---
    std::copy(unambig_totals, unambig_totals + n_components, theta_init_out);
    for (const auto& ec : ec_data) {
        const int n = ec.n;
        const int k = ec.k;
        const int32_t* cidx = ec.comp_idx.data();
        const double*  wt   = ec.wt_flat.data();
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

    // --- Asymmetric Dirichlet prior: gDNA-only (calibration v6) ---
    // Mode-aware baseline: VBEM needs +0.5 (Jeffreys) to cancel digamma's
    // sparsifying bias; MAP has no bias, baseline 0.0.
    const double baseline = use_vbem ? 0.5 : 0.0;
    for (int i = 0; i < n_components; ++i) {
        if (eligible[i] <= 0.0) {
            prior_out[i] = 0.0;
        } else if (i == gdna_idx) {
            prior_out[i] = baseline + std::max(alpha_gdna, EM_LOG_EPSILON);
        } else {
            prior_out[i] = baseline;
        }
    }
}
```

Empty-locus branch in `run_locus_em_native` collapses to:

```cpp
if (n_units == 0 || n_candidates == 0) {
    std::vector<double> alpha(nc);
    double total = 0.0;
    const double baseline = use_vbem ? 0.5 : 0.0;
    for (size_t i = 0; i < nc; ++i) {
        const double p = (pe_ptr[i] > 0.0) ? baseline : 0.0;
        alpha[i] = unambig_totals[i] + p;
        total   += alpha[i];
    }
    // ... same normalization as before ...
}
```

(Behaviour change: empty locus under MAP previously returned a `1e-300` floor on every eligible component; this becomes `0.0`. After M-step normalisation the difference is below double-precision representability, so still bit-safe in practice. Worth a unit test of the empty-locus path before/after.)

### Verification protocol (must complete before merge)

1. **Build:** `pip install --no-build-isolation -e .`
2. **Goldens unchanged:** `pytest tests/test_golden_output.py -v` — must pass without `--update-golden`.
3. **Full suite:** `pytest tests/ -v` minus the deleted tests — must pass.
4. **Lint:** `ruff check src/ tests/`.
5. **Bit-identity smoke test:** Run [scripts/benchmarking](../../scripts/benchmarking) on one VCAP condition (`gdna_high_ss_0.90_nrna_none`) before and after PR-A. Diff `quant.feather` row-by-row — must match exactly.
6. **VCAP profile re-run:** `peak_rss_mb` should drop by ~ (n_loci × max_n_components × 4 bytes) ≈ a few hundred MB (one float32 array per locus eliminated). `compute_eb_gdna_priors` time should drop ~ 10–20 % (`build_prior_weight_rna` + per-locus float32 alloc gone).

### Definition of done

- All bullet points in §2 verified.
- The string `alpha_rna` appears nowhere in `src/rigel/`.
- The string `prior_weight_rna` appears nowhere in `src/rigel/`.
- `nrna_weight` appears only in `src/rigel/sim/reads.py` (and the simulator usage is a different concept; rename if necessary).

### Risks

- **Test fixtures pass non-zero `α_rna`.** Several tests (`test_em_impl.py`, `test_em_prior_weight.py`) deliberately exercise the now-removed code path. Their assertions about how non-zero `α_rna` shapes `θ` lose their target. Recommendation: delete; the production code path no longer offers this lever.
- **Downstream `summary.json` consumers.** The `"alpha_rna"` per-locus column disappears from outputs. Audit `scripts/benchmarking/analysis.py` and any notebooks that read it.

---

## 3. PR-B — Cleanup C + D: rename + documentation

### Goal
Make the surviving code self-explain. **Bit-identical EM output.**

### Renames (sequence with `vscode_renameSymbol` so all references update atomically)

| Old | New | Rationale |
|---|---|---|
| `compute_ovr_prior_and_warm_start` | `compute_gdna_prior_and_warm_start` | "OVR" (overlap-resolution) was the v3 redesign's name for the now-removed RNA-allocation prior. With α_rna gone the function only does two things: build the gDNA-asymmetric prior and the coverage-weighted warm start. |
| `MultiLocusPrior.gdna_prior_count` ↔ `PriorTable.alpha_gdna` | Pick one. Use `alpha_gdna` everywhere. | They are bit-identical aliases today and the redundancy is confusing. |
| `assemble_priors` | unchanged (good name) | Keep. |
| `enable_gdna_for_multilocus` | unchanged | Keep. |
| `prior_eligible` (C++) | unchanged | Keep. |

### Documentation rewrites

Files to update:

- [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py) module docstring (lines 1-22): replace the `(alpha_gdna, alpha_rna, prior_weight_rna)` triple with the single-argument story.
- [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py) `MultiLocusPrior` docstring (lines ~141-167): drop the `eta_r` paragraph; restate the asymmetric prior in two sentences.
- [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) header comment (lines ~735-742): replace the "Coverage-weighted warm start + unified OVR prior" section with the actual three-line prior model below.
- [docs/calibration/](../../docs/calibration/) anywhere that references the OVR prior or `α_rna`.
- [docs/bayesian_prior/bayesian_prior_plan_v3.md](../../docs/bayesian_prior/bayesian_prior_plan_v3.md): annotate the document with a "Status: shipped, see source-of-truth in this folder" header so a future reader does not mistake it for live design.
- [README.md](../../README.md), [CLAUDE.md](../../CLAUDE.md), [.github/copilot-instructions.md](../../.github/copilot-instructions.md): one-line edits if any of them mention `nrna_weight` or `alpha_rna`.

### Canonical 3-line prior model (paste into module docstring + a new `docs/calibration/prior_model.md`)

> Rigel uses a per-multi-locus asymmetric Dirichlet prior with **only the gDNA component carrying calibration mass**:
>
> 1. **gDNA component:** `prior[gdna] = baseline + max(η_g(ℓ), ε)` where `η_g(ℓ)` is the calibration-derived global expected gDNA pseudocount and `baseline = 0.5` under VBEM (Jeffreys, cancels the digamma sparsification bias) or `0.0` under MAP.
> 2. **RNA components:** `prior[t] = baseline` (the same Jeffreys/zero baseline, no per-component allocation).
> 3. **Eligibility:** an independent per-component `0/1` flag from `enable_gdna_for_multilocus` controls whether the gDNA component is even instantiated; when ineligible, `prior[gdna] = 0` and `prior[t] = baseline` for all RNA components.

### Verification

- `pytest tests/ -v` passes (renames are syntactic).
- `pytest tests/test_golden_output.py` passes (no math changed).
- Read the new docs end-to-end: any reader unfamiliar with the codebase can derive the prior table from the documentation alone.

### Risks

- Renames touch many files. Use the language server's rename, not search-and-replace.
- Docs in `docs/calibration/` and `docs/bayesian_prior/` are voluminous; budget effort accordingly.

---

## 4. PR-C — Cleanup E: MAP / VBEM mode discipline

### Goal
Make the (radically different) MAP and VBEM regularisations explicit at the configuration boundary. **Bit-identical EM output for the production default (VBEM).**

### Background

After PR-A, MAP gets `prior[t] = 0` for every RNA component while VBEM gets `prior[t] = 0.5` (Jeffreys). These are not knobs we should expose blindly: MAP without any RNA prior degenerates on isoform-rich loci (the M-step zero-variance singularity that Jeffreys exists to prevent), and most users will want VBEM regardless.

### Implementation

| File | Change |
|---|---|
| [src/rigel/config.py](../../src/rigel/config.py) | Default `EMConfig.mode = "vbem"` (verify; if already, leave). Add a `Literal["vbem", "map"]` annotation. |
| [src/rigel/cli.py](../../src/rigel/cli.py) | Promote `--em-mode` to `--em-mode {vbem,map}` with explicit help text: `vbem` (default; recommended), `map` (research/debug only — no RNA prior, may degenerate on isoform-rich loci). |
| [src/rigel/pipeline.py](../../src/rigel/pipeline.py) or [src/rigel/estimator.py](../../src/rigel/estimator.py) | When `mode == "map"`, emit a single `logger.warning(...)` once per run: `"EM mode 'map' selected: RNA prior is uniform zero. Quantitative results may differ materially from the recommended 'vbem' mode."` |
| [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) | No code change; the existing `use_vbem` flag already gates the baseline. |
| [docs/parameters.md](../../docs/parameters.md) | Document the mode trade-off in one paragraph, not five. |

### Optional follow-up (not in this PR)

Consider a `mode = "vbem_jeffreys" | "vbem_uniform" | "map"` triple and let users dial the baseline. Out of scope here.

### Verification

- `pytest tests/ -v` passes.
- Run the CLI with `--em-mode map` on the conftest mini fixtures; confirm the warning fires and output is produced.
- Goldens unchanged (default mode is unchanged).

### Risks

- If any internal script or test passes `mode="map"`, the new warning will surface in CI logs. Acceptable.

---

## 5. PR-D — Cleanup F: `coverage_weights` audit (now scoped)

### Goal
Decide whether `coverage_weights` survives at all. **This is the only PR that may change EM output.**

### What's left of `coverage_weights` after PRs A–C

After PR-A, `coverage_weights` has exactly **two** consumers in C++:

1. **Warm-start `θ_init`** in `compute_gdna_prior_and_warm_start` ([em_solver.cpp:765-787 post-PR-A](../../src/rigel/native/em_solver.cpp#L765-L787)). Uses per-row coverage-normalised shares to seed `θ`.
2. **Tie-break secondary sort key** in `build_equiv_classes` ([em_solver.cpp:281-296](../../src/rigel/native/em_solver.cpp#L281-L296)). Used only when `ll_flat` rows are exactly equal. Pure determinism.

The 2.08 GB float64 array in `ScoredFragments` exists solely for these two purposes.

### Experiment

A single C++ macro toggles between three configurations:

| Config | scoring.cpp behaviour | Expected outcome |
|---|---|---|
| **A — baseline** | `cov_wt = compute_fragment_weight(...)` (current) | reference |
| **B — float32 cast** | `cov_wt = static_cast<float>(compute_fragment_weight(...))` (and float32 throughout) | should be bit-equivalent within EM tolerance |
| **C — uniform** | `cov_wt = 1.0` always | tells us whether the trapezoid weighting matters at all |

Implementation:

- Add `RIGEL_COVERAGE_WEIGHT_MODE` build-time macro to [src/rigel/native/scoring.cpp](../../src/rigel/native/scoring.cpp) that switches between the three. Default = A. CMake exposes a `-DRIGEL_COVERAGE_WEIGHT_MODE={trapezoid,trapezoid_f32,uniform}` option for the build.
- Add a tie-break test toggle so config C still produces deterministic output (fall back to lexicographic on `ll_flat` only, then `comp_idx`, then `t_indices`).

### Metrics (for each config × {VCAP, two sim conditions})

- Per-transcript abundance Spearman ρ vs A. Target: B > 0.9999; C > 0.999.
- Per-locus EM iteration count distribution. Target: C median ≤ 110 % of A median.
- Total mRNA / nRNA / gDNA mass fractions. Target: within 0.5 % of A.
- Wall time + peak RSS deltas (these are the prize).

Expected savings:

| Config | Peak RSS savings | Wall savings |
|---|---|---|
| B (float32) | ~1 GB on VCAP | ~5 % in scoring + EM build |
| C (uniform) | ~2 GB on VCAP | ~10 % in scoring + EM build + warm start |

### Decision rule

- **C ≈ A** ⇒ ship C. Drop the `coverage_weights` column from `ScoredFragments`, scoring.cpp, CSR, EM ABI, and tests. Replace the warm-start coverage shares with `1/k`. Replace the sort tie-break with `comp_idx` ordering.
- **C diverges, B ≈ A** ⇒ ship B. Float32 cast end-to-end. ABI changes per [§1.2 of the Tier-1 plan](../performance/tier1_plan.md#12).
- **Both diverge** ⇒ ship neither. Document the structural dependency in the EM and revisit only if profiling later demands it.

### Code-archaeology hypothesis

The trapezoid-coverage warm start was a port of salmon's ECC weighting. salmon's EM started cold (uniform `θ`) and benefitted; rigel's EM has a calibration-derived gDNA prior plus the coverage-share warm start does its work in **one** pass before SQUAREM. SQUAREM then makes O(10) accelerated steps on the M-step posterior, dwarfing any warm-start asymmetry. **Predicted outcome: C ≈ A within tolerance**, allowing the array to be removed entirely.

### Verification (post-decision)

- Re-run all golden tests with `--update-golden` if the chosen config produces ≥ 1e-5 relative drift; document in [tests/golden/CHANGELOG.md](../../tests/golden/) (create if absent).
- VCAP profile re-run: report new `peak_rss_mb` and `quant_from_buffer` wall.
- Sim-suite benchmark via [scripts/benchmarking](../../scripts/benchmarking) on at least three conditions (no/low/high gDNA × ss=0.90).

### Risks

- **Determinism regression** if the tie-break collapses to a non-stable comparator. Add a regression test that exercises duplicate `ll_flat` rows and asserts deterministic `θ` across repeated runs.
- **nRNA / mRNA tradeoff drift.** The warm start affects which side of the EM basin the solver converges to in pathological loci. The benchmark conditions in scripts/benchmarking/configs/default.yaml should cover this; if they don't, add one with `nrna_rand` + isoform-rich locus.

---

## 6. PR sequencing & expected payoffs

| PR | Items | Math change? | LOC removed | Expected RSS Δ | Expected wall Δ |
|---|---|---|---|---|---|
| **PR-A** | §2 (Cleanups A + B) | No (bit-safe) | ~250 | ~300–500 MB | −10 to −20 % on `compute_eb_gdna_priors` |
| **PR-B** | §3 (Cleanups C + D) | No (rename + docs) | 0 (net) | 0 | 0 |
| **PR-C** | §4 (Cleanup E) | No (default unchanged) | ~10 net adds | 0 | 0 |
| **PR-D** | §5 (Cleanup F) | **Yes — gated on benchmark** | up to ~400 if config C ships | up to **−2 GB** on VCAP | up to −10 % in scoring + EM build |

Combined wall-time payoff if all four merge with config C: estimated **~4 s on VCAP quant** plus the entire `prior_weight_rna` allocation overhead per locus.

---

## 7. Definition of done (cleanup as a whole)

A cleanup PR is done when:

1. The full test suite passes (`pytest tests/ -v`), excluding the documented pre-existing strand-LLR failure.
2. `ruff check src/ tests/` is clean.
3. For PRs A, B, C: golden outputs are unchanged and a VCAP smoke run produces a `quant.feather` byte-identical to the pre-PR run on at least one condition.
4. For PR-D: the benchmark report in §5 is published (small markdown drop in `docs/performance/`) and the chosen configuration is justified against the decision rule.
5. After all four PRs, the canonical 3-line prior model in §3 is the only prior-related explanation a reader needs.
