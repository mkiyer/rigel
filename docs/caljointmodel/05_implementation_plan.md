# 05 — Implementation Plan

Seven phases. Each is independently committable. After Phase 2 the
pipeline does not run end-to-end. After Phase 5 it does. Phase 4 is
broken into three sub-phases, one per calibration goal G2–G4 (G1 was
dropped — see [`00_overview.md`](00_overview.md) §1).

```
Phase 1 — Archive            (move salvageable algorithms to archive/)
Phase 2 — Burn               (delete legacy, leave gaping hole)
Phase 3 — Scaffold           (native substrate extension, public types, validation)
Phase 4 — Implement          (G2 + G3 unified machinery, G4)
Phase 5 — Integrate          (locus prior rewrite, pipeline wiring)
Phase 6 — Validate           (synthetic, paralog, hybrid-capture, armis2)
Phase 7 — Cleanup            (regenerate goldens, magic-number audit)
```

---

## Phase 1 — Archive

**Goal.** Preserve the SRD gDNA FL recovery + empirical-Bayes FL
smoothing algorithms for reference (they are no longer needed by the
calibrator, but may be useful to downstream EM scoring in the future).
Move everything else into the archive tree.

**Tasks.**

1. `mkdir -p archive/calibration_legacy_2026_05/src/rigel/calibration archive/calibration_legacy_2026_05/tests`.
2. `git mv` ARCHIVE files per [`02_failure_audit.md`](02_failure_audit.md) §4 into the archive tree, preserving relative paths.
3. `git mv` DELETE files into the archive tree as well, but rename their extension to `DELETED.md` so they cannot be re-imported.
4. `git mv` calibration tests targeting deleted modules into `archive/.../tests/`.
5. Write `archive/calibration_legacy_2026_05/README.md`: SHA, ARCHIVE vs DELETE distinction, no-resurrection warning citing this design.

**Done.**
- [ ] `ls src/rigel/calibration/` shows only the KEEP list.
- [ ] `ls archive/calibration_legacy_2026_05/` contains all moved files.
- [ ] Archive README exists.

**State.** Repo does not import; tests don't collect. Expected.

---

## Phase 2 — Burn

**Goal.** Excise every call to deleted modules from the live codebase. Replace the calibration call site in `pipeline.py` with a `NotImplementedError` stub.

**Tasks.**

1. Create `src/rigel/calibration/_stub.py`:
   - Empty `CalibrationConfig` and `CalibrationResult` dataclasses with the [`04_interface_contract.md`](04_interface_contract.md) §3 / §5 field signatures, defaulting to zeros / empty arrays.
   - `calibrate(...)` raises `NotImplementedError("Calibration burn-down; see docs/caljointmodel/")`.
2. Rewrite `src/rigel/calibration/__init__.py` to re-export only the three symbols.
3. Strip `quant_from_buffer` in `pipeline.py` to call only the stub. Delete the entire legacy calibration block (~150 lines).
4. Delete calibration imports from `src/rigel/estimator.py`, `src/rigel/locus.py`. Replace the locus prior call with a stub that raises `NotImplementedError` until Phase 5.
5. Delete legacy calibration sub-configs from `src/rigel/config.py`. Add a single stub `calibration: CalibrationConfig` field on `PipelineConfig`.
6. Strip calibration-related CLI flags from `cli.py`.
7. `CHANGELOG.md` entry: `2026-05-29: Calibration system burn-down (replacement: docs/caljointmodel/)`.

**Done.**
- [ ] `python -c "import rigel"` succeeds.
- [ ] `pytest --collect-only tests/` is clean.
- [ ] `rigel quant ...` reaches calibration call and raises `NotImplementedError`.
- [ ] `git grep -nE 'background_model|fusion|latent_states|boundary_sweep|strand_deconv|adaptive_prior|p_unexpressed' src/` returns zero hits.
- [ ] `wc -l src/rigel/calibration/*.py` < 200 lines.

---

## Phase 3 — Scaffold

**Goal.** Land the native substrate extension, the public types, and the input-validation layer. No inference yet.

**Tasks.**

1. **Native: add `CalibrationAggregates` struct** per [`04_interface_contract.md`](04_interface_contract.md) §2.1.
   - In `bam_scanner.cpp::observe_calibration_fragment`, increment the new region counters (and the appropriate boundary counters when the fragment crosses a boundary).
   - Add nanobind bindings in `_bam_impl`.
   - Rebuild: `pip install --no-build-isolation -e .`.
2. **Python: substrate adapter.**
   - Add `CalibrationSubstrate` dataclass in `scan_payload.py` per §2.3.
   - Add `build_calibration_substrate(buffer, index)` helper in `pipeline.py`.
3. **Public types in `calibrate.py`.**
   - Real `CalibrationConfig` per §3.
   - Real `CalibrationResult` per §5 with `__post_init__` invariants per §5.1.
   - Stub `calibrate(...)` that validates substrate (raises `CalibrationSubstrateError` per §4.1) and returns a placeholder result with masses = 0, $\omega = 1$, hyperparameters at the initialization values from [`03_inference.md`](03_inference.md) §7.
   - Delete `_stub.py`; update `__init__.py`.
4. Add exception types `CalibrationSubstrateError`, `CalibrationConvergenceError`.
5. Tests in `tests/calibration/`:
   - `test_native_aggregates.py` — round-trip a synthetic BAM through native; verify all 9 counters.
   - `test_substrate_invariants.py` — each invariant violation raises.
   - `test_config_defaults.py` — defaults match doc 04 §3.
   - `test_result_invariants.py` — each `__post_init__` check fires.

**Done.**
- [ ] Native build succeeds.
- [ ] All scaffold tests pass.
- [ ] `rigel quant ...` runs through calibration stub returning placeholder result; aborts at locus prior consumer.

---

## Phase 4 — Implement

Three sub-phases. Each independently committable with its own test set.

### Phase 4a — G2 / G3 unified: per-region/boundary deconvolution

**Goal.** Implement the per-region (and per-boundary) soft allocation using count + strand channels and the splice-deterministic-RNA rule.

The function is one piece of code, applied three times (contained, left, right) sharing global hyperparameters.

**Tasks.**

1. Implement `_per_region_estep(substrate_set, pi_g_prior, rho_d_bb, rho_r_bb, kappa_rna, eps_s, omega, rho_0, L_eff, phi, M_d_unspl_prev) -> (M_g, M_d, k_plus_g_hat, k_plus_d_hat)` per [`03_inference.md`](03_inference.md) §3. `kappa_d = 0.5` is a module-level constant, not an argument.
   - Vectorized over $|S|$.
   - Returns mass arrays consistent with the §3.6 totals.
2. Tests in `tests/calibration/test_g2_g3_deconvolution.py`:
   - On synthetic substrate with known $\pi_r^{(g)}$, $M_r^{(g)}$ recovered within 5% relative for regions with $\geq 20$ fragments.
   - **Hybrid-capture sanity**: captured exon with 90% capture, depleted intronic flank → contained gDNA mass low, boundary gDNA mass high.
   - **Paralog sanity**: pair-anchored gDNA straddlers with 126/26 strand split → $M_r^{(d, \text{cont})} < 10$ via the three-leg rescue (count NB + BB strand + zero spliced; see [`02_failure_audit.md`](02_failure_audit.md) §2.5). **Risk-flagged**: weaker than the prior FL-bearing design; benchmarks in Phase 6 must verify.
   - Mass aggregation conserves the substrate counts: $M_r^{(g)} + M_r^{(d)}$ equals total fragment count to $10^{-9}$.
   - Vectorized output matches scalar reference on 1000 random regions.
3. Boundary half-split helper:
   - `_boundary_half_split(M_g_L, M_g_R) -> M_g_boundary_contribution`. Trivial; covered by an explicit conservation test.

**Done.**
- [ ] All G2/G3 tests pass.

### Phase 4b — G4: per-region exposure (closed-form Gamma)

**Goal.** Implement the closed-form Gamma posterior on $\omega_r$ per [`03_inference.md`](03_inference.md) §4.

**Tasks.**

1. Implement `_update_exposure(M_g_tot, rho_0, L_eff, phi) -> (omega, log_omega_var)`.
2. Tests in `tests/calibration/test_g4_exposure.py`:
   - Synthetic dataset with known $\omega_r$: recovered values within 10% relative for regions with $\geq 10$ expected gDNA fragments.
   - Empty region: $\omega_r = 1$, $\log\sigma_r^2 = \phi$.
   - Variance scales as $1/(1/\phi + M)$.

**Done.**
- [ ] All G4 tests pass.

### Phase 4c — Global M-step + outer loop

**Goal.** Implement the five global hyperparameter M-steps and the outer-loop driver.

**Tasks.**

1. Implement `_m_step_rho_0` (closed form per [`03_inference.md`](03_inference.md) §5.1). No `_m_step_kappa_d` exists — $\kappa_d = 0.5$ is fixed.
2. Implement `_m_step_eps_s` (§5.5).
3. Implement `_m_step_phi` — `scipy.optimize.minimize_scalar` bounded on `(_PHI_FLOOR, 1e2)` with moment estimator warm start (§5.4).
4. Implement `_m_step_rho_d_bb` — `minimize_scalar` bounded on `(_BB_FLOOR, 1 - _BB_FLOOR)` (§5.2). Symmetric BB ($\alpha = \beta$, mean 0.5 fixed).
5. Implement `_m_step_rho_r_bb` — same pattern, asymmetric BB driven by pre-computed per-region $\kappa_r^{\text{rna}}$ (§5.3).
6. Implement `_update_pi_g_prior` (§5.6).
7. Implement the `calibrate(...)` outer loop per [`03_inference.md`](03_inference.md) §1.
8. Tests in `tests/calibration/test_m_step_and_outer.py`:
   - Each M-step recovers its true value within tolerance on synthetic data (per §1.4 of [`06_validation_plan.md`](06_validation_plan.md)).
   - Mass-change diagnostic monotone decreasing over outer iterations.
   - Outer loop converges in ≤ 25 iterations on synthetic.
9. End-to-end test `tests/calibration/test_e2e_synthetic.py`:
   - $R = 1000$, simulated unspliced + spliced + boundary substrates with known $(\rho_0, \phi, \rho_d^{\text{BB}}, \rho_r^{\text{BB}}, \omega_r)$. Recover hyperparameters within tolerances; recover per-region exposure within median log-error 0.10.

**Done.**
- [ ] All M-step and outer-loop tests pass.
- [ ] E2E synthetic test passes.
- [ ] No Tier-1 magic numbers in `calibrate.py` (self-audit per [`02_failure_audit.md`](02_failure_audit.md) §6).
- [ ] Mass-change is monotone decreasing on every synthetic run.
- [ ] `rigel quant ...` runs through calibration successfully; still aborts at locus prior consumer (Phase 5).

---

## Phase 5 — Integrate

**Goal.** Wire the calibrator into the locus prior consumer.

**Tasks.**

1. Rewrite `src/rigel/locus.py::compute_locus_priors_from_partitions` to consume `CalibrationResult` per [`04_interface_contract.md`](04_interface_contract.md) §6. Implement the §6.2 pseudocount formulas verbatim.
2. Delete the Phase 2 `NotImplementedError` stub.
3. Tests in `tests/calibration/test_locus_prior_consumer.py`:
   - Synthetic locus with known masses → $\alpha_t^{(d)}, \alpha_t^{(g)}$ match the §6.2 formulas to $10^{-10}$.
   - Mass-conservation invariant (§6.3) holds.
   - Symmetric paralog locus → symmetric pseudocounts.
   - Empty-region case → only gDNA prior contribution.
4. Smoke test: `rigel quant` on a minimal scenario BAM runs end-to-end, produces non-empty `quant.feather`.

**Done.**
- [ ] `rigel quant` runs end-to-end on at least one scenario BAM.
- [ ] All consumer tests pass.
- [ ] `pytest tests/test_pipeline_smoke.py` passes.

---

## Phase 6 — Validate

**Goal.** Prove the rewrite passes the regressions the old system passed AND fixes the paralog + hybrid-capture cases.

**Tasks** (full details in [`06_validation_plan.md`](06_validation_plan.md)):

1. **Paralog regression** — `tests/scenarios_aligned/test_paralogs.py`. Expected: t1 ≈ t2 ≈ 500.
2. **Hybrid-capture regression** — new scenario per [`06_validation_plan.md`](06_validation_plan.md) §3.2.
3. **Synthetic benchmark sweeps** — `scripts/benchmark/configs/locus_simple_*.yaml`.
4. **Armis2 real-data smoke** — corner conditions per [`06_validation_plan.md`](06_validation_plan.md) §5.
5. **Numerical-stability stress** — tiny region sets, all-zero regions, degenerate strand.

**Done.**
- [ ] Paralog test passes naturally (no tolerance loosening).
- [ ] Hybrid-capture scenario test passes.
- [ ] Synthetic sweeps document any regression with root cause.
- [ ] All `CalibrationResult` fields finite across all runs.
- [ ] Mass-change history monotone in every run.

---

## Phase 7 — Cleanup

**Tasks.**

1. Regenerate `tests/golden/*` with `pytest tests/ --update-golden`.
2. `ruff check src/ tests/` and `ruff format src/ tests/`.
3. Magic-number audit: ≤ 8 in calibration, each annotated.
4. Update `CLAUDE.md` and `.github/copilot-instructions.md` for new architecture.
5. Write `docs/caljointmodel/07_postmortem.md`.
6. Decide on `archive/calibration_legacy_2026_05/` retention. Recommendation: keep until next minor release.

**Done.**
- [ ] `pytest tests/` all green.
- [ ] `ruff check` clean.
- [ ] Magic-number audit ≤ 8, each annotated.
- [ ] Docs updated.

---

## Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Three-leg paralog rescue (count NB + BB strand + zero spliced) empirically weaker than the prior FL-bearing design | Medium-High | High | Phase 6 paralog regression and synthetic sweeps must verify rescue. If it fails, options: (a) reintroduce FL only for paralog-flagged regions, (b) tighten BB priors, (c) accept partial regression and document. See [`02_failure_audit.md`](02_failure_audit.md) §2.5 risk flag. |
| Per-region BB strand likelihood is unstable at very small $\rho_d^{\text{BB}}$ | Low | Medium | scipy's `betabinom.logpmf` handles the limit; Newton bound on $\rho_d^{\text{BB}}$ prevents drift to 0. |
| Boundary half-split mis-attributes mass at very large boundaries | Low | Medium | `boundary_split_factor` is a config knob. Sweep in Phase 6. |
| Mass-change diagnostic increases between iterations (EM violation) | Very low | High | Indicates bug. Raise `CalibrationConvergenceError`. |
| Locus prior pseudocount underflows or overflows on extreme regions | Medium | High | Consumer unit tests in Phase 5 cover extreme-mass cases; clip $\alpha_t \in [\alpha_{\min}, \alpha_{\max}]$. |

## Schedule

No calendar estimate. Sequence is the contract; throughput is not.
Phases 1, 2, 7 are mechanical and small. Phase 3 includes a native
rebuild. Phase 4 (especially 4b and 4d) and Phase 6 are where the
substance is.
