# New Calibration Implementation Log

Date started: 2026-05-26
Plan: `docs/newcal/new_cal_plan_v6.md`

## Phase I - Completed Baseline

- Added additive v6 primitives without rewiring the current calibration orchestrator.
- Added `calibration.boundaries.BoundaryTable` derived from current region-side ledger slots.
- Added compartment-aware strand deconvolution helpers sharing one library-level `kappa_d`.
- Added `calibration.background_model.BackgroundModel` and `fit_background_model()`.
- Verified focused Phase I tests plus existing strand deconvolution tests.

Deviation from plan:

- Did not remove the old density precision-cap/exposure path during Phase I. Kept it intact to preserve current behavior until v6 downstream wiring is gated and benchmarked.

Issues encountered:

- Terminal boundary slots needed explicit zero-fill semantics. The implementation scatters only internal region sides, so terminal-compatible mass does not enter the boundary table.
- TS_NONE/TS_AMBIG regions need to preserve raw totals in the compartment path; otherwise background fitting would silently lose intergenic mass. The new compartment path treats those regions conservatively as all gDNA when strand deconvolution is ineligible.

## Phase II - Boundary Imputation And Sweeps

### 2026-05-26 Checkpoint 1

Progress:

- Began Phase II implementation.
- Inspected existing Gamma/NB predictive math in `density_model.py` and Phase I boundary/background primitives.
- Chosen scope remains additive: add `boundary_model.py` and `boundary_sweep.py`, export the API, and test directly without orchestrator rewiring.

Deviations from plan:

- None yet.

Issues encountered:

- None yet.

### 2026-05-26 Checkpoint 2

Progress:

- Added `calibration.boundary_model` with `BoundaryLocalPosterior`, local Gamma-excess construction, and NB contained-count prediction.
- Added `calibration.boundary_sweep` with evidence-reliability transfer weights and simple sequential forward/reverse scans.
- Exported Phase II API from `rigel.calibration`.

Deviations from plan:

- Added a small shared helper, `predict_contained_gdna_from_excess()`, in `boundary_model.py` so local imputation and swept evidence use identical NB prediction math.

Issues encountered:

- Need to be careful that only excess Gamma evidence is swept. The background Gamma posterior is used only at prediction time; it is not stored in `alpha_excess`/`beta_excess` and is not propagated across boundaries.
- While implementing transfer weights, noticed `BoundaryTable.left_region_index()` / `right_region_index()` returned adjacent regions for terminal boundaries despite documenting `-1`. Fixed the helper contract and added a test. Terminal mass/opportunity behavior was already zero-filled, so this was a helper/debug correctness issue rather than a sweep-output issue.

### 2026-05-26 Checkpoint 3

Progress:

- Added `tests/test_boundary_model.py` for local background shrinkage, strong boundary prediction, and strand-deconvolved boundary counts.
- Added `tests/test_boundary_sweep.py` for zero-weight identity, attenuated propagation, transfer-weight behavior, and an internal-exon-like propagation case.
- Re-ran focused Phase I/II tests plus calibration/pipeline smoke tests.

Validation:

- `ruff check` on edited calibration modules and new tests passed.
- `pytest tests/test_boundary_model.py tests/test_boundary_sweep.py tests/test_boundaries.py tests/test_compartment_strand_deconv.py tests/test_background_model.py tests/test_strand_deconv.py tests/test_calibrate.py tests/test_calibration_integration.py tests/test_pipeline_wiring.py -v` passed: 44 tests.

Deviations from plan:

- Phase II remains additive and is not wired into `_orchestrator.py`. This matches migration safety, but means benchmark behavior is unchanged until Phase III/IV wiring.

Issues encountered:

- No unresolved Phase II implementation issues at this checkpoint.

## Phase III - Four-State Tensor Solver

### 2026-05-26 Design Checkpoint 1

Progress:

- Published the Phase III design document at `docs/newcal/phase3_four_state_solver.md`.
- Captured the intended four biological states, annotation-derived log priors, log Bayes factor tensor assembly, calibration E/M loop, `RegionCalibration` schema, and verification strategy.
- Confirmed Phase III remains design-only so far; no solver code has been implemented yet.

Deviations from plan:

- None implemented yet. The design document is intended to specialize the Phase III section of `docs/newcal/new_cal_plan_v6.md` before code changes begin.

Issues encountered:

- The first attempt to append this checkpoint through the terminal corrupted the log tail. Repaired the checkpoint and will use patch/file-edit tools for log updates going forward.
- The design document needs review for consistency with v6 before implementation, especially `gamma_r` in `RegionCalibration`, exact field names for region effective lengths, and whether `region_arrays` should be explicit in the E-step signature.

### 2026-05-26 Design Checkpoint 2

Progress:

- Reviewed the detailed Phase III design against current Phase I/II APIs: `DensityObservation`, `RegionArrays`, `BackgroundModel`, `BoundaryLocalPosterior`, `BoundarySweepResult`, and compartment strand deconvolution.
- Updated `docs/newcal/phase3_four_state_solver.md` so the pseudocode uses actual fields (`observation.contained_leff`, `observation.spliced_count`) and passes `RegionArrays` explicitly to the E-step.
- Added `gamma_r` to the `RegionCalibration` schema and E-step output derivation, matching the v6 downstream diagnostics/exposure contract.
- Fixed tensor assembly pseudocode to use explicit singleton slices for NumPy broadcasting.
- Clarified that fractional routed counts are rounded for the NB density likelihood term.

Deviations from plan:

- The detailed design now treats `capture_enrichment_target` as a scalar carried on `RegionCalibration` for traceability, in addition to pass diagnostics. It remains an internal stabilizer, not a user-facing knob.
- `DensityObservation.spliced_count` is promoted from diagnostic-only to expression evidence in the v6 tensor path. This is an intentional Phase III contract extension and should not affect the current v5 density model until Phase IV wiring.
- The intronic prior was corrected to give `background` a positive soft preference, matching the v6 text that intron-only regions should favor background before data are considered.

Issues encountered:

- Exact strand evidence is still specified as a wrapper around existing `strand_log_likelihood_d_grid()`, but the precise state-specific log-likelihood mapping needs careful implementation and tests.
- The E-step currently defines `mu_gdna` and `A_r` as posterior mixtures of off-target and captured predictions. This matches the smooth-state contract, but it should be stress-tested on ambiguous unstranded capture fixtures.
- The current `DensityObservation` docstring says `spliced_count` is diagnostic-only; Phase III implementation should update or qualify that comment when `latent_states.py` starts consuming it.

### 2026-05-26 Audit Checkpoint - Phases 0/I/II And Phase III Priors

Progress:

- Audited `docs/newcal/new_cal_plan_v6.md` against implemented Phase I/II modules.
- Confirmed additive Phase I/II primitives exist and are exported: `BoundaryTable`, background fitting, compartment strand deconvolution, local boundary imputation, and boundary sweeps.
- Ran focused Phase I/II primitive tests: `pytest tests/test_boundaries.py tests/test_compartment_strand_deconv.py tests/test_background_model.py tests/test_boundary_model.py tests/test_boundary_sweep.py -v` passed 17 tests.
- Reviewed Phase 0 outside-review dispositions as implementation constraints. Phase I/II satisfy the accepted decisions that apply to them: simple sequential sweeps, evidence-reliability transfer weights, exact boundary table contract, one library-level `kappa_d`, and no new C++/probe BED/native boundary storage.
- Revised Phase III annotation priors so intergenic and intron-only regions favor `background` rather than `gdna_only_capture` before evidence.

Deviations from plan:

- Phase I in v6 says to remove the old `400` precision-cap path and old density exposure source. Current implementation has intentionally stayed additive: `_orchestrator.py` still uses `fit_density_evidence()` and `RegionExposure.from_density()`, and `density_model.py` still defines `DENSITY_PRIOR_MAX_PRECISION = 400`. This is compatible with migration safety but not a full implementation of the Phase I removal item.
- Phase I/II primitives are not wired into the live calibration orchestrator yet. They are available and tested, but production quantification still follows the current v5 density/fusion path until Phase IV gating/wiring.
- Phase III detailed plan originally gave intergenic regions a modest capture preference. That is now corrected: annotation alone should not prefer capture; boundary/density likelihoods should earn captured-state probability.

Issues encountered:

- No Phase I/II test failures in the focused suite.
- Audit finding: the current implementation state should be described as “Phase I/II primitives implemented and validated, but not yet live in the orchestrator,” not “v6 Phase I/II fully wired.”

### 2026-05-26 Phase III Implementation Checkpoint 1

Progress:

- Added `src/rigel/calibration/latent_states.py` with four-state constants, annotation-derived priors, log Bayes factor builders, log-tensor assembly, and probability normalization.
- Added `src/rigel/calibration/calibration_iteration.py` with `RegionCalibration`, `CalibrationStepResult`, `calibration_e_step()`, `calibration_m_step()`, and `run_calibration_iteration()`.
- Exported the Phase III API from `rigel.calibration`.
- Updated `DensityObservation` documentation to note that `spliced_count` remains diagnostic-only for the v5 density model but is consumed by the additive v6 tensor path.
- Added `tests/test_latent_states.py` and `tests/test_calibration_iteration.py`.

Validation:

- `ruff check src/rigel/calibration/latent_states.py src/rigel/calibration/calibration_iteration.py src/rigel/calibration/__init__.py src/rigel/calibration/density_observation.py tests/test_latent_states.py tests/test_calibration_iteration.py` passed.
- `pytest tests/test_latent_states.py tests/test_calibration_iteration.py -v` passed: 6 tests.
- Combined focused suite passed: `pytest tests/test_boundaries.py tests/test_compartment_strand_deconv.py tests/test_background_model.py tests/test_boundary_model.py tests/test_boundary_sweep.py tests/test_latent_states.py tests/test_calibration_iteration.py -v` passed: 23 tests.
- VS Code diagnostics reported no errors in edited Phase III files/tests.

Deviations from plan:

- `calibration_e_step()` returns a `CalibrationStepResult` dataclass rather than the tuple sketched in the design. The dataclass carries the sweep, local posterior, log tensor, flags, and `sum_log_evidence`, which are needed for diagnostics and iteration.
- The capture logBF uses boundary excess over the off-target expectation (`max(sweep.mu_sweep - rho_off * contained_leff, 0)`) rather than raw `sweep.mu_sweep`. This avoids assigning no-excess regions `p_captured ~= 0.5` just because `sweep.mu_sweep` includes the background Gamma prior.
- The strand logBF is a conservative summary term based on `RegionGdnaChannelEstimate.contained_rna_lower` and `contained_precision`; exact `strand_log_likelihood_d_grid()` tensor wiring is deferred because the tensor layer does not currently receive raw folded strand counts.
- The scalar M-step refits `rho_off` and updates `capture_enrichment_target` diagnostically, but it carries `kappa_d` through from `strand_channels` instead of refitting it. A live `kappa_d` refit requires payload/folded strand counts and is better completed during Phase IV wiring.
- `RegionCalibration` is not yet added to `CalibrationResult`; this remains Phase IV migration work so existing orchestrator outputs stay unchanged.

Issues encountered:

- No lint, test, or editor diagnostic failures after implementation.
- The exact strand tensor remains the main modeling gap in Phase III. The current summary term is intentionally conservative and tested only for directionality, not exact likelihood calibration.

## Phase IV - Clean Cutover Plan

### 2026-05-26 Design Checkpoint 1

Progress:

- Published the detailed Phase IV plan at `docs/newcal/phase4_clean_cutover_plan.md`.
- Inspected the current live cutover surfaces: `_orchestrator.py`, `_result.py`, `prior.py`, `exposure.py`, `pipeline.py`, and calibration config.
- Identified the old live dependency chain that Phase IV must replace: `fit_density_evidence()` -> `fuse_density_and_strand()` -> `RegionExposure.from_density()` -> `calibration.fused_region_gdna` -> `assemble_priors()`.
- Updated `docs/newcal/new_cal_plan_v6.md` to remove stale `use_v6_calibration` / migration-gate language and describe Phase IV as a clean replacement.

Deviations from plan:

- The earlier v6 plan described a cautious migration with compatibility fields and an internal `use_v6_calibration` flag. This is now superseded by the clean-cutover directive: the new calibration path becomes the only live path, with no old density/fusion/exposure state preserved on `CalibrationResult`.

Issues encountered:

- The native EM still accepts an `enable_gdna` array. The Phase IV plan resolves this by keeping it as a technical likelihood-eligibility bit only: true when a locus has at least one unspliced unit with finite gDNA likelihood, independent of `upper_gdna` or any calibration posterior threshold.
- `density_model.py`, `integration.py`, and `exposure.py` still have tests and may have reusable low-level math. The plan requires removing them from live production state first; physical deletion can be a cleanup checkpoint after production references and tests are migrated.

### 2026-05-26 Implementation Checkpoint 1

Progress:

- Rewrote the live calibration handoff around `RegionCalibration`. `CalibrationResult` now carries `region_calibration`, final-pass `background_model`, `boundary_local`, `boundary_sweep`, `strand_channels`, and the EM `prior_table` diagnostics; the old live `density_evidence`, `region_gdna`, `region_exposure`, and `fused_region_gdna` state was removed from the result contract.
- Rebuilt `_orchestrator.py` around the v6 flow: build FL models once, construct region arrays/ledger/observation, estimate compartment strand channels, fit `BackgroundModel`, build `BoundaryTable`, run `run_calibration_iteration()`, and return the v6-only result.
- Cut `assemble_priors()` over to `RegionCalibration.mu_gdna`, `upper_gdna`, `rna_lower`, and `A_r`. Native `enable_gdna` remains a technical candidate-eligibility bit rather than a calibration posterior threshold.
- Updated pipeline/config/CLI wiring for the clean cutover: removed the old density exposure knob, added `background_trim_fraction` and `max_calibration_passes`, required `calibration.region_calibration` before quantification, and logged v6 calibration summaries.
- Updated simulation/benchmark helper reporting to read v6 `region_calibration.rho_off`.
- Updated focused tests and intentionally refreshed golden outputs after the clean cutover behavior change.

Deviations from plan:

- Kept `density_model.py`, `integration.py`, and the old exposure module available for low-level tests/imports, but removed them from the live calibration result/API exports and production downstream path. Physical deletion is deferred cleanup rather than part of this checkpoint.
- Added a smooth RNA-lower-bound likelihood penalty in `assemble_priors()`: loci with strong region-level RNA lower-bound evidence and little plausible gDNA inflate the gDNA EM effective length denominator. This is not a gDNA on/off gate; `enable_gdna` still reflects only native technical eligibility, and likelihood-supported gDNA remains available.
- Reduced the v6 background Gamma alpha floor from `1.0` to `1.0e-3` so sparse no-evidence background regions do not manufacture a whole-count gDNA prior.
- Adjusted native EM gDNA initialization so the gDNA component does not receive the RNA VBEM Jeffreys baseline. With zero calibration prior, gDNA warm-start mass is assigned only when its likelihood weight beats RNA; otherwise it remains technically eligible but unsupported.

Issues encountered:

- Full-suite validation after the initial clean cutover found a real no-gDNA regression: 10 scenario tests failed because pure RNA fragments leaked large mass into gDNA when `enable_gdna` became technical eligibility only.
- Diagnostic scripts showed two contributing causes: the background alpha floor produced positive `mu_gdna` in no-evidence regions, and many pure-RNA unspliced fragments had finite or locally favorable gDNA likelihoods, allowing EM to siphon mass even with tiny priors.
- The fix preserves the clean cutover and avoids restoring old `fused_region_gdna` / posterior-threshold gating by moving suppression into continuous prior/effective-length scaling.

Validation:

- Native rebuild after C++ changes: `pip install --no-build-isolation -e .` succeeded.
- Focused leakage and true-gDNA checks passed: `tests/test_native_gdna_eligibility.py`, representative no-gDNA scenario tests, `tests/scenarios/test_gdna_diagnosis.py`, `tests/scenarios/test_single_exon.py::TestSingleExon::test_gdna_sweep`, and `tests/scenarios/test_wide_intron.py`.
- Focused calibration/native suite passed after final changes: `tests/test_calibration_prior.py`, `tests/test_calibrate.py`, `tests/test_calibration_result.py`, `tests/test_pipeline_wiring.py`, `tests/test_rna_lower_confidence.py`, and `tests/test_native_gdna_eligibility.py`.
- Ruff passed on `src/` and `tests/`.
- Golden outputs were intentionally updated and reverified: `tests/test_golden_output.py` passed.
- Final full suite passed: `pytest tests/ -v` reported `1087 passed in 26.95s`.

### 2026-05-26 Implementation Checkpoint 2 - Cleanup To Baseline

Progress:

- Reverted the provisional no-gDNA leakage suppression added in Checkpoint 1.
- Restored native EM initialization to the clean baseline: gDNA again receives the mode-aware baseline prior when technically eligible, and warm-start mass is not filtered by a gDNA-vs-RNA likelihood heuristic.
- Removed the hardcoded RNA-lower-bound gDNA likelihood penalty from `assemble_priors()`. `gdna_eff_len` is again only the FL/geometric gDNA effective length multiplied by `RegionCalibration.A_r`.
- Restored the background Gamma `alpha_floor` default from `1.0e-3` to `1.0`, matching the pre-suppression v6 baseline.
- Kept the Phase IV cutover itself intact: `CalibrationResult` remains v6-only, `assemble_priors()` consumes `RegionCalibration`, and native `enable_gdna` remains technical eligibility only (`candidate_enable`) rather than `candidate_enable & (gdna_upper > eps)`.
- Refreshed golden outputs to this clean baseline so they no longer encode the discarded suppression patch.

Deviations from plan:

- The provisional smooth RNA-lower penalty, background pseudo-count reduction, and native EM warm-start changes from Checkpoint 1 are explicitly abandoned. They solved tests but were not adopted as modeling policy.
- No replacement gDNA aggressiveness knob was added in this checkpoint. That work should be designed separately as a bidirectional Bayesian prior problem, likely with explicit RNA and gDNA prior gravity.

Issues encountered:

- The clean baseline still overestimates gDNA in zero-gDNA sentinel cases. This is now treated as an open operating-point/modeling issue rather than an implementation bug to patch around.
- Representative zero-gDNA severity after cleanup:
	- `test_pure_mrna_baseline`: expected T1 = 1000, observed T1 = 774, gDNA = 226.
	- `test_single_exon::test_baseline`: expected t1 = 339, observed t1 = 80, diff = 259.
	- `test_wide_intron::test_baseline_no_false_nrna`: expected t1 = 500, observed t1 = 244, diff = 256.
- The likely underlying issue remains asymmetric prior gravity: RNA components have no comparable prior mass while any positive gDNA prior and competitive gDNA likelihood can pull ambiguous unspliced fragments into gDNA.

Validation:

- Native rebuild after reversion: `pip install --no-build-isolation -e .` succeeded.
- Ruff passed on `src/` and `tests/`.
- Focused clean-cutover suite passed: `tests/test_calibration_prior.py`, `tests/test_native_gdna_eligibility.py`, `tests/test_calibrate.py`, `tests/test_calibration_result.py`, and `tests/test_pipeline_wiring.py` (`21 passed`).
- Clean-baseline goldens were updated and reverified: `tests/test_golden_output.py` passed (`21 passed`).
- Representative zero-gDNA sentinel scenarios intentionally fail under this baseline and should remain visible until the prior/aggressiveness design is settled.

### 2026-05-26 Documentation Checkpoint - Current gDNA Prior Construction

Progress:

- Published `docs/newcal/gdna_prior_current_construction.md` as a descriptive audit of the current calibration-to-EM gDNA prior path.
- Documented that the live per-locus EM gDNA prior count is the geometrically allocated `RegionCalibration.mu_gdna` count, not simply an arbitrary Python-side pseudocount.
- Documented the separate denominator/exposure path: FL-marginal locus gDNA effective length multiplied by the bp-weighted regional `A_r` exposure multiplier and floored at `1.0`.
- Documented where the independent floors and baselines actually enter: background Gamma `alpha_floor = 1.0` / `beta_floor = 1.0`, exposure minimum `1.0e-4`, native VBEM baseline `0.5`, and numerical guards.
- Documented the native warm-start and final prior vector so future prior design can distinguish initialization behavior from modeled prior mass.

Deviations from plan:

- None. This is documentation only and does not change the clean baseline.

Issues encountered:

- The audit reinforces the same open modeling issue from Checkpoint 2: current calibration adds explicit extra prior count only to gDNA, while RNA has no matched calibrated prior gravity at the projection layer.
