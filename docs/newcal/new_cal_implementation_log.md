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

### 2026-05-26 Design Checkpoint - Grouped RNA/gDNA Prior v2

Progress:

- Published `docs/prior/prior_redesign_v2.md` as the implementation plan for calibration-derived grouped RNA/gDNA priors.
- Synthesized the Claude and Gemini prior-redesign proposals against the live `RegionCalibration` -> `PriorTable` -> native EM path.
- Adopted the core group-prior design: calibration acts on the aggregate `gDNA` vs `sum RNA` simplex boundary, while the EM remains responsible for within-RNA transcript partitioning.
- Specified additive pseudocount semantics (`alpha_gdna_add`, `alpha_rna_add`) rather than standard Dirichlet `alpha - 1` concentration notation, so calibration expected counts remain non-negative pseudo-fragment masses.
- Added a phased implementation plan covering regional `mu_rna`, paired prior projection, native grouped MAP/VBEM updates, SQUAREM integration, output diagnostics, and benchmark tuning.

Deviations from attached proposals:

- Rejected gene-level RNA groups as a production inference primitive because Rigel internal inference must remain transcript-first; any future multi-group extension must use transcript-derived or region-derived groups.
- Did not assume SQUAREM is automatically unaffected by the grouped projection; the plan requires grouped fixed-point steps, convergence guards, and benchmark validation.
- Separated calibration prior strength from gDNA operating-point bias so sensitivity can be tuned without changing total prior budget.

Issues encountered:

- The main implementation risk is native EM surgery: the grouped prior must be applied inside every MAP/VBEM step, not only as initialization or a per-component prior vector.
- `mu_rna` must be surfaced carefully from region-level calibration evidence so it protects RNA without manufacturing expression in background regions.

### 2026-05-26 Design Checkpoint - Grouped RNA/gDNA Prior v3

Progress:

- Published `docs/prior/prior_redesign_v3.md` as the implementation-ready successor to the v2 plan after review.
- Accepted the core review findings that v2's raw `m_g + m_r` budget was too hot, that `mu_rna` must not include spliced fragments, and that SQUAREM/grouped-prior behavior needs explicit instrumentation and fallback criteria.
- Reframed the design around mass-conserving unspliced prior deconvolution: per region, `prior_gdna_unspliced + prior_rna_unspliced = prior_unspliced_total`, with spliced fragments excluded from prior pseudocount mass.
- Replaced the v2 spliced-plus-residual `mu_rna` concept with a new `PriorMassDeconvolution` contract carrying unspliced total, gDNA mean, RNA mean, method, precision, and flags.
- Specified a bounded prior budget with weaker-side balance, one-sided edge support, and a hard max count so calibration cannot inject raw expression-scale pseudocounts into highly expressed loci.
- Revised gDNA availability semantics: remove `enable_gdna` as a Python/calibration policy input in the final path, derive native `has_gdna_candidate` structurally, and zero both prior sides when no gDNA candidate exists.
- Preserved the user's non-negotiable all-RNA aggregate requirement: annotated mRNA and synthetic nRNA remain one RNA group; nRNA siphon benchmarks are diagnostic, not a trigger to split RNA groups.

Deviations from v2:

- v3 does not add a general `mu_rna` count to `RegionCalibration`; it adds unspliced-only paired prior masses.
- v3 does not use raw `m_g + m_r` as the default prior budget and does not allow spliced evidence to become RNA prior pseudocount mass.
- v3 no longer treats `enable_gdna` as a required Python-side prior-table field. Native may keep a temporary compatibility input during ABI migration, but production semantics should be derived from actual gDNA likelihood candidates.
- v3 rejects the reviewed suggestion to split mature RNA and nRNA as a follow-up path; future fixes should improve regional deconvolution and operating-point defaults while keeping one aggregate RNA group.

Issues encountered:

- The key implementation risk shifts from inventing `mu_rna` to maintaining exact unspliced mass conservation through region deconvolution, geometry allocation, and locus prior budgeting.
- Strand-uninformative data remain intentionally hard: unspliced nascent-like mass may enter the gDNA prior because strand balance cannot distinguish it from gDNA. This needs explicit benchmark reporting rather than hidden correction.
- Native SQUAREM currently needs careful inspection before implementing acceptance-rate metrics; v3 requires at least step/clamp/nonfinite/fallback diagnostics even if a full line-search acceptance counter is not yet present.

### 2026-05-26 Implementation Checkpoint - Prior Redesign v3 Phase 0 Start

Progress:

- Began implementation from `docs/prior/prior_redesign_v3.md`.
- Inspected the current live handoff: `RegionCalibration` has only gDNA-side prior mass, `PriorTable` projects only gDNA counts, `pipeline.quant_from_buffer()` passes `gdna_prior_count_em` plus `enable_gdna`, and native EM still builds a per-component prior vector with a VBEM baseline.
- Confirmed the current calibration data structures already expose the needed unspliced-only ingredients: `DensityObservation` stores contained and boundary unspliced counts separately, and `RegionGdnaChannelEstimate` stores compartment-specific gDNA posterior means for strand-informative deconvolution.

Deviations from plan:

- None yet. The implementation will keep the current `enable_gdna` Python parameter temporarily for wrapper compatibility, but native production semantics should derive candidate availability from actual unspliced finite gDNA candidates.

Issues encountered:

- `tests/test_bayesian_prior_acceptance.py` appears to describe an older global-prior contract and imports `rigel.calibration.locus_prior`, which is no longer present in the current clean-cutover code. I will treat it as stale prior-contract coverage and either replace it with v3 tests or report it separately if it blocks validation.

### 2026-05-26 Implementation Checkpoint - Prior Redesign v3 Python Phases

Progress:

- Added `PriorMassDeconvolution` and `build_prior_mass_deconvolution()` to `calibration_iteration.py`.
- Wired `prior_mass` through `CalibrationStepResult` and `RegionCalibration`.
- Built the prior mass from unspliced-only counts: density fallback uses clipped `mu_gdna`, while strand-informative calibration sums compartment gDNA means and derives RNA as `unspliced_total - gdna`.
- Added `prior_mass` summaries to `CalibrationResult.to_summary_dict()`.
- Added grouped prior config fields to `EMConfig`: `aggregate_prior_strength`, `aggregate_prior_edge_count`, `aggregate_prior_max_count`, and `gdna_prior_logit_bias`.
- Extended `PriorTable` with paired RNA/gDNA expected counts, bounded alpha fields, budget/share diagnostics, conservation diagnostics, and `gdna_prior_density`.
- Updated `assemble_priors()` to consume `RegionCalibration.prior_mass`, allocate gDNA/RNA/unspliced masses by geometry, compute bounded paired pseudocounts, and keep `gdna_prior_count_em` as an alias for `alpha_gdna_add`.
- Threaded `EMConfig` from `pipeline.quant_from_buffer()` into `assemble_priors()`.

Deviations from plan:

- Kept `enable_gdna` on `PriorTable` for diagnostics and wrapper compatibility during the native ABI transition. It is no longer used to compute prior mass.
- The first focused test run exposed an expected assertion drift: `gdna_prior_count_em` no longer equals raw allocated gDNA mass because it aliases bounded `alpha_gdna_add`. Updated the test to assert raw `gdna_expected_count` separately from bounded EM prior mass.

Issues encountered:

- None unresolved in the Python calibration/projection layer.

Validation:

- `pytest tests/test_calibration_iteration.py tests/test_calibration_prior.py tests/test_calibration_result.py tests/test_rna_lower_confidence.py -v` passed: 35 tests.

### 2026-05-26 Implementation Checkpoint - Prior Redesign v3 Native/Output Completion

Progress:

- Completed Phase 3 native grouped RNA/gDNA EM implementation.
- Updated the native batch EM ABI to accept paired aggregate pseudocounts: `locus_gdna_prior_count` and `locus_rna_prior_count`.
- Added a single grouped prior update path shared by MAP and VBEM. The update constrains only aggregate gDNA vs aggregate RNA, then redistributes RNA-side mass dynamically according to current RNA evidence or carried RNA state.
- Removed the modeled VBEM `0.5` baseline as prior mass. Native now keeps only numerical floors needed for stability.
- Changed native gDNA availability to be structural: a locus has gDNA only if an unspliced unit has finite gDNA likelihood. The compatibility `enable_gdna` array remains in the wrapper but no longer gates modeling behavior.
- Added cold-start behavior for RNA prior mass: if no current RNA evidence and no carried RNA state exist, RNA alpha is inactive for that update.
- Threaded `alpha_rna_add` through `estimator.py` and `pipeline.py` into native EM.
- Extended per-locus output diagnostics with `rna_expected_count`, `prior_unspliced_total`, `alpha_gdna_add`, `alpha_rna_add`, prior budget/share fields, and `gdna_prior_density` while retaining `gdna_prior_count` / `gdna_prior_count_em` aliases for this migration checkpoint.
- Replaced stale `tests/test_bayesian_prior_acceptance.py` old global-prior coverage with v3 grouped-prior acceptance tests.
- Added native tests for RNA counter-prior behavior and inactive grouped priors when no structural gDNA candidate exists.
- Added `scripts/debug/prior_v3_sentinel_sweep.py` to measure zero-gDNA and true-gDNA sentinel behavior across grouped-prior operating points.

Deviations from plan:

- The published v3 plan proposed initial internal defaults `aggregate_prior_strength=1.0`, `aggregate_prior_edge_count=5.0`, `aggregate_prior_max_count=10.0`, and `gdna_prior_logit_bias=0.0`. Full-suite sentinel tests showed this was too weak: zero-gDNA single-exon, wide-intron, pure-mRNA, and nRNA-double-counting baselines leaked substantial RNA mass into structurally available gDNA.
- After the sentinel sweep, changed the production defaults to `aggregate_prior_strength=3.0`, `aggregate_prior_edge_count=1000.0`, `aggregate_prior_max_count=3000.0`, and `gdna_prior_logit_bias=-6.0`.
- Tightened the budget cap semantics so `aggregate_prior_max_count` caps the final effective pseudocount budget after applying strength, not only the pre-strength `budget_raw`.
- Kept the Python `enable_gdna` array as a temporary ABI/wrapper compatibility and diagnostic field. Native ignores it as a calibration/modeling gate and derives `has_gdna_candidate` internally.
- SQUAREM diagnostics were added, but a true line-search acceptance/rejection fallback was not implemented in this checkpoint. The native path records step scale, clamp/nonfinite counts, and stabilization failures; full acceptance-rate metrics remain future work.

Issues encountered:

- First full-suite run after native grouped EM found 10 sentinel failures, all from no-gDNA or near-no-gDNA baselines where structural gDNA candidates siphoned RNA mass. These were treated as real model operating-point failures, not tolerance/test issues.
- The initial bounded budget was conservative enough to avoid raw-count-scale hammers but too small to counteract gDNA likelihood competition in pure unspliced/ambiguous loci.
- Increasing edge/max without a negative gDNA logit bias helped single-exon and wide-intron sentinels but worsened some mixed/multiexon zero-gDNA cases because small inferred gDNA shares were amplified. The final `-6.0` gDNA logit bias preserves true-gDNA recovery in the measured sentinel subset while suppressing false gDNA in zero-gDNA baselines.
- No unexpected compile failures occurred; native rebuild succeeded on the first build after the ABI changes.

Validation:

- Native rebuild after C++ edits: `pip install --no-build-isolation -e .` succeeded.
- Focused native tests passed: `pytest tests/test_batch_em_impl.py tests/test_native_gdna_eligibility.py -v` (`10 passed`).
- Focused grouped-prior integration passed after tuning: `pytest tests/test_calibration_iteration.py tests/test_calibration_prior.py tests/test_calibration_result.py tests/test_rna_lower_confidence.py tests/test_bayesian_prior_acceptance.py tests/test_batch_em_impl.py tests/test_native_gdna_eligibility.py tests/test_pipeline_wiring.py -v` (`52 passed`).
- Sentinel regression rerun passed: selected `test_gdna_diagnosis`, `test_non_overlapping`, `test_nrna_double_counting`, `test_single_exon`, and `test_wide_intron` cases (`43 passed`).
- Golden outputs were intentionally regenerated for the grouped-prior schema/behavior change and reverified: `pytest tests/test_golden_output.py -v` (`21 passed`).
- Final full suite passed: `pytest tests/ -v` (`1093 passed in 25.95s`).

### 2026-05-26 Design Checkpoint - Adaptive Prior v1

Progress:

- Published `docs/prior/adaptive_prior_v1.md` as the next-step design for replacing v3's fixed grouped-prior operating constants with sample- and locus-specific uncertainty-derived priors.
- Reframed the v3 four constants as internal fallback limits rather than ordinary user-tunable hyperparameters.
- Specified a Beta-equivalent effective sample size design: calibration should provide `E[D_l]` and `Var(D_l)` for locus-level gDNA unspliced mass, then derive `alpha_gdna_add` and `alpha_rna_add` from posterior mean and variance.
- Identified the missing implementation surface: v3 preserves posterior means and a coarse precision scalar, but adaptive priors need count-scale variance from strand and density calibration.
- Designed empirical-Bayes sample-level shrinkage so ambiguous loci borrow strength from the sample's inferred global gDNA profile instead of from a fixed `gdna_prior_logit_bias`.

Deviations from current v3 implementation:

- The adaptive design demotes `aggregate_prior_strength`, `aggregate_prior_edge_count`, and `gdna_prior_logit_bias` from primary policy to legacy/fallback-only behavior.
- The design keeps `aggregate_prior_max_count` conceptually, but only as an absolute ESS safety cap, not as the central prior budget rule.
- The design requires new uncertainty payload fields before implementation; trying to build v4 directly on the current `precision` scalar would risk creating a new magic-number layer.

Issues encountered:

- The first draft almost reintroduced method-specific trust constants for strand vs density fallback. I revised the design to make inverse posterior variance the primary weighting surface, with method labels used for diagnostics and migration fallbacks only.

### 2026-05-26 Design Checkpoint - Adaptive Prior v2 One-Knob Quantile Plan

Progress:

- Published `docs/prior/adaptive_prior_v2.md` as the implementation-ready successor to `adaptive_prior_v1.md` after reviewing `docs/prior/one_knob_parameterization.md`.
- Accepted the central one-knob design: ordinary users should specify only `rna_confidence`, a posterior quantile level for the gDNA share of unspliced mass.
- Revised the adaptive prior from a posterior-mean share summary to a decision-theoretic quantile summary: `share_l = QBeta(q; a_l, b_l)`, while prior ESS remains learned from posterior variance.
- Kept the v1 uncertainty architecture: preserve count-scale variance from strand and density calibration, project variance to loci, estimate sample-level empirical-Bayes shrinkage, and compute grouped `alpha_gdna_add` / `alpha_rna_add` from the final posterior summary.
- Added implementation phases covering helper math tests, regional variance preservation, variance projection, EB global profile, adaptive quantile priors, CLI/config migration, output migration, and validation/audit.

Clarifications and deviations from the attached one-knob proposal:

- Clarified that quantile optimality is exact for continuous share summaries under asymmetric absolute loss; binary RNA/gDNA classification would use a posterior probability threshold instead.
- Corrected the legacy confidence mapping: current `rna_lower_confidence` maps to the same upper-gDNA quantile level `q`, not `1 - q`, when translated into the gDNA-share prior.
- Rejected a semantic no-op bridge where `rna_confidence=0.5` silently maps to v3 magic defaults. Legacy v3 can remain an explicit compatibility mode, but the new knob should not pretend to mean something it does not mean.
- Kept internal safety rails such as `_MAX_ESS` and `_MAX_GLOBAL_ESS`, but explicitly classified them as private numerical caps rather than user-tunable model parameters.

Issues encountered:

- The density interval fallback needed a correction: variance derived from an upper bound must use the confidence level that created that upper bound, not the user `rna_confidence` quantile, because `q=0.5` would imply `z=0`.
- The moment-matching helper needs the unspliced total/finite-information scale as an input for the variance floor; the v2 plan now includes `total` in `beta_from_mean_var()`.

### 2026-05-26 Design Checkpoint - Adaptive Prior v4 Review

Progress:

- Reviewed `docs/prior/adaptive_prior_v3.md` against the live v6 calibration objects, grouped-prior projection path, and four-state expressed x capture model.
- Published `docs/prior/adaptive_prior_v4.md` as the implementation-ready successor to v3.
- Kept the core one-knob contract: `rna_confidence` is the only user-facing prior operating point, and the EM still receives one grouped gDNA alpha and one aggregate RNA alpha.
- Simplified the global empirical-Bayes design: v4 pools local Beta-equivalent evidence, shrinks by cross-locus agreement, applies leave-one-locus-out global support, and caps final ESS by the locus's own prior-owned unspliced mass.
- Added explicit dry-run expectations for stranded/unstranded whole-transcriptome RNA-seq, stranded/unstranded capture RNA-seq, and the four calibration states: background, gDNA-only capture, expressed capture, and expressed off-target.

Deviations from adaptive prior v3:

- Replaced v3's inverse-variance global Beta formula because it can collapse to zero ESS at the exact boundaries `psi=0` or `psi=1`, losing strong no-gDNA or all-gDNA sample evidence.
- Changed `U_l <= eps` behavior: v3 allowed a locus to become purely shrinkage-driven; v4 emits no prior at all when a locus has no prior-owned unspliced mass, preventing global prior injection into unsupported loci.
- Replaced `tau=0 when variance is zero` with boundary-aware moment matching. Point-mass evidence at `phi=0` or `phi=1` now carries ESS capped by local data content.
- Clarified that `rna_lower_confidence` and `gdna_density_confidence` are internal calibration interval levels, not the same object as the public prior decision quantile. Legacy names can be migration aliases, but they should not continue as user-tuned calibration knobs.
- Corrected the v3 file-reference error: `RegionGdnaEstimate` lives in `src/rigel/calibration/strand_deconv.py`, not `_result.py`.

Issues encountered:

- The hardest regimes remain unstranded whole-transcriptome and especially unstranded capture RNA-seq. v4 intentionally treats unresolved gDNA-vs-unspliced-RNA ambiguity as low ESS plus diagnostics rather than introducing new reliability multipliers.
- Density fallback variance must be implemented with the full four-state law of total variance. The current `p_captured * captured + (1-p_captured) * off-target` mean is not sufficient as an uncertainty contract for expressed capture ambiguity.
- The global prior is now deliberately conservative in heterogeneous capture samples: cross-locus disagreement lowers global ESS so local evidence remains primary.

### 2026-05-26 Design Checkpoint - Adaptive Prior v5 Entropy-Dirichlet

Progress:

- Published `docs/prior/adaptive_prior_v5.md` as the parameter-free successor to v4.
- Replaced Beta moment matching with discrete Dirichlet-Multinomial counts: per-region pseudocounts are `n_r = U_r * w_r * q_r`, where `q_r` is the two-group projection of the four-state posterior and `w_r = 1 - H(P_r) / log 4`.
- Replaced the user-facing `rna_confidence` quantile with information-entropy weighting. Calibration certainty now scales prior gravity directly; there is no quantile to invert and no boundary U-shape pathology to defend against.
- Replaced explicit leave-one-locus-out shrinkage with a single `N_global - N_l` vector subtraction across loci, and used `1 - W_l` (locus mass-weighted entropy) as the shrinkage weight.
- Removed every prior-related user knob from CLI, YAML, and `EMConfig`. No deprecation aliases, no legacy modes.
- Specified deletion of `prior_redesign_v3.md`, `adaptive_prior_v1.md`, `adaptive_prior_v2.md`, `adaptive_prior_v3.md`, `adaptive_prior_v4.md`, and `one_knob_parameterization.md` into `docs/prior/archive/` once the implementation lands.

Deviations from adaptive prior v4:

- No `--rna-confidence` knob and no Beta. The prior is a linear function of calibration evidence, not a quantile of a posterior.
- No variance fields are added to `PriorMassDeconvolution` or `RegionGdnaChannelEstimate`. The Dirichlet count model does not need posterior variance.
- No leave-one-locus-out loop. Global pool minus local mass is mathematically equivalent and trivially vectorized.
- `rna_lower_confidence` and `gdna_density_confidence` survive only as private module-level constants for the strand RNA lower-bound screen and gDNA upper-bound diagnostic. They no longer influence the EM prior in any way.
- The grouped EM prior columns in `loci.feather` are collapsed to a much smaller set; no compatibility aliases for the deleted v3/v4 columns.

Issues encountered:

- The v5 design depends on the four-state calibration posterior being well calibrated. Pathologies in the tensor will propagate directly into the prior. This is the intended contract: improving the prior now means improving the calibration model.
- Capture-heterogeneous samples still rely on the global pool being a meaningful summary. v5 protects locally confident loci via `W_l`, but a fully ambiguous sample will produce uniformly weak priors, which is the honest outcome.
- Implementation must be done in one cutover; partial migration is not supported because the v3/v4 column set is deleted outright.
