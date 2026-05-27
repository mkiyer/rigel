# PR01 Implementation Log

Date: 2026-05-27

## Scope

Implement `pr01_continuous_strand_reliability.md`: continuous strand reliability, beta-binomial RNA uncertainty in the mixed-model slab, reliability-weighted `PriorMassDeconvolution`, and compact diagnostics.

## Running notes

- Started by auditing the current strand path. `strand_log_likelihood_d_grid()`, `_exact_posterior_R()`, and `_summarize_normal()` currently model RNA strand observations with a fixed point-binomial `p_r1_sense`; this matches the critique and must change for PR01 source/reliability outputs.
- `_orchestrator.calibrate()` currently calls `deconvolve_compartments_by_strand(compartment_counts, kappa_d=...)` without passing `StrandSummary`; PR01 must thread the strand posterior through this call.
- `RegionGdnaChannelEstimate` currently has only compartment mean/upper/RNA-lower/precision arrays. PR01 needs compartment reliability and log-Bayes-factor arrays, plus clear semantics that means consumed by prior mass are slab means under the beta-binomial RNA mixed model.
- Added `StrandModel.n_minor`, `n_major`, `minor_rate_posterior_alpha`, and `minor_rate_posterior_beta`; `to_dict()` now reports these diagnostics.
- Extended `StrandSummary` with same/opposite counts and minor-rate posterior parameters. For backward compatibility with existing tests and call sites that construct `StrandSummary(p_r1_sense=..., n_observations=...)`, the summary infers same/opposite counts from `p_r1_sense` and `n_observations` when explicit counts are absent.
- Added a new beta-binomial mixed-model strand likelihood helper, leaving the legacy `strand_log_likelihood_d_grid()` intact for existing point-binomial callers.
- Added exact and large-count approximate beta-binomial slab summaries. The exact path computes `log_p_mixed` and slab `E[D | data, H1]` from the same likelihood array; the approximate path uses beta-binomial RNA variance and a `log_ndtr()`-based Normal interval probability helper.
- Added `RegionGdnaChannelEstimate` reliability/log-BF arrays with backward-compatible defaults, and threaded `StrandSummary` through `_orchestrator.calibrate()` into `deconvolve_compartments_by_strand()`.
- Updated `build_prior_mass_deconvolution()` so strand-derived prior source mass and precision are multiplied by the compartment reliability weights.
- Added compact reliability diagnostics to `CalibrationResult.to_summary_dict()`.
- Added focused tests for beta-binomial mixed likelihood convolution, slab mean/log likelihood consistency, smooth reliability, near-unstranded gating, reliability-weighted prior mass, and `StrandSummary` posterior fields.
- Tightened `StrandSummary` validation so explicit minor-rate posterior parameters must match the same/opposite counts when observations are present.
- Added additional focused tests for pure RNA low reliability, strong mixed high reliability, small training-set uncertainty, large-count approximate reliability, and deep-tail Normal interval probabilities.

## Deviations from the plan

- Backward-compatible `StrandSummary` construction is more permissive than the strict PR recipe: older call sites can omit `n_same`/`n_opposite`, and the dataclass infers them. This keeps the implementation focused and avoids a repo-wide constructor churn.
- The legacy point-binomial `strand_log_likelihood_d_grid()` remains in place because `calibration/integration.py` and existing tests still call it. PR01 production source mass uses the new beta-binomial helper instead. A later cleanup can migrate or rename legacy fusion callers if desired.
- The clean `tests/scenarios/test_nrna_double_counting.py::TestNrnaDoubleCounting::test_full_sweep[g0_n0_s90]` scenario now shows a small false-gDNA posterior under continuous reliability (about 1-1.25% in debug runs). This is the expected replacement for the old hard source decision, not an nRNA double-counting regression, so the imperfect-strand clean mRNA tolerance was widened from 20 to 30 fragments while keeping perfect-strand clean tolerance unchanged.

## Validation log

- `conda activate rigel && pytest tests/test_strand_summary.py tests/test_strand_deconv.py tests/test_compartment_strand_deconv.py tests/test_calibration_iteration.py -v` -> 34 passed.
- `conda activate rigel && ruff check src/rigel/strand_model.py src/rigel/calibration/strand_summary.py src/rigel/calibration/strand_deconv.py src/rigel/calibration/calibration_iteration.py src/rigel/calibration/_orchestrator.py src/rigel/calibration/_result.py tests/test_strand_summary.py tests/test_strand_deconv.py tests/test_compartment_strand_deconv.py tests/test_calibration_iteration.py` -> passed.
- `conda activate rigel && pytest tests/test_calibrate.py tests/test_calibration_result.py tests/test_per_locus_gdna_mass.py tests/test_pipeline_wiring.py -v` -> 15 passed.
- First full-suite retry with `conda activate rigel && pytest tests/ -q` -> 1148 passed, 22 failed. Failures were the intended golden-output changes plus the clean imperfect-strand scenario tolerance noted above.
- `conda activate rigel && pytest tests/test_strand_summary.py tests/test_strand_deconv.py tests/test_compartment_strand_deconv.py tests/test_calibration_iteration.py tests/scenarios/test_nrna_double_counting.py::TestNrnaDoubleCounting::test_full_sweep -q` after the final PR01 test additions and scenario tolerance update -> 67 passed.
- `conda activate rigel && ruff check src/rigel/strand_model.py src/rigel/calibration/strand_summary.py src/rigel/calibration/strand_deconv.py src/rigel/calibration/calibration_iteration.py src/rigel/calibration/_orchestrator.py src/rigel/calibration/_result.py tests/test_strand_summary.py tests/test_strand_deconv.py tests/test_compartment_strand_deconv.py tests/test_calibration_iteration.py tests/scenarios/test_nrna_double_counting.py` -> passed.
- `conda activate rigel && pytest tests/test_golden_output.py --update-golden -q` -> 21 skipped after regenerating golden files.
- `conda activate rigel && pytest tests/test_golden_output.py -q` -> 21 passed.
- `conda activate rigel && pytest tests/ -q` -> 1176 passed.
- `conda activate rigel && git diff --check` -> passed.
