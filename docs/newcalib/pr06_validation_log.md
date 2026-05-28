# PR06 Validation Log

## 2026-05-28: Zero-gDNA Calibration Mass Audit

### Trigger

The `anti_intron_ss065_nrna50` no-gDNA scenario inferred nonzero calibration gDNA mass:

- Simulated gDNA: `0`
- `RegionUnsplicedMass`: `325.517 / 1720 = 18.93%` gDNA
- Exposure learned from that mass: `omega_p50=8.742`, `omega_p95=40.077`, `omega_max=50.934`

### Root Cause

The mass was produced upstream of EM in `build_region_unspliced_mass()` tier 1. The production PR03 bridge was using the strand-only mixed slab from `deconvolve_compartments_by_strand()` as unconditional regional gDNA mass. In low-strand, nRNA-heavy zero-gDNA regions, normal RNA strand fluctuations created positive mixed-slab gDNA expectations. Because the handoff did not honor deterministic-zero density evidence, those local strand fluctuations became `RegionUnsplicedMass.gdna_mass`.

Diagnostic comparison confirmed this was a bridge regression, not fundamental evidence for gDNA: the existing density/strand fusion path returned `rho_ref=0`, `rho_ref_source=ZERO`, `mean_sum=0.000`, `upper_sum=0.000` for the same case.

### Fix

- `calibrate()` now computes `fit_density_evidence()` after `build_density_observation()`.
- If density evidence returns `rho_ref_source == "ZERO"`, `run_calibration_iteration()` receives `force_zero_gdna_mass=True`.
- `build_region_unspliced_mass()` now honors that flag by returning `gdna_mass=0`, `rna_mass=total_mass`, `METHOD_BACKGROUND_FALLBACK`, zero precision, and `FLAG_M_IMPUTED_FROM_BACKGROUND`.
- Exposure fitting now only allows high-confidence unexpressed rows (`p_unexpressed >= 0.80`) into the tau2 pool and neutralizes omega for all non-pool rows. This prevents expressed/nRNA rows from learning EM exposure even when raw diagnostics carry ratios.

### Post-Fix Diagnostic

For `anti_intron_ss065_nrna50`:

- `RegionUnsplicedMass`: `0 / 1720 = 0.00%` gDNA
- Exposure: `tau2=0`, `method=no_pool_neutral`, `omega=1` for all regions
- The raw strand deconvolution still reports noisy local slab means, which is useful diagnostically, but they no longer become calibration mass without density support.

### Validation

- `pytest tests/test_region_unspliced_mass.py tests/test_exposure.py tests/test_calibration_integration.py -q` -> `40 passed`
- `pytest tests/scenarios/test_antisense_intronic.py tests/scenarios_aligned/test_paralogs.py -q` -> `31 passed`
- Refreshed golden outputs intentionally with `pytest tests/test_golden_output.py --update-golden -q`
- `pytest tests/test_golden_output.py -q` -> `21 passed`
- Focused calibration/output/scenario suite -> `135 passed`
- Ruff on touched files -> passed
- Full suite: `pytest tests/ -q` -> `1206 passed`
- Benchmark smoke:
  - `python -m scripts.benchmarking status -c scripts/benchmarking/configs/default.yaml` reached config parsing but could not enumerate conditions because `/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks/runs` is absent on this machine.
  - `python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml --dry-run` failed for the same missing `runs/` directory before command generation.

### Notes

The old PR05 conditional xfails in the antisense and aligned paralog scenarios were removed after the strict assertions passed. The same diagnostic still reports some final EM gDNA assignment in the no-gDNA case; that is downstream likelihood competition and is separate from the fixed calibration-mass/exposure bug.