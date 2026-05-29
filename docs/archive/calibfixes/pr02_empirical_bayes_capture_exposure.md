# PR 2: One-Class Source-Reliable Capture Exposure

> **Superseded by** `pr02_capture_exposure_v3.md`. This file is retained as
> the v1/v2 planning trail and should not be implemented directly.

## Goal

Fix `A_r` by estimating local capture exposure from source-reliable gDNA mass, not from latent calibration states. PR02 should work for both ordinary RNA-seq and hybrid-capture RNA-seq without requiring any external BED or probe file.

The first implementation should be deliberately simple:

```text
unit exposure = source-reliable gDNA mass / off-target gDNA expectation
A_r           = 1 + local target weight * positive unit excess exposure
```

No panel file, no seed hyperparameters, and no mandatory captured/off-panel mixture model are needed to fix the current bug.

## What Is Broken Today

The root failure is in `calibration_iteration.py`:

```text
mu_gdna = p_captured * captured_mu + (1 - p_captured) * off_target_mu
A_r     = mu_gdna / (rho_off * contained_leff)
gamma_r = captured_mu / (rho_off * contained_leff)
```

This is wrong because `p_captured` is a latent stratum posterior, not a gDNA source posterior. The four latent states describe expression/capture context:

```text
unexpressed_offtarget
unexpressed_capture
expressed_capture
expressed_offtarget
```

They do not say whether the molecules are RNA or gDNA. When expressed captured RNA pushes a region toward a capture stratum, the current code can turn that into higher `mu_gdna`, then into higher `A_r`, then into larger gDNA EM effective length. That is a source/exposure feedback loop built on the wrong variable.

The concrete fixes are:

1. Stop deriving `A_r` from latent-state-weighted `mu_gdna`.
2. Stop using density fallback `mu_gdna` as exposure evidence when no source-reliable channel exists.
3. Estimate exposure only from PR01 source-reliable gDNA mass.
4. Pool exposure at MultiLocus/gene scale, not fine-region scale.
5. Keep `A_r = 1` anywhere source exposure is not identifiable.

Everything else is secondary.

## Parent Contracts

This PR keeps the v5 invariants:

- latent states are expression/capture strata, not RNA/gDNA source labels;
- source mass is explicit and comes from `PriorMassDeconvolution.gdna_unspliced_mean` only when that mass is source-reliable;
- `A_r` is a sampling-opportunity multiplier, not a local abundance ratio from latent states;
- RNA prior mass remains grouped, with no transcript-level RNA floor;
- reliability and exposure diagnostics must travel with the fitted exposure.

## Design Decisions

### No external panel files

Remove all external BED/probe concepts from PR02 production behavior. Do not add these config fields or CLI flags:

```text
capture_probes
capture_probe_format
capture_gdna_split_penalty
capture_min_overlap
capture_validation_panel
capture_validation_panel_format
```

Benchmark comparison against a synthetic truth panel belongs in benchmarking analysis, not in `rigel quant` calibration input.

### No seed hyperparameters

Do not add seed thresholds such as:

```text
capture_seed_quantile
capture_seed_min_excess
capture_seed_min_enrichment
capture_min_unit_opportunity
```

Those knobs are symptoms of the two-class design trying to discover a panel before it has a stable exposure estimator. PR02 should not need to seed captured loci. It should estimate a continuous exposure weight for each unit from source-reliable gDNA/background ratios.

### No `ann` variable name

Use `target_weight` for the annotation-derived surface:

```text
T_{r,c} = target_weight for region r, compartment c
```

`T` does not mean "captured" and does not predict the panel. It only says this compartment is targetable enough for local capture exposure to matter. Annotation is a proxy for targetability; source-reliable gDNA decides exposure.

### One-class before two-class

Use a one-class exposure model in PR02. A true two-class captured/off-panel mixture is not required for the root fix.

A two-class model answers: "which loci belong to an inferred panel?"

PR02 only needs to answer: "given source-reliable gDNA mass, how much larger is local gDNA sampling opportunity than the off-target background?"

Those are different problems. The second one is enough to repair `A_r` and works naturally for both capture and non-capture samples:

- non-capture RNA-seq: source-reliable gDNA tracks background, so `A_r ~= 1`;
- capture RNA-seq with gDNA: target units have source-reliable gDNA excess, so `A_r > 1`;
- RNA-rich regions with no source-reliable gDNA: `A_r = 1`;
- unstranded capture-on: source split is still deferred, so PR02 should not invent exposure from density states.

A two-class model can be revisited later if benchmarks show that continuous one-class weights are too noisy or need stronger target-panel reporting. It should not block the root fix.

## Configuration Contract

Keep production configuration minimal:

```python
capture_exposure_mode: Literal["auto", "off"] = "auto"
```

Meanings:

- `auto`: compute source-reliable one-class exposure where the source channel is identifiable; otherwise return neutral `A_r = 1`.
- `off`: force neutral `A_r = 1`.

CLI flags:

```text
--capture-exposure-mode {auto,off}
```

No other PR02 tuning knobs are required.

## Exposure Unit

Fit exposure at a locus-level unit `u`:

1. MultiLocus if the region/locus machinery already groups transcripts into a shared EM neighborhood;
2. otherwise isolated annotated gene;
3. otherwise no capture unit.

Fine regions are observations, not exposure units. A fine region can be noisy or RNA-dominated. A locus/gene unit has enough source-reliable mass for a stable background-normalized ratio.

Each region maps to one exposure unit for PR02. Regions with no unit get neutral exposure.

## Target Weight

Use `T_{r,c}` as a deterministic targetability weight in `[0, 1]`:

```text
contained exon compartment:       T_{r,contained} = 1 for targetable exonic unit regions
boundary compartments:            T_{r,boundary} in [0, 1] when boundary evidence belongs to the unit
intronic/intergenic compartments: T_{r,c} = 0 unless boundary-compatible with the unit
```

This is not a prediction that the region was captured. It only gates where local capture exposure is allowed to modify gDNA effective length.

For the first implementation, `T` can be simple:

- contained exonic parts of a MultiLocus/gene unit get `1`;
- boundary evidence adjacent to that unit contributes to fitting the unit exposure;
- unrelated intronic/intergenic regions get `0` for region-level `A_r`.

## Evidence Aggregation

For each region `r`, compartment `c`, and unit `u`:

```text
B_{r,c} = rho_off * leff_{r,c}                  # off-target gDNA expectation
G_{r,c} = source-reliable gDNA mean from PR01   # not latent density fallback
T_{r,c} = target_weight
```

Aggregate to unit `u`:

```text
B_u = sum_{r,c in u} T_{r,c} * B_{r,c}
G_u = sum_{r,c in u} T_{r,c} * G_{r,c}
```

Then compute the one-class excess exposure:

```text
lambda_u_raw = G_u / B_u

theta_u = max(lambda_u_raw - 1, 0)
lambda_u = 1 + theta_u
```

If `B_u <= 0`, if the unit has no source-reliable gDNA channel, or if `capture_exposure_mode == "off"`:

```text
theta_u = 0
lambda_u = 1
```

There is no seed set and no fitted captured class.

## Source-Reliable Means Only

This point is critical. `G_{r,c}` must come from PR01 source reliability, not from latent density fallback.

Allowed evidence:

- stranded source deconvolution means weighted by PR01 reliability;
- boundary source means when PR01 marks them source-informative;
- future strand-free FL/boundary source split, once implemented.

Disallowed evidence:

- `p_captured`;
- `p_unexpressed_capture`;
- `p_expressed_capture`;
- `sweep.mu_sweep` by itself;
- `mu_gdna` when it was derived from latent-state mixing;
- raw observed count excess in RNA-rich exons.

This is why PR02 does not solve unstranded capture-on source splitting. Without a source-reliable `G`, exposure must stay neutral.

## Region Exposure

Build region exposure from the unit ratio:

```text
A_r = 1 + T_{r,contained} * theta_{unit(r)}
gamma_r = A_r
```

Boundary compartments contribute to fitting `theta_u`, but the region-level `A_r` used for EM gDNA effective length should initially use contained target weight only. This keeps the first implementation easy to reason about.

No unit, no target weight, or no source reliability gives:

```text
A_r = 1
gamma_r = 1
```

## One-Class Versus Two-Class

Do not implement a captured/off-panel mixture in PR02 unless the one-class model fails benchmarks.

The one-class model is enough because the current failure is not panel classification. The current failure is that exposure is being estimated from the wrong object. Once `A_r` is a source-reliable background ratio, non-capture samples naturally return neutral exposure.

A two-class model adds complexity and new failure modes:

- it needs seeding or sparse-mixture initialization;
- it can hallucinate a captured class in ordinary RNA-seq if residuals are structured;
- it introduces extra diagnostics before the basic exposure ratio is proven;
- it answers panel discovery, not the minimum `A_r` repair.

Panel discovery should be a report derived from the continuous one-class exposures, not an internal state that determines `A_r`.

## Learned Panel Output

Still emit a diagnostic file:

```text
inferred_capture_panel.bed
+ inferred_capture_panel.tsv
```

This output is derived after exposure fitting. It must not feed back into calibration.

Recommended TSV columns:

```text
unit_id
unit_name
unit_kind
chrom
start
end
strand
B_u
G_u
lambda_u
theta_u
source_reliability_mass
n_regions
```

For the BED score, use a deterministic continuous score such as:

```text
score = scaled log1p(theta_u) with zero score when theta_u = 0
```

No threshold is needed for correctness. The report can include all units with scores, or only units with `theta_u > 0` if we want a compact BED.

## New Module: `capture_exposure.py`

Add flags:

```python
CAPTURE_EXPOSURE_DISABLED = np.uint16(0x1)
CAPTURE_EXPOSURE_NO_UNIT = np.uint16(0x2)
CAPTURE_EXPOSURE_NO_TARGET_WEIGHT = np.uint16(0x4)
CAPTURE_EXPOSURE_NO_SOURCE_RELIABILITY = np.uint16(0x8)
CAPTURE_EXPOSURE_NO_BACKGROUND = np.uint16(0x10)
CAPTURE_EXPOSURE_NEUTRAL = np.uint16(0x20)
```

Add dataclasses:

```python
@dataclass(frozen=True, slots=True)
class CaptureUnitMap:
    unit_id: np.ndarray                    # int32[R], -1 when no unit
    unit_names: tuple[str, ...]
    unit_kind: np.ndarray                  # uint8[U], multilocus/gene/etc.
    contained_target_weight: np.ndarray    # float32[R], in [0, 1]
    boundary_left_target_weight: np.ndarray
    boundary_right_target_weight: np.ndarray

@dataclass(frozen=True, slots=True)
class CaptureExposureFit:
    A_r: np.ndarray                        # float32[R]
    gamma_r: np.ndarray                    # float32[R]
    region_theta: np.ndarray               # float32[R]
    region_target_weight: np.ndarray       # float32[R]
    unit_theta: np.ndarray                 # float64[U]
    unit_lambda: np.ndarray                # float64[U]
    unit_background: np.ndarray            # float64[U], B_u
    unit_gdna_source: np.ndarray           # float64[U], G_u
    unit_source_reliability_mass: np.ndarray
    flags: np.ndarray                      # uint16[R]
    unit_flags: np.ndarray                 # uint16[U]
```

Summary diagnostics:

```text
n_units
n_units_neutral
n_units_no_source_reliability
n_units_exposure_gt_1
A_r summary
unit_lambda summary
unit_background summary
unit_gdna_source summary
source_reliability_mass summary
```

## Integration Recipe

Current code derives exposure after `prior_mass` is built:

```python
denominator = rho_off * contained_leff
exposure = _safe_exposure(mu_gdna, denominator)
captured_exposure = _safe_exposure(captured_mu, denominator)
```

Replace that with:

```python
capture_fit = fit_capture_exposure(
    observation=observation,
    background=background,
    prior_mass=prior_mass,
    strand_channels=strand_channels,
    unit_map=capture_unit_map,
    mode=config.capture_exposure_mode,
)

A_r = capture_fit.A_r
gamma_r = capture_fit.gamma_r
```

Important implementation details:

1. Build `CaptureUnitMap` once after `region_arrays`, loci, and boundaries are available.
2. Call `fit_capture_exposure()` after `build_prior_mass_deconvolution()`.
3. Inside `fit_capture_exposure()`, ignore `PriorMassDeconvolution` rows whose method/source does not represent source-reliable evidence.
4. Keep `PriorMassDeconvolution` itself unchanged unless PR02 needs clearer method flags.
5. Store `capture_fit` on `CalibrationStepResult` and final `RegionCalibration` for diagnostics.
6. Delete `_safe_exposure()` if no other call site remains.

`prior.py` can keep consuming `region_calibration.A_r` via `bp_weighted_mean_exposure_over_blocks()`. The key change is the provenance of `A_r`.

## Files To Edit

| File | Purpose |
|---|---|
| `src/rigel/config.py` | Add `capture_exposure_mode` only. |
| `src/rigel/cli.py` | Add `--capture-exposure-mode {auto,off}`. |
| `src/rigel/calibration/capture_exposure.py` | New unit mapping, one-class exposure fit, diagnostics, learned panel writer helpers. |
| `src/rigel/calibration/calibration_iteration.py` | Replace latent-state ratio exposure with `fit_capture_exposure()`. |
| `src/rigel/calibration/_orchestrator.py` | Build/pass `CaptureUnitMap`. |
| `src/rigel/calibration/_result.py` | Add compact capture-exposure diagnostics. |
| `src/rigel/pipeline.py` or output writer | Emit `inferred_capture_panel.bed` and TSV diagnostics. |
| `src/rigel/calibration/prior.py` | Keep consuming `A_r` only as exposure. |
| `scripts/benchmarking/runner.py` | Do not pass probe files to `rigel quant`. |
| `tests/test_capture_exposure.py` | Unit tests for one-class ratios and source-reliability gating. |
| `tests/test_per_locus_gdna_mass.py` | Verify `A_r` affects gDNA effective length as exposure only. |

## Tests

Add focused tests for:

1. ordinary capture-off sample: `G_u ~= B_u` gives `A_r ~= 1`;
2. no source-reliable channel: `A_r == 1` even if latent capture probability is high;
3. expressed captured RNA with no PR01 gDNA source evidence does not inflate `A_r`;
4. strong source-reliable local gDNA excess gives `lambda_u = G_u / B_u`;
5. negative/local depleted ratios are floored at neutral exposure, not below `1`;
6. boundary source evidence contributes to unit `G_u` and `B_u` but does not create independent units;
7. fine-region splitting does not change unit exposure when summed to the same locus;
8. `capture_exposure_mode="off"` returns all-one exposure;
9. `prior.py` multiplies gDNA effective length by the new one-class `A_r`;
10. learned panel BED/TSV is derived from `theta_u` and does not affect `A_r`.

## Validation Commands

```bash
conda activate rigel && ruff check src/rigel/config.py src/rigel/cli.py src/rigel/calibration tests/test_capture_exposure.py tests/test_per_locus_gdna_mass.py
conda activate rigel && pytest tests/test_capture_exposure.py tests/test_per_locus_gdna_mass.py tests/test_calibration_iteration.py tests/test_cli.py -v
```

Broader checks before merge:

```bash
conda activate rigel && pytest tests/test_calibrate.py tests/test_calibration_result.py tests/test_pipeline_wiring.py tests/test_golden_output.py -v
```

Golden outputs should be updated only when changed diagnostics or quantification are intentional and explained.

## Benchmark Gate

After PR01 + PR02, rerun the hybrid-capture suite without providing probe files as calibration inputs.

Required outcomes:

- capture-off conditions: `A_r` remains effectively `1` except for source-reliable local gDNA excess;
- no-gDNA capture-on: RNA depth does not inflate exposure;
- high-gDNA stranded capture-on: target loci with source-reliable excess get elevated `A_r`;
- high-gDNA capture-off: no broad false exposure inflation;
- unstranded capture-on: no material regression; source split remains deferred;
- learned panel diagnostics recover known synthetic captured genes/loci when compared externally by benchmark analysis.

If these fail because one-class ratios are too noisy, then design a follow-up two-class or shrinkage model using the observed failure mode. Do not preemptively add seed knobs.

## Review Checklist

- [ ] No external BED/probe input remains in PR02 production calibration.
- [ ] No seed hyperparameters were introduced.
- [ ] No `ann` variable naming remains; use `target_weight`/`T`.
- [ ] `A_r` no longer comes from `mu_gdna / denominator`.
- [ ] `A_r` no longer uses latent capture-state probabilities.
- [ ] Density fallback source mass does not fit exposure.
- [ ] Exposure unit is MultiLocus/gene, not fine region.
- [ ] Off-target/no-source units return neutral `A_r = 1`.
- [ ] Boundary source evidence contributes to unit exposure fitting.
- [ ] `prior.py` consumes `A_r` only as exposure, not source mass.
- [ ] Learned panel output is diagnostic only and does not affect calibration.
