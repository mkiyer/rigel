# Strand-Based gDNA / RNA Deconvolution - Implementation Plan v2

Date: 2026-05-23
Status: ready for implementation
Supersedes: `strand_model_impl_plan_v1.md`
Design basis: `rnaseq_mode_aware_gdna_density_plan.md`
Primary use case: strand-specific hybrid-capture RNA-seq

## 0. Why v2 Exists

v1 had the right statistical direction, but the implementation shape was too
large. It proposed a new `strand_deconv/` subpackage with seven modules,
several new dataclasses, per-source kappa objects, exon self-training, exact
posterior grids, FFT convolutions, and multiple confidence levels. That is not
the right first implementation.

The calibration package is already too fragmented. We should not solve a
complexity problem by adding more files. The first implementation should be:

1. one new calibration module;
2. one new region-level output object;
3. one strand-count builder;
4. one kappa estimate reused from existing code;
5. one deconvolution routine;
6. one integration point in `calibrate(...)`;
7. one prior-assembly consumer.

This is enough to unlock the strand-specific capture use case and remove the
current fractional-cutover blocker without turning calibration into a maze.

## 1. Big-Picture Pipeline

The intended mental model is simple:

```text
BAM
  -> C++ scanner
  -> fractional region accumulator: float32[R, 12]
  -> calibration
       - FL models
       - global density diagnostics
       - strand-based region gDNA estimate
  -> locus prior assembly
  -> EM
```

Calibration is allowed to be statistically careful. It is not allowed to become
architecturally ornate. The strand model is one calibration method that consumes
the same region evidence table everyone else consumes.

## 2. Critique Of v1

### 2.1 Too many modules

The v1 module plan:

```text
calibration/strand_deconv/
  observations.py
  kappa_d.py
  seed_training.py
  exon_screen.py
  posterior.py
  tail_approx.py
  result.py
region_gdna.py
```

is overfit to the whiteboard. Most of those names describe helper functions,
not ownership boundaries. They would make future agents search across eight
files to understand one algorithm.

v2 adds exactly one implementation file:

```text
src/rigel/calibration/strand_deconv.py
```

Everything else is additive wiring to existing modules.

### 2.2 Too many parameters

v1 exposed or implied many knobs: multiple confidence levels, exact-grid
thresholds, saddlepoint switches, kappa clipping bounds, minimum source counts,
no-RNA exon thresholds, self-training iteration count, source-specific
conservative overrides, and posterior degeneracy flags.

v2 uses four internal constants and no new CLI/config surface in the first PR:

```python
STRAND_ALPHA = 0.95
MAX_EXACT_POSTERIOR_N = 128
MIN_KAPPA_TRAINING_REGIONS = 20
CONSERVATIVE_KAPPA_FALLBACK = 20.0
```

These are implementation constants, not user-facing parameters. We can tune
them after benchmark evidence. Do not add more until a failing validation case
requires one.

### 2.3 Exon self-training is deferred

The v1 no-RNA exon screen is statistically sensible, but it is not required for
the first working strand model. It adds a loop, new thresholds, spliced-evidence
policy, source-specific diagnostics, and a failure mode where the model trains
on its own mistakes.

v2 does not train on exons. It trains `kappa_d` from high-purity intergenic and
intron-only regions. If those are sparse, it uses a conservative overdispersed
fallback, not an overconfident binomial fallback.

This is the key simplification: missing kappa evidence should make RNA lower
bounds weaker, not make the model more certain.

### 2.4 Exact posterior everywhere is unnecessary

Exact convolution is useful for low counts. It is not worth designing an FFT
batching subsystem before we have a working model. v2 uses:

- exact discrete posterior only for `n <= 128`;
- a moment-matched normal lower-bound approximation for larger `n`.

This is easier to implement, easier to test, and likely fast enough. If real
benchmarks show a calibration problem at high counts, we improve the numerical
method inside the same module.

## 3. Naming Rules

All variables remain transcript-relative:

| Name | Meaning |
|---|---|
| `k_sense` | observed unspliced count on transcript-sense strand |
| `k_antisense` | observed unspliced count on transcript-antisense strand |
| `n` | `k_sense + k_antisense` |
| `p_r1_sense` | P(R1 aligns transcript-sense), from `StrandModel` |
| `R` | unknown RNA-derived unspliced count |
| `D` | unknown gDNA-derived unspliced count, `D = n - R` |
| `kappa_d` | beta-binomial concentration for gDNA strand imbalance |

No `major`, `minor`, `primary`, `secondary`, or protocol-rotated channel names
in code. When a signed scalar is needed, use:

```python
strand_sign = 1.0 if p_r1_sense >= 0.5 else -1.0
signed_excess = strand_sign * (k_sense - 0.5 * n)
strand_contrast = abs(p_r1_sense - 0.5)
```

This preserves sense/antisense names while still handling R1-sense and
R1-antisense protocols correctly.

## 4. Minimal Code Changes

### 4.1 Add one module

Create:

```text
src/rigel/calibration/strand_deconv.py
```

This module owns the full v2 strand algorithm:

```python
@dataclass(frozen=True, slots=True)
class StrandRegionCounts:
    k_sense: np.ndarray
    k_antisense: np.ndarray
    n_total: np.ndarray
    eligible: np.ndarray
    p_r1_sense: float

@dataclass(frozen=True, slots=True)
class RegionGdnaEstimate:
    mean_count: np.ndarray
    upper_count: np.ndarray
    rna_lower_count: np.ndarray
    precision: np.ndarray
    flags: np.ndarray
    kappa_d: float
    p_r1_sense: float
    alpha: float
```

That is the complete new schema. No separate `StrandKappaFit`, no per-source
stat dataclasses, no nested result package.

### 4.2 Reuse existing modules

Use existing code instead of adding new owners:

| Need | Existing owner |
|---|---|
| transcript-strand class | `fractional_evidence.transcript_strand_class` via `RegionArrays.ts_class` |
| sense/antisense projection | `fractional_evidence.sense_antisense_split` |
| region-count sorted view | `PayloadArrays.region_counts_sorted` |
| gDNA strand beta-binomial MoM | `density_global.estimate_strand_balance` |
| strand protocol estimate | `strand_summary.StrandSummary` / scanner `StrandModel` |
| result hand-off | `_result.CalibrationResult` |
| prior consumer | `locus_prior.py` replacement path |

Do not extend `PayloadArrays` with spliced fields for v2. The existing
`region_counts_sorted` already contains all 12 channels. The strand module can
call `sense_antisense_split` directly on that matrix.

### 4.3 Add only two `CalibrationResult` fields

Add optional fields to `_result.CalibrationResult`:

```python
region_gdna: RegionGdnaEstimate | None = None
strand_region_counts: StrandRegionCounts | None = None
```

`strand_region_counts` is diagnostic/debug surface and can be dropped later if
it is not useful. `region_gdna` is the production hand-off.

Do not change output file schemas until prior assembly consumes the estimate.
The first PR should preserve current golden outputs except for `summary.json`.

## 5. v2 Algorithm

### 5.1 Build strand-region counts

Function:

```python
def build_strand_region_counts(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    *,
    p_r1_sense: float,
) -> StrandRegionCounts:
    ...
```

Implementation:

1. For each unspliced compartment (`contained`, `boundary_left`,
   `boundary_right`), call `sense_antisense_split(...)` on
   `payload_arrays.region_counts_sorted` and `region_arrays.ts_class`.
2. Sum `sense` across the three compartments into `k_sense`.
3. Sum `antisense` across the three compartments into `k_antisense`.
4. Set `n_total = k_sense + k_antisense`.
5. Set `eligible = (ts_class == TS_POS or ts_class == TS_NEG) and n_total > 0`.

Boundary evidence is not given a special policy in v2. The fractional
accumulator already assigned boundary mass to region compartments. The strand
model consumes that accounting as-is.

### 5.2 Train gDNA strand overdispersion

Function:

```python
def estimate_kappa_d(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
) -> StrandBalanceEstimate:
    ...
```

Training set:

```text
intergenic regions OR intron-only regions
contained unspliced POS/NEG counts only
```

Rationale:

- intergenic and intron-only contained regions are the cleanest available gDNA
  sources;
- using contained counts avoids boundary-policy debates in the first version;
- exon self-training is deferred;
- capture sparsity is handled by a conservative fallback.

Implementation:

1. Build mask: `is_intergenic(signature) | is_intron_only(signature)`.
2. Call existing `estimate_strand_balance(...)` with
   `contained_unspliced_pos`, `contained_unspliced_neg`, and the mask.
3. If the returned estimate has fewer than `MIN_KAPPA_TRAINING_REGIONS`, use
   `CONSERVATIVE_KAPPA_FALLBACK` and mark the summary as fallback.
4. If the estimate falls back because residual variance is below binomial
   expectation, keep the high-kappa result. Clean data may legitimately look
   binomial.
5. If training is sparse, do not use high-kappa fallback. Sparse evidence must
   weaken RNA lower bounds.

This reuses the existing statistical implementation and fixes the main v1
conservatism issue.

### 5.3 Deconvolve one region

Model:

```text
D = n - R
DNA_sense | D ~ BetaBinomial(D, mean=0.5, concentration=kappa_d)
RNA_sense | R ~ Binomial(R, p_r1_sense)
k_sense = DNA_sense + RNA_sense
```

Moment estimate:

```text
rna_mean = clip((k_sense - 0.5 * n) / (p_r1_sense - 0.5), 0, n)
gdna_mean = n - rna_mean
```

Identifiability guard:

```text
if abs(p_r1_sense - 0.5) < STRAND_CONTRAST_NUMERICAL_FLOOR:
    rna_mean = 0
    rna_lower = 0
    gdna_mean = n
    gdna_upper = n
    precision = 0
```

### 5.4 Lower-bound method

Function:

```python
def deconvolve_regions_by_strand(
    counts: StrandRegionCounts,
    *,
    kappa_d: float,
    alpha: float = STRAND_ALPHA,
) -> RegionGdnaEstimate:
    ...
```

For `n <= MAX_EXACT_POSTERIOR_N`, compute the exact discrete posterior over
`R = 0..n` by direct convolution. Keep this implementation simple and readable;
no FFT, no batching, no native code.

For `n > MAX_EXACT_POSTERIOR_N`, use a normal approximation with the variance
from the model:

```text
Var(DNA_sense | D) = 0.25 * D * (D + kappa_d) / (1 + kappa_d)
Var(RNA_sense | R) = R * p_r1_sense * (1 - p_r1_sense)
```

The normal path finds the smallest RNA count whose model can explain the
observed sense/antisense imbalance at confidence `alpha`. It should be
conservative: when in doubt, lower `rna_lower_count`, raise `gdna_upper_count`,
and lower `precision`.

Output:

```text
rna_lower_count = lower alpha-confidence bound on RNA-derived count
gdna_upper_count = n - rna_lower_count
mean_count = gdna_mean
upper_count = gdna_upper_count
```

v2 emits only 95% bounds. Add 99% later only if downstream prior calibration
needs it.

### 5.5 Precision

Use one simple precision proxy:

```text
precision = strand_contrast^2 * n / (1 + overdispersion_factor)
```

where:

```text
strand_contrast = abs(p_r1_sense - 0.5)
overdispersion_factor ~= 1 + n / (1 + kappa_d)
```

Then normalize to `[0, 1]`:

```text
precision = precision / (1 + precision)
```

This is intentionally not a second statistical model. It is a monotone weight:
more fragments and stronger strandedness increase confidence; higher DNA
overdispersion decreases confidence.

### 5.6 Flags

Keep only flags that change downstream behavior or debugging:

```python
FLAG_INELIGIBLE = 1 << 0
FLAG_NEAR_UNSTRANDED = 1 << 1
FLAG_KAPPA_FALLBACK = 1 << 2
FLAG_APPROX_NORMAL = 1 << 3
```

No low-count, posterior-degenerate, training-seed, no-RNA-exon, or source flags
in v2. If the exact posterior underflows, the code should switch to the normal
path and set `FLAG_APPROX_NORMAL`.

## 6. Integration Plan

### M0 - Calibration cleanup guardrail

Before adding the algorithm, make a small cleanup PR that adds a package-level
comment or docstring rule in `src/rigel/calibration/__init__.py`:

```text
New calibration algorithms must start in one module. Split only after a second
independent algorithm needs the same helper.
```

This is small but important. It stops the current pattern of module-per-idea.

No broad renames in M0. Renaming 19 files before landing functionality would
create risk without helping users.

### M1 - Add `strand_deconv.py`

Add:

- `StrandRegionCounts`
- `RegionGdnaEstimate`
- `build_strand_region_counts(...)`
- `estimate_kappa_d(...)`
- `deconvolve_regions_by_strand(...)`

Keep all helper functions private in the same file:

- `_exact_region_posterior(...)`
- `_normal_region_bound(...)`
- `_beta_binomial_pmf(...)`
- `_binomial_pmf(...)`

Acceptance:

- no new package directory;
- no new config object;
- one test file covers the public functions.

### M2 - Wire calibration result

Modify `_orchestrator.calibrate(...)`:

```text
region_arrays = RegionArrays.from_region_df(...)
payload_arrays = PayloadArrays.from_payload(...)
strand_summary = provided or StrandSummary.from_model(scan_trained.strand_model)
counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=...)
kappa = estimate_kappa_d(region_arrays, payload_arrays)
region_gdna = deconvolve_regions_by_strand(counts, kappa_d=kappa.kappa)
return build_calibration_result(..., region_gdna=region_gdna, strand_region_counts=counts)
```

Modify `_result.CalibrationResult` and `build_calibration_result(...)` to carry
these optional fields and summarize:

```json
"strand_deconv": {
  "enabled": true,
  "alpha": 0.95,
  "p_r1_sense": 0.99,
  "kappa_d": 123.4,
  "n_eligible_regions": 12345,
  "n_approx_normal": 678,
  "kappa_fallback": false
}
```

Acceptance:

- `calibrate(...)` returns a result with `region_gdna` populated for
  strand-specific runs;
- existing tests pass except any test intentionally pinned to the old summary
  schema.

### M3 - Consume in prior assembly

Replace the current fractional-cutover fail-fast in `locus_prior.py` with the
smallest useful consumer:

1. for each locus footprint, find overlapping region rows;
2. sum `region_gdna.mean_count` into the locus gDNA prior count;
3. sum `region_gdna.upper_count` for diagnostics;
4. use `precision` to set prior strength;
5. leave density/exposure fallback unchanged or conservative for ineligible
   regions.

Do not redesign all locus-prior machinery in this PR. The goal is to remove the
blocker and let EM run with strand-informed gDNA priors.

Acceptance:

- `FractionalCutoverPending` no longer blocks strand-specific calibration;
- EM can run end-to-end on a small strand-specific fixture;
- prior diagnostics show the strand-derived gDNA contribution by locus.

### M4 - Validate before adding features

Run focused tests and one benchmark before implementing exon self-training,
capture exposure, or source-specific kappa.

Validation cases:

1. R1-sense synthetic: high `k_sense`, low `k_antisense` -> positive RNA lower
   bound.
2. R1-antisense synthetic: low `k_sense`, high `k_antisense` -> same RNA lower
   bound after reflection.
3. gDNA-only synthetic: 50/50 split -> RNA lower bound usually zero.
4. sparse kappa training: conservative fallback -> weaker RNA lower bounds than
   high-kappa fit.
5. near-unstranded: all RNA lower bounds zero.

Acceptance:

- protocol symmetry holds;
- true gDNA count is below `upper_count` at least 95% of the time in simulation;
- true RNA count is above `rna_lower_count` at least 95% of the time;
- no runtime explosion on a representative payload.

## 7. Package Cleanup Roadmap

The cleanup should be real, but staged after the first working strand path.
Do not combine large file moves with statistical implementation.

### C1 - Stop adding private one-off modules

Effective immediately:

- no new `_foo.py` for a helper;
- no new subpackage for one algorithm;
- no dataclass unless it crosses a module boundary or is serialized.

### C2 - Collapse dead fractional-cutover stubs

After M3 removes the prior blocker, delete or fold the modules that only exist
to preserve shapes during the cutover. Candidates to audit:

```text
_locus_n_obs.py
_region_index_py.py
_regional_exposure.py
_exposure.py
density_loco.py
regions.py
```

Each deletion must be test-driven: remove one module, run the smallest related
tests, then continue.

### C3 - Target package shape

Long term, calibration should read like this:

```text
calibration/
  __init__.py
  orchestrator.py       # top-level calibrate()
  result.py             # CalibrationResult only
  payload.py            # scan payload + sorted arrays
  evidence.py           # signatures, masks, sense/antisense helpers
  fl.py                 # fragment-length calibration
  gdna.py               # global density + strand deconv + region estimates
  prior.py              # locus/multilocus prior assembly
```

That is eight conceptual files, not nineteen-plus. The transition can happen
with compatibility shims after functionality is restored.

## 8. Tests

Use current repo style: tests live directly under `tests/`, not under a new
nested test package.

Add one main test file:

```text
tests/test_strand_deconv.py
```

Required tests:

1. `test_build_strand_region_counts_ts_pos_and_ts_neg`
   - verifies `sense_antisense_split` semantics are preserved.
2. `test_deconvolution_protocol_symmetry`
   - `p` vs `1 - p` with reflected counts gives identical outputs.
3. `test_gdna_only_has_zero_or_small_rna_lower_bound`
   - 50/50 gDNA regions do not produce confident RNA.
4. `test_rna_excess_produces_positive_lower_bound`
   - strong protocol-consistent imbalance produces RNA lower bound.
5. `test_sparse_kappa_training_uses_conservative_fallback`
   - sparse training does not become high-confidence binomial.
6. `test_near_unstranded_returns_conservative_output`
   - `p ~= 0.5` returns zero RNA lower bound and all-gDNA upper bound.
7. `test_exact_and_normal_paths_are_monotone`
   - larger protocol-consistent imbalance never lowers `rna_lower_count`.

Add one integration test only after M2:

```text
tests/test_calibration_strand_deconv_integration.py
```

It should build a tiny synthetic payload, call `calibrate(...)`, and assert that
`CalibrationResult.region_gdna` is populated with finite arrays.

## 9. What v2 Does Not Do

These are intentionally out of the first implementation:

- no exon self-training;
- no capture BED ingestion;
- no capture exposure model;
- no source-specific kappa dataclasses;
- no 99% or 99.9% bounds;
- no FFT posterior engine;
- no native C++;
- no broad package renames in the same PR as the algorithm.

Each can be added later if benchmarks justify it. The first implementation has
to be small enough that a human can hold it in their head.

## 10. Immediate Implementation Checklist

1. Add `src/rigel/calibration/strand_deconv.py`.
2. Add `tests/test_strand_deconv.py`.
3. Wire optional `region_gdna` into `_result.CalibrationResult`.
4. Wire `calibrate(...)` to build counts, estimate kappa, and deconvolve.
5. Add integration test.
6. Run the focused tests:

```bash
conda activate rigel && pytest tests/test_strand_deconv.py tests/test_calibration_strand_deconv_integration.py -v
```

7. Only after those pass, connect `RegionGdnaEstimate` to the smallest working
   fractional prior path.

## 11. Guiding Rule

If an implementation step wants a new module, a new public dataclass, or a new
parameter, ask one question first:

```text
Can the first working strand-specific capture implementation be correct without it?
```

If yes, defer it. This model is not conceptually complicated. The code should
make that obvious.
