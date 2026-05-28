# PR 05 v3 - Fair Component Exposure for EM

## Status

This plan supersedes [pr05_impl_plan_v2.md](pr05_impl_plan_v2.md).

The v2 plan applied `RegionExposure.omega` only to the gDNA denominator. That is the wrong model for
capture libraries. Exposure is learned from gDNA because gDNA is a stable molecular source for
estimating capture/accessibility/mappability variation, but the exposure process itself is not
specific to gDNA. A bait sees sequence, not molecule class. If a locus has 100x exposure, fragments
from gDNA, mature RNA, and nascent RNA sequence that overlaps that exposed interval are all more
visible.

PR 05 therefore applies exposure fairly to every EM component: each transcript component gets a
component-specific exposure-adjusted EM effective length, and each locus gDNA component gets its own
exposure-adjusted EM effective length. The additive Bayesian prior mass split remains separate.

User review resolved the main open choices: keep TPM's existing meaning, add exposure-adjusted TPM
only as a diagnostic, use bp-weighted component exposure factors for PR05, require region-table
coverage for all indexed sequence, and keep transcript exposure factors stored separately from
effective lengths until the EM input packing step.

## 1. Goal

Implement a clean, inspectable path from calibration to native EM inputs:

```text
RegionUnsplicedMass.gdna_mass/rna_mass  ->  additive Dirichlet mass split
RegionExposure.omega                    ->  per-component EM effective lengths
```

The outputs to native EM should be:

```text
alpha_gdna_add[locus]       from regional gDNA mass
alpha_rna_add[locus]        from regional RNA mass
transcript_eff_len_em[t]    transcript effective length adjusted by transcript exposure
gdna_eff_len_em[locus]      locus gDNA effective length adjusted by locus exposure
```

No native C++ ABI change is expected: the native solver already receives one `t_eff_lens` array and
one `gdna_eff_len` array.

## 2. Corrected Model

### 2.1 Exposure Is Source-Agnostic

`omega_r` is a regional observation-scale factor:

```text
omega_r = regional sampling visibility / global sampling visibility
```

It is estimated from gDNA mass because gDNA abundance is comparatively stable. Once estimated,
`omega_r` applies to every component whose molecule geometry overlaps region `r`.

Do not describe `omega` as a gDNA-only factor in downstream code, comments, tests, or output.

### 2.2 Exposure Adjusts EM Effective Lengths Symmetrically

Native EM uses component effective lengths as denominators:

```text
log weight_c = log(theta_c) + log_likelihood_c - log(eff_len_em_c)
```

For PR 05, the component EM effective length is:

```text
eff_len_em_c = max(eff_len_unweighted_c * exposure_factor_c, 1.0)
```

This applies to:

```text
transcript components: c = transcript t
gDNA component:        c = gDNA shadow for multi-locus l
```

High exposure increases the EM effective length of every component that overlaps the exposed
sequence. If every component in a locus has the same exposure factor, posterior competition is
unchanged because the same scale is applied to all competitors. That is the desired fairness
property. Exposure should not give gDNA a private boost.

### 2.3 What PR 05 Is Not Doing

PR 05 does not add a per-fragment numerator term such as `+log(omega_at_fragment)`. A fully local
exposure likelihood would use both a numerator and a denominator:

```text
log_lik_fragment_component += log(local_fragment_exposure)
log_lik_fragment_component -= log(integral_exposure_over_component)
```

That is a larger model and would require fragment-to-region exposure lookup for both RNA and gDNA.
PR 05 instead implements a component-level denominator approximation that is symmetric across source
classes and cheap enough for production.

## 3. Current Post-PR04 State

Production after PR 04:

- `RegionCalibration.region_exposure` contains per-region `omega`.
- `RegionUnsplicedMass` contains regional `gdna_mass`, `rna_mass`, and `total_mass`.
- `compute_adaptive_prior()` already projects regional mass to per-locus additive alpha arrays.
- `AbundanceEstimator` already has separate transcript effective length arrays:
  - `effective_lengths`: public/output effective lengths.
  - `em_effective_lengths`: EM-only effective lengths passed to native EM.
- Native EM already receives:
  - `t_eff_lens_arr`: per-transcript EM effective lengths.
  - `locus_gdna_eff_lens`: per-locus gDNA EM effective lengths.
- Current PR04 wiring keeps all exposure-neutral by setting transcript EM lengths equal to public
  lengths and gDNA EM length equal to unweighted gDNA length.

## 4. Target Contracts

### 4.1 Mass Prior Contract

Keep the adaptive prior mass split unchanged:

```text
alpha_gdna_add <- RegionUnsplicedMass.gdna_mass
alpha_rna_add  <- RegionUnsplicedMass.rna_mass
```

The existing `p_unexpressed` gate, overlap projection, ESS cap, structural gDNA gate, RNA call bias,
and locus ESS shrinkage remain the prior machinery.

`RegionExposure.omega` must not be passed into `compute_adaptive_prior()` in PR 05.

### 4.2 Component Exposure Contract

Introduce a component exposure/effective-length table, separate from adaptive prior math:

```python
@dataclass(frozen=True, slots=True)
class ComponentExposureTable:
    transcript_eff_len_unweighted: np.ndarray
    transcript_exposure_factor: np.ndarray
    transcript_eff_len_em: np.ndarray

    gdna_eff_len_unweighted: np.ndarray
    gdna_exposure_factor: np.ndarray
    gdna_eff_len_em: np.ndarray
    gdna_eff_len_adjustment_ratio: np.ndarray
```

Keep the exposure factors and effective lengths unraveled. `transcript_exposure_factor` and
`gdna_exposure_factor` are first-class arrays, not temporary intermediates. The EM effective-length
arrays are derived from them only after validation and numerical clipping.

Recommended module placement:

```text
src/rigel/calibration/_exposure.py      low-level region-overlap exposure helpers
src/rigel/calibration/prior.py          adaptive alpha mass + EM input assembly
```

If `prior.py` grows too mixed, add `src/rigel/calibration/em_inputs.py` and make
`assemble_priors()` a compatibility wrapper around the new helper during PR05.

### 4.3 PriorTable Contract

Update `PriorTable` to carry only locus-level denominator diagnostics needed by downstream output,
plus the existing prior diagnostics:

```python
@dataclass(frozen=True, slots=True)
class PriorTable:
    alpha_gdna_add: np.ndarray
    alpha_rna_add: np.ndarray
    enable_gdna: np.ndarray

    gdna_eff_len_em: np.ndarray
    gdna_eff_len_unweighted: np.ndarray
    gdna_exposure_factor: np.ndarray
    gdna_eff_len_adjustment_ratio: np.ndarray

    # Existing adaptive-prior diagnostics remain.
```

The transcript arrays should be passed through `ComponentExposureTable` or directly installed into
`AbundanceEstimator` before EM. They do not need to live in `PriorTable` unless that is the smallest
implementation path.

Delete from final PR05 production contracts:

```text
gdna_em_exposure_weight
gdna_eff_len_weight_ratio
```

### 4.4 Output Columns

Transcript and nRNA outputs already expose EM-only effective lengths. Rename or repurpose the ratio
column so the sign is obvious:

```text
effective_length                 unexposed/public effective length
em_effective_length              exposure-adjusted EM effective length
em_exposure_factor               em_effective_length / effective_length
```

Replace the current `em_exposure_weight` output name with `em_exposure_factor` if schema churn is
acceptable in PR05. No backward compatibility is required by the roadmap.

Keep conventional TPM columns exposure-neutral, but add side-by-side diagnostic exposure-adjusted
TPM columns for benchmarking/debugging. These diagnostics use `em_effective_length` as the length
denominator and the same count pool as the matching regular TPM column:

```text
tpm                             existing unexposed TPM
tpm_em_exposure                 diagnostic TPM using em_effective_length
tpm_total_rna                   existing total-RNA-normalized TPM, where present
tpm_total_rna_em_exposure       diagnostic total-RNA TPM using em_effective_length
```

The diagnostic TPM columns are not the primary abundance contract in PR05.

Locus output columns:

```text
gdna_eff_len_unweighted
gdna_exposure_factor
gdna_eff_len_em
gdna_eff_len_adjustment_ratio
```

where:

```text
gdna_eff_len_adjustment_ratio = gdna_eff_len_em / gdna_eff_len_unweighted
```

For high exposure this ratio is greater than 1, not less than 1.

### 4.5 Summary JSON

Update the CLI summary block to expose the corrected sign. The current CLI has only a
`_locus_series_summary()` helper, so PR05 should either add a small generic dataframe-series helper
or add a transcript-specific helper for `quant_df` before using the transcript summary below:

```python
"em_exposure": {
    "transcript_factor": _quant_series_summary("em_exposure_factor"),
    "gdna_factor": _locus_series_summary("gdna_exposure_factor"),
    "gdna_eff_len_em": _locus_series_summary("gdna_eff_len_em"),
    "gdna_eff_len_unweighted": _locus_series_summary("gdna_eff_len_unweighted"),
    "gdna_eff_len_adjustment_ratio": _locus_series_summary(
        "gdna_eff_len_adjustment_ratio"
    ),
}
```

Keep the exact JSON shape modest; the important point is that summaries no longer describe exposure
as gDNA-only.

## 5. Exposure Aggregation Algorithm

### 5.1 Use Base-Pair Weighted Region Overlap In PR 05

Do not implement the v2 FL start-window exposure helper in PR05. It is more complex than needed and
risks a Python-side performance regression in mega-loci.

For a component with genomic blocks `B_c`, compute:

```text
covered_bp_c = sum_r overlap_bp(B_c, region_r)
weighted_bp_c = sum_r overlap_bp(B_c, region_r) * omega_r
exposure_factor_c = weighted_bp_c / covered_bp_c
```

Then:

```text
eff_len_em_c = max(eff_len_unweighted_c * exposure_factor_c, 1.0)
```

This is a deliberate PR05 approximation: exposure is averaged across the component geometry by base
pairs. It avoids fragment-to-region mapping and avoids looping over possible fragment starts.

### 5.2 Component Geometry

Use the molecule geometry for each component:

- Mature transcript components use their exon intervals.
- Annotated single-exon/nRNA transcripts use their single exon interval.
- Synthetic nRNA transcripts use their synthetic single-exon span already reconstructed in
  `TranscriptIndex._t_exon_intervals`.
- gDNA components use the merged `MultiLocus.loci` intervals.

This means mature mRNA is not charged for intronic exposure unless the intronic sequence is part of
that transcript component. Nascent RNAs need no special case because they are represented as
transcripts with their own blocks.

### 5.3 Region Coverage Is Required

Do not hard-code missing regions to `omega = 1.0` for normal production paths.

The region table is expected to tile every indexed reference, including intergenic sequence. A locus
or transcript block that is not covered by `RegionArrays` indicates an index/calibration contract
problem. PR05 should validate coverage before native EM and fail with a controlled Python error in
strict production mode. It should not silently substitute neutral exposure for intergenic sequence.

This addresses the v2 gap: intergenic/off-target regions should already have estimated exposure
factors learned from the calibration model.

The strict coverage rule is not permission to send unsafe numbers into native EM. At the final Python
boundary where exposure factors are converted to C++ effective lengths, enforce finite positive
values and apply a defensive floor:

```text
MIN_EXPOSURE_FACTOR = 1.0e-4
exposure_factor_c = max(exposure_factor_c, MIN_EXPOSURE_FACTOR)
eff_len_em_c = max(eff_len_unweighted_c * exposure_factor_c, 1.0)
```

If a narrow non-strict recovery path is ever needed for malformed indexes, it must be explicit,
diagnostic-flagged, and finite (for example neutral factor `1.0` for uncovered bases). The default
PR05 production path should raise before native EM rather than recover silently.

### 5.4 Efficient Sweep Implementation

Add a reusable helper in `_exposure.py`:

```python
@dataclass(frozen=True, slots=True)
class RegionWeightedExposure:
    exposure_factor: np.ndarray
    covered_bp: np.ndarray
    weighted_bp: np.ndarray
    flags: np.ndarray


def component_bp_weighted_exposure(
    *,
    block_ref_ids: np.ndarray,
    block_starts: np.ndarray,
    block_ends: np.ndarray,
    component_ids: np.ndarray,
    n_components: int,
    region_arrays: RegionArrays,
    omega: np.ndarray,
    strict_coverage: bool = True,
    min_exposure_factor: float = 1.0e-4,
) -> RegionWeightedExposure:
    ...
```

The helper should process sorted blocks and sorted regions by reference. Complexity should scale with
actual block-region overlaps, not `n_components * n_regions`.

Implementation sketch:

1. Validate `omega.shape == region_arrays.start.shape`, all values finite and positive.
2. Validate blocks have `end > start` and valid ref ids.
3. Sort blocks by `(ref_id, start, end, component_id)`.
4. For each reference, sweep region intervals and component blocks.
5. For every overlap, accumulate:

```python
np.add.at(covered_bp, component_ids, overlap_bp)
np.add.at(weighted_bp, component_ids, overlap_bp * max(omega_region, min_exposure_factor))
```

6. Track uncovered block bases. In strict mode, raise `ValueError` if any block base is uncovered.
7. Return `exposure_factor = max(weighted_bp / covered_bp, min_exposure_factor)`; for genuinely
  zero-width components, return neutral factor `1.0` with a diagnostic flag.

### 5.5 Transcript Blocks

Build flat transcript block arrays from `TranscriptIndex._t_exon_intervals`:

```text
component_id = transcript index
block_ref_id = transcript ref id
block_start/block_end = exon coordinates
```

Use this once per quantification run to compute:

```text
transcript_exposure_factor[t]
transcript_eff_len_em[t] = max(transcript_eff_len_unweighted[t] * factor[t], 1.0)
```

Store `transcript_exposure_factor` separately from `transcript_eff_len_em` so diagnostics can show
whether a transcript changed because of exposure or because of the unweighted FL effective length.

Avoid a per-locus transcript loop. The flat block sweep over exon intervals is the intended efficient
algorithm.

If `_t_exon_intervals` is not available for a loaded index, expose a small `TranscriptIndex` method
to return exon CSR/block arrays from the loaded interval table rather than reconstructing from the
GTF.

### 5.6 gDNA Locus Blocks

Build flat gDNA block arrays from `multi_loci`:

```text
component_id = multi_locus_id
block_ref_id/start/end = each Locus in MultiLocus.loci
```

Compute:

```text
gdna_exposure_factor[locus]
gdna_eff_len_unweighted[locus] = existing gdna_eff_len_for_loci(...)
gdna_eff_len_em[locus] = max(gdna_eff_len_unweighted[locus] * factor[locus], 1.0)
```

Keep `gdna_eff_len_for_loci()` as the unweighted FL-marginal gDNA opportunity helper. Do not create
an FL-window exposure helper in PR05.

## 6. EM Input Assembly

Add a focused helper that assembles exposure-adjusted EM lengths independently of prior alpha mass:

```python
def compute_component_exposure_table(
    *,
    index: TranscriptIndex,
    multi_loci: list[MultiLocus],
    transcript_eff_len_unweighted: np.ndarray,
    gdna_fl: FragmentLengthModel,
    region_arrays: RegionArrays,
    region_exposure: RegionExposure,
) -> ComponentExposureTable:
    ...
```

Refactor `assemble_priors()` or its successor into clear phases:

```text
1. Build RegionArrays once.
2. Compute native gDNA eligibility.
3. Compute adaptive alpha mass from RegionUnsplicedMass.
4. Compute ComponentExposureTable from RegionExposure.
5. Return prior table plus exposure table to pipeline.
```

Recommended minimal API for PR05:

```python
@dataclass(frozen=True, slots=True)
class EMInputTable:
    prior: PriorTable
    exposure: ComponentExposureTable
```

If changing callers to `assemble_em_inputs()` is too broad for one PR, keep `assemble_priors()` as
the public function but have it return enough fields for pipeline to install transcript EM lengths
and pass gDNA EM lengths. The implementation should still keep the mass-prior helper and exposure
helper separate.

## 7. Pipeline Wiring

### 7.1 Geometry Setup

Current `_setup_geometry_and_estimator()` computes unexposed RNA effective lengths and constructs
`AbundanceEstimator` with `effective_lengths_em=None`.

PR05 should either:

1. Compute transcript exposure before estimator construction and pass `effective_lengths_em`, or
2. Construct the estimator as today, then install `estimator._t_eff_len_em` from the exposure table
   before calling native EM.

Preferred implementation:

```text
compute unexposed transcript effective lengths
compute exposure table once multi_loci are known
set estimator._t_eff_len_em = exposure.transcript_eff_len_em before EM
pass exposure.gdna_eff_len_em to _run_locus_em_partitioned()
```

This keeps scoring and fragment routing unchanged and changes only EM normalization.

### 7.2 `_run_locus_em_partitioned()`

Rename Python-level parameters to make the native input explicit:

```python
gdna_eff_len_em: np.ndarray
gdna_eff_len_unweighted: np.ndarray | None = None
gdna_exposure_factor: np.ndarray | None = None
gdna_eff_len_adjustment_ratio: np.ndarray | None = None
```

Inside `_call_batch_em()`, continue passing to the native wrapper as:

```python
gdna_eff_len=batch_gdna_eff_len_em
```

The native wrapper name can remain `gdna_eff_len` because it is the native ABI name.

### 7.3 Locus Metadata

Update `_build_locus_meta()` to write:

```text
gdna_eff_len_unweighted
gdna_exposure_factor
gdna_eff_len_em
gdna_eff_len_adjustment_ratio
```

Delete `gdna_em_exposure_weight` and `gdna_eff_len_weight_ratio`.

### 7.4 Transcript and nRNA Output

`AbundanceEstimator.get_counts_df()` and `get_nrna_counts_df()` already output both unexposed and
EM effective lengths. Update the ratio name:

```text
em_exposure_factor = em_effective_length / effective_length
```

Keep public `effective_length` and conventional TPM denominators unexposed in PR05 unless a separate
analysis decision is made. `em_effective_length` is the EM normalization surface.

## 8. Tests

### 8.1 Exposure Helper Tests

Add focused tests for `component_bp_weighted_exposure()`:

1. Uniform `omega = 1` returns factor 1 for every component.
2. Whole-component `omega = 50` returns factor 50.
3. Mixed component:

```text
component blocks: [0, 1000)
regions: [0, 100) omega=10, [100, 1000) omega=1
expected factor = (100*10 + 900*1) / 1000 = 1.9
```

4. Multi-block component accumulates all blocks.
5. Mature transcript exons ignore high-exposure intronic regions not in the exon blocks.
6. Synthetic nRNA single-exon span includes intronic regions because its transcript block spans them.
7. Intergenic region exposure is used when gDNA locus blocks overlap intergenic regions.
8. Missing region coverage raises a controlled Python `ValueError` in strict mode; no `KeyError`, no
  C++ call, and no silent neutral fallback.

### 8.2 Component Effective Length Tests

Add tests for `compute_component_exposure_table()`:

1. Transcript EM effective length equals unweighted length when exposure factor is 1.
2. Transcript EM effective length is `unweighted * factor` for high/low factors.
3. gDNA EM effective length is `unweighted * gdna_factor`, not `unweighted / gdna_factor`.
4. `gdna_eff_len_adjustment_ratio` is greater than 1 for high exposure and less than 1 for low
   exposure unless the floor is active.
5. Alpha arrays from `compute_adaptive_prior()` are unchanged when only `omega` changes.

### 8.3 Fairness / Sign Sentinel Tests

Add a native or wrapper-level EM test that catches the v2 bug:

1. Build a locus with two transcript components and one gDNA component, all compatible with the same
  unspliced evidence.
2. Run EM with unweighted lengths.
3. Run EM again after multiplying both transcript and gDNA EM effective lengths by the same factor
  (for example `omega = 500`).
4. Assert all posterior fragment assignments and the transcript/gDNA split are unchanged within
  tolerance.

This test fails if exposure is applied only to gDNA.

Add a second test where transcript and gDNA exposure factors differ because their component
geometries differ. Assert the posterior shift follows the component denominator ratio.

Add a mixed-baiting simulation-style test: exons have high exposure, introns have low exposure, and
the locus contains mature transcript, synthetic/annotated nRNA, and gDNA candidates. The test should
confirm that mature transcript exposure is computed from exon blocks while nRNA/gDNA exposure tracks
the broader genomic span. This documents the intended PR05 behavior for exon-targeted capture.

Add a heterogeneous-exposure diagnostic test: a long multi-exon transcript where only one exon is
highly exposed. Assert the component factor is the bp-weighted average and document that PR05
intentionally smooths local exposure spikes until a future fragment-local likelihood PR.

### 8.4 Pipeline / Output Tests

Update `tests/test_pipeline_wiring.py` so the mock prior/exposure table includes transcript and gDNA
EM effective lengths.

Update output tests:

- `quant.feather` and `nrna_quant.feather` include `em_exposure_factor`.
- `quant.feather` includes diagnostic exposure-adjusted TPM columns side by side with regular TPM.
- `nrna_quant.feather` includes the corresponding diagnostic exposure-adjusted TPM column.
- `loci.feather` includes `gdna_exposure_factor`, `gdna_eff_len_em`, and
  `gdna_eff_len_adjustment_ratio`.
- Removed old columns are absent.

### 8.5 Golden Outputs

Golden output headers will change. Regenerate only after unit tests confirm:

```text
uniform exposure applied to all components leaves posterior competition unchanged
high exposure increases all affected component EM effective lengths
no path divides gDNA effective length by exposure
missing region coverage fails before native EM in strict mode
```

## 9. Validation Commands

Targeted:

```bash
conda activate rigel && pytest \
    tests/test_exposure.py \
    tests/test_calibration_prior.py \
    tests/test_per_locus_gdna_mass.py \
    tests/test_pipeline_wiring.py \
    tests/test_batch_em_impl.py \
    -v
```

Golden/output contract:

```bash
conda activate rigel && pytest tests/test_golden_output.py -v
```

Full suite:

```bash
conda activate rigel && pytest tests/ -q
```

Lint touched files:

```bash
conda activate rigel && ruff check \
    src/rigel/calibration/_exposure.py \
    src/rigel/calibration/prior.py \
    src/rigel/pipeline.py \
    src/rigel/estimator.py \
    src/rigel/cli.py \
    tests/test_exposure.py \
    tests/test_calibration_prior.py \
    tests/test_per_locus_gdna_mass.py \
    tests/test_pipeline_wiring.py \
    tests/test_batch_em_impl.py
```

Final grep sentinel:

```bash
grep -R --exclude-dir='__pycache__' -n -E \
  'gdna_eff_len.*\/.*exposure|\/.*gdna_exposure|gdna_em_exposure_weight|gdna_eff_len_weight_ratio' \
  src/rigel tests
```

Expected: no production path divides gDNA effective length by exposure and no output uses the old
asymmetric names.

## 10. Implementation Order

1. Mark PR05 v2 as superseded by this v3 plan.
2. Add `component_bp_weighted_exposure()` with strict coverage validation, finite positive exposure
  checks, and a defensive exposure-factor floor.
3. Add transcript block extraction helper from `TranscriptIndex._t_exon_intervals`.
4. Add gDNA locus block extraction helper from `MultiLocus.loci`.
5. Add `ComponentExposureTable` and `compute_component_exposure_table()`.
6. Refactor prior assembly so adaptive mass and component exposure are separate phases.
7. Install `transcript_eff_len_em` into `AbundanceEstimator` before native EM.
8. Pass `gdna_eff_len_em` to `_run_locus_em_partitioned()`.
9. Update locus metadata and transcript/nRNA output schemas, including diagnostic exposure-adjusted
  TPM columns.
10. Update CLI summary.
11. Add fairness/sign sentinel EM tests plus mixed-baiting, heterogeneous-exposure, and missing
  coverage robustness tests.
12. Update golden outputs intentionally.
13. Run targeted tests, golden test, full suite, and final grep.

## 11. Resolved Decisions And Known Limitations

### Decision 1 - TPM Reporting

Primary TPM remains exposure-neutral in PR05. Keep public `effective_length` and existing TPM
denominators unexposed. Add diagnostic exposure-adjusted TPM columns using `em_effective_length` so
benchmarking can inspect how much exposure changes the reported scale without redefining the main
abundance contract.

A later reporting PR can decide whether to add exposure-normalized expression outputs.

### Decision 2 - Component-Level Approximation Versus Local Exposure Likelihood

Use bp-weighted component exposure factors. Do not implement fragment-local exposure numerator terms
yet.

The full local likelihood is more exact but requires fragment-to-region exposure lookup for every
candidate and should be a separate PR with performance profiling. PR05 must document the smoothing
limitation near the code that computes `transcript_eff_len_em`: highly heterogeneous capture across a
long transcript is collapsed to one component-level average.

### Decision 3 - Missing Region Handling

Strict by default. The region table should cover intergenic regions and every indexed sequence used
by EM. Missing coverage is an index/calibration bug and should raise a controlled Python error before
native EM. Exposure factors are still numerically floored and validated at the packing boundary so
the native solver never ingests zero, negative, NaN, or infinite effective-length adjustments.

### Decision 4 - Where To Store Transcript Exposure Arrays

Store one exposure factor per transcript separately from effective length. Install the derived
`transcript_eff_len_em` on `AbundanceEstimator` for native EM, and report both `em_effective_length`
and `em_exposure_factor` in transcript/nRNA output. Avoid adding large transcript arrays to
`PriorTable` unless it substantially simplifies the first implementation.

## 12. Done Means

- Exposure is described and implemented as source-agnostic.
- `omega` is learned from gDNA but applied to all overlapping EM components.
- Transcript EM effective lengths are exposure-adjusted.
- Transcript exposure factors are stored separately from the derived EM effective lengths.
- gDNA EM effective lengths are exposure-adjusted in the same direction.
- Uniform exposure applied to every component in a locus does not change gDNA/RNA posterior
  competition.
- Mature transcript exposure is computed from transcript exon blocks; synthetic nRNA uses its own
  transcript block.
- Intergenic region exposure comes from the region table; there is no production neutral fallback for
  missing regions.
- Missing region coverage fails safely in Python before native EM, and exposure factors are finite
  and floored before packing C++ inputs.
- Regular TPM remains exposure-neutral, with exposure-adjusted TPM provided only as a diagnostic.
- No production path divides gDNA effective length by exposure.
