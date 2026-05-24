# Fine Region Phase 4 Calibration Implementation Plan

Date: 2026-05-22
Status: planning document
Prerequisites: Phase 0-3 region partition, native fractional accumulator, Python payload cutover
Companion plan: `rnaseq_mode_aware_gdna_density_plan.md` covers the
mode-aware density architecture across no-capture, hybrid-capture,
unstranded, and strand-specific RNA-seq libraries.

## 1. Purpose

Phase 4 restores Rigel to a fully functional `rigel quant` after the fractional
accumulator cutover. The new region partition and 12-channel fractional payload
are now the foundation. Calibration must consume that foundation and produce the
objects needed by downstream locus EM.

Calibration has four concrete goals:

1. Estimate RNA fragment length (FL) from spliced RNA fragments. This is already
   done by `calibration.fl` from scanner-trained `SPLICED_ANNOT` counts.
2. Estimate gDNA FL from gDNA-compatible unspliced FL pools. This is already
   structurally ready: `INTERGENIC_{CONTAINED,BOUNDARY}` plus
   `INTRONIC_{CONTAINED,BOUNDARY}` aggregate into `FLModels.gdna`.
3. Estimate gDNA mass and density in every region of the genome partition.
4. Aggregate per-region gDNA estimates into a global gDNA estimate, then project
   the same estimates into per-locus gDNA prior counts for EM.

The outcome of Phase 4 is not just a density table. It is a restored pipeline:
scan -> calibration -> scoring -> locus construction -> prior assembly -> EM ->
quant outputs.

## 2. Current Foundation

Already implemented:

- `CalibrationScanPayload` carries `region_counts[R,12]`, `channel_mass[12]`,
  `signature_mass[16]`, six FL pools, and explicit exclusion counters.
- `FractionalEvidenceView` exposes signature masks, strand-relative splits, FL
  pool helpers, and per-channel views.
- `RegionArrays` and `PayloadArrays` give sorted, per-reference arrays for
  region geometry and hot unspliced channel views.
- `fractional_boundary_side_exposure(lengths_bp, gdna_fl)` implements the exact
  per-side fractional exposure formula.
- `calibrate(...)` builds `FLModels` and diagnostics, then returns a result with
  no global density, uniform regional exposure, and an empty prior table.
- `quant_from_buffer(...)` raises `FractionalCutoverPending` after calibration.
  The CLI exits 64 and writes a partial `summary.json`.

Phase 4 must replace the fail-fast boundary with real calibration and prior
assembly.

## 3. Core Design Decisions

### 3.1 Region signatures are the semantic contract

Every estimator consumes the 4-bit region `signature` and the 12 fractional
channels. Do not reconstruct the old 8-mask payload. Do not reintroduce the
old orient axis.

Use these region classes for estimator policy:

| Class | Signatures | Main use |
| --- | --- | --- |
| `INTERGENIC` | `0x0` | Highest-purity gDNA anchor |
| `INTRON_POS` | `0x8` | gDNA anchor with possible pos-strand nRNA |
| `INTRON_NEG` | `0x4` | gDNA anchor with possible neg-strand nRNA |
| `INTRON_BISTRAND` | `0xC` | Lower-confidence intronic anchor |
| `EXON_POS` | `0x2`, `0xA` | RNA-contaminated; usually impute/project gDNA |
| `EXON_NEG` | `0x1`, `0x5` | RNA-contaminated; usually impute/project gDNA |
| `EXON_BISTRAND` | `0x3` | RNA-contaminated; low direct reliability |
| `MIXED_EXON_INTRON` | all remaining signatures with exon and intron bits | RNA-contaminated; use neighbor/projection plus diagnostics |

The names `coarse` and `fine` should not appear in new implementation APIs.
Those were migration terms, not runtime semantics.

### 3.2 Do not depend on `boundary_kind`

Phase 4 should derive boundary relationships from `signature`,
`left_signature`, and `right_signature`. `boundary_kind_{left,right}` can remain
on disk for now, but the estimator should not read it. After Phase 4 is green,
we can remove boundary-kind instrumentation in a separate cleanup if it remains
unused.

### 3.3 Calibration no longer uses splicing-anchor tolerance

`splicing_anchor_tolerance` remains resolver provenance only. Calibration must
not use it for boundary filtering or exposure denominators. Tiny overhangs are
naturally down-weighted by the fractional accumulator when the estimator uses
the intron/intergenic side of the boundary.

### 3.4 gDNA is unspliced and unstranded

The estimator uses unspliced channels for gDNA. Spliced channels are RNA-only
evidence and should be used to diagnose or downweight RNA contamination, not as
gDNA numerator mass.

gDNA should be symmetric across POS/NEG channels. Strand-specific RNA creates
sense/antisense imbalance in transcript-bearing regions. `sense_antisense_split`
is the only API for that decomposition.

### 3.5 Estimate mass and uncertainty, not just point densities

Every region-level estimate should carry:

- observed candidate mass;
- exposure/opportunity;
- raw density;
- shrunk density;
- expected gDNA mass;
- precision or uncertainty proxy;
- source/flags explaining how the estimate was obtained.

This is necessary because per-locus priors should not treat a tiny, noisy region
and a long, well-observed region as equally certain.

## 4. New Calibration Products

### 4.1 `RegionGdnaTable`

Add a new dataclass, likely in a new module `calibration/region_gdna.py` or in
`density_global.py` if keeping the existing module boundary is preferred.

Required arrays, all in `RegionArrays` sorted-position order:

```python
@dataclass(frozen=True, slots=True)
class RegionGdnaTable:
    ref_id: np.ndarray              # int32[R]
    start: np.ndarray               # int64[R]
    end: np.ndarray                 # int64[R]
    length: np.ndarray              # float64[R]
    signature: np.ndarray           # uint8[R]
    left_signature: np.ndarray      # uint8[R], 0xFF sentinel allowed
    right_signature: np.ndarray     # uint8[R], 0xFF sentinel allowed

    observed_unspliced: np.ndarray  # float64[R]
    observed_spliced: np.ndarray    # float64[R]
    observed_gdna_candidate: np.ndarray

    contained_exposure: np.ndarray  # float64[R]
    left_exposure: np.ndarray       # float64[R]
    right_exposure: np.ndarray      # float64[R]
    total_exposure: np.ndarray      # float64[R]

    rho_raw: np.ndarray             # float64[R]
    rho_hat: np.ndarray             # float64[R]
    gdna_mass: np.ndarray           # float64[R] = rho_hat * total_exposure
    precision: np.ndarray           # float64[R]
    source: np.ndarray              # uint8[R]
    flags: np.ndarray               # uint32[R]
```

Source enum examples:

- `INTERGENIC_DIRECT`
- `INTRON_STRAND_DECONV`
- `INTRON_TOTAL_WEAK`
- `EXON_NEIGHBOR_IMPUTED`
- `MIXED_NEIGHBOR_IMPUTED`
- `GLOBAL_FALLBACK`
- `NO_EXPOSURE`

Flags should include at least:

- `FLAG_LOW_EXPOSURE`
- `FLAG_STRAND_UNIDENTIFIABLE`
- `FLAG_RNA_CONTAMINATION_RISK`
- `FLAG_NEIGHBOR_IMPUTED`
- `FLAG_GLOBAL_FALLBACK`
- `FLAG_ZERO_OBSERVED`
- `FLAG_REF_EDGE`

### 4.2 Revised `GlobalDensityTable`

The old `GlobalDensityTable` represented a few class-level densities. Phase 4
should make it a summary of the region estimator:

```python
@dataclass(frozen=True, slots=True)
class GlobalDensityTable:
    gdna_fl: FragmentLengthModel
    region_gdna: RegionGdnaTable
    n_gdna_global: float
    total_exposure: float
    rho_global: float
    by_region_class: dict[str, RegionClassSummary]
    estimator_version: str = "fractional_v1"
```

`to_summary_dict()` must emit JSON-safe diagnostics:

- `n_gdna_global`
- `rho_global`
- `total_exposure`
- per-class observed mass, estimated mass, exposure, density, number of regions,
  and fraction imputed;
- strand-deconvolution usage counts;
- fallback counts;
- FL model quality and gDNA FL mean.

Keep compatibility properties only if downstream code still expects them. Do not
preserve misleading fields with old semantics.

### 4.3 Real `RegionalGdnaExposure`

`RegionalGdnaExposure.build(...)` should consume `RegionGdnaTable`, not raw
legacy arrays. It should produce the denominator exposure weights used by EM:

```text
rho_ref = weighted_quantile(rho_hat, weights=total_exposure, q=reference_quantile)
A_r = clip(rho_hat / rho_ref, floor, 1.0)
```

Initial defaults:

- floor: `1e-4` or existing minimum exposure weight;
- mode: `kappa_eb` or rename to `fractional`;
- if disabled in config, return `uniform`.

The regional exposure object must retain region lookup arrays (`ref_offsets`,
`ref_id`, `start`, `end`) because `footprint_exposure_weight` depends on them.
It also needs an efficient `weighted_length_on_ref(...)` implementation if that
method is not already present in the current stub.

### 4.4 Restored `PriorTable`

`PriorTable` remains the EM hand-off. Phase 4 should keep these arrays:

- `gdna_prior_count`: physical expected gDNA count per multi-locus before EM
  denominator weighting;
- `gdna_prior_count_em`: count actually passed to EM after exposure/effective
  length adjustments;
- `gdna_eff_len`: FL-marginal gDNA effective length per multi-locus;
- `enable_gdna`: explicit eligibility flag;
- `gdna_eff_len_unweighted`;
- `gdna_em_exposure_weight`.

## 5. Region-Level Estimator

### 5.1 Evidence extraction

Build vectorized arrays from `FractionalEvidenceView` and `PayloadArrays`:

```text
U_contained = contained_unspliced_pos + contained_unspliced_neg
U_left      = boundary_left_unspliced_pos + boundary_left_unspliced_neg
U_right     = boundary_right_unspliced_pos + boundary_right_unspliced_neg
U_total     = U_contained + U_left + U_right
S_total     = contained_spliced + boundary_left_spliced + boundary_right_spliced
```

`PayloadArrays` currently materializes only unspliced channels. Phase 4 either:

1. extends it with spliced totals, or
2. reads spliced totals through `FractionalEvidenceView` when diagnostics need
   them.

For performance, prefer extending `PayloadArrays` once Phase 4 starts using
spliced totals frequently.

### 5.2 Exposure extraction

For each region span `S = end - start`:

```text
E_contained = l_eff_contained(S, gdna_fl)
E_side      = fractional_boundary_side_exposure(S, gdna_fl)
E_left      = E_side if left_signature != 0xFF else 0
E_right     = E_side if right_signature != 0xFF else 0
E_total     = E_contained + E_left + E_right
```

For internal regions, `E_total` should be close to region length. This is a
useful invariant: under uniform gDNA, the fractional accumulator allocates one
fragment's mass across regions without length bias.

Unit tests must cover:

- large-region degenerate FL: `E_side = (L - 1) / 2`;
- short-region degenerate FL: `E_side = S / 2`;
- internal total exposure approaches `S`;
- reference-edge regions lose the missing-side exposure.

### 5.3 Direct high-confidence gDNA anchors

#### Intergenic regions

For `signature == 0x0`:

```text
Y_g = U_total
E_g = E_total
rho_raw = Y_g / E_g
```

This is the highest-purity gDNA channel.

#### Intron-only regions

For `signature in {0x4, 0x8, 0xC}`:

- Use `sense_antisense_split(..., splice=UNSPLICED)` for contained and boundary
  channels, summed across compartments.
- If the spliced strand model is identifiable, estimate and subtract
  strand-specific nRNA contamination.
- If strand is not identifiable, use intron-only total as a weak gDNA anchor and
  mark `FLAG_STRAND_UNIDENTIFIABLE | FLAG_RNA_CONTAMINATION_RISK`.

For identifiable strand, use the RNA/gDNA two-component moment:

```text
q = strand_summary.strand_specificity
contrast = 2*q - 1
sense = unspliced mass on transcript-sense channel
antisense = unspliced mass on transcript-antisense channel
rna_hat = max((sense - antisense) / contrast, 0)
gdna_hat = clamp(sense + antisense - rna_hat, 0, sense + antisense)
```

For bistrand intron signatures (`0xC`), strand-relative decomposition is
ambiguous. Treat them as weak anchors unless Phase 4 adds a bistrand-specific
split.

### 5.4 Boundary-side use

Boundary evidence should be used on the receiving region side. Do not sum both
sides of an exon-intron pair to estimate gDNA. The important case remains:

```text
100 bp unspliced RNA-like overhang: 98 bp exon, 2 bp intron
exon-side mass   = 0.98  # RNA-contaminated; do not use as gDNA anchor
intron-side mass = 0.02  # small noise contribution
```

Therefore:

- Intronic or intergenic side mass is a gDNA-compatible numerator.
- Exonic side mass is RNA-contaminated and should not be a primary gDNA anchor.
- Mixed signatures are reported separately and either downweighted or imputed.

### 5.5 Exon and mixed regions

Exon and mixed regions need a gDNA estimate for EM, but their direct unspliced
mass is often dominated by mRNA/nRNA. Phase 4 should not pretend those counts
are clean.

Initial implementation policy:

1. Compute a weak direct exonic candidate only when strand is identifiable:
   use the same sense/antisense moment as introns, but assign lower precision
   and mark RNA-contamination risk.
2. Compute an imputed density from neighboring high-confidence regions on the
   same reference.
3. Precision-weight the weak direct candidate and the imputed candidate.
4. If neither is available, fall back to the global density.

Neighbor imputation should use `left_signature` / `right_signature` and sorted
reference arrays, not `boundary_kind`.

Recommended first imputation algorithm:

- Build an anchor mask from intergenic and confident intron-only regions.
- For every non-anchor region, find the nearest left and right anchors on the
  same reference.
- Combine anchors with inverse-distance weights capped by a maximum window
  (for example 100 kb initially, configurable later only if needed).
- Weight each anchor by both distance and anchor precision.
- If no anchor is found within the window, use the global density.

This is intentionally simple, vectorizable, and debuggable. It can later be
replaced with a smoother or capture-aware model without changing downstream
interfaces.

### 5.6 Empirical-Bayes shrinkage

Raw region densities are sparse. Every direct estimate should be shrunk toward a
class-level or global prior:

```text
rho_hat = (Y_g + alpha * rho_prior) / (E_g + alpha)
precision = E_g + alpha
```

Initial priors:

- intergenic -> global intergenic density;
- confident intron-only -> global intron direct density;
- weak intron-only -> mixture of global intron and global intergenic;
- exon/mixed -> neighbor-imputed density, then global fallback.

Use existing `KappaEstimate` machinery where it still makes sense, but rewrite
its inputs around fractional region evidence. Do not import old integer payload
helpers.

## 6. Global gDNA Estimate

Once every region has `rho_hat` and `total_exposure`:

```text
gdna_mass_r = rho_hat_r * total_exposure_r
n_gdna_global = sum_r gdna_mass_r
rho_global = n_gdna_global / sum_r total_exposure_r
```

This is the primary global gDNA estimate. It should replace the removed
`estimate_global_gdna_fragments(...)` call.

Diagnostics must also report anchor-only totals:

- intergenic direct observed mass and exposure;
- intron direct observed mass and exposure;
- exon/mixed imputed mass;
- number and fraction of regions using global fallback;
- estimated global gDNA contamination rate when total fragment count is known.

Do not use gDNA FL pool total as the global gDNA count. FL pool mass is training
evidence for the gDNA FL shape, not an unbiased total gDNA estimator.

## 7. Per-Locus Prior Assembly

`assemble_priors(...)` must project region-level gDNA estimates into each
`MultiLocus` built from EM units.

### 7.1 Locus-local expected gDNA count

For each locus interval, find overlapping regions with `RegionIndexPy`. For each
overlapping region:

1. Prorate contained mass using `contained_exposure_clipped(...)`:

   ```text
   contained_fraction = E_contained_clip / max(E_contained_full, eps)
   ```

2. Include left/right boundary-side mass only when that side boundary is inside
   the locus window using `boundary_side_in_window(...)`.
3. Combine contained and boundary components into a locus-local expected gDNA
   mass.

The result is a `LocusGdnaEstimate` with:

- intergenic contribution;
- intron contribution;
- boundary contribution;
- exon/mixed imputed contribution;
- total expected gDNA;
- local density tuple diagnostics;
- fallback flags.

### 7.2 Multi-locus aggregation

For each `MultiLocus`:

```text
n_obs = number of EM units in the multi-locus
n_gdna = sum per-locus expected gDNA
n_rna = max(n_obs - n_gdna, 0)
pi_gdna = n_gdna / max(n_obs, 1)
gdna_prior_count = n_gdna
```

Then compute EM denominator fields:

```text
gdna_eff_len_unweighted = gdna_eff_len_for_loci(...)
gdna_em_exposure_weight = footprint_exposure_weight(..., regional_exposure)
gdna_eff_len = gdna_eff_len_unweighted * gdna_em_exposure_weight
gdna_prior_count_em = gdna_prior_count * gdna_em_exposure_weight
```

The exact scaling of `gdna_prior_count_em` should be verified against existing
EM behavior. The important invariant is that the prior count and the effective
length denominator are transformed consistently.

### 7.3 gDNA eligibility

`enable_gdna` should be independent of prior size. A locus should enable gDNA
when:

- at least one unspliced EM unit has finite gDNA likelihood; and
- `gdna_eff_len > 0`; and
- the locus is not explicitly disallowed by annotation or config.

Do not disable gDNA solely because the expected prior count is tiny. A small
prior can still regularize EM without removing the gDNA component.

## 8. Pipeline Restoration

Phase 4 must restore the body of `quant_from_buffer(...)` intentionally, not by
blindly pasting the deleted legacy block.

Required steps:

1. Import the restored `assemble_priors`, `build_multi_loci`, and
   `partition_and_free` dependencies.
2. Create geometry and estimator using `calibration.regional_exposure`.
3. Score fragments using `calibration.fl_models.rna` and
   `calibration.fl_models.gdna`.
4. Build multi-loci from `em_data`.
5. Call the new fractional `assemble_priors(...)` before `partition_and_free`.
6. Replace `calibration = calibration.with_priors(prior_table)` with a real
   frozen dataclass update.
7. Run partitioned EM with `gdna_prior_count_em`, `gdna_eff_len`, and
   `enable_gdna`.
8. Set `stats.n_gdna_em` from estimator EM result.
9. Set `stats.n_gdna_global` from `calibration.global_densities.n_gdna_global`.
10. Return `PipelineResult` normally so the CLI writes all quantification output
    files.

After this is complete, `FractionalCutoverPending` should not be raised in the
normal quant path.

## 9. Module-by-Module Missing Pieces

### 9.1 `fractional_evidence.py`

Missing:

- signature-class helpers for the estimator (`region_class_from_signature`,
  class-name table, direct-anchor masks);
- boundary-side evidence builder that joins left/right channel mass with
  `left_signature` and `right_signature`;
- optional spliced channel helpers if not added to `PayloadArrays`.

### 9.2 `_arrays.py`

Missing:

- `left_signature` and `right_signature` in `RegionArrays`;
- optionally `length` to avoid recomputing `end - start`;
- spliced channel totals if used frequently;
- a memory-conscious path that avoids unnecessary full `region_counts[order, :]`
  copies if large genomes make this expensive.

### 9.3 `density_global.py` or new `region_gdna.py`

Missing:

- real `compute_global_densities(...)`;
- region-level gDNA estimator;
- class summaries;
- global aggregation;
- replacement for `estimate_global_gdna_fragments(...)` or a new explicit
  property on `GlobalDensityTable`;
- JSON summary schema.

### 9.4 `_regional_exposure.py`

Missing:

- real `RegionalGdnaExposure.build(...)` from `RegionGdnaTable`;
- `weighted_length_on_ref(...)` if absent;
- per-class summaries based on the new region classes;
- removal of integer-payload language from docstrings.

### 9.5 `locus_prior.py`

Missing:

- real `estimate_locus_gdna(...)` using region estimates and clipped exposure;
- real `expected_gdna_count_global(...)` or replacement with region-table
  projection;
- real `enable_gdna_for_multilocus(...)`;
- real `assemble_multilocus_prior(...)`;
- real `assemble_priors(...)` returning populated arrays;
- updated diagnostics and fallback flags for region-estimator sources.

### 9.6 `_result.py`

Missing:

- `CalibrationResult.with_priors(...)` should become `dataclasses.replace(...)`
  again;
- dataframe schemas should include new region-estimator diagnostics where useful
  (`n_gdna_exon_imputed`, `n_gdna_global_fallback`, uncertainty/precision
  fields);
- `to_summary_dict()` should expose the new global density and regional exposure
  schema.

### 9.7 `_orchestrator.py`

Missing:

- build `RegionArrays`, `PayloadArrays`, and `FractionalEvidenceView` once;
- build `FLModels`;
- call `compute_global_densities(...)`;
- call `RegionalGdnaExposure.build(...)` when enabled;
- return `CalibrationResult` with real `global_densities` and
  `regional_exposure`.

### 9.8 `pipeline.py`

Missing:

- restore `quant_from_buffer(...)` after `assemble_priors` exists;
- remove the fail-fast raise;
- set `stats.n_gdna_global` from the new global estimate;
- remove outdated cutover docstrings/comments.

### 9.9 `cli.py`

Missing:

- remove the special exit-64 cutover path once normal quant works;
- keep any useful handling only for genuinely exceptional calibration failures;
- ensure `summary.json` contains final calibration and EM outputs.

### 9.10 `tests/conftest.py`

Missing:

- remove `collect_ignore` entries for calibration tests;
- remove pytest hooks that convert `FractionalCutoverPending` to skips;
- rewrite ignored tests against the fractional payload.

### 9.11 Test factories

Add `tests/_fractional_factories.py` with helpers:

- `empty_payload(n_regions, signatures=None)`;
- `payload_with_region_mass(...)`;
- `payload_with_pool(...)`;
- `region_df_from_signatures(...)`;
- `simple_region_arrays(...)`;
- `fragment_length_model_at(length)`;
- `strand_summary(...)`.

### 9.12 Scripts and docs

Missing:

- port `scripts/profiling/profiler.py` off deleted `_orient` imports;
- port `scripts/debug/diagnose_vcap_capture_exposure.py` or mark it obsolete;
- update `docs/MANUAL.md` and `docs/parameters.md` calibration schema;
- refresh `CHANGELOG.md` with Phase 4 restored quant behavior.

## 10. Tests Required

### 10.1 Unit tests

New or rewritten tests:

- `test_fractional_evidence.py`
  - signature class masks over all 16 signatures;
  - boundary-side evidence uses neighbor signatures, not boundary_kind;
  - no legacy orient symbols.

- `test_region_gdna_estimator.py`
  - intergenic direct estimate;
  - intron strand-deconvolution formula;
  - unidentifiable strand fallback;
  - bistrand intron weak anchor;
  - exon/mixed neighbor imputation;
  - global fallback;
  - precision increases with exposure;
  - region mass sums are nonnegative and finite.

- `test_density_global.py`
  - `compute_global_densities` returns populated `GlobalDensityTable`;
  - `n_gdna_global == sum(region_gdna.gdna_mass)`;
  - class summaries match region table aggregates;
  - no `FractionalCutoverPending`.

- `test_regional_exposure.py`
  - real `RegionalGdnaExposure.build`;
  - uniform fallback when disabled;
  - `weighted_length_on_ref` integrates region weights correctly;
  - reference-edge and empty-reference behavior.

- `test_locus_prior.py` / `test_assemble_priors.py`
  - region-to-locus clipping;
  - boundary side inclusion only when side lies in locus window;
  - multi-locus aggregation;
  - `enable_gdna` independent of prior magnitude;
  - arrays have correct dtype/shape and no NaNs.

- `test_calibrate_orchestrator.py`
  - calibration result has FL models, global densities, regional exposure,
    diagnostics, and empty priors before locus assembly.

- `test_calibration_result.py`
  - `with_priors` updates frozen result;
  - summary schema includes new calibration fields.

### 10.2 Integration tests

- Restore `test_pipeline_wiring.py` against the current scan tuple and payload.
- Restore `test_pipeline_smoke.py` so at least one quant scenario exits 0 and
  writes all output files.
- Restore `test_pipeline_integration_v6.py` expectations around populated prior
  arrays.
- Add an end-to-end synthetic case with known gDNA-only intergenic/intronic
  fragments where global gDNA estimate is close to truth.
- Add a micro-overhang case showing intron-side boundary mass contributes `o/L`,
  not 1.0.

### 10.3 Golden tests

Goldens will move when quant resumes. Update only after unit and integration
behavior is understood.

Protocol:

1. Run targeted integration tests.
2. Run `pytest tests/ --update-golden` only when changes are intentional.
3. Document changed golden files and why.
4. Verify `summary.json` contains new calibration schema and no removed cutover
   status.

## 11. Performance Requirements

Phase 4 must be vectorized and memory-conscious.

Rules:

- Use NumPy arrays, not per-region Python loops, for genome-wide estimation.
- Keep payload `region_counts` as float32; use float64 for density/exposure
  products.
- Avoid extra full `R x 12` copies beyond the sorted view unless profiling shows
  they are acceptable.
- Use per-reference CSR offsets for neighbor imputation and locus projection.
- Keep region-estimator arrays O(R), not O(R x classes x bins).
- Do not build pandas DataFrames in hot paths; only build them for output or
  tests.

Expected scale guardrails:

- Region estimator should be linear in number of regions.
- Locus projection should be linear in total number of overlapping regions
  visited by all loci, not genome-wide per locus.
- The estimator should be benchmarked on large real indexes before enabling
  expensive smoothing variants.

## 12. Implementation Order

### P4.0 Planning lock

- Land this document.
- Decide whether the first implementation uses a new `region_gdna.py` module or
  keeps the estimator in `density_global.py`.
- Confirm `boundary_kind` is not consumed by Phase 4.

### P4.1 Fractional factories and tests

- Add `tests/_fractional_factories.py`.
- Add active tests for region signatures, exposures, and simple payloads.
- Start removing `collect_ignore` one module at a time.

### P4.2 Region evidence and exposure arrays

- Extend `RegionArrays` with neighbor signatures and length.
- Extend `PayloadArrays` or `FractionalEvidenceView` with spliced totals.
- Add boundary-side evidence builder.

### P4.3 Region gDNA estimator

- Implement direct intergenic and intron estimates.
- Implement strand-deconvolution and unidentifiable fallback.
- Implement EB shrinkage.
- Implement exon/mixed imputation.
- Produce `RegionGdnaTable`.

### P4.4 Global density table

- Implement `compute_global_densities(...)`.
- Aggregate `n_gdna_global` from the region table.
- Implement summary output.

### P4.5 Regional exposure

- Implement `RegionalGdnaExposure.build(...)` from `RegionGdnaTable`.
- Implement/verify `weighted_length_on_ref(...)`.
- Add tests for footprint weighting.

### P4.6 Locus prior assembly

- Implement `estimate_locus_gdna(...)`.
- Implement `assemble_multilocus_prior(...)` and `assemble_priors(...)`.
- Verify all `PriorTable` arrays are populated and consistent.

### P4.7 Restore quant pipeline

- Restore `quant_from_buffer(...)` scoring, locus construction, prior assembly,
  partitioning, and EM.
- Re-enable normal CLI output.
- Remove the exit-64 cutover path from normal operation.

### P4.8 Remove cutover gates

- Remove `FractionalCutoverPending` from normal calibration modules.
- Remove pytest skip hooks and `collect_ignore` entries.
- Port or retire legacy scripts.

### P4.9 Validation and benchmarking

- Run `ruff check src/ tests/`.
- Run the full pytest suite.
- Refresh goldens intentionally.
- Run synthetic locus sweep.
- Run the full-scale benchmark subset that previously exposed gDNA -> RNA
  leakage.

## 13. Acceptance Criteria

Phase 4 is complete when all of the following are true:

- `rigel quant` exits 0 on standard synthetic and scenario tests.
- `quant.feather`, `gene_quant.feather`, `nrna_quant.feather`, `loci.feather`,
  and `summary.json` are written normally.
- `CalibrationResult.global_densities` is populated.
- `CalibrationResult.regional_exposure` is populated or intentionally uniform
  when disabled.
- `CalibrationResult.prior_table` is populated after `quant_from_buffer`.
- `stats.n_gdna_global` comes from the per-region gDNA aggregate.
- `summary.json` reports region-density diagnostics and no
  `fractional_cutover_pending` status.
- The normal code path does not raise `FractionalCutoverPending`.
- `tests/conftest.py` no longer ignores the calibration/prior test modules.
- Full pytest and ruff are green.
- Benchmarks do not show a material regression in gDNA -> RNA leakage or scan
  performance.

## 14. Open Decisions

These should be resolved before or during P4.3:

1. Exact exon/mixed direct-evidence blending: how much weight should the
   strand-deconvolved exonic candidate get relative to neighbor imputation?
2. Anchor imputation window: fixed bp radius, nearest-neighbor only, or adaptive
   by anchor exposure?
3. EB shrinkage hyperparameters: reuse `KappaEstimate`, introduce explicit
   region-class priors, or start with a simple exposure pseudocount?
4. Regional exposure floor and reference quantile defaults.
5. Whether `boundary_kind_{left,right}` should be deleted immediately after
   Phase 4 or retained until benchmark review.
6. Whether to expose `RegionGdnaTable` as an optional output file for debugging.

Recommendation for the first implementation: choose the simplest vectorized
version that restores a correct working pipeline, surfaces uncertainty and
fallback flags, and keeps richer smoothing models as follow-up improvements.
