# Strand-Based gDNA / RNA Deconvolution — Implementation Plan v5

Date: 2026-05-23
Status: implementation-ready
Supersedes: `strand_model_impl_plan_v4.md`
Design basis: `rnaseq_mode_aware_gdna_density_plan.md`, `gdna_exposure_plan_v4.3.md`
Primary use case: strand-specific hybrid-capture RNA-seq

## 0. What Changed From v4

v4's critique:

1. **Two parallel exposure arrays.** v4 had both `CaptureProfile.region_weight`
   and `RegionGdnaEstimate.exposure_weight` as `float32[R]` but never said how
   they relate or where each is consumed. v5 keeps exactly one per-region
   exposure array, `A_r`, owned by a single `RegionExposure` object that is
   written once and read many times.
2. **The "inferred" exposure path was a wishlist.** The user has committed to
   unsupervised capture exposure as the path forward. v5 specifies the
   algorithm concretely (Phase 7) and includes a measurable acceptance test.
3. **BED capture support is deferred.** v5 does not ship `--capture-targets`.
   The schema leaves room for it later, but Phase 7 owns the exposure story.
4. **The ordered implementation in v4 was a checklist, not phases.** v5
   organises the work as eight phases, each with explicit inputs, outputs,
   tests, and acceptance criteria. The pipeline becomes runnable end-to-end
   in Phase 6, with exposure improving quality in Phase 7.
5. **Prior semantics were ambiguous.** v4 said "weight by
   `precision * exposure_weight`" without specifying whether `mean_count` is
   exposure-aware. v5 fixes a clean separation: deconvolution is
   exposure-blind (it estimates `D` from strand counts only); exposure enters
   the prior on the denominator side as a gDNA effective-length scaling.
6. **A floor constant survived implicitly.** v5 reuses
   `STRAND_CONTRAST_NUMERICAL_FLOOR` already defined in `strand_summary.py`
   rather than introducing a new one. The kappa-not-identified case is named
   `DNA_STRAND_UNIFORM_PRIOR` and reuses the existing `STRAND_KAPPA_MIN /
   STRAND_KAPPA_MAX` bounds from `density_global.py`.

## 1. System Shape

```text
BAM
  -> C++ scanner
  -> FragmentBuffer + StrandModel + fractional region accumulator
  -> calibration.calibrate(...)
       1. FL models
       2. strand kappa_d seed from high-purity regions
       3. exon self-training screen using the rna_lower_confidence rule
       4. strand kappa_d refit
       5. exposure-blind RegionGdnaEstimate (per-region D, R lower bound)
       6. unsupervised RegionExposure (per-region A_r)
  -> prior assembly (numerator from RegionGdnaEstimate, denominator scaled by A_r)
  -> EM
```

Two artifacts cross the calibration boundary:

- `FLModels` (existing, unchanged).
- `RegionGdnaEstimate` (new) and `RegionExposure` (new), bundled into
  `CalibrationResult`.

## 2. The Single Statistical User Knob

```python
@dataclass(frozen=True)
class CalibrationConfig:
    rna_lower_confidence: float = 0.95
```

CLI:

```bash
rigel quant ... --rna-lower-confidence 0.95
```

Semantics:

- Controls the RNA lower bound `P(R >= rna_lower_count | data) >=
  rna_lower_confidence`.
- Controls the exon self-training screen via the dual rule
  `P(R > 0 | data, kappa_d_seed) <= 1 - rna_lower_confidence`.
- Reported in `summary.json` and on the `RegionGdnaEstimate`.

No `STRAND_ALPHA`. No exon n-min. No spliced threshold parameter. The user
has exactly one knob to set, and it controls all per-region calls and the
self-training screen together.

## 3. Schema Contracts

### 3.1 Per-region strand-deconvolution output

```python
@dataclass(frozen=True, slots=True)
class StrandRegionCounts:
    k_sense: np.ndarray       # float32[R]
    k_antisense: np.ndarray   # float32[R]
    n_total: np.ndarray       # float32[R] = k_sense + k_antisense
    eligible: np.ndarray      # bool[R]
    p_r1_sense: float

@dataclass(frozen=True, slots=True)
class RegionGdnaEstimate:
    n_total: np.ndarray          # float32[R]
    mean_count: np.ndarray       # float32[R]   E[D | data]
    upper_count: np.ndarray      # float32[R]   conservative upper bound on D
    rna_lower_count: np.ndarray  # float32[R]   confidence bound on R
    precision: np.ndarray        # float32[R]   in [0, 1]
    flags: np.ndarray            # uint8[R]
    kappa_d: float
    kappa_d_n_seed_regions: int
    kappa_d_n_exon_self_training: int
    p_r1_sense: float
    rna_lower_confidence: float
```

`RegionGdnaEstimate` carries no exposure field. The deconvolution is
exposure-blind by design: it answers "given the observed sense/antisense
split, how much of `n_total` is gDNA?" and nothing else.

### 3.2 Per-region exposure output

```python
@dataclass(frozen=True, slots=True)
class RegionExposure:
    mode: Literal["uniform", "unsupervised"]
    A_r: np.ndarray              # float32[R] in (0, 1], 1.0 means fully exposed
    rho_r: np.ndarray            # float32[R] per-region gDNA density (frag / bp)
    rho_ref: float               # reference density used for normalization
    reference_quantile: float    # weighted quantile setting used
    eligible: np.ndarray         # bool[R]    rows that contributed to rho_ref
    flags: np.ndarray            # uint8[R]   FLAG_AT_FLOOR, FLAG_NO_DENSITY, ...
```

`RegionExposure.uniform(R)` returns `A_r = 1.0` for all rows. This is the
only shape used until Phase 7.

### 3.3 Per-locus prior output

```python
@dataclass(frozen=True, slots=True)
class PriorTable:
    gdna_prior_count: np.ndarray   # float64[L]  sum of mean_count over reliable rows
    rna_prior_count: np.ndarray    # float64[L]  sum of (n_total - mean_count)
    gdna_eff_len: np.ndarray       # float64[L]  exposure-scaled denominator
    enable_gdna: np.ndarray        # uint8[L]
```

`gdna_eff_len = gdna_eff_len_raw * A_locus`, where `A_locus` is a
bp-weighted average of `A_r` over the locus footprint (the rule from
`gdna_exposure_plan_v4.3.md §5`). Numerator counts are not multiplied by
`A_r`; only the denominator is.

## 4. Phases

Each phase is one PR. Do not bundle.

### Phase 1 — Cleanup

**Inputs:** current tree.

**Outputs:** trimmed `src/rigel/calibration/`, no `FractionalCutoverPending`,
no stale CLI flags.

**Source deletions:**

```text
src/rigel/calibration/_regional_exposure.py
src/rigel/calibration/_locus_n_obs.py
src/rigel/calibration/_region_index_py.py
src/rigel/calibration/density_loco.py
src/rigel/calibration/_kappa.py
src/rigel/calibration/locus_prior.py
src/rigel/calibration/errors.py
```

**Module renames / replacements (later phases):** `locus_prior.py` →
`prior.py` (Phase 5), new files `strand_deconv.py` (Phase 3) and
`exposure.py` (Phase 7).

**Trim, do not delete:** `_exposure.py` keeps `gdna_eff_len_for_loci`,
`contained_exposure_clipped`, `fractional_boundary_side_exposure`,
`_merged_blocks`, and a new private helper used by §5.3 below
(`bp_weighted_mean_exposure_over_blocks(...)`). Delete
`transcript_exposure_weights` and `footprint_exposure_weight` plus their
dependency on `RegionalGdnaExposure`.

**Trim `density_global.py`:** remove `KappaEstimate` (deleted), the `kappa`
field on `GlobalGdnaDensity`, and `ExonCompositeDensity`. Keep
`StrandBalanceEstimate`, `estimate_strand_balance`, `compute_global_densities`,
`l_eff_contained`.

**Trim `_result.py`:** drop `regional_exposure`, `regional_weighting_stats`,
`prior_table`, `with_priors`, `with_regional_weighting_stats`,
`multi_locus_prior_df`, `per_locus_gdna_df`. Replace with placeholder
`region_gdna: RegionGdnaEstimate | None = None` and
`region_exposure: RegionExposure | None = None` until later phases populate
them.

**Strip orchestrator:** `calibrate(...)` loses
`regional_exposure_enabled`, `regional_exposure_reference_quantile`,
`resolver_splicing_anchor_tolerance`. Adds `rna_lower_confidence` (used in
Phase 3, accepted now to avoid a config churn).

**Pipeline and CLI:**

- Delete `_apply_unit_gdna_weights` in `pipeline.py`.
- Keep the `FractionalCutoverPending` raise in `quant_from_buffer(...)` for
  now (Phase 6 removes it). Replace the raise text with a clearer
  `NotImplementedError("rigel quant: locus EM lands in Phase 6")`.
- Delete `errors.py`, the CLI `except FractionalCutoverPending` clause, the
  `_write_cutover_summary` helper, and the
  `--regional-exposure` / `--regional-exposure-reference-quantile` flags.
- Delete `CalibrationConfig.regional_exposure_enabled` and
  `regional_exposure_reference_quantile`.

**Tests:**

- Delete `tests/test_regional_exposure.py`,
  `tests/test_assemble_priors.py`,
  `tests/test_locus_prior_fused.py`,
  `tests/test_fractional_cutover_cleanup.py`,
  `tests/test_gdna_eff_len.py`,
  `tests/test_pipeline_integration_v6.py`,
  `tests/test_calibrate_orchestrator.py`.
- Other tests stay green.

**Acceptance:**

- `pytest tests/` green (with deleted files removed).
- `grep -R FractionalCutoverPending src/` returns nothing.
- `grep -R RegionalGdnaExposure src/` returns nothing.
- `wc -l src/rigel/calibration/*.py` under ~1500 lines.

### Phase 2 — Config knob

**Inputs:** trimmed tree from Phase 1.

**Outputs:** one new CLI flag, one new config field, both validated.

- `CalibrationConfig.rna_lower_confidence: float = 0.95`.
- `__post_init__` validates `0.5 <= rna_lower_confidence < 1.0`.
- CLI `--rna-lower-confidence FLOAT`.
- The orchestrator threads the value into `calibrate(...)` and writes it
  unchanged into `summary.json` under `"calibration_config"`.

**Tests:**

- `test_rna_lower_confidence_default_is_095`.
- `test_rna_lower_confidence_validation_rejects_out_of_range`.
- `test_cli_rna_lower_confidence_round_trips_to_config`.

**Acceptance:** flag and config field exist, are validated, are echoed in
`summary.json`. No behavior change yet (still raises in Phase 6 territory).

### Phase 3 — Strand deconvolution module

**Inputs:** Phase 2 done.

**Outputs:** new file `src/rigel/calibration/strand_deconv.py` and tests.

**Public surface:**

```python
def build_strand_region_counts(region_arrays, payload_arrays, *, p_r1_sense)
    -> StrandRegionCounts: ...

def estimate_kappa_d(region_arrays, payload_arrays, counts, strand_summary,
                    *, rna_lower_confidence)
    -> StrandBalanceEstimate: ...   # includes seed + exon refit

def screen_no_rna_exons(counts, region_arrays, payload_arrays,
                        kappa_d_seed, *, rna_lower_confidence)
    -> np.ndarray: ...

def deconvolve_regions_by_strand(counts, *, kappa_d, rna_lower_confidence)
    -> RegionGdnaEstimate: ...
```

**Generative model:**

```text
D = n - R
DNA_sense | D ~ BetaBinomial(D, mean=0.5, concentration=kappa_d)
RNA_sense | R ~ Binomial(R, p_r1_sense)
k_sense   = DNA_sense + RNA_sense
```

**Numerical paths** (private, not a user knob):

- exact discrete posterior over `R in {0, ..., round(n)}` when
  `n_total <= MAX_EXACT_POSTERIOR_N` (constant lives in
  `strand_deconv.py`; this is numerical, not policy).
- normal approximation otherwise, with model variance
  `Var(DNA_sense | D) = 0.25 D (D + kappa_d) / (1 + kappa_d)` and
  `Var(RNA_sense | R) = R p_r1_sense (1 - p_r1_sense)`. Use the standard
  Gaussian quantile derived from `rna_lower_confidence`.

**Identifiability guard:** reuse `STRAND_CONTRAST_NUMERICAL_FLOOR` from
`strand_summary.py`. When `|p_r1_sense - 0.5| <` floor, set
`rna_lower_count = 0`, `mean_count = n_total`, `upper_count = n_total`,
`precision = 0`, raise `FLAG_NEAR_UNSTRANDED`.

**Exon self-training:** the rule from §4 below.

**Kappa fallback:** when MoM fails (existing `estimate_strand_balance` already
flags this), use `pi_d ~ Beta(1, 1)` i.e. `kappa_d = DNA_STRAND_UNIFORM_PRIOR
= 2.0`. Defined in this module, named, justified by the comment.

**Flags (uint8):**

```python
FLAG_INELIGIBLE      = 1 << 0
FLAG_NEAR_UNSTRANDED = 1 << 1
FLAG_KAPPA_FALLBACK  = 1 << 2
FLAG_APPROX_NORMAL   = 1 << 3
FLAG_EXON_SELF_TRAIN = 1 << 4
```

**Tests** (`tests/test_strand_deconv.py`):

- `test_observations_ts_pos_and_ts_neg_basic`.
- `test_protocol_symmetry` (R1-sense vs R1-antisense produce reflected outputs).
- `test_near_unstranded_is_conservative`.
- `test_rna_lower_confidence_monotone_in_bound` (higher confidence => smaller
  `rna_lower_count`).
- `test_exact_and_normal_agree_at_boundary` (n at exact-vs-normal boundary,
  agreement within 1%).
- `test_exon_screen_accepts_balanced_no_spliced_exon`.
- `test_exon_screen_rejects_clear_rna_exon`.
- `test_exon_screen_rejects_spliced_exon`.
- `test_exon_screen_requires_both_strand_channels`.
- `test_kappa_refit_uses_accepted_exons` (refit `n_regions` > seed
  `n_regions` and final `kappa_d` is the pooled MoM estimate, not the seed
  value).
- `test_kappa_fallback_when_unidentifiable`.

**Acceptance:** `pytest tests/test_strand_deconv.py -v` green. Module is
self-contained — no orchestration wiring yet.

### Phase 4 — Wire into `CalibrationResult`

**Inputs:** Phase 3 done.

**Outputs:** `calibrate(...)` populates `region_gdna` and uniform
`region_exposure` on `CalibrationResult`; `summary.json` gains a
`strand_deconv` block. Pipeline still raises in `quant_from_buffer`.

**Orchestrator additions** (in `_orchestrator.py`):

```python
counts  = build_strand_region_counts(...)
balance = estimate_kappa_d(..., rna_lower_confidence=...)
rge     = deconvolve_regions_by_strand(counts, kappa_d=balance.kappa,
                                       rna_lower_confidence=...)
exp     = RegionExposure.uniform(R=int(rge.n_total.size))
return build_calibration_result(..., region_gdna=rge, region_exposure=exp,
                                strand_region_counts=counts)
```

`CalibrationResult` fields become required (no `None`).

**Summary block:**

```json
"strand_deconv": {
  "rna_lower_confidence": 0.95,
  "p_r1_sense": 0.99,
  "kappa_d": 123.4,
  "kappa_d_n_seed_regions": 11873,
  "kappa_d_n_exon_self_training": 8412,
  "kappa_d_fallback_used": false,
  "n_regions_eligible": 612944,
  "n_regions_approx_normal": 671,
  "n_regions_near_unstranded": 12
}
```

**Tests:** `tests/test_calibrate.py` (replaces deleted orchestrator tests).
Cover: result carries non-None `region_gdna`, non-None
`region_exposure.uniform`, forbidden legacy fields not present, summary
block present.

**Acceptance:** `calibrate(...)` runs on a synthetic calibration payload
and returns a populated result. `pytest tests/` green.

### Phase 5 — Minimal prior assembly

**Inputs:** Phase 4 done.

**Outputs:** new file `src/rigel/calibration/prior.py` with
`assemble_priors(...)` returning a `PriorTable`.

**Algorithm (per locus):**

1. Find region rows overlapping the locus footprint.
2. Filter out rows flagged `FLAG_INELIGIBLE` or `FLAG_NEAR_UNSTRANDED`.
3. Sum `mean_count` over surviving rows -> `gdna_prior_count`.
4. Sum `(n_total - mean_count)` over surviving rows -> `rna_prior_count`.
5. Compute `gdna_eff_len_raw` from `gdna_eff_len_for_loci` (existing).
6. Compute `A_locus = bp_weighted_mean_exposure_over_blocks(blocks,
   region_exposure)`. In Phase 5 with uniform exposure this is exactly
   `1.0`.
7. `gdna_eff_len = max(gdna_eff_len_raw * A_locus, 1.0)`.
8. `enable_gdna = gdna_prior_count > 0`.

`precision` is **not** used to scale counts. It is reported as a per-locus
diagnostic (`gdna_prior_precision_mean`) but the prior is a sum of
counts. This preserves the physical interpretation that
`gdna_prior_count` is "expected gDNA fragments observed in the locus
footprint during calibration."

**Tests** (`tests/test_prior.py`):

- `test_locus_with_strand_pure_intron_recovers_gdna_count`.
- `test_locus_with_strong_rna_exon_gets_small_gdna_prior`.
- `test_locus_excludes_unstranded_rows`.
- `test_uniform_exposure_leaves_eff_len_unchanged`.
- `test_overlapping_loci_share_region_contributions`.

**Acceptance:** prior file under 200 lines. Tests green.

### Phase 6 — Unblock `rigel quant`

**Inputs:** Phase 5 done.

**Outputs:** `pipeline.quant_from_buffer(...)` runs end-to-end with the new
`PriorTable`; CLI no longer special-cases anything.

**Pipeline edits:**

1. Build loci as today (no change).
2. `prior_table = assemble_priors(loci=loci, region_arrays=...,
   region_gdna=cal.region_gdna, region_exposure=cal.region_exposure)`.
3. Pass `prior_table.gdna_prior_count`, `prior_table.rna_prior_count`,
   `prior_table.gdna_eff_len`, `prior_table.enable_gdna` to the EM exactly
   where the old prior path consumed them.
4. Remove the `NotImplementedError` from Phase 1.

**Tests:**

- `tests/test_pipeline_end_to_end_strand.py` — tiny synthetic BAM, run
  `run_pipeline(...)`, assert `AbundanceEstimator` is produced and EM
  converges.

**Acceptance:**

- `rigel quant` runs end-to-end on a strand-specific fixture with uniform
  exposure.
- `grep -R FractionalCutoverPending src/ tests/` is empty.
- All existing tests still green.

### Phase 7 — Unsupervised capture exposure

**Inputs:** Phase 6 done. This is the phase where the user-committed
unsupervised exposure model lands.

**Outputs:** new file `src/rigel/calibration/exposure.py` implementing
`compute_region_exposure(...) -> RegionExposure`.

**Algorithm (concrete and committed):**

For every region with `eligible == True`:

1. `D_r = region_gdna.mean_count[r]` (exposure-blind expected gDNA count).
2. `L_r = max(end[r] - start[r], 1.0)` (region length in bp).
3. `rho_r = D_r / L_r` (per-region gDNA density, frag/bp).

Cross-region reference:

```text
rho_ref = weighted_quantile(rho_r[eligible], weights=L_r[eligible],
                            q=0.95)
```

The reference quantile is **not** a user knob. It is a numerical
normalisation choice and lives in `exposure.py` as
`EXPOSURE_REFERENCE_QUANTILE = 0.95`. (If validation later shows it
matters, we promote it; today it does not.)

Per-region exposure:

```text
A_r = clip(rho_r / rho_ref, A_FLOOR, 1.0)
```

`A_FLOOR = 1e-3` is a numerical guard against zero gDNA effective length
in the prior. It is private to `exposure.py`.

**No empirical-Bayes shrinkage.** The strand deconvolution has already
absorbed sampling variance into `mean_count`. The exposure is a
deterministic transform of `mean_count` and `L_r`. If the user wants more
shrinkage later, that is a research item, not a config knob.

**Where it plugs in:**

`_orchestrator.calibrate(...)` replaces the uniform call from Phase 4:

```python
exp = compute_region_exposure(region_arrays, region_gdna=rge)
```

`prior.assemble_priors(...)` from Phase 5 already consumes
`region_exposure` correctly via `A_locus`. No prior code changes.

**Summary additions:**

```json
"region_exposure": {
  "mode": "unsupervised",
  "rho_ref": 0.0312,
  "reference_quantile": 0.95,
  "n_regions_eligible": 612944,
  "n_regions_at_floor": 17321,
  "rho_q05": 1.1e-5,
  "rho_q50": 5.4e-4,
  "rho_q95": 0.0298,
  "A_q05": 3.5e-4,
  "A_q50": 0.017,
  "A_q95": 0.95
}
```

**Tests** (`tests/test_exposure.py`):

- `test_uniform_when_all_densities_equal` (all `A_r ≈ 1`).
- `test_zero_density_regions_get_floor`.
- `test_reference_quantile_anchors_top_5_percent_to_one`.
- `test_dense_regions_have_higher_A_than_sparse`.
- `test_exposure_propagates_into_gdna_eff_len`.

**Acceptance:**

- New unit tests green.
- On the lab's strand-specific capture fixture, the bp-weighted mean of
  `A_r` over annotated capture targets is at least 4x the bp-weighted mean
  over non-target regions. This is the measurable commitment to
  unsupervised exposure quality. If this fails, the phase is not done —
  fix the algorithm, not the assertion.
- `rigel quant` results on the same fixture show a measurable improvement
  in gDNA -> RNA misassignment on captured loci relative to Phase 6
  uniform exposure. Quantify in PR description; no regression on
  synthetic-24.

### Phase 8 — Validation

**Inputs:** Phases 1-7 done.

**Outputs:** a validation report in
`docs/benchmarks/strand_deconv_v5_validation_<date>.md`.

**Validation cases:**

1. R1-sense (`p_r1_sense = 0.99`) pure-RNA region at varied `n_total`:
   `rna_lower_count` covers true `R` with empirical frequency >=
   `rna_lower_confidence`.
2. R1-antisense (`p_r1_sense = 0.01`) reflected counts: identical coverage
   to (1). Protocol symmetry.
3. Pure-gDNA region at varied `n_total`: `rna_lower_count == 0` with
   frequency >= `rna_lower_confidence`.
4. Mixed region with known `R`, varied `kappa_d` and `p_r1_sense`:
   `mean_count` close to `D_true`.
5. Sparse intergenic data: kappa fallback widens RNA lower bounds;
   coverage holds.
6. Exposure recovery: on a synthetic hybrid-capture dataset, `A_r` over
   targets is >=4x `A_r` off-target (matches Phase 7 acceptance, but
   reported with full quantile distributions).
7. End-to-end on a real strand-specific capture truth dataset: per-locus
   gDNA fraction error vs uniform-exposure baseline, vs Phase 7
   unsupervised exposure.

If anything fails: fix the model, not the test.

## 5. Why The Implementation Is Safe In This Order

- **Phase 1** is pure deletion. If it breaks anything that wasn't already
  raising `FractionalCutoverPending`, that's a real bug uncovered, not a
  v5 regression.
- **Phase 2** adds one inert knob.
- **Phase 3** ships the math in isolation. Wrong math is caught by unit
  tests before any pipeline plumbing depends on it.
- **Phase 4** wires the result into `CalibrationResult`. The pipeline still
  raises, so no downstream consumer sees the new fields unless they're a
  test.
- **Phase 5** writes the prior against the new schema, still without EM.
- **Phase 6** is the first phase where production `rigel quant` runs. It
  uses uniform exposure, so any deviation from a known-good baseline is
  attributable to the strand model, not exposure.
- **Phase 7** is the first phase where exposure shapes EM. By isolating
  the change to one new file, regressions are easy to bisect: if the
  measurable acceptance test in §4 fails, the bug is in `exposure.py`
  alone.
- **Phase 8** is validation, not implementation. It reports, it doesn't
  decide.

## 6. Out Of Scope For v5

- `--capture-targets BED` (deferred; the unsupervised path in Phase 7 is
  the committed approach).
- nRNA-vs-mRNA per-locus decomposition.
- Mappability correction.
- FFT-based posterior batching.
- Native (C++) strand-deconv kernel.
- 99% and 99.9% bounds.
- Iterative exon self-training (more than one screen pass).
- Per-source kappa diagnostics dataclasses.
- Predictive feature-based exposure (the simple density-ratio path is
  enough until validated against a real dataset).

## 7. Definition Of Done

- No `FractionalCutoverPending` references in production code.
- `rigel quant` runs on a strand-specific fixture end-to-end (Phase 6).
- `RegionExposure` is populated by the unsupervised algorithm (Phase 7).
- The measurable exposure acceptance (>=4x target/off-target ratio) holds
  on the validation fixture.
- The user has exactly one statistical knob, `--rna-lower-confidence`, and
  it controls both RNA lower bounds and the exon self-training screen.
- `summary.json` reports the strand and exposure blocks shown above.
- All eight phases have green tests.
