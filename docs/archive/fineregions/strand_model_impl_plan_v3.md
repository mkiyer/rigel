# Strand-Based gDNA / RNA Deconvolution — Implementation Plan v3

Date: 2026-05-23
Status: ready for implementation (replaces v2)
Supersedes: `strand_model_impl_plan_v2.md`, `strand_model_impl_plan_v1.md`
Design basis: `rnaseq_mode_aware_gdna_density_plan.md`
Primary use case: strand-specific hybrid-capture RNA-seq

## 0. What Changed From v2

Three deltas, all driven by reading the existing code:

1. **Cleanup is now a prerequisite, not a roadmap item.** The calibration
   package today carries six modules whose entire job is to raise
   `FractionalCutoverPending` or hold uniform-mode placeholders. Those go
   first, with their tests, before any new strand code is written.
2. **Exon self-training is in v3.** A trained `kappa_d` from intergenic +
   intron-only regions is used to identify exon regions whose split is
   plausibly gDNA-only, and `kappa_d` is then refit on the union. This is
   one extra private function in the same module and one extra test — no
   new dataclasses, no new files.
3. **The pipeline-side `FractionalCutoverPending` raise is also a
   prerequisite for the strand model to be useful.** v3 schedules its
   removal explicitly, so `rigel quant` returns an EM result on a strand-
   specific BAM at the end of the milestone chain.

## 1. End-To-End Mental Model

This is the entire system v3 commits to. Nothing else lives in calibration.

```text
BAM
  -> C++ scanner (rigel.native)               # one pass, single source of truth
  -> FragmentBuffer + StrandModel + FL counts
  -> fractional region accumulator           # float32[R, 12] per region
  -> calibration.calibrate(...)
       1. FL models (RNA / gDNA / global)
       2. strand kappa_d seed (intergenic + intron-only)
       3. exon self-training screen
       4. strand kappa_d refit
       5. per-region gDNA/RNA deconvolution -> RegionGdnaEstimate
  -> locus prior assembly                    # consumes RegionGdnaEstimate only
  -> EM
```

Two and only two artifacts cross the calibration boundary:

- `FLModels` — already exists, unchanged.
- `RegionGdnaEstimate` — new, defined below.

Anything else currently flowing through (regional exposure, locus n_obs
helpers, `KappaEstimate`, locoregional shrinkage, `RegionalGdnaExposure`,
`transcript_exposure_weights`, the legacy unit-gDNA-weights helper, the
PriorTable / MultiLocusPrior / LocusGdnaEstimate fail-fast stubs) is
deleted in M1 before the new code lands.

## 2. Naming Rules (unchanged from v2)

Transcript-relative throughout, no `major` / `minor` / `primary` /
`secondary` / rotated-frame names anywhere in code, tests, or summaries.

| Name | Meaning |
|---|---|
| `k_sense` | observed unspliced count on the transcript-sense strand |
| `k_antisense` | observed unspliced count on the transcript-antisense strand |
| `n` | `k_sense + k_antisense` |
| `p_r1_sense` | `P(R1 aligns transcript-sense)`, from `StrandModel` |
| `R` | unknown RNA-derived count for a region |
| `D` | unknown gDNA-derived count for a region, `D = n − R` |
| `kappa_d` | beta-binomial concentration for gDNA strand imbalance |

Protocol direction enters the model only through `p_r1_sense`. The signed
moment estimator `R_hat = (k_sense − 0.5·n) / (p_r1_sense − 0.5)` flips
itself between R1-sense and R1-antisense preps; no rotation logic.

## 3. M1 — Calibration Cleanup (Delete First)

Goal: make the package small enough that the strand work lands into a
tree a human can hold in their head. Each deletion is independently
test-driven.

### 3.1 Delete these source modules

| File | Reason |
|---|---|
| `src/rigel/calibration/_regional_exposure.py` | All non-`uniform()` methods raise `FractionalCutoverPending`; field is unread by EM. |
| `src/rigel/calibration/_exposure.py` | `transcript_exposure_weights`, `boundary_crossing_exposure`, `fractional_boundary_side_exposure`, `gdna_eff_len_for_loci` — feed only the regional-exposure / locus-prior fail-fast paths. |
| `src/rigel/calibration/_locus_n_obs.py` | Used by the deleted legacy prior path. |
| `src/rigel/calibration/_region_index_py.py` | Used only by the deleted legacy prior path. |
| `src/rigel/calibration/density_loco.py` | `shrink_to_loco` is exported but not imported anywhere outside `__init__.py`. |
| `src/rigel/calibration/_kappa.py` | NB MoM `KappaEstimate` only feeds `density_global`'s per-density-type kappa diagnostic, which is a summary field nothing consumes. Strand kappa (`estimate_strand_balance`) is unrelated and stays. |
| `src/rigel/calibration/locus_prior.py` | Entire file is fail-fast stubs + dataclass shells. Replaced in M3 by a real, minimal `prior.py`. |
| `src/rigel/calibration/_diagnostics.py` | Verify after the above: if `Diagnostics` is only the now-removed `n_observed` breakdown, drop. Otherwise reduce to the surviving fields. |

### 3.2 Trim `density_global.py`

After `_kappa.py` is gone:

- Remove `KappaEstimate` import and the `kappa` field on `GlobalGdnaDensity`.
- Remove `ExonCompositeDensity` (placeholder; never populated downstream).
- Keep `compute_global_densities(...)`, `estimate_strand_balance(...)`,
  `l_eff_contained(...)`, `GlobalDensityTable`, `GlobalGdnaDensity`,
  `StrandBalanceEstimate`. These are the only pieces other code touches.
- Reduce `to_summary_dict()` to fields that still have meaning.

After this, `density_global.py` should be roughly half its current size.

### 3.3 Trim `_result.py`

- Drop `regional_exposure`, `regional_weighting_stats`,
  `multi_locus_prior_df`, `per_locus_gdna_df`, `prior_table`,
  `with_priors`, `with_regional_weighting_stats`,
  `resolver_splicing_anchor_tolerance`.
- Keep `fl_models`, `global_densities`, `diagnostics` (or whatever
  survives from §3.1), and add the two new fields from §4.
- Update `build_calibration_result(...)` to match.

### 3.4 Strip the orchestrator

`src/rigel/calibration/_orchestrator.py` shrinks to:

```python
def calibrate(
    *,
    index: TranscriptIndex,
    payload: CalibrationScanPayload,
    scan_trained: FragmentLengthModels,
    strand_summary: StrandSummary,
    fl_prior_ess: float = POOL_EB_PRIOR_ESS,
    pool_quality_good: int = POOL_QUALITY_GOOD_THRESHOLD,
    pool_quality_weak: int = POOL_QUALITY_WEAK_THRESHOLD,
) -> CalibrationResult:
    ...
```

`regional_exposure_enabled`, `regional_exposure_reference_quantile`,
`resolver_splicing_anchor_tolerance` parameters are removed. `strand_summary`
becomes required (callers always have one in `pipeline.run_pipeline`).

### 3.5 Strip the public surface

`src/rigel/calibration/__init__.py` is reduced to the names that survive,
in alphabetical order:

```text
CalibrationResult, FLModels, GlobalDensityTable, GlobalGdnaDensity,
RegionArrays, PayloadArrays, RegionGdnaEstimate, StrandBalanceEstimate,
StrandSummary, build_fl_models, calibrate, compute_global_densities,
estimate_strand_balance, l_eff_contained
```

Everything else either becomes private or stops being a calibration
public API.

### 3.6 Clean the pipeline and CLI

`src/rigel/pipeline.py`:

- Delete `_apply_unit_gdna_weights(...)` (≈90 lines, regional-exposure
  legacy).
- Delete the `regional_exposure` parameter and branch in `quant_from_buffer`
  callers; drop the import of `transcript_exposure_weights`.
- In `quant_from_buffer(...)`, remove the `FractionalCutoverPending` raise.
  The function will be the real EM driver again in M4; in M1 it can return
  a `NotImplementedError` with a precise message so the test surface is
  honest about the staged work.
- Drop the calibration TYPE_CHECKING import block referencing
  `_regional_exposure`.

`src/rigel/cli.py`:

- Delete the `FractionalCutoverPending` import, except-clause,
  `_write_cutover_summary(...)`, and the CLI flags
  `--regional-exposure`, `--no-regional-exposure`,
  `--regional-exposure-reference-quantile`.
- Drop the corresponding entries from the config-transform table.

`src/rigel/config.py`:

- Delete `regional_exposure_enabled` and
  `regional_exposure_reference_quantile` from `CalibrationConfig`.

`src/rigel/calibration/errors.py`:

- Delete the file. `FractionalCutoverPending` has no surviving consumer.

### 3.7 Delete legacy tests

Remove (whole files):

```text
tests/test_regional_exposure.py
tests/test_assemble_priors.py
tests/test_locus_prior_fused.py
tests/test_fractional_cutover_cleanup.py
tests/test_gdna_eff_len.py
tests/test_pipeline_integration_v6.py        # rewrite in M4
tests/test_calibrate_orchestrator.py         # rewrite to match trimmed calibrate()
```

The rewrites in M4 will be much smaller than the originals.

### 3.8 Acceptance for M1

- `pytest tests/` green after the deletions, except for tests deliberately
  removed in §3.7.
- `grep -R FractionalCutoverPending src/ tests/` returns nothing.
- `wc -l src/rigel/calibration/*.py` is under ~1500 total (was 4210).

## 4. M2 — Add `strand_deconv.py`

Goal: introduce one new module that owns the entire strand algorithm
(observations, two-pass kappa training, deconvolution, flags), and wire
its output through `CalibrationResult`. No other files added.

### 4.1 File and public surface

Create `src/rigel/calibration/strand_deconv.py` with exactly:

```python
@dataclass(frozen=True, slots=True)
class StrandRegionCounts:
    k_sense: np.ndarray         # float32[R]
    k_antisense: np.ndarray     # float32[R]
    n_total: np.ndarray         # float32[R]
    eligible: np.ndarray        # bool[R]
    p_r1_sense: float

@dataclass(frozen=True, slots=True)
class RegionGdnaEstimate:
    mean_count: np.ndarray      # float32[R]  ~ posterior E[D]
    upper_count: np.ndarray     # float32[R]  ~ 95% upper bound on D
    rna_lower_count: np.ndarray # float32[R]  ~ 95% lower bound on R
    precision: np.ndarray       # float32[R]  in [0, 1]
    flags: np.ndarray           # uint8[R]
    kappa_d: float
    kappa_d_n_training_regions: int
    kappa_d_n_exon_self_training: int
    p_r1_sense: float
    alpha: float                # 0.95 for v3

def build_strand_region_counts(...) -> StrandRegionCounts: ...
def estimate_kappa_d(...) -> StrandBalanceEstimate: ...
def screen_no_rna_exons(...) -> np.ndarray: ...        # bool mask
def deconvolve_regions_by_strand(...) -> RegionGdnaEstimate: ...
```

All private helpers are leading-underscore functions in the same file:
`_exact_region_posterior`, `_normal_region_bound`, `_beta_binomial_pmf`,
`_binomial_pmf`, `_apply_identifiability_guard`, `_compute_precision`,
`_pool_training_counts`.

### 4.2 Constants

Implementation-only constants (no CLI / config surface in M2):

```python
STRAND_ALPHA = 0.95
MAX_EXACT_POSTERIOR_N = 128

# kappa_d training
MIN_KAPPA_SEED_REGIONS = 20         # below this, fall back conservatively
CONSERVATIVE_KAPPA_FALLBACK = 20.0  # weak; widens RNA lower bound

# exon self-training screen
EXON_SCREEN_N_MIN = 10              # need this many fragments to test an exon
EXON_SCREEN_SPLICED_MAX = 1.0       # spliced evidence rules an exon out
EXON_SCREEN_R_FRAC = 0.05           # threshold on the RNA fraction
EXON_SCREEN_P_REJECT = 0.01         # accept exon only if P(R/n > R_FRAC) < this
```

Tunable later, never at the API surface.

### 4.3 Step 1 — observations

`build_strand_region_counts(region_arrays, payload_arrays, *, p_r1_sense)`:

1. For each of `contained`, `boundary_left`, `boundary_right`, call
   `sense_antisense_split(payload_arrays.region_counts_sorted,
   region_arrays.ts_class, compartment=..., splice=SPLICE_UNSPLICED)`.
2. Sum the three `sense` arrays into `k_sense`; same for antisense.
3. `n_total = k_sense + k_antisense`.
4. `eligible = ((ts_class == TS_POS) | (ts_class == TS_NEG)) & (n_total > 0)`.
5. Return `StrandRegionCounts`.

Spliced fields are not stored on the output (used only inside the exon
screen; computed on demand there).

### 4.4 Step 2 — kappa_d seed (intergenic + intron-only)

`_seed_kappa_d(region_arrays, payload_arrays)`:

1. `mask = is_intergenic(signature) | is_intron_only(signature)`.
2. Call `estimate_strand_balance(contained_unspliced_pos,
   contained_unspliced_neg, mask)` (kept from §3.2).
3. If `n_regions < MIN_KAPPA_SEED_REGIONS`, return the
   `CONSERVATIVE_KAPPA_FALLBACK` value tagged
   `fallback_reason="insufficient seed regions"`. Sparse evidence must
   *not* default to high-confidence binomial.

The returned object is a `StrandBalanceEstimate` (existing dataclass).

### 4.5 Step 3 — exon self-training screen

`screen_no_rna_exons(counts, region_arrays, payload_arrays,
strand_summary, kappa_d_seed) -> np.ndarray (bool[R])`:

For each region with `is_exon_any(signature) == True` and `ts_class ∈
{TS_POS, TS_NEG}`:

1. **Hard filters** (each on its own line of reasoning):
   - `n_total >= EXON_SCREEN_N_MIN`;
   - spliced sense+antisense mass `<= EXON_SCREEN_SPLICED_MAX`
     (computed once per call with one extra `sense_antisense_split` for
     `SPLICE_SPLICED`);
   - `has_both_strands(signature) == False` (strand-pure exon only).
2. **Probabilistic filter:** using the seed `kappa_d` and the normal-
   approximation tail (§4.7), compute
   `P(R / n_total >= EXON_SCREEN_R_FRAC | k_sense, n_total, p_r1_sense,
   kappa_d_seed)`. Accept the exon as no-RNA training only if this tail
   probability is below `EXON_SCREEN_P_REJECT`.

Return the boolean mask. Do not iterate. The screen runs once.

The screen reuses the same normal-approximation routine used by the
production deconvolution path, so this is one shared `_normal_region_bound`
implementation tested in two ways.

### 4.6 Step 4 — refit `kappa_d`

`estimate_kappa_d(...)`:

1. Compute the seed (§4.4).
2. Run the exon screen (§4.5).
3. Pool sufficient statistics across (intergenic ∪ intron-only ∪
   accepted exons). For the symmetric beta-binomial MoM,
   `S = Σ (k_sense − 0.5 n)²`, `B = Σ 0.25 n`,
   `A = Σ 0.25 n²` are additive across regions, so pooling is just three
   sums — exactly what `estimate_strand_balance` already does internally.
4. Return one `StrandBalanceEstimate` carrying:
   - the refit `kappa`;
   - `n_regions` = seed regions + accepted exons;
   - `fallback_used` only if both the seed was sparse *and* the exon
     screen accepted fewer than `MIN_KAPPA_SEED_REGIONS` exons.

If the screen accepts a `kappa_d` higher than the seed (more data,
tighter estimate), accept it. If the screen accepts a `kappa_d` lower
than the seed (less consistency), accept it; lower `kappa_d` widens RNA
lower bounds, which is the conservative direction.

Record `kappa_d_n_exon_self_training` (count of accepted exons) on the
final `RegionGdnaEstimate` so it appears in the JSON summary and tests.

### 4.7 Step 5 — per-region deconvolution

`deconvolve_regions_by_strand(counts, *, kappa_d, alpha=STRAND_ALPHA)`:

Model (already stated in v2 §5.3):

```text
D = n − R
DNA_sense | D ~ BetaBinomial(D, mean=0.5, concentration=kappa_d)
RNA_sense | R ~ Binomial(R, p_r1_sense)
k_sense = DNA_sense + RNA_sense
```

Two numerical paths, chosen per region:

- `n_total <= MAX_EXACT_POSTERIOR_N` → exact discrete posterior over `R ∈
  {0, …, round(n)}` by direct convolution. Uniform prior on `R`.
- otherwise → normal approximation with the model variance
  `Var(DNA_sense | D) = 0.25 D (D + kappa_d) / (1 + kappa_d)` and
  `Var(RNA_sense | R) = R p_r1_sense (1 − p_r1_sense)`, find the smallest
  `R` whose model can plausibly explain the observed `k_sense` at
  confidence `alpha`.

Identifiability guard: if `|p_r1_sense − 0.5| <
STRAND_CONTRAST_NUMERICAL_FLOOR`, set `rna_lower_count = 0`,
`mean_count = n_total`, `upper_count = n_total`, `precision = 0`,
`FLAG_NEAR_UNSTRANDED`.

### 4.8 Precision and flags

```python
strand_contrast = abs(p_r1_sense - 0.5)
overdispersion = 1.0 + n_total / (1.0 + kappa_d)
raw = (strand_contrast ** 2) * n_total / overdispersion
precision = raw / (1.0 + raw)             # in [0, 1]
```

Flags (uint8):

```python
FLAG_INELIGIBLE       = 1 << 0
FLAG_NEAR_UNSTRANDED  = 1 << 1
FLAG_KAPPA_FALLBACK   = 1 << 2
FLAG_APPROX_NORMAL    = 1 << 3
FLAG_EXON_SELF_TRAIN  = 1 << 4   # this region contributed to refit
```

### 4.9 Wire through the orchestrator

`calibrate(...)` body, after `compute_global_densities(...)`:

```python
counts = build_strand_region_counts(
    region_arrays, payload_arrays, p_r1_sense=strand_summary.p_r1_sense
)
kappa = estimate_kappa_d(region_arrays, payload_arrays, counts, strand_summary)
region_gdna = deconvolve_regions_by_strand(counts, kappa_d=kappa.kappa)
return build_calibration_result(
    payload=payload,
    scan_trained=scan_trained,
    global_densities=global_densities,
    fl_models=fl_models,
    region_gdna=region_gdna,
    strand_region_counts=counts,
)
```

`CalibrationResult` gains:

```python
region_gdna: RegionGdnaEstimate
strand_region_counts: StrandRegionCounts
```

(both required — there is no longer a no-strand-information code path
in v3; `StrandSummary.uninformative()` is still legal and will simply
produce all-`FLAG_NEAR_UNSTRANDED` regions).

### 4.10 Summary block

`CalibrationResult.to_summary_dict()` adds:

```json
"strand_deconv": {
  "alpha": 0.95,
  "p_r1_sense": 0.99,
  "kappa_d": 123.4,
  "kappa_d_n_seed_regions": 11873,
  "kappa_d_n_exon_self_training": 8412,
  "kappa_d_fallback": false,
  "n_regions_eligible": 612944,
  "n_regions_approx_normal": 671
}
```

### 4.11 Tests (added in M2)

One new file, `tests/test_strand_deconv.py`. Sections:

- `test_observations_ts_pos_and_ts_neg_basic` — sense/antisense
  semantics.
- `test_observations_ambiguous_and_intergenic_are_ineligible`.
- `test_seed_kappa_recovers_known_value_within_20pct`.
- `test_seed_kappa_falls_back_when_sparse`.
- `test_exon_screen_rejects_clear_rna` —
  `p_r1_sense=0.99`, 90/10 sense/antisense exon → rejected.
- `test_exon_screen_accepts_balanced_exon` — 50/50 with spliced=0 →
  accepted.
- `test_exon_screen_rejects_when_spliced_evidence_present` —
  even 50/50 unspliced + 2 spliced → rejected.
- `test_refit_kappa_uses_exon_self_training_when_available` —
  refit `n_regions` > seed `n_regions`.
- `test_deconvolution_protocol_symmetry` — reflected counts with
  `p_r1_sense ↔ 1 − p_r1_sense` produce the same outputs.
- `test_deconvolution_near_unstranded_is_conservative`.
- `test_deconvolution_exact_and_normal_agree_at_boundary` —
  evaluate at `n = MAX_EXACT_POSTERIOR_N ± 1`, agreement within 1%.
- `test_deconvolution_monotone_in_observed_imbalance`.

### 4.12 Acceptance for M2

- `pytest tests/test_strand_deconv.py -v` green.
- `pytest tests/` green except for tests intentionally rewritten in M4.
- `summary.json` carries the new `strand_deconv` block.
- No new modules under `src/rigel/calibration/` besides
  `strand_deconv.py`.

## 5. M3 — Minimal Prior Assembly

Goal: a single, small `prior.py` that consumes `RegionGdnaEstimate` and
emits the prior table batch EM needs. Replaces the deleted
`locus_prior.py`.

### 5.1 New file

`src/rigel/calibration/prior.py` with exactly:

```python
@dataclass(frozen=True, slots=True)
class LocusPrior:
    locus_id: int
    n_obs: float
    gdna_prior_count: float
    rna_prior_count: float
    enable_gdna: bool

@dataclass(frozen=True, slots=True)
class PriorTable:
    per_locus: tuple[LocusPrior, ...]
    gdna_prior_count: np.ndarray        # float64[L]
    rna_prior_count: np.ndarray         # float64[L]
    enable_gdna: np.ndarray             # uint8[L]

def assemble_priors(
    *,
    loci: Sequence[Locus],
    region_arrays: RegionArrays,
    region_gdna: RegionGdnaEstimate,
) -> PriorTable: ...
```

### 5.2 Algorithm

For each locus footprint `(ref, start, end)`:

1. Find region rows overlapping the footprint using the existing
   sorted/per-ref CSR in `RegionArrays`.
2. Sum `region_gdna.mean_count[overlapping]` → `gdna_prior_count`.
3. Sum `(region_gdna.n_total − region_gdna.mean_count)[overlapping]` →
   `rna_prior_count`.
4. `enable_gdna = gdna_prior_count > 0`.
5. Weight the contribution of each region by
   `region_gdna.precision[overlap_row]` when summing. Regions with
   `FLAG_NEAR_UNSTRANDED` or `FLAG_INELIGIBLE` contribute zero.

That is the whole prior. No multi-locus aggregation, no boundary policies,
no exon-vs-intron weighting. The strand model already accounts for the
hardest part (RNA vs gDNA per region); the prior is just a sum.

### 5.3 Tests

`tests/test_prior.py`:

- `test_locus_with_one_strand_pure_intron_recovers_gdna_count`.
- `test_locus_with_strong_rna_exon_gets_small_gdna_prior`.
- `test_locus_with_only_unstranded_regions_disables_gdna`.
- `test_overlapping_loci_share_region_contributions_correctly`.

### 5.4 Acceptance for M3

- Tests above pass.
- `prior.py` is < 200 lines.
- Nothing under `src/rigel/calibration/` references the old `PriorTable`
  shape.

## 6. M4 — Unblock `rigel quant`

Goal: remove the `FractionalCutoverPending` fail-fast in
`pipeline.quant_from_buffer(...)` and run EM end-to-end with the new
prior table.

### 6.1 Pipeline integration

`quant_from_buffer(...)`:

1. Build loci from `FragmentBuffer` (existing path; do not change).
2. `prior_table = assemble_priors(loci=loci, region_arrays=...,
   region_gdna=calibration.region_gdna)`.
3. Pass `prior_table.gdna_prior_count` and `prior_table.rna_prior_count`
   to the EM exactly where the old `PriorTable.gdna_prior_count` was
   consumed, before this file was stubbed.
4. Return a populated `AbundanceEstimator`.

The CLI no longer special-cases a cutover pending exception. Standard
error paths are sufficient.

### 6.2 Rewrites of the tests removed in §3.7

Add (replacements for the previously deleted tests, much smaller):

- `tests/test_calibrate.py` — covers the trimmed `calibrate(...)` API,
  asserts `region_gdna` and `fl_models` are populated, asserts no
  forbidden fields are present.
- `tests/test_pipeline_end_to_end_strand.py` — synthesizes a tiny BAM
  fixture, calls `run_pipeline(...)`, checks that an
  `AbundanceEstimator` is produced and that EM iterations reach
  convergence.

### 6.3 Acceptance for M4

- `rigel quant` runs end-to-end on the strand-specific capture fixture
  without raising.
- New end-to-end test green.
- `grep -R FractionalCutoverPending src/ tests/` returns nothing.
- All other existing tests still green.

## 7. M5 — Validation

After M4 is merged, run focused statistical validation before adding any
new features.

Validation cases (synthesized; cheap):

1. R1-sense (`p_r1_sense = 0.99`), pure-RNA exon at varied `n_total`:
   `rna_lower_count` covers true `R` with frequency ≥ 0.95.
2. R1-antisense (`p_r1_sense = 0.01`), reflected counts: identical
   coverage to (1). Protocol symmetry check.
3. Pure-gDNA region at varied `n_total`: `rna_lower_count == 0` with
   frequency ≥ 0.95.
4. Mixed region with known `R`, varied `kappa_d` and `p_r1_sense`:
   `mean_count` close to `D_true`; coverage holds.
5. Sparse intergenic data: conservative fallback widens
   `rna_lower_count`; coverage remains ≥ 0.95.
6. End-to-end on the lab's strand-specific capture truth data:
   per-locus gDNA fraction error reported.

These are validation, not unit tests; if any of them fails, fix the
model, not the test.

## 8. Final Calibration Package Shape

After M1–M4, `src/rigel/calibration/` is exactly these files:

```text
__init__.py            # public surface only
_arrays.py             # RegionArrays + PayloadArrays
_orchestrator.py       # calibrate(...)
_result.py             # CalibrationResult + builder
density_global.py      # trimmed: contained gDNA densities + strand balance
fl.py                  # fragment-length models (unchanged)
fractional_evidence.py # masks, ts_class, sense_antisense_split
prior.py               # NEW in M3: assemble_priors + PriorTable
regions.py             # index-time region partition (unchanged; loaded by index.py)
scan_payload.py        # CalibrationScanPayload (unchanged)
signature.py           # bit / channel constants (unchanged)
strand_deconv.py       # NEW in M2: the strand model
strand_summary.py      # StrandSummary (unchanged)
```

Twelve files. The four still-leading-underscore ones are genuine
package-internal modules (orchestrator, arrays, result, prior dataclass
shells), not the "private stub" pattern. The package is small enough that
a new contributor can read every file in a single sitting.

## 9. Ordered Implementation Checklist

Each step is a self-contained PR. Do not bundle.

1. **M1.a Source deletions** — remove the modules in §3.1. Trim
   `density_global.py`, `_result.py`, `__init__.py`, `_orchestrator.py`,
   `pipeline.py`, `cli.py`, `config.py`. Delete `errors.py`. Tests for
   the surviving APIs stay green.
2. **M1.b Test deletions** — remove the test files in §3.7.
3. **M2 Strand deconv** — add `strand_deconv.py` and
   `tests/test_strand_deconv.py`. Wire `region_gdna` and
   `strand_region_counts` into `CalibrationResult`. `summary.json`
   gains the `strand_deconv` block.
4. **M3 Prior assembly** — add `prior.py` and `tests/test_prior.py`.
5. **M4 Pipeline** — remove the `FractionalCutoverPending` raise from
   `pipeline.quant_from_buffer(...)`. Add the two replacement tests.
   `rigel quant` runs end-to-end.
6. **M5 Validation** — run §7. If anything fails, fix in the same PR
   that diagnoses it.

## 10. Out Of Scope For v3

- Capture BED ingestion and capture exposure modeling.
- nRNA-vs-mRNA per-locus decomposition (orthogonal; v3 only delivers
  gDNA vs RNA at the region level).
- Mappability correction.
- FFT-based posterior batching.
- Native (C++) strand-deconv kernel.
- Per-source kappa diagnostics dataclasses.
- 99% and 99.9% bounds (added later if prior calibration asks for them).
- Iterative exon self-training (more than one screen pass).

Anything on this list is a feature flag away after v3 lands; none of it
is required for the first working strand-specific capture quant.

## 11. The One-Question Filter

Before adding any new module, dataclass, public function, or parameter
during implementation, ask:

```text
Can the first working strand-specific capture quant be correct without it?
```

If yes, defer. The point of v3 is not the strand model — that is the
easy part. The point is to leave behind a calibration package a human
can read and trust.
