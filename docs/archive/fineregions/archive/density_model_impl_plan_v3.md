# gDNA Density Model Implementation Plan v3

Date: 2026-05-24
Status: implementation ready
Companion: `calibration_roadmap_v3.md`
Supersedes: `density_model_impl_plan_v1.md`, `density_model_impl_plan_v2.md`

## 0. Synthesis Summary

This plan turns the v1 architectural design and the v2 statistical refinements
into a step-by-step implementation. Every section here names concrete files,
functions, dataclasses, and tests. It also lists the existing code that is
stale, half-stale, or in the wrong module, and what to do with it.

The model in one paragraph: fit a Gamma background prior on contained
DNA-dominant anchors using exposure-adjusted method of moments. Cap the prior
opportunity so a captured boundary can move the posterior. For every region,
update the chosen prior with local unspliced boundary evidence. Predict the
held-out contained background with the Negative Binomial predictive
distribution and emit unbounded mean, upper, variance, tail tension, and a
relative-exposure scalar. Density evidence is independent of strand
deconvolution and EM until later phases consume it.

The implementation lands in six phases:

```text
Phase 0 — Cleanup of stale density and strand scaffolding
Phase 1 — RegionCountLedger memory-light evidence bundle
Phase 2 — DensityObservation: counts and FL-aware opportunities
Phase 3 — DensityModel: MoM Gamma prior, local update, NB predictive
Phase 4 — CalibrationResult / orchestrator wiring and summary
Phase 5 — RegionExposure becomes the density-derived surface
Phase 6 — Config, CLI, docs, and dead-script removal
```

No phase changes EM behavior. `quant_from_buffer` stays blocked at its current
`NotImplementedError`. Locus prior assembly and EM integration are out of
scope for this plan and belong to roadmap v3 Phase 5/6.

## 1. Decisions Locked from v1/v2

The following are locked and not re-litigated in v3:

- Background prior is fit only from contained anchors (INTERGENIC, INTRON-only).
- Boundary counts are local Bayesian evidence, never prior training data.
- Only unspliced channels enter the density likelihood.
- Estimator is exposure-adjusted **robust** method of moments — squared
  Pearson residuals with an upper-tail trim — so a handful of nRNA-bearing
  introns cannot blow up `phi` and collapse `beta`.
- Prior precision is bounded by a coefficient-of-variation floor
  (`beta_cap = 1 / (rho * CV_min^2)`), a depth-invariant statistical scale,
  not a magic boundary-count multiplier.
- Density evidence is computed without using `strand_summary`, `kappa_d`, or
  `RegionGdnaEstimate`.
- `RegionExposure` produced by density is a relative scalar that may exceed 1.
- v3 does not assemble locus priors or unblock `quant_from_buffer`.

The notation used throughout:

```text
C_r   = contained_unspliced_pos + contained_unspliced_neg
B_l   = boundary_left_unspliced_pos + boundary_left_unspliced_neg
B_r   = boundary_right_unspliced_pos + boundary_right_unspliced_neg
B_tot = B_l + B_r

L_c   = contained effective length under gDNA FL
L_l   = left boundary-side opportunity under gDNA FL
L_r   = right boundary-side opportunity under gDNA FL
L_b   = L_l + L_r
```

## 2. Inventory: What Exists and What Should Change

### 2.1 Keep as-is

- `calibration/_arrays.py` — `RegionArrays` and `PayloadArrays` are exactly
  what we need.
- `calibration/_exposure.fractional_boundary_side_exposure` — correct per-side
  boundary opportunity for any region length.
- `calibration/density_global.l_eff_contained` — keep but move (see §2.3).
- `calibration/fractional_evidence.is_intergenic`, `is_intron_only`,
  `transcript_strand_class`.
- `calibration/strand_deconv.py` — independent strand path; do not touch.
- `calibration/_diagnostics.py` and `calibration/_fl_sources.py`.
- `calibration/scan_payload.py`.

### 2.2 Half-stale, replace within Phase 0

`calibration/density_global.py` mixes three unrelated things:

1. `compute_global_densities` produces an `EXON-INTRON` slot that is hardcoded
   to `_empty_density("EXON-INTRON")` and an `exon_contained` slot that is
   never populated. These are legacy.
2. `GlobalGdnaDensity` carries dead fields: `strand_active`, `rho_uncorrected`,
   `strand_correction_fragments`, `n_strand_informative_regions`,
   `strand_informative_exposure_fraction`, `n_uninf_observed`,
   `n_fragments_estimated`, `n_rows_eligible`.
3. `estimate_strand_balance` and `StrandBalanceEstimate` model **strand**
   overdispersion, not density; they belong with `strand_deconv.py`.

Action in Phase 0:

- Move `estimate_strand_balance` and `StrandBalanceEstimate` to
  `calibration/strand_balance.py`. Re-export from `calibration/__init__.py`
  to preserve any existing imports.
- Delete `EXON-INTRON` and `exon_contained` paths from `GlobalDensityTable`
  and `compute_global_densities`. The remaining function returns a slim
  `ContainedRhoSummary` with INTERGENIC and INTRON rho, n_fragments,
  eff_length, and n_regions_used. This stays in `CalibrationResult` only as a
  short-term summary block alongside the new `density_evidence`.
- Mark `GlobalDensityTable`, `GlobalGdnaDensity`, `DensityType`, and the
  `_contained_density_for_mask`/`_empty_density` helpers as Phase 0
  replacement candidates. Phase 4 deletes them once `DensityEvidence` is the
  source of truth for the summary.

### 2.3 Move

- `l_eff_contained` should live next to `fractional_boundary_side_exposure` in
  `calibration/_exposure.py`. Re-export from the new `density_model.py` only
  if convenient. Update `density_global.py` (or its successor) and tests.

### 2.4 Replace

- `CalibrationResult.global_densities` becomes optional and is removed once
  Phase 4 lands.
- `RegionExposure.uniform` stops being the only producer; Phase 5 introduces a
  `RegionExposure.from_density(...)` constructor.

### 2.5 Delete

These debug scripts reference removed fields (`lambda_gdna`, `kappa`,
`exon_intron`, `EXON-INTRON`, `strand_active`) and are already broken against
the current Phase 4 API:

- `scripts/debug/diagnose_vcap_capture_exposure.py`
- `scripts/debug/analyze_vcap_exon_deconv_validation.py`
- `scripts/debug/inspect_zero_nrna_candidate_sets.py`
- `scripts/debug/diagnose_cap_at_1_exposure.py`
- `scripts/debug/analyze_density_computation.py`
- `scripts/debug/analyze_density_root_cause.py`
- `scripts/debug/test_gdna_calibration_sweep.py`
- `scripts/debug/analyze_synthetic_24_deep.py`
- `scripts/debug/analyze_vcap_gdna_false_rna_by_region.py`
- `scripts/debug/oracle_boundary_splice_tag_probe.py`

Confirm with the owner before deleting. If any are still wanted, port them to
the new `density_evidence` summary in the same PR as Phase 4. Do not preserve
broken scripts.

`src/rigel/sim/locus_sweep.py` references
`cal.global_densities.intergenic.lambda_gdna` (lines ~895/902/930). Port to
the new `density_evidence` summary or remove the block in Phase 4.

`src/rigel/sim/analysis.py` reads `EXON-INTRON` from summary dicts; replace
with `density_evidence.priors["INTRON"]` summary in Phase 4 or drop the
metric.

### 2.6 Test inventory

Currently in `tests/conftest.py::collect_ignore`:

- `tests/test_density_global.py` — replace with new tests in Phase 0; do not
  unlock the old file.
- `tests/test_exposure.py` — rewrite for `RegionExposure.from_density` in
  Phase 5.
- `tests/test_per_locus_gdna_mass.py`, `tests/test_ndarray_util.py`,
  `tests/test_calibration_result.py`, `tests/test_calibration_accumulator.py`,
  `tests/test_bayesian_prior_acceptance.py`, `tests/test_pipeline_wiring.py`,
  `tests/test_profiler.py` — out of scope for v3; explicitly leave skipped.

Active and useful:

- `tests/test_strand_deconv.py` — untouched.
- `tests/test_fractional_density_global.py` — must continue to pass against
  whatever survives in `compute_global_densities` after Phase 0, or be
  retired in Phase 4 when that function is removed.
- `tests/test_calibrate.py`, `tests/test_fl_eff_len_cache.py`.

## 3. Module Layout (Target)

After all six phases land, `calibration/` looks like:

```text
calibration/
  __init__.py
  _arrays.py
  _diagnostics.py
  _exposure.py             # contained_exposure_clipped,
                           # fractional_boundary_side_exposure,
                           # l_eff_contained (moved here)
  _fl_sources.py
  _orchestrator.py         # FL -> ledger -> obs -> evidence -> result
  _result.py               # density_evidence is the gDNA channel
  evidence.py              # RegionCountLedger (Phase 1)
  density_observation.py   # DensityObservation builder (Phase 2)
  density_model.py         # GammaRatePrior, DensityEvidence, fit + predict (Phase 3)
  exposure.py              # RegionExposure incl. from_density (Phase 5)
  fl.py
  fractional_evidence.py
  regions.py
  scan_payload.py
  signature.py
  strand_balance.py        # StrandBalanceEstimate, estimate_strand_balance (moved)
  strand_deconv.py
  strand_summary.py
```

`density_global.py` is removed.

## 4. Phase 0 — Cleanup (mergeable on its own)

Goal: make the surface ready for the density model without changing observable
calibration outputs.

Steps:

1. Add `calibration/strand_balance.py`. Move `StrandBalanceEstimate` and
   `estimate_strand_balance` from `density_global.py`. Re-export from
   `calibration/__init__.py`.
2. Move `l_eff_contained` to `_exposure.py`. Keep a deprecation re-export in
   `density_global.py` for one PR cycle.
3. Trim `GlobalGdnaDensity` to fields still consumed: `type`, `rho`,
   `n_fragments`, `eff_length_bp`, `n_regions_used`. Remove strand-prefixed
   and `rho_uncorrected`-style fields. Update `to_summary_dict`.
4. Remove the `exon_intron` and `exon_contained` slots from
   `GlobalDensityTable`. `compute_global_densities` returns `intergenic`,
   `intron`, `gdna_fl`, `strand_balance` only.
5. Audit imports of moved/removed names across `src/rigel` and `tests/`. The
   only legitimate consumers are the orchestrator and `_result.py`.
6. Confirm `scripts/debug/*` deletion list (see §2.5) with the owner. Delete
   in this PR or list as a follow-up.
7. Run `pytest tests/ -v`. The pre-existing `test_calibration.py::TestStrandLLR`
   failure is allowed; nothing else should regress.

New tests in Phase 0:

- `tests/test_strand_balance.py` — coverage moved from any latent
  `test_density_global` strand-balance assertions.
- No new density tests yet.

Exit criteria:

- `density_global.py` contains only `compute_global_densities` and the trimmed
  `GlobalGdnaDensity`/`GlobalDensityTable`.
- All non-debug tests pass.

## 5. Phase 1 — RegionCountLedger

Goal: one memory-light bundle of fractional counts so density, strand, and
future fusion all read from the same source.

New file: `calibration/evidence.py`.

```python
@dataclass(frozen=True, slots=True)
class RegionCountLedger:
    contained_unspliced_pos: np.ndarray   # float32[R]
    contained_unspliced_neg: np.ndarray
    boundary_left_unspliced_pos: np.ndarray
    boundary_left_unspliced_neg: np.ndarray
    boundary_right_unspliced_pos: np.ndarray
    boundary_right_unspliced_neg: np.ndarray
    contained_spliced_pos: np.ndarray
    contained_spliced_neg: np.ndarray
    boundary_left_spliced_pos: np.ndarray
    boundary_left_spliced_neg: np.ndarray
    boundary_right_spliced_pos: np.ndarray
    boundary_right_spliced_neg: np.ndarray
```

Built from `PayloadArrays.region_counts_sorted` by `signature.channel_index`.
The 12 channels already exist; this dataclass is a typed view, not a copy.

Derived totals as methods:

```python
def contained_unspliced_total(self, out=None) -> np.ndarray
def boundary_left_unspliced_total(self, out=None) -> np.ndarray
def boundary_right_unspliced_total(self, out=None) -> np.ndarray
def boundary_unspliced_total(self, out=None) -> np.ndarray
def spliced_total(self, out=None) -> np.ndarray
```

Notes:

- Each array is a view into `region_counts_sorted`. Do not allocate new
  storage in the constructor.
- Reductions return `float64`.

Builder:

```python
def build_region_count_ledger(payload_arrays: PayloadArrays) -> RegionCountLedger
```

`PayloadArrays.contained_unspliced_total` etc. become thin wrappers around
`RegionCountLedger.contained_unspliced_total`. The duplicated unspliced-total
materialization in `PayloadArrays.from_payload` can be deleted at this point
or in Phase 4.

Tests: `tests/test_region_count_ledger.py`.

- Channels are correctly indexed for every signature.
- Totals are bitwise equal to direct slice sums.
- The ledger does not copy `region_counts_sorted`.
- `TS_NONE` / `TS_AMBIG` rows retain observed counts.

Exit criteria: orchestrator and strand path both go through the ledger.

## 6. Phase 2 — DensityObservation

Goal: a sorted, frozen, per-region table of gDNA-relevant counts and FL-aware
opportunities. No model.

New file: `calibration/density_observation.py`.

```python
@dataclass(frozen=True, slots=True)
class DensityObservation:
    # counts (float32[R])
    contained_count: np.ndarray
    boundary_left_count: np.ndarray
    boundary_right_count: np.ndarray
    boundary_count: np.ndarray
    observed_compatible_count: np.ndarray

    # opportunities (float64[R])
    contained_leff: np.ndarray
    boundary_left_leff: np.ndarray
    boundary_right_leff: np.ndarray
    boundary_leff: np.ndarray

    # classification (bool[R])
    anchor_intergenic: np.ndarray
    anchor_intron: np.ndarray
    is_anchor: np.ndarray

    # diagnostics
    spliced_count: np.ndarray
    region_length: np.ndarray
```

Builder:

```python
def build_density_observation(
    region_arrays: RegionArrays,
    ledger: RegionCountLedger,
    gdna_fl: FragmentLengthModel,
) -> DensityObservation
```

Implementation:

- Counts come from the ledger's unspliced totals (POS + NEG).
- `contained_leff = l_eff_contained(end - start, gdna_fl)`.
- `boundary_left_leff = boundary_right_leff =
  fractional_boundary_side_exposure(end - start, gdna_fl)`. Both sides share
  the same per-region side opportunity since the helper depends only on the
  receiving-region length.
- `anchor_intergenic = is_intergenic(signature)`,
  `anchor_intron = is_intron_only(signature)`, `is_anchor` is their union.
- `observed_compatible_count = contained_count + boundary_count`.
- `spliced_count = contained_spliced_total + boundary_spliced_total` from the
  ledger. Stored as a diagnostic; never used by the model.

Tests: `tests/test_density_observation.py`.

- Anchor masks are exact for known signatures.
- Contained `L_c` reproduces `l_eff_contained` for hand-computed FL PMFs.
- Boundary `L_side(S)` reproduces `sum_ell h(ell) * min((ell - 1) / 2, S / 2)`.
- Short regions get small but finite boundary opportunity and zero contained
  opportunity when `S < gdna_fl.min_positive`.
- Spliced counts have no effect on any non-spliced field.

Exit criteria: `DensityObservation` builds successfully on every existing
test fixture and on the calibration scenario BAMs.

## 7. Phase 3 — DensityModel

Goal: turn `DensityObservation` into `DensityEvidence` via MoM Gamma prior,
local update, and Negative Binomial predictive.

New file: `calibration/density_model.py`.

### 7.1 Dataclasses

```python
@dataclass(frozen=True, slots=True)
class GammaRatePrior:
    family: str               # "INTERGENIC" | "INTRON" | "ALL"
    alpha: float
    beta: float
    mean_density: float       # alpha / beta
    phi: float                # 1 / alpha
    beta_raw: float
    beta_cap: float
    cap_applied: bool
    n_regions: int
    n_fragments: float
    eff_length: float
    residual_sum: float                  # S_trim
    poisson_variance_sum: float          # B_trim
    extra_variance_basis_sum: float      # A_trim
    trim_upper: float                    # fraction of residuals discarded
    n_trimmed: int                       # rows excluded from the variance fit
    trimmed_mu_fraction: float           # sum(mu_discarded) / sum(mu_all)
    fallback_depth: int
    fit_status: str
```

```python
@dataclass(frozen=True, slots=True)
class DensityEvidence:
    # predictive
    mean_unbounded: np.ndarray            # float64[R]
    upper_unbounded: np.ndarray
    variance_unbounded: np.ndarray
    precision_unbounded: np.ndarray

    # posterior
    alpha_post: np.ndarray
    beta_post: np.ndarray
    rho_post: np.ndarray
    relative_exposure: np.ndarray         # rho_post / rho_ref

    # local information
    w_boundary_opportunity: np.ndarray
    w_boundary_count: np.ndarray

    # tension diagnostics
    tail_probability: np.ndarray
    expected_tail_count: np.ndarray
    density_over_observed_ratio: np.ndarray

    # provenance and quality
    prior_family: np.ndarray              # uint8[R]: 0=INTERGENIC, 1=INTRON, 2=ALL
    fallback_depth: np.ndarray            # uint8[R]
    flags: np.ndarray                     # uint8[R]
    confidence: float

    # priors and reference
    priors: dict[str, GammaRatePrior]
    rho_ref: float
```

### 7.2 Background prior fit

Estimator:

```python
def fit_gamma_prior(
    counts: np.ndarray,
    opportunities: np.ndarray,
    *,
    family: str,
    beta_cap: float,
    phi_floor: float = DENSITY_PHI_FLOOR,
    trim_upper: float = DENSITY_PHI_TRIM_UPPER,
    min_regions: int = DENSITY_MIN_PRIOR_REGIONS,
    min_eff_length: float = DENSITY_MIN_EFF_LENGTH,
    fallback_depth: int = 0,
) -> GammaRatePrior
```

Rationale: the INTRON anchor mask is signature-only and does not screen for
expression. A small minority of transcribed introns will carry nascent-RNA
counts that look like extreme upper-tail outliers under the Gamma-Poisson
null. A naive sum of squared residuals lets a handful of those rows inflate
`phi_hat`, drive `alpha_raw -> 0`, and collapse `beta_raw -> 0`, leaving the
prior with no influence at all. v3 uses an upper-trimmed estimator so a few
outliers cannot kill the background.

Robust algorithm (upper-trimmed Pearson MoM):

```text
eligible = (opportunities >= min_eff_length) & isfinite(counts)
if eligible.sum() < min_regions:
    raise InsufficientAnchors

C   = counts[eligible].astype(float64)
L   = opportunities[eligible].astype(float64)

rho = C.sum() / L.sum()
mu  = rho * L

# Pearson-style squared residuals (variance-stabilising scale).
r2  = (C - mu) ** 2

# Discard the top `trim_upper` fraction of squared residuals. The retained
# rows define a contiguous lower-tail subset; sum mu and mu^2 over the same
# subset so the (S - B) / A identity stays exposure-adjusted.
k_keep = max(min_regions, int(ceil((1.0 - trim_upper) * r2.size)))
keep   = argpartition(r2, k_keep - 1)[:k_keep]

S_trim = r2[keep].sum()
B_trim = mu[keep].sum()
A_trim = (mu[keep] ** 2).sum()

phi_hat   = max(0.0, (S_trim - B_trim) / max(A_trim, eps))
alpha_raw = 1.0 / max(phi_hat, phi_floor)
beta_raw  = alpha_raw / max(rho, eps)

beta  = min(beta_raw, beta_cap)
alpha = rho * beta
```

Notes:

- Trimming is one-sided (upper). Under Gamma-Poisson the lower tail rarely
  produces leverage; under nRNA contamination only the upper tail is
  inflated. Two-sided trimming would discard legitimate quiet regions.
- The Pearson scaling `r2 / mu` would be more standard, but plain `(C - mu)^2`
  preserves the closed-form `(S - B) / A` identity. Using the same trimming
  index for `S`, `B`, and `A` keeps the estimator unbiased to first order
  under the Gamma-Poisson null.
- Default `DENSITY_PHI_TRIM_UPPER = 0.05` (drop the top 5% of squared
  residuals). MAD on `(C - mu) / sqrt(mu)` is a documented alternative; we
  prefer trimming because it preserves the existing variance identity and is
  trivial to vectorise.
- `rho`, `n_fragments = C.sum()`, and `eff_length = L.sum()` are computed on
  the full eligible set; only the variance estimator is trimmed.
- The prior also records `n_trimmed`, the trimmed top-residual sum, and the
  fraction of total `mu` discarded so dataset-level nRNA contamination is
  observable post-hoc.

Return a `GammaRatePrior` with diagnostic moments preserved (full and trimmed).

`beta_cap` is a depth-invariant CV floor on the prior, not a multiplier of
boundary opportunities:

```python
def compute_beta_cap(rho: float) -> float:
    """Cap prior precision so the prior CV cannot fall below DENSITY_PRIOR_MIN_CV.

    Gamma(alpha, beta) has CV = 1 / sqrt(alpha) = 1 / sqrt(rho * beta).
    Requiring CV >= cv_min is equivalent to beta <= 1 / (rho * cv_min^2).
    With cv_min = 0.05 this gives beta_cap = 400 / rho — a depth-invariant
    scale that is identical for shallow and ultra-deep libraries.
    """
    if rho <= 0.0 or not np.isfinite(rho):
        return float("inf")
    return float(DENSITY_PRIOR_MAX_PRECISION / rho)
```

`beta_cap` is computed per family using that family's `rho`. The all-anchor
fallback prior uses its own `rho` to derive its own cap.

### 7.3 Fallback selection

For every region:

```text
0: family prior for the anchor class (INTERGENIC or INTRON)
1: all-anchor prior fit across INTERGENIC ∪ INTRON
2: broad weak fallback = Gamma(rho_ref * beta_cap, beta_cap)
```

Each region records the depth used. Non-anchor regions (exon, mixed) always
start at depth 1 because they have no native family.

### 7.4 Local update and predictive

```text
alpha_post = alpha_prior + B_tot
beta_post  = beta_prior  + L_b
rho_post   = alpha_post / beta_post

mean_c     = rho_post * L_c
var_c      = mean_c * (1 + L_c / beta_post)
p_nb       = beta_post / (beta_post + L_c)
upper_c    = scipy.stats.nbinom.ppf(confidence, alpha_post, p_nb)

mean_unbounded  = B_tot + mean_c
upper_unbounded = B_tot + upper_c
var_unbounded   = var_c
```

Quality and tension:

```text
w_boundary_opportunity = L_b / (beta_prior + L_b)
w_boundary_count       = where(B_tot > 0, B_tot / (alpha_prior + B_tot), 0.0)

tail_probability       = 1 - scipy.stats.nbinom.cdf(C_r, alpha_post, p_nb)
expected_tail_count    = E[(N_c^gdna - C_r)+ | boundary]
density_over_observed_ratio = mean_unbounded / max(B_tot + C_r, eps)

cv                  = sqrt(var_unbounded) / max(mean_unbounded, eps)
precision_unbounded = 1.0 / (1.0 + cv)
```

`expected_tail_count` initial implementation: stable truncated sum for small
means and a normal-approximation closed form for large means. Choose the
threshold by mean magnitude, not by `alpha_post`.

### 7.5 Flags

```python
FLAG_LOW_CONTAINED_OPPORTUNITY = 0x01
FLAG_LOW_BOUNDARY_OPPORTUNITY  = 0x02
FLAG_PRIOR_DOMINATED           = 0x04
FLAG_HIGH_TAIL_TENSION         = 0x08
FLAG_NON_ANCHOR                = 0x10
FLAG_FALLBACK_USED             = 0x20
```

Rules:

```text
FLAG_LOW_CONTAINED_OPPORTUNITY <- L_c < DENSITY_MIN_EFF_LENGTH
FLAG_LOW_BOUNDARY_OPPORTUNITY  <- L_b < DENSITY_MIN_BOUNDARY_OPPORTUNITY
FLAG_PRIOR_DOMINATED           <- w_boundary_opportunity < DENSITY_MIN_BOUNDARY_INFO
FLAG_HIGH_TAIL_TENSION         <- tail_probability > DENSITY_TAIL_PROBABILITY_WARN
FLAG_NON_ANCHOR                <- not is_anchor
FLAG_FALLBACK_USED             <- fallback_depth > 0
```

### 7.6 Public entry point

```python
def fit_density_evidence(
    observation: DensityObservation,
    *,
    confidence: float = 0.95,
    min_eff_length: float = DENSITY_MIN_EFF_LENGTH,
) -> DensityEvidence
```

Internal constants in `density_model.py`:

```python
# Prior precision floor: CV of the Gamma prior cannot fall below 5%.
DENSITY_PRIOR_MIN_CV             = 0.05
DENSITY_PRIOR_MAX_PRECISION      = 1.0 / (DENSITY_PRIOR_MIN_CV ** 2)  # = 400.0

# Robust MoM controls.
DENSITY_PHI_FLOOR                = 1.0e-6
DENSITY_PHI_TRIM_UPPER           = 0.05    # drop top 5% squared residuals

# Opportunity / information thresholds.
DENSITY_MIN_EFF_LENGTH           = 1.0
DENSITY_MIN_BOUNDARY_OPPORTUNITY = 1.0
DENSITY_MIN_BOUNDARY_INFO        = 0.05
DENSITY_MIN_PRIOR_REGIONS        = 20

# Tail tension: under a well-calibrated model the per-region tail
# probability is uniformly distributed, so the median region sits near 0.5.
# Flagging is reserved for genuinely surprising regions: those where the
# model predicts substantially more contained background than was observed
# (suggesting unmodelled depletion or capture-specific masking).
DENSITY_TAIL_PROBABILITY_WARN    = 0.95
```

### 7.7 Tests

`tests/test_density_model.py`:

- `fit_gamma_prior` recovers mean density on Poisson anchors.
- `fit_gamma_prior` reports finite `phi` on overdispersed simulated anchors.
- `phi_hat == 0` falls back to `phi_floor` instead of producing `inf` alpha.
- **Robust trim (nRNA simulation):** inject 1–5% of anchors with synthetic
  high-count outliers (`C_i = 50 * mu_i`) into an otherwise Poisson anchor
  set. The trimmed estimator returns `phi_hat` and `beta_raw` within a small
  tolerance of the untrimmed Poisson fit. A standard SSE estimator on the
  same input collapses `beta_raw` toward zero; assert that v3 does not.
- `compute_beta_cap(rho)` returns `400 / rho` for `rho > 0` and `inf` for
  `rho == 0`. Scaling counts and opportunities by the same constant leaves
  `beta_cap` unchanged (depth invariance).
- `beta_raw > beta_cap` preserves `alpha / beta == rho` after capping and
  yields prior CV `>= DENSITY_PRIOR_MIN_CV` to within float tolerance.
- Sparse anchors trigger fallback depth 1 then depth 2 deterministically.
- Local update is correct: zero boundary leaves posterior equal to prior.
- High boundary count produces `rho_post > rho_ref`.
- NB predictive mean and variance match `scipy.stats.nbinom` for integer cases.
- `precision_unbounded` is monotone in `var_unbounded` at fixed mean.
- **Tension threshold:** a region whose contained count matches its mean
  prediction has `tail_probability` near 0.5 and is NOT flagged with
  `FLAG_HIGH_TAIL_TENSION`. A region with `mean_unbounded >> observed` does
  trigger the flag.
- Flag bits are set correctly on hand-crafted regions.

Exit criteria: `fit_density_evidence` produces finite output on the existing
calibration scenarios, with no exceptions on production-equivalent data.

## 8. Phase 4 — CalibrationResult and Orchestrator

Goal: density evidence becomes a first-class field in `CalibrationResult`,
the orchestrator builds it, and the legacy `global_densities` block is
removed.

### 8.1 Changes to `_result.py`

```python
@dataclass(frozen=True, slots=True)
class CalibrationResult:
    fl_models: FLModels
    diagnostics: Diagnostics
    region_gdna: "RegionGdnaEstimate"
    region_exposure: "RegionExposure"
    density_evidence: "DensityEvidence"
    density_observation: "DensityObservation"
    n_multi_loci: int = 0
    rna_lower_confidence: float = 0.95
```

- Remove `global_densities`.
- `to_summary_dict()` no longer emits a `global_densities` block. It emits a
  new `density_evidence` block with prior families, cap usage, flag counts,
  and tension summaries.
- `build_calibration_result` requires `density_evidence` and
  `density_observation` arguments.

### 8.2 Changes to `_orchestrator.py`

```text
1. Build RegionArrays and PayloadArrays as today.
2. Build RegionCountLedger from PayloadArrays.
3. Build FL models as today.
4. Run strand path as today (independent).
5. Build DensityObservation from region_arrays + ledger + fl_models.gdna.
6. Fit DensityEvidence with confidence = config.gdna_density_confidence.
7. Build RegionExposure from density (see Phase 5).
8. Assemble CalibrationResult.
```

The orchestrator no longer calls `compute_global_densities`. The function and
the file `density_global.py` are deleted in this phase.

### 8.3 Downstream consumers

Replace every reading of `cal.global_densities` with `cal.density_evidence`.
Confirmed sites and ports:

- `src/rigel/sim/locus_sweep.py` — port the `intergenic.lambda_gdna` printout
  to `density_evidence.priors["INTERGENIC"].mean_density`.
- `src/rigel/sim/analysis.py` — port `EXON-INTRON` reads to the new tension
  summaries; if the metric is not reproducible, document and drop.
- Test files in `tests/golden/` may include `global_densities` JSON. Regenerate
  with `pytest tests/ --update-golden`.

### 8.4 Summary schema

```text
density_evidence:
  confidence: float
  n_regions: int
  n_anchor_intergenic: int
  n_anchor_intron: int
  rho_ref: float
  priors:
    INTERGENIC: { alpha, beta, mean_density, phi, beta_raw, beta_cap,
                  cap_applied, n_regions, fit_status }
    INTRON:     { ... }
    ALL:        { ... }   # all-anchor fallback prior, may be null
  flags:
    n_low_contained_opportunity: int
    n_low_boundary_opportunity: int
    n_prior_dominated: int
    n_high_tail_tension: int
    n_fallback_used: int
  relative_exposure: { min, p50, p95, max, mean }
  tail_probability:   { p50, p95, max }
  mean_unbounded:     { p50, p95, max }
```

Avoid serializing per-region arrays.

### 8.5 Tests

- `tests/test_calibrate.py` — assert `density_evidence` is populated;
  strand-only changes do not move density evidence on the same payload.
- Golden output regeneration if needed.

Exit criteria:

- `density_global.py` is deleted.
- `compute_global_densities` is gone.
- All non-debug tests pass.

## 9. Phase 5 — RegionExposure from Density

Goal: replace the `RegionExposure.uniform` scaffold with a density-derived
surface that is wired into `CalibrationResult.region_exposure`. EM still does
not consume this; only the summary and downstream work do.

### 9.1 Constructor

`calibration/exposure.py` gains:

```python
@classmethod
def from_density(
    cls,
    density_evidence: "DensityEvidence",
    *,
    max_exposure: float | None = None,
) -> "RegionExposure":
    """Build a relative-exposure surface from posterior density."""
```

Implementation:

```text
A_r       = density_evidence.relative_exposure  # may exceed 1
rho_r     = density_evidence.rho_post
eligible  = density_evidence.flags & FLAG_LOW_BOUNDARY_OPPORTUNITY == 0
if max_exposure is not None:
    A_r = minimum(A_r, max_exposure)
rho_ref   = density_evidence.rho_ref
```

`RegionExposure` docstring and dataclass are updated to make explicit that
`A_r` can exceed 1. Remove the `(0, 1]` contract from the current docstring.

`RegionExposure.uniform(R)` remains as a fallback for code paths that have no
density evidence (synthetic tests). It is not called by the orchestrator.

### 9.2 Tests

- `tests/test_region_exposure_from_density.py` — verify alignment with
  `density_evidence.relative_exposure`, eligibility wiring, and clipping.
- Rewrite `tests/test_exposure.py` for the new contract; remove from
  `collect_ignore` in `conftest.py`.

Exit criteria:

- `RegionExposure.from_density` is the orchestrator's exposure source.
- `region_exposure.to_summary_dict()` reports density-derived statistics.

## 10. Phase 6 — Config, CLI, Docs

### 10.1 Config

`src/rigel/config.py` gains:

```python
@dataclass(frozen=True)
class CalibrationConfig:
    ...
    gdna_density_confidence: float = 0.95
    density_min_eff_length: float = 1.0
```

Validate `0.5 <= gdna_density_confidence < 1.0`.

Internal constants stay in `density_model.py`.

### 10.2 CLI

`src/rigel/cli.py`:

```text
--gdna-density-confidence FLOAT (default 0.95)
    Posterior confidence for the per-region gDNA NB upper bound.
```

Wire via `_ParamSpec("gdna_density_confidence",
"calibration.gdna_density_confidence")`.

The other thresholds (`density_min_eff_length`, internal flag thresholds) are
not exposed yet.

### 10.3 Docs

- `docs/MANUAL.md` — short subsection on the new density confidence flag.
- `docs/fineregions/calibration_roadmap_v3.md` — mark Phase 1–4 (this plan)
  as complete on landing; keep Phase 5 (fusion) and beyond as open.
- Move `density_model_impl_plan_v1.md` and `density_model_impl_plan_v2.md`
  into `docs/fineregions/archive/` once v3 lands so v3 is the obvious entry
  point.

## 11. Removal Checklist (After Phase 6)

Things that should no longer exist anywhere:

- `calibration/density_global.py`
- `GlobalGdnaDensity`, `GlobalDensityTable`, `DensityType`,
  `_contained_density_for_mask`, `_empty_density`
- `DENSITY_PRIOR_EQUIV_BOUNDARIES` (never landed; replaced by
  `DENSITY_PRIOR_MIN_CV` / `DENSITY_PRIOR_MAX_PRECISION`)
- `CalibrationResult.global_densities`
- `to_summary_dict()` keys `"global_densities"`, `"EXON-INTRON"`,
  `"exon_contained"`
- Legacy fields `lambda_gdna`, `strand_active`, `rho_uncorrected`,
  `strand_correction_fragments`, `n_strand_informative_regions`,
  `strand_informative_exposure_fraction`, `n_uninf_observed`,
  `n_fragments_estimated`, `n_rows_eligible`
- The `scripts/debug/*` files listed in §2.5 (or their ported successors)

Verification grep:

```text
grep -rEn "global_densities|EXON-INTRON|exon_contained|lambda_gdna" src/ tests/
```

After Phase 6 this should return zero non-archive hits.

## 12. Risks and Unexpected Aspects

1. Golden outputs change. Mitigation: explicit `--update-golden` step in the
   Phase 4 PR with reviewer-visible diff.
2. Memory increase from per-region float64 evidence arrays. Mitigation:
   `slots=True`; the evidence dataclass has ~10 per-region arrays. At 10^7
   regions this is on the order of 1 GB. If too large, switch
   `mean_unbounded`, `upper_unbounded`, and `variance_unbounded` to float32
   after measurement.
3. NB quantile cost for very large `R`. Mitigation: vectorize
   `scipy.stats.nbinom.ppf`; if it becomes a hotspot, add a saddlepoint or
   normal approximation path keyed off `alpha_post` and `L_c / beta_post`.
4. `RegionExposure.A_r` can exceed 1 once it is density-derived. This breaks
   the current docstring contract and any caller assuming a probability.
   Audit every reader of `RegionExposure.A_r` before Phase 5 lands; this is
   the highest-risk silent-behavior change in this plan.
5. Beta cap behavior varies between datasets. Mitigation: cap is a CV floor
   (`beta_cap = 1 / (rho * CV_min^2)`), which is depth invariant — shallow
   and ultra-deep libraries get the same maximum prior precision relative to
   their own density. `cap_applied` is exposed in the prior diagnostics so
   dataset behavior is visible.
6. Intronic anchors contaminated by nascent RNA. Mitigation: the trimmed
   variance estimator absorbs sparse upper-tail outliers without collapsing
   the prior. Per-family diagnostics record the trimmed-mass fraction so
   systematic contamination is observable, and downstream fusion / robust
   anchor filtering can still act on it. v3 does not silently exclude
   contaminated anchors.
7. Capture-mode boundary leakage into INTRON-only anchors could elevate the
   background prior. The `beta_cap` keeps the prior movable, but the global
   `rho` will still drift. Record the per-family `n_fragments` and
   `eff_length` so post-hoc analysis can detect this.
8. Fractional counts in MoM. Sums of squared residuals remain well-defined
   for fractional `C_i`; the bias is small and the cap absorbs the worst
   cases.

## 13. Things v3 Intentionally Does Not Do

- Numerical NB MLE for the background prior.
- Reference-specific priors.
- Capture-mode hard classifier.
- Spliced channels feeding the density numerator.
- Density-strand probabilistic fusion (separate roadmap phase).
- Locus-level prior assembly or `quant_from_buffer` unblocking.
- Mappability-corrected denominators.

These belong to later roadmap phases and have their own design surface.

## 14. Done Definition

v3 is done when:

- All six phases are merged.
- `CalibrationResult.density_evidence` is populated for every calibration
  run.
- `RegionExposure.from_density` is the only exposure source in the
  orchestrator path.
- `density_global.py` and the listed legacy fields/scripts no longer exist.
- `pytest tests/ -v` passes except for the pre-existing
  `tests/test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna`
  failure.
