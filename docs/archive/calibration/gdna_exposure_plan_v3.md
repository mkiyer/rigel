# gDNA Regional Exposure - Implementation Plan v3

**Status**: implemented (R1-R6 review fixes folded in; 13-step plan landed, all tests passing).
**Date**: 2026-05-18
**Supersedes**: `gdna_exposure_plan_v2.md`.
**Companion review**: `gdna_exposure_plan_v2_review.md`.

---

## 1. Summary

Hybrid-capture gDNA is not uniformly exposed across the genome. The v6
model currently normalizes the gDNA component by an unweighted genomic
opportunity length, so large off-target introns and mega-loci can make a
true gDNA fragment look too cheap to explain as gDNA and too easy to route
to RNA.

v3 learns a per-region relative gDNA exposure weight `A_r in (0, 1]` from
the conservative calibration evidence already collected by v6, then applies
that same weight to both sides of the gDNA likelihood:

```text
gdna_log_liks[u] += log A(ref_u, midpoint_u)      # numerator
L_g[locus]       = sum over valid gDNA starts of A(midpoint)  # denominator
```

The EM solver is unchanged. The gDNA prior count used by EM is unchanged in
v3; a regional prior alternative is emitted only as a diagnostic.

---

## 2. Locked Decisions

1. **Numerator and denominator ship together.** There is no enabled
   denominator-only production mode.
2. **The prior is unchanged.** `gdna_prior_count` remains the existing v6
   global-only expected count. Regional prior counts are diagnostic only.
3. **One CLI flag.** `--regional-exposure {auto,off}` is the only user-facing
   surface.
4. **Flat config field.** Use `CalibrationConfig.regional_exposure_enabled`
   rather than a nested config dataclass; the current CLI config mapper is
   one level deep.
5. **No midpoint in per-locus partitions.** `genomic_midpoint` is a global
   `ScoredFragments` field used before `partition_and_free()`; it is not
   scattered into `LocusPartition`.
6. **Cross-reference units are skipped for the numerator.** If one EM unit's
   candidate transcripts span multiple references, per-unit `log A` is not
   applied to that unit and the skip is counted.
7. **Auto-uniform attenuation is deterministic.** It uses a closed-form null
   spread approximation, no random Poisson fallback in v3.

---

## 3. Files Touched

| Path | Change |
|---|---|
| `src/rigel/calibration/_regional_exposure.py` | New module: exposure table builder, weighted quantiles, vectorized lookups, regional diagnostic prior helper. |
| `src/rigel/calibration/_arrays.py` | Extend `RegionArrays` with strand and `PayloadArrays` with orientation-resolved intron/boundary arrays. |
| `src/rigel/calibration/_exposure.py` | Add `weighted_gdna_eff_len_for_loci()`. |
| `src/rigel/calibration/_orchestrator.py` | Build `RegionalGdnaExposure` after `compute_global_densities()`. |
| `src/rigel/calibration/_result.py` | Add `regional_exposure` and application stats to `CalibrationResult` summary. |
| `src/rigel/calibration/locus_prior.py` | Thread exposure into `assemble_priors()`, compute weighted `L_g`, keep v6 prior, emit regional prior diagnostic. |
| `src/rigel/config.py` | Add flat `CalibrationConfig.regional_exposure_enabled: bool = True`. |
| `src/rigel/cli.py` | Add `--regional-exposure {auto,off}` and config transform. |
| `src/rigel/native/scoring.cpp` | Append per-unit `genomic_midpoint: int64[N_units]` to `StreamingScorer.finish()`. |
| `src/rigel/scan.py` | Unpack `genomic_midpoint` and populate `ScoredFragments`. |
| `src/rigel/scored_fragments.py` | Add `genomic_midpoint` to global `ScoredFragments` only. |
| `src/rigel/pipeline.py` | Apply per-unit `log A` after prior assembly and before partitioning; pass new diagnostics into locus output. |
| `src/rigel/estimator.py` | Add locus output columns for regional diagnostics. |
| `tests/test_regional_exposure.py` | New exposure-builder and lookup tests. |
| `tests/test_weighted_eff_len.py` | New weighted effective-length tests. |
| Existing routing/prior tests | Extend for no-op parity and per-unit numerator behavior. |
| `scripts/debug/gdna_regional_exposure_diag.py` | Pre-merge diagnostic script for VCaP hotspots. |

After editing `src/rigel/native/scoring.cpp`, recompile:

```bash
conda activate rigel && pip install --no-build-isolation -e .
```

---

## 4. Regional Exposure Module

### 4.1 Public Dataclasses

Create `src/rigel/calibration/_regional_exposure.py`:

```python
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

ExposureMode = Literal["uniform", "regional"]

REFERENCE_QUANTILE = 0.95
LOG_A_FLOOR = float(np.log(1.0e-4))
LOG_RHO_FLOOR = float(np.log(np.finfo(np.float64).tiny))
LOG_RHO_CLIP_NATS = 8.0   # spread metric: clip log rho_hat at log(rho_global_k) - 8
SPREAD_EPS = 1.0e-12
Z_Q95 = 1.6448536269514722  # standard normal 95th quantile


@dataclass(frozen=True, slots=True)
class RegionalWeightApplicationStats:
    n_units_seen: int = 0
    n_units_weighted: int = 0
    n_units_no_gdna: int = 0
    n_units_missing_midpoint: int = 0
    n_units_cross_ref_skipped: int = 0


@dataclass(frozen=True, slots=True)
class RegionalGdnaExposure:
    rho_hat: np.ndarray       # float64[R], aligned to RegionArrays sorted order
    log_weight: np.ndarray    # float64[R], log A_r in [LOG_A_FLOOR, 0]
    weight: np.ndarray        # float64[R], exp(log_weight)
    mode: ExposureMode
    rho_ref: float
    n_at_floor: int
    per_class: dict[str, dict[str, float]] = field(default_factory=dict)

    # Region geometry copied from RegionArrays so lookup does not need a
    # separate RegionIndexPy object.
    ref_offsets: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int32))
    ref_id: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int32))
    start: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    end: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
```

`RegionalGdnaExposure` owns the lookup table. Its arrays must be aligned to
`RegionArrays` sorted-position order, because `locus_prior.py` scratch arrays
also use that order.

### 4.2 Constructor API

```python
@classmethod
def uniform(cls, region_arrays: RegionArrays) -> RegionalGdnaExposure:
    """Return identity exposure, A_r == 1 and log A_r == 0."""

@classmethod
def build(
    cls,
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    global_densities: GlobalDensityTable,
    gdna_fl: FragmentLengthModel,
    *,
    strand_summary: StrandSummary | None = None,
    splicing_anchor_tolerance: int = 0,
    enabled: bool = True,
) -> RegionalGdnaExposure:
    """Build regional exposure or return uniform when disabled."""
```

If `enabled` is false, `build()` returns `uniform(region_arrays)`.

### 4.3 Lookup API

```python
def log_weights_for_positions(
    self,
    ref_ids: np.ndarray,      # int32[N]
    positions: np.ndarray,    # int64[N]
) -> np.ndarray:
    """Vectorized log A(ref, pos). Uniform mode returns zeros."""

def weighted_length_on_ref(self, ref_id: int, start: int, end: int) -> float:
    """Return integral of A(x) over [start, end) on ref_id."""
```

Lookup implementation uses one Python loop over references and
`np.searchsorted` within each ref's region slice. Positions outside any
region return `log A = 0` for per-unit lookup (the per-unit midpoint is
always a single point, and the region partition does not cover unmapped
contigs).

For `weighted_length_on_ref`, the region partition is gap-free by
construction (`regions.py` emits exactly one `INTERGENIC` between genic
spans). A gap inside `[start, end)` on a covered contig therefore signals a
corrupted `RegionArrays` and **must raise** `RuntimeError`; silently filling
the gap with identity weight would mask a partition bug. The one explicit
exception is when the entire query interval lies on a contig that has no
regions at all (e.g., reference present in the FASTA but absent from the
GTF); that case returns `end - start` (identity exposure) and is reported
once per contig in calibration diagnostics.

---

## 5. Exposure Mathematics

### 5.1 Per-Region Exposure `E_r`

For each region `r` in class `k`:

| Class | Exposure `E_r` |
|---|---|
| `INTERGENIC` | `l_eff_contained(end_r - start_r, gdna_fl)` |
| `INTRON` | `l_eff_contained(end_r - start_r, gdna_fl)` |
| `EXON` | `(eligible_left + eligible_right) * B_cross(gdna_fl, K)` |

`B_cross` is `boundary_crossing_exposure(gdna_fl,
splicing_anchor_tolerance=K)` and `K` must match the scanner payload.

### 5.2 Per-Region Conservative Count `Y_r`

Counts are computed in `RegionArrays` sorted order. The strand-correction
formula `total - (same - opp) / signed_strand_contrast` matches v6 exactly
(see `density_global._gdna_count_moment`). To avoid two sources of truth,
**extract and reuse**:

- `_strand_identifiable_rows(strands)` from `density_global.py` for the
  per-region identifiability mask.
- A new shared helper `_strand_correction_active(strand_summary)` that
  encapsulates the v6 gate: `abs(signed_strand_contrast) >=
  max(STRAND_CONTRAST_NUMERICAL_FLOOR, contrast_margin(0.99))`. This must
  be promoted out of `_compute_density()` so both the global density code
  and `_regional_exposure.py` import it.
- `_gdna_count_moment(n_same, n_opp, signed_strand_contrast=...)` for the
  bias-corrected count moment.

Intergenic:

```python
Y_ig = payload_arrays.intergenic_per_region.astype(np.float64)
```

Intergenic regions are strand-uninformative, so no strand correction is
applied.

Intron:

```python
raw = payload_arrays.intron_by_orient.sum(axis=1).astype(np.float64)
if _strand_correction_active(strand_summary):
    identifiable = _strand_identifiable_rows(region_arrays.strand)
    same = payload_arrays.intron_by_orient[:, ORIENT_SAME]
    opp  = payload_arrays.intron_by_orient[:, ORIENT_OPP]
    corrected = _gdna_count_moment(same, opp, signed_strand_contrast=ssc)
    Y_in = np.where(identifiable, np.maximum(corrected, 0.0), raw)
else:
    Y_in = raw
```

Exon boundary (eligibility identical to v6 `_channel_exon_intron`):

```python
eligible_left  = bf_left.astype(np.int64)
eligible_right = bf_right.astype(np.int64)
raw = (eligible_left * payload_arrays.u_left_by_orient.sum(axis=1) +
       eligible_right * payload_arrays.u_right_by_orient.sum(axis=1)).astype(np.float64)
if _strand_correction_active(strand_summary):
    identifiable = _strand_identifiable_rows(region_arrays.strand)
    same = (eligible_left * payload_arrays.u_left_by_orient[:, ORIENT_SAME] +
            eligible_right * payload_arrays.u_right_by_orient[:, ORIENT_SAME])
    opp  = (eligible_left * payload_arrays.u_left_by_orient[:, ORIENT_OPP] +
            eligible_right * payload_arrays.u_right_by_orient[:, ORIENT_OPP])
    corrected = _gdna_count_moment(same, opp, signed_strand_contrast=ssc)
    Y_ex = np.where(identifiable, np.maximum(corrected, 0.0), raw)
else:
    Y_ex = raw
```

`PayloadArrays` therefore needs the orientation-resolved arrays from
`CalibrationScanPayload`, and `RegionArrays` needs a per-region `strand`
field.

### 5.3 EB-Shrunk Density

For class `k`:

```text
rho_hat_r = (Y_r + kappa_k * rho_global_k) / (E_r + kappa_k)
```

When `E_r == 0`, this returns `rho_global_k` as long as `kappa_k > 0`, which
is the desired no-local-evidence behavior.

### 5.4 Exposure-Weighted Quantiles

Define `_weighted_quantile(values, weights, q)` precisely:

1. Keep rows where `values` is finite and `weights > 0`.
2. Sort by `values` ascending.
3. Compute cumulative weights.
4. Return the first value whose cumulative weight is at least `q * total_weight`.
5. If total weight is zero, return the class fallback value.

For each class:

```text
rho_ref_k = weighted_q95(rho_hat_r, weights=E_r)
```

If a class has zero total exposure, fall back to its `rho_global_k`. The global
reference is:

```text
rho_ref = max(rho_ref_intergenic, rho_ref_intron, rho_ref_exon)
```

If `rho_ref <= 0`, return uniform exposure.

### 5.5 Auto-Uniform Signal

For each class, compute observed log-density spread. We clip `log rho_hat`
at a *class-scaled* floor rather than `LOG_RHO_FLOOR ~ -708` so that
float64-tiny shrinkage residuals do not dominate the spread metric:

```text
log_rho_floor_k = log(rho_global_k) - LOG_RHO_CLIP_NATS    # LOG_RHO_CLIP_NATS = 8.0
log_rho_r = log(max(rho_hat_r, exp(log_rho_floor_k)))
observed_spread_k = weighted_q95(log_rho_r, E_r) - weighted_q50(log_rho_r, E_r)
```

`LOG_RHO_CLIP_NATS = 8.0` corresponds to a ~3000x dynamic range below the
class-global density, far wider than any real per-region shrinkage residual
but tight enough to prevent the q50 from collapsing onto a representation
floor. `LOG_RHO_FLOOR` is still used as a defensive last-resort cap inside
the vector before quantile sorting.

Compute a deterministic null spread using the delta-method approximation under
`Y_r ~ Poisson(rho_global_k * E_r)`:

```text
Var(log rho_hat_r) ~= E_r / (rho_global_k * (E_r + kappa_k)^2)
sigma2_pooled = weighted_mean(Var(log rho_hat_r), weights=E_r)
null_spread_k = Z_Q95 * sqrt(max(sigma2_pooled, 0))
```

If `rho_global_k <= 0` or the class has no positive exposure, set
`null_spread_k = 0` and `signal_k = 0`.

Then:

```text
signal_k = clamp((observed_spread_k - null_spread_k) / max(observed_spread_k, SPREAD_EPS), 0, 1)
```

There is no random fallback in v3. The deterministic approximation is emitted
as diagnostics and can be replaced later if benchmarks show it is inadequate.

### 5.6 Region Weight

For each region in class `k`:

```text
raw_log_ratio_r = min(log(rho_hat_r) - log(rho_ref), 0)
log_A_r = signal_k * raw_log_ratio_r
log_A_r = max(log_A_r, LOG_A_FLOOR)
A_r = exp(log_A_r)
```

If every `signal_k == 0`, set `mode = "uniform"` and return exact identity
arrays (`log_weight == 0`, `weight == 1`). Otherwise set `mode = "regional"`.

`LOG_A_FLOOR` is a numerical safety constant, not a user-facing parameter.
Emit `n_at_floor` and per-class `signal`, `observed_spread`, `null_spread`,
and weighted quantiles in `summary.json`.

---

## 6. Weighted gDNA Effective Length

Add to `src/rigel/calibration/_exposure.py`:

```python
def weighted_gdna_eff_len_for_loci(
    loci: tuple | list,
    ref_lengths: Mapping[str | int, int] | Sequence[int],
    fl: FragmentLengthModel,
    exposure: RegionalGdnaExposure,
    *,
    min_value: float = 1.0,
) -> float:
    ...
```

Use `TYPE_CHECKING` or local imports to avoid a circular import between
`_exposure.py` and `_regional_exposure.py`.

Uniform fast path:

```python
if exposure.mode == "uniform":
    return gdna_eff_len_for_loci(loci, ref_lengths, fl, min_value=min_value)
```

Regional path mirrors `gdna_eff_len_for_loci()` exactly:

1. Iterate over `positive_ell = np.flatnonzero(fl.pmf[1:] > 0) + 1`.
2. For each `ell`, build valid start windows `[a - ell + 1, b)` for every
   `Locus`, clipped to contig-valid starts.
3. Merge overlapping start windows per reference.
4. Shift each merged start window to a midpoint window using `ell // 2`.
5. Add `pmf[ell] * exposure.weighted_length_on_ref(ref_id, mlo, mhi)`.
6. Return `max(total, min_value)`.

When `A_r == 1`, the regional path is mathematically identical to the existing
unweighted function. The production uniform branch preserves bit-exactness.

---

## 7. Native Scoring Change

Append one element to `StreamingScorer.finish()`:

```text
genomic_midpoint: int64[N_units]
```

Append it at the end of the returned tuple so existing positional consumers of
the first 23 elements are not shifted.

Implementation in `scoring.cpp`:

- Add `std::vector<int64_t>* v_gmid` to `FillState` / `StreamingScorer`.
- Non-multimap EM units:
  - If the unit has a finite gDNA hypothesis and `g_sta >= 0` and `g_fp > 0`,
    push `int64_t(g_sta) + int64_t(g_fp) / 2`.
  - Otherwise push `INT64_MIN`.
- Multimapper EM units:
  - Track the first valid unspliced member midpoint in the MM group.
  - If the group emits a finite gDNA likelihood, push that midpoint.
  - Otherwise push `INT64_MIN`.
- Deterministic-unambiguous and chimeric outputs are not EM units and do not
  get entries in `genomic_midpoint`.

Python changes:

- `scan.py` unpacks `result[-1]` as `genomic_midpoint`.
- `ScoredFragments` gains `genomic_midpoint: np.ndarray`.
- Do not add `genomic_midpoint` to `LocusPartition`.

---

## 8. Pipeline Order

Inside `quant_from_buffer`, the regional path must run in this order:

1. Score fragments into global `ScoredFragments`.
2. Build `multi_loci`.
3. Assign locus IDs.
4. Assemble priors using the regional exposure table. This computes weighted
   `gdna_eff_len`, but `enable_gdna` still sees the original finite/non-finite
   `gdna_log_liks`.
5. Apply per-unit `log A` in place on `em_data.gdna_log_liks`; capture the
   returned `RegionalWeightApplicationStats` and persist it onto
   `calibration_result` via
   `calibration_result.with_regional_weighting_stats(stats)` so it reaches
   `summary.json`.
6. Partition and free.
7. Run EM.

This preserves the existing invariant that `assemble_priors()` runs before
`partition_and_free()` nulls global arrays.

**Apply-order invariant.** `assemble_priors()` calls
`enable_gdna_for_multilocus()`, which routes on `np.isfinite(gdna_log_liks)`.
The per-unit weight added in step 5 is always finite (clipped to
`LOG_A_FLOOR ~ -9.21`), so applying weights after `assemble_priors()` is
behaviorally identical to applying them before for the purpose of
`enable_gdna` eligibility. We apply after purely for code clarity (single
ownership of the global `gdna_log_liks` array).

### 8.1 Per-Unit Weight Function

Add to `pipeline.py`:

```python
def _apply_unit_gdna_weights(
    em_data: ScoredFragments,
    exposure: RegionalGdnaExposure,
    index: TranscriptIndex,
) -> RegionalWeightApplicationStats:
    ...
```

Algorithm:

```python
if exposure.mode == "uniform":
    return stats with n_units_seen=em_data.n_units

# Invariant guard: every EM unit must have >=1 candidate. np.*.reduceat
# silently misbehaves on empty groups, so assert before reducing.
assert (np.diff(em_data.offsets) > 0).all(), "empty EM unit detected"

finite = np.isfinite(em_data.gdna_log_liks)
has_mid = em_data.genomic_midpoint != np.iinfo(np.int64).min
base_mask = finite & has_mid

candidate_refs = index.t_to_ref_arr[em_data.t_indices]
unit_starts = em_data.offsets[:-1]
unit_min_ref = np.minimum.reduceat(candidate_refs, unit_starts)
unit_max_ref = np.maximum.reduceat(candidate_refs, unit_starts)
same_ref = unit_min_ref == unit_max_ref

mask = base_mask & same_ref
log_w = exposure.log_weights_for_positions(unit_min_ref[mask], em_data.genomic_midpoint[mask])
em_data.gdna_log_liks[mask] += log_w.astype(em_data.gdna_log_liks.dtype, copy=False)
```

Units with finite gDNA likelihood but candidate transcripts on multiple refs are
left unchanged (`A=1`) and counted as `n_units_cross_ref_skipped`.

`index.t_to_ref_arr` does not currently exist. Add it in
`TranscriptIndex.load()` using the same categorical-to-canonical mapping
already reconstructed in `locus.py` and `locus_prior.py`.

---

## 9. Prior Assembly and Locus Diagnostics

### 9.1 `PriorTable`

Extend `PriorTable`:

```python
gdna_eff_len_unweighted: np.ndarray        # float64[n_loci]
gdna_prior_count_regional: np.ndarray      # float64[n_loci], diagnostic only
```

Update `PriorTable.empty()` accordingly.

### 9.2 `assemble_priors()`

Signature:

```python
def assemble_priors(
    ...,
    regional_exposure: RegionalGdnaExposure | None = None,
) -> PriorTable:
```

If `regional_exposure is None`, create `RegionalGdnaExposure.uniform(region_arrays)`.

For each `MultiLocus`:

```python
unweighted = gdna_eff_len_for_loci(ml.loci, ref_lengths_arr, gdna_fl, min_value=1.0)
weighted = weighted_gdna_eff_len_for_loci(
    ml.loci,
    ref_lengths_arr,
    gdna_fl,
    regional_exposure,
    min_value=1.0,
)
gdna_eff_len_arr[idx] = weighted
gdna_eff_len_unweighted_arr[idx] = unweighted
```

`gdna_prior_count_arr[idx]` remains the existing v6 `eta_g` from
`expected_gdna_count_global()`.

Compute `gdna_prior_count_regional` as a diagnostic expected count under the
regional density field:

```text
eta_regional(Locus) =
    sum_over_overlapping_regions rho_hat_r * eff_clip_core_r
  + sum_over_eligible_exon_sides rho_hat_r * B_cross
```

The first term covers contained fragments. The second mirrors the existing
boundary-crossing exposure term. The diagnostic is summed across constituent
`Locus` intervals of a `MultiLocus`.

### 9.3 `loci.feather` Data Path

`loci.feather` is built from `estimator.locus_results`, so diagnostics must
flow through `_run_locus_em_partitioned()`:

1. Add optional arrays to `_run_locus_em_partitioned()`:
   - `gdna_eff_len_unweighted`
   - `gdna_prior_count_regional`
2. Add fields to `_build_locus_meta()`:
   - `gdna_eff_len_unweighted`
   - `gdna_eff_len_weight_ratio`
   - `gdna_prior_count_regional`
3. Extend `AbundanceEstimator.get_loci_df()` with those columns.

Do not mutate `loci_df` in `cli.py`; the estimator owns the table schema.

---

## 10. Calibration Result and Summary

### 10.1 Building the Exposure Table

In `_orchestrator.calibrate()` after `compute_global_densities()`:

```python
region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
regional_exposure = RegionalGdnaExposure.build(
    region_arrays,
    payload_arrays,
    global_densities,
    fl_models.gdna,
    strand_summary=strand_summary,
    splicing_anchor_tolerance=int(getattr(payload, "splicing_anchor_tolerance", 0)),
    enabled=regional_exposure_enabled,
)
```

`calibrate()` gains a keyword:

```python
regional_exposure_enabled: bool = True
```

`pipeline.run_pipeline()` passes `config.calibration.regional_exposure_enabled`.

### 10.2 `CalibrationResult`

Add fields:

```python
regional_exposure: RegionalGdnaExposure | None = None
regional_weighting_stats: RegionalWeightApplicationStats | None = None
```

Production `calibrate()` always supplies `regional_exposure`. Direct unit-test
construction may leave it as `None`; `quant_from_buffer()` should defensively
fall back to uniform exposure in that case.

`with_priors()` must preserve both fields via `dataclasses.replace()`.

Add:

```python
def with_regional_weighting_stats(self, stats: RegionalWeightApplicationStats) -> CalibrationResult:
    return dataclasses.replace(self, regional_weighting_stats=stats)
```

`to_summary_dict()` emits:

```json
"regional_exposure": {
  "mode": "regional",
  "rho_ref": 1.81e-4,
  "n_regions": 296039,
  "n_at_floor": 142,
  "log_a_floor": -9.21034,
  "per_class": {
    "INTERGENIC": {
      "rho_q05": 4.0e-6,
      "rho_q50": 7.9e-6,
      "rho_q95": 1.2e-5,
      "rho_ref_class": 1.2e-5,
      "observed_log_spread": 0.42,
      "null_log_spread": 0.39,
      "signal": 0.07
    }
  },
  "application": {
    "n_units_seen": 1000000,
    "n_units_weighted": 980000,
    "n_units_no_gdna": 10000,
    "n_units_missing_midpoint": 1000,
    "n_units_cross_ref_skipped": 9000
  }
}
```

---

## 11. Config and CLI

### 11.1 Config

In `CalibrationConfig`:

```python
regional_exposure_enabled: bool = True
```

No nested config dataclass in v3.

### 11.2 CLI

Add:

```text
--regional-exposure {auto,off}
```

Mapping:

```text
auto -> calibration.regional_exposure_enabled = True
off  -> calibration.regional_exposure_enabled = False
```

Add a `_ParamSpec` transform for this enum so config files and
`summary.json["command"]` record the user-facing value.

---

## 12. Test Plan

### 12.1 Regional Exposure Tests

New `tests/test_regional_exposure.py`:

| Test | Assertion |
|---|---|
| Uniform constructor | `mode == "uniform"`, `weight == 1`, `log_weight == 0`. |
| Deterministic equal-density input | `signal == 0`, mode uniform, bit-exact identity weights. |
| Bimodal density input | High-density regions have weight near 1; low-density regions have lower weights; no public parameters needed. |
| Zero-exposure exon | `rho_hat == rho_global_exon`; weight finite. |
| Strand-corrected intron count | Per-region formula matches the existing global moment on a small synthetic table. |
| Vector lookup | `log_weights_for_positions()` matches scalar `RegionIndexPy.overlap()` for random positions. |
| Outside/sentinel positions | Return identity log weight `0`. |
| Floor clipping | `n_at_floor` increments and `log_weight == LOG_A_FLOOR`. |

### 12.2 Weighted Effective Length Tests

New `tests/test_weighted_eff_len.py`:

| Test | Assertion |
|---|---|
| Uniform fast path | Bit-exact equality with `gdna_eff_len_for_loci()`. |
| Fixed-weight partition invariance | Splitting a same-weight region leaves the result unchanged. |
| Two-state geometry | Far from contig edges, 10 kb at `A=1` plus 90 kb at `A=0.01` gives approximately `10_900` for long fixed FL. |
| Overlapping `Locus` intervals | Matches the merged-window expectation; no double counting. |
| Floor-weight geometry | All regions at `exp(LOG_A_FLOOR)` produce a tiny positive length, clamped only if below `min_value`. |

### 12.3 Prior Assembly Tests

Extend `tests/test_assemble_priors.py`:

| Test | Assertion |
|---|---|
| `regional_exposure=None` | Existing v6 `gdna_eff_len`, `gdna_prior_count`, and `enable_gdna` are bit-exact. |
| Uniform exposure | Same bit-exact parity as no exposure. |
| Regional exposure | `gdna_eff_len <= gdna_eff_len_unweighted`; `gdna_prior_count` unchanged; regional diagnostic finite. |
| Locus output diagnostics | New arrays flow through to `loci.feather` columns. |

### 12.4 Per-Unit Numerator Tests

Add tests near routing/pipeline tests:

| Test | Assertion |
|---|---|
| Uniform exposure | `_apply_unit_gdna_weights()` leaves `gdna_log_liks` unchanged. |
| Regional exposure | Finite same-ref units receive exactly `log A(ref, midpoint)`. |
| `-inf` gDNA units | Unchanged and counted as no-gDNA. |
| Missing midpoint | Unchanged and counted. |
| Cross-reference unit | Unchanged and counted as skipped. |

### 12.5 Native Tests

Extend scanner/scoring tests:

- `StreamingScorer.finish()` returns one appended `genomic_midpoint` array.
- Non-spliced finite-gDNA units have midpoint `g_sta + g_fp // 2`.
- Spliced/no-gDNA units have `INT64_MIN`.
- Multimapper groups use the first valid unspliced member midpoint.

### 12.6 Regression and Benchmarks

Required before merge:

```bash
conda activate rigel
pytest tests/ -v
ruff check src/ tests/ scripts/debug/gdna_regional_exposure_diag.py
```

Benchmark gates:

- Synthetic uniform sweeps: `--regional-exposure auto` must keep mRNA/nRNA/gDNA
  relative-error drift <= 0.5% versus current golden outputs.
- VCaP no-mm: gDNA-source -> RNA false positives should drop by at least 70%
  versus the current baseline, or the diagnostic report must explain why the
  predicted exposure shift did not materialize.

---

## 13. Stage-0 Diagnostic Script

Add `scripts/debug/gdna_regional_exposure_diag.py`.

Inputs:

- BAM path.
- Index directory.
- Optional BED/TSV of hotspot intervals. Default to the five VCaP no-mm
  hotspots from the prior investigation.

The script reruns `scan_and_buffer()` and `calibrate()` to obtain the payload;
saved `summary.json` is insufficient because it does not contain per-region
counts or boundary flux.

Outputs under `results/vcap_regional_exposure_diag_<date>/`:

- `regional_exposure_summary.tsv`: class quantiles, `rho_ref`, signal, floor hits.
- `hotspot_eff_len.tsv`: unweighted and weighted `L_g`, ratio, mean/median weight.
- `hotspot_regions.tsv`: top and bottom contributing regions by weighted exposure.
- `fragment_examples.tsv`: sampled false-RNA gDNA-source fragments with midpoint,
  region weight, old/new gDNA log likelihood, and predicted posterior shift.
- Markdown report in `docs/benchmarks/`.

Greenlight before enabling the production patch by default:

- Mean predicted combined gDNA log-posterior gain on the five VCaP hotspots is
  at least 2.0 nats.
- Uniform synthetic diagnostics produce mode `uniform` or negligible weights
  consistent with the <= 0.5% benchmark gate.

---

## 14. Rollback

`--regional-exposure off` sets `CalibrationConfig.regional_exposure_enabled = False`,
which builds `RegionalGdnaExposure.uniform(region_arrays)`. The uniform fast paths
then restore v6 behavior:

- Per-unit numerator adds zero.
- Weighted `L_g` delegates to `gdna_eff_len_for_loci()` bit-exactly.
- `gdna_prior_count` is unchanged.

The only native ABI change is the appended `genomic_midpoint` return value from
`StreamingScorer.finish()`. Reverting the C++ and Python unpacking changes fully
restores the pre-v3 binary interface.

---

## 15. Deferred to Later Versions

- Per-hit gDNA exposure weighting inside the C++ multimapper LSE.
- Replacing the EM `gdna_prior_count` with a regional expected-count prior.
- Applying exposure weights to mRNA/nRNA effective lengths.
- Footprint-average per-fragment weights.
- External BED target intervals.
- Random or exact Poisson null simulation for the auto-uniform signal.

---

## 16. Implementation Checklist

1. Add `t_to_ref_arr` to `TranscriptIndex.load()`.
2. Extend `RegionArrays` and `PayloadArrays` with strand/orientation fields.
3. Implement `_regional_exposure.py` and tests.
4. Implement `weighted_gdna_eff_len_for_loci()` and tests.
5. Build exposure in `calibrate()` and store on `CalibrationResult`.
6. Extend `PriorTable` and `assemble_priors()` diagnostics.
7. Add native `genomic_midpoint` and recompile.
8. Thread midpoint into global `ScoredFragments`.
9. Apply per-unit weights after prior assembly and before partitioning.
10. Extend locus output diagnostics.
11. Add CLI/config plumbing.
12. Run tests, synthetic benchmarks, and VCaP no-mm diagnostic.