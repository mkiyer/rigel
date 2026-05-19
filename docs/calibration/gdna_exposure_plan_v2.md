# gDNA Regional Exposure — Implementation Plan v2

**Status**: implementation-ready, supersedes `gdna_exposure_plan_v1.md`.
**Date**: 2026-05-18
**Companion**: `gdna_exposure_plan_v1_review.md` (the review that drove this rewrite).

---

## 1. One-paragraph summary

Learn a per-region gDNA-exposure weight `A_r ∈ (0, 1]` from the conservative
gDNA evidence v6 calibration already collects, then apply it to **both** the
per-fragment gDNA log-likelihood (`+ log A(midpoint)`) and the locus-wide
gDNA effective length used as the EM denominator. Ship denominator and
numerator together in one production patch. Keep the EM prior unchanged in
v2; only emit a diagnostic alternative. Add no user-facing knobs except a
single `--regional-exposure {auto,off}` flag. Estimate every threshold from
the data using exposure-weighted statistics and EB shrinkage that already
exist in `calibration/`.

---

## 2. What changed from v1

The review made seven structural corrections; v2 adopts every one and
removes the multi-stage rollout that motivated several of them.

1. **One shipping unit**, not three stages. The denominator term and the
   per-unit numerator term are mass-conserving siblings; shipping the
   denominator alone is wrong, and the v1 plan said so on one page and
   contradicted itself on another. v2 ships both together; "Stage 0" is a
   diagnostic that has no production code path.
2. **`ref_id` is required for per-unit lookup.** The native ABI gains one
   per-unit `genomic_midpoint` field; the unit's reference is inferred in
   Python from its candidate transcripts (`em_data.locus_t_indices →
   index.t_to_ref`). No new column on `FragmentBuffer`/`ResolvedFragment`.
3. **Weighted `L_g` weights midpoints of expanded start windows**, not
   regions clipped to `[a, b)`. Equivalent re-derivation of
   `gdna_eff_len_for_loci`'s per-`ell` window logic with one extra
   intersection step.
4. **EB split-invariance is not a real invariant** and is dropped from the
   test plan. The geometric invariance that *is* real — fixed weights,
   geometry unchanged under sub-division of a region with the same
   `A` — stays.
5. **`gdna_prior_count` is unchanged in v2.** Mixing global and per-region
   densities in `ρ_global · A_r` is not a clean expected-count expression.
   We emit a diagnostic alternative `Σ_r ρ̂_r · L_r` and revisit the
   prior in v3 once VCaP and synthetic results agree.
6. **Parameter removal.** `c_max`, `min_exposure_bp`,
   `min_class_regions_for_quantile`, `min_boundary_events_for_exon_refine`,
   `capture_index_threshold` are all gone from the public surface. Where a
   threshold is still mathematically necessary, it is learned from the
   data or it is a numerical safety floor that is reported, not tuned.
7. **Lookup is vectorized** (`np.searchsorted` per ref) — not a Python
   loop over units.

A second-order change: the auto-uniform mode uses a continuous attenuation
(`signal_k ∈ [0, 1]` derived from observed vs. null log-density spread)
instead of a hard threshold. Uniform libraries get `signal ≈ 0`, capture
libraries get `signal ≈ 1`, no knob needed.

---

## 3. Locked mathematics

Notation:
- Region `r` belongs to class `k ∈ {ig, in, ex}` (intergenic, intron, exon).
- `Y_r`, `E_r`: per-region gDNA-proxy count and FL-weighted exposure
  (defined in §3.1).
- `ρ_global_k`, `κ_k`: existing fields on `GlobalGdnaDensity`.
- `gdna_fl`: `FragmentLengthModel` already produced by calibration.

### 3.1 Per-region count and exposure (no new estimator)

These are the v6 calibration numerators/denominators, applied per region
instead of summed per class.

| Class | `Y_r` | `E_r` |
|---|---|---|
| `INTERGENIC` | `payload.per_region_counts[r, MASK_INTERGENIC]` (strand-corrected via the existing `compute_global_densities` correction when active) | `l_eff_contained(end_r − start_r, gdna_fl)` |
| `INTRON` | analogous, with the same strand correction | analogous |
| `EXON` | `(u_left[r] + u_right[r])` restricted to eligible sides (`bf_left`/`bf_right`) | `n_eligible_sides_r × B_cross(gdna_fl)` |

For `EXON` regions with `n_eligible_sides_r == 0`, we set `E_r = 0` and
the shrinkage in §3.2 returns `ρ_global_ex` (the prior wins). No
`min_boundary_events` threshold.

### 3.2 Per-region shrunk density

Use the existing `density_loco.shrink_to_loco` formula on a per-region
basis:

$$
\hat{\rho}_r \;=\; \frac{Y_r + \kappa_k \cdot \rho^{\text{global}}_k}{E_r + \kappa_k}
$$

Limits:
- `E_r = 0` ⇒ `ρ̂_r = ρ_global_k` (no local evidence).
- `Y_r = 0`, `E_r → ∞` ⇒ `ρ̂_r → 0` (genuinely empty region).
- `κ_k → 0` ⇒ pure local; `κ_k → ∞` ⇒ pure global.

### 3.3 Reference density (exposure-weighted, no quantile threshold)

For each class `k`, the per-class reference is the exposure-weighted Q95
of `ρ̂_r`:

$$
\rho^{\text{ref}}_k \;=\; Q_{0.95}^{w=E_r}\bigl(\hat{\rho}_r : r \in k\bigr)
$$

Tiny regions contribute proportionally to their exposure rather than being
excluded by a hard `min_exposure_bp` threshold. If a class has zero total
exposure (structurally impossible in v6 — there is always at least one
intergenic region — but defensively), we fall back to `ρ_global_k`.

The global reference is

$$
\rho^{\text{ref}} \;=\; \max_k \rho^{\text{ref}}_k
$$

### 3.4 Auto-uniform attenuation (no threshold)

For each class compute observed and null log-density spreads:

- `observed_spread_k = Q95_w(log ρ̂_r) − Q50_w(log ρ̂_r)`
- `null_spread_k = expected value of the same quantity if Y_r ~ Poisson(ρ_global_k · E_r) for every region`

`null_spread_k` is computed analytically from the per-region exposures
under the Poisson + EB-shrinkage null (no fitting, no simulation needed
in the common case; see Appendix A for the closed-form approximation,
otherwise we draw 256 Poisson samples per region — still parameter-free).

The per-class signal weight is:

$$
\text{signal}_k \;=\; \max\!\left(0,\; \frac{\text{observed\_spread}_k - \text{null\_spread}_k}{\max(\text{observed\_spread}_k,\; \varepsilon)}\right)
$$

A region's log-weight is attenuated by its class signal:

$$
\log A_r \;=\; \text{signal}_k(r) \cdot \min\!\left(\log \hat{\rho}_r - \log \rho^{\text{ref}},\; 0\right)
$$

so `A_r ∈ (0, 1]`. A numerical log-floor `log A_r ≥ LOG_A_FLOOR` (default
`log(1e-4) ≈ −9.2`) prevents `-inf` posterior catastrophes; this floor is
not a tuning knob — it is reported in `summary.json` and triggered counts
are emitted so we can monitor whether real data is hitting it.

**Properties**:
- Uniform library: `observed_spread_k ≈ null_spread_k` ⇒ `signal_k ≈ 0`
  ⇒ `log A_r ≈ 0` ⇒ `A_r ≈ 1`. Bit-exact `A_r = 1` is achieved by branching
  on `signal_k == 0` (we test the equality before exponentiation).
- Capture library: `observed_spread_k ≫ null_spread_k` ⇒ `signal_k ≈ 1` ⇒
  `A_r ≈ ρ̂_r / ρ^ref`.
- No `capture_index_threshold` knob.

### 3.5 Weighted gDNA effective length for a `MultiLocus`

For a `MultiLocus` `M = ⋃_i [a_i, b_i)` on (per-`i`) refs:

$$
\tilde{L}_g(M) \;=\; \sum_{\ell : h_g(\ell) > 0} h_g(\ell) \;\cdot\; \sum_{\text{ref } r} \int_{\text{starts}(\ell, r)} A_r\!\bigl(s + \lfloor \ell / 2 \rfloor\bigr) \, ds
$$

Algorithm — **mirrors `gdna_eff_len_for_loci` exactly, with one extra
intersection step per `ell`**:

```
for ell in positive_ell:                   # same iteration as the existing impl
    for ref in refs_in_locus:
        windows_start = []                 # build the merged start windows
        for [a, b) in loci on this ref:
            lo = max(a - ell + 1, 0)
            hi = min(b, ref_len - ell + 1)
            if hi > lo: windows_start.append((lo, hi))
        windows_start = merge(windows_start)        # exactly as today
        # Shift each start window to a midpoint window:
        shift = ell // 2
        midpoint_windows = [(lo + shift, hi + shift) for lo, hi in windows_start]
        # New step: intersect with region partition on this ref and weight
        for (mlo, mhi) in midpoint_windows:
            for (region_idx, overlap_bp) in regions_overlapping([mlo, mhi)):
                total += pmf[ell] * overlap_bp * A_r[region_idx]
```

The midpoint windows are a 1-bp-precise shift of the existing merged
start windows, so when `A_r ≡ 1` the per-`ell` sum
`Σ overlap_bp = (hi − lo)` reproduces `gdna_eff_len_for_loci` exactly.

**Uniform fast path**: when `exposure.mode == "uniform"`, we return
`gdna_eff_len_for_loci(...)` directly — bit-exact no-op, no new arithmetic.

### 3.6 Per-unit gDNA log-likelihood update

For each unit `u` with finite `gdna_log_liks[u]`:

$$
\text{gdna\_log\_liks}[u] \mathrel{+}= \max\bigl(\log A(\text{ref}_u,\; \text{midpoint}_u),\; \texttt{LOG\_A\_FLOOR}\bigr)
$$

- `midpoint_u`: from the new per-unit `genomic_midpoint` field (§4.4).
- `ref_u`: inferred in Python from `em_data.locus_t_indices[u]` via the
  index's `t_to_ref` table (the first candidate's reference; all
  candidates of one unit live in the same locus, hence — in practice —
  the same ref).
- Units with `gdna_log_liks[u] = −∞` (no gDNA hypothesis) and units with
  `midpoint_u == INT64_MIN` (sentinel) are skipped.

Lookups are vectorized: group unit indices by `ref_u`, then use
`np.searchsorted(region_starts_for_ref, midpoint_u)` to find the
containing region in O(log R) per unit with one Python-side loop over
references.

### 3.7 EM prior: unchanged in v2

`gdna_prior_count` keeps the existing v6 `η_g` formula. The diagnostic
alternative `η_g^reg = Σ_{r ∈ Locus} ρ̂_r · L_core_r` is computed and
emitted to `loci.feather` as `gdna_prior_count_regional` but not fed to
the EM. The rationale is that the well-defined regional expected-count
formula uses `ρ̂_r` directly, not `ρ_global · A_r` (the latter is not a
proper density), and we want the VCaP/synthetic benchmark to isolate the
likelihood-side effect before touching the prior.

---

## 4. Code integration

### 4.1 Files touched

| Path | Change |
|---|---|
| `src/rigel/calibration/_regional_exposure.py` | new module — table builder + vectorized lookup |
| `src/rigel/calibration/_exposure.py` | add `weighted_gdna_eff_len_for_loci` with exact uniform branch |
| `src/rigel/calibration/_orchestrator.py` | call the builder; pass result to `CalibrationResult` |
| `src/rigel/calibration/_result.py` | add `regional_exposure` field + summary dict block |
| `src/rigel/calibration/locus_prior.py` | route `regional_exposure` to `assemble_priors`; populate `gdna_eff_len` via the weighted helper; emit diagnostic `gdna_prior_count_regional` (not used by EM) |
| `src/rigel/config.py` | `RegionalExposureConfig(enabled: bool = True)` only |
| `src/rigel/cli.py` | one flag: `--regional-exposure {auto,off}` |
| `src/rigel/native/scoring.cpp` | emit per-unit `genomic_midpoint: int64[N_units]` (sentinel `INT64_MIN`) |
| `src/rigel/scan.py` | thread the new array out of `StreamingScorer.finish()` |
| `src/rigel/scored_fragments.py` | add `genomic_midpoint: np.ndarray` field |
| `src/rigel/pipeline.py` | call `_apply_unit_gdna_weights` after scoring, before partitioning; emit `loci.feather` diagnostics |
| `tests/test_regional_exposure.py` | new — see §6 |
| `tests/test_weighted_eff_len.py` | new — see §6 |
| `tests/test_assemble_priors.py` | extend with `regional_exposure=None`/`uniform` parity tests |
| `scripts/debug/gdna_regional_exposure_diag.py` | Stage-0 diagnostic |

`em_solver.cpp` is unchanged; the weighted `L_g` value goes into the
existing `locus_gdna_eff_lens` parameter.

### 4.2 `_regional_exposure.py` — public surface

```python
"""Per-region gDNA exposure for capture-aware likelihood and effective length.

Single owner of A_r. Reuses v6 EB shrinkage from density_loco.
"""

from dataclasses import dataclass
import numpy as np
from typing import Literal
from ._arrays import RegionArrays, PayloadArrays
from .density_global import GlobalDensityTable
from ._exposure import boundary_crossing_exposure
from ..frag_length_model import FragmentLengthModel

ExposureMode = Literal["uniform", "regional"]

# Numerical safety floor on log A — emitted, not tuned. Region weights
# below exp(LOG_A_FLOOR) ≈ 1e-4 are clipped here so per-fragment scores
# never become -inf from this term alone.
LOG_A_FLOOR: float = -9.21


@dataclass(frozen=True)
class RegionalGdnaExposure:
    """Per-region exposure aligned to RegionArrays sort order.

    All arrays have length R == number of regions. Shape and length
    are validated at construction.
    """

    rho_hat: np.ndarray              # float64, (R,)
    log_weight: np.ndarray           # float64, (R,) — log A_r
    weight: np.ndarray               # float64, (R,) — A_r
    mode: ExposureMode               # "uniform" → weight ≡ 1, log_weight ≡ 0
    rho_ref: float                   # max over classes of per-class Q95
    n_at_floor: int                  # count of log_weight == LOG_A_FLOOR
    per_class: dict[str, dict[str, float]]  # for summary.json

    # Sorted-per-ref CSR for fast searchsorted lookup
    ref_offsets: np.ndarray          # int32, (n_refs + 1,)
    region_starts_by_ref: np.ndarray # int64, (R,) — region.start, grouped per ref
    region_ends_by_ref: np.ndarray   # int64, (R,)
    sorted_region_idx: np.ndarray    # int64, (R,) — back to RegionArrays index

    @classmethod
    def build(
        cls,
        region_arrays: RegionArrays,
        payload_arrays: PayloadArrays,
        global_densities: GlobalDensityTable,
        gdna_fl: FragmentLengthModel,
        *,
        splicing_anchor_tolerance: int = 0,
    ) -> "RegionalGdnaExposure":
        """Build the per-region exposure from existing v6 inputs."""
        ...

    @classmethod
    def uniform(cls, region_arrays: RegionArrays) -> "RegionalGdnaExposure":
        """Identity exposure: weight ≡ 1, log_weight ≡ 0, mode='uniform'."""
        ...

    def weights_for_positions(
        self,
        ref_ids: np.ndarray,    # int32, (N,)
        positions: np.ndarray,  # int64, (N,)
    ) -> np.ndarray:
        """Vectorized A(ref, pos) lookup. Returns float64 (N,) in (0, 1].

        Positions outside any region return 1.0 (treated as no
        information; identity weight). Sentinel positions (INT64_MIN)
        also return 1.0; callers must mask separately if they need
        them flagged.
        """
        ...

    def log_weights_for_positions(
        self,
        ref_ids: np.ndarray,
        positions: np.ndarray,
    ) -> np.ndarray:
        """Same as weights_for_positions but returns log A (already clipped
        to LOG_A_FLOOR). Uniform mode returns zeros without a lookup."""
        ...
```

Implementation notes:

- `build` computes per-region `Y_r`, `E_r` (§3.1), shrinks (§3.2),
  computes exposure-weighted Q95 per class (§3.3), runs the
  signal/null calculation (§3.4), and emits the table. All arrays are
  aligned to `RegionArrays` sort order so a `RegionIndexPy.overlap`
  result indexes into them directly.
- `uniform` is a structural no-op; production code must call this when
  `enabled=False` so the rest of the pipeline can rely on the same
  table type.
- `weights_for_positions` uses one `np.searchsorted` per reference
  (small Python loop bounded by the number of references — typically
  25 for hg38). For each ref `r`: `idx = np.searchsorted(starts_r,
  positions_r, side='right') - 1`; positions where `idx < 0` or
  `positions_r >= ends_r[idx]` fall outside any region.

### 4.3 `weighted_gdna_eff_len_for_loci` in `_exposure.py`

```python
def weighted_gdna_eff_len_for_loci(
    loci: tuple | list,
    ref_lengths,
    fl: FragmentLengthModel,
    exposure: RegionalGdnaExposure,
    *,
    min_value: float = 1.0,
) -> float:
    """Per-MultiLocus weighted gDNA overlap effective length.

    Uniform fast path: when exposure.mode == 'uniform', delegates to
    gdna_eff_len_for_loci for bit-exact backward compatibility.

    Algorithm (§3.5): mirrors gdna_eff_len_for_loci's per-ell merged
    start windows, then shifts to midpoint windows and intersects with
    the region partition, weighting each overlap_bp by A_r.
    """
    if exposure.mode == "uniform":
        return gdna_eff_len_for_loci(loci, ref_lengths, fl, min_value=min_value)
    ...
```

The intersection-with-region-partition step uses the per-ref CSR of
region starts stored on `RegionalGdnaExposure`; one `np.searchsorted`
locates the first region overlapping the midpoint window, then we
walk forward.

### 4.4 Native ABI change: per-unit `genomic_midpoint`

`StreamingScorer.finish()` in `scoring.cpp` currently returns a tuple of
23 arrays/scalars (see `scan.py` line ~145). v2 appends one element:

- `genomic_midpoint: int64_t[N_units]`. For each unit:
  - **Unambiguous / strand-ambig units (NH = 1)**: `g_sta + g_fp / 2`
    where `g_sta` and `g_fp` are the read's genomic start and gDNA
    footprint (already tracked in `score_mm_alignment` and the non-MM
    path).
  - **Multimapper units (NH > 1)**: the midpoint of the *primary*
    alignment (the first member in the qname group). Per-hit
    weighting is a future-stage refinement.
  - **Units without a gDNA hypothesis** (spliced-only, chimeric):
    sentinel `INT64_MIN`.

Python adds `genomic_midpoint: np.ndarray` (`int64`) to
`ScoredFragments`. `locus_partition.partition_and_free` scatters it to
per-locus partitions alongside `gdna_log_liks`.

Recompile required (`pip install --no-build-isolation -e .`).

### 4.5 Pipeline integration

`pipeline.quant_from_buffer` change, preserving the
"assemble_priors before partition_and_free" invariant:

```python
scored_fragments = ...                                          # existing
# 1. Build / propagate exposure table
exposure = calibration_result.regional_exposure                 # always present

# 2. Apply per-unit numerator (in-place, no-op in uniform mode)
if exposure.mode == "regional":
    _apply_unit_gdna_weights(
        em_data=scored_fragments,
        exposure=exposure,
        index=index,
    )

# 3. Assemble priors (uses weighted L_g; prior count unchanged in v2)
prior_table = assemble_priors(
    multi_loci, scored_fragments, index, payload, global_densities,
    regional_exposure=exposure,
    ...
)

# 4. partition_and_free, run EM (unchanged paths)
```

`_apply_unit_gdna_weights` is a single vectorized function:

```python
def _apply_unit_gdna_weights(em_data, exposure, index):
    finite = np.isfinite(em_data.gdna_log_liks)
    has_mid = em_data.genomic_midpoint != np.iinfo(np.int64).min
    mask = finite & has_mid
    if not mask.any():
        return
    # Infer ref from first candidate of each unit:
    first_t = em_data.t_indices[em_data.offsets[:-1]]
    ref_ids = index.t_to_ref_arr[first_t]  # int32
    log_w = exposure.log_weights_for_positions(
        ref_ids[mask], em_data.genomic_midpoint[mask]
    )
    em_data.gdna_log_liks[mask] += log_w
```

Multimapper note: the `first_t`-by-`offsets` trick assumes every
unit's candidate transcripts share a reference, which is true for
locus-local fragments by construction. We add an assert in v2 that
checks this on a sample, gated by a debug flag.

### 4.6 Config and CLI

```python
# config.py
@dataclass(frozen=True)
class RegionalExposureConfig:
    enabled: bool = True
```

```text
# cli.py
--regional-exposure {auto,off}   # default 'auto'; 'off' → enabled=False
```

No other knobs. Q95 constant, `LOG_A_FLOOR`, and the null-spread Poisson
sample count are module-level constants in `_regional_exposure.py`,
visible to anyone reading the source, not surfaced as configuration.

### 4.7 `loci.feather` and `summary.json` additions

`loci.feather` adds three columns (existing columns unchanged):

- `gdna_eff_len_unweighted` — what v6 would have produced.
- `gdna_eff_len_weight_ratio` — `gdna_eff_len / gdna_eff_len_unweighted` ∈ (0, 1].
- `gdna_prior_count_regional` — diagnostic alternative prior (see §3.7),
  not consumed by EM.

`summary.json["calibration"]` gets:

```json
"regional_exposure": {
  "mode": "regional",
  "rho_ref": 1.81e-4,
  "n_regions": 296039,
  "n_at_floor": 142,
  "log_a_floor": -9.21,
  "per_class": {
    "INTERGENIC": {
      "exposure_weighted_rho_q05/q50/q95": [...],
      "observed_log_spread": 2.4,
      "null_log_spread": 0.3,
      "signal": 0.875,
      "rho_ref_class": 1.2e-5
    },
    "INTRON":  { ... },
    "EXON":    { ... }
  }
}
```

---

## 5. Stage-0 diagnostic (no production code change)

**Path**: `scripts/debug/gdna_regional_exposure_diag.py`.

**Inputs**: the index dir and the input BAM. The script **reruns**
`scan_and_buffer` + `calibrate(...)` to obtain `payload` (it cannot
reconstruct per-region counts from `summary.json` alone — the v1 plan was
wrong on this point).

**Outputs** under `results/vcap_regional_exposure_diag_<date>/`:
- Global: `capture_index` per class, `signal_k`, `ρ_ref_k`, `ρ_ref`,
  number of regions by weight decile.
- Per hotspot locus: `gdna_eff_len_unweighted`, `gdna_eff_len_weighted`,
  `weight_ratio`, top/bottom 5 contributing regions, mean/median `A`
  inside the locus.
- Per-fragment (sample 500 hotspot fragments tagged as gDNA-source in
  the oracle): current and predicted gDNA log-likelihood with and
  without the per-unit term; predicted posterior swing.

**Greenlight criterion** for merging the production patch:
On the five known VCaP hotspots, the predicted *combined*
(numerator + denominator) log-posterior gain for true gDNA fragments must
be ≥ 2.0 nats on average. If smaller, do not ship the patch as-is —
investigate the per-fragment derivation first.

---

## 6. Test plan

All tests are required to pass before merge. Numeric tolerances are
inline with each.

### 6.1 `tests/test_regional_exposure.py` (new)

| ID | Test | Assertion |
|---|---|---|
| E1 | `uniform(region_arrays)` returns identity | `mode == 'uniform'`, `weight ≡ 1.0`, `log_weight ≡ 0.0` |
| E2 | `build()` on synthetic-uniform Y_r (Poisson with class-uniform ρ) | `signal_k ≈ 0` for every class; `mode == 'uniform'`; `weight ≡ 1.0` bit-exact |
| E3 | `build()` on bimodal Y_r (10 high-density + 90 low-density regions, 100× ratio) | `signal_k ≈ 1`; `weight` ≈ {1.0, 0.01}; `n_at_floor == 0` |
| E4 | `build()` shrinkage on a tiny region with 0 counts | `ρ̂_r ≈ ρ_global_k`; `weight_r ≈ ρ_global / ρ_ref` |
| E5 | `weights_for_positions` correctness against per-call `RegionIndexPy.overlap` | exact match for 10k random (ref, pos) pairs |
| E6 | `weights_for_positions` on positions outside any region | returns 1.0 |
| E7 | `weights_for_positions` on sentinel `INT64_MIN` | returns 1.0 (caller is responsible for masking) |
| E8 | `log_weights_for_positions` uniform-mode short circuit | returns `np.zeros(N)` without calling `searchsorted` (verified via mock) |
| E9 | LOG_A_FLOOR clipping triggers when `ρ̂_r / ρ_ref < exp(LOG_A_FLOOR)` | `n_at_floor` increments; `log_weight == LOG_A_FLOOR` |

### 6.2 `tests/test_weighted_eff_len.py` (new)

| ID | Test | Assertion |
|---|---|---|
| W1 | **Uniform identity**: `exposure.mode == 'uniform'` ⇒ `weighted_gdna_eff_len_for_loci == gdna_eff_len_for_loci` | bit-exact for 20 randomly generated loci |
| W2 | **Geometric partition invariance** (fixed weights, no re-build): split a region of constant `A` into two adjacent sub-regions with the same `A` and rerun | result unchanged to float64 precision |
| W3 | **Two-state geometry**: 100 kb locus, first 10 kb at `A=1`, last 90 kb at `A=0.01`, far from contig ends and long FL | result ≈ `10_000 + 90_000 × 0.01 = 10_900` within `1%` (FL-boundary tail) |
| W4 | **Multi-Locus merge**: two overlapping `Locus` intervals on one ref | matches the result for the merged interval, identical `A` per region |
| W5 | **No regions overlap the midpoint window**: degenerate exposure (all `A_r = 0`) | result clamps to `min_value` |

### 6.3 `tests/test_assemble_priors.py` (extend existing)

| ID | Test | Assertion |
|---|---|---|
| AP1 | `regional_exposure=None` ⇒ behavior unchanged from v6 | bit-exact `gdna_eff_len`, `gdna_prior_count`, `enable_gdna` |
| AP2 | `regional_exposure=uniform(R)` ⇒ behavior unchanged from v6 | bit-exact |
| AP3 | `regional_exposure=build()` on bimodal synthetic | `gdna_eff_len` strictly less than unweighted; `gdna_prior_count` unchanged (v2 prior decision) |
| AP4 | `gdna_prior_count_regional` is populated and != `gdna_prior_count` on regional data | column exists and is finite |

### 6.4 Per-unit numerator (`tests/test_pipeline_routing.py` or new)

| ID | Test | Assertion |
|---|---|---|
| U1 | `_apply_unit_gdna_weights` no-op in uniform mode (or when not called) | `gdna_log_liks` unchanged |
| U2 | Bimodal synthetic: per-unit delta equals `log A(midpoint)` for each unit | exact for each unit |
| U3 | `-inf` units untouched | `gdna_log_liks[is_-inf]` unchanged |
| U4 | Sentinel `midpoint == INT64_MIN` units untouched | unchanged |

### 6.5 Integration / benchmark gates

- All `scripts/benchmark/configs/locus_simple_*.yaml`: relative-error
  drift on mRNA/nRNA/gDNA `≤ 0.5%` vs. golden, with `--regional-exposure
  auto` (auto should detect uniform on synthetic).
- VCaP no_mm: FP gDNA→RNA misclassification count drops by `≥ 70%` vs.
  current main baseline (the deep-dive math predicts a swing from 2.64M
  to under 800K).
- Existing 1074+ tests: all green.

---

## 7. Risk register and rollback

| Risk | Mitigation |
|---|---|
| Auto-detector under-attenuates a real capture sample | `signal_k` is continuous; even partial attenuation moves in the right direction. The diagnostic script reports `observed_spread − null_spread` per class so misdetection is visible. |
| Auto-detector over-attenuates a uniform sample | `signal_k ≈ 0` produces `log A ≈ 0` ⇒ `A ≈ 1`. Worst case is a vanishing per-fragment additive term and `L̃_g ≈ L_g` (within float64). Synthetic suite catches > 0.5% drift. |
| Multimapper `ref_u` inference takes the first candidate's ref but other candidates are on a different ref | In practice, units in a v6 `MultiLocus` are locus-local. We add a debug-gated assert that all candidates of a unit share a ref on a sampled subset. |
| `LOG_A_FLOOR` triggers too often (region weight collapse) | `n_at_floor` is emitted in `summary.json`; flag if `> 1%` of regions hit the floor. |
| Native ABI churn breaks `scan_payload`-only consumers | Only `StreamingScorer.finish()` tuple changes; existing positional access reads forward (the new element is appended at the end). |

**Rollback**: `--regional-exposure off` ⇒ `RegionalExposureConfig(enabled=False)`
⇒ `RegionalGdnaExposure.uniform(region_arrays)` ⇒ every code path
short-circuits to v6 behavior. The native ABI change is the only
irreversible part; reverting the C++ patch restores pre-v2 binaries.

---

## 8. Out of scope (deferred to v3)

- Replacing `gdna_prior_count` with a regional formula (requires the
  v2 VCaP+synthetic benchmark to land first so we can attribute
  changes correctly).
- Per-hit weighting of multimapper gDNA likelihoods inside the C++
  scorer.
- Applying `A_r` to mRNA / nRNA effective lengths.
- Footprint-average per-fragment weights (midpoint is the v2 default;
  the diagnostic script reports midpoint-vs-footprint divergence so we
  can decide if v3 needs the upgrade).
- BED-file target intervals.

---

## 9. Appendix A — null log-spread closed form

For a single class with regions indexed by `r`:

- Under the null `Y_r ~ Poisson(ρ_g · E_r)`, the EB-shrunk estimate is
  `ρ̂_r = (Y_r + κ · ρ_g) / (E_r + κ)`.
- The conditional mean is `ρ_g` and the variance is
  `Var(ρ̂_r) = ρ_g · E_r / (E_r + κ)²`.
- The exposure-weighted distribution of `log ρ̂_r` under the null is
  approximately normal in the high-`E` tail (where `E · ρ_g ≫ 1`) with
  variance `(ρ_g · E_r) / [(E_r + κ)² · ρ_g²] = 1 / [(E_r + κ)² / E_r · ρ_g]`
  (delta method).
- For exposure-weighted Q95 vs. Q50, the spread is `≈ z_0.95 · σ_pooled`
  where `σ_pooled² = (Σ_r E_r · Var(log ρ̂_r)) / Σ_r E_r`.

This is an `O(R)` closed-form computation. When the high-`E` regime is
not dominant (rare; `Σ_r 1[ρ_g E_r ≥ 1] / R < 0.5`) we fall back to
256 Poisson samples per region — still parameter-free, ~`O(R)` time.

---

## 10. Appendix B — concrete formula summary

```
# Per-region (k = class of r):
Y_r = conservative count (intergenic / intron / boundary by class; §3.1)
E_r = FL-weighted exposure (contained-l_eff / n_sides × B_cross; §3.1)
rho_hat_r = (Y_r + kappa_k * rho_global_k) / (E_r + kappa_k)

# Per-class reference + auto attenuation:
rho_ref_k = Q95(rho_hat_r, weights = E_r, in class k)
rho_ref   = max_k rho_ref_k
observed_spread_k = Q95_w(log rho_hat_r) - Q50_w(log rho_hat_r)
null_spread_k     = analytical or Poisson-resampled under H0 (Appendix A)
signal_k = max(0, observed_spread_k - null_spread_k) / max(observed_spread_k, eps)

# Per-region weight:
log_A_r = signal_k(r) * min(log(rho_hat_r) - log(rho_ref), 0.0)
log_A_r = max(log_A_r, LOG_A_FLOOR)
A_r     = exp(log_A_r)

# Weighted MultiLocus gDNA effective length (uniform fast path elided):
L_tilde_g(M) = sum_ell h_g(ell) * sum over midpoint windows of overlap_bp * A_r

# Per-unit gDNA log-likelihood:
gdna_log_liks[u] += log_A_r(ref_u, midpoint_u)    # vectorized via searchsorted

# Prior in v2: unchanged. Diagnostic only:
gdna_prior_count_regional = sum over r in locus of rho_hat_r * L_core_r
```

Seven lines, plus the existing v6 calibration plumbing. That is the
entire v2 model.
