# gDNA Regional Exposure — Implementation Plan v1

**Status**: implementation-ready plan, not yet implemented.
**Date**: 2026-05-18
**Supersedes**: `gdna_exposure_model_v3.md` (model concept), incorporates v3 decisions.

---

## 0. Critique and revisions applied to v3

The v3 concept document is correct in shape but ambiguous in three places that
matter for code. This plan locks the decisions; readers of v3 should treat this
document as the authoritative source going forward.

| v3 statement | Issue | Decision in v1 |
|---|---|---|
| "Exon-region gDNA density from boundary evidence" (§ Exon Regions) | Underspecified: a single boundary side's density is noisy, and boundary flux already feeds the EXON-INTRON channel. Need a precise formula. | Use the per-region **EXON-INTRON shrunk density** from the existing boundary channel as the *baseline* for exonic positions. Per-exon refinement is shrunk to this baseline using a separate, locally-evaluated exon-boundary count when ≥ `min_boundary_events_for_exon_refine`. See §3.4. |
| "Rho_ref = robust high quantile, e.g. Q95" (§ Region Weights) | Q95 is unstable on a multi-channel mixture (intergenic / intron / boundary-projected-exon have different scales). | Compute the reference per **density channel** (intergenic, intron, boundary-projected-exon), then take the **mass-weighted maximum** of the per-channel Q95 values. The capture index reports each channel's spread. See §3.5. |
| "Midpoint convention for A(start, ell)" (§ Weighted Effective Length) | Midpoint is fine for short fragments but biased at sharp on/off boundaries when ell > 200 bp. | Keep midpoint convention in v1 production for performance, but verify on diagnostics with a footprint-average reference implementation. Flag in `summary.json` if footprint-average and midpoint diverge by > 5%. See §3.6. |
| "lambda_floor = 0.005" (capture cap c=200) | Hard cap, not justified empirically. | Keep `lambda_floor = 1/c_max` configurable with default `c_max = 200`. Emit a clipping diagnostic in `summary.json` (`n_regions_at_floor`). See §3.5. |
| "Do not propagate variance into Dirichlet strength" (§ Uncertainty) | Correct for v1 scope but contradicts the long-term direction discussed with the user. | Confirm v1 will not modify `gdna_prior_count` *strength* — but the prior magnitude (`η_g`) **does** change because it is computed against weighted exposure. Document this distinction clearly. See §3.7. |
| "Stage 1 denominator-only is acceptable" (§ Per-Fragment Score) | Tempting but creates a non-mass-conserving regime that we then have to remove. | Stage 1 includes both denominator (weighted `L_g`) **and** the per-unit numerator term `log A(midpoint)`. We will not ship a denominator-only intermediate to production. The pure denominator variant exists only inside the Stage 0 diagnostic script for ablation. |

Two additions to v3 that are not optional in v1:

- **Partition-invariance test** (v3 §Tests #2) is mandatory before Stage 1
  enters production; it is the canonical guard against the density² bug that
  the earlier "weight ∝ mass" proposal contained.
- **gDNA-proxy fragment categories must be enumerated explicitly** in code, not
  described in prose. See §3.3 for the exact mask combinations.

---

## 1. Scope and explicit non-goals

### In scope (v1)

1. Per-region shrunk gDNA-proxy density `ρ̂_r` derived from existing v6
   calibration counts (no new scanner work).
2. Per-region relative exposure weight `A_r ∈ [λ_floor, 1]` with global
   reference density.
3. Weighted `L_eff_gdna` for every `MultiLocus`, replacing the unweighted
   `gdna_eff_len_for_loci` result in the prior table.
4. Weighted expected gDNA pseudocount `η̃_g` in `expected_gdna_count_global`.
5. Per-unit additive gDNA log-likelihood term `+ log A(midpoint)` applied in
   Python after scoring (no C++ scorer changes).
6. New diagnostics in `summary.json` and `loci.feather`.
7. A configuration flag to fully disable the regional model (force
   `A_r = 1`) for regression safety.

### Explicit non-goals (deferred)

- Per-locus enrichment ratios `c_L`.
- HMM/change-point target segmentation.
- EM-posterior-driven exposure refinement (circular learning).
- Footprint-average per-fragment weights (midpoint only in v1).
- Applying the same weights to mRNA or nRNA target effective lengths.
- Native-side fast lookup for exposure weights (Python-side per-unit is enough
  for v1; native path is a Stage 3 follow-up if profiling demands it).
- Variance-weighted Dirichlet prior strength.
- Sequence-feature regression (GC, mappability) on `ρ̂_r`.
- External BED target file (will be added as a follow-on once Stage 1/2 land).

---

## 2. Locked mathematical specification

### 2.1 Region partition (unchanged from v6)

The genome is partitioned into disjoint regions of three types (`RegionType`):

- `INTERGENIC` — between genic spans.
- `INTRON` — interstice of a genic span not covered by any exon.
- `EXON` — union of all-strand exons inside a genic span.

`region_df` is provided by the index; the calibration scanner emits per-region
fragment counts and per-region boundary flux. See `src/rigel/calibration/regions.py`
and `src/rigel/calibration/scan_payload.py`.

### 2.2 Per-region gDNA-proxy count `Y_r` and exposure `E_r`

For each region `r`:

| Region type | `Y_r` (numerator) | `E_r` (denominator) | Source |
|---|---|---|---|
| `INTERGENIC` | `payload.per_region_counts[r, MASK_INTERGENIC]` (strand-corrected when strand model active; same formula as v6 `compute_global_densities`) | `l_eff_contained(end_r − start_r, gdna_fl)` | `density_global.py` already emits both |
| `INTRON` | `payload.per_region_counts[r, MASK_INTRON]` (strand-corrected analogously) | `l_eff_contained(end_r − start_r, gdna_fl)` | same |
| `EXON` (boundary-projected) | `u_left[r] + u_right[r]` restricted to eligible sides (`boundary_flux_left/right` flags) | `n_eligible_sides_r × B_cross` (one per eligible side) | `density_global.py` already emits both |

This re-uses the *exact* numerators and denominators that v6 calibration
already computes; v1 only changes whether we use a single per-type density or a
per-region shrunk density.

### 2.3 Empirical-Bayes shrinkage

Per region in class `k ∈ {intergenic, intron, exon}`:

$$
\hat{\rho}_r \;=\; \frac{Y_r + \kappa_k \cdot \rho^{\text{global}}_k}{E_r + \kappa_k}
$$

where `κ_k` and `ρ^global_k` are the existing `kappa` and `rho` fields on
`GlobalGdnaDensity` (one per class). This is **exactly the existing
`density_loco.shrink_to_loco` formula** applied per region. No new estimator.

Properties:

- `Y_r = 0` and `E_r = 0` (impossible region) ⇒ `ρ̂_r = ρ^global_k` (prior wins).
- `E_r → ∞` ⇒ `ρ̂_r → Y_r / E_r` (data wins).
- Coordinate-invariant: splitting a region into two adjacent equal-density
  sub-regions leaves the shrunk densities of both sub-regions equal to the
  parent's shrunk density (proof: `Y` and `E` both split proportionally, the
  prior term shifts identically).

### 2.4 Reference density and per-region weight

Per class `k`, compute a robust per-class reference:

$$
\rho^{\text{ref}}_k \;=\; Q_{0.95}\bigl(\{\hat{\rho}_r : r \in \text{class } k,\; E_r \geq E_{\min}\}\bigr)
$$

with `E_min` (a small absolute exposure threshold, e.g. 100 bp of effective
length) to exclude estimator-noisy tiny regions from the quantile.

The global reference is the mass-weighted max across classes:

$$
\rho^{\text{ref}} \;=\; \max_k \rho^{\text{ref}}_k
$$

Per-region weight:

$$
A_r \;=\; \mathrm{clamp}\!\left(\frac{\hat{\rho}_r}{\rho^{\text{ref}}},\; \lambda_{\text{floor}},\; 1\right)
$$

with `λ_floor = 1 / c_max`. The clamp's upper bound at 1 is essential: it
guarantees `A_r ≤ 1` so the weighted likelihood is always ≤ the unweighted
likelihood, and on uniform libraries we recover the current model exactly
(every region's density is near `ρ^ref`).

**Uniform-library detection**. If the capture index

$$
\text{capture\_index} \;=\; \max_k \left[ Q_{0.95}(\log \hat{\rho}_r^{(k)}) - Q_{0.50}(\log \hat{\rho}_r^{(k)}) \right]
$$

is below a threshold (default `0.6 nats ≈ 1.8×`), set `A_r ≡ 1` for all
regions and tag `mode = "uniform"` in diagnostics. This is the no-op
behavior for polyA / ribominus libraries.

### 2.5 Weighted gDNA effective length for a `MultiLocus`

For each `MultiLocus` `M` (a tuple of half-open `Locus` intervals on
references):

$$
\tilde{L}_{g}(M) \;=\; \sum_{\ell \geq 1} h_g(\ell) \sum_{(p,p+\ell) \cap M \neq \emptyset} A\bigl(\text{midpoint}(p, \ell)\bigr)
$$

where `A(x)` is the piecewise-constant function defined by the region
partition at genomic position `x`. The fragment-length marginal and the
start-window logic (`gdna_eff_len_for_loci`) are unchanged; we only weight
each valid start position by `A` of its midpoint.

Equivalent algorithmic statement (the form we implement):

1. For each constituent `Locus` `[a, b)` on ref `r`:
   1. Restrict `region_df` to regions overlapping `[a, b)`.
   2. For each `(sub_a, sub_b)` = region clipped to `[a, b)`:
      - Compute `eff = l_eff_contained(sub_b - sub_a, gdna_fl)` (existing call)
        + boundary-overlap corrections from the existing
        `gdna_eff_len_for_loci` per-`ell` sum.
      - Multiply each `ell`'s contribution by `A_r` (this is the only new
        operation).
2. Sum across regions and across `Locus` intervals; merge windows by
   reference to avoid double counting (already done in
   `gdna_eff_len_for_loci`).

**Invariant (must be a unit test)**: when all `A_r = 1`, the weighted
function returns the same value as the unweighted
`gdna_eff_len_for_loci` to full float64 precision (we use the same
arithmetic; the weight is a no-op multiplication).

### 2.6 Weighted expected gDNA pseudocount

`expected_gdna_count_global` currently computes `η_g` as a sum of four
mass terms (intergenic-contained, intron-contained, boundary-crossing,
exon-contained). v1 changes each contained term to use the weighted
per-region exposure:

$$
\tilde{\eta}_g(\ell)
\;=\; \rho^{\text{ig}} \sum_{r \in \text{IG}\cap L} A_r \cdot \ell^{\text{core}}_r
\;+\; \rho^{\text{in}} \sum_{r \in \text{IN}\cap L} A_r \cdot \ell^{\text{core}}_r
\;+\; \rho^{\text{b}} \!\left(s_\ell B_{\text{cross}} + \sum_{r \in \text{EX}\cap L} A_r \cdot \ell^{\text{core}}_r\right)
$$

where `ℓ^core_r` is the FL-PMF-weighted exposure of region `r` clipped to
the locus (already produced by `_compute_locus_scratch`). The boundary
crossing term is left unweighted (boundary events are tied to a single
position with eligibility flags, not to a span).

Note: `ρ^ig`, `ρ^in`, `ρ^b` here are the **existing global densities**
(not the per-region shrunk values). The per-region `A_r` is the *only*
new factor. This preserves the EB prior's role as a global prior and
locates capture awareness solely in the exposure layer.

### 2.7 Per-fragment additive log-weight (mass-conservation)

For a fragment `f` with genomic midpoint `m_f`:

$$
\log p_g(f) \;\mathrel{+}=\; \log\max\bigl(A(m_f),\; \lambda_{\text{floor}}\bigr)
$$

This is added in Python to the per-unit `gdna_log_liks` array
produced by the C++ scorer. Mass conservation: when integrated against
the locus opportunity space, the additive term and the weighted
denominator `\tilde{L}_g` cancel out for fragments distributed in
proportion to `A`, leaving the per-fragment posterior balanced.

**Proof sketch.** The gDNA component's contribution to the equivalence-class
log-likelihood is `log p_g(f) - log L̃_g`. With weights:
`log p_g^old(f) + log A(m_f) - log[\sum_p A(p) \cdot h_g(ell)]`. For a fragment
at a position where `A = A_r`, this equals `log p_g^old(f) + log A_r - log L_g
- log(\bar{A})` where `\bar A` is the locus-average weight. The two `log A`
terms produce the correct positional posterior boost in high-weight regions and
demotion in low-weight regions, with the locus-marginal `\bar A` ensuring total
mass is preserved.

---

## 3. Detailed integration with current code

### 3.1 Files touched (summary)

| Path | Change | Type |
|---|---|---|
| `src/rigel/calibration/_regional_exposure.py` | **new module**: `RegionalGdnaExposure` table builder + helpers | add |
| `src/rigel/calibration/_exposure.py` | add `weighted_gdna_eff_len_for_loci()`; keep existing `gdna_eff_len_for_loci` as the `A_r ≡ 1` reference | edit |
| `src/rigel/calibration/locus_prior.py` | thread `RegionalGdnaExposure` through `assemble_priors`; new branch in `expected_gdna_count_global`; new `gdna_eff_len` path | edit |
| `src/rigel/calibration/_result.py` | extend `CalibrationResult` with `regional_exposure: RegionalGdnaExposure | None`; extend summary dict | edit |
| `src/rigel/calibration/_orchestrator.py` | build `RegionalGdnaExposure` after `compute_global_densities`; pass through `build_calibration_result` | edit |
| `src/rigel/config.py` | new `CalibrationConfig.regional_exposure: RegionalExposureConfig` (or extend existing calibration block); CLI surface | edit |
| `src/rigel/cli.py` | expose disable / max-c knobs (or hide behind experimental flag for v1.0) | edit |
| `src/rigel/pipeline.py` | per-unit gDNA log-weight addition; pass `regional_exposure` into `assemble_priors`; surface diagnostics into `loci.feather` | edit |
| `src/rigel/scan.py` | populate per-unit `genomic_midpoint` array (one int64 per unit) for the post-scoring weight lookup | edit |
| `src/rigel/scored_fragments.py` | add `genomic_midpoint: np.ndarray` field on `ScoredFragments` | edit |
| `src/rigel/native/scoring.cpp` | minimal: emit a per-unit `genomic_midpoint` int64 alongside `gdna_log_liks` (we already track `g_sta` and `g_fp` per hit; we just project to a per-unit representative midpoint) | edit, recompile |
| `src/rigel/native/em_solver.cpp` | none for v1 (denominator change reuses existing `locus_gdna_eff_lens` parameter; numerator change is applied in Python before the EM call) | unchanged |
| `tests/test_regional_exposure.py` | **new test module** | add |
| `tests/test_assemble_priors.py` | extend with regional-exposure no-op regression | edit |
| `tests/test_em_impl.py` | partition-invariance integration test | edit |
| `tests/golden/` | regenerate after Stage 1+2 ship | data |
| `scripts/debug/gdna_regional_exposure_diag.py` | Stage 0 diagnostic | add |
| `scripts/benchmark/configs/vcap_no_mm_regional.yaml` | Stage 1/2 benchmark | add |
| `docs/calibration/gdna_exposure_plan_v1.md` | this doc | add |

### 3.2 New module: `src/rigel/calibration/_regional_exposure.py`

Public surface (single source of truth for the exposure layer):

```python
"""Regional gDNA exposure for capture-aware effective length."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from typing import Literal
from ._arrays import RegionArrays, PayloadArrays
from .density_global import GlobalDensityTable
from ..frag_length_model import FragmentLengthModel
from ._exposure import boundary_crossing_exposure

ExposureMode = Literal["uniform", "regional"]


@dataclass(frozen=True)
class RegionalExposureConfig:
    """Configuration for the regional gDNA exposure model."""

    enabled: bool = True
    c_max: float = 200.0                # caps the dynamic range; lambda_floor = 1/c_max
    min_exposure_bp: float = 100.0      # E_min for quantile inclusion
    capture_index_threshold: float = 0.6  # nats; below → mode='uniform'
    min_class_regions_for_quantile: int = 20  # below → fall back to global density
    min_boundary_events_for_exon_refine: int = 5  # for per-exon density refinement


@dataclass(frozen=True)
class RegionalGdnaExposure:
    """Per-region shrunk gDNA density and relative exposure weight.

    All arrays are aligned to ``RegionArrays`` sort order (sorted by
    ``(ref_id, start)``), so a single sorted-position index from
    ``RegionIndexPy`` accesses them safely. Length R.
    """

    rho_shrunk: np.ndarray         # float64, (R,) — shrunk per-region density
    exposure_weight: np.ndarray    # float64, (R,) — A_r in [lambda_floor, 1]
    uncertainty_proxy: np.ndarray  # float64, (R,) — 1/sqrt(E_r + kappa_class)
    mode: ExposureMode             # "uniform" → all A_r = 1
    rho_ref: float                 # the global reference density
    capture_index: float           # nats
    lambda_floor: float
    n_regions_at_floor: int
    n_regions_at_ceiling: int      # count of A_r == 1.0 (informational)
    # Per-class diagnostics, for summary.json
    per_class: dict[str, dict[str, float]]

    @classmethod
    def build(
        cls,
        region_arrays: RegionArrays,
        payload_arrays: PayloadArrays,
        global_densities: GlobalDensityTable,
        gdna_fl: FragmentLengthModel,
        *,
        config: RegionalExposureConfig,
        splicing_anchor_tolerance: int = 0,
    ) -> "RegionalGdnaExposure":
        ...

    @classmethod
    def uniform(cls, n_regions: int) -> "RegionalGdnaExposure":
        """Return the no-op exposure (A_r ≡ 1) for ``n_regions`` regions."""
        ...

    def weight_for_position(self, ref_id: int, pos: int,
                            region_index: "RegionIndexPy") -> float:
        """Look up A(pos) on ref_id via region overlap. O(log R) per call."""
        ...
```

Implementation outline of `RegionalGdnaExposure.build`:

1. Compute per-region `E_r` for INTERGENIC and INTRON via
   `l_eff_contained(end - start, gdna_fl)`. For EXON, use the exon's own
   span only as a fallback `E_r` (the boundary-flux denominator is per
   eligible side, not per region); the per-region density for exon regions
   comes from a separate code path described below.
2. Compute per-region `Y_r` (strand-corrected when the existing strand
   model is active, mirroring `compute_global_densities`).
3. Compute `ρ̂_r = (Y_r + κ_k · ρ_global_k) / (E_r + κ_k)` using
   `density_loco.shrink_to_loco`.
4. For EXON regions: assign `ρ̂_r = ρ_global_exon` (the boundary-derived
   global density), then **optionally refine** with local boundary events
   when `n_boundary_events_in_region ≥ min_boundary_events_for_exon_refine`,
   using the same `shrink_to_loco` formula with `Y = boundary events`,
   `E = n_eligible_sides × B_cross`. This is the v3 §Exon Regions clause
   made precise.
5. Compute per-class `ρ^ref_k = Q_{0.95}({ρ̂_r : E_r ≥ E_min})` when the
   class has ≥ `min_class_regions_for_quantile`; otherwise fall back to the
   global density `ρ_global_k`.
6. `ρ^ref = max_k ρ^ref_k`. Compute `capture_index = max_k IQR_log` (see §2.4).
7. If `capture_index < threshold`, return `cls.uniform(R)` and tag
   `mode='uniform'`.
8. Otherwise compute `A_r = clamp(ρ̂_r / ρ^ref, λ_floor, 1)`.
9. Record per-class diagnostics in `per_class`.

Lookup helper (`weight_for_position`) wraps `RegionIndexPy.overlap` for a
single position; for batched use, callers should iterate the units
themselves (cheap; see §3.8).

### 3.3 Strand-corrected proxy counts: explicit mask

To eliminate ambiguity in "conservative gDNA proxy", v1 uses exactly these
numerators (all already available in `PayloadArrays` / `CalibrationScanPayload`):

- **INTERGENIC**: `payload.per_region_counts[:, MASK_INTERGENIC]`.
- **INTRON**: `payload.per_region_counts[:, MASK_INTRON]`, optionally
  subtracted by the estimated nRNA-source fragment count using the existing
  strand-model correction in `density_global.compute_global_densities`
  (`strand_correction_fragments`). v1 reuses the same correction; if
  `strand_active=False`, no subtraction is applied.
- **EXON (boundary side)**: `payload.u_left + payload.u_right` restricted
  to sides where `region_arrays.bf_left/bf_right` is true (capture-tile
  eligibility from the index). These are the same eligible sides
  `_boundary_term_prorated` uses.

We do **not** use exon-contained unspliced counts (`MASK_EXON_*`) in v1
because the mature RNA contamination problem there is exactly what we are
trying to fix downstream. This matches the v3 decision.

### 3.4 EXON region density: explicit hierarchy

For each EXON region the assigned density is:

1. If `n_eligible_boundary_sides_r ≥ 1` and
   `n_boundary_events_r ≥ min_boundary_events_for_exon_refine`:
   `ρ̂_r = shrink_to_loco(N_b_r, n_eligible × B_cross, ρ_global_exon, κ_exon)`.
2. Else: `ρ̂_r = ρ_global_exon` (no local exon-side evidence; trust the global).

This is the only piece of v3 that was prose-only; codifying it here removes
implementation ambiguity. Step 1 reuses the per-region boundary count we
already compute in `_boundary_term_prorated` (factor out as
`boundary_events_per_region` helper in v1).

### 3.5 Reference density and clipping

- Per-class quantile uses NumPy's `np.quantile(..., method='linear')` after
  filtering on `E_r ≥ E_min`.
- Eligibility check `n_class_regions ≥ min_class_regions_for_quantile`
  (default 20) protects against degenerate samples.
- The cross-class max is recorded as `rho_ref` and emitted in `summary.json`.
- `n_regions_at_floor`, `n_regions_at_ceiling`, and per-class
  density-distribution percentiles are emitted for QC.

### 3.6 Midpoint vs. footprint-average

v1 production: midpoint only.
`midpoint(start, ell) = start + ell // 2`. Lookup is O(log R) via
`RegionIndexPy.overlap(ref_id, midpoint, midpoint + 1)`.

Stage 0 diagnostic script also computes the footprint-average weight
`Ā(start, ell) = mean_{p ∈ [start, start + ell)} A(p)` analytically by
intersecting `[start, start+ell)` with the region partition. We will
sample 1000 representative fragments per hotspot locus and check
`|log A_mid - log Ā| < 0.05 nats` (5%). If the divergence is larger, we
will pull the footprint-average path forward into Stage 2.

### 3.7 `assemble_priors` integration

The single new parameter:

```python
def assemble_priors(
    multi_loci: list[MultiLocus],
    em_data: ScoredFragments,
    index: TranscriptIndex,
    payload: CalibrationScanPayload,
    global_densities: GlobalDensityTable,
    *,
    gdna_fl: FragmentLengthModel | None = None,
    intergenic_flank_bp: int = INTERGENIC_FLANK_BP_DEFAULT,
    splicing_anchor_tolerance: int = 0,
    regional_exposure: RegionalGdnaExposure | None = None,  # <-- NEW
) -> PriorTable:
```

If `regional_exposure is None`, we construct `RegionalGdnaExposure.uniform(R)`
internally (identity transformation); existing behavior is preserved.

Two call sites change inside the per-locus loop:

1. `gdna_eff_len_arr[idx] = weighted_gdna_eff_len_for_loci(ml.loci, ref_lengths_arr, gdna_fl, regional_exposure, region_index, min_value=1.0)`.
2. `expected_gdna_count_global(...)` becomes
   `expected_gdna_count_global(..., regional_exposure=regional_exposure)`
   and applies `A_r` in the contained-term sums.

The boundary-crossing term (`s_ell · B_cross`) is **not** weighted in v1:
boundary events are tied to a single eligible side, not to a spatial span,
and the eligibility flag itself already encodes "is this side a captured
boundary".

Implementation note: `_compute_locus_scratch` already computes
`eff_clip_core` per region. We add a per-locus weighted version:
`eff_clip_core_weighted = eff_clip_core * A_r[region_ids]`. The two
arrays live side by side in `_LocusScratch`.

### 3.8 Per-unit gDNA log-weight (the numerator term)

In `pipeline.quant_from_buffer`, after `scan.FragmentRouter` returns
`ScoredFragments` and before `_run_locus_em_partitioned` is called:

```python
if regional_exposure.mode == "regional":
    # One numpy operation per unit; lookups are O(log R) per unit.
    # Total cost ≈ N_units × log(R) — negligible vs. EM cost.
    weights = _apply_unit_gdna_weights(
        em_data=scored_fragments,
        region_index=region_index,
        exposure=regional_exposure,
        lambda_floor_log=np.log(regional_exposure.lambda_floor),
    )
    # In-place: scored_fragments.gdna_log_liks += weights
```

Where `genomic_midpoint` is a new per-unit `int64` array we propagate from
the scorer (see §3.10). The function adds `log A(midpoint_unit)` to
`gdna_log_liks[unit]`; `-inf` entries (units with no gDNA hypothesis)
remain `-inf`.

For multimappers (NH > 1, the units where `gdna_log_liks` is an LSE over
several hits), v1 approximates with the *primary alignment's* midpoint.
The error is bounded because all hits of a multimapper share a candidate
neighborhood and tend to be in similar exposure regimes. Profiling on
`include_multimap=true` runs is needed before promoting this to a
per-hit weighted LSE inside the scorer; that promotion is Stage 3.

### 3.9 `gdna_eff_len` invariant when `mode == 'uniform'`

The integration test is:

> Build a `RegionalGdnaExposure` with `mode='uniform'` (synthetic or via
> the auto-detector on a uniform input). Run `assemble_priors` twice:
> once with `regional_exposure=None` and once with the uniform exposure
> table. Assert `np.array_equal(gdna_eff_len_a, gdna_eff_len_b)`
> bit-exactly and `np.array_equal(gdna_prior_count_a, gdna_prior_count_b)`
> bit-exactly.

This is the regression guard against any unintended drift in the no-op
path. It must pass before Stage 1 enters main.

### 3.10 Native ABI change: per-unit `genomic_midpoint`

Today the C++ scorer's `score_chunk` returns an LSE'd `gdna_log_liks`
per unit but no per-unit position. We need a single new per-unit array:

- **What**: `genomic_midpoint: int64[N_units]`, one representative
  genomic midpoint per unit. For unambiguous and ambig-strand fragments,
  this is the unit's only fragment's midpoint. For multimappers, the
  primary alignment's midpoint.
- **Where**: emitted by `score_chunk` alongside `gdna_log_liks` and
  threaded into `ScoredFragments`.
- **Negative sentinel**: `INT64_MIN` for units with no gDNA hypothesis
  (spliced-only or chimeric units). Python ignores these in the lookup.

This is a minimal-surface change. The C++ scorer already computes per-hit
genomic start (`cp.g_sta`) and footprint (`cp.g_fp`); the midpoint is
`g_sta + g_fp / 2`.

Two compatibility paths:

- If a user runs an old `scored_fragments.feather` without the new column,
  we fall back to `mode='uniform'` and warn.
- The new column is appended to the existing tuple returned by
  `StreamingScorer.finish()`; the Python side adds it positionally.

### 3.11 `summary.json` additions

```json
"calibration": {
  ...
  "regional_exposure": {
    "mode": "regional",
    "enabled": true,
    "capture_index": 2.43,
    "rho_ref": 1.81e-04,
    "lambda_floor": 0.005,
    "n_regions": 296039,
    "n_regions_at_floor": 142385,
    "n_regions_at_ceiling": 3120,
    "per_class": {
      "INTERGENIC": {
        "n_regions_used": 33120,
        "rho_quantiles": {"q05": 4.0e-6, "q50": 7.9e-6, "q95": 1.2e-5},
        "rho_ref_class": 1.2e-5
      },
      "INTRON":      {"n_regions_used": 259622, "rho_quantiles": {...}, "rho_ref_class": 8.5e-5},
      "EXON":        {"n_regions_used": 309202, "rho_quantiles": {...}, "rho_ref_class": 1.81e-4}
    }
  }
}
```

### 3.12 `loci.feather` additions

Six new columns (informational; do not alter existing schemas in place):

| Column | Type | Meaning |
|---|---|---|
| `gdna_eff_len_unweighted` | float64 | Same formula the v6 pipeline would emit today (call `gdna_eff_len_for_loci` with `A_r ≡ 1`). |
| `gdna_eff_len` | float64 | **Reused name**; now refers to weighted. |
| `gdna_eff_len_weight_ratio` | float64 | `gdna_eff_len / gdna_eff_len_unweighted`. ∈ (0, 1]. |
| `gdna_prior_count_unweighted` | float64 | Same as today; informational. |
| `gdna_prior_count` | float64 | **Reused name**; now weighted. |
| `n_regions_in_locus` | int64 | Diagnostic. |

The `gdna_eff_len` and `gdna_prior_count` column names keep their meaning
(the EM-active values) but their values change when `mode='regional'`. The
`*_unweighted` siblings allow direct measurement of the weight ratio per
locus, which is the actionable diagnostic for the VCaP hotspots.

### 3.13 Config and CLI

In `src/rigel/config.py`:

```python
@dataclass(frozen=True)
class CalibrationConfig:
    ...
    regional_exposure: "RegionalExposureConfig" = field(
        default_factory=RegionalExposureConfig
    )
```

CLI exposure (in `src/rigel/cli.py`):

- `--regional-exposure {auto,off}` (default `auto`) — `off` forces
  `enabled=False` and the model is the v6 unweighted model.
- `--regional-exposure-c-max FLOAT` (default 200.0; advanced).

That is the entire user-facing surface. All other tuning knobs
(`min_exposure_bp`, `capture_index_threshold`,
`min_class_regions_for_quantile`, `min_boundary_events_for_exon_refine`) are
config-only and not advertised as CLI flags in v1.

---

## 4. Stage-by-stage execution plan

### Stage 0 — Diagnostic script (no production code changes)

**Path**: `scripts/debug/gdna_regional_exposure_diag.py`.

**Purpose**: validate the hypothesis on the five VCaP hotspot loci
before any production code change.

**Inputs**: index dir, calibration outputs (we can either rerun
`calibrate(...)` on the no_mm payload or load the existing
`summary.json` plus the index `region_df` to reconstruct per-region
counts; rerunning calibrate is preferred and is cheap).

**Outputs** (`results/vcap_regional_exposure_diag_2026-05-18/`):
- Per-locus table for the five hotspots:
  - `gdna_eff_len_unweighted`
  - `gdna_eff_len_weighted` (midpoint convention)
  - `gdna_eff_len_weighted_footprint` (footprint-average reference)
  - `predicted_log_post_gain = log(unweighted / weighted)`
  - `top_5_regions_by_weight`, `bottom_5_regions_by_weight`
  - `n_regions_in_locus`, `mean_A`, `median_A`, `Q05_A`, `Q95_A`
- Per-fragment table for 500 sampled hotspot FP fragments:
  - `A(midpoint)`, `A_footprint`, `predicted_post_gain_per_fragment`
  - Predicted gDNA posterior share assuming the EM converged to the
    same θ ratios.
- Global summary: `capture_index`, `rho_ref` per class, `mode`.

**Acceptance criterion before promoting to Stage 1**:
On the five hotspot loci, the predicted gDNA log-posterior gain must be
≥ 2.0 nats on average (equivalent to ≥ 7× swing in posterior odds).
If smaller, stop and investigate before changing production code.

### Stage 1 — Production exposure layer (no per-fragment numerator)

**Goal**: ship the weighted denominator + weighted `η_g` with strict
backward compatibility on uniform libraries.

**Tasks**:

1. Implement `_regional_exposure.py` per §3.2.
2. Implement `weighted_gdna_eff_len_for_loci` in `_exposure.py`. Internally
   reuse the existing `gdna_eff_len_for_loci` per-`ell` window logic; add
   a per-window weight multiplication using `RegionIndexPy.overlap` for
   each merged window's midpoint. Add unit tests:
   - Identity test (A ≡ 1) per-locus.
   - Partition-invariance test: split a region; weighted length unchanged.
   - Two-state test: 100 kb locus, 10 kb at A=1 and 90 kb at A=0.01,
     `c_max=200`, asserts `L̃_g ≈ 10 + 90/200 = 10.45 kb` to ≤ 1%
     tolerance (boundary FL terms).
3. Extend `expected_gdna_count_global` with `regional_exposure` kwarg.
4. Extend `assemble_priors` to build / accept `RegionalGdnaExposure`.
5. Extend `CalibrationResult` and `_orchestrator.calibrate` to surface
   the table.
6. Add `regional_exposure` summary block to `to_summary_dict`.
7. Add new `loci.feather` columns; do not remove old columns.
8. CLI/config plumbing per §3.13.
9. Run `pytest tests/ -v` (must pass with all existing 1074+ tests
   green; new tests added).
10. Run synthetic benchmarks (`scripts/benchmark/configs/locus_simple_*.yaml`)
    and confirm zero drift on uniform inputs.
11. Run VCaP no_mm benchmark with `--regional-exposure auto` and report
    `gdna_eff_len_weight_ratio` per hotspot locus.

**Ship gate for Stage 1**:
- All tests pass.
- Uniform-input regression delta < 1e-9 on per-locus `gdna_eff_len` and
  `gdna_prior_count`.
- VCaP `gdna_eff_len_weight_ratio` at the five hotspot loci ≤ 0.1
  (i.e. ≥ 10× shrinkage), as predicted by Stage 0.

### Stage 2 — Per-fragment numerator term

**Goal**: complete the mass-conserving model.

**Tasks**:

1. Add `genomic_midpoint: int64[N_units]` to the
   `StreamingScorer.finish()` return tuple (`scoring.cpp` change;
   recompile per `pip install --no-build-isolation -e .`).
2. Plumb through `scan.py` → `ScoredFragments` (`scored_fragments.py`).
3. Add `_apply_unit_gdna_weights(em_data, region_index, exposure,
   lambda_floor_log) -> None` in `pipeline.py` (in-place modification of
   `gdna_log_liks`).
4. Call it in `quant_from_buffer` after `assemble_priors` and before
   `_run_locus_em_partitioned`. Guard with `if exposure.mode == 'regional'`.
5. New tests:
   - Unit: `_apply_unit_gdna_weights` adds the correct log-weight per
     unit; `-inf` units stay `-inf`; out-of-range midpoints (sentinel)
     are no-ops.
   - Integration: VCaP-style synthetic with bimodal A; assert that
     Stage 2 reduces gDNA-routed posterior FP rate vs. Stage 1.
6. Regenerate `tests/golden/` outputs that depend on per-locus EM mass
   (these are the ones that legitimately change under Stage 2).
7. Re-run synthetic benchmarks; confirm zero drift on uniform inputs
   (the additive per-unit term is `log A = log 1 = 0` for `mode='uniform'`).
8. Re-run VCaP no_mm; report gDNA-source FP count vs Stage 1 and vs.
   pre-v1 baseline (target: < 500K FPs from the current 2.64M, per
   yesterday's deep-dive math).

**Ship gate for Stage 2**:
- All tests pass.
- Uniform-input regression delta < 1e-9 on `gdna_log_liks`.
- VCaP no_mm FP count drops by ≥ 70% from baseline.

### Stage 3 — Deferred items (only if Stage 2 results warrant)

- Footprint-average per-fragment weight (if the midpoint-vs-footprint
  divergence in Stage 0 exceeded 5%).
- Per-hit weighting inside the scorer for multimappers.
- Coverage-segmentation-based exon weight refinement.
- Apply weights to mRNA / nRNA effective lengths.
- Variance-weighted Dirichlet prior strength.
- BED-file optional input.

Each item has its own micro-plan; none should be tackled while Stage 2 is
unmerged.

---

## 5. Test plan (mandatory before each stage ships)

### 5.1 Unit tests (`tests/test_regional_exposure.py` — new)

| ID | Test | Assertion |
|---|---|---|
| U1 | `uniform()` constructor on R=10 | All `A_r == 1`, mode='uniform' |
| U2 | `build()` on a tiny synthetic payload with all-equal densities | `capture_index < threshold`, `mode='uniform'` |
| U3 | `build()` on a bimodal density (10 high-density regions + 90 low-density), `c=100` | `A_r ≈ {1.0, 0.01}`; `capture_index ≈ log(100) ≈ 4.6` |
| U4 | Shrinkage: a tiny region with 0 counts and 1 with 1000 counts but small `E` | Tiny region's `ρ̂_r ≈ ρ_global`; large region trusts data |
| U5 | Reference: per-class quantile with `n_class_regions < min` | Falls back to `ρ_global_k` |
| U6 | `weight_for_position` correctness on a synthetic partition | Returns `A_r` for `r` containing `pos` |

### 5.2 `weighted_gdna_eff_len_for_loci` tests

| ID | Test | Assertion |
|---|---|---|
| W1 | Identity: `A ≡ 1` | Bit-exact match with `gdna_eff_len_for_loci` |
| W2 | Partition invariance: split one region into two adjacent equal-A pieces | Result unchanged |
| W3 | Two-state: 100 kb locus, 10 kb at A=1, 90 kb at A=0.01 | `L̃ ≈ 10 + 90/100` within 1% (FL boundary terms) |
| W4 | Multi-`Locus` `MultiLocus` with overlapping windows | Sum equals single-Locus call on the merged span |

### 5.3 `assemble_priors` regression (`tests/test_assemble_priors.py` extension)

| ID | Test | Assertion |
|---|---|---|
| A1 | `regional_exposure=None` reproduces v6 prior table | Bit-exact `gdna_eff_len`, `gdna_prior_count` |
| A2 | `regional_exposure=uniform(R)` reproduces v6 prior table | Bit-exact |
| A3 | `regional_exposure=build(...)` on uniform-input synthetic | Bit-exact (because auto-detect → uniform) |
| A4 | `regional_exposure=build(...)` on bimodal-input synthetic | `gdna_eff_len < unweighted`; `gdna_prior_count` reduced proportionally |

### 5.4 Stage 2 integration

| ID | Test | Assertion |
|---|---|---|
| S1 | Per-unit weight no-op when `mode='uniform'` | `gdna_log_liks` unchanged |
| S2 | Per-unit weight on bimodal synthetic | Per-unit deltas equal `log A(midpoint)` |
| S3 | `-inf` units untouched | `gdna_log_liks[-inf_units]` still `-inf` after call |
| S4 | End-to-end: bimodal synthetic, capture-style | gDNA posterior mass on truly-gDNA fragments increases by ≥ 50% vs. Stage 1 |

### 5.5 Benchmarks (must pass before main)

- `scripts/benchmark/configs/locus_simple_baseline.yaml` and all other
  `locus_simple_*` configs: zero drift (within Δ < 0.5% on mRNA/nRNA/gDNA
  relative error) on uniform inputs.
- VCaP no_mm: FP count and per-hotspot ZW distribution emit in a new
  benchmark report `docs/benchmarks/vcap_no_mm_regional_2026-05-XX.md`.
- VCaP with_mm: confirm multimapper-included case also improves (or at
  least does not regress).

---

## 6. Risk register and rollback

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Auto-detect mis-classifies a capture sample as uniform | Medium | Loss of correction in real use | Capture index emitted in `summary.json`; threshold is config-tunable. User can force `--regional-exposure auto` regardless; users can set explicit mode via env var if needed. |
| Auto-detect mis-classifies a uniform sample as capture | Medium | Mild posterior drift on a fraction of loci | `A_r` is clipped to `[lambda_floor, 1]`; the absolute drift on uniform-like samples is bounded by `(1 - lambda_floor)` per region. Test on synthetic suite catches > 0.5% drift. |
| Per-region exon density misestimated when boundary evidence is absent | High in panels | Wrong A on exon regions of unannotated/untargeted genes | Fallback to `ρ_global_exon`; the dominant case (well-covered captured exons) has boundary support and refines correctly. |
| `genomic_midpoint` plumbing breaks an old `scored_fragments` consumer | Low | Test failure | Add field at end of tuple; positional access only. |
| Stage 1 ships denominator-only by accident, breaking mass conservation in real data | Low (we explicitly chose not to ship Stage 1 standalone to production) | Posterior drift in favor of gDNA in low-A regions | CI gate: do not enable `regional_exposure` in default config until Stage 2 lands. |
| Per-unit Python lookup is too slow on huge buffers | Low | Wall-time regression | Profiled at design time: `N_units × log(R)` ≈ 50M × 25 = 1.25G ops, ≈ 5 s in pure Python; vectorizable with sorted-unit ref grouping if needed. |

**Rollback path**: any of the three layers can be disabled independently:

- **Stage 2 rollback**: set `regional_exposure.enabled = False` in config or
  `--regional-exposure off` on CLI; per-unit weights are not applied and
  the prior table reverts to unweighted.
- **Stage 1 rollback** (post-merge bug): same flag also disables weighted
  `L_g` / `η_g` and we fall back to the v6 baseline exactly.
- **ABI rollback**: the only ABI change is the new per-unit
  `genomic_midpoint` field. Reverting to `None` everywhere also restores
  the pre-v1 behavior bit-exactly.

---

## 7. Effort estimate (informational)

| Stage | Stage-level effort | Why |
|---|---|---|
| Stage 0 (diagnostic) | small | Pure Python; reuses calibration outputs |
| Stage 1 (denominator + prior) | medium | New module, two integration points, ~6 tests, regression harness |
| Stage 2 (per-fragment) | small-to-medium | One ABI change, one new pipeline call, two new tests |
| Stage 3 (deferred) | varies | Each item is its own micro-project |

Total new test files: 1. Modified test files: 3-4. Modified source files: 10.
Recompile required at Stage 2 (native ABI change).

---

## 8. Open questions to answer in Stage 0

These are answered by running the Stage 0 diagnostic; the answers feed
config defaults.

1. What is the empirical `capture_index` for the VCaP no_mm BAM?
   (Hypothesis: 2.5–4.5 nats. If < 0.6, the threshold needs lowering.)
2. What is the empirical `ρ^ref` and how does it compare across the
   three classes?
3. What fraction of regions land at `lambda_floor`? (If > 80%, the floor
   may be too low for real loci.)
4. On the five hotspot loci, does the midpoint convention give the same
   weights as the footprint-average within 5%?
5. What is the predicted per-fragment gDNA log-posterior gain on the five
   hotspots, and does it exceed the 2-nat ship gate?

The Stage 0 script emits a written report answering each of these so
Stage 1 can be tuned and started with empirical defaults.

---

## 9. Appendix — formula consolidation

For ease of reference during implementation, the locked v1 formulas:

```
# 1. Per-region shrunk density (per class k):
rho_hat_r = (Y_r + kappa_k * rho_global_k) / (E_r + kappa_k)

# 2. Per-class reference (filter on E_r >= E_min, require >= min_class regions):
rho_ref_k = quantile_95(rho_hat_r for r in class k)
rho_ref   = max(rho_ref_k for k in {ig, in, ex})

# 3. Capture index (detect uniform vs. regional):
capture_index = max over k of [Q95(log rho_hat_r) - Q50(log rho_hat_r)]
mode = "uniform" if capture_index < threshold else "regional"

# 4. Per-region weight:
A_r = clamp(rho_hat_r / rho_ref, lambda_floor, 1.0)

# 5. Weighted MultiLocus gDNA effective length:
L_tilde_g(M) = sum_ell h_g(ell)
               * sum over merged start windows w of
                   (w_end - w_start) * A_r_containing_midpoint(w, ell)

# 6. Weighted expected gDNA pseudocount (per Locus L):
eta_tilde_g(L) =
    rho_global_ig * sum_{r in IG cap L} A_r * l_core_r
  + rho_global_in * sum_{r in IN cap L} A_r * l_core_r
  + rho_global_ex * sum_{r in EX cap L} A_r * l_core_r
  + rho_global_ex * s_ell * B_cross   # boundary term unweighted

# 7. Per-unit gDNA log-likelihood update (applied in Python post-scoring):
gdna_log_liks[u] += max(log A(midpoint(u)), log lambda_floor)
   if mode == "regional" and gdna_log_liks[u] is finite
```

These are the only seven equations the implementation needs to realise.
Everything else is plumbing.
