# M6 — Per-`Locus` gDNA mass + `MultiLocus` priors

**Predecessor:** [calibration_v6_plan.md §2.6.1, §2.7, §2.8](calibration_v6_plan.md).
**Depends on (all green):** M2 (`region_df` on `TranscriptIndex`), M3
(`CalibrationScanPayload`), M4 (`GlobalDensityTable`, `KappaEstimate`,
`l_eff_overlap`), M5 (`MultiLocus`/`Locus`, `prior_weight_rna` ABI).

**Goal:** the numerical core of v6 calibration. Convert the global
density block + per-region payload into per-`Locus` gDNA mass
estimates, aggregate to `MultiLocus`-level Dirichlet priors
(`α_gdna`, `α_rna`), and assemble the `PriorTable` consumed by the
batch EM. Pure-Python pass; no new C++.

---

## 1. Surface contract

After M6, calibration callers can write:

```python
from rigel.calibration import (
    assemble_priors, PriorTable, MultiLocusPrior, LocusGdnaEstimate,
    estimate_locus_gdna, assemble_multilocus_prior,
    shrink_to_loco, partition_units_to_loci,
    C_BASE_DEFAULT,
    FLAG_INTERGENIC_ZERO_LEFF, FLAG_INTRON_ZERO_LEFF,
    FLAG_EXON_INTRON_NO_ELIGIBLE, FLAG_PI_CLIPPED,
)

prior_table = assemble_priors(
    multi_loci=multi_loci,
    em_data=em_data,
    index=index,
    payload=payload,
    global_densities=global_densities,
    gdna_fl_mean=global_densities.gdna_fl_mean,
    c_base=C_BASE_DEFAULT,
)
# prior_table.alpha_gdna, .alpha_rna, .prior_weight_rna ⇒ batch EM kwargs
```

The legacy `compute_locus_priors_from_partitions(...)` path stays alive
until M8b switches the pipeline over. **No call-site migration in M6.**

---

## 2. Module layout

Five new files, one extension to an existing one:

```
src/rigel/calibration/
├── density_loco.py            (new) shrink_to_loco — closed-form NB EB
├── _arrays.py                 (new) RegionArrays, PayloadArrays — pre-flattened views
├── _region_index_py.py        (new) per-ref searchsorted region overlap
├── _locus_n_obs.py            (new) partition_units_to_loci (anchor routing)
└── locus_prior.py             (extend) LocusGdnaEstimate, MultiLocusPrior,
                                        PriorTable, estimate_locus_gdna,
                                        assemble_multilocus_prior, assemble_priors
```

**Why a Python `_region_index_py.py` instead of reusing the M3 native
`RegionIndex`?** The M3 index lives behind the `BamScanner` ABI and
carries scanner-side state. The locus-prior pass runs in Python, on a
per-locus loop dominated by O(`len(ml.loci)`) calls (≈1.0 average,
≤3 in the long tail). A 30-line `numpy.searchsorted`-keyed lookup is
faster than the round-trip cost of a nanobind call per query and
keeps the locus-prior assembly entirely in pure Python — easier to
test, easier to vectorize later. We keep the option to swap for the
native `RegionIndex` once profiling justifies it.

---

## 3. Algorithms

### 3.1 `shrink_to_loco` — closed-form NB EB

**File:** `src/rigel/calibration/density_loco.py`

Per [§2.6.1](calibration_v6_plan.md#261-two-level-eb-shrinkage):

```python
def shrink_to_loco(
    n_loco: float,
    leff_loco: float,
    rho_global: float,
    kappa: float,
) -> float:
    """ρ̂_loco = (N + κ·ρ̂_global) / (L_eff + κ)."""
    denom = leff_loco + kappa
    if denom <= 0.0:
        return rho_global   # degenerate locus → fall back to global
    return (n_loco + kappa * rho_global) / denom
```

Pure scalar; no numpy import needed. Vectorization is unwarranted
(the per-locus pass calls this 3 times per `Locus`, avg `len(ml.loci)
≈ 1`). One-line tests pin the (κ→∞ ⇒ ρ_global), (κ→0 ⇒ N/L_eff),
(L_eff=0 ⇒ ρ_global) limits.

### 3.2 `RegionArrays` / `PayloadArrays` — pre-flattened views

**File:** `src/rigel/calibration/_arrays.py`

Per-call `region_df["start"].to_numpy()` is the dominant overhead at
locus-prior assembly time on real indexes (~1.5M regions). We
materialize the five hot columns once:

```python
@dataclass(frozen=True, slots=True)
class RegionArrays:
    ref_id:         np.ndarray   # int32, (R,)
    start:          np.ndarray   # int64, (R,)
    end:            np.ndarray   # int64, (R,)
    type:           np.ndarray   # uint8, (R,)
    bf_left:        np.ndarray   # bool,  (R,)
    bf_right:       np.ndarray   # bool,  (R,)
    # Per-ref CSR: (n_refs + 1,) int32 offsets so that
    # ref_id[ref_offsets[r]:ref_offsets[r+1]] is sorted by start.
    ref_offsets:    np.ndarray   # int32, (n_refs + 1,)

    @classmethod
    def from_region_df(cls, region_df: pd.DataFrame, ref_name_to_id: dict[str,int]) -> Self: ...

@dataclass(frozen=True, slots=True)
class PayloadArrays:
    intergenic_per_region:  np.ndarray  # int64, (R,) — per_region_counts[:, MASK_INTERGENIC]
    intron_per_region:      np.ndarray  # int64, (R,)
    u_left:                 np.ndarray  # int64, (R,)
    u_right:                np.ndarray  # int64, (R,)

    @classmethod
    def from_payload(cls, payload: CalibrationScanPayload) -> Self: ...
```

`RegionArrays.from_region_df` translates `ref_name → ref_id` once
(via the index's `ref_name_to_id`) and validates per-ref sortedness
in DEBUG mode only (the M2 invariant guarantees it; we trust and
move on in production).

### 3.3 `_region_index_py.RegionIndexPy` — per-ref overlap

**File:** `src/rigel/calibration/_region_index_py.py`

Per-reference CSR over sorted region intervals; binary search +
linear scan. No external deps.

```python
@dataclass(frozen=True, slots=True)
class RegionIndexPy:
    arrays: RegionArrays

    def overlap(self, ref_id: int, qstart: int, qend: int) -> np.ndarray:
        """Return int32[] of region_ids overlapping [qstart, qend) on ref_id.

        Half-open interval; uses np.searchsorted on the per-ref slice
        of `arrays.start` (right) and `arrays.end` (left). Result is
        sorted ascending by region_id (== ascending by start within a
        ref by M2 invariant).
        """
```

Implementation (≈ 25 LoC):

```python
def overlap(self, ref_id, qstart, qend):
    a = self.arrays
    lo = a.ref_offsets[ref_id]
    hi = a.ref_offsets[ref_id + 1]
    if lo == hi:
        return np.empty(0, dtype=np.int32)
    starts = a.start[lo:hi]
    ends   = a.end  [lo:hi]
    # First region with start >= qend : exclusive upper bound.
    j_hi = lo + np.searchsorted(starts, qend, side="left")
    # Linear walk back from j_hi to find the first row with end > qstart.
    # Per-ref regions are non-overlapping ⇒ the overlapping window is contiguous.
    # Using searchsorted on `ends` directly is OK because per-ref ends are
    # also sorted (regions don't overlap), so:
    j_lo = lo + np.searchsorted(ends, qstart, side="right")
    if j_lo >= j_hi:
        return np.empty(0, dtype=np.int32)
    return np.arange(j_lo, j_hi, dtype=np.int32)
```

Two `searchsorted` calls per `overlap` invocation. The "ends sorted
because regions are non-overlapping" trick is the cleanup win — no
linear scan in the hot path.

### 3.4 `partition_units_to_loci` — anchor-transcript routing

**File:** `src/rigel/calibration/_locus_n_obs.py`

Per [§2.8](calibration_v6_plan.md#28-n_obsl--anchor-transcript-routing).
Routes EM units within a `MultiLocus` to its constituent `Locus`
intervals.

```python
def partition_units_to_loci(
    ml: MultiLocus,
    em_data: ScoredFragments,
    t_to_locus_idx: np.ndarray,   # int8, (N_T,) — local Locus index per global tid
) -> tuple[np.ndarray, ...]:
    """Return one int32[] of unit indices per Locus in ml.loci.

    Fast path (len(ml.loci) == 1, ≈99% of MultiLoci): return
    (ml.unit_indices,) directly — zero allocation.

    Slow path: gather anchor transcripts, then bin into local Locus
    indices via the precomputed t_to_locus_idx.
    """
    if len(ml.loci) == 1:
        return (ml.unit_indices,)

    anchor_t = em_data.locus_t_indices[ml.unit_indices]    # int32, (n_units,)
    bins     = t_to_locus_idx[anchor_t]                    # int8, (n_units,)
    # Sanity: every anchor must map to a local Locus.
    if (bins < 0).any():
        raise RuntimeError(
            "partition_units_to_loci: anchor transcript not in any Locus "
            "of MultiLocus(multi_locus_id={ml.multi_locus_id}). "
            "build_multi_loci invariant violation."
        )
    return tuple(
        ml.unit_indices[bins == j].astype(np.int32, copy=False)
        for j in range(len(ml.loci))
    )
```

`t_to_locus_idx` is built once per `MultiLocus` by `assemble_priors`
(below) — it is `int8` because `len(ml.loci) ≤ 127` is overwhelmingly
true (the largest paralog cluster in human GENCODE is 19 loci).

**Anti-pattern explicitly forbidden:** materializing per-fragment
genomic coordinates from the buffer to bin units geometrically. See
`/memories/repo/multilocus-partition-design-2026-04.md`.

### 3.5 `estimate_locus_gdna` — per-`Locus` core

**File:** `src/rigel/calibration/locus_prior.py` (extend the M5 stub).

Per [§2.7.1](calibration_v6_plan.md#271-algorithm-per-locus). One
`Locus` in, one `LocusGdnaEstimate` out:

```python
def estimate_locus_gdna(
    locus:            Locus,
    n_obs:            int,
    region_index:     RegionIndexPy,
    region_arrays:    RegionArrays,
    payload_arrays:   PayloadArrays,
    global_densities: GlobalDensityTable,
    gdna_fl_mean:     float,
) -> LocusGdnaEstimate:
    region_ids = region_index.overlap(locus.ref_id, locus.start, locus.end)
    if region_ids.size == 0:
        raise RuntimeError(
            f"Locus(ref={locus.ref}, start={locus.start}, end={locus.end}) "
            f"overlaps no regions. BAM reference does not match index — "
            f"rebuild the index or check the BAM header."
        )

    types   = region_arrays.type[region_ids]
    starts  = region_arrays.start[region_ids]
    ends    = region_arrays.end[region_ids]

    # Clip each region to the locus interval.
    cl_lo = np.maximum(starts, locus.start)
    cl_hi = np.minimum(ends, locus.end)
    cl_len = (cl_hi - cl_lo).astype(np.float64)            # > 0 by overlap query
    leff   = cl_len + (gdna_fl_mean - 1.0)                 # l_eff_overlap

    # --- Per-type contributions ---
    fallback_flags = 0

    # INTERGENIC
    m = (types == int(RegionType.INTERGENIC))
    n_ig, leff_ig, rho_loco_ig = _shrink_one(
        m, region_ids, payload_arrays.intergenic_per_region, leff,
        global_densities.intergenic,
    )
    n_gdna_ig = rho_loco_ig * leff_ig
    if leff_ig <= 0.0: fallback_flags |= FLAG_INTERGENIC_ZERO_LEFF

    # INTRON — same pattern
    m = (types == int(RegionType.INTRON))
    n_in, leff_in, rho_loco_in = _shrink_one(
        m, region_ids, payload_arrays.intron_per_region, leff,
        global_densities.intron,
    )
    n_gdna_in = rho_loco_in * leff_in
    if leff_in <= 0.0: fallback_flags |= FLAG_INTRON_ZERO_LEFF

    # EXON-INTRON — boundary-flux density × full exonic L_eff in locus.
    n_ei, leff_ei, n_eligible_boundaries, rho_loco_ei = _shrink_exon_intron(
        types, region_ids, region_arrays, payload_arrays, leff,
        global_densities.exon_intron, gdna_fl_mean,
    )
    n_gdna_ei = rho_loco_ei * leff_ei
    if n_eligible_boundaries == 0: fallback_flags |= FLAG_EXON_INTRON_NO_ELIGIBLE

    n_gdna_total = n_gdna_ig + n_gdna_in + n_gdna_ei
    pi_unclipped = n_gdna_total / n_obs if n_obs > 0 else 0.0
    pi_gdna = float(np.clip(pi_unclipped, 0.0, 1.0))
    if pi_unclipped > 1.0: fallback_flags |= FLAG_PI_CLIPPED

    return LocusGdnaEstimate(
        locus=locus, n_obs=n_obs,
        n_gdna_intergenic=n_gdna_ig,
        n_gdna_intron=n_gdna_in,
        n_gdna_exon_intron=n_gdna_ei,
        n_gdna=n_gdna_total,
        pi_gdna=pi_gdna,
        rho_loco=(rho_loco_ig, rho_loco_in, rho_loco_ei),
        leff_loco=(leff_ig, leff_in, leff_ei),
        n_eligible_boundaries=int(n_eligible_boundaries),
        fallback_flags=fallback_flags,
    )
```

Helpers (file-local):

```python
def _shrink_one(
    type_mask, region_ids, count_col, leff, global_d,
):
    rids = region_ids[type_mask]
    if rids.size == 0:
        return 0.0, 0.0, global_d.rho
    n_loco    = float(count_col[rids].sum())
    leff_loco = float(leff[type_mask].sum())
    rho_loco  = shrink_to_loco(n_loco, leff_loco, global_d.rho, global_d.kappa.value)
    return n_loco, leff_loco, rho_loco

def _shrink_exon_intron(
    types, region_ids, region_arrays, payload_arrays, leff, global_d, gdna_fl_mean,
):
    m = (types == int(RegionType.EXON))
    if not m.any():
        return 0.0, 0.0, 0, global_d.rho
    rids = region_ids[m]
    bf_l = region_arrays.bf_left [rids]
    bf_r = region_arrays.bf_right[rids]
    n_eligible = int(bf_l.sum() + bf_r.sum())
    leff_full  = float(leff[m].sum())   # full exonic L_eff inside locus

    if n_eligible == 0:
        # No flag-eligible boundaries in this locus → use global density,
        # but still multiply by full exonic L_eff (mass we expect, not
        # per-eligible-boundary mass).
        return 0.0, leff_full, 0, global_d.rho

    n_loco_eligible    = int(
        (payload_arrays.u_left [rids] * bf_l).sum() +
        (payload_arrays.u_right[rids] * bf_r).sum()
    )
    leff_loco_eligible = float(n_eligible) * gdna_fl_mean
    rho_loco = shrink_to_loco(
        n_loco_eligible, leff_loco_eligible, global_d.rho, global_d.kappa.value,
    )
    return n_loco_eligible, leff_full, n_eligible, rho_loco
```

**Critical:** for EXON-INTRON the **density** uses the eligible
boundary slice (capture-aware), but the **predicted count** uses
`leff_full` (the full exonic L_eff inside the locus). This is the
locked plan §2.7.1 step 3.

### 3.6 `assemble_multilocus_prior` — `MultiLocus` aggregation

```python
def assemble_multilocus_prior(
    ml:                    MultiLocus,
    per_locus_estimates:   tuple[LocusGdnaEstimate, ...],
) -> MultiLocusPrior:
    n_obs  = sum(e.n_obs  for e in per_locus_estimates)
    n_gdna = sum(e.n_gdna for e in per_locus_estimates)
    n_rna  = max(0.0, n_obs - n_gdna)
    pi_gdna = (n_gdna / n_obs) if n_obs > 0 else 0.0
    return MultiLocusPrior(
        multi_locus_id=ml.multi_locus_id,
        n_obs=n_obs, n_gdna=n_gdna, n_rna=n_rna,
        pi_gdna=float(np.clip(pi_gdna, 0.0, 1.0)),
        per_locus=per_locus_estimates,
    )
```

### 3.7 `assemble_priors` — orchestrator

```python
def assemble_priors(
    multi_loci:       list[MultiLocus],
    em_data:          ScoredFragments,
    index:            TranscriptIndex,
    payload:          CalibrationScanPayload,
    global_densities: GlobalDensityTable,
    *,
    gdna_fl_mean:     float | None = None,
    c_base:           float = C_BASE_DEFAULT,
) -> PriorTable:
    if gdna_fl_mean is None:
        gdna_fl_mean = global_densities.gdna_fl_mean

    region_arrays  = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(payload)
    region_index   = RegionIndexPy(arrays=region_arrays)

    n_loci = len(multi_loci)
    alpha_gdna       = np.zeros(n_loci, dtype=np.float64)
    alpha_rna        = np.zeros(n_loci, dtype=np.float64)
    multi_locus_priors: list[MultiLocusPrior] = [None] * n_loci

    # Pre-bind the per-transcript Locus map per MultiLocus (slow-path only).
    for ml in multi_loci:
        # Build a tiny t_to_locus_idx[N_T_local] once, then route units.
        t_to_locus_idx = _build_t_to_locus_idx(ml, index)
        units_per_locus = partition_units_to_loci(ml, em_data, t_to_locus_idx)

        per_locus_est = tuple(
            estimate_locus_gdna(
                locus, n_obs=int(units_per_locus[j].size),
                region_index=region_index,
                region_arrays=region_arrays,
                payload_arrays=payload_arrays,
                global_densities=global_densities,
                gdna_fl_mean=gdna_fl_mean,
            )
            for j, locus in enumerate(ml.loci)
        )

        ml_prior = assemble_multilocus_prior(ml, per_locus_est)
        multi_locus_priors[ml.multi_locus_id] = ml_prior

        # Dirichlet pseudocount at strength c_base, scaled by π̂.
        alpha_gdna[ml.multi_locus_id] = c_base * ml_prior.pi_gdna
        alpha_rna [ml.multi_locus_id] = c_base * (1.0 - ml_prior.pi_gdna)

    # M5 stub: per-component RNA weight defaults to all-ones.
    # M8b will swap this for build_prior_weight_rna(ml, em_data, nrna_weight=0.0).
    prior_weight_rna = [build_prior_weight_rna(ml, em_data) for ml in multi_loci]

    return PriorTable(
        multi_locus_priors=tuple(multi_locus_priors),
        alpha_gdna=alpha_gdna,
        alpha_rna=alpha_rna,
        prior_weight_rna=prior_weight_rna,
        c_base_value=c_base,
    )
```

`_build_t_to_locus_idx(ml, index)` — file-local helper:

```python
def _build_t_to_locus_idx(ml: MultiLocus, index: TranscriptIndex) -> np.ndarray:
    """Return int8[N_T_local] mapping each ml.transcript_indices[i] to its
    local Locus index in ml.loci (or -1 if a violation, which estimate_locus_gdna
    will then catch)."""
    n_loci_in_ml = len(ml.loci)
    if n_loci_in_ml == 1:
        # All transcripts trivially belong to the only Locus.
        # partition_units_to_loci's fast path doesn't read this, so we
        # return a 1-element marker for cheap dispatch.
        return np.zeros(0, dtype=np.int8)

    # Slow path: each transcript's (ref_id, span) must match exactly one
    # ml.loci entry. We use the transcript table's ref_id + start/end.
    t_idx   = ml.transcript_indices
    t_ref   = index.t_df["ref_id"].values[t_idx]   # int32
    t_start = index.t_df["start"].values[t_idx]    # int64
    out = np.full(t_idx.size, -1, dtype=np.int8)
    for j, loc in enumerate(ml.loci):
        m = (t_ref == loc.ref_id) & (t_start >= loc.start) & (t_start < loc.end)
        out[m] = j
    return out
```

Note: `ml.transcript_indices` are global tids. `t_to_locus_idx[bins]`
in `partition_units_to_loci` indexes into this `int8` array using
**global anchor tids that have already been remapped to local
indices**. Practically, we keep `t_to_locus_idx` indexed by the
**local** position of each anchor transcript in `ml.transcript_indices`
and provide a `searchsorted` lookup. Concretely:

```python
def partition_units_to_loci(ml, em_data, t_to_local_locus):
    if len(ml.loci) == 1:
        return (ml.unit_indices,)
    anchor_t = em_data.locus_t_indices[ml.unit_indices]
    # ml.transcript_indices is sorted ascending by build_multi_loci.
    local_pos = np.searchsorted(ml.transcript_indices, anchor_t)
    # Bounds: anchor_t must lie in ml.transcript_indices (locus invariant).
    bins = t_to_local_locus[local_pos]
    if (bins < 0).any():
        raise RuntimeError(...)
    return tuple(
        ml.unit_indices[bins == j].astype(np.int32, copy=False)
        for j in range(len(ml.loci))
    )
```

This makes the slow path O(`n_units`) with two array gathers and one
`searchsorted`. The `int8` `t_to_local_locus` is `len(ml.transcript_indices)`
long, not `N_T` — both space-efficient and cache-friendly.

### 3.8 Result schemas (final)

```python
# locus_prior.py

C_BASE_DEFAULT = 10.0

FLAG_INTERGENIC_ZERO_LEFF      = 1 << 0
FLAG_INTRON_ZERO_LEFF          = 1 << 1
FLAG_EXON_INTRON_NO_ELIGIBLE   = 1 << 2
FLAG_PI_CLIPPED                = 1 << 3

@dataclass(frozen=True, slots=True)
class LocusGdnaEstimate:
    locus:                 Locus
    n_obs:                 int
    n_gdna_intergenic:     float
    n_gdna_intron:         float
    n_gdna_exon_intron:    float
    n_gdna:                float
    pi_gdna:               float                 # clipped [0, 1]
    rho_loco:              tuple[float, float, float]   # (intergenic, intron, exon-intron)
    leff_loco:             tuple[float, float, float]
    n_eligible_boundaries: int
    fallback_flags:        int                   # bitmask of FLAG_*

@dataclass(frozen=True, slots=True)
class MultiLocusPrior:
    multi_locus_id:        int
    n_obs:                 int
    n_gdna:                float
    n_rna:                 float
    pi_gdna:               float
    per_locus:             tuple[LocusGdnaEstimate, ...]

@dataclass(frozen=True, slots=True)
class PriorTable:
    multi_locus_priors:    tuple[MultiLocusPrior, ...]
    alpha_gdna:            np.ndarray            # float64, (n_loci,)
    alpha_rna:             np.ndarray            # float64, (n_loci,)
    prior_weight_rna:      list[np.ndarray]      # float32, [n_loci][n_components_i]
    c_base_value:          float
```

Indexing: `multi_locus_priors[ml.multi_locus_id]` returns the prior
for that MultiLocus. `alpha_gdna[ml.multi_locus_id]` etc. (these are
the same shape and convention as today's
`compute_locus_priors_from_partitions` outputs, so M8b's pipeline
swap is a one-liner).

---

## 4. Cleanup wins (folded into M6)

While reading the surfaces I noted several improvements that fall
naturally out of this work — each touches the same files so cost is
incremental:

1. **Single source of `L_eff` arithmetic.** `density_loco.py` and
   `locus_prior.py` both use `l_eff_overlap` (already in
   `density_global.py`). Re-export from `density_loco.py` so callers
   never import the *_global module just for `l_eff_overlap`.
2. **`MASK_INTERGENIC` / `MASK_INTRON` constants** are currently
   duplicated as `_MASK_*` private names in `density_global.py`. M6
   introduces `_arrays.py` which needs them; promote them to public
   module attrs in `scan_payload.py` alongside `MASK_N_STATES` so both
   `density_global` and `_arrays` import from the canonical site. Drop
   the private duplicates.
3. **`Locus.span` is already a property** (M5); use it everywhere
   instead of `loc.end - loc.start`. One-line cleanup, removes a
   subtraction at every clipped-region computation that doesn't depend
   on the locus boundary.
4. **`MultiLocus` already exposes `gdna_span` as a precomputed field**
   — the locus-prior pass should NOT recompute it. Audit `assemble_priors`
   to ensure the only span arithmetic is the per-region clip.
5. **`partition_units_to_loci` deletion target.** The legacy
   `compute_locus_priors_from_partitions` reads partitioned log-liks
   to compute γ-posteriors; M6 does NOT need any partition data. M8c
   deletes that function. M6 leaves it in place but adds a
   module-level docstring note: "DEPRECATED post-M8b: use
   `assemble_priors`."
6. **Skip-list for the transcript-table column read.** `_build_t_to_locus_idx`
   reads `index.t_df["start"].values` and `["ref_id"].values` once per
   slow-path MultiLocus. Hoist these into module-level numpy arrays
   on the first call (or reach into the existing `index.t_df` cached
   numpy arrays if any exist — verify during implementation).

---

## 5. Tests

### 5.1 Unit tests

`tests/test_density_loco.py` — ≥ 6 cases:
* (κ=0) ⇒ n_loco / leff_loco.
* (κ=∞ proxy: 1e9) ⇒ rho_global within rtol=1e-6.
* leff_loco=0 ⇒ rho_global.
* n_loco=0, leff_loco>0 ⇒ rho_global · κ / (leff + κ) (correct shrinkage to global at zero local count).
* Monotonicity: doubling n_loco at fixed leff strictly increases output.
* Non-finite input ⇒ raises (degenerate global).

`tests/test_partition_units.py` — ≥ 5 cases:
* Single-Locus MultiLocus ⇒ fast-path returns identity tuple, zero allocation (assert `result[0] is ml.unit_indices`).
* Two-`Locus` paralog MultiLocus, all units anchored to Locus 0 ⇒ result lengths (n, 0).
* Mixed anchors ⇒ correct binning per anchor.
* Anchor outside `ml.transcript_indices` ⇒ `RuntimeError` with the canonical message.
* Slow-path output dtype is `int32` and is a view (or contiguous array) over `ml.unit_indices`.

`tests/test_per_locus_gdna_mass.py` — ≥ 8 cases pinning §2.7.1 end-to-end with hand-built fakes:
* Single INTERGENIC region locus, n=10 fragments at known L_eff ⇒ exact n_gdna, pi_gdna=1.0.
* Single INTRON region, n=0 ⇒ pi_gdna=0.0; rho_loco shrinks to global.
* EXON-only locus with no eligible boundaries ⇒ density falls back to global, n_gdna = rho_global · L_eff_full.
* EXON-only locus with both boundaries eligible, hand-counted u_L/u_R ⇒ exact density and exact predicted count.
* **Short-locus overlap-L_eff regression:** Locus shorter than gdna_fl_mean. Containment formula would give L_eff=0 here; overlap formula gives ≈ gdna_fl_mean. Assert n_gdna > 0. (Master regression for `/memories/repo/calibration-overlap-leff-2026-04.md`.)
* Locus with all three region types ⇒ n_gdna = sum of per-type contributions.
* π̂ > 1 case (capture leakage) ⇒ pi_gdna == 1.0 and `FLAG_PI_CLIPPED` set.
* No-region overlap ⇒ `RuntimeError` with rebuild-index hint.

`tests/test_assemble_priors.py` — ≥ 5 cases:
* MultiLocus with single Locus ⇒ MultiLocusPrior fields match LocusGdnaEstimate.
* Two-Locus MultiLocus ⇒ aggregated n_obs, n_gdna are sums; pi_gdna is the weighted average.
* PriorTable.alpha_gdna == c_base · pi_gdna; PriorTable.alpha_rna == c_base · (1 - pi_gdna).
* PriorTable.prior_weight_rna is a list of float32 ndarrays of length 2 + 1 = 3 (transcript count = 2, +1 gDNA) for a 2-transcript MultiLocus.
* End-to-end sanity: synthetic mini_index + synthetic payload (50/50 gDNA/RNA) ⇒ `mean(pi_gdna)` ≈ 0.5 within ±0.1.

### 5.2 Performance budget

`tests/test_assemble_priors.py::test_perf_5k_multi_loci_under_1s` —
construct 5,000 synthetic MultiLoci (each with 1–3 Loci), run
`assemble_priors`, assert wall < 1.0 s on Apple M1 / M3 dev hardware.
Skipped on CI machines via `pytest.mark.skipif("CI" in os.environ)`.

---

## 6. Exit gate

* All M6 unit tests green (≥ 24 cases per plan).
* Protected suite green.
* `grep -rn 'merged_intervals' src/ tests/` ⇒ 0 hits (already true post-M5).
* `grep -rn '_MASK_INTERGENIC\|_MASK_INTRON' src/rigel/calibration/` ⇒ only the (now public) `MASK_*` constants in `scan_payload.py`.
* No new C++ code; no nanobind ABI change.
* Performance budget satisfied on local hardware.

---

## 7. Commit ordering

Single commit `M6: per-Locus gDNA mass + MultiLocus priors`. The
diff is wide (5 new files, 1 extension) but mechanically scoped: every
new file is independently reviewable and unit-tested.

If the diff is too wide for one PR, split as:

1. `M6-a: density_loco + arrays + region_index_py` (+ `test_density_loco`, `test_region_index_py`).
2. `M6-b: partition_units_to_loci` (+ `test_partition_units`).
3. `M6-c: estimate_locus_gdna + assemble_multilocus_prior` (+ `test_per_locus_gdna_mass`).
4. `M6-d: assemble_priors orchestrator + PriorTable` (+ `test_assemble_priors`).

Each sub-commit is independently testable; each leaves `pipeline.py`
untouched (the new surface is dormant until M8b switches it on).

---

## 8. Risks

| # | Risk | Mitigation |
|---|---|---|
| 1 | `_build_t_to_locus_idx` mismatch with `build_multi_loci` invariants | The function does a strict containment check `start >= loc.start AND start < loc.end`; estimate_locus_gdna asserts no `-1` bins survive. Master test exercises the failure mode. |
| 2 | Per-call `region_arrays.from_region_df` overhead | Materialized once at `assemble_priors` entry; not called per-MultiLocus. |
| 3 | EXON-INTRON branch with `n_eligible == 0` returning `rho=global` and counting full L_eff inflates gDNA on capture-naive loci | Documented; `FLAG_EXON_INTRON_NO_ELIGIBLE` set on the estimate; M9 measures false-positive rate. |
| 4 | `PriorTable.prior_weight_rna` is a list (not a contiguous array) — slow C++ unmarshalling | The `batch_locus_em_partitioned` ABI already accepts `list[ndarray] \| None`; per-locus arrays are tiny (≤ a few KB each); no measurable overhead at 200K loci. |
| 5 | Short-locus overlap-L_eff regression returning if a future refactor flips the formula | Explicit golden test pins behaviour; `l_eff_overlap` has a docstring forbidding the alternative. |
