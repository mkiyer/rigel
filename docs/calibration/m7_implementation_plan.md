# M7 — Pool FL models + `CalibrationResult` v6 schema

**Predecessor:** [calibration_v6_plan.md §M7, §2.9, §2.10](calibration_v6_plan.md).
**Depends on (all green):** M3 (`CalibrationScanPayload` with
`fl_hist`), M4 (`GlobalDensityTable`), M6 (`PriorTable`,
`MultiLocusPrior`, `LocusGdnaEstimate`).

**Goal:** publish the final calibration result. Convert the M3 payload
+ scan-trained `frag_length_models` into the four pool FL models
required by downstream scoring/EM, then assemble the immutable v6
`CalibrationResult` carrying global densities, pool models, the
M6 prior table, and per-MultiLocus / per-Locus diagnostic dataframes.

---

## 1. Surface contract

After M7 lands, calibration callers can write:

```python
from rigel.calibration import (
    CalibrationResult,
    PoolFLModels,
    compute_pool_fl_models,
    build_calibration_result,
)

pool = compute_pool_fl_models(
    payload=payload,
    scan_trained=frag_length_models,        # from BamScanner
    global_density_block=global_densities,  # for n_pool_annotation_gap
)

calibration = build_calibration_result(
    payload=payload,
    global_densities=global_densities,
    pool=pool,
    prior_table=prior_table,
    multi_loci=multi_loci,
)

# Immutable. To attach a new prior_table later (e.g., refit), use:
calibration2 = calibration.with_priors(new_prior_table)
```

The legacy `_simple.CalibrationResult` (v1) coexists. M7 publishes the
v6 schema in **a new file** `_result_v6.py` to avoid collision. M8c
deletes `_simple.py` and renames `_result_v6.py` → `result.py`.

---

## 2. Module layout

Two new files, two extensions:

```
src/rigel/calibration/
├── _fl_pool.py           (new) PoolFLModels, compute_pool_fl_models
├── _result_v6.py         (new) CalibrationResult, build_calibration_result,
│                                _build_multi_locus_prior_df,
│                                _build_per_locus_gdna_df
└── __init__.py           (extend) re-export PoolFLModels, CalibrationResult,
                                    compute_pool_fl_models, build_calibration_result
```

The legacy `_result.py` is left untouched; M8c deletes it.

---

## 3. Pool FL models

### 3.1 The pool taxonomy

Per [§2.9](calibration_v6_plan.md#29-pool-fl-models). Four pools:

| Pool | Source mask(s) | Source model |
|---|---|---|
| `gdna` | `fl_hist[2] ∪ fl_hist[3] ∪ fl_hist[4]` (all unspliced gDNA-bearing) | derived from `payload.fl_hist` |
| `rna_nrna_unspliced` | scan-trained nRNA model | `scan_trained.nrna_model` (passthrough) |
| `rna_mrna_spliced` | scan-trained mRNA model | `scan_trained.mrna_model` (passthrough) |
| `global` | scan-trained generic FL model | `scan_trained.global_model` (passthrough) |

**Why are RNA / global passthrough?** The scanner already trains
high-quality FL models from unique mappers using the full BAM, with
proper genomic-vs-spliced span normalization. Re-deriving them from
`fl_hist` would (a) discard the scanner's per-strand stratification
and (b) fall into the **spliced-genomic-span trap**: `fl_hist[1]`
(EXON-only mask) contains predominantly *spliced* mRNA fragments
whose **genomic** span equals the genomic distance between the outer
read ends — NOT the fragment length. Histogramming this column as a
fragment-length distribution is the textbook way to overestimate FL
by 100s of bp on real data. M7 explicitly rejects mask 1 as a FL
source.

The gDNA pool is the *only* one we re-derive from `fl_hist` because
(a) the scanner has no way to identify gDNA fragments at scan time
(SRD calibration runs *after* the scan), and (b) gDNA fragments are
unspliced ⇒ genomic span == FL ⇒ `fl_hist[2,3,4]` is a valid FL
distribution.

### 3.2 Quality classifier

Per [§2.9.2](calibration_v6_plan.md#292-quality-classifier). Each pool
gets one of three labels:

| Label | When | Action |
|---|---|---|
| `good` | `n_pool >= 5000` AND model is well-conditioned | use as-is |
| `weak` | `200 <= n_pool < 5000` (gDNA only — scanner pools never hit this branch) | EB-shrink the empirical FL toward `pool.global` with ESS=1000 |
| `fallback` | `n_pool < 200` (gDNA only) | identity copy of `pool.global`, flagged |

Constants:

```python
POOL_QUALITY_GOOD_THRESHOLD = 5000
POOL_QUALITY_WEAK_THRESHOLD = 200
POOL_EB_PRIOR_ESS           = 1000.0
```

### 3.3 `PoolFLModels` schema

```python
# src/rigel/calibration/_fl_pool.py

@dataclass(frozen=True, slots=True)
class PoolFLModels:
    # The four FL distributions, each a FragmentLengthModel.
    gdna:                      FragmentLengthModel
    rna_nrna_unspliced:        FragmentLengthModel
    rna_mrna_spliced:          FragmentLengthModel
    global_:                   FragmentLengthModel    # trailing _ to avoid keyword

    # Quality + provenance.
    gdna_quality:              Literal["good", "weak", "fallback"]
    gdna_n_fragments:          int
    gdna_fl_mean:              float
    gdna_eb_ess:               float          # 0.0 if not shrunk
    gdna_used_global_fallback: bool

    # Per-pool annotation gap counters: mask-N counts that landed in
    # `fl_hist[mask]` but did not contribute to any pool. Useful for
    # diagnosing "why is gDNA n smaller than I expected?".
    n_pool_annotation_gap:     dict[str, int]

    def to_summary_dict(self) -> dict[str, object]: ...
```

### 3.4 `compute_pool_fl_models` algorithm

```python
def compute_pool_fl_models(
    payload:               CalibrationScanPayload,
    scan_trained:          FragLengthModels,            # from BamScanner
    *,
    global_density_block:  GlobalDensityTable,           # for diagnostics only
) -> PoolFLModels:
    # 1. Pool gDNA fl_hist over masks 2 (INTRON_ONLY), 3 (EXON+INTRON),
    #    4 (INTERGENIC_ONLY). All three are unspliced ⇒ genomic span == FL.
    gdna_hist = (
        payload.fl_hist[mask.INTRON_ONLY] +
        payload.fl_hist[mask.EXON_INTRON] +
        payload.fl_hist[mask.INTERGENIC_ONLY]
    )                                                    # int64, (1024,)
    n_gdna = int(gdna_hist.sum())

    if n_gdna >= POOL_QUALITY_GOOD_THRESHOLD:
        gdna_model = FragmentLengthModel.from_histogram(gdna_hist)
        quality, eb_ess, fallback = "good", 0.0, False

    elif n_gdna >= POOL_QUALITY_WEAK_THRESHOLD:
        # EB-shrink toward global with ESS=1000.
        gdna_model = FragmentLengthModel.from_histogram(gdna_hist).finalize_with_prior(
            prior=scan_trained.global_model,
            prior_ess=POOL_EB_PRIOR_ESS,
        )
        quality, eb_ess, fallback = "weak", POOL_EB_PRIOR_ESS, False

    else:
        # Below 200 fragments: copy the scan-trained global model.
        gdna_model = scan_trained.global_model.copy()
        quality, eb_ess, fallback = "fallback", 0.0, True

    # 2. Annotation-gap diagnostics — fragments that landed in fl_hist
    #    but did not contribute to gdna pool (i.e. mask 0, 1, 5, 6, 7).
    n_pool_annotation_gap = {
        f"mask_{m}": int(payload.fl_hist[m].sum())
        for m in (0, 1, 5, 6, 7)
    }

    return PoolFLModels(
        gdna=gdna_model,
        rna_nrna_unspliced=scan_trained.nrna_model,         # passthrough
        rna_mrna_spliced=scan_trained.mrna_model,           # passthrough
        global_=scan_trained.global_model,                  # passthrough
        gdna_quality=quality,
        gdna_n_fragments=n_gdna,
        gdna_fl_mean=float(gdna_model.mean),
        gdna_eb_ess=eb_ess,
        gdna_used_global_fallback=fallback,
        n_pool_annotation_gap=n_pool_annotation_gap,
    )
```

`FragmentLengthModel.finalize_with_prior` may not exist yet — verify
in `src/rigel/frag_length_model.py` during implementation. If absent,
add it as a tiny method:

```python
def finalize_with_prior(self, prior: "FragmentLengthModel", prior_ess: float) -> "FragmentLengthModel":
    """EB shrinkage: posterior ∝ self.counts + prior_ess · prior.pmf."""
    n_self = float(self.counts.sum())
    if n_self <= 0.0:
        return prior.copy()
    pseudo_counts = prior_ess * prior.pmf       # (1024,) float64
    posterior_counts = self.counts.astype(np.float64) + pseudo_counts
    return FragmentLengthModel.from_pmf(posterior_counts / posterior_counts.sum())
```

### 3.5 Mask name constants

Define alongside `MASK_N_STATES`:

```python
# src/rigel/calibration/scan_payload.py
class mask:
    """Symbolic names for the 8 fragment-class bitmasks (3-bit code: EXON|INTRON|INTERGENIC)."""
    NONE             = 0b000   # 0 — fragment hit no annotated region
    EXON_ONLY        = 0b001   # 1 — spliced mRNA-bearing
    INTRON_ONLY      = 0b010   # 2 — unspliced, every block intronic
    EXON_INTRON      = 0b011   # 3 — unspliced, mixed exon/intron blocks
    INTERGENIC_ONLY  = 0b100   # 4 — every block intergenic
    EXON_INTERGENIC  = 0b101   # 5 — readthrough into intergenic
    INTRON_INTERGENIC= 0b110   # 6 — intron + intergenic
    ALL              = 0b111   # 7 — exon + intron + intergenic
```

This is a **cleanup win**: today the three `_MASK_*` constants in
`density_global.py` are private with no symbolic friends for the
remaining five masks. M7 promotes them to a single public class
that both M4's `density_global` and M7's `_fl_pool` import from.
Drop the `_MASK_*` private duplicates.

---

## 4. `CalibrationResult` v6 schema

### 4.1 The dataclass

```python
# src/rigel/calibration/_result_v6.py

@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Immutable result of v6 calibration. Carries everything quant_from_buffer needs."""

    # --- Block 1: global densities (M4) ---
    global_densities:       GlobalDensityTable

    # --- Block 2: pool FL models (M7) ---
    pool:                   PoolFLModels

    # --- Block 3: per-MultiLocus priors (M6) ---
    prior_table:            PriorTable

    # --- Block 4: provenance ---
    payload_summary:        dict[str, int]    # n_observed, n_excluded_*, n_unobserved
    n_multi_loci:           int

    # --- Block 5: derived diagnostic dataframes (lazy in spirit, eager in code) ---
    multi_locus_prior_df:   pd.DataFrame      # one row per MultiLocus
    per_locus_gdna_df:      pd.DataFrame      # one row per Locus

    # --- Convenience accessors ---
    @property
    def alpha_gdna(self) -> np.ndarray: return self.prior_table.alpha_gdna
    @property
    def alpha_rna(self) -> np.ndarray: return self.prior_table.alpha_rna
    @property
    def prior_weight_rna(self) -> list[np.ndarray]: return self.prior_table.prior_weight_rna
    @property
    def gdna_fl_mean(self) -> float: return self.pool.gdna_fl_mean

    # --- Mutator-style helpers (return new instance; frozen-safe) ---
    def with_priors(self, prior_table: PriorTable) -> "CalibrationResult":
        """Return a copy with prior_table replaced (for refit / posthoc adjustment)."""
        return dataclasses.replace(
            self,
            prior_table=prior_table,
            multi_locus_prior_df=_build_multi_locus_prior_df(prior_table.multi_locus_priors),
            # per_locus_gdna_df is unchanged (depends only on per_locus estimates,
            # which prior_table.multi_locus_priors carries forward).
            per_locus_gdna_df=_build_per_locus_gdna_df(prior_table.multi_locus_priors),
        )

    def to_summary_dict(self) -> dict[str, object]:
        """Flatten to JSON-serializable dict for summary.json."""
        return {
            "global_densities": self.global_densities.to_summary_dict(),
            "pool":             self.pool.to_summary_dict(),
            "n_multi_loci":     self.n_multi_loci,
            "payload_summary":  self.payload_summary,
            "c_base":           self.prior_table.c_base_value,
            "mean_pi_gdna":     float(np.mean([
                ml.pi_gdna for ml in self.prior_table.multi_locus_priors
            ])) if self.prior_table.multi_locus_priors else 0.0,
        }
```

**Explicit non-fields** (per [§2.10.4](calibration_v6_plan.md#2104-explicit-non-fields)):
* No `version` field. The class is the version.
* No `active` / `enabled` flag. Existence implies use.
* No mutable cache slots. Re-call `with_priors(...)` to update.
* No `payload` reference. The result is decoupled from the raw scan
  data; downstream code reads diagnostic dataframes, not `fl_hist`.

### 4.2 `build_calibration_result` orchestrator

```python
def build_calibration_result(
    *,
    payload:           CalibrationScanPayload,
    global_densities:  GlobalDensityTable,
    pool:              PoolFLModels,
    prior_table:       PriorTable,
    multi_loci:        list[MultiLocus],
) -> CalibrationResult:
    payload_summary = {
        "n_observed":           payload.n_observed,
        "n_excluded_multimap":  payload.n_excluded_multimap,
        "n_excluded_chimera":   payload.n_excluded_chimera,
        "n_excluded_artifact":  payload.n_excluded_artifact,
        "n_unobserved":         payload.n_unobserved,
        "n_unannotated_ref":    payload.n_unannotated_ref,
    }

    return CalibrationResult(
        global_densities      = global_densities,
        pool                  = pool,
        prior_table           = prior_table,
        payload_summary       = payload_summary,
        n_multi_loci          = len(multi_loci),
        multi_locus_prior_df  = _build_multi_locus_prior_df(prior_table.multi_locus_priors),
        per_locus_gdna_df     = _build_per_locus_gdna_df(prior_table.multi_locus_priors),
    )
```

### 4.3 Diagnostic dataframes

```python
def _build_multi_locus_prior_df(
    mlps: tuple[MultiLocusPrior, ...]
) -> pd.DataFrame:
    return pd.DataFrame({
        "multi_locus_id": [m.multi_locus_id for m in mlps],
        "n_obs":          [m.n_obs          for m in mlps],
        "n_gdna":         [m.n_gdna         for m in mlps],
        "n_rna":          [m.n_rna          for m in mlps],
        "pi_gdna":        [m.pi_gdna        for m in mlps],
        "n_loci":         [len(m.per_locus) for m in mlps],
    })

def _build_per_locus_gdna_df(
    mlps: tuple[MultiLocusPrior, ...]
) -> pd.DataFrame:
    rows = []
    for ml in mlps:
        for est in ml.per_locus:
            rows.append({
                "multi_locus_id":        ml.multi_locus_id,
                "ref":                   est.locus.ref,
                "start":                 est.locus.start,
                "end":                   est.locus.end,
                "span":                  est.locus.span,
                "n_obs":                 est.n_obs,
                "n_gdna":                est.n_gdna,
                "n_gdna_intergenic":     est.n_gdna_intergenic,
                "n_gdna_intron":         est.n_gdna_intron,
                "n_gdna_exon_intron":    est.n_gdna_exon_intron,
                "pi_gdna":               est.pi_gdna,
                "n_eligible_boundaries": est.n_eligible_boundaries,
                "fallback_flags":        est.fallback_flags,
            })
    return pd.DataFrame(rows)
```

Both dataframes have **locked column orders** — tests pin them.

---

## 5. Cleanup wins (folded into M7)

1. **Symbolic mask names** (§3.5): one public `mask` class replaces
   four private `_MASK_*` ints across `density_global.py` and
   `_fl_pool.py`. Improves readability site-wide; no behaviour
   change.
2. **`with_priors` via `dataclasses.replace`** (idiomatic, not a
   custom mutator): preserves frozen contract; one line at every call
   site that today copies fields manually.
3. **Decouple result from payload.** Today's `_simple.CalibrationResult`
   stores a reference to the raw payload — invites diagnostic code
   to reach back into `fl_hist` and recompute things. M7's
   `CalibrationResult` carries diagnostic dataframes only. M8c
   deletion of `_simple.py` then mechanical.
3. **Eager dataframes vs lazy properties.** Both diagnostic frames are
   built eagerly at `build_calibration_result` time. They are at most
   `O(n_multi_loci)` and `O(n_loci)` rows respectively, ≤ a few
   hundred KB on real indexes. Eager construction removes a class of
   "frozen dataclass + lazy cache" anti-patterns.
4. **No `version` / `active` flags.** v6 is the only schema as of
   M8c; the legacy v1 is deleted, not deprecated-then-deleted. This
   removes ~10 conditional branches in pipeline.py.
5. **Unify `to_summary_dict` discipline.** `GlobalDensityTable`
   already has one; M7 adds the same to `PoolFLModels` and
   `CalibrationResult`. `summary.json` writer becomes a single
   `result.to_summary_dict()` call.

---

## 6. Tests

### 6.1 Pool-FL tests

`tests/test_pool_fl.py` — ≥ 8 cases:
* `good` branch: pool with 10,000 fragments → quality="good", `gdna_eb_ess == 0.0`, model exactly matches `from_histogram` output.
* `weak` branch: pool with 1,000 fragments → quality="weak", `gdna_eb_ess == 1000.0`, model is a weighted blend (verify mean lies between empirical and global means).
* `fallback` branch: pool with 50 fragments → quality="fallback", `gdna_used_global_fallback == True`, model PMF is bit-identical to `scan_trained.global_model`.
* `fallback` branch with empty pool: `n_gdna == 0` → fallback, no division-by-zero.
* RNA passthrough integrity: `rna_nrna_unspliced is scan_trained.nrna_model` (object identity), same for `rna_mrna_spliced`, `global_`.
* Spliced-genomic-span trap regression: assert `fl_hist[1]` (EXON-only) does NOT contribute to `gdna_n_fragments`. Construct a payload where mask-1 has 100,000 entries with mean genomic span 5000 bp; assert `pool.gdna_fl_mean` is unaffected.
* Annotation gap counters: build a payload with known per-mask totals; assert `n_pool_annotation_gap` matches mask-{0,1,5,6,7} sums exactly.
* `to_summary_dict()` round-trips through `json.dumps` without error.

### 6.2 `CalibrationResult` tests

`tests/test_calibration_result_v6.py` — ≥ 10 cases:
* Construction round-trip: `build_calibration_result` returns a frozen instance; `dataclasses.is_dataclass` and `result.__dataclass_fields__` checks.
* Frozen contract: `result.gdna_fl_mean = 0.0` raises `FrozenInstanceError`.
* `with_priors` returns a new instance; original is unmodified (test field-by-field).
* `with_priors` rebuilds `multi_locus_prior_df` and `per_locus_gdna_df` (verify rows match the new prior_table).
* Convenience accessors: `result.alpha_gdna is result.prior_table.alpha_gdna` (object identity, no copy).
* `multi_locus_prior_df` schema lock: exact column order `["multi_locus_id", "n_obs", "n_gdna", "n_rna", "pi_gdna", "n_loci"]`.
* `per_locus_gdna_df` schema lock: 14 columns in locked order (per §4.3).
* Empty MultiLoci list: `n_multi_loci == 0`, both dataframes are empty with correct columns.
* `to_summary_dict()` is JSON-serializable (test via `json.dumps(default=str)`).
* `payload_summary` reflects all six payload counters.

### 6.3 Integration

`tests/test_m7_integration.py` — 1 case end-to-end:
* Build a synthetic `CalibrationScanPayload` + `region_df` + `multi_loci` of moderate size (≈20 MultiLoci, ≈40 Loci, ≈5K fragments). Run the full M4→M6→M7 chain. Assert: `result.n_multi_loci == 20`, `len(result.per_locus_gdna_df) == 40`, `0.0 <= result.alpha_gdna.mean() <= result.prior_table.c_base_value`.

### 6.4 Protected suite

* All tests in `tests/test_calibration_simple.py` continue to pass (legacy v1 path untouched).
* All M3/M4/M5/M6 tests continue to pass.

---

## 7. Exit gate

* All M7 unit + integration tests green (≥ 19 cases).
* Protected suite green.
* `grep -rn '_MASK_INTERGENIC\|_MASK_INTRON\|_MASK_EXON' src/rigel/` ⇒ 0 hits (all migrated to `mask.*`).
* `grep -rn 'from .calibration import CalibrationResult' src/rigel/` ⇒ still resolves to v1 (the M8b switch comes later); v6 is imported only by tests.
* No new C++ code; no nanobind ABI change.

---

## 8. Commit ordering

Single commit `M7: pool FL models + CalibrationResult v6 schema`.

If split is preferred:

1. `M7-a: mask names + PoolFLModels + compute_pool_fl_models` (+ `test_pool_fl`).
2. `M7-b: CalibrationResult v6 + diagnostic dataframes` (+ `test_calibration_result_v6`).
3. `M7-c: M4→M6→M7 integration test`.

---

## 9. Risks

| # | Risk | Mitigation |
|---|---|---|
| 1 | `FragmentLengthModel.finalize_with_prior` may not exist; signature mismatch | First implementation step is to grep for the symbol; if absent, add per §3.4 with a 5-line method. Tested in isolation before `_fl_pool`. |
| 2 | `scan_trained.{nrna,mrna,global}_model` API drift across the FragLengthModels surface | Pin via `tests/test_pool_fl.py::test_passthrough_identity` — the test fails loudly on any rename. |
| 3 | Two `CalibrationResult` symbols in flight (v1 in `_simple.py`, v6 in `_result_v6.py`) confuse callers | M7 does NOT export the v6 name from `__init__.py` under the bare `CalibrationResult` alias. Tests import from `rigel.calibration._result_v6`. M8c renames + re-exports. |
| 4 | Eager dataframe construction inflates memory at scale (millions of Loci) | Per-Locus row is 14 fields × ~100B = ~1.4 KB. 1M Loci ⇒ ~1.4 GB; for indexes that large the dataframes become a real cost. **Defer to M9 if profiling flags it** — easy retrofit via lazy properties + a `materialize()` method. Not a M7 blocker (typical run: ≤ 500K Loci ⇒ ~700 MB peak, comparable to existing index footprint). |
| 5 | `with_priors` rebuilds `per_locus_gdna_df` even though it depends only on `prior_table.multi_locus_priors[i].per_locus` (which `with_priors` accepts wholesale, so it's stale or not) | The new `prior_table` carries new `MultiLocusPrior` instances ⇒ new `per_locus` ⇒ correct rebuild. Test pins the behaviour. |
