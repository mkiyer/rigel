# M7 — Fragment-Length Models + CalibrationResult v6 (Revised v2)

**Status**: revised plan (v2).  Supersedes the v1 draft and the M7
commit `b1badaf`.

**Scope**: pure-Python; no C++ ABI change in M7.  Calibration becomes
the sole owner of the finalized FL distribution surface.

---

## 1. The One Policy

> **Three FL distributions exist in v6: `global`, `rna`, `gdna`.
> They are owned by exactly one module: `rigel.calibration.fl`.
> They are built by exactly one function: `build_fl_models(...)`.
> They are shrunk by exactly one rule: Empirical-Bayes Dirichlet
> shrinkage of `rna` and `gdna` toward `global`, with a single ESS.**

Everything else in this document follows from that sentence.

---

## 2. Producers vs Owner — the data flow

```
                    ┌─────────────────────────────┐
   BAM scanner ─►   │  raw int64 count vectors    │
   (C++)            │  • global_counts            │
                    │  • rna_counts (SPLICED)     │
                    └─────────────────────────────┘
                                  │
   Calibration accumulator ─►     │
   (C++, gates on UNSPLICED ∩     │  • gdna_counts
    {intron, exon-intron,         │
     intergenic})                 │
                                  ▼
                    ┌─────────────────────────────┐
                    │  rigel.calibration.fl       │
                    │  build_fl_models(...)       │  ← single owner,
                    │                             │    single policy
                    │  global = from raw          │
                    │  rna    = EB-shrunk(global) │
                    │  gdna   = EB-shrunk(global) │
                    └─────────────────────────────┘
                                  │
                                  ▼
                          FLModels (immutable)
                                  │
                                  ▼
                       CalibrationResult (M7)
                                  │
                                  ▼
                            EM solver (M8)
```

**Producers** of raw histograms emit `np.ndarray[int64]` count vectors
and stop there.  They never call `finalize()`, never apply a prior,
never wrap counts in a `FragmentLengthModel`.

**One owner** (`rigel.calibration.fl`) takes the raw counts, applies
the one policy, and returns three finalized `FragmentLengthModel`
instances inside an immutable `FLModels` struct.

**Consumers** (the EM scorer, downstream reporting, `summary.json`)
read finalized models from `FLModels` only.  They never see raw
counts, never re-finalize, never know that "mask 2 + 3 + 4" exists.

---

## 3. Implications for existing code

### 3.1 `FragmentLengthModels` (the legacy class)

The `FragmentLengthModels` container in `frag_length_model.py` is a
**legacy artefact of SRD-v1**:

* It owns three roles today: (1) raw histogram accumulator at scan
  time, (2) post-scan finaliser via `build_scoring_models()` +
  `finalize(prior_ess=...)`, (3) the container the scorer reads
  `rna_model` / `gdna_model` from.
* Roles (2) and (3) are duplicates of v6's `FLModels`.

**M7 policy**: the v6 path does **not** call `build_scoring_models()`
or `finalize(prior_ess=...)` on `FragmentLengthModels`.  We treat the
scanner-side instance purely as a typed view onto raw count vectors
(`global_model.counts`, `category_models[SPLICED_ANNOT].counts`).
The legacy `rna_model` / `gdna_model` slots are never read on the v6
path.

**M8c**: `FragmentLengthModels` is deleted.  The scanner is updated to
return raw count arrays directly (a one-line C++ change at the
finalize site).  Until then, we tolerate the duplicate scaffolding —
calling `.counts` on the scanner-trained sub-models is the cleanest
way to get raw histograms without touching C++.

### 3.2 `_simple.py` (SRD-v1 calibration)

Continues to use `FragmentLengthModels` exactly as today.  Untouched
in M7.  Deleted in M8c.

### 3.3 `pipeline.py`

Continues to construct `FragmentLengthModels` for SRD-v1 today.  In
M8 it will additionally surface the v6 `CalibrationResult.fl_models`
and wire those into the scorer.  Untouched in M7.

### 3.4 The 3-bit region mask

The mask is an internal abstraction that lives in:
* C++ accumulator (`accumulator.{h,cpp}`).
* `scan_payload.py` `MASK_*` constants.
* `density_global.py` (legitimately uses per-mask region counts).

It is **never** named in:
* Public field names.
* `summary.json` keys.
* Diagnostic dataframe columns.

Where M7 needs to sum mask rows to get a histogram (e.g. gDNA =
`fl_hist[2] + [3] + [4]`), the sum lives **inside** `build_fl_models`
as a private helper.  No other module touches that arithmetic.

### 3.5 EXON_ONLY fragments

Counted, not modelled.  EXON_ONLY = mask 0b001.  It contributes to
`n_observed` and to `Diagnostics.n_exon_only`.  It does **not** go into
any FL histogram on the v6 path (the SPLICED ones do via the scanner
training site; the unspliced-exon-only ones are a small, low-signal
bucket that has no calibration role).

---

## 4. Public API

### 4.1 `rigel.calibration.fl` (the FL submodule)

```python
# rigel/calibration/fl.py

POOL_EB_PRIOR_ESS              = 1000.0
POOL_QUALITY_GOOD_THRESHOLD    = 5000   # n ≥ this → no shrinkage
POOL_QUALITY_WEAK_THRESHOLD    = 200    # n < this → identity-fallback to global

Quality = Literal["good", "weak", "fallback"]


@dataclass(frozen=True, slots=True)
class FLModels:
    """The three finalized FL distributions.  Sole product of the
    calibration FL pipeline; sole input the EM scorer needs."""

    global_: FragmentLengthModel    # unconditional anchor (no prior)
    rna:     FragmentLengthModel    # SPLICED, EB-shrunk to global
    gdna:    FragmentLengthModel    # UNSPLICED ∩ {intron,intergenic,exon-intron}

    rna_quality:  Quality
    gdna_quality: Quality
    n_global:     int
    n_rna:        int
    n_gdna:       int

    def to_summary_dict(self) -> dict[str, object]: ...


def build_fl_models(
    *,
    global_counts: np.ndarray,    # int64 (max_size+1,) — every observed fragment
    rna_counts:    np.ndarray,    # int64 (max_size+1,) — SPLICED only
    gdna_counts:   np.ndarray,    # int64 (max_size+1,) — gDNA pool only
    max_size:      int,
    prior_ess:     float = POOL_EB_PRIOR_ESS,
    good_threshold: int = POOL_QUALITY_GOOD_THRESHOLD,
    weak_threshold: int = POOL_QUALITY_WEAK_THRESHOLD,
) -> FLModels:
    """Build the three FL distributions under one EB policy."""
```

The function is the single owner.  It takes raw count vectors as
function arguments — not payload objects, not scanner objects — so
the data flow is explicit and the function is trivially testable.

### 4.2 `rigel.calibration.scan_payload` — sourcing raw counts

The producers' raw counts reach `build_fl_models` via two thin
helpers (still in the calibration module — the FL submodule does
not import the scanner):

```python
# rigel/calibration/_fl_sources.py  (private)

def extract_global_counts(scan_trained: FragmentLengthModels) -> np.ndarray:
    """Raw global FL histogram from the scanner."""
    return scan_trained.global_model.counts.astype(np.int64, copy=False)

def extract_rna_counts(scan_trained: FragmentLengthModels) -> np.ndarray:
    """Raw SPLICED-ANNOT FL histogram from the scanner."""
    from rigel.splice import SpliceType
    return scan_trained.category_models[SpliceType.SPLICED_ANNOT].counts.astype(
        np.int64, copy=False)

def extract_gdna_counts(payload: CalibrationScanPayload) -> np.ndarray:
    """Raw gDNA pool FL histogram from the calibration accumulator.

    Sum of mask rows {INTRON_ONLY, EXON_INTRON, INTERGENIC_ONLY}.
    Mask 0b001 (EXON_ONLY) contains spliced-fragment GENOMIC SPANS,
    not fragment lengths — explicitly excluded.
    """
    h = payload.fl_hist
    return h[MASK_INTRON] + (h[MASK_EXON | MASK_INTRON]) + h[MASK_INTERGENIC]
```

These three functions are the **only** code that knows where raw
counts live.  After M7+a (C++ collapse), only their bodies change;
`build_fl_models` is unaffected.

### 4.3 `Diagnostics` — named breakdown of `n_observed`

A typed struct replaces the prior `dict[str, int]`:

```python
@dataclass(frozen=True, slots=True)
class Diagnostics:
    """Named decomposition of payload.n_observed (sums to n_observed)."""
    n_exon_only:              int
    n_intron_only:            int
    n_intergenic_only:        int
    n_exon_intron:            int
    n_exon_intergenic:        int
    n_intron_intergenic:      int
    n_exon_intron_intergenic: int
    n_unannotated:            int

    def total(self) -> int:
        return (self.n_exon_only + self.n_intron_only +
                self.n_intergenic_only + self.n_exon_intron +
                self.n_exon_intergenic + self.n_intron_intergenic +
                self.n_exon_intron_intergenic + self.n_unannotated)

    def to_summary_dict(self) -> dict[str, int]: ...

    @classmethod
    def from_payload(cls, payload: CalibrationScanPayload) -> "Diagnostics": ...
```

### 4.4 `rigel.calibration._result_v6.CalibrationResult`

```python
@dataclass(frozen=True, slots=True)
class CalibrationResult:
    global_densities: GlobalDensityTable    # M4
    fl_models:        FLModels              # M7  (single FL surface)
    prior_table:      PriorTable            # M6

    diagnostics:   Diagnostics              # named breakdown of n_observed
    n_multi_loci:  int

    multi_locus_prior_df: pd.DataFrame
    per_locus_gdna_df:    pd.DataFrame

    # Convenience zero-copy aliases
    @property
    def alpha_gdna(self) -> np.ndarray:        ...
    @property
    def alpha_rna(self) -> np.ndarray:         ...
    @property
    def prior_weight_rna(self) -> list[np.ndarray]: ...

    @property
    def gdna_fl_mean(self) -> float:   return self.fl_models.gdna.mean
    @property
    def rna_fl_mean(self) -> float:    return self.fl_models.rna.mean
    @property
    def global_fl_mean(self) -> float: return self.fl_models.global_.mean

    def with_priors(self, prior_table: PriorTable) -> "CalibrationResult": ...
    def to_summary_dict(self) -> dict[str, object]: ...
```

### 4.5 Orchestrator

```python
def build_calibration_result(
    *,
    payload:          CalibrationScanPayload,
    scan_trained:     FragmentLengthModels,    # raw-counts source only
    global_densities: GlobalDensityTable,
    prior_table:      PriorTable,
    fl_prior_ess:     float = POOL_EB_PRIOR_ESS,
) -> CalibrationResult:
    fl_models = build_fl_models(
        global_counts = extract_global_counts(scan_trained),
        rna_counts    = extract_rna_counts(scan_trained),
        gdna_counts   = extract_gdna_counts(payload),
        max_size      = scan_trained.max_size,
        prior_ess     = fl_prior_ess,
    )
    diagnostics = Diagnostics.from_payload(payload)
    return CalibrationResult(
        global_densities = global_densities,
        fl_models        = fl_models,
        prior_table      = prior_table,
        diagnostics      = diagnostics,
        n_multi_loci     = len(prior_table.multi_locus_priors),
        multi_locus_prior_df = build_multi_locus_prior_df(
            prior_table.multi_locus_priors),
        per_locus_gdna_df    = build_per_locus_gdna_df(
            prior_table.multi_locus_priors),
    )
```

That is the full v6 calibration entry point — six explicit lines and
zero ambiguity about who owns what.

---

## 5. The single EB policy (concrete)

```python
def _finalize(counts: np.ndarray, max_size: int,
              global_counts: np.ndarray | None,
              prior_ess: float | None) -> FragmentLengthModel:
    """One primitive, used three times in build_fl_models."""
    fl = FragmentLengthModel(max_size=max_size)
    fl.counts[: counts.size] = counts.astype(np.float64, copy=False)
    fl._total_weight = float(fl.counts.sum())
    if global_counts is None:
        fl.finalize()                          # global anchor
    else:
        fl.finalize(prior_counts=global_counts.astype(np.float64),
                    prior_ess=prior_ess)
    return fl


def build_fl_models(*, global_counts, rna_counts, gdna_counts,
                    max_size, prior_ess, good_threshold, weak_threshold) -> FLModels:
    n_g, n_r, n_d = int(global_counts.sum()), int(rna_counts.sum()), int(gdna_counts.sum())

    global_ = _finalize(global_counts, max_size, None, None)

    rna,  rna_q  = _classify_and_build(rna_counts,  global_counts, max_size,
                                       prior_ess, good_threshold, weak_threshold,
                                       global_)
    gdna, gdna_q = _classify_and_build(gdna_counts, global_counts, max_size,
                                       prior_ess, good_threshold, weak_threshold,
                                       global_)
    return FLModels(global_=global_, rna=rna, gdna=gdna,
                    rna_quality=rna_q, gdna_quality=gdna_q,
                    n_global=n_g, n_rna=n_r, n_gdna=n_d)


def _classify_and_build(counts, global_counts, max_size, prior_ess,
                        good_threshold, weak_threshold,
                        global_model) -> tuple[FragmentLengthModel, Quality]:
    n = int(counts.sum())
    if n >= good_threshold:                     # pure empirical
        return _finalize(counts, max_size, None, None), "good"
    if n >= weak_threshold:                     # EB-shrunk
        return _finalize(counts, max_size, global_counts, prior_ess), "weak"
    return global_model, "fallback"             # identity-share with global
```

Three branches.  One primitive.  Symmetric for RNA and gDNA.  No
duplicate logic anywhere.

---

## 6. Migration from commit `b1badaf` (one commit)

1. **Rename module**: `_fl_pool.py` → `fl.py` (public submodule, no
   leading underscore — it's the user-facing FL surface for v6).
2. **Rename type**: `PoolFLModels` → `FLModels`; rename fields per
   §4.1 (`gdna` / `rna` / `global_`).
3. **Switch RNA source**: previously a passthrough of the scanner's
   already-finalised `rna_model`.  Now reads raw counts via
   `extract_rna_counts(scan_trained)` and applies the M7 EB policy.
   The scanner-side `rna_model.finalize(...)` becomes unused on the
   v6 path (see §3.1).
4. **Add** `Diagnostics` (typed struct).  Delete `n_pool_annotation_gap`
   (`dict[str, int]`).  Delete `gdna_eb_ess` and
   `gdna_used_global_fallback` from `FLModels` — both are derivable
   from `gdna_quality`.
5. **Update** `_result_v6.py`:
   * field rename `pool` → `fl_models`
   * field type change `payload_summary: dict` →
     `diagnostics: Diagnostics`
   * add `rna_fl_mean`, `global_fl_mean` properties
6. **Add** `_fl_sources.py` (three private extractor functions, §4.2).
7. **Update** `__init__.py` exports:
   * remove `PoolFLModels`, `compute_pool_fl_models`
   * add `FLModels`, `build_fl_models`, `Diagnostics`,
     `POOL_EB_PRIOR_ESS`, `POOL_QUALITY_GOOD_THRESHOLD`,
     `POOL_QUALITY_WEAK_THRESHOLD`, `Quality`
   * `CalibrationResultV6` alias kept until M8c
8. **Refactor tests**:
   * `tests/test_pool_fl.py` → `tests/test_fl_models.py`
   * remove the "RNA passthrough object identity" test (it asserted
     a leaky behaviour we are removing)
   * add **EB symmetry** test: identical raw RNA & gDNA counts
     ⇒ identical `log_likelihood` curves
   * add **EXON_ONLY isolation** test: 1e6 EXON_ONLY fragments must
     not shift `gdna.log_likelihood` at any FL bin
   * add **`Diagnostics.total() == payload.n_observed`** sentinel
   * add **no `mask_*` keys in JSON** regression guard
   * `tests/test_calibration_result_v6.py`: rename `pool` → `fl_models`;
     add `rna_fl_mean` / `global_fl_mean` checks; type-check
     `diagnostics` is a `Diagnostics`, not a dict.

---

## 7. M1–M6 audit (deltas vs v1 of this plan)

Same conclusions as v1.  One addition:

* **`FragmentLengthModels.gdna_model` slot is dead on the v6 path**
  (the v1 SRD path injected calibration's gDNA model into it; v6 does
  not).  Marked for deletion in M8c alongside the rest of the class.

---

## 8. Why this is "simple, elegant, beautiful"

* **One owner.** `rigel.calibration.fl.build_fl_models` is the sole
  function that turns raw counts into finalized `FragmentLengthModel`
  objects on the v6 path.
* **One policy.** `_finalize(...)` is called three times (once for
  global, once for RNA, once for gDNA).  The classifier
  (`_classify_and_build`) is called twice (once for RNA, once for
  gDNA) with identical signature.  Total FL logic: ~30 lines.
* **One product.** `FLModels` is the only finalized FL surface
  downstream code sees.  No duplicate `rna_model` slot, no
  scanner-side `finalize(prior_ess=...)` competing for ownership.
* **Explicit data flow.** `build_fl_models` takes three count
  vectors as kwargs.  No hidden dependency on a payload object, no
  pulling fields off a scanner — the function signature *is* the
  contract.
* **Zero recomputation.** Each raw histogram is accumulated once (in
  C++), extracted once (one-liner per source), and finalized once
  (one EB call per pool).  Three `np.ndarray` → three
  `FragmentLengthModel`.  Done.
* **Mask abstraction stays sealed.**  The 3-bit mask is referenced
  in exactly one Python site after M7: `extract_gdna_counts(...)`.
  No public field, no JSON key, no dataframe column ever names it.

---

## 9. Risks & open questions

1. **Symmetric-ESS regression vs SRD-v1 weights.**  The scanner's
   `FragmentLengthModels.finalize(prior_ess=...)` was tuned for SRD-v1.
   Switching v6 to a single `POOL_EB_PRIOR_ESS = 1000` is a deliberate
   policy unification; downstream EM behaviour will be re-validated
   in M8 against the benchmark suite.

2. **`SpliceType.SPLICED_ANNOT` is the sole RNA source.**
   `build_scoring_models()` historically copied SPLICED_ANNOT counts
   into `rna_model` and *also* logged a warning when SPLICED_ANNOT was
   empty.  v6 routes around `build_scoring_models()` entirely; if
   SPLICED_ANNOT is empty, RNA falls into the `fallback` branch
   (identity-shares `global_`).  Acceptable behaviour, identical
   outcome to the legacy "no spliced data" path.

3. **C++ collapse (M7+a) deferred.**  Today the scanner emits an
   8×1024 `fl_hist` of which only 3 rows feed v6.  This is wasteful
   but correct.  A standalone follow-up commit will reduce it to
   3×1024 and rename the C++ accumulator slots; M7 keeps `fl_hist`
   shape stable so the change is purely cosmetic.

---

## 10. Out-of-scope follow-ups

* **M7+a**: collapse C++ `fl_hist` 8×N → 3×N; emit named arrays
  (`global_fl_hist`, `rna_fl_hist`, `gdna_fl_hist`).
  `_fl_sources.py` becomes a one-liner-per-source.
* **M8**: wire `CalibrationResult` v6 into the scorer
  (`scoring.py`) and `pipeline.quant_from_buffer`.  The scorer reads
  `fl_models.rna.log_likelihood(...)` /
  `fl_models.gdna.log_likelihood(...)` and never mentions
  `FragmentLengthModels`.
* **M8c**: delete `FragmentLengthModels`, `_simple.py`, `_result.py`
  (legacy v1).  Rename `_result_v6.py` → `result.py`,
  `CalibrationResultV6` → `CalibrationResult`.
