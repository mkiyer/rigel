# M8 — `calibrate()` orchestrator + pipeline integration

**Status**: detailed plan, ready for implementation.  Refines and (in
several places) supersedes `docs/calibration/calibration_v6_plan.md`
§M8 to reflect the actual state of the M7 v2 surface that shipped in
commit `db25fad`.

**Scope**: pure-Python; no C++ ABI change.  Replace the SRD-v1 orchestrator
(`calibrate_gdna`) and the v1 `CalibrationResult` schema with the v6
surface composed from M3–M7.  Three focused commits (M8a / M8b / M8c).

---

## 1. The One Sentence

> **One module owns calibration orchestration: `rigel.calibration._orchestrator`.
> One function composes the M3–M7 building blocks: `calibrate(...)`.
> One result type flows through the pipeline: `CalibrationResult` (v6 schema).
> The scorer reads finalized FL distributions directly from `result.fl_models`;
> there is no backflow mutation of `FragmentLengthModels`.**

Everything in this plan follows from that sentence.

---

## 2. Goals

1. Replace `calibrate_gdna(...)` with a single `calibrate(...)` that
   composes the already-shipped M4 (global densities), M7 (FL models),
   and M6 (priors) primitives — no recomputation, no parallel paths.
2. Plumb the v6 `CalibrationResult` through `pipeline.run_pipeline` and
   `pipeline.quant_from_buffer` end-to-end.  Eliminate the v1
   `frag_length_models.gdna_model = gdna_copy` mutation by having the
   scorer consume `FLModels` directly.
3. Delete the SRD-v1 surface: `_simple.py`, `_categorize.py`,
   `_fl_mixture.py`, `_fl_empirical_bayes.py`, the v1 `_result.py`,
   their CLI flags, their config fields, and their tests.

## 3. Non-goals

- No C++ change.  `bam_scanner.cpp::set_regions` is wired (M3); the
  scanner-side `FragmentLengthModels` instance becomes a typed view that
  carries `global_model` + `category_models[SPLICED_ANNOT]` raw counts
  only — its v1 `rna_model` / `gdna_model` slots disappear in M8c.
- No new statistical primitives.  Every per-sample quantity needed by
  the new orchestrator already exists in M3–M7.
- No benchmark validation.  M9 owns the validation matrix.
- No deprecation grace period.  Master plan §M8b paragraph allowed
  "deprecation-warning shims" for v1 CLI flags; we drop that — v1 flags
  are deleted in M8c, not warned about.  `nrna_weight` and `c_base` are
  already configurable via the M5/M6 surfaces.

---

## 4. Issues to resolve / improvement opportunities discovered post-M7

These are the deltas between the master plan §M8 and the actual code
state.  All are fixed in this plan.

### 4.1 The master-plan `calibrate(...)` signature is stale

Master plan §M8 prescribes:

```python
def calibrate(buffer, index, payload, frag_length_models, strand_models,
              *, pool_quality_thresholds=(5_000, 200)) -> CalibrationResult:
    pool = compute_pool_fl_models(payload.fl_hist, ...)
    global_dens = compute_global_densities(...)
    return build_calibration_result(pool_models=pool, ..., region_df=...,
                                    strand_specificity=...)
```

But M7 v2 shipped a different surface: there is no `compute_pool_fl_models`,
no `PoolFLModels`, no `_fl_pool.py`.  The owner is `rigel.calibration.fl`
with `build_fl_models(*, global_counts, rna_counts, gdna_counts, max_size,
prior_ess, good_threshold, weak_threshold) -> FLModels`, and counts are
extracted by the three pure functions in `_fl_sources.py`.

`build_calibration_result(...)` already exists in `_result_v6.py` and has
the signature `(*, payload, scan_trained, global_densities, prior_table,
fl_prior_ess)`.  It internally calls `extract_*_counts` + `build_fl_models`.

**Resolution**: M8 `calibrate(...)` is a six-line orchestrator that calls
`compute_global_densities(...)` and `build_calibration_result(...)`.  No
new FL-related code.  No `pool_quality_thresholds` parameter — the
M7 v2 thresholds are constants in `fl.py` (`POOL_QUALITY_GOOD_THRESHOLD`,
`POOL_QUALITY_WEAK_THRESHOLD`); we add a single `fl_prior_ess` knob in
`CalibrationConfig` for the EB shrinkage strength.

### 4.2 `build_calibration_result` requires a `prior_table` upfront — but priors don't exist until `build_loci`

The master plan's flow is "`calibrate(...)` returns a `CalibrationResult`,
then later `result.with_priors(...)` backfills".  Today `build_calibration_result`
requires `prior_table` as a kwarg.

**Resolution**: make `prior_table` default to `PriorTable.empty()` (a new
classmethod returning a degenerate table with zero-length arrays).
`calibrate(...)` passes that empty table; `quant_from_buffer` later calls
`assemble_priors(...)` and `result.with_priors(...)` to swap in the real
table.  `CalibrationResult.alpha_gdna` etc. on the empty result are
zero-length arrays, never indexed because no locus exists yet.

### 4.3 The `pipeline.py:778` mutation is a v1 backflow bridge

```python
frag_length_models.gdna_model = gdna_copy  # SRD-v1 backflow
```

After `calibrate_gdna(...)` produces a `gdna_fl_model`, v1 code copies
it into the scanner-trained `FragmentLengthModels` so the scorer can read
it via `frag_length_models.gdna_model`.  This is the only reason
`FragmentLengthModels` has a mutable `gdna_model` slot.

**Resolution**: `FragmentScorer` is refactored in M8b to take
`fl_models: FLModels` instead of `frag_length_models: FragmentLengthModels`.
The mutation goes away.  Likewise `pipeline.py:915`
`frag_length_models.build_scoring_models()` becomes a no-op and is removed.

### 4.4 `pipeline.py:384–385` reaches into `frag_length_models.rna_model.compute_all_transcript_eff_lens(...)`

This works because `FLModels.rna` IS a `FragmentLengthModel`.  The fix is
a one-line substitution: `fl_models.rna.compute_all_transcript_eff_lens(...)`.

### 4.5 `_run_locus_em_partitioned` is missing the M5 `prior_weight_rna_per_locus` parameter

M5's commit shipped `build_prior_weight_rna(...)` and the C++ ABI for
`prior_weight_rna`, but the Python entrypoint `_run_locus_em_partitioned`
was never updated to plumb the per-locus arrays into `run_em_batch(...)`.

**Resolution**: M8b adds the `prior_weight_rna_per_locus: list[np.ndarray]
| None = None` parameter and forwards it to `run_em_batch`.  Pipeline
constructs the list from `result.prior_table.prior_weight_rna` (which is
already populated by `assemble_priors`).

### 4.6 `n_oor` vs `n_unannotated_ref` naming drift

Master plan §2.10 calls the out-of-region counter `n_oor`; the shipped
`CalibrationScanPayload` field is `n_unannotated_ref`.  No code change
required — `Diagnostics.from_payload` and `to_summary_dict` already
handle the actual field name.  This plan flags the discrepancy as
documented-only; M8 does not rename.

### 4.7 `boundary_flux_gdna_summary` and `kappa_diagnostics` keys

Master plan §M7 prescribed these as fields on `CalibrationResult`, but
M7 v2 omitted them (single-owner discipline: `Diagnostics` is the typed
container, not a free-form `dict`).  The current `CalibrationResult` v6
exposes `kappa` per-density via `global_densities.{intergenic,intron,
exon_intron}.kappa` and `boundary_flux` totals via
`global_densities.exon_intron.{n_fragments, eff_length_bp,
n_regions_used}`.  M8b's `to_summary_dict()` flattens those into the
documented summary keys without adding new dataclass fields.

### 4.8 Master plan's "v1 keys emitted in parallel" is rejected

We do a clean cut.  Pre-existing external dashboards consume
`summary.json` from rigel runs that no one else ships; double-emit costs
nothing but adds a sunset commitment.  `summary.json` v6 keys land in
M8b; v1 keys are removed in the same commit.

---

## 5. Audit — current state at HEAD `db25fad` (post-M7 v2)

| Surface | Location | Status |
|---|---|---|
| `calibrate_gdna(...)` (v1) | `src/rigel/calibration/_simple.py` | live, replaced in M8b |
| `CalibrationResult` (v1) | `src/rigel/calibration/_result.py` | live, deleted in M8c |
| `CalibrationResultV6` alias | `src/rigel/calibration/__init__.py:24` | re-exports `_result_v6.CalibrationResult`; renamed in M8c |
| `build_calibration_result` | `src/rigel/calibration/_result_v6.py:180` | shipped (M7 v2); kwarg `prior_table` made optional in M8a |
| `FLModels`, `build_fl_models` | `src/rigel/calibration/fl.py` | shipped (M7 v2) |
| `compute_global_densities` | `src/rigel/calibration/density_global.py:100` | shipped (M4) |
| `assemble_priors` | `src/rigel/calibration/locus_prior.py:352` | shipped (M6) |
| `compute_locus_priors_from_partitions` | `src/rigel/locus.py:202` | live, replaced in M8b by `assemble_priors` |
| `FragmentScorer.from_models(... frag_length_models, ...)` | `src/rigel/scoring.py:94` | refactored in M8b to take `fl_models: FLModels` |
| `_run_locus_em_partitioned` | `src/rigel/pipeline.py:514` | gains `prior_weight_rna_per_locus` in M8b |
| `frag_length_models.gdna_model = gdna_copy` mutation | `src/rigel/pipeline.py:778` | deleted in M8b |
| `frag_length_models.build_scoring_models()` | `src/rigel/pipeline.py:915` | deleted in M8b |
| `set_regions(...)` C++ wiring | `bam_scanner.cpp:990` | shipped (M3); no change |
| `CalibrationConfig` v1 fields | `src/rigel/config.py:154` | replaced in M8a (additive) and cleaned in M8c |
| CLI `--calibration-*` flags | `src/rigel/cli.py:1077` | replaced in M8a (additive) and cleaned in M8c |
| Tests `test_calibration_simple.py`, `test_categorize.py`, `test_gdna.py`, `test_gdna_harmonic_length.py` | `tests/` | deleted in M8c |

---

## 6. Design

### 6.1 `_orchestrator.py` — the six-line orchestrator

```python
# src/rigel/calibration/_orchestrator.py
from .density_global import compute_global_densities
from ._result_v6   import build_calibration_result, CalibrationResult
from .scan_payload import CalibrationScanPayload
from ..frag_length_model import FragmentLengthModels  # type-only at runtime
from ..index import TranscriptIndex
from .fl import POOL_EB_PRIOR_ESS

def calibrate(
    *,
    index:        TranscriptIndex,
    payload:      CalibrationScanPayload,
    scan_trained: FragmentLengthModels,
    fl_prior_ess: float = POOL_EB_PRIOR_ESS,
) -> CalibrationResult:
    if index.region_df is None:
        raise RuntimeError(
            "Index has no region table. Rebuild the index "
            "(rigel index --fasta ... --gtf ...). Older indexes "
            "are not supported."
        )
    # Build FL models first so we can pass gdna_fl_mean to global densities.
    # build_calibration_result internally calls build_fl_models, so we
    # call it once with an empty PriorTable; densities go in unchanged.
    global_densities = compute_global_densities(
        index.region_df, payload, gdna_fl_mean=None,  # filled lazily below
    )
    return build_calibration_result(
        payload          = payload,
        scan_trained     = scan_trained,
        global_densities = global_densities,
        prior_table      = None,           # empty default — see §6.2
        fl_prior_ess     = fl_prior_ess,
    )
```

The `gdna_fl_mean=None` path in `compute_global_densities` is the M4
contract — the function uses the global-mean fallback when the gdna
mean is unknown at call time.  After `build_calibration_result` runs we
have `result.fl_models.gdna.mean`; the per-`Locus` code in
`assemble_priors` re-runs the locoregional shrink with that value, so
re-calling `compute_global_densities` with the now-known mean is not
necessary (priors get the right value in §6.4).

> **Improvement opportunity (resolved here):** the master plan called
> for `calibrate(...)` to know the gDNA-FL mean before computing
> densities.  In the shipped M4 the densities only use the mean for the
> `EXON-INTRON` denominator; `assemble_priors` recomputes loco
> densities with the up-to-date mean.  We don't double-compute global
> densities.

### 6.2 `PriorTable.empty()` — the lazy-backfill seed

Add a classmethod to `locus_prior.PriorTable`:

```python
@classmethod
def empty(cls) -> "PriorTable":
    return cls(
        multi_locus_priors = (),
        alpha_gdna         = np.empty(0, dtype=np.float64),
        alpha_rna          = np.empty(0, dtype=np.float64),
        prior_weight_rna   = (),
        c_base_value       = C_BASE_DEFAULT,
    )
```

`build_calibration_result(..., prior_table=None)` substitutes
`PriorTable.empty()`.  The eager `multi_locus_prior_df` /
`per_locus_gdna_df` builders return empty DataFrames with the documented
schema.  `with_priors(real_table)` rebuilds them.

### 6.3 `CalibrationConfig` (v6, additive in M8a, sole schema in M8c)

```python
@dataclass(frozen=True, slots=True)
class CalibrationConfig:
    fl_prior_ess:      float = 1000.0   # POOL_EB_PRIOR_ESS
    nrna_weight:       float = 0.0      # M5 default
    c_base:            float = 10.0     # C_BASE_DEFAULT
    # v1 fields kept ONLY in M8a/M8b (deleted in M8c):
    exon_fit_tolerance_bp: int = 0
    fl_prior_ess_v1:        float = 500.0   # renamed to avoid collision
    max_iter:               int = 1000
    tol:                    float = 1e-4
```

Wait — there's a name collision: M8a needs both v6 `fl_prior_ess` (M7
default 1000) and v1 `fl_prior_ess` (default 500).  We resolve it by
renaming the v1 field to `fl_prior_ess_v1` for the M8a/M8b interim and
deleting both v1 fields in M8c (`fl_prior_ess` becomes the sole knob,
default 1000).  Pipeline reads `config.calibration.fl_prior_ess` in M8b.

### 6.4 `pipeline.py` — the new shape

`scan_and_buffer` (signature unchanged in M8b; the calibration payload
is already the 5th tuple element).

`run_pipeline` (M8b):

```python
result = scan_and_buffer(...)
stats, strand_models, frag_length_models, buffer, payload = result
calibration = calibrate(
    index        = index,
    payload      = payload,
    scan_trained = frag_length_models,
    fl_prior_ess = config.calibration.fl_prior_ess,
)
quant_result = quant_from_buffer(
    buffer       = buffer,
    index        = index,
    config       = config,
    strand_models= strand_models,
    calibration  = calibration,
)
```

`quant_from_buffer` (M8b):

```python
em_data, multi_loci = build_loci(buffer, index, calibration.fl_models, ...)
prior_table = assemble_priors(
    multi_loci       = multi_loci,
    em_data          = em_data,
    index            = index,
    payload          = calibration_payload,   # carried through; see note
    global_densities = calibration.global_densities,
    gdna_fl_mean     = calibration.fl_models.gdna.mean,
    c_base           = config.calibration.c_base,
)
calibration = calibration.with_priors(prior_table)

prior_weight_rna_per_locus = [
    build_prior_weight_rna(ml, em_data, nrna_weight=config.calibration.nrna_weight)
    for ml in multi_loci
]

scorer = FragmentScorer(fl_models=calibration.fl_models, ...)
em_results = _run_locus_em_partitioned(
    em_data, multi_loci,
    alpha_gdna                 = calibration.alpha_gdna,
    alpha_rna                  = calibration.alpha_rna,
    prior_weight_rna_per_locus = prior_weight_rna_per_locus,
    ...
)
```

Note: `assemble_priors` needs the raw `payload`, not just the
`CalibrationResult`.  We carry the payload through `quant_from_buffer`
as a parameter (unchanged end-to-end signature: pipeline already passes
it to `quant_from_buffer` after `scan_and_buffer`).

### 6.5 `FragmentScorer` ABI change

Old:
```python
FragmentScorer.from_models(strand_models, frag_length_models, index, ...)
# reads frag_length_models.rna_model, frag_length_models.gdna_model
```

New (M8b):
```python
FragmentScorer.from_models(strand_models, fl_models: FLModels, index, ...)
# reads fl_models.rna, fl_models.gdna
```

This is a hard break, no compat shim.  Two call sites in
`pipeline.py` and one in `tests/`; tests adjust along with the change.

### 6.6 `summary.json` (M8b — v6 keys only)

`CalibrationResult.to_summary_dict()` already produces v6-shaped output
(via M7 v2).  CLI's `cal_dict = result.calibration.to_summary_dict()`
keeps working.  The legacy keys (`pi_pool` at top-level instead of
nested, `srd.*`, etc.) are dropped.  Schema:

```json
{
  "calibration": {
    "pi_pool":        <float>,
    "fl_models":      { "global": {...}, "rna": {...}, "gdna": {...},
                        "rna_quality": "...", "gdna_quality": "..." },
    "global_densities": { "intergenic": {...}, "intron": {...},
                          "exon_intron": {...} },
    "kappa": { "intergenic": {...}, "intron": {...}, "exon_intron": {...} },
    "diagnostics":    { "n_unannotated": ..., "n_intron_only": ..., ... }
  }
}
```

### 6.7 CLI flags (v6)

| New flag (M8a, sole in M8c) | Maps to | Default |
|---|---|---|
| `--cal-fl-prior-ess` | `calibration.fl_prior_ess` | `1000.0` |
| `--cal-nrna-weight` | `calibration.nrna_weight` | `0.0` |
| `--cal-c-base` | `calibration.c_base` | `10.0` |

Master plan listed `--cal-quality-good` / `--cal-quality-weak`; we omit
them (the thresholds are constants in `fl.py`; if a future user needs
to tune them, that's a one-line `CalibrationConfig` addition then).

V1 flags (`--calibration-exon-fit-tolerance-bp`,
`--calibration-fl-prior-ess`) are kept in M8a additively (with renamed
underlying field), and deleted in M8c.

---

## 7. Migration — three commits

### M8a — Additive (new surface alongside SRD-v1)

Files added:

- `src/rigel/calibration/_orchestrator.py` — `calibrate(...)` (§6.1).

Files modified:

- `src/rigel/calibration/locus_prior.py` — add `PriorTable.empty()`.
- `src/rigel/calibration/_result_v6.py` — `prior_table` kwarg defaults
  to `None` → `PriorTable.empty()`; helpers handle empty table.
- `src/rigel/calibration/__init__.py` — export `calibrate`.
- `src/rigel/config.py` — rename existing `fl_prior_ess` to
  `fl_prior_ess_v1`; add v6 fields `fl_prior_ess=1000.0`,
  `nrna_weight=0.0`, `c_base=10.0`.
- `src/rigel/cli.py` — add `--cal-fl-prior-ess`, `--cal-nrna-weight`,
  `--cal-c-base`; rename existing `--calibration-fl-prior-ess` →
  internally maps to `fl_prior_ess_v1`.

Files **not** touched: `pipeline.py`, `scoring.py`, `_simple.py`.

Tests added (target ≥ 6):

- `tests/test_calibrate_orchestrator.py`:
  - `test_calibrate_returns_v6_result` — happy path on synthetic payload.
  - `test_calibrate_raises_on_legacy_index` — `region_df=None` → RuntimeError.
  - `test_calibrate_with_priors_roundtrip` — empty seed → `with_priors` swap.
  - `test_calibrate_fl_models_owns_finalize` — no `FragmentLengthModels`
    mutation observed.
  - `test_calibrate_fl_prior_ess_propagates` — different ESS → different
    EB-shrunk densities.
  - `test_prior_table_empty_invariants` — empty table dataframe schemas
    match v6 schema with zero rows.

**Exit gate**: ≥ 6 new tests green; protected suite green; v1 path still
wired and unchanged.

### M8b — Switchover (pipeline replaces v1 wiring)

Files modified:

- `src/rigel/scoring.py` — `FragmentScorer.from_models(...)` accepts
  `fl_models: FLModels` instead of `frag_length_models: FragmentLengthModels`.
  Reads from `fl_models.rna` / `fl_models.gdna`.
- `src/rigel/pipeline.py`:
  - `run_pipeline` calls `calibrate(...)` instead of `calibrate_gdna(...)`.
  - Drop the `frag_length_models.gdna_model = gdna_copy` mutation
    (line 778) and the `frag_length_models.build_scoring_models()`
    call (line 915).
  - `quant_from_buffer` accepts the v6 `CalibrationResult`; calls
    `assemble_priors(...)` then `calibration.with_priors(...)` after
    `build_loci`.
  - Switch line 384–385 to `fl_models.rna.compute_all_transcript_eff_lens(...)`.
  - `_run_locus_em_partitioned` gains `prior_weight_rna_per_locus`
    parameter; forwards to `run_em_batch`.
  - Build the per-locus `prior_weight_rna_per_locus` list from
    `calibration.prior_table.prior_weight_rna` (per-locus arrays already
    built by `assemble_priors` via `build_prior_weight_rna`).
- `src/rigel/cli.py` — `cal_dict = result.calibration.to_summary_dict()`
  keeps working; verify keys are v6-only.
- `scripts/profiling/profiler.py:70` — switch from `calibrate_gdna(...)`
  to `calibrate(...)`.

Files **not** touched: `_simple.py`, `_categorize.py`, `_fl_mixture.py`,
`_fl_empirical_bayes.py`, the v1 `_result.py`, the v1 tests.  They
become dead code in this commit, deleted in M8c.

Golden outputs regenerated; diffs documented in commit message.

Tests added (target ≥ 8):

- `tests/test_pipeline_integration_v6.py`:
  - `test_pipeline_pristine_rna` — `gdna_none` synthetic, `pi_pool ≈ 0`.
  - `test_pipeline_with_gdna` — contaminated synthetic, `pi_pool > 0.05`.
  - `test_pipeline_no_v1_mutation` — assert
    `frag_length_models.gdna_model is None` after `run_pipeline`.
  - `test_scorer_consumes_fl_models` — FragmentScorer constructed with
    `FLModels`, reads `.rna` / `.gdna`.
  - `test_prior_weight_rna_per_locus_plumbed` — non-zero `nrna_weight`
    propagates from CLI through to EM.
  - `test_with_priors_dataframes_populated` — post-`build_loci`,
    `result.multi_locus_prior_df` has `n_multi_loci` rows.
  - `test_assemble_priors_called_once` — counter on the call site.
  - `test_summary_json_v6_schema` — keys present, no `srd.*` leakage.

Existing tests adjusted: `tests/test_pipeline_smoke.py`,
`tests/test_oracle_bam.py`, `tests/test_estimator.py` (FragmentScorer
constructor change), `tests/test_calibration_simple.py` (left in place
for M8b, deleted in M8c), `tests/test_golden_output.py` (regenerated).

**Exit gate**: full suite green (incl. golden regeneration); benchmark
matrix not regressed vs M8a (smoke check on `gdna_none_ss_1.00_nrna_none`
+ `gdna_high_ss_0.90_nrna_none`).

### M8c — Subtractive (legacy deletion)

Files deleted:

- `src/rigel/calibration/_simple.py`
- `src/rigel/calibration/_categorize.py`
- `src/rigel/calibration/_fl_mixture.py`
- `src/rigel/calibration/_fl_empirical_bayes.py`
- `src/rigel/calibration/_result.py` (v1)
- `tests/test_calibration_simple.py`
- `tests/test_categorize.py`
- `tests/test_gdna.py`
- `tests/test_gdna_harmonic_length.py`

Files renamed:

- `src/rigel/calibration/_result_v6.py` → `_result.py`.
- (No file renames for tests; `tests/test_calibration_result_v6.py`
  keeps its name.)

Symbols renamed:

- `CalibrationResultV6` (alias in `__init__.py`) → `CalibrationResult`
  (canonical export).  Internal class name was already `CalibrationResult`
  inside `_result_v6.py`; only the alias goes away.

Files modified:

- `src/rigel/calibration/__init__.py` — re-export from `._result`,
  drop `CalibrationResultV6` alias.
- `src/rigel/config.py` — drop `fl_prior_ess_v1`, `exon_fit_tolerance_bp`,
  `max_iter`, `tol`.  `fl_prior_ess` becomes the sole field (already
  added in M8a).
- `src/rigel/cli.py` — drop `--calibration-exon-fit-tolerance-bp`,
  `--calibration-fl-prior-ess`, drop the corresponding `_PARAM_SPECS`
  entries.
- `CHANGELOG.md` — entry under "Breaking changes".

Optional cleanup (deferred if scope creeps): move general FL utilities
that survived to `src/rigel/frag_length_eb.py` /
`src/rigel/frag_length_mixture.py`.  M7 v2's single-owner discipline
already eliminated `_fl_empirical_bayes.py`'s call sites except
`_simple.py`; once `_simple.py` is gone, the file is dead and is
deleted (no move needed).  Same for `_fl_mixture.py`.

**Exit gate**: full suite green; `git grep -r 'calibrate_gdna\|_simple\|_categorize\|_fl_mixture\|_fl_empirical_bayes\|CalibrationResultV6'` returns nothing in `src/` or `tests/`; CHANGELOG entry written.

---

## 8. Test plan summary

| Commit | Tests added | Tests modified | Tests deleted |
|---|---|---|---|
| M8a | 6 (`test_calibrate_orchestrator.py`) | 0 | 0 |
| M8b | 8 (`test_pipeline_integration_v6.py`) | ~5 (scorer ABI, golden regen) | 0 |
| M8c | 0 | 1 (`test_calibration_result_v6.py` import paths) | 4 (v1 tests) |

Combined ≥ 14 new tests green; full suite green at every commit boundary.

Protected suite (`tests/test_index.py`, `test_buffer.py`, `test_resolution.py`,
`test_em_impl.py`, `test_estimator.py`, `test_pipeline_smoke.py`,
`test_oracle_bam.py`, `test_cli.py`, `test_index_integrity.py`) green at
each commit.

---

## 9. Risks

| # | Risk | Mitigation |
|---|---|---|
| 1 | `FragmentScorer` ABI change breaks downstream callers (notebooks, scripts) | Single-commit cut; grep `from rigel.scoring import` in `scripts/`; update in M8b. |
| 2 | Golden output regeneration hides a real regression | M8b commit message documents the diff (per-locus alpha_gdna / alpha_rna magnitudes); M9 benchmark matrix catches numerical regressions. |
| 3 | `PriorTable.empty()` introduces a path that crashes if `with_priors` is forgotten | `quant_from_buffer` always calls `with_priors`; assertion `len(alpha_gdna) == n_multi_loci` in `_run_locus_em_partitioned` catches the omission. |
| 4 | `assemble_priors` needs `payload`, but `CalibrationResult` doesn't carry it | Pipeline plumbs the raw `payload` from `scan_and_buffer` into `quant_from_buffer` directly; this is already its existing shape. |
| 5 | Profiler script (`scripts/profiling/profiler.py:70`) breaks at M8b | Updated in M8b alongside the pipeline change. |
| 6 | Scanner-side `FragmentLengthModels` still has `rna_model`/`gdna_model` slots after M8c | Acceptable; those slots become optional v1 vestiges that are simply never set on the v6 path.  A follow-up plan can prune them. |

---

## 10. Out of scope

- C++ ABI change (no recompilation in M8).
- Benchmark validation (M9).
- Mappability-adjusted `L_eff` (future plan).
- `c_base(ℓ)` formula (constant 10.0 stays).
- Renaming scanner-side `FragmentLengthModels` slots.

---

## 11. References

- `docs/calibration/calibration_v6_plan.md` §M8 — master spec
  (this document refines and partly supersedes).
- `docs/calibration/m7_implementation_plan.md` (v2) — style + ownership
  template; M7 v2 commit `db25fad`.
- `docs/calibration/m6_implementation_plan.md` — `assemble_priors` API.
- `docs/calibration/m5_implementation_plan.md` — `prior_weight_rna` ABI.
- Memory: `/memories/repo/rigel-public-interface-2026-03-09.md`,
  `/memories/repo/rigel-tripartite-prior-2026-03-14.md`.
