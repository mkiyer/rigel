# Fractional Accumulator — Python Cutover Plan

Date: 2026-05-22
Status: planning document, native cutover complete; Python cutover in progress
Parent design: [fractional_accumulator_phase3_plan.md](fractional_accumulator_phase3_plan.md) (§9.3+)
Companion: [fine_region_migration_phase0_log.html](fine_region_migration_phase0_log.html)

Locked decisions (updated 2026-05-22):

- **§3 orient cleanup**: delete `_orient.py`; replace with a single
  `sense_antisense_split(...)` function in `fractional_evidence.py`. Move
  `StrandSummary` to `strand_summary.py` next to `strand_model.py`.
- **§5 density / prior policy**: Option A (fail fast). `quant_from_buffer`
  refuses to run with a typed error until the Phase 4 estimator lands.
- **§6.3 FL EB / mixture**: end-to-end float64 fractional weights (option
  6.3b).

## 1. Purpose And Scope

The Phase 3 native cutover is complete:

- `CalibrationAccumulator` and `CalibrationPayload` carry fractional mass
  (`region_counts[R, 12]` float32, `fl_pool_mass[6, 1024]` + `fl_pool_total[6]`
  float64, `channel_mass[12]`, `signature_mass[16]`, and the new counters
  `n_excluded_strand_ambig` / `n_fl_unavailable`).
- `BamScanner.set_regions(...)` no longer takes `splicing_anchor_tolerance`;
  the resolver still consumes it for implicit-splice slack and the scanner
  surfaces it on the payload as `resolver_splicing_anchor_tolerance`.
- Native rebuild succeeds; the new payload dict is exported from `scan()`.

Every Python consumer of the legacy schema (`global_counts`,
`per_region_counts`, `fl_hist`, `u_left`, `u_right`, `*_by_orient`,
`n_below_tolerance`, `splicing_anchor_tolerance` as a calibration field) is
now broken at import time. This document is the implementation plan for the
hard cutover that brings the Python layer back in line with the fractional
payload.

In scope:

- Rewrite `scan_payload.py` and its validation.
- Add `fractional_evidence.py` (the `FractionalEvidenceView`) as the single
  Python interface to the new evidence.
- Rewrite `_fl_sources.py`, `_diagnostics.py`, `_arrays.py`.
- Decide and implement the policy for `density_global.compute_global_densities`,
  `_regional_exposure.RegionalGdnaExposure`, and `locus_prior.*`.
- Update `_result.CalibrationResult` and `_exposure` helpers.
- Update `BamScanConfig.splicing_anchor_tolerance` docstring and CLI metadata.
- Rewrite test factories that hand-build calibration payloads.
- Run pytest + ruff and update goldens / log progress.

Out of scope (deferred to Phase 4 or its own change set):

- The new gDNA density and per-locus prior estimator. Phase 3 either fails
  fast or ships a clearly labelled minimal current-compatible view (§5
  decision).
- RNA-side FL extraction (`extract_rna_counts` reads
  `scan_trained.category_models[SPLICED_ANNOT]`, which is unchanged).
- Index file schema. Existing v4 indexes already carry `signature`,
  `neighbor_signature`, and `boundary_kind_{left,right}` and require no
  reindex.

## 2. Decision Already Locked

These came out of the §9.1 planning checkpoint, the §11 resolutions in
the parent plan, and the explicit user direction recorded in repo memory
`fine-region-fractional-boundary-design-2026-05-22.md`.

1. **Payload schema is the contract**. `CalibrationScanPayload` carries
   the native arrays as-is. No compatibility synthesis of legacy fields.
2. **Six FL pools**, indexed by `fl_pool_index(signature, compartment)`,
   are the only FL evidence surface. gDNA FL aggregate is computed in
   Python as `INTERGENIC_CONTAINED + INTERGENIC_BOUNDARY +
   INTRONIC_CONTAINED + INTRONIC_BOUNDARY`. EXON pools are diagnostic.
3. **FL pools are unspliced-only and strand-collapsed.** `region_counts`
   keeps the strand axis (pos / neg). Strand-ambiguous fragments emit no
   mass and are counted in `n_excluded_strand_ambig`.
4. **Resolver K is provenance only.** `BamScanConfig.splicing_anchor_tolerance`
   stays public and CLI-exposed, but it controls resolver implicit-splice
   slack only. Calibration accumulator and payload validation must not
   read it as a calibration tolerance.
5. **`boundary_kind_{left,right}` columns stay in `regions.feather`** as
   index-only metadata. The Phase 3 accumulator and Python consumers do
   not read them. Future estimators may join them back in via the region
   table if needed.
6. **Genomic-strand storage is canonical.** The accumulator stores fragment
   mass on POS / NEG channels (genomic strand). Sense / antisense relative
   to a region is derived at consumer time from
   `(region.signature, channel.strand)`. See §4.

## 3. Orient Cleanup — Clean Redesign

The legacy `_orient.py` mixes two unrelated concerns that the fractional
cutover lets us separate cleanly. We delete the module and split it
into two single-purpose pieces.

### 3.1 What `_orient.py` currently contains

1. `ORIENT_SAME` / `ORIENT_OPP` / `ORIENT_UNINF` / `ORIENT_N` integer
   constants — the third axis of the legacy `*_by_orient` per-region
   arrays. After the fractional cutover the accumulator no longer
   classifies fragments by orient, so these constants describe a
   storage layout that no longer exists.
2. `classify_orient(region_strand, fragment_strand)` — the per-fragment
   routing helper used at accumulation time. Native code already does
   the equivalent in C++ for `n_excluded_strand_ambig`, and the new
   accumulator stores fragment mass by genomic strand. This helper
   has no remaining call site.
3. `StrandSummary` — a summary of `StrandModel` (the read-1 sense
   probability and its standard error). Used by the calibration
   pipeline to warn when strand is unidentifiable. Completely
   independent of (1) and (2).

### 3.2 Target state

- **Delete** `src/rigel/calibration/_orient.py`.
- **Delete** `ORIENT_*` constants and `classify_orient` entirely. They
  describe an accumulator that no longer exists.
- **Move** `StrandSummary` to a new single-purpose module
  `src/rigel/calibration/strand_summary.py`. It is the
  `StrandModel`-derived public summary; nothing about it is
  orient-axis-shaped.
- **Add** one elegant API in `fractional_evidence.py` that subsumes
  every legitimate downstream use of orient (see §4). No integer
  enum, no lookup table, no `*_capability` field — just a function
  that takes a channel slice of `region_counts` and returns the
  `(sense, antisense, strand_ambiguous)` mass triple per region.

### 3.3 Why this is cleaner than keeping the enums

- The only legitimate downstream question is *"of this region's mass,
  how much agrees with the transcript strand, how much disagrees, and
  how much can't be classified?"* That is exactly the
  `sense_antisense_split` return signature. No consumer needs the
  intermediate `SAME` / `OPP` / `UNINF` per-fragment label.
- Removing the enum removes a whole axis of accidental complexity:
  callers stop carrying around `(R, 3)` integer arrays and stop
  reasoning about an `ORIENT_N == 3` magic number.
- Strand-ambiguous fragments are now accounted for explicitly at scan
  time (`n_excluded_strand_ambig` in the payload), so consumers do not
  need a per-fragment "uninformative" bucket anymore. The only
  remaining ambiguity is *per-region* (signatures with both
  transcript strands), and that is a property of the region table
  exclusively.

### 3.4 Imports to delete

Grep target after the cutover:

```text
rg 'from .* import .*ORIENT_(SAME|OPP|UNINF|N)\b' src/rigel/ tests/
rg 'from .* import .*classify_orient\b' src/rigel/ tests/
rg '\b_orient\b' src/rigel/ tests/
```

All three must return empty. `tests/test_orient_routing.py` is
deleted with the module.

## 4. Sense / Antisense Split — The Only Orient API

### 4.1 Per-region transcript-strand class

Given the signature bit layout

```text
bit 3 (0x8) intron_pos    bit 2 (0x4) intron_neg
bit 1 (0x2) exon_pos      bit 0 (0x1) exon_neg
```

define, vectorized over the region table, a single int8 column:

```python
TS_NONE  =  0   # no transcript strand (intergenic)
TS_POS   = +1   # only pos transcript strand bits set
TS_NEG   = -1   # only neg transcript strand bits set
TS_AMBIG =  2   # both pos and neg transcript strand bits set

def transcript_strand_class(signature: np.ndarray) -> np.ndarray:
    has_pos = (signature & 0b1010) != 0   # BIT_INTRON_POS | BIT_EXON_POS
    has_neg = (signature & 0b0101) != 0   # BIT_INTRON_NEG | BIT_EXON_NEG
    out = np.zeros(signature.shape, dtype=np.int8)
    out[ has_pos & ~has_neg] = TS_POS
    out[~has_pos &  has_neg] = TS_NEG
    out[ has_pos &  has_neg] = TS_AMBIG
    return out
```

This is the *only* orient-derived per-region datum the cutover
introduces. It is cached on `RegionArrays` (one int8 column, 1 B per
region).

### 4.2 The split function

`fractional_evidence.sense_antisense_split` is the entire public API
for strand-relative analysis of any channel slice:

```python
def sense_antisense_split(
    region_counts: np.ndarray,        # float32[R, 12]
    ts_class:      np.ndarray,        # int8[R], from transcript_strand_class
    *,
    compartment:   int,               # CONTAINED | BOUNDARY_LEFT | BOUNDARY_RIGHT
    splice:        int,               # UNSPLICED | SPLICED
) -> SenseSplit:
    """Return per-region (sense, antisense, ambiguous) float32 mass.

    For each region:
    * TS_POS : sense = mass on POS channel, antisense = mass on NEG channel
    * TS_NEG : sense = mass on NEG channel, antisense = mass on POS channel
    * TS_AMBIG or TS_NONE : both channels go to ambiguous; sense and
      antisense are zero.

    The function is a pure numpy expression — no Python loops, no
    intermediate label arrays — and returns float32 views with shape
    (R,) each.
    """
```

`SenseSplit` is a tiny frozen dataclass `(sense, antisense, ambiguous)`
of three float32[R] arrays.

### 4.3 What replaces every old orient use

| Old code                                               | New expression                                                   |
| ------------------------------------------------------ | ---------------------------------------------------------------- |
| `intron_by_orient[:, ORIENT_OPP]` (gDNA inside intron) | `split(CONTAINED, UNSPLICED).antisense` on intron-only regions   |
| `intron_by_orient[:, ORIENT_SAME]` (likely nRNA)       | `split(CONTAINED, UNSPLICED).sense` on intron-only regions       |
| `exon_contained_counts_by_orient[:, ORIENT_OPP]`       | `split(CONTAINED, UNSPLICED).antisense` on exon-only regions     |
| `u_left_by_orient[:, ORIENT_OPP]`                      | `split(BOUNDARY_LEFT, UNSPLICED).antisense`                      |
| `intron_by_orient[:, ORIENT_UNINF]`                    | `split(CONTAINED, UNSPLICED).ambiguous` on intron-only regions   |

No consumer ever needs a per-fragment SAME/OPP/UNINF label after the
cutover. The split function is sufficient and self-explanatory.

## 5. Density / Prior Policy — Locked: Option A (Fail Fast)

The Phase 3 cutover ships the fractional payload, evidence view, FL
models, and diagnostics. The density / regional-exposure / locus-prior
layer is rewritten under Phase 4. During the interim:

- `density_global.compute_global_densities`,
  `_regional_exposure.RegionalGdnaExposure.build`, and
  `locus_prior.assemble_priors` raise a typed
  `FractionalCutoverPending` exception with a single shared message
  pointing at the parent plan §11 and this document.
- `calibrate(...)` still returns a `CalibrationResult` carrying the
  fractional payload, FL models, diagnostics, and provenance. The
  `global_densities`, `regional_exposure`, and `priors` fields are
  `None`.
- `pipeline.quant_from_buffer` refuses to run end-to-end: it raises
  `FractionalCutoverPending` immediately after calibration completes,
  with the same message. There is no "uniform priors" fallback flag.
- `rigel quant` therefore prints the typed error and exits non-zero
  after the scan / calibrate stages finish. `rigel index` is
  unaffected.

`FractionalCutoverPending` lives in `rigel.calibration.errors` (new
module, one class, one shared message constant). Every fail-fast site
imports the same exception and message; Phase 4 will delete them all
in a single grep-and-replace.

What the user can still do on `main` during Phase 3:

- `rigel index ...` — unchanged.
- `rigel quant ...` — runs through scan + calibration, emits the
  fractional payload, FL models, and the new `summary.json`
  diagnostics, then exits with `FractionalCutoverPending`. Useful for
  validating the new diagnostics on real data.
- Direct API access to the payload, `FractionalEvidenceView`, and
  `Diagnostics` for ad-hoc analysis.

What blocks until Phase 4:

- Per-locus priors, batch EM, and any quant `.feather` output.
- The benchmarking workflows in `scripts/benchmarking/`.

This is the explicit, deliberate scope for the Phase 3 cutover.

## 6. Module-By-Module Rewrite

Each subsection is one todo from the active list.

### 6.1 `scan_payload.py` (todo 7)

Replace the dataclass to match the parent plan §3.3:

```python
@dataclass(frozen=True, slots=True)
class CalibrationScanPayload:
    region_counts: np.ndarray         # float32[R, 12]
    channel_mass: np.ndarray          # float64[12]
    signature_mass: np.ndarray        # float64[16]
    fl_pool_mass: np.ndarray          # float64[6, 1024]
    fl_pool_total: np.ndarray         # float64[6]

    n_observed: int
    n_excluded_multimap: int
    n_excluded_chimera: int
    n_excluded_artifact: int
    n_excluded_strand_ambig: int
    n_unobserved: int
    n_unannotated_ref: int
    n_fl_unavailable: int

    resolver_splicing_anchor_tolerance: int
    n_regions: int
```

`from_scan_dict` validates:

- Shapes / dtypes / contiguity exactly per the parent plan §3.3.
- All arrays finite and non-negative.
- `region_counts.sum()` ≈ `n_observed` (float32 tolerance ≈
  `1e-6 * n_observed + 1`).
- `channel_mass.sum() ≈ region_counts.sum()`.
- `signature_mass.sum() ≈ region_counts.sum()`.
- `fl_pool_mass.sum(axis=1)` matches `fl_pool_total` exactly (both
  built natively in the same pass).
- `fl_pool_mass.sum()` ≤ `n_observed - n_excluded_strand_ambig`
  (spliced fragments and `n_fl_unavailable` cases reduce it; equality
  is not required).
- Balance identity:
  `n_observed + n_excluded_multimap + n_excluded_chimera +
  n_excluded_artifact + n_excluded_strand_ambig + n_unobserved ==
  n_total` (when `n_total` is supplied; `n_unannotated_ref` is a
  *subset* of `n_observed` and must not double-count).

Remove every legacy mask constant from this module
(`MASK_EXON`, `MASK_INTRON`, `MASK_INTERGENIC`, `MASK_N_STATES`,
`FL_HIST_N_BINS`, `ORIENT_N`). Re-export `N_FL_POOLS`, `FL_POOL_*`
indices, `CHANNEL_*` indices, and `channel_index(...)` from
`signature.py` for ergonomic imports.

Add channel-view convenience properties on the payload (or as free
functions in `fractional_evidence.py`):

- `contained_unspliced(strand=None)` → float32[R] or float32[R, 2]
- `boundary_left_unspliced(strand=None)`
- `boundary_right_unspliced(strand=None)`
- `contained_spliced(strand=None)`
- `region_total_mass()` → float32[R]

### 6.2 `fractional_evidence.py` (todo 8 — new module)

Single public interface for fractional evidence. Owns:

- Signature-derived per-region masks: `is_intergenic`,
  `is_intron_only`, `is_exon_any`, `has_both_strands`.
- The `transcript_strand_class(signature)` int8 column from §4.1.
- `sense_antisense_split(region_counts, ts_class, *, compartment, splice)`
  from §4.2 and the `SenseSplit` dataclass.
- `FractionalEvidenceView(payload, region_arrays)` — caches
  `ts_class`, signature masks, and exposes per-channel helpers
  (`contained_unspliced`, `boundary_left_unspliced`, etc.) and
  `.split(compartment, splice) -> SenseSplit`.
- FL pool helpers: `gdna_fl_mass(payload)`,
  `rna_candidate_fl_mass(payload)`, `pool(payload, name)`,
  `pool_total(payload, name)`.
- Optional `BoundarySideEvidence` dataclass joining
  `region_counts[:, boundary_*]` to `signature` /
  `neighbor_signature`. Phase 4 will consume it; Phase 3 ships the
  builder so the interface is frozen.

No `ORIENT_*` constants are exported. No `_orient` re-export. No
backward-compatibility shim — every consumer either calls
`sense_antisense_split` or reads channel slices directly.

### 6.3 `_fl_sources.py` (todo 9)

`extract_gdna_counts(payload)` becomes a one-liner:

```python
def extract_gdna_counts(payload: CalibrationScanPayload) -> np.ndarray:
    return gdna_fl_mass(payload)   # float64[1024]
```

`extract_global_counts` and `extract_rna_counts` are unchanged
(they read `scan_trained`, not the calibration payload).

**Open**: the FL mixture (`_fl_mixture.py`) and EB Dirichlet
shrinkage (`_fl_empirical_bayes.py`) currently consume `int64[1024]`
histograms. The new aggregate is float64 fractional mass. Two
options:

- **6.3a** Round to `int64` at the boundary (cheap, preserves the
  existing pipeline contract, loses fractional precision near low
  counts).
- **6.3b** Update the FL mixture / EB layer to accept float64
  weights directly. This is the right answer long-term and is
  small in scope, but it is a behavior change that needs its own
  green tests.

Decision lock pending. Default for the cutover: 6.3b (float64 weights
end-to-end). If the FL EB changes get hairy, fall back to 6.3a and
file a follow-up.

### 6.4 `_diagnostics.py` (todo 10)

Replace the 8-mask `Diagnostics` with a fractional one:

```python
@dataclass(frozen=True, slots=True)
class Diagnostics:
    # Fragment-level counters (exact, integers; sum to n_total).
    n_observed: int
    n_excluded_multimap: int
    n_excluded_chimera: int
    n_excluded_artifact: int
    n_excluded_strand_ambig: int
    n_unobserved: int
    n_unannotated_ref: int          # subset of n_observed
    n_fl_unavailable: int           # subset of n_observed

    # Fractional mass summaries (float64).
    total_region_mass: float
    mass_by_coarse_class: Mapping[str, float]   # INTERGENIC, INTRON, EXON
    mass_by_compartment: Mapping[str, float]    # CONTAINED, BOUNDARY_LEFT, ...
    mass_by_splice: Mapping[str, float]
    mass_by_strand: Mapping[str, float]
    mass_by_signature: np.ndarray               # float64[16]

    # FL pool totals (float64).
    fl_pool_total: Mapping[str, float]          # 6 named pools
```

`total()` returns `n_observed + n_excluded_*` and is the accountability
sentinel. `to_summary_dict()` emits everything as nested JSON-safe
primitives for `summary.json`.

### 6.5 `_arrays.py` (todo 11)

`RegionArrays` keeps `ref_id`, `start`, `end`, `signature`,
`ref_offsets`, `order`, `n_refs`. Drop `type`, `strand`, `bf_left`,
`bf_right` from the *required* surface — they are no longer hot
inputs. If specific consumers still want them as joined columns,
they read them directly off `region_df` (which is small per locus).

`PayloadArrays` is replaced by a thin sorted-channel view:

```python
@dataclass(frozen=True, slots=True)
class PayloadArrays:
    region_counts_sorted: np.ndarray   # float32[R, 12], reordered by RegionArrays.order

    # Channel slices precomputed for hot loops.
    contained_unspliced_pos: np.ndarray  # float32[R]
    contained_unspliced_neg: np.ndarray
    boundary_left_unspliced_pos: np.ndarray
    boundary_left_unspliced_neg: np.ndarray
    boundary_right_unspliced_pos: np.ndarray
    boundary_right_unspliced_neg: np.ndarray

    # Aggregates used by regional/locus consumers.
    contained_unspliced_total: np.ndarray
    boundary_left_unspliced_total: np.ndarray
    boundary_right_unspliced_total: np.ndarray
```

Open: do we expose spliced channels at this layer or push that to
`FractionalEvidenceView`? Default for cutover: not in `PayloadArrays`
(keep this layer minimal and unspliced-focused, mirroring current
density consumers); add them later if a consumer needs them.

The temporary `_derive_signature_from_coarse(...)` bridge for old
unit tests is removed. Tests must construct v4-signature region
frames directly (helpers will live in a new `tests/_fractional_factories.py`,
see §8).

### 6.6 `density_global.py` / `_regional_exposure.py` / `locus_prior.py` (todo 12)

Per §5 (locked Option A), these three modules are reduced to skeletons:

- All top-level entry points (`compute_global_densities`,
  `RegionalGdnaExposure.build`, `estimate_locus_gdna`,
  `expected_gdna_count_global`, `assemble_priors`,
  `assemble_multilocus_prior`) raise `FractionalCutoverPending`
  immediately.
- The dataclass *shapes* are preserved (`GlobalGdnaDensity`,
  `GlobalDensityTable`, `RegionalGdnaExposure`, `LocusGdnaEstimate`,
  `MultiLocusPrior`, `PriorTable`) so `_result.CalibrationResult`
  type signatures and downstream type hints keep importing.
- All internal helpers that depended on the legacy payload schema
  (`_gdna_count_moment`, `_strand_identifiable_rows`,
  `kappa_opportunity_bp`, `precision_opportunity`,
  `strand_reliability_power`, `strand_correction_usable`, the
  `_orient`-based reweighting) are deleted from the public surface.
  If Phase 4 needs them, it will re-derive them against the
  fractional evidence interface.
- `_exposure.boundary_crossing_exposure(K)` is deleted outright (no
  `_legacy` namespace). `fractional_boundary_side_exposure(lengths_bp,
  gdna_fl)` from the parent plan §6 lands here even though no Phase 3
  consumer calls it — placing it now freezes the geometry interface
  for Phase 4.

This is intentionally aggressive: the smaller the surface that
survives the cutover, the less Phase 4 has to disentangle.

### 6.7 `_result.py` (todo 13)

- Drop `n_below_tolerance` and the K calibration semantics from
  `CalibrationResult`.
- Add `resolver_splicing_anchor_tolerance` (provenance only).
- Replace the `Diagnostics` instance with the §6.4 version.
- Per-locus / multi-locus diagnostic frame columns
  (`_PER_LOCUS_COLUMNS`, `_MULTI_LOCUS_COLUMNS`) are written by
  `locus_prior` consumers. Under Option A they will be unreachable
  until Phase 4; the column names themselves stay defined here so
  test imports keep working, but they will move when the new
  estimator lands.

### 6.8 `_exposure.py` updates (still todo 13)

- Delete `boundary_crossing_exposure(K)` from the public surface.
  Anything still importing it must move to
  `fractional_boundary_side_exposure(...)`.
- `contained_exposure_clipped`, `boundary_side_in_window`,
  `footprint_exposure_weight`, `transcript_exposure_weights`,
  `gdna_eff_len_for_loci` are unchanged in shape but their callers
  in `locus_prior` raise during Option A. Leave the helpers in
  place — they are correct and Phase 4 will need them.

### 6.9 `_orchestrator.py` (still todo 13)

`calibrate(...)` keeps its signature. Under the locked Option A:

1. Builds FL models (unchanged inputs, float64 fractional weights per
   §6.3b).
2. Builds the `Diagnostics` summary from the fractional payload.
3. Returns a `CalibrationResult` carrying payload, diagnostics, FL
   models, `resolver_splicing_anchor_tolerance` provenance, and
   `global_densities = None`, `regional_exposure = None`,
   `priors = None`.
4. `with_priors(...)` raises `FractionalCutoverPending`.

### 6.10 `pipeline.py` — `quant_from_buffer` refusal (still todo 13)

`quant_from_buffer` is updated to:

1. Run `scan_and_buffer` (unchanged).
2. Run `calibrate(...)` end-to-end (produces the fractional payload,
   FL models, and diagnostics).
3. Persist `summary.json` from the new diagnostics so users can
   inspect Phase 3 output on real data.
4. Raise `FractionalCutoverPending` before any EM / per-locus prior
   work, with a single typed message: *"gDNA density and per-locus
   priors are pending the Phase 4 estimator (see
   docs/fineregions/fractional_accumulator_python_cutover.md §5)."*

`rigel quant` catches `FractionalCutoverPending` at the CLI boundary
and exits with code 64 (`EX_USAGE`) and the typed message — no
stack trace.

### 6.11 `config.py` + `cli.py` (todo 14)

- Update `BamScanConfig.splicing_anchor_tolerance` docstring to say
  resolver-only.
- Update CLI help in `cli.py` for `--splicing-anchor-tolerance`.
- Keep the default (3) and the validation (`>= 0`).
- No new CLI flags. (The eventual Phase 4 estimator may add an
  `--estimator` switch; out of scope.)

## 7. Pipeline / `scan_payload.from_scan_dict` Wiring

After §6.1, the pipeline call path becomes:

```python
scan_dict = scanner.scan(...)
cal_dict = scan_dict["calibration"]    # new fractional keys
payload = CalibrationScanPayload.from_scan_dict(
    cal_dict, n_total=stats.n_read_names
)
```

`stats.n_read_names` (already produced by the scanner) is the
denominator for the balance assertion. The pipeline must continue to
pass it.

`install_calibration_regions(...)` in `pipeline.py` already drops K
(done as part of the native cutover). Verify no other call site still
passes it positionally.

## 8. Test Strategy

Run order (each step must be green before the next):

1. **Schema tests for `CalibrationScanPayload`**
   - Shape / dtype / contiguity validation.
   - Mass-identity tolerances.
   - Balance assertion.
   - `fl_pool_mass.sum(axis=1) == fl_pool_total`.
2. **FractionalEvidenceView**
   - Mask correctness over all 16 signatures.
   - Sense / antisense derivation table (§4) — exhaustive over
     `(sig, strand)`.
   - `gdna_fl_mass` equals the sum of the four non-exon pools.
3. **Native + Python integration**
   - End-to-end scan against an existing synthetic scenario, validate
     payload, view, and diagnostics dump.
   - Worker-merge equivalence: single-thread vs. multi-thread within
     float32 tolerance.
4. **Test factories**
   - Replace all legacy payload factories in:
     - `tests/test_calibrate_orchestrator.py`
     - `tests/test_regional_exposure.py`
     - `tests/test_assemble_priors.py`
     - `tests/test_locus_prior_fused.py`
     - `tests/test_calibration_result.py`
     - `tests/test_per_locus_gdna_mass.py`
     - `tests/test_density_global.py`
     - `tests/test_calibration_accumulator.py`
     - `tests/test_bayesian_prior_acceptance.py`
   - Land a shared helper module
     `tests/_fractional_factories.py` exporting
     `empty_payload(R)`, `payload_with_region_mass(...)`,
     `payload_with_pool(...)`, etc. against the new schema.
5. **Density / prior tests under Option A**
   - Convert tests that exercise `compute_global_densities`,
     `RegionalGdnaExposure.build`, `locus_prior.*` to assert
     `FractionalCutoverPending` is raised. Mark them
     `@pytest.mark.xfail(reason="Phase 4 estimator")` only if Option B
     is selected instead.
6. **CLI / config**
   - Help text shows resolver-only wording.
   - `BamScanConfig` round-trips the field unchanged.
7. **Goldens**
   - Goldens consuming `summary.json` fields that disappear
     (`n_below_tolerance`, integer mask counts) must be regenerated
     with `pytest --update-golden` and committed in the same change
     set; document each renamed key.

## 9. Ordering And Stop-Points

The implementation order locks in dependency order:

1. `scan_payload.py` (todo 7) — every other module imports it; nothing
   else compiles until this lands.
2. `signature.py` re-exports (already done).
3. `fractional_evidence.py` (todo 8).
4. `_fl_sources.py` (todo 9).
5. `_diagnostics.py` (todo 10).
6. `_arrays.py` (todo 11) and the `Option A` fail-fast wiring in
   `density_global` / `_regional_exposure` / `locus_prior` (todo 12).
7. `_result.py` + `_exposure.py` cleanup (todo 13).
8. `config.py` + `cli.py` docstring / help updates (todo 14).
9. Test factories + assertion rewrites (todo 15).
10. Full pytest + ruff and golden refresh (todo 16).

Useful stop-points to land independently if scope or time is tight:

- **S1**: 1–4 land green with their unit tests. Density/prior layer
  still raises `ModuleNotFoundError` for missing imports; mark as
  expected-fail. Buys us the new payload and FL extraction path.
- **S2**: through step 8. End-to-end scan + calibrate runs under
  Option A, producing payload + FL models + diagnostics. Quant runs
  block until Phase 4.
- **S3**: full green test suite with all factories rewritten.

## 10. Open Questions (After 2026-05-22 Locks)

The three items the user locked on 2026-05-22 are removed from this
list (orient cleanup → §3, Option A → §5, float64 FL weights → §6.3).
The remaining open items shape Phase 4, not this cutover:

1. **Mixed-signature EXON_INTRON boundary handling.** Parent plan §11
   Q3. Phase 3 surfaces them in diagnostics only.
2. **Spliced FL coverage.** Parent plan §11 Q4. Phase 3 leaves spliced
   fragments out of FL pools.
3. **First downstream estimator interface.** Parent plan §11 Q5.
   Recommendation: `FractionalEvidenceView` is the public interface;
   the raw payload is serialization-shaped only.
4. **Spliced channels in `PayloadArrays`.** §6.5. Default: exclude
   until Phase 4 needs them.
5. **Per-locus diagnostic columns.** §6.7. Recommendation: rename to
   `mass_*` and refresh goldens. Final names land with the Phase 4
   estimator.
6. **`extract_rna_counts` location.** §6.3. Default: keep co-located
   in `_fl_sources` for now.
7. **`boundary_kind_{left,right}` index columns.** Phase 3 does not
   read them. Default: keep through Phase 4; revisit once the new
   estimator is locked.
8. **`tests/test_v5_*` source files.** Confirm none remain as `.py`
   sources; delete leftover `__pycache__` entries during the cutover.
9. **`PayloadArrays.from_payload` reordering.** With float32
   `region_counts` the `[order, :]` copy is ≈ 48 B/region, so ≈ 480
   MB once at R ≈ 10⁷. Watch item; if it bites, reorder *queries*
   instead of the payload.

## 11. Acceptance Criteria

Phase 3 Python cutover is done when:

- All Python imports of legacy payload symbols are gone
  (`grep -RIn 'fl_hist\|per_region_counts\|u_left\|u_right\|by_orient\|MASK_EXON\|MASK_INTRON\|MASK_INTERGENIC\|n_below_tolerance' src/rigel/` returns nothing).
- `pytest tests/ -v` is green under the chosen Option (A or B). Tests
  that document Phase 4 dependency are clearly named and `xfail`ed
  with a reference to this document.
- `ruff check src/ tests/` and `ruff format --check src/ tests/` are
  clean.
- `summary.json` from a calibration run contains the new diagnostics
  surface and no legacy keys.
- `docs/fineregions/fine_region_migration_phase0_log.html` records
  the completion row, including the chosen §5 option and §6.3 option.
- This document gets a final "Status: complete" change in its header.
