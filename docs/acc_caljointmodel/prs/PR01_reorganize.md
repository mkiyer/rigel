# PR 1 — Reorganize + region geometry & strand annotation *(Step 2: code cleanup and organization)*

**Parent plan:** [`../00_implementation_plan.md`](../00_implementation_plan.md) §7.
**Type:** Python-only (index rebuild required — schema bump on `regions.feather`).
**Build required:** no (C++ untouched).
**Status:** Decisions resolved (§7). Magic-number budget: **0 new
constants** (the signature bits + strand-class values are an encoding, not
tunables). Ready to implement on your go-ahead.

## Goal

Put the calibration package and the region-geometry substrate into the
shape the calibrator needs, *before* any inference is written. Three
things land:

1. **Package layout** — break the calibrator surface out of the
   `__init__.py` stub into purpose files (`result.py`, `substrate.py`,
   `calibrate.py`, `priors.py` as skeletons), so PRs 2–6 fill bodies
   without reshuffling.
2. **Region geometry + strand annotation** — restore the per-region
   `signature` (4-bit exon/intron × strand) and the `RegionArrays`
   sorted-CSR geometry view + boundary↔region index mapping. This is the
   **D2 prerequisite**: the strand-balance model needs a per-region
   strand class (`+ / − / ambiguous`), which the current minimal
   `regions.py` dropped during the burn.
3. **Seam repair** — wire `calibrate(...)` with its real inputs (still
   raising), and add an index-alignment guard so the accumulator's region
   order is provably identical to `region_df`.

No inference; `rigel quant` still aborts at the `calibrate()` stub. This
is the "organization" half of the two-step cleanup.

> All recovery is from `fc96902` ("new calib system implemented") per the
> [Recovered-Code Map](../00_implementation_plan.md) §5, scrubbed of the
> obsolete pieces noted below.

## Background — what exists vs what we recover

- **Live & kept:** `AccumulatorPayload` ([scan_payload.py](../../../src/rigel/scan_payload.py)),
  the native accumulator, the minimal `regions.py` partition builder,
  `index.region_df` / `index.ref_name_to_id` / `index.ref_names`.
- **Recover (adapt):** `RegionArrays` (geometry) and the region
  `signature` computation; `transcript_strand_class` (signature → strand
  class); the boundary↔region index mapping.
- **Recover but SCRUB:** the pre-burn `signature.py` also defines a
  *12-channel* layout (`N_COMPARTMENTS×N_SPLICE×N_STRAND`) and 6 FL pools.
  Those are **obsolete** — the accumulator emits a 4-channel
  `region_contained` + separate boundary arrays, and there is no FL. We
  recover **only** the region-signature half (`BIT_*`, `pack_signature`,
  `coarse_strand_from_signature`, `is_ambiguous_signature`,
  `transcript_strand_class`).
- **Deferred to PR 2:** the count/mass substrate views (the calibrator-
  facing `CalibrationSubstrate` with channel reductions + boundary→region
  mass attribution). PR 1 provides only the geometry + index mapping they
  build on.

## Target package layout (end of PR 1)

```
src/rigel/calibration/
  __init__.py        # re-export CalibrationConfig, CalibrationResult, calibrate,
                     #   CalibrationSubstrate(skel), CalibrationSubstrateError,
                     #   CalibrationConvergenceError, assemble_priors(skel)
  regions.py         # EDIT: add per-region `signature` (+ strand class) to the
                     #   partition (see Q1/Q2)
  signature.py       # NEW (recovered, scrubbed): 4-bit region signature bits,
                     #   pack_signature, coarse_strand/type, is_ambiguous,
                     #   transcript_strand_class
  region_arrays.py   # NEW (recovered): RegionArrays.from_region_df + boundary↔
                     #   region index mapping (left_region_index/right_region_index)
  strand_summary.py  # KEEP
  result.py          # NEW skeleton: CalibrationResult (empty→real schema in PR 2)
  substrate.py       # NEW skeleton: CalibrationSubstrate (real adapter in PR 2)
  calibrate.py       # NEW skeleton: calibrate() (placeholder in PR 2)
  priors.py          # NEW skeleton: assemble_priors() (real consumer in PR 6)
  errors.py          # NEW: CalibrationSubstrateError, CalibrationConvergenceError
```

The empty placeholder `CalibrationResult` and the `calibrate` stub move
out of `__init__.py` into `result.py` / `calibrate.py`; `__init__.py`
becomes pure re-exports.

## Tasks

### T1 — Package skeletons (organization)

- Create `result.py`, `substrate.py`, `calibrate.py`, `priors.py`,
  `errors.py` with correct signatures + docstrings and
  `NotImplementedError` bodies (or empty dataclasses). Move the inline
  `CalibrationResult`/`calibrate` out of `__init__.py`.
- `errors.py`: `CalibrationSubstrateError(ValueError)`,
  `CalibrationConvergenceError(RuntimeError)` (doc 04 §9).
- `__init__.py`: re-export the public surface (no logic).
- `calibrate.py` skeleton signature (firmed in PR 2):
  `def calibrate(*, payload, region_arrays, strand_model, config) -> CalibrationResult: raise NotImplementedError(...)`.

### T2 — Recover `signature.py` (region signature only)

`git show fc96902:src/rigel/calibration/signature.py` → keep:
`BIT_INTRON_POS/NEG`, `BIT_EXON_POS/NEG`, `pack_signature`,
`validate_signature`, `coarse_type_from_signature`,
`coarse_strand_from_signature`, `is_ambiguous_signature`. Drop the
12-channel `channel_index`/compartment/FL-pool constants (obsolete under
the 4-channel accumulator). Add `transcript_strand_class` (from
`fc96902:fractional_evidence.py`; TS_NONE=0/TS_POS=1/TS_NEG=−1/TS_AMBIG=2).

### T3 — The merged-signature region partition

**Region definition (locked, Q1).** Every genome position carries a 4-bit
signature {`exon_pos`, `exon_neg`, `intron_pos`, `intron_neg`}. A *region*
is a maximal non-overlapping interval over which the signature is constant;
by construction **a region and each of its neighbours have different
signatures** (adjacent equal-signature segments MUST be merged). Boundaries
therefore sit exactly at signature transitions.

Recover the event-sweep builder from `fc96902:regions.py`
(`_add_interval_events` → sort events → sweep emitting one segment per
inter-event gap, with the merge rule
`rows[-1].end == start and rows[-1].signature == sig → extend`).
Per-transcript events add `BIT_EXON_{POS,NEG}` over exons and
`BIT_INTRON_{POS,NEG}` over introns.

**Minimal schema (Q2):**
`REGION_COLUMNS = [region_id, ref_name, start, end, length, signature]` —
persist `signature` only; derive `type` / `strand` / `ts_class` on load
from the signature. Drop the fc96902 `left/right_signature`,
`boundary_kind`, `boundary_flux_*`, and per-bit bool columns (not needed).
Update the dtype map and `validate_against_ref_lengths` (add the
neighbour-differs invariant).

- `index.py`: regions.feather gains `signature`; **force index rebuild,
  no backward compat (Q3).**
- The merged partition changes the region set, so re-baseline
  `tests/test_scanner_accumulator_integration.py` expected values.

### T4 — `RegionArrays` + boundary↔region mapping

`git show fc96902:src/rigel/calibration/_arrays.py` → recover
`RegionArrays.from_region_df(region_df, ref_name_to_id)` into
`region_arrays.py` (geometry: `ref_id/start/end/signature/ts_class/
region_size_bp/ref_offsets/order/n_refs`). Drop `PayloadArrays` here —
it reads the old 12-channel payload and is superseded by PR 2's
`CalibrationSubstrate`.

Add the **boundary↔region index mapping** (recovered logic from
`fc96902:boundaries.py` `left_region_index()`/`right_region_index()`),
computed from the `AccumulatorPayload` topology offsets
(`ref_region_offsets`, `ref_boundary_offsets`): for region `r` in ref `f`
with local index `i`, `left_boundary(r) = ref_boundary_offsets[f] + i`,
`right_boundary(r) = + i + 1`; and the inverse (`left_region_of(b)`,
`right_region_of(b)`). This is the substrate the PR-2 `CalibrationSubstrate`
uses to attribute `boundary[r].mass_right` + `boundary[r+1].mass_left`
to region `r` (per [§4 D1](../00_implementation_plan.md)).

### T5 — Seam repair + index-alignment guard

- `pipeline.py::run_pipeline`: build
  `region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)`
  and wire
  `calibrate(payload=calibration_payload, region_arrays=region_arrays, strand_model=strand_models, config=config.calibration)`
  (still raises). All four inputs are live locals at the call site today.
- **Index-alignment guard:** assert `len(region_arrays.start) == payload.r_total`
  and per-ref `ref_region_offsets` agree, so calibration arrays (ref-major)
  line up 1:1 with `region_df`. Raise `CalibrationSubstrateError` on mismatch.

### T6 — Tests (`tests/calibration/`)

- `test_signature.py` — `pack_signature` / `coarse_strand` /
  `transcript_strand_class` on hand-built exon/intron geometry incl.
  overlap (AMBIG) and intergenic (NONE). (Recover asserts from the
  fc96902 signature/region tests where useful.)
- `test_region_arrays.py` — `from_region_df` CSR ordering, `ts_class`,
  `ref_offsets`; boundary↔region mapping round-trips across ref seams.
- `test_region_index_alignment.py` — `R_total == len(region_df)` and
  per-ref offset equality on a synthetic index (parent §8.3).
- `test_package_surface.py` — `rigel.calibration` exposes exactly the
  intended names; skeletons raise `NotImplementedError`, not `AttributeError`.

## Acceptance gate

- [ ] `import rigel`, `pytest --collect-only tests/` clean.
- [ ] `tests/calibration/` (new) pass; substrate suite still green.
- [ ] A freshly built index's `regions.feather` carries `signature`;
      `RegionArrays.from_region_df` consumes it; `ts_class` matches the
      transcript geometry.
- [ ] `rigel quant` reaches `calibrate(payload, region_arrays, strand_model, config)`
      (real args) and aborts there with `NotImplementedError` — alignment
      guard passes en route.
- [ ] `ruff check src/ tests/` clean.
- [ ] No new full-suite failure modes (still only `calibrate()` aborts).

## Rollback

Revert the PR. Index rebuild reverts to the minimal `regions.feather`.

## 7. Resolved decisions (your critique, 2026-05-29)

- **Q1 — Partition granularity: MERGE (locked).** Region = maximal
  constant-4-bit-signature interval; a region and its neighbours must have
  different signatures (adjacent equal-signature segments merged). See T3.
- **Q2 — Intergenic: include in BOTH density and strand models.**
  Intergenic = signature `0` = strand `NONE`. For the strand model,
  `NONE`-strand regions get an **arbitrary (fixed) sense/antisense
  assignment** (fair: gDNA is unstranded, ~50/50). Caveat: unannotated
  transcription (incl. antisense inside annotated regions) injects
  strand/density outliers → **overdispersion**, which the model must
  absorb. The robustness mechanism is a **PR 3/PR 5 design discussion**
  (note below). PR 1 only tags `NONE` correctly and does not exclude
  intergenic; it persists `signature` only and derives `type`/`strand`/
  `ts_class` on load.
- **Q3 — Index schema bump: yes, force rebuild, no backward compat.**
- **Q4 — PR1/PR2 split: yes** (geometry in PR 1; the count/mass
  `CalibrationSubstrate` in PR 2).
- **Q5 — `calibrate()` signature firms up across PR1→PR2: agree.**
- **Q6 — Clean / concise / organized; ≤25 modules, ≤~8 constants.**
  **Hard rule: pause and discuss before introducing ANY new magic number /
  heuristic / parameter.** PR 1 introduces **zero** new tunables. Recorded
  in agent memory and the master-plan conventions.

### Robustness vs outlier-clipping (Q2 follow-up — PR 3 / PR 5 discussion)

The legacy seeding clipped the top-X% density / strand-imbalanced regions —
a magic-number cliff the burn removed. The v6 **principled substitute is
fitted overdispersion**: the NB count dispersion `φ` and the Beta-Binomial
strand dispersions `ρ_d` / `ρ_r` give heavy-tailed likelihoods that
down-weight outliers *without* a percentile cliff, and there is **no
two-state expressed/unexpressed latent** to seed wrongly. Whether
overdispersion alone suffices — vs a principled robustification such as
median / Huber moment warm-starts (still constant-free) — is to be settled
**empirically at PR 5** with benchmarks. Per Q6, **no clip constant will be
introduced without discussion.** Flagged here so it is not forgotten.

## 8. Execution notes (implemented 2026-05-29)

Implemented on branch `pr1-reorganize`. **Magic-number budget: 0 new tunables**
(signature bits + strand-class integers are an encoding). Deviations and
decisions worth recording:

- **Boundary terminal mapping — natural & round-tripping (refinement of the
  recovered logic).** `fc96902:boundaries.py` set *both* sides of every
  reference-terminal boundary to `-1`. The recovered
  `boundary_region_indices` instead attributes the one real adjacent region
  to a terminal (only the off-edge side is `-1`): the start terminal's right
  region is region 0, the end terminal's left region is the last region. This
  makes the forward/inverse maps exact inverses on internal seams and a clean
  one-sided attribution at terminals, so the round-trip test holds. It is
  numerically harmless (terminal boundaries carry ~zero crossing mass) and is
  the **more correct** input for the PR 3 D2 strand model, which orients each
  boundary side to its adjacent region's strand.
- **`RegionType` / `RegionStrand` moved into `signature.py`.** At `fc96902`
  these lived in `regions.py` and `signature.py` imported them back (a
  circular dependency worked around with function-local imports). They are an
  encoding, so they now live in `signature.py`, giving a clean one-way
  dependency `regions.py → signature.py` (leaf, imports only numpy + enum).
- **Dropped `has_intron_flag` / `has_exon_flag`** from `signature.py` (unused,
  not in the T2 keep-list).
- **Index rebuild forced via `INDEX_FORMAT_VERSION` 4 → 5** (Q3). Every test
  builds its index fresh, so the bump introduced **no** new failure mode.
- **`test_scanner_accumulator_integration.py` needed no re-baseline.** Its
  expectations are derived from the index at runtime (not hard-coded), so the
  merged partition is absorbed automatically; it passes unchanged.
- **One extra test file beyond the T6 list:** `tests/calibration/test_regions.py`
  directly covers the merge invariant + signature assignment of the new
  builder (the riskiest new logic), with worked exon/intron geometry.

### Acceptance gate — results

- `import rigel` OK; `pytest --collect-only tests/` → 891 collected, clean.
- `tests/calibration/` (24) + integration (3) → **27 passed**.
- Fresh index `regions.feather` columns =
  `[region_id, ref_name, start, end, length, signature]` (uint8);
  `RegionArrays.from_region_df` consumes it; `ts_class` matches geometry.
- `rigel quant` reaches
  `calibrate(payload=…, region_arrays=…, strand_model=…, config=…)` with live
  args and aborts there (`NotImplementedError`); the alignment guard passes
  en route.
- `ruff check src/ tests/` clean.
- Full suite: **248 failures/errors, every one `NotImplementedError`** from the
  `calibrate` / `quant_from_buffer` stubs — identical failure mode and count
  to the PR 0 baseline. No new modes (no `RuntimeError` from the version bump,
  no import/attribute errors).
