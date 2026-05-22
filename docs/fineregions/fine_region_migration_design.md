# Fine-Grained Region Migration Design

Date: 2026-05-21

## Purpose

Rigel currently persists a coarse calibration region partition with three region types:
`INTERGENIC`, `INTRON`, and `EXON`, plus a collapsed strand code. This loses important
local structure in loci with overlapping genes, overlapping isoforms, long single-exon
transcripts, and read-through transcripts. The goal of this migration is to restore a
fine-grained, strand-aware, non-overlapping genome partition and to make BAM scan
calibration evidence granular enough to improve gDNA, nRNA, and mRNA deconvolution.

This migration is intentionally not backward compatible. Index format, scan payloads,
calibration data structures, tests, and golden outputs should all move together.

## Historical Audit

The old fine-region code is in git commit `14c307f`.

Useful pieces from that commit:

- `14c307f:src/rigel/index.py`, `build_region_table(...)`
  - Built a reference-complete non-overlapping region table.
  - Collected boundaries from transcript starts, transcript ends, exon starts, exon ends,
    and reference boundaries.
  - Formed atomic bins between sorted boundaries, assigned boolean flags, then merged
    adjacent bins with identical signatures.
  - Persisted columns: `region_id`, `ref`, `start`, `end`, `length`, `exon_pos`,
    `exon_neg`, `tx_pos`, `tx_neg`.

- `14c307f:src/rigel/native/resolve_context.h`, `RegionAccumulator`
  - Accumulated float64 fractional per-region counts in shape `[n_regions, 4]`.
  - Columns were `unspliced_pos`, `unspliced_neg`, `spliced_pos`, `spliced_neg`.
  - Fractional weighting used `overlap_bp / total_aligned_bp / num_hits`, which kept a
    fragment's total contribution at one unit across all blocks and regions.
  - Fragment length observations were emitted only for unique, unspliced, single-region
    fragments.

- `14c307f:tests/test_regions.py`
  - Good tests to resurrect: reference-complete tiling, empty references, globally
    sequential `region_id`, no gaps or overlaps, adjacent identical-signature merging,
    opposite-strand overlap flags.

- `14c307f:tests/test_region_acc_nh_weighting.py`
  - Useful fractional-count invariant tests, but current calibration excludes
    multimappers. Keep the one-unit fractional-count invariant; do not reintroduce
    multimapper calibration evidence in this migration unless we deliberately design it.

Pieces not to reuse directly:

- The historical schema used `tx_pos` and `tx_neg`, not explicit intron flags. That is
  insufficient for the new requested 16-state model if a same-strand exon from one
  transcript overlaps an intron from another transcript. Persist explicit intron flags.
- The historical accumulator did not separate contained support from left/right boundary
  flux.
- The historical accumulator used a cgranges-backed region index. Current main has a
  better native per-reference CSR `RegionIndex` and a better scanner accountability
  model in `src/rigel/native/calibration/`.
- The old calibration EM is not the target architecture. Treat it as background, not as
  code to revive.

## Current-Code Audit

Current region and calibration flow:

- [src/rigel/calibration/regions.py](../../src/rigel/calibration/regions.py)
  defines `RegionType`, `RegionStrand`, `RegionRecord`, and `emit_regions(...)`. This is
  the core coarse-region implementation to replace.

- [src/rigel/index.py](../../src/rigel/index.py)
  builds both `intervals.feather` and `regions.feather` in `build_index_artifacts(...)`.
  The interval table is for transcript resolution and can stay mostly separate; the
  region table should move to fine signatures.

- [src/rigel/pipeline.py](../../src/rigel/pipeline.py)
  wires regions into the native scanner in `_wire_calibration_regions(...)`. It currently
  derives a 3-bit coarse type mask from `RegionType` and passes collapsed strand codes to
  C++.

- [src/rigel/native/calibration/region_index.h](../../src/rigel/native/calibration/region_index.h)
  stores starts, ends, type masks, and strand codes. This should become the native holder
  for fine signatures and any derived coarse fields.

- [src/rigel/native/calibration/accumulator.cpp](../../src/rigel/native/calibration/accumulator.cpp)
  currently emits integer per-region mask counts, integer boundary flux, and orientation
  counters. This should be replaced or substantially refactored into a fractional fine
  accumulator.

- [src/rigel/calibration/scan_payload.py](../../src/rigel/calibration/scan_payload.py)
  validates the current integer payload. This should become the schema lock for the new
  fine-region scan payload.

- [src/rigel/calibration/density_global.py](../../src/rigel/calibration/density_global.py)
  computes current global gDNA densities from coarse classes. The first fine-region
  calibration implementation should be able to reproduce equivalent coarse densities via
  a compatibility adapter, then incrementally use the richer evidence.

- [src/rigel/calibration/_regional_exposure.py](../../src/rigel/calibration/_regional_exposure.py)
  learns per-region exposure weights from coarse evidence. This is the main downstream
  consumer that should benefit from fine signatures and boundary flux.

- [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py)
  consumes `RegionArrays`, `PayloadArrays`, local region overlaps, regional exposure, and
  gDNA effective lengths. It should not know about old coarse schema details after this
  migration.

Cleanup opportunities before implementation:

- Remove `set_regions_legacy(...)` from native scanner bindings.
- Remove optional legacy-region behavior from `scan_and_buffer(...)`; a missing fine
  region table should be an index error.
- Stop using `tx_pos`/`tx_neg` as a proxy for introns in calibration code. If a transcript
  span flag is needed, derive `tx_pos = intron_pos | exon_pos` and
  `tx_neg = intron_neg | exon_neg`.
- Consolidate region column naming around `ref_name`, not historical `ref`.
- Replace `RegionType`/`RegionStrand` as persisted semantics with fine `signature` plus
  explicit derived coarse fields.
- Ensure region IDs are stable globally in FASTA reference order, then start order.
- Audit `_wire_calibration_regions(...)` for references that have no transcripts. Fine
  regions are reference-complete, so scanner region ref IDs should follow the index/BAM
  reference mapping rather than silently dropping intergenic-only references.

## Canonical Region Signature

The persisted fine-region signature is a 4-bit unsigned integer:

```text
bit 3: intron_pos
bit 2: intron_neg
bit 1: exon_pos
bit 0: exon_neg
```

With this canonical layout:

```text
signature = (intron_pos << 3) | (intron_neg << 2) | (exon_pos << 1) | exon_neg
```

The TODO text contains one example that labels `b0100` as `intron_pos`. Under the
canonical bit layout above, `b0100` is `intron_neg` and `b1000` is `intron_pos`.
The coarse-state table in the TODO is consistent with the canonical layout, so this
design adopts the canonical layout and treats the prose example as a swapped-bit typo.

Persist explicit booleans as well as the packed signature for readability and testability:

- `intron_pos: bool`
- `intron_neg: bool`
- `exon_pos: bool`
- `exon_neg: bool`
- `signature: uint8`

Derived helpers:

```text
has_pos = intron_pos | exon_pos
has_neg = intron_neg | exon_neg
tx_pos = has_pos
tx_neg = has_neg
coarse_type = EXON if (exon_pos | exon_neg), else INTRON if (intron_pos | intron_neg), else INTERGENIC
coarse_ambig = has_pos & has_neg
coarse_strand = NONE, POS, NEG, or AMBIG from has_pos/has_neg
```

Coarse derivation table:

| Signature | Fine state | Coarse type | Ambig |
| --- | --- | --- | --- |
| `b0000` | intergenic | INTERGENIC | no |
| `b0001` | exon_neg | EXON | no |
| `b0010` | exon_pos | EXON | no |
| `b0011` | exon_pos + exon_neg | EXON | yes |
| `b0100` | intron_neg | INTRON | no |
| `b0101` | intron_neg + exon_neg | EXON | no |
| `b0110` | intron_neg + exon_pos | EXON | yes |
| `b0111` | intron_neg + exon_pos + exon_neg | EXON | yes |
| `b1000` | intron_pos | INTRON | no |
| `b1001` | intron_pos + exon_neg | EXON | yes |
| `b1010` | intron_pos + exon_pos | EXON | no |
| `b1011` | intron_pos + exon_pos + exon_neg | EXON | yes |
| `b1100` | intron_pos + intron_neg | INTRON | yes |
| `b1101` | intron_pos + intron_neg + exon_neg | EXON | yes |
| `b1110` | intron_pos + intron_neg + exon_pos | EXON | yes |
| `b1111` | all four flags | EXON | yes |

## Fine Region Table Schema

Bump `INDEX_FORMAT_VERSION` and require the new schema at load time.

Recommended `regions.feather` columns:

| Column | Type | Notes |
| --- | --- | --- |
| `region_id` | int64 | Globally sequential in FASTA reference order, then start order. |
| `ref_name` | string | Reference name from FASTA. |
| `start` | int64 | 0-based inclusive. |
| `end` | int64 | 0-based exclusive. |
| `length` | int64 | `end - start`; persisted to make diagnostics and tests simple. |
| `signature` | uint8 | Canonical 4-bit fine state. |
| `intron_pos` | bool | Explicit flag. |
| `intron_neg` | bool | Explicit flag. |
| `exon_pos` | bool | Explicit flag. |
| `exon_neg` | bool | Explicit flag. |
| `coarse_type` | uint8 | Derived convenience, not compatibility. |
| `coarse_strand` | uint8 | Derived convenience, not compatibility. |
| `coarse_ambig` | bool | Derived convenience. |

Potential optional columns for later calibration design:

- `left_signature`, `right_signature`: adjacent region signatures within the same ref,
  with sentinel values at reference ends. Useful for boundary-class filters without
  repeated neighbor lookups.
- `left_boundary_kind`, `right_boundary_kind`: compact flags for exon/intron/TSS/TES
  boundary provenance. This may help distinguish terminal transcript boundaries from
  internal exon-intron boundaries.

Do not persist old `tx_pos_bp`, `tx_neg_bp`, `exon_pos_bp`, `exon_neg_bp`, or old
`boundary_flux_left/right` as primary schema. If bp summaries are needed, derive them
from fine flags and `length`.

## Phase 1: Index Builder

Implement `build_fine_region_table(transcripts, ref_lengths) -> pd.DataFrame` in
[src/rigel/calibration/regions.py](../../src/rigel/calibration/regions.py) or a new
small helper module imported by [src/rigel/index.py](../../src/rigel/index.py).

Algorithm per reference:

1. Group real, non-synthetic transcripts by reference.
2. For each transcript, sort exons by coordinate.
3. Add exon interval events for each exon:
   - positive transcript: `exon_pos += 1` at exon start, `exon_pos -= 1` at exon end.
   - negative transcript: `exon_neg += 1` at exon start, `exon_neg -= 1` at exon end.
4. Add intron interval events for each gap between consecutive exons of the same
   transcript:
   - positive transcript: `intron_pos += 1` at intron start, `intron_pos -= 1` at intron end.
   - negative transcript: `intron_neg += 1` at intron start, `intron_neg -= 1` at intron end.
5. Include reference boundaries `0` and `ref_length`.
6. Sweep boundaries in ascending order. Apply events at boundary `b`, then the state to
   the right of `b` applies to `[b, next_boundary)`.
7. Emit a region only when the 4-bit signature changes from the previous emitted segment.
   This preserves the old adjacent-identical-bin merge behavior.
8. Empty references emit exactly one `b0000` intergenic region covering the whole ref.
9. Validate: positive length, full coverage, no gaps, no overlaps, region IDs equal row
   index, every signature in `[0, 15]`, and derived coarse fields match the signature.

Important detail: explicit intron interval events are required. `tx_pos && !exon_pos` is
not sufficient when a same-strand exon from one transcript overlaps an intron from another
transcript.

Tests to add or revive:

- Empty-reference single intergenic region.
- Multi-reference full coverage.
- Adjacent same-signature regions merge.
- Opposite-strand exons create `b0011`.
- Same-strand exon over intron creates `b1010` or `b0101`.
- Opposite-strand exon over intron creates `b1001` or `b0110` and `coarse_ambig=True`.
- Read-through transcript spanning two genes does not collapse all exons into one coarse
  region unless signatures are truly identical.

## Phase 2: BAM Scan Fine Accumulator

Replace the current integer mask accumulator with a fractional fine accumulator.

The new scan payload should store one float64 count unit per eligible fragment, split over
region compartments:

```text
compartment: contained, boundary_left, boundary_right
splice:      unspliced, spliced
strand:      pos, neg
```

Recommended native layout:

```text
region_counts: float64[R, 3, 2, 2]
```

or a flattened equivalent with named accessors in Python:

```text
contained_unspliced_pos
contained_unspliced_neg
contained_spliced_pos
contained_spliced_neg
boundary_left_unspliced_pos
boundary_left_unspliced_neg
boundary_left_spliced_pos
boundary_left_spliced_neg
boundary_right_unspliced_pos
boundary_right_unspliced_neg
boundary_right_spliced_pos
boundary_right_spliced_neg
```

Fragment observation policy:

- Keep current calibration exclusion gates for the first implementation:
  - exclude multimappers from calibration evidence;
  - exclude chimeric fragments;
  - exclude splice artifacts;
  - observe at most one eligible non-chimeric hit per molecule.
- Continue emitting accountability counters for excluded and unobserved fragments.
- Keep fragment length model training separate and conservative.

Fractional accumulation algorithm:

1. Input is the chosen fragment's compatible aligned blocks `(ref_id, start, end)` and
   fragment metadata: splice status and strand.
2. Reject mixed-reference or mixed-strand observations for calibration.
3. Let `total_aligned_bp = sum(end - start for all compatible blocks)`. This denominator
   is required for the invariant that a fragment contributes exactly one unit across all
   regions and compartments. Using per-block denominators would make a multi-block spliced
   fragment contribute one unit per block, which conflicts with the one-fragment-one-unit
   invariant.
4. For each aligned block, query overlapping regions in the native CSR `RegionIndex`.
5. If the block overlaps one region, add `block_len / total_aligned_bp` to that region's
   `contained` compartment.
6. If the block overlaps multiple regions, for each overlapping region:
   - `overlap_bp = min(block_end, region_end) - max(block_start, region_start)`.
   - `base_weight = overlap_bp / total_aligned_bp`.
   - `cross_left = block_start < region_start && block_end > region_start`.
   - `cross_right = block_start < region_end && block_end > region_end`.
   - `n_crossed = cross_left + cross_right`.
   - Add `base_weight / n_crossed` to each crossed side.
7. Sum of all emitted contained and boundary contributions for an observed fragment must be
   `1.0` within floating-point tolerance.

The examples in the TODO follow from this algorithm:

- Block `[80, 140)` over `[0,100)`, `[100,150)` contributes `20/60` to R1 right and
  `40/60` to R2 left.
- Block `[50,250)` over `[0,100)`, `[100,150)`, `[150,300)` contributes `50/200` to R1
  right, `50/200/2` to R2 left, `50/200/2` to R2 right, and `100/200` to R3 left.

Native implementation notes:

- Keep `RegionIndex::overlap_into(...)`; add `signature(rid)`, `coarse_type(rid)`, and
  maybe neighbor signature accessors.
- Store fine signatures in native region index, not only derived coarse masks.
- Replace current `int64_t` count vectors with `double` vectors for fractional evidence.
- Keep integer accountability counters separate from fractional region evidence.
- Consider tracking `n_fractional_unit_violations` or debug-only assertions for the
  one-unit invariant.
- Keep all arrays row-major and contiguous for nanobind export.

Python payload updates:

- Replace `CalibrationScanPayload.per_region_counts[:, MASK_*]` with a typed
  `FineCalibrationScanPayload` or a revised `CalibrationScanPayload` with explicit
  compartment arrays.
- `PayloadArrays` should become a thin sorted-order view of the 12 fine count columns.
- Remove old 8-state mask validation as the primary schema. If fragment-level masks are
  retained for FL histograms, keep them in a separately named block.

Tests to add:

- Single contained block contributes exactly one unit to contained.
- Boundary crossing examples from TODO.
- Two-block spliced fragment contributes exactly one unit total, not two.
- Left/right flux splits correctly for middle regions crossing both boundaries.
- Strand and splice channels route to the expected column.
- Excluded multimapper/chimera/artifact accountability still balances `n_read_names`.

## Phase 3: Fragment Length Model

Keep the first implementation conservative:

- Continue training RNA FL from scanner-trained resolved RNA observations.
- Continue training gDNA FL from gDNA-compatible unique observations.
- Define gDNA-compatible using fine-region evidence, not old coarse persisted types.

Initial gDNA-compatible categories:

- Pure intergenic contained evidence: signature `b0000`.
- Pure intronic contained evidence: signatures with intron bits and no exon bits.
- Exon-intron boundary evidence: boundary flux touching a transition where one side has
  exon evidence and the adjacent side has intron evidence.

Implementation approach:

1. Keep a fragment-level FL histogram block in the scan payload if it remains the simplest
   way to build `fl_models.gdna`.
2. Derive its categories from fine signatures and boundary transitions.
3. Only after Phase 2 tests pass, evaluate whether the older v4-style region-weighted FL
   model is worth reviving.

Do not let the FL model migration block the region schema and accumulation migration. The
first target is parity with current conservative gDNA FL behavior, then improvement.

## Phase 4: Calibration Design Direction

Phase 4 is not fully specified yet, but Phase 2 should produce enough evidence for these
models:

1. Coarse-parity adapter
   - Reconstruct current-style `INTERGENIC`, `INTRON`, `EXON-INTRON`, and
     `EXON-CONTAINED` numerators from fine counts.
   - Use this adapter as a regression bridge while removing old persisted coarse schema.

2. Per-signature density table
   - Estimate local/global gDNA density by fine signature and compartment.
   - Separate contained evidence from boundary flux evidence.
   - Treat gDNA as strand-symmetric and use strand imbalance to downweight RNA-like
     contained exonic evidence.

3. Boundary-transition model
   - Use left/right flux with adjacent signatures to estimate how much gDNA density should
     project across exon-intron and exon-intergenic transitions.
   - This can recover today's boundary-flux behavior and extend it to finer transitions.

4. Uncertainty-aware regional exposure
   - Tiny regions and boundary-only regions need stronger shrinkage.
   - Use effective opportunity and overdispersion to compute a per-region reliability or
     posterior precision.
   - Feed this into regional exposure weights and later into the gDNA prior strength.

5. Prior and EM integration
   - `expected_gdna_count_global(...)` and `RegionalGdnaExposure.build(...)` should consume
     fine payload arrays directly.
   - `gdna_eff_len_for_loci(...)` can remain geometric at first, but exposure weighting
     should use fine regions rather than coarse exon unions.
   - The gDNA prior should eventually carry uncertainty, not just a hard pseudocount.

## Implementation Sequence

### Phase 0: Cleanup

Deliverables:

- Remove native `set_regions_legacy(...)` and any missing-region fallback.
- Rename current coarse-region concepts so future diffs are obvious:
  - `RegionType` -> `CoarseRegionType` if retained as derived convenience.
  - `RegionStrand` -> `CoarseRegionStrand` if retained as derived convenience.
- Add a small pure-Python signature helper module or functions:
  - `pack_signature(...)`
  - `coarse_type_from_signature(...)`
  - `coarse_strand_from_signature(...)`
  - `coarse_ambig_from_signature(...)`
- Add tests for all 16 signature derivations before touching native code.
- Decide whether region ref IDs should use `index.ref_name_to_id` everywhere in the
  scanner. Do this before Phase 2 to avoid losing intergenic-only references.

### Phase 1: Fine Index

Deliverables:

- `INDEX_FORMAT_VERSION = 4` or higher.
- New fine `regions.feather` schema.
- Strict load-time schema validation.
- Region builder tests revived from `14c307f` and extended for explicit intron flags.
- No legacy load path.

### Phase 2: Fine BAM Accumulation

Deliverables:

- Native region index stores fine signatures.
- Native accumulator emits fractional 12-channel region evidence.
- Python payload dataclass validates float64 shapes and accountability counters.
- Unit tests for contained and boundary flux examples.
- Native rebuild required after C++ changes:
  `conda activate rigel && pip install --no-build-isolation -e .`

### Phase 3: gDNA FL Parity

Deliverables:

- Define fine-signature-derived gDNA FL source categories.
- Keep current FL quality thresholds and summaries unless tests show they must change.
- Tests proving current gDNA-compatible scenarios still train a usable gDNA FL model.

### Phase 4: Calibration Parity, Then Extension

Deliverables:

- Coarse-parity adapter from fine evidence to current density estimators.
- Replace direct old coarse payload reads in `density_global.py`, `_regional_exposure.py`,
  and `locus_prior.py`.
- Once parity passes, introduce fine-signature and boundary-transition density estimates.
- Add uncertainty fields to regional exposure summaries.

### Phase 5: Validation and Benchmarking

Deliverables:

- Unit tests for signature, index, scan payload, density parity, regional exposure, and
  locus prior.
- Golden output updates only after expected behavioral changes are understood.
- Synthetic benchmark comparison focused on:
  - pool-level gDNA/nRNA/mRNA error;
  - gDNA false RNA in large exonic regions;
  - high-density capture hotspots;
  - unstranded and weak-strand-identifiable libraries;
  - runtime and memory on large connected components.

## Open Questions

1. Should Phase 2 observe multimappers with `1/NH` fractional region evidence, as the old
   accumulator did, or keep current conservative exclusion? Recommendation: keep exclusion
   for this migration and revisit later.
2. Should terminal exon/intergenic boundaries contribute to boundary flux, or should the
   first calibration pass restrict boundary flux to exon-intron transitions? Recommendation:
   collect all boundary flux, then filter by adjacent fine signatures in calibration.
3. Should `signature` be the only persisted state, with booleans derived on load, or should
   booleans be persisted too? Recommendation: persist both during this major refactor for
   human auditability and robust tests.
4. Should spliced fragments contribute to boundary flux? The requested payload includes
   spliced left/right channels, so collect them. Calibration can choose whether to use them.
5. How should uncertainty enter the EM prior? Recommendation: first expose per-region
   evidence opportunity and posterior precision in calibration summaries, then downweight
   `gdna_prior_count_em` rather than changing likelihoods immediately.

## First PR Checklist

- [ ] Add signature helper functions and tests for all 16 states.
- [ ] Replace coarse region emission with fine region emission.
- [ ] Bump index format and fail hard on old indexes.
- [ ] Update `RegionArrays` to carry fine signatures and derived coarse fields.
- [ ] Update `_wire_calibration_regions(...)` to pass signatures to native code.
- [ ] Replace native calibration payload arrays with fractional fine evidence.
- [ ] Rebuild native extensions.
- [ ] Add scan accumulator tests for contained and boundary examples.
- [ ] Add coarse-parity adapter tests before changing calibration math.
