# Workstream A: Region Partition

This is the cleanest place to begin coding because the repository already has
most of the underlying machinery.

## 1. What Already Exists

`src/rigel/index.py` already constructs a genome tiling and writes
`intervals.feather`.

Relevant existing concepts:

- exon intervals
- transcript-span intervals
- unambiguous intronic intervals
- intergenic intervals
- splice junction table
- transcript and nRNA tables

The redesign should avoid replacing that machinery. Instead it should add a
calibration-oriented flattened region layer alongside it.

One important correction is needed up front: `intervals.feather` is not yet a
true non-overlapping partition. It is a layered interval table containing EXON,
TRANSCRIPT, INTERGENIC, and UNAMBIG_INTRON rows that may share boundaries and
overlap each other. Workstream A therefore needs its own flattened region index
rather than a direct reuse of `intervals.feather`.

## 2. Implementation Goal

Produce a calibration-region table with constant annotation context per region.

More concretely, the output should be a true flattened non-overlapping region
index of the genome:

- non-overlapping half-open intervals per reference
- one `region_id` per atomic interval
- persisted metadata in `regions.feather` with optional `regions.tsv`
- a cgranges lookup keyed by `region_id`
- region metadata sufficient for later purity modeling and fragment overlap
  summarization

The first-pass table should support:

- interval identity: `region_id`, `ref`, `start`, `end`, `length`
- four stored annotation flags: `exon_pos`, `exon_neg`, `tx_pos`, `tx_neg`
- effective lengths for region-level density calculations

## 3. Recommended First-Pass Design

### 3.1 Add a calibration-region builder on top of the index

Recommended shape:

- add a pure helper in `src/rigel/index.py` or a new nearby module such as
  `src/rigel/calibration_regions.py`
- build flattened regions directly from transcript and exon coordinates already
  present during index build
- keep the output as a pandas DataFrame initially for ease of inspection and
  testing
- write `regions.feather` at index build time, with optional `regions.tsv`
- load `regions.feather` into a dedicated cgranges structure at index load time

### 3.2 Build flattened non-overlapping bins from boundary sweeps

The transcript and exon annotations are the right substrate, but not the right
final object.

Recommended first implementation:

1. collect all transcript-span start and end coordinates per reference
2. collect all exon start and end coordinates per reference
3. union and sort those coordinates into one boundary array per reference
4. form adjacent half-open bins between successive boundaries
5. discard zero-length bins
6. assign the four boolean flags to each bin by testing whether that bin is
  spanned by a pos/neg transcript span and by a pos/neg exon
7. assign `region_id` after the flattened bins are finalized

This gives a deterministic non-overlapping partition of the annotation geometry
itself. It is the flattened index.

`UNAMBIG_INTRON` should remain available as a validation or auxiliary signal,
but it should not be the primitive geometry for the calibration partition.

### 3.3 Use exactly four stored flags

The current interval representation is useful but not yet the right public shape
for calibration. The first revision should store exactly these four booleans:

- `exon_pos`
- `exon_neg`
- `tx_pos`
- `tx_neg`

Everything else is inferred:

- strand pos if `exon_pos or tx_pos`
- strand neg if `exon_neg or tx_neg`
- strand ambiguous if `(exon_pos or tx_pos) and (exon_neg or tx_neg)`
- strand none if all four flags are false
- intergenic if all four flags are false
- genic if any flag is true
- intronic if both exon flags are false and either transcript-span flag is true

This is the main design refinement relative to the current docs.

The implementation should store only these four booleans in `regions.feather`.
No `_ambig` columns should be written. Ambiguity is an interpretation of the
four stored flags, not persisted state.

Use `tx_*` rather than `t_*` in code and docs. `t_*` is too easily confused
with transcript indices like `t_index` and with transcript-ID mappings.

It should be built at index-build time from the same transcript objects already
used to produce the existing index artifacts.

### 3.4 Keep transcript resolution separate in the first pass

The region partition may eventually support sidecar mappings such as
`region_id -> transcript_id set`, but that should not be part of the first
implementation contract.

Reason:

- transcript resolution already has a dedicated index and merge path
- mixing calibration-region lookup with transcript resolution would duplicate
  hot-path logic prematurely
- the immediate goal is calibration and purity modeling, not replacing the main
  fragment resolver

So the first-pass contract should be:

- calibration regions have their own `region_id`, `regions.feather`, and
  cgranges index
- transcript resolution continues to use the existing `TranscriptIndex` path
- any region-to-transcript mapping is an optional sidecar for diagnostics or a
  later simplification pass

### 3.5 Keep tiny-region merging optional

Do not merge adjacent intervals in the first implementation unless diagnostics
show that the raw partition is too fragmented for stable evidence estimation.

The simplest initial rule is:

- preserve exact atomic tiling from the boundary-swept partition
- postpone smoothing or merging to a later revision

### 3.6 Define the fragment-facing overlap contract now

Fragments will often overlap multiple atomic regions. That is expected and
should not be treated as an immediate ambiguity.

The right first contract is a merged fragment-region summary with fields such as:

- `region_ids`
- `overlap_bp_total`
- `exon_pos_bp`
- `exon_neg_bp`
- `tx_pos_bp`
- `tx_neg_bp`

Important caveat:

- exon and transcript-span base-pair totals are not additive with each other
- a base that lies in an exon of a transcript contributes to both exon and
  transcript-span totals

That is acceptable as long as these are documented as layer-specific overlap
totals rather than disjoint coverage bins.

## 4. Step-by-Step Coding Plan

1. identify the current loaded interval table shape in `TranscriptIndex.load()`
2. define the `regions.feather` schema with the four stored flags
3. implement a deterministic flattened-index builder from transcript and exon
  boundaries at index build time
4. write `regions.feather` and optional `regions.tsv`
5. load `regions.feather` and build a cgranges index keyed by `region_id`
5. add tests for interval coverage, non-overlap, and reference completeness
6. add tests for representative annotation contexts:
  intergenic, exon-only, transcript-only, intronic-only, opposite-strand
  overlap, and strand-ambiguous overlap
7. expose a read path so downstream calibration code can access the region table
  without re-deriving it repeatedly
8. keep any optional `region_id -> transcript_id` sidecar explicitly separate
  from the core partition contract

## 5. Concrete Code Touchpoints

Most likely touchpoints:

- `src/rigel/index.py`
- possibly a new module: `src/rigel/calibration_regions.py`
- `regions.feather` / optional `regions.tsv` build and load paths
- possibly a small loader/index helper for `region_id` cgranges construction
- tests near index integrity and annotation behavior

Likely test additions:

- `tests/test_index.py`
- `tests/test_index_integrity.py`
- a new dedicated file such as `tests/test_calibration_regions.py`

## 6. Exit Criteria

Workstream A is complete when:

1. the genome can be loaded as a deterministic atomic calibration-region table
2. every region has stable `exon_pos`, `exon_neg`, `tx_pos`, and `tx_neg` flags
3. the region table is non-overlapping and reference-complete
4. `regions.feather` is written at index build time and loaded successfully
5. a `region_id` cgranges lookup exists for downstream fragment overlap
6. the table can be consumed downstream without revisiting annotation parsing