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

The redesign should avoid replacing that machinery. Instead it should expose a
calibration-oriented region layer derived from it.

One important correction is needed up front: `intervals.feather` is not yet a
true non-overlapping partition. It is a layered interval table containing EXON,
TRANSCRIPT, INTERGENIC, and UNAMBIG_INTRON rows that may share boundaries and
overlap each other. Workstream A therefore needs a deterministic remapping from
that layered representation to atomic calibration regions.

## 2. Implementation Goal

Produce a calibration-region table with constant annotation context per region.

More concretely, the output should be a true atomic partition of the genome:

- non-overlapping half-open intervals per reference
- one `region_id` per atomic interval
- a cgranges lookup keyed by `region_id`
- region metadata sufficient for later purity modeling and fragment overlap
  summarization

The first-pass table should support:

- interval identity: `region_id`, `ref`, `start`, `end`, `length`
- an implementation-stable annotation schema for exon and transcript-span
  presence by strand
- context ambiguity flags
- effective lengths for region-level density calculations

## 3. Recommended First-Pass Design

### 3.1 Add a calibration-region builder on top of the index

Recommended shape:

- add a pure helper in `src/rigel/index.py` or a new nearby module such as
  `src/rigel/calibration_regions.py`
- derive calibration regions from existing interval outputs rather than from
  raw GTF reprocessing
- keep the output as a pandas DataFrame initially for ease of inspection and
  testing

### 3.2 Build a true atomic partition from boundary sweeps

The existing interval table is the right substrate, but not the right final
object.

Recommended first implementation:

1. restrict the calibration partition to EXON, TRANSCRIPT, and INTERGENIC
   layers from `intervals.feather`
2. collect every unique start and end coordinate per reference from those rows
3. sort the boundaries and form adjacent half-open atomic segments
4. evaluate annotation membership on each atomic segment
5. discard zero-length segments

This gives a deterministic non-overlapping partition without requiring a second
annotation parser.

`UNAMBIG_INTRON` should remain available as a validation or auxiliary signal,
but it should not be the primitive geometry for the calibration partition.

### 3.3 Use four primitive booleans and six reported flags

The current interval representation is useful but not yet the right public shape
for calibration. The first revision should normalize each region into fields
such as:

- primitive source-of-truth booleans:
  - `exon_pos`
  - `exon_neg`
  - `tx_pos`
  - `tx_neg`
- derived reported flags:
  - `exon_ambig = exon_pos and exon_neg`
  - `tx_ambig = tx_pos and tx_neg`
- derived region properties:
  - `is_genic = exon_pos or exon_neg or tx_pos or tx_neg`
  - `is_intergenic = not is_genic`
  - `is_intronic_only = (not exon_pos) and (not exon_neg) and (tx_pos or tx_neg)`
  - `has_pos = exon_pos or tx_pos`
  - `has_neg = exon_neg or tx_neg`
  - `is_context_ambiguous = exon_ambig or tx_ambig`

This is the main design refinement relative to the current docs.

The implementation should treat the four primitive booleans as authoritative
state and materialize the two `*_ambig` flags as derived columns for reporting.
That avoids duplicated state and prevents impossible combinations from being
stored.

Use `tx_*` rather than `t_*` in code and docs. `t_*` is too easily confused
with transcript indices like `t_index` and with transcript-ID mappings.

It should be built as a calibration-specific refinement of the existing index
tiling rather than a full refactor of the index and resolver stack.

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

- calibration regions have their own `region_id` and cgranges index
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
- `exon_ambig_bp`
- `tx_pos_bp`
- `tx_neg_bp`
- `tx_ambig_bp`

Important caveat:

- exon and transcript-span base-pair totals are not additive with each other
- a base that lies in an exon of a transcript contributes to both exon and
  transcript-span totals

That is acceptable as long as these are documented as layer-specific overlap
totals rather than disjoint coverage bins.

## 4. Step-by-Step Coding Plan

1. identify the current loaded interval table shape in `TranscriptIndex.load()`
2. define the calibration-region schema with four primitive booleans and two
  derived ambiguity flags
3. implement a deterministic boundary-sweep conversion from index intervals to
  atomic calibration regions
4. build a cgranges index keyed by `region_id`
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
- possibly a small loader/index helper for `region_id` cgranges construction
- tests near index integrity and annotation behavior

Likely test additions:

- `tests/test_index.py`
- `tests/test_index_integrity.py`
- a new dedicated file such as `tests/test_calibration_regions.py`

## 6. Exit Criteria

Workstream A is complete when:

1. the genome can be loaded as a deterministic atomic calibration-region table
2. every region has stable primitive flags and invariant-derived ambiguity flags
3. the region table is non-overlapping and reference-complete
4. a `region_id` cgranges lookup exists for downstream fragment overlap
5. the table can be consumed downstream without revisiting annotation parsing