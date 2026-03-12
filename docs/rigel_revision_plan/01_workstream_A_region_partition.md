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

## 2. Implementation Goal

Produce a calibration-region table with constant annotation context per region.

The first-pass table should support:

- interval identity: `region_id`, `ref`, `start`, `end`, `length`
- a four-flag annotation summary
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

### 3.2 Normalize interval context into a four-flag schema

The current interval representation is useful but not yet the right public shape
for calibration. The first revision should normalize each region into fields
such as:

- `is_genic`
- `has_exon_pos`
- `has_exon_neg`
- `has_intron_pos`
- `has_intron_neg`
- `is_intergenic`
- `is_context_ambiguous`

This four-flag scheme is simpler and more composable than a large context enum,
and it handles opposite-strand overlaps directly.

It should be built as a calibration-specific refinement of the existing index
tiling rather than a full refactor of the index and resolver stack.

### 3.3 Keep tiny-region merging optional

Do not merge adjacent intervals in the first implementation unless diagnostics
show that the raw partition is too fragmented for stable evidence estimation.

The simplest initial rule is:

- preserve exact tiling from the index-derived partition
- postpone smoothing or merging to a later revision

## 4. Step-by-Step Coding Plan

1. identify the current loaded interval table shape in `TranscriptIndex.load()`
2. define a calibration-region schema in code and docs
3. implement a deterministic conversion from index intervals to calibration
   regions
4. add tests for interval coverage, disjointness, and reference completeness
5. add tests for representative annotation contexts:
   intergenic, exon-only, unambiguous intron, mixed or ambiguous
6. expose a read path so downstream calibration code can access the region table
   without re-deriving it repeatedly

## 5. Concrete Code Touchpoints

Most likely touchpoints:

- `src/rigel/index.py`
- possibly a new module: `src/rigel/calibration_regions.py`
- tests near index integrity and annotation behavior

Likely test additions:

- `tests/test_index.py`
- `tests/test_index_integrity.py`
- a new dedicated file such as `tests/test_calibration_regions.py`

## 6. Exit Criteria

Workstream A is complete when:

1. the genome can be loaded as a deterministic calibration-region table
2. every region has a stable context label and ambiguity flag
3. the table can be consumed downstream without revisiting annotation parsing