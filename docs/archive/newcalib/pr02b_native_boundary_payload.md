# PR 02B - Native Boundary Payload (deferred)

## Status

Split out of the original [pr02_native_support_and_boundaries.md](pr02_native_support_and_boundaries.md).
This is a separate goal from the per-region support counts in
[pr02a_v2.md](pr02a_v2.md) and lands after PR 02A. It can be sequenced
flexibly relative to PR 03/04 because nothing in the EB exposure model
strictly requires native boundary tables.

## Goal

Replace the Python-side transitional `BoundaryTable` (built in
`src/rigel/calibration/boundaries.py` by walking sorted region arrays) with
a native boundary payload emitted by the BAM scanner. The native payload
becomes the single source of truth for boundary geometry and per-boundary
support; the Python builder is deleted once the cutover is verified.

Motivation:

- The Python builder is reconstructing information the native scan already
  has on hand.
- Per-boundary fragment support is needed for boundary-local exposure /
  deconvolution diagnostics and is awkward to compute Python-side from
  fractional region channels.
- Reducing the boundary code path to a single native producer simplifies
  future calibration changes.

## Payload Shape

Add a `boundary` sub-dict to the scanner's `calibration` payload:

```text
boundary_pos[B]              int64    1-based or 0-based per the convention used
                                      by RegionIndex (mirror RegionIndex exactly)
boundary_ref_id[B]           int32
boundary_is_terminal[B]      uint8    1 if this boundary is the outer edge of
                                      a reference contig's region span (no
                                      neighbor on one side), else 0
left_region_idx[B]           int32    -1 if no region on the left side
right_region_idx[B]          int32    -1 if no region on the right side
boundary_unspliced_support[B] uint64  fragments that crossed this boundary
                                      with splice_type == UNSPLICED
boundary_spliced_support[B]   uint64  fragments that crossed this boundary
                                      with splice_type == SPLICED
boundary_leff_left[B]        float64  effective length of the left region as
                                      seen by the FL distribution at this
                                      boundary (matches existing semantics)
boundary_leff_right[B]       float64  symmetric
```

`B` is enumerated once at `set_regions()` time by walking the sorted region
index. Terminal boundaries (start of the first region or end of the last
region per contig) get `left_region_idx == -1` or `right_region_idx == -1`
respectively.

A fragment "crosses" boundary `b` iff its mapped extent has
`frag_start < boundary_pos[b] < frag_end`. This is computed once per
fragment per boundary it crosses; same once-per-(fragment, boundary)
invariant as PR 02A.

## Native Implementation Sketch

1. **Boundary enumeration.** In `RegionIndex` (or a new `BoundaryIndex`),
   after regions are loaded, do a linear pass per contig and emit the
   union of region starts and ends. Adjacent regions sharing a coordinate
   produce one boundary with both `left_region_idx` and `right_region_idx`
   set. Store boundaries sorted by `(ref_id, pos)` for fast lookup.

2. **Per-fragment crossing detection.** In `CalibrationAccumulator::observe`,
   after the existing per-block region-overlap loop, perform a single range
   query on the boundary index for boundaries in `(frag_start, frag_end)`.
   For each returned boundary, increment the appropriate support counter.

3. **Effective lengths.** Computed once at index build time; constant
   per boundary.

4. **Merging.** Element-wise `add_into` over the support vectors only;
   geometry vectors are identical across workers.

## Python Wiring

- `src/rigel/calibration/scan_payload.py`: add a `BoundaryScanPayload` nested
  dataclass with the same validation pattern as the region arrays.
- `src/rigel/calibration/boundaries.py`: keep the Python builder for one
  release as an oracle, behind a feature flag, used only for the exact
  reconstruction test below.
- `src/rigel/calibration/_arrays.py`: no reordering needed; boundary indices
  reference the native (unsorted) region order. The Python sort permutation
  is applied only when consumers ask for boundary-keyed-by-sorted-region.

## Cutover Plan

1. Land the native payload alongside the existing Python `BoundaryTable`.
2. Add an exact equality test: for every benchmark BAM, the native payload
   reproduces the geometry of the Python `BoundaryTable` row-for-row, and
   the per-boundary support counts derived from a Python reference
   computation agree with the native counts.
3. Switch every production consumer to the native payload.
4. Delete `build_boundary_table` and its tests.

Do **not** do steps 3-4 until the exact reconstruction test is green on
all calibration scenarios.

## Tests

- Geometry parity: native payload equals `build_boundary_table` output on
  the standard synthetic fixtures.
- Crossing count: a hand-crafted fragment crossing exactly two boundaries
  contributes `+1` to both, regardless of how many region overlaps it has.
- Terminal flags: the first and last boundary per contig have one
  `*_region_idx == -1` and `boundary_is_terminal == 1`.
- Spliced vs unspliced separation: same as PR 02A's case 3, lifted to the
  boundary level.

## Out of Scope

- Any redesign of how exposure / deconvolution consume boundaries; this PR
  only changes who produces the boundary table.
- Changing the boundary effective-length formula; it must reproduce the
  current Python definition exactly.

## Done Means

- BAM scanner emits a `boundary` sub-dict with the schema above.
- `BoundaryScanPayload` validates it and is the new public source of
  boundary information.
- The Python `build_boundary_table` is removed from the production path
  (kept only behind a test-only flag if useful).
- All boundary-related calibration tests pass against the native payload.
