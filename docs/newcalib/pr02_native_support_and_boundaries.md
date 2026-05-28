# PR 02 - Native Observation Support and Boundary Payload

## Goal

Extend the native fractional accumulator so calibration can distinguish fractional mass from
physical fragment support. This is required for exposure shrinkage.

This PR can land before the two-state calibration rewrite if it is easier to validate independently.
Because it touches C++ under `src/rigel/native/`, rebuild with editable install before testing.

## Current State

`src/rigel/native/calibration/accumulator.h` exposes `CalibrationPayload` with:

- `region_counts[R, 12]` fractional mass
- `channel_mass[12]`
- `signature_mass[16]`
- `fl_pool_mass[6, 1024]`
- `fl_pool_total[6]`
- global counters such as `n_observed` and `n_unannotated_ref`

`src/rigel/native/calibration/accumulator.cpp` adds fractional mass inside `observe()` for every
overlapping `(fragment block, region)` pair. It does not currently expose a per-region count of
unique physical fragments touching the region.

## New Payload Fields

Add per-region support arrays to `CalibrationPayload`:

```text
region_unspliced_support[R]        uint32 or uint64
region_total_support[R]            uint32 or uint64, diagnostic optional
region_unspliced_mass_sq[R]        float64, optional but strongly recommended
```

`region_unspliced_support[r]` is the core field. It increments once for an accepted unspliced
physical fragment that contributes positive overlap mass to region `r`.

`region_unspliced_mass_sq[r]` stores the sum over fragments of squared total fractional mass into
that region. This makes it possible to compute a fractional effective sample size:

```text
ESS_r = M_r^2 / sum_i w_ir^2
```

This is better than a raw count when a region has many tiny overlap contributions. If implementation
risk is too high, add `region_unspliced_support` first and add `region_unspliced_mass_sq` in the same
PR only after the count invariant is green.

## Counting Rules

Count support at the physical fragment level:

- A fragment is eligible only after the existing gates: unique mapper, non-chimeric, non-artifact,
  routable strand, and positive aligned length.
- Count unspliced fragments for exposure support.
- Count each fragment once per region, even if it overlaps multiple exon blocks in that region.
- Count each fragment once per region, even if the fragment contributes to both left and right
  boundary channels of that region.
- Do not count unannotated-reference observations for any region.
- Multimappers remain excluded from calibration support, matching current fractional mass behavior.

The native loop can still add fractional channel mass as it does today. The support count needs a
second per-fragment reduction step that sees the unique set of region IDs touched by the fragment.

## Implementation Sketch

In `CalibrationAccumulator`:

- Add scratch vectors such as `touched_regions_`, `touched_region_mass_`.
- During `observe()`, whenever `overlap_bp > 0`, append `(rid, w)` to scratch in addition to adding
  channel mass.
- After all blocks are processed, sort the scratch pairs by `rid`, reduce duplicate region IDs, and
  update support arrays once per region.
- For `region_unspliced_mass_sq`, square the reduced per-fragment total mass for the region, not the
  individual channel pieces.

The fully-spanned case currently splits mass between boundary-left and boundary-right channels. The
support count must still increment once for that region.

## Python Wiring

Update these wrappers and validators:

- `src/rigel/native/bam_scanner.cpp`: include the new arrays in the `calibration` dict.
- `src/rigel/calibration/scan_payload.py`: validate shape, dtype, non-negativity, and expose the
  fields on `CalibrationScanPayload`.
- `src/rigel/calibration/_arrays.py`: sort the support arrays with `RegionArrays.order`.
- `src/rigel/calibration/region_count_ledger.py`: expose support views or a companion support table.
- `src/rigel/calibration/_diagnostics.py`: include support summaries and invariants.

## Native Boundary Payload

The current Python `BoundaryTable` in `src/rigel/calibration/boundaries.py` is transitional. It
derives boundary objects from region-side slots after scanning.

If this PR stays small, keep `BoundaryTable` for now and only add support arrays. If native work is
already stable, extend the payload with a native boundary table in the same general area:

```text
boundary_pos[B]
boundary_ref_id[B]
boundary_is_terminal[B]
left_region_idx[B]
right_region_idx[B]
boundary_counts[B, channel-like fields]
boundary_unspliced_support[B]
boundary_leff_left[B]
boundary_leff_right[B]
```

The migration path should be:

1. Add native boundary payload and Python wrapper.
2. Validate it exactly against the existing Python-derived `BoundaryTable`.
3. Switch calibration to consume the native payload.
4. Delete `build_boundary_table` production usage.

Do not do steps 3 and 4 until the exact reconstruction test is green.

## Tests

- Native accumulator unit test for one fragment spanning both boundaries of one region: mass splits
  across two boundary channels, support increments once.
- Native accumulator unit test for a spliced multi-block fragment hitting the same region twice:
  total support increments once if `region_total_support` is implemented; unspliced support does not
  increment.
- Payload validation test for shape and dtype errors.
- Python sorting test to ensure support arrays follow `RegionArrays.order`.
- If native boundary payload lands, exact equality test against the current `BoundaryTable` builder.

## Build and Test

After C++ edits:

```bash
conda activate rigel && pip install --no-build-isolation -e .
conda activate rigel && pytest tests/test_pipeline_wiring.py tests/test_region_count_ledger.py -v
```

Add the actual native accumulator test file to this command when the PR creates it.

## Done Means

- The calibration payload exposes physical support counts.
- Fractional mass identities remain unchanged.
- Existing scan balance assertions still pass.
- Support counts are not inflated by boundary-left/right splitting.
