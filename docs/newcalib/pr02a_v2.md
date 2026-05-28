# PR 02A - Native Per-Region Physical Support Counts (v2)

## Status

Supersedes the "Native Observation Support" section of
[pr02_native_support_and_boundaries.md](pr02_native_support_and_boundaries.md).
The boundary-payload portion of the original PR 02 is split off into
[pr02b_native_boundary_payload.md](pr02b_native_boundary_payload.md).

## Goal

Add two per-region integer support vectors to the native calibration payload:

```text
region_unspliced_support[R]   uint32
region_spliced_support[R]     uint32
```

Each entry is the number of distinct physical fragments that contributed
**any** overlap mass to the region, partitioned by the fragment's splice class.
Each fragment is counted at most once per region, regardless of how many of
its blocks touch the region or whether it is contained, boundary-left,
boundary-right, or fully spanning. These counts are the effective sample
sizes the EB exposure model will use for shrinkage (see PR 04). They are
**not** related to fractional mass; the existing fractional `region_counts`
contract is untouched.

## Critique of Original PR 02 (resolved here)

1. **No fractional ESS.** Original called for `region_unspliced_mass_sq[r]` to
   compute `M_r^2 / sum_i w_ir^2`. This is rejected. ESS for shrinkage is the
   integer count of fragments touching the region. A fragment that contributes
   a tiny overlap fraction is still one independent sample of the regional
   sampling rate; squared-mass weighting confuses sample identity with mass
   weight. Drop the field entirely.
2. **Missing spliced counts.** Original tracked only unspliced support. We need
   both classes because spliced support is a useful diagnostic and is required
   by downstream calibration audits (e.g. distinguishing intronic regions with
   real RNA support from those with only crossing spliced reads).
3. **Once-per-region invariant must be explicit.** Original spent words on the
   fully-spanned case adding to two channels. The support count is
   compartment-agnostic: one increment per `(fragment, region)` pair, full stop.
4. **Keep it simple.** No `region_total_support` synthetic field, no float
   accumulators, no ESS arithmetic. Two `uint32` vectors that we increment.

## What Changes in Native Code

### `src/rigel/native/calibration/accumulator.h`

Add to `CalibrationPayload`:

```cpp
std::vector<std::uint32_t> region_unspliced_support;  // size = n_regions
std::vector<std::uint32_t> region_spliced_support;    // size = n_regions
```

Constructor allocates both with size `n_regions`, zero-initialized. No new
counters; no new constants.

Add one small scratch vector to `CalibrationAccumulator`:

```cpp
std::vector<std::int32_t> touched_regions_;  // reused per observe()
```

Reserve a small warm capacity (~16) alongside `block_hits_`.

### `src/rigel/native/calibration/accumulator.cpp`

Inside `observe()`, before the per-block loop begins, clear
`touched_regions_`. Inside the inner loop where we already test
`overlap_bp > 0` and add fractional mass, also push `rid` onto
`touched_regions_`. After the per-block loops complete, reduce duplicates
and increment the right support vector exactly once per unique region:

```cpp
std::sort(touched_regions_.begin(), touched_regions_.end());
auto last = std::unique(touched_regions_.begin(), touched_regions_.end());
auto& support = (splice_idx == region_signature::kSpliceUnspliced)
                    ? payload_.region_unspliced_support
                    : payload_.region_spliced_support;
for (auto it = touched_regions_.begin(); it != last; ++it) {
    ++support[static_cast<size_t>(*it)];
}
```

Rules baked into the placement:

- Push happens only when `overlap_bp > 0`, so empty intersections never count.
- Push happens inside the same gate that adds fractional mass, so the invariant
  "support counted iff fractional mass added" holds by construction.
- The sort+unique reduction makes the per-region increment exactly once,
  independent of how many blocks of the fragment hit the region and independent
  of compartment (contained vs left vs right vs fully spanning).
- The splice class is the fragment's splice class, decided once per fragment
  (matches how `region_counts` channels are routed today). Mixed-block
  fragments cannot occur: `splice_type` is a per-fragment label.
- `n_unannotated_ref` already detects fragments with no region hit; if
  `touched_regions_` is empty after both loops, the sort/unique loop is a
  no-op. No extra branch needed.

### `merge_from`

Add two more `add_into(...)` calls for the new vectors. They are `uint32`;
overflow is theoretically possible only at extremely large depth. To keep
the merge robust, perform addition in `uint64` accumulator inside `add_into`
or change the storage type to `uint64`. Recommendation: store as `uint64_t`
to match the cost of one cache line per ~8 regions and avoid any future
overflow audit. The Python payload is converted to `np.uint64` at the
wrapper boundary.

Update the suggested types accordingly:

```cpp
std::vector<std::uint64_t> region_unspliced_support;
std::vector<std::uint64_t> region_spliced_support;
```

### `src/rigel/native/bam_scanner.cpp`

In the block that builds `cal_dict`, add (right after the `region_counts`
emission):

```cpp
cal_dict["region_unspliced_support"] =
    vec_to_ndarray(std::move(payload.region_unspliced_support));
cal_dict["region_spliced_support"] =
    vec_to_ndarray(std::move(payload.region_spliced_support));
```

Use the existing 1-D `vec_to_ndarray` helper (these are shape `(R,)`).

## Python Wiring

### `src/rigel/calibration/scan_payload.py`

Extend `CalibrationScanPayload`:

```python
region_unspliced_support: np.ndarray  # uint64[R]
region_spliced_support:   np.ndarray  # uint64[R]
```

In `from_scan_dict`, validate both with `_check_array(name, arr, np.uint64,
(n_regions,))`. The existing helper already enforces non-negativity (vacuous
for unsigned), C-contiguity, and shape.

Add three lightweight invariants alongside the existing mass identities:

1. Per-region support never exceeds `n_observed`:
   `region_unspliced_support.max(initial=0) <= n_observed`
   `region_spliced_support.max(initial=0)  <= n_observed`
2. Sum of unspliced support is bounded by `n_observed - n_unannotated_ref`
   times the maximum number of regions a single fragment can touch. We cannot
   tighten this without knowing fragment fan-out, so the only universal bound
   we check is the per-region cap above.
3. A region with positive unspliced fractional mass must have positive
   unspliced support:
   `support == 0` implies `unspliced fractional mass == 0`. This is the
   contract that downstream code will rely on, so it is worth a one-shot
   assertion in the validator. Use a relative tolerance equal to the float32
   noise floor of `region_counts`.

### `src/rigel/calibration/_arrays.py`

`RegionArrays.from_payload` already reorders fractional arrays by `order`.
Add the two support vectors to `PayloadArrays` (or whatever the existing
sorted view dataclass is called) and reorder identically:

```python
unspliced_support_sorted = payload.region_unspliced_support[order]
spliced_support_sorted   = payload.region_spliced_support[order]
```

Store both as `np.uint64` and assert C-contiguity after the gather.

### `src/rigel/calibration/region_count_ledger.py`

`RegionCountLedger` currently exposes fractional totals only. Add two
read-only properties that return the sorted support vectors:

```python
@property
def unspliced_support(self) -> np.ndarray:  # uint64[R]
    return self._payload_arrays.region_unspliced_support_sorted

@property
def spliced_support(self) -> np.ndarray:    # uint64[R]
    return self._payload_arrays.region_spliced_support_sorted
```

No new derived quantities. Downstream PRs (the EB exposure model) will
consume the property directly.

### `src/rigel/calibration/_diagnostics.py`

Add summary lines to the existing diagnostics block:

- total unspliced support, total spliced support
- count of regions with `unspliced_support == 0` vs total
- 50th / 90th / 99th percentile of nonzero unspliced support

These are observational only; nothing gates on them.

## Tests

New native-level unit test (Python through the BAM scanner is easiest; a
direct gtest is acceptable if one already exists for the accumulator):

1. **Once-per-region invariant, fully spanning.** One unspliced fragment whose
   single block fully spans one region. Expect:
   - fractional mass split 50/50 between `BOUNDARY_LEFT` and
     `BOUNDARY_RIGHT` channels (unchanged behavior).
   - `region_unspliced_support[r] == 1`, `region_spliced_support[r] == 0`.

2. **Once-per-region invariant, multi-block contained.** One spliced fragment
   with two exon blocks both falling inside the same region. Expect:
   - `region_spliced_support[r] == 1`, `region_unspliced_support[r] == 0`.

3. **Spliced and unspliced separation.** Two fragments touching the same
   region, one spliced and one unspliced. Expect both support vectors to
   read `1` at that region.

4. **Multi-region fanout.** One unspliced fragment overlapping three distinct
   regions. Expect `region_unspliced_support` to read `1` at each of the
   three, and `sum(region_unspliced_support) == 3` for the fragment.

5. **Strand-ambiguous fragment.** Already returns early without mass; assert
   support vectors remain zero.

6. **No region hits.** A fragment with `n_unannotated_ref += 1` adds zero
   support.

7. **Payload validator.** Construct a scan dict with an inconsistent shape
   for `region_unspliced_support` and assert `from_scan_dict` raises.

8. **Sort invariance.** Round-trip through `RegionArrays.from_payload` and
   verify the sorted support vector indexed by the inverse permutation
   equals the original payload vector.

Existing scan-balance tests (mass identities, `region_counts` sums) must
continue to pass unchanged.

## Build and Test

After C++ edits:

```bash
conda activate rigel && pip install --no-build-isolation -e .
conda activate rigel && pytest \
    tests/test_calibration_accumulator.py \
    tests/test_pipeline_wiring.py \
    tests/test_region_count_ledger.py \
    -v
```

If a dedicated accumulator test file does not exist yet, add
`tests/test_calibration_accumulator.py` housing the eight cases above.

## Out of Scope

- Native boundary payload (`boundary_pos`, `left_region_idx`, etc.). Tracked
  in [pr02b_native_boundary_payload.md](pr02b_native_boundary_payload.md).
- Any consumption of the new support vectors by calibration math. Tracked in
  PR 03 (regional gDNA mass contract) and PR 04 (EB exposure factor).
- Fractional ESS / squared-mass accumulation. Explicitly rejected; see the
  critique above.

## Done Means

- Native payload exposes `region_unspliced_support[R]` and
  `region_spliced_support[R]` as `uint64`.
- Both vectors are incremented exactly once per `(fragment, region)` pair,
  gated identically to fractional mass.
- Python wrapper validates shape and dtype, exposes sorted views, and asserts
  the "positive mass implies positive support" invariant.
- All existing calibration tests pass; the new eight-case test file is green.
- No changes to existing fractional mass identities, channels, signatures, or
  FL pools.
