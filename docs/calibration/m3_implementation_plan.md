# M3 — Calibration Accumulator: Refactor & Cleanup Plan

**Status:** M3 was implemented and committed in the milestone-3 work (the
accumulator, `RegionIndex`, scanner integration, and `CalibrationScanPayload`
all exist and pass 24 contract tests).  This plan is **not** an
implementation plan; it is a **simplification plan** to remove
over-engineered scaffolding that survived from earlier design iterations
and is now dead complexity.  After this refactor lands, M4–M9 can build on
a clean accumulator surface with no legacy.

**Rule for this milestone:** no behavioural change.  The
`CalibrationScanPayload` produced by the post-refactor scanner must be
**byte-identical** to the pre-refactor payload on every existing test
fixture (verified by `tests/test_calibration_accumulator.py` and the full
`pytest tests/` run).  This is a pure code-quality pass.

---

## 1. Motivation — what's wrong today

### 1.1 The inline + spill + dedup machinery in `observe()` is dead complexity

[src/rigel/native/calibration/accumulator.cpp](src/rigel/native/calibration/accumulator.cpp) contains ~70 lines of three-way dedup
(new-inline vs committed-prefix vs pre-existing-spill, then new-spill vs
inline-set vs older-spill).  This was written when `RegionIndex` could
return overlapping regions per call.  M2 locked the invariant that
**regions on a single reference are non-overlapping**, so each
`overlap_into` call now returns **sorted, distinct** region IDs by
construction.  The only remaining source of duplicates is multi-block
(spliced) fragments where two exon blocks both straddle the same
intronic gap and therefore both touch the boundary exon region.

That residual case can be handled in **one of two ways**, both of which
are an order of magnitude simpler than the current code:

* **Option A — Sorted merge:** maintain a single sorted `std::vector<int32_t>
  hits_` per fragment.  Each `overlap_into` call returns a sorted run; merge
  it into `hits_` with a two-pointer pass that skips duplicates.  ~15
  LoC, allocation-free after warm-up (`hits_.clear()` retains capacity).
* **Option B — Visited bitmap:** maintain a `std::vector<int32_t> epoch_`
  (size `n_regions`) and a fragment-scoped `current_epoch_`.  A region is
  "already added" iff `epoch_[rid] == current_epoch_`.  `++current_epoch_`
  per fragment.  O(1) dedup, but costs `4 * n_regions` bytes of resident
  memory per worker (~4–40 MB for human references).

**Decision:** Option A.  Per-fragment hit counts are ~1–10; sorted merge
of two ≤10-element arrays is faster than a bitmap touch and uses zero
additional memory.  The bitmap approach is only better at
`hits_per_fragment ≫ log(n_regions)` which never happens.

### 1.2 The inline-array vs spill-vector split is premature optimization

The current accumulator carries `std::array<int32_t, 16> inline_buf_` and
`std::vector<int32_t> spill_buf_`, with a custom split protocol in
`overlap_into<N>`.  The intent was to avoid heap traffic for the common
case (≤16 hits).  But:

* `std::vector<int32_t>::clear()` does not free memory; after warm-up
  the vector's backing storage is reused exactly like the inline array.
* The cost of the inline/spill split is paid in code complexity at
  every read site — the dedup loop in `observe()`, the bump-region
  fan-out, and the boundary-flux fan-out all have to fork on
  `i < n_inline ? inline_buf_[i] : spill_buf_[...]`.

**Decision:** delete the inline buffer.  Use a single
`std::vector<int32_t> hits_` reserved to capacity 16 in the constructor.
The first ~10 fragments grow it to its working size; from then on
`clear()` + `push_back()` is allocation-free.

### 1.3 `region_index != nullptr` guard pollutes the hot path

[src/rigel/native/bam_scanner.cpp](src/rigel/native/bam_scanner.cpp) currently treats the calibration path as
optional: `if (region_index != nullptr) ws.cal_acc.note_multimap();` etc.
Once the M8b switchover lands, calibration is mandatory and `set_regions`
will always have been called.  Until then we keep the guard, but we can
**collapse the four guarded sites into one method** on the scanner —
`maybe_observe(qname_group_state&)` — so the guard appears exactly once.

**Decision:** introduce `BamScanner::maybe_observe(...)` containing all
four observation-policy decisions (multimap, chimera, artifact, observe,
unobserved).  The four call sites in `process_qname_group_threaded`
become one.

### 1.4 `n_oor` is poorly named and partly dead

`n_oor` increments when `mask == 0` after the overlap pass — i.e., a
fragment whose every exon block landed entirely outside any region.
Under the M2 tiling invariant this is **impossible** for a fragment
whose `ref_id` is one of the references we partitioned.  It can only
happen on an unannotated decoy / contig.

The name "out of region" is opaque.  Rename to `n_unannotated_ref` and
add a one-line comment explaining the only mechanism that produces it.
The counter is still useful as a QC signal (large value ⇒ user pointed
the BAM at a different reference assembly than the index).

### 1.5 `kExonBit` is buried in an anonymous namespace

The mask bits are documented in three different places:
- `accumulator.h` block comment ("bit 0 = EXON, bit 1 = INTRON, ...")
- `accumulator.cpp` `constexpr uint8_t kExonBit = 0b001;` (anon ns)
- `regions.py` `RegionType` IntEnum + the comment above
  `type_mask = 1 << (2 - type)`

**Decision:** lift the mask constants to public namespace
`rigel::calibration::mask` in `region_index.h`:

```cpp
namespace rigel::calibration::mask {
constexpr uint8_t EXON       = 0b001;
constexpr uint8_t INTRON     = 0b010;
constexpr uint8_t INTERGENIC = 0b100;
constexpr int     N_STATES   = 8;
}  // namespace rigel::calibration::mask
```

Use them at every read site.  Python side already has `RegionType`;
no change needed there.

### 1.6 Worker merge is hand-rolled five times

`CalibrationAccumulator::merge_from` element-wise-adds:
- `global_counts` (8-element array)
- `per_region_counts` (R*8 vector)
- `fl_hist` (8*1024 vector)
- `u_left` (R vector)
- `u_right` (R vector)
- six scalar counters

Five separate `for (size_t i = 0; ...)` loops differ only by which
buffer.  Replace with a single `add_into(dst, src)` helper template, or
use `std::transform` with `std::plus<>`.  Cosmetic but reduces the
chance of asymmetric merges in future edits.

### 1.7 `CalibrationScanPayload.from_scan_dict` re-indirects through `d[...]`

[src/rigel/calibration/scan_payload.py](src/rigel/calibration/scan_payload.py) accepts a `Mapping[str, Any]` and pulls
each field with `d["..."]` plus shape/dtype/balance asserts.  This
indirection is fine for the validator pattern, but the caller side
(`pipeline.py`) does:

```python
calibration_payload = (
    CalibrationScanPayload.from_scan_dict(scan_result.get("calibration"),
                                          n_total=stats.n_read_names)
    if scan_result.get("calibration") is not None
    else None
)
```

The `None` branch exists because the scanner returns `None` when
`set_regions` was never called.  After M3 cleanup makes calibration
mandatory in the scan path (single point of guarding in `bam_scanner.cpp`,
not the Python caller), the `None` branch can be deleted; the validator
will then refuse `None` with the existing `_ERR_SUFFIX` message.

But mandatory calibration is an M8b decision; for M3 cleanup we leave
the `None` branch and just simplify the call site.

---

## 2. Concrete edits

Six commits, each focused and independently revertable.

### 2.1 Commit `m3-cleanup-1`: lift mask constants

* Add `rigel::calibration::mask::{EXON,INTRON,INTERGENIC,N_STATES}` to
  [src/rigel/native/calibration/region_index.h](src/rigel/native/calibration/region_index.h).
* Replace `kExonBit`, raw `0b001` / `0b010` / `0b100` literals, and
  `kMaskCard` references in `accumulator.{h,cpp}` and `bam_scanner.cpp`
  with the named constants.
* Delete the anonymous-namespace `constexpr uint8_t kExonBit` in
  `accumulator.cpp`.

### 2.2 Commit `m3-cleanup-2`: collapse inline + spill into a single sorted hit vector

* In `accumulator.h`:
  * Delete `std::array<int32_t, kInlineCap> inline_buf_`,
    `std::vector<int32_t> spill_buf_`, and `kInlineCap`.
  * Add `std::vector<int32_t> hits_` (constructor reserves 16).
* In `region_index.h`:
  * Replace the templated `overlap_into<N>(... out_inline, out_spill,
    already_inline)` with a single non-template
    `overlap_into(int32_t ref_id, int64_t qstart, int64_t qend,
    std::vector<int32_t>& out)` that **appends** sorted region IDs.
  * Caller is responsible for calling `clear()` between fragments; this
    matches the `hits_` semantics.
* In `accumulator.cpp`:
  * Replace the dedup body of `observe()` with a sorted-merge over
    per-block hit lists.  Concretely: keep one `scratch_block_` vector
    that `overlap_into` fills for the current exon block, and merge it
    into the persistent fragment-scope `hits_` vector with a two-pointer
    pass that skips duplicates.  The mask is OR-ed exactly when a region
    is **first added** to `hits_`.
  * Per-region count fan-out and boundary-flux fan-out simplify to a
    single `for (int32_t rid : hits_)` loop each.
* Net diff: ~70 LoC deleted, ~25 LoC added.

### 2.3 Commit `m3-cleanup-3`: rename `n_oor` → `n_unannotated_ref`

* `accumulator.h`: rename the field on `CalibrationPayload`.
* `accumulator.cpp`: rename the increment.
* `bam_scanner.cpp`: rename the `nb::dict` key.
* `scan_payload.py`: rename the field on `CalibrationScanPayload`,
  validator, and tests.
* Add a one-line comment at the increment site explaining the mechanism
  (decoy / contig / index-mismatch).
* Update [tests/test_calibration_accumulator.py](tests/test_calibration_accumulator.py) field reference and
  the `_good_payload_dict` fixture key.

### 2.4 Commit `m3-cleanup-4`: collapse worker-merge into one helper

* In `accumulator.cpp`, replace the five buffer-add loops in
  `merge_from` with a single helper:

  ```cpp
  template <typename T>
  static void add_into(std::vector<T>& dst, const std::vector<T>& src) {
      for (size_t i = 0; i < dst.size(); ++i) dst[i] += src[i];
  }
  ```

  (`std::array` overload likewise.)
* Apply.  Six scalar-counter additions stay inline.

### 2.5 Commit `m3-cleanup-5`: introduce `BamScanner::maybe_observe`

* Extract the four observation-policy sites from
  `process_qname_group_threaded` into a single private method
  `BamScanner::maybe_observe(WorkerState& ws, /* group state */)`.
* Move the `if (region_index != nullptr)` guard into that method so
  every other site is unconditional.
* Verify `tests/test_calibration_accumulator.py::TestObservationPolicy`
  still pins multimap / chimera; verify
  `TestPipelinePayload::test_run_pipeline_attaches_payload` still passes.

### 2.6 Commit `m3-cleanup-6`: documentation pass

* Drop the references to `m3_implementation_plan.md` from header
  comments now that this document supersedes it.
* Add an "M3 invariants" block comment to `accumulator.h` enumerating
  the four invariants the refactor relies on:
  1. Regions are per-ref non-overlapping (M2).
  2. Per-fragment hit counts are O(10).
  3. The scanner gates non-observable fragments before calling
     `observe()` (handled by `maybe_observe`).
  4. `payload.fl_hist` records `frag_end - frag_start` (genomic span,
     not FL); see plan §M7 for the consumer-side semantics.

---

## 3. Test impact

* No new tests required.  The 24 existing `test_calibration_accumulator.py`
  tests pin the contract.
* `TestWorkerMergeEquality::test_one_vs_four_workers_identical` is the
  master regression guard for commits 2 and 4.
* `TestMaskCorrectness` (5 cases) and `TestBoundaryFlux` (4 cases) are
  the master regression guards for commit 2.
* `TestObservationPolicy` is the master regression guard for commit 5.
* `TestPayloadValidation` covers the field renames in commit 3.

If any of these go red after a commit, **fix the implementation, not the
test** (per `.github/copilot-instructions.md`).

---

## 4. Exit gate

* All 6 commits land sequentially with each commit's `pytest tests/`
  green.
* `tests/test_calibration_accumulator.py`: 24 pass + 1 skip (unchanged).
* Net LoC change: ~−80 (the codebase shrinks).
* `accumulator.cpp` body of `observe()` is ≤ 50 LoC (currently ~95).
* The `inline + spill + already_inline` vocabulary is gone from the
  codebase (`grep` returns zero hits).

---

## 5. What this milestone deliberately does NOT do

* **Does not** widen `fl_hist` to `(8, 2, kFlBins)` (spliced/unspliced
  strata) — tracked as an M7 follow-up; out of scope here.
* **Does not** make calibration mandatory in the scanner — that is M8b.
* **Does not** delete the SRD-v1 calibration code paths — that is M8c.
* **Does not** change the on-disk schema, the nanobind ABI, or the
  Python public surface.  This is purely an internal cleanup.
