# M3 — C++ Calibration Scanner: Implementation Plan

**Status:** Implementation contract for milestone M3 of
`docs/calibration/calibration_v6_plan.md`. This plan supersedes the M3 sketch
in §4 of the parent plan with concrete file layouts, function signatures,
and a simpler-than-spec'd architecture that takes advantage of what M1
and M2 already give us.

**Predecessor work (lost):** `src/rigel/native/calibration/` from the
2026-05-01 incident contained `calibration_sink.h` + a per-worker
accumulator that interleaved scanning and observation logic. The lesson
from re-deriving it: most of that complexity came from defending against
inconsistent region partitions. M2 closed those failure modes via
`validate_against_ref_lengths`, so we can build a much smaller scanner.

---

## 0. Goals (in priority order)

1. **Per-fragment observation in O(log R + k)** where R is the region
   table size and k is the number of regions a fragment touches. Real
   libraries hit k ≤ 4 for ~99% of fragments.
2. **Allocation-free hot path.** No heap traffic in the per-fragment
   `observe()` after the first warm-up fragment per worker.
3. **Lossless worker merge.** N-worker output bit-identical to
   1-worker output. This is the only thing parallelism can silently
   corrupt; the merge-equality test is the gate.
4. **Single observation predicate** at exactly one site in
   `process_qname_group_threaded`. The spec text talks about "three
   observation sites" — that conflates *counting* with *observing*.
   In our design only unique-mapper, non-chimeric, non-artifact
   fragments are *observed*; the other classes are *counted only*.

---

## 1. Simplifications vs. the parent plan §4-M3

| Topic | Parent plan | This plan |
|---|---|---|
| Observation sites | "three sites: intergenic-resolved, chimeric, genic-resolved" | **one** site at the end of `process_qname_group_threaded`; chimeric/multimap/artifact/oor are counter-only branches |
| Region overlap engine | "binary-search keyed by ref_id; cgranges-like" | Per-ref CSR (`ref_starts[ref_id]` → `[lo, hi)` slice into a flat region array). Non-overlapping by M2 invariant ⇒ a single `lower_bound` on `ends` + linear scan to `start ≥ query.end` |
| `boundary_tol` | parameter on `observe()` | compile-time constant `kBoundaryTol = 0`; can be promoted to a knob later if real data needs it |
| `SmallRegionSet` | "inline-16 + spill vector with linear-scan dedup" | `std::array<int32_t,16>` + `std::vector<int32_t>` spill — both members of `WorkerState`, reused per fragment via `clear()` |
| `set_regions` API surface | `(ref_names, starts, ends, type_masks)` | Four typed numpy arrays `(ref_ids: int32, starts: int64, ends: int64, type_masks: uint8)`. The Python caller (M2 already gave us `ref_name_to_id`) does ref-name resolution; the C++ never sees strings |
| Region/payload schema in C++ | spread across `region_index.h` + `accumulator.h` | one `CalibrationPayload` POD struct in `accumulator.h` is the ONLY thing exported to Python; numpy ownership transferred via `vec_to_ndarray` move |
| `Pass D` boundary-flux loop | separate code path | folded into the same overlap walk; two side effects in one pass |

What we keep from the parent plan verbatim: the 8-state mask layout
(bit 0 = EXON, bit 1 = INTRON, bit 2 = INTERGENIC), the `fl_hist`
shape `(8, 1024)`, the `u_L`/`u_R` counters, the per-fragment `+1`
counting unit, and the four exclusion counters.

---

## 2. File layout (new + modified)

```
src/rigel/native/calibration/                       (new directory)
├── region_index.h        ~120 LoC   header-only
├── small_region_set.h    ~60  LoC   header-only
├── accumulator.h         ~120 LoC   header + inline observe() body
└── accumulator.cpp       ~150 LoC   merge_from + numpy export helpers

src/rigel/native/bam_scanner.cpp                    (modified)
  + include "calibration/accumulator.h"
  + WorkerState gains: CalibrationAccumulator cal_acc;
  + BamScanner gains: std::unique_ptr<RegionIndex> region_index_;
                       CalibrationAccumulator cal_acc_merged_;
  + nb::def("set_regions", ...)
  + one observation call inside process_qname_group_threaded
  + merge step inside scan() after workers.join()
  + build_result() emits result["calibration"] dict

src/rigel/calibration/scan_payload.py               (new)
  ~80 LoC frozen dataclass + from_scan_dict validator

CMakeLists.txt                                       (modified)
  + add the calibration/*.cpp source to _bam_impl target

tests/test_calibration_accumulator.py               (new)
  ~9 tests, ~250 LoC
```

Total new C++: **~450 LoC** (header-heavy by design, hot path is
inlinable). Total modified C++: ~80 LoC of additions to bam_scanner.cpp.

---

## 3. Detailed designs

### 3.1 `RegionIndex` (header-only, immutable after `set`)

```cpp
namespace rigel::calibration {

class RegionIndex {
public:
    // Single setup call. Arrays are copied into private buffers.
    // Inputs MUST be sorted by (ref_id, start) and per-ref
    // contiguous + non-overlapping. M2's validate_against_ref_lengths
    // guarantees this; we additionally `assert` it in DEBUG builds.
    void set(const int32_t* ref_ids,
             const int64_t* starts,
             const int64_t* ends,
             const uint8_t* type_masks,
             int64_t n_regions,
             int32_t n_refs);

    // O(log n + k) overlap. Appends region indices into `out_inline`,
    // spilling to `out_spill` only when inline is full.  Returns total
    // count.  Output is in start order.
    template <std::size_t N>
    int32_t overlap_into(int32_t ref_id,
                         int64_t start,
                         int64_t end,
                         std::array<int32_t, N>& out_inline,
                         std::vector<int32_t>& out_spill) const;

    // Trivial accessors (used by accumulator's boundary-flux pass)
    inline int64_t  start(int32_t rid)     const { return starts_[rid]; }
    inline int64_t  end(int32_t rid)       const { return ends_[rid]; }
    inline uint8_t  type_mask(int32_t rid) const { return type_masks_[rid]; }
    inline int64_t  n_regions()            const { return starts_.size(); }
    inline int32_t  n_refs()               const { return n_refs_; }

private:
    // Flat per-region columns (parallel arrays, region_id is the index).
    std::vector<int64_t> starts_;
    std::vector<int64_t> ends_;
    std::vector<uint8_t> type_masks_;

    // Per-ref CSR offsets: ref_offsets_[r], ref_offsets_[r+1] is
    // [lo, hi) into the flat arrays for ref_id r.
    std::vector<int64_t> ref_offsets_;
    int32_t n_refs_ = 0;
};

}  // namespace rigel::calibration
```

**Overlap implementation:**

```cpp
template <std::size_t N>
int32_t RegionIndex::overlap_into(
    int32_t ref_id, int64_t qstart, int64_t qend,
    std::array<int32_t, N>& out_inline,
    std::vector<int32_t>& out_spill) const
{
    if (ref_id < 0 || ref_id >= n_refs_) return 0;
    const int64_t lo = ref_offsets_[ref_id];
    const int64_t hi = ref_offsets_[ref_id + 1];
    if (lo == hi) return 0;

    // First region whose end > qstart (strictly).
    auto it = std::upper_bound(
        starts_.data() + lo, starts_.data() + hi, qstart,
        [](int64_t v, int64_t s){ return v < s; });
    // upper_bound on starts gives the first start > qstart; we want
    // the region BEFORE that if its end > qstart.
    int64_t i = (it - starts_.data());
    if (i > lo && ends_[i - 1] > qstart) i = i - 1;

    int32_t n = 0;
    while (i < hi && starts_[i] < qend) {
        if (n < static_cast<int32_t>(N))      out_inline[n] = static_cast<int32_t>(i);
        else                                   out_spill.push_back(static_cast<int32_t>(i));
        ++n; ++i;
    }
    return n;
}
```

The non-overlap invariant (M2) means we never need a second-pass dedup.
The query never returns the same region twice. This is what makes the
algorithm O(log n + k) instead of O(log n + R) for a general interval
tree.

### 3.2 `SmallRegionSet` (header-only)

This is just a thin wrapper around the inline + spill arrays held by
`RegionIndex::overlap_into`. We don't need a class — the wrapper would
add no value over a `std::array<int32_t,16> inline_; std::vector<int32_t> spill_;`
pair. Plan-text "SmallRegionSet" maps to: **two named member vars on
`CalibrationAccumulator`** that are `clear()`'d at the top of every
`observe()`. No template, no class, no header file needed.

(If a future caller needs to dedup ACROSS multiple `overlap_into` calls
per fragment — e.g., one per exon block — we'd want a real
`SmallRegionSet`. For now, the union is computed by calling
`overlap_into` once per exon block and feeding all results into the
SAME inline+spill buffers. Linear-scan dedup across the combined list
is O(k²) but k ≤ 4 in practice; we measure this in the perf test.)

### 3.3 `CalibrationAccumulator`

```cpp
namespace rigel::calibration {

struct CalibrationPayload {
    // Fixed-size global counters
    std::array<int64_t, 8> global_counts {};

    // Per-region counts: per_region_counts_[8 * region_id + mask]
    std::vector<int64_t> per_region_counts;   // size = 8 * n_regions

    // Per-mask FL histogram: fl_hist_[1024 * mask + length]
    std::vector<int64_t> fl_hist;             // size = 8 * 1024

    // Boundary-flux counters (EXON regions only, but indexed by region_id)
    std::vector<int64_t> u_left;              // size = n_regions
    std::vector<int64_t> u_right;             // size = n_regions

    // Exclusion / observation counters
    int64_t n_observed             = 0;
    int64_t n_excluded_multimap    = 0;
    int64_t n_excluded_chimera     = 0;
    int64_t n_excluded_artifact    = 0;
    int64_t n_oor                  = 0;
};

class CalibrationAccumulator {
public:
    explicit CalibrationAccumulator(int64_t n_regions);

    // The hot-path observation. exons is the fragment's R1∪R2 exon
    // blocks (sorted, merged) in genomic coords on a SINGLE ref.
    // Caller has already gated on:
    //   - unique mapper
    //   - non-chimeric (so all exons share ref_id)
    //   - non-artifact splice_type
    //
    // Increments global_counts[mask], per_region_counts[region][mask],
    // fl_hist[mask][min(fl, 1023)], u_left[]/u_right[] (only when
    // splice_type == SPLICE_UNSPLICED), and n_observed. Sets the
    // n_oor flag if mask == 0.
    void observe(int8_t splice_type,
                 int32_t ref_id,
                 int64_t frag_start,
                 int64_t frag_end,
                 const ExonBlock* exons,
                 int32_t n_exons,
                 const RegionIndex& regions);

    // Counter-only branches (called from process_qname_group_threaded
    // for fragments that should not enter observe()).
    inline void note_multimap()  { payload_.n_excluded_multimap++; }
    inline void note_chimera()   { payload_.n_excluded_chimera++; }
    inline void note_artifact()  { payload_.n_excluded_artifact++; }

    // Worker-merge: payload_-additive.
    void merge_from(CalibrationAccumulator& other);

    // Move out the payload (called once after merge).
    CalibrationPayload take() && { return std::move(payload_); }

private:
    CalibrationPayload payload_;
    // Reused per-fragment scratch.
    std::array<int32_t, 16> inline_buf_ {};
    std::vector<int32_t> spill_buf_;
};

}  // namespace rigel::calibration
```

**`observe()` algorithm (folded Pass A + Pass D):**

```cpp
void CalibrationAccumulator::observe(...) {
    // Pass A: walk every exon block, query overlap into shared
    // inline+spill buffers. Build mask along the way.
    spill_buf_.clear();
    int32_t n_inline = 0;
    uint8_t mask = 0;

    // Loop over exon blocks; merge their overlaps in-place.
    // Linear-scan dedup against what we have so far (k ≤ 4 in practice).
    for (int32_t e = 0; e < n_exons; ++e) {
        const ExonBlock& blk = exons[e];
        // Inline temp for this block
        std::array<int32_t, 16> blk_inline;
        std::vector<int32_t> blk_spill;
        int32_t n_blk = regions.overlap_into(
            blk.ref_id, blk.start, blk.end, blk_inline, blk_spill);
        for (int32_t i = 0; i < n_blk; ++i) {
            int32_t rid = (i < 16) ? blk_inline[i] : blk_spill[i - 16];
            // Linear-scan dedup against current set
            bool dup = false;
            for (int32_t j = 0; j < n_inline; ++j)
                if (inline_buf_[j] == rid) { dup = true; break; }
            if (!dup) {
                for (int32_t rs : spill_buf_)
                    if (rs == rid) { dup = true; break; }
            }
            if (dup) continue;
            // Append + update mask via type_mask
            mask |= (1u << (2 - regions.type_mask(rid)));  // see §2.3
            if (n_inline < 16) inline_buf_[n_inline++] = rid;
            else                spill_buf_.push_back(rid);
        }
    }

    // Always: bump global mask counter, FL histogram.
    payload_.global_counts[mask]++;
    int64_t fl = frag_end - frag_start;
    int64_t fl_idx = std::min<int64_t>(std::max<int64_t>(fl, 0), 1023);
    payload_.fl_hist[mask * 1024 + fl_idx]++;
    payload_.n_observed++;

    if (mask == 0) { payload_.n_oor++; return; }

    // Per-region counts.
    auto bump_region = [&](int32_t rid){
        payload_.per_region_counts[rid * 8 + mask]++;
    };
    for (int32_t i = 0; i < n_inline; ++i) bump_region(inline_buf_[i]);
    for (int32_t rid : spill_buf_)         bump_region(rid);

    // Pass D: boundary-flux for unspliced single-ref fragments.
    if (splice_type == SPLICE_UNSPLICED) {
        constexpr uint8_t kExonBit = 0b001;
        auto check_one = [&](int32_t rid){
            if ((regions.type_mask(rid) & kExonBit) == 0) return;
            const int64_t rs = regions.start(rid);
            const int64_t re = regions.end(rid);
            if (frag_start <  rs && frag_end > rs) payload_.u_left[rid]++;
            if (frag_start <  re && frag_end > re) payload_.u_right[rid]++;
        };
        for (int32_t i = 0; i < n_inline; ++i) check_one(inline_buf_[i]);
        for (int32_t rid : spill_buf_)         check_one(rid);
    }
}
```

The whole hot path is one `std::upper_bound` per exon block + a constant
amount of arithmetic per overlapped region. No allocations after the
spill buffer reaches steady state (which it almost never needs to).

### 3.4 `set_regions` Python-callable

```cpp
// In BamScanner public:
void set_regions(nb::ndarray<int32_t, nb::ndim<1>, nb::c_contig> ref_ids,
                 nb::ndarray<int64_t, nb::ndim<1>, nb::c_contig> starts,
                 nb::ndarray<int64_t, nb::ndim<1>, nb::c_contig> ends,
                 nb::ndarray<uint8_t, nb::ndim<1>, nb::c_contig> type_masks,
                 int32_t n_refs)
{
    // Length agreement
    const int64_t n = static_cast<int64_t>(ref_ids.shape(0));
    if (starts.shape(0) != n || ends.shape(0) != n || type_masks.shape(0) != n)
        throw std::invalid_argument("set_regions: array length mismatch");
    if (region_index_) throw std::runtime_error(
        "set_regions: regions already set; create a new BamScanner");

    region_index_ = std::make_unique<RegionIndex>();
    region_index_->set(ref_ids.data(), starts.data(), ends.data(),
                       type_masks.data(), n, n_refs);
}
```

**Caller-side (in `pipeline.py`, M3 wiring):**

```python
region_df = index.region_df
ref_ids = region_df["ref_name"].map(index.ref_name_to_id).to_numpy(np.int32)
starts  = region_df["start"].to_numpy(np.int64)
ends    = region_df["end"].to_numpy(np.int64)
# Plan §2.3: type_mask = 1 << (2 - type)
type_masks = (1 << (2 - region_df["type"].to_numpy())).astype(np.uint8)
scanner.set_regions(ref_ids, starts, ends, type_masks,
                    n_refs=len(index.ref_names))
```

### 3.5 Worker integration (`bam_scanner.cpp` deltas)

Three changes to `WorkerState`:
```cpp
struct WorkerState {
    ResolverScratch scratch;
    FragmentAccumulator accumulator;
    BamScanStats stats;
    StrandObservations strand_obs;
    FragLenObservations fraglen_obs;
    rigel::calibration::CalibrationAccumulator cal_acc;  // NEW

    WorkerState(int32_t n_transcripts, int64_t n_regions)
        : scratch(n_transcripts), cal_acc(n_regions) {}
};
```

The single observation site goes at the **end** of
`process_qname_group_threaded`, after all hits and all counters have
been resolved. We pass `region_index` by reference into the function
signature (already shares `ctx` the same way).

Decision tree at the observation site (executed once per qname group):

```cpp
// Observation policy (§2.10) — exactly one branch fires per molecule.
if (is_multimap) {
    cal_acc.note_multimap();
} else if (any_hit_chimeric && !any_non_chimeric_resolved) {
    cal_acc.note_chimera();
} else {
    // Find the fragment to observe: the resolved non-chimeric hit if
    // any, else the (single) intergenic fragment.
    const AssembledFragment* obs_frag = ...;       // tracked in loop
    int8_t obs_splice = ...;                        // tracked in loop
    if (obs_splice == SPLICE_ARTIFACT) {
        cal_acc.note_artifact();
    } else if (obs_frag != nullptr) {
        const int32_t ref_id = obs_frag->exons[0].ref_id;
        const int64_t fs = obs_frag->exons.front().start;
        const int64_t fe = obs_frag->exons.back().end;
        cal_acc.observe(obs_splice, ref_id, fs, fe,
                        obs_frag->exons.data(),
                        static_cast<int32_t>(obs_frag->exons.size()),
                        *region_index);
    }
}
```

We track `obs_frag` + `obs_splice` inside the existing `for (...) all_hits`
loop. Since unique mappers have at most one hit, this is a single
assignment, not a "best of" reduction. The `obs_frag` pointer is
backed by the loop-scoped `AssembledFragment frag` — care is taken to
hold a copy or to defer observation until after the loop body retains
the right object. **Implementation note:** because `frag` is constructed
fresh each iteration, we copy the `exons` vector once into a
`WorkerState` scratch field (`obs_exons_`) when we see the chosen hit.
This is the only allocation in the observation path; for unique mappers
it happens at most once per molecule.

(Alternatively, since `is_unique_mapper ⇒ all_hits.size() ≤ 1`, we
break out of the loop after handling that single hit, do the
observation, then continue with the post-loop counters. This avoids
the copy entirely. We'll measure both in the perf test and keep the
simpler one.)

### 3.6 Worker merge (after `workers.join()`)

In the existing reader thread, immediately after the per-worker
`stats`/`strand_obs`/`fraglen_obs` merge:

```cpp
for (auto& ws_ptr : worker_states) {
    cal_acc_merged_.merge_from(ws_ptr->cal_acc);
}
```

`merge_from` is element-wise additive over the eight `payload_` fields.
The two parallel arrays (per_region_counts, fl_hist) and the two
counter vectors (u_left, u_right) are sized identically across all
workers (we passed `n_regions` to every WorkerState constructor), so
this is straight `+=` in tight loops.

### 3.7 `build_result` — payload export

```cpp
nb::dict cal_dict;
cal_dict["global_counts"]      = vec_to_ndarray(  // (8,) int64
    std::vector<int64_t>(cal_acc_merged_.payload().global_counts.begin(),
                         cal_acc_merged_.payload().global_counts.end()));
cal_dict["per_region_counts"]  = vec_to_ndarray2d(  // (n_regions, 8) int64
    std::move(cal_acc_merged_.payload().per_region_counts),
    n_regions, 8);
cal_dict["fl_hist"]            = vec_to_ndarray2d(  // (8, 1024) int64
    std::move(cal_acc_merged_.payload().fl_hist), 8, 1024);
cal_dict["u_left"]             = vec_to_ndarray(std::move(payload.u_left));
cal_dict["u_right"]            = vec_to_ndarray(std::move(payload.u_right));
cal_dict["n_observed"]         = payload.n_observed;
cal_dict["n_excluded_multimap"]= payload.n_excluded_multimap;
cal_dict["n_excluded_chimera"] = payload.n_excluded_chimera;
cal_dict["n_excluded_artifact"]= payload.n_excluded_artifact;
cal_dict["n_oor"]              = payload.n_oor;
result["calibration"] = cal_dict;
```

We already have `vec_to_ndarray` (1D); we add a one-line `vec_to_ndarray2d`
helper to `ndarray_util.h` that wraps the same buffer with shape `(M,N)`.

If `set_regions` was never called, `result["calibration"] = nb::none()`.
The Python side treats this as "no calibration payload available" and
short-circuits downstream consumers. (This branch only matters for
direct-from-Python `BamScanner` users; the `pipeline.py` path always
calls `set_regions`.)

### 3.8 `CalibrationScanPayload` (Python)

```python
# src/rigel/calibration/scan_payload.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping
import numpy as np


@dataclass(frozen=True, slots=True)
class CalibrationScanPayload:
    global_counts: np.ndarray            # shape (8,),       int64
    per_region_counts: np.ndarray        # shape (R, 8),     int64
    fl_hist: np.ndarray                  # shape (8, 1024),  int64
    u_left: np.ndarray                   # shape (R,),       int64
    u_right: np.ndarray                  # shape (R,),       int64
    n_observed: int
    n_excluded_multimap: int
    n_excluded_chimera: int
    n_excluded_artifact: int
    n_oor: int

    @classmethod
    def from_scan_dict(
        cls,
        d: Mapping,
        *,
        n_total: int | None = None,
        n_unresolved: int | None = None,
    ) -> "CalibrationScanPayload":
        """Validate shapes + dtypes, then construct.

        If ``n_total`` is provided, asserts the balance:
            n_observed + n_excluded_multimap + n_excluded_chimera
            + n_excluded_artifact + (n_unresolved or 0) == n_total

        On any shape/dtype/balance violation: ValueError with the
        canonical "Rebuild the index or rerun the scan" message.
        """
        # Normalize numpy + check dtypes/shapes...
        # Build instance, then balance check.
```

The balance assertion is the single most important downstream contract.
Any future regression that double-counts a fragment in the per-fragment
loop is caught at payload construction time, with one human-readable
error message — not 30 commits later in the EM convergence behaviour.

---

## 4. CMakeLists wiring

```cmake
# In src/rigel/native/CMakeLists.txt or root CMakeLists:
nanobind_add_module(_bam_impl STABLE_ABI ...
    bam_scanner.cpp
    calibration/accumulator.cpp        # NEW
)
target_include_directories(_bam_impl PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}        # existing
    # calibration/*.h are picked up via "calibration/foo.h" relative
)
```

No extra link dependencies. No new compile flags — we keep the existing
`-O3` + LTO from the parent target.

---

## 5. Tests (`tests/test_calibration_accumulator.py`)

Nine focused tests that hit the unique failure modes:

1. **`test_set_regions_basic`** — `BamScanner.set_regions` accepts the
   four expected dtypes, rejects mismatched lengths, rejects a second
   call on the same scanner instance.

2. **`test_observation_policy_excludes_multimapper`** — synthetic
   2-record qname group with `NH=2`; assert `n_excluded_multimap == 1`,
   `n_observed == 0`, all per-region counts zero.

3. **`test_observation_policy_excludes_chimera`** — qname group whose
   only hit resolves chimeric (trans). Assert `n_excluded_chimera == 1`.

4. **`test_observation_policy_excludes_artifact`** — qname group with
   a SJ on the splice blacklist (becomes `SPLICE_ARTIFACT`). Assert
   `n_excluded_artifact == 1`.

5. **`test_8_state_mask_correctness`** — fixture: 4 regions tiling
   a single 2 kb reference (INTERGENIC[0,500), EXON[500,800),
   INTRON[800,1200), EXON[1200,2000)). Synthesize one fragment per
   mask using hand-built BAM records; assert
   `payload.global_counts == expected_per_mask_vector`.

6. **`test_boundary_flux_counters`** — single unspliced fragment that
   straddles the EXON[500,800) left boundary. Assert
   `u_left[1] == 1, u_right[1] == 0`. Repeat for right-straddle and
   for full-region-spanning (both flags `==1`).

7. **`test_per_mask_fl_histogram`** — three fragments of length 150,
   200, 250 mapping to EXON_ONLY. Assert `fl_hist[1, 150] == 1`,
   `fl_hist[1, 200] == 1`, `fl_hist[1, 250] == 1`, all other rows zero.

8. **`test_worker_merge_equality`** — same BAM, run scan with
   `n_workers ∈ {1, 4}`. Assert all numpy arrays in
   `payload` are exactly equal across the two runs (no
   `assert_allclose` — exact byte equality on int64 counters).

9. **`test_payload_balance_assertion`** — tampered scan dict where
   `n_observed + excluded counters ≠ n_total - n_unresolved`;
   `CalibrationScanPayload.from_scan_dict` raises `ValueError`.

A 10th test, **`test_pipeline_smoke_attaches_payload`**, lives in
`tests/test_pipeline_smoke.py` (or new `test_pipeline_calibration.py`
if the smoke test would grow too unwieldy): runs the existing mini
pipeline end-to-end and asserts
`pipeline_result.calibration_payload is not None`.

---

## 6. Pipeline integration (minimal at M3)

The plan §4-M3 says "`PipelineResult.calibration_payload` populated
end-to-end". We implement the smallest hookup that satisfies this:

1. Add `calibration_payload: CalibrationScanPayload | None` slot to
   the existing `PipelineResult` (or whatever `scan_and_buffer`
   returns — TBD by code reading; either a new tuple element or a
   field).
2. In `scan_and_buffer`, after building `scanner` and *before*
   `scanner.scan(...)`:
   ```python
   _wire_calibration_regions(scanner, index)   # the snippet from §3.4
   ```
3. After `scanner.scan(...)` returns, build the payload:
   ```python
   cal_dict = scan_result.get("calibration")
   payload = (
       CalibrationScanPayload.from_scan_dict(
           cal_dict, n_total=stats.total, n_unresolved=...)
       if cal_dict is not None else None
   )
   ```

No M3 consumer reads `calibration_payload` yet — M4 will. The hookup
exists at M3 only so the smoke test can pin that the pipeline glue
works end-to-end.

---

## 7. What we explicitly DON'T do at M3

- **No global-density computation.** That's M4.
- **No per-locus prior assembly.** That's M6.
- **No `pool_fl_models`.** That's M7.
- **No `boundary_tol` parameter.** Hard-coded `0`; revisited if a real
  benchmark shows it matters.
- **No multimapper observation modes.** Multimappers are excluded
  unconditionally per §2.10.
- **No mappability adjustments to `L_eff(R)`.** Documented future work.
- **No SmallRegionSet class.** Two scratch buffers on
  `CalibrationAccumulator` do the same job with less code.

---

## 8. Risks specific to M3

| # | Risk | Mitigation |
|---|---|---|
| 1 | `set_regions` ABI drifts before M4 reads it | Single caller (`pipeline.py`), single test (`test_set_regions_basic`); no third-party call sites until M8 |
| 2 | Per-fragment dedup is O(k²) and k blows up on real loci | Test 5 verifies on 4-region fixture; if real data shows k > 16 we promote to a hash set; instrumented via a `max_k_observed` debug counter |
| 3 | Worker merge silently truncates | `test_worker_merge_equality` is exact byte equality, not approximate |
| 4 | Observation site fires on chimeric/artifact fragments | Three dedicated exclusion tests (#2, #3, #4) cover every ineligible class one fragment at a time |
| 5 | `obs_frag` pointer dangles after loop exits | Implementation note in §3.5: copy exons into a WorkerState scratch field, OR observe inside the loop for unique mappers (which have ≤ 1 hit). Code review checks both options before commit |
| 6 | Synthetic BAM records in tests don't exercise the splice-artifact path because the blacklist is empty | Test 4 builds a tiny `splice_blacklist.feather` with one entry and confirms an artifact-classified fragment is counted, not observed |

---

## 9. Exit gate (recap from parent plan §4-M3)

- New tests #1–#9 green.
- `test_pipeline_smoke.py::test_pipeline_calibration_payload_populated` green.
- Protected suite (parent §6) green.
- `intervals.feather`, `regions.feather`, EM golden outputs unchanged.
- `pip install --no-build-isolation -e .` succeeds without warnings on
  macOS and Linux (the project's two CI targets).
- Single commit. Diff size estimate: ~+800 LoC C++, ~+250 LoC tests,
  ~+30 LoC Python wiring.

---

## 10. Commit discipline reminders (from parent plan §3)

- Commit M3 as exactly one PR-sized commit.
- Stage all new files (including `accumulator.cpp`,
  `scan_payload.py`, the test file, and the `CMakeLists.txt` change)
  immediately after the test suite goes green.
- If the worker-merge test fails after a code edit, do NOT widen its
  tolerance — investigate. It is the only sentinel for parallelism
  bugs in the calibration path.
- No bulk regex on `src/`. No exception.
