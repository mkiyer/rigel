/**
 * accumulator.h — Per-worker calibration observation accumulator.
 *
 * The hot path (``observe()``) is allocation-free after warm-up: it
 * uses two reusable scratch ``std::vector<int32_t>`` buffers (per-block
 * hits + persistent fragment-scope hits) whose ``clear()`` retains
 * capacity.  Per-fragment hit counts are O(10), so working set stays
 * well within L1.
 *
 * M3 invariants this code relies on:
 *   1. Regions are per-ref non-overlapping (M2): each ``overlap_into``
 *      call returns strictly increasing region IDs with no internal
 *      duplicates.
 *   2. Per-fragment hit counts are O(10).
 *   3. The scanner gates non-observable fragments (multimappers,
 *      chimeras, splice artifacts) before calling ``observe()``.
 *   4. ``payload.fl_hist`` records ``frag_end - frag_start`` (genomic
 *      span, not fragment length); see calibration plan §M7 for the
 *      consumer-side semantics.
 *
 * Mask layout (single source of truth lives in ``region_index.h``):
 *   bit 0 = EXON, bit 1 = INTRON, bit 2 = INTERGENIC.
 */

#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "constants.h"
#include "calibration/region_index.h"

namespace rigel::calibration {

struct CalibrationPayload {
    static constexpr int32_t kFlBins = 1024;

    // Global per-mask counters (size mask::N_STATES).
    std::array<int64_t, mask::N_STATES> global_counts {};

    // Per-region per-mask counts: shape (n_regions, mask::N_STATES) row-major.
    std::vector<int64_t> per_region_counts;

    // Per-mask FL histogram: shape (mask::N_STATES, kFlBins) row-major.
    std::vector<int64_t> fl_hist;

    // Boundary-flux counters (only EXON regions get nonzero values,
    // but indexed by region_id for O(1) lookup).
    std::vector<int64_t> u_left;
    std::vector<int64_t> u_right;

    // Counters
    int64_t n_observed             = 0;
    int64_t n_excluded_multimap    = 0;
    int64_t n_excluded_chimera     = 0;
    int64_t n_excluded_artifact    = 0;
    int64_t n_unobserved           = 0;  // non-mm group, no eligible hit
    // Observed fragment whose every block fell outside any region.
    // Under the M2 tiling invariant this is only possible for
    // alignments to references the calibration partition didn't cover
    // (e.g. decoys / contigs / index-mismatched BAMs).  Useful as a QC
    // signal — large values mean the BAM was aligned to a different
    // reference than the index was built from.
    int64_t n_unannotated_ref      = 0;
    // Observed fragment with raw region hits whose every aligned-block
    // overlap was below the boundary-tolerance threshold q(K) = max(K,1)
    // (so no hit contributed a type bit to ``obs_mask``). These
    // fragments still increment ``global_counts[0]``, ``fl_hist[0]``,
    // and ``n_observed`` — this is an auxiliary QC subcounter.
    // Always zero when boundary_tolerance == 0.
    int64_t n_below_tolerance      = 0;
};

class CalibrationAccumulator {
public:
    explicit CalibrationAccumulator(int64_t n_regions,
                                    int32_t boundary_tolerance = 0) {
        if (boundary_tolerance < 0) {
            throw std::invalid_argument(
                "CalibrationAccumulator: boundary_tolerance must be >= 0");
        }
        boundary_tolerance_ = boundary_tolerance;
        payload_.per_region_counts.assign(
            static_cast<size_t>(n_regions) * mask::N_STATES, 0);
        payload_.fl_hist.assign(
            static_cast<size_t>(mask::N_STATES) *
                CalibrationPayload::kFlBins,
            0);
        payload_.u_left.assign(static_cast<size_t>(n_regions), 0);
        payload_.u_right.assign(static_cast<size_t>(n_regions), 0);
        n_regions_ = n_regions;
        // Warm reserve to amortize the first ~10 fragments' allocations.
        hits_.reserve(16);
        hits_overlap_bp_.reserve(16);
        block_hits_.reserve(16);
        block_overlap_bp_.reserve(16);
        merge_buf_.reserve(16);
        merge_buf_bp_.reserve(16);
        qualified_hits_.reserve(16);
    }

    /**
     * Hot-path observation. ``exons`` is the chosen fragment's
     * R1∪R2 exon blocks (sorted, single-ref). The caller must have
     * gated on: unique mapper, non-chimeric, non-artifact splice.
     */
    void observe(int8_t splice_type,
                 int32_t ref_id,
                 int64_t frag_start,
                 int64_t frag_end,
                 const ExonBlock* exons,
                 int32_t n_exons,
                 const RegionIndex& regions);

    inline void note_multimap()  { ++payload_.n_excluded_multimap; }
    inline void note_chimera()   { ++payload_.n_excluded_chimera;  }
    inline void note_artifact()  { ++payload_.n_excluded_artifact; }
    inline void note_unobserved(){ ++payload_.n_unobserved;        }

    /// Worker merge: element-wise additive into *this.
    void merge_from(const CalibrationAccumulator& other);

    const CalibrationPayload& payload() const { return payload_; }
    CalibrationPayload&       payload()       { return payload_; }

    int64_t n_regions() const { return n_regions_; }
    int32_t boundary_tolerance() const { return boundary_tolerance_; }

private:
    CalibrationPayload payload_;
    int64_t n_regions_ = 0;
    int32_t boundary_tolerance_ = 0;

    // Per-fragment scratch (reused across observe() calls).  After
    // warm-up these never reallocate: clear() keeps capacity.
    //
    // ``hits_overlap_bp_`` tracks per-fragment exact aligned-block
    // overlap (summed across exon blocks) parallel to ``hits_``; it
    // is used to gate ``obs_mask`` and per-region fan-out on the
    // boundary-tolerance threshold ``q(K) = max(K, 1)``.
    // ``qualified_hits_`` holds the subset of ``hits_`` whose
    // overlap_bp >= q(K) — only these contribute mask bits, per-region
    // counts, and boundary-flux events.  At K=0 (q=1) every hit
    // returned by ``RegionIndex::overlap_into`` already has overlap
    // >= 1, so ``qualified_hits_ == hits_`` and behavior collapses to
    // the pre-tolerance code path bit-for-bit.
    std::vector<int32_t> hits_;             // sorted, deduped fragment-scope hits
    std::vector<int64_t> hits_overlap_bp_;  // per-rid summed overlap (parallel to hits_)
    std::vector<int32_t> block_hits_;       // per-exon-block scratch
    std::vector<int64_t> block_overlap_bp_; // per-block overlap (parallel to block_hits_)
    std::vector<int32_t> merge_buf_;        // sorted-merge output buffer (rids)
    std::vector<int64_t> merge_buf_bp_;     // sorted-merge output buffer (overlap)
    std::vector<int32_t> qualified_hits_;   // hits_ filtered by overlap >= q(K)
};

}  // namespace rigel::calibration
