/**
 * accumulator.h — Per-worker calibration observation accumulator.
 *
 * The hot path (``observe()``) is allocation-free after warm-up: it
 * uses two reusable scratch buffers (inline + spill) for region IDs.
 *
 * 8-state mask bit layout (matches Python plan §2.3):
 *   bit 0 = EXON
 *   bit 1 = INTRON
 *   bit 2 = INTERGENIC
 * The ``type_mask`` stored in :class:`RegionIndex` is already the
 * single-bit mask for that region's type, so OR-ing them produces
 * the union mask directly.
 *
 * See ``docs/calibration/m3_implementation_plan.md`` §3.3.
 */

#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "constants.h"
#include "calibration/region_index.h"

namespace rigel::calibration {

struct CalibrationPayload {
    static constexpr int32_t kFlBins = 1024;
    static constexpr int32_t kMaskCard = 8;

    // Global per-mask counters (size 8).
    std::array<int64_t, kMaskCard> global_counts {};

    // Per-region per-mask counts: shape (n_regions, 8) row-major.
    std::vector<int64_t> per_region_counts;

    // Per-mask FL histogram: shape (8, kFlBins) row-major.
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
    int64_t n_oor                  = 0;
};

class CalibrationAccumulator {
public:
    explicit CalibrationAccumulator(int64_t n_regions) {
        payload_.per_region_counts.assign(
            static_cast<size_t>(n_regions) * CalibrationPayload::kMaskCard, 0);
        payload_.fl_hist.assign(
            static_cast<size_t>(CalibrationPayload::kMaskCard) *
                CalibrationPayload::kFlBins,
            0);
        payload_.u_left.assign(static_cast<size_t>(n_regions), 0);
        payload_.u_right.assign(static_cast<size_t>(n_regions), 0);
        n_regions_ = n_regions;
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

private:
    CalibrationPayload payload_;
    int64_t n_regions_ = 0;

    // Per-fragment scratch (reused across observe() calls)
    static constexpr std::size_t kInlineCap = 16;
    std::array<int32_t, kInlineCap> inline_buf_ {};
    std::vector<int32_t> spill_buf_;
};

}  // namespace rigel::calibration
