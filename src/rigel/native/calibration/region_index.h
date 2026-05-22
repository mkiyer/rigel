/**
 * region_index.h — Per-reference CSR region overlap engine.
 *
 * Given a region partition that is per-reference contiguous and
 * non-overlapping (the M2 invariant enforced by
 * :func:`rigel.calibration.regions.validate_against_ref_lengths`),
 * answer overlap queries in O(log n + k) where n is the number of
 * regions on the queried reference and k is the number of hit regions.
 *
 * The non-overlap invariant guarantees that a single ``overlap_into``
 * call returns at most one hit per region, and the returned region
 * IDs are strictly increasing (regions are stored in start order).
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "calibration/region_signature.h"

namespace rigel::calibration {

// Public mask layout — single source of truth for the 8-state coding.
// Used by RegionIndex (region type masks) and CalibrationAccumulator
// (per-fragment OR'd masks, FL histogram strata).
namespace mask {
constexpr uint8_t EXON       = 0b001;
constexpr uint8_t INTRON     = 0b010;
constexpr uint8_t INTERGENIC = 0b100;
constexpr int     N_STATES   = 8;
}  // namespace mask

class RegionIndex {
public:
    static constexpr uint8_t kSignatureSentinel = 0xFF;
    static constexpr uint8_t kMaxBoundaryKind = 5;

    RegionIndex() = default;

    /**
     * Single setup call. Arrays are copied into private buffers.
     *
     * Inputs MUST be sorted by (ref_id, start) with regions per-ref
         * contiguous and non-overlapping. Fine signatures and boundary
         * metadata are carried for the fractional accumulator migration;
         * ``type_masks`` and ``strands`` remain the Phase 2 bridge used by
         * the integer accumulator.
     */
    void set(const int32_t* ref_ids,
             const int64_t* starts,
             const int64_t* ends,
             const uint8_t* signatures,
             const uint8_t* left_signatures,
             const uint8_t* right_signatures,
             const uint8_t* boundary_kind_left,
             const uint8_t* boundary_kind_right,
             const uint8_t* type_masks,
             const uint8_t* strands,
             int64_t n_regions,
             int32_t n_refs)
    {
        if (n_refs < 0)
            throw std::invalid_argument("RegionIndex::set: n_refs < 0");
        if (n_regions < 0)
            throw std::invalid_argument("RegionIndex::set: n_regions < 0");

        starts_.assign(starts, starts + n_regions);
        ends_.assign(ends, ends + n_regions);
    signatures_.assign(signatures, signatures + n_regions);
    left_signatures_.assign(left_signatures, left_signatures + n_regions);
    right_signatures_.assign(right_signatures, right_signatures + n_regions);
    boundary_kind_left_.assign(boundary_kind_left, boundary_kind_left + n_regions);
    boundary_kind_right_.assign(boundary_kind_right, boundary_kind_right + n_regions);
        type_masks_.assign(type_masks, type_masks + n_regions);
        strands_.assign(strands, strands + n_regions);
        n_refs_ = n_refs;

        // Build per-ref CSR offsets.  Inputs are sorted by (ref_id, start),
        // so we can scan once.
        ref_offsets_.assign(static_cast<size_t>(n_refs) + 1, 0);
        int32_t cur_ref = 0;
        for (int64_t i = 0; i < n_regions; ++i) {
            int32_t r = ref_ids[i];
            if (r < 0 || r >= n_refs) {
                throw std::invalid_argument(
                    "RegionIndex::set: ref_id out of range");
            }
            if (r < cur_ref) {
                throw std::invalid_argument(
                    "RegionIndex::set: regions not sorted by ref_id");
            }
            if (ends[i] <= starts[i]) {
                throw std::invalid_argument(
                    "RegionIndex::set: region end must be greater than start");
            }
            if (i > 0 && r == ref_ids[i - 1]) {
                if (starts[i] < starts[i - 1]) {
                    throw std::invalid_argument(
                        "RegionIndex::set: regions not sorted by start within ref");
                }
                if (starts[i] < ends[i - 1]) {
                    throw std::invalid_argument(
                        "RegionIndex::set: regions overlap within ref");
                }
            }
            if (!valid_signature(signatures[i])) {
                throw std::invalid_argument(
                    "RegionIndex::set: signature out of range");
            }
            if (!valid_neighbor_signature(left_signatures[i]) ||
                !valid_neighbor_signature(right_signatures[i])) {
                throw std::invalid_argument(
                    "RegionIndex::set: neighbor signature out of range");
            }
            if (boundary_kind_left[i] > kMaxBoundaryKind ||
                boundary_kind_right[i] > kMaxBoundaryKind) {
                throw std::invalid_argument(
                    "RegionIndex::set: boundary kind out of range");
            }
            while (cur_ref < r) {
                ref_offsets_[cur_ref + 1] = i;
                ++cur_ref;
            }
        }
        while (cur_ref < n_refs) {
            ref_offsets_[cur_ref + 1] = n_regions;
            ++cur_ref;
        }
    }

    /**
     * O(log n + k) overlap query. Region IDs that overlap [qstart, qend)
     * on ``ref_id`` are *appended* to ``out`` in strictly increasing
     * order.  The caller owns ``out`` and is responsible for clearing
     * (or otherwise managing) it between calls.
     */
    void overlap_into(int32_t ref_id,
                      int64_t qstart,
                      int64_t qend,
                      std::vector<int32_t>& out) const
    {
        if (ref_id < 0 || ref_id >= n_refs_) return;
        const int64_t lo = ref_offsets_[ref_id];
        const int64_t hi = ref_offsets_[ref_id + 1];
        if (lo == hi || qend <= qstart) return;

        // Find first region whose start > qstart.  The region just
        // before it is a candidate iff its end > qstart.
        const int64_t* sp = starts_.data();
        auto it = std::upper_bound(sp + lo, sp + hi, qstart);
        int64_t i = static_cast<int64_t>(it - sp);
        if (i > lo && ends_[i - 1] > qstart) i = i - 1;

        while (i < hi && starts_[i] < qend) {
            out.push_back(static_cast<int32_t>(i));
            ++i;
        }
    }

    inline int64_t  start(int32_t rid)     const { return starts_[rid]; }
    inline int64_t  end(int32_t rid)       const { return ends_[rid]; }
    inline uint8_t  signature(int32_t rid) const { return signatures_[rid]; }
    inline uint8_t  left_signature(int32_t rid) const { return left_signatures_[rid]; }
    inline uint8_t  right_signature(int32_t rid) const { return right_signatures_[rid]; }
    inline uint8_t  boundary_kind_left(int32_t rid) const { return boundary_kind_left_[rid]; }
    inline uint8_t  boundary_kind_right(int32_t rid) const { return boundary_kind_right_[rid]; }
    inline uint8_t  type_mask(int32_t rid) const { return type_masks_[rid]; }
    inline uint8_t  strand(int32_t rid)    const { return strands_[rid]; }
    inline int64_t  n_regions()            const { return static_cast<int64_t>(starts_.size()); }
    inline int32_t  n_refs()               const { return n_refs_; }

private:
    static constexpr bool valid_signature(uint8_t value) {
        return value < region_signature::kNSignatures;
    }

    static constexpr bool valid_neighbor_signature(uint8_t value) {
        return valid_signature(value) || value == kSignatureSentinel;
    }

    std::vector<int64_t> starts_;
    std::vector<int64_t> ends_;
    std::vector<uint8_t> signatures_;
    std::vector<uint8_t> left_signatures_;
    std::vector<uint8_t> right_signatures_;
    std::vector<uint8_t> boundary_kind_left_;
    std::vector<uint8_t> boundary_kind_right_;
    std::vector<uint8_t> type_masks_;
    std::vector<uint8_t> strands_;
    std::vector<int64_t> ref_offsets_;   // size = n_refs_ + 1
    int32_t n_refs_ = 0;
};

}  // namespace rigel::calibration
