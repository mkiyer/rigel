/**
 * region_index.h — Per-reference CSR region overlap engine.
 *
 * Given a region partition that is per-reference contiguous and
 * non-overlapping (the M2 invariant enforced by
 * :func:`rigel.calibration.regions.validate_against_ref_lengths`),
 * answer overlap queries in O(log n + k) where n is the number of
 * regions on the queried reference and k is the number of hit regions.
 *
 * The non-overlap invariant is what makes this so cheap — a query
 * never returns the same region twice, so no dedup pass is needed.
 *
 * See ``docs/calibration/m3_implementation_plan.md`` §3.1.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace rigel::calibration {

class RegionIndex {
public:
    RegionIndex() = default;

    /**
     * Single setup call. Arrays are copied into private buffers.
     *
     * Inputs MUST be sorted by (ref_id, start) with regions per-ref
     * contiguous and non-overlapping. ``type_masks[i]`` is the
     * pre-computed bit mask for region ``i`` (bit 0 = EXON,
     * bit 1 = INTRON, bit 2 = INTERGENIC).
     */
    void set(const int32_t* ref_ids,
             const int64_t* starts,
             const int64_t* ends,
             const uint8_t* type_masks,
             int64_t n_regions,
             int32_t n_refs)
    {
        if (n_refs < 0)
            throw std::invalid_argument("RegionIndex::set: n_refs < 0");
        if (n_regions < 0)
            throw std::invalid_argument("RegionIndex::set: n_regions < 0");

        starts_.assign(starts, starts + n_regions);
        ends_.assign(ends, ends + n_regions);
        type_masks_.assign(type_masks, type_masks + n_regions);
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
     * O(log n + k) overlap query. Region indices that overlap
     * [qstart, qend) on ``ref_id`` are appended into ``out_inline``
     * up to capacity N, then to ``out_spill``. Returns total count.
     *
     * Output is in start order (regions are stored sorted).
     */
    template <std::size_t N>
    int32_t overlap_into(int32_t ref_id,
                         int64_t qstart,
                         int64_t qend,
                         std::array<int32_t, N>& out_inline,
                         std::vector<int32_t>& out_spill,
                         int32_t already_inline) const
    {
        if (ref_id < 0 || ref_id >= n_refs_) return 0;
        const int64_t lo = ref_offsets_[ref_id];
        const int64_t hi = ref_offsets_[ref_id + 1];
        if (lo == hi || qend <= qstart) return 0;

        // Find first region whose start > qstart.  The region just
        // before it is a candidate iff its end > qstart.
        const int64_t* sp = starts_.data();
        auto it = std::upper_bound(sp + lo, sp + hi, qstart);
        int64_t i = static_cast<int64_t>(it - sp);
        if (i > lo && ends_[i - 1] > qstart) i = i - 1;

        int32_t n = 0;
        int32_t inline_pos = already_inline;
        while (i < hi && starts_[i] < qend) {
            const int32_t rid = static_cast<int32_t>(i);
            if (inline_pos < static_cast<int32_t>(N)) {
                out_inline[inline_pos++] = rid;
            } else {
                out_spill.push_back(rid);
            }
            ++n;
            ++i;
        }
        return n;
    }

    inline int64_t  start(int32_t rid)     const { return starts_[rid]; }
    inline int64_t  end(int32_t rid)       const { return ends_[rid]; }
    inline uint8_t  type_mask(int32_t rid) const { return type_masks_[rid]; }
    inline int64_t  n_regions()            const { return static_cast<int64_t>(starts_.size()); }
    inline int32_t  n_refs()               const { return n_refs_; }

private:
    std::vector<int64_t> starts_;
    std::vector<int64_t> ends_;
    std::vector<uint8_t> type_masks_;
    std::vector<int64_t> ref_offsets_;   // size = n_refs_ + 1
    int32_t n_refs_ = 0;
};

}  // namespace rigel::calibration
