/**
 * accumulator.h — Fractional accumulator (per-reference region/boundary split).
 *
 * Canonical spec: docs/accumulator/00_design.md
 * Python reference: tests/native/_accumulator_reference.py
 *
 * One Accumulator describes ONE reference. Construction takes the sorted
 * boundary-position array of length N+1 (int64 genomic coordinates), giving
 * N contiguous regions and N+1 boundaries:
 *
 *     regions[i]    = [boundaries[i], boundaries[i+1])
 *     boundaries[i] is positioned at the i-th coordinate
 *
 * Storage:
 *   Region:   uint32[4] contained                                  (16 B)
 *   Boundary: float32[4] mass_left, float32[4] mass_right,
 *             uint32[4]  flux                                       (48 B)
 *
 * Channel encoding (4 channels):
 *   ch = (spliced ? 2 : 0) + (strand_pos ? 0 : 1)
 *
 * The native implementation must match the Python reference byte-for-byte.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rigel::accumulator {

inline constexpr std::size_t kNChannels = 4;

inline int channel_idx(bool spliced, bool strand_pos) noexcept {
    return (spliced ? 2 : 0) + (strand_pos ? 0 : 1);
}

struct Region {
    std::uint32_t contained[kNChannels];  // 16 B
};
static_assert(sizeof(Region) == 16, "Region must be 16 bytes");

struct Boundary {
    float          mass_left[kNChannels];   // 16 B
    float          mass_right[kNChannels];  // 16 B
    std::uint32_t  flux[kNChannels];        // 16 B
};
static_assert(sizeof(Boundary) == 48, "Boundary must be 48 bytes");

class Accumulator {
public:
    /// Construct with a sorted, strictly increasing array of boundary
    /// positions. `boundaries` is moved into the accumulator. Length must
    /// be >= 1; n_regions is `boundaries.size() - 1`.
    explicit Accumulator(std::vector<std::int64_t> boundaries);

    std::size_t n_regions()    const noexcept { return regions_.size(); }
    std::size_t n_boundaries() const noexcept { return boundaries_.size(); }

    const std::int64_t* boundary_positions() const noexcept {
        return boundary_positions_.data();
    }
    std::size_t n_boundary_positions() const noexcept {
        return boundary_positions_.size();
    }

    Region*         regions_data()        noexcept { return regions_.data(); }
    Boundary*       boundaries_data()     noexcept { return boundaries_.data(); }
    const Region*   regions_data()   const noexcept { return regions_.data(); }
    const Boundary* boundaries_data()const noexcept { return boundaries_.data(); }

    /// Region index containing `pos`, or -1 if `pos` is outside
    /// [boundaries.front(), boundaries.back()).
    std::int64_t region_of_pos(std::int64_t pos) const noexcept;

    /// Deposit a single fragment's evidence per docs/accumulator/00_design.md.
    /// `block_starts`/`block_ends` have length `n_blocks`. Empty / fully
    /// out-of-range fragments are no-ops.
    void deposit(const std::int64_t* block_starts,
                 const std::int64_t* block_ends,
                 std::size_t n_blocks,
                 bool spliced,
                 bool strand_pos);

    /// Element-wise sum of `other` into this accumulator. Requires identical
    /// boundary positions (asserts at start). Used to merge per-worker
    /// accumulators after a parallel scan.
    void merge_from(const Accumulator& other);

private:
    std::vector<std::int64_t> boundary_positions_;  // size = n_regions + 1
    std::vector<Region>       regions_;             // size = n_regions
    std::vector<Boundary>     boundaries_;          // size = n_regions + 1
};

// ============================================================================
// AccumulatorSet — one Accumulator per reference, flat shared partition.
// ============================================================================
//
// Inputs:
//   - `boundary_positions`: flat int64 array of size B_pos_total holding the
//     concatenated boundary positions for all references.
//   - `ref_pos_offsets`: int64 array of size n_refs + 1; the boundary
//     positions for reference f live in
//         boundary_positions[ref_pos_offsets[f] .. ref_pos_offsets[f+1])
//     and define n_regions_f = (ref_pos_offsets[f+1] - ref_pos_offsets[f]) - 1.
//   - References with fewer than 2 positions (n_regions == 0) are valid and
//     produce an empty Accumulator with no regions/boundaries.
//
// Total counts:
//   R_total      = sum over refs of n_regions_f
//   B_obj_total  = sum over refs of (n_regions_f + 1) = R_total + n_refs
//                  for refs with >=1 region (refs with 0 regions contribute 0).
//
class AccumulatorSet {
public:
    AccumulatorSet(const std::int64_t* boundary_positions,
                   std::size_t n_positions,
                   const std::int64_t* ref_pos_offsets,
                   std::size_t n_refs);

    /// Number of references managed.
    std::size_t n_refs() const noexcept { return accs_.size(); }

    /// Reference accumulator at index `ref_id` (0-based).
    Accumulator&       at(std::int32_t ref_id);
    const Accumulator& at(std::int32_t ref_id) const;

    /// Element-wise merge of `other` into `this`. Requires the per-ref
    /// boundary positions to match exactly.
    void merge_from(const AccumulatorSet& other);

private:
    std::vector<Accumulator> accs_;
};

}  // namespace rigel::accumulator
