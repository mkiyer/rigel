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
 *             uint32[4]  flux_left, uint32[4] flux_right            (64 B)
 *
 * Channel encoding (4 channels):  ch = (spliced ? 2 : 0) + (primary ? 0 : 1)
 *   ch0 = unspliced & primary    ch1 = unspliced & !primary
 *   ch2 = spliced   & primary    ch3 = spliced   & !primary
 * where the scanner sets `primary` to:
 *   unspliced: read aligned to the '+' genome strand (genome strand);
 *   spliced:   read is SENSE to its splice-motif strand
 *              (align_strand == sj_strand) — transcript-relative.
 * So ch0/ch1 are unspliced genome +/−; ch2/ch3 are spliced sense/antisense.
 * The accumulator is channel-agnostic — the scanner assigns `primary`.
 *
 * Boundary flux is PER SIDE: `flux_left[ch]` / `flux_right[ch]` count
 * fragment-events touching the left / right side of the boundary. An
 * unspliced contiguous crossing credits both sides of its one boundary; a
 * spliced intron-skip credits one side of each flanking boundary (so no false
 * exon-intron flux). Mass (`mass_left`/`mass_right`) is per side as before.
 *
 * The native implementation must match the Python reference byte-for-byte.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rigel::accumulator {

inline constexpr std::size_t kNChannels = 4;

inline int channel_idx(bool spliced, bool primary) noexcept {
    return (spliced ? 2 : 0) + (primary ? 0 : 1);
}

// gDNA fragment-length (FL) pools (PR 4c). Each UNSPLICED fragment's FL mass is
// binned into one of 6 pools = 3 region-types {0=intergenic, 1=intronic,
// 2=exonic} × 2 compartments {0=contained, 1=boundary-crossing}. The gDNA FL is
// later aggregated from the intergenic+intronic pools (both compartments); see
// rigel.calibration.fl. FL pooling is OPTIONAL — an Accumulator built with empty
// region_types or max_fl <= 0 skips it entirely (existing behaviour unchanged).
inline constexpr std::size_t kNRegionTypes    = 3;
inline constexpr std::size_t kNFlCompartments = 2;
inline constexpr std::size_t kNFlPools        = kNRegionTypes * kNFlCompartments;  // 6

inline std::size_t fl_pool_idx(std::uint8_t region_type, bool boundary) noexcept {
    return static_cast<std::size_t>(region_type) * kNFlCompartments + (boundary ? 1u : 0u);
}

struct Region {
    std::uint32_t contained[kNChannels];  // 16 B
};
static_assert(sizeof(Region) == 16, "Region must be 16 bytes");

struct Boundary {
    float          mass_left[kNChannels];   // 16 B
    float          mass_right[kNChannels];  // 16 B
    std::uint32_t  flux_left[kNChannels];   // 16 B
    std::uint32_t  flux_right[kNChannels];  // 16 B
};
static_assert(sizeof(Boundary) == 64, "Boundary must be 64 bytes");

class Accumulator {
public:
    /// Construct with a sorted, strictly increasing array of boundary
    /// positions. `boundaries` is moved into the accumulator. Length must
    /// be >= 1; n_regions is `boundaries.size() - 1`.
    ///
    /// FL pooling (PR 4c) is enabled iff `region_types` is non-empty (its
    /// length must then equal n_regions) AND `max_fl > 0`; otherwise the FL
    /// pools stay empty and deposit() does no FL binning.
    explicit Accumulator(std::vector<std::int64_t> boundaries,
                         std::vector<std::uint8_t> region_types = {},
                         int max_fl = 0);

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

    /// gDNA FL pools (PR 4c). Flat float64, pool-major: pool `p` occupies
    /// `[p*(max_fl+1), (p+1)*(max_fl+1))`, FL bin `min(footprint, max_fl)`.
    /// Empty when FL pooling is disabled.
    const double* fl_pool_data() const noexcept { return fl_pool_mass_.data(); }
    std::size_t   fl_pool_size() const noexcept { return fl_pool_mass_.size(); }
    int           max_fl()       const noexcept { return max_fl_; }

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
                 bool primary);

    /// Element-wise sum of `other` into this accumulator. Requires identical
    /// boundary positions (asserts at start). Used to merge per-worker
    /// accumulators after a parallel scan.
    void merge_from(const Accumulator& other);

private:
    std::vector<std::int64_t> boundary_positions_;  // size = n_regions + 1
    std::vector<Region>       regions_;             // size = n_regions
    std::vector<Boundary>     boundaries_;          // size = n_regions + 1
    std::vector<std::uint8_t> region_types_;        // size = n_regions, or 0 (FL off)
    int                       max_fl_ = 0;          // FL pooling off iff <= 0
    std::vector<double>       fl_pool_mass_;         // kNFlPools*(max_fl+1), or 0
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
    /// FL pooling (PR 4c): `region_types` is the flat ref-major per-region
    /// type array (length R_total = sum of n_regions over refs); it is sliced
    /// per reference and forwarded to each Accumulator with `max_fl`. Pass
    /// `region_types == nullptr` / `max_fl == 0` to disable FL pooling.
    AccumulatorSet(const std::int64_t* boundary_positions,
                   std::size_t n_positions,
                   const std::int64_t* ref_pos_offsets,
                   std::size_t n_refs,
                   const std::uint8_t* region_types = nullptr,
                   std::size_t n_region_types = 0,
                   int max_fl = 0);

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
