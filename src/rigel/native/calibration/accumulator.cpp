/**
 * accumulator.cpp — Fractional accumulator implementation.
 *
 * Spec: docs/accumulator/00_design.md
 * Python reference: tests/native/_accumulator_reference.py
 *
 * Must match the Python reference byte-for-byte (float32 masses,
 * uint32 counts).
 */
#include "accumulator.h"

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <utility>

namespace rigel::accumulator {

namespace {

// Per-slice scratch entry: (region_idx, start, end).
struct Slice {
    std::int64_t region_idx;
    std::int64_t start;
    std::int64_t end;
};

}  // anonymous namespace

Accumulator::Accumulator(std::vector<std::int64_t> boundaries,
                         std::vector<std::uint8_t> region_types,
                         int max_fl)
    : boundary_positions_(std::move(boundaries)),
      region_types_(std::move(region_types)),
      max_fl_(max_fl)
{
    // An empty position array is permitted and yields a no-op accumulator
    // (n_regions == 0, n_boundaries == 0). This lets AccumulatorSet allocate
    // a placeholder for references that have no annotated partition.
    if (boundary_positions_.empty()) {
        regions_.clear();
        boundaries_.clear();
        region_types_.clear();
        max_fl_ = 0;
        return;
    }
    // Strict monotonicity check.
    for (std::size_t i = 1; i < boundary_positions_.size(); ++i) {
        if (boundary_positions_[i] <= boundary_positions_[i - 1]) {
            throw std::invalid_argument(
                "Accumulator: boundaries must be strictly increasing");
        }
    }
    const std::size_t n = boundary_positions_.size() - 1;
    regions_.assign(n, Region{});
    boundaries_.assign(n + 1, Boundary{});
    // POD zero-init via value-init in assign above; be explicit for clarity.
    if (n > 0) std::memset(regions_.data(), 0, n * sizeof(Region));
    std::memset(boundaries_.data(), 0, (n + 1) * sizeof(Boundary));

    // FL pooling (PR 4c) is enabled iff a per-region type array of the right
    // length is supplied AND max_fl > 0; otherwise it is a no-op.
    if (max_fl_ > 0 && region_types_.size() == n) {
        fl_pool_mass_.assign(kNFlPools * (static_cast<std::size_t>(max_fl_) + 1), 0.0);
    } else {
        region_types_.clear();
        max_fl_ = 0;
        fl_pool_mass_.clear();
    }
}

std::int64_t Accumulator::region_of_pos(std::int64_t pos) const noexcept {
    if (boundary_positions_.size() < 2) return -1;
    const std::int64_t lo = boundary_positions_.front();
    const std::int64_t hi = boundary_positions_.back();
    if (pos < lo || pos >= hi) return -1;
    // np.searchsorted(boundaries, pos, side='right') - 1
    auto it = std::upper_bound(boundary_positions_.begin(),
                               boundary_positions_.end(), pos);
    return static_cast<std::int64_t>(
               std::distance(boundary_positions_.begin(), it)) - 1;
}

void Accumulator::deposit(const std::int64_t* block_starts,
                          const std::int64_t* block_ends,
                          std::size_t n_blocks,
                          bool spliced,
                          bool primary)
{
    if (n_blocks == 0 || regions_.empty()) return;

    const int ch = channel_idx(spliced, primary);
    const std::int64_t edge_lo = boundary_positions_.front();
    const std::int64_t edge_hi = boundary_positions_.back();

    // 1. Expand each block into per-region slices.
    //    Worst case: each block straddles every region. Reserve a small
    //    initial capacity; vector growth amortizes.
    std::vector<Slice> slices;
    slices.reserve(n_blocks * 2);

    for (std::size_t b = 0; b < n_blocks; ++b) {
        std::int64_t blk_start = block_starts[b];
        std::int64_t blk_end   = block_ends[b];
        if (blk_end <= blk_start) continue;
        std::int64_t s = std::max(blk_start, edge_lo);
        std::int64_t e = std::min(blk_end,   edge_hi);
        if (e <= s) continue;
        std::int64_t cur = s;
        std::int64_t r = region_of_pos(cur);
        while (cur < e && r != -1 &&
               r < static_cast<std::int64_t>(regions_.size()))
        {
            const std::int64_t region_end = boundary_positions_[r + 1];
            const std::int64_t slice_end  = std::min(e, region_end);
            slices.push_back(Slice{r, cur, slice_end});
            cur = slice_end;
            ++r;
        }
    }

    if (slices.empty()) return;

    // 2. Compute L = sum of slice lengths.
    std::int64_t L = 0;
    for (const auto& sl : slices) {
        L += (sl.end - sl.start);
    }
    if (L <= 0) return;
    const double inv_L = 1.0 / static_cast<double>(L);

    // gDNA FL pooling (PR 4c): bin this UNSPLICED fragment's footprint (L) into
    // the FL pool(s) for the region-type(s) it touches — CONTAINED for a
    // single-region fragment, BOUNDARY (fractional per side) for a crossing.
    // Spliced fragments are excluded (their genomic span is not the FL; the RNA
    // FL is the scanner's SPLICED-ANNOT channel). No-op when FL pooling is off.
    const bool fl_on = !spliced && !fl_pool_mass_.empty();
    const std::size_t fl_row = static_cast<std::size_t>(max_fl_) + 1;
    const std::size_t fl_bin =
        fl_on ? static_cast<std::size_t>(std::min<std::int64_t>(L, max_fl_)) : 0;

    // 3. Single-region (all slices in same region) → contained.
    bool all_same = true;
    const std::int64_t r0 = slices.front().region_idx;
    for (std::size_t i = 1; i < slices.size(); ++i) {
        if (slices[i].region_idx != r0) { all_same = false; break; }
    }
    if (all_same) {
        regions_[static_cast<std::size_t>(r0)].contained[ch] += 1u;
        if (fl_on) {
            const std::size_t pool = fl_pool_idx(
                region_types_[static_cast<std::size_t>(r0)], /*boundary=*/false);
            fl_pool_mass_[pool * fl_row + fl_bin] += 1.0;
        }
        return;
    }

    // 4. Crossing path: distribute each slice's mass across the boundaries it
    //    crosses, conserving fragment mass (§4.3 of 00_design.md). A slice
    //    crosses its LEFT boundary iff it is not the first slice, and its RIGHT
    //    boundary iff it is not the last. A region the fragment *encompasses* —
    //    a fully-traversed interior slice overlapping BOTH its boundaries —
    //    splits its mass 50/50: half to the RIGHT side of its left boundary
    //    (mass_right) and half to the LEFT side of its right boundary
    //    (mass_left). End slices keep full mass on their single crossed side.
    //
    //    Flux is PER SIDE and integer (NOT split): the left region's slice
    //    credits flux_left of its right boundary; the right region's slice
    //    credits flux_right of its left boundary. A contiguous crossing credits
    //    both sides of its one boundary; a spliced jump credits one side of each
    //    flanking boundary, leaving the intron-facing sides at zero (no false
    //    exon-intron flux). Slices are monotonic, so each side is credited at
    //    most once per fragment.
    const std::size_t n_sl = slices.size();
    for (std::size_t i = 0; i < n_sl; ++i) {
        const Slice& sl = slices[i];
        const bool crosses_left  = (i > 0);
        const bool crosses_right = (i + 1 < n_sl);
        const int n_cross = (crosses_left ? 1 : 0) + (crosses_right ? 1 : 0);
        if (n_cross == 0) continue;  // defensive; single-region handled above
        const double share = static_cast<double>(sl.end - sl.start) * inv_L
                             / static_cast<double>(n_cross);
        if (crosses_right) {
            Boundary& bo = boundaries_[static_cast<std::size_t>(sl.region_idx + 1)];
            bo.mass_left[ch] += static_cast<float>(share);
            bo.flux_left[ch] += 1u;
        }
        if (crosses_left) {
            Boundary& bi = boundaries_[static_cast<std::size_t>(sl.region_idx)];
            bi.mass_right[ch] += static_cast<float>(share);
            bi.flux_right[ch] += 1u;
        }
    }

    // gDNA FL (crossing): each slice's fractional mass → the BOUNDARY-compartment
    // pool of its region-type, at FL bin = min(footprint, max_fl).
    if (fl_on) {
        for (const auto& sl : slices) {
            const std::size_t pool = fl_pool_idx(
                region_types_[static_cast<std::size_t>(sl.region_idx)], /*boundary=*/true);
            fl_pool_mass_[pool * fl_row + fl_bin] +=
                static_cast<double>(sl.end - sl.start) * inv_L;
        }
    }
}

void Accumulator::merge_from(const Accumulator& other) {
    if (boundary_positions_.size() != other.boundary_positions_.size() ||
        !std::equal(boundary_positions_.begin(), boundary_positions_.end(),
                    other.boundary_positions_.begin()))
    {
        throw std::invalid_argument(
            "Accumulator::merge_from: boundary positions differ");
    }
    const std::size_t n_r = regions_.size();
    const std::size_t n_b = boundaries_.size();
    for (std::size_t i = 0; i < n_r; ++i) {
        for (std::size_t c = 0; c < kNChannels; ++c) {
            regions_[i].contained[c] += other.regions_[i].contained[c];
        }
    }
    for (std::size_t i = 0; i < n_b; ++i) {
        for (std::size_t c = 0; c < kNChannels; ++c) {
            boundaries_[i].mass_left[c]  += other.boundaries_[i].mass_left[c];
            boundaries_[i].mass_right[c] += other.boundaries_[i].mass_right[c];
            boundaries_[i].flux_left[c]  += other.boundaries_[i].flux_left[c];
            boundaries_[i].flux_right[c] += other.boundaries_[i].flux_right[c];
        }
    }
    // gDNA FL pools (PR 4c): element-wise sum (identical partition ⇒ same size).
    if (fl_pool_mass_.size() == other.fl_pool_mass_.size()) {
        for (std::size_t i = 0; i < fl_pool_mass_.size(); ++i) {
            fl_pool_mass_[i] += other.fl_pool_mass_[i];
        }
    }
}

// ============================================================================
// AccumulatorSet
// ============================================================================

AccumulatorSet::AccumulatorSet(const std::int64_t* boundary_positions,
                               std::size_t n_positions,
                               const std::int64_t* ref_pos_offsets,
                               std::size_t n_refs,
                               const std::uint8_t* region_types,
                               std::size_t n_region_types,
                               int max_fl)
{
    const bool fl_on = (region_types != nullptr && max_fl > 0);
    accs_.reserve(n_refs);
    std::size_t region_off = 0;  // cumulative region offset into region_types
    for (std::size_t f = 0; f < n_refs; ++f) {
        const std::int64_t lo = ref_pos_offsets[f];
        const std::int64_t hi = ref_pos_offsets[f + 1];
        if (lo < 0 || hi < lo ||
            static_cast<std::size_t>(hi) > n_positions)
        {
            throw std::invalid_argument(
                "AccumulatorSet: ref_pos_offsets out of range");
        }
        // hi - lo == 0 yields an empty (no-op) Accumulator placeholder for
        // references with no annotated partition. n_regions = positions - 1.
        std::vector<std::int64_t> slice(
            boundary_positions + lo, boundary_positions + hi);
        const std::size_t pos_count = static_cast<std::size_t>(hi - lo);
        const std::size_t n_regions_f = (pos_count >= 1) ? (pos_count - 1) : 0;
        std::vector<std::uint8_t> types_slice;
        if (fl_on && n_regions_f > 0) {
            if (region_off + n_regions_f > n_region_types) {
                throw std::invalid_argument(
                    "AccumulatorSet: region_types shorter than total regions");
            }
            types_slice.assign(region_types + region_off,
                               region_types + region_off + n_regions_f);
        }
        accs_.emplace_back(std::move(slice), std::move(types_slice),
                           fl_on ? max_fl : 0);
        region_off += n_regions_f;
    }
    // Validate one beyond-end offset.
    if (n_refs > 0 &&
        static_cast<std::size_t>(ref_pos_offsets[n_refs]) != n_positions)
    {
        throw std::invalid_argument(
            "AccumulatorSet: ref_pos_offsets[n_refs] must equal n_positions");
    }
    if (fl_on && region_off != n_region_types) {
        throw std::invalid_argument(
            "AccumulatorSet: region_types length != total region count");
    }
    (void)n_positions;  // unused when n_refs == 0
}

Accumulator& AccumulatorSet::at(std::int32_t ref_id) {
    if (ref_id < 0 ||
        static_cast<std::size_t>(ref_id) >= accs_.size())
    {
        throw std::out_of_range("AccumulatorSet::at: ref_id out of range");
    }
    return accs_[static_cast<std::size_t>(ref_id)];
}

const Accumulator& AccumulatorSet::at(std::int32_t ref_id) const {
    if (ref_id < 0 ||
        static_cast<std::size_t>(ref_id) >= accs_.size())
    {
        throw std::out_of_range("AccumulatorSet::at: ref_id out of range");
    }
    return accs_[static_cast<std::size_t>(ref_id)];
}

void AccumulatorSet::merge_from(const AccumulatorSet& other) {
    if (accs_.size() != other.accs_.size()) {
        throw std::invalid_argument(
            "AccumulatorSet::merge_from: n_refs differs");
    }
    for (std::size_t f = 0; f < accs_.size(); ++f) {
        accs_[f].merge_from(other.accs_[f]);
    }
}

}  // namespace rigel::accumulator
