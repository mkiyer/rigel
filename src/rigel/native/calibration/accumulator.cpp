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

Accumulator::Accumulator(std::vector<std::int64_t> boundaries)
    : boundary_positions_(std::move(boundaries))
{
    // An empty position array is permitted and yields a no-op accumulator
    // (n_regions == 0, n_boundaries == 0). This lets AccumulatorSet allocate
    // a placeholder for references that have no annotated partition.
    if (boundary_positions_.empty()) {
        regions_.clear();
        boundaries_.clear();
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
                          bool strand_pos)
{
    if (n_blocks == 0 || regions_.empty()) return;

    const int ch = channel_idx(spliced, strand_pos);
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

    // 3. Single-region (all slices in same region) → contained.
    bool all_same = true;
    const std::int64_t r0 = slices.front().region_idx;
    for (std::size_t i = 1; i < slices.size(); ++i) {
        if (slices[i].region_idx != r0) { all_same = false; break; }
    }
    if (all_same) {
        regions_[static_cast<std::size_t>(r0)].contained[ch] += 1u;
        return;
    }

    // 4. Crossing path: for each consecutive slice pair, deposit on the
    //    boundary between r_left's right-boundary and r_right's left-
    //    boundary (these are the same boundary iff adjacent; different
    //    iff non-adjacent).
    //
    //    Collect touched boundary indices into a small scratch vector
    //    and dedupe before incrementing flux.
    std::vector<std::int64_t> touched;
    touched.reserve(slices.size() * 2);

    for (std::size_t i = 0; i + 1 < slices.size(); ++i) {
        const Slice& a = slices[i];
        const Slice& b = slices[i + 1];
        const std::int64_t len_left  = a.end - a.start;
        const std::int64_t len_right = b.end - b.start;
        const std::int64_t b_out = a.region_idx + 1;  // right boundary of r_left
        const std::int64_t b_in  = b.region_idx;      // left  boundary of r_right

        boundaries_[static_cast<std::size_t>(b_out)].mass_left[ch] +=
            static_cast<float>(static_cast<double>(len_left) * inv_L);
        boundaries_[static_cast<std::size_t>(b_in)].mass_right[ch] +=
            static_cast<float>(static_cast<double>(len_right) * inv_L);

        touched.push_back(b_out);
        touched.push_back(b_in);
    }

    // Dedupe touched and increment flux.
    std::sort(touched.begin(), touched.end());
    touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
    for (std::int64_t b : touched) {
        boundaries_[static_cast<std::size_t>(b)].flux[ch] += 1u;
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
            boundaries_[i].flux[c]       += other.boundaries_[i].flux[c];
        }
    }
}

// ============================================================================
// AccumulatorSet
// ============================================================================

AccumulatorSet::AccumulatorSet(const std::int64_t* boundary_positions,
                               std::size_t n_positions,
                               const std::int64_t* ref_pos_offsets,
                               std::size_t n_refs)
{
    accs_.reserve(n_refs);
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
        // references with no annotated partition.
        std::vector<std::int64_t> slice(
            boundary_positions + lo, boundary_positions + hi);
        accs_.emplace_back(std::move(slice));
    }
    // Validate one beyond-end offset.
    if (n_refs > 0 &&
        static_cast<std::size_t>(ref_pos_offsets[n_refs]) != n_positions)
    {
        throw std::invalid_argument(
            "AccumulatorSet: ref_pos_offsets[n_refs] must equal n_positions");
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
