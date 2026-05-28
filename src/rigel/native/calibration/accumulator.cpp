/**
 * accumulator.cpp - Out-of-line bodies for fractional CalibrationAccumulator.
 *
 * Per-fragment fractional mass routing:
 *   1. Compute total_aligned_bp across exon blocks.
 *   2. Gate strand-ambiguous fragments to n_excluded_strand_ambig.
 *   3. Compute fl_idx once per fragment (unspliced only).
 *   4. For each block, query region overlaps; for each (block, region)
 *      pair, compute overlap_bp, w = overlap_bp / total_aligned_bp,
 *      classify compartment from fragment-vs-region bounds, and add
 *      w to region_counts and the appropriate channel/signature/FL pool.
 *
 * The fully-spanned case (cross_left && cross_right) splits w in half
 * between BOUNDARY_LEFT and BOUNDARY_RIGHT region_counts channels but
 * adds the total w to the BOUNDARY FL pool exactly once (left/right
 * map to the same pool).
 */

#include "calibration/accumulator.h"

#include <algorithm>
#include <stdexcept>

namespace rigel::calibration {

namespace {

inline int64_t block_region_overlap(int64_t blk_start, int64_t blk_end,
                                    int64_t region_start, int64_t region_end) {
    const int64_t lo = std::max<int64_t>(blk_start, region_start);
    const int64_t hi = std::min<int64_t>(blk_end, region_end);
    return hi > lo ? (hi - lo) : 0;
}

}  // namespace

void CalibrationAccumulator::observe(
    int8_t splice_type,
    int32_t /*ref_id*/,
    int64_t frag_start,
    int64_t frag_end,
    const ExonBlock* exons,
    int32_t n_exons,
    int8_t fragment_strand,
    const RegionIndex& regions)
{
    // Strand gate: ambiguous/unknown strand is not routable into the
    // strand-resolved 12-channel layout. Count and return.
    int strand_idx;
    if (fragment_strand == STRAND_POS) {
        strand_idx = region_signature::kChannelStrandPos;
    } else if (fragment_strand == STRAND_NEG) {
        strand_idx = region_signature::kChannelStrandNeg;
    } else {
        ++payload_.n_excluded_strand_ambig;
        return;
    }

    // Total aligned base pairs (sum of M/D/=/X spans, per block).
    int64_t total_aligned_bp = 0;
    for (int32_t e = 0; e < n_exons; ++e) {
        const ExonBlock& blk = exons[e];
        const int64_t len = static_cast<int64_t>(blk.end) -
                            static_cast<int64_t>(blk.start);
        if (len > 0) total_aligned_bp += len;
    }
    if (total_aligned_bp <= 0) {
        ++payload_.n_unobserved;
        return;
    }

    // From here on the fragment counts as observed.
    ++payload_.n_observed;

    const int splice_idx = (splice_type == SPLICE_UNSPLICED)
                               ? region_signature::kSpliceUnspliced
                               : region_signature::kSpliceSpliced;

    // FL eligibility (unspliced only; strand-collapsed in pools).
    int32_t fl_idx = -1;
    if (splice_idx == region_signature::kSpliceUnspliced) {
        const int64_t fl_raw = frag_end - frag_start;
        if (fl_raw >= 0 && fl_raw < CalibrationPayload::kFlBins) {
            fl_idx = static_cast<int32_t>(fl_raw);
        } else if (fl_raw >= CalibrationPayload::kFlBins) {
            // Clip to last bin so very long fragments still contribute
            // mass to the FL pool (consistent with the legacy cap).
            fl_idx = CalibrationPayload::kFlBins - 1;
        } else {
            // fl_raw < 0 should not happen for well-formed unspliced
            // fragments, but count it for audit and skip FL contribution.
            ++payload_.n_fl_unavailable;
        }
    }

    const double inv_total =
        1.0 / static_cast<double>(total_aligned_bp);

    bool any_region_hit = false;
    touched_regions_.clear();

    for (int32_t e = 0; e < n_exons; ++e) {
        const ExonBlock& blk = exons[e];

        block_hits_.clear();
        regions.overlap_into(blk.ref_id,
                             static_cast<int64_t>(blk.start),
                             static_cast<int64_t>(blk.end),
                             block_hits_);
        if (block_hits_.empty()) continue;
        any_region_hit = true;

        for (int32_t rid : block_hits_) {
            const int64_t rs = regions.start(rid);
            const int64_t re = regions.end(rid);
            const int64_t overlap_bp = block_region_overlap(
                static_cast<int64_t>(blk.start),
                static_cast<int64_t>(blk.end),
                rs, re);
            if (overlap_bp <= 0) continue;

            // Record this region as touched by the current fragment.
            // Duplicates (multi-block fragments hitting the same region)
            // are reduced to a single support increment below.
            touched_regions_.push_back(rid);

            const double w = static_cast<double>(overlap_bp) * inv_total;

            // Fragment-vs-region boundary classification (per fragment,
            // not per block; spliced fragments use end-to-end extent).
            const bool cross_left  = frag_start < rs;
            const bool cross_right = frag_end   > re;

            const std::uint8_t sig = regions.signature(rid);
            const size_t reg_base =
                static_cast<size_t>(rid) * CalibrationPayload::kChannels;

            auto add_channel = [&](int compartment, double mass) {
                const int ch = region_signature::channel_index(
                    compartment, splice_idx, strand_idx);
                payload_.region_counts[reg_base + static_cast<size_t>(ch)] +=
                    static_cast<float>(mass);
                payload_.channel_mass[static_cast<size_t>(ch)] += mass;
                payload_.signature_mass[static_cast<size_t>(sig)] += mass;

                if (fl_idx >= 0) {
                    const int pool = region_signature::fl_pool_index(
                        sig, compartment);
                    payload_.fl_pool_mass[
                        static_cast<size_t>(pool) *
                            CalibrationPayload::kFlBins +
                        static_cast<size_t>(fl_idx)] += mass;
                    payload_.fl_pool_total[static_cast<size_t>(pool)] += mass;
                }
            };

            if (cross_left && cross_right) {
                // Fully spans the region: split region_counts mass
                // 50/50 between BOUNDARY_LEFT and BOUNDARY_RIGHT.
                // FL pool: left and right boundary slots collapse to
                // the same pool, so add total w exactly once via the
                // left compartment to avoid double-counting.
                const double half = 0.5 * w;
                const int ch_l = region_signature::channel_index(
                    region_signature::kCompartmentBoundaryLeft,
                    splice_idx, strand_idx);
                const int ch_r = region_signature::channel_index(
                    region_signature::kCompartmentBoundaryRight,
                    splice_idx, strand_idx);
                payload_.region_counts[reg_base + static_cast<size_t>(ch_l)] +=
                    static_cast<float>(half);
                payload_.region_counts[reg_base + static_cast<size_t>(ch_r)] +=
                    static_cast<float>(half);
                payload_.channel_mass[static_cast<size_t>(ch_l)] += half;
                payload_.channel_mass[static_cast<size_t>(ch_r)] += half;
                payload_.signature_mass[static_cast<size_t>(sig)] += w;

                if (fl_idx >= 0) {
                    const int pool = region_signature::fl_pool_index(
                        sig, region_signature::kCompartmentBoundaryLeft);
                    payload_.fl_pool_mass[
                        static_cast<size_t>(pool) *
                            CalibrationPayload::kFlBins +
                        static_cast<size_t>(fl_idx)] += w;
                    payload_.fl_pool_total[static_cast<size_t>(pool)] += w;
                }
            } else if (cross_left) {
                add_channel(region_signature::kCompartmentBoundaryLeft, w);
            } else if (cross_right) {
                add_channel(region_signature::kCompartmentBoundaryRight, w);
            } else {
                add_channel(region_signature::kCompartmentContained, w);
            }
        }
    }

    if (!any_region_hit) {
        ++payload_.n_unannotated_ref;
    }

    // Per-region support: increment exactly once per unique region the
    // current fragment touched. Compartment-agnostic (contained, left,
    // right, and fully-spanning all contribute one increment per region).
    if (!touched_regions_.empty()) {
        std::sort(touched_regions_.begin(), touched_regions_.end());
        auto last = std::unique(touched_regions_.begin(),
                                touched_regions_.end());
        auto& support = (splice_idx == region_signature::kSpliceUnspliced)
                            ? payload_.region_unspliced_support
                            : payload_.region_spliced_support;
        for (auto it = touched_regions_.begin(); it != last; ++it) {
            ++support[static_cast<size_t>(*it)];
        }
    }
}

namespace {

template <typename T>
inline void add_into(std::vector<T>& dst, const std::vector<T>& src) {
    const size_t n = dst.size();
    for (size_t i = 0; i < n; ++i) dst[i] += src[i];
}

template <typename T, std::size_t N>
inline void add_into(std::array<T, N>& dst, const std::array<T, N>& src) {
    for (std::size_t i = 0; i < N; ++i) dst[i] += src[i];
}

}  // namespace

void CalibrationAccumulator::merge_from(const CalibrationAccumulator& other) {
    if (other.n_regions_ != n_regions_) {
        throw std::runtime_error(
            "CalibrationAccumulator::merge_from: n_regions mismatch");
    }
    add_into(payload_.region_counts,    other.payload_.region_counts);
    add_into(payload_.channel_mass,     other.payload_.channel_mass);
    add_into(payload_.signature_mass,   other.payload_.signature_mass);
    add_into(payload_.fl_pool_mass,     other.payload_.fl_pool_mass);
    add_into(payload_.fl_pool_total,    other.payload_.fl_pool_total);
    add_into(payload_.region_unspliced_support,
             other.payload_.region_unspliced_support);
    add_into(payload_.region_spliced_support,
             other.payload_.region_spliced_support);
    payload_.n_observed              += other.payload_.n_observed;
    payload_.n_excluded_multimap     += other.payload_.n_excluded_multimap;
    payload_.n_excluded_chimera      += other.payload_.n_excluded_chimera;
    payload_.n_excluded_artifact     += other.payload_.n_excluded_artifact;
    payload_.n_excluded_strand_ambig += other.payload_.n_excluded_strand_ambig;
    payload_.n_unobserved            += other.payload_.n_unobserved;
    payload_.n_unannotated_ref       += other.payload_.n_unannotated_ref;
    payload_.n_fl_unavailable        += other.payload_.n_fl_unavailable;
}

}  // namespace rigel::calibration
