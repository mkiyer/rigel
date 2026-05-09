/**
 * accumulator.cpp — out-of-line bodies for CalibrationAccumulator.
 *
 * The dedup strategy relies on the M2 invariant: each ``overlap_into``
 * call returns strictly increasing region IDs.  Across exon blocks
 * the per-block lists may overlap (a spliced fragment whose two
 * exons both touch the same boundary exon region), so we sorted-merge
 * each new block list into a persistent fragment-scope ``hits_``
 * vector with a two-pointer pass that skips duplicates.
 */

#include "calibration/accumulator.h"

#include <algorithm>
#include <stdexcept>

namespace rigel::calibration {

namespace {

// Compute exact aligned-block overlap of [blk_start, blk_end) with the
// region [region_start, region_end). Returns >= 1 because RegionIndex
// only returns hits with positive overlap on at least one block.
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
    const RegionIndex& regions)
{
    hits_.clear();
    hits_overlap_bp_.clear();
    qualified_hits_.clear();

    for (int32_t e = 0; e < n_exons; ++e) {
        const ExonBlock& blk = exons[e];

        block_hits_.clear();
        regions.overlap_into(blk.ref_id,
                             static_cast<int64_t>(blk.start),
                             static_cast<int64_t>(blk.end),
                             block_hits_);
        if (block_hits_.empty()) continue;

        // Compute exact per-region overlap_bp for this block.
        block_overlap_bp_.clear();
        block_overlap_bp_.reserve(block_hits_.size());
        for (int32_t rid : block_hits_) {
            block_overlap_bp_.push_back(block_region_overlap(
                static_cast<int64_t>(blk.start),
                static_cast<int64_t>(blk.end),
                regions.start(rid),
                regions.end(rid)));
        }

        if (hits_.empty()) {
            // Fast path: first non-empty block — adopt its hits wholesale.
            hits_.swap(block_hits_);
            hits_overlap_bp_.swap(block_overlap_bp_);
            continue;
        }

        // Sorted-merge block_hits_ into hits_, summing per-region
        // overlap_bp where the same region appears in both lists
        // (e.g. a spliced fragment whose two exons both touch the
        // same boundary region).
        merge_buf_.clear();
        merge_buf_bp_.clear();
        size_t i = 0;
        size_t j = 0;
        const size_t na = hits_.size();
        const size_t nb = block_hits_.size();
        while (i < na && j < nb) {
            const int32_t a = hits_[i];
            const int32_t b = block_hits_[j];
            if (a < b) {
                merge_buf_.push_back(a);
                merge_buf_bp_.push_back(hits_overlap_bp_[i]);
                ++i;
            } else if (a > b) {
                merge_buf_.push_back(b);
                merge_buf_bp_.push_back(block_overlap_bp_[j]);
                ++j;
            } else {
                merge_buf_.push_back(a);
                merge_buf_bp_.push_back(hits_overlap_bp_[i] + block_overlap_bp_[j]);
                ++i;
                ++j;
            }
        }
        while (i < na) {
            merge_buf_.push_back(hits_[i]);
            merge_buf_bp_.push_back(hits_overlap_bp_[i]);
            ++i;
        }
        while (j < nb) {
            merge_buf_.push_back(block_hits_[j]);
            merge_buf_bp_.push_back(block_overlap_bp_[j]);
            ++j;
        }
        hits_.swap(merge_buf_);
        hits_overlap_bp_.swap(merge_buf_bp_);
    }

    // Boundary-tolerance qualification: a hit contributes its type bit
    // (and gets per-region / boundary-flux fan-out) only if its exact
    // aligned-block overlap is >= q = max(K, 1). At K=0 this is a
    // tautology because RegionIndex::overlap_into already only returns
    // hits with overlap >= 1, so qualified_hits_ == hits_ and behavior
    // is bit-identical to the pre-tolerance code path.
    const int64_t q = static_cast<int64_t>(std::max(splicing_anchor_tolerance_, 1));
    uint8_t obs_mask = 0;
    for (size_t k = 0; k < hits_.size(); ++k) {
        if (hits_overlap_bp_[k] >= q) {
            qualified_hits_.push_back(hits_[k]);
            obs_mask |= regions.type_mask(hits_[k]);
        }
    }

    // Always: bump global mask counter, FL histogram, n_observed.
    payload_.global_counts[obs_mask]++;
    const int64_t fl_raw = frag_end - frag_start;
    int64_t fl_idx = fl_raw < 0 ? 0 : fl_raw;
    if (fl_idx >= CalibrationPayload::kFlBins) {
        fl_idx = CalibrationPayload::kFlBins - 1;
    }
    payload_.fl_hist[static_cast<size_t>(obs_mask) * CalibrationPayload::kFlBins
                     + static_cast<size_t>(fl_idx)]++;
    payload_.n_observed++;

    if (obs_mask == 0) {
        // Split: raw hits existed but all sub-tolerance => below-tolerance;
        // otherwise (no raw hits at all, or all hits had type_mask 0) =>
        // unannotated_ref. At K=0 (q=1) every hit is qualified, so the
        // first branch can never fire and the legacy semantics hold.
        if (!hits_.empty() && qualified_hits_.empty()) {
            payload_.n_below_tolerance++;
        } else {
            payload_.n_unannotated_ref++;
        }
        return;
    }

    // Per-region count fan-out and (for unspliced single-ref fragments)
    // boundary-flux fan-out.  Only qualified hits participate.  The
    // boundary-flux endpoint predicate uses ``frag_start + q <= b &&
    // frag_end >= b + q`` which collapses to the strict ``< b && > b``
    // crossing predicate at K=0 (q=1) on integer half-open coordinates.
    const bool flux_eligible = (splice_type == SPLICE_UNSPLICED);
    for (int32_t rid : qualified_hits_) {
        payload_.per_region_counts[
            static_cast<size_t>(rid) * mask::N_STATES
            + static_cast<size_t>(obs_mask)]++;
        if (flux_eligible && (regions.type_mask(rid) & mask::EXON) != 0) {
            const int64_t rs = regions.start(rid);
            const int64_t re = regions.end(rid);
            if (frag_start + q <= rs && frag_end >= rs + q) payload_.u_left[rid]++;
            if (frag_start + q <= re && frag_end >= re + q) payload_.u_right[rid]++;
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
    add_into(payload_.global_counts,     other.payload_.global_counts);
    add_into(payload_.per_region_counts, other.payload_.per_region_counts);
    add_into(payload_.fl_hist,           other.payload_.fl_hist);
    add_into(payload_.u_left,            other.payload_.u_left);
    add_into(payload_.u_right,           other.payload_.u_right);
    payload_.n_observed          += other.payload_.n_observed;
    payload_.n_excluded_multimap += other.payload_.n_excluded_multimap;
    payload_.n_excluded_chimera  += other.payload_.n_excluded_chimera;
    payload_.n_excluded_artifact += other.payload_.n_excluded_artifact;
    payload_.n_unobserved        += other.payload_.n_unobserved;
    payload_.n_unannotated_ref   += other.payload_.n_unannotated_ref;
    payload_.n_below_tolerance   += other.payload_.n_below_tolerance;
}

}  // namespace rigel::calibration
