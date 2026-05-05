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

#include <stdexcept>

namespace rigel::calibration {

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
    uint8_t obs_mask = 0;

    for (int32_t e = 0; e < n_exons; ++e) {
        const ExonBlock& blk = exons[e];

        block_hits_.clear();
        regions.overlap_into(blk.ref_id,
                             static_cast<int64_t>(blk.start),
                             static_cast<int64_t>(blk.end),
                             block_hits_);
        if (block_hits_.empty()) continue;

        if (hits_.empty()) {
            // Fast path: first non-empty block — adopt its hits wholesale.
            for (int32_t rid : block_hits_) obs_mask |= regions.type_mask(rid);
            hits_.swap(block_hits_);
            continue;
        }

        // Sorted-merge block_hits_ into hits_, skipping duplicates and
        // OR-ing type masks for each newly added region.
        merge_buf_.clear();
        size_t i = 0;
        size_t j = 0;
        const size_t na = hits_.size();
        const size_t nb = block_hits_.size();
        while (i < na && j < nb) {
            const int32_t a = hits_[i];
            const int32_t b = block_hits_[j];
            if (a < b) {
                merge_buf_.push_back(a);
                ++i;
            } else if (a > b) {
                obs_mask |= regions.type_mask(b);
                merge_buf_.push_back(b);
                ++j;
            } else {
                merge_buf_.push_back(a);
                ++i;
                ++j;
            }
        }
        while (i < na) merge_buf_.push_back(hits_[i++]);
        while (j < nb) {
            const int32_t b = block_hits_[j++];
            obs_mask |= regions.type_mask(b);
            merge_buf_.push_back(b);
        }
        hits_.swap(merge_buf_);
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
        payload_.n_unannotated_ref++;
        return;
    }

    // Per-region count fan-out and (for unspliced single-ref fragments)
    // boundary-flux fan-out.  Single linear pass over the deduped hits.
    const bool flux_eligible = (splice_type == SPLICE_UNSPLICED);
    for (int32_t rid : hits_) {
        payload_.per_region_counts[
            static_cast<size_t>(rid) * mask::N_STATES
            + static_cast<size_t>(obs_mask)]++;
        if (flux_eligible && (regions.type_mask(rid) & mask::EXON) != 0) {
            const int64_t rs = regions.start(rid);
            const int64_t re = regions.end(rid);
            if (frag_start < rs && frag_end > rs) payload_.u_left[rid]++;
            if (frag_start < re && frag_end > re) payload_.u_right[rid]++;
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
}

}  // namespace rigel::calibration
