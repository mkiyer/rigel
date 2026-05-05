/**
 * accumulator.cpp — out-of-line bodies for CalibrationAccumulator.
 */

#include "calibration/accumulator.h"

#include <algorithm>
#include <stdexcept>

namespace rigel::calibration {

namespace {
constexpr uint8_t kExonBit = 0b001;
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
    spill_buf_.clear();
    int32_t n_inline = 0;
    uint8_t mask = 0;

    // Pass A: walk every exon block, query overlap into shared
    // inline+spill buffers, dedup against what is already there.
    for (int32_t e = 0; e < n_exons; ++e) {
        const ExonBlock& blk = exons[e];
        const int32_t pre_inline = n_inline;
        const std::size_t pre_spill = spill_buf_.size();

        const int32_t k = regions.overlap_into(
            blk.ref_id,
            static_cast<int64_t>(blk.start),
            static_cast<int64_t>(blk.end),
            inline_buf_, spill_buf_, n_inline);
        if (k == 0) continue;

        // overlap_into appends to inline first (greedily up to capacity),
        // then spills the remainder.  Compute the split.
        const int32_t cap_remaining =
            static_cast<int32_t>(kInlineCap) - pre_inline;
        const int32_t k_inline =
            std::min<int32_t>(k, std::max<int32_t>(0, cap_remaining));

        // Dedup new inline entries [pre_inline, pre_inline + k_inline)
        // against committed prefix [0, pre_inline) and pre-existing
        // spill (size pre_spill). Compact in place; update mask.
        int32_t write = pre_inline;
        for (int32_t j = 0; j < k_inline; ++j) {
            const int32_t rid = inline_buf_[pre_inline + j];
            bool dup = false;
            for (int32_t q = 0; q < write; ++q) {
                if (inline_buf_[q] == rid) { dup = true; break; }
            }
            if (!dup) {
                for (std::size_t q = 0; q < pre_spill; ++q) {
                    if (spill_buf_[q] == rid) { dup = true; break; }
                }
            }
            if (dup) continue;
            inline_buf_[write++] = rid;
            mask |= regions.type_mask(rid);
        }
        n_inline = write;

        // Dedup new spill entries against [0, n_inline) and
        // pre-existing spill prefix.
        std::size_t spill_write = pre_spill;
        for (std::size_t j = pre_spill; j < spill_buf_.size(); ++j) {
            const int32_t rid = spill_buf_[j];
            bool dup = false;
            for (int32_t q = 0; q < n_inline; ++q) {
                if (inline_buf_[q] == rid) { dup = true; break; }
            }
            if (!dup) {
                for (std::size_t q = 0; q < spill_write; ++q) {
                    if (spill_buf_[q] == rid) { dup = true; break; }
                }
            }
            if (dup) continue;
            spill_buf_[spill_write++] = rid;
            mask |= regions.type_mask(rid);
        }
        spill_buf_.resize(spill_write);
    }

    // Always: bump global mask counter, FL histogram, n_observed.
    payload_.global_counts[mask]++;
    const int64_t fl_raw = frag_end - frag_start;
    int64_t fl_idx = fl_raw < 0 ? 0 : fl_raw;
    if (fl_idx >= CalibrationPayload::kFlBins) {
        fl_idx = CalibrationPayload::kFlBins - 1;
    }
    payload_.fl_hist[static_cast<size_t>(mask) * CalibrationPayload::kFlBins
                     + static_cast<size_t>(fl_idx)]++;
    payload_.n_observed++;

    if (mask == 0) {
        payload_.n_oor++;
        return;
    }

    auto bump_region = [&](int32_t rid) {
        payload_.per_region_counts[
            static_cast<size_t>(rid) * CalibrationPayload::kMaskCard
            + static_cast<size_t>(mask)]++;
    };
    for (int32_t i = 0; i < n_inline; ++i) bump_region(inline_buf_[i]);
    for (int32_t rid : spill_buf_)         bump_region(rid);

    // Pass D: boundary-flux for unspliced single-ref fragments.
    if (splice_type == SPLICE_UNSPLICED) {
        auto check_one = [&](int32_t rid) {
            if ((regions.type_mask(rid) & kExonBit) == 0) return;
            const int64_t rs = regions.start(rid);
            const int64_t re = regions.end(rid);
            if (frag_start < rs && frag_end > rs) payload_.u_left[rid]++;
            if (frag_start < re && frag_end > re) payload_.u_right[rid]++;
        };
        for (int32_t i = 0; i < n_inline; ++i) check_one(inline_buf_[i]);
        for (int32_t rid : spill_buf_)         check_one(rid);
    }
}

void CalibrationAccumulator::merge_from(const CalibrationAccumulator& other) {
    if (other.n_regions_ != n_regions_) {
        throw std::runtime_error(
            "CalibrationAccumulator::merge_from: n_regions mismatch");
    }
    for (size_t m = 0; m < CalibrationPayload::kMaskCard; ++m) {
        payload_.global_counts[m] += other.payload_.global_counts[m];
    }
    const auto& src_prc = other.payload_.per_region_counts;
    for (size_t i = 0; i < src_prc.size(); ++i) {
        payload_.per_region_counts[i] += src_prc[i];
    }
    const auto& src_fl = other.payload_.fl_hist;
    for (size_t i = 0; i < src_fl.size(); ++i) {
        payload_.fl_hist[i] += src_fl[i];
    }
    const auto& src_ul = other.payload_.u_left;
    for (size_t i = 0; i < src_ul.size(); ++i) {
        payload_.u_left[i] += src_ul[i];
    }
    const auto& src_ur = other.payload_.u_right;
    for (size_t i = 0; i < src_ur.size(); ++i) {
        payload_.u_right[i] += src_ur[i];
    }
    payload_.n_observed          += other.payload_.n_observed;
    payload_.n_excluded_multimap += other.payload_.n_excluded_multimap;
    payload_.n_excluded_chimera  += other.payload_.n_excluded_chimera;
    payload_.n_excluded_artifact += other.payload_.n_excluded_artifact;
    payload_.n_unobserved        += other.payload_.n_unobserved;
    payload_.n_oor               += other.payload_.n_oor;
}

}  // namespace rigel::calibration
