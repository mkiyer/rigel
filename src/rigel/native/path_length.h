// path_length.h — L, the length of a fragment's path, defined once for the resolver, which scores each candidate
// transcript by it in the EM, and the accumulator, which deposits by it in calibration (`calibration/accumulator.cpp`).
// A path is the fragment's genomic extent [start, end) with introns cut out: the junctions the reads observed and the
// introns a hypothesis implies in the unsequenced gaps. L is the total length of its segments. The specification is
// `tests/native/_accumulator_reference.py`; `tests/test_transcript_space_fl.py` holds the resolver's length to the
// accumulator's for every candidate.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "constants.h"

namespace rigel {

/// The introns as a sorted, DISJOINT set clipped to [start, end).
///
/// ⚠ Overlapping AND ABUTTING introns both merge. Overlapping ones are contradictory observations of one
/// molecule -- a real BAM produces them when the mates disagree about an acceptor. Abutting ones imply a
/// zero-length exon, which is physically impossible, so no single molecule can legitimately use both.
///
/// ⚠ Sort the RAW pairs and clip inside the loop, in that order, to match the spec exactly.
inline void normalise_introns(const IntronBlock* observed,
                              std::size_t n_observed,
                              const IntronBlock* implied,
                              std::size_t n_implied,
                              std::int64_t start,
                              std::int64_t end,
                              std::vector<std::pair<std::int64_t, std::int64_t>>& out)
{
    // ⚠ Both ranges are emplaced BEFORE the sort, so the union needs no second scratch vector and the
    // caller never has to think about which list came first.
    out.clear();
    out.reserve(n_observed + n_implied);
    for (std::size_t i = 0; i < n_observed; ++i) {
        out.emplace_back(static_cast<std::int64_t>(observed[i].start),
                         static_cast<std::int64_t>(observed[i].end));
    }
    for (std::size_t i = 0; i < n_implied; ++i) {
        out.emplace_back(static_cast<std::int64_t>(implied[i].start),
                         static_cast<std::int64_t>(implied[i].end));
    }
    std::sort(out.begin(), out.end());

    std::size_t kept = 0;
    for (std::size_t i = 0; i < out.size(); ++i) {
        const std::int64_t intron_start = std::max(out[i].first, start);
        const std::int64_t intron_end   = std::min(out[i].second, end);
        if (intron_end <= intron_start) continue;  // zero-length, or entirely outside the fragment
        if (kept > 0 && intron_start <= out[kept - 1].second) {
            out[kept - 1].second = std::max(out[kept - 1].second, intron_end);
            continue;
        }
        out[kept++] = {intron_start, intron_end};  // kept <= i, so this never clobbers an unread slot
    }
    out.resize(kept);
}

/// The path's contiguous genomic segments: [start, end) with the (normalised) introns cut out.
inline void build_segments(std::int64_t start,
                           std::int64_t end,
                           const std::vector<std::pair<std::int64_t, std::int64_t>>& introns,
                           std::vector<std::pair<std::int64_t, std::int64_t>>& out)
{
    out.clear();
    std::int64_t cursor = start;
    for (const auto& [intron_start, intron_end] : introns) {
        if (intron_start > cursor) out.emplace_back(cursor, intron_start);
        cursor = intron_end;
    }
    if (end > cursor) out.emplace_back(cursor, end);
}

/// L: the total length of the path's segments, [start, end) less the observed and the implied introns. ``introns``
/// and ``segments`` are the caller's scratch and hold the normalised introns and the segments on return, which the
/// deposit walks.
///
/// ⚠ ONE definition of L. Deriving it any other way invites two formulas for one quantity, and the obvious second
/// formula -- (end - start) - sum(intron) -- disagrees by up to 1.5x on overlapping introns, goes NEGATIVE on a wide
/// overlap, and miscounts an intron that runs past the fragment's end.
inline std::int64_t path_length(const IntronBlock* observed,
                                std::size_t n_observed,
                                const IntronBlock* implied,
                                std::size_t n_implied,
                                std::int64_t start,
                                std::int64_t end,
                                std::vector<std::pair<std::int64_t, std::int64_t>>& introns,
                                std::vector<std::pair<std::int64_t, std::int64_t>>& segments)
{
    normalise_introns(observed, n_observed, implied, n_implied, start, end, introns);
    build_segments(start, end, introns, segments);
    std::int64_t length = 0;
    for (const auto& [a, b] : segments) length += b - a;
    return length;
}

}  // namespace rigel
