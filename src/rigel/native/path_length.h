// path_length.h — L, the length of a fragment's path, defined once for the resolver (which scores each candidate
// transcript by it in the EM) and the calibration accumulator (which deposits by it). A path is the fragment's genomic
// extent [start, end) with introns cut out: the junctions the reads observed and the introns a hypothesis implies in
// the unsequenced gaps. L is the total length of its segments. Specification: `tests/native/_accumulator_reference.py`.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "constants.h"

namespace rigel {

/// The introns as a sorted, disjoint set clipped to [start, end). Overlapping and abutting introns both merge:
/// overlapping ones are contradictory observations of one molecule (mates that disagree about an acceptor), and
/// abutting ones imply a zero-length exon, which no molecule has. The raw pairs are sorted first and clipped inside the
/// loop, the order the specification uses.
inline void normalise_introns(const IntronBlock* observed,
                              std::size_t n_observed,
                              const IntronBlock* implied,
                              std::size_t n_implied,
                              std::int64_t start,
                              std::int64_t end,
                              std::vector<std::pair<std::int64_t, std::int64_t>>& out)
{
    // Both lists go in before the sort, so the union needs no second scratch vector and their order does not matter.
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
/// This is the one definition of L. The shortcut (end - start) - sum(intron) is wrong on overlapping introns and on an
/// intron that runs past the fragment's end.
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
