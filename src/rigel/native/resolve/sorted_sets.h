// sorted_sets.h — transcript sets as sorted, unique int32 vectors, and how one fragment's per-block sets merge into its
// candidates.
#pragma once

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <vector>

namespace rigel {

inline void sort_unique(std::vector<int32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

inline std::vector<int32_t> sorted_intersection(const std::vector<int32_t>& a, const std::vector<int32_t>& b) {
    std::vector<int32_t> out;
    out.reserve(std::min(a.size(), b.size()));
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));
    return out;
}

inline std::vector<int32_t> sorted_union(const std::vector<int32_t>& a, const std::vector<int32_t>& b) {
    std::vector<int32_t> out;
    out.reserve(a.size() + b.size());
    std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));
    return out;
}

/// One fragment's candidates from its per-block sets: the intersection of the non-empty sets, or their union when that
/// is empty. An empty set is a block that overlaps no transcript, so it constrains nothing; no non-empty set gives an
/// empty result.
inline std::vector<int32_t> merge_sets(const std::vector<std::vector<int32_t>>& sets) {
    std::vector<const std::vector<int32_t>*> non_empty;
    for (const auto& s : sets)
        if (!s.empty()) non_empty.push_back(&s);
    if (non_empty.empty()) return {};

    std::vector<int32_t> isect = *non_empty[0];
    for (size_t i = 1; i < non_empty.size() && !isect.empty(); i++) isect = sorted_intersection(isect, *non_empty[i]);
    if (!isect.empty()) return isect;

    std::vector<int32_t> all = *non_empty[0];
    for (size_t i = 1; i < non_empty.size(); i++) all = sorted_union(all, *non_empty[i]);
    return all;
}

}  // namespace rigel
