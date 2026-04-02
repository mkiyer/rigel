/**
 * constants.h — Shared constants and helper types for rigel C++ extensions.
 *
 * Must match the Python IntEnum values exactly (rigel.types,
 * rigel.categories).  Included by resolve.cpp, bam_scanner.cpp, and
 * any future native modules.
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace rigel {

// ================================================================
// Enum-mirror constants
// ================================================================

// IntervalType (rigel.types.IntervalType)
static constexpr int8_t ITYPE_EXON           = 0;
static constexpr int8_t ITYPE_TRANSCRIPT      = 1;

// SpliceType (rigel.categories.SpliceType)
static constexpr int32_t SPLICE_UNSPLICED       = 0;
static constexpr int32_t SPLICE_SPLICED_UNANNOT = 1;
static constexpr int32_t SPLICE_SPLICED_ANNOT   = 2;

// MergeOutcome (rigel.types.MergeOutcome)
static constexpr int32_t MC_INTERSECTION          = 0;
static constexpr int32_t MC_INTERSECTION_NONEMPTY = 1;
static constexpr int32_t MC_UNION                 = 2;
static constexpr int32_t MC_EMPTY                 = 3;

// ChimeraType (rigel.types.ChimeraType)
static constexpr int32_t CHIMERA_NONE           = 0;
static constexpr int32_t CHIMERA_TRANS           = 1;
static constexpr int32_t CHIMERA_CIS_STRAND_SAME = 2;
static constexpr int32_t CHIMERA_CIS_STRAND_DIFF = 3;

// Strand (rigel.types.Strand)
static constexpr int32_t STRAND_NONE      = 0;
static constexpr int32_t STRAND_POS       = 1;
static constexpr int32_t STRAND_NEG       = 2;
static constexpr int32_t STRAND_AMBIGUOUS = 3;  // POS | NEG

// FragmentClass (rigel.buffer — fragment classification codes)
static constexpr int32_t FRAG_UNAMBIG           = 0;
static constexpr int32_t FRAG_AMBIG_SAME_STRAND = 1;
static constexpr int32_t FRAG_AMBIG_OPP_STRAND  = 2;
static constexpr int32_t FRAG_MULTIMAPPER       = 3;
static constexpr int32_t FRAG_CHIMERIC          = 4;

// ZF assignment flags bitfield (rigel.annotate — BAM ZF:i tag)
static constexpr int32_t ZF_RESOLVED  = 0x1;   // bit 0: scored and assigned
static constexpr int32_t ZF_GDNA      = 0x2;   // bit 1: gDNA component won
static constexpr int32_t ZF_NRNA      = 0x4;   // bit 2: single-exon transcript
static constexpr int32_t ZF_SYNTHETIC = 0x8;   // bit 3: rigel-generated nRNA span

// ================================================================
// Scoring constants (shared by scoring.cpp and Python side)
// ================================================================
static constexpr double LOG_HALF       = -0.6931471805599453;       // log(0.5)
static constexpr double TAIL_DECAY_LP  = -0.01005033585350145;      // log(0.99)

// ================================================================
// Geometric helper types shared across modules
// ================================================================

struct ExonBlock {
    int32_t ref_id;
    int32_t start;
    int32_t end;
    int32_t strand;
};

struct IntronBlock {
    int32_t ref_id;
    int32_t start;
    int32_t end;
    int32_t strand;
};

// Merge result: sorted vector of transcript indices + criteria
struct MergeResult {
    std::vector<int32_t> t_inds;
    int32_t criteria;
    bool is_empty() const { return t_inds.empty(); }
};

// Chimera detection result
struct ChimeraResult {
    int32_t type;
    int32_t gap;
};

// SJ exact-match map key
struct SJKey {
    int32_t ref_id;
    int32_t start;
    int32_t end;
    int32_t strand;

    bool operator==(const SJKey& o) const {
        return ref_id == o.ref_id && start == o.start &&
               end == o.end && strand == o.strand;
    }
};

struct SJKeyHash {
    size_t operator()(const SJKey& k) const {
        size_t h = std::hash<int32_t>()(k.ref_id);
        h ^= std::hash<int32_t>()(k.start)  + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>()(k.end)    + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>()(k.strand) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

// ================================================================
// RawResolveResult — internal result from _resolve_core()
// ================================================================

struct RawResolveResult {
    std::vector<int32_t> t_inds;
    int32_t ambig_strand = 0;
    int32_t splice_type = 0;
    int32_t exon_strand = 0;
    int32_t sj_strand = 0;
    int32_t genomic_footprint = -1;
    int32_t genomic_start = -1;
    int32_t merge_criteria = MC_EMPTY;
    int32_t read_length = 0;
    int32_t chimera_type = CHIMERA_NONE;
    int32_t chimera_gap = -1;
    std::unordered_map<int32_t, int32_t> frag_length_map;
    // Parallel arrays to t_inds
    std::vector<int32_t> t_exon_bp;
    std::vector<int32_t> t_intron_bp;
};

// ================================================================
// Set operations
// ================================================================

// Check whether two sorted vectors share at least one element.
inline bool has_intersection(const std::vector<int32_t>& a,
                             const std::vector<int32_t>& b) {
    auto ai = a.begin(), bi = b.begin();
    while (ai != a.end() && bi != b.end()) {
        if (*ai < *bi) ++ai;
        else if (*ai > *bi) ++bi;
        else return true;
    }
    return false;
}

// Sorted-vector intersection.
inline std::vector<int32_t> vec_intersect(const std::vector<int32_t>& a,
                                          const std::vector<int32_t>& b) {
    std::vector<int32_t> out;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                          std::back_inserter(out));
    return out;
}

// Progressive set merging identical to Python merge_sets_with_criteria().
inline MergeResult merge_sets(const std::vector<std::vector<int32_t>>& sets) {
    if (sets.empty()) return {{}, MC_EMPTY};

    // Separate non-empty sets
    std::vector<const std::vector<int32_t>*> non_empty;
    for (const auto& s : sets)
        if (!s.empty()) non_empty.push_back(&s);
    if (non_empty.empty()) return {{}, MC_EMPTY};

    bool has_empty_set = (non_empty.size() < sets.size());

    // 1. Intersection of ALL sets (only meaningful if none are empty)
    if (!has_empty_set) {
        std::vector<int32_t> isect = *non_empty[0];
        for (size_t i = 1; i < non_empty.size() && !isect.empty(); i++)
            isect = vec_intersect(isect, *non_empty[i]);
        if (!isect.empty())
            return {std::move(isect), MC_INTERSECTION};
    } else {
        // 2. Intersection of non-empty sets only
        if (non_empty.size() == 1) {
            return {*non_empty[0], MC_INTERSECTION_NONEMPTY};
        }
        std::vector<int32_t> isect = *non_empty[0];
        for (size_t i = 1; i < non_empty.size() && !isect.empty(); i++)
            isect = vec_intersect(isect, *non_empty[i]);
        if (!isect.empty())
            return {std::move(isect), MC_INTERSECTION_NONEMPTY};
    }

    // 3. Union of all sets
    std::unordered_set<int32_t> all;
    for (const auto& s : sets)
        for (int32_t v : s) all.insert(v);
    if (!all.empty()) {
        std::vector<int32_t> sorted_union(all.begin(), all.end());
        std::sort(sorted_union.begin(), sorted_union.end());
        return {std::move(sorted_union), MC_UNION};
    }
    return {{}, MC_EMPTY};
}

// Intrachromosomal chimera detection
inline ChimeraResult detect_chimera(
    const std::vector<ExonBlock>& exons,
    const std::vector<std::vector<int32_t>>& exon_t_sets)
{
    // Filter to blocks with non-empty transcript sets
    std::vector<int> item_idx;
    for (size_t i = 0; i < exon_t_sets.size(); i++)
        if (!exon_t_sets[i].empty()) item_idx.push_back(static_cast<int>(i));

    int n = static_cast<int>(item_idx.size());
    if (n <= 1) return {CHIMERA_NONE, -1};

    // Union-find
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);

    auto find = [&](int x) -> int {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra != rb) parent[ra] = rb;
    };

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (has_intersection(exon_t_sets[item_idx[i]],
                                 exon_t_sets[item_idx[j]]))
                unite(i, j);

    // Group into connected components
    std::unordered_map<int, std::vector<int>> components;
    for (int i = 0; i < n; i++)
        components[find(i)].push_back(i);
    if (components.size() <= 1) return {CHIMERA_NONE, -1};

    // Strand characterisation
    std::unordered_set<int32_t> unique_strands;
    for (const auto& [root, members] : components) {
        int32_t strand = STRAND_NONE;
        for (int idx : members) strand |= exons[item_idx[idx]].strand;
        unique_strands.insert(strand);
    }
    int32_t chimera_type = (unique_strands.size() == 1)
        ? CHIMERA_CIS_STRAND_SAME : CHIMERA_CIS_STRAND_DIFF;

    // Minimum gap between components
    std::vector<std::vector<int>> comp_list;
    comp_list.reserve(components.size());
    for (auto& [root, members] : components)
        comp_list.push_back(std::move(members));

    int32_t min_gap = std::numeric_limits<int32_t>::max();
    for (size_t ci = 0; ci < comp_list.size(); ci++) {
        for (size_t cj = ci + 1; cj < comp_list.size(); cj++) {
            for (int bi : comp_list[ci]) {
                for (int bj : comp_list[cj]) {
                    const auto& blki = exons[item_idx[bi]];
                    const auto& blkj = exons[item_idx[bj]];
                    int32_t gap;
                    if (blki.end <= blkj.start) gap = blkj.start - blki.end;
                    else if (blkj.end <= blki.start) gap = blki.start - blkj.end;
                    else gap = 0;
                    min_gap = std::min(min_gap, gap);
                }
            }
        }
    }
    return {chimera_type, min_gap};
}

}  // namespace rigel
