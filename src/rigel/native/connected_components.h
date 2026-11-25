// connected_components.h — the loci: connected components of the fragment-to-transcript graph.
//
// Every candidate transcript of a unit is united with the unit's first (union-find with path halving and union by
// rank, O(N_candidates · α(N_transcripts))); each component's transcripts and units are then listed as sorted CSRs.
#pragma once

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "ndarray_util.h"

namespace nb = nanobind;

namespace rigel {

// Disjoint-set forest with path halving and union by rank.
struct UnionFind {
    std::vector<int32_t> parent;
    std::vector<int32_t> rank;

    explicit UnionFind(int32_t n) : parent(n), rank(n, 0) { std::iota(parent.begin(), parent.end(), 0); }

    int32_t find(int32_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    void unite(int32_t a, int32_t b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) ++rank[a];
    }
};

// Each transcript's component, numbered from 0 in order of each component's lowest transcript; -1 for a transcript no
// unit names.
inline std::vector<int32_t> transcript_components(const int64_t* offsets, const int32_t* t_idx, int64_t n_units,
                                                  int32_t n_transcripts, int32_t& n_comp) {
    UnionFind uf(n_transcripts);
    std::vector<bool> active(n_transcripts, false);
    for (int64_t u = 0; u < n_units; ++u) {
        int32_t first_t = -1;
        for (int64_t j = offsets[u]; j < offsets[u + 1]; ++j) {
            const int32_t t = t_idx[j];
            active[t] = true;
            if (first_t < 0) first_t = t;
            else uf.unite(first_t, t);
        }
    }
    std::vector<int32_t> labels(n_transcripts, -1);
    std::unordered_map<int32_t, int32_t> root_to_label;
    n_comp = 0;
    for (int32_t t = 0; t < n_transcripts; ++t) {
        if (!active[t]) continue;
        const auto [it, added] = root_to_label.try_emplace(uf.find(t), n_comp);
        if (added) ++n_comp;
        labels[t] = it->second;
    }
    return labels;
}

// Each unit's component: that of its first candidate, since every candidate is in one; -1 for a unit with none.
inline std::vector<int32_t> unit_components(const int64_t* offsets, const int32_t* t_idx, int64_t n_units,
                                            const std::vector<int32_t>& t_labels) {
    std::vector<int32_t> labels(n_units, -1);
    for (int64_t u = 0; u < n_units; ++u)
        if (offsets[u + 1] > offsets[u]) labels[u] = t_labels[t_idx[offsets[u]]];
    return labels;
}

// The indices grouped by label as a CSR, ascending within each group; a negative label is in no group. Returns the
// offsets (int64, n_groups + 1) and the flat indices (int32).
inline std::pair<nb::object, nb::object> group_by_label(const std::vector<int32_t>& labels, int32_t n_groups) {
    std::vector<int64_t> offsets(static_cast<size_t>(n_groups) + 1, 0);
    for (const int32_t c : labels)
        if (c >= 0) ++offsets[c + 1];
    for (int32_t c = 0; c < n_groups; ++c) offsets[c + 1] += offsets[c];
    std::vector<int32_t> flat(static_cast<size_t>(offsets[n_groups]));
    std::vector<int64_t> cursor(offsets.begin(), offsets.end() - 1);
    for (size_t i = 0; i < labels.size(); ++i)
        if (labels[i] >= 0) flat[cursor[labels[i]]++] = static_cast<int32_t>(i);
    return {vec_to_ndarray(std::move(offsets)), vec_to_ndarray(std::move(flat))};
}

// The loci of the global CSR (offsets int64[n_units + 1], t_indices int32[n_candidates]): (n_components,
// comp_t_offsets, comp_t_flat, comp_u_offsets, comp_u_flat), each component's transcript and unit indices as a sorted
// CSR.
// The CSR is checked here once, so the builders read every candidate unguarded: its offsets start at 0, never fall and
// end at its candidates' count, and every candidate is one of the `n_transcripts`.
inline nb::tuple connected_components(Array1<int64_t> offsets_arr, Array1<int32_t> t_indices_arr,
                                      int32_t n_transcripts) {
    const int64_t* offsets = offsets_arr.data();
    const int32_t* t_idx = t_indices_arr.data();
    const int64_t n_units = static_cast<int64_t>(offsets_arr.shape(0)) - 1;
    if (n_units < 0 || n_transcripts < 0 || offsets[0] != 0 ||
        offsets[n_units] != static_cast<int64_t>(t_indices_arr.shape(0)))
        throw std::invalid_argument("connected_components: the offsets must start at 0 and end at the candidates' count");
    for (int64_t u = 0; u < n_units; ++u)
        if (offsets[u + 1] < offsets[u]) throw std::invalid_argument("connected_components: the offsets fall");
    for (size_t j = 0; j < t_indices_arr.shape(0); ++j)
        if (t_idx[j] < 0 || t_idx[j] >= n_transcripts)
            throw std::invalid_argument("connected_components: candidate transcript " + std::to_string(t_idx[j]) +
                                        " is outside the " + std::to_string(n_transcripts));
    int32_t n_comp = 0;
    const std::vector<int32_t> t_labels = transcript_components(offsets, t_idx, n_units, n_transcripts, n_comp);
    const auto [t_offsets, t_flat] = group_by_label(t_labels, n_comp);
    const auto [u_offsets, u_flat] = group_by_label(unit_components(offsets, t_idx, n_units, t_labels), n_comp);
    return nb::make_tuple(n_comp, t_offsets, t_flat, u_offsets, u_flat);
}

}  // namespace rigel
