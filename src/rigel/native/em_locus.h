// em_locus.h — one locus's EM sub-problem: its components and its units' candidates.
//
// A locus has n_t transcript components (annotated transcripts and synthetic spans) in [0, n_t), one gDNA component
// at n_t, and, where some unit has a place outside every locus, the intergenic term at n_t + 1. The sub-problem is
// read from the locus's partition (the per-locus arrays `rigel.locus_partition` scatters) with every candidate
// remapped to its local component.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "constants.h"

namespace rigel {

// The intergenic term: a fixed rate no locus parameter carries. A multimapper's place with no transcript candidate is
// in no locus, so its gDNA comes at the library's intergenic rate `λ` (fragments per start, in the unit the capture
// efficiencies are read in), beside the locus's components, whose absolute rates are `M θ_c / L̃_c`. `M` is the
// locus's own data count, re-estimated with θ: the maximiser of `Σ_u log(Σ_c M θ_c q_uc / L̃_c + λ b_u) + Σ_d log(M θ_d)
// − M` under the grouped prior on θ gives `M = Σ_u Σ_{c in locus} r_uc + N_det`. In the locus's frame the term is one
// more component, last, of weight `log λ − log M`; it takes no share of the prior and none of θ's normalisation, and
// its slot of the EM's state holds `M` so SQUAREM accelerates the two together: `M / N` under MAP, where the state is
// a share, and `M` under VBEM, where it is a count. A locus with no intergenic place has no such component.
struct IntergenicTerm {
    int    idx = -1;        // its component index, the last one; -1: the locus has none
    double log_rate = -std::numeric_limits<double>::infinity();
    double n_all = 0.0;     // the locus's fragments: its units with a candidate and its deterministic ones
};

// The intergenic column's weight, from the locus's data count `M`.
inline double intergenic_log_weight(const IntergenicTerm& ig, double data_count) {
    return ig.log_rate - std::log(std::max(data_count, EM_LOG_EPSILON));
}

// The weights at a share theta, as MAP's E-step reads them: log(theta_c + eps) − log L̃_c per locus component, and the
// intergenic term's at the data count `data_count`.
inline void map_log_weights(const double* theta, const double* log_eff_len, int n_locus, const IntergenicTerm& ig,
                            double data_count, double* out) {
    for (int c = 0; c < n_locus; ++c) out[c] = std::log(theta[c] + EM_LOG_EPSILON) - log_eff_len[c];
    if (ig.idx >= 0) out[ig.idx] = intergenic_log_weight(ig, data_count);
}

// The weights at the converged theta, which the assignment reads: MAP's E-step's, but a component at a share of
// exactly 0 (the EM left it at its floor) is -inf, so it takes no posterior. At log(0 + eps) it would take about eps
// of every unit it is a candidate of: a count no fragment supports, and a candidate a whole-count draw can be handed.
inline void converged_log_weights(const double* theta, const double* log_eff_len, int n_locus,
                                  const IntergenicTerm& ig, double data_count, double* out) {
    map_log_weights(theta, log_eff_len, n_locus, ig, data_count, out);
    for (int c = 0; c < n_locus; ++c)
        if (c != ig.idx && !(theta[c] > 0.0)) out[c] = -std::numeric_limits<double>::infinity();
}

// One locus's partition: its candidate CSR and per-unit arrays, and its transcripts.
struct PartitionView {
    const int64_t* offsets;
    const int32_t* t_indices;
    const float*   log_liks;
    const float*   coverage_wts;
    const uint8_t* count_cols;
    const uint8_t* is_spliced;
    const float*   gdna_log_liks;
    const float*   intergenic_log_liks;  // per unit; -inf where a unit has no intergenic place
    int     n_units;
    int64_t n_candidates;
    const int32_t* transcript_indices;   // the locus's transcripts (`Locus.transcript_indices`)
    int n_transcripts;
};

struct LocalCandidate {
    int32_t local_comp;
    double  log_lik;
    double  cov_wt;
    uint8_t count_col;
};

// One locus's sub-problem, reused across the loci a thread solves.
struct LocusSubProblem {
    int n_t;           // number of transcripts in locus
    int n_components;  // n_t + 1, or n_t + 2 with the intergenic term
    int gdna_idx;      // = n_t     (single gDNA component)
    int ig_idx = -1;   // = n_t + 1 when a unit has an intergenic place (the intergenic term), else -1
    int n_local_units;
    bool has_gdna_candidate = false;

    // Local CSR
    std::vector<int64_t>  offsets;      // [n_local_units + 1]
    std::vector<int32_t>  t_indices;    // local component indices
    std::vector<double>   log_liks;
    std::vector<double>   coverage_wts;
    std::vector<uint8_t>  count_cols;

    // Per-component
    std::vector<double>   unambig_totals;   // [n_components]
    std::vector<double>   log_eff_len;      // [n_components] log L̃ per component
    // [n_components] the RNA prior's allocation weight, remapped from the per-transcript lane; empty when the
    // caller supplies none, and the prior then follows the evidence.
    std::vector<double>   component_rna_prior_weight;

    std::vector<int32_t>  local_to_global_t; // [n_t]
};

// Unit ui's candidates in buf, as local components: its transcripts that belong to the locus, then the gDNA component
// (n_t) when the unit is unspliced with a finite gDNA term, and the intergenic term (n_t + 1) when it has an
// intergenic place, sorted by component as the equivalence-class builder requires. The scorer supplies log h_G(ell_f)
// and the non-length terms for gDNA; the EM applies the locus's -log L̃_gDNA. Returns how many there are.
inline int unit_candidates(const LocusSubProblem& sub, const PartitionView& pv, int ui,
                           const std::vector<int32_t>& local_map, std::vector<LocalCandidate>& buf, bool& has_gdna) {
    const int64_t p_start = pv.offsets[ui], p_end = pv.offsets[ui + 1];
    const double gdna_ll = static_cast<double>(pv.gdna_log_liks[ui]);
    has_gdna = pv.is_spliced[ui] == 0 && std::isfinite(gdna_ll);
    const double ig_ll = sub.ig_idx >= 0 ? static_cast<double>(pv.intergenic_log_liks[ui])
                                         : -std::numeric_limits<double>::infinity();
    const bool has_ig = std::isfinite(ig_ll);
    buf.resize(static_cast<size_t>(p_end - p_start) + (has_gdna ? 1 : 0) + (has_ig ? 1 : 0));
    int k = 0;
    for (int64_t j = p_start; j < p_end; ++j)  // every candidate is one of the locus's transcripts (read_partitions)
        buf[k++] = {local_map[pv.t_indices[j]], static_cast<double>(pv.log_liks[j]),
                    static_cast<double>(pv.coverage_wts[j]), pv.count_cols[j]};
    if (has_gdna) buf[k++] = {sub.gdna_idx, gdna_ll, 1.0, 0};
    if (has_ig) buf[k++] = {sub.ig_idx, ig_ll, 1.0, 0};
    if (k > 1) {
        std::sort(buf.begin(), buf.begin() + k,
                  [](const LocalCandidate& a, const LocalCandidate& b) { return a.local_comp < b.local_comp; });
    }
    return k;
}

// The sub-problem's local CSR, every unit's candidates in turn. local_map maps a global transcript to its local
// component (-1 outside the locus).
inline void gather_candidates(LocusSubProblem& sub, const PartitionView& pv, const std::vector<int32_t>& local_map) {
    const int n_u = pv.n_units;
    // Sized to the worst case: the input candidates plus one gDNA and one intergenic candidate per unit.
    const size_t max_out = static_cast<size_t>(pv.n_candidates) + static_cast<size_t>(n_u) * (sub.ig_idx >= 0 ? 2 : 1);
    sub.t_indices.resize(max_out);
    sub.log_liks.resize(max_out);
    sub.coverage_wts.resize(max_out);
    sub.count_cols.resize(max_out);
    sub.offsets.resize(n_u + 1);
    sub.offsets[0] = 0;
    std::vector<LocalCandidate> buf;
    size_t cursor = 0;
    for (int ui = 0; ui < n_u; ++ui) {
        bool has_gdna;
        const int k = unit_candidates(sub, pv, ui, local_map, buf, has_gdna);
        sub.has_gdna_candidate = sub.has_gdna_candidate || has_gdna;
        for (int i = 0; i < k; ++i, ++cursor) {
            sub.t_indices[cursor] = buf[i].local_comp;
            sub.log_liks[cursor] = buf[i].log_lik;
            sub.coverage_wts[cursor] = buf[i].cov_wt;
            sub.count_cols[cursor] = buf[i].count_col;
        }
        sub.offsets[ui + 1] = static_cast<int64_t>(cursor);
    }
    sub.t_indices.resize(cursor);
    sub.log_liks.resize(cursor);
    sub.coverage_wts.resize(cursor);
    sub.count_cols.resize(cursor);
}

// The per-component arrays: each transcript's global index, unambiguous total, RNA prior weight (empty when the caller
// supplies none; the gDNA component is no recipient and keeps 0) and log effective length. A component's yield enters
// unfloored: a yield of 0 is a component with no start position, which cannot emit, so +inf here is -inf in every
// E-step weight, and a row no component can emit is skipped.
inline void fill_components(LocusSubProblem& sub, const PartitionView& pv, double gdna_eff_len,
                            const double* all_unambig_row_sums, const double* all_t_eff_lens,
                            const double* all_t_rna_prior_weight) {
    const int n_t = sub.n_t, nc = sub.n_components;
    const int32_t* t_arr = pv.transcript_indices;
    const double inf = std::numeric_limits<double>::infinity();
    sub.local_to_global_t.assign(t_arr, t_arr + n_t);
    sub.unambig_totals.assign(nc, 0.0);
    sub.log_eff_len.assign(nc, 0.0);
    sub.component_rna_prior_weight.clear();
    if (all_t_rna_prior_weight != nullptr) sub.component_rna_prior_weight.assign(nc, 0.0);
    for (int i = 0; i < n_t; ++i) {
        sub.unambig_totals[i] = all_unambig_row_sums[t_arr[i]];
        if (all_t_rna_prior_weight != nullptr) sub.component_rna_prior_weight[i] = all_t_rna_prior_weight[t_arr[i]];
        const double eff_len = all_t_eff_lens[t_arr[i]];
        sub.log_eff_len[i] = eff_len > 0.0 ? std::log(eff_len) : inf;
    }
    // gDNA: the locus's capture-contracted effective length (uncontracted off capture), as for the RNA ones.
    sub.log_eff_len[sub.gdna_idx] = gdna_eff_len > 0.0 ? std::log(gdna_eff_len) : inf;
}

// Extract one locus's sub-problem from its partition. `all_t_eff_lens` is the EM's capture-contracted effective length
// per transcript (computed in Python) and `gdna_eff_len` the locus's gDNA counterpart. `local_map` is the caller's
// scratch over the global transcripts, -1 everywhere on entry and on return.
inline void extract_locus_sub_problem(LocusSubProblem& sub, const PartitionView& pv, double gdna_eff_len,
                                      const double* all_unambig_row_sums, const double* all_t_eff_lens,
                                      const double* all_t_rna_prior_weight, std::vector<int32_t>& local_map) {
    const int n_t = pv.n_transcripts;
    bool any_intergenic = false;
    for (int ui = 0; ui < pv.n_units && !any_intergenic; ++ui)
        any_intergenic = std::isfinite(static_cast<double>(pv.intergenic_log_liks[ui]));
    sub.n_t = n_t;
    sub.n_local_units = pv.n_units;
    sub.gdna_idx = n_t;
    sub.ig_idx = any_intergenic ? n_t + 1 : -1;
    sub.n_components = n_t + (any_intergenic ? 2 : 1);
    sub.has_gdna_candidate = false;

    const int map_size = static_cast<int>(local_map.size());
    auto map_transcripts = [&](bool to_local) {
        for (int i = 0; i < n_t; ++i) {
            const int32_t gt = pv.transcript_indices[i];
            if (gt >= 0 && gt < map_size) local_map[gt] = to_local ? i : -1;
        }
    };
    map_transcripts(true);
    gather_candidates(sub, pv, local_map);
    map_transcripts(false);
    fill_components(sub, pv, gdna_eff_len, all_unambig_row_sums, all_t_eff_lens, all_t_rna_prior_weight);
}

// The locus's intergenic term, when a unit has an intergenic place: its component, the library's log rate, and N.
inline IntergenicTerm intergenic_term(const LocusSubProblem& sub, double log_rate) {
    IntergenicTerm ig;
    if (sub.ig_idx < 0) return ig;
    ig.idx = sub.ig_idx;
    ig.log_rate = log_rate;
    double n_all = 0.0;
    for (int u = 0; u < sub.n_local_units; ++u)
        if (sub.offsets[u + 1] > sub.offsets[u]) n_all += 1.0;
    for (int c = 0; c < sub.n_components; ++c) n_all += sub.unambig_totals[c];
    ig.n_all = n_all;
    return ig;
}

}  // namespace rigel
