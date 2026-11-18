// em_assign.h — the locus EM's posterior assignment at the converged theta.
//
// FRACTIONAL scatters every unit's posterior as it stands. SAMPLE gives each unit to one component, count first: each
// component's fractional count in the locus (transcripts and gDNA alike) is rounded by largest remainder, then one pass
// draws each unit from its posterior re-weighted per candidate by (count still owed / posterior mass still to come). A
// unit whose candidates are all full keeps its most probable one (stranded) and a repair along augmenting paths puts
// every count back on target, or within one of its fractional count where the exact targets are unreachable. The
// result is feasible, not optimal: no candidate is dropped per fragment, and the repair ignores posteriors.
// Gate: tests/test_estimator.py (TestWholeCounts).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

#include "constants.h"
#include "em_locus.h"

namespace rigel {

// SplitMix64: a small deterministic PRNG, one stream per locus.
struct SplitMix64 {
    static constexpr uint64_t GAMMA = 0x9e3779b97f4a7c15ULL;  // 2^64/φ, the golden-ratio increment
    uint64_t state;

    explicit SplitMix64(uint64_t seed) : state(seed) {}

    // Stream `index` of a run's seed: each locus draws from its own.
    static SplitMix64 stream(uint64_t seed, uint64_t index) { return SplitMix64(seed ^ (index * GAMMA)); }

    uint64_t next() {
        state += GAMMA;
        uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    // A uniform double in [0, 1).
    double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
};

// Every unit's posterior at the components' log-weights, aligned with sub.t_indices. A unit no component can emit
// keeps all zeros and is assigned nothing, in either mode.
inline std::vector<double> unit_posteriors(const LocusSubProblem& sub, const std::vector<double>& log_weights) {
    std::vector<double> post(static_cast<size_t>(sub.offsets[sub.n_local_units]), 0.0);
    for (int ui = 0; ui < sub.n_local_units; ++ui) {
        const int64_t s = sub.offsets[ui], e = sub.offsets[ui + 1];
        double max_val = -std::numeric_limits<double>::infinity();
        for (int64_t k = s; k < e; ++k) max_val = std::max(max_val, sub.log_liks[k] + log_weights[sub.t_indices[k]]);
        if (max_val == -std::numeric_limits<double>::infinity()) continue;
        double sum_exp = 0.0;
        for (int64_t k = s; k < e; ++k) {
            post[k] = std::exp(sub.log_liks[k] + log_weights[sub.t_indices[k]] - max_val);
            sum_exp += post[k];
        }
        if (sum_exp > 0.0 && std::isfinite(sum_exp)) {
            double inv = 1.0 / sum_exp;
            for (int64_t k = s; k < e; ++k) post[k] *= inv;
        } else {
            for (int64_t k = s; k < e; ++k) post[k] = 0.0;
        }
    }
    return post;
}

// SAMPLE's whole counts and the draw that meets them. `choice[ui]` is the unit's candidate position, -1 for a unit
// with no posterior.
class WholeCountDraw {
public:
    WholeCountDraw(const LocusSubProblem& sub, const std::vector<double>& post)
        : sub_(sub), post_(post), nc_(sub.n_components), to_come_(nc_, 0.0), quota_(nc_, 0) {
        for (size_t k = 0; k < post.size(); ++k) to_come_[sub.t_indices[k]] += post[k];
        n_frac_ = to_come_;
        round_counts();
    }

    // The draw, then the repair of any count the stranded units leave off target.
    std::vector<int> run(SplitMix64& rng, int& stranded, int& unrepaired) {
        draw(rng, stranded);
        unrepaired = 0;
        if (stranded > 0) repair(unrepaired);
        return std::move(choice_);
    }

private:
    int comp(int ui, int pos) const { return sub_.t_indices[sub_.offsets[ui] + pos]; }
    int64_t count(int c) const { return static_cast<int64_t>(members_[c].size()); }

    // Each component's fractional count, rounded within the locus by largest remainder: the fractional counts sum to
    // the number of units with a posterior, an integer up to rounding.
    void round_counts() {
        double total = 0.0;
        int64_t floors = 0;
        for (int c = 0; c < nc_; ++c) {
            total += to_come_[c];
            quota_[c] = static_cast<int64_t>(std::floor(to_come_[c]));
            floors += quota_[c];
        }
        const int64_t short_by = std::llround(total) - floors;
        if (short_by <= 0) return;
        std::vector<int> order(nc_);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            double ra = to_come_[a] - std::floor(to_come_[a]);
            double rb = to_come_[b] - std::floor(to_come_[b]);
            if (ra != rb) return ra > rb;                                        // the largest remainder first
            if (to_come_[a] != to_come_[b]) return to_come_[a] > to_come_[b];    // a tie to the larger count
            return a < b;
        });
        for (int64_t i = 0; i < short_by && i < nc_; ++i) quota_[order[i]] += 1;
    }

    // Units in order, each from its own posterior re-weighted by (count owed / mass to come). A unit whose candidates
    // all hold their counts keeps its most probable one, for now, and is stranded.
    void draw(SplitMix64& rng, int& stranded) {
        std::vector<int64_t> owed(quota_);
        auto drawable_at = [&](int64_t k) { return owed[sub_.t_indices[k]] > 0 && post_[k] > 0.0; };
        auto pull = [&](int32_t c) { return static_cast<double>(owed[c]) / std::max(to_come_[c], EM_LOG_EPSILON); };
        stranded = 0;
        choice_.assign(sub_.n_local_units, -1);
        for (int ui = 0; ui < sub_.n_local_units; ++ui) {
            const int64_t s = sub_.offsets[ui], e = sub_.offsets[ui + 1];
            bool has_posterior = false;
            double drawable = 0.0;
            for (int64_t k = s; k < e; ++k) {
                if (post_[k] > 0.0) has_posterior = true;
                if (drawable_at(k)) drawable += post_[k] * pull(sub_.t_indices[k]);
            }
            if (!has_posterior) continue;
            int pick = -1;
            if (drawable > 0.0) {
                double target = rng.uniform() * drawable, acc = 0.0;
                for (int64_t k = s; k < e; ++k) {
                    if (!drawable_at(k)) continue;
                    acc += post_[k] * pull(sub_.t_indices[k]);
                    pick = static_cast<int>(k - s);   // the last drawable one, if rounding leaves target >= acc
                    if (target < acc) break;
                }
            } else {
                double best = -1.0;
                for (int64_t k = s; k < e; ++k) {
                    if (post_[k] > best) { best = post_[k]; pick = static_cast<int>(k - s); }
                }
                ++stranded;
            }
            choice_[ui] = pick;
            owed[sub_.t_indices[s + pick]] -= 1;
            for (int64_t k = s; k < e; ++k) to_come_[sub_.t_indices[k]] -= post_[k];
        }
    }

    // Puts the counts back on their rounded targets along augmenting paths. A path exists whenever the targets are
    // reachable at all; the units no path places are counted in `unrepaired`. Counts within one of every fractional
    // count are always reachable: the posterior is a fractional assignment inside [floor(n_c), ceil(n_c)], so by flow
    // integrality an integral one exists, and the same paths find it: counts above their ceiling shed, then counts
    // below their floor draw; each path fixes one unit of violation and creates none.
    void repair(int& unrepaired) {
        members_.assign(nc_, {});
        slot_.assign(sub_.n_local_units, -1);
        for (int ui = 0; ui < sub_.n_local_units; ++ui) {
            if (choice_[ui] < 0) continue;
            const int c = comp(ui, choice_[ui]);
            slot_[ui] = static_cast<int>(members_[c].size());
            members_[c].push_back(ui);
        }
        seen_.assign(nc_, -1);
        via_unit_.assign(nc_, -1);
        via_pos_.assign(nc_, -1);
        from_comp_.assign(nc_, -1);
        shed_to(quota_);
        for (int c = 0; c < nc_; ++c) unrepaired += static_cast<int>(std::max<int64_t>(count(c) - quota_[c], 0));
        if (unrepaired == 0) return;
        std::vector<int64_t> lo(nc_), hi(nc_);
        for (int c = 0; c < nc_; ++c) {
            lo[c] = static_cast<int64_t>(std::floor(n_frac_[c]));
            hi[c] = static_cast<int64_t>(std::ceil(n_frac_[c]));
        }
        shed_to(hi);
        fill_to(lo);
    }

    // Unit ui moves to its candidate at position to_pos.
    void move(int ui, int to_pos) {
        const int from = comp(ui, choice_[ui]), to = comp(ui, to_pos);
        const int last = members_[from].back();
        members_[from][slot_[ui]] = last;
        slot_[last] = slot_[ui];
        members_[from].pop_back();
        slot_[ui] = static_cast<int>(members_[to].size());
        members_[to].push_back(ui);
        choice_[ui] = to_pos;
    }

    // While a component holds more than its target, a breadth-first search links it to one below its target through
    // units that can move (each only to another of its own candidates with a posterior above zero), and moving the
    // chain puts one unit back.
    void shed_to(const std::vector<int64_t>& target) {
        for (int c0 = 0; c0 < nc_; ++c0) {
            while (count(c0) > target[c0]) {
                const int found = forward_path(c0, target);
                if (found < 0) break;   // no component under its target is reachable from c0
                for (int y = found; y != c0;) {
                    const int x = from_comp_[y];
                    move(via_unit_[y], via_pos_[y]);
                    y = x;
                }
            }
        }
    }

    int forward_path(int c0, const std::vector<int64_t>& target) {
        ++search_;
        queue_.assign(1, c0);
        seen_[c0] = search_;
        for (size_t qi = 0; qi < queue_.size(); ++qi) {
            const int x = queue_[qi];
            for (int v : members_[x]) {
                const int64_t s = sub_.offsets[v], e = sub_.offsets[v + 1];
                for (int64_t k = s; k < e; ++k) {
                    const int y = sub_.t_indices[k];
                    if (post_[k] <= 0.0 || seen_[y] == search_) continue;
                    seen_[y] = search_;
                    from_comp_[y] = x;
                    via_unit_[y] = v;
                    via_pos_[y] = static_cast<int>(k - s);
                    if (count(y) < target[y]) return y;
                    queue_.push_back(y);
                }
            }
        }
        return -1;
    }

    // While a component holds fewer than its floor, a backward search through the units that could move into it finds
    // a component above its floor, and moving the chain brings one unit in.
    void fill_to(const std::vector<int64_t>& lo) {
        std::vector<std::vector<std::pair<int, int>>> into(nc_);   // component -> (unit, candidate position)
        for (int ui = 0; ui < sub_.n_local_units; ++ui) {
            const int64_t s = sub_.offsets[ui];
            for (int64_t k = s; k < sub_.offsets[ui + 1]; ++k)
                if (post_[k] > 0.0) into[sub_.t_indices[k]].emplace_back(ui, static_cast<int>(k - s));
        }
        for (int c0 = 0; c0 < nc_; ++c0) {
            while (count(c0) < lo[c0]) {
                const int found = backward_path(c0, lo, into);
                if (found < 0) break;
                for (int x = found; x != c0;) {
                    const int y = from_comp_[x];
                    move(via_unit_[x], via_pos_[x]);
                    x = y;
                }
            }
        }
    }

    int backward_path(int c0, const std::vector<int64_t>& lo, const std::vector<std::vector<std::pair<int, int>>>& into) {
        ++search_;
        queue_.assign(1, c0);
        seen_[c0] = search_;
        for (size_t qi = 0; qi < queue_.size(); ++qi) {
            const int y = queue_[qi];
            for (const auto& [v, pos] : into[y]) {
                if (choice_[v] < 0) continue;
                const int x = comp(v, choice_[v]);
                if (seen_[x] == search_) continue;
                seen_[x] = search_;
                from_comp_[x] = y;        // v moves x -> y
                via_unit_[x] = v;
                via_pos_[x] = pos;
                if (count(x) > lo[x]) return x;
                queue_.push_back(x);
            }
        }
        return -1;
    }

    const LocusSubProblem& sub_;
    const std::vector<double>& post_;
    const int nc_;
    std::vector<double> to_come_;   // each component's posterior mass among the units not yet drawn
    std::vector<double> n_frac_;    // the fractional counts themselves
    std::vector<int64_t> quota_;
    std::vector<int> choice_;
    // the repair's state: each component's units, each unit's place among them, and the searches' marks
    std::vector<std::vector<int>> members_;
    std::vector<int> slot_, seen_, via_unit_, via_pos_, from_comp_, queue_;
    int search_ = 0;
};

// The accumulators every locus's assignment adds into: the per-transcript count matrix (n_transcripts, n_cols), the
// posterior sums the confidence metrics read, and the optional per-unit annotation (nullptr when not requested).
struct AssignmentSinks {
    double* em_counts;
    double* posterior_sum;
    double* n_assigned;
    int n_transcripts;
    int n_cols;
    int32_t* winner_tid;
    float*   winner_post;
    int32_t* n_candidates;
};

// One locus's assigned mass: its transcript components (annotated transcripts and synthetic spans alike; the caller
// separates the two), its gDNA, and the intergenic term outside the locus; and SAMPLE's units the draw stranded and
// the repair could not place.
struct LocusAssignment {
    double rna = 0.0, gdna = 0.0, intergenic = 0.0;
    int stranded = 0, unrepaired = 0;
};

// A unit's annotation: its candidate count, and its winner (SAMPLE's draw; FRACTIONAL's most probable candidate) as a
// transcript or an EM_WINNER_* code with the winner's posterior.
inline void annotate_unit(const LocusSubProblem& sub, const IntergenicTerm& ig, const double* posteriors, int64_t s,
                          int seg_len, int winner, bool fractional, const AssignmentSinks& out, int at) {
    if (fractional) {
        double best_post = -1.0;
        for (int j = 0; j < seg_len; ++j) {
            if (posteriors[j] > best_post) {
                best_post = posteriors[j];
                winner = j;
            }
        }
    }
    out.n_candidates[at] = seg_len;
    if (winner < 0) {
        out.winner_tid[at] = EM_WINNER_NONE;
        out.winner_post[at] = 0.0f;
        return;
    }
    const int32_t comp = sub.t_indices[s + winner];
    out.winner_tid[at] = comp < sub.n_t ? sub.local_to_global_t[comp]
                       : comp == ig.idx ? EM_WINNER_INTERGENIC
                       : EM_WINNER_GDNA;
    out.winner_post[at] = static_cast<float>(posteriors[winner]);
}

// Assign each unit's posterior at the converged theta and scatter it into the accumulators. The weights are MAP's
// E-step's at theta, whichever EM ran, with the intergenic term's at the locus's converged data count. rng is this
// locus's stream (SAMPLE only); `unit_offset` is the locus's first unit in the per-unit annotation.
inline LocusAssignment assign_posteriors(const LocusSubProblem& sub, const double* theta, const IntergenicTerm& ig,
                                         double data_count, int assignment_mode, SplitMix64& rng,
                                         const AssignmentSinks& out, int unit_offset) {
    std::vector<double> log_weights(sub.n_components);
    map_log_weights(theta, sub.log_eff_len.data(), sub.n_components, ig, data_count, log_weights.data());
    const std::vector<double> post = unit_posteriors(sub, log_weights);
    const bool fractional = assignment_mode == ASSIGN_FRACTIONAL;
    LocusAssignment res;
    std::vector<int> choice;
    if (!fractional) choice = WholeCountDraw(sub, post).run(rng, res.stranded, res.unrepaired);

    for (int ui = 0; ui < sub.n_local_units; ++ui) {
        const int64_t s = sub.offsets[ui];
        const int seg_len = static_cast<int>(sub.offsets[ui + 1] - s);
        if (seg_len == 0) continue;
        const double* posteriors = post.data() + s;
        const int winner = fractional ? -1 : choice[ui];   // SAMPLE's candidate position within the unit
        if (out.winner_tid != nullptr)
            annotate_unit(sub, ig, posteriors, s, seg_len, winner, fractional, out, unit_offset + ui);
        for (int j = 0; j < seg_len; ++j) {
            const double p = fractional ? posteriors[j] : (j == winner ? 1.0 : 0.0);
            if (p == 0.0) continue;
            const int32_t comp = sub.t_indices[s + j];
            if (comp < sub.n_t) {
                const int32_t global_t = sub.local_to_global_t[comp];
                const uint8_t col = sub.count_cols[s + j];
                if (global_t < 0 || global_t >= out.n_transcripts || col >= out.n_cols) continue;
                out.em_counts[global_t * out.n_cols + col] += p;
                res.rna += p;
                out.posterior_sum[global_t] += posteriors[j] * posteriors[j];
                out.n_assigned[global_t] += posteriors[j];
            } else if (comp == ig.idx) {
                res.intergenic += p;
            } else {
                res.gdna += p;
            }
        }
    }
    return res;
}

}  // namespace rigel
