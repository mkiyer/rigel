// em_squarem.h — the locus EM under SQUAREM: MAP-EM or VBEM with the grouped calibration prior.
//
// The two variants share every step but three: MAP's state is the normalised share theta and VBEM's the Dirichlet
// alpha, so they differ in the E-step's weights (log theta against digamma(alpha)), the floor a component keeps (0
// against VBEM_SQUAREM_PRIOR_FLOOR) and the normalisation (MAP renormalises its state after every step).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "constants.h"
#include "em_estep.h"
#include "em_locus.h"
#include "em_prior.h"
#include "thread_pool.h"

namespace rigel {

// digamma by its asymptotic series, after shifting x up to 6 where the series is accurate.
inline double digamma(double x) {
    double result = 0.0;
    while (x < 6.0) {
        result -= 1.0 / x;
        x += 1.0;
    }
    // ψ(x) ≈ ln(x) - 1/(2x) - Σ B_{2k}/(2k·x^{2k})
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    result += std::log(x) - 0.5 * inv_x
        - inv_x2 * (1.0/12.0
        - inv_x2 * (1.0/120.0
        - inv_x2 * (1.0/252.0
        - inv_x2 * (1.0/240.0
        - inv_x2 * (1.0/132.0
        - inv_x2 * (691.0/32760.0
        - inv_x2 * (1.0/12.0)))))));
    return result;
}

// The warm start's counts, projected by the grouped prior: each component's deterministic count plus each ambiguous
// unit's coverage-weighted shares. The intergenic slot, which SQUAREM seeds itself, starts at 0 and takes no share.
inline void compute_grouped_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const AggregatePrior& aggregate_prior,
    double*       warm_counts_out,
    int           n_components_all,
    int           ig_idx)
{
    const int n_components = ig_idx >= 0 ? ig_idx : n_components_all;
    std::fill(warm_counts_out, warm_counts_out + n_components_all, 0.0);
    std::vector<double> warm_raw(static_cast<size_t>(n_components_all), 0.0);
    std::copy(unambig_totals, unambig_totals + n_components, warm_raw.begin());
    // the coverage shares, each unit's weights already normalised within its candidate set (`build_equiv_classes`)
    for (const auto& ec : ec_data) {
        const int k = ec.k;
        const int32_t* cidx = ec.comp_idx.data();
        const double* wt = ec.wt_flat.data();
        for (int i = 0; i < ec.n; ++i)
            for (int j = 0; j < k; ++j)
                if (cidx[j] != ig_idx) warm_raw[cidx[j]] += wt[i * k + j];
    }
    apply_grouped_prior_update(warm_raw.data(), nullptr, aggregate_prior, warm_counts_out, n_components);
}

// SQUAREM step length: backtracked, never clamped. SQUAREM jumps to `state0 + 2·step·r + step²·v`, which can carry a
// shrinking component below the floor, where it never returns. Clamping there would let the accelerator, not the EM,
// decide which components live. So the step halves its excess over 1 (the plain double step `state2`) until no
// component `state2` keeps above the floor lands below it; a component `state2` itself takes below the floor does not
// hold the step back. `backtracks` counts the halvings.
inline double backtracked_squarem_step(
    double step,
    const std::vector<double>& state0,
    const std::vector<double>& r_vec,
    const std::vector<double>& v_vec,
    const std::vector<double>& state2,
    double floor,
    int& backtracks)
{
    const size_t nc = state0.size();
    while (step > 1.0) {
        bool feasible = true;
        for (size_t i = 0; i < nc; ++i) {
            if (!(state2[i] > floor)) continue;
            const double e = state0[i] + 2.0 * step * r_vec[i] + step * step * v_vec[i];
            if (!(e >= floor)) {
                feasible = false;
                break;
            }
        }
        if (feasible) return step;
        ++backtracks;
        const double shorter = 1.0 + 0.5 * (step - 1.0);
        step = (shorter < step) ? shorter : 1.0;  // exactly 1 once the excess no longer halves
    }
    return 1.0;
}

// SQUAREM's step statistics over one run.
struct SquaremStats {
    int iterations = 0;            // SQUAREM iterations completed
    double step_scale_sum = 0.0, step_scale_max = 0.0;
    int step_scale_count = 0;
    int backtracks = 0;            // halvings of an extrapolation step (backtracked_squarem_step)
    int nonfinite = 0;
    int stabilization_fails = 0;   // non-finite entries of a stabilising step

    double step_scale_mean() const {
        return step_scale_count > 0 ? step_scale_sum / static_cast<double>(step_scale_count) : 0.0;
    }
};

struct EMResult {
    std::vector<double> theta;
    double data_count = 0.0;  // the locus's own data count M at convergence, when it has an intergenic term
    SquaremStats stats;
};

// One locus's EM: its E-step, its prior, and the variant. The state's locus components are [0, n_locus); the
// intergenic slot, if any, is the last and holds M (VBEM) or M / N (MAP).
struct LocusEM {
    EStep& estep;
    const double* log_eff_len;
    const double* unambig_totals;
    const AggregatePrior& prior;
    const IntergenicTerm& ig;
    int n_components;
    bool vbem;

    int n_locus() const { return ig.idx >= 0 ? ig.idx : n_components; }
    double floor() const { return vbem ? VBEM_SQUAREM_PRIOR_FLOOR : 0.0; }

    // The E-step's weights at a state: VBEM's digamma(alpha_c) − digamma(Σ alpha) − log L̃_c, MAP's log(theta_c + eps)
    // − log L̃_c, and the intergenic term's at the state's data count.
    void log_weights(const double* state, double* out) const {
        const int n = n_locus();
        if (!vbem) {
            map_log_weights(state, log_eff_len, n, ig, ig.idx >= 0 ? state[ig.idx] * ig.n_all : 0.0, out);
            return;
        }
        double alpha_sum = 0.0;
        for (int i = 0; i < n; ++i) alpha_sum += state[i];
        double dg_sum = digamma(std::max(alpha_sum, EM_LOG_EPSILON));
        for (int i = 0; i < n; ++i) out[i] = digamma(std::max(state[i], EM_LOG_EPSILON)) - dg_sum - log_eff_len[i];
        if (ig.idx >= 0) out[ig.idx] = intergenic_log_weight(ig, state[ig.idx]);
    }

    // One EM step, state → next: the E-step at the state's weights, then the grouped M-step on the raw counts (the
    // deterministic plus the posterior totals), MAP's renormalised to a share. The intergenic slot takes the locus's
    // data count.
    void step(const double* state, double* em_totals, double* next) const {
        const int n = n_locus();
        std::vector<double> weights(static_cast<size_t>(n_components));
        log_weights(state, weights.data());
        std::fill(em_totals, em_totals + n_components, 0.0);
        estep.run(weights.data(), em_totals);
        std::vector<double> raw_counts(static_cast<size_t>(n_components));
        double data = 0.0;
        for (int i = 0; i < n_components; ++i) raw_counts[i] = unambig_totals[i] + em_totals[i];
        for (int i = 0; i < n; ++i) data += raw_counts[i];
        apply_grouped_prior_update(raw_counts.data(), state, prior, next, n);
        if (!vbem) {
            double total = 0.0;
            for (int i = 0; i < n; ++i) total += next[i];
            if (total > 0.0) {
                double inv_total = 1.0 / total;
                for (int i = 0; i < n; ++i) next[i] *= inv_total;
            }
        }
        if (ig.idx >= 0) next[ig.idx] = vbem ? data : data / ig.n_all;
    }

    // The starting state from the warm start's counts: VBEM's floored alpha, MAP's normalised share; the intergenic
    // slot starts with every fragment the locus's (M = N).
    std::vector<double> seed(const double* init_counts) const {
        const int n = n_locus();
        std::vector<double> state(static_cast<size_t>(n_components));
        if (vbem) {
            for (int i = 0; i < n_components; ++i) state[i] = std::max(init_counts[i], VBEM_SQUAREM_PRIOR_FLOOR);
            if (ig.idx >= 0) state[ig.idx] = ig.n_all;
            return state;
        }
        double total = 0.0;
        for (int i = 0; i < n; ++i) {
            state[i] = nonnegative_finite(init_counts[i]);
            total += state[i];
        }
        if (total > 0.0) {
            double inv = 1.0 / total;
            for (int i = 0; i < n; ++i) state[i] *= inv;
        }
        if (ig.idx >= 0) state[ig.idx] = 1.0;
        return state;
    }

    // MAP's extrapolated state, renormalised to a share; a state with no mass left falls back to the plain double step.
    void normalise(std::vector<double>& extrap, const std::vector<double>& state2) const {
        if (vbem) return;
        const int n = n_locus();
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += extrap[i];
        if (s > 0.0) {
            double inv = 1.0 / s;
            for (int i = 0; i < n; ++i) extrap[i] *= inv;
        } else {
            std::copy(state2.begin(), state2.end(), extrap.begin());
        }
    }

    // The change between two states the convergence test reads: the L1 distance of the normalised shares (MAP's
    // state is one already), plus the intergenic slot's change in M / N.
    double change(const std::vector<double>& before, const std::vector<double>& after) const {
        double delta = 0.0;
        if (!vbem) {
            for (int i = 0; i < n_components; ++i) delta += std::abs(after[i] - before[i]);
            return delta;
        }
        const int n = n_locus();
        double sum_old = 0.0, sum_new = 0.0;
        for (int i = 0; i < n; ++i) {
            sum_old += before[i];
            sum_new += after[i];
        }
        if (sum_old > 0.0 && sum_new > 0.0) {
            double inv_old = 1.0 / sum_old;
            double inv_new = 1.0 / sum_new;
            for (int i = 0; i < n; ++i) delta += std::abs(after[i] * inv_new - before[i] * inv_old);
        }
        if (ig.idx >= 0 && ig.n_all > 0.0) delta += std::abs(after[ig.idx] - before[ig.idx]) / ig.n_all;
        return delta;
    }

    // theta at the converged state, the intergenic slot zeroed, and the locus's data count M beside it. A VBEM
    // component at the floor is dead (the floor stands in for an alpha of 0), so its share is exactly 0, as a MAP
    // component's at its floor is.
    void read_out(const std::vector<double>& state, EMResult& out) const {
        const int n = n_locus();
        out.theta.assign(state.begin(), state.end());
        if (vbem) {
            double total = 0.0;
            for (int i = 0; i < n; ++i) total += state[i];
            if (total > 0.0) {
                double inv = 1.0 / total;
                for (int i = 0; i < n; ++i) out.theta[i] = state[i] > VBEM_SQUAREM_PRIOR_FLOOR ? state[i] * inv : 0.0;
            }
        }
        if (ig.idx >= 0) {
            out.data_count = vbem ? state[ig.idx] : state[ig.idx] * ig.n_all;
            out.theta[ig.idx] = 0.0;
        }
    }
};

// A state's entries floored: a non-finite entry is set to the floor and counted, as is one below it.
inline int clamp_to_floor(std::vector<double>& state, double floor) {
    int nonfinite = 0;
    for (double& x : state) {
        if (!std::isfinite(x)) {
            x = floor;
            ++nonfinite;
        } else if (x < floor) {
            x = floor;
        }
    }
    return nonfinite;
}

// SQUAREM's extrapolation from the two plain steps state0 → state1 → state2, floored and (MAP) renormalised.
inline void extrapolate(const LocusEM& em, const std::vector<double>& state0, const std::vector<double>& state1,
                        const std::vector<double>& state2, std::vector<double>& r_vec, std::vector<double>& v_vec,
                        std::vector<double>& extrap, SquaremStats& stats) {
    const size_t nc = state0.size();
    double sv2 = 0.0, srv = 0.0;
    for (size_t i = 0; i < nc; ++i) {
        r_vec[i] = state1[i] - state0[i];
        v_vec[i] = (state2[i] - state1[i]) - r_vec[i];
        sv2 += v_vec[i] * v_vec[i];
        srv += r_vec[i] * v_vec[i];
    }
    if (sv2 == 0.0) {
        std::copy(state2.begin(), state2.end(), extrap.begin());
        return;
    }
    double step = std::max(-srv / sv2, 1.0);
    if (!std::isfinite(step)) {
        step = 1.0;
        ++stats.nonfinite;
    }
    step = backtracked_squarem_step(step, state0, r_vec, v_vec, state2, em.floor(), stats.backtracks);
    stats.step_scale_sum += step;
    stats.step_scale_max = std::max(stats.step_scale_max, step);
    ++stats.step_scale_count;
    for (size_t i = 0; i < nc; ++i) {
        // At step 1 the jump is the plain double step: copied, because recomputing it rounds on the scale of state0
        // and could carry a small live component below the floor.
        extrap[i] = (step == 1.0) ? state2[i] : state0[i] + 2.0 * step * r_vec[i] + step * step * v_vec[i];
    }
    // Only a component the plain step itself carried to the floor can land below it here.
    stats.nonfinite += clamp_to_floor(extrap, em.floor());
    em.normalise(extrap, state2);
}

// The EM from the warm start's counts to convergence under SQUAREM: two plain steps, the extrapolation, and a
// stabilising step per iteration, until the state changes by less than convergence_delta.
inline EMResult run_squarem(
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,
    const double* unambig_totals,
    const AggregatePrior& aggregate_prior,
    const double* init_counts,
    int           n_components,
    int           max_iterations,
    double        convergence_delta,
    bool          use_vbem,
    ThreadPool*   pool,
    const IntergenicTerm& ig)
{
    EStep estep(ec_data, pool);
    const LocusEM em{estep, log_eff_len, unambig_totals, aggregate_prior, ig, n_components, use_vbem};
    const size_t nc = static_cast<size_t>(n_components);
    std::vector<double> state0 = em.seed(init_counts);
    std::vector<double> state1(nc), state2(nc), extrap(nc), next(nc), r_vec(nc), v_vec(nc), em_totals(nc, 0.0);
    EMResult out;
    SquaremStats& stats = out.stats;
    const int max_sq_iters = std::max(max_iterations / SQUAREM_BUDGET_DIVISOR, 1);
    for (int iter = 0; iter < max_sq_iters; ++iter) {
        em.step(state0.data(), em_totals.data(), state1.data());
        em.step(state1.data(), em_totals.data(), state2.data());
        extrapolate(em, state0, state1, state2, r_vec, v_vec, extrap, stats);
        em.step(extrap.data(), em_totals.data(), next.data());
        const int nonfinite = clamp_to_floor(next, em.floor());
        stats.nonfinite += nonfinite;
        stats.stabilization_fails += nonfinite;
        const double delta = em.change(state0, next);
        std::swap(state0, next);
        stats.iterations = iter + 1;
        if (delta < convergence_delta) break;
    }
    em.read_out(state0, out);
    return out;
}

}  // namespace rigel
