/**
 * em_solver.cpp — C++ EM solver for hulkrna locus-level abundance estimation.
 *
 * Replaces the Python hot path: _em_step, _vbem_step, _build_equiv_classes,
 * _apply_bias_correction_uniform, and the SQUAREM acceleration loop.
 *
 * Module: hulkrna._em_impl
 *
 * Build:
 *   Part of the hulkrna scikit-build-core build — see CMakeLists.txt.
 *   Pure C++17 + nanobind + numpy.  No external dependencies.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <unordered_map>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/tuple.h>

namespace nb = nanobind;

// ================================================================
// Constants (must match estimator.py exactly)
// ================================================================

static constexpr double EM_LOG_EPSILON = 1e-300;
static constexpr int    MAX_FRAG_LEN  = 1000000;
static constexpr int    SQUAREM_BUDGET_DIVISOR = 3;

// ================================================================
// Array type aliases
// ================================================================

using i32_1d = nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig>;
using i64_1d = nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig>;
using f64_1d = nb::ndarray<const double,  nb::ndim<1>, nb::c_contig>;

// Mutable variants for in-place modification
using f64_1d_mut = nb::ndarray<double, nb::ndim<1>, nb::c_contig>;

// ================================================================
// digamma — self-contained asymptotic series implementation
// ================================================================

static inline double digamma(double x) {
    // Handle non-positive values that can arise from EM_LOG_EPSILON
    if (x <= 0.0) {
        if (x == 0.0) return -1e300;
        // For very small positive values after clamping
        return -1e300;
    }

    double result = 0.0;
    // Shift x up so asymptotic series is accurate (need x >= 6)
    while (x < 6.0) {
        result -= 1.0 / x;
        x += 1.0;
    }
    // Asymptotic expansion: ψ(x) ≈ ln(x) - 1/(2x) - Σ B_{2k}/(2k·x^{2k})
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

// ================================================================
// Equivalence class — groups units with identical candidate sets
// ================================================================

struct EquivClass {
    std::vector<int32_t> comp_idx;  // k component indices
    std::vector<double>  ll_flat;   // n*k log-likelihoods (row-major)
    std::vector<double>  wt_flat;   // n*k coverage weights (row-major)
    mutable std::vector<double> scratch;  // n*k workspace (reused each iteration)
    int n;  // number of units in this class
    int k;  // number of components per unit
};

// ================================================================
// Equivalence class builder
// ================================================================
//
// Replaces Python _build_equiv_classes(). Groups CSR units by their
// candidate component set (the ordered tuple of t_indices per unit)
// into dense matrices for efficient batch processing.

// Hash for vector<int32_t> keys
struct VecHash {
    size_t operator()(const std::vector<int32_t>& v) const noexcept {
        // FNV-1a hash
        size_t h = 14695981039346656037ULL;
        for (int32_t val : v) {
            h ^= static_cast<size_t>(static_cast<uint32_t>(val));
            h *= 1099511628211ULL;
        }
        return h;
    }
};

static std::vector<EquivClass> build_equiv_classes(
    const int64_t* offsets,
    const int32_t* t_indices,
    const double*  log_liks,
    const double*  coverage_wts,
    int n_units)
{
    if (n_units == 0) return {};

    // Group units by candidate component set
    std::unordered_map<std::vector<int32_t>, std::vector<int>, VecHash> class_map;
    class_map.reserve(static_cast<size_t>(n_units));

    for (int u = 0; u < n_units; ++u) {
        auto start = static_cast<size_t>(offsets[u]);
        auto end   = static_cast<size_t>(offsets[u + 1]);
        if (start == end) continue;

        std::vector<int32_t> key(t_indices + start, t_indices + end);
        class_map[std::move(key)].push_back(static_cast<int>(start));
    }

    // Build dense matrices per class
    std::vector<EquivClass> result;
    result.reserve(class_map.size());

    for (auto& [key, start_list] : class_map) {
        int k = static_cast<int>(key.size());
        int n = static_cast<int>(start_list.size());

        EquivClass ec;
        ec.comp_idx = std::move(key);
        ec.n = n;
        ec.k = k;
        ec.ll_flat.resize(static_cast<size_t>(n) * k);
        ec.wt_flat.resize(static_cast<size_t>(n) * k);
        ec.scratch.resize(static_cast<size_t>(n) * k);

        for (int i = 0; i < n; ++i) {
            int s = start_list[i];
            for (int j = 0; j < k; ++j) {
                ec.ll_flat[static_cast<size_t>(i) * k + j] =
                    log_liks[static_cast<size_t>(s) + j];
                ec.wt_flat[static_cast<size_t>(i) * k + j] =
                    coverage_wts[static_cast<size_t>(s) + j];
            }
        }

        result.push_back(std::move(ec));
    }

    return result;
}

// ================================================================
// Bias correction (uniform fast-path)
// ================================================================
//
// Applies -log(max(L - frag_len + 1, 1)) to each candidate's log_lik.
// Replaces Python _apply_bias_correction_uniform().

static void apply_bias_correction_uniform(
    double*        log_liks,        // mutated in-place
    const int32_t* t_indices,
    const int32_t* tx_starts,
    const int32_t* tx_ends,
    const int64_t* profile_lengths,
    size_t         n_candidates)
{
    for (size_t i = 0; i < n_candidates; ++i) {
        int64_t frag_len = static_cast<int64_t>(tx_ends[i]) -
                           static_cast<int64_t>(tx_starts[i]);
        if (frag_len < 0) frag_len = 0;
        if (frag_len >= MAX_FRAG_LEN) frag_len = MAX_FRAG_LEN - 1;

        int64_t prof_len = profile_lengths[t_indices[i]];
        int64_t eff_len = prof_len - frag_len + 1;
        if (eff_len < 1) eff_len = 1;

        log_liks[i] -= std::log(static_cast<double>(eff_len));
    }
}

// ================================================================
// EM step kernel — the hot inner loop
// ================================================================
//
// For one equivalence class: computes posteriors from log_weights,
// normalizes rows via log-sum-exp, accumulates column sums into
// em_totals.

static inline void em_step_kernel(
    const EquivClass& ec,
    const double*     log_weights,   // [n_components]
    double*           em_totals)     // [n_components], accumulated
{
    const int n = ec.n;
    const int k = ec.k;
    const double* ll = ec.ll_flat.data();
    double* scratch = ec.scratch.data();
    const int32_t* cidx = ec.comp_idx.data();

    // 1. Add log_weights to log-likelihoods
    for (int i = 0; i < n; ++i) {
        const int offset = i * k;
        for (int j = 0; j < k; ++j) {
            scratch[offset + j] = ll[offset + j] + log_weights[cidx[j]];
        }
    }

    // 2. Log-sum-exp normalize (row-wise) + accumulate column sums
    for (int i = 0; i < n; ++i) {
        double* row = scratch + i * k;

        // Find row max
        double max_val = row[0];
        for (int j = 1; j < k; ++j) {
            if (row[j] > max_val) max_val = row[j];
        }

        // exp and sum
        double row_sum = 0.0;
        for (int j = 0; j < k; ++j) {
            row[j] = std::exp(row[j] - max_val);
            row_sum += row[j];
        }

        // Normalize
        if (row_sum > 0.0 && std::isfinite(row_sum)) {
            double inv_sum = 1.0 / row_sum;
            for (int j = 0; j < k; ++j) {
                row[j] *= inv_sum;
            }
        } else {
            for (int j = 0; j < k; ++j) {
                row[j] = 0.0;
            }
        }
    }

    // 3. Accumulate column sums into em_totals
    for (int j = 0; j < k; ++j) {
        double col_sum = 0.0;
        for (int i = 0; i < n; ++i) {
            col_sum += scratch[i * k + j];
        }
        em_totals[cidx[j]] += col_sum;
    }
}

// ================================================================
// MAP-EM step: theta → theta_new
// ================================================================

static void map_em_step(
    const double* theta,
    const std::vector<EquivClass>& ec_data,
    const double* log_eff_len,
    const double* unique_totals,
    const double* prior,
    double*       em_totals,    // zeroed then accumulated
    double*       theta_new,    // output: normalized
    int           n_components)
{
    // Compute log_weights = log(theta + epsilon) - log_eff_len
    std::vector<double> log_weights(static_cast<size_t>(n_components));
    for (int i = 0; i < n_components; ++i) {
        log_weights[i] = std::log(theta[i] + EM_LOG_EPSILON) - log_eff_len[i];
    }

    // Zero em_totals
    std::fill(em_totals, em_totals + n_components, 0.0);

    // E-step: accumulate posteriors
    for (const auto& ec : ec_data) {
        em_step_kernel(ec, log_weights.data(), em_totals);
    }

    // M-step: theta_new = (unique_totals + em_totals + prior), normalized
    double total = 0.0;
    for (int i = 0; i < n_components; ++i) {
        theta_new[i] = unique_totals[i] + em_totals[i] + prior[i];
        total += theta_new[i];
    }
    if (total > 0.0) {
        double inv_total = 1.0 / total;
        for (int i = 0; i < n_components; ++i) {
            theta_new[i] *= inv_total;
        }
    }
}

// ================================================================
// VBEM step: alpha → alpha_new
// ================================================================

static void vbem_step(
    const double* alpha,
    const std::vector<EquivClass>& ec_data,
    const double* log_eff_len,
    const double* unique_totals,
    const double* prior,
    double*       em_totals,
    double*       alpha_new,
    int           n_components)
{
    // Compute alpha_sum
    double alpha_sum = 0.0;
    for (int i = 0; i < n_components; ++i) {
        alpha_sum += alpha[i];
    }

    // Compute log_weights = digamma(max(alpha, eps)) - digamma(max(alpha_sum, eps)) - log_eff_len
    double dg_sum = digamma(std::max(alpha_sum, EM_LOG_EPSILON));
    std::vector<double> log_weights(static_cast<size_t>(n_components));
    for (int i = 0; i < n_components; ++i) {
        log_weights[i] = digamma(std::max(alpha[i], EM_LOG_EPSILON))
                        - dg_sum - log_eff_len[i];
    }

    // Zero em_totals
    std::fill(em_totals, em_totals + n_components, 0.0);

    // E-step
    for (const auto& ec : ec_data) {
        em_step_kernel(ec, log_weights.data(), em_totals);
    }

    // M-step: alpha_new = unique_totals + em_totals + prior (unnormalized)
    for (int i = 0; i < n_components; ++i) {
        alpha_new[i] = unique_totals[i] + em_totals[i] + prior[i];
    }
}

// ================================================================
// Coverage-weighted warm start + OVR prior
// ================================================================

static void compute_ovr_prior_and_warm_start(
    const std::vector<EquivClass>& ec_data,
    const double* unique_totals,
    const double* eligible,    // [n_components] 1.0 if eligible, 0.0 otherwise
    double        alpha_flat,
    double        gamma,
    double*       prior_out,       // [n_components] output
    double*       theta_init_out,  // [n_components] output
    int           n_components)
{
    // Initialize theta_init from unique_totals
    std::copy(unique_totals, unique_totals + n_components, theta_init_out);

    // Accumulate coverage totals from all equivalence classes
    std::vector<double> coverage_totals(static_cast<size_t>(n_components), 0.0);
    int n_ambiguous = 0;

    for (const auto& ec : ec_data) {
        const int n = ec.n;
        const int k = ec.k;
        n_ambiguous += n;
        const int32_t* cidx = ec.comp_idx.data();
        const double* wt = ec.wt_flat.data();

        for (int i = 0; i < n; ++i) {
            // Mask weights by eligibility, normalize row
            double row_sum = 0.0;
            for (int j = 0; j < k; ++j) {
                double w = wt[i * k + j] * eligible[cidx[j]];
                row_sum += w;
            }
            if (row_sum == 0.0) row_sum = 1.0;
            double inv_row_sum = 1.0 / row_sum;

            for (int j = 0; j < k; ++j) {
                double share = (wt[i * k + j] * eligible[cidx[j]]) * inv_row_sum;
                theta_init_out[cidx[j]] += share;
                coverage_totals[cidx[j]] += share;
            }
        }
    }

    // Build prior: alpha_flat + gamma-scaled OVR coverage weights
    if (n_ambiguous > 0 && gamma > 0.0) {
        double inv_n = gamma / static_cast<double>(n_ambiguous);
        for (int i = 0; i < n_components; ++i) {
            double ovr_i = coverage_totals[i] * inv_n;
            prior_out[i] = (eligible[i] > 0.0) ? (alpha_flat + ovr_i) : 0.0;
        }
    } else {
        for (int i = 0; i < n_components; ++i) {
            prior_out[i] = (eligible[i] > 0.0) ? alpha_flat : 0.0;
        }
    }
}

// ================================================================
// SQUAREM acceleration wrapper
// ================================================================

struct EMResult {
    std::vector<double> theta;
    std::vector<double> alpha;
    std::vector<double> em_totals;
};

static EMResult run_squarem(
    const std::vector<EquivClass>& ec_data,
    const double* log_eff_len,
    const double* unique_totals,
    double*       prior,
    const double* theta_init,
    int           n_components,
    int           max_iterations,
    double        convergence_delta,
    bool          use_vbem,
    double        prune_threshold,
    int           post_prune_iters)
{
    int max_sq_iters = std::max(max_iterations / SQUAREM_BUDGET_DIVISOR, 1);
    size_t nc = static_cast<size_t>(n_components);

    std::vector<double> em_totals(nc, 0.0);

    // Temporary vectors for SQUAREM
    std::vector<double> state0(nc);
    std::vector<double> state1(nc);
    std::vector<double> state2(nc);
    std::vector<double> state_extrap(nc);
    std::vector<double> state_new(nc);
    std::vector<double> r_vec(nc);
    std::vector<double> v_vec(nc);

    std::vector<double> theta(nc);
    std::vector<double> alpha_out(nc);

    if (use_vbem) {
        // ---- VBEM with SQUAREM acceleration ----
        // state = Dirichlet parameters alpha
        for (size_t i = 0; i < nc; ++i) {
            state0[i] = theta_init[i] + prior[i];  // initial alpha
        }

        for (int iter = 0; iter < max_sq_iters; ++iter) {
            // Two plain VBEM steps
            vbem_step(state0.data(), ec_data, log_eff_len,
                      unique_totals, prior, em_totals.data(),
                      state1.data(), n_components);

            vbem_step(state1.data(), ec_data, log_eff_len,
                      unique_totals, prior, em_totals.data(),
                      state2.data(), n_components);

            // SQUAREM extrapolation
            double sv2 = 0.0, srv = 0.0;
            for (size_t i = 0; i < nc; ++i) {
                r_vec[i] = state1[i] - state0[i];
                v_vec[i] = (state2[i] - state1[i]) - r_vec[i];
                sv2 += v_vec[i] * v_vec[i];
                srv += r_vec[i] * v_vec[i];
            }

            if (sv2 == 0.0) {
                std::copy(state2.begin(), state2.end(), state_extrap.begin());
            } else {
                double step = std::max(-srv / sv2, 1.0);
                for (size_t i = 0; i < nc; ++i) {
                    state_extrap[i] = state0[i] + 2.0 * step * r_vec[i]
                                    + step * step * v_vec[i];
                    // Alpha must stay positive
                    if (state_extrap[i] < EM_LOG_EPSILON)
                        state_extrap[i] = EM_LOG_EPSILON;
                }
            }

            // Stabilisation step
            vbem_step(state_extrap.data(), ec_data, log_eff_len,
                      unique_totals, prior, em_totals.data(),
                      state_new.data(), n_components);

            // Convergence check on normalized theta
            double sum_old = 0.0, sum_new = 0.0;
            for (size_t i = 0; i < nc; ++i) {
                sum_old += state0[i];
                sum_new += state_new[i];
            }
            double delta = 0.0;
            if (sum_old > 0.0 && sum_new > 0.0) {
                double inv_old = 1.0 / sum_old;
                double inv_new = 1.0 / sum_new;
                for (size_t i = 0; i < nc; ++i) {
                    delta += std::abs(state_new[i] * inv_new
                                    - state0[i] * inv_old);
                }
            }

            std::swap(state0, state_new);

            if (delta < convergence_delta) break;
        }

        // alpha_out = state0 (converged Dirichlet params)
        std::copy(state0.begin(), state0.end(), alpha_out.begin());

        // theta = normalized alpha
        double total = 0.0;
        for (size_t i = 0; i < nc; ++i) total += alpha_out[i];
        if (total > 0.0) {
            double inv = 1.0 / total;
            for (size_t i = 0; i < nc; ++i) theta[i] = alpha_out[i] * inv;
        } else {
            std::copy(alpha_out.begin(), alpha_out.end(), theta.begin());
        }

    } else {
        // ---- MAP-EM with SQUAREM acceleration ----
        // state = normalized theta
        double total = 0.0;
        for (size_t i = 0; i < nc; ++i) {
            state0[i] = theta_init[i] + prior[i];
            total += state0[i];
        }
        if (total > 0.0) {
            double inv = 1.0 / total;
            for (size_t i = 0; i < nc; ++i) state0[i] *= inv;
        }

        for (int iter = 0; iter < max_sq_iters; ++iter) {
            // Two plain EM steps
            map_em_step(state0.data(), ec_data, log_eff_len,
                        unique_totals, prior, em_totals.data(),
                        state1.data(), n_components);

            map_em_step(state1.data(), ec_data, log_eff_len,
                        unique_totals, prior, em_totals.data(),
                        state2.data(), n_components);

            // SQUAREM extrapolation
            double sv2 = 0.0, srv = 0.0;
            for (size_t i = 0; i < nc; ++i) {
                r_vec[i] = state1[i] - state0[i];
                v_vec[i] = (state2[i] - state1[i]) - r_vec[i];
                sv2 += v_vec[i] * v_vec[i];
                srv += r_vec[i] * v_vec[i];
            }

            if (sv2 == 0.0) {
                std::copy(state2.begin(), state2.end(), state_extrap.begin());
            } else {
                double alpha_step = std::max(-srv / sv2, 1.0);
                for (size_t i = 0; i < nc; ++i) {
                    state_extrap[i] = state0[i]
                        + 2.0 * alpha_step * r_vec[i]
                        + alpha_step * alpha_step * v_vec[i];
                    if (state_extrap[i] < 0.0) state_extrap[i] = 0.0;
                }
                double s = 0.0;
                for (size_t i = 0; i < nc; ++i) s += state_extrap[i];
                if (s > 0.0) {
                    double inv = 1.0 / s;
                    for (size_t i = 0; i < nc; ++i) state_extrap[i] *= inv;
                } else {
                    std::copy(state2.begin(), state2.end(),
                              state_extrap.begin());
                }
            }

            // Stabilisation step
            map_em_step(state_extrap.data(), ec_data, log_eff_len,
                        unique_totals, prior, em_totals.data(),
                        state_new.data(), n_components);

            // Convergence
            double delta = 0.0;
            for (size_t i = 0; i < nc; ++i) {
                delta += std::abs(state_new[i] - state0[i]);
            }

            std::swap(state0, state_new);

            if (delta < convergence_delta) break;
        }

        // theta = state0 (converged normalized theta)
        std::copy(state0.begin(), state0.end(), theta.begin());

        // alpha_out = unique_totals + em_totals + prior
        for (size_t i = 0; i < nc; ++i) {
            alpha_out[i] = unique_totals[i] + em_totals[i] + prior[i];
        }
    }

    // ----------------------------------------------------------
    // Post-EM pruning
    // ----------------------------------------------------------
    if (prune_threshold >= 0.0) {
        bool any_pruned = false;
        std::vector<bool> prune_mask(nc, false);

        for (size_t i = 0; i < nc; ++i) {
            double data_count = unique_totals[i] + em_totals[i];
            double denom = std::max(alpha_out[i], EM_LOG_EPSILON);
            double evidence_ratio = data_count / denom;
            if (unique_totals[i] == 0.0
                && evidence_ratio < prune_threshold) {
                prune_mask[i] = true;
                any_pruned = true;
            }
        }
        // Never prune gDNA (last component)
        prune_mask[nc - 1] = false;

        if (any_pruned) {
            for (size_t i = 0; i < nc; ++i) {
                if (prune_mask[i]) prior[i] = 0.0;
            }

            if (use_vbem) {
                for (size_t i = 0; i < nc; ++i) {
                    if (prune_mask[i]) alpha_out[i] = EM_LOG_EPSILON;
                }
                for (int p = 0; p < post_prune_iters; ++p) {
                    std::vector<double> a_new(nc);
                    vbem_step(alpha_out.data(), ec_data, log_eff_len,
                              unique_totals, prior, em_totals.data(),
                              a_new.data(), n_components);
                    std::swap(alpha_out, a_new);
                }
                double total = 0.0;
                for (size_t i = 0; i < nc; ++i) total += alpha_out[i];
                if (total > 0.0) {
                    double inv = 1.0 / total;
                    for (size_t i = 0; i < nc; ++i) theta[i] = alpha_out[i] * inv;
                } else {
                    std::copy(alpha_out.begin(), alpha_out.end(), theta.begin());
                }
            } else {
                for (size_t i = 0; i < nc; ++i) {
                    if (prune_mask[i]) theta[i] = 0.0;
                }
                double total = 0.0;
                for (size_t i = 0; i < nc; ++i) total += theta[i];
                if (total > 0.0) {
                    double inv = 1.0 / total;
                    for (size_t i = 0; i < nc; ++i) theta[i] *= inv;
                }
                for (int p = 0; p < post_prune_iters; ++p) {
                    std::vector<double> t_new(nc);
                    map_em_step(theta.data(), ec_data, log_eff_len,
                                unique_totals, prior, em_totals.data(),
                                t_new.data(), n_components);
                    std::swap(theta, t_new);
                }
                for (size_t i = 0; i < nc; ++i) {
                    alpha_out[i] = unique_totals[i] + em_totals[i] + prior[i];
                }
            }
        }
    }

    return { std::move(theta), std::move(alpha_out), std::move(em_totals) };
}

// ================================================================
// Top-level entry point: run_locus_em_native()
// ================================================================
//
// Takes CSR per-locus data + config, returns (theta, alpha, em_totals).
// Replaces the entire body of AbundanceEstimator.run_locus_em() from
// bias correction through post-prune.

static std::tuple<nb::ndarray<nb::numpy, double, nb::ndim<1>>,
                  nb::ndarray<nb::numpy, double, nb::ndim<1>>,
                  nb::ndarray<nb::numpy, double, nb::ndim<1>>>
run_locus_em_native(
    // CSR data
    i64_1d offsets,
    i32_1d t_indices,
    f64_1d_mut log_liks,        // mutated in-place by bias correction
    f64_1d coverage_wts,
    i32_1d tx_starts,
    i32_1d tx_ends,
    i64_1d bias_profiles,
    // Per-component vectors
    f64_1d unique_totals_arr,
    f64_1d effective_lens,
    f64_1d prior_eligible,
    // Scalar config
    int    n_components,
    double prior_alpha,
    double prior_gamma,
    int    max_iterations,
    double convergence_delta,
    bool   use_vbem,
    double prune_threshold,
    int    post_prune_iters)
{
    size_t nc = static_cast<size_t>(n_components);
    size_t n_candidates = t_indices.shape(0);
    int n_units = static_cast<int>(offsets.shape(0)) - 1;

    // Get raw pointers
    const int64_t*  off_ptr  = offsets.data();
    const int32_t*  ti_ptr   = t_indices.data();
    double*         ll_ptr   = log_liks.data();
    const double*   cw_ptr   = coverage_wts.data();
    const int32_t*  txs_ptr  = tx_starts.data();
    const int32_t*  txe_ptr  = tx_ends.data();
    const int64_t*  bp_ptr   = bias_profiles.data();
    const double*   ut_ptr   = unique_totals_arr.data();
    const double*   el_ptr   = effective_lens.data();
    const double*   pe_ptr   = prior_eligible.data();

    // Copy unique_totals (we need a mutable copy)
    std::vector<double> unique_totals(ut_ptr, ut_ptr + nc);

    // 1. Apply bias correction (uniform fast-path)
    if (n_candidates > 0) {
        apply_bias_correction_uniform(
            ll_ptr, ti_ptr, txs_ptr, txe_ptr, bp_ptr, n_candidates);
    }

    // 2. Handle empty locus
    if (n_units == 0 || n_candidates == 0) {
        // theta = (unique_totals + prior) / sum, alpha = unique_totals + prior
        std::vector<double> alpha(nc);
        double total = 0.0;
        for (size_t i = 0; i < nc; ++i) {
            double p = (pe_ptr[i] > 0.0) ? prior_alpha : 0.0;
            alpha[i] = unique_totals[i] + p;
            total += alpha[i];
        }
        std::vector<double> theta(nc);
        if (total > 0.0) {
            double inv = 1.0 / total;
            for (size_t i = 0; i < nc; ++i) theta[i] = alpha[i] * inv;
        } else {
            std::copy(alpha.begin(), alpha.end(), theta.begin());
        }
        std::vector<double> em_totals(nc, 0.0);

        // Return as numpy arrays
        auto* theta_data = new double[nc];
        auto* alpha_data = new double[nc];
        auto* em_data    = new double[nc];
        std::copy(theta.begin(), theta.end(), theta_data);
        std::copy(alpha.begin(), alpha.end(), alpha_data);
        std::copy(em_totals.begin(), em_totals.end(), em_data);

        size_t shape[1] = { nc };
        nb::capsule theta_owner(theta_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
        nb::capsule alpha_owner(alpha_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
        nb::capsule em_owner(em_data,       [](void* p) noexcept { delete[] static_cast<double*>(p); });

        return std::make_tuple(
            nb::ndarray<nb::numpy, double, nb::ndim<1>>(theta_data, 1, shape, std::move(theta_owner)),
            nb::ndarray<nb::numpy, double, nb::ndim<1>>(alpha_data, 1, shape, std::move(alpha_owner)),
            nb::ndarray<nb::numpy, double, nb::ndim<1>>(em_data,    1, shape, std::move(em_owner))
        );
    }

    // 3. Compute log(effective_lengths)
    std::vector<double> log_eff_len(nc);
    for (size_t i = 0; i < nc; ++i) {
        log_eff_len[i] = std::log(el_ptr[i]);
    }

    // 4. Build equivalence classes
    auto ec_data = build_equiv_classes(
        off_ptr, ti_ptr, ll_ptr, cw_ptr, n_units);

    // 5. Coverage-weighted warm start + OVR prior
    std::vector<double> prior(nc);
    std::vector<double> theta_init(nc);
    compute_ovr_prior_and_warm_start(
        ec_data, unique_totals.data(), pe_ptr,
        prior_alpha, prior_gamma,
        prior.data(), theta_init.data(), n_components);

    // 6. Run SQUAREM
    auto result = run_squarem(
        ec_data, log_eff_len.data(), unique_totals.data(),
        prior.data(), theta_init.data(),
        n_components, max_iterations, convergence_delta,
        use_vbem, prune_threshold, post_prune_iters);

    // 7. Return as numpy arrays
    auto* theta_out = new double[nc];
    auto* alpha_out = new double[nc];
    auto* em_out    = new double[nc];
    std::copy(result.theta.begin(), result.theta.end(), theta_out);
    std::copy(result.alpha.begin(), result.alpha.end(), alpha_out);
    std::copy(result.em_totals.begin(), result.em_totals.end(), em_out);

    size_t shape[1] = { nc };
    nb::capsule theta_owner(theta_out, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    nb::capsule alpha_owner(alpha_out, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    nb::capsule em_owner(em_out,       [](void* p) noexcept { delete[] static_cast<double*>(p); });

    return std::make_tuple(
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(theta_out, 1, shape, std::move(theta_owner)),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(alpha_out, 1, shape, std::move(alpha_owner)),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(em_out,    1, shape, std::move(em_owner))
    );
}

// ================================================================
// nanobind module definition
// ================================================================

NB_MODULE(_em_impl, m) {
    m.doc() = "C++ EM solver for hulkrna locus-level abundance estimation.\n\n"
              "Provides run_locus_em_native() which replaces the Python EM hot path:\n"
              "_em_step, _vbem_step, _build_equiv_classes, SQUAREM loop.";

    m.def("run_locus_em_native", &run_locus_em_native,
          nb::arg("offsets"),
          nb::arg("t_indices"),
          nb::arg("log_liks"),
          nb::arg("coverage_wts"),
          nb::arg("tx_starts"),
          nb::arg("tx_ends"),
          nb::arg("bias_profiles"),
          nb::arg("unique_totals"),
          nb::arg("effective_lens"),
          nb::arg("prior_eligible"),
          nb::arg("n_components"),
          nb::arg("prior_alpha"),
          nb::arg("prior_gamma"),
          nb::arg("max_iterations"),
          nb::arg("convergence_delta"),
          nb::arg("use_vbem"),
          nb::arg("prune_threshold"),
          nb::arg("post_prune_iters"),
          "Run EM for a single locus sub-problem.\n\n"
          "Takes CSR per-locus data + config, returns (theta, alpha, em_totals).\n"
          "Replaces the Python EM hot path with a single C++ call.");
}
