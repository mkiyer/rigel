/**
 * em_solver.cpp — C++ EM solver for rigel locus-level abundance estimation.
 *
 * Replaces the Python hot path: _em_step, _vbem_step, _build_equiv_classes,
 * _apply_bias_correction_uniform, and the SQUAREM acceleration loop.
 *
 * Module: rigel._em_impl
 *
 * Build:
 *   Part of the rigel scikit-build-core build — see CMakeLists.txt.
 *   Pure C++17 + nanobind + numpy.  No external dependencies.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <unordered_map>
#include <vector>

#include <atomic>
#include <thread>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/tuple.h>

namespace nb = nanobind;

// ================================================================
// Constants — single source of truth, exported to Python via module attrs.
// Python imports these from rigel._em_impl instead of redefining.
// ================================================================

static constexpr double EM_LOG_EPSILON = 1e-300;
static constexpr int    MAX_FRAG_LEN  = 1000000;
static constexpr int    SQUAREM_BUDGET_DIVISOR = 3;
static constexpr double NRNA_FRAC_CLAMP_EPS = 1e-8;

// ================================================================
// Array type aliases
// ================================================================

using i32_1d = nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig>;
using i64_1d = nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig>;
using f32_1d = nb::ndarray<const float,   nb::ndim<1>, nb::c_contig>;
using f64_1d = nb::ndarray<const double,  nb::ndim<1>, nb::c_contig>;
using u8_1d  = nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig>;

// Mutable variants for in-place modification
using f64_1d_mut = nb::ndarray<double, nb::ndim<1>, nb::c_contig>;
using f64_2d_mut = nb::ndarray<double, nb::ndim<2>, nb::c_contig>;
using f64_2d     = nb::ndarray<const double, nb::ndim<2>, nb::c_contig>;

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

struct EmEquivClass {
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

static std::vector<EmEquivClass> build_equiv_classes(
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
    std::vector<EmEquivClass> result;
    result.reserve(class_map.size());

    for (auto& [key, start_list] : class_map) {
        int k = static_cast<int>(key.size());
        int n = static_cast<int>(start_list.size());

        EmEquivClass ec;
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

    // ---- Deterministic ordering ----
    // Multi-threaded BAM scanning produces fragments in non-deterministic
    // order.  The unordered_map above inherits that non-determinism in both
    // (a) the iteration order of equiv classes, and (b) the row order of
    // units within each class.  Since the EM E-step accumulates column sums
    // over rows, and FP addition is non-associative, different row orders
    // produce ULP-level differences that SQUAREM amplifies across iterations
    // until post-EM pruning flips borderline decisions, causing large
    // cascading output differences.
    //
    // Fix: sort equiv classes by comp_idx, and sort rows within each class
    // by their log-likelihood fingerprint.  This makes the EM iteration
    // fully deterministic regardless of input fragment order.

    // Sort equiv classes by comp_idx (lexicographic)
    std::sort(result.begin(), result.end(),
              [](const EmEquivClass& a, const EmEquivClass& b) {
                  return a.comp_idx < b.comp_idx;
              });

    // Sort rows within each equiv class by log-lik (lexicographic)
    for (auto& ec : result) {
        if (ec.n <= 1) continue;
        int k = ec.k;
        int n = ec.n;

        // Build sort index
        std::vector<int> idx(static_cast<size_t>(n));
        std::iota(idx.begin(), idx.end(), 0);

        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            for (int j = 0; j < k; ++j) {
                double va = ec.ll_flat[static_cast<size_t>(a) * k + j];
                double vb = ec.ll_flat[static_cast<size_t>(b) * k + j];
                if (va != vb) return va < vb;
            }
            for (int j = 0; j < k; ++j) {
                double va = ec.wt_flat[static_cast<size_t>(a) * k + j];
                double vb = ec.wt_flat[static_cast<size_t>(b) * k + j];
                if (va != vb) return va < vb;
            }
            return false;
        });

        // Check if already sorted (common case for single-unit classes)
        bool already_sorted = true;
        for (int i = 0; i < n; ++i) {
            if (idx[i] != i) { already_sorted = false; break; }
        }
        if (already_sorted) continue;

        // Reorder ll_flat and wt_flat
        std::vector<double> new_ll(static_cast<size_t>(n) * k);
        std::vector<double> new_wt(static_cast<size_t>(n) * k);
        for (int i = 0; i < n; ++i) {
            auto src = static_cast<size_t>(idx[i]) * k;
            auto dst = static_cast<size_t>(i) * k;
            std::copy(ec.ll_flat.data() + src, ec.ll_flat.data() + src + k,
                      new_ll.data() + dst);
            std::copy(ec.wt_flat.data() + src, ec.wt_flat.data() + src + k,
                      new_wt.data() + dst);
        }
        ec.ll_flat = std::move(new_ll);
        ec.wt_flat = std::move(new_wt);
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
    const EmEquivClass& ec,
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
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,
    const double* unambig_totals,
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

    // M-step: theta_new = (unambig_totals + em_totals + prior), normalized
    double total = 0.0;
    for (int i = 0; i < n_components; ++i) {
        theta_new[i] = unambig_totals[i] + em_totals[i] + prior[i];
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
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,
    const double* unambig_totals,
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

    // M-step: alpha_new = unambig_totals + em_totals + prior (unnormalized)
    for (int i = 0; i < n_components; ++i) {
        alpha_new[i] = unambig_totals[i] + em_totals[i] + prior[i];
    }
}

// ================================================================
// Linked MAP-EM step: (theta_t, nrna_frac) → (theta_t_new, nrna_frac_new)
// ================================================================
//
// Uses the linked mRNA-nRNA model where:
//   log_weight[mRNA i]  = log(θ_t[i]) + log(1-nrna_frac[i]) - log_eff_len[i]
//   log_weight[nRNA i]  = log(θ_t[i]) + log(nrna_frac[i])   - log_eff_len[n_t+i]
//   log_weight[gDNA]    = log(γ)                     - log_eff_len[2*n_t]
//
// M-step: θ_t normalised jointly with γ; nrna_frac from MAP Beta.

static void linked_map_em_step(
    const double* theta_t,       // [n_t + 1]: θ_1..θ_T, γ
    const double* nrna_frac,           // [n_t]: nrna_frac_1..nrna_frac_T
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,   // [2*n_t + 1]
    const double* unambig_totals, // [2*n_t + 1]
    const double* prior,         // [2*n_t + 1]
    const double* nrna_frac_alpha,     // [n_t]
    const double* nrna_frac_beta,      // [n_t]
    double*       em_totals,     // [2*n_t + 1]: zeroed then accumulated
    double*       theta_t_new,   // [n_t + 1]: output (normalised)
    double*       nrna_frac_new,       // [n_t]: output
    int           n_transcripts,
    int           n_components)
{
    int n_t = n_transcripts;
    int gdna_idx = 2 * n_t;

    // Build log_weights from linked parameterisation
    std::vector<double> log_weights(static_cast<size_t>(n_components));
    for (int i = 0; i < n_t; ++i) {
        double log_theta = std::log(theta_t[i] + EM_LOG_EPSILON);
        // mRNA component i
        log_weights[i] = log_theta
            + std::log(1.0 - nrna_frac[i] + EM_LOG_EPSILON)
            - log_eff_len[i];
        // nRNA component n_t + i
        log_weights[n_t + i] = log_theta
            + std::log(nrna_frac[i] + EM_LOG_EPSILON)
            - log_eff_len[n_t + i];
    }
    // gDNA component
    log_weights[gdna_idx] = std::log(theta_t[n_t] + EM_LOG_EPSILON)
                           - log_eff_len[gdna_idx];

    // Zero em_totals
    std::fill(em_totals, em_totals + n_components, 0.0);

    // E-step: accumulate posteriors (kernel unchanged)
    for (const auto& ec : ec_data) {
        em_step_kernel(ec, log_weights.data(), em_totals);
    }

    // M-step: aggregate into theta_t and update nrna_frac
    double total = 0.0;
    for (int i = 0; i < n_t; ++i) {
        double mrna_count = unambig_totals[i] + em_totals[i];
        double nrna_count = unambig_totals[n_t + i] + em_totals[n_t + i];
        double total_count = mrna_count + nrna_count;

        // θ_t ∝ total_data + sum of mRNA & nRNA priors
        theta_t_new[i] = total_count + prior[i] + prior[n_t + i];
        total += theta_t_new[i];

        // nrna_frac MAP Beta: (nrna + α - 1) / (total + α + β - 2)
        double a = nrna_frac_alpha[i];
        double b = nrna_frac_beta[i];
        double nrna_frac_den = total_count + a + b - 2.0;
        if (nrna_frac_den > 0.0) {
            double nrna_frac_val = (nrna_count + a - 1.0) / nrna_frac_den;
            if (nrna_frac_val < NRNA_FRAC_CLAMP_EPS) nrna_frac_val = NRNA_FRAC_CLAMP_EPS;
            else if (nrna_frac_val > 1.0 - NRNA_FRAC_CLAMP_EPS) nrna_frac_val = 1.0 - NRNA_FRAC_CLAMP_EPS;
            nrna_frac_new[i] = nrna_frac_val;
        } else {
            nrna_frac_new[i] = nrna_frac[i];  // keep previous
        }
    }

    // gDNA
    theta_t_new[n_t] = unambig_totals[gdna_idx]
                     + em_totals[gdna_idx] + prior[gdna_idx];
    total += theta_t_new[n_t];

    // Normalise θ_t so Σ θ_t + γ = 1
    if (total > 0.0) {
        double inv = 1.0 / total;
        for (int i = 0; i <= n_t; ++i) {
            theta_t_new[i] *= inv;
        }
    }
}

// ================================================================
// Coverage-weighted warm start + OVR prior
// ================================================================

static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,    // [n_components] 1.0 if eligible, 0.0 otherwise
    double        alpha_flat,
    double        gamma,
    double*       prior_out,       // [n_components] output
    double*       theta_init_out,  // [n_components] output
    int           n_components)
{
    // Initialize theta_init from unambig_totals
    std::copy(unambig_totals, unambig_totals + n_components, theta_init_out);

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
    std::vector<double> nrna_frac;       // [n_transcripts], empty for classic mode
};

static EMResult run_squarem(
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,
    const double* unambig_totals,
    double*       prior,
    const double* theta_init,
    int           n_components,
    int           max_iterations,
    double        convergence_delta,
    bool          use_vbem,
    double        prune_threshold)
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
                      unambig_totals, prior, em_totals.data(),
                      state1.data(), n_components);

            vbem_step(state1.data(), ec_data, log_eff_len,
                      unambig_totals, prior, em_totals.data(),
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
                      unambig_totals, prior, em_totals.data(),
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
                        unambig_totals, prior, em_totals.data(),
                        state1.data(), n_components);

            map_em_step(state1.data(), ec_data, log_eff_len,
                        unambig_totals, prior, em_totals.data(),
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
                        unambig_totals, prior, em_totals.data(),
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

        // alpha_out = unambig_totals + em_totals + prior
        for (size_t i = 0; i < nc; ++i) {
            alpha_out[i] = unambig_totals[i] + em_totals[i] + prior[i];
        }
    }

    // ----------------------------------------------------------
    // Post-EM pruning
    // ----------------------------------------------------------
    if (prune_threshold >= 0.0) {
        bool any_pruned = false;
        std::vector<bool> prune_mask(nc, false);

        for (size_t i = 0; i < nc; ++i) {
            double data_count = unambig_totals[i] + em_totals[i];
            double denom = std::max(alpha_out[i], EM_LOG_EPSILON);
            double evidence_ratio = data_count / denom;
            if (unambig_totals[i] == 0.0
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
                {
                    std::vector<double> a_new(nc);
                    vbem_step(alpha_out.data(), ec_data, log_eff_len,
                              unambig_totals, prior, em_totals.data(),
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
                // Single post-prune EM step to redistribute mass
                {
                    std::vector<double> t_new(nc);
                    map_em_step(theta.data(), ec_data, log_eff_len,
                                unambig_totals, prior, em_totals.data(),
                                t_new.data(), n_components);
                    std::swap(theta, t_new);
                }
                for (size_t i = 0; i < nc; ++i) {
                    alpha_out[i] = unambig_totals[i] + em_totals[i] + prior[i];
                }
            }
        }
    }

    return { std::move(theta), std::move(alpha_out), std::move(em_totals), {} };
}

// ================================================================
// Linked SQUAREM acceleration
// ================================================================
//
// SQUAREM on theta_t (n_t+1 elements).  nrna_frac is derived from the
// M-step after each EM iteration — NOT extrapolated.

static EMResult linked_run_squarem(
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,     // [2*n_t + 1]
    const double* unambig_totals,   // [2*n_t + 1]
    double*       prior,           // [2*n_t + 1], mutable for pruning
    const double* theta_t_init,    // [n_t + 1]
    const double* nrna_frac_init,        // [n_t]
    const double* nrna_frac_alpha,       // [n_t]
    const double* nrna_frac_beta,        // [n_t]
    int           n_transcripts,
    int           n_components,    // 2*n_t + 1
    int           max_iterations,
    double        convergence_delta,
    double        prune_threshold)
{
    int n_t = n_transcripts;
    int nc  = n_components;
    int ns  = n_t + 1;  // SQUAREM state size (theta_t + gamma)
    int max_sq_iters = std::max(max_iterations / SQUAREM_BUDGET_DIVISOR, 1);

    std::vector<double> em_totals(static_cast<size_t>(nc), 0.0);
    std::vector<double> nrna_frac(nrna_frac_init, nrna_frac_init + n_t);
    std::vector<double> nrna_frac_tmp(static_cast<size_t>(n_t));

    // SQUAREM state vectors (theta_t space)
    std::vector<double> state0(static_cast<size_t>(ns));
    std::vector<double> state1(static_cast<size_t>(ns));
    std::vector<double> state2(static_cast<size_t>(ns));
    std::vector<double> state_extrap(static_cast<size_t>(ns));
    std::vector<double> state_new(static_cast<size_t>(ns));
    std::vector<double> r_vec(static_cast<size_t>(ns));
    std::vector<double> v_vec(static_cast<size_t>(ns));

    // Initialize state0 = theta_t_init + collapsed prior, normalised
    for (int i = 0; i < n_t; ++i) {
        state0[i] = theta_t_init[i] + prior[i] + prior[n_t + i];
    }
    state0[n_t] = theta_t_init[n_t] + prior[2 * n_t];
    {
        double s = 0.0;
        for (int i = 0; i < ns; ++i) s += state0[i];
        if (s > 0.0) {
            double inv = 1.0 / s;
            for (int i = 0; i < ns; ++i) state0[i] *= inv;
        }
    }

    for (int iter = 0; iter < max_sq_iters; ++iter) {
        // Two linked EM steps
        linked_map_em_step(
            state0.data(), nrna_frac.data(), ec_data,
            log_eff_len, unambig_totals, prior,
            nrna_frac_alpha, nrna_frac_beta,
            em_totals.data(), state1.data(), nrna_frac_tmp.data(),
            n_transcripts, n_components);
        std::copy(nrna_frac_tmp.begin(), nrna_frac_tmp.end(), nrna_frac.begin());

        linked_map_em_step(
            state1.data(), nrna_frac.data(), ec_data,
            log_eff_len, unambig_totals, prior,
            nrna_frac_alpha, nrna_frac_beta,
            em_totals.data(), state2.data(), nrna_frac_tmp.data(),
            n_transcripts, n_components);
        std::copy(nrna_frac_tmp.begin(), nrna_frac_tmp.end(), nrna_frac.begin());

        // SQUAREM extrapolation on theta_t
        double sv2 = 0.0, srv = 0.0;
        for (int i = 0; i < ns; ++i) {
            r_vec[i] = state1[i] - state0[i];
            v_vec[i] = (state2[i] - state1[i]) - r_vec[i];
            sv2 += v_vec[i] * v_vec[i];
            srv += r_vec[i] * v_vec[i];
        }

        if (sv2 == 0.0) {
            std::copy(state2.begin(), state2.end(), state_extrap.begin());
        } else {
            double alpha_step = std::max(-srv / sv2, 1.0);
            for (int i = 0; i < ns; ++i) {
                state_extrap[i] = state0[i]
                    + 2.0 * alpha_step * r_vec[i]
                    + alpha_step * alpha_step * v_vec[i];
                if (state_extrap[i] < 0.0) state_extrap[i] = 0.0;
            }
            double s = 0.0;
            for (int i = 0; i < ns; ++i) s += state_extrap[i];
            if (s > 0.0) {
                double inv = 1.0 / s;
                for (int i = 0; i < ns; ++i) state_extrap[i] *= inv;
            } else {
                std::copy(state2.begin(), state2.end(),
                          state_extrap.begin());
            }
        }

        // Stabilisation step
        linked_map_em_step(
            state_extrap.data(), nrna_frac.data(), ec_data,
            log_eff_len, unambig_totals, prior,
            nrna_frac_alpha, nrna_frac_beta,
            em_totals.data(), state_new.data(), nrna_frac_tmp.data(),
            n_transcripts, n_components);
        std::copy(nrna_frac_tmp.begin(), nrna_frac_tmp.end(), nrna_frac.begin());

        // Convergence check
        double delta = 0.0;
        for (int i = 0; i < ns; ++i) {
            delta += std::abs(state_new[i] - state0[i]);
        }
        std::swap(state0, state_new);
        if (delta < convergence_delta) break;
    }

    // Decompose theta_t + nrna_frac → full theta[2*n_t+1]
    std::vector<double> theta(static_cast<size_t>(nc));
    for (int i = 0; i < n_t; ++i) {
        theta[i]       = state0[i] * (1.0 - nrna_frac[i]);
        theta[n_t + i] = state0[i] * nrna_frac[i];
    }
    theta[2 * n_t] = state0[n_t];

    // alpha_out = unambig_totals + em_totals + prior
    std::vector<double> alpha_out(static_cast<size_t>(nc));
    for (int i = 0; i < nc; ++i) {
        alpha_out[i] = unambig_totals[i] + em_totals[i] + prior[i];
    }

    // ----------------------------------------------------------
    // Post-EM pruning at transcript level
    // ----------------------------------------------------------
    if (prune_threshold >= 0.0) {
        bool any_pruned = false;
        std::vector<bool> prune_mask(static_cast<size_t>(n_t), false);

        for (int i = 0; i < n_t; ++i) {
            double mrna_data = unambig_totals[i] + em_totals[i];
            double nrna_data = unambig_totals[n_t + i] + em_totals[n_t + i];
            double total_data = mrna_data + nrna_data;
            double total_alpha = alpha_out[i] + alpha_out[n_t + i];
            double evidence = total_data
                / std::max(total_alpha, EM_LOG_EPSILON);
            if (unambig_totals[i] == 0.0
                && evidence < prune_threshold) {
                prune_mask[i] = true;
                any_pruned = true;
            }
        }

        if (any_pruned) {
            // Zero prior and state for pruned transcripts
            for (int i = 0; i < n_t; ++i) {
                if (prune_mask[i]) {
                    prior[i] = 0.0;
                    prior[n_t + i] = 0.0;
                    state0[i] = 0.0;
                }
            }
            // Re-normalise state
            double s = 0.0;
            for (int i = 0; i < ns; ++i) s += state0[i];
            if (s > 0.0) {
                double inv = 1.0 / s;
                for (int i = 0; i < ns; ++i) state0[i] *= inv;
            }

            // Post-prune EM iteration (single step)
            {
                std::vector<double> t_new(static_cast<size_t>(ns));
                linked_map_em_step(
                    state0.data(), nrna_frac.data(), ec_data,
                    log_eff_len, unambig_totals, prior,
                    nrna_frac_alpha, nrna_frac_beta,
                    em_totals.data(), t_new.data(), nrna_frac_tmp.data(),
                    n_transcripts, n_components);
                std::swap(state0, t_new);
                std::copy(nrna_frac_tmp.begin(), nrna_frac_tmp.end(), nrna_frac.begin());
            }

            // Rebuild full theta + alpha
            for (int i = 0; i < n_t; ++i) {
                theta[i]       = state0[i] * (1.0 - nrna_frac[i]);
                theta[n_t + i] = state0[i] * nrna_frac[i];
            }
            theta[2 * n_t] = state0[n_t];

            for (int i = 0; i < nc; ++i) {
                alpha_out[i] = unambig_totals[i] + em_totals[i] + prior[i];
            }
        }
    }

    return { std::move(theta), std::move(alpha_out),
             std::move(em_totals), std::move(nrna_frac) };
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
    f64_1d unambig_totals_arr,
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
    // Linked model parameters
    int    n_transcripts,
    f64_1d nrna_frac_alpha_arr,
    f64_1d nrna_frac_beta_arr)
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
    const double*   ut_ptr   = unambig_totals_arr.data();
    const double*   el_ptr   = effective_lens.data();
    const double*   pe_ptr   = prior_eligible.data();

    bool linked = (n_transcripts > 0);
    size_t n_t = linked ? static_cast<size_t>(n_transcripts) : 0;
    const double* ea_ptr = linked ? nrna_frac_alpha_arr.data() : nullptr;
    const double* eb_ptr = linked ? nrna_frac_beta_arr.data() : nullptr;

    // Copy unambig_totals (we need a mutable copy)
    std::vector<double> unambig_totals(ut_ptr, ut_ptr + nc);

    // 1. Apply bias correction (uniform fast-path)
    if (n_candidates > 0) {
        apply_bias_correction_uniform(
            ll_ptr, ti_ptr, txs_ptr, txe_ptr, bp_ptr, n_candidates);
    }

    // 2. Handle empty locus
    if (n_units == 0 || n_candidates == 0) {
        // theta = (unambig_totals + prior) / sum, alpha = unambig_totals + prior
        std::vector<double> alpha(nc);
        double total = 0.0;
        for (size_t i = 0; i < nc; ++i) {
            double p = (pe_ptr[i] > 0.0) ? prior_alpha : 0.0;
            alpha[i] = unambig_totals[i] + p;
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

        // Compute nrna_frac (prior mean) for linked model
        size_t nrna_frac_alloc = (n_t > 0) ? n_t : 1;
        auto* nrna_frac_data = new double[nrna_frac_alloc];
        if (linked) {
            for (size_t i = 0; i < n_t; ++i) {
                double sum_ab = ea_ptr[i] + eb_ptr[i];
                nrna_frac_data[i] = (sum_ab > 0.0) ? ea_ptr[i] / sum_ab : 0.5;
            }
        }

        // Return as numpy arrays
        auto* theta_data = new double[nc];
        auto* alpha_data = new double[nc];
        auto* em_data    = new double[nc];
        std::copy(theta.begin(), theta.end(), theta_data);
        std::copy(alpha.begin(), alpha.end(), alpha_data);
        std::copy(em_totals.begin(), em_totals.end(), em_data);

        size_t shape[1]     = { nc };
        size_t shape_nrna_frac[1] = { n_t };
        nb::capsule theta_owner(theta_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
        nb::capsule alpha_owner(alpha_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
        nb::capsule em_owner(em_data,       [](void* p) noexcept { delete[] static_cast<double*>(p); });
        nb::capsule nrna_frac_owner(nrna_frac_data,     [](void* p) noexcept { delete[] static_cast<double*>(p); });

        return std::make_tuple(
            nb::ndarray<nb::numpy, double, nb::ndim<1>>(theta_data, 1, shape, std::move(theta_owner)),
            nb::ndarray<nb::numpy, double, nb::ndim<1>>(alpha_data, 1, shape, std::move(alpha_owner)),
            nb::ndarray<nb::numpy, double, nb::ndim<1>>(em_data,    1, shape, std::move(em_owner)),
            nb::ndarray<nb::numpy, double, nb::ndim<1>>(nrna_frac_data,   1, shape_nrna_frac, std::move(nrna_frac_owner))
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
        ec_data, unambig_totals.data(), pe_ptr,
        prior_alpha, prior_gamma,
        prior.data(), theta_init.data(), n_components);

    // 6. Run SQUAREM (linked or classic)
    EMResult result;
    if (linked) {
        // Collapse warm start to theta_t space
        int nt_i = n_transcripts;
        std::vector<double> theta_t_init(static_cast<size_t>(nt_i + 1));
        for (int i = 0; i < nt_i; ++i) {
            theta_t_init[i] = theta_init[i] + theta_init[nt_i + i];
        }
        theta_t_init[nt_i] = theta_init[2 * nt_i];

        // nrna_frac_init = prior mean, clamped
        std::vector<double> nrna_frac_init(n_t);
        for (size_t i = 0; i < n_t; ++i) {
            double sum_ab = ea_ptr[i] + eb_ptr[i];
            nrna_frac_init[i] = (sum_ab > 0.0) ? ea_ptr[i] / sum_ab : 0.5;
            if (nrna_frac_init[i] < NRNA_FRAC_CLAMP_EPS) nrna_frac_init[i] = NRNA_FRAC_CLAMP_EPS;
            else if (nrna_frac_init[i] > 1.0 - NRNA_FRAC_CLAMP_EPS) nrna_frac_init[i] = 1.0 - NRNA_FRAC_CLAMP_EPS;
        }

        result = linked_run_squarem(
            ec_data, log_eff_len.data(), unambig_totals.data(),
            prior.data(), theta_t_init.data(), nrna_frac_init.data(),
            ea_ptr, eb_ptr,
            n_transcripts, n_components,
            max_iterations, convergence_delta,
            prune_threshold);
    } else {
        result = run_squarem(
            ec_data, log_eff_len.data(), unambig_totals.data(),
            prior.data(), theta_init.data(),
            n_components, max_iterations, convergence_delta,
            use_vbem, prune_threshold);
    }

    // 7. Return as numpy arrays
    size_t nrna_frac_alloc = (n_t > 0) ? n_t : 1;
    auto* theta_out = new double[nc];
    auto* alpha_out = new double[nc];
    auto* em_out    = new double[nc];
    auto* nrna_frac_out   = new double[nrna_frac_alloc];
    std::copy(result.theta.begin(), result.theta.end(), theta_out);
    std::copy(result.alpha.begin(), result.alpha.end(), alpha_out);
    std::copy(result.em_totals.begin(), result.em_totals.end(), em_out);
    if (n_t > 0) {
        std::copy(result.nrna_frac.begin(), result.nrna_frac.end(), nrna_frac_out);
    }

    size_t shape[1]     = { nc };
    size_t shape_nrna_frac[1] = { n_t };
    nb::capsule theta_owner(theta_out, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    nb::capsule alpha_owner(alpha_out, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    nb::capsule em_owner(em_out,       [](void* p) noexcept { delete[] static_cast<double*>(p); });
    nb::capsule nrna_frac_owner(nrna_frac_out,     [](void* p) noexcept { delete[] static_cast<double*>(p); });

    return std::make_tuple(
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(theta_out, 1, shape, std::move(theta_owner)),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(alpha_out, 1, shape, std::move(alpha_owner)),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(em_out,    1, shape, std::move(em_owner)),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(nrna_frac_out,   1, shape_nrna_frac, std::move(nrna_frac_owner))
    );
}

// ================================================================
// Batch locus EM — Phase 2A: single C++ call for all loci
// ================================================================
//
// Replaces the Python per-locus for-loop:
//   for locus in loci:
//       build_locus_em_data → run_locus_em → assign_locus_ambiguous
//
// Processes all loci in a single C++ call, eliminating 29K Python→C++
// round-trips and all numpy/pandas per-locus overhead.

// Numerical-stability floor for prior components (single source — also
// exported via module attr; see NB_MODULE block at bottom).
static constexpr double EM_PRIOR_EPSILON = 1e-10;

// Per-locus candidate record (used during sub-problem extraction)
struct LocalCandidate {
    int32_t local_comp;
    double  log_lik;
    double  cov_wt;
    int32_t tx_start;
    int32_t tx_end;
    uint8_t count_col;
};

// Per-locus sub-problem (stack-allocated, reused across loci)
struct LocusSubProblem {
    int n_t;           // number of transcripts in locus
    int n_components;  // 2*n_t + 1
    int gdna_idx;      // = 2*n_t
    int n_local_units;

    // Local CSR
    std::vector<int64_t>  offsets;      // [n_local_units + 1]
    std::vector<int32_t>  t_indices;    // local component indices
    std::vector<double>   log_liks;
    std::vector<double>   coverage_wts;
    std::vector<int32_t>  tx_starts;
    std::vector<int32_t>  tx_ends;
    std::vector<uint8_t>  count_cols;

    // Per-unit metadata
    std::vector<int32_t>  locus_t_arr;   // best transcript (global) per unit
    std::vector<uint8_t>  locus_ct_arr;  // count col for best transcript

    // Per-component
    std::vector<double>   unambig_totals;   // [n_components]
    std::vector<double>   prior;            // [n_components]
    std::vector<int64_t>  bias_profiles;    // [n_components]
    std::vector<double>   nrna_frac_alpha;  // [n_t]
    std::vector<double>   nrna_frac_beta;   // [n_t]
    std::vector<double>   eligible;         // [n_components]

    // Local→global transcript mapping
    std::vector<int32_t>  local_to_global_t; // [n_t]
};

// Extract per-locus sub-problem from global CSR data.
// Reimplements Python build_locus_em_data() entirely in C++.
static void extract_locus_sub_problem(
    LocusSubProblem& sub,
    // Locus definition
    const int32_t* t_arr, int n_t,
    const int32_t* u_arr, int n_u,
    double gdna_init,
    // Global CSR
    const int64_t* g_offsets,
    const int32_t* g_t_indices,
    const double*  g_log_liks,
    const double*  g_coverage_wts,
    const int32_t* g_tx_starts,
    const int32_t* g_tx_ends,
    const uint8_t* g_count_cols,
    const uint8_t* g_is_spliced,
    const double*  g_gdna_log_liks,
    const int32_t* g_genomic_footprints,
    const int32_t* g_locus_t_indices,
    const uint8_t* g_locus_count_cols,
    int32_t nrna_base,
    // Per-transcript data (global arrays, all N_T_TOTAL elements)
    const double*  all_unambig_row_sums,
    const double*  all_nrna_init,
    const double*  all_nrna_frac_alpha,
    const double*  all_nrna_frac_beta,
    const int64_t* all_t_starts,
    const int64_t* all_t_ends,
    const int64_t* all_t_lengths,
    const double*  all_transcript_spans,
    const double*  all_exonic_lengths,
    double mean_frag,
    // Scratch: global→local mapping, caller-owned, reset between calls
    int32_t* local_map, int local_map_size)
{
    sub.n_t = n_t;
    sub.n_components = 2 * n_t + 1;
    sub.gdna_idx = 2 * n_t;
    sub.n_local_units = n_u;

    // --- Build global→local mapping ---
    int max_global = 0;
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        int nrna_gt = nrna_base + gt;
        if (gt + 1 > max_global) max_global = gt + 1;
        if (nrna_gt + 1 > max_global) max_global = nrna_gt + 1;
    }
    // Ensure we don't exceed scratch buffer
    if (max_global > local_map_size) max_global = local_map_size;

    // Reset scratch region to -1
    for (int i = 0; i < max_global; ++i) local_map[i] = -1;

    // mRNA: global_t → local_i
    // nRNA: nrna_base + global_t → n_t + local_i
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        if (gt >= 0 && gt < local_map_size) local_map[gt] = i;
        int ngt = nrna_base + gt;
        if (ngt >= 0 && ngt < local_map_size) local_map[ngt] = n_t + i;
    }

    // --- Gather + dedup candidates from global CSR ---
    // For each unit, extract candidates, map to local, dedup by
    // (unit, local_comp) keeping best log_lik.

    int nc = sub.n_components;

    // Temporary per-unit dedup buffers (indexed by local_comp)
    struct BestCandidate {
        double  log_lik;
        double  cov_wt;
        int32_t tx_start;
        int32_t tx_end;
        uint8_t count_col;
    };
    std::vector<BestCandidate> best_buf(nc);
    std::vector<bool> seen(nc, false);

    // Pre-size output vectors (estimate: total global candidates for these units)
    size_t total_est = 0;
    for (int ui = 0; ui < n_u; ++ui) {
        int u = u_arr[ui];
        total_est += static_cast<size_t>(g_offsets[u + 1] - g_offsets[u]);
    }
    // Add space for gDNA candidates (at most n_u)
    total_est += static_cast<size_t>(n_u);

    sub.t_indices.clear();
    sub.log_liks.clear();
    sub.coverage_wts.clear();
    sub.tx_starts.clear();
    sub.tx_ends.clear();
    sub.count_cols.clear();
    sub.t_indices.reserve(total_est);
    sub.log_liks.reserve(total_est);
    sub.coverage_wts.reserve(total_est);
    sub.tx_starts.reserve(total_est);
    sub.tx_ends.reserve(total_est);
    sub.count_cols.reserve(total_est);

    sub.offsets.resize(n_u + 1);
    sub.offsets[0] = 0;

    sub.locus_t_arr.resize(n_u);
    sub.locus_ct_arr.resize(n_u);

    // Compute locus span for gDNA bias profile
    int64_t locus_start = all_t_starts[t_arr[0]];
    int64_t locus_end = all_t_ends[t_arr[0]];
    for (int i = 1; i < n_t; ++i) {
        int64_t ts = all_t_starts[t_arr[i]];
        int64_t te = all_t_ends[t_arr[i]];
        if (ts < locus_start) locus_start = ts;
        if (te > locus_end) locus_end = te;
    }
    int64_t locus_span = locus_end - locus_start;

    for (int ui = 0; ui < n_u; ++ui) {
        int u = u_arr[ui];
        auto g_start = g_offsets[u];
        auto g_end = g_offsets[u + 1];

        // Per-unit metadata
        sub.locus_t_arr[ui] = g_locus_t_indices[u];
        sub.locus_ct_arr[ui] = g_locus_count_cols[u];

        // Reset seen flags for this unit
        for (int c = 0; c < nc; ++c) seen[c] = false;

        // Gather and dedup RNA candidates
        for (auto j = g_start; j < g_end; ++j) {
            int32_t global_t = g_t_indices[j];
            if (global_t < 0 || global_t >= local_map_size) continue;
            int32_t local = local_map[global_t];
            if (local < 0 || local >= nc) continue;

            double ll = g_log_liks[j];
            if (!seen[local] || ll > best_buf[local].log_lik) {
                best_buf[local] = {ll, g_coverage_wts[j],
                                   g_tx_starts[j], g_tx_ends[j],
                                   g_count_cols[j]};
                seen[local] = true;
            }
        }

        // Add gDNA candidate for unspliced units
        bool is_spliced = (g_is_spliced[u] != 0);
        double gdna_ll = g_gdna_log_liks[u];
        if (!is_spliced && std::isfinite(gdna_ll)) {
            int32_t gdna_comp = sub.gdna_idx;
            int32_t footprint = g_genomic_footprints[u];
            if (!seen[gdna_comp] || gdna_ll > best_buf[gdna_comp].log_lik) {
                best_buf[gdna_comp] = {gdna_ll, 1.0, 0, footprint, 0};
                seen[gdna_comp] = true;
            }
        }

        // Collect deduped candidates for this unit (sorted by local_comp
        // for stable equivalence class keys)
        for (int c = 0; c < nc; ++c) {
            if (seen[c]) {
                sub.t_indices.push_back(c);
                sub.log_liks.push_back(best_buf[c].log_lik);
                sub.coverage_wts.push_back(best_buf[c].cov_wt);
                sub.tx_starts.push_back(best_buf[c].tx_start);
                sub.tx_ends.push_back(best_buf[c].tx_end);
                sub.count_cols.push_back(best_buf[c].count_col);
            }
        }

        sub.offsets[ui + 1] = static_cast<int64_t>(sub.t_indices.size());
    }

    // --- Build per-component arrays ---
    sub.local_to_global_t.resize(n_t);
    for (int i = 0; i < n_t; ++i) sub.local_to_global_t[i] = t_arr[i];

    // Unambig totals: mRNA slot = row sum of unambig_counts[t_arr[i], :]
    sub.unambig_totals.assign(nc, 0.0);
    for (int i = 0; i < n_t; ++i) {
        sub.unambig_totals[i] = all_unambig_row_sums[t_arr[i]];
    }
    // nRNA and gDNA slots are 0 (no double-counting)

    // Bias profiles: [mRNA lengths, nRNA lengths, locus_span]
    sub.bias_profiles.resize(nc);
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        sub.bias_profiles[i] = all_t_lengths[gt];           // mRNA
        sub.bias_profiles[n_t + i] = all_t_ends[gt] - all_t_starts[gt]; // nRNA
    }
    sub.bias_profiles[sub.gdna_idx] = locus_span;

    // Prior: EM_PRIOR_EPSILON for eligible, 0.0 for ineligible
    sub.prior.assign(nc, EM_PRIOR_EPSILON);

    // Zero gDNA prior when gdna_init == 0
    if (gdna_init == 0.0) {
        sub.prior[sub.gdna_idx] = 0.0;
    }

    // Zero nRNA prior for single-exon transcripts
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        double span = all_transcript_spans[gt];
        double exon_len = all_exonic_lengths[gt];
        if (span <= exon_len) {
            sub.prior[n_t + i] = 0.0;  // nRNA slot
        }
    }

    // Zero nRNA prior when nrna_init is zero
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        if (all_nrna_init[gt] == 0.0) {
            sub.prior[n_t + i] = 0.0;
        }
    }

    // Zero unambig_totals for dead components
    for (int c = 0; c < nc; ++c) {
        if (sub.prior[c] == 0.0) sub.unambig_totals[c] = 0.0;
    }

    // Eligible = (prior > 0)
    sub.eligible.resize(nc);
    for (int c = 0; c < nc; ++c) {
        sub.eligible[c] = (sub.prior[c] > 0.0) ? 1.0 : 0.0;
    }

    // nrna_frac alpha/beta
    sub.nrna_frac_alpha.resize(n_t);
    sub.nrna_frac_beta.resize(n_t);
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        sub.nrna_frac_alpha[i] = all_nrna_frac_alpha[gt];
        sub.nrna_frac_beta[i] = all_nrna_frac_beta[gt];
    }

    // Clean up local_map scratch for next call
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        if (gt >= 0 && gt < local_map_size) local_map[gt] = -1;
        int ngt = nrna_base + gt;
        if (ngt >= 0 && ngt < local_map_size) local_map[ngt] = -1;
    }
}

// Assign posteriors after EM convergence.
// Reimplements Python assign_locus_ambiguous() entirely in C++.
// Scatters results into the provided accumulator arrays.
static void assign_posteriors(
    const LocusSubProblem& sub,
    const double* theta,
    double confidence_threshold,
    // Output accumulators (accumulated across loci)
    double* em_counts_2d,          // [N_T, n_cols], row-major
    double* em_high_conf_2d,       // [N_T, n_cols]
    double* nrna_em_counts,        // [N_T]
    double* gdna_locus_counts_2d,  // [N_T, n_cols]
    double* posterior_sum,         // [N_T]
    double* n_assigned,            // [N_T]
    // Per-locus accumulation
    double& mrna_total,
    double& nrna_total,
    double& gdna_total,
    int N_T_TOTAL,  // total transcripts for bounds checking
    int n_cols)     // number of splice-strand columns (actual 2D stride)
{
    int n_t = sub.n_t;
    int nc  = sub.n_components;
    int gdna_idx = sub.gdna_idx;
    int n_units = sub.n_local_units;
    const int32_t* local_to_global = sub.local_to_global_t.data();

    // Effective lengths are all 1.0, so log_eff_len = 0.
    // log_weights = log(theta + eps)
    std::vector<double> log_weights(nc);
    for (int c = 0; c < nc; ++c) {
        log_weights[c] = std::log(theta[c] + EM_LOG_EPSILON);
    }

    mrna_total = 0.0;
    nrna_total = 0.0;
    gdna_total = 0.0;

    // Process each unit
    for (int ui = 0; ui < n_units; ++ui) {
        auto s = sub.offsets[ui];
        auto e = sub.offsets[ui + 1];
        int seg_len = static_cast<int>(e - s);
        if (seg_len == 0) continue;

        // Compute log posteriors
        // log_posterior[j] = log_lik[j] + log_weights[t_indices[j]]
        double max_val = -1e300;
        for (int j = 0; j < seg_len; ++j) {
            int32_t comp = sub.t_indices[s + j];
            double lp = sub.log_liks[s + j] + log_weights[comp];
            if (lp > max_val) max_val = lp;
        }

        // Log-sum-exp normalization
        double sum_exp = 0.0;
        std::vector<double> posteriors(seg_len);
        for (int j = 0; j < seg_len; ++j) {
            int32_t comp = sub.t_indices[s + j];
            double lp = sub.log_liks[s + j] + log_weights[comp];
            posteriors[j] = std::exp(lp - max_val);
            sum_exp += posteriors[j];
        }
        if (sum_exp > 0.0 && std::isfinite(sum_exp)) {
            double inv = 1.0 / sum_exp;
            for (int j = 0; j < seg_len; ++j) posteriors[j] *= inv;
        } else {
            for (int j = 0; j < seg_len; ++j) posteriors[j] = 0.0;
        }
        // Track max mRNA posterior for high-confidence gating
        double max_mrna_posterior = 0.0;

        // Scatter posteriors
        for (int j = 0; j < seg_len; ++j) {
            int32_t comp = sub.t_indices[s + j];
            double p = posteriors[j];
            if (p == 0.0) continue;

            if (comp < n_t) {
                // mRNA
                int32_t global_t = local_to_global[comp];
                uint8_t col = sub.count_cols[s + j];
                if (global_t < 0 || global_t >= N_T_TOTAL || col >= n_cols) continue;
                em_counts_2d[global_t * n_cols + col] += p;
                mrna_total += p;
                if (p > max_mrna_posterior) max_mrna_posterior = p;

                // Confidence tracking
                posterior_sum[global_t] += p * p;
                n_assigned[global_t] += p;
            } else if (comp < gdna_idx) {
                // nRNA
                int32_t local_t = comp - n_t;
                int32_t global_t = local_to_global[local_t];
                nrna_em_counts[global_t] += p;
                nrna_total += p;
            } else {
                // gDNA
                gdna_total += p;
            }
        }

        // High-confidence mRNA: if max mRNA posterior >= threshold,
        // scatter all mRNA posteriors in this unit to hc counts
        if (max_mrna_posterior >= confidence_threshold) {
            for (int j = 0; j < seg_len; ++j) {
                int32_t comp = sub.t_indices[s + j];
                double p = posteriors[j];
                if (comp < n_t && p > 0.0) {
                    int32_t global_t = local_to_global[comp];
                    uint8_t col = sub.count_cols[s + j];
                    if (global_t >= 0 && global_t < N_T_TOTAL && col < n_cols) {
                        em_high_conf_2d[global_t * n_cols + col] += p;
                    }
                }
            }
        }

        // gDNA locus attribution
        double gdna_unit_sum = 0.0;
        for (int j = 0; j < seg_len; ++j) {
            if (sub.t_indices[s + j] == gdna_idx) {
                gdna_unit_sum += posteriors[j];
            }
        }
        if (gdna_unit_sum > 0.0) {
            int32_t lt = sub.locus_t_arr[ui];
            uint8_t lct = sub.locus_ct_arr[ui];
            if (lt >= 0 && lt < N_T_TOTAL && lct < n_cols) {
                gdna_locus_counts_2d[lt * n_cols + lct] += gdna_unit_sum;
            }
        }
    }
}

// ----------------------------------------------------------------
// Top-level batch entry point
// ----------------------------------------------------------------

static std::tuple<
    double,  // total_gdna_em
    nb::ndarray<nb::numpy, double, nb::ndim<1>>,  // locus_mrna[n_loci]
    nb::ndarray<nb::numpy, double, nb::ndim<1>>,  // locus_nrna[n_loci]
    nb::ndarray<nb::numpy, double, nb::ndim<1>>   // locus_gdna[n_loci]
>
batch_locus_em(
    // --- Global CSR (ScoredFragments) ---
    i64_1d g_offsets,
    i32_1d g_t_indices,
    f64_1d g_log_liks,
    f64_1d g_coverage_wts,
    i32_1d g_tx_starts,
    i32_1d g_tx_ends,
    u8_1d  g_count_cols,
    u8_1d  g_is_spliced,
    f64_1d g_gdna_log_liks,
    i32_1d g_genomic_footprints,
    i32_1d g_locus_t_indices,
    u8_1d  g_locus_count_cols,
    int32_t nrna_base_index,
    // --- Locus definitions (flattened CSR) ---
    i64_1d locus_t_offsets,    // [n_loci + 1]
    i32_1d locus_t_flat,       // concatenated transcript indices
    i64_1d locus_u_offsets,    // [n_loci + 1]
    i32_1d locus_u_flat,       // concatenated unit indices
    f64_1d gdna_inits,         // [n_loci]
    // --- Per-transcript data ---
    f64_2d unambig_counts,     // [N_T, n_cols]
    f64_1d nrna_init_arr,
    f64_1d nrna_frac_alpha_arr,
    f64_1d nrna_frac_beta_arr,
    i64_1d t_starts_arr,
    i64_1d t_ends_arr,
    i64_1d t_lengths_arr,
    f64_1d transcript_spans_arr,
    f64_1d exonic_lengths_arr,
    // --- Mutable output accumulators ---
    f64_2d_mut em_counts_out,          // [N_T, n_cols]
    f64_2d_mut em_high_conf_out,       // [N_T, n_cols]
    f64_1d_mut nrna_em_counts_out,     // [N_T]
    f64_2d_mut gdna_locus_counts_out,  // [N_T, n_cols]
    f64_1d_mut posterior_sum_out,       // [N_T]
    f64_1d_mut n_assigned_out,         // [N_T]
    // --- EM config ---
    double mean_frag,
    int    max_iterations,
    double convergence_delta,
    double prior_alpha,
    double prior_gamma,
    bool   use_vbem,
    double prune_threshold,
    double confidence_threshold,
    int    n_transcripts_total,
    int    n_splice_strand_cols,
    int    n_threads)
{
    int n_loci = static_cast<int>(locus_t_offsets.shape(0)) - 1;
    int N_T = n_transcripts_total;
    int N_COLS = n_splice_strand_cols;

    // Get raw pointers to all input arrays
    const int64_t*  goff  = g_offsets.data();
    const int32_t*  gti   = g_t_indices.data();
    const double*   gll   = g_log_liks.data();
    const double*   gcw   = g_coverage_wts.data();
    const int32_t*  gtxs  = g_tx_starts.data();
    const int32_t*  gtxe  = g_tx_ends.data();
    const uint8_t*  gcc   = g_count_cols.data();
    const uint8_t*  gis   = g_is_spliced.data();
    const double*   ggll  = g_gdna_log_liks.data();
    const int32_t*  ggfp  = g_genomic_footprints.data();
    const int32_t*  glti  = g_locus_t_indices.data();
    const uint8_t*  glcc  = g_locus_count_cols.data();

    const int64_t*  lt_off = locus_t_offsets.data();
    const int32_t*  lt_fl  = locus_t_flat.data();
    const int64_t*  lu_off = locus_u_offsets.data();
    const int32_t*  lu_fl  = locus_u_flat.data();
    const double*   gi_ptr = gdna_inits.data();

    const double*   uac    = unambig_counts.data();  // [N_T, N_COLS] row-major
    const double*   nri    = nrna_init_arr.data();
    const double*   nfa    = nrna_frac_alpha_arr.data();
    const double*   nfb    = nrna_frac_beta_arr.data();
    const int64_t*  ts_ptr = t_starts_arr.data();
    const int64_t*  te_ptr = t_ends_arr.data();
    const int64_t*  tl_ptr = t_lengths_arr.data();
    const double*   sp_ptr = transcript_spans_arr.data();
    const double*   el_ptr = exonic_lengths_arr.data();

    // Mutable output pointers
    double* em_out    = em_counts_out.data();
    double* hc_out    = em_high_conf_out.data();
    double* nrna_out  = nrna_em_counts_out.data();
    double* gdna_out  = gdna_locus_counts_out.data();
    double* psum_out  = posterior_sum_out.data();
    double* nass_out  = n_assigned_out.data();

    // Pre-compute per-transcript unambig row sums
    std::vector<double> unambig_row_sums(N_T, 0.0);
    for (int t = 0; t < N_T; ++t) {
        double s = 0.0;
        for (int c = 0; c < N_COLS; ++c) s += uac[t * N_COLS + c];
        unambig_row_sums[t] = s;
    }

    // Allocate scratch buffer for global-to-local mapping
    int local_map_size = nrna_base_index + N_T + 1;

    // Allocate per-locus result arrays
    std::vector<double> locus_mrna_vec(n_loci, 0.0);
    std::vector<double> locus_nrna_vec(n_loci, 0.0);
    std::vector<double> locus_gdna_vec(n_loci, 0.0);
    double* locus_mrna_data = locus_mrna_vec.data();
    double* locus_nrna_data = locus_nrna_vec.data();
    double* locus_gdna_data = locus_gdna_vec.data();

    std::atomic<double> total_gdna_em{0.0};

    // Determine actual thread count
    int actual_threads = n_threads;
    if (actual_threads <= 0) {
        int hw = static_cast<int>(std::thread::hardware_concurrency());
        actual_threads = (hw > 0) ? hw : 1;
    }
    if (actual_threads < 1) actual_threads = 1;

    // Release the GIL for the compute-intensive parallel section.
    // All numpy data pointers have been extracted above; pure C++ from here
    // until the scoped block ends and GIL is re-acquired for numpy return
    // array creation.
    {
        nb::gil_scoped_release release;

        // --- Main locus loop (thread-parallel via std::thread) ---
        //
        // Key safety property: loci are connected components of the
        // transcript-unit bipartite graph, so each transcript and each unit
        // belongs to exactly one locus.  This means different loci write to
        // DISJOINT indices in the output arrays (em_out, hc_out, nrna_out,
        // gdna_out, psum_out, nass_out).  No data races on those arrays.
        //
        // The only shared scalar accumulator is total_gdna_em, handled via
        // std::atomic<double>.
        //
        // Each thread gets its own LocusSubProblem and local_map scratch.
        // Work is distributed dynamically via an atomic counter (chunks of 16).
        constexpr int CHUNK_SIZE = 16;
        std::atomic<int> next_locus{0};

        auto worker_fn = [&]() {
            // Thread-private scratch state (allocated once per thread, reused
            // across all loci assigned to this thread by dynamic scheduling).
            LocusSubProblem sub;
            std::vector<int32_t> local_map(local_map_size, -1);
            double local_gdna = 0.0;

            for (;;) {
                int chunk_start = next_locus.fetch_add(CHUNK_SIZE,
                    std::memory_order_relaxed);
                if (chunk_start >= n_loci) break;
                int chunk_end = std::min(chunk_start + CHUNK_SIZE, n_loci);
                for (int li = chunk_start; li < chunk_end; ++li) {
            // Get locus transcript and unit index ranges
            auto t_start = lt_off[li];
            auto t_end   = lt_off[li + 1];
            auto u_start = lu_off[li];
            auto u_end   = lu_off[li + 1];
            int n_t = static_cast<int>(t_end - t_start);
            int n_u = static_cast<int>(u_end - u_start);

            if (n_u == 0) {
                locus_mrna_data[li] = 0.0;
                locus_nrna_data[li] = 0.0;
                locus_gdna_data[li] = 0.0;
                continue;
            }

            const int32_t* t_arr = lt_fl + t_start;
            const int32_t* u_arr = lu_fl + u_start;
            double gdna_init = gi_ptr[li];

            // 1. Extract sub-problem
            extract_locus_sub_problem(
                sub, t_arr, n_t, u_arr, n_u, gdna_init,
                goff, gti, gll, gcw, gtxs, gtxe, gcc,
                gis, ggll, ggfp, glti, glcc,
                nrna_base_index,
                unambig_row_sums.data(), nri, nfa, nfb,
                ts_ptr, te_ptr, tl_ptr, sp_ptr, el_ptr,
                mean_frag,
                local_map.data(), local_map_size);

            int nc = sub.n_components;
            int sub_n_t = sub.n_t;
            size_t n_candidates = sub.t_indices.size();
            int n_local_units = sub.n_local_units;

            // 2. Apply bias correction (mutates log_liks in-place)
            if (n_candidates > 0) {
                apply_bias_correction_uniform(
                    sub.log_liks.data(),
                    sub.t_indices.data(),
                    sub.tx_starts.data(),
                    sub.tx_ends.data(),
                    sub.bias_profiles.data(),
                    n_candidates);
            }

            // 3. Handle empty sub-problem
            if (n_local_units == 0 || n_candidates == 0) {
                locus_mrna_data[li] = 0.0;
                locus_nrna_data[li] = 0.0;
                locus_gdna_data[li] = 0.0;
                continue;
            }

            // 4. Compute log(effective_lengths) — all 1.0, so log = 0.0
            std::vector<double> log_eff_len(nc, 0.0);

            // 5. Build equivalence classes
            auto ec_data = build_equiv_classes(
                sub.offsets.data(),
                sub.t_indices.data(),
                sub.log_liks.data(),
                sub.coverage_wts.data(),
                n_local_units);

            // 6. Coverage-weighted warm start + OVR prior
            std::vector<double> prior(nc);
            std::vector<double> theta_init(nc);
            compute_ovr_prior_and_warm_start(
                ec_data,
                sub.unambig_totals.data(),
                sub.eligible.data(),
                prior_alpha, prior_gamma,
                prior.data(), theta_init.data(), nc);

            // 7. Run SQUAREM (linked or classic)
            bool linked = (sub_n_t > 0);
            EMResult result;
            if (linked) {
                // Collapse warm start to theta_t space
                int ns = sub_n_t + 1;
                std::vector<double> theta_t_init(ns);
                for (int i = 0; i < sub_n_t; ++i) {
                    theta_t_init[i] = theta_init[i] + theta_init[sub_n_t + i];
                }
                theta_t_init[sub_n_t] = theta_init[2 * sub_n_t];

                // nrna_frac_init = prior mean, clamped
                std::vector<double> nrna_frac_init(sub_n_t);
                for (int i = 0; i < sub_n_t; ++i) {
                    double sa = sub.nrna_frac_alpha[i] + sub.nrna_frac_beta[i];
                    nrna_frac_init[i] = (sa > 0.0)
                        ? sub.nrna_frac_alpha[i] / sa : 0.5;
                    if (nrna_frac_init[i] < NRNA_FRAC_CLAMP_EPS)
                        nrna_frac_init[i] = NRNA_FRAC_CLAMP_EPS;
                    else if (nrna_frac_init[i] > 1.0 - NRNA_FRAC_CLAMP_EPS)
                        nrna_frac_init[i] = 1.0 - NRNA_FRAC_CLAMP_EPS;
                }

                result = linked_run_squarem(
                    ec_data, log_eff_len.data(),
                    sub.unambig_totals.data(),
                    prior.data(),
                    theta_t_init.data(),
                    nrna_frac_init.data(),
                    sub.nrna_frac_alpha.data(),
                    sub.nrna_frac_beta.data(),
                    sub_n_t, nc,
                    max_iterations, convergence_delta,
                    prune_threshold);
            } else {
                result = run_squarem(
                    ec_data, log_eff_len.data(),
                    sub.unambig_totals.data(),
                    prior.data(),
                    theta_init.data(),
                    nc, max_iterations, convergence_delta,
                    use_vbem, prune_threshold);
            }

            // 8. Assign posteriors (writes to disjoint transcript
            //    indices — safe across threads, no atomics needed)
            double locus_mrna = 0.0, locus_nrna = 0.0, locus_gdna = 0.0;
            assign_posteriors(
                sub, result.theta.data(), confidence_threshold,
                em_out, hc_out, nrna_out, gdna_out,
                psum_out, nass_out,
                locus_mrna, locus_nrna, locus_gdna,
                N_T, N_COLS);

            locus_mrna_data[li] = locus_mrna;
            locus_nrna_data[li] = locus_nrna;
            locus_gdna_data[li] = locus_gdna;

            local_gdna += locus_gdna;
                } // end chunk loop
            } // end work-stealing loop

            // Flush thread-local gdna to shared atomic
            double prev = total_gdna_em.load(std::memory_order_relaxed);
            while (!total_gdna_em.compare_exchange_weak(
                prev, prev + local_gdna,
                std::memory_order_relaxed, std::memory_order_relaxed)) {}
        }; // end worker_fn

        if (actual_threads <= 1) {
            // Single-threaded: run directly, no thread overhead
            worker_fn();
        } else {
            std::vector<std::thread> threads;
            threads.reserve(actual_threads);
            for (int t = 0; t < actual_threads; ++t)
                threads.emplace_back(worker_fn);
            for (auto& th : threads)
                th.join();
        }
    } // end gil_scoped_release — GIL re-acquired here

    double total_gdna_em_val = total_gdna_em.load(std::memory_order_relaxed);

    // Package per-locus results as numpy arrays (copy from vector)
    size_t shape[1] = { static_cast<size_t>(n_loci) };

    // Allocate new arrays, copy data, and transfer ownership via capsules
    auto* mrna_copy = new double[n_loci];
    auto* nrna_copy = new double[n_loci];
    auto* gdna_copy = new double[n_loci];
    std::memcpy(mrna_copy, locus_mrna_vec.data(), n_loci * sizeof(double));
    std::memcpy(nrna_copy, locus_nrna_vec.data(), n_loci * sizeof(double));
    std::memcpy(gdna_copy, locus_gdna_vec.data(), n_loci * sizeof(double));

    nb::capsule mrna_owner(mrna_copy, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    nb::capsule nrna_owner(nrna_copy, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    nb::capsule gdna_owner(gdna_copy, [](void* p) noexcept { delete[] static_cast<double*>(p); });

    return std::make_tuple(
        total_gdna_em_val,
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(
            mrna_copy, 1, shape, std::move(mrna_owner)),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(
            nrna_copy, 1, shape, std::move(nrna_owner)),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(
            gdna_copy, 1, shape, std::move(gdna_owner))
    );
}


// ================================================================
// Phase 3 — C++ Union-Find Connected Components
// ================================================================
//
// Replaces scipy.sparse.csgraph.connected_components for locus building.
// Uses disjoint-set (union-find) with path compression and union by rank.
// Time: O(N_candidates * α(N_transcripts)) ≈ O(N_candidates).
//
// For each EM unit, all candidate transcripts in that unit are connected.
// We union the first transcript in each unit with every other transcript
// in that unit.  After processing all units, find() on each transcript
// gives its component root, which we relabel to sequential 0-based IDs.

namespace {

/// Disjoint-set forest with path compression and union by rank.
struct UnionFind {
    std::vector<int32_t> parent;
    std::vector<int32_t> rank;

    explicit UnionFind(int32_t n) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int32_t find(int32_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];  // path halving
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

}  // anonymous namespace


/// Build loci (connected components) from the global CSR fragment data.
///
/// Arguments:
///   offsets      — int64[n_units + 1]: CSR row pointers
///   t_indices    — int32[n_candidates]: candidate transcript/nRNA indices
///   nrna_base    — int32: first nRNA shadow index (= num_transcripts)
///   n_transcripts — int32: number of real transcripts
///
/// Returns a tuple of:
///   labels       — int32[n_transcripts]: component label per transcript
///                  (sequential 0-based, only active transcripts get real
///                  labels; inactive transcripts get -1)
///   n_components — int32: number of connected components
static nb::tuple connected_components_native(
    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig>  offsets_arr,
    nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig>  t_indices_arr,
    int32_t nrna_base,
    int32_t n_transcripts)
{
    const int64_t* offsets = offsets_arr.data();
    const int32_t* t_idx   = t_indices_arr.data();
    const int64_t  n_units = static_cast<int64_t>(offsets_arr.shape(0)) - 1;

    if (n_units <= 0 || n_transcripts <= 0) {
        // Return all -1 labels, 0 components
        auto* labels = new int32_t[n_transcripts];
        std::fill(labels, labels + n_transcripts, -1);
        size_t shape[1] = {static_cast<size_t>(n_transcripts)};
        nb::capsule owner(labels, [](void* p) noexcept { delete[] static_cast<int32_t*>(p); });
        return nb::make_tuple(
            nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(labels, 1, shape, std::move(owner)),
            nb::int_(0)
        );
    }

    UnionFind uf(n_transcripts);

    // Track which transcripts are actually referenced
    std::vector<bool> active(n_transcripts, false);

    for (int64_t u = 0; u < n_units; ++u) {
        int64_t start = offsets[u];
        int64_t end   = offsets[u + 1];
        if (start >= end) continue;

        // Find first valid transcript in this unit
        int32_t first_t = -1;
        for (int64_t j = start; j < end; ++j) {
            int32_t t = t_idx[j];
            if (t >= nrna_base) t -= nrna_base;  // map nRNA → base transcript
            if (t >= 0 && t < n_transcripts) {
                first_t = t;
                active[t] = true;
                break;
            }
        }
        if (first_t < 0) continue;

        // Union all other transcripts in this unit with the first
        for (int64_t j = start + 1; j < end; ++j) {
            int32_t t = t_idx[j];
            if (t >= nrna_base) t -= nrna_base;
            if (t >= 0 && t < n_transcripts && t != first_t) {
                active[t] = true;
                uf.unite(first_t, t);
            }
        }
    }

    // Assign sequential component labels to active transcripts.
    // labels[t] = component index (0-based) for active transcripts, -1 otherwise.
    std::vector<int32_t> labels(n_transcripts, -1);
    std::unordered_map<int32_t, int32_t> root_to_label;
    int32_t n_comp = 0;

    for (int32_t t = 0; t < n_transcripts; ++t) {
        if (!active[t]) continue;
        int32_t root = uf.find(t);
        auto it = root_to_label.find(root);
        if (it == root_to_label.end()) {
            root_to_label[root] = n_comp;
            labels[t] = n_comp;
            ++n_comp;
        } else {
            labels[t] = it->second;
        }
    }

    // --- Build per-component transcript and unit lists (CSR form) ---

    // 1. Count transcripts per component
    std::vector<int64_t> comp_t_counts(n_comp, 0);
    for (int32_t t = 0; t < n_transcripts; ++t) {
        if (labels[t] >= 0) comp_t_counts[labels[t]]++;
    }

    // 2. Assign each unit to a component via its first transcript
    std::vector<int32_t> unit_label(n_units, -1);
    for (int64_t u = 0; u < n_units; ++u) {
        int64_t start = offsets[u];
        int64_t end   = offsets[u + 1];
        for (int64_t j = start; j < end; ++j) {
            int32_t t = t_idx[j];
            if (t >= nrna_base) t -= nrna_base;
            if (t >= 0 && t < n_transcripts && labels[t] >= 0) {
                unit_label[u] = labels[t];
                break;
            }
        }
    }

    // 3. Count units per component
    std::vector<int64_t> comp_u_counts(n_comp, 0);
    for (int64_t u = 0; u < n_units; ++u) {
        if (unit_label[u] >= 0) comp_u_counts[unit_label[u]]++;
    }

    // 4. Build CSR offsets via prefix sum
    auto* ct_off = new int64_t[n_comp + 1];
    auto* cu_off = new int64_t[n_comp + 1];
    ct_off[0] = 0;
    cu_off[0] = 0;
    for (int32_t c = 0; c < n_comp; ++c) {
        ct_off[c + 1] = ct_off[c] + comp_t_counts[c];
        cu_off[c + 1] = cu_off[c] + comp_u_counts[c];
    }
    int64_t total_t = ct_off[n_comp];
    int64_t total_u = cu_off[n_comp];

    // 5. Fill flat arrays (iterate ascending → output is sorted)
    auto* ct_flat = new int32_t[std::max(total_t, int64_t(1))];
    auto* cu_flat = new int32_t[std::max(total_u, int64_t(1))];

    // Reuse counts as write cursors
    std::fill(comp_t_counts.begin(), comp_t_counts.end(), 0);
    std::fill(comp_u_counts.begin(), comp_u_counts.end(), 0);

    for (int32_t t = 0; t < n_transcripts; ++t) {
        int32_t c = labels[t];
        if (c < 0) continue;
        ct_flat[ct_off[c] + comp_t_counts[c]++] = t;
    }
    for (int64_t u = 0; u < n_units; ++u) {
        int32_t c = unit_label[u];
        if (c < 0) continue;
        cu_flat[cu_off[c] + comp_u_counts[c]++] = static_cast<int32_t>(u);
    }

    // Wrap in numpy arrays with capsule ownership
    size_t ct_off_shape[1] = {static_cast<size_t>(n_comp + 1)};
    size_t ct_flat_shape[1] = {static_cast<size_t>(total_t)};
    size_t cu_off_shape[1] = {static_cast<size_t>(n_comp + 1)};
    size_t cu_flat_shape[1] = {static_cast<size_t>(total_u)};

    nb::capsule own_ct_off(ct_off, [](void* p) noexcept { delete[] static_cast<int64_t*>(p); });
    nb::capsule own_ct_flat(ct_flat, [](void* p) noexcept { delete[] static_cast<int32_t*>(p); });
    nb::capsule own_cu_off(cu_off, [](void* p) noexcept { delete[] static_cast<int64_t*>(p); });
    nb::capsule own_cu_flat(cu_flat, [](void* p) noexcept { delete[] static_cast<int32_t*>(p); });

    return nb::make_tuple(
        nb::int_(n_comp),
        nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(ct_off, 1, ct_off_shape, std::move(own_ct_off)),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(ct_flat, 1, ct_flat_shape, std::move(own_ct_flat)),
        nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(cu_off, 1, cu_off_shape, std::move(own_cu_off)),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(cu_flat, 1, cu_flat_shape, std::move(own_cu_flat))
    );
}


// ================================================================
// nanobind module definition
// ================================================================

NB_MODULE(_em_impl, m) {
    m.doc() = "C++ EM solver for rigel locus-level abundance estimation.\n\n"
              "Provides run_locus_em_native() which replaces the Python EM hot path:\n"
              "_em_step, _vbem_step, _build_equiv_classes, SQUAREM loop.\n\n"
              "Also provides batch_locus_em() which replaces the entire per-locus\n"
              "Python for-loop with a single C++ call.";

    m.def("run_locus_em_native", &run_locus_em_native,
          nb::arg("offsets"),
          nb::arg("t_indices"),
          nb::arg("log_liks"),
          nb::arg("coverage_wts"),
          nb::arg("tx_starts"),
          nb::arg("tx_ends"),
          nb::arg("bias_profiles"),
          nb::arg("unambig_totals"),
          nb::arg("effective_lens"),
          nb::arg("prior_eligible"),
          nb::arg("n_components"),
          nb::arg("prior_alpha"),
          nb::arg("prior_gamma"),
          nb::arg("max_iterations"),
          nb::arg("convergence_delta"),
          nb::arg("use_vbem"),
          nb::arg("prune_threshold"),
          nb::arg("n_transcripts"),
          nb::arg("nrna_frac_alpha"),
          nb::arg("nrna_frac_beta"),
          "Run EM for a single locus sub-problem.\n\n"
          "Takes CSR per-locus data + config, returns (theta, alpha, em_totals, nrna_frac).\n"
          "When n_transcripts > 0, uses the linked mRNA-nRNA model.\n"
          "Replaces the Python EM hot path with a single C++ call.");

    m.def("batch_locus_em", &batch_locus_em,
          nb::arg("g_offsets"),
          nb::arg("g_t_indices"),
          nb::arg("g_log_liks"),
          nb::arg("g_coverage_wts"),
          nb::arg("g_tx_starts"),
          nb::arg("g_tx_ends"),
          nb::arg("g_count_cols"),
          nb::arg("g_is_spliced"),
          nb::arg("g_gdna_log_liks"),
          nb::arg("g_genomic_footprints"),
          nb::arg("g_locus_t_indices"),
          nb::arg("g_locus_count_cols"),
          nb::arg("nrna_base_index"),
          nb::arg("locus_t_offsets"),
          nb::arg("locus_t_flat"),
          nb::arg("locus_u_offsets"),
          nb::arg("locus_u_flat"),
          nb::arg("gdna_inits"),
          nb::arg("unambig_counts"),
          nb::arg("nrna_init"),
          nb::arg("nrna_frac_alpha"),
          nb::arg("nrna_frac_beta"),
          nb::arg("t_starts"),
          nb::arg("t_ends"),
          nb::arg("t_lengths"),
          nb::arg("transcript_spans"),
          nb::arg("exonic_lengths"),
          nb::arg("em_counts_out"),
          nb::arg("em_high_conf_out"),
          nb::arg("nrna_em_counts_out"),
          nb::arg("gdna_locus_counts_out"),
          nb::arg("posterior_sum_out"),
          nb::arg("n_assigned_out"),
          nb::arg("mean_frag"),
          nb::arg("max_iterations"),
          nb::arg("convergence_delta"),
          nb::arg("prior_alpha"),
          nb::arg("prior_gamma"),
          nb::arg("use_vbem"),
          nb::arg("prune_threshold"),
          nb::arg("confidence_threshold"),
          nb::arg("n_transcripts_total"),
          nb::arg("n_splice_strand_cols"),
          nb::arg("n_threads") = 0,
          "Run locus EM for ALL loci in a single C++ call.\n\n"
          "Replaces the Python per-locus for-loop:\n"
          "  build_locus_em_data -> run_locus_em -> assign_locus_ambiguous\n"
          "Returns (total_gdna_em, locus_mrna, locus_nrna, locus_gdna).\n\n"
          "n_threads: 0 = all cores, 1 = sequential, N = cap at N threads.");

    m.def("connected_components", &connected_components_native,
          nb::arg("offsets"),
          nb::arg("t_indices"),
          nb::arg("nrna_base"),
          nb::arg("n_transcripts"),
          "Find connected components of the fragment→transcript overlap graph.\n\n"
          "Uses union-find with path compression and union by rank.\n"
          "Returns (n_comp, comp_t_offsets, comp_t_flat, comp_u_offsets,\n"
          "comp_u_flat) where the CSR pairs (offsets, flat) give sorted\n"
          "transcript indices and unit indices for each component.");

    // ----------------------------------------------------------------
    // Export constants so Python imports from this single source of truth.
    // ----------------------------------------------------------------
    m.attr("EM_LOG_EPSILON")         = EM_LOG_EPSILON;
    m.attr("MAX_FRAG_LEN")          = MAX_FRAG_LEN;
    m.attr("SQUAREM_BUDGET_DIVISOR") = SQUAREM_BUDGET_DIVISOR;
    m.attr("NRNA_FRAC_CLAMP_EPS")   = NRNA_FRAC_CLAMP_EPS;
    m.attr("EM_PRIOR_EPSILON")       = EM_PRIOR_EPSILON;
}
