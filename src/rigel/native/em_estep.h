// em_estep.h — the locus EM's equivalence classes and its E-step over them.
//
// The units sharing one candidate component set form an equivalence class, held as dense (n, k) matrices. The E-step
// turns the components' log-weights into every unit's posterior and adds the posteriors' column sums into the
// components' totals, in an order that is the same at every thread count.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "fast_exp.h"
#include "thread_pool.h"

namespace rigel {

struct EmEquivClass {
    std::vector<int32_t> comp_idx;  // k component indices
    std::vector<double>  ll_flat;   // n*k log-likelihoods (row-major)
    std::vector<double>  wt_flat;   // n*k coverage weights (row-major)
    int n;  // number of units in this class
    int k;  // number of components per unit
};

// FNV-1a hash of a candidate component set.
struct VecHash {
    size_t operator()(const std::vector<int32_t>& v) const noexcept {
        size_t h = 14695981039346656037ULL;
        for (int32_t val : v) {
            h ^= static_cast<size_t>(static_cast<uint32_t>(val));
            h *= 1099511628211ULL;
        }
        return h;
    }
};

// The units of a candidate CSR grouped by their (sorted) candidate set, in comp_idx order.
inline std::vector<EmEquivClass> build_equiv_classes(
    const int64_t* offsets,
    const int32_t* t_indices,
    const double*  log_liks,
    const double*  coverage_wts,
    int n_units)
{
    if (n_units == 0) return {};

    std::unordered_map<std::vector<int32_t>, std::vector<int>, VecHash> class_map;
    class_map.reserve(static_cast<size_t>(n_units));
    for (int u = 0; u < n_units; ++u) {
        auto start = static_cast<size_t>(offsets[u]);
        auto end   = static_cast<size_t>(offsets[u + 1]);
        if (start == end) continue;
        std::vector<int32_t> key(t_indices + start, t_indices + end);
        class_map[std::move(key)].push_back(u);
    }

    std::vector<EmEquivClass> result;
    result.reserve(class_map.size());
    for (auto& [key, unit_list] : class_map) {
        int k = static_cast<int>(key.size());
        int n = static_cast<int>(unit_list.size());
        EmEquivClass ec;
        ec.comp_idx = key;
        ec.n = n;
        ec.k = k;
        ec.ll_flat.resize(static_cast<size_t>(n) * k);
        ec.wt_flat.resize(static_cast<size_t>(n) * k);
        for (int i = 0; i < n; ++i) {
            auto s = static_cast<size_t>(offsets[unit_list[i]]);
            for (int j = 0; j < k; ++j) {
                ec.ll_flat[static_cast<size_t>(i) * k + j] = log_liks[s + static_cast<size_t>(j)];
                ec.wt_flat[static_cast<size_t>(i) * k + j] = coverage_wts[s + static_cast<size_t>(j)];
            }
        }
        result.push_back(std::move(ec));
    }

    // Sort the classes by comp_idx. `unordered_map` iterates in no fixed order, and the E-step's column sums
    // depend on addition order in floating point, so an unstable order gives ULP differences that SQUAREM
    // amplifies. Rows within a class need no sort: `build_multi_loci` orders each locus's units by `frag_id`
    // (BAM order), so `unit_list` is already canonical.
    std::sort(result.begin(), result.end(),
              [](const EmEquivClass& a, const EmEquivClass& b) { return a.comp_idx < b.comp_idx; });
    return result;
}

// Kahan summation: O(eps) error rather than O(n·eps), for the E-step's column sums over large classes.
struct KahanAccumulator {
    double sum = 0.0;
    double c = 0.0;

    inline void add(double val) {
        double y = val - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
};

// row[j] ← exp(row[j] − max_val) for j < k, and their sum. A whole NEON block below fast_exp's cutoff is written as
// zeros without evaluating it.
inline double exp_shifted_sum(double* row, int k, double max_val) {
    double row_sum = 0.0;
    int j = 0;
#if RIGEL_HAS_NEON
    {
        float64x2_t sum_v = vdupq_n_f64(0.0);
        const float64x2_t max_v = vdupq_n_f64(max_val);
        const float64x2_t cutoff_v = vdupq_n_f64(rigel::detail::EXP_CUTOFF);
        const float64x2_t zero_v = vdupq_n_f64(0.0);
        for (; j + 2 <= k; j += 2) {
            float64x2_t v = vld1q_f64(row + j);
            v = vsubq_f64(v, max_v);
            uint64x2_t mask = vcltq_f64(v, cutoff_v);
            if (vgetq_lane_u64(mask, 0) & vgetq_lane_u64(mask, 1)) {
                vst1q_f64(row + j, zero_v);
                continue;
            }
            v = rigel::fast_exp_neon(v);
            vst1q_f64(row + j, v);
            sum_v = vaddq_f64(sum_v, v);
        }
        row_sum = vgetq_lane_f64(sum_v, 0) + vgetq_lane_f64(sum_v, 1);
    }
#endif
    for (; j < k; ++j) {
        double e = rigel::fast_exp_scalar(row[j] - max_val);
        row[j] = e;
        row_sum += e;
    }
    return row_sum;
}

// For every row of one equivalence class: the posteriors from log_weights, normalised by log-sum-exp, with the k
// column sums Kahan-summed in row order into col_sums. Thread-safe: the class is read-only, each row's posteriors
// live in a stack-local buffer, and col_sums is the caller's slot. A row no component can emit (every weight -inf)
// is outside the model and skipped.
inline void em_step_kernel(
    const EmEquivClass& ec,
    const double*     log_weights,   // [n_components]
    double*           col_sums)      // [k], overwritten
{
    const int k = ec.k;
    const double* ll = ec.ll_flat.data();
    const int32_t* cidx = ec.comp_idx.data();

    // Row buffer: on the stack for k <= MAX_K_STACK, on the heap above it.
    constexpr int MAX_K_STACK = 512;
    double stack_row[MAX_K_STACK];
    std::vector<double> heap_row;
    double* row = stack_row;
    if (k > MAX_K_STACK) {
        heap_row.resize(k);
        row = heap_row.data();
    }

    // One Kahan accumulator per column, fed in row order.
    std::vector<KahanAccumulator> col_acc(k);
    for (int i = 0; i < ec.n; ++i) {
        double max_val = ll[i * k] + log_weights[cidx[0]];
        row[0] = max_val;
        for (int j = 1; j < k; ++j) {
            double val = ll[i * k + j] + log_weights[cidx[j]];
            row[j] = val;
            if (val > max_val) max_val = val;
        }
        if (max_val == -std::numeric_limits<double>::infinity()) continue;
        const double row_sum = exp_shifted_sum(row, k, max_val);
        if (row_sum > 0.0 && std::isfinite(row_sum)) {
            double inv_sum = 1.0 / row_sum;
            for (int c = 0; c < k; ++c) {
                row[c] *= inv_sum;
                col_acc[c].add(row[c]);
            }
        }
    }
    for (int j = 0; j < k; ++j) col_sums[j] = col_acc[j].sum;
}

// The E-step: per-class column sums, added in class order. Each equivalence class writes its column sums to its own
// slot, and the slots are added into em_totals in class order, so the totals are bit-identical for any thread count.
// A class is never split across threads, since a split would re-associate its row sums. Built once per locus and
// reused by every step.
class EStep {
public:
    // pool: the locus's threads, or nullptr to run on the caller's thread alone.
    EStep(const std::vector<EmEquivClass>& ec_data, ThreadPool* pool)
        : ec_data_(ec_data), pool_(pool), slot_(ec_data.size() + 1, 0)
    {
        const int n_ec = static_cast<int>(ec_data.size());
        for (int i = 0; i < n_ec; ++i) slot_[i + 1] = slot_[i] + static_cast<size_t>(ec_data[i].k);
        sums_.resize(slot_[n_ec]);
        if (pool_) deal_classes();
    }

    // Adds the posteriors' column sums into em_totals [n_components], zeroed by the caller.
    void run(const double* log_weights, double* em_totals) {
        const int n_ec = static_cast<int>(ec_data_.size());
        auto run_class = [&](int i) { em_step_kernel(ec_data_[i], log_weights, sums_.data() + slot_[i]); };
        auto add_class = [&](int i) {
            const int32_t* cidx = ec_data_[i].comp_idx.data();
            const double* s = sums_.data() + slot_[i];
            for (int j = 0; j < ec_data_[i].k; ++j) em_totals[cidx[j]] += s[j];
        };
        if (pool_) {
            pool_->run_parallel([&](int tid) {
                for (int c = thread_start_[tid]; c < thread_start_[tid + 1]; ++c) run_class(thread_classes_[c]);
            });
            for (int i = 0; i < n_ec; ++i) add_class(i);
        } else {
            // The same additions in the same order, each class's while its slot is still in cache.
            for (int i = 0; i < n_ec; ++i) {
                run_class(i);
                add_class(i);
            }
        }
    }

private:
    // The classes dealt heaviest first, each to the thread with the least work so far.
    void deal_classes() {
        const int n_ec = static_cast<int>(ec_data_.size());
        auto work = [&](int i) { return static_cast<int64_t>(ec_data_[i].n) * ec_data_[i].k; };
        std::vector<int> order(static_cast<size_t>(n_ec));
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return work(a) != work(b) ? work(a) > work(b) : a < b; });
        const int n_threads = pool_->n_threads();
        std::vector<int64_t> load(static_cast<size_t>(n_threads), 0);
        std::vector<std::vector<int>> dealt(static_cast<size_t>(n_threads));
        for (int i : order) {
            const int t = static_cast<int>(std::min_element(load.begin(), load.end()) - load.begin());
            dealt[t].push_back(i);
            load[t] += work(i);
        }
        thread_start_.assign(static_cast<size_t>(n_threads) + 1, 0);
        for (int t = 0; t < n_threads; ++t) {
            thread_start_[t + 1] = thread_start_[t] + static_cast<int>(dealt[t].size());
            thread_classes_.insert(thread_classes_.end(), dealt[t].begin(), dealt[t].end());
        }
    }

    const std::vector<EmEquivClass>& ec_data_;
    ThreadPool* pool_;
    std::vector<size_t> slot_;   // [n_ec + 1]: class i's k sums are sums_[slot_[i], slot_[i + 1])
    std::vector<double> sums_;
    std::vector<int> thread_start_;    // threaded only: thread t runs thread_classes_[thread_start_[t], thread_start_[t + 1])
    std::vector<int> thread_classes_;
};

// Σ_units logsumexp_c(log_weights[c] + log_lik_{u,c}): the EM objective's data term, which ranks two converged
// solutions by likelihood.
inline double marginal_data_loglik(const std::vector<EmEquivClass>& ec_data, const double* log_weights) {
    double total = 0.0;
    for (const auto& ec : ec_data) {
        const int k = ec.k;
        const double* ll = ec.ll_flat.data();
        const int32_t* cidx = ec.comp_idx.data();
        for (int i = 0; i < ec.n; ++i) {
            double mx = ll[i * k] + log_weights[cidx[0]];
            for (int j = 1; j < k; ++j) {
                double v = ll[i * k + j] + log_weights[cidx[j]];
                if (v > mx) mx = v;
            }
            if (!(mx > -std::numeric_limits<double>::infinity())) continue;  // no component can emit it
            double s = 0.0;
            for (int j = 0; j < k; ++j) s += std::exp((ll[i * k + j] + log_weights[cidx[j]]) - mx);
            total += mx + std::log(s);
        }
    }
    return total;
}

}  // namespace rigel
