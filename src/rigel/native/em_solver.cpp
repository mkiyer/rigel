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
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <vector>

#include <atomic>
#include <thread>

#include "fast_exp.h"
#include "thread_pool.h"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

// ================================================================
// Constants — single source of truth, exported to Python via module attrs.
// Python imports these from rigel._em_impl instead of redefining.
// ================================================================

static constexpr double EM_LOG_EPSILON = 1e-300;
static constexpr int    MAX_FRAG_LEN  = 1000000;
static constexpr int    SQUAREM_BUDGET_DIVISOR = 3;

// Target number of element-operations per E-step parallel task.
// Each equivalence class row with k components costs O(k); tasks are
// sized to ~ESTEP_TASK_WORK_TARGET / k rows for load-balanced threading.
static constexpr int    ESTEP_TASK_WORK_TARGET = 4096;

// VBEM SQUAREM clamp floor: minimum alpha value after SQUAREM extrapolation
// or stabilization.  Prevents components from entering the digamma absorbing
// regime (psi(a) ~ -1/a for small a) where recovery is impossible.
// At 0.1, psi(0.1) ~ -10.4, giving weight ~3e-5 — small enough not to steal
// mass, but large enough that a component with genuine read support can
// accumulate evidence and recover.  Any value >= 0.01 avoids the absorbing
// barrier; values below ~0.01 are effectively dead in double precision.
static constexpr double VBEM_CLAMP_FLOOR = 0.1;

// Assignment mode constants (must match Python _ASSIGNMENT_MODE_MAP)
static constexpr int ASSIGN_FRACTIONAL = 0;
static constexpr int ASSIGN_MAP        = 1;
static constexpr int ASSIGN_SAMPLE     = 2;

// ================================================================
// Profiling instrumentation — per-locus and aggregate statistics
// ================================================================

using hrclock = std::chrono::steady_clock;

/// Per-locus profiling statistics collected during batch_locus_em_partitioned.
struct LocusProfile {
    int locus_idx = -1;
    int n_transcripts = 0;
    int n_units = 0;
    int n_components = 0;        // n_t + 1
    int n_equiv_classes = 0;
    int64_t ec_total_elements = 0; // sum of n*k across all ECs
    int max_ec_width = 0;        // max k across ECs
    int max_ec_depth = 0;        // max n across ECs
    int squarem_iterations = 0;
    int estep_threads_used = 0;
    bool is_mega_locus = false;

    // Sub-phase wall times in microseconds
    double extract_us = 0.0;
    double bias_us = 0.0;
    double build_ec_us = 0.0;
    double warm_start_us = 0.0;
    double squarem_us = 0.0;
    double assign_us = 0.0;
    double total_us = 0.0;

    // VBEM-specific: digamma calls per E-step
    int64_t digamma_calls_per_estep = 0;
};

// ================================================================
// SplitMix64 — lightweight, deterministic, thread-local PRNG
// ================================================================

struct SplitMix64 {
    uint64_t state;

    explicit SplitMix64(uint64_t seed) : state(seed) {}

    uint64_t next() {
        state += 0x9e3779b97f4a7c15ULL;
        uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    /// Return a uniform double in [0, 1).
    double uniform() {
        return static_cast<double>(next() >> 11) * 0x1.0p-53;
    }
};

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

    // Group units by candidate component set — store unit index u
    std::unordered_map<std::vector<int32_t>, std::vector<int>, VecHash> class_map;
    class_map.reserve(static_cast<size_t>(n_units));

    for (int u = 0; u < n_units; ++u) {
        auto start = static_cast<size_t>(offsets[u]);
        auto end   = static_cast<size_t>(offsets[u + 1]);
        if (start == end) continue;

        std::vector<int32_t> key(t_indices + start, t_indices + end);
        class_map[std::move(key)].push_back(u);
    }

    // Build dense matrices per class
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
        // scratch is no longer used by em_step_kernel_range (uses
        // stack-local row buffer instead), so skip allocation.

        for (int i = 0; i < n; ++i) {
            int u = unit_list[i];
            auto s = static_cast<size_t>(offsets[u]);
            for (int j = 0; j < k; ++j) {
                ec.ll_flat[static_cast<size_t>(i) * k + j] =
                    log_liks[s + static_cast<size_t>(j)];
                ec.wt_flat[static_cast<size_t>(i) * k + j] =
                    coverage_wts[s + static_cast<size_t>(j)];
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
    // potentially causing large cascading output differences.
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
// Kahan summation helper for numerical stability
// ================================================================
//
// Reduces FP accumulation error from O(n*eps) to O(eps).  Critical
// for column-sum accumulators in the E-step where n can be 500K+.

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

// ================================================================
// EM step kernel — the hot inner loop (range-based)
// ================================================================
//
// Processes rows [row_start, row_end) of one equivalence class:
// computes posteriors from log_weights, normalizes rows via
// log-sum-exp, accumulates column sums into em_totals with
// Kahan summation for numerical stability.
//
// Thread safety: multiple threads may process non-overlapping row
// ranges of the same EC concurrently.  Each writes to disjoint
// regions of ec.scratch (which is mutable).  em_totals must be
// thread-private when called from parallel_estep.

static inline void em_step_kernel_range(
    const EmEquivClass& ec,
    int               row_start,
    int               row_end,
    const double*     log_weights,   // [n_components]
    double*           em_totals)     // [n_components], accumulated
{
    const int k = ec.k;
    const double* ll = ec.ll_flat.data();
    const int32_t* cidx = ec.comp_idx.data();

    // Stack-local row buffer: avoids writing to the heap-allocated
    // ec.scratch buffer entirely.  ec.scratch is never read after
    // the E-step (assign_posteriors recomputes from theta + CSR),
    // so eliminating the writes reduces memory traffic.
    // For k <= 512 (covers all practical ECs), use stack allocation.
    constexpr int MAX_K_STACK = 512;
    double stack_row[MAX_K_STACK];
    std::vector<double> heap_row;
    double* row;
    if (k <= MAX_K_STACK) {
        row = stack_row;
    } else {
        heap_row.resize(k);
        row = heap_row.data();
    }

    // Per-column Kahan accumulators for fused column-sum accumulation.
    // Fusing the column sums into the row-processing loop eliminates a
    // separate column-stride pass over scratch, improving cache locality.
    // Each column's accumulator sees values in the same order as before
    // (row_start to row_end), so results are bit-for-bit identical.
    std::vector<KahanAccumulator> col_acc(k);

    for (int i = row_start; i < row_end; ++i) {

        // Compute row values and find max
        double max_val = ll[i * k] + log_weights[cidx[0]];
        row[0] = max_val;
        for (int j = 1; j < k; ++j) {
            double val = ll[i * k + j] + log_weights[cidx[j]];
            row[j] = val;
            if (val > max_val) max_val = val;
        }

        // Exp and sum — vectorized with early-zero skip
        double row_sum = 0.0;
        int j = 0;

#if RIGEL_HAS_AVX512F
        {
            __m512d sum_v = _mm512_setzero_pd();
            const __m512d max_v = _mm512_set1_pd(max_val);
            const __m512d cutoff_v = _mm512_set1_pd(rigel::detail::EXP_CUTOFF);

            for (; j + 8 <= k; j += 8) {
                __m512d v = _mm512_loadu_pd(row + j);
                v = _mm512_sub_pd(v, max_v);

                // Whole-vector early-zero skip: if ALL 8 lanes < cutoff
                __mmask8 mask = _mm512_cmp_pd_mask(v, cutoff_v, _CMP_LT_OQ);
                if (mask == 0xFF) {
                    _mm512_storeu_pd(row + j, _mm512_setzero_pd());
                    continue;
                }

                v = rigel::fast_exp_avx512(v);
                _mm512_storeu_pd(row + j, v);
                sum_v = _mm512_add_pd(sum_v, v);
            }
            row_sum = _mm512_reduce_add_pd(sum_v);
        }
        // AVX-512 scalar tail
        for (; j < k; ++j) {
            double x = row[j] - max_val;
            double e = rigel::fast_exp_scalar(x);
            row[j] = e;
            row_sum += e;
        }
#elif RIGEL_HAS_AVX2 && RIGEL_HAS_FMA
        {
            __m256d sum_v = _mm256_setzero_pd();
            const __m256d max_v = _mm256_set1_pd(max_val);
            const __m256d cutoff_v = _mm256_set1_pd(rigel::detail::EXP_CUTOFF);

            for (; j + 4 <= k; j += 4) {
                __m256d v = _mm256_loadu_pd(row + j);
                v = _mm256_sub_pd(v, max_v);

                // Whole-vector early-zero skip: if ALL 4 lanes < cutoff
                __m256d cmp = _mm256_cmp_pd(v, cutoff_v, _CMP_LT_OQ);
                if (_mm256_movemask_pd(cmp) == 0xF) {
                    _mm256_storeu_pd(row + j, _mm256_setzero_pd());
                    continue;
                }

                v = rigel::fast_exp_avx2(v);
                _mm256_storeu_pd(row + j, v);
                sum_v = _mm256_add_pd(sum_v, v);
            }
            // Horizontal sum of 4 lanes
            __m128d lo = _mm256_castpd256_pd128(sum_v);
            __m128d hi = _mm256_extractf128_pd(sum_v, 1);
            lo = _mm_add_pd(lo, hi);
            row_sum = _mm_cvtsd_f64(lo) + _mm_cvtsd_f64(_mm_unpackhi_pd(lo, lo));
        }
        // AVX2 scalar tail
        for (; j < k; ++j) {
            double x = row[j] - max_val;
            double e = rigel::fast_exp_scalar(x);
            row[j] = e;
            row_sum += e;
        }
#elif RIGEL_HAS_NEON
        {
            float64x2_t sum_v = vdupq_n_f64(0.0);
            const float64x2_t max_v = vdupq_n_f64(max_val);
            const float64x2_t cutoff_v = vdupq_n_f64(rigel::detail::EXP_CUTOFF);
            const float64x2_t zero_v = vdupq_n_f64(0.0);

            for (; j + 2 <= k; j += 2) {
                float64x2_t v = vld1q_f64(row + j);
                v = vsubq_f64(v, max_v);

                // Whole-vector early-zero skip: if BOTH lanes < cutoff, store zeros
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
        // NEON scalar tail
        for (; j < k; ++j) {
            double x = row[j] - max_val;
            double e = rigel::fast_exp_scalar(x);
            row[j] = e;
            row_sum += e;
        }
#else
        // Pure scalar path (no SIMD)
        for (; j < k; ++j) {
            double x = row[j] - max_val;
            double e = rigel::fast_exp_scalar(x);
            row[j] = e;
            row_sum += e;
        }
#endif

        // Normalize and accumulate column sums (fused)
        if (row_sum > 0.0 && std::isfinite(row_sum)) {
            double inv_sum = 1.0 / row_sum;
            for (int j = 0; j < k; ++j) {
                row[j] *= inv_sum;
                col_acc[j].add(row[j]);
            }
        } else {
            // col_acc: adding 0.0 is a no-op, skip
        }
    }

    // Flush column sums to em_totals
    for (int j = 0; j < k; ++j) {
        em_totals[cidx[j]] += col_acc[j].sum;
    }
}

// ================================================================
// Parallel E-step: task-based load balancing with Kahan reduction
// ================================================================
//
// Breaks large ECs into row-range tasks of bounded size
// (~ESTEP_TASK_WORK_TARGET/k rows), partitions tasks across threads
// by actual computational cost (n*k), and reduces with Kahan summation.
//
// Determinism: tasks are created and assigned in a fixed order;
// each thread processes its partition sequentially; cross-thread
// reduction sums in thread-index order with Kahan.  Result is
// fully deterministic.

struct EStepTask {
    int ec_idx;
    int row_start;
    int row_end;
    int64_t cost;
};

static void parallel_estep(
    const std::vector<EmEquivClass>& ec_data,
    const double*  log_weights,
    double*        em_totals,      // [n_components], zeroed by caller
    int            n_components,
    int            n_threads,
    rigel::EStepThreadPool* pool = nullptr)
{
    // 1. Break ECs into granular tasks for load balance
    std::vector<EStepTask> tasks;
    int64_t total_cost = 0;

    for (int i = 0; i < static_cast<int>(ec_data.size()); ++i) {
        int n = ec_data[i].n;
        int k = ec_data[i].k;
        int chunk_rows = std::max(1, ESTEP_TASK_WORK_TARGET / std::max(k, 1));

        for (int r = 0; r < n; r += chunk_rows) {
            int r_end = std::min(r + chunk_rows, n);
            int64_t cost = static_cast<int64_t>(r_end - r) * k;
            tasks.push_back({i, r, r_end, cost});
            total_cost += cost;
        }
    }

    // 2. Partition tasks by cumulative cost (greedy sequential)
    std::vector<int> thread_bounds(n_threads + 1, 0);
    {
        int cur_thread = 1;
        int64_t cur_cost = 0;
        int64_t target = total_cost / n_threads;

        for (size_t t = 0; t < tasks.size(); ++t) {
            cur_cost += tasks[t].cost;
            if (cur_cost >= target && cur_thread < n_threads) {
                thread_bounds[cur_thread++] = static_cast<int>(t + 1);
                cur_cost = 0;
            }
        }
        while (cur_thread <= n_threads) {
            thread_bounds[cur_thread++] = static_cast<int>(tasks.size());
        }
    }

    // 3. Thread-local em_totals buffers
    std::vector<std::vector<double>> local_totals(n_threads);
    for (int t = 0; t < n_threads; ++t) {
        local_totals[t].assign(static_cast<size_t>(n_components), 0.0);
    }

    // 4. Worker: process assigned tasks
    auto worker = [&](int tid) {
        int start_task = thread_bounds[tid];
        int end_task   = thread_bounds[tid + 1];
        double* my_totals = local_totals[tid].data();

        for (int t = start_task; t < end_task; ++t) {
            em_step_kernel_range(ec_data[tasks[t].ec_idx],
                                 tasks[t].row_start,
                                 tasks[t].row_end,
                                 log_weights,
                                 my_totals);
        }
    };

    // 5. Launch threads (use pool if available, else spawn/join)
    if (pool) {
        pool->run_parallel(worker);
    } else {
        std::vector<std::thread> threads;
        threads.reserve(n_threads - 1);
        for (int t = 1; t < n_threads; ++t) {
            threads.emplace_back(worker, t);
        }
        worker(0);
        for (auto& th : threads) th.join();
    }

    // 6. Deterministic Kahan reduction across threads
    for (int i = 0; i < n_components; ++i) {
        KahanAccumulator acc;
        for (int t = 0; t < n_threads; ++t) {
            acc.add(local_totals[t][i]);
        }
        em_totals[i] += acc.sum;
    }
}

// ================================================================
// MAP-EM step: theta → theta_new
// ================================================================
//
// Plain MAP-EM with per-component Dirichlet prior.
// gDNA strand symmetry is enforced via per-fragment log(0.5) in the likelihood
// (added at scoring time), so no M-step coupling is needed.

static void map_em_step(
    const double* theta,
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,
    const double* unambig_totals,
    const double* prior,
    double*       em_totals,    // zeroed then accumulated
    double*       theta_new,    // output: normalized
    int           n_components,
    int           estep_threads = 1,
    rigel::EStepThreadPool* pool = nullptr)
{
    // Compute log_weights = log(theta + epsilon) - log_eff_len
    std::vector<double> log_weights(static_cast<size_t>(n_components));
    for (int i = 0; i < n_components; ++i) {
        log_weights[i] = std::log(theta[i] + EM_LOG_EPSILON) - log_eff_len[i];
    }

    // Zero em_totals
    std::fill(em_totals, em_totals + n_components, 0.0);

    // E-step: accumulate posteriors
    if (estep_threads > 1) {
        parallel_estep(ec_data, log_weights.data(), em_totals,
                       n_components, estep_threads, pool);
    } else {
        for (const auto& ec : ec_data) {
            em_step_kernel_range(ec, 0, ec.n, log_weights.data(), em_totals);
        }
    }

    // Plain MAP-EM: theta_new = (unambig + em + prior)
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
    int           n_components,
    int           estep_threads = 1,
    rigel::EStepThreadPool* pool = nullptr)
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
    if (estep_threads > 1) {
        parallel_estep(ec_data, log_weights.data(), em_totals,
                       n_components, estep_threads, pool);
    } else {
        for (const auto& ec : ec_data) {
            em_step_kernel_range(ec, 0, ec.n, log_weights.data(), em_totals);
        }
    }

    // M-step: alpha_new = unambig_totals + em_totals + prior (unnormalized)
    for (int i = 0; i < n_components; ++i) {
        alpha_new[i] = unambig_totals[i] + em_totals[i] + prior[i];
    }
}

// ================================================================
// Coverage-weighted warm start + unified OVR prior
// ================================================================

static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,    // [n_components] 1.0 if eligible, 0.0 otherwise
    double        locus_gamma,       // calibration gDNA fraction ∈ [0,1]
    double        total_pseudocount, // C (total prior budget, default 1.0)
    int           gdna_idx,          // index of gDNA component (-1 if none)
    double*       prior_out,       // [n_components] output
    double*       theta_init_out,  // [n_components] output
    int           n_components)
{
    // Initialize theta_init from unambig_totals
    std::copy(unambig_totals, unambig_totals + n_components, theta_init_out);

    // Accumulate coverage totals from all equivalence classes
    std::vector<double> coverage_totals(static_cast<size_t>(n_components), 0.0);

    for (const auto& ec : ec_data) {
        const int n = ec.n;
        const int k = ec.k;
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

    // --- Unified OVR prior: budget-constrained distribution ---

    // Compute total RNA coverage (excluding gDNA and ineligible)
    double total_rna_coverage = 0.0;
    int n_rna_eligible = 0;
    for (int i = 0; i < n_components; ++i) {
        if (eligible[i] > 0.0 && i != gdna_idx) {
            total_rna_coverage += coverage_totals[i];
            ++n_rna_eligible;
        }
    }

    double rna_budget = (1.0 - locus_gamma) * total_pseudocount;
    double gdna_alpha = std::max(locus_gamma * total_pseudocount,
                                 EM_LOG_EPSILON);

    for (int i = 0; i < n_components; ++i) {
        if (eligible[i] <= 0.0) {
            prior_out[i] = 0.0;
        } else if (i == gdna_idx) {
            prior_out[i] = gdna_alpha;
        } else if (total_rna_coverage > 0.0) {
            prior_out[i] = std::max(
                rna_budget * coverage_totals[i] / total_rna_coverage,
                EM_LOG_EPSILON);
        } else if (n_rna_eligible > 0) {
            // No coverage data — distribute uniformly
            prior_out[i] = std::max(rna_budget / n_rna_eligible,
                                    EM_LOG_EPSILON);
        } else {
            prior_out[i] = EM_LOG_EPSILON;
        }
    }

    // --- gDNA warm-start override ---
    if (gdna_idx >= 0 && gdna_idx < n_components
        && eligible[gdna_idx] > 0.0)
    {
        double total_theta = 0.0;
        for (int i = 0; i < n_components; ++i)
            total_theta += theta_init_out[i];

        double others = total_theta - theta_init_out[gdna_idx];
        if (others > 0.0 && locus_gamma < 1.0 - 1e-10) {
            theta_init_out[gdna_idx] =
                (locus_gamma / std::max(1.0 - locus_gamma, 1e-10))
                * others;
        } else {
            theta_init_out[gdna_idx] =
                std::max(locus_gamma * total_theta, EM_LOG_EPSILON);
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
    int squarem_iterations = 0;  // number of SQUAREM iterations completed
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
    int           estep_threads = 1,
    rigel::EStepThreadPool* pool = nullptr)
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
    int completed_iterations = 0;

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
                      state1.data(), n_components, estep_threads, pool);

            vbem_step(state1.data(), ec_data, log_eff_len,
                      unambig_totals, prior, em_totals.data(),
                      state2.data(), n_components, estep_threads, pool);

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
                    // Clamp to recoverable floor: keeps components out of
                    // the digamma absorbing regime (prior alone is ~1e-6
                    // in mega-loci, deep in the death zone).
                    double floor_i = std::max(prior[i], VBEM_CLAMP_FLOOR);
                    if (state_extrap[i] < floor_i)
                        state_extrap[i] = floor_i;
                }
            }

            // Stabilisation step
            vbem_step(state_extrap.data(), ec_data, log_eff_len,
                      unambig_totals, prior, em_totals.data(),
                      state_new.data(), n_components, estep_threads, pool);

            // Floor-clamp + convergence check on normalized theta
            double sum_old = 0.0, sum_new = 0.0;
            for (size_t i = 0; i < nc; ++i) {
                // Clamp stabilisation output to recoverable floor —
                // prevents the absorbing barrier at tiny alpha values.
                double floor_i = std::max(prior[i], VBEM_CLAMP_FLOOR);
                if (state_new[i] < floor_i) {
                    state_new[i] = floor_i;
                }
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

            if (delta < convergence_delta) {
                completed_iterations = iter + 1;
                break;
            }
            completed_iterations = iter + 1;
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
            // Two EM steps
            map_em_step(state0.data(), ec_data, log_eff_len,
                        unambig_totals, prior, em_totals.data(),
                        state1.data(), n_components,
                        estep_threads, pool);

            map_em_step(state1.data(), ec_data, log_eff_len,
                        unambig_totals, prior, em_totals.data(),
                        state2.data(), n_components,
                        estep_threads, pool);

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
                        state_new.data(), n_components,
                        estep_threads, pool);

            // Convergence
            double delta = 0.0;
            for (size_t i = 0; i < nc; ++i) {
                delta += std::abs(state_new[i] - state0[i]);
            }

            std::swap(state0, state_new);

            if (delta < convergence_delta) {
                completed_iterations = iter + 1;
                break;
            }
            completed_iterations = iter + 1;
        }

        // theta = state0 (converged normalized theta)
        std::copy(state0.begin(), state0.end(), theta.begin());

        // alpha_out = unambig_totals + em_totals + prior
        for (size_t i = 0; i < nc; ++i) {
            alpha_out[i] = unambig_totals[i] + em_totals[i] + prior[i];
        }
    }

    return { std::move(theta), std::move(alpha_out), std::move(em_totals),
             completed_iterations };
}

// ================================================================
// Top-level entry point: run_locus_em_native()
// ================================================================
//
// Takes CSR per-locus data + config, returns (theta, alpha, em_totals).
// Replaces the entire body of AbundanceEstimator.run_locus_em() from
// bias correction through EM convergence.

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
    f64_1d unambig_totals_arr,
    f64_1d effective_lens,
    f64_1d prior_eligible,
    // Scalar config
    int    n_components,
    double total_pseudocount,
    int    max_iterations,
    double convergence_delta,
    bool   use_vbem,
    // Transcript count for prior warm-start
    int    n_transcripts)
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

    // Copy unambig_totals (we need a mutable copy)
    std::vector<double> unambig_totals(ut_ptr, ut_ptr + nc);

    // 1. Apply bias correction (uniform fast-path)
    if (n_candidates > 0) {
        apply_bias_correction_uniform(
            ll_ptr, ti_ptr, txs_ptr, txe_ptr, bp_ptr, n_candidates);
    }

    // 2. Handle empty locus
    if (n_units == 0 || n_candidates == 0) {
        // Count eligible components for uniform distribution
        int n_eligible = 0;
        for (size_t i = 0; i < nc; ++i) {
            if (pe_ptr[i] > 0.0) ++n_eligible;
        }
        double uniform_alpha = (n_eligible > 0)
            ? std::max(total_pseudocount / n_eligible, EM_LOG_EPSILON)
            : EM_LOG_EPSILON;

        std::vector<double> alpha(nc);
        double total = 0.0;
        for (size_t i = 0; i < nc; ++i) {
            double p = (pe_ptr[i] > 0.0) ? uniform_alpha : 0.0;
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

    // 5. Coverage-weighted warm start + unified OVR prior
    //    (run_locus_em_native: no gDNA component, all budget to RNA)
    std::vector<double> prior(nc);
    std::vector<double> theta_init(nc);
    compute_ovr_prior_and_warm_start(
        ec_data, unambig_totals.data(), pe_ptr,
        0.0,               // locus_gamma = 0 (no gDNA)
        total_pseudocount,
        -1,                // gdna_idx = -1 (no gDNA)
        prior.data(), theta_init.data(), n_components);

    // 6. Run SQUAREM
    EMResult result = run_squarem(
        ec_data, log_eff_len.data(), unambig_totals.data(),
        prior.data(), theta_init.data(),
        n_components, max_iterations, convergence_delta,
        use_vbem);

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
// Batch locus EM — single C++ call for all loci
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
    int n_components;  // n_t + 1
    int gdna_idx;      // = n_t     (single gDNA component)
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
    std::vector<double>   eligible;         // [n_components]

    // Local→global transcript mapping
    std::vector<int32_t>  local_to_global_t; // [n_t]
};

// Assign posteriors after EM convergence.
// Reimplements Python assign_locus_ambiguous() entirely in C++.
// Scatters results into the provided accumulator arrays.
//
// assignment_mode: ASSIGN_FRACTIONAL (0), ASSIGN_MAP (1), ASSIGN_SAMPLE (2)
// min_posterior: components with posterior < min_posterior are zeroed before
//   discrete (MAP/sample) assignment.  Ignored for fractional mode.
// rng: thread-local SplitMix64 instance (only used for sample mode).
static void assign_posteriors(
    const LocusSubProblem& sub,
    const double* theta,
    int assignment_mode,
    double min_posterior,
    SplitMix64& rng,
    // Output accumulators (accumulated across loci)
    double* em_counts_2d,          // [N_T, n_cols], row-major
    double* gdna_locus_counts_2d,  // [N_T, n_cols]
    double* posterior_sum,         // [N_T]
    double* n_assigned,            // [N_T]
    // Per-locus accumulation
    double& mrna_total,
    double& gdna_total,
    int N_T_TOTAL,  // total transcripts for bounds checking
    int n_cols,     // number of splice-strand columns (actual 2D stride)
    // --- Per-unit assignment output (nullable) ---
    int32_t* out_winner_tid,   // [total_units] or nullptr
    float*   out_winner_post,  // [total_units] or nullptr
    int16_t* out_n_candidates, // [total_units] or nullptr
    int      out_offset)       // write offset for this locus
{
    int n_t = sub.n_t;
    int nc  = sub.n_components;
    int gdna = sub.gdna_idx;
    int n_units = sub.n_local_units;
    const int32_t* local_to_global = sub.local_to_global_t.data();

    // Effective lengths are all 1.0, so log_eff_len = 0.
    // log_weights = log(theta + eps)
    std::vector<double> log_weights(nc);
    for (int c = 0; c < nc; ++c) {
        log_weights[c] = std::log(theta[c] + EM_LOG_EPSILON);
    }

    mrna_total = 0.0;
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

        // ---- Discrete assignment dispatch ----
        // For MAP/sample modes: threshold, renormalize, then select winner.
        // For fractional mode: use raw posteriors as-is.
        // The "weights" vector holds the final assignment weights (sum to 1).
        // In discrete modes, exactly one entry is 1.0 and the rest are 0.0.
        std::vector<double> weights(seg_len);

        int winner = -1;  // MAP/sample winner index (used by annotation output)

        if (assignment_mode == ASSIGN_FRACTIONAL) {
            // Traditional EM: scatter fractional posteriors
            for (int j = 0; j < seg_len; ++j) weights[j] = posteriors[j];
        } else {
            // MAP or sample: threshold then renormalize
            double renorm_sum = 0.0;
            for (int j = 0; j < seg_len; ++j) {
                if (posteriors[j] >= min_posterior) {
                    weights[j] = posteriors[j];
                    renorm_sum += posteriors[j];
                } else {
                    weights[j] = 0.0;
                }
            }
            if (renorm_sum > 0.0) {
                double inv = 1.0 / renorm_sum;
                for (int j = 0; j < seg_len; ++j) weights[j] *= inv;
            }

            // Find winner index
            winner = -1;
            if (assignment_mode == ASSIGN_MAP) {
                // Maximum a posteriori: pick the highest posterior
                double best = -1.0;
                for (int j = 0; j < seg_len; ++j) {
                    if (weights[j] > best) {
                        best = weights[j];
                        winner = j;
                    }
                }
            } else {
                // Sample: categorical draw from renormalized posteriors
                double u = rng.uniform();
                double cumulative = 0.0;
                for (int j = 0; j < seg_len; ++j) {
                    cumulative += weights[j];
                    if (u < cumulative) {
                        winner = j;
                        break;
                    }
                }
                // Edge case: rounding — assign to last non-zero
                if (winner < 0) {
                    for (int j = seg_len - 1; j >= 0; --j) {
                        if (weights[j] > 0.0) { winner = j; break; }
                    }
                }
            }

            // Zero everything, set winner to 1.0
            for (int j = 0; j < seg_len; ++j) weights[j] = 0.0;
            if (winner >= 0) weights[winner] = 1.0;
        }

        // Track max mRNA posterior (for diagnostics)

        // --- Per-unit annotation output ---
        if (out_winner_tid != nullptr) {
            int ann_winner = -1;
            if (assignment_mode == ASSIGN_FRACTIONAL) {
                // Fractional mode: find MAP winner for annotation display
                double best_post = -1.0;
                for (int j = 0; j < seg_len; ++j) {
                    if (posteriors[j] > best_post) {
                        best_post = posteriors[j];
                        ann_winner = j;
                    }
                }
            } else {
                ann_winner = winner;
            }

            int write_idx = out_offset + ui;
            out_n_candidates[write_idx] = static_cast<int16_t>(
                std::min(seg_len, static_cast<int>(INT16_MAX)));

            if (ann_winner >= 0) {
                int32_t comp = sub.t_indices[s + ann_winner];
                if (comp < n_t) {
                    out_winner_tid[write_idx] = local_to_global[comp];
                } else {
                    out_winner_tid[write_idx] = -2;  // gDNA
                }
                out_winner_post[write_idx] = static_cast<float>(
                    posteriors[ann_winner]);
            } else {
                out_winner_tid[write_idx] = -1;
                out_winner_post[write_idx] = 0.0f;
            }
        }

        // Scatter assignment weights
        for (int j = 0; j < seg_len; ++j) {
            int32_t comp = sub.t_indices[s + j];
            double p = weights[j];
            if (p == 0.0) continue;

            if (comp < n_t) {
                // Transcript (mRNA or synthetic nRNA — both handled the same way)
                int32_t global_t = local_to_global[comp];
                uint8_t col = sub.count_cols[s + j];
                if (global_t < 0 || global_t >= N_T_TOTAL || col >= n_cols) continue;
                em_counts_2d[global_t * n_cols + col] += p;
                mrna_total += p;

                // Confidence tracking (use original posteriors)
                posterior_sum[global_t] += posteriors[j] * posteriors[j];
                n_assigned[global_t] += posteriors[j];
            } else {
                // gDNA
                gdna_total += p;
            }
        }

        // gDNA locus attribution
        double gdna_unit_sum = 0.0;
        for (int j = 0; j < seg_len; ++j) {
            int32_t c = sub.t_indices[s + j];
            if (c == gdna) {
                gdna_unit_sum += weights[j];
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

// ================================================================
// Partition scatter functions — array-by-array global→per-locus scatter
// ================================================================

/// Build per-locus CSR offsets from global offsets and per-locus unit lists.
///
/// For each locus li with units [u0, u1, ...]:
///   partition_offsets[li][0] = 0
///   partition_offsets[li][k+1] = partition_offsets[li][k]
///                                + (g_offsets[u+1] - g_offsets[u])
///
/// Returns a Python list of int64 numpy arrays, one per locus.
static nb::list build_partition_offsets(
    i64_1d g_offsets,
    nb::list locus_units,
    int n_loci)
{
    const int64_t* goff = g_offsets.data();

    // Extract raw pointers from locus_units under GIL
    struct UnitInfo {
        const int32_t* data;
        int n;
    };
    std::vector<UnitInfo> unit_infos(n_loci);
    for (int li = 0; li < n_loci; ++li) {
        auto arr = nb::cast<i32_1d>(locus_units[li]);
        unit_infos[li] = {arr.data(), static_cast<int>(arr.shape(0))};
    }

    // Compute offsets (sequential — I/O bound)
    struct PartOffsets {
        int64_t* data;
        int n;  // n_units for this locus
    };
    std::vector<PartOffsets> results(n_loci);
    for (int li = 0; li < n_loci; ++li) {
        int n_u = unit_infos[li].n;
        const int32_t* u_arr = unit_infos[li].data;
        int64_t* out = new int64_t[n_u + 1];
        out[0] = 0;
        for (int k = 0; k < n_u; ++k) {
            int u = u_arr[k];
            out[k + 1] = out[k] + (goff[u + 1] - goff[u]);
        }
        results[li] = {out, n_u};
    }

    // Wrap results as numpy arrays
    nb::list result_list;
    for (int li = 0; li < n_loci; ++li) {
        int64_t* ptr = results[li].data;
        size_t shape[1] = {static_cast<size_t>(results[li].n + 1)};
        nb::capsule owner(ptr, [](void* p) noexcept {
            delete[] static_cast<int64_t*>(p);
        });
        result_list.append(
            nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(ptr, 1, shape, std::move(owner)));
    }
    return result_list;
}

/// Scatter per-candidate data from global CSR into per-locus arrays.
///
/// For each locus, copies candidate segments:
///   For unit k with global index u = locus_units[li][k]:
///     memcpy(dst + p_off[k], src + g_off[u], (g_off[u+1]-g_off[u]) * sizeof(T))
///
/// Returns a Python list of numpy arrays, one per locus.
template <typename T>
static nb::list scatter_candidates_impl(
    nb::ndarray<const T, nb::ndim<1>, nb::c_contig> global_arr,
    i64_1d g_offsets,
    nb::list locus_units,
    nb::list partition_offsets,
    int n_loci)
{
    const T* src = global_arr.data();
    const int64_t* goff = g_offsets.data();

    // Extract pointers under GIL
    struct LInfo {
        const int32_t* units;
        const int64_t* p_off;
        int n_units;
        int64_t n_candidates;
    };
    std::vector<LInfo> infos(n_loci);
    for (int li = 0; li < n_loci; ++li) {
        auto u_arr = nb::cast<i32_1d>(locus_units[li]);
        auto p_arr = nb::cast<i64_1d>(partition_offsets[li]);
        int n_u = static_cast<int>(u_arr.shape(0));
        const int64_t* poff = p_arr.data();
        infos[li] = {u_arr.data(), poff, n_u, poff[n_u]};
    }

    // Scatter (sequential)
    struct Result { T* data; int64_t n; };
    std::vector<Result> results(n_loci);
    for (int li = 0; li < n_loci; ++li) {
        const auto& info = infos[li];
        T* dst = new T[info.n_candidates > 0 ? info.n_candidates : 1];
        for (int k = 0; k < info.n_units; ++k) {
            int u = info.units[k];
            int64_t g_start = goff[u];
            int64_t seg_len = goff[u + 1] - g_start;
            if (seg_len > 0) {
                std::memcpy(dst + info.p_off[k],
                            src + g_start,
                            static_cast<size_t>(seg_len) * sizeof(T));
            }
        }
        results[li] = {dst, info.n_candidates};
    }

    // Wrap as numpy arrays
    nb::list result_list;
    for (int li = 0; li < n_loci; ++li) {
        T* ptr = results[li].data;
        size_t shape[1] = {static_cast<size_t>(results[li].n)};
        nb::capsule owner(ptr, [](void* p) noexcept {
            delete[] static_cast<T*>(p);
        });
        result_list.append(
            nb::ndarray<nb::numpy, T, nb::ndim<1>>(ptr, 1, shape, std::move(owner)));
    }
    return result_list;
}

/// Scatter per-unit data from global array into per-locus arrays.
///
/// For each locus, gathers elements: dst[k] = src[units[k]]
///
/// Returns a Python list of numpy arrays, one per locus.
template <typename T>
static nb::list scatter_units_impl(
    nb::ndarray<const T, nb::ndim<1>, nb::c_contig> global_arr,
    nb::list locus_units,
    int n_loci)
{
    const T* src = global_arr.data();

    // Extract pointers under GIL
    struct UInfo { const int32_t* data; int n; };
    std::vector<UInfo> infos(n_loci);
    for (int li = 0; li < n_loci; ++li) {
        auto arr = nb::cast<i32_1d>(locus_units[li]);
        infos[li] = {arr.data(), static_cast<int>(arr.shape(0))};
    }

    // Gather (sequential)
    struct Result { T* data; int n; };
    std::vector<Result> results(n_loci);
    for (int li = 0; li < n_loci; ++li) {
        int n_u = infos[li].n;
        const int32_t* u_arr = infos[li].data;
        T* dst = new T[n_u > 0 ? n_u : 1];
        for (int k = 0; k < n_u; ++k) {
            dst[k] = src[u_arr[k]];
        }
        results[li] = {dst, n_u};
    }

    // Wrap as numpy arrays
    nb::list result_list;
    for (int li = 0; li < n_loci; ++li) {
        T* ptr = results[li].data;
        size_t shape[1] = {static_cast<size_t>(results[li].n)};
        nb::capsule owner(ptr, [](void* p) noexcept {
            delete[] static_cast<T*>(p);
        });
        result_list.append(
            nb::ndarray<nb::numpy, T, nb::ndim<1>>(ptr, 1, shape, std::move(owner)));
    }
    return result_list;
}

// ================================================================
// PartitionView — per-locus CSR data for partition-native EM
// ================================================================

struct PartitionView {
    // Per-locus CSR data (contiguous, 0-indexed)
    const int64_t* offsets;
    const int32_t* t_indices;
    const double*  log_liks;
    const double*  coverage_wts;
    const int32_t* tx_starts;
    const int32_t* tx_ends;
    const uint8_t* count_cols;
    const uint8_t* is_spliced;
    const double*  gdna_log_liks;
    const int32_t* genomic_footprints;
    const int32_t* locus_t_indices;
    const uint8_t* locus_count_cols;
    int     n_units;
    int64_t n_candidates;

    // Locus transcript membership (from Locus.transcript_indices)
    const int32_t* transcript_indices;
    int n_transcripts;
};

// Extract per-locus sub-problem from a PartitionView.
//
// Remaps global transcript indices to local component indices [0, n_t),
// appends a gDNA component (index n_t) for unspliced units, and sorts
// candidates within each unit by component index (required by the
// equivalence-class builder downstream).
//
// Candidates are written to pre-allocated output arrays via a write cursor
// to avoid dynamic allocation in the inner loop. A reusable sort buffer
// (std::vector<LocalCandidate>) is shared across all units — resize() is
// a no-op when existing capacity is sufficient, so after the first unit
// there are effectively zero allocations.
static void extract_locus_sub_problem_from_partition(
    LocusSubProblem& sub,
    const PartitionView& pv,
    double locus_gamma,
    int64_t gdna_span,
    const double*  all_unambig_row_sums,
    const int64_t* all_t_lengths,
    int32_t* local_map, int local_map_size)
{
    int n_t = pv.n_transcripts;
    int n_u = pv.n_units;
    const int32_t* t_arr = pv.transcript_indices;

    sub.n_t = n_t;
    sub.n_local_units = n_u;
    sub.n_components = n_t + 1;
    sub.gdna_idx = n_t;
    int nc = sub.n_components;

    // --- Build global→local mapping ---
    int max_global = 0;
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        if (gt + 1 > max_global) max_global = gt + 1;
    }
    if (max_global > local_map_size) max_global = local_map_size;

    for (int i = 0; i < max_global; ++i) local_map[i] = -1;
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        if (gt >= 0 && gt < local_map_size) local_map[gt] = i;
    }

    // --- Pre-allocate output arrays to worst-case size ---
    // Maximum output candidates = input RNA candidates + one gDNA per unit.
    size_t max_out = static_cast<size_t>(pv.n_candidates)
                   + static_cast<size_t>(n_u);

    sub.t_indices.resize(max_out);
    sub.log_liks.resize(max_out);
    sub.coverage_wts.resize(max_out);
    sub.tx_starts.resize(max_out);
    sub.tx_ends.resize(max_out);
    sub.count_cols.resize(max_out);

    sub.offsets.resize(n_u + 1);
    sub.offsets[0] = 0;

    sub.locus_t_arr.resize(n_u);
    sub.locus_ct_arr.resize(n_u);

    // Reusable sort buffer — persists across units, only grows.
    std::vector<LocalCandidate> sort_buf;

    size_t cursor = 0;

    for (int ui = 0; ui < n_u; ++ui) {
        auto p_start = pv.offsets[ui];
        auto p_end   = pv.offsets[ui + 1];
        int width_in = static_cast<int>(p_end - p_start);

        sub.locus_t_arr[ui] = pv.locus_t_indices[ui];
        sub.locus_ct_arr[ui] = pv.locus_count_cols[ui];

        // Determine if this unit gets a gDNA candidate
        bool is_spliced = (pv.is_spliced[ui] != 0);
        double gdna_ll = pv.gdna_log_liks[ui];
        bool has_gdna = (!is_spliced && std::isfinite(gdna_ll));
        int width_out = width_in + (has_gdna ? 1 : 0);

        // Fill sort buffer with remapped candidates
        sort_buf.resize(width_out);
        int k = 0;

        for (auto j = p_start; j < p_end; ++j) {
            int32_t global_t = pv.t_indices[j];
            if (global_t < 0 || global_t >= local_map_size) continue;
            int32_t local = local_map[global_t];
            if (local < 0 || local >= nc) continue;

            sort_buf[k++] = {local, pv.log_liks[j], pv.coverage_wts[j],
                             pv.tx_starts[j], pv.tx_ends[j],
                             pv.count_cols[j]};
        }

        // Append gDNA candidate (component = n_t, always the largest index)
        if (has_gdna) {
            int32_t footprint = pv.genomic_footprints[ui];
            sort_buf[k++] = {sub.gdna_idx, gdna_ll, 1.0,
                             0, footprint, 0};
        }

        // Trim to actual count (candidates may have been skipped above)
        int actual = k;

        // Sort by local component index (required by equiv-class builder).
        // For n <= 1 this is a no-op. std::sort handles small n efficiently
        // via insertion sort internally.
        if (actual > 1) {
            std::sort(sort_buf.begin(), sort_buf.begin() + actual,
                      [](const LocalCandidate& a, const LocalCandidate& b) {
                          return a.local_comp < b.local_comp;
                      });
        }

        // Write sorted candidates directly to pre-allocated output
        for (int i = 0; i < actual; ++i) {
            const auto& c = sort_buf[i];
            sub.t_indices[cursor]    = c.local_comp;
            sub.log_liks[cursor]     = c.log_lik;
            sub.coverage_wts[cursor] = c.cov_wt;
            sub.tx_starts[cursor]    = c.tx_start;
            sub.tx_ends[cursor]      = c.tx_end;
            sub.count_cols[cursor]   = c.count_col;
            ++cursor;
        }

        sub.offsets[ui + 1] = static_cast<int64_t>(cursor);
    }

    // Trim output arrays to actual size
    sub.t_indices.resize(cursor);
    sub.log_liks.resize(cursor);
    sub.coverage_wts.resize(cursor);
    sub.tx_starts.resize(cursor);
    sub.tx_ends.resize(cursor);
    sub.count_cols.resize(cursor);

    // --- Build per-component arrays ---
    sub.local_to_global_t.resize(n_t);
    for (int i = 0; i < n_t; ++i) sub.local_to_global_t[i] = t_arr[i];

    sub.unambig_totals.assign(nc, 0.0);
    for (int i = 0; i < n_t; ++i) {
        sub.unambig_totals[i] = all_unambig_row_sums[t_arr[i]];
    }

    sub.bias_profiles.resize(nc);
    for (int i = 0; i < n_t; ++i) {
        sub.bias_profiles[i] = all_t_lengths[t_arr[i]];
    }
    sub.bias_profiles[sub.gdna_idx] = gdna_span;

    sub.prior.assign(nc, EM_PRIOR_EPSILON);
    if (locus_gamma == 0.0) {
        sub.prior[sub.gdna_idx] = 0.0;
    }

    for (int c = 0; c < nc; ++c) {
        if (sub.prior[c] == 0.0) sub.unambig_totals[c] = 0.0;
    }

    sub.eligible.resize(nc);
    for (int c = 0; c < nc; ++c) {
        sub.eligible[c] = (sub.prior[c] > 0.0) ? 1.0 : 0.0;
    }

    // Clean up local_map scratch for next call
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        if (gt >= 0 && gt < local_map_size) local_map[gt] = -1;
    }
}

// ----------------------------------------------------------------
// Partition-native batch EM entry point
// ----------------------------------------------------------------

static std::tuple<
    double,  // total_gdna_em
    nb::ndarray<nb::numpy, double, nb::ndim<1>>,  // locus_mrna[n_loci]
    nb::ndarray<nb::numpy, double, nb::ndim<1>>,  // locus_gdna[n_loci]
    nb::list,  // locus_stats
    nb::object,  // out_winner_tid (ndarray or None)
    nb::object,  // out_winner_post (ndarray or None)
    nb::object   // out_n_candidates (ndarray or None)
>
batch_locus_em_partitioned(
    // Per-locus partition data (list of 12-tuples)
    nb::list partition_tuples,
    // Per-locus transcript membership (list of int32[])
    nb::list locus_transcript_indices,
    // Per-locus scalars
    f64_1d   locus_gammas,
    i64_1d   gdna_spans,
    // Per-transcript globals
    f64_2d   unambig_counts,
    i64_1d   t_lengths_arr,
    // Mutable output accumulators
    f64_2d_mut em_counts_out,
    f64_2d_mut gdna_locus_counts_out,
    f64_1d_mut posterior_sum_out,
    f64_1d_mut n_assigned_out,
    // EM config
    int    max_iterations,
    double convergence_delta,
    double total_pseudocount,
    bool   use_vbem,
    int    assignment_mode,
    double assignment_min_posterior,
    uint64_t rng_seed,
    int    n_transcripts_total,
    int    n_splice_strand_cols,
    int    n_threads,
    bool   emit_locus_stats,
    bool   emit_assignments)
{
    int n_loci = static_cast<int>(nb::len(partition_tuples));
    int N_T = n_transcripts_total;
    int N_COLS = n_splice_strand_cols;

    // --- Under GIL: extract PartitionViews from tuples ---
    std::vector<PartitionView> views(n_loci);
    for (int i = 0; i < n_loci; ++i) {
        nb::tuple tup = nb::borrow<nb::tuple>(partition_tuples[i]);
        auto& v = views[i];
        auto off_arr = nb::cast<i64_1d>(tup[0]);
        v.offsets          = off_arr.data();
        v.t_indices        = nb::cast<i32_1d>(tup[1]).data();
        v.log_liks         = nb::cast<f64_1d>(tup[2]).data();
        v.coverage_wts     = nb::cast<f64_1d>(tup[3]).data();
        v.tx_starts        = nb::cast<i32_1d>(tup[4]).data();
        v.tx_ends          = nb::cast<i32_1d>(tup[5]).data();
        v.count_cols       = nb::cast<u8_1d>(tup[6]).data();
        v.is_spliced       = nb::cast<u8_1d>(tup[7]).data();
        v.gdna_log_liks    = nb::cast<f64_1d>(tup[8]).data();
        v.genomic_footprints = nb::cast<i32_1d>(tup[9]).data();
        v.locus_t_indices  = nb::cast<i32_1d>(tup[10]).data();
        v.locus_count_cols = nb::cast<u8_1d>(tup[11]).data();
        v.n_units = static_cast<int>(off_arr.shape(0)) - 1;
        v.n_candidates = v.offsets[v.n_units];

        auto t_arr = nb::cast<i32_1d>(locus_transcript_indices[i]);
        v.transcript_indices = t_arr.data();
        v.n_transcripts = static_cast<int>(t_arr.shape(0));
    }

    const double*   lg_ptr = locus_gammas.data();
    const int64_t*  gs_ptr = gdna_spans.data();
    const double*   uac    = unambig_counts.data();
    const int64_t*  tl_ptr = t_lengths_arr.data();

    double* em_out    = em_counts_out.data();
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

    int local_map_size = N_T + 1;

    // --- Per-unit assignment output arrays ---
    int total_units = 0;
    std::vector<int> locus_write_offsets(n_loci, 0);
    for (int i = 0; i < n_loci; ++i) {
        locus_write_offsets[i] = total_units;
        total_units += views[i].n_units;
    }

    int32_t* out_tid_ptr   = nullptr;
    float*   out_post_ptr  = nullptr;
    int16_t* out_ncand_ptr = nullptr;

    std::vector<int32_t> out_tid_vec;
    std::vector<float>   out_post_vec;
    std::vector<int16_t> out_ncand_vec;

    if (emit_assignments && total_units > 0) {
        out_tid_vec.assign(total_units, -1);
        out_post_vec.assign(total_units, 0.0f);
        out_ncand_vec.assign(total_units, 0);
        out_tid_ptr   = out_tid_vec.data();
        out_post_ptr  = out_post_vec.data();
        out_ncand_ptr = out_ncand_vec.data();
    }

    std::vector<double> locus_mrna_vec(n_loci, 0.0);
    std::vector<double> locus_gdna_vec(n_loci, 0.0);
    double* locus_mrna_data = locus_mrna_vec.data();
    double* locus_gdna_data = locus_gdna_vec.data();

    std::vector<LocusProfile> locus_profiles(
        emit_locus_stats ? static_cast<size_t>(n_loci) : 0);

    std::atomic<double> total_gdna_em{0.0};

    int actual_threads = n_threads;
    if (actual_threads <= 0) {
        int hw = static_cast<int>(std::thread::hardware_concurrency());
        actual_threads = (hw > 0) ? hw : 1;
    }
    if (actual_threads < 1) actual_threads = 1;

    {
        nb::gil_scoped_release release;

        // --- Scheduling ---
        std::vector<int64_t> locus_work(static_cast<size_t>(n_loci));
        int64_t total_work = 0;
        for (int li = 0; li < n_loci; ++li) {
            locus_work[li] = static_cast<int64_t>(views[li].n_transcripts) *
                             static_cast<int64_t>(views[li].n_units);
            total_work += locus_work[li];
        }

        std::vector<int> locus_order(static_cast<size_t>(n_loci));
        std::iota(locus_order.begin(), locus_order.end(), 0);
        std::sort(locus_order.begin(), locus_order.end(),
                  [&](int a, int b) { return locus_work[a] > locus_work[b]; });

        int64_t fair_share = (actual_threads > 1 && total_work > 0)
            ? total_work / actual_threads
            : total_work + 1;

        int mega_end = 0;
        for (int i = 0; i < n_loci; ++i) {
            if (locus_work[locus_order[i]] >= fair_share)
                ++mega_end;
            else
                break;
        }

        // Lambda: process one locus from partition data
        auto process_locus = [&](int li, int estep_thr,
                                 LocusSubProblem& sub,
                                 std::vector<int32_t>& local_map_vec,
                                 bool mega,
                                 rigel::EStepThreadPool* pool = nullptr) -> double {
            auto locus_t0 = hrclock::now();
            const auto& pv = views[li];
            int n_t = pv.n_transcripts;
            int n_u = pv.n_units;

            if (n_u == 0) {
                locus_mrna_data[li] = 0.0;
                locus_gdna_data[li] = 0.0;
                return 0.0;
            }

            // 1. Extract sub-problem from partition
            auto t1 = hrclock::now();
            extract_locus_sub_problem_from_partition(
                sub, pv, lg_ptr[li], gs_ptr[li],
                unambig_row_sums.data(), tl_ptr,
                local_map_vec.data(), local_map_size);
            auto t2 = hrclock::now();

            int nc = sub.n_components;
            size_t n_candidates = sub.t_indices.size();
            int n_local_units = sub.n_local_units;

            // 2. Apply bias correction
            if (n_candidates > 0) {
                apply_bias_correction_uniform(
                    sub.log_liks.data(),
                    sub.t_indices.data(),
                    sub.tx_starts.data(),
                    sub.tx_ends.data(),
                    sub.bias_profiles.data(),
                    n_candidates);
            }
            auto t3 = hrclock::now();

            // 3. Handle empty sub-problem
            if (n_local_units == 0 || n_candidates == 0) {
                locus_mrna_data[li] = 0.0;
                locus_gdna_data[li] = 0.0;
                return 0.0;
            }

            // 4. Log effective lengths (all 1.0 → log = 0.0)
            std::vector<double> log_eff_len(nc, 0.0);

            // 5. Build equivalence classes
            auto ec_data = build_equiv_classes(
                sub.offsets.data(),
                sub.t_indices.data(),
                sub.log_liks.data(),
                sub.coverage_wts.data(),
                n_local_units);
            auto t4 = hrclock::now();

            // 6. OVR prior + coverage-weighted warm start
            std::vector<double> prior(nc);
            std::vector<double> theta_init(nc);
            compute_ovr_prior_and_warm_start(
                ec_data,
                sub.unambig_totals.data(),
                sub.eligible.data(),
                lg_ptr[li], total_pseudocount, sub.gdna_idx,
                prior.data(), theta_init.data(), nc);
            auto t5 = hrclock::now();

            // 7. SQUAREM
            EMResult result = run_squarem(
                ec_data, log_eff_len.data(),
                sub.unambig_totals.data(),
                prior.data(),
                theta_init.data(),
                nc, max_iterations, convergence_delta,
                use_vbem,
                estep_thr, pool);
            auto t6 = hrclock::now();

            // 8. Assign posteriors
            SplitMix64 locus_rng(rng_seed ^ (static_cast<uint64_t>(li) * 0x9e3779b97f4a7c15ULL));
            double locus_mrna = 0.0, locus_gdna = 0.0;
            assign_posteriors(
                sub, result.theta.data(),
                assignment_mode, assignment_min_posterior, locus_rng,
                em_out, gdna_out,
                psum_out, nass_out,
                locus_mrna, locus_gdna,
                N_T, N_COLS,
                out_tid_ptr, out_post_ptr, out_ncand_ptr,
                locus_write_offsets[li]);
            auto t7 = hrclock::now();

            locus_mrna_data[li] = locus_mrna;
            locus_gdna_data[li] = locus_gdna;

            if (emit_locus_stats) {
                auto us = [](auto a, auto b) {
                    return std::chrono::duration<double, std::micro>(b - a).count();
                };
                LocusProfile& prof = locus_profiles[li];
                prof.locus_idx = li;
                prof.n_transcripts = n_t;
                prof.n_units = n_u;
                prof.n_components = nc;
                prof.n_equiv_classes = static_cast<int>(ec_data.size());
                prof.squarem_iterations = result.squarem_iterations;
                prof.estep_threads_used = estep_thr;
                prof.is_mega_locus = mega;

                int64_t total_elems = 0;
                int max_k = 0, max_n = 0;
                for (const auto& ec : ec_data) {
                    total_elems += static_cast<int64_t>(ec.n) * ec.k;
                    if (ec.k > max_k) max_k = ec.k;
                    if (ec.n > max_n) max_n = ec.n;
                }
                prof.ec_total_elements = total_elems;
                prof.max_ec_width = max_k;
                prof.max_ec_depth = max_n;
                prof.digamma_calls_per_estep = use_vbem ? nc : 0;
                prof.extract_us = us(locus_t0, t2);
                prof.bias_us = us(t2, t3);
                prof.build_ec_us = us(t3, t4);
                prof.warm_start_us = us(t4, t5);
                prof.squarem_us = us(t5, t6);
                prof.assign_us = us(t6, t7);
                prof.total_us = us(locus_t0, t7);
            }

            return locus_gdna;
        }; // end process_locus

        // ---- Phase 1: mega-loci ----
        {
            LocusSubProblem sub;
            std::vector<int32_t> local_map_vec(local_map_size, -1);

            std::unique_ptr<rigel::EStepThreadPool> pool;
            if (actual_threads > 1 && mega_end > 0) {
                pool = std::make_unique<rigel::EStepThreadPool>(actual_threads);
            }

            for (int i = 0; i < mega_end; ++i) {
                int li = locus_order[i];
                double gdna = process_locus(li, actual_threads, sub, local_map_vec,
                                            true, pool.get());
                double prev = total_gdna_em.load(std::memory_order_relaxed);
                while (!total_gdna_em.compare_exchange_weak(
                    prev, prev + gdna,
                    std::memory_order_relaxed, std::memory_order_relaxed)) {}
            }
        }

        // ---- Phase 2: work-steal normal loci ----
        int n_phase2 = n_loci - mega_end;
        if (n_phase2 > 0) {
            constexpr int CHUNK_SIZE = 16;
            std::atomic<int> next_idx{0};

            auto worker_fn = [&]() {
                LocusSubProblem sub;
                std::vector<int32_t> local_map_vec(local_map_size, -1);
                double local_gdna = 0.0;

                for (;;) {
                    int chunk_start = next_idx.fetch_add(CHUNK_SIZE,
                        std::memory_order_relaxed);
                    if (chunk_start >= n_phase2) break;
                    int chunk_end = std::min(chunk_start + CHUNK_SIZE, n_phase2);
                    for (int idx = chunk_start; idx < chunk_end; ++idx) {
                        int li = locus_order[mega_end + idx];
                        local_gdna += process_locus(li, 1, sub, local_map_vec,
                                                    false);
                    }
                }

                double prev = total_gdna_em.load(std::memory_order_relaxed);
                while (!total_gdna_em.compare_exchange_weak(
                    prev, prev + local_gdna,
                    std::memory_order_relaxed, std::memory_order_relaxed)) {}
            };

            if (actual_threads <= 1) {
                worker_fn();
            } else {
                std::vector<std::thread> threads;
                threads.reserve(actual_threads);
                for (int t = 0; t < actual_threads; ++t)
                    threads.emplace_back(worker_fn);
                for (auto& th : threads)
                    th.join();
            }
        }
    } // end gil_scoped_release

    double total_gdna_em_val = total_gdna_em.load(std::memory_order_relaxed);

    size_t shape[1] = {static_cast<size_t>(n_loci)};
    auto* mrna_copy = new double[n_loci];
    auto* gdna_copy = new double[n_loci];
    std::memcpy(mrna_copy, locus_mrna_vec.data(), n_loci * sizeof(double));
    std::memcpy(gdna_copy, locus_gdna_vec.data(), n_loci * sizeof(double));

    nb::capsule mrna_owner(mrna_copy, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    nb::capsule gdna_owner(gdna_copy, [](void* p) noexcept { delete[] static_cast<double*>(p); });

    nb::list stats_list;
    if (emit_locus_stats) {
        std::vector<size_t> sorted_idx(locus_profiles.size());
        std::iota(sorted_idx.begin(), sorted_idx.end(), 0);
        std::sort(sorted_idx.begin(), sorted_idx.end(),
                  [&](size_t a, size_t b) {
                      return locus_profiles[a].total_us >
                             locus_profiles[b].total_us;
                  });

        for (size_t si : sorted_idx) {
            const auto& p = locus_profiles[si];
            if (p.total_us == 0.0 && p.n_units == 0) continue;

            nb::dict d;
            d["locus_idx"] = p.locus_idx;
            d["n_transcripts"] = p.n_transcripts;
            d["n_units"] = p.n_units;
            d["n_components"] = p.n_components;
            d["n_equiv_classes"] = p.n_equiv_classes;
            d["ec_total_elements"] = p.ec_total_elements;
            d["max_ec_width"] = p.max_ec_width;
            d["max_ec_depth"] = p.max_ec_depth;
            d["squarem_iterations"] = p.squarem_iterations;
            d["estep_threads_used"] = p.estep_threads_used;
            d["is_mega_locus"] = p.is_mega_locus;
            d["digamma_calls_per_estep"] = p.digamma_calls_per_estep;
            d["extract_us"] = p.extract_us;
            d["bias_us"] = p.bias_us;
            d["build_ec_us"] = p.build_ec_us;
            d["warm_start_us"] = p.warm_start_us;
            d["squarem_us"] = p.squarem_us;
            d["assign_us"] = p.assign_us;
            d["total_us"] = p.total_us;
            stats_list.append(d);
        }
    }

    return std::make_tuple(
        total_gdna_em_val,
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(
            mrna_copy, 1, shape, std::move(mrna_owner)),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(
            gdna_copy, 1, shape, std::move(gdna_owner)),
        stats_list,
        [&]() -> nb::object {
            if (!emit_assignments || total_units <= 0) return nb::none();
            size_t u_shape[1] = {static_cast<size_t>(total_units)};
            auto* p = new int32_t[total_units];
            std::memcpy(p, out_tid_vec.data(), total_units * sizeof(int32_t));
            nb::capsule own(p, [](void* x) noexcept { delete[] static_cast<int32_t*>(x); });
            return nb::cast(nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(p, 1, u_shape, std::move(own)));
        }(),
        [&]() -> nb::object {
            if (!emit_assignments || total_units <= 0) return nb::none();
            size_t u_shape[1] = {static_cast<size_t>(total_units)};
            auto* p = new float[total_units];
            std::memcpy(p, out_post_vec.data(), total_units * sizeof(float));
            nb::capsule own(p, [](void* x) noexcept { delete[] static_cast<float*>(x); });
            return nb::cast(nb::ndarray<nb::numpy, float, nb::ndim<1>>(p, 1, u_shape, std::move(own)));
        }(),
        [&]() -> nb::object {
            if (!emit_assignments || total_units <= 0) return nb::none();
            size_t u_shape[1] = {static_cast<size_t>(total_units)};
            auto* p = new int16_t[total_units];
            std::memcpy(p, out_ncand_vec.data(), total_units * sizeof(int16_t));
            nb::capsule own(p, [](void* x) noexcept { delete[] static_cast<int16_t*>(x); });
            return nb::cast(nb::ndarray<nb::numpy, int16_t, nb::ndim<1>>(p, 1, u_shape, std::move(own)));
        }()
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
///   t_indices    — int32[n_candidates]: candidate transcript indices
///   n_transcripts — int32: number of transcripts
///
/// Returns a tuple of:
///   n_components — int32: number of connected components
///   comp_t_offsets, comp_t_flat — CSR of transcript indices per component
///   comp_u_offsets, comp_u_flat — CSR of unit indices per component
static nb::tuple connected_components_native(
    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig>  offsets_arr,
    nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig>  t_indices_arr,
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

        // Collect all transcript indices from this unit.
        // All candidates are now direct transcript indices.
        int32_t first_t = -1;
        for (int64_t j = start; j < end; ++j) {
            int32_t t = t_idx[j];
            if (t >= 0 && t < n_transcripts) {
                active[t] = true;
                if (first_t < 0) {
                    first_t = t;
                } else {
                    uf.unite(first_t, t);
                }
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
              "Also provides batch_locus_em_partitioned() which replaces the entire\n"
              "per-locus Python for-loop with a single C++ call.";

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
          nb::arg("total_pseudocount"),
          nb::arg("max_iterations"),
          nb::arg("convergence_delta"),
          nb::arg("use_vbem"),
          nb::arg("n_transcripts"),
          "Run EM for a single locus sub-problem.\n\n"
          "Takes CSR per-locus data + config, returns (theta, alpha, em_totals).\n"
          "Replaces the Python EM hot path with a single C++ call.");

    // ---- Partition scatter functions ----
    m.def("build_partition_offsets", &build_partition_offsets,
          nb::arg("g_offsets"),
          nb::arg("locus_units"),
          nb::arg("n_loci"),
          "Build per-locus CSR offsets from global offsets and locus unit lists.");

    m.def("scatter_candidates_f64",
          &scatter_candidates_impl<double>,
          nb::arg("global_arr"), nb::arg("g_offsets"),
          nb::arg("locus_units"), nb::arg("partition_offsets"),
          nb::arg("n_loci"),
          "Scatter per-candidate float64 array into per-locus arrays.");
    m.def("scatter_candidates_i32",
          &scatter_candidates_impl<int32_t>,
          nb::arg("global_arr"), nb::arg("g_offsets"),
          nb::arg("locus_units"), nb::arg("partition_offsets"),
          nb::arg("n_loci"),
          "Scatter per-candidate int32 array into per-locus arrays.");
    m.def("scatter_candidates_u8",
          &scatter_candidates_impl<uint8_t>,
          nb::arg("global_arr"), nb::arg("g_offsets"),
          nb::arg("locus_units"), nb::arg("partition_offsets"),
          nb::arg("n_loci"),
          "Scatter per-candidate uint8 array into per-locus arrays.");

    m.def("scatter_units_f64",
          &scatter_units_impl<double>,
          nb::arg("global_arr"), nb::arg("locus_units"), nb::arg("n_loci"),
          "Scatter per-unit float64 array into per-locus arrays.");
    m.def("scatter_units_i32",
          &scatter_units_impl<int32_t>,
          nb::arg("global_arr"), nb::arg("locus_units"), nb::arg("n_loci"),
          "Scatter per-unit int32 array into per-locus arrays.");
    m.def("scatter_units_u8",
          &scatter_units_impl<uint8_t>,
          nb::arg("global_arr"), nb::arg("locus_units"), nb::arg("n_loci"),
          "Scatter per-unit uint8 array into per-locus arrays.");
    m.def("scatter_units_i64",
          &scatter_units_impl<int64_t>,
          nb::arg("global_arr"), nb::arg("locus_units"), nb::arg("n_loci"),
          "Scatter per-unit int64 array into per-locus arrays.");

    // ---- Partition-native batch EM ----
    m.def("batch_locus_em_partitioned", &batch_locus_em_partitioned,
          nb::arg("partition_tuples"),
          nb::arg("locus_transcript_indices"),
          nb::arg("locus_gammas"),
          nb::arg("gdna_spans"),
          nb::arg("unambig_counts"),
          nb::arg("t_lengths"),
          nb::arg("em_counts_out"),
          nb::arg("gdna_locus_counts_out"),
          nb::arg("posterior_sum_out"),
          nb::arg("n_assigned_out"),
          nb::arg("max_iterations"),
          nb::arg("convergence_delta"),
          nb::arg("total_pseudocount"),
          nb::arg("use_vbem"),
          nb::arg("assignment_mode"),
          nb::arg("assignment_min_posterior"),
          nb::arg("rng_seed"),
          nb::arg("n_transcripts_total"),
          nb::arg("n_splice_strand_cols"),
          nb::arg("n_threads") = 0,
          nb::arg("emit_locus_stats") = false,
          nb::arg("emit_assignments") = false,
          "Run locus EM from per-locus partition data.\n\n"
          "Accepts a list of 12-tuples (one per locus) containing partition\n"
          "arrays, plus a list of transcript index arrays per locus.\n"
          "Returns (total_gdna_em, locus_mrna, locus_gdna, locus_stats,\n"
          " out_winner_tid, out_winner_post, out_n_candidates).");

    m.def("connected_components", &connected_components_native,
          nb::arg("offsets"),
          nb::arg("t_indices"),
          nb::arg("n_transcripts"),
          "Find connected components of the fragment→transcript overlap graph.\n\n"
          "Uses union-find with path compression and union by rank.\n"
          "All candidates are direct transcript indices.\n"
          "Returns (n_comp, comp_t_offsets, comp_t_flat, comp_u_offsets,\n"
          "comp_u_flat) where the CSR pairs (offsets, flat) give sorted\n"
          "transcript indices and unit indices for each component.");

    // ----------------------------------------------------------------
    // Export constants so Python imports from this single source of truth.
    // ----------------------------------------------------------------
    m.attr("EM_LOG_EPSILON")         = EM_LOG_EPSILON;
    m.attr("MAX_FRAG_LEN")          = MAX_FRAG_LEN;
    m.attr("SQUAREM_BUDGET_DIVISOR") = SQUAREM_BUDGET_DIVISOR;
    m.attr("EM_PRIOR_EPSILON")       = EM_PRIOR_EPSILON;
    m.attr("ESTEP_TASK_WORK_TARGET") = ESTEP_TASK_WORK_TARGET;

    // ----------------------------------------------------------------
    // Test-only: expose fast_exp for accuracy validation from Python
    // ----------------------------------------------------------------
    m.def("_fast_exp_test_array",
          [](nb::ndarray<double, nb::ndim<1>, nb::c_contig> inputs) {
              const size_t n = inputs.shape(0);
              const double* src = inputs.data();
              // Allocate output array
              double* out = new double[n];
              size_t i = 0;
#if RIGEL_HAS_AVX512F
              for (; i + 8 <= n; i += 8) {
                  __m512d v = _mm512_loadu_pd(src + i);
                  v = rigel::fast_exp_avx512(v);
                  _mm512_storeu_pd(out + i, v);
              }
#elif RIGEL_HAS_AVX2 && RIGEL_HAS_FMA
              for (; i + 4 <= n; i += 4) {
                  __m256d v = _mm256_loadu_pd(src + i);
                  v = rigel::fast_exp_avx2(v);
                  _mm256_storeu_pd(out + i, v);
              }
#elif RIGEL_HAS_NEON
              for (; i + 2 <= n; i += 2) {
                  float64x2_t v = vld1q_f64(src + i);
                  v = rigel::fast_exp_neon(v);
                  vst1q_f64(out + i, v);
              }
#endif
              for (; i < n; ++i) {
                  out[i] = rigel::fast_exp_scalar(src[i]);
              }
              // Return as numpy array, transfer ownership
              nb::capsule owner(out, [](void* p) noexcept {
                  delete[] static_cast<double*>(p);
              });
              return nb::ndarray<nb::numpy, double, nb::ndim<1>>(
                  out, {n}, owner);
          },
          nb::arg("inputs"),
          "Apply fast_exp to an array (test-only). Returns numpy array.");
}
