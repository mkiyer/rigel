/**
 * em_solver.cpp — the per-locus EM (module rigel._em_impl).
 *
 * Builds each locus's sub-problem and equivalence classes, runs MAP-EM or VBEM under SQUAREM with the grouped
 * gDNA/RNA calibration prior, assigns the posteriors, and finds the loci as connected components.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <cstdint>
#include <cstring>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <atomic>
#include <thread>

#include "fast_exp.h"
#include "thread_pool.h"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/tuple.h>

namespace nb = nanobind;

// ================================================================
// Constants
// ================================================================

static constexpr double EM_LOG_EPSILON = 1e-300;
// One SQUAREM iteration runs three VBEM steps (two plain, one stabilising), so `iterations / 3` SQUAREM
// iterations spend the same VBEM budget `EMConfig.iterations` names.
static constexpr int    SQUAREM_BUDGET_DIVISOR = 3;

// VBEM floor: the smallest alpha a component keeps. digamma(1e-300) ≈ −1e300, so after max-subtraction the E-step
// gives a component at the floor exactly zero responsibility and the evidence-proportional prior gives it nothing: it
// stays dead unless it holds deterministic fragments. Only the EM's own step may put a component there, never
// SQUAREM's extrapolation (`backtracked_squarem_step`).
static constexpr double VBEM_SQUAREM_PRIOR_FLOOR = EM_LOG_EPSILON;

// Mode codes are labels, not quantities. The module exports them for rigel.estimator, and
// `batch_locus_em_partitioned` refuses any other value.
//
// The assignment (`assign_posteriors`)
static constexpr int ASSIGN_FRACTIONAL = 0;  // every candidate keeps its posterior share
static constexpr int ASSIGN_SAMPLE     = 1;  // count first, then one component per unit
// The warm start (`compute_grouped_warm_start`)
static constexpr int WARM_START_COVERAGE = 0;  // the default coverage-weighted seed
static constexpr int WARM_START_PRIOR    = 1;  // the prior alone: the seed is zeroed
static constexpr int WARM_START_UNIFORM  = 2;  // every component equal

// The winner code of a unit the intergenic term wins; -2 is the locus's own gDNA.
static constexpr int32_t EM_WINNER_INTERGENIC = -3;

// ================================================================
// The intergenic term: a fixed rate no locus parameter carries
// ================================================================
// A multimapper's place with no transcript candidate is in no locus, so its gDNA comes at the library's intergenic
// rate `λ` (fragments per start, in the unit the capture efficiencies are read in), beside the locus's components,
// whose absolute rates are `M θ_c / L̃_c`. `M` is the locus's own data count, re-estimated with θ: the maximiser of
// `Σ_u log(Σ_c M θ_c q_uc / L̃_c + λ b_u) + Σ_d log(M θ_d) − M` under the grouped prior on θ gives
// `M = Σ_u Σ_{c in locus} r_uc + N_det`. In the locus's frame the term is one more component, last, of weight
// `log λ − log M`; it takes no share of the prior and none of θ's normalisation, and its slot of the EM's state holds
// `M` so SQUAREM accelerates the two together: `M / N` under MAP, where the state is a share, and `M` under VBEM,
// where it is a count. A locus with no intergenic place has no such component.
struct IntergenicTerm {
    int    idx = -1;        // its component index, the last one; -1: the locus has none
    double log_rate = -std::numeric_limits<double>::infinity();
    double n_all = 0.0;     // the locus's fragments: its units with a candidate and its deterministic ones
};

// The intergenic column's weight, from the locus's data count `M`.
static inline double intergenic_log_weight(const IntergenicTerm& ig, double data_count) {
    return ig.log_rate - std::log(std::max(data_count, EM_LOG_EPSILON));
}

// ================================================================
// Per-locus profiling statistics (batch_locus_em_partitioned with emit_locus_stats)
// ================================================================

using hrclock = std::chrono::steady_clock;

struct LocusProfile {
    int locus_idx = -1;
    int n_transcripts = 0;
    int n_units = 0;
    int n_components = 0;        // n_t + 1, or n_t + 2 with the intergenic term
    int n_equiv_classes = 0;
    int64_t ec_total_elements = 0; // sum of n*k across all ECs
    int max_ec_width = 0;        // max k across ECs
    int max_ec_depth = 0;        // max n across ECs
    int squarem_iterations = 0;
    int estep_threads_used = 0;
    bool is_mega_locus = false;
    double squarem_step_scale_mean = 0.0;
    double squarem_step_scale_max = 0.0;
    int squarem_backtrack_count = 0;
    int squarem_nonfinite_count = 0;
    int squarem_grouped_stabilization_fail_count = 0;
    int assignment_stranded = 0;    // SAMPLE: units the draw left with every candidate full (assign_posteriors)
    int assignment_unrepaired = 0;  // SAMPLE: of those, the ones no repair path could place
    double gdna_eff_len = 1.0;

    // Final marginal data log-likelihood at the converged theta (diagnostic; only when emit_locus_stats).
    double final_data_loglik = 0.0;

    // Sub-phase wall times in microseconds
    double extract_us = 0.0;
    double build_ec_us = 0.0;
    double warm_start_us = 0.0;
    double squarem_us = 0.0;
    double assign_us = 0.0;
    double total_us = 0.0;
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

// Array type aliases

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
// Equivalence classes: the units sharing one candidate component set, as dense matrices
// ================================================================

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

    // Sort the classes by comp_idx. `unordered_map` iterates in no fixed order, and the E-step's column sums
    // depend on addition order in floating point, so an unstable order gives ULP differences that SQUAREM
    // amplifies. Rows within a class need no sort: `build_multi_loci` orders each locus's units by `frag_id`
    // (BAM order), so `unit_list` is already canonical.
    std::sort(result.begin(), result.end(),
              [](const EmEquivClass& a, const EmEquivClass& b) {
                  return a.comp_idx < b.comp_idx;
              });

    return result;
}

// ================================================================
// Kahan summation: O(eps) error rather than O(n·eps), for the E-step's column sums over large classes
// ================================================================

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
// EM step kernel — the hot inner loop
// ================================================================
// For every row of one equivalence class: the posteriors from log_weights, normalised by log-sum-exp, with the k
// column sums Kahan-summed in row order into col_sums. Thread-safe: the class is read-only, each row's posteriors
// live in a stack-local buffer, and col_sums is the caller's slot.

static inline void em_step_kernel(
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
    double* row;
    if (k <= MAX_K_STACK) {
        row = stack_row;
    } else {
        heap_row.resize(k);
        row = heap_row.data();
    }

    // One Kahan accumulator per column, fed in row order.
    std::vector<KahanAccumulator> col_acc(k);

    for (int i = 0; i < ec.n; ++i) {

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

                // Whole-vector early-zero skip: all 8 lanes below the cutoff
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

                // Whole-vector early-zero skip: all 4 lanes below the cutoff
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

                // Whole-vector early-zero skip: both lanes below the cutoff
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
            for (int c = 0; c < k; ++c) {
                row[c] *= inv_sum;
                col_acc[c].add(row[c]);
            }
        } else {
            // a row no component can emit (every weight -inf) lands here: outside the model, skipped
        }
    }

    for (int j = 0; j < k; ++j) {
        col_sums[j] = col_acc[j].sum;
    }
}

// ================================================================
// The E-step: per-class column sums, added in class order
// ================================================================
// Each equivalence class writes its column sums to its own slot, and the slots are added into em_totals in class
// order, so the totals are bit-identical for any thread count. A class is never split across threads, since a split
// would re-associate its row sums. Built once per locus and reused by every step.

class EStep {
public:
    // pool: the locus's threads, or nullptr to run on the caller's thread alone.
    EStep(const std::vector<EmEquivClass>& ec_data, rigel::ThreadPool* pool)
        : ec_data_(ec_data), pool_(pool), slot_(ec_data.size() + 1, 0)
    {
        const int n_ec = static_cast<int>(ec_data.size());
        for (int i = 0; i < n_ec; ++i) {
            slot_[i + 1] = slot_[i] + static_cast<size_t>(ec_data[i].k);
        }
        sums_.resize(slot_[n_ec]);
        if (pool_) {
            // The classes dealt heaviest first, each to the thread with the least work so far.
            auto work = [&](int i) { return static_cast<int64_t>(ec_data[i].n) * ec_data[i].k; };
            std::vector<int> order(static_cast<size_t>(n_ec));
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                return work(a) != work(b) ? work(a) > work(b) : a < b;
            });
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
    }

    // Adds the posteriors' column sums into em_totals [n_components], zeroed by the caller.
    void run(const double* log_weights, double* em_totals) {
        const int n_ec = static_cast<int>(ec_data_.size());
        auto run_class = [&](int i) {
            em_step_kernel(ec_data_[i], log_weights, sums_.data() + slot_[i]);
        };
        auto add_class = [&](int i) {
            const int32_t* cidx = ec_data_[i].comp_idx.data();
            const double* s = sums_.data() + slot_[i];
            for (int j = 0; j < ec_data_[i].k; ++j) {
                em_totals[cidx[j]] += s[j];
            }
        };
        if (pool_) {
            pool_->run_parallel([&](int tid) {
                for (int c = thread_start_[tid]; c < thread_start_[tid + 1]; ++c) {
                    run_class(thread_classes_[c]);
                }
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
    const std::vector<EmEquivClass>& ec_data_;
    rigel::ThreadPool* pool_;
    std::vector<size_t> slot_;   // [n_ec + 1]: class i's k sums are sums_[slot_[i], slot_[i + 1])
    std::vector<double> sums_;
    std::vector<int> thread_start_;    // threaded only: thread t runs thread_classes_[thread_start_[t], thread_start_[t + 1])
    std::vector<int> thread_classes_;
};

// ================================================================
// Grouped RNA/gDNA prior update
// ================================================================

// The calibration prior for one locus: pseudo-fragments of gDNA and of RNA, and how the RNA share is allocated.
struct AggregatePrior {
    double gdna_prior_fragments = 0.0;
    double rna_prior_fragments  = 0.0;
    int    gdna_index = -1;
    bool   has_gdna_candidate = false;

    // Per-component allocation weight for the RNA prior; nullptr allocates in proportion to the evidence each
    // component already carries. Both are the additive pseudocount `a_i = rna_prior · w_i / Σw`; the evidence rule
    //     out[i] = raw[i] · (1 + rna_prior/rna_count) = raw[i] + rna_prior · raw[i]/rna_count
    // is the case `w_i = raw[i]`. They differ at `raw[i] == 0`: evidence weights make `out[i] = 0` an absorbing
    // state, so a component with no evidence is never revived, and a strictly positive weight has no such state.
    const double* component_rna_prior_weight = nullptr;
};

static inline double nonnegative_finite(double x) {
    return (std::isfinite(x) && x > 0.0) ? x : 0.0;
}

// Split the locus into gDNA and RNA, add the calibration prior to each, and share the RNA prior over the RNA
// components in proportion to the evidence each already carries (or to the caller's allocation weight).
//
// Every RNA component, annotated transcript or synthetic span, receives its share: under evidence weights the prior
// is the common factor `(1 + rna_prior/rna_count)` and moves only the gDNA:RNA split, while a withheld share would
// move RNA between components the data cannot tell apart. Within one call the RNA components sum to
// `rna_count + rna_prior` and gDNA to `gdna_count + gdna_prior`, exactly.
static void apply_grouped_prior_update(
    const double* raw_counts,
    const double* carried_state,
    const AggregatePrior& aggregate_prior,
    double* out_counts,
    int n_components)
{
    const int gdna_index = aggregate_prior.gdna_index;
    const bool has_gdna = aggregate_prior.has_gdna_candidate
        && gdna_index >= 0 && gdna_index < n_components;

    // The RNA pool's evidence and carried totals; every RNA component is a recipient.
    double rna_count = 0.0, rna_carried = 0.0;
    for (int i = 0; i < n_components; ++i) {
        if (i == gdna_index) continue;
        rna_count += nonnegative_finite(raw_counts[i]);
        if (carried_state != nullptr) {
            rna_carried += nonnegative_finite(carried_state[i]);
        }
    }

    // The RNA components' total allocation weight, when the caller supplies one.
    const double* rna_prior_weight = aggregate_prior.component_rna_prior_weight;
    double rna_weight = 0.0;
    if (rna_prior_weight != nullptr) {
        for (int i = 0; i < n_components; ++i) {
            if (i == gdna_index) continue;
            rna_weight += nonnegative_finite(rna_prior_weight[i]);
        }
    }
    const bool weighted = rna_weight > EM_LOG_EPSILON;

    // The gDNA pseudocount needs a gDNA component to land on.
    double gdna_prior = has_gdna ? nonnegative_finite(aggregate_prior.gdna_prior_fragments) : 0.0;
    // The RNA pseudocount lands on the RNA components, so it is not gated on a gDNA component (an all-spliced locus
    // has none). It is dropped only when nothing can receive it, tested on the denominator the chosen branch divides
    // by (the weight total, else the count, else the carried state): any other total could keep a live prior at a
    // locus with zero count but carried alpha (VBEM), multiply it by zero and move the gDNA:RNA split.
    double rna_prior  = nonnegative_finite(aggregate_prior.rna_prior_fragments);
    const double prior_recipients = weighted
        ? rna_weight
        : ((rna_count > EM_LOG_EPSILON) ? rna_count : rna_carried);
    if (prior_recipients <= EM_LOG_EPSILON) {
        rna_prior = 0.0;
    }

    const double gdna_count = has_gdna ? nonnegative_finite(raw_counts[gdna_index]) : 0.0;
    const double gdna_total = gdna_count + gdna_prior;
    const double rna_total  = rna_count + rna_prior;

    std::fill(out_counts, out_counts + n_components, 0.0);
    if (has_gdna) {
        out_counts[gdna_index] = gdna_total;
    }

    // Weighted allocation: out[i] = raw[i] + rna_prior · w_i / Σw, so the pool sums to rna_count + rna_prior. A zero
    // weight total is not `weighted` and takes the evidence rule below, which keeps the prior and the split.
    if (weighted) {
        const double scale = rna_prior / rna_weight;
        for (int i = 0; i < n_components; ++i) {
            if (i == gdna_index) continue;
            out_counts[i] = nonnegative_finite(raw_counts[i])
                + nonnegative_finite(rna_prior_weight[i]) * scale;
        }
        return;
    }

    // Evidence weights: each RNA component's evidence scaled by rna_total / rna_count.
    if (rna_count > EM_LOG_EPSILON) {
        const double inv = 1.0 / rna_count;
        for (int i = 0; i < n_components; ++i) {
            if (i == gdna_index) continue;
            out_counts[i] = rna_total * nonnegative_finite(raw_counts[i]) * inv;
        }
    } else if (rna_carried > EM_LOG_EPSILON && carried_state != nullptr) {
        // Zero evidence, reachable under VBEM (which passes `alpha` as the carried state): the prior-only pool is
        // shared in proportion to the carried alpha. The warm start passes nullptr and never reaches this.
        const double inv = 1.0 / rna_carried;
        for (int i = 0; i < n_components; ++i) {
            if (i == gdna_index) continue;
            out_counts[i] = rna_total * nonnegative_finite(carried_state[i]) * inv;
        }
    }
}

// ================================================================
// MAP-EM step: theta → theta_new, under the grouped prior (`apply_grouped_prior_update`)
// ================================================================

static void map_em_step(
    const double* theta,
    EStep&        estep,
    const double* log_eff_len,
    const double* unambig_totals,
    const AggregatePrior& aggregate_prior,
    double*       em_totals,    // zeroed then accumulated
    double*       theta_new,    // output: normalized
    int           n_components,
    const IntergenicTerm& ig = IntergenicTerm{})
{
    // The locus's own components; the intergenic term, if any, is the last and its slot holds M / N.
    const int n_locus = ig.idx >= 0 ? ig.idx : n_components;

    // Compute log_weights = log(theta + epsilon) - log_eff_len
    std::vector<double> log_weights(static_cast<size_t>(n_components));
    for (int i = 0; i < n_locus; ++i) {
        log_weights[i] = std::log(theta[i] + EM_LOG_EPSILON) - log_eff_len[i];
    }
    if (ig.idx >= 0) log_weights[ig.idx] = intergenic_log_weight(ig, theta[ig.idx] * ig.n_all);

    std::fill(em_totals, em_totals + n_components, 0.0);

    estep.run(log_weights.data(), em_totals);

    // Grouped MAP-EM: raw counts first, then aggregate prior projection.
    std::vector<double> raw_counts(static_cast<size_t>(n_components));
    double total = 0.0, data = 0.0;
    for (int i = 0; i < n_components; ++i) {
        raw_counts[i] = unambig_totals[i] + em_totals[i];
    }
    for (int i = 0; i < n_locus; ++i) data += raw_counts[i];
    apply_grouped_prior_update(
        raw_counts.data(), theta, aggregate_prior, theta_new, n_locus);
    for (int i = 0; i < n_locus; ++i) {
        total += theta_new[i];
    }

    if (total > 0.0) {
        double inv_total = 1.0 / total;
        for (int i = 0; i < n_locus; ++i) {
            theta_new[i] *= inv_total;
        }
    }
    if (ig.idx >= 0) theta_new[ig.idx] = data / ig.n_all;
}

// ================================================================
// VBEM step: alpha → alpha_new
// ================================================================

static void vbem_step(
    const double* alpha,
    EStep&        estep,
    const double* log_eff_len,
    const double* unambig_totals,
    const AggregatePrior& aggregate_prior,
    double*       em_totals,
    double*       alpha_new,
    int           n_components,
    const IntergenicTerm& ig = IntergenicTerm{})
{
    // The locus's own components; the intergenic term, if any, is the last and its slot holds M.
    const int n_locus = ig.idx >= 0 ? ig.idx : n_components;

    double alpha_sum = 0.0;
    for (int i = 0; i < n_locus; ++i) {
        alpha_sum += alpha[i];
    }

    // Compute log_weights = digamma(max(alpha, eps)) - digamma(max(alpha_sum, eps)) - log_eff_len
    double dg_sum = digamma(std::max(alpha_sum, EM_LOG_EPSILON));
    std::vector<double> log_weights(static_cast<size_t>(n_components));
    for (int i = 0; i < n_locus; ++i) {
        log_weights[i] = digamma(std::max(alpha[i], EM_LOG_EPSILON))
                        - dg_sum - log_eff_len[i];
    }
    if (ig.idx >= 0) log_weights[ig.idx] = intergenic_log_weight(ig, alpha[ig.idx]);

    std::fill(em_totals, em_totals + n_components, 0.0);

    estep.run(log_weights.data(), em_totals);

    // M-step: grouped count-space update (unnormalized alpha-like mass).
    std::vector<double> raw_counts(static_cast<size_t>(n_components));
    double data = 0.0;
    for (int i = 0; i < n_components; ++i) {
        raw_counts[i] = unambig_totals[i] + em_totals[i];
    }
    for (int i = 0; i < n_locus; ++i) data += raw_counts[i];
    apply_grouped_prior_update(
        raw_counts.data(), alpha, aggregate_prior, alpha_new, n_locus);
    if (ig.idx >= 0) alpha_new[ig.idx] = data;
}

// ================================================================
// Marginal data log-likelihood at a given theta (diagnostic)
// ================================================================
// Σ_units logsumexp_c(log_weights[c] + log_lik_{u,c}): the EM objective's data term, which ranks two converged
// solutions by likelihood.
static double marginal_data_loglik(
    const std::vector<EmEquivClass>& ec_data, const double* log_weights)
{
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
            for (int j = 0; j < k; ++j) {
                s += std::exp((ll[i * k + j] + log_weights[cidx[j]]) - mx);
            }
            total += mx + std::log(s);
        }
    }
    return total;
}

// ================================================================
// Warm start, projected by the grouped prior
// ================================================================

static void compute_grouped_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const AggregatePrior& aggregate_prior,
    double*       warm_counts_out,
    int           n_components_all,
    int           warm_mode,   // a WARM_START_* code
    int           ig_idx = -1) // the intergenic term's slot, which SQUAREM seeds itself
{
    // The seed and the prior are the locus's own components'; the intergenic slot starts at 0 here.
    const int n_components = ig_idx >= 0 ? ig_idx : n_components_all;
    std::fill(warm_counts_out, warm_counts_out + n_components_all, 0.0);
    // WARM_START_COVERAGE seeds the deterministic counts plus each ambiguous unit's coverage-weighted shares, then
    // applies the grouped prior. WARM_START_PRIOR zeroes the seed, so theta starts from the prior alone; it needs an
    // allocation weight, since under evidence weights an all-zero seed gives an all-zero RNA pool. WARM_START_UNIFORM
    // starts every component equal, so where the EM lands is a property of the likelihood, not of the seed.
    std::vector<double> warm_raw(static_cast<size_t>(n_components_all), 0.0);
    if (warm_mode == WARM_START_PRIOR) {
        apply_grouped_prior_update(
            warm_raw.data(), nullptr, aggregate_prior, warm_counts_out, n_components);
        return;
    }
    if (warm_mode == WARM_START_UNIFORM) {
        std::fill(warm_raw.begin(), warm_raw.end(), 1.0);
        apply_grouped_prior_update(
            warm_raw.data(), nullptr, aggregate_prior, warm_counts_out, n_components);
        return;
    }
    std::copy(unambig_totals, unambig_totals + n_components, warm_raw.begin());

    for (const auto& ec : ec_data) {
        const int n = ec.n;
        const int k = ec.k;
        const int32_t* cidx = ec.comp_idx.data();
        const double* wt = ec.wt_flat.data();

        for (int i = 0; i < n; ++i) {
            // Normalize coverage weights within the actual candidate set.
            double row_sum = 0.0;
            for (int j = 0; j < k; ++j) {
                if (cidx[j] == ig_idx) continue;  // the intergenic term takes no warm-start share
                double w = wt[i * k + j];
                row_sum += w;
            }
            if (row_sum == 0.0) row_sum = 1.0;
            double inv_row_sum = 1.0 / row_sum;

            for (int j = 0; j < k; ++j) {
                if (cidx[j] == ig_idx) continue;
                double share = wt[i * k + j] * inv_row_sum;
                warm_raw[cidx[j]] += share;
            }
        }
    }
    apply_grouped_prior_update(
        warm_raw.data(), nullptr, aggregate_prior, warm_counts_out, n_components);
}

// ================================================================
// SQUAREM step length: backtracked, never clamped
// ================================================================
// SQUAREM jumps to `state0 + 2·step·r + step²·v`, which can carry a shrinking component below the floor, where it
// never returns. Clamping there would let the accelerator, not the EM, decide which components live. So the step
// halves its excess over 1 (the plain double step `state2`) until no component `state2` keeps above the floor lands
// below it; a component `state2` itself takes below the floor does not hold the step back. `backtracks` counts the
// halvings.
static double backtracked_squarem_step(
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

// ================================================================
// SQUAREM acceleration wrapper
// ================================================================

struct EMResult {
    std::vector<double> theta;
    double data_count = 0.0;  // the locus's own data count M at convergence, when it has an intergenic term
    int squarem_iterations = 0;  // number of SQUAREM iterations completed
    double squarem_step_scale_mean = 0.0;
    double squarem_step_scale_max = 0.0;
    int squarem_backtrack_count = 0;       // halvings of an extrapolation step (backtracked_squarem_step)
    int squarem_nonfinite_count = 0;
    int squarem_grouped_stabilization_fail_count = 0;
};

static EMResult run_squarem(
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,
    const double* unambig_totals,
    const AggregatePrior& aggregate_prior,
    const double* init_counts,
    int           n_components,
    int           max_iterations,
    double        convergence_delta,
    bool          use_vbem,
    rigel::ThreadPool* pool = nullptr,
    const IntergenicTerm& ig = IntergenicTerm{})
{
    int max_sq_iters = std::max(max_iterations / SQUAREM_BUDGET_DIVISOR, 1);
    size_t nc = static_cast<size_t>(n_components);
    // The locus's own components; the intergenic slot, if any, is the last and holds M (VBEM) or M / N (MAP).
    const size_t n_locus = ig.idx >= 0 ? static_cast<size_t>(ig.idx) : nc;
    EStep estep(ec_data, pool);

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
    double data_count_out = 0.0;
    int completed_iterations = 0;
    double step_scale_sum = 0.0;
    double step_scale_max = 0.0;
    int step_scale_count = 0;
    int backtrack_count = 0;
    int nonfinite_count = 0;
    int stabilization_fail_count = 0;

    if (use_vbem) {
        // ---- VBEM under SQUAREM; the state is the Dirichlet alpha ----
        for (size_t i = 0; i < nc; ++i) {
            state0[i] = std::max(init_counts[i], VBEM_SQUAREM_PRIOR_FLOOR);
        }
        if (ig.idx >= 0) state0[ig.idx] = ig.n_all;  // every fragment the locus's, to start

        for (int iter = 0; iter < max_sq_iters; ++iter) {
            // Two plain VBEM steps
            vbem_step(state0.data(), estep, log_eff_len,
                      unambig_totals, aggregate_prior, em_totals.data(),
                      state1.data(), n_components, ig);

            vbem_step(state1.data(), estep, log_eff_len,
                      unambig_totals, aggregate_prior, em_totals.data(),
                      state2.data(), n_components, ig);

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
                if (!std::isfinite(step)) {
                    step = 1.0;
                    ++nonfinite_count;
                }
                step = backtracked_squarem_step(step, state0, r_vec, v_vec, state2,
                                                VBEM_SQUAREM_PRIOR_FLOOR, backtrack_count);
                step_scale_sum += step;
                step_scale_max = std::max(step_scale_max, step);
                ++step_scale_count;
                for (size_t i = 0; i < nc; ++i) {
                    // At step 1 the jump is the plain double step: copied, because recomputing it rounds
                    // on the scale of state0 and could carry a small live component below the floor.
                    state_extrap[i] = (step == 1.0)
                        ? state2[i]
                        : state0[i] + 2.0 * step * r_vec[i] + step * step * v_vec[i];
                    // Only a component the plain step itself carried to the floor can land below it here.
                    double floor_i = VBEM_SQUAREM_PRIOR_FLOOR;
                    if (!std::isfinite(state_extrap[i])) {
                        state_extrap[i] = floor_i;
                        ++nonfinite_count;
                    } else if (state_extrap[i] < floor_i) {
                        state_extrap[i] = floor_i;
                    }
                }
            }

            // Stabilisation step
            vbem_step(state_extrap.data(), estep, log_eff_len,
                      unambig_totals, aggregate_prior, em_totals.data(),
                      state_new.data(), n_components, ig);

            // Clamp to the floor, then test convergence on normalised theta.
            double sum_old = 0.0, sum_new = 0.0;
            for (size_t i = 0; i < nc; ++i) {
                double floor_i = VBEM_SQUAREM_PRIOR_FLOOR;
                if (!std::isfinite(state_new[i])) {
                    state_new[i] = floor_i;
                    ++nonfinite_count;
                    ++stabilization_fail_count;
                } else if (state_new[i] < floor_i) {
                    state_new[i] = floor_i;
                }
                if (i < n_locus) {
                    sum_old += state0[i];
                    sum_new += state_new[i];
                }
            }
            double delta = 0.0;
            if (sum_old > 0.0 && sum_new > 0.0) {
                double inv_old = 1.0 / sum_old;
                double inv_new = 1.0 / sum_new;
                for (size_t i = 0; i < n_locus; ++i) {
                    delta += std::abs(state_new[i] * inv_new
                                    - state0[i] * inv_old);
                }
            }
            if (ig.idx >= 0 && ig.n_all > 0.0) delta += std::abs(state_new[ig.idx] - state0[ig.idx]) / ig.n_all;

            std::swap(state0, state_new);

            if (delta < convergence_delta) {
                completed_iterations = iter + 1;
                break;
            }
            completed_iterations = iter + 1;
        }

        std::copy(state0.begin(), state0.end(), alpha_out.begin());

        // theta = normalized alpha, over the locus's own components
        double total = 0.0;
        for (size_t i = 0; i < n_locus; ++i) total += alpha_out[i];
        if (total > 0.0) {
            double inv = 1.0 / total;
            for (size_t i = 0; i < n_locus; ++i) theta[i] = alpha_out[i] * inv;
        } else {
            std::copy(alpha_out.begin(), alpha_out.begin() + n_locus, theta.begin());
        }
        if (ig.idx >= 0) {
            theta[ig.idx] = 0.0;
            data_count_out = alpha_out[ig.idx];
        }

    } else {
        // ---- MAP-EM under SQUAREM; the state is the normalised theta ----
        double total = 0.0;
        for (size_t i = 0; i < n_locus; ++i) {
            state0[i] = nonnegative_finite(init_counts[i]);
            total += state0[i];
        }
        if (total > 0.0) {
            double inv = 1.0 / total;
            for (size_t i = 0; i < n_locus; ++i) state0[i] *= inv;
        }
        if (ig.idx >= 0) state0[ig.idx] = 1.0;  // M = N: every fragment the locus's, to start

        for (int iter = 0; iter < max_sq_iters; ++iter) {
            // Two EM steps
            map_em_step(state0.data(), estep, log_eff_len,
                        unambig_totals, aggregate_prior, em_totals.data(),
                        state1.data(), n_components, ig);

            map_em_step(state1.data(), estep, log_eff_len,
                        unambig_totals, aggregate_prior, em_totals.data(),
                        state2.data(), n_components, ig);

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
                if (!std::isfinite(alpha_step)) {
                    alpha_step = 1.0;
                    ++nonfinite_count;
                }
                alpha_step = backtracked_squarem_step(alpha_step, state0, r_vec, v_vec, state2,
                                                      0.0, backtrack_count);
                step_scale_sum += alpha_step;
                step_scale_max = std::max(step_scale_max, alpha_step);
                ++step_scale_count;
                for (size_t i = 0; i < nc; ++i) {
                    state_extrap[i] = (alpha_step == 1.0)
                        ? state2[i]
                        : state0[i] + 2.0 * alpha_step * r_vec[i]
                          + alpha_step * alpha_step * v_vec[i];
                    if (!std::isfinite(state_extrap[i])) {
                        state_extrap[i] = 0.0;
                        ++nonfinite_count;
                    } else if (state_extrap[i] < 0.0) {
                        state_extrap[i] = 0.0;
                    }
                }
                double s = 0.0;
                for (size_t i = 0; i < n_locus; ++i) s += state_extrap[i];
                if (s > 0.0) {
                    double inv = 1.0 / s;
                    for (size_t i = 0; i < n_locus; ++i) state_extrap[i] *= inv;
                } else {
                    std::copy(state2.begin(), state2.end(),
                              state_extrap.begin());
                }
            }

            // Stabilisation step
            map_em_step(state_extrap.data(), estep, log_eff_len,
                        unambig_totals, aggregate_prior, em_totals.data(),
                        state_new.data(), n_components, ig);

            for (size_t i = 0; i < nc; ++i) {
                if (!std::isfinite(state_new[i])) {
                    state_new[i] = 0.0;
                    ++nonfinite_count;
                    ++stabilization_fail_count;
                }
            }

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

        std::copy(state0.begin(), state0.end(), theta.begin());
        if (ig.idx >= 0) {
            data_count_out = state0[ig.idx] * ig.n_all;
            theta[ig.idx] = 0.0;
        }
    }

    EMResult out{ std::move(theta), data_count_out, completed_iterations };
    out.squarem_step_scale_mean = step_scale_count > 0
        ? step_scale_sum / static_cast<double>(step_scale_count)
        : 0.0;
    out.squarem_step_scale_max = step_scale_max;
    out.squarem_backtrack_count = backtrack_count;
    out.squarem_nonfinite_count = nonfinite_count;
    out.squarem_grouped_stabilization_fail_count = stabilization_fail_count;
    return out;
}

// ================================================================
// Batch locus EM — single C++ call for all loci
// ================================================================

// Per-locus candidate record (used during sub-problem extraction)
struct LocalCandidate {
    int32_t local_comp;
    double  log_lik;
    double  cov_wt;
    uint8_t count_col;
};

// Per-locus sub-problem (stack-allocated, reused across loci)
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

    // Local→global transcript mapping
    std::vector<int32_t>  local_to_global_t; // [n_t]
};

// Assign each unit's posterior at the converged theta and scatter it into the accumulators. FRACTIONAL scatters
// every posterior as it stands. SAMPLE gives each unit to one component, count first: each component's fractional
// count in the locus (transcripts and gDNA alike) is rounded by largest remainder, then one pass draws each unit from
// its posterior re-weighted per candidate by (count still owed / posterior mass still to come). A unit whose
// candidates are all full keeps its most probable one (stranded) and a repair along augmenting paths puts every
// count back on target, or within one of its fractional count where the exact targets are unreachable. The result
// is feasible, not optimal: no candidate is dropped per fragment, and the repair ignores posteriors.
// Gate: tests/test_estimator.py (TestWholeCounts). rng: this locus's SplitMix64 stream (SAMPLE only).
static void assign_posteriors(
    const LocusSubProblem& sub,
    const double* theta,
    const double* log_eff_len,
    int assignment_mode,
    SplitMix64& rng,
    int& stranded,     // SAMPLE: units the draw left with every candidate full
    int& unrepaired,   // SAMPLE: of those, the ones no repair path could place (0 whenever the counts are reachable)
    // Output accumulators (accumulated across loci)
    double* em_counts_2d,          // [N_T, n_cols], row-major
    double* posterior_sum,         // [N_T]
    double* n_assigned,            // [N_T]
    // Per-locus totals. rna_total is the mass on every transcript component, annotated transcripts and synthetic
    // spans alike; the caller separates the two. rna_total + gdna_total + ig_total is the assigned mass.
    double& rna_total,
    double& gdna_total,
    int N_T_TOTAL,  // total transcripts for bounds checking
    int n_cols,     // number of splice-strand columns (actual 2D stride)
    // --- Per-unit assignment output (nullable) ---
    int32_t* out_winner_tid,   // [total_units] or nullptr
    float*   out_winner_post,  // [total_units] or nullptr
    int16_t* out_n_candidates, // [total_units] or nullptr
    int      out_offset,       // write offset for this locus
    // The intergenic term: its weight from the locus's converged data count, and its mass, outside the locus
    const IntergenicTerm& ig,
    double   data_count,
    double&  ig_total)
{
    int n_t = sub.n_t;
    int nc  = sub.n_components;
    int n_units = sub.n_local_units;
    const int32_t* local_to_global = sub.local_to_global_t.data();
    stranded = 0;
    unrepaired = 0;

    // log(theta[c] + eps) - log L̃_c: the MAP E-step's weights, at the converged theta whichever EM ran.
    std::vector<double> log_weights(nc);
    for (int c = 0; c < nc; ++c) {
        log_weights[c] = std::log(theta[c] + EM_LOG_EPSILON) - log_eff_len[c];
    }
    if (ig.idx >= 0) log_weights[ig.idx] = intergenic_log_weight(ig, data_count);

    // Every unit's posterior, candidate-aligned with sub.t_indices. A unit whose log-sum-exp is not finite keeps
    // all zeros and is assigned nothing, in either mode.
    const size_t n_cand = static_cast<size_t>(sub.offsets[n_units]);
    std::vector<double> post(n_cand, 0.0);
    for (int ui = 0; ui < n_units; ++ui) {
        auto s = sub.offsets[ui];
        auto e = sub.offsets[ui + 1];
        if (e == s) continue;
        double max_val = -1e300;
        for (auto k = s; k < e; ++k) {
            double lp = sub.log_liks[k] + log_weights[sub.t_indices[k]];
            if (lp > max_val) max_val = lp;
        }
        double sum_exp = 0.0;
        for (auto k = s; k < e; ++k) {
            post[k] = std::exp(sub.log_liks[k] + log_weights[sub.t_indices[k]] - max_val);
            sum_exp += post[k];
        }
        if (sum_exp > 0.0 && std::isfinite(sum_exp)) {
            double inv = 1.0 / sum_exp;
            for (auto k = s; k < e; ++k) post[k] *= inv;
        } else {
            for (auto k = s; k < e; ++k) post[k] = 0.0;
        }
    }

    // SAMPLE's decision: one candidate position per unit (-1: the unit has no posterior and gets nothing).
    std::vector<int> choice;
    if (assignment_mode == ASSIGN_SAMPLE) {
        // Count first: each component's fractional count, rounded within the locus by largest remainder.
        std::vector<double> to_come(nc, 0.0);  // each component's posterior mass among the units not yet drawn
        for (size_t k = 0; k < n_cand; ++k) to_come[sub.t_indices[k]] += post[k];
        const std::vector<double> n_frac(to_come);   // the fractional counts themselves
        std::vector<int64_t> quota(nc, 0);
        double total = 0.0;
        int64_t floors = 0;
        for (int c = 0; c < nc; ++c) {
            total += to_come[c];
            quota[c] = static_cast<int64_t>(std::floor(to_come[c]));
            floors += quota[c];
        }
        // The fractional counts sum to the number of units with a posterior: an integer, up to rounding.
        int64_t short_by = std::llround(total) - floors;
        if (short_by > 0) {
            std::vector<int> order(nc);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                double ra = to_come[a] - std::floor(to_come[a]);
                double rb = to_come[b] - std::floor(to_come[b]);
                if (ra != rb) return ra > rb;                                    // the largest remainder first
                if (to_come[a] != to_come[b]) return to_come[a] > to_come[b];    // a tie to the larger count
                return a < b;
            });
            for (int64_t i = 0; i < short_by && i < nc; ++i) quota[order[i]] += 1;
        }

        // The draw: units in order, each from its own posterior re-weighted by (count owed / mass to come).
        std::vector<int64_t> owed(quota);
        choice.assign(n_units, -1);
        for (int ui = 0; ui < n_units; ++ui) {
            auto s = sub.offsets[ui];
            auto e = sub.offsets[ui + 1];
            bool has_posterior = false;
            double drawable = 0.0;
            for (auto k = s; k < e; ++k) {
                if (post[k] > 0.0) has_posterior = true;
                int32_t c = sub.t_indices[k];
                if (owed[c] > 0 && post[k] > 0.0) {
                    drawable += post[k] * (static_cast<double>(owed[c]) / std::max(to_come[c], 1e-300));
                }
            }
            if (!has_posterior) continue;
            int pick = -1;
            if (drawable > 0.0) {
                double target = rng.uniform() * drawable, acc = 0.0;
                for (auto k = s; k < e; ++k) {
                    int32_t c = sub.t_indices[k];
                    if (owed[c] > 0 && post[k] > 0.0) {
                        acc += post[k] * (static_cast<double>(owed[c]) / std::max(to_come[c], 1e-300));
                        pick = static_cast<int>(k - s);   // the last drawable one, if rounding leaves target >= acc
                        if (target < acc) break;
                    }
                }
            } else {
                // Every candidate already holds its count: the unit keeps its most probable one, for now.
                double best = -1.0;
                for (auto k = s; k < e; ++k) {
                    if (post[k] > best) { best = post[k]; pick = static_cast<int>(k - s); }
                }
                ++stranded;
            }
            choice[ui] = pick;
            owed[sub.t_indices[s + pick]] -= 1;
            for (auto k = s; k < e; ++k) to_come[sub.t_indices[k]] -= post[k];
        }

        // The repair: a breadth-first search links each over-full component to an under-full one through units that
        // can move (each only to another of its own candidates with a posterior above zero); moving the chain puts
        // both on their counts. A path exists whenever the rounded counts are reachable at all.
        if (stranded > 0) {
            std::vector<std::vector<int>> members(nc);
            std::vector<int> slot(n_units, -1);
            for (int ui = 0; ui < n_units; ++ui) {
                if (choice[ui] < 0) continue;
                int c = sub.t_indices[sub.offsets[ui] + choice[ui]];
                slot[ui] = static_cast<int>(members[c].size());
                members[c].push_back(ui);
            }
            auto move = [&](int ui, int to_pos) {
                int from = sub.t_indices[sub.offsets[ui] + choice[ui]];
                int to = sub.t_indices[sub.offsets[ui] + to_pos];
                int last = members[from].back();
                members[from][slot[ui]] = last;
                slot[last] = slot[ui];
                members[from].pop_back();
                slot[ui] = static_cast<int>(members[to].size());
                members[to].push_back(ui);
                choice[ui] = to_pos;
            };
            std::vector<int64_t> over(nc);
            for (int c = 0; c < nc; ++c) over[c] = static_cast<int64_t>(members[c].size()) - quota[c];
            std::vector<int> seen(nc, -1), via_unit(nc, -1), via_pos(nc, -1), from_comp(nc, -1);
            std::vector<int> queue;
            int search = 0;
            for (int c0 = 0; c0 < nc; ++c0) {
                while (over[c0] > 0) {
                    ++search;
                    queue.assign(1, c0);
                    seen[c0] = search;
                    int found = -1;
                    for (size_t qi = 0; qi < queue.size() && found < 0; ++qi) {
                        int x = queue[qi];
                        for (int v : members[x]) {
                            auto s = sub.offsets[v];
                            auto e = sub.offsets[v + 1];
                            for (auto k = s; k < e; ++k) {
                                int y = sub.t_indices[k];
                                if (post[k] <= 0.0 || seen[y] == search) continue;
                                seen[y] = search;
                                from_comp[y] = x;
                                via_unit[y] = v;
                                via_pos[y] = static_cast<int>(k - s);
                                if (over[y] < 0) { found = y; break; }
                                queue.push_back(y);
                            }
                            if (found >= 0) break;
                        }
                    }
                    if (found < 0) break;          // no component under its count is reachable from c0
                    over[found] += 1;
                    over[c0] -= 1;
                    for (int y = found; y != c0; ) {
                        int x = from_comp[y];
                        move(via_unit[y], via_pos[y]);
                        y = x;
                    }
                }
            }
            unrepaired = 0;
            for (int c = 0; c < nc; ++c) if (over[c] > 0) unrepaired += static_cast<int>(over[c]);

            // Fallback when the exact targets are unreachable. Counts within one of every fractional count always are
            // reachable: the posterior is a fractional assignment inside [floor(n_c), ceil(n_c)], so by flow
            // integrality an integral one exists. The same augmenting paths find it: counts above their ceiling shed,
            // then counts below their floor draw; each path fixes one unit of violation and creates none.
            if (unrepaired > 0) {
                std::vector<int64_t> lo(nc), hi(nc);
                for (int c = 0; c < nc; ++c) {
                    lo[c] = static_cast<int64_t>(std::floor(n_frac[c]));
                    hi[c] = static_cast<int64_t>(std::ceil(n_frac[c]));
                }
                auto count = [&](int c) { return static_cast<int64_t>(members[c].size()); };
                // shed: forward from a component above its ceiling to one below its ceiling
                for (int c0 = 0; c0 < nc; ++c0) {
                    while (count(c0) > hi[c0]) {
                        ++search;
                        queue.assign(1, c0);
                        seen[c0] = search;
                        int found = -1;
                        for (size_t qi = 0; qi < queue.size() && found < 0; ++qi) {
                            int x = queue[qi];
                            for (int v : members[x]) {
                                auto s = sub.offsets[v];
                                auto e = sub.offsets[v + 1];
                                for (auto k = s; k < e; ++k) {
                                    int y = sub.t_indices[k];
                                    if (post[k] <= 0.0 || seen[y] == search) continue;
                                    seen[y] = search;
                                    from_comp[y] = x;
                                    via_unit[y] = v;
                                    via_pos[y] = static_cast<int>(k - s);
                                    if (count(y) < hi[y]) { found = y; break; }
                                    queue.push_back(y);
                                }
                                if (found >= 0) break;
                            }
                        }
                        if (found < 0) break;
                        for (int y = found; y != c0; ) {
                            int x = from_comp[y];
                            move(via_unit[y], via_pos[y]);
                            y = x;
                        }
                    }
                }
                // fill: backward from a component below its floor, through the units that could move into it
                std::vector<std::vector<std::pair<int, int>>> into(nc);   // component -> (unit, candidate position)
                for (int ui = 0; ui < n_units; ++ui) {
                    auto s = sub.offsets[ui];
                    for (auto k = s; k < sub.offsets[ui + 1]; ++k) {
                        if (post[k] > 0.0) into[sub.t_indices[k]].emplace_back(ui, static_cast<int>(k - s));
                    }
                }
                for (int c0 = 0; c0 < nc; ++c0) {
                    while (count(c0) < lo[c0]) {
                        ++search;
                        queue.assign(1, c0);
                        seen[c0] = search;
                        int found = -1;
                        for (size_t qi = 0; qi < queue.size() && found < 0; ++qi) {
                            int y = queue[qi];
                            for (const auto& [v, pos] : into[y]) {
                                if (choice[v] < 0) continue;
                                int x = sub.t_indices[sub.offsets[v] + choice[v]];
                                if (seen[x] == search) continue;
                                seen[x] = search;
                                from_comp[x] = y;        // v moves x -> y
                                via_unit[x] = v;
                                via_pos[x] = pos;
                                if (count(x) > lo[x]) { found = x; break; }
                                queue.push_back(x);
                            }
                        }
                        if (found < 0) break;
                        for (int x = found; x != c0; ) {
                            int y = from_comp[x];
                            move(via_unit[x], via_pos[x]);
                            x = y;
                        }
                    }
                }
            }
        }
    }

    rna_total = 0.0;
    gdna_total = 0.0;
    ig_total = 0.0;
    std::vector<double> weights;

    for (int ui = 0; ui < n_units; ++ui) {
        auto s = sub.offsets[ui];
        auto e = sub.offsets[ui + 1];
        int seg_len = static_cast<int>(e - s);
        if (seg_len == 0) continue;
        const double* posteriors = post.data() + s;
        weights.assign(seg_len, 0.0);
        int winner = -1;  // SAMPLE's component index within the unit (also the annotation output)
        if (assignment_mode == ASSIGN_FRACTIONAL) {
            for (int j = 0; j < seg_len; ++j) weights[j] = posteriors[j];
        } else if (choice[ui] >= 0) {
            winner = choice[ui];
            weights[winner] = 1.0;
        }

        // --- Per-unit annotation output ---
        if (out_winner_tid != nullptr) {
            int ann_winner = winner;
            if (assignment_mode == ASSIGN_FRACTIONAL) {
                // Fractional mode: the most probable component, for annotation display
                double best_post = -1.0;
                for (int j = 0; j < seg_len; ++j) {
                    if (posteriors[j] > best_post) {
                        best_post = posteriors[j];
                        ann_winner = j;
                    }
                }
            }

            int write_idx = out_offset + ui;
            out_n_candidates[write_idx] = static_cast<int16_t>(
                std::min(seg_len, static_cast<int>(INT16_MAX)));

            if (ann_winner >= 0) {
                int32_t comp = sub.t_indices[s + ann_winner];
                if (comp < n_t) {
                    out_winner_tid[write_idx] = local_to_global[comp];
                } else if (comp == ig.idx) {
                    out_winner_tid[write_idx] = EM_WINNER_INTERGENIC;
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
                // A transcript component, annotated transcript or synthetic span: RNA, into rna_total
                // (the Python caller separates the two).
                int32_t global_t = local_to_global[comp];
                uint8_t col = sub.count_cols[s + j];
                if (global_t < 0 || global_t >= N_T_TOTAL || col >= n_cols) continue;
                em_counts_2d[global_t * n_cols + col] += p;
                rna_total += p;

                // Confidence tracking (use original posteriors)
                posterior_sum[global_t] += posteriors[j] * posteriors[j];
                n_assigned[global_t] += posteriors[j];
            } else if (comp == ig.idx) {
                // the intergenic term: gDNA outside every locus
                ig_total += p;
            } else {
                // gDNA
                gdna_total += p;
            }
        }
    }
}

// ================================================================
// Partition scatter functions — array-by-array global→per-locus scatter
// ================================================================

/// Per-locus CSR offsets from the global offsets and each locus's unit list: offsets[li][0] = 0 and
/// offsets[li][k+1] = offsets[li][k] + (g_offsets[u+1] - g_offsets[u]) for its k-th unit u. One int64 array per locus.
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

    // Compute offsets (sequential)
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

/// Scatter a per-candidate global CSR array into per-locus arrays: the segment of the locus's k-th unit u,
/// src[g_off[u], g_off[u+1]), is copied to dst + p_off[k]. One array per locus.
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

/// Gather a per-unit global array into per-locus arrays: dst[k] = src[units[k]]. One array per locus.
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
    const float*   log_liks;
    const float*   coverage_wts;
    const uint8_t* count_cols;
    const uint8_t* is_spliced;
    const float*   gdna_log_liks;
    const float*   intergenic_log_liks = nullptr;  // per unit; -inf where a unit has no intergenic place
    int     n_units;
    int64_t n_candidates;

    // Locus transcript membership (from Locus.transcript_indices)
    const int32_t* transcript_indices;
    int n_transcripts;

};

// Extract one locus's sub-problem from its PartitionView: remap global transcript indices to local components
// [0, n_t), append the gDNA component (n_t) to unspliced units and the intergenic term (n_t + 1) to units with an
// intergenic place, and sort each unit's candidates by component, as the equivalence-class builder requires.
// `all_t_eff_lens` is the EM's capture-contracted effective length per transcript (computed in Python) and
// `gdna_eff_len` the locus's gDNA counterpart. The inner loop writes through a cursor and allocates nothing.
static void extract_locus_sub_problem_from_partition(
    LocusSubProblem& sub,
    const PartitionView& pv,
    double gdna_eff_len,
    const double*  all_unambig_row_sums,
    const double*  all_t_eff_lens,
    const double*  all_t_rna_prior_weight,
    int32_t* local_map, int local_map_size)
{
    int n_t = pv.n_transcripts;
    int n_u = pv.n_units;
    const int32_t* t_arr = pv.transcript_indices;

    sub.n_t = n_t;
    sub.n_local_units = n_u;
    sub.gdna_idx = n_t;
    sub.has_gdna_candidate = false;
    // The intergenic term is a component only where some unit has an intergenic place.
    bool any_intergenic = false;
    for (int ui = 0; ui < n_u && !any_intergenic; ++ui)
        any_intergenic = std::isfinite(static_cast<double>(pv.intergenic_log_liks[ui]));
    sub.ig_idx = any_intergenic ? n_t + 1 : -1;
    sub.n_components = n_t + (any_intergenic ? 2 : 1);
    int nc = sub.n_components;

    // --- Build global→local mapping --- (every entry arrives -1: the scratch is allocated so and reset below)
    for (int i = 0; i < n_t; ++i) {
        int gt = t_arr[i];
        if (gt >= 0 && gt < local_map_size) local_map[gt] = i;
    }

    // Pre-size to the worst case: the input candidates plus one gDNA and one intergenic candidate per unit.
    size_t max_out = static_cast<size_t>(pv.n_candidates)
                   + static_cast<size_t>(n_u) * (any_intergenic ? 2 : 1);

    sub.t_indices.resize(max_out);
    sub.log_liks.resize(max_out);
    sub.coverage_wts.resize(max_out);
    sub.count_cols.resize(max_out);

    sub.offsets.resize(n_u + 1);
    sub.offsets[0] = 0;

    // Reusable sort buffer — persists across units, only grows.
    std::vector<LocalCandidate> sort_buf;

    size_t cursor = 0;

    for (int ui = 0; ui < n_u; ++ui) {
        auto p_start = pv.offsets[ui];
        auto p_end   = pv.offsets[ui + 1];
        int width_in = static_cast<int>(p_end - p_start);

        // Determine if this unit gets a gDNA candidate
        bool is_spliced = (pv.is_spliced[ui] != 0);
        double gdna_ll = static_cast<double>(pv.gdna_log_liks[ui]);
        bool has_gdna = (!is_spliced && std::isfinite(gdna_ll));
        double ig_ll = any_intergenic ? static_cast<double>(pv.intergenic_log_liks[ui])
                                      : -std::numeric_limits<double>::infinity();
        bool has_ig = std::isfinite(ig_ll);
        int width_out = width_in + (has_gdna ? 1 : 0) + (has_ig ? 1 : 0);

        // Fill sort buffer with remapped candidates
        sort_buf.resize(width_out);
        int k = 0;

        for (auto j = p_start; j < p_end; ++j) {
            int32_t global_t = pv.t_indices[j];
            if (global_t < 0 || global_t >= local_map_size) continue;
            int32_t local = local_map[global_t];
            if (local < 0 || local >= nc) continue;

            double log_lik = static_cast<double>(pv.log_liks[j]);
            double coverage_wt = static_cast<double>(pv.coverage_wts[j]);
            sort_buf[k++] = {local, log_lik, coverage_wt, pv.count_cols[j]};
        }

        // The gDNA candidate (component n_t): the scorer supplies log h_G(ell_f) and the non-length terms; the EM
        // applies the locus's -log L̃_gDNA.
        if (has_gdna) {
            sort_buf[k++] = {sub.gdna_idx, gdna_ll, 1.0, 0};
            sub.has_gdna_candidate = true;
        }
        // The intergenic term (component n_t + 1, the largest): the unit's gDNA terms at its places outside every
        // locus, which the EM reads at the intergenic rate.
        if (has_ig) {
            sort_buf[k++] = {sub.ig_idx, ig_ll, 1.0, 0};
        }

        // Trim to actual count (candidates may have been skipped above)
        int actual = k;

        // Sort by local component index (required by the equivalence-class builder).
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
            sub.count_cols[cursor]   = c.count_col;
            ++cursor;
        }

        sub.offsets[ui + 1] = static_cast<int64_t>(cursor);
    }

    // Trim output arrays to actual size
    sub.t_indices.resize(cursor);
    sub.log_liks.resize(cursor);
    sub.coverage_wts.resize(cursor);
    sub.count_cols.resize(cursor);

    // --- Build per-component arrays ---
    sub.local_to_global_t.resize(n_t);
    for (int i = 0; i < n_t; ++i) sub.local_to_global_t[i] = t_arr[i];

    sub.unambig_totals.assign(nc, 0.0);
    for (int i = 0; i < n_t; ++i) {
        sub.unambig_totals[i] = all_unambig_row_sums[t_arr[i]];
    }

    // The RNA prior's allocation weight, left empty when the caller supplies none (every site then tests the same
    // nullptr). The gDNA component is not a recipient and keeps 0.
    sub.component_rna_prior_weight.clear();
    if (all_t_rna_prior_weight != nullptr) {
        sub.component_rna_prior_weight.assign(nc, 0.0);
        for (int i = 0; i < n_t; ++i) {
            sub.component_rna_prior_weight[i] = all_t_rna_prior_weight[t_arr[i]];
        }
    }

    // A component's yield enters unfloored. A yield of 0 is a component with no start position, which
    // cannot emit: +inf here is -inf in every E-step weight, and a row no component can emit is skipped.
    const double kInf = std::numeric_limits<double>::infinity();
    sub.log_eff_len.assign(nc, 0.0);
    for (int i = 0; i < n_t; ++i) {
        double Le = all_t_eff_lens[t_arr[i]];
        sub.log_eff_len[i] = Le > 0.0 ? std::log(Le) : kInf;
    }
    // gDNA: this MultiLocus's capture-contracted effective length (uncontracted off capture), as for the RNA ones.
    double GLe = gdna_eff_len;
    sub.log_eff_len[sub.gdna_idx] = (GLe > 0.0 ? std::log(GLe) : kInf);

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
    nb::ndarray<nb::numpy, double, nb::ndim<1>>,  // locus_rna_total[n_loci] (annotated transcripts + synthetic spans)
    nb::ndarray<nb::numpy, double, nb::ndim<1>>,  // locus_gdna[n_loci]
    nb::list,  // locus_stats
    nb::object,  // out_winner_tid (ndarray or None)
    nb::object,  // out_winner_post (ndarray or None)
    nb::object,  // out_n_candidates (ndarray or None)
    double       // total_intergenic_em: the intergenic term's mass, outside every locus
>
batch_locus_em_partitioned(
    // Per-locus partition data (list of 8-tuples)
    nb::list partition_tuples,
    // Per-locus transcript membership (list of int32[])
    nb::list locus_transcript_indices,
    // Per-locus grouped additive calibration priors
    f64_1d   locus_gdna_prior_count,
    f64_1d   locus_rna_prior_count,
    // Per-locus capture-contracted effective length of the gDNA component (uncontracted off capture).
    f64_1d   locus_gdna_eff_lens,
    // Per-transcript globals
    f64_2d   unambig_counts,
    f64_1d   t_eff_lens_arr,
    // The RNA prior's per-transcript allocation weight; empty allocates in proportion to current evidence.
    f64_1d   t_rna_prior_weight_arr,
    // A WARM_START_* code: COVERAGE (the coverage-weighted seed), PRIOR (seed zeroed: theta starts from the
    // prior alone) or UNIFORM (every component equal).
    int      warm_start_mode,
    // Mutable output accumulators
    f64_2d_mut em_counts_out,
    f64_1d_mut posterior_sum_out,
    f64_1d_mut n_assigned_out,
    // EM config
    int    max_iterations,
    double convergence_delta,
    bool   use_vbem,
    int    assignment_mode,
    uint64_t rng_seed,
    int    n_transcripts_total,
    int    n_splice_strand_cols,
    int    n_threads,
    bool   emit_locus_stats,
    bool   emit_assignments,
    // The intergenic gDNA rate, fragments per start in the efficiencies' unit; 0 reads no intergenic term.
    double intergenic_rate)
{
    // An undefined code is refused: an unknown assignment code would index the per-unit draw only ASSIGN_SAMPLE
    // fills, and an unknown warm start would silently run the coverage seed.
    if (assignment_mode != ASSIGN_FRACTIONAL && assignment_mode != ASSIGN_SAMPLE) {
        throw std::invalid_argument(
            "batch_locus_em_partitioned: assignment_mode must be ASSIGN_FRACTIONAL or ASSIGN_SAMPLE; got " +
            std::to_string(assignment_mode));
    }
    if (warm_start_mode != WARM_START_COVERAGE && warm_start_mode != WARM_START_PRIOR &&
        warm_start_mode != WARM_START_UNIFORM) {
        throw std::invalid_argument(
            "batch_locus_em_partitioned: warm_start_mode must be WARM_START_COVERAGE, WARM_START_PRIOR or "
            "WARM_START_UNIFORM; got " + std::to_string(warm_start_mode));
    }

    int n_loci = static_cast<int>(nb::len(partition_tuples));
    int N_T = n_transcripts_total;
    int N_COLS = n_splice_strand_cols;

    // --- Under GIL: extract PartitionViews from tuples ---
    std::vector<PartitionView> views(n_loci);
    for (int i = 0; i < n_loci; ++i) {
        nb::tuple tup = nb::borrow<nb::tuple>(partition_tuples[i]);
        if (nb::len(tup) != 8)
            throw std::invalid_argument("batch_locus_em_partitioned: a locus's partition is an 8-tuple");
        auto& v = views[i];
        auto off_arr = nb::cast<i64_1d>(tup[0]);
        v.offsets          = off_arr.data();
        v.t_indices        = nb::cast<i32_1d>(tup[1]).data();
        // Float payloads are cast without conversion: a float64 array raises rather than being copied.
        v.log_liks         = nb::cast<f32_1d>(tup[2], false).data();
        v.coverage_wts     = nb::cast<f32_1d>(tup[3], false).data();
        v.count_cols       = nb::cast<u8_1d>(tup[4]).data();
        v.is_spliced       = nb::cast<u8_1d>(tup[5]).data();
        v.gdna_log_liks    = nb::cast<f32_1d>(tup[6], false).data();
        v.intergenic_log_liks = nb::cast<f32_1d>(tup[7], false).data();
        v.n_units = static_cast<int>(off_arr.shape(0)) - 1;
        v.n_candidates = v.offsets[v.n_units];

        auto t_arr = nb::cast<i32_1d>(locus_transcript_indices[i]);
        v.transcript_indices = t_arr.data();
        v.n_transcripts = static_cast<int>(t_arr.shape(0));
    }

    if (static_cast<int>(locus_gdna_eff_lens.shape(0)) != n_loci) {
        throw std::runtime_error(
            "batch_locus_em_partitioned: locus_gdna_eff_lens length must equal n_loci");
    }
    if (static_cast<int>(locus_gdna_prior_count.shape(0)) != n_loci) {
        throw std::runtime_error(
            "batch_locus_em_partitioned: locus_gdna_prior_count length must equal n_loci");
    }
    if (static_cast<int>(locus_rna_prior_count.shape(0)) != n_loci) {
        throw std::runtime_error(
            "batch_locus_em_partitioned: locus_rna_prior_count length must equal n_loci");
    }

    const double*   gp_ptr = locus_gdna_prior_count.data();
    const double*   rp_ptr = locus_rna_prior_count.data();
    const double*   gel_ptr = locus_gdna_eff_lens.data();
    const double*   uac    = unambig_counts.data();
    const double*   tel_ptr = t_eff_lens_arr.data();
    // Checked against the transcript axis, not the locus one: a wrong-axis array would index out of bounds,
    // so it is refused rather than resized.
    if (t_rna_prior_weight_arr.size() != 0 &&
        static_cast<int>(t_rna_prior_weight_arr.shape(0)) != n_transcripts_total) {
        throw std::runtime_error(
            "batch_locus_em_partitioned: t_rna_prior_weight must be empty or of length "
            "n_transcripts_total");
    }
    const double*   trpw_ptr = (t_rna_prior_weight_arr.size() == 0)
        ? nullptr : t_rna_prior_weight_arr.data();

    double* em_out    = em_counts_out.data();
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

    std::vector<double> locus_rna_vec(n_loci, 0.0);
    std::vector<double> locus_gdna_vec(n_loci, 0.0);
    std::vector<double> locus_ig_vec(n_loci, 0.0);
    double* locus_rna_data = locus_rna_vec.data();
    double* locus_gdna_data = locus_gdna_vec.data();
    double* locus_ig_data = locus_ig_vec.data();
    const double log_intergenic_rate = intergenic_rate > 0.0
        ? std::log(intergenic_rate) : -std::numeric_limits<double>::infinity();

    std::vector<LocusProfile> locus_profiles(
        emit_locus_stats ? static_cast<size_t>(n_loci) : 0);

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
                                 rigel::ThreadPool* pool = nullptr) {
            auto locus_t0 = hrclock::now();
            const auto& pv = views[li];
            int n_t = pv.n_transcripts;
            int n_u = pv.n_units;

            if (n_u == 0) {
                locus_rna_data[li] = 0.0;
                locus_gdna_data[li] = 0.0;
                return;
            }

            // 1. Extract sub-problem from partition
            extract_locus_sub_problem_from_partition(
                sub, pv,
                gel_ptr[li],
                unambig_row_sums.data(), tel_ptr, trpw_ptr,
                local_map_vec.data(), local_map_size);
            auto t2 = hrclock::now();

            int nc = sub.n_components;
            size_t n_candidates = sub.t_indices.size();
            int n_local_units = sub.n_local_units;

            // 2. Handle empty sub-problem
            if (n_local_units == 0 || n_candidates == 0) {
                locus_rna_data[li] = 0.0;
                locus_gdna_data[li] = 0.0;
                return;
            }

            // 3. log effective length per component (RNA L̃_t plus gDNA L̃_M).
            const double* log_eff_len_ptr = sub.log_eff_len.data();

            // 4. Build equivalence classes
            auto ec_data = build_equiv_classes(
                sub.offsets.data(),
                sub.t_indices.data(),
                sub.log_liks.data(),
                sub.coverage_wts.data(),
                n_local_units);
            auto t4 = hrclock::now();

            // The intergenic term, when a unit of the locus has an intergenic place.
            IntergenicTerm ig;
            if (sub.ig_idx >= 0) {
                ig.idx = sub.ig_idx;
                ig.log_rate = log_intergenic_rate;
                double n_all = 0.0;
                for (int u = 0; u < n_local_units; ++u)
                    if (sub.offsets[u + 1] > sub.offsets[u]) n_all += 1.0;
                for (int c = 0; c < nc; ++c) n_all += sub.unambig_totals[c];
                ig.n_all = n_all;
            }

            // 5. Grouped aggregate prior and warm start.
            AggregatePrior aggregate_prior{
                nonnegative_finite(gp_ptr[li]),
                nonnegative_finite(rp_ptr[li]),
                sub.gdna_idx,
                sub.has_gdna_candidate,
                sub.component_rna_prior_weight.empty()
                    ? nullptr : sub.component_rna_prior_weight.data(),
            };
            std::vector<double> init_counts(nc);
            compute_grouped_warm_start(
                ec_data,
                sub.unambig_totals.data(),
                aggregate_prior,
                init_counts.data(), nc, warm_start_mode, ig.idx);
            auto t5 = hrclock::now();

            // 6. SQUAREM
            EMResult result = run_squarem(
                ec_data, log_eff_len_ptr,
                sub.unambig_totals.data(),
                aggregate_prior,
                init_counts.data(),
                nc, max_iterations, convergence_delta,
                use_vbem,
                pool, ig);
            auto t6 = hrclock::now();

            // 7. Assign posteriors
            SplitMix64 locus_rng(rng_seed ^ (static_cast<uint64_t>(li) * 0x9e3779b97f4a7c15ULL));
            double locus_rna = 0.0, locus_gdna = 0.0, locus_ig = 0.0;
            int stranded = 0, unrepaired = 0;
            assign_posteriors(
                sub, result.theta.data(),
                log_eff_len_ptr,
                assignment_mode, locus_rng, stranded, unrepaired,
                em_out,
                psum_out, nass_out,
                locus_rna, locus_gdna,
                N_T, N_COLS,
                out_tid_ptr, out_post_ptr, out_ncand_ptr,
                locus_write_offsets[li],
                ig, result.data_count, locus_ig);
            auto t7 = hrclock::now();

            locus_rna_data[li] = locus_rna;
            locus_gdna_data[li] = locus_gdna;
            locus_ig_data[li] = locus_ig;

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
                prof.squarem_step_scale_mean = result.squarem_step_scale_mean;
                prof.squarem_step_scale_max = result.squarem_step_scale_max;
                prof.squarem_backtrack_count = result.squarem_backtrack_count;
                prof.squarem_nonfinite_count = result.squarem_nonfinite_count;
                prof.squarem_grouped_stabilization_fail_count =
                    result.squarem_grouped_stabilization_fail_count;
                prof.assignment_stranded = stranded;
                prof.assignment_unrepaired = unrepaired;
                prof.gdna_eff_len = gel_ptr[li];
                {
                    // Marginal data log-lik at the converged theta (which fixed point is the MLE). With an intergenic
                    // term it is the conditional one, without the Poisson count terms or the prior: not that objective.
                    std::vector<double> lw(static_cast<size_t>(nc));
                    for (int c = 0; c < nc; ++c) {
                        lw[c] = std::log(result.theta[c] + EM_LOG_EPSILON) - log_eff_len_ptr[c];
                    }
                    if (ig.idx >= 0) lw[ig.idx] = intergenic_log_weight(ig, result.data_count);
                    prof.final_data_loglik = marginal_data_loglik(ec_data, lw.data());
                }

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
                prof.extract_us = us(locus_t0, t2);
                prof.build_ec_us = us(t2, t4);
                prof.warm_start_us = us(t4, t5);
                prof.squarem_us = us(t5, t6);
                prof.assign_us = us(t6, t7);
                prof.total_us = us(locus_t0, t7);
            }
        }; // end process_locus

        // ---- Phase 1: mega-loci ----
        {
            LocusSubProblem sub;
            std::vector<int32_t> local_map_vec(local_map_size, -1);

            std::unique_ptr<rigel::ThreadPool> pool;
            if (actual_threads > 1 && mega_end > 0) {
                pool = std::make_unique<rigel::ThreadPool>(actual_threads);
            }

            for (int i = 0; i < mega_end; ++i) {
                int li = locus_order[i];
                process_locus(li, actual_threads, sub, local_map_vec, true, pool.get());
            }
        }

        // ---- Phase 2: work-steal normal loci, one at a time ----
        // One at a time because the order is heaviest first: a chunk would put the heaviest on one worker. Which
        // worker runs a locus changes no number: each locus has its own stream and totals are summed in locus order.
        int n_phase2 = n_loci - mega_end;
        if (n_phase2 > 0) {
            std::atomic<int> next_idx{0};

            auto worker_fn = [&]() {
                LocusSubProblem sub;
                std::vector<int32_t> local_map_vec(local_map_size, -1);

                for (int idx; (idx = next_idx.fetch_add(1, std::memory_order_relaxed)) < n_phase2;) {
                    process_locus(locus_order[mega_end + idx], 1, sub, local_map_vec, false);
                }
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

    // Summed in locus order once the solve is done, never as the workers finish: an arrival-order sum
    // re-associates from run to run, so the reported total would not be one number.
    double total_gdna_em_val = 0.0, total_ig_em_val = 0.0;
    for (int li = 0; li < n_loci; ++li) total_gdna_em_val += locus_gdna_data[li];
    for (int li = 0; li < n_loci; ++li) total_ig_em_val += locus_ig_data[li];

    size_t shape[1] = {static_cast<size_t>(n_loci)};
    auto* rna_copy = new double[n_loci];
    auto* gdna_copy = new double[n_loci];
    std::memcpy(rna_copy, locus_rna_vec.data(), n_loci * sizeof(double));
    std::memcpy(gdna_copy, locus_gdna_vec.data(), n_loci * sizeof(double));

    nb::capsule rna_owner(rna_copy, [](void* p) noexcept { delete[] static_cast<double*>(p); });
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
            d["n_rna_components"] = p.n_transcripts;
            d["n_em_fragments"] = p.n_units;
            d["n_components"] = p.n_components;
            d["n_equiv_classes"] = p.n_equiv_classes;
            d["ec_total_elements"] = p.ec_total_elements;
            d["max_ec_width"] = p.max_ec_width;
            d["max_ec_depth"] = p.max_ec_depth;
            d["squarem_iterations"] = p.squarem_iterations;
            d["estep_threads_used"] = p.estep_threads_used;
            d["is_mega_locus"] = p.is_mega_locus;
            d["squarem_step_scale_mean"] = p.squarem_step_scale_mean;
            d["squarem_step_scale_max"] = p.squarem_step_scale_max;
            d["squarem_backtrack_count"] = p.squarem_backtrack_count;
            d["squarem_nonfinite_count"] = p.squarem_nonfinite_count;
            d["squarem_grouped_stabilization_fail_count"] =
                p.squarem_grouped_stabilization_fail_count;
            d["assignment_stranded"] = p.assignment_stranded;
            d["assignment_unrepaired"] = p.assignment_unrepaired;
            d["gdna_em_effective_length"] = p.gdna_eff_len;
            d["final_data_loglik"] = p.final_data_loglik;
            d["extract_us"] = p.extract_us;
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
            rna_copy, 1, shape, std::move(rna_owner)),
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
        }(),
        total_ig_em_val
    );
}
// ================================================================
// Union-Find connected components
// ================================================================
// Every candidate transcript of a unit is united with the unit's first; find() then gives each transcript its
// component root, relabelled to sequential 0-based IDs. O(N_candidates · α(N_transcripts)).

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


/// The loci: connected components of the global CSR (offsets int64[n_units + 1], t_indices int32[n_candidates]).
/// Returns (n_components, comp_t_offsets, comp_t_flat, comp_u_offsets, comp_u_flat), each component's transcript
/// and unit indices as a sorted CSR.
static nb::tuple connected_components_native(
    nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig>  offsets_arr,
    nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig>  t_indices_arr,
    int32_t n_transcripts)
{
    const int64_t* offsets = offsets_arr.data();
    const int32_t* t_idx   = t_indices_arr.data();
    const int64_t  n_units = static_cast<int64_t>(offsets_arr.shape(0)) - 1;

    if (n_units <= 0 || n_transcripts <= 0) {
        // No component: the same five-part CSR result as any other input, every array empty but the two offsets
        // ([0]), so a caller unpacks one shape
        auto* ct_off = new int64_t[1]{0};
        auto* cu_off = new int64_t[1]{0};
        auto* ct_flat = new int32_t[1];
        auto* cu_flat = new int32_t[1];
        size_t one[1] = {1}, none[1] = {0};
        nb::capsule own_ct_off(ct_off, [](void* p) noexcept { delete[] static_cast<int64_t*>(p); });
        nb::capsule own_ct_flat(ct_flat, [](void* p) noexcept { delete[] static_cast<int32_t*>(p); });
        nb::capsule own_cu_off(cu_off, [](void* p) noexcept { delete[] static_cast<int64_t*>(p); });
        nb::capsule own_cu_flat(cu_flat, [](void* p) noexcept { delete[] static_cast<int32_t*>(p); });
        return nb::make_tuple(
            nb::int_(0),
            nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(ct_off, 1, one, std::move(own_ct_off)),
            nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(ct_flat, 1, none, std::move(own_ct_flat)),
            nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(cu_off, 1, one, std::move(own_cu_off)),
            nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(cu_flat, 1, none, std::move(own_cu_flat))
        );
    }

    UnionFind uf(n_transcripts);

    // Track which transcripts are actually referenced
    std::vector<bool> active(n_transcripts, false);

    for (int64_t u = 0; u < n_units; ++u) {
        int64_t start = offsets[u];
        int64_t end   = offsets[u + 1];
        if (start >= end) continue;

        // Unite every candidate transcript of the unit with its first.
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

    // labels[t]: the 0-based component of an active transcript, -1 otherwise.
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
              "Provides batch_locus_em_partitioned(), the production locus EM\n"
              "entry point used by the Python quantification pipeline.";

    // ---- Partition scatter functions ----
    m.def("build_partition_offsets", &build_partition_offsets,
          nb::arg("g_offsets"),
          nb::arg("locus_units"),
          nb::arg("n_loci"),
          "Build per-locus CSR offsets from global offsets and locus unit lists.");

    m.def("scatter_candidates_f32",
          &scatter_candidates_impl<float>,
          nb::arg("global_arr"), nb::arg("g_offsets"),
          nb::arg("locus_units"), nb::arg("partition_offsets"),
          nb::arg("n_loci"),
          "Scatter per-candidate float32 array into per-locus arrays.");
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

    m.def("scatter_units_f32",
          &scatter_units_impl<float>,
          nb::arg("global_arr"), nb::arg("locus_units"), nb::arg("n_loci"),
          "Scatter per-unit float32 array into per-locus arrays.");
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
          nb::arg("locus_gdna_prior_count"),
          nb::arg("locus_rna_prior_count"),
          nb::arg("locus_gdna_eff_lens"),
          nb::arg("unambig_counts"),
          nb::arg("t_eff_lens"),
          nb::arg("t_rna_prior_weight"),
          nb::arg("warm_start_mode"),
          nb::arg("em_counts_out"),
          nb::arg("posterior_sum_out"),
          nb::arg("n_assigned_out"),
          nb::arg("max_iterations"),
          nb::arg("convergence_delta"),
          nb::arg("use_vbem"),
          nb::arg("assignment_mode"),
          nb::arg("rng_seed"),
          nb::arg("n_transcripts_total"),
          nb::arg("n_splice_strand_cols"),
          nb::arg("n_threads") = 0,
          nb::arg("emit_locus_stats") = false,
          nb::arg("emit_assignments") = false,
          nb::arg("intergenic_rate") = 0.0,
          "Run locus EM from per-locus partition data.\n\n"
          "Accepts a list of 8-tuples (one per locus) of partition arrays.\n"
          "Returns (total_gdna_em, locus_rna_total, locus_gdna, locus_stats,\n"
          " out_winner_tid, out_winner_post, out_n_candidates, total_intergenic_em).\n"
          "\n"
          "locus_rna_total is the sum of posteriors assigned to any\n"
          "transcript-like component (annotated transcripts and synthetic spans);\n"
          "splitting into annotated transcripts vs synthetic spans is done by\n"
          "the Python caller.");

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

    // The log floor, exported for the grouped-prior-update gate.
    m.attr("EM_LOG_EPSILON")         = EM_LOG_EPSILON;

    // The mode codes, the one definition rigel.estimator reads; batch_locus_em_partitioned refuses any other.
    m.attr("ASSIGN_FRACTIONAL")      = ASSIGN_FRACTIONAL;
    m.attr("ASSIGN_SAMPLE")          = ASSIGN_SAMPLE;
    m.attr("WARM_START_COVERAGE")    = WARM_START_COVERAGE;
    m.attr("WARM_START_PRIOR")       = WARM_START_PRIOR;
    m.attr("WARM_START_UNIFORM")     = WARM_START_UNIFORM;
    m.attr("EM_WINNER_INTERGENIC")   = EM_WINNER_INTERGENIC;

    // Test-only: the static `apply_grouped_prior_update`, exposed so its conservation identity can be tested (for
    // given raw_counts, the RNA components sum to rna_count + rna_prior). The identity holds per call, not end to
    // end: a different rna_prior changes theta, hence the next E-step's raw_counts and the converged gDNA total.
    // An empty carried_state / rna_prior_weight means nullptr, as on batch_locus_em_partitioned's lane.
    m.def("_apply_grouped_prior_update_test",
          [](f64_1d raw_counts,
             f64_1d carried_state,
             f64_1d rna_prior_weight,
             double gdna_prior_fragments,
             double rna_prior_fragments,
             int    gdna_index,
             bool   has_gdna_candidate) {
              const size_t n = raw_counts.shape(0);
              if (carried_state.size() != 0 && carried_state.shape(0) != n) {
                  throw std::runtime_error(
                      "_apply_grouped_prior_update_test: carried_state must be empty or "
                      "the same length as raw_counts");
              }
              if (rna_prior_weight.size() != 0 && rna_prior_weight.shape(0) != n) {
                  throw std::runtime_error(
                      "_apply_grouped_prior_update_test: rna_prior_weight must be empty or "
                      "the same length as raw_counts");
              }
              AggregatePrior aggregate_prior;
              aggregate_prior.gdna_prior_fragments = gdna_prior_fragments;
              aggregate_prior.rna_prior_fragments  = rna_prior_fragments;
              aggregate_prior.gdna_index           = gdna_index;
              aggregate_prior.has_gdna_candidate   = has_gdna_candidate;
              aggregate_prior.component_rna_prior_weight =
                  (rna_prior_weight.size() == 0) ? nullptr : rna_prior_weight.data();

              double* out = new double[n];
              apply_grouped_prior_update(
                  raw_counts.data(),
                  (carried_state.size() == 0) ? nullptr : carried_state.data(),
                  aggregate_prior,
                  out,
                  static_cast<int>(n));
              nb::capsule owner(out, [](void* p) noexcept {
                  delete[] static_cast<double*>(p);
              });
              return nb::ndarray<nb::numpy, double, nb::ndim<1>>(out, {n}, owner);
          },
          nb::arg("raw_counts"),
          nb::arg("carried_state"),
          nb::arg("rna_prior_weight"),
          nb::arg("gdna_prior_fragments"),
          nb::arg("rna_prior_fragments"),
          nb::arg("gdna_index"),
          nb::arg("has_gdna_candidate"),
          "Run one grouped prior update and return out_counts (test-only).\n\n"
          "Empty carried_state / rna_prior_weight mean nullptr.");

    // Test-only: fast_exp over an array, for accuracy tests.
    m.def("_fast_exp_test_array",
          [](nb::ndarray<double, nb::ndim<1>, nb::c_contig> inputs) {
              const size_t n = inputs.shape(0);
              const double* src = inputs.data();
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
