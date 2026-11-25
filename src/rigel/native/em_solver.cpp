// em_solver.cpp — the locus EM's native module, `_em_impl`.
//
// One call solves every locus: its sub-problem (`em_locus.h`), its equivalence classes and E-step (`em_estep.h`),
// MAP-EM or VBEM under SQUAREM with the grouped calibration prior (`em_squarem.h`, `em_prior.h`), and the posterior
// assignment (`em_assign.h`). Beside it, the scatter that builds the loci's partitions (`em_scatter.h`) and the
// connected components that find the loci (`connected_components.h`).

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "connected_components.h"
#include "constants.h"
#include "em_assign.h"
#include "em_estep.h"
#include "em_locus.h"
#include "em_prior.h"
#include "em_scatter.h"
#include "em_squarem.h"
#include "ndarray_util.h"
#include "thread_pool.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;
using namespace rigel;

namespace {

using f64_1d = Array1<double>;
using f64_2d = nb::ndarray<const double, nb::ndim<2>, nb::c_contig>;
using f64_1d_out = nb::ndarray<double, nb::ndim<1>, nb::c_contig>;
using f64_2d_out = nb::ndarray<double, nb::ndim<2>, nb::c_contig>;
using u8_1d = Array1<uint8_t>;

using hrclock = std::chrono::steady_clock;

// One locus's statistics (`emit_locus_stats`); a locus with nothing to solve keeps locus_idx -1 and has no row.
struct LocusProfile {
    int locus_idx = -1;
    int n_transcripts = 0;
    int n_units = 0;
    int n_components = 0;        // n_t + 1, or n_t + 2 with the intergenic term
    int n_equiv_classes = 0;
    int64_t ec_total_elements = 0; // sum of n*k across all ECs: the E-step's entries, each distinct row once
    int max_ec_width = 0;        // max k across ECs
    int max_ec_depth = 0;        // max n (distinct rows) across ECs
    int estep_threads_used = 0;
    bool is_mega_locus = false;
    SquaremStats squarem;
    LocusAssignment assignment;
    double gdna_eff_len = 1.0;
    // The marginal data log-likelihood at the converged theta (which fixed point is the MLE). With an intergenic term
    // it is the conditional one, without the Poisson count terms or the prior: not that objective.
    double final_data_loglik = 0.0;
    // Sub-phase wall times in microseconds
    double extract_us = 0.0, build_ec_us = 0.0, warm_start_us = 0.0, squarem_us = 0.0, assign_us = 0.0, total_us = 0.0;

    void describe_classes(const std::vector<EmEquivClass>& ec_data) {
        n_equiv_classes = static_cast<int>(ec_data.size());
        for (const auto& ec : ec_data) {
            ec_total_elements += static_cast<int64_t>(ec.n) * ec.k;
            max_ec_width = std::max(max_ec_width, ec.k);
            max_ec_depth = std::max(max_ec_depth, ec.n);
        }
    }

    nb::dict to_dict() const {
        nb::dict d;
        d["locus_idx"] = locus_idx;
        d["n_rna_components"] = n_transcripts;
        d["n_em_fragments"] = n_units;
        d["n_components"] = n_components;
        d["n_equiv_classes"] = n_equiv_classes;
        d["ec_total_elements"] = ec_total_elements;
        d["max_ec_width"] = max_ec_width;
        d["max_ec_depth"] = max_ec_depth;
        d["squarem_iterations"] = squarem.iterations;
        d["estep_threads_used"] = estep_threads_used;
        d["is_mega_locus"] = is_mega_locus;
        d["squarem_step_scale_mean"] = squarem.step_scale_mean();
        d["squarem_step_scale_max"] = squarem.step_scale_max;
        d["squarem_backtrack_count"] = squarem.backtracks;
        d["squarem_nonfinite_count"] = squarem.nonfinite;
        d["squarem_grouped_stabilization_fail_count"] = squarem.stabilization_fails;
        d["assignment_stranded"] = assignment.stranded;
        d["assignment_unrepaired"] = assignment.unrepaired;
        d["gdna_em_effective_length"] = gdna_eff_len;
        d["final_data_loglik"] = final_data_loglik;
        d["extract_us"] = extract_us;
        d["build_ec_us"] = build_ec_us;
        d["warm_start_us"] = warm_start_us;
        d["squarem_us"] = squarem_us;
        d["assign_us"] = assign_us;
        d["total_us"] = total_us;
        return d;
    }
};

// One thread's scratch, reused across the loci it solves.
struct LocusScratch {
    LocusSubProblem sub;
    std::vector<int32_t> local_map;   // global transcript -> local component, -1 between loci
    explicit LocusScratch(int n_transcripts) : local_map(static_cast<size_t>(n_transcripts), -1) {}
};

// The loci's partitions (each a `LocusPartition`, its arrays read by name) and transcript lists, as views, each
// checked here once so the
// EM reads every candidate unguarded: it holds a unit; its offsets start at 0 and never fall; its candidate and unit
// arrays are as long as its offsets say; its transcripts are among the `n_transcripts_total`; every candidate is one
// of its transcripts; and every count column is one of the `n_cols`. Every array must be of its exact dtype: the views
// outlive the casts that checked them.
std::vector<PartitionView> read_partitions(nb::list partitions, nb::list locus_transcript_indices,
                                           int n_transcripts_total, int n_cols) {
    using i64 = Array1<int64_t>;
    using i32 = Array1<int32_t>;
    using f32 = Array1<float>;
    using u8 = Array1<uint8_t>;
    std::vector<PartitionView> views(nb::len(partitions));
    if (nb::len(locus_transcript_indices) != views.size())
        throw std::invalid_argument("batch_locus_em_partitioned: one transcript list per locus partition");
    // the locus that holds each transcript: loci are disjoint, and a candidate is checked against its own locus's list
    std::vector<int64_t> named_by(static_cast<size_t>(std::max(n_transcripts_total, 0)), -1);
    for (size_t i = 0; i < views.size(); ++i) {
        const auto refuse = [i](const std::string& what) {
            throw std::invalid_argument("batch_locus_em_partitioned: locus " + std::to_string(i) + " " + what);
        };
        // by name: the four float arrays are two pairs of one shape and dtype, which a position would not tell apart
        const nb::object partition = partitions[i];
        const auto field = [&](const char* name) {
            nb::object value = nb::getattr(partition, name, nb::none());
            if (value.is_none()) refuse(std::string("has no array ") + name);
            return value;
        };
        const auto offsets = nb::cast<i64>(field("offsets"), false);
        const auto t_indices = nb::cast<i32>(field("t_indices"), false);
        const auto log_liks = nb::cast<f32>(field("log_liks"), false);
        const auto coverage_wts = nb::cast<f32>(field("coverage_weights"), false);
        const auto count_cols = nb::cast<u8>(field("count_cols"), false);
        const auto gdna_excluded = nb::cast<u8>(field("gdna_excluded"), false);
        const auto gdna_log_liks = nb::cast<f32>(field("gdna_log_liks"), false);
        const auto intergenic_log_liks = nb::cast<f32>(field("intergenic_log_liks"), false);
        const auto transcripts = nb::cast<i32>(locus_transcript_indices[i], false);
        const int64_t n_units = static_cast<int64_t>(offsets.shape(0)) - 1;
        if (n_units < 1) refuse("holds no unit");
        const int64_t* off = offsets.data();
        if (off[0] != 0) refuse("has offsets that do not start at 0");
        for (int64_t u = 0; u < n_units; ++u)
            if (off[u + 1] < off[u]) refuse("has offsets that fall");
        const size_t n_candidates = static_cast<size_t>(off[n_units]);
        if (t_indices.shape(0) != n_candidates || log_liks.shape(0) != n_candidates ||
            coverage_wts.shape(0) != n_candidates || count_cols.shape(0) != n_candidates)
            refuse("has a candidate array whose length is not its offsets' last");
        const size_t units = static_cast<size_t>(n_units);
        if (gdna_excluded.shape(0) != units || gdna_log_liks.shape(0) != units || intergenic_log_liks.shape(0) != units)
            refuse("has a unit array whose length is not its unit count");
        for (size_t k = 0; k < transcripts.shape(0); ++k) {
            const int32_t t = transcripts.data()[k];
            if (t < 0 || t >= n_transcripts_total)
                refuse("names transcript " + std::to_string(t) + ", outside the " + std::to_string(n_transcripts_total));
            const int64_t owner = named_by[static_cast<size_t>(t)];
            if (owner == static_cast<int64_t>(i)) refuse("names transcript " + std::to_string(t) + " twice");
            if (owner >= 0)
                refuse("names transcript " + std::to_string(t) + ", which locus " + std::to_string(owner) + " holds");
            named_by[static_cast<size_t>(t)] = static_cast<int64_t>(i);
        }
        for (size_t j = 0; j < n_candidates; ++j) {
            const int32_t t = t_indices.data()[j];
            if (t < 0 || t >= n_transcripts_total || named_by[static_cast<size_t>(t)] != static_cast<int64_t>(i))
                refuse("has candidate transcript " + std::to_string(t) + ", which is not one of its transcripts");
            if (count_cols.data()[j] >= n_cols)
                refuse("has count column " + std::to_string(count_cols.data()[j]) + ", past its " +
                       std::to_string(n_cols) + " columns");
        }
        PartitionView& v = views[i];
        v.offsets = off;
        v.t_indices = t_indices.data();
        v.log_liks = log_liks.data();
        v.coverage_wts = coverage_wts.data();
        v.count_cols = count_cols.data();
        v.gdna_excluded = gdna_excluded.data();
        v.gdna_log_liks = gdna_log_liks.data();
        v.intergenic_log_liks = intergenic_log_liks.data();
        v.n_units = static_cast<int>(n_units);
        v.n_candidates = off[n_units];
        v.transcript_indices = transcripts.data();
        v.n_transcripts = static_cast<int>(transcripts.shape(0));
    }
    return views;
}

// An undefined assignment code is refused: it would index the per-unit draw only ASSIGN_SAMPLE fills.
void check_assignment_mode(int assignment_mode) {
    if (assignment_mode != ASSIGN_FRACTIONAL && assignment_mode != ASSIGN_SAMPLE) {
        throw std::invalid_argument(
            "batch_locus_em_partitioned: assignment_mode must be ASSIGN_FRACTIONAL or ASSIGN_SAMPLE; got " +
            std::to_string(assignment_mode));
    }
}

void check_length(const f64_1d& a, size_t n, const char* what) {
    if (a.shape(0) != n)
        throw std::runtime_error(std::string("batch_locus_em_partitioned: ") + what + " length must equal n_loci");
}

// The loci heaviest first (work = transcripts × units), and how many lead the order as mega-loci: those with at least
// a thread's fair share of the work, solved one at a time with the E-step on every thread.
std::vector<int> schedule_loci(const std::vector<PartitionView>& views, int n_threads, int& mega_end) {
    const int n_loci = static_cast<int>(views.size());
    std::vector<int64_t> work(static_cast<size_t>(n_loci));
    int64_t total_work = 0;
    for (int li = 0; li < n_loci; ++li) {
        work[li] = static_cast<int64_t>(views[li].n_transcripts) * static_cast<int64_t>(views[li].n_units);
        total_work += work[li];
    }
    std::vector<int> order(static_cast<size_t>(n_loci));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return work[a] > work[b]; });
    const int64_t fair_share = (n_threads > 1 && total_work > 0) ? total_work / n_threads : total_work + 1;
    mega_end = 0;
    while (mega_end < n_loci && work[order[mega_end]] >= fair_share) ++mega_end;
    return order;
}

// The batch: every locus's inputs, the accumulators the assignment adds into, and each locus's totals.
struct BatchEM {
    const std::vector<PartitionView>& views;
    const double* gdna_prior;      // per locus, pseudo-fragments
    const double* rna_prior;
    const double* gdna_eff_len;    // per locus, the gDNA component's capture-contracted effective length
    const double* t_eff_len;       // per transcript
    const double* t_rna_prior_weight;  // per transcript, or nullptr: the RNA prior follows the evidence
    std::vector<double> t_unambig_totals;   // per transcript, over its count columns
    int max_iterations, assignment_mode;
    double convergence_delta;
    bool use_vbem;
    uint64_t rng_seed;
    double log_intergenic_rate;
    bool emit_locus_stats;
    AssignmentSinks sinks;
    std::vector<int> unit_offsets;   // each locus's first unit in the per-unit annotation
    std::vector<double> locus_rna, locus_gdna, locus_ig;
    std::vector<LocusProfile> profiles;

    // One locus, end to end: extract, classes, warm start, SQUAREM, assignment. `pool` runs the E-step on every
    // thread (a mega-locus), or is nullptr.
    void solve(int li, LocusScratch& scratch, ThreadPool* pool, bool mega) {
        auto now = [&] { return emit_locus_stats ? hrclock::now() : hrclock::time_point{}; };
        const auto t0 = now();
        const PartitionView& pv = views[li];
        LocusSubProblem& sub = scratch.sub;
        extract_locus_sub_problem(sub, pv, gdna_eff_len[li], t_unambig_totals.data(), t_eff_len, t_rna_prior_weight,
                                  scratch.local_map);
        const auto t_extract = now();
        const IntergenicTerm ig = intergenic_term(sub, log_intergenic_rate);
        const auto ec_data = build_equiv_classes(sub.offsets.data(), sub.t_indices.data(), sub.log_liks.data(),
                                                 sub.coverage_wts.data(), sub.n_local_units, ig.idx);
        const auto t_classes = now();

        const AggregatePrior prior{gdna_prior[li], rna_prior[li], sub.gdna_idx,
                                   sub.has_gdna_candidate,
                                   sub.component_rna_prior_weight.empty() ? nullptr
                                                                          : sub.component_rna_prior_weight.data()};
        std::vector<double> init_counts(sub.n_components);
        compute_grouped_warm_start(ec_data, sub.unambig_totals.data(), prior, init_counts.data(), sub.n_components,
                                   ig.idx);
        const auto t_warm = now();
        const EMResult em = run_squarem(ec_data, sub.log_eff_len.data(), sub.unambig_totals.data(), prior,
                                        init_counts.data(), sub.n_components, max_iterations, convergence_delta,
                                        use_vbem, pool, ig);
        const auto t_em = now();
        SplitMix64 rng = SplitMix64::stream(rng_seed, static_cast<uint64_t>(li));
        const LocusAssignment assigned = assign_posteriors(sub, em.theta.data(), ig, em.data_count, assignment_mode,
                                                           rng, sinks, unit_offsets[li]);
        const auto t_assign = now();
        locus_rna[li] = assigned.rna;
        locus_gdna[li] = assigned.gdna;
        locus_ig[li] = assigned.intergenic;
        if (!emit_locus_stats) return;

        auto us = [](auto a, auto b) { return std::chrono::duration<double, std::micro>(b - a).count(); };
        LocusProfile& prof = profiles[li];
        prof.locus_idx = li;
        prof.n_transcripts = pv.n_transcripts;
        prof.n_units = pv.n_units;
        prof.n_components = sub.n_components;
        prof.describe_classes(ec_data);
        prof.estep_threads_used = pool != nullptr ? pool->n_threads() : 1;
        prof.is_mega_locus = mega;
        prof.squarem = em.stats;
        prof.assignment = assigned;
        prof.gdna_eff_len = gdna_eff_len[li];
        std::vector<double> lw(static_cast<size_t>(sub.n_components));
        converged_log_weights(em.theta.data(), sub.log_eff_len.data(), sub.n_components, ig, em.data_count,
                              lw.data());
        prof.final_data_loglik = marginal_data_loglik(ec_data, lw.data());
        prof.extract_us = us(t0, t_extract);
        prof.build_ec_us = us(t_extract, t_classes);
        prof.warm_start_us = us(t_classes, t_warm);
        prof.squarem_us = us(t_warm, t_em);
        prof.assign_us = us(t_em, t_assign);
        prof.total_us = us(t0, t_assign);
    }

    // Every locus: the mega-loci first, one at a time with the E-step on the pool, then the rest pulled one at a time
    // by the pool's threads. One at a time because the order is heaviest first: a chunk would put the heaviest on one
    // worker. Which worker runs a locus changes no number: each locus has its own stream and its own totals.
    void solve_all(int n_threads) {
        int mega_end = 0;
        const std::vector<int> order = schedule_loci(views, n_threads, mega_end);
        const int n_loci = static_cast<int>(order.size());
        const int n_transcripts = sinks.n_transcripts;
        std::unique_ptr<ThreadPool> pool = n_threads > 1 ? std::make_unique<ThreadPool>(n_threads) : nullptr;
        {
            LocusScratch scratch(n_transcripts);
            for (int i = 0; i < mega_end; ++i) solve(order[i], scratch, pool.get(), true);
        }
        std::atomic<int> next{mega_end};
        auto worker = [&](int) {
            LocusScratch scratch(n_transcripts);
            for (int i; (i = next.fetch_add(1, std::memory_order_relaxed)) < n_loci;)
                solve(order[i], scratch, nullptr, false);
        };
        if (pool) pool->run_parallel(worker);
        else worker(0);
    }
};

// Each transcript's unambiguous total over its count columns.
std::vector<double> unambig_totals(const f64_2d& unambig_counts, int n_transcripts, int n_cols) {
    std::vector<double> sums(n_transcripts, 0.0);
    for (int t = 0; t < n_transcripts; ++t) {
        double s = 0.0;
        for (int c = 0; c < n_cols; ++c) s += unambig_counts.data()[t * n_cols + c];
        sums[t] = s;
    }
    return sums;
}

// Each locus's first unit in the per-unit annotation, and the total.
std::vector<int> first_units(const std::vector<PartitionView>& views, int& total_units) {
    std::vector<int> first(views.size(), 0);
    total_units = 0;
    for (size_t li = 0; li < views.size(); ++li) {
        first[li] = total_units;
        total_units += views[li].n_units;
    }
    return first;
}

// The per-unit annotation arrays, or None when not requested or empty.
template <typename T>
nb::object annotation(std::vector<T>& v) {
    if (v.empty()) return nb::none();
    return vec_to_ndarray(std::move(v));
}

nb::tuple batch_locus_em_partitioned(
    nb::list partitions,                // per locus, its partition (`LocusPartition`), its arrays read by name
    nb::list locus_transcript_indices,  // per locus, int32 transcript indices
    f64_1d locus_gdna_prior_count,      // per locus, the grouped calibration prior's pseudo-fragments
    f64_1d locus_rna_prior_count,
    f64_1d locus_gdna_eff_lens,         // per locus, the gDNA component's capture-contracted effective length
    f64_2d unambig_counts,              // per transcript and count column
    f64_1d t_eff_lens_arr,              // per transcript, the EM's capture-contracted effective length
    f64_1d t_rna_prior_weight_arr,      // per transcript, or empty: the RNA prior follows the current evidence
    f64_2d_out em_counts_out,           // the accumulators, added into
    f64_1d_out posterior_sq_sum_out,
    f64_1d_out posterior_mass_out,
    int max_iterations,
    double convergence_delta,
    bool use_vbem,
    int assignment_mode,
    uint64_t rng_seed,
    int n_threads,
    bool emit_locus_stats,
    bool emit_assignments,
    double intergenic_rate,             // fragments per start in the efficiencies' unit; 0 reads no intergenic term
    u8_1d t_is_synthetic_arr)           // per transcript, or empty: with it, every unit's class masses are returned
{
    check_assignment_mode(assignment_mode);
    // The count matrix is the transcript axis and the column axis; every per-transcript array is checked against it
    // here, once, so the EM reads and writes every one unguarded.
    const int n_transcripts_total = static_cast<int>(unambig_counts.shape(0));
    const int n_cols = static_cast<int>(unambig_counts.shape(1));
    const auto per_transcript = [&](size_t len, const char* name) {
        if (static_cast<int>(len) != n_transcripts_total)
            throw std::invalid_argument(std::string("batch_locus_em_partitioned: ") + name +
                                        " is not one per transcript");
    };
    per_transcript(t_eff_lens_arr.shape(0), "t_eff_lens");
    per_transcript(posterior_sq_sum_out.shape(0), "posterior_sq_sum_out");
    per_transcript(posterior_mass_out.shape(0), "posterior_mass_out");
    per_transcript(em_counts_out.shape(0), "em_counts_out");
    if (static_cast<int>(em_counts_out.shape(1)) != n_cols)
        throw std::invalid_argument("batch_locus_em_partitioned: em_counts_out's columns are not unambig_counts'");
    const std::vector<PartitionView> views =
        read_partitions(partitions, locus_transcript_indices, n_transcripts_total, n_cols);
    const size_t n_loci = views.size();
    check_length(locus_gdna_eff_lens, n_loci, "locus_gdna_eff_lens");
    check_length(locus_gdna_prior_count, n_loci, "locus_gdna_prior_count");
    check_length(locus_rna_prior_count, n_loci, "locus_rna_prior_count");
    // A prior count is a count of pseudo-fragments: one that is not finite and nonnegative is the caller's error.
    for (const f64_1d* prior : {&locus_gdna_prior_count, &locus_rna_prior_count})
        for (size_t li = 0; li < n_loci; ++li)
            if (!(std::isfinite(prior->data()[li]) && prior->data()[li] >= 0.0))
                throw std::invalid_argument(
                    "batch_locus_em_partitioned: a locus prior count is negative or not finite");
    // Checked against the transcript axis, not the locus one: a wrong-axis array would index out of bounds, so it is
    // refused rather than resized.
    if (t_rna_prior_weight_arr.size() != 0 &&
        static_cast<int>(t_rna_prior_weight_arr.shape(0)) != n_transcripts_total) {
        throw std::runtime_error(
            "batch_locus_em_partitioned: t_rna_prior_weight must be empty or of length n_transcripts_total");
    }
    if (t_is_synthetic_arr.size() != 0 && static_cast<int>(t_is_synthetic_arr.shape(0)) != n_transcripts_total) {
        throw std::runtime_error(
            "batch_locus_em_partitioned: t_is_synthetic must be empty or of length n_transcripts_total");
    }

    int total_units = 0;
    std::vector<int> unit_offsets = first_units(views, total_units);
    const bool annotate = emit_assignments && total_units > 0;
    std::vector<int32_t> winner_tid(annotate ? total_units : 0, EM_WINNER_NONE);
    std::vector<float> winner_post(annotate ? total_units : 0, 0.0f);
    std::vector<int32_t> n_candidates(annotate ? total_units : 0, 0);
    // Every unit's class masses, in the per-unit annotation's row order, when the caller sorts the transcripts.
    const bool classes = t_is_synthetic_arr.size() != 0 && total_units > 0;
    std::vector<double> class_mass(classes ? static_cast<size_t>(N_EM_CLASSES) * total_units : 0, 0.0);

    BatchEM batch{views,
                  locus_gdna_prior_count.data(),
                  locus_rna_prior_count.data(),
                  locus_gdna_eff_lens.data(),
                  t_eff_lens_arr.data(),
                  t_rna_prior_weight_arr.size() == 0 ? nullptr : t_rna_prior_weight_arr.data(),
                  unambig_totals(unambig_counts, n_transcripts_total, n_cols),
                  max_iterations,
                  assignment_mode,
                  convergence_delta,
                  use_vbem,
                  rng_seed,
                  intergenic_rate > 0.0 ? std::log(intergenic_rate) : NEG_INF,
                  emit_locus_stats,
                  AssignmentSinks{em_counts_out.data(), posterior_sq_sum_out.data(), posterior_mass_out.data(),
                                  n_transcripts_total, n_cols, annotate ? winner_tid.data() : nullptr,
                                  annotate ? winner_post.data() : nullptr, annotate ? n_candidates.data() : nullptr,
                                  classes ? class_mass.data() : nullptr,
                                  classes ? t_is_synthetic_arr.data() : nullptr},
                  std::move(unit_offsets),
                  std::vector<double>(n_loci, 0.0),
                  std::vector<double>(n_loci, 0.0),
                  std::vector<double>(n_loci, 0.0),
                  std::vector<LocusProfile>(emit_locus_stats ? n_loci : 0)};
    {
        nb::gil_scoped_release release;
        batch.solve_all(thread_budget(n_threads));
    }

    // Summed in locus order once the solve is done, never as the workers finish: an arrival-order sum re-associates
    // from run to run, so the reported total would not be one number.
    double total_gdna = 0.0, total_ig = 0.0;
    for (size_t li = 0; li < n_loci; ++li) total_gdna += batch.locus_gdna[li];
    for (size_t li = 0; li < n_loci; ++li) total_ig += batch.locus_ig[li];
    // In locus order, so two runs write the same rows in the same order.
    nb::list stats;
    for (const LocusProfile& p : batch.profiles)
        if (p.locus_idx >= 0) stats.append(p.to_dict());
    nb::object class_masses = nb::none();
    if (classes) {
        class_masses = vec_to_ndarray(std::move(class_mass),
                                      {static_cast<size_t>(total_units), static_cast<size_t>(N_EM_CLASSES)});
    }
    return nb::make_tuple(total_gdna, vec_to_ndarray(std::move(batch.locus_rna)),
                          vec_to_ndarray(std::move(batch.locus_gdna)), stats, annotation(winner_tid),
                          annotation(winner_post), annotation(n_candidates), total_ig, class_masses);
}

}  // namespace

NB_MODULE(_em_impl, m) {
    m.doc() = "The locus EM: every locus solved in one call, the partition scatter that builds the loci's arrays, and "
              "the connected components that find the loci.";

    m.def("build_partition_offsets", &build_partition_offsets, nb::arg("g_offsets"), nb::arg("locus_units"),
          nb::arg("n_loci"), "Build per-locus CSR offsets from global offsets and locus unit lists.");
    m.def("scatter_candidates_f32", &scatter_candidates<float>, nb::arg("global_arr"), nb::arg("g_offsets"),
          nb::arg("locus_units"), nb::arg("partition_offsets"), nb::arg("n_loci"),
          "Scatter per-candidate float32 array into per-locus arrays.");
    m.def("scatter_candidates_i32", &scatter_candidates<int32_t>, nb::arg("global_arr"), nb::arg("g_offsets"),
          nb::arg("locus_units"), nb::arg("partition_offsets"), nb::arg("n_loci"),
          "Scatter per-candidate int32 array into per-locus arrays.");
    m.def("scatter_candidates_u8", &scatter_candidates<uint8_t>, nb::arg("global_arr"), nb::arg("g_offsets"),
          nb::arg("locus_units"), nb::arg("partition_offsets"), nb::arg("n_loci"),
          "Scatter per-candidate uint8 array into per-locus arrays.");
    m.def("scatter_units_f32", &scatter_units<float>, nb::arg("global_arr"), nb::arg("locus_units"), nb::arg("n_loci"),
          "Scatter per-unit float32 array into per-locus arrays.");
    m.def("scatter_units_u8", &scatter_units<uint8_t>, nb::arg("global_arr"), nb::arg("locus_units"),
          nb::arg("n_loci"), "Scatter per-unit uint8 array into per-locus arrays.");
    m.def("scatter_units_i64", &scatter_units<int64_t>, nb::arg("global_arr"), nb::arg("locus_units"),
          nb::arg("n_loci"), "Scatter per-unit int64 array into per-locus arrays.");

    m.def("batch_locus_em_partitioned", &batch_locus_em_partitioned,
          nb::arg("partitions"),
          nb::arg("locus_transcript_indices"),
          nb::arg("locus_gdna_prior_count"),
          nb::arg("locus_rna_prior_count"),
          nb::arg("locus_gdna_eff_lens"),
          nb::arg("unambig_counts"),
          nb::arg("t_eff_lens"),
          nb::arg("t_rna_prior_weight"),
          nb::arg("em_counts_out").noconvert(),
          nb::arg("posterior_sq_sum_out").noconvert(),
          nb::arg("posterior_mass_out").noconvert(),
          nb::arg("max_iterations"),
          nb::arg("convergence_delta"),
          nb::arg("use_vbem"),
          nb::arg("assignment_mode"),
          nb::arg("rng_seed"),
          nb::arg("n_threads") = 0,
          nb::arg("emit_locus_stats") = false,
          nb::arg("emit_assignments") = false,
          nb::arg("intergenic_rate") = 0.0,
          nb::arg("t_is_synthetic"),
          "Run locus EM from per-locus partition data.\n\n"
          "Accepts one partition per locus (a LocusPartition), each array read by its name.\n"
          "Returns (total_gdna_em, locus_rna, locus_gdna, locus_stats,\n"
          " out_winner_tid, out_winner_post, out_n_candidates, total_intergenic_em, class_masses).\n"
          "\n"
          "locus_rna is the sum of posteriors assigned to any\n"
          "transcript-like component (annotated transcripts and synthetic spans);\n"
          "splitting into annotated transcripts vs synthetic spans is done by\n"
          "the Python caller.");

    m.def("connected_components", &connected_components, nb::arg("offsets"), nb::arg("t_indices"),
          nb::arg("n_transcripts"),
          "Find connected components of the fragment→transcript overlap graph.\n\n"
          "Uses union-find with path halving and union by rank.\n"
          "All candidates are direct transcript indices.\n"
          "Returns (n_comp, comp_t_offsets, comp_t_flat, comp_u_offsets,\n"
          "comp_u_flat) where the CSR pairs (offsets, flat) give sorted\n"
          "transcript indices and unit indices for each component.");

    // The assignment codes, the one definition rigel.estimator reads; batch_locus_em_partitioned refuses any other.
    m.attr("ASSIGN_FRACTIONAL") = ASSIGN_FRACTIONAL;
    m.attr("ASSIGN_SAMPLE") = ASSIGN_SAMPLE;
    m.attr("EM_WINNER_GDNA") = EM_WINNER_GDNA;
    m.attr("EM_WINNER_INTERGENIC") = EM_WINNER_INTERGENIC;
    nb::list class_names;
    for (int c = 0; c < N_EM_CLASSES; ++c) class_names.append(nb::none());
    class_names[EM_CLASS_ANNOTATED] = "annotated";
    class_names[EM_CLASS_NASCENT] = "nascent";
    class_names[EM_CLASS_GDNA] = "gdna";
    class_names[EM_CLASS_INTERGENIC] = "intergenic";
    m.attr("EM_CLASSES") = nb::tuple(class_names);

    rigel::bind_testing(m.def_submodule("testing", "The EM's kernels bound for the tests."));
}
