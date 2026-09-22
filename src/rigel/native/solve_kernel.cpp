// solve_kernel.cpp — THE SOLVE's one native module (`_solve_impl`): the calibration sweep's block pipeline in ONE call
// per sweep, a pool of threads over the locus blocks, with ψ's kernel (`psi_kernel.h`) and the composition transfer's
// (`transfer_kernel.h`) as its pieces, and every binding the gates read the same code through.
//
//   solve_blocks     THE SWEEP: every locus block of the chain — the prior rows (the fitted gDNA arm read from the
//                    landscape's curve, the intron factory's rows), the SELF-SOLVE ψ, the own-evidence precision, the
//                    message LAYER (the builders, the two passes, the solve — unless the policy is silent), the FINAL ψ,
//                    the write-back on the owned slots, `has_composition`, the counts of
//                    the backbone's assertions — each block on one thread of a pool, on that thread's arena; the belief
//                    written in place, the diagnostics capture on request.
//                    Bit-identical at every thread count: a block's arithmetic is the same whatever thread takes it, and
//                    nothing is reduced across blocks but integer counts (`sweep.solve_chain`).
//   psi_solve        ψ over a slot list, the slots pulled one at a time by the same pool — the pre-sweep solve
//                    (`region_geometry.init_beliefs`) and the gates (`simplex_logodds._solve_regions_logodds_all`).
//   psi_cube, posterior_median, compose, gdna_arm — ψ's pieces for the gates (`simplex_logodds`).
//   transfer_prepare, transfer_pass, transfer_solve — the policy's three kernels on tables the call allocates and
//                    returns, for the transfer gates (`tests/calibration/_transfer_harness.py`).
//   rows             the row constructors of `transfer_rows.h`, the builders' flag predicates, and the constructions the
//                    block pipeline absorbed — the factory rows, the factor precision, the strand evidence, the
//                    log-gamma — bound for the gates. Nothing in `src/` reads them.
//
// One implementation of every arithmetic, in one module: the gates and the production path cannot drift apart.

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>

#include "psi_kernel.h"
#include "thread_pool.h"
#include "transfer_kernel.h"

namespace nb = nanobind;
using namespace transfer_kernel;
namespace P = psi_kernel;

namespace {

using Vec = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using Mat = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using Cube = nb::ndarray<double, nb::ndim<3>, nb::c_contig, nb::device::cpu>;
using BoolVec = nb::ndarray<bool, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using BoolMat = nb::ndarray<bool, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using IdxVec = nb::ndarray<int64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using IdxMat = nb::ndarray<int32_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using KindMat = nb::ndarray<int8_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using FlagVec = nb::ndarray<uint16_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

static_assert(sizeof(bool) == 1, "the bit tables are byte arrays");

// ---- numpy arrays owned by the module: a vector handed to Python without a copy --------------------------

template <class T>
nb::ndarray<nb::numpy, T, nb::ndim<1>> own1(std::vector<T>&& v) {
    auto* p = new std::vector<T>(std::move(v));
    nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<T>*>(q); });
    const size_t shape[1] = {p->size()};
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(p->data(), 1, shape, std::move(owner));
}
template <class T>
nb::ndarray<nb::numpy, T, nb::ndim<2>> own2(std::vector<T>&& v, size_t rows, size_t cols) {
    if (v.size() != rows * cols) v.resize(rows * cols);
    auto* p = new std::vector<T>(std::move(v));
    nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<T>*>(q); });
    const size_t shape[2] = {rows, cols};
    return nb::ndarray<nb::numpy, T, nb::ndim<2>>(p->data(), 2, shape, std::move(owner));
}
// a byte vector as a numpy bool array (bool is one byte)
nb::ndarray<nb::numpy, bool, nb::ndim<1>> own_bool1(std::vector<uint8_t>&& v) {
    auto* p = new std::vector<uint8_t>(std::move(v));
    nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<uint8_t>*>(q); });
    const size_t shape[1] = {p->size()};
    return nb::ndarray<nb::numpy, bool, nb::ndim<1>>(reinterpret_cast<bool*>(p->data()), 1, shape, std::move(owner));
}
nb::ndarray<nb::numpy, bool, nb::ndim<2>> own_bool2(std::vector<uint8_t>&& v, size_t rows, size_t cols) {
    auto* p = new std::vector<uint8_t>(std::move(v));
    nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<uint8_t>*>(q); });
    const size_t shape[2] = {rows, cols};
    return nb::ndarray<nb::numpy, bool, nb::ndim<2>>(reinterpret_cast<bool*>(p->data()), 2, shape, std::move(owner));
}
inline bool* as_bool(std::vector<uint8_t>& v) { return reinterpret_cast<bool*>(v.data()); }

// ═══ THE THREAD POOL — the locus EM's, one for the module, shared by the block pool and ψ's slot pool ═══════════
// Persistent across the calls of a run and rebuilt only when the budget changes; the calls are serialised on it.
// `n_threads` is the budget: 0 is every core, as the locus EM reads it.

rigel::EStepThreadPool& pool_of(int n_threads) {
    static std::unique_ptr<rigel::EStepThreadPool> pool;
    if (!pool || pool->n_threads() != n_threads) pool = std::make_unique<rigel::EStepThreadPool>(n_threads);
    return *pool;
}
std::mutex& pool_mutex() {
    static std::mutex m;
    return m;
}
int resolve_threads(int n_threads, int64_t n_tasks) {
    int t = n_threads > 0 ? n_threads : static_cast<int>(std::thread::hardware_concurrency());
    if (t < 1) t = 1;
    if (static_cast<int64_t>(t) > n_tasks) t = static_cast<int>(std::max<int64_t>(n_tasks, 1));
    return t;
}

// ═══ THE CONSTRUCTIONS THE BLOCK PIPELINE ABSORBED ═════════════════════════════════════════════════════════════
// Each is the Python it replaces term for term and in its operation order — the factory rows and the factor precision
// of `density_deconv`, the strand evidence of `region_init` — so the block is bit-identical to the loop it replaces.

// numpy's pairwise summation of a contiguous float64 row (`pairwise_sum` in its loops: eight accumulators over blocks
// of at most 128, halves above), the summation order `np.sum(axis=1)` uses — the one reduction the factor precision
// makes, reproduced so that the precision is the same bits
constexpr int PW_BLOCKSIZE = 128;
double pairwise_sum(const double* a, int n) {
    if (n < 8) {
        double res = 0.0;
        for (int i = 0; i < n; ++i) res += a[i];
        return res;
    }
    if (n <= PW_BLOCKSIZE) {
        double r[8];
        for (int j = 0; j < 8; ++j) r[j] = a[j];
        int i = 8;
        for (; i < n - (n % 8); i += 8)
            for (int j = 0; j < 8; ++j) r[j] += a[i + j];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) res += a[i];
        return res;
    }
    int n2 = n / 2;
    n2 -= n2 % 8;
    return pairwise_sum(a, n2) + pairwise_sum(a + n2, n - n2);
}

constexpr double DENSITY_EPS = 1.0e-12;  // density_deconv._EPS
constexpr double OWN_EPS = 1.0e-9;       // region_init._EPS: the own-evidence predicate's guard, and the clips

// the fitted gDNA background the factory scores each intron against (density_deconv.GdnaBackground)
struct Background {
    double log_mu_bg, alpha, size;
    bool informative;
    // effective size: 1/α_eff = 1/α + 1/size — the per-region over-dispersion ⊕ the posterior's own width
    double alpha_eff() const {
        double inv_alpha = std::isfinite(alpha) ? 1.0 / std::max(alpha, DENSITY_EPS) : 0.0;
        inv_alpha += 1.0 / std::max(size, DENSITY_EPS);
        return inv_alpha <= DENSITY_EPS ? std::numeric_limits<double>::infinity() : 1.0 / inv_alpha;
    }
};

// The intron factory's λ-factor for one slot: log NegBinom(f_g·C; ρ_bg·E, α_eff) over the σ(λ) grid, the
// row offset so its max is 0 (an f_g-independent constant is irrelevant to ψ) — via `_log_negbinom`, the mean/size
// parameterisation, r → ∞ the exact Poisson limit
// density_deconv._log_negbinom: log NegBinom(g; mean μ, size r) for continuous g ≥ 0 in the mean/size
// parameterisation, Γln(g+r) − Γln(r) − Γln(g+1) + r·log(r/(r+μ)) + g·log(μ/(r+μ)); r → ∞ is the exact Poisson
// limit g·log μ − μ − Γln(g+1) (taken directly: no Γln(∞))
// ⛔ every product is a named temporary: numpy rounds each operation, and a product left inside the sum lets the
// compiler fuse it into a multiply-add — one rounding fewer, and a row a bit off numpy's
double log_negbinom(double g, double mu, double size) {
    mu = std::max(mu, DENSITY_EPS);
    const double lg_g1 = std::lgamma(g + 1.0);
    if (!std::isfinite(size)) {
        const double t = g * std::log(mu);
        return t - mu - lg_g1;
    }
    const double r = std::max(size, DENSITY_EPS);
    const double rpm = r + mu;
    const double t_r = r * (std::log(r) - std::log(rpm));
    const double t_g = g * (std::log(mu) - std::log(rpm));
    return std::lgamma(g + r) - std::lgamma(r) - lg_g1 + t_r + t_g;
}

void factory_row(const Background& bg, double count, double eff, const double* fg, int K, double* out) {
    if (!bg.informative) { std::fill(out, out + K, 0.0); return; }
    const double Eg = std::max(eff, DENSITY_EPS);
    const double alpha_eff = bg.alpha_eff();
    const double mu = std::exp(bg.log_mu_bg) * Eg;  // the background gDNA count location
    for (int k = 0; k < K; ++k) {
        const double g = std::clamp(fg[k], DENSITY_EPS, 1.0 - DENSITY_EPS) * count;
        out[k] = log_negbinom(g, mu, alpha_eff);
    }
    const double m = vmax(out, K);
    for (int k = 0; k < K; ++k) out[k] -= m;
}

// The factor's precision for one row: the composition evidence a λ-factor carries, read off its own
// curvature — τ = 1/Var_λ under the normalised factor; a flat row carries none (never the grid's own width)
double factor_precision_row(const double* row, const double* lam, const double* lam2, int K, std::vector<double>& w) {
    if (!(ptp(row, K) > DENSITY_EPS)) return 0.0;
    w.resize(static_cast<size_t>(K) * 3);
    double* e = w.data();
    double* wl = e + K;
    double* wl2 = wl + K;
    const double m = vmax(row, K);
    for (int k = 0; k < K; ++k) e[k] = std::exp(row[k] - m);
    const double z = std::max(pairwise_sum(e, K), DENSITY_EPS);
    for (int k = 0; k < K; ++k) e[k] /= z;
    for (int k = 0; k < K; ++k) { wl[k] = e[k] * lam[k]; wl2[k] = e[k] * lam2[k]; }
    const double mu = pairwise_sum(wl, K);
    const double mu2 = mu * mu;  // a named product: no fused multiply-add against numpy's two roundings
    const double var = pairwise_sum(wl2, K) - mu2;
    return var > DENSITY_EPS ? 1.0 / std::max(var, DENSITY_EPS) : 0.0;
}

// region_init.strand_evidence for one slot: the reference-free strand composition evidence I_strand at the
// message-free fg_loc — N_eff·disc·[f_g(1−f_g)]² / (4 p(1−p)), p = κ + f_g(½ − κ), the count overdispersed
double strand_evidence(double u_pos, double u_neg, double fg_loc, double kappa, double od_r, double disc) {
    const double n_raw = u_pos + u_neg;
    const double n_str = n_raw / (1.0 + std::max(n_raw - 1.0, 0.0) * od_r);
    const double fgl = std::clamp(fg_loc, OWN_EPS, 1.0 - OWN_EPS);
    const double pmix = std::clamp(kappa + fgl * (0.5 - kappa), OWN_EPS, 1.0 - OWN_EPS);
    const double sq = fgl * (1.0 - fgl);
    return n_str * disc * (sq * sq) / (4.0 * pmix * (1.0 - pmix));
}

// ═══ ψ FOR THE PRE-SWEEP SOLVE AND THE GATES ══════════════════════════════════════════════════════════════════

// the fitted gDNA arm's arrays as a call receives them: None, or (log_rho, logP, mass, eff)
struct ArmArrays {
    Vec log_rho, logP, mass, eff;
    P::Arm arm;
    ArmArrays(nb::object gdna, int m) {
        if (gdna.is_none()) return;
        nb::tuple t = nb::cast<nb::tuple>(gdna);
        log_rho = nb::cast<Vec>(t[0]); logP = nb::cast<Vec>(t[1]); mass = nb::cast<Vec>(t[2]); eff = nb::cast<Vec>(t[3]);
        if (log_rho.shape(0) < 1 || logP.shape(0) != log_rho.shape(0))
            throw std::invalid_argument("psi: the gDNA arm's curve is not two (G,) arrays");
        if (static_cast<int>(mass.shape(0)) != m || static_cast<int>(eff.shape(0)) != m)
            throw std::invalid_argument("psi: the gDNA arm's support is not one (mass, eff) per slot");
        arm = P::Arm{log_rho.data(), logP.data(), static_cast<int>(log_rho.shape(0)), mass.data(), eff.data()};
    }
};

// the two prior rows (the λ-factor's, the delivered composition's) as a call receives them: None or (m, K)
struct RowPriors {
    bool has_l, has_r;
    Mat l, r;
    RowPriors(nb::object lam_logprior, nb::object row_logprior, int K)
        : has_l(!lam_logprior.is_none()), has_r(!row_logprior.is_none()) {
        if (has_l) l = nb::cast<Mat>(lam_logprior);
        if (has_r) r = nb::cast<Mat>(row_logprior);
        if ((has_l && static_cast<int>(l.shape(1)) != K) || (has_r && static_cast<int>(r.shape(1)) != K))
            throw std::invalid_argument("psi: a prior row is not on the solve grid");
    }
    const double* row_of(const Mat& m, bool present, int i) const {
        return present ? m.data() + static_cast<size_t>(i) * m.shape(1) : nullptr;
    }
    P::SlotInputs inputs(const P::Arm& arm, int i, double u_pos, double u_neg, double fg_ref, double fpos_ref,
                         double fneg_ref, bool ap, bool an, const P::Delivered* row) const {
        return P::SlotInputs{u_pos, u_neg, fg_ref, fpos_ref, fneg_ref, ap, an, arm.built(),
                             arm.built() ? arm.shift(i) : 0.0, row_of(l, has_l, i), row_of(r, has_r, i), row};
    }
};

P::Rows unpack_rows(int m, int K, IdxVec& cube_slot, Mat& cube_pos, BoolVec& has_pos, Mat& cube_neg, BoolVec& has_neg,
                    Vec& cube_u, Vec& cube_total, Vec& cube_opportunity, Vec& cube_rho) {
    const int d = static_cast<int>(cube_slot.shape(0));
    if (d && (static_cast<int>(cube_pos.shape(1)) != K || static_cast<int>(cube_u.shape(0)) != K))
        throw std::invalid_argument("psi: a delivered row is not on the solve grid");
    return P::unpack_rows(m, K, d, cube_slot.data(), cube_pos.data(), has_pos.data(), cube_neg.data(), has_neg.data(),
                          cube_u.data(), cube_total.data(), cube_opportunity.data(), cube_rho.data());
}

// THE SLOT POOL: every slot is solved on its own cube with nothing shared but the read-only grid and the delivered
// rows, and writes its own four outputs, so the slot list is pulled one slot at a time by the pool — the same
// arithmetic per slot in the same order within the slot, whatever thread takes it: BIT-IDENTICAL to the serial loop
// at every thread count. One slot at a time (no chunk constant) balances the load on asymmetric cores.
void psi_solve(IdxVec slots, Vec u_pos, Vec u_neg, BoolVec allow_pos, BoolVec allow_neg, Vec fg_ref, Vec fpos_ref,
               Vec fneg_ref, double kappa, double od_g, double od_r, Vec lam, nb::object gdna, nb::object lam_logprior,
               nb::object row_logprior, IdxVec cube_slot, Mat cube_pos, BoolVec cube_has_pos, Mat cube_neg,
               BoolVec cube_has_neg, Vec cube_u, Vec cube_total, Vec cube_opportunity, Vec cube_rho, int n_tilt,
               Vec out_fg, Vec out_fpos, Vec out_fneg, Vec out_var, int n_threads) {
    const int m = static_cast<int>(u_pos.shape(0)), K = static_cast<int>(lam.shape(0));
    if (n_tilt < 2) throw std::invalid_argument("psi_solve: the tilt needs at least two nodes");
    ArmArrays A(gdna, m);
    RowPriors R(lam_logprior, row_logprior, K);
    P::Rows D = unpack_rows(m, K, cube_slot, cube_pos, cube_has_pos, cube_neg, cube_has_neg, cube_u, cube_total,
                            cube_opportunity, cube_rho);
    const P::Grid g(lam.data(), K, kappa, od_g, od_r, n_tilt, &A.arm);
    const int64_t n_sel = static_cast<int64_t>(slots.shape(0));
    const int64_t* sl = slots.data();
    const double *up = u_pos.data(), *un = u_neg.data(), *fr = fg_ref.data(), *pr = fpos_ref.data(), *nr = fneg_ref.data();
    const bool *ap = allow_pos.data(), *an = allow_neg.data();
    double *o_fg = out_fg.data(), *o_fp = out_fpos.data(), *o_fn = out_fneg.data(), *o_v = out_var.data();
    auto solve_one = [&](int64_t q, P::Scratch& S) {
        const int i = static_cast<int>(sl[q]);
        const P::SlotInputs s = R.inputs(A.arm, i, up[i], un[i], fr[i], pr[i], nr[i], ap[i], an[i],
                                         D.index[i] >= 0 ? &D.rows[D.index[i]] : nullptr);
        P::solve_slot(g, s, S, o_fg[i], o_fp[i], o_fn[i], o_v[i]);
    };
    const int T = resolve_threads(n_threads, n_sel);
    if (T == 1) {
        P::Scratch S;
        for (int64_t q = 0; q < n_sel; ++q) solve_one(q, S);
        return;
    }
    std::vector<P::Scratch> scratch(T);
    std::atomic<int64_t> next{0};
    auto worker = [&](int tid) {
        P::Scratch& S = scratch[tid];
        for (int64_t q; (q = next.fetch_add(1, std::memory_order_relaxed)) < n_sel;) solve_one(q, S);
    };
    nb::gil_scoped_release release;
    std::lock_guard<std::mutex> lock(pool_mutex());
    pool_of(T).run_parallel(worker);
}

// ψ itself, for the gates: the cube (m, K, C) with the two strand-fraction grids and the tilt beside it,
// every slot of one class (ambig: n_tilt + 2 columns; else one)
void psi_cube(Vec u_pos, Vec u_neg, BoolVec allow_pos, BoolVec allow_neg, Vec fg_ref, Vec fpos_ref, Vec fneg_ref,
              double kappa, double od_g, double od_r, Vec lam, nb::object gdna, nb::object lam_logprior,
              nb::object row_logprior, IdxVec cube_slot, Mat cube_pos, BoolVec cube_has_pos, Mat cube_neg,
              BoolVec cube_has_neg, Vec cube_u, Vec cube_total, Vec cube_opportunity, Vec cube_rho, int n_tilt,
              bool ambig, Cube out_psi, Cube out_fpos, Cube out_fneg, Cube out_tau) {
    const int m = static_cast<int>(u_pos.shape(0)), K = static_cast<int>(lam.shape(0));
    ArmArrays A(gdna, m);
    RowPriors R(lam_logprior, row_logprior, K);
    const P::Grid g(lam.data(), K, kappa, od_g, od_r, n_tilt, &A.arm);
    const int C = g.columns(ambig);
    if (static_cast<int>(out_psi.shape(1)) != K || static_cast<int>(out_psi.shape(2)) != C)
        throw std::invalid_argument("psi_cube: the output is not (m, K, columns)");
    P::Rows D = unpack_rows(m, K, cube_slot, cube_pos, cube_has_pos, cube_neg, cube_has_neg, cube_u, cube_total,
                            cube_opportunity, cube_rho);
    std::vector<double> scratch;
    for (int i = 0; i < m; ++i) {
        const P::SlotInputs s = R.inputs(A.arm, i, u_pos.data()[i], u_neg.data()[i], fg_ref.data()[i], fpos_ref.data()[i],
                                         fneg_ref.data()[i], allow_pos.data()[i], allow_neg.data()[i],
                                         D.index[i] >= 0 ? &D.rows[D.index[i]] : nullptr);
        const size_t off = static_cast<size_t>(i) * K * C;
        P::slot_cube(g, s, ambig, out_psi.data() + off, out_fpos.data() + off, out_fneg.data() + off, out_tau.data() + off,
                     scratch);
    }
}

void posterior_median_rows(Mat post, Vec lam, Vec out) {
    const int m = static_cast<int>(post.shape(0)), K = static_cast<int>(post.shape(1));
    std::vector<double> cdf;
    for (int i = 0; i < m; ++i)
        out.data()[i] = P::posterior_median(post.data() + static_cast<size_t>(i) * K, lam.data(), K, cdf);
}

// the fitted gDNA arm as an (m, K) matrix, for the gates: the kernel's own construction at every (slot, cell)
void gdna_arm_rows(Vec log_rho, Vec logP, Vec lam, Vec mass, Vec eff, Mat out) {
    const int m = static_cast<int>(mass.shape(0)), K = static_cast<int>(lam.shape(0));
    if (static_cast<int>(out.shape(0)) != m || static_cast<int>(out.shape(1)) != K)
        throw std::invalid_argument("gdna_arm: the output is not (m, K)");
    if (logP.shape(0) != log_rho.shape(0) || log_rho.shape(0) < 1 || static_cast<int>(eff.shape(0)) != m)
        throw std::invalid_argument("gdna_arm: the curve is two (G,) arrays and the support one (mass, eff) per slot");
    const P::Arm arm{log_rho.data(), logP.data(), static_cast<int>(log_rho.shape(0)), mass.data(), eff.data()};
    std::vector<double> log_frac(K);
    for (int k = 0; k < K; ++k) log_frac[k] = P::arm_abscissa(sigmoid(lam.data()[k]));
    for (int i = 0; i < m; ++i) {
        const double sh = arm.shift(i);
        double* o = out.data() + static_cast<size_t>(i) * K;
        for (int k = 0; k < K; ++k) o[k] = arm.at(log_frac[k], sh);
    }
}

void compose_rows(Vec f_g, Vec w_pos, BoolVec allow_pos, BoolVec allow_neg, Vec out_fpos, Vec out_fneg) {
    const int m = static_cast<int>(f_g.shape(0));
    for (int i = 0; i < m; ++i)
        P::compose(f_g.data()[i], w_pos.data()[i], allow_pos.data()[i], allow_neg.data()[i], out_fpos.data()[i],
                   out_fneg.data()[i]);
}

// ═══ THE BLOCK PIPELINE ═══════════════════════════════════════════════════════════════════════════════════════

enum Policy { SILENT = 0, TRANSFER = 1 };

// the backbone's assertions, counted per block on the owned slots as (violations, eligible); the backbone raises
enum Assertion { POP_AT_MOST_THREE, POP_REACHES_THREE, LAM_ROWS_FINITE, CUBE_ROWS_FINITE, WRITEBACK_ONLY_SOLVABLE, N_ASSERT };
const char* const ASSERTION_NAMES[N_ASSERT] = {"population_at_most_three", "population_reaches_three", "lam_rows_finite",
                                               "cube_rows_finite", "writeback_only_solvable"};

struct ChainArrays {  // the whole chain, read-only
    int64_t n;
    const bool *is_bnd, *is_exon, *fp, *fn, *exon_pos, *exon_neg, *terminal;
    const int64_t *left, *right;
    const uint16_t* flags;
    const double *cnt, *spl, *sj, *sj_lo, *sj_hi, *route_lo, *route_hi;  // (n, 2)
    const double *a_g, *a_r;
    const double *bel_fpos, *bel_fneg, *bel_fg;  // the incoming belief
};

struct Params {
    double kappa, od_g, od_r, disc;  // ψ's strand model; the protocol's discriminability (0: the channel is dead)
    const double* lam; int K; int n_tilt;
    // the policy: its kind, and the strand model ITS own claims read (a policy's, not the sweep's — the cache gates
    // hold the two apart), the library's coordinates and strand witness
    int policy; bool has_strand; double pol_kappa, pol_od_g, pol_od_r; double rho_gdna, rho_rna; bool split_live;
};

struct Landscape { bool built = false; const double* log_rho = nullptr; const double* logP = nullptr; int G = 0; };

struct Factory {  // the intron factory: its inputs (the rows built per block), rows given as one array, or none
    enum { NONE = 0, INPUTS = 1, ROWS = 2 } mode = NONE;
    const bool* is_intron = nullptr; const double* count = nullptr; const double* eff = nullptr; Background bg{};
    const double* rows = nullptr;  // (n, K) chain-wide, ROWS only
};

struct Outputs {  // written in place, the owned slots of each block
    double *fpos, *fneg, *fg, *var; bool* has_comp;
    // the diagnostics capture, or nullptr
    double *fg_loc = nullptr, *fg_strand = nullptr, *tau_lam = nullptr, *tau_fac = nullptr, *lam_rows = nullptr;
};

struct Lane3 {  // one received lane as the capture returns it
    std::vector<uint8_t> present, has_witness; std::vector<double> profile, count, opp, rna_count, rna_var;
};
struct RecvCopy { std::vector<uint8_t> has_nbr, has_comp; std::vector<double> comp; Lane3 lanes[3]; };

struct BlockResult {
    int64_t counts[N_ASSERT][2] = {};
    bool rows_delivered = false, cube_delivered = false, captured = false;
    // the capture: the cube rows the solve delivered (in the block's local slots) and the two received tables on
    // the owned slots
    std::vector<int64_t> cslot; std::vector<double> cpos, cneg, ctotal, copp, crho; std::vector<uint8_t> chas_pos, chas_neg;
    RecvCopy recv[2];
};

// one thread's tables for one block, resized per block, the capacity kept across blocks
struct Arena {
    int n = 0, K = 0;
    std::vector<double> n_u, n_s, flux, cnt_col[2], factory, own_rows, lane_own[3], comp[2], prof[2][3], rows, cube_pos,
        cube_neg, cube_total, cube_opp, cube_rho, f_par[6], lane_wit[3], face_store, flux_store[3], recv_count[2][3],
        recv_opp[2][3], recv_rna[2][3], recv_rnav[2][3], fg_loc, fp_loc, fn_loc, tau_lam, tau_fac, out, bound, row, w,
        lam2;
    std::vector<char> is_intron, is_intergenic;
    std::vector<int64_t> left, right, seq_f, seq_b, cube_slot;
    std::vector<int8_t> f_kind;
    std::vector<int32_t> f_row, f_row2, flux_index[3];
    std::vector<uint8_t> own_mask, lane_face[3], lane_two[3], lane_own_mask[3], lane_wit_mask[3], empty_g, empty_r,
        has_nbr[2], has_comp[2], present[2][3], has_wit[2][3], written, cube_has_pos, cube_has_neg, held, solvable,
        has_own, own_ev;
    std::vector<const double*> factory_ptr, row_ptr;
    std::unique_ptr<Scratch> S;  // transfer_rows' scratch, built for K
    P::Scratch psi;

    void size(int n_, int K_) {
        if (K_ != K || !S) S = std::make_unique<Scratch>(K_);
        n = n_; K = K_;
        const size_t nk = static_cast<size_t>(n) * K;
        auto grow = [](auto& v, size_t m) { if (v.size() < m) v.resize(m); };
        for (auto* v : {&n_u, &n_s, &flux, &cnt_col[0], &cnt_col[1], &fg_loc, &fp_loc, &fn_loc, &tau_lam, &tau_fac,
                        &cube_total, &cube_opp, &cube_rho})
            grow(*v, n);
        for (auto* v : {&factory, &own_rows, &lane_own[0], &lane_own[1], &lane_own[2], &comp[0], &comp[1], &rows,
                        &cube_pos, &cube_neg})
            grow(*v, nk);
        for (int s = 0; s < 2; ++s)
            for (int l = 0; l < 3; ++l) {
                grow(prof[s][l], nk);
                grow(recv_count[s][l], n); grow(recv_opp[s][l], n); grow(recv_rna[s][l], n); grow(recv_rnav[s][l], n);
                grow(present[s][l], n); grow(has_wit[s][l], n);
            }
        for (int l = 0; l < 3; ++l) {
            grow(lane_wit[l], 2 * n); grow(lane_face[l], 2 * n); grow(lane_two[l], 2 * n); grow(lane_own_mask[l], n);
            grow(lane_wit_mask[l], n); grow(flux_index[l], 2 * n);
        }
        for (int t = 0; t < 6; ++t) grow(f_par[t], 2 * n);
        grow(f_kind, 2 * n); grow(f_row, 2 * n); grow(f_row2, 2 * n);
        grow(is_intron, n); grow(is_intergenic, n);
        grow(left, n); grow(right, n); grow(seq_f, n); grow(seq_b, n); grow(cube_slot, n);
        for (auto* v : {&own_mask, &empty_g, &empty_r, &has_nbr[0], &has_nbr[1], &has_comp[0], &has_comp[1], &written,
                        &cube_has_pos, &cube_has_neg, &held, &solvable, &has_own, &own_ev})
            grow(*v, n);
        grow(factory_ptr, n); grow(row_ptr, n);
        grow(lam2, K); grow(w, 3 * static_cast<size_t>(K));
    }
};

// one locus block, end to end, on one thread
void solve_one_block(const ChainArrays& C, const Params& p, const P::Grid& g, const Landscape& L, const Factory& F,
                     int64_t start, int64_t stop, int64_t end, bool want_capture, Outputs& out, Arena& A,
                     BlockResult& res) {
    const int n = static_cast<int>(end - start), n_owned = static_cast<int>(stop - start), K = p.K;
    A.size(n, K);
    Scratch& S = *A.S;
    for (int k = 0; k < K; ++k) A.lam2[k] = p.lam[k] * p.lam[k];
    // the block's slices of the chain, its links re-based (a neighbour outside the block is no neighbour)
    const bool *is_bnd = C.is_bnd + start, *is_exon = C.is_exon + start, *fp = C.fp + start, *fn = C.fn + start;
    const bool *exon_pos = C.exon_pos + start, *exon_neg = C.exon_neg + start, *term = C.terminal + start;
    const uint16_t* flags = C.flags + start;
    const double *cnt = C.cnt + 2 * start, *spl = C.spl + 2 * start, *sj = C.sj + 2 * start;
    const double *sj_lo = C.sj_lo + 2 * start, *sj_hi = C.sj_hi + 2 * start, *route_lo = C.route_lo + 2 * start;
    const double *route_hi = C.route_hi + 2 * start, *a_g = C.a_g + start, *a_r = C.a_r + start;
    const double *bel_fpos = C.bel_fpos + start, *bel_fneg = C.bel_fneg + start, *bel_fg = C.bel_fg + start;
    for (int i = 0; i < n; ++i) {
        const int64_t l = C.left[start + i] - start, r = C.right[start + i] - start;
        A.left[i] = (l >= 0 && l < n) ? l : -1;
        A.right[i] = (r >= 0 && r < n) ? r : -1;
        A.n_u[i] = cnt[2 * i] + cnt[2 * i + 1];
        A.n_s[i] = spl[2 * i] + spl[2 * i + 1];
        A.flux[i] = sj[2 * i] + sj[2 * i + 1];
        A.cnt_col[0][i] = cnt[2 * i]; A.cnt_col[1][i] = cnt[2 * i + 1];
        A.seq_f[i] = i; A.seq_b[i] = n - 1 - i;
    }
    Chain::classify(n, is_bnd, is_exon, fp, fn, A.is_intron.data(), A.is_intergenic.data());
    // THE PRIOR ROWS: the factory's per intron (its inputs, or the rows given), the arm's shift per slot
    const double* fg = g.fg.data();
    for (int i = 0; i < n; ++i) {
        const double* r = nullptr;
        if (F.mode == Factory::INPUTS) {
            if (F.is_intron[start + i]) {
                double* dest = A.factory.data() + static_cast<size_t>(i) * K;
                factory_row(F.bg, F.count[start + i], F.eff[start + i], fg, K, dest);
                r = dest;
            }
        } else if (F.mode == Factory::ROWS) {
            r = F.rows + static_cast<size_t>(start + i) * K;
        }
        A.factory_ptr[i] = r;
    }
    // THE SELF-SOLVE: every slot with a live strand and a fragment, on its own cube, the incoming belief the freeze
    auto in_slots = [&](int i) {
        const double signal = cnt[2 * i] + cnt[2 * i + 1] + A.n_u[i] + A.n_s[i];
        return (fp[i] || fn[i]) && signal > 0.0;
    };
    for (int i = 0; i < n; ++i) {
        as_bool(A.solvable)[i] = (fp[i] || fn[i]) && A.n_u[i] > 0.0;
        A.fg_loc[i] = bel_fg[i];
        if (!in_slots(i)) continue;
        const P::SlotInputs s{cnt[2 * i], cnt[2 * i + 1], bel_fg[i], bel_fpos[i], bel_fneg[i], fp[i], fn[i], L.built,
                              L.built ? P::Arm::shift_of(A.n_u[i], a_g[i]) : 0.0, A.factory_ptr[i], nullptr, nullptr};
        double f_g, f_p, f_n, v_g;
        P::solve_slot(g, s, A.psi, f_g, f_p, f_n, v_g);
        if (as_bool(A.solvable)[i]) A.fg_loc[i] = f_g;
    }
    // THE OWN EVIDENCE: the single-strand strand precision at fg_loc + the factory row's curvature
    for (int i = 0; i < n; ++i) {
        const bool single = fp[i] != fn[i];
        const double i_strand = strand_evidence(cnt[2 * i], cnt[2 * i + 1], A.fg_loc[i], p.kappa, p.od_r, p.disc);
        double tau = single ? i_strand : 0.0;
        double fac = 0.0;
        if (F.mode != Factory::NONE) {
            fac = A.factory_ptr[i] ? factor_precision_row(A.factory_ptr[i], p.lam, A.lam2.data(), K, A.w) : 0.0;
            tau = tau + fac;
        }
        A.tau_lam[i] = tau; A.tau_fac[i] = fac;
        as_bool(A.has_own)[i] = tau > 0.0;
        as_bool(A.own_ev)[i] = tau > OWN_EPS;
    }
    // THE LAYER
    bool has_rows = false;
    int d = 0;
    P::Rows delivered;
    std::fill_n(as_bool(A.held), n, false);
    const bool run_layer = p.policy == TRANSFER;
    if (run_layer) {
        Chain c;
        c.n = n; c.K = K; c.lam = p.lam;
        c.is_bnd = is_bnd; c.is_exon = is_exon; c.fp = fp; c.fn = fn; c.exon_pos = exon_pos; c.exon_neg = exon_neg;
        c.has_own = as_bool(A.has_own); c.left = A.left.data(); c.right = A.right.data(); c.flags = flags;
        c.n_u = A.n_u.data(); c.n_s = A.n_s.data(); c.a_g = a_g; c.a_r = a_r; c.belief = bel_fg; c.cnt = cnt;
        c.route_lo = route_lo; c.route_hi = route_hi; c.sj_lo = sj_lo; c.sj_hi = sj_hi; c.flux = A.flux.data();
        c.src = A.factory_ptr.data(); c.is_intron = A.is_intron.data(); c.is_intergenic = A.is_intergenic.data();
        c.has_strand = p.has_strand; c.kappa = p.pol_kappa; c.od_g = p.pol_od_g; c.od_r = p.pol_od_r;
        // the tables, their bits cleared, their matrices unfilled
        std::fill_n(as_bool(A.own_mask), n, false);
        RowsOut own{A.own_rows.data(), as_bool(A.own_mask), K};
        std::fill_n(A.f_kind.data(), 2 * n, static_cast<int8_t>(NONE));
        FacesOut Fc{&c, A.f_kind.data(), A.f_row.data(), A.f_row2.data(), A.f_par[0].data(), A.f_par[1].data(),
                    A.f_par[2].data(), A.f_par[3].data(), A.f_par[4].data(), A.f_par[5].data(), RowStore{&A.face_store, K}};
        LaneOut lanes_out[3];
        for (int l = 0; l < 3; ++l) {
            std::fill_n(as_bool(A.lane_own_mask[l]), n, false);
            std::fill_n(as_bool(A.lane_wit_mask[l]), n, false);
            std::fill_n(A.flux_index[l].data(), 2 * n, -1);
            lanes_out[l] = LaneOut{as_bool(A.lane_face[l]), as_bool(A.lane_two[l]),
                                   RowsOut{A.lane_own[l].data(), as_bool(A.lane_own_mask[l]), K},
                                   RowStore{&A.flux_store[l], K}, A.flux_index[l].data(), A.lane_wit[l].data(),
                                   as_bool(A.lane_wit_mask[l])};
        }
        const bool gdna_built = p.rho_gdna > 0.0;
        prepare_block(c, own, Fc, gdna_built ? &lanes_out[0] : nullptr, lanes_out[1], lanes_out[2], p.rho_gdna,
                      p.rho_rna, S);
        // the lanes as the pass reads them; the RNA lanes' witness columns under the protocol
        for (int i = 0; i < n; ++i) {
            as_bool(A.empty_g)[i] = !(A.n_u[i] > 0.0) || !(a_g[i] > 0.0);
            as_bool(A.empty_r)[i] = !(A.n_u[i] > 0.0) || !(a_r[i] > 0.0);
        }
        std::vector<LaneView> lanes;
        if (gdna_built)
            lanes.push_back(LaneView{0, as_bool(A.lane_face[0]), as_bool(A.lane_two[0]), as_bool(A.empty_g),
                                     A.lane_own[0].data(), as_bool(A.lane_own_mask[0]), A.n_u.data(), a_g, nullptr,
                                     A.lane_wit[0].data(), as_bool(A.lane_wit_mask[0])});
        for (int s = 0; s < 2; ++s) {
            const int col_read = read_column(s, p.has_strand, p.pol_kappa);
            lanes.push_back(LaneView{1 + s, as_bool(A.lane_face[1 + s]), as_bool(A.lane_two[1 + s]), as_bool(A.empty_r),
                                     A.lane_own[1 + s].data(), as_bool(A.lane_own_mask[1 + s]), A.cnt_col[col_read].data(),
                                     a_r, p.split_live ? A.cnt_col[1 - col_read].data() : nullptr, A.lane_wit[1 + s].data(),
                                     as_bool(A.lane_wit_mask[1 + s])});
        }
        const FacesView Fv{A.f_kind.data(), A.f_row.data(), A.f_row2.data(), A.f_par[0].data(), A.f_par[1].data(),
                           A.f_par[2].data(), A.f_par[3].data(), A.f_par[4].data(), A.f_par[5].data(),
                           A.face_store.data(), K};
        // THE TWO PASSES on the two received tables — the forward pass reads each node's low neighbour, the backward
        // its high one; the backbone's bit `has_neighbour` says the side exists
        ReceivedView R[2];
        for (int s = 0; s < 2; ++s) {
            std::fill_n(as_bool(A.has_comp[s]), n, false);
            R[s].has_neighbour = as_bool(A.has_nbr[s]); R[s].has_composition = as_bool(A.has_comp[s]);
            R[s].composition = A.comp[s].data(); R[s].K = K;
            for (int l = 0; l < 3; ++l) {
                std::fill_n(as_bool(A.present[s][l]), n, false);
                R[s].lanes[l] = LevelsView{as_bool(A.present[s][l]), A.prof[s][l].data(), A.recv_count[s][l].data(),
                                           A.recv_opp[s][l].data(), as_bool(A.has_wit[s][l]), A.recv_rna[s][l].data(),
                                           A.recv_rnav[s][l].data(), K};
            }
            const int64_t* nbr = s == 0 ? A.left.data() : A.right.data();
            for (int i = 0; i < n; ++i) R[s].has_neighbour[i] = nbr[i] >= 0;
            pass_block(p.lam, K, (s == 0 ? A.seq_f : A.seq_b).data(), n, nbr, term, A.own_rows.data(), as_bool(A.own_mask),
                       Fv, lanes, R[s], S, A.out);
        }
        // THE SOLVE: the two held tables into ψ's two channels
        SolveLane G, rna[2];
        if (gdna_built) {
            G.built = true; G.field = 0; G.empty = as_bool(A.empty_g); G.total = A.n_u.data(); G.a = a_g;
            G.rho_ref = p.rho_gdna; G.own_rows = A.lane_own[0].data(); G.own_mask = as_bool(A.lane_own_mask[0]);
            G.flux_store = A.flux_store[0].data(); G.flux_index = A.flux_index[0].data();
        }
        for (int s = 0; s < 2; ++s) {
            rna[s].built = true; rna[s].field = 1 + s; rna[s].empty = as_bool(A.empty_r); rna[s].total = A.n_u.data();
            rna[s].a = a_r; rna[s].rho_ref = p.rho_rna; rna[s].own_rows = A.lane_own[1 + s].data();
            rna[s].own_mask = as_bool(A.lane_own_mask[1 + s]); rna[s].flux_store = A.flux_store[1 + s].data();
            rna[s].flux_index = A.flux_index[1 + s].data();
        }
        std::fill_n(A.rows.data(), static_cast<size_t>(n) * K, 0.0);
        CubeOut cube{A.cube_slot.data(), A.cube_pos.data(), as_bool(A.cube_has_pos), A.cube_neg.data(),
                     as_bool(A.cube_has_neg), A.cube_total.data(), A.cube_opp.data(), A.cube_rho.data(), n};
        const auto [live, n_cube] = solve_block(p.lam, K, n, R, G, rna, fp, fn, A.rows.data(), as_bool(A.written), cube,
                                                S, A.bound, A.row);
        has_rows = live; d = n_cube;
        for (int i = 0; i < n; ++i) as_bool(A.held)[i] = R[0].has_composition[i] || R[1].has_composition[i];
        if (d)
            delivered = P::unpack_rows(n, K, d, A.cube_slot.data(), A.cube_pos.data(), as_bool(A.cube_has_pos),
                                       A.cube_neg.data(), as_bool(A.cube_has_neg), p.lam, A.cube_total.data(),
                                       A.cube_opp.data(), A.cube_rho.data());
        if (want_capture) {
            res.captured = true;
            for (int q = 0; q < d; ++q) {
                res.cslot.push_back(A.cube_slot[q]);
                const double *cp = A.cube_pos.data() + static_cast<size_t>(q) * K, *cn = A.cube_neg.data() + static_cast<size_t>(q) * K;
                res.cpos.insert(res.cpos.end(), cp, cp + K); res.cneg.insert(res.cneg.end(), cn, cn + K);
                res.chas_pos.push_back(as_bool(A.cube_has_pos)[q]); res.chas_neg.push_back(as_bool(A.cube_has_neg)[q]);
                res.ctotal.push_back(A.cube_total[q]); res.copp.push_back(A.cube_opp[q]); res.crho.push_back(A.cube_rho[q]);
            }
            for (int s = 0; s < 2; ++s) {
                RecvCopy& rc = res.recv[s];
                rc.has_nbr.assign(A.has_nbr[s].begin(), A.has_nbr[s].begin() + n_owned);
                rc.has_comp.assign(A.has_comp[s].begin(), A.has_comp[s].begin() + n_owned);
                rc.comp.assign(A.comp[s].begin(), A.comp[s].begin() + static_cast<size_t>(n_owned) * K);
                for (int l = 0; l < 3; ++l) {
                    Lane3& ln = rc.lanes[l];
                    ln.present.assign(A.present[s][l].begin(), A.present[s][l].begin() + n_owned);
                    ln.has_witness.assign(A.has_wit[s][l].begin(), A.has_wit[s][l].begin() + n_owned);
                    ln.profile.assign(A.prof[s][l].begin(), A.prof[s][l].begin() + static_cast<size_t>(n_owned) * K);
                    ln.count.assign(A.recv_count[s][l].begin(), A.recv_count[s][l].begin() + n_owned);
                    ln.opp.assign(A.recv_opp[s][l].begin(), A.recv_opp[s][l].begin() + n_owned);
                    ln.rna_count.assign(A.recv_rna[s][l].begin(), A.recv_rna[s][l].begin() + n_owned);
                    ln.rna_var.assign(A.recv_rnav[s][l].begin(), A.recv_rnav[s][l].begin() + n_owned);
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) A.row_ptr[i] = has_rows ? A.rows.data() + static_cast<size_t>(i) * K : nullptr;
    // THE CHECKS on what was delivered, over the owned slots
    for (int i = 0; i < n_owned; ++i) {
        const int pop = 1 + (fp[i] ? 1 : 0) + (fn[i] ? 1 : 0);
        res.counts[POP_AT_MOST_THREE][0] += pop > 3;
        res.counts[POP_AT_MOST_THREE][1] += 1;
        res.counts[POP_REACHES_THREE][1] += pop >= 3;
    }
    res.rows_delivered = has_rows;
    if (has_rows) {
        for (int i = 0; i < n_owned; ++i) {
            const double* r = A.rows.data() + static_cast<size_t>(i) * K;
            bool finite = true;
            for (int k = 0; k < K; ++k) finite = finite && std::isfinite(r[k]);
            res.counts[LAM_ROWS_FINITE][0] += !finite;
        }
        res.counts[LAM_ROWS_FINITE][1] += n_owned;
    }
    res.cube_delivered = d > 0;
    if (d > 0) {
        for (int q = 0; q < d; ++q) {
            const P::Delivered& row = delivered.rows[q];
            const int64_t i = A.cube_slot[q];
            if (i >= n_owned) continue;
            bool bad = false;
            for (const double* prof : {row.pos, row.neg}) {
                if (prof == nullptr) continue;
                for (int k = 0; k < K; ++k) bad = bad || !std::isfinite(prof[k]);
            }
            res.counts[CUBE_ROWS_FINITE][0] += bad;
            res.counts[CUBE_ROWS_FINITE][1] += 1;
        }
    }
    // THE FINAL SOLVE: the arm, the factory row and the delivered row per cell, the cube rows at the AMBIG slots
    for (int i = 0; i < n_owned; ++i) {
        const bool solvable = as_bool(A.solvable)[i];
        res.counts[WRITEBACK_ONLY_SOLVABLE][1] += !solvable;
        if (want_capture) {
            out.fg_loc[start + i] = A.fg_loc[i]; out.tau_lam[start + i] = A.tau_lam[i]; out.tau_fac[start + i] = A.tau_fac[i];
            out.fg_strand[start + i] = 0.0;
            if (out.lam_rows) {
                double* dest = out.lam_rows + static_cast<size_t>(start + i) * K;
                if (has_rows) std::copy(A.rows.data() + static_cast<size_t>(i) * K, A.rows.data() + static_cast<size_t>(i + 1) * K, dest);
                else std::fill(dest, dest + K, 0.0);
            }
        }
        if (!in_slots(i)) continue;
        const P::SlotInputs s{cnt[2 * i], cnt[2 * i + 1], bel_fg[i], bel_fpos[i], bel_fneg[i], fp[i], fn[i], L.built,
                              L.built ? P::Arm::shift_of(A.n_u[i], a_g[i]) : 0.0, A.factory_ptr[i], A.row_ptr[i],
                              delivered.index.empty() || delivered.index[i] < 0 ? nullptr : &delivered.rows[delivered.index[i]]};
        double f_g, f_p, f_n, v_g;
        P::solve_slot(g, s, A.psi, f_g, f_p, f_n, v_g);
        if (solvable) {  // THE WRITE-BACK: only a solvable slot; a locked or empty one keeps the incoming belief
            out.fpos[start + i] = std::clamp(f_p, 0.0, 1.0);
            out.fneg[start + i] = std::clamp(f_n, 0.0, 1.0);
            out.fg[start + i] = std::clamp(f_g, 0.0, 1.0);
            out.var[start + i] = v_g;
        }
        if (want_capture) {  // the strand-only solve: no prior, no messages
            const P::SlotInputs bare{cnt[2 * i], cnt[2 * i + 1], bel_fg[i], bel_fpos[i], bel_fneg[i], fp[i], fn[i], false, 0.0,
                                     nullptr, nullptr, nullptr};
            double b_g, b_p, b_n, b_v;
            P::solve_slot(g, bare, A.psi, b_g, b_p, b_n, b_v);
            out.fg_strand[start + i] = b_g;
        }
    }
    // has_composition: an own composition channel, structural certainty, or a composition held from either side
    for (int i = 0; i < n_owned; ++i)
        out.has_comp[start + i] = as_bool(A.own_ev)[i] || (!fp[i] && !fn[i]) || as_bool(A.held)[i];
}

// ---- the call ---------------------------------------------------------------------------------------------

// an optional array argument: its element pointer, the array kept alive for the call; nullptr for None
template <class A> const typename A::Scalar* opt(nb::handle h, std::vector<nb::object>& keep) {
    if (h.is_none()) return nullptr;
    keep.push_back(nb::borrow<nb::object>(h));
    return nb::cast<A>(h).data();
}

nb::dict lane_capture(Lane3&& ln, int n_owned, int K) {
    nb::dict d;
    d["present"] = own_bool1(std::move(ln.present));
    d["profile"] = own2(std::move(ln.profile), n_owned, K);
    d["count"] = own1(std::move(ln.count));
    d["opportunity"] = own1(std::move(ln.opp));
    d["has_witness"] = own_bool1(std::move(ln.has_witness));
    d["rna_count"] = own1(std::move(ln.rna_count));
    d["rna_count_var"] = own1(std::move(ln.rna_var));
    return d;
}

// THE SWEEP'S CALL. `blocks` is (B, 3) int64: start, stop, end per block. `gdna` is None or (log_rho, logP); `factory`
// None, ("inputs", is_intron, count, eff, log_mu_bg, alpha, size, informative) or ("rows", rows). The belief arrays are
// written in place on the owned slots (they arrive as copies of the incoming belief); `diagnostics` is None or a dict of
// the capture's arrays to fill. Returns a dict: `counts` (B, 5, 2) int64 per assertion (violations, eligible),
// `rows_delivered` / `cube_delivered` (B,) bool, and under a capture `cubes` — a list per block: None, or the cube rows the
// solve delivered as (slot, pos, has_pos, neg, has_neg, total, opportunity, rho) in the block's local slots — and
// `received`, a list per block (None, or the two received tables as dicts).
nb::dict solve_blocks(BoolVec is_boundary, BoolVec is_exon, BoolVec free_pos, BoolVec free_neg, BoolVec exon_pos, BoolVec exon_neg,
                      BoolVec terminal, IdxVec left, IdxVec right, FlagVec flags, Mat cnt, Mat spliced, Mat sj_count,
                      Mat sj_count_lo, Mat sj_count_hi, Mat route_rate_lo, Mat route_rate_hi, Vec eff_gdna, Vec eff_rna,
                      Vec belief_fpos, Vec belief_fneg, Vec belief_fg, nb::ndarray<int64_t, nb::ndim<2>, nb::c_contig> blocks,
                      double kappa, double od_g, double od_r, double disc, Vec lam, int n_tilt, int policy, bool has_strand,
                      double policy_kappa, double policy_od_g, double policy_od_r, double rho_gdna, double rho_rna,
                      bool split_live, nb::object gdna, nb::object factory,
                      Vec out_fpos, Vec out_fneg, Vec out_fg, Vec out_var, BoolVec out_has_composition,
                      nb::object diagnostics, int n_threads) {
    const int64_t n = static_cast<int64_t>(free_pos.shape(0));
    const int K = static_cast<int>(lam.shape(0)), B = static_cast<int>(blocks.shape(0));
    if (n_tilt < 2) throw std::invalid_argument("solve_blocks: the tilt needs at least two nodes");
    if (static_cast<int>(cnt.shape(1)) != 2 || static_cast<int64_t>(cnt.shape(0)) != n)
        throw std::invalid_argument("solve_blocks: the counts are (n, 2)");
    std::vector<nb::object> keep;  // the optional arrays, alive for the call
    ChainArrays C{n, is_boundary.data(), is_exon.data(), free_pos.data(), free_neg.data(), exon_pos.data(), exon_neg.data(),
                  terminal.data(), left.data(), right.data(), flags.data(), cnt.data(), spliced.data(), sj_count.data(),
                  sj_count_lo.data(), sj_count_hi.data(), route_rate_lo.data(), route_rate_hi.data(), eff_gdna.data(),
                  eff_rna.data(), belief_fpos.data(), belief_fneg.data(), belief_fg.data()};
    Params p{kappa, od_g, od_r, disc, lam.data(), K, n_tilt, policy, has_strand, policy_kappa, policy_od_g, policy_od_r,
             rho_gdna, rho_rna, split_live};
    Landscape L;
    if (!gdna.is_none()) {
        nb::tuple t = nb::cast<nb::tuple>(gdna);
        L.log_rho = opt<Vec>(t[0], keep); L.logP = opt<Vec>(t[1], keep);
        L.G = static_cast<int>(nb::cast<Vec>(t[0]).shape(0));
        if (L.G < 1 || static_cast<int>(nb::cast<Vec>(t[1]).shape(0)) != L.G)
            throw std::invalid_argument("solve_blocks: the gDNA arm's curve is two (G,) arrays");
        L.built = true;
    }
    Factory F;
    if (!factory.is_none()) {
        nb::tuple t = nb::cast<nb::tuple>(factory);
        const std::string mode = nb::cast<std::string>(t[0]);
        if (mode == "inputs") {
            F.mode = Factory::INPUTS;
            F.is_intron = opt<BoolVec>(t[1], keep); F.count = opt<Vec>(t[2], keep); F.eff = opt<Vec>(t[3], keep);
            F.bg = Background{nb::cast<double>(t[4]), nb::cast<double>(t[5]), nb::cast<double>(t[6]), nb::cast<bool>(t[7])};
        } else if (mode == "rows") {
            F.mode = Factory::ROWS;
            Mat rows = nb::cast<Mat>(t[1]);
            if (static_cast<int64_t>(rows.shape(0)) != n || static_cast<int>(rows.shape(1)) != K)
                throw std::invalid_argument("solve_blocks: the factory rows are (n, K)");
            F.rows = opt<Mat>(t[1], keep);
        } else {
            throw std::invalid_argument("solve_blocks: the factory is None, ('inputs', ...) or ('rows', rows)");
        }
    }
    Outputs out{out_fpos.data(), out_fneg.data(), out_fg.data(), out_var.data(), out_has_composition.data()};
    const bool want_capture = !diagnostics.is_none();
    if (want_capture) {
        nb::dict cd = nb::cast<nb::dict>(diagnostics);
        out.fg_loc = nb::cast<Vec>(cd["fg_loc"]).data(); out.fg_strand = nb::cast<Vec>(cd["fg_strand"]).data();
        out.tau_lam = nb::cast<Vec>(cd["tau_lam"]).data(); out.tau_fac = nb::cast<Vec>(cd["tau_fac"]).data();
        out.lam_rows = cd["lam_rows"].is_none() ? nullptr : nb::cast<Mat>(cd["lam_rows"]).data();
    }
    const int64_t* bl = blocks.data();
    P::Arm curve{L.log_rho, L.logP, L.G, nullptr, nullptr};
    const P::Grid g(lam.data(), K, kappa, od_g, od_r, n_tilt, &curve);
    std::vector<BlockResult> results(B);
    const int T = resolve_threads(n_threads, B);
    {
        std::vector<Arena> arenas(T);
        std::atomic<int64_t> next{0};
        std::mutex err_mutex;
        std::exception_ptr err;
        auto worker = [&](int tid) {
            Arena& A = arenas[tid];
            for (int64_t b; (b = next.fetch_add(1, std::memory_order_relaxed)) < B;) {
                try {
                    solve_one_block(C, p, g, L, F, bl[3 * b], bl[3 * b + 1], bl[3 * b + 2], want_capture, out, A,
                                    results[b]);
                } catch (...) {
                    std::lock_guard<std::mutex> lk(err_mutex);
                    if (!err) err = std::current_exception();
                }
            }
        };
        nb::gil_scoped_release release;
        if (T == 1) {
            worker(0);
        } else {
            std::lock_guard<std::mutex> lock(pool_mutex());
            pool_of(T).run_parallel(worker);
        }
        if (err) std::rethrow_exception(err);
    }
    // the results, under the GIL: the counts, the capture
    std::vector<int64_t> counts(static_cast<size_t>(B) * N_ASSERT * 2);
    std::vector<uint8_t> rows_del(B), cube_del(B);
    nb::list cubes, received;
    for (int b = 0; b < B; ++b) {
        BlockResult& r = results[b];
        for (int a = 0; a < N_ASSERT; ++a) {
            counts[(static_cast<size_t>(b) * N_ASSERT + a) * 2] = r.counts[a][0];
            counts[(static_cast<size_t>(b) * N_ASSERT + a) * 2 + 1] = r.counts[a][1];
        }
        rows_del[b] = r.rows_delivered; cube_del[b] = r.cube_delivered;
        if (!r.captured) { cubes.append(nb::none()); received.append(nb::none()); }
        else {
            const size_t d = r.cslot.size();
            if (!d) cubes.append(nb::none());
            else
                cubes.append(nb::make_tuple(own1(std::move(r.cslot)), own2(std::move(r.cpos), d, K),
                                            own_bool1(std::move(r.chas_pos)), own2(std::move(r.cneg), d, K),
                                            own_bool1(std::move(r.chas_neg)), own1(std::move(r.ctotal)),
                                            own1(std::move(r.copp)), own1(std::move(r.crho))));
            const int n_owned = static_cast<int>(bl[3 * b + 1] - bl[3 * b]);
            nb::list sides;
            for (int s = 0; s < 2; ++s) {
                RecvCopy& rc = r.recv[s];
                nb::dict d;
                d["has_neighbour"] = own_bool1(std::move(rc.has_nbr));
                d["has_composition"] = own_bool1(std::move(rc.has_comp));
                d["composition"] = own2(std::move(rc.comp), n_owned, K);
                d["level_gdna"] = lane_capture(std::move(rc.lanes[0]), n_owned, K);
                d["level_rna_pos"] = lane_capture(std::move(rc.lanes[1]), n_owned, K);
                d["level_rna_neg"] = lane_capture(std::move(rc.lanes[2]), n_owned, K);
                sides.append(d);
            }
            received.append(sides);
        }
    }
    nb::dict result;
    const size_t shape[3] = {static_cast<size_t>(B), static_cast<size_t>(N_ASSERT), 2};
    auto* cp = new std::vector<int64_t>(std::move(counts));
    nb::capsule owner(cp, [](void* q) noexcept { delete static_cast<std::vector<int64_t>*>(q); });
    result["counts"] = nb::ndarray<nb::numpy, int64_t, nb::ndim<3>>(cp->data(), 3, shape, std::move(owner));
    result["rows_delivered"] = own_bool1(std::move(rows_del));
    result["cube_delivered"] = own_bool1(std::move(cube_del));
    result["cubes"] = cubes;
    result["received"] = received;
    nb::list names;
    for (int a = 0; a < N_ASSERT; ++a) names.append(nb::str(ASSERTION_NAMES[a]));
    result["assertions"] = names;
    return result;
}

// ═══ THE TRANSFER KERNELS FOR THE GATES — on tables the call allocates and returns ═══════════════════════════════

// the tables of one block as the gates read them: a dict of fresh arrays (`tests/calibration/_transfer_harness.py`)
nb::dict transfer_prepare(Vec lam, BoolVec is_boundary, BoolVec is_exon, BoolVec free_pos, BoolVec free_neg, BoolVec exon_pos,
                          BoolVec exon_neg, IdxVec left, IdxVec right, FlagVec flags, Mat cnt, Mat spliced, Mat sj_count,
                          Mat sj_count_lo, Mat sj_count_hi, Mat route_rate_lo, Mat route_rate_hi, Vec eff_gdna, Vec eff_rna,
                          Vec belief_fg, BoolVec has_own_composition, nb::object factory_rows, bool has_strand, double kappa,
                          double od_g, double od_r, double rho_gdna, double rho_rna, bool split_live) {
    const int n = static_cast<int>(free_pos.shape(0)), K = static_cast<int>(lam.shape(0));
    std::vector<double> n_u(n), n_s(n), flux(n), cnt_col[2] = {std::vector<double>(n), std::vector<double>(n)};
    for (int i = 0; i < n; ++i) {
        n_u[i] = cnt.data()[2 * i] + cnt.data()[2 * i + 1];
        n_s[i] = spliced.data()[2 * i] + spliced.data()[2 * i + 1];
        flux[i] = sj_count.data()[2 * i] + sj_count.data()[2 * i + 1];
        cnt_col[0][i] = cnt.data()[2 * i]; cnt_col[1][i] = cnt.data()[2 * i + 1];
    }
    std::vector<char> is_intron(n), is_intergenic(n);
    Chain::classify(n, is_boundary.data(), is_exon.data(), free_pos.data(), free_neg.data(), is_intron.data(), is_intergenic.data());
    std::vector<const double*> src(n, nullptr);
    Mat rows_given;
    if (!factory_rows.is_none()) {
        rows_given = nb::cast<Mat>(factory_rows);
        if (static_cast<int>(rows_given.shape(0)) != n || static_cast<int>(rows_given.shape(1)) != K)
            throw std::invalid_argument("transfer_prepare: the factory rows are (n, K)");
        for (int i = 0; i < n; ++i) src[i] = rows_given.data() + static_cast<size_t>(i) * K;
    }
    Chain c;
    c.n = n; c.K = K; c.lam = lam.data();
    c.is_bnd = is_boundary.data(); c.is_exon = is_exon.data(); c.fp = free_pos.data(); c.fn = free_neg.data();
    c.exon_pos = exon_pos.data(); c.exon_neg = exon_neg.data(); c.has_own = has_own_composition.data();
    c.left = left.data(); c.right = right.data(); c.flags = flags.data();
    c.n_u = n_u.data(); c.n_s = n_s.data(); c.a_g = eff_gdna.data(); c.a_r = eff_rna.data(); c.belief = belief_fg.data();
    c.cnt = cnt.data(); c.route_lo = route_rate_lo.data(); c.route_hi = route_rate_hi.data();
    c.sj_lo = sj_count_lo.data(); c.sj_hi = sj_count_hi.data(); c.flux = flux.data();
    c.src = src.data(); c.is_intron = is_intron.data(); c.is_intergenic = is_intergenic.data();
    c.has_strand = has_strand; c.kappa = kappa; c.od_g = od_g; c.od_r = od_r;
    const size_t nk = static_cast<size_t>(n) * K;
    std::vector<double> own_rows(nk), face_store, f_par[6];
    std::vector<uint8_t> own_mask(n);
    std::vector<int8_t> f_kind(2 * n, static_cast<int8_t>(NONE));
    std::vector<int32_t> f_row(2 * n, -1), f_row2(2 * n, -1);
    for (auto& v : f_par) v.assign(2 * n, 0.0);
    RowsOut own{own_rows.data(), as_bool(own_mask), K};
    FacesOut F{&c, f_kind.data(), f_row.data(), f_row2.data(), f_par[0].data(), f_par[1].data(), f_par[2].data(),
               f_par[3].data(), f_par[4].data(), f_par[5].data(), RowStore{&face_store, K}};
    struct LaneBuf {
        std::vector<uint8_t> face, two, own_mask, wit_mask; std::vector<double> own_rows, wit, store; std::vector<int32_t> flux_index;
    } lb[3];
    LaneOut lanes[3];
    for (int l = 0; l < 3; ++l) {
        lb[l].face.assign(2 * n, 0); lb[l].two.assign(2 * n, 0); lb[l].own_mask.assign(n, 0); lb[l].wit_mask.assign(n, 0);
        lb[l].own_rows.assign(nk, 0.0); lb[l].wit.assign(2 * n, 0.0); lb[l].flux_index.assign(2 * n, -1);
        lanes[l] = LaneOut{as_bool(lb[l].face), as_bool(lb[l].two), RowsOut{lb[l].own_rows.data(), as_bool(lb[l].own_mask), K},
                           RowStore{&lb[l].store, K}, lb[l].flux_index.data(), lb[l].wit.data(), as_bool(lb[l].wit_mask)};
    }
    const bool gdna_built = rho_gdna > 0.0;
    Scratch S(K);
    prepare_block(c, own, F, gdna_built ? &lanes[0] : nullptr, lanes[1], lanes[2], rho_gdna, rho_rna, S);
    nb::dict out;
    out["lam"] = own1(std::vector<double>(lam.data(), lam.data() + K));
    out["n_u"] = own1(std::move(n_u));
    nb::dict o;
    o["rows"] = own2(std::move(own_rows), n, K);
    o["mask"] = own_bool1(std::move(own_mask));
    out["own"] = o;
    nb::dict f;
    const int n_rows = F.store.n_rows;
    f["kind"] = own2(std::move(f_kind), n, 2);
    f["row"] = own2(std::move(f_row), n, 2);
    f["row2"] = own2(std::move(f_row2), n, 2);
    const char* par_names[6] = {"n_u", "n_s", "a_b", "a_x", "width", "var"};
    for (int t = 0; t < 6; ++t) f[par_names[t]] = own2(std::move(f_par[t]), n, 2);
    face_store.resize(static_cast<size_t>(n_rows) * K);
    f["rows"] = own2(std::move(face_store), n_rows, K);
    out["faces"] = f;
    nb::dict ls;
    const char* lane_names[3] = {"gdna", "pos", "neg"};
    for (int l = 0; l < 3; ++l) {
        if (l == 0 && !gdna_built) continue;
        nb::dict d;
        d["field"] = l;
        d["face"] = own_bool2(std::move(lb[l].face), n, 2);
        d["two_sided"] = own_bool2(std::move(lb[l].two), n, 2);
        std::vector<uint8_t> empty(n);
        for (int i = 0; i < n; ++i)
            empty[i] = !(c.n_u[i] > 0.0) || !((l == 0 ? c.a_g : c.a_r)[i] > 0.0);
        d["empty"] = own_bool1(std::move(empty));
        d["own_rows"] = own2(std::move(lb[l].own_rows), n, K);
        d["own_mask"] = own_bool1(std::move(lb[l].own_mask));
        d["rho_ref"] = l == 0 ? rho_gdna : rho_rna;
        d["total"] = out["n_u"];
        if (l == 0) {
            d["count"] = out["n_u"];
            d["a"] = own1(std::vector<double>(eff_gdna.data(), eff_gdna.data() + n));
            d["other"] = nb::none();
        } else {
            const int col_read = read_column(l - 1, has_strand, kappa);
            d["count"] = own1(std::vector<double>(cnt_col[col_read]));
            d["a"] = own1(std::vector<double>(eff_rna.data(), eff_rna.data() + n));
            nb::object other = nb::none();
            if (split_live) other = nb::cast(own1(std::vector<double>(cnt_col[1 - col_read])));
            d["other"] = other;
        }
        const int n_flux = lanes[l].flux.n_rows;
        lb[l].store.resize(static_cast<size_t>(n_flux) * K);
        d["flux_rows"] = own2(std::move(lb[l].store), n_flux, K);
        d["flux_index"] = own2(std::move(lb[l].flux_index), n, 2);
        d["witness"] = own2(std::move(lb[l].wit), n, 2);
        d["witness_mask"] = own_bool1(std::move(lb[l].wit_mask));
        ls[lane_names[l]] = d;
    }
    out["lanes"] = ls;
    return out;
}

// a received table's dict (the harness's `Received.arrays()`) as the pass writes it
struct ReceivedArrays {
    BoolVec has_nbr, has_comp; Mat comp;
    struct L { BoolVec present, has_wit; Mat profile; Vec count, opp, rna, rnav; } lanes[3];
    ReceivedView view;
    explicit ReceivedArrays(nb::dict d) {
        has_nbr = nb::cast<BoolVec>(d["has_neighbour"]); has_comp = nb::cast<BoolVec>(d["has_composition"]);
        comp = nb::cast<Mat>(d["composition"]);
        const char* names[3] = {"level_gdna", "level_rna_pos", "level_rna_neg"};
        view.has_neighbour = has_nbr.data(); view.has_composition = has_comp.data(); view.composition = comp.data();
        view.K = static_cast<int>(comp.shape(1));
        for (int l = 0; l < 3; ++l) {
            nb::dict ld = nb::cast<nb::dict>(d[names[l]]);
            L& x = lanes[l];
            x.present = nb::cast<BoolVec>(ld["present"]); x.profile = nb::cast<Mat>(ld["profile"]);
            x.count = nb::cast<Vec>(ld["count"]); x.opp = nb::cast<Vec>(ld["opportunity"]);
            x.has_wit = nb::cast<BoolVec>(ld["has_witness"]); x.rna = nb::cast<Vec>(ld["rna_count"]);
            x.rnav = nb::cast<Vec>(ld["rna_count_var"]);
            view.lanes[l] = LevelsView{x.present.data(), x.profile.data(), x.count.data(), x.opp.data(), x.has_wit.data(),
                                       x.rna.data(), x.rnav.data(), view.K};
        }
    }
};

// the tables dict `transfer_prepare` returns (or the harness builds by hand), as the pass and the solve read them
struct TablesArrays {
    Vec lam; Mat own_rows; BoolVec own_mask;
    KindMat f_kind; IdxMat f_row, f_row2; Mat f_par[6]; Mat f_rows;
    struct L {
        bool built = false; int field = 0; BoolMat face, two; BoolVec empty, own_mask, wit_mask; Mat own_rows, flux_rows, wit;
        IdxMat flux_index; Vec count, a, other, total; bool has_other = false; double rho_ref = 0.0;
    } lanes[3];
    int n = 0, K = 0;
    explicit TablesArrays(nb::dict t) {
        lam = nb::cast<Vec>(t["lam"]); K = static_cast<int>(lam.shape(0));
        nb::dict o = nb::cast<nb::dict>(t["own"]);
        own_rows = nb::cast<Mat>(o["rows"]); own_mask = nb::cast<BoolVec>(o["mask"]); n = static_cast<int>(own_mask.shape(0));
        nb::dict f = nb::cast<nb::dict>(t["faces"]);
        f_kind = nb::cast<KindMat>(f["kind"]); f_row = nb::cast<IdxMat>(f["row"]); f_row2 = nb::cast<IdxMat>(f["row2"]);
        const char* par_names[6] = {"n_u", "n_s", "a_b", "a_x", "width", "var"};
        for (int q = 0; q < 6; ++q) f_par[q] = nb::cast<Mat>(f[par_names[q]]);
        f_rows = nb::cast<Mat>(f["rows"]);
        nb::dict ls = nb::cast<nb::dict>(t["lanes"]);
        const char* lane_names[3] = {"gdna", "pos", "neg"};
        for (int l = 0; l < 3; ++l) {
            if (!ls.contains(lane_names[l])) continue;
            nb::dict d = nb::cast<nb::dict>(ls[lane_names[l]]);
            L& x = lanes[l];
            x.built = true; x.field = nb::cast<int>(d["field"]);
            x.face = nb::cast<BoolMat>(d["face"]); x.two = nb::cast<BoolMat>(d["two_sided"]); x.empty = nb::cast<BoolVec>(d["empty"]);
            x.own_rows = nb::cast<Mat>(d["own_rows"]); x.own_mask = nb::cast<BoolVec>(d["own_mask"]);
            x.count = nb::cast<Vec>(d["count"]); x.a = nb::cast<Vec>(d["a"]); x.total = nb::cast<Vec>(d["total"]);
            x.has_other = !d["other"].is_none();
            if (x.has_other) x.other = nb::cast<Vec>(d["other"]);
            x.flux_rows = nb::cast<Mat>(d["flux_rows"]); x.flux_index = nb::cast<IdxMat>(d["flux_index"]);
            x.wit = nb::cast<Mat>(d["witness"]); x.wit_mask = nb::cast<BoolVec>(d["witness_mask"]);
            x.rho_ref = nb::cast<double>(d["rho_ref"]);
        }
    }
    FacesView faces() const {
        return FacesView{f_kind.data(), f_row.data(), f_row2.data(), f_par[0].data(), f_par[1].data(), f_par[2].data(),
                         f_par[3].data(), f_par[4].data(), f_par[5].data(), f_rows.data(), K};
    }
    std::vector<LaneView> lane_views() const {
        std::vector<LaneView> v;
        for (int l = 0; l < 3; ++l) {
            const L& x = lanes[l];
            if (!x.built) continue;
            v.push_back(LaneView{x.field, x.face.data(), x.two.data(), x.empty.data(), x.own_rows.data(), x.own_mask.data(),
                                 x.count.data(), x.a.data(), x.has_other ? x.other.data() : nullptr, x.wit.data(), x.wit_mask.data()});
        }
        return v;
    }
    SolveLane solve_lane(int l) const {
        const L& x = lanes[l];
        SolveLane s;
        if (!x.built) return s;
        s.built = true; s.field = x.field; s.empty = x.empty.data(); s.total = x.total.data(); s.a = x.a.data();
        s.rho_ref = x.rho_ref; s.own_rows = x.own_rows.data(); s.own_mask = x.own_mask.data();
        s.flux_store = x.flux_rows.data(); s.flux_index = x.flux_index.data();
        return s;
    }
};

// one directional pass on the tables, the received table written in place — the gates' whole passes and single hops
void transfer_pass(nb::dict tables, nb::dict received, IdxVec seq, IdxVec nbr, BoolVec terminal) {
    TablesArrays T(tables);
    ReceivedArrays R(received);
    if (R.view.K != T.K) throw std::invalid_argument("transfer_pass: the tables and the received table disagree on the grid");
    const std::vector<LaneView> lanes = T.lane_views();
    Scratch S(T.K);
    std::vector<double> out;
    pass_block(T.lam.data(), T.K, seq.data(), static_cast<int>(seq.shape(0)), nbr.data(), terminal.data(), T.own_rows.data(),
               T.own_mask.data(), T.faces(), lanes, R.view, S, out);
}

// the two received tables into ψ's channels: (live, rows (n, K), cube dict or None)
nb::tuple transfer_solve(nb::dict tables, nb::dict from_left, nb::dict from_right, BoolVec free_pos, BoolVec free_neg) {
    TablesArrays T(tables);
    ReceivedArrays Lft(from_left), Rgt(from_right);
    const int n = T.n, K = T.K;
    const ReceivedView sides[2] = {Lft.view, Rgt.view};
    const SolveLane G = T.solve_lane(0), rna[2] = {T.solve_lane(1), T.solve_lane(2)};
    std::vector<double> rows(static_cast<size_t>(n) * K, 0.0), cpos(static_cast<size_t>(n) * K), cneg(static_cast<size_t>(n) * K),
        ctotal(n), copp(n), crho(n), bound, row;
    std::vector<uint8_t> written(n), chas_pos(n), chas_neg(n);
    std::vector<int64_t> cslot(n);
    CubeOut cube{cslot.data(), cpos.data(), as_bool(chas_pos), cneg.data(), as_bool(chas_neg), ctotal.data(), copp.data(),
                 crho.data(), n};
    Scratch S(K);
    const auto [live, d] = solve_block(T.lam.data(), K, n, sides, G, rna, free_pos.data(), free_neg.data(), rows.data(),
                                       as_bool(written), cube, S, bound, row);
    nb::object cube_out = nb::none();
    if (d) {
        cslot.resize(d); cpos.resize(static_cast<size_t>(d) * K); cneg.resize(static_cast<size_t>(d) * K);
        chas_pos.resize(d); chas_neg.resize(d); ctotal.resize(d); copp.resize(d); crho.resize(d);
        nb::dict c;
        c["slot"] = own1(std::move(cslot));
        c["profile_pos"] = own2(std::move(cpos), d, K);
        c["has_pos"] = own_bool1(std::move(chas_pos));
        c["profile_neg"] = own2(std::move(cneg), d, K);
        c["has_neg"] = own_bool1(std::move(chas_neg));
        c["u"] = own1(std::vector<double>(T.lam.data(), T.lam.data() + K));
        c["total"] = own1(std::move(ctotal));
        c["opportunity"] = own1(std::move(copp));
        c["rho_ref"] = own1(std::move(crho));
        cube_out = c;
    }
    return nb::make_tuple(live, own2(std::move(rows), n, K), cube_out);
}

// ═══ THE ROW CONSTRUCTORS, THE FLAG PREDICATES AND THE ABSORBED CONSTRUCTIONS, BOUND FOR THE GATES ═════════════
// The one implementation of every row constructor (`transfer_rows.h`), flag predicate (`transfer_kernel.h`) and
// absorbed construction (above) is read by the gates through these bindings — as ψ's is through `psi_cube` — each a
// fresh array from its own arguments, on a scratch of its own. Nothing in `src/` reads them.

using Row = nb::ndarray<nb::numpy, double, nb::ndim<1>>;

Row make_row(std::vector<double>&& v) { return own1(std::move(v)); }

int K_of(const Vec& v) { return static_cast<int>(v.shape(0)); }

nb::object flank(int64_t x) { return x < 0 ? nb::none() : nb::cast(x); }

void bind_rows(nb::module_& m) {
    nb::module_ r = m.def_submodule(
        "rows", "The row constructors of transfer_rows.h, the builders' flag predicates and the constructions the block "
                "pipeline absorbed, bound for the gates.");
    r.attr("EPS") = EPS;
    r.attr("MARGINAL_NODES") = own1(std::vector<double>(MARGINAL_NODES, MARGINAL_NODES + N_MARGINAL_NODES));
    r.def("trigamma", &trigamma, nb::arg("x"), "zeta(2, x), the counting variance's one home.");
    r.def("lgamma", [](Vec x) {
        const int n = K_of(x);
        std::vector<double> out(n);
        for (int i = 0; i < n; ++i) out[i] = std::lgamma(x.data()[i]);
        return make_row(std::move(out));
    }, nb::arg("x"), "log Γ(x) per element — libm's, the kernel's own log-gamma: the factory rows' one arithmetic "
       "wherever they are built (scipy's gammaln, cephes, differs in the last bits).");
    r.def("count_logvar", &count_logvar, nb::arg("n"));
    r.def("hop_price", &hop_price, nb::arg("n_s"), nb::arg("a_s"), nb::arg("n_x"), nb::arg("a_x"));
    r.def("blur_row", [](Vec row, Vec lam, double v) {
        const int K = K_of(row);
        std::vector<double> out(K), scratch;
        blur_row(row.data(), K, lam.data(), v, out.data(), scratch);
        return make_row(std::move(out));
    }, nb::arg("row"), nb::arg("lam"), nb::arg("v"));
    r.def("lower_side", [](Vec p) {
        const int K = K_of(p);
        std::vector<double> out(K);
        lower_side(p.data(), K, out.data());
        return make_row(std::move(out));
    }, nb::arg("profile"));
    r.def("face_map_lambda", [](Vec lam, double n_u, double a_g_b, double a_r_b, double e_g_e, double e_r_e, double s) {
        const int K = K_of(lam);
        std::vector<double> out(K);
        face_map_lambda(lam.data(), K, n_u, a_g_b, a_r_b, e_g_e, e_r_e, s, out.data());
        return make_row(std::move(out));
    }, nb::arg("lam"), nb::arg("n_u"), nb::arg("a_g_b"), nb::arg("a_r_b"), nb::arg("e_g_e"), nb::arg("e_r_e"), nb::arg("s"));
    r.def("transport_row", [](Vec row, Vec lam, Vec map, double n_u, double n_s) {
        const int K = K_of(lam);
        Scratch S(K);
        std::vector<double> out(K);
        transport_row(row.data(), lam.data(), K, map.data(), n_u, n_s, S, out.data());
        return make_row(std::move(out));
    }, nb::arg("row"), nb::arg("lam"), nb::arg("lam_e_of_u"), nb::arg("n_u"), nb::arg("n_s"));
    r.def("splice_out_row", [](Vec row_e, Vec lam, double n_u, double n_s, double a_g_b, double a_g_e, Vec nodes) {
        const int K = K_of(lam);
        Scratch S(K);
        std::vector<double> out(K);
        splice_out_row(row_e.data(), lam.data(), K, n_u, n_s, a_g_b, a_g_e, nodes.data(), K_of(nodes), S, out.data());
        return make_row(std::move(out));
    }, nb::arg("row_e"), nb::arg("lam"), nb::arg("n_u"), nb::arg("n_s"), nb::arg("a_g_b"), nb::arg("a_g_e"), nb::arg("nodes"));
    r.def("level_map_lambda", [](Vec lam, double density_b, double opportunity_i, double total_i) {
        const int K = K_of(lam);
        std::vector<double> out(K);
        level_map_lambda(lam.data(), K, density_b, opportunity_i, total_i, out.data());
        return make_row(std::move(out));
    }, nb::arg("lam"), nb::arg("density_b"), nb::arg("opportunity_i"), nb::arg("total_i"));
    r.def("level_row", [](Vec row_b, Vec lam, Vec map, double v) {
        const int K = K_of(lam);
        Scratch S(K);
        std::vector<double> out(K);
        level_row(row_b.data(), lam.data(), K, map.data(), v, S, out.data());
        return make_row(std::move(out));
    }, nb::arg("row_b"), nb::arg("lam"), nb::arg("lam_i_of_b"), nb::arg("v"));
    r.def("level_bound_row", [](Vec lam, double density_b, double opportunity_i, double total_i, double v) {
        const int K = K_of(lam);
        std::vector<double> out(K);
        level_bound_row(lam.data(), K, density_b, opportunity_i, total_i, v, out.data());
        return make_row(std::move(out));
    }, nb::arg("lam"), nb::arg("density_b"), nb::arg("opportunity_i"), nb::arg("total_i"), nb::arg("v"));
    r.def("edge_level_row", [](Vec lam, double n_b, double n_e, double a_g_b, double a_g_e) {
        const int K = K_of(lam);
        std::vector<double> out(K);
        edge_level_row(lam.data(), K, n_b, n_e, a_g_b, a_g_e, out.data());
        return make_row(std::move(out));
    }, nb::arg("lam"), nb::arg("n_b"), nb::arg("n_e"), nb::arg("a_g_b"), nb::arg("a_g_e"));
    r.def("poisson_level", [](Vec u, double n, double a, double rho_ref) {
        const int K = K_of(u);
        std::vector<double> out(K);
        poisson_level(u.data(), K, n, a, rho_ref, out.data());
        return make_row(std::move(out));
    }, nb::arg("u"), nb::arg("n"), nb::arg("a"), nb::arg("rho_ref"));
    r.def("level_of_profile", [](Vec row, Vec lam, Vec u, double n, double a, double rho_ref) {
        const int K = K_of(lam);
        Scratch S(K);
        std::vector<double> out(K);
        level_of_profile(row.data(), lam.data(), u.data(), K, n, a, rho_ref, S, out.data());
        return make_row(std::move(out));
    }, nb::arg("row"), nb::arg("lam"), nb::arg("u"), nb::arg("n"), nb::arg("a"), nb::arg("rho_ref"));
    r.def("rna_level_of_profile", [](Vec row, Vec lam, Vec u, double n, double a_r, double rho_ref) {
        const int K = K_of(lam);
        Scratch S(K);
        std::vector<double> out(K);
        rna_level_of_profile(row.data(), lam.data(), u.data(), K, n, a_r, rho_ref, S, out.data());
        return make_row(std::move(out));
    }, nb::arg("row"), nb::arg("lam"), nb::arg("u"), nb::arg("n"), nb::arg("a_r"), nb::arg("rho_ref"));
    r.def("profile_of_level", [](Vec p, Vec u, Vec lam, double n, double a, double rho_ref) {
        const int K = K_of(lam);
        Scratch S(K);
        std::vector<double> out(K);
        profile_of_level(p.data(), u.data(), lam.data(), K, n, a, rho_ref, S, out.data());
        return make_row(std::move(out));
    }, nb::arg("profile"), nb::arg("u"), nb::arg("lam"), nb::arg("n"), nb::arg("a"), nb::arg("rho_ref"));
    r.def("rna_row_of_level", [](Vec p, Vec u, Vec lam, double n, double a_r, double rho_ref) {
        const int K = K_of(lam);
        Scratch S(K);
        std::vector<double> out(K);
        rna_row_of_level(p.data(), u.data(), lam.data(), K, n, a_r, rho_ref, S, out.data());
        return make_row(std::move(out));
    }, nb::arg("profile"), nb::arg("u"), nb::arg("lam"), nb::arg("n"), nb::arg("a_r"), nb::arg("rho_ref"));
    r.def("flux_level", [](Vec u, double count, double rate, double rho_ref, double v) -> nb::object {
        const int K = K_of(u);
        Scratch S(K);
        std::vector<double> out(K);
        if (!flux_level(u.data(), K, count, rate, rho_ref, v, S, out.data())) return nb::none();
        return nb::cast(make_row(std::move(out)));
    }, nb::arg("u"), nb::arg("count"), nb::arg("rate"), nb::arg("rho_ref"), nb::arg("v") = 0.0,
       "The certified flux as a lower-sided RNA level, or None for a zero count.");
    // the builders' flag predicates, on the boundary flags
    r.def("face_is_licensed", [](int64_t f, bool fp_e, bool fn_e, bool fp_i, bool fn_i) {
        return face_is_licensed(static_cast<int>(f), fp_e, fn_e, fp_i, fn_i);
    }, nb::arg("flags"), nb::arg("fp_e"), nb::arg("fn_e"), nb::arg("fp_i"), nb::arg("fn_i"));
    r.def("boundary_shares_strand", &boundary_shares_strand, nb::arg("fp_b"), nb::arg("fn_b"), nb::arg("fp_i"), nb::arg("fn_i"));
    r.def("outside_flank", [](int64_t f, int64_t left, int64_t right) {
        int64_t o, i;
        outside_flank(static_cast<int>(f), left, right, o, i);
        return nb::make_tuple(flank(o), flank(i));
    }, nb::arg("flags"), nb::arg("left"), nb::arg("right"), "(outside, inside) of a terminus boundary, or (None, None).");
    r.def("junction_exon_side", [](int64_t f, int64_t left, int64_t right) {
        return flank(junction_exon_side(static_cast<int>(f), left, right));
    }, nb::arg("flags"), nb::arg("left"), nb::arg("right"), "The flank on a junction's exon side, or None.");
    r.def("junction_flanks", [](int64_t f, int64_t left, int64_t right) {
        int64_t c, e;
        junction_flanks(static_cast<int>(f), left, right, c, e);
        return nb::make_tuple(flank(c), flank(e));
    }, nb::arg("flags"), nb::arg("left"), nb::arg("right"), "(C, E) at an exon|exon junction with no terminus, or (None, None).");
    r.def("read_column", [](int64_t col, bool has_strand, double kappa) { return read_column(static_cast<int>(col), has_strand, kappa); },
          nb::arg("col"), nb::arg("has_strand"), nb::arg("kappa"),
          "The genome-strand column strand `col`'s RNA reads on: its own when the library reads sense, else the other.");
    // the constructions the block pipeline absorbed
    r.def("factory_rows", [](BoolVec is_intron, Vec count, Vec eff, double log_mu_bg, double alpha, double size, bool informative,
                             Vec lam) {
        const int n = K_of(count), K = K_of(lam);
        std::vector<double> fg(K), out(static_cast<size_t>(n) * K, 0.0);
        for (int k = 0; k < K; ++k) fg[k] = sigmoid(lam.data()[k]);
        const Background bg{log_mu_bg, alpha, size, informative};
        for (int i = 0; i < n; ++i)
            if (is_intron.data()[i]) factory_row(bg, count.data()[i], eff.data()[i], fg.data(), K, out.data() + static_cast<size_t>(i) * K);
        return own2(std::move(out), n, K);
    }, nb::arg("is_intron"), nb::arg("count"), nb::arg("eff"), nb::arg("log_mu_bg"), nb::arg("alpha"), nb::arg("size"),
       nb::arg("informative"), nb::arg("lam"),
       "The intron factory's λ-factor rows (n, K) — log NegBinom(f_g·C; ρ_bg·E, α_eff) at the intron slots, max-normalised, "
       "zero elsewhere — the block pipeline's own construction.");
    r.def("log_negbinom", [](Vec g, Vec mu, double size) {
        const int n = K_of(g);
        std::vector<double> out(n);
        for (int i = 0; i < n; ++i) out[i] = log_negbinom(g.data()[i], mu.data()[i], size);
        return make_row(std::move(out));
    }, nb::arg("g"), nb::arg("mu"), nb::arg("size"),
       "log NegBinom(g; mean mu, size) per element, continuous g, the mean/size parameterisation; size = inf is the "
       "exact Poisson limit.");
    r.def("factor_precision", [](Mat rows, Vec lam) {
        const int m = static_cast<int>(rows.shape(0)), K = static_cast<int>(rows.shape(1));
        if (K_of(lam) != K) throw std::invalid_argument("factor_precision: the rows and the grid disagree");
        std::vector<double> lam2(K), w, out(m);
        for (int k = 0; k < K; ++k) lam2[k] = lam.data()[k] * lam.data()[k];
        for (int i = 0; i < m; ++i) out[i] = factor_precision_row(rows.data() + static_cast<size_t>(i) * K, lam.data(), lam2.data(), K, w);
        return make_row(std::move(out));
    }, nb::arg("rows"), nb::arg("lam"),
       "The composition evidence each λ-factor row carries, read off its own curvature: 1/Var_λ under the normalised row, "
       "0 for a flat one.");
    r.def("strand_evidence", [](Vec u_pos, Vec u_neg, Vec fg_loc, double kappa, double od_r, double disc) {
        const int m = K_of(u_pos);
        std::vector<double> out(m);
        for (int i = 0; i < m; ++i) out[i] = strand_evidence(u_pos.data()[i], u_neg.data()[i], fg_loc.data()[i], kappa, od_r, disc);
        return make_row(std::move(out));
    }, nb::arg("u_pos"), nb::arg("u_neg"), nb::arg("fg_loc"), nb::arg("kappa"), nb::arg("od_r"), nb::arg("disc"),
       "The reference-free strand composition evidence I_strand at fg_loc, the count overdispersed, `disc` the protocol's "
       "discriminability (0: the channel is dead).");
    r.def("pairwise_sum", [](Vec a) { return pairwise_sum(a.data(), K_of(a)); }, nb::arg("a"),
          "numpy's pairwise summation of a contiguous float64 row — the order `np.sum(axis=1)` uses.");
    r.attr("OWN_EVIDENCE_EPS") = OWN_EPS;
}

}  // namespace

NB_MODULE(_solve_impl, m) {
    m.doc() = "The calibration solve's kernels: the block pipeline in one call per sweep, ψ, the composition transfer's "
              "builders, pass and solve, and the row constructors — one module, one implementation.";
    m.def("solve_blocks", &solve_blocks, nb::arg("is_boundary"), nb::arg("is_exon"), nb::arg("free_pos"), nb::arg("free_neg"),
          nb::arg("exon_pos"), nb::arg("exon_neg"), nb::arg("terminal"), nb::arg("left"), nb::arg("right"), nb::arg("flags"),
          nb::arg("cnt"), nb::arg("spliced"), nb::arg("sj_count"), nb::arg("sj_count_lo"), nb::arg("sj_count_hi"),
          nb::arg("route_rate_lo"), nb::arg("route_rate_hi"), nb::arg("eff_gdna"), nb::arg("eff_rna"), nb::arg("belief_fpos"),
          nb::arg("belief_fneg"), nb::arg("belief_fg"), nb::arg("blocks"), nb::arg("kappa"), nb::arg("od_g"), nb::arg("od_r"),
          nb::arg("disc"), nb::arg("lam"), nb::arg("n_tilt"), nb::arg("policy"), nb::arg("has_strand"),
          nb::arg("policy_kappa"), nb::arg("policy_od_g"), nb::arg("policy_od_r"), nb::arg("rho_gdna"),
          nb::arg("rho_rna"), nb::arg("split_live"), nb::arg("gdna").none(), nb::arg("factory").none(),
          nb::arg("out_fpos"), nb::arg("out_fneg"), nb::arg("out_fg"), nb::arg("out_var"), nb::arg("out_has_composition"),
          nb::arg("diagnostics").none(), nb::arg("n_threads"),
          "Solve every locus block of the chain — the prior rows, the self-solve, the own evidence, the message layer, the "
          "final solve, the write-back, the assertions' counts — on a pool of threads, one block at a time; the belief and "
          "`has_composition` written in place on the owned slots; the counts and the capture returned. "
          "Bit-identical at every thread count.");
    m.def("psi_solve", &psi_solve, nb::arg("slots"), nb::arg("u_pos"), nb::arg("u_neg"), nb::arg("allow_pos"),
          nb::arg("allow_neg"), nb::arg("fg_ref"), nb::arg("fpos_ref"), nb::arg("fneg_ref"), nb::arg("kappa"),
          nb::arg("od_g"), nb::arg("od_r"), nb::arg("lam"), nb::arg("gdna").none(),
          nb::arg("lam_logprior").none(), nb::arg("row_logprior").none(), nb::arg("cube_slot"), nb::arg("cube_pos"),
          nb::arg("cube_has_pos"), nb::arg("cube_neg"), nb::arg("cube_has_neg"), nb::arg("cube_u"), nb::arg("cube_total"),
          nb::arg("cube_opportunity"), nb::arg("cube_rho"), nb::arg("n_tilt"), nb::arg("out_fg"), nb::arg("out_fpos"),
          nb::arg("out_fneg"), nb::arg("out_var"), nb::arg("n_threads"),
          "Solve every slot in `slots` on its own cube and write f_g, f_pos, f_neg and Var(log f_g) in place, on "
          "`n_threads` threads (0: every core) — bit-identical at every thread count. `gdna` is the fitted arm as "
          "(log_rho, logP, mass, eff) or None; `lam_logprior` / `row_logprior` the λ-factor and delivered rows or None.");
    m.def("psi_cube", &psi_cube, nb::arg("u_pos"), nb::arg("u_neg"), nb::arg("allow_pos"), nb::arg("allow_neg"),
          nb::arg("fg_ref"), nb::arg("fpos_ref"), nb::arg("fneg_ref"), nb::arg("kappa"), nb::arg("od_g"),
          nb::arg("od_r"), nb::arg("lam"), nb::arg("gdna").none(), nb::arg("lam_logprior").none(),
          nb::arg("row_logprior").none(), nb::arg("cube_slot"), nb::arg("cube_pos"), nb::arg("cube_has_pos"), nb::arg("cube_neg"),
          nb::arg("cube_has_neg"), nb::arg("cube_u"), nb::arg("cube_total"), nb::arg("cube_opportunity"),
          nb::arg("cube_rho"), nb::arg("n_tilt"), nb::arg("ambig"), nb::arg("out_psi"), nb::arg("out_fpos"),
          nb::arg("out_fneg"), nb::arg("out_tau"),
          "ψ over the (m, K, columns) cube for slots of one class, with f_pos, f_neg and the tilt beside it.");
    m.def("posterior_median", &posterior_median_rows, nb::arg("post"), nb::arg("lam"), nb::arg("out"),
          "The continuous ½-quantile of each row's λ posterior, read on λ's midpoint edges and mapped through σ.");
    m.def("compose", &compose_rows, nb::arg("f_g"), nb::arg("w_pos"), nb::arg("allow_pos"), nb::arg("allow_neg"),
          nb::arg("out_fpos"), nb::arg("out_fneg"),
          "The composition as the image of (f_g, w_pos) on the admissible strands.");
    m.def("gdna_arm", &gdna_arm_rows, nb::arg("log_rho"), nb::arg("logP"), nb::arg("lam"), nb::arg("mass"), nb::arg("eff"),
          nb::arg("out"),
          "The fitted gDNA arm as (m, K): the landscape's curve at log σ(λ) + log M − log E per slot and cell, numpy's "
          "interpolation with the ends held — the kernel's own construction, written in place for the gates.");
    m.def("transfer_prepare", &transfer_prepare, nb::arg("lam"), nb::arg("is_boundary"), nb::arg("is_exon"),
          nb::arg("free_pos"), nb::arg("free_neg"), nb::arg("exon_pos"), nb::arg("exon_neg"), nb::arg("left"),
          nb::arg("right"), nb::arg("flags"), nb::arg("cnt"), nb::arg("spliced"), nb::arg("sj_count"), nb::arg("sj_count_lo"),
          nb::arg("sj_count_hi"), nb::arg("route_rate_lo"), nb::arg("route_rate_hi"), nb::arg("eff_gdna"), nb::arg("eff_rna"),
          nb::arg("belief_fg"), nb::arg("has_own_composition"), nb::arg("factory_rows").none(), nb::arg("has_strand"),
          nb::arg("kappa"), nb::arg("od_g"), nb::arg("od_r"), nb::arg("rho_gdna"), nb::arg("rho_rna"), nb::arg("split_live"),
          "Build one block's own claims, face rules and level lanes into fresh tables and return them as a dict — `own` "
          "(rows, mask), `faces` (kind, row, row2, n_u, n_s, a_b, a_x, width, var, rows), `lanes` (gdna where the library "
          "has a gDNA coordinate, pos, neg: each its faces, two-sided faces, emptiness, own rows and mask, count, a, other, "
          "flux rows and index, witnesses), `lam`, `n_u`. For the gates.");
    m.def("transfer_pass", &transfer_pass, nb::arg("tables"), nb::arg("received"), nb::arg("seq"), nb::arg("nbr"),
          nb::arg("terminal"),
          "Run one pass in chain order on the tables: for every destination in `seq` with a neighbour `nbr[i] >= 0` that "
          "is not a terminal, apply the face's composition rule and carry each lane's level, writing the received table's "
          "arrays in place. For the gates.");
    m.def("transfer_solve", &transfer_solve, nb::arg("tables"), nb::arg("from_left"), nb::arg("from_right"), nb::arg("free_pos"),
          nb::arg("free_neg"),
          "The policy's solve for one block: the two received tables into the fused λ rows and the cube delivery at the "
          "AMBIG nodes. Returns (live, rows, cube dict or None). For the gates.");
    bind_rows(m);
}
