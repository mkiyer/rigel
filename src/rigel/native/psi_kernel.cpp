// psi_kernel.cpp — ψ, THE per-slot solve of the calibration sweep, for one block in one native call.
//
// `simplex_logodds._solve_regions_logodds_all` (the dispatcher: the reference defaults, the signal mask, the
// delivered cube rows packed) makes one call here for every slot with a live strand and a fragment; each
// slot is solved on its own (λ, θ) cube in one pass, so the read-out is chunk-exact by construction and the
// slots are independent. Per slot: ψ = the strand term (`transfer_rows.h`) + the two Jeffreys arms (½ log f_g
// + ½ log(1 − f_g)) + the fitted gDNA prior + the λ-factor row + the delivered RNA level row + the θ
// quadrature's log-weights, over K λ cells and — at an AMBIG slot — `n_tilt` θ nodes across the strand
// term's peak plus the two tilt atoms (τ = ±1), a single-strand slot's one column being its live strand;
// then the read-out: `f_g` the continuous ½-quantile of the θ-marginal on λ (through σ), `Var(log f_g)` its
// grid moment, `w_+` the RNA-mass-weighted share of the + strand over every column, the composition their
// image on the admissible strands. Gates: `tests/calibration/test_vertex_reference.py` (the cube, the
// quadrature, the atoms, the read-out, through `psi_cube`, `posterior_median` and `compose` — the same code),
// `test_sweep.py` (chunk-exactness), `profiling/sweep_replay.py replay --tolerance` on a captured sweep.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "thread_pool.h"
#include "transfer_rows.h"

namespace nb = nanobind;
using transfer_rows::EPS;
using transfer_rows::interp;
using transfer_rows::log_expit;
using transfer_rows::sigmoid;
using transfer_rows::strand_term;
using transfer_rows::strand_variance;

namespace {

using Vec = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using Mat = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using Cube = nb::ndarray<double, nb::ndim<3>, nb::c_contig, nb::device::cpu>;
using BoolVec = nb::ndarray<bool, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using IdxVec = nb::ndarray<int64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

constexpr double JEFFREYS_REF = 0.5;  // simplex_logodds._JEFFREYS_REF: the reference exponent of both arms
constexpr double PI = 3.14159265358979323846;
const double LOG_PI = std::log(PI), LOG_2 = std::log(2.0);
const double NEG_INF = -std::numeric_limits<double>::infinity();

// the θ quadrature's truncation: the strand term's mass outside the window is below double precision —
// erfc(√T) < ε₆₄ at T = −log ε₆₄, with ε₆₄ = 2⁻⁵² (numpy's finfo(float64).eps); the node count that
// resolves the peak inside that window is derived from it in the dispatcher (simplex_logodds._TILT_NODES)
const double T_NATS = -std::log(std::numeric_limits<double>::epsilon());

// a delivered RNA level row (one row of simplex_logodds.CubeRows): the held profile per strand over
// u = log(ρ/ρ_ref) — the one grid `u` every row is on — the slot's total and RNA opportunity, the lanes'
// reference density
struct Delivered {
    const double* pos; const double* neg; const double* u;  // nullptr where a strand's profile is absent
    double total, opportunity, rho_ref;
};

// ---- the θ window for one (slot, λ) --------------------------------------------------------------------
// EQUATIONS.md §9e: at fixed λ the strand term is a Gaussian in τ with centre τ̂ = d/a and width σ_τ = σ_p/|a|
// (a = (1 − f_g)(κ − ½), d = u₊/n − ½); the window is where the term lies within T of its maximum ON the
// domain — [τ̂ − ρ, τ̂ + ρ] ∩ [−1, 1] with ρ = √((clip(τ̂) − τ̂)² + 2σ_τ²T) — and the whole domain where the
// slot has no strand information.
inline void tilt_window(double u_pos, double n, double var, double fg, double kappa, double& th_lo, double& th_hi) {
    const bool counted = n > 0.0;
    const double a = (1.0 - fg) * (kappa - 0.5);
    double tau_lo = -1.0, tau_hi = 1.0;
    if (counted && a != 0.0) {
        const double d = u_pos / n - 0.5;
        const double sig_p = std::sqrt(var) / n;
        const double tau_hat = d / a;
        const double sig_tau = sig_p / std::fabs(a);
        const double c = std::clamp(tau_hat, -1.0, 1.0) - tau_hat;
        const double rho = std::sqrt(c * c + 2.0 * (sig_tau * sig_tau) * T_NATS);
        tau_lo = std::max(-1.0, tau_hat - rho);
        tau_hi = std::min(1.0, tau_hat + rho);
    }
    th_lo = std::asin(tau_lo);
    th_hi = std::asin(tau_hi);
}

// ---- one slot's cube ------------------------------------------------------------------------------------

struct SlotInputs {
    double u_pos, u_neg, fg_ref, fpos_ref, fneg_ref;
    bool ap, an;
    const double* gdna_prior;  // (K,) or nullptr
    const double* lam_prior;   // (K,) or nullptr
    const Delivered* row;      // or nullptr
};

struct Grid {
    const double* lam; int K;
    std::vector<double> fg, arm;  // σ(λ), and ½ log f_g + ½ log(1 − f_g) per cell
    double kappa, od_g, od_r;
    int n_tilt;
    explicit Grid(const double* lam_, int K_, double kappa_, double od_g_, double od_r_, int n_tilt_)
        : lam(lam_), K(K_), fg(K_), arm(K_), kappa(kappa_), od_g(od_g_), od_r(od_r_), n_tilt(n_tilt_) {
        for (int k = 0; k < K; ++k) {
            fg[k] = sigmoid(lam[k]);
            arm[k] = JEFFREYS_REF * log_expit(lam[k]) + JEFFREYS_REF * log_expit(-lam[k]);
        }
    }
    int columns(bool ambig) const { return ambig ? n_tilt + 2 : 1; }
};

// ψ over one slot's (K, C) cube, with the tilt and the two strand-fraction grids beside it; C = 1 for a
// single-strand slot (its live strand's τ = ±1), n_tilt + 2 at an AMBIG one (the windowed θ nodes with their
// trapezoid log-weights − log π, then the two atoms at log-weight 0). Every buffer is K·C long.
void slot_cube(const Grid& g, const SlotInputs& s, bool ambig, double* psi, double* f_pos, double* f_neg,
               double* tau, std::vector<double>& scratch) {
    const int K = g.K, C = g.columns(ambig);
    const double n = s.u_pos + s.u_neg;
    const double var = strand_variance(n, s.fg_ref, s.fpos_ref, s.fneg_ref, g.kappa, g.od_g, g.od_r);
    const double half_log_var = 0.5 * std::log(var);
    scratch.assign(static_cast<size_t>(K) * C, 0.0);
    double* log_w = scratch.data();
    // the tilt grid and its weights
    if (!ambig) {
        const double t = s.ap && !s.an ? 1.0 : -1.0;
        for (int k = 0; k < K; ++k) tau[k] = t;
    } else {
        const int Kt = g.n_tilt;
        for (int k = 0; k < K; ++k) {
            double th_lo, th_hi;
            tilt_window(s.u_pos, n, var, g.fg[k], g.kappa, th_lo, th_hi);
            const double h = (th_hi - th_lo) / (Kt - 1);
            const double log_h = std::log(h);
            double* tk = tau + static_cast<size_t>(k) * C;
            double* wk = log_w + static_cast<size_t>(k) * C;
            // τ = sin θ at the nodes by the rotation recurrence from (sin, cos) of θ_lo and of h — two
            // transcendentals per (slot, λ) instead of one per node, a few ulps over the window; clamped
            // to the sine's range, which the rounding can leave by an ulp at a domain end (a share
            // (1 ∓ τ)/2 below zero has no logarithm)
            double sn = std::sin(th_lo), cs = std::cos(th_lo);
            const double sh = std::sin(h), ch = std::cos(h);
            for (int t = 0; t < Kt; ++t) {
                tk[t] = std::clamp(sn, -1.0, 1.0);
                wk[t] = ((t == 0 || t == Kt - 1) ? log_h - LOG_2 : log_h) - LOG_PI;
                const double sn2 = sn * ch + cs * sh;
                cs = cs * ch - sn * sh;
                sn = sn2;
            }
            tk[Kt] = 1.0; tk[Kt + 1] = -1.0;  // the tilt atom: pure +, pure −, at log-weight 0
            wk[Kt] = 0.0; wk[Kt + 1] = 0.0;
        }
    }
    // ψ before the delivered row and the weights: strand + arms + the fitted prior + the λ-factor row
    for (int k = 0; k < K; ++k) {
        const double fg = g.fg[k], f_act = 1.0 - fg;
        const double base = g.arm[k] + (s.gdna_prior ? s.gdna_prior[k] : 0.0) + (s.lam_prior ? s.lam_prior[k] : 0.0);
        for (int t = 0; t < C; ++t) {
            const size_t at = static_cast<size_t>(k) * C + t;
            const double tv = tau[at];
            f_pos[at] = f_act * (1.0 + tv) / 2.0;
            f_neg[at] = f_act * (1.0 - tv) / 2.0;
            const double p = 0.5 * fg + g.kappa * f_pos[at] + (1.0 - g.kappa) * f_neg[at];
            psi[at] = strand_term(s.u_pos, n, p, var, half_log_var) + base;
        }
    }
    // the delivered row's map: each held profile read at the density every cell implies,
    // f_s·n/a_r relative to ρ_ref, summed over the strands and max-normalised over the whole row
    if (s.row != nullptr && (s.row->pos != nullptr || s.row->neg != nullptr)) {
        std::vector<double> rv(static_cast<size_t>(K) * C, 0.0);
        const double scale = s.row->total / s.row->opportunity, log_rho = std::log(s.row->rho_ref);
        double rmax = NEG_INF;
        for (int k = 0; k < K; ++k) {
            const double f_act = 1.0 - g.fg[k];
            for (int t = 0; t < C; ++t) {
                const size_t at = static_cast<size_t>(k) * C + t;
                double v = 0.0;
                for (int strand = 0; strand < 2; ++strand) {
                    const double* prof = strand == 0 ? s.row->pos : s.row->neg;
                    if (prof == nullptr) continue;
                    const double sign = strand == 0 ? 1.0 : -1.0;
                    const double u_s = std::log(f_act * (1.0 + sign * tau[at]) / 2.0 * scale) - log_rho;
                    double r;
                    interp(&u_s, 1, s.row->u, prof, K, prof[0], prof[K - 1], &r);
                    v += r;
                }
                rv[at] = v;
                if (v > rmax) rmax = v;
            }
        }
        for (size_t at = 0; at < rv.size(); ++at) psi[at] += rv[at] - rmax;
    }
    // the weights, and THE WITNESS: a held level on a strand certifies that strand carries RNA, so the atom
    // that puts all the RNA on the OTHER strand is out
    if (ambig) {
        const int Kt = g.n_tilt;
        for (size_t at = 0; at < static_cast<size_t>(K) * C; ++at) psi[at] += log_w[at];
        if (s.row != nullptr) {
            for (int k = 0; k < K; ++k) {
                if (s.row->neg != nullptr) psi[static_cast<size_t>(k) * C + Kt] = NEG_INF;
                if (s.row->pos != nullptr) psi[static_cast<size_t>(k) * C + Kt + 1] = NEG_INF;
            }
        }
    }
}

// ---- the read-out ---------------------------------------------------------------------------------------

// the continuous ½-quantile of a λ posterior, read on λ's midpoint edges and mapped through σ
// (simplex_logodds._posterior_median_fg; DESIGN.md §6c): the grid mass as a histogram, the crossing bin
// interpolated; a posterior with no mass reads the bin the empty CDF leaves it in
double posterior_median(const double* post, const double* lam, int K, std::vector<double>& cdf) {
    double tot = 0.0;
    for (int k = 0; k < K; ++k) tot += post[k];
    const double scale = tot > 0.0 ? 1.0 / tot : 1.0;
    cdf.resize(K + 1);
    cdf[0] = 0.0;
    double acc = 0.0;
    for (int k = 0; k < K; ++k) {
        acc += post[k] * scale;
        cdf[k + 1] = acc;
    }
    int below = 0;
    for (int k = 0; k <= K; ++k) if (cdf[k] < 0.5) ++below;
    const int j = std::clamp(below, 1, K);
    const double lo = cdf[j - 1], hi = cdf[j], span = hi - lo;
    const double t = span > 0.0 ? (0.5 - lo) / span : 0.5;
    auto edge = [&](int e) {  // the histogram edges on the uniform λ lattice, the outer half-bins mirrored
        if (e == 0) return lam[0] - 0.5 * (lam[1] - lam[0]);
        if (e == K) return lam[K - 1] + 0.5 * (lam[K - 1] - lam[K - 2]);
        return 0.5 * (lam[e - 1] + lam[e]);
    };
    const double e0 = edge(j - 1), e1 = edge(j);
    return sigmoid(e0 + t * (e1 - e0));
}

// the composition as the image of ψ's two parameters (simplex_logodds._compose; DESIGN.md §6c): the RNA
// total is 1 − f_g exactly, the share w_+ is clamped and restricted to the admissible strands
inline void compose(double f_g, double w_pos, bool ap, bool an, double& f_pos, double& f_neg) {
    const double fr = 1.0 - f_g;
    double w = std::clamp(w_pos, 0.0, 1.0);
    if (!(ap && an)) w = ap ? 1.0 : 0.0;
    f_pos = ap ? fr * w : 0.0;
    f_neg = an ? fr * (1.0 - w) : 0.0;
}

struct Scratch {
    std::vector<double> psi, f_pos, f_neg, tau, log_w, post_lam, cdf;
};

// one slot, solved: the cube, one posterior over it, the read-out
void solve_slot(const Grid& g, const SlotInputs& s, Scratch& S, double& fg_out, double& fp_out, double& fn_out,
                double& vg_out) {
    const bool ambig = s.ap && s.an;
    const int K = g.K, C = g.columns(ambig);
    const size_t cells = static_cast<size_t>(K) * C;
    S.psi.resize(cells); S.f_pos.resize(cells); S.f_neg.resize(cells); S.tau.resize(cells);
    slot_cube(g, s, ambig, S.psi.data(), S.f_pos.data(), S.f_neg.data(), S.tau.data(), S.log_w);
    // one posterior over the cube (max-shifted), its θ-marginal per λ, the RNA-mass moments
    double m = NEG_INF;
    for (size_t at = 0; at < cells; ++at) if (S.psi[at] > m) m = S.psi[at];
    if (!std::isfinite(m)) m = 0.0;
    S.post_lam.assign(K, 0.0);
    double z = 0.0, m_pos = 0.0, m_neg = 0.0;
    for (int k = 0; k < K; ++k) {
        double row = 0.0;
        for (int t = 0; t < C; ++t) {
            const size_t at = static_cast<size_t>(k) * C + t;
            const double p = std::exp(S.psi[at] - m);
            row += p;
            m_pos += p * S.f_pos[at];
            m_neg += p * S.f_neg[at];
        }
        S.post_lam[k] = row;
        z += row;
    }
    const double inv_z = z > 0.0 ? 1.0 / z : 0.0;
    for (int k = 0; k < K; ++k) S.post_lam[k] *= inv_z;
    m_pos *= inv_z; m_neg *= inv_z;
    // the λ read-out on the θ-marginal: the ½-quantile and the log-variance moment
    const double f_g = posterior_median(S.post_lam.data(), g.lam, K, S.cdf);
    double e1 = 0.0, e2 = 0.0;
    for (int k = 0; k < K; ++k) {
        const double lf = log_expit(g.lam[k]);
        e1 += S.post_lam[k] * lf;
        e2 += S.post_lam[k] * (lf * lf);
    }
    const double var_g = std::max(e2 - e1 * e1, 0.0);
    // the tilt share, and the composition as the image of the two parameters; a slot with no counts reports 0
    const double rna = m_pos + m_neg;
    const double w_pos = rna > 0.0 ? m_pos / rna : 0.5;
    const bool active = (s.u_pos + s.u_neg) > 0.0;
    if (!active) { fg_out = fp_out = fn_out = vg_out = 0.0; return; }
    fg_out = std::clamp(f_g, 0.0, 1.0);
    compose(fg_out, w_pos, s.ap, s.an, fp_out, fn_out);
    vg_out = var_g;
}

// ---- the calls -----------------------------------------------------------------------------------------

struct Rows {  // the delivered cube rows, packed by the dispatcher: parallel arrays over the delivered slots
    std::vector<int> index;  // per slot of the block: its row, or -1
    std::vector<Delivered> rows;
};

Rows unpack_rows(int m, int K, IdxVec& cube_slot, Mat& cube_pos, BoolVec& has_pos, Mat& cube_neg, BoolVec& has_neg,
                 Vec& cube_u, Vec& cube_total, Vec& cube_opportunity, Vec& cube_rho) {
    Rows R;
    R.index.assign(m, -1);
    const int d = static_cast<int>(cube_slot.shape(0));
    if (d && (static_cast<int>(cube_pos.shape(1)) != K || static_cast<int>(cube_u.shape(0)) != K))
        throw std::invalid_argument("psi: a delivered row is not on the solve grid");
    for (int r = 0; r < d; ++r) {
        const int slot = static_cast<int>(cube_slot.data()[r]);
        if (slot < 0 || slot >= m) throw std::invalid_argument("psi: a delivered row names a slot outside the block");
        R.index[slot] = r;
        R.rows.push_back(Delivered{has_pos.data()[r] ? cube_pos.data() + static_cast<size_t>(r) * K : nullptr,
                                   has_neg.data()[r] ? cube_neg.data() + static_cast<size_t>(r) * K : nullptr,
                                   cube_u.data(), cube_total.data()[r],
                                   cube_opportunity.data()[r], cube_rho.data()[r]});
    }
    return R;
}

inline const double* prior_row(Mat& holder, bool present, int i) {
    return present ? holder.data() + static_cast<size_t>(i) * holder.shape(1) : nullptr;
}

// THE THREADS. Every slot is solved on its own cube with nothing shared but the read-only grid and the
// delivered rows, and writes its own four outputs, so the slot list is pulled one slot at a time by a pool
// of threads: the same arithmetic per slot in the same order within the slot, whatever thread takes it and
// in whatever order — BIT-IDENTICAL to the serial loop at every thread count. One slot at a time (no
// chunk, no granularity constant) balances the load on asymmetric cores; an AMBIG slot is n_tilt + 2
// columns against a single-strand slot's one. The pool persists across calls (852 a sweep) and is rebuilt
// only when the budget changes; the calls are serialised on it. `n_threads` is the budget: 0 is every
// core, as the locus EM reads it.
rigel::EStepThreadPool& pool_of(int n_threads) {
    static std::unique_ptr<rigel::EStepThreadPool> pool;
    if (!pool || pool->n_threads() != n_threads) pool = std::make_unique<rigel::EStepThreadPool>(n_threads);
    return *pool;
}
std::mutex& pool_mutex() {
    static std::mutex m;
    return m;
}
int resolve_threads(int n_threads, int64_t n_slots) {
    int t = n_threads > 0 ? n_threads : static_cast<int>(std::thread::hardware_concurrency());
    if (t < 1) t = 1;
    if (static_cast<int64_t>(t) > n_slots) t = static_cast<int>(std::max<int64_t>(n_slots, 1));
    return t;
}

void psi_solve(IdxVec slots, Vec u_pos, Vec u_neg, BoolVec allow_pos, BoolVec allow_neg, Vec fg_ref, Vec fpos_ref,
               Vec fneg_ref, double kappa, double od_g, double od_r, Vec lam, nb::object gdna_logprior,
               nb::object lam_logprior, IdxVec cube_slot, Mat cube_pos, BoolVec cube_has_pos, Mat cube_neg,
               BoolVec cube_has_neg, Vec cube_u, Vec cube_total, Vec cube_opportunity, Vec cube_rho, int n_tilt,
               Vec out_fg, Vec out_fpos, Vec out_fneg, Vec out_var, int n_threads) {
    const int m = static_cast<int>(u_pos.shape(0)), K = static_cast<int>(lam.shape(0));
    if (n_tilt < 2) throw std::invalid_argument("psi_solve: the tilt needs at least two nodes");
    const bool has_g = !gdna_logprior.is_none(), has_l = !lam_logprior.is_none();
    Mat g_prior = has_g ? nb::cast<Mat>(gdna_logprior) : Mat();
    Mat l_prior = has_l ? nb::cast<Mat>(lam_logprior) : Mat();
    if ((has_g && static_cast<int>(g_prior.shape(1)) != K) || (has_l && static_cast<int>(l_prior.shape(1)) != K))
        throw std::invalid_argument("psi_solve: a prior is not on the solve grid");
    Rows R = unpack_rows(m, K, cube_slot, cube_pos, cube_has_pos, cube_neg, cube_has_neg, cube_u, cube_total,
                         cube_opportunity, cube_rho);
    const Grid g(lam.data(), K, kappa, od_g, od_r, n_tilt);
    const int64_t n_sel = static_cast<int64_t>(slots.shape(0));
    const int64_t* sl = slots.data();
    const double *up = u_pos.data(), *un = u_neg.data(), *fr = fg_ref.data(), *pr = fpos_ref.data(), *nr = fneg_ref.data();
    const bool *ap = allow_pos.data(), *an = allow_neg.data();
    double *o_fg = out_fg.data(), *o_fp = out_fpos.data(), *o_fn = out_fneg.data(), *o_v = out_var.data();
    auto solve_one = [&](int64_t q, Scratch& S) {
        const int i = static_cast<int>(sl[q]);
        SlotInputs s{up[i], un[i], fr[i], pr[i], nr[i], ap[i], an[i], prior_row(g_prior, has_g, i),
                     prior_row(l_prior, has_l, i), R.index[i] >= 0 ? &R.rows[R.index[i]] : nullptr};
        solve_slot(g, s, S, o_fg[i], o_fp[i], o_fn[i], o_v[i]);
    };
    const int T = resolve_threads(n_threads, n_sel);
    if (T == 1) {
        Scratch S;
        for (int64_t q = 0; q < n_sel; ++q) solve_one(q, S);
        return;
    }
    std::vector<Scratch> scratch(T);
    std::atomic<int64_t> next{0};
    auto worker = [&](int tid) {
        Scratch& S = scratch[tid];
        for (int64_t q; (q = next.fetch_add(1, std::memory_order_relaxed)) < n_sel;) solve_one(q, S);
    };
    nb::gil_scoped_release release;
    std::lock_guard<std::mutex> lock(pool_mutex());
    pool_of(T).run_parallel(worker);
}

// ψ itself, for the gates: the cube (m, K, C) with the two strand-fraction grids and the tilt beside it,
// every slot of one class (ambig: n_tilt + 2 columns; else one)
void psi_cube(Vec u_pos, Vec u_neg, BoolVec allow_pos, BoolVec allow_neg, Vec fg_ref, Vec fpos_ref, Vec fneg_ref,
              double kappa, double od_g, double od_r, Vec lam, nb::object gdna_logprior, nb::object lam_logprior,
              IdxVec cube_slot, Mat cube_pos, BoolVec cube_has_pos, Mat cube_neg, BoolVec cube_has_neg, Vec cube_u,
              Vec cube_total, Vec cube_opportunity, Vec cube_rho, int n_tilt, bool ambig, Cube out_psi, Cube out_fpos,
              Cube out_fneg, Cube out_tau) {
    const int m = static_cast<int>(u_pos.shape(0)), K = static_cast<int>(lam.shape(0));
    const bool has_g = !gdna_logprior.is_none(), has_l = !lam_logprior.is_none();
    Mat g_prior = has_g ? nb::cast<Mat>(gdna_logprior) : Mat();
    Mat l_prior = has_l ? nb::cast<Mat>(lam_logprior) : Mat();
    const Grid g(lam.data(), K, kappa, od_g, od_r, n_tilt);
    const int C = g.columns(ambig);
    if (static_cast<int>(out_psi.shape(1)) != K || static_cast<int>(out_psi.shape(2)) != C)
        throw std::invalid_argument("psi_cube: the output is not (m, K, columns)");
    Rows R = unpack_rows(m, K, cube_slot, cube_pos, cube_has_pos, cube_neg, cube_has_neg, cube_u, cube_total,
                         cube_opportunity, cube_rho);
    std::vector<double> scratch;
    for (int i = 0; i < m; ++i) {
        SlotInputs s{u_pos.data()[i], u_neg.data()[i], fg_ref.data()[i], fpos_ref.data()[i], fneg_ref.data()[i],
                     allow_pos.data()[i], allow_neg.data()[i], prior_row(g_prior, has_g, i),
                     prior_row(l_prior, has_l, i), R.index[i] >= 0 ? &R.rows[R.index[i]] : nullptr};
        const size_t off = static_cast<size_t>(i) * K * C;
        slot_cube(g, s, ambig, out_psi.data() + off, out_fpos.data() + off, out_fneg.data() + off, out_tau.data() + off,
                  scratch);
    }
}

void posterior_median_rows(Mat post, Vec lam, Vec out) {
    const int m = static_cast<int>(post.shape(0)), K = static_cast<int>(post.shape(1));
    std::vector<double> cdf;
    for (int i = 0; i < m; ++i) out.data()[i] = posterior_median(post.data() + static_cast<size_t>(i) * K, lam.data(), K, cdf);
}

void compose_rows(Vec f_g, Vec w_pos, BoolVec allow_pos, BoolVec allow_neg, Vec out_fpos, Vec out_fneg) {
    const int m = static_cast<int>(f_g.shape(0));
    for (int i = 0; i < m; ++i)
        compose(f_g.data()[i], w_pos.data()[i], allow_pos.data()[i], allow_neg.data()[i], out_fpos.data()[i],
                out_fneg.data()[i]);
}

}  // namespace

NB_MODULE(_psi_impl, m) {
    m.doc() = "ψ, the calibration sweep's per-slot solve on the (λ, θ) cube, and its pieces for the gates.";
    m.def("psi_solve", &psi_solve, nb::arg("slots"), nb::arg("u_pos"), nb::arg("u_neg"), nb::arg("allow_pos"),
          nb::arg("allow_neg"), nb::arg("fg_ref"), nb::arg("fpos_ref"), nb::arg("fneg_ref"), nb::arg("kappa"),
          nb::arg("od_g"), nb::arg("od_r"), nb::arg("lam"), nb::arg("gdna_logprior").none(),
          nb::arg("lam_logprior").none(), nb::arg("cube_slot"), nb::arg("cube_pos"), nb::arg("cube_has_pos"),
          nb::arg("cube_neg"), nb::arg("cube_has_neg"), nb::arg("cube_u"), nb::arg("cube_total"),
          nb::arg("cube_opportunity"), nb::arg("cube_rho"), nb::arg("n_tilt"), nb::arg("out_fg"), nb::arg("out_fpos"),
          nb::arg("out_fneg"), nb::arg("out_var"), nb::arg("n_threads"),
          "Solve every slot in `slots` on its own cube and write f_g, f_pos, f_neg and Var(log f_g) in place, on "
          "`n_threads` threads (0: every core) — bit-identical at every thread count.");
    m.def("psi_cube", &psi_cube, nb::arg("u_pos"), nb::arg("u_neg"), nb::arg("allow_pos"), nb::arg("allow_neg"),
          nb::arg("fg_ref"), nb::arg("fpos_ref"), nb::arg("fneg_ref"), nb::arg("kappa"), nb::arg("od_g"),
          nb::arg("od_r"), nb::arg("lam"), nb::arg("gdna_logprior").none(), nb::arg("lam_logprior").none(),
          nb::arg("cube_slot"), nb::arg("cube_pos"), nb::arg("cube_has_pos"), nb::arg("cube_neg"),
          nb::arg("cube_has_neg"), nb::arg("cube_u"), nb::arg("cube_total"), nb::arg("cube_opportunity"),
          nb::arg("cube_rho"), nb::arg("n_tilt"), nb::arg("ambig"), nb::arg("out_psi"), nb::arg("out_fpos"),
          nb::arg("out_fneg"), nb::arg("out_tau"),
          "ψ over the (m, K, columns) cube for slots of one class, with f_pos, f_neg and the tilt beside it.");
    m.def("posterior_median", &posterior_median_rows, nb::arg("post"), nb::arg("lam"), nb::arg("out"),
          "The continuous ½-quantile of each row's λ posterior, read on λ's midpoint edges and mapped through σ.");
    m.def("compose", &compose_rows, nb::arg("f_g"), nb::arg("w_pos"), nb::arg("allow_pos"), nb::arg("allow_neg"),
          nb::arg("out_fpos"), nb::arg("out_fneg"),
          "The composition as the image of (f_g, w_pos) on the admissible strands.");
}
