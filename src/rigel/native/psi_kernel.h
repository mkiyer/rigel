// psi_kernel.h — ψ, THE per-slot solve of the calibration sweep, on one slot's own (λ, θ) cube.
//
// Every slot with a live strand and a fragment is solved on its own cube in one pass, so the read-out is chunk-exact
// by construction and the slots are independent — of each other and of the thread that takes them. Per slot: ψ = the
// strand term (`transfer_rows.h`) + the two Jeffreys arms (½ log f_g + ½ log(1 − f_g)) + the fitted gDNA arm — the
// landscape's curve read at the density each cell implies (`Arm`; no (n, K) matrix of it ever exists) — + the λ-factor
// row and the delivered composition row, the two added per cell, + the delivered RNA level row + the θ quadrature's
// log-weights, over K λ cells and — at an AMBIG slot — `n_tilt` θ nodes across the strand term's peak plus the two tilt
// atoms (τ = ±1), a single-strand slot's one column being its live strand; then the read-out: `f_g` the continuous
// ½-quantile of the θ-marginal on λ (through σ), `Var(log f_g)` its grid moment, `w_+` the RNA-mass-weighted share of
// the + strand over every column, the composition their image on the admissible strands. The block pipeline
// (`solve_kernel.cpp`) calls `solve_slot` for every slot of a block; the bindings for the gates
// (`psi_cube`, `posterior_median`, `compose`, `gdna_arm`) read the same code. Gates:
// `tests/calibration/test_vertex_reference.py`, `test_sweep.py` (chunk- and thread-exactness),
// `profiling/sweep_replay.py replay --tolerance` on a captured sweep.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "transfer_rows.h"

namespace psi_kernel {

using transfer_rows::interp;
using transfer_rows::log_expit;
using transfer_rows::sigmoid;
using transfer_rows::strand_term;
using transfer_rows::strand_variance;

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

// THE FITTED gDNA ARM: the landscape's curve — `logP` over a natural-log density grid `log_rho`
// (`landscape.DensityLandscape`) — read at the density every cell implies, log ρ_c = log f_g + log M − log E (the
// slot's unspliced count on its gDNA opportunity), with numpy's interpolation and its end rules: the ends held
// constant off the grid ("no more information out here", never a linear extension of the last slope). Bare — no
// reference, no measure term, no Jacobian: the curve is a density in log-rate, so its conversion to a linear-rate
// density cancels the change of variable exactly, per component.
constexpr double ARM_EPS = 1.0e-12;  // landscape._EPS: the clips on the fraction, the mass and the opportunity

inline double arm_abscissa(double fg) { return std::log(std::clamp(fg, ARM_EPS, 1.0 - ARM_EPS)); }

struct Arm {
    const double* log_rho = nullptr; const double* logP = nullptr; int G = 0;
    const double* mass = nullptr; const double* eff = nullptr;  // the per-slot gDNA support
    bool built() const { return log_rho != nullptr; }
    // the slot's shift of a cell's log-fraction onto the density axis: log M − log E
    static double shift_of(double mass_i, double eff_i) {
        return std::log(std::max(mass_i, ARM_EPS)) - std::log(std::max(eff_i, ARM_EPS));
    }
    double shift(int i) const { return shift_of(mass[i], eff[i]); }
    double at(double log_frac, double shift) const {
        const double x = log_frac + shift;
        double out;
        interp(&x, 1, log_rho, logP, G, logP[0], logP[G - 1], &out);
        return out;
    }
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
    bool has_arm; double arm_shift;  // the fitted gDNA arm at this slot — its log M − log E — or no arm
    const double* lam_prior;   // (K,) the λ-factor row, or nullptr
    const double* row_prior;   // (K,) the delivered composition row, or nullptr
    const Delivered* row;      // the delivered RNA level row, or nullptr
};

struct Grid {
    const double* lam; int K;
    // σ(λ); the two Jeffreys arms ½ log f_g + ½ log(1 − f_g); log of the clipped σ(λ), the fitted arm's abscissa
    std::vector<double> fg, jeffreys, log_frac;
    double kappa, od_g, od_r;
    int n_tilt;
    const Arm* arm;  // the fitted gDNA arm, built or not
    explicit Grid(const double* lam_, int K_, double kappa_, double od_g_, double od_r_, int n_tilt_, const Arm* arm_)
        : lam(lam_), K(K_), fg(K_), jeffreys(K_), log_frac(K_), kappa(kappa_), od_g(od_g_), od_r(od_r_),
          n_tilt(n_tilt_), arm(arm_) {
        for (int k = 0; k < K; ++k) {
            fg[k] = sigmoid(lam[k]);
            jeffreys[k] = JEFFREYS_REF * log_expit(lam[k]) + JEFFREYS_REF * log_expit(-lam[k]);
            log_frac[k] = arm_abscissa(fg[k]);
        }
    }
    int columns(bool ambig) const { return ambig ? n_tilt + 2 : 1; }
};

// ψ over one slot's (K, C) cube, with the tilt and the two strand-fraction grids beside it; C = 1 for a
// single-strand slot (its live strand's τ = ±1), n_tilt + 2 at an AMBIG one (the windowed θ nodes with their
// trapezoid log-weights − log π, then the two atoms at log-weight 0). Every buffer is K·C long.
inline void slot_cube(const Grid& g, const SlotInputs& s, bool ambig, double* psi, double* f_pos, double* f_neg,
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
    // ψ before the delivered level row and the weights: strand + the Jeffreys arms + the fitted arm at the cell's
    // density + (the λ-factor row + the delivered composition row)
    for (int k = 0; k < K; ++k) {
        const double fg = g.fg[k], f_act = 1.0 - fg;
        const double base = g.jeffreys[k] + (s.has_arm ? g.arm->at(g.log_frac[k], s.arm_shift) : 0.0) +
                            ((s.lam_prior ? s.lam_prior[k] : 0.0) + (s.row_prior ? s.row_prior[k] : 0.0));
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
inline double posterior_median(const double* post, const double* lam, int K, std::vector<double>& cdf) {
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
inline void solve_slot(const Grid& g, const SlotInputs& s, Scratch& S, double& fg_out, double& fp_out, double& fn_out,
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

// ---- the delivered cube rows, unpacked from the table ψ's dispatcher and the block pipeline hand over -------

struct Rows {  // per slot its row, or -1; the rows themselves
    std::vector<int> index;
    std::vector<Delivered> rows;
};

// `d` delivered rows over `m` slots: the slot of each, the two profiles with their presence bits, the one grid, the
// totals, the opportunities, the reference densities — the arrays of simplex_logodds.CubeRows
inline Rows unpack_rows(int m, int K, int d, const int64_t* slot, const double* pos, const bool* has_pos, const double* neg,
                        const bool* has_neg, const double* u, const double* total, const double* opportunity,
                        const double* rho) {
    Rows R;
    R.index.assign(m, -1);
    for (int r = 0; r < d; ++r) {
        const int i = static_cast<int>(slot[r]);
        if (i < 0 || i >= m) throw std::invalid_argument("psi: a delivered row names a slot outside the block");
        R.index[i] = r;
        R.rows.push_back(Delivered{has_pos[r] ? pos + static_cast<size_t>(r) * K : nullptr,
                                   has_neg[r] ? neg + static_cast<size_t>(r) * K : nullptr, u, total[r], opportunity[r],
                                   rho[r]});
    }
    return R;
}

}  // namespace psi_kernel
