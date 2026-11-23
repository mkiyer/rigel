// transfer_rows.h — the row constructors of the composition transfer, shared by the transfer's builders, pass and
// solve (`transfer_kernel.h`) and ψ's kernel (`psi_kernel.h`), and bound for the gates in `_solve_impl.testing`.
// Each function is a pure function of one face's or one node's numbers over the solve grid `lam`
// (`f_g = sigma(lam)`, K points) and writes a max-normalised log-row.
// Gate: tests/calibration/test_transfer_*.py, which recompute with these functions.
#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "fast_exp.h"

namespace transfer_rows {

constexpr double EPS = 1.0e-9;
constexpr double TINY = 1.0e-300;
// NONE, or one of the five composition rules a directed face can carry; bound by name as `testing.RULE_CODES`
constexpr int NONE = 0, FORWARD = 1, TRANSPORT = 2, SPLICE_OUT = 3, EDGE = 4, LEVEL = 5;

// ---- scalar pieces ---------------------------------------------------------------------------------

// trigamma(x) = zeta(2, x): the recurrence up to x >= 12, then the asymptotic series through
// B_14 / x^15, whose truncation there is below 1e-17 relative.
inline double trigamma(double x) {
    double acc = 0.0;
    while (x < 12.0) {
        acc += 1.0 / (x * x);
        x += 1.0;
    }
    const double x2 = 1.0 / (x * x);
    double s = 1.0 / x + 0.5 * x2;
    double p = x2 / x;  // 1/x^3
    s += p / 6.0;
    p *= x2; s -= p / 30.0;
    p *= x2; s += p / 42.0;
    p *= x2; s -= p / 30.0;
    p *= x2; s += p * 5.0 / 66.0;
    p *= x2; s -= p * 691.0 / 2730.0;
    p *= x2; s += p * 7.0 / 6.0;
    return acc + s;
}

inline double count_logvar(double n) { return trigamma(n + 0.5); }

// hop_price(n_s, a_s, n_x, a_x): both counts' counting plus the discrepancy beyond it.
inline double hop_price(double n_s, double a_s, double n_x, double a_x) {
    double v = count_logvar(n_s) + count_logvar(n_x);
    if (n_s > 0.0 && n_x > 0.0) {
        const double l = std::log((n_x / a_x) / (n_s / a_s));
        v += std::max(0.0, l * l - (1.0 / n_s + 1.0 / n_x));
    }
    return v;
}

inline double sigmoid(double x) { return 1.0 / (1.0 + std::exp(-x)); }  // scipy.special.expit
// scipy.special.log_expit: log sigma(x), exact in the depleted tail (never forms 1 - sigma)
inline double log_expit(double x) { return x < 0.0 ? x - std::log1p(std::exp(x)) : -std::log1p(std::exp(-x)); }

// The strand variance, frozen at the reference composition (f_ref, ref_pos, ref_neg): the count sets the
// precision, never the composition.
inline double strand_variance(double n, double f_ref, double ref_pos, double ref_neg, double kappa) {
    const double p_ref = 0.5 * f_ref + kappa * ref_pos + (1.0 - kappa) * ref_neg;
    return std::max(n * p_ref * (1.0 - p_ref), EPS);
}
// The strand term at one cell: the Gaussian log-likelihood of the + column count u_pos at the mean n·p,
// p = ½ f_g + κ f_+ + (1 − κ) f_−, at the frozen variance.
inline double strand_term(double u_pos, double n, double p, double var, double half_log_var) {
    const double d = u_pos - n * p;
    return -0.5 * (d * d) / var - half_log_var;
}

// ---- row pieces (length K) ---------------------------------------------------------------------------

inline double vmax(const double* r, int K) {
    double m = r[0];
    for (int j = 1; j < K; ++j) if (r[j] > m) m = r[j];
    return m;
}
inline double vmin(const double* r, int K) {
    double m = r[0];
    for (int j = 1; j < K; ++j) if (r[j] < m) m = r[j];
    return m;
}
inline double ptp(const double* r, int K) { return vmax(r, K) - vmin(r, K); }
inline void norm_inplace(double* r, int K) {
    const double m = vmax(r, K);
    for (int j = 0; j < K; ++j) r[j] -= m;
}

// numpy.interp(x, xp, fp, left, right) for increasing (not necessarily strictly) xp of length K.
inline void interp(const double* x, int n, const double* xp, const double* fp, int K, double left,
                   double right, double* out) {
    int j = -1;  // the last query's interval, which an ascending query walks forward from instead of searching
    double last = 0.0;
    for (int i = 0; i < n; ++i) {
        const double xv = x[i];
        if (std::isnan(xv)) { out[i] = xv; continue; }
        if (xv < xp[0]) { out[i] = left; continue; }
        if (xv > xp[K - 1]) { out[i] = right; continue; }
        if (j >= 0 && xv >= last) {
            while (j + 1 < K && xp[j + 1] <= xv) ++j;
        } else {
            int lo = 0, hi = K;  // upper_bound: the largest j with xp[j] <= xv
            while (lo < hi) {
                const int mid = (lo + hi) >> 1;
                if (xp[mid] <= xv) lo = mid + 1; else hi = mid;
            }
            j = lo - 1;
        }
        last = xv;
        if (j >= K - 1 || xv == xp[j]) { out[i] = fp[j]; continue; }
        out[i] = fp[j] + (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]) * (xv - xp[j]);
    }
}

// blur_row: the delta-method counting width — a Gaussian blur of variance v along lam on a
// max-normalised log-row, edge-padded; a non-positive v only re-normalises. The kernel reaches
// BLUR_HALF_WIDTH_SD standard deviations each side: the Gaussian mass beyond them is erfc(4/√2) ≈ 6.3e-5.
constexpr double BLUR_HALF_WIDTH_SD = 4.0;
// The kernel's half-width in grid steps, at least one.
inline int blur_radius(double v, double step) {
    return std::max(static_cast<int>(std::ceil(BLUR_HALF_WIDTH_SD * std::sqrt(v) / step)), 1);
}

inline void blur_row_at_step(const double* row, int K, double dlam, double v, double* out,
                             std::vector<double>& scratch) {
    const double m = vmax(row, K);
    for (int j = 0; j < K; ++j) out[j] = row[j] - m;
    if (v > 0.0 && K > 1) {
        const int half = blur_radius(v, dlam);
        const int W = 2 * half + 1;
        scratch.resize(2u * static_cast<size_t>(W) + static_cast<size_t>(K + 2 * half) + 2u * static_cast<size_t>(K));
        double* kern = scratch.data();
        double* logkern = kern + W;
        double* edged = logkern + W;
        double* pr = edged + (K + 2 * half);
        double* sum = pr + K;
        for (int t = 0; t < W; ++t) {
            const double x = (t - half) * dlam;
            logkern[t] = kern[t] = -0.5 * x * x / v;
        }
        rigel::fast_exp_inplace(kern, W);
        double ksum = 0.0;
        for (int t = 0; t < W; ++t) ksum += kern[t];
        for (int t = 0; t < W; ++t) kern[t] /= ksum;
        // the row's probabilities, edge-padded
        std::copy(out, out + K, edged + half);
        rigel::fast_exp_inplace(edged + half, K);
        std::fill(edged, edged + half, edged[half]);
        std::fill(edged + half + K, edged + K + 2 * half, edged[half + K - 1]);
        // the convolution, every output at once per kernel weight (the kernel is symmetric)
        std::fill(sum, sum + K, 0.0);
        for (int t = 0; t < W; ++t) {
            const double kt = kern[t];
            const double* e = edged + t;
            for (int j = 0; j < K; ++j) sum[j] += e[j] * kt;
        }
        for (int j = 0; j < K; ++j) {
            const double s = sum[j];
            if (s >= std::numeric_limits<double>::min()) {
                pr[j] = std::log(s);
            } else {
                // Preserve the same convolution below the normal probability
                // range; a probability floor would erase relative evidence.
                double peak = -std::numeric_limits<double>::infinity();
                for (int t = 0; t < W; ++t) {
                    const int i = std::clamp(j + t - half, 0, K - 1);
                    peak = std::max(peak, out[i] + logkern[t]);
                }
                if (peak == -std::numeric_limits<double>::infinity()) {
                    pr[j] = peak;
                    continue;
                }
                double scaled = 0.0;
                for (int t = 0; t < W; ++t) {
                    const int i = std::clamp(j + t - half, 0, K - 1);
                    scaled += std::exp(out[i] + logkern[t] - peak);
                }
                pr[j] = peak + std::log(scaled) - std::log(ksum);
            }
        }
        std::copy(pr, pr + K, out);
    }
    norm_inplace(out, K);
}

inline void blur_row(const double* row, int K, const double* lam, double v, double* out,
                     std::vector<double>& scratch) {
    blur_row_at_step(row, K, K > 1 ? lam[1] - lam[0] : 0.0, v, out, scratch);
}

inline void lower_side(const double* p, int K, double* out) {
    double m = -std::numeric_limits<double>::infinity();
    for (int j = 0; j < K; ++j) {
        if (p[j] > m) m = p[j];
        out[j] = m;
    }
    norm_inplace(out, K);
}

// intersect: the pointwise minimum of the parts already accumulated in out (the first part copied
// in, the rest folded in with fold_min), then max-normalised — bounds intersect, they do not multiply.
inline void fold_min(const double* p, int K, double* out) {
    for (int j = 0; j < K; ++j) out[j] = std::min(out[j], p[j]);
}

// face_map_lambda(lam, n_u, a_g_b, a_r_b, e_g_e, e_r_e, s)
inline void face_map_lambda(const double* lam, int K, double n_u, double a_g_b, double a_r_b, double e_g_e,
                            double e_r_e, double s, double* out) {
    for (int j = 0; j < K; ++j) {
        const double sig = sigmoid(lam[j]);
        const double g_arm = n_u * sig / a_g_b * e_g_e;
        const double r_arm = (n_u * (1.0 - sig) / a_r_b + s) * e_r_e;
        out[j] = std::log(std::max(g_arm, TINY)) - std::log(std::max(r_arm, TINY));
    }
}

struct Scratch {
    std::vector<double> a, b, c, d, e, f, g, h, blur;
    std::vector<double> extended_axis, extended_source, extended_blur;
    std::vector<double> level_count;  // a level row's expected count at each grid point
    explicit Scratch(int K) : a(K), b(K), c(K), d(K), e(K), f(K), g(K), h(K), level_count(K) {}
};

// transport_row(row, lam, lam_e_of_u, n_u, n_s)
inline void transport_row(const double* row, const double* lam, int K, const double* map, double n_u,
                          double n_s, Scratch& S, double* out) {
    const double m = vmax(row, K);
    for (int j = 0; j < K; ++j) S.a[j] = row[j] - m;
    interp(lam, K, map, lam, K, lam[0], lam[K - 1], S.b.data());  // lam_u_of_x
    interp(S.b.data(), K, lam, S.a.data(), K, S.a[0], S.a[K - 1], S.c.data());
    blur_row(S.c.data(), K, lam, trigamma(n_u + 0.5) + trigamma(n_s + 0.5), out, S.blur);
}

// The reverse face uses the same component opportunities and central RNA rate as the forward map.
inline void splice_out_row(const double* row_e, const double* lam, int K, double n_u, double n_s,
                           double a_g_b, double a_g_e, double a_r_b, double a_r_e, double splice_rate,
                           const double* nodes, int n_nodes, Scratch& S,
                           double* out) {
    if (!(n_u > 0.0 && a_g_b > 0.0 && a_g_e > 0.0 && a_r_b > 0.0 && a_r_e > 0.0) || ptp(row_e, K) <= EPS) {
        std::fill(out, out + K, 0.0);
        return;
    }
    const double m = vmax(row_e, K);
    for (int j = 0; j < K; ++j) S.a[j] = row_e[j] - m;
    const double sd = std::sqrt(trigamma(n_s + 0.5) + trigamma(n_u + 0.5));
    // The forward face map read in reverse, with each component's own opportunity and the route-summed
    // certified rate; only that rate varies across the marginal nodes. The map at each node is
    // face_map_lambda(lam, n_u, a_g_b, a_r_b, a_g_e, a_r_e, s) with its node-independent half (the gDNA arm's
    // log, the RNA arm's unspliced density) computed once per face, in the same operation order, so the bits
    // match.
    double* g_log = S.g.data();
    double* r_unspl = S.h.data();
    for (int j = 0; j < K; ++j) {
        const double sig = sigmoid(lam[j]);
        g_log[j] = std::log(std::max(n_u * sig / a_g_b * a_g_e, TINY));
        r_unspl[j] = n_u * (1.0 - sig) / a_r_b;
    }
    std::fill(S.d.begin(), S.d.end(), 0.0);
    for (int t = 0; t < n_nodes; ++t) {
        const double s = splice_rate * std::exp(nodes[t] * sd);
        for (int j = 0; j < K; ++j) S.b[j] = g_log[j] - std::log(std::max((r_unspl[j] + s) * a_r_e, TINY));
        interp(S.b.data(), K, lam, S.a.data(), K, S.a[0], S.a[K - 1], S.c.data());
        rigel::fast_exp_inplace(S.c.data(), K);
        for (int j = 0; j < K; ++j) S.d[j] += S.c[j];
    }
    for (int j = 0; j < K; ++j) out[j] = std::log(std::max(S.d[j] / n_nodes, TINY));
    norm_inplace(out, K);
    if (ptp(out, K) <= EPS) std::fill(out, out + K, 0.0);
}

// level_row(row_b, lam, lam_i_of_b, v)
inline void level_row(const double* row_b, const double* lam, int K, const double* map, double v, Scratch& S,
                      double* out) {
    if (ptp(row_b, K) <= EPS) {
        std::fill(out, out + K, 0.0);
        return;
    }
    const double m = vmax(row_b, K);
    for (int j = 0; j < K; ++j) S.a[j] = row_b[j] - m;
    interp(lam, K, map, lam, K, lam[0], lam[K - 1], S.b.data());  // the preimage
    interp(S.b.data(), K, lam, S.a.data(), K, S.a[0], S.a[K - 1], S.c.data());
    blur_row(S.c.data(), K, lam, v, out, S.blur);
    if (ptp(out, K) <= EPS) std::fill(out, out + K, 0.0);
}

// ---- the builders' rows: the maps and levels a node's own numbers make -----------------------------

// level_map_lambda(lam, density_b, opportunity_i, total_i): the level-kept map lam_i(lam_b).
inline void level_map_lambda(const double* lam, int K, double density_b, double opportunity_i, double total_i,
                             double* out) {
    for (int j = 0; j < K; ++j) {
        const double f = std::clamp(sigmoid(lam[j]) * density_b * opportunity_i / total_i, EPS, 1.0 - EPS);
        out[j] = std::log(f / (1.0 - f));
    }
}

// level_bound_row(lam, density_b, opportunity_i, total_i, v): the crossing total's one-sided bound. v is a sum of
// counting variances, so positive.
inline void level_bound_row(const double* lam, int K, double density_b, double opportunity_i, double total_i,
                            double v, double* out) {
    assert(v > 0.0);
    const double log_d = std::log(density_b);
    for (int j = 0; j < K; ++j) {
        const double x = std::log(std::max(sigmoid(lam[j]) * total_i / opportunity_i, TINY)) - log_d;
        const double xp = std::max(x, 0.0);
        out[j] = -0.5 * (xp * xp) / v;
    }
    norm_inplace(out, K);
}

// edge_level_row(lam, n_b, n_e, a_g_b, a_g_e): the edge's Poisson level below, nothing above.
inline void edge_level_row(const double* lam, int K, double n_b, double n_e, double a_g_b, double a_g_e,
                           double* out) {
    if (!(n_b > 0.0)) {
        std::fill(out, out + K, 0.0);
        return;
    }
    for (int j = 0; j < K; ++j) {
        const double c = sigmoid(lam[j]) * n_e * a_g_b / a_g_e;
        out[j] = c >= n_b ? 0.0 : n_b * std::log(std::max(c, TINY) / n_b) - (c - n_b);
    }
}

// poisson_level(u, n, a, rho_ref): a structurally pure gDNA count as a level profile.
inline void poisson_level(const double* u, int K, double n, double a, double rho_ref, double* out) {
    for (int j = 0; j < K; ++j) {
        const double c = rho_ref * std::exp(u[j]) * a;
        out[j] = n > 0.0 ? n * std::log(std::max(c, TINY)) - c : -c;
    }
    norm_inplace(out, K);
}

// the Poisson tail every level's own total adds above itself: n log(c/n) - (c - n) where c >= n
inline double total_tail(double c, double n) {
    return c >= n ? n * std::log(std::max(c, TINY) / n) - (c - n) : 0.0;
}

// level_of_profile(row, lam, u, n, a, rho_ref): a composition profile read as a gDNA level.
inline void level_of_profile(const double* row, const double* lam, const double* u, int K, double n, double a,
                             double rho_ref, Scratch& S, double* out) {
    const double m = vmax(row, K);
    for (int j = 0; j < K; ++j) S.a[j] = row[j] - m;
    for (int j = 0; j < K; ++j) {
        const double c = rho_ref * std::exp(u[j]) * a;
        S.level_count[j] = c;
        const double f = std::clamp(c / n, EPS, 1.0 - EPS);
        S.b[j] = std::log(f / (1.0 - f));
    }
    interp(S.b.data(), K, lam, S.a.data(), K, S.a[0], S.a[K - 1], out);
    for (int j = 0; j < K; ++j) out[j] += total_tail(S.level_count[j], n);
    norm_inplace(out, K);
}

// rna_level_of_profile(row, lam, u, n, a_r, rho_ref): a single-strand profile read as its RNA level.
inline void rna_level_of_profile(const double* row, const double* lam, const double* u, int K, double n,
                                 double a_r, double rho_ref, Scratch& S, double* out) {
    const double m = vmax(row, K);
    for (int j = 0; j < K; ++j) S.a[j] = row[j] - m;
    for (int j = 0; j < K; ++j) {
        const double c = rho_ref * std::exp(u[j]) * a_r;
        S.level_count[j] = c;
        const double f_r = std::clamp(c / n, EPS, 1.0 - EPS);
        S.b[j] = std::log((1.0 - f_r) / f_r);
    }
    interp(S.b.data(), K, lam, S.a.data(), K, S.a[0], S.a[K - 1], out);
    for (int j = 0; j < K; ++j) out[j] += total_tail(S.level_count[j], n);
    norm_inplace(out, K);
}

// profile_of_level(profile, u, lam, n, a, rho_ref): a held gDNA level read as this node's composition row
// through its own total — a pure coordinate change, u(lam) = log(sigma(lam) n / (a rho_ref)).
inline void profile_of_level(const double* p, const double* u, const double* lam, int K, double n, double a,
                             double rho_ref, Scratch& S, double* out) {
    for (int j = 0; j < K; ++j) S.a[j] = std::log(sigmoid(lam[j]) * n / (a * rho_ref));
    interp(S.a.data(), K, u, p, K, p[0], p[K - 1], out);
    norm_inplace(out, K);
}

// rna_row_of_level(profile, u, lam, n, a_r, rho_ref): a held RNA level read as a single-strand node's
// composition row — u_s(lam) = log((1 − sigma(lam)) n / (a_r rho_ref)); a lower-only level (non-decreasing
// in u) is non-increasing in lam: "at least this much RNA" is "at most this much gDNA".
inline void rna_row_of_level(const double* p, const double* u, const double* lam, int K, double n, double a_r,
                             double rho_ref, Scratch& S, double* out) {
    for (int j = 0; j < K; ++j) S.a[j] = std::log(1.0 / (1.0 + std::exp(lam[j])) * n / (a_r * rho_ref));
    interp(S.a.data(), K, u, p, K, p[0], p[K - 1], out);
    norm_inplace(out, K);
}

// flux_level(u, count, rate, rho_ref, v): the certified flux as a lower-sided RNA level; false = no claim.
inline bool flux_level(const double* u, int K, double count, double rate, double rho_ref, double v, Scratch& S,
                       double* out) {
    if (!(count > 0.0 && rate > 0.0)) return false;
    if (v > 0.0 && K > 1) {
        // Evaluate the known likelihood across the Gaussian footprint instead
        // of repeating its endpoint values outside the retained table.
        const double step = u[1] - u[0];
        const int half = blur_radius(v, step);
        const unsigned M = static_cast<unsigned>(K) + 2u * static_cast<unsigned>(half);
        if (M < static_cast<unsigned>(K) || M > static_cast<unsigned>(std::numeric_limits<int>::max()))
            throw std::length_error("flux level grid is too large");
        S.extended_axis.resize(M);
        S.extended_source.resize(M);
        S.extended_blur.resize(M);
        for (int j = 0; j < M; ++j) S.extended_axis[j] = u[0] + (j - half) * step;
        poisson_level(S.extended_axis.data(), M, count, count / rate, rho_ref, S.extended_source.data());
        blur_row_at_step(S.extended_source.data(), M, step, v, S.extended_blur.data(), S.blur);
        std::copy(S.extended_blur.begin() + half, S.extended_blur.begin() + half + K, S.b.begin());
        lower_side(S.b.data(), K, out);
    } else {
        poisson_level(u, K, count, count / rate, rho_ref, S.a.data());
        lower_side(S.a.data(), K, out);
    }
    return true;
}

}  // namespace transfer_rows
