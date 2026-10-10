// transfer_rows.h — the ROW CONSTRUCTORS of the composition transfer, shared by the transfer's builders, pass
// and solve (`transfer_kernel.h`) and ψ's kernel (`psi_kernel.h`), and bound for the gates as `native.transfer_rows`.
// Every function here is a pure function of one face's or one node's numbers over the solve grid `lam`
// (`f_g = sigma(lam)`, K points), writing a max-normalised log-row. The gates are the transfer gates
// (`tests/calibration/test_transfer_*.py`), which recompute with these very functions.
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace transfer_rows {

constexpr double EPS = 1.0e-9;
constexpr double TINY = 1.0e-300;
//: the five kinds of composition rule a directed face can carry
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

// THE FROZEN STRAND VARIANCE (simplex_logodds: the count-zero-information freeze): the mixture's variance at
// the REFERENCE composition (f_ref, ref_pos, ref_neg) — the count sets precision, never composition.
inline double strand_variance(double n, double f_ref, double ref_pos, double ref_neg, double kappa, double od_g,
                              double od_r) {
    const double rscale = kappa * (1.0 - kappa);
    const double p_ref = 0.5 * f_ref + kappa * ref_pos + (1.0 - kappa) * ref_neg;
    const double nf = n * f_ref, np = n * ref_pos, nn = n * ref_neg;
    const double var = n * p_ref * (1.0 - p_ref) + (nf * nf) * 0.25 * od_g + (np * np) * rscale * od_r +
                       (nn * nn) * rscale * od_r;
    return std::max(var, EPS);
}
// THE STRAND TERM at one cell: the Gaussian log-likelihood of the + column count u_pos at the mean n·p,
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
    for (int i = 0; i < n; ++i) {
        const double xv = x[i];
        if (std::isnan(xv)) { out[i] = xv; continue; }
        if (xv < xp[0]) { out[i] = left; continue; }
        if (xv > xp[K - 1]) { out[i] = right; continue; }
        int lo = 0, hi = K;  // upper_bound: the largest j with xp[j] <= xv
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (xp[mid] <= xv) lo = mid + 1; else hi = mid;
        }
        const int j = lo - 1;
        if (j >= K - 1 || xv == xp[j]) { out[i] = fp[j]; continue; }
        out[i] = fp[j] + (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]) * (xv - xp[j]);
    }
}

// blur_row: the delta-method counting width — a Gaussian blur of variance v along lam on a
// max-normalised log-row, edge-padded; a non-positive v only re-normalises.
inline int blur_radius(double v, double step) {
    return std::max(static_cast<int>(std::ceil(4.0 * std::sqrt(v) / step)), 1);
}

inline void blur_row_at_step(const double* row, int K, double dlam, double v, double* out,
                             std::vector<double>& scratch) {
    const double m = vmax(row, K);
    for (int j = 0; j < K; ++j) out[j] = row[j] - m;
    if (v > 0.0 && K > 1) {
        const int half = blur_radius(v, dlam);
        const int W = 2 * half + 1;
        scratch.resize(2u * static_cast<size_t>(W) + static_cast<size_t>(K + 2 * half) + static_cast<size_t>(K));
        double* kern = scratch.data();
        double* logkern = kern + W;
        double* edged = logkern + W;
        double* pr = edged + (K + 2 * half);
        double ksum = 0.0;
        for (int t = 0; t < W; ++t) {
            const double x = (t - half) * dlam;
            logkern[t] = -0.5 * x * x / v;
            kern[t] = std::exp(logkern[t]);
            ksum += kern[t];
        }
        for (int t = 0; t < W; ++t) kern[t] /= ksum;
        const double e0 = std::exp(out[0]), e1 = std::exp(out[K - 1]);
        for (int t = 0; t < half; ++t) edged[t] = e0;
        for (int j = 0; j < K; ++j) edged[half + j] = std::exp(out[j]);
        for (int t = 0; t < half; ++t) edged[half + K + t] = e1;
        // numpy.convolve(edged, kern, "valid") reads the kernel reversed; it is symmetric term for
        // term (kern[t] and kern[W-1-t] are the same expression of the same |x|), so read it forward
        for (int j = 0; j < K; ++j) {
            double s = 0.0;
            const double* e = edged + j;
            for (int t = 0; t < W; ++t) s += e[t] * kern[t];
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
    explicit Scratch(int K) : a(K), b(K), c(K), d(K), e(K), f(K), g(K), h(K) {}
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
    // The map is the forward face map read in reverse, with EACH component's own opportunity and
    // the route-summed certified rate. Only that rate varies across the existing marginal nodes.
    // The map at every node is face_map_lambda(lam, n_u, a_g_b, a_r_b, a_g_e, a_r_e, s), with
    // its node-independent half hoisted: the gDNA arm's log and the RNA arm's unspliced density are the
    // same at every node (only the node's spliced density s joins the RNA arm), so they are computed once
    // per face — the same operations in the same order per cell, so the same bits, at half the node loop's
    // transcendentals
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
        for (int j = 0; j < K; ++j) S.d[j] += std::exp(S.c[j]);
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

// level_bound_row(lam, density_b, opportunity_i, total_i, v): the crossing total's one-sided bound.
inline void level_bound_row(const double* lam, int K, double density_b, double opportunity_i, double total_i,
                            double v, double* out) {
    const double den = std::max(v, 1e-12), log_d = std::log(density_b);
    for (int j = 0; j < K; ++j) {
        const double x = std::log(std::max(sigmoid(lam[j]) * total_i / opportunity_i, TINY)) - log_d;
        const double xp = std::max(x, 0.0);
        out[j] = -0.5 * (xp * xp) / den;
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
        const double f = std::clamp(c / n, EPS, 1.0 - EPS);
        S.b[j] = std::log(f / (1.0 - f));
    }
    interp(S.b.data(), K, lam, S.a.data(), K, S.a[0], S.a[K - 1], out);
    for (int j = 0; j < K; ++j) out[j] += total_tail(rho_ref * std::exp(u[j]) * a, n);
    norm_inplace(out, K);
}

// rna_level_of_profile(row, lam, u, n, a_r, rho_ref): a single-strand profile read as its RNA level.
inline void rna_level_of_profile(const double* row, const double* lam, const double* u, int K, double n,
                                 double a_r, double rho_ref, Scratch& S, double* out) {
    const double m = vmax(row, K);
    for (int j = 0; j < K; ++j) S.a[j] = row[j] - m;
    for (int j = 0; j < K; ++j) {
        const double c = rho_ref * std::exp(u[j]) * a_r;
        const double f_r = std::clamp(c / n, EPS, 1.0 - EPS);
        S.b[j] = std::log((1.0 - f_r) / f_r);
    }
    interp(S.b.data(), K, lam, S.a.data(), K, S.a[0], S.a[K - 1], out);
    for (int j = 0; j < K; ++j) out[j] += total_tail(rho_ref * std::exp(u[j]) * a_r, n);
    norm_inplace(out, K);
}

// profile_of_level(profile, u, lam, n, a, rho_ref): a held gDNA level read as THIS node's composition row
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
