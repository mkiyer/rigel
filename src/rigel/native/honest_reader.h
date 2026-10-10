/**
 * honest_reader.h — the honest capture reader's evaluator: a slot's log evidence over gDNA density from its own
 * strand columns with the RNA amount integrated out and its delivered local factors (the composition row, the
 * held DNA level, the held RNA levels) on fixed lattices, and the posterior MODE under the landscape.
 *
 *   log L(x) = log ∫ Pois(u; ρEg/2 + q r) Pois(v; ρEg/2 + (1−q) r) C(log(ρEg/r)) H_pos H_neg r^(−1/2) dr + D(log ρ)
 *   x* = argmax_x [log L(x) + logP(x)]  on the landscape's grid, parabola-refined
 *
 * Shared by the standalone module (honest_reader.cpp, the gates' entry point) and the block solve
 * (solve_kernel.cpp: the reader streams with the blocks after the last refit's final psi). The executable
 * specification is tests/native/_honest_reader_reference.py, held curve by curve and mode by mode in
 * tests/native/test_honest_reader.py; the lattice rules are parameterised in `Rules`, and the shipped rules
 * (the coarse lattice at every slot, the own-term window centred under the slot's q, the zero-read search) are
 * the specification's. The derived reductions behind flags are each gated separately and ship off.
 * Pure C++17, no Python.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace rigel::honest {

constexpr double kPi = 3.14159265358979323846;
constexpr double kNegInf = -std::numeric_limits<double>::infinity();

// ── the lattice rules (the prototype's constants, each a derivation already in the tree) ─────────────────────
struct Rules {
    double coarse_step;   // the count solver's lattice, 0.2 nat; the coarse RNA lattice's step
    double window_sd;     // sqrt(2T), T = -log eps64: psi's window rule, in standard deviations
    int tilt_nodes;       // psi's K_t = 24 midpoint nodes on the strand share
    bool use_coarse;      // the replica keeps the coarse lattice at every slot; production keeps it only where the
                          // window step equals the coarse step (a shallow slot's integrand spans the whole support)
    bool ystar_slot_q;    // the replica centres the window at the own-term mode under the slot's q (the prototype
                          // computed it once per slot); the alternative uses the tilt node's own assignment
    double refine;        // divides every step (the convergence check)
    bool derived_tilt;    // the share's node count from its own standard deviation (never more than tilt_nodes)
    bool skip_far_tilt;   // base share nodes beyond sqrt(2T) standard deviations of the observed share are not evaluated
    int search_stride;    // the coarse-to-fine density search for ZERO-READ slots (0: the full grid for every slot)
};

inline Rules replica_rules(double coarse_step, double window_sd, int tilt_nodes, int search_stride) {
    return Rules{coarse_step, window_sd, tilt_nodes, true, true, 1.0, false, false, search_stride};
}

struct Level {
    const double* grid = nullptr;
    const double* vals = nullptr;
    int n = 0;
    double origin = 1.0;
    bool present = false;
};

struct Slot {
    double u, v, eg, er, q;
    bool both, admits;
    const double* comp;  // n_lam values on the lambda lattice
    Level dna, pos, neg;
};

// np.interp: linear inside, held ends outside, knots increasing
inline double interp_held(double x, const double* knots, const double* vals, int n) {
    if (n == 0) return 0.0;
    if (x <= knots[0]) return vals[0];
    if (x >= knots[n - 1]) return vals[n - 1];
    int hi = int(std::upper_bound(knots, knots + n, x) - knots);
    int lo = hi - 1;
    double t = (x - knots[lo]) / (knots[hi] - knots[lo]);
    return vals[lo] + t * (vals[hi] - vals[lo]);
}

inline double logsumexp_v(const std::vector<double>& a) {
    double m = kNegInf;
    for (double x : a) if (x > m) m = x;
    if (m == kNegInf) return kNegInf;
    double s = 0.0;
    for (double x : a) s += std::exp(x - m);
    return m + std::log(s);
}

// the RNA amount (log) maximising the own term with its r^(-1/2) measure: the right-hand zero of r d/dr log f
inline double rna_mode_y(double a, double u, double v, double q, double y_lo, double y_hi) {
    double lo = y_lo, hi = y_hi;
    for (int it = 0; it < 50; ++it) {
        double mid = 0.5 * (lo + hi);
        double r = std::exp(mid);
        double g = -r + 0.5;
        if (u > 0) g += u * q * r / (a / 2 + q * r);
        if (v > 0) g += v * (1 - q) * r / (a / 2 + (1 - q) * r);
        if (g > 0) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// the per-slot evaluator; its vectors are the scratch, so no allocation happens per node once warmed
struct Evaluator {
    const Slot& s;
    const double* lam;
    int n_lam;
    const Rules& R;
    double n, sd, h, coarse_step, lg, y_hi, y_lo;
    bool comp_live;
    std::vector<double> coarse, offs;
    std::vector<double> y, logw, f, parts, hyps, phi, K;
    std::vector<double> r_coarse, comp_coarse;
    double comp_key = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> yA, fA, yB, fB, yC, fC;

    Evaluator(const Slot& slot, const double* lam_, int n_lam_, const Rules& rules)
        : s(slot), lam(lam_), n_lam(n_lam_), R(rules) {
        n = s.u + s.v;
        sd = 1.0 / std::sqrt(n + 1.0) / R.refine;
        coarse_step = R.coarse_step / R.refine;
        h = std::min(coarse_step, sd);
        lg = std::lgamma(s.u + 1.0) + std::lgamma(s.v + 1.0);
        double lo = s.comp[0], hi = s.comp[0];
        for (int i = 1; i < n_lam; ++i) { lo = std::min(lo, s.comp[i]); hi = std::max(hi, s.comp[i]); }
        comp_live = (hi - lo) > 0.0;
        y_hi = std::log(n + 12.0 * std::sqrt(n) + 12.0);
        y_lo = y_hi - 20.0;
        if (R.use_coarse || !(h < coarse_step)) {
            // np.arange(y_lo, y_hi + step, step): ceil((stop - start) / step) nodes. A shallow slot (step not below
            // the coarse step) has an integrand as broad as the support: the coarse lattice IS its window.
            int cnt = int(std::ceil((y_hi + coarse_step - y_lo) / coarse_step));
            coarse.resize(cnt);
            for (int i = 0; i < cnt; ++i) coarse[i] = y_lo + i * coarse_step;
        }
        int k = int(std::ceil(R.window_sd));
        for (int i = -k; i <= k; ++i) offs.push_back(i * h);
        r_coarse.resize(coarse.size());
        for (size_t i = 0; i < coarse.size(); ++i) r_coarse[i] = std::exp(coarse[i]);
        comp_coarse.assign(coarse.size(), 0.0);
    }

    double dna_term(double rho) const {
        if (!s.dna.present) return 0.0;
        return interp_held(std::log(rho) - std::log(s.dna.origin), s.dna.grid, s.dna.vals, s.dna.n);
    }

    static double interp_shifted(double x, const double* knots, const double* vals, int n, double shift) {
        return interp_held(x - shift, knots, vals, n);
    }

    double level_terms(double yy, double fpos, double fneg) const {
        double out = 0.0;
        if (s.pos.present) {
            if (fpos == 0.0) out += s.pos.vals[0];
            else out += interp_shifted(yy, s.pos.grid, s.pos.vals, s.pos.n, std::log(s.pos.origin) + std::log(s.er) - std::log(fpos));
        }
        if (s.neg.present) {
            if (fneg == 0.0) out += s.neg.vals[0];
            else out += interp_shifted(yy, s.neg.grid, s.neg.vals, s.neg.n, std::log(s.neg.origin) + std::log(s.er) - std::log(fneg));
        }
        return out;
    }

    // log f at every node of `y` for density dna = rho*Eg and strand assignment (qq, fpos, fneg)
    void log_f(double dna, double log_dna, double qq, double fpos, double fneg) {
        int m = int(y.size());
        f.resize(m);
        double a2 = dna / 2.0;
        for (int i = 0; i < m; ++i) {
            double yy = y[i];
            double r = std::exp(yy);
            double out = -dna - r - lg + 0.5 * yy;
            if (s.u > 0) out += s.u * std::log(a2 + qq * r);
            if (s.v > 0) out += s.v * std::log(a2 + (1.0 - qq) * r);
            if (comp_live) out += interp_held(log_dna - yy, lam, s.comp, n_lam);
            out += level_terms(yy, fpos, fneg);
            f[i] = out;
        }
    }

    // log f on the coarse nodes with the cached r and composition values (identical values, fewer evaluations)
    void log_f_coarse(double dna, double log_dna, double qq, double fpos, double fneg) {
        int m = int(coarse.size());
        f.resize(m);
        if (comp_live && log_dna != comp_key) {
            for (int i = 0; i < m; ++i) comp_coarse[i] = interp_held(log_dna - coarse[i], lam, s.comp, n_lam);
            comp_key = log_dna;
        }
        double a2 = dna / 2.0;
        for (int i = 0; i < m; ++i) {
            double yy = coarse[i];
            double r = r_coarse[i];
            double out = -dna - r - lg + 0.5 * yy;
            if (s.u > 0) out += s.u * std::log(a2 + qq * r);
            if (s.v > 0) out += s.v * std::log(a2 + (1.0 - qq) * r);
            if (comp_live) out += comp_coarse[i];
            out += level_terms(yy, fpos, fneg);
            f[i] = out;
        }
    }

    void eval_window(double centre, double dna, double log_dna, double qq, double fpos, double fneg, std::vector<double>& yy, std::vector<double>& ff) {
        yy.clear(); ff.clear();
        for (double o : offs) yy.push_back(std::min(std::max(centre + o, y_lo - 5.0), y_hi + 1.0));
        y.swap(yy); log_f(dna, log_dna, qq, fpos, fneg); y.swap(yy);
        ff = f;
    }

    static void merge3(const std::vector<double>& ya, const std::vector<double>& fa, const std::vector<double>& yb, const std::vector<double>& fb,
                       const std::vector<double>& yc, const std::vector<double>& fc, std::vector<double>& y, std::vector<double>& f) {
        y.clear(); f.clear();
        size_t ia = 0, ib = 0, ic = 0;
        const double inf = std::numeric_limits<double>::infinity();
        while (ia < ya.size() || ib < yb.size() || ic < yc.size()) {
            double va = ia < ya.size() ? ya[ia] : inf, vb = ib < yb.size() ? yb[ib] : inf, vc = ic < yc.size() ? yc[ic] : inf;
            if (va <= vb && va <= vc) { y.push_back(va); f.push_back(fa[ia++]); }
            else if (vb <= vc) { y.push_back(vb); f.push_back(fb[ib++]); }
            else { y.push_back(vc); f.push_back(fc[ic++]); }
        }
    }

    // trapezoid log-weights from the node spacings (a duplicated node, from clipping, gets zero weight)
    void weights_of(const std::vector<double>& yy) {
        int m = int(yy.size());
        logw.assign(m, kNegInf);
        if (m == 1) { logw[0] = 0.0; return; }
        auto lg0 = [](double d) { return d > 0.0 ? std::log(d) : kNegInf; };
        logw[0] = lg0(0.5 * (yy[1] - yy[0]));
        logw[m - 1] = lg0(0.5 * (yy[m - 1] - yy[m - 2]));
        for (int i = 1; i < m - 1; ++i) logw[i] = lg0(0.5 * ((yy[i] - yy[i - 1]) + (yy[i + 1] - yy[i])));
    }

    // the inner integral over y for one density node and one strand assignment: the coarse nodes, the window
    // around the own term's mode, then the window around the peak found on those; every node evaluated once
    double inner(double dna, double log_dna, double qq, double fpos, double fneg) {
        bool windowed = h < coarse_step;
        yA = coarse;
        log_f_coarse(dna, log_dna, qq, fpos, fneg);
        fA = f;
        if (!windowed) {
            weights_of(yA);
            parts.resize(yA.size());
            for (size_t i = 0; i < yA.size(); ++i) parts[i] = fA[i] + logw[i];
            return logsumexp_v(parts);
        }
        double y_star = rna_mode_y(dna, s.u, s.v, R.ystar_slot_q ? s.q : qq, y_lo - 5.0, y_hi + 1.0);
        eval_window(y_star, dna, log_dna, qq, fpos, fneg, yB, fB);
        double best_y = yA.empty() ? yB[0] : yA[0], best_f = yA.empty() ? fB[0] : fA[0];
        for (size_t i = 0; i < yA.size(); ++i) if (fA[i] > best_f) { best_f = fA[i]; best_y = yA[i]; }
        for (size_t i = 0; i < yB.size(); ++i) if (fB[i] > best_f) { best_f = fB[i]; best_y = yB[i]; }
        eval_window(best_y, dna, log_dna, qq, fpos, fneg, yC, fC);
        merge3(yA, fA, yB, fB, yC, fC, y, f);
        weights_of(y);
        parts.resize(y.size());
        for (size_t i = 0; i < y.size(); ++i) parts[i] = f[i] + logw[i];
        return logsumexp_v(parts);
    }

    void s_window(double centre, double sd_s, std::vector<double>& out) const {
        int k = int(std::ceil(R.window_sd));
        for (int i = -k; i <= k; ++i) {
            double sw = centre + i * sd_s;
            if (sw > 0.0 && sw < 1.0) out.push_back(sw);
        }
    }

    // log L at every density node of `rho`, written to `out`
    void log_L(const double* rho, int nx, double* out) {
        if (s.er <= 0.0 || !s.admits) {
            for (int ix = 0; ix < nx; ++ix) {
                double dna = rho[ix] * s.eg;
                double own = -dna - lg;
                if (s.u > 0) own += s.u * std::log(dna / 2);
                if (s.v > 0) own += s.v * std::log(dna / 2);
                out[ix] = own + (comp_live ? s.comp[n_lam - 1] : 0.0) + dna_term(rho[ix]);
            }
            return;
        }
        if (!s.both) {
            for (int ix = 0; ix < nx; ++ix) {
                double dna = rho[ix] * s.eg;
                out[ix] = inner(dna, std::log(dna), s.q, 1.0, 0.0) + dna_term(rho[ix]);
            }
            return;
        }
        // both-strand with no strand reads and no RNA level: the share integrand is a constant, the continuum
        // integrates to that constant and so does each admitted atom, and the -log 3 cancels the count of three
        if (s.u == 0.0 && s.v == 0.0 && !s.pos.present && !s.neg.present) {
            for (int ix = 0; ix < nx; ++ix) {
                double dna = rho[ix] * s.eg;
                out[ix] = inner(dna, std::log(dna), s.q, 1.0, 0.0) + dna_term(rho[ix]);
            }
            return;
        }
        // both-strand: the arcsine continuum over the share plus the pure atoms the witness bits admit
        double sd_s = 1.0 / (std::sqrt(n + 1.0) * std::max(std::fabs(2 * s.q - 1), 1.0 / std::sqrt(n + 1.0))) / R.refine;
        int m = int(R.tilt_nodes * R.refine);
        if (R.derived_tilt && !s.pos.present && !s.neg.present) {
            int md = int(std::ceil((kPi / 2) / sd_s));
            m = std::max(1, std::min(m, md));
        }
        std::vector<double> base;
        bool located = n > 0 && std::fabs(2 * s.q - 1) * std::sqrt(n + 1) > 1.0;
        double s_star = located ? (s.u / n - (1 - s.q)) / (2 * s.q - 1) : 0.5;
        for (int i = 0; i < m; ++i) {
            double ph = (i + 0.5) * (kPi / 2) / m;
            if (R.skip_far_tilt && located && !s.pos.present && !s.neg.present) {
                double sh = std::cos(ph); sh *= sh;
                if (std::fabs(sh - s_star) > R.window_sd * sd_s) continue;
            }
            base.push_back(ph);
        }
        if (located) {
            std::vector<double> win;
            s_window(s_star, sd_s, win);
            for (double sw : win) base.push_back(std::acos(std::sqrt(sw)));
        }
        std::vector<double> mixed(nx);
        std::vector<double> phi_prev, K_prev;
        auto tilt = [&](std::vector<double> ph) {
            ph.push_back(0.0); ph.push_back(kPi / 2);
            std::sort(ph.begin(), ph.end());
            ph.erase(std::unique(ph.begin(), ph.end()), ph.end());
            int np_ = int(ph.size());
            K.assign(size_t(nx) * np_, 0.0);
            int np_prev = int(phi_prev.size());
            for (int j = 0; j < np_; ++j) {
                int jp = -1;
                if (np_prev) {
                    auto it = std::lower_bound(phi_prev.begin(), phi_prev.end(), ph[j]);
                    if (it != phi_prev.end() && *it == ph[j]) jp = int(it - phi_prev.begin());
                }
                if (jp >= 0) {
                    for (int ix = 0; ix < nx; ++ix) K[size_t(ix) * np_ + j] = K_prev[size_t(ix) * np_prev + jp];
                    continue;
                }
                double sh = std::cos(ph[j]); sh *= sh;
                double qq = s.q * sh + (1.0 - s.q) * (1.0 - sh);
                for (int ix = 0; ix < nx; ++ix) {
                    double dna = rho[ix] * s.eg;
                    K[size_t(ix) * np_ + j] = inner(dna, std::log(dna), qq, sh, 1.0 - sh);
                }
            }
            phi_prev = ph; K_prev = K;
            std::vector<double> lw(np_);
            auto lg0 = [](double d) { return d > 0.0 ? std::log(d) : kNegInf; };
            lw[0] = lg0(0.5 * (ph[1] - ph[0]));
            lw[np_ - 1] = lg0(0.5 * (ph[np_ - 1] - ph[np_ - 2]));
            for (int j = 1; j < np_ - 1; ++j) lw[j] = lg0(0.5 * ((ph[j] - ph[j - 1]) + (ph[j + 1] - ph[j])));
            std::vector<double> row(np_);
            for (int ix = 0; ix < nx; ++ix) {
                for (int j = 0; j < np_; ++j) row[j] = K[size_t(ix) * np_ + j] + lw[j];
                mixed[ix] = logsumexp_v(row) + std::log(2.0 / kPi);
            }
            phi = ph;
        };
        tilt(base);
        {
            int np_ = int(phi.size());
            std::vector<double> rounded;
            for (int ix = 0; ix < nx; ++ix) {
                int best = 0;
                for (int j = 1; j < np_; ++j) if (K[size_t(ix) * np_ + j] > K[size_t(ix) * np_ + best]) best = j;
                double c = std::cos(phi[best]);
                rounded.push_back(std::round(c * c / sd_s) * sd_s);
            }
            std::sort(rounded.begin(), rounded.end());
            rounded.erase(std::unique(rounded.begin(), rounded.end()), rounded.end());
            std::vector<double> extra;
            for (double c : rounded) s_window(c, sd_s, extra);
            std::sort(extra.begin(), extra.end());
            extra.erase(std::unique(extra.begin(), extra.end()), extra.end());
            if (!extra.empty()) {
                std::vector<double> ph = base;
                for (double e : extra) ph.push_back(std::acos(std::sqrt(e)));
                tilt(ph);
            }
        }
        for (int ix = 0; ix < nx; ++ix) {
            hyps.clear();
            hyps.push_back(mixed[ix]);
            double dna = rho[ix] * s.eg;
            if (!s.neg.present) hyps.push_back(inner(dna, std::log(dna), s.q, 1.0, 0.0));
            if (!s.pos.present) hyps.push_back(inner(dna, std::log(dna), 1.0 - s.q, 0.0, 1.0));
            out[ix] = logsumexp_v(hyps) - std::log(3.0) + dna_term(rho[ix]);
        }
    }
};

// the mode of log L + logP on the grid, parabola-refined
inline double mode_of(const double* logL, const double* x, const double* logP, int nx, double* max_post) {
    int best = 0;
    double bv = kNegInf;
    for (int i = 0; i < nx; ++i) {
        double p = logL[i] + logP[i];
        if (p > bv) { bv = p; best = i; }
    }
    if (max_post) *max_post = bv;
    if (best <= 0 || best >= nx - 1) return x[best];
    double a = logL[best - 1] + logP[best - 1], b = bv, c = logL[best + 1] + logP[best + 1];
    double den = a - 2.0 * b + c;
    double off = den < 0.0 ? 0.5 * (a - c) / den : 0.0;
    return x[best] + off * (x[best + 1] - x[best]);
}

// log L at single grid nodes on demand (the search evaluates a subset)
struct Searcher {
    Evaluator& ev;
    const double* rho;
    const double* x;
    const double* logP;
    int nx;
    std::vector<double> logL;
    std::vector<uint8_t> done;
    Searcher(Evaluator& e, const double* rho_, const double* x_, const double* logP_, int nx_)
        : ev(e), rho(rho_), x(x_), logP(logP_), nx(nx_), logL(nx_, std::numeric_limits<double>::quiet_NaN()), done(nx_, 0) {}
    double at(int i) {
        if (!done[i]) { ev.log_L(rho + i, 1, &logL[i]); done[i] = 1; }
        return logL[i];
    }
    double post(int i) { return at(i) + logP[i]; }
};

// the local maxima of the landscape, shared by every slot's search
inline std::vector<int> prior_peaks_of(const double* logP, int nx) {
    std::vector<int> out;
    for (int i = 0; i < nx; ++i) {
        bool l = (i == 0) || logP[i] >= logP[i - 1], r = (i == nx - 1) || logP[i] >= logP[i + 1];
        if (l && r) out.push_back(i);
    }
    return out;
}

// the mode by a coarse-to-fine search (ZERO-READ slots only: their log L is the log of a convolution with the
// RNA-amount kernel and is smooth at the kernel's width): every `stride`-th node, the landscape's own peaks and
// the slot's own DNA witnesses as candidates, a stride-wide window around every local maximum of the candidate
// posterior, the argmax over everything evaluated, widened until it has evaluated neighbours
inline double mode_by_search(Evaluator& ev, const Slot& s, const double* rho, const double* x, const double* logP, int nx, int stride,
                             const std::vector<int>& prior_peaks, double* max_post) {
    Searcher S(ev, rho, x, logP, nx);
    std::vector<int> cand;
    for (int i = 0; i < nx; i += stride) cand.push_back(i);
    if (cand.back() != nx - 1) cand.push_back(nx - 1);
    for (int i : prior_peaks) cand.push_back(i);
    auto nearest = [&](double xv) {
        if (!(xv > -std::numeric_limits<double>::infinity())) return -1;
        int i = int(std::lower_bound(x, x + nx, xv) - x);
        if (i >= nx) i = nx - 1;
        if (i > 0 && std::fabs(x[i - 1] - xv) < std::fabs(x[i] - xv)) i = i - 1;
        return i;
    };
    if (s.eg > 0.0) {
        double dna_col = s.q >= 0.5 ? s.v : s.u;
        if (dna_col > 0.0) { int i = nearest(std::log(2.0 * dna_col / s.eg)); if (i >= 0) cand.push_back(i); }
        if (s.u + s.v > 0.0) { int i = nearest(std::log((s.u + s.v) / s.eg)); if (i >= 0) cand.push_back(i); }
    }
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
    std::vector<double> pc(cand.size());
    for (size_t k = 0; k < cand.size(); ++k) pc[k] = S.post(cand[k]);
    for (size_t k = 0; k < cand.size(); ++k) {
        bool left_ok = (k == 0) || pc[k] >= pc[k - 1];
        bool right_ok = (k + 1 == cand.size()) || pc[k] >= pc[k + 1];
        if (!(left_ok && right_ok)) continue;
        int lo = std::max(0, cand[k] - stride), hi = std::min(nx - 1, cand[k] + stride);
        for (int i = lo; i <= hi; ++i) S.at(i);
    }
    int best = -1; double bv = kNegInf;
    for (;;) {
        for (int i = 0; i < nx; ++i) if (S.done[i]) { double p = S.logL[i] + logP[i]; if (p > bv) { bv = p; best = i; } }
        bool grew = false;
        if (best > 0 && !S.done[best - 1]) { S.at(best - 1); grew = true; }
        if (best < nx - 1 && !S.done[best + 1]) { S.at(best + 1); grew = true; }
        if (!grew) break;
        bv = kNegInf; best = -1;
    }
    if (max_post) *max_post = bv;
    if (best <= 0 || best >= nx - 1) return x[best];
    double a = S.logL[best - 1] + logP[best - 1], b = bv, c = S.logL[best + 1] + logP[best + 1];
    double den = a - 2.0 * b + c;
    double off = den < 0.0 ? 0.5 * (a - c) / den : 0.0;
    return x[best] + off * (x[best + 1] - x[best]);
}

// one slot's mode under the rules: the full grid, or the search for a zero-read slot when the rules allow it
inline double slot_mode(const Slot& s, const double* lam, int n_lam, const Rules& R, const double* rho, const double* x,
                        const double* logP, int nx, const std::vector<int>& prior_peaks, std::vector<double>& logL_scratch,
                        double* max_post) {
    Evaluator ev(s, lam, n_lam, R);
    if (R.search_stride > 0 && s.u == 0.0 && s.v == 0.0)
        return mode_by_search(ev, s, rho, x, logP, nx, R.search_stride, prior_peaks, max_post);
    logL_scratch.resize(nx);
    ev.log_L(rho, nx, logL_scratch.data());
    return mode_of(logL_scratch.data(), x, logP, nx, max_post);
}

}  // namespace rigel::honest
