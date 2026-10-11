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
 * tests/native/test_honest_reader.py; the lattice rules (`Rules`) are the specification's. The evaluation is
 * organised for speed without changing a value: the density-node loop is
 * outermost and everything that depends on the density node alone (the composition row, the own term's mode and
 * its window) or on the strand assignment alone (the held levels on the coarse lattice) is computed once; the
 * Poisson logs and the trapezoid sums' exponentials run over the lattice through vmath.h.
 * Pure C++17, no Python.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "vmath.h"

namespace rigel::honest {

constexpr double kPi = 3.14159265358979323846;
constexpr double kNegInf = -std::numeric_limits<double>::infinity();

// ── the lattice rules (each a derivation already in the tree) ────────────────────────────────────────────────
struct Rules {
    double coarse_step;   // the count solver's lattice, 0.2 nat; the coarse RNA lattice's step, kept at every slot
    double window_sd;     // sqrt(2T), T = -log eps64: psi's window rule, in standard deviations
    int tilt_nodes;       // psi's K_t = 24 midpoint nodes on the strand share
    int search_stride;    // the coarse-to-fine density search for ZERO-READ slots (0: the full grid for every slot)
};

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

// np.interp as above, with the search started at `hint` (the interval found for the previous argument of an ordered
// run) and moved to the argument's own interval: the same interval `interp_held` finds, so the same value bit for bit
inline double interp_held_from(double x, const double* knots, const double* vals, int n, int& hint) {
    if (n == 0) return 0.0;
    if (x <= knots[0]) { hint = 0; return vals[0]; }
    if (x >= knots[n - 1]) { hint = n - 2; return vals[n - 1]; }
    int lo = hint < 0 ? 0 : (hint > n - 2 ? n - 2 : hint);
    while (lo > 0 && knots[lo] > x) --lo;
    while (lo < n - 2 && knots[lo + 1] <= x) ++lo;
    hint = lo;
    double t = (x - knots[lo]) / (knots[lo + 1] - knots[lo]);
    return vals[lo] + t * (vals[lo + 1] - vals[lo]);
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

// a strand assignment of the RNA amount: a node of the share's arcsine lattice (`phi`), or a pure atom. Its held
// RNA-level terms on the coarse lattice do not depend on the density node, so they are evaluated once per slot
// and read at every density node.
struct Assignment {
    double phi;  // the lattice angle (an atom carries NaN)
    double qq, fpos, fneg;
    double shift_pos = 0.0, shift_neg = 0.0;  // log(origin) + log(Er) - log(f), the level's argument shift
    std::vector<double> lev;                  // level terms on the coarse lattice
    bool lev_ready = false;
};

// the per-slot evaluator; its vectors are the scratch, so no allocation happens per node once warmed. The
// density-node loop is outermost: the composition row on the coarse lattice, the own term's mode and the window
// around it (its nodes, their r and composition values) are the node's context, set once and read by every
// strand assignment evaluated there; a value that depends on neither (the level terms on the coarse lattice,
// the coarse lattice's weights) is computed once per slot. Every expression is the specification's, term by
// term in the same order, so the result is bit-identical to evaluating each node from scratch.
struct Evaluator {
    const Slot& s;
    const double* lam;
    int n_lam;
    const Rules& R;
    double n, sd, h, coarse_step, lg, y_hi, y_lo;
    bool comp_live, windowed;
    double c_pos = 0.0, c_neg = 0.0;  // log(origin) + log(Er) per delivered RNA level
    std::vector<double> coarse, offs, r_coarse, logw_coarse;  // logw_coarse: the coarse lattice's trapezoid weights
    // the density node's context
    double ctx_dna = std::numeric_limits<double>::quiet_NaN(), ctx_log_dna = 0.0, ctx_a2 = 0.0;
    std::vector<double> comp_row, yB, rB, compB;
    // scratch
    std::vector<double> y, logw, f, parts, hyps, fA, fB, yC, fC;
    std::vector<double> argA, argB, la, lb, ex, rC, compC;  // the Poisson pair's rates and logs, the exponentials, window C's values
    std::vector<double> K;
    // the share lattice at integer multiples of its resolution (the second pass's windows), cached
    std::vector<Assignment> share_cache;
    std::vector<uint8_t> share_ready;

    Evaluator(const Slot& slot, const double* lam_, int n_lam_, const Rules& rules)
        : s(slot), lam(lam_), n_lam(n_lam_), R(rules) {
        n = s.u + s.v;
        sd = 1.0 / std::sqrt(n + 1.0);
        coarse_step = R.coarse_step;
        h = std::min(coarse_step, sd);
        windowed = h < coarse_step;
        lg = std::lgamma(s.u + 1.0) + std::lgamma(s.v + 1.0);
        double lo = s.comp[0], hi = s.comp[0];
        for (int i = 1; i < n_lam; ++i) { lo = std::min(lo, s.comp[i]); hi = std::max(hi, s.comp[i]); }
        comp_live = (hi - lo) > 0.0;
        y_hi = std::log(n + 12.0 * std::sqrt(n) + 12.0);
        y_lo = y_hi - 20.0;
        {
            // np.arange(y_lo, y_hi + step, step): ceil((stop - start) / step) nodes, at every slot; a shallow slot
            // (window step not below the coarse step) has an integrand as broad as the support, and the coarse
            // lattice IS its window
            int cnt = int(std::ceil((y_hi + coarse_step - y_lo) / coarse_step));
            coarse.resize(cnt);
            for (int i = 0; i < cnt; ++i) coarse[i] = y_lo + i * coarse_step;
        }
        int k = int(std::ceil(R.window_sd));
        for (int i = -k; i <= k; ++i) offs.push_back(i * h);
        r_coarse.resize(coarse.size());
        for (size_t i = 0; i < coarse.size(); ++i) r_coarse[i] = std::exp(coarse[i]);
        weights_of(coarse);
        logw_coarse = logw;
        comp_row.assign(coarse.size(), 0.0);
        yB.resize(offs.size()); rB.resize(offs.size()); compB.assign(offs.size(), 0.0);
        if (s.pos.present) c_pos = std::log(s.pos.origin) + std::log(s.er);
        if (s.neg.present) c_neg = std::log(s.neg.origin) + std::log(s.er);
    }

    double dna_term(double rho) const {
        if (!s.dna.present) return 0.0;
        return interp_held(std::log(rho) - std::log(s.dna.origin), s.dna.grid, s.dna.vals, s.dna.n);
    }

    // the held RNA levels' terms at one RNA amount under an assignment (its shifts already carry -log f)
    double level_terms(double yy, const Assignment& a) const {
        double out = 0.0;
        if (s.pos.present) {
            if (a.fpos == 0.0) out += s.pos.vals[0];
            else out += interp_held(yy - a.shift_pos, s.pos.grid, s.pos.vals, s.pos.n);
        }
        if (s.neg.present) {
            if (a.fneg == 0.0) out += s.neg.vals[0];
            else out += interp_held(yy - a.shift_neg, s.neg.grid, s.neg.vals, s.neg.n);
        }
        return out;
    }

    Assignment assignment(double phi, double qq, double fpos, double fneg) const {
        Assignment a;
        a.phi = phi; a.qq = qq; a.fpos = fpos; a.fneg = fneg;
        if (s.pos.present && fpos != 0.0) a.shift_pos = c_pos - std::log(fpos);
        if (s.neg.present && fneg != 0.0) a.shift_neg = c_neg - std::log(fneg);
        return a;
    }

    Assignment node_of(double phi) const {
        double sh = std::cos(phi); sh *= sh;
        double qq = s.q * sh + (1.0 - s.q) * (1.0 - sh);
        return assignment(phi, qq, sh, 1.0 - sh);
    }

    void ready(Assignment& a) {
        if (a.lev_ready) return;
        a.lev.resize(coarse.size());
        int hp = 0, hn = 0;
        for (size_t i = 0; i < coarse.size(); ++i) a.lev[i] = level_terms_from(coarse[i], a, hp, hn);
        a.lev_ready = true;
    }

    // the window's nodes around a centre: positions clipped to the support's margins, their r and composition
    void prepare_window(double centre) {
        const int w = int(offs.size());
        for (int k = 0; k < w; ++k) yB[k] = std::min(std::max(centre + offs[k], y_lo - 5.0), y_hi + 1.0);
        vmath::vexp(yB.data(), rB.data(), w);
        if (comp_live) {
            int hint = n_lam - 2;
            for (int k = 0; k < w; ++k) compB[k] = interp_held_from(ctx_log_dna - yB[k], lam, s.comp, n_lam, hint);
        }
    }

    // the density node's context: dna = rho * Eg; the own term's window is centred under the slot's q
    void set_node(double dna, double log_dna) {
        ctx_dna = dna; ctx_log_dna = log_dna; ctx_a2 = dna / 2.0;
        if (comp_live) {
            int hint = n_lam - 2;
            for (size_t i = 0; i < coarse.size(); ++i) comp_row[i] = interp_held_from(log_dna - coarse[i], lam, s.comp, n_lam, hint);
        }
        if (windowed) prepare_window(rna_mode_y(dna, s.u, s.v, s.q, y_lo - 5.0, y_hi + 1.0));
    }

    // log-sum-exp with the exponentials over the array (vmath::vexp)
    double lse(const std::vector<double>& a) {
        double m = kNegInf;
        for (double v : a) if (v > m) m = v;
        if (m == kNegInf) return kNegInf;
        const int n_ = int(a.size());
        ex.resize(n_);
        for (int i = 0; i < n_; ++i) ex[i] = a[i] - m;
        vmath::vexp(ex.data(), ex.data(), n_);
        double sum = 0.0;
        for (int i = 0; i < n_; ++i) sum += ex[i];
        return m + std::log(sum);
    }

    // log Σ w·e^f, the trapezoid weights applied to the exponentials (one log for the sum, none per node)
    double lse_weighted(const std::vector<double>& fv, const std::vector<double>& wv) {
        double m = kNegInf;
        for (double v : fv) if (v > m) m = v;
        if (m == kNegInf) return kNegInf;
        const int n_ = int(fv.size());
        ex.resize(n_);
        for (int i = 0; i < n_; ++i) ex[i] = fv[i] - m;
        vmath::vexp(ex.data(), ex.data(), n_);
        double sum = 0.0;
        for (int i = 0; i < n_; ++i) sum += wv[i] * ex[i];
        return m + std::log(sum);
    }

    // the Poisson pair's log terms at the nodes with rates `r` (n_ of them), vectorised: la = log(a2 + qq r),
    // lb = log(a2 + (1 - qq) r), each only when its column has reads
    void poisson_logs(const double* r, int n_, double qq) {
        argA.resize(n_); argB.resize(n_); la.resize(n_); lb.resize(n_);
        const double a2 = ctx_a2;
        for (int i = 0; i < n_; ++i) { argA[i] = a2 + qq * r[i]; argB[i] = a2 + (1.0 - qq) * r[i]; }
        if (s.u > 0) vmath::vlog(argA.data(), la.data(), n_);
        if (s.v > 0) vmath::vlog(argB.data(), lb.data(), n_);
    }

    // the held RNA levels' terms along an ascending run of nodes, the searches moving with the run
    double level_terms_from(double yy, const Assignment& a, int& hp, int& hn) const {
        double out = 0.0;
        if (s.pos.present) {
            if (a.fpos == 0.0) out += s.pos.vals[0];
            else out += interp_held_from(yy - a.shift_pos, s.pos.grid, s.pos.vals, s.pos.n, hp);
        }
        if (s.neg.present) {
            if (a.fneg == 0.0) out += s.neg.vals[0];
            else out += interp_held_from(yy - a.shift_neg, s.neg.grid, s.neg.vals, s.neg.n, hn);
        }
        return out;
    }

    // log f at the ascending nodes `yy` with rates `r` and composition values `comp_vals` under assignment `a`, into `out`
    void log_f_nodes(const double* yy, const double* r, const double* comp_vals, int n_, const Assignment& a, double* out) {
        poisson_logs(r, n_, a.qq);
        const double dna = ctx_dna;
        int hp = 0, hn = 0;
        for (int i = 0; i < n_; ++i) {
            double v = -dna - r[i] - lg + 0.5 * yy[i];
            if (s.u > 0) v += s.u * la[i];
            if (s.v > 0) v += s.v * lb[i];
            if (comp_live) v += comp_vals[i];
            v += level_terms_from(yy[i], a, hp, hn);
            out[i] = v;
        }
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

    // trapezoid weights from the node spacings (a duplicated node, from clipping, gets zero weight)
    void weights_of(const std::vector<double>& yy) {
        int m = int(yy.size());
        logw.assign(m, 0.0);
        if (m <= 1) { if (m == 1) logw[0] = 1.0; return; }
        logw[0] = 0.5 * (yy[1] - yy[0]);
        logw[m - 1] = 0.5 * (yy[m - 1] - yy[m - 2]);
        for (int i = 1; i < m - 1; ++i) logw[i] = 0.5 * ((yy[i] - yy[i - 1]) + (yy[i + 1] - yy[i]));
    }

    // the inner integral over y at the current density node for one strand assignment: the coarse nodes, the
    // window around the own term's mode, then the window around the peak found on those; every node evaluated once
    double inner(Assignment& a) {
        ready(a);
        const double dna = ctx_dna;
        const int m = int(coarse.size());
        fA.resize(m);
        poisson_logs(r_coarse.data(), m, a.qq);
        for (int i = 0; i < m; ++i) {
            double out = -dna - r_coarse[i] - lg + 0.5 * coarse[i];
            if (s.u > 0) out += s.u * la[i];
            if (s.v > 0) out += s.v * lb[i];
            if (comp_live) out += comp_row[i];
            out += a.lev[i];
            fA[i] = out;
        }
        if (!windowed) return lse_weighted(fA, logw_coarse);
        const int w = int(offs.size());
        fB.resize(w);
        log_f_nodes(yB.data(), rB.data(), compB.data(), w, a, fB.data());
        double best_y = m ? coarse[0] : yB[0], best_f = m ? fA[0] : fB[0];
        for (int i = 0; i < m; ++i) if (fA[i] > best_f) { best_f = fA[i]; best_y = coarse[i]; }
        for (int k = 0; k < w; ++k) if (fB[k] > best_f) { best_f = fB[k]; best_y = yB[k]; }
        yC.resize(w); fC.resize(w); rC.resize(w); compC.resize(w);
        for (int k = 0; k < w; ++k) yC[k] = std::min(std::max(best_y + offs[k], y_lo - 5.0), y_hi + 1.0);
        vmath::vexp(yC.data(), rC.data(), w);
        if (comp_live) {
            int hint = n_lam - 2;
            for (int k = 0; k < w; ++k) compC[k] = interp_held_from(ctx_log_dna - yC[k], lam, s.comp, n_lam, hint);
        } else {
            for (int k = 0; k < w; ++k) compC[k] = 0.0;
        }
        log_f_nodes(yC.data(), rC.data(), compC.data(), w, a, fC.data());
        merge3(coarse, fA, yB, fB, yC, fC, y, f);
        weights_of(y);
        return lse_weighted(f, logw);
    }

    void s_window(double centre, double sd_s, std::vector<double>& out) const {
        int k = int(std::ceil(R.window_sd));
        for (int i = -k; i <= k; ++i) {
            double sw = centre + i * sd_s;
            if (sw > 0.0 && sw < 1.0) out.push_back(sw);
        }
    }

    // the share lattice's node at k times the share's resolution (the per-node second pass's window), cached
    Assignment& share_node(int k, double sd_s) {
        if (int(share_cache.size()) <= k) { share_cache.resize(k + 1); share_ready.resize(k + 1, 0); }
        if (!share_ready[k]) { share_cache[k] = node_of(std::acos(std::sqrt(k * sd_s))); share_ready[k] = 1; }
        return share_cache[k];
    }

    // the share's base lattice from a list of angles: the endpoints added, sorted, duplicates removed
    std::vector<Assignment> lattice_of(std::vector<double> ph) const {
        ph.push_back(0.0); ph.push_back(kPi / 2);
        std::sort(ph.begin(), ph.end());
        ph.erase(std::unique(ph.begin(), ph.end()), ph.end());
        std::vector<Assignment> nodes;
        nodes.reserve(ph.size());
        for (double p : ph) nodes.push_back(node_of(p));
        return nodes;
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
        // single-strand, and both-strand with no strand reads and no RNA level: the share integrand is a
        // constant, the continuum integrates to that constant and so does each admitted atom, and the -log 3
        // cancels the count of three
        if (!s.both || (s.u == 0.0 && s.v == 0.0 && !s.pos.present && !s.neg.present)) {
            Assignment single = assignment(std::numeric_limits<double>::quiet_NaN(), s.q, 1.0, 0.0);
            for (int ix = 0; ix < nx; ++ix) {
                double dna = rho[ix] * s.eg;
                set_node(dna, std::log(dna));
                out[ix] = inner(single) + dna_term(rho[ix]);
            }
            return;
        }
        // both-strand: the arcsine continuum over the share plus the pure atoms the witness bits admit. The base
        // lattice (the midpoint nodes, the endpoints, the window around the observed share where the columns
        // resolve it) is shared by every density node; each density node then adds its own window at integer
        // multiples of the share's resolution around its best base node, so its log L depends on nothing but
        // itself.
        double sd_s = 1.0 / (std::sqrt(n + 1.0) * std::max(std::fabs(2 * s.q - 1), 1.0 / std::sqrt(n + 1.0)));
        int m = R.tilt_nodes;
        std::vector<double> base;
        bool located = n > 0 && std::fabs(2 * s.q - 1) * std::sqrt(n + 1) > 1.0;
        double s_star = located ? (s.u / n - (1 - s.q)) / (2 * s.q - 1) : 0.5;
        for (int i = 0; i < m; ++i) base.push_back((i + 0.5) * (kPi / 2) / m);
        if (located) {
            std::vector<double> win;
            s_window(s_star, sd_s, win);
            for (double sw : win) base.push_back(std::acos(std::sqrt(sw)));
        }
        std::vector<Assignment> nodes = lattice_of(base);
        const int np1 = int(nodes.size());
        Assignment atom_pos = assignment(std::numeric_limits<double>::quiet_NaN(), s.q, 1.0, 0.0);
        Assignment atom_neg = assignment(std::numeric_limits<double>::quiet_NaN(), 1.0 - s.q, 0.0, 1.0);
        std::vector<double> col(np1);
        std::vector<std::pair<double, int>> order;  // (angle, index: base j, or -(k + 1) for the share cache's k)
        std::vector<double> row;
        const int kw = int(std::ceil(R.window_sd));
        auto lg0 = [](double d) { return d > 0.0 ? std::log(d) : kNegInf; };
        for (int ix = 0; ix < nx; ++ix) {
            double dna = rho[ix] * s.eg;
            set_node(dna, std::log(dna));
            for (int j = 0; j < np1; ++j) col[j] = inner(nodes[j]);
            order.clear();
            for (int j = 0; j < np1; ++j) order.emplace_back(nodes[j].phi, j);
            int best = 0;
            for (int j = 1; j < np1; ++j) if (col[j] > col[best]) best = j;
            double c = std::cos(nodes[best].phi);
            int kc = int(std::round(c * c / sd_s));
            for (int i = -kw; i <= kw; ++i) {
                int k = kc + i;
                double sw = k * sd_s;
                if (!(sw > 0.0 && sw < 1.0)) continue;
                order.emplace_back(share_node(k, sd_s).phi, -(k + 1));
            }
            std::sort(order.begin(), order.end());
            order.erase(std::unique(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first == b.first; }), order.end());
            const int np_ = int(order.size());
            row.resize(np_);
            for (int j = 0; j < np_; ++j) {
                double lw = j == 0 ? lg0(0.5 * (order[1].first - order[0].first))
                          : j == np_ - 1 ? lg0(0.5 * (order[np_ - 1].first - order[np_ - 2].first))
                          : lg0(0.5 * ((order[j].first - order[j - 1].first) + (order[j + 1].first - order[j].first)));
                int id = order[j].second;
                row[j] = (id >= 0 ? col[id] : inner(share_node(-id - 1, sd_s))) + lw;
            }
            hyps.clear();
            hyps.push_back(lse(row) + std::log(2.0 / kPi));
            if (!s.neg.present) hyps.push_back(inner(atom_pos));
            if (!s.pos.present) hyps.push_back(inner(atom_neg));
            out[ix] = lse(hyps) - std::log(3.0) + dna_term(rho[ix]);
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
