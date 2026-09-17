// pass_kernel.cpp — ONE DIRECTIONAL PASS of the composition-transfer policy, as one native call.
//
// The backbone (`calibration/sweep.py::_pass`) runs a pass in chain order: for every destination with
// a neighbour on the pass's side that is not a terminal, the recipient receives what its neighbour
// sends. What a hop does is the transfer policy's (`calibration/messages/transfer.py`,
// `faces.py`, `lanes.py`, `transfer_rows.py`): the composition rule at the face — FORWARD, TRANSPORT,
// SPLICE_OUT, EDGE, LEVEL, or none — applied to the sender's own claim and what it holds from its far
// side, and each level lane whose faces include this one carrying its population's level (the sender
// emits; a full recipient prices the hop and takes the lower side). This file is that arithmetic on
// the same tables the Python reads (`Faces`, `LevelLane`, `Received`), row by row over the solve grid,
// with the per-hop Python call removed. The executable specification is the Python; the gates are
// `tests/calibration/test_pass_kernel.py` (both passes on real captured blocks, every table compared)
// and `profiling/sweep_replay.py replay --tolerance` on a captured sweep.
//
// Numerics: every row operation is the Python one term for term (max-normalisation, the edge-padded
// Gaussian blur at half-width ceil(4 sqrt(v) / dlam), linear interpolation with numpy's end rules,
// trigamma for the counting variance). Summation orders differ from numpy's, so a port is held to
// the replay's derived budget, not to bits.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

namespace {

using Vec = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using Mat = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using BoolVec = nb::ndarray<bool, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using BoolMat = nb::ndarray<bool, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using IdxVec = nb::ndarray<int64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using IdxMat = nb::ndarray<int32_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using KindMat = nb::ndarray<int8_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

constexpr double EPS = 1.0e-9;  // transfer_rows.EPS
constexpr double TINY = 1.0e-300;
constexpr int NONE = 0, FORWARD = 1, TRANSPORT = 2, SPLICE_OUT = 3, EDGE = 4, LEVEL = 5;

// ---- scalar pieces ---------------------------------------------------------------------------------

// trigamma(x) = zeta(2, x): the recurrence up to x >= 12, then the asymptotic series through
// B_14 / x^15, whose truncation there is below 1e-17 relative.
double trigamma(double x) {
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
double hop_price(double n_s, double a_s, double n_x, double a_x) {
    double v = count_logvar(n_s) + count_logvar(n_x);
    if (n_s > 0.0 && n_x > 0.0) {
        const double l = std::log((n_x / a_x) / (n_s / a_s));
        v += std::max(0.0, l * l - (1.0 / n_s + 1.0 / n_x));
    }
    return v;
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
void interp(const double* x, int n, const double* xp, const double* fp, int K, double left,
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
void blur_row(const double* row, int K, const double* lam, double v, double* out,
              std::vector<double>& scratch) {
    const double m = vmax(row, K);
    for (int j = 0; j < K; ++j) out[j] = row[j] - m;
    if (v > 0.0 && K > 1) {
        const double dlam = lam[1] - lam[0];
        const int half = std::max(static_cast<int>(std::ceil(4.0 * std::sqrt(v) / dlam)), 1);
        const int W = 2 * half + 1;
        scratch.resize(static_cast<size_t>(W) + static_cast<size_t>(K + 2 * half) + static_cast<size_t>(K));
        double* kern = scratch.data();
        double* edged = kern + W;
        double* pr = edged + (K + 2 * half);
        double ksum = 0.0;
        for (int t = 0; t < W; ++t) {
            const double x = (t - half) * dlam;
            kern[t] = std::exp(-0.5 * x * x / v);
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
            pr[j] = s;
        }
        for (int j = 0; j < K; ++j) out[j] = std::log(std::max(pr[j], TINY));
    }
    norm_inplace(out, K);
}

inline void lower_side(const double* p, int K, double* out) {
    double m = -std::numeric_limits<double>::infinity();
    for (int j = 0; j < K; ++j) {
        if (p[j] > m) m = p[j];
        out[j] = m;
    }
    norm_inplace(out, K);
}

// face_map_lambda(lam, n_u, a_g_b, a_r_b, e_g_e, e_r_e, s)
void face_map_lambda(const double* lam, int K, double n_u, double a_g_b, double a_r_b, double e_g_e,
                     double e_r_e, double s, double* out) {
    for (int j = 0; j < K; ++j) {
        const double sig = 1.0 / (1.0 + std::exp(-lam[j]));
        const double g_arm = n_u * sig / a_g_b * e_g_e;
        const double r_arm = (n_u * (1.0 - sig) / a_r_b + s) * e_r_e;
        out[j] = std::log(std::max(g_arm, TINY)) - std::log(std::max(r_arm, TINY));
    }
}

struct Scratch {
    std::vector<double> a, b, c, d, e, f, blur;
    explicit Scratch(int K) : a(K), b(K), c(K), d(K), e(K), f(K) {}
};

// transport_row(row, lam, lam_e_of_u, n_u, n_s)
void transport_row(const double* row, const double* lam, int K, const double* map, double n_u,
                   double n_s, Scratch& S, double* out) {
    const double m = vmax(row, K);
    for (int j = 0; j < K; ++j) S.a[j] = row[j] - m;
    interp(lam, K, map, lam, K, lam[0], lam[K - 1], S.b.data());  // lam_u_of_x
    interp(S.b.data(), K, lam, S.a.data(), K, S.a[0], S.a[K - 1], S.c.data());
    blur_row(S.c.data(), K, lam, trigamma(n_u + 0.5) + trigamma(n_s + 0.5), out, S.blur);
}

// splice_out_row(row_e, lam, n_u, n_s, a_g_b, a_g_e)
void splice_out_row(const double* row_e, const double* lam, int K, double n_u, double n_s,
                    double a_g_b, double a_g_e, const double* nodes, int n_nodes, Scratch& S,
                    double* out) {
    if (!(n_u > 0.0 && a_g_b > 0.0 && a_g_e > 0.0) || ptp(row_e, K) <= EPS) {
        std::fill(out, out + K, 0.0);
        return;
    }
    const double m = vmax(row_e, K);
    for (int j = 0; j < K; ++j) S.a[j] = row_e[j] - m;
    const double sd = std::sqrt(trigamma(n_s + 0.5) + trigamma(n_u + 0.5));
    std::fill(S.d.begin(), S.d.end(), 0.0);
    for (int t = 0; t < n_nodes; ++t) {
        const double s = n_s / a_g_b * std::exp(nodes[t] * sd);
        face_map_lambda(lam, K, n_u, a_g_b, a_g_b, a_g_e, a_g_e, s, S.b.data());
        interp(S.b.data(), K, lam, S.a.data(), K, S.a[0], S.a[K - 1], S.c.data());
        for (int j = 0; j < K; ++j) S.d[j] += std::exp(S.c[j]);
    }
    for (int j = 0; j < K; ++j) out[j] = std::log(std::max(S.d[j] / n_nodes, TINY));
    norm_inplace(out, K);
    if (ptp(out, K) <= EPS) std::fill(out, out + K, 0.0);
}

// level_row(row_b, lam, lam_i_of_b, v)
void level_row(const double* row_b, const double* lam, int K, const double* map, double v, Scratch& S,
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

// ---- the tables ---------------------------------------------------------------------------------------

struct FacesView {
    KindMat kind; IdxMat row; IdxMat row2;
    Mat n_u, n_s, a_b, a_x, width, var;
    Mat rows;
};

struct LevelsView {  // one lane's Received table, written in place
    BoolVec present; Mat profile; Vec count; Vec opportunity; BoolVec has_witness; Vec rna_count; Vec rna_count_var;
};

struct LaneView {
    int field;  // 0 gdna, 1 rna_pos, 2 rna_neg — which Levels table of the Received
    BoolMat face; BoolMat two_sided; BoolVec empty;
    Mat own_level; BoolVec own_mask;
    Vec count; Vec a;
    bool has_other; Vec other;
    Mat flux_witness; BoolVec flux_mask;  // (n, 2): count, opportunity
};

inline const double* row_of(const Mat& m, int i) { return m.data() + static_cast<size_t>(i) * m.shape(1); }
inline double* row_of(Mat& m, int i) { return m.data() + static_cast<size_t>(i) * m.shape(1); }
inline bool at2(const BoolMat& m, int i, int side) { return m.data()[static_cast<size_t>(i) * 2 + side]; }

inline void levels_write(LevelsView& L, int K, int i, const double* profile, double count,
                         double opportunity, bool witness, double rna_count, double rna_var) {
    L.present.data()[i] = true;
    double* dest = row_of(L.profile, i);
    if (dest != profile) std::copy(profile, profile + K, dest);
    L.count.data()[i] = count;
    L.opportunity.data()[i] = opportunity;
    L.has_witness.data()[i] = witness;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    L.rna_count.data()[i] = witness ? rna_count : nan;
    L.rna_count_var.data()[i] = witness ? rna_var : nan;
}

inline void levels_forward(LevelsView& L, int K, int s, int i) {
    L.present.data()[i] = L.present.data()[s];
    std::copy(row_of(L.profile, s), row_of(L.profile, s) + K, row_of(L.profile, i));
    L.count.data()[i] = L.count.data()[s];
    L.opportunity.data()[i] = L.opportunity.data()[s];
    L.has_witness.data()[i] = L.has_witness.data()[s];
    L.rna_count.data()[i] = L.rna_count.data()[s];
    L.rna_count_var.data()[i] = L.rna_count_var.data()[s];
}

// intersect the present parts (pointwise minimum), max-normalised, into out; whether any part
bool intersect2(const double* p, const double* q, int K, double* out) {
    if (p == nullptr && q == nullptr) return false;
    if (p != nullptr && q != nullptr) {
        for (int j = 0; j < K; ++j) out[j] = std::min(p[j], q[j]);
    } else {
        const double* r = p != nullptr ? p : q;
        std::copy(r, r + K, out);
    }
    norm_inplace(out, K);
    return true;
}

// LevelLane.emit(s, x, levels): what s sends toward x, written into row x. Whether anything was sent.
bool lane_emit(const LaneView& ln, LevelsView& L, int K, int s, int x, int side, Scratch& S) {
    const double* own = ln.own_mask.data()[s] ? row_of(ln.own_level, s) : nullptr;
    const bool far = L.present.data()[s];
    if (ln.empty.data()[s] && own == nullptr) {
        if (far) levels_forward(L, K, s, x);
        return far;
    }
    const double* own_used = own;
    if (own != nullptr && !at2(ln.two_sided, x, side)) {
        lower_side(own, K, S.a.data());
        own_used = S.a.data();
    }
    const double* held = far ? row_of(L.profile, s) : nullptr;
    if (!intersect2(own_used, held, K, S.b.data())) return false;
    if (ln.empty.data()[s]) {
        if (!ln.flux_mask.data()[s]) throw std::runtime_error("an empty node with an own level has no flux witness");
        const double* w = row_of(ln.flux_witness, s);
        levels_write(L, K, x, S.b.data(), w[0], w[1], false, 0.0, 0.0);
        return true;
    }
    if (ln.has_other) {
        const double c_r = ln.count.data()[s], c_o = ln.other.data()[s];
        levels_write(L, K, x, S.b.data(), c_r, ln.a.data()[s], true, c_r - c_o, c_r + c_o);
    } else {
        levels_write(L, K, x, S.b.data(), ln.count.data()[s], ln.a.data()[s], false, 0.0, 0.0);
    }
    return true;
}

// LevelLane.receive(levels, s, x): a FULL recipient re-prices row x in place.
void lane_receive(const LaneView& ln, LevelsView& L, int K, const double* lam, int x, int side,
                  Scratch& S) {
    const double n_sent = L.count.data()[x], a_sent = L.opportunity.data()[x];
    const double cnt = ln.count.data()[x], a_x = ln.a.data()[x];
    double v;
    if (!ln.has_other || !L.has_witness.data()[x]) {
        v = hop_price(n_sent, a_sent, cnt, a_x);
    } else {
        v = count_logvar(n_sent) + count_logvar(cnt);
        const double n_s = L.rna_count.data()[x], v_s = L.rna_count_var.data()[x];
        const double n_x = cnt - ln.other.data()[x], v_x = cnt + ln.other.data()[x];
        if (n_s > 0.0 && n_x > 0.0) {
            const double l = std::log((n_x / a_x) / (n_s / a_sent));
            v += std::max(0.0, l * l - (v_s / (n_s * n_s) + v_x / (n_x * n_x)));
        }
    }
    const double* held = row_of(L.profile, x);
    const double* p = held;
    if (!at2(ln.two_sided, x, side)) {
        lower_side(held, K, S.a.data());
        p = S.a.data();
    }
    if (v > 0.0) {
        blur_row(p, K, lam, v, S.b.data(), S.blur);
        p = S.b.data();
    }
    if (ln.has_other) {
        const double c_o = ln.other.data()[x];
        levels_write(L, K, x, p, cnt, a_x, true, cnt - c_o, cnt + c_o);
    } else {
        levels_write(L, K, x, p, cnt, a_x, false, 0.0, 0.0);
    }
}

// Faces.apply(s, i, own, held): the composition rule at the face into i from side; the row for i in
// out, or nullptr for no claim.
const double* faces_apply(const FacesView& F, int i, int side, const double* own, const double* held,
                          const double* lam, int K, const double* nodes, int n_nodes, Scratch& S,
                          double* out) {
    const size_t at = static_cast<size_t>(i) * 2 + side;
    const int k = F.kind.data()[at];
    if (k == FORWARD) {
        if (own == nullptr && held == nullptr) return nullptr;
        for (int j = 0; j < K; ++j) out[j] = (own ? own[j] : 0.0) + (held ? held[j] : 0.0);
        norm_inplace(out, K);
        return out;
    }
    if (k == EDGE) {
        const int r = F.row.data()[at];
        if (r < 0) return nullptr;
        const double* row = row_of(F.rows, r);
        if (ptp(row, K) <= EPS) return nullptr;
        std::copy(row, row + K, out);
        return out;
    }
    if (k == LEVEL) {
        if (own == nullptr) {
            const int r2 = F.row2.data()[at];
            if (r2 < 0) return nullptr;
            std::copy(row_of(F.rows, r2), row_of(F.rows, r2) + K, out);
            return out;
        }
        level_row(own, lam, K, row_of(F.rows, F.row.data()[at]), F.var.data()[at], S, out);
        return out;
    }
    if (k != TRANSPORT && k != SPLICE_OUT) return nullptr;
    if (own == nullptr && held == nullptr) return nullptr;
    double* sending = S.e.data();
    for (int j = 0; j < K; ++j) sending[j] = (own ? own[j] : 0.0) + (held ? held[j] : 0.0);
    norm_inplace(sending, K);
    const double n_u = F.n_u.data()[at], n_s = F.n_s.data()[at];
    double* tmp = S.f.data();
    if (k == TRANSPORT) {
        transport_row(sending, lam, K, row_of(F.rows, F.row.data()[at]), n_u, n_s, S, tmp);
    } else {
        splice_out_row(sending, lam, K, n_u, n_s, F.a_b.data()[at], F.a_x.data()[at], nodes, n_nodes, S, tmp);
    }
    const double w = F.width.data()[at];
    if (w > 0.0) {
        blur_row(tmp, K, lam, w, out, S.blur);
    } else {
        std::copy(tmp, tmp + K, out);
    }
    return out;
}

// ---- the pass -----------------------------------------------------------------------------------------

void transfer_pass(Vec lam, IdxVec seq, IdxVec nbr, BoolVec terminal, Mat own, BoolVec own_mask, KindMat f_kind, IdxMat f_row,
                   IdxMat f_row2, Mat f_n_u, Mat f_n_s, Mat f_a_b, Mat f_a_x, Mat f_width, Mat f_var, Mat f_rows,
                   Vec marginal_nodes, nb::list lanes, BoolVec has_composition, Mat composition, nb::list levels) {
    const int K = static_cast<int>(lam.shape(0));
    const int n = static_cast<int>(seq.shape(0));
    if (static_cast<int>(composition.shape(1)) != K || static_cast<int>(own.shape(1)) != K)
        throw std::invalid_argument("transfer_pass: the rows and the grid disagree");
    FacesView F{f_kind, f_row, f_row2, f_n_u, f_n_s, f_a_b, f_a_x, f_width, f_var, f_rows};
    std::vector<LevelsView> tables;
    for (nb::handle h : levels) {
        nb::tuple t = nb::cast<nb::tuple>(h);
        tables.push_back(LevelsView{nb::cast<BoolVec>(t[0]), nb::cast<Mat>(t[1]), nb::cast<Vec>(t[2]),
                                    nb::cast<Vec>(t[3]), nb::cast<BoolVec>(t[4]), nb::cast<Vec>(t[5]),
                                    nb::cast<Vec>(t[6])});
    }
    std::vector<LaneView> lns;
    for (nb::handle h : lanes) {
        nb::tuple t = nb::cast<nb::tuple>(h);
        LaneView ln{nb::cast<int>(t[0]), nb::cast<BoolMat>(t[1]), nb::cast<BoolMat>(t[2]), nb::cast<BoolVec>(t[3]),
                    nb::cast<Mat>(t[4]), nb::cast<BoolVec>(t[5]), nb::cast<Vec>(t[6]), nb::cast<Vec>(t[7]),
                    !t[8].is_none(), Vec(), nb::cast<Mat>(t[9]), nb::cast<BoolVec>(t[10])};
        if (ln.has_other) ln.other = nb::cast<Vec>(t[8]);
        if (ln.field < 0 || ln.field >= static_cast<int>(tables.size()))
            throw std::invalid_argument("transfer_pass: a lane names a table the Received does not hold");
        lns.push_back(std::move(ln));
    }
    const double* L = lam.data();
    const double* nodes = marginal_nodes.data();
    const int n_nodes = static_cast<int>(marginal_nodes.shape(0));
    Scratch S(K);
    std::vector<double> out(K);
    const int64_t* seqp = seq.data();
    const int64_t* nbrp = nbr.data();
    const bool* term = terminal.data();
    bool* has_comp = has_composition.data();
    for (int idx = 0; idx < n; ++idx) {
        const int i = static_cast<int>(seqp[idx]);
        const int s = static_cast<int>(nbrp[i]);
        if (s < 0 || term[i]) continue;
        const int side = s < i ? 0 : 1;
        if (F.kind.data()[static_cast<size_t>(i) * 2 + side] != NONE) {
            const double* own_s = own_mask.data()[s] ? row_of(own, s) : nullptr;
            const double* held = has_comp[s] ? row_of(composition, s) : nullptr;
            const double* r = faces_apply(F, i, side, own_s, held, L, K, nodes, n_nodes, S, out.data());
            if (r != nullptr && ptp(r, K) > EPS) {
                double* dest = row_of(composition, i);
                std::copy(r, r + K, dest);
                norm_inplace(dest, K);
                has_comp[i] = true;
            }
        }
        for (const LaneView& ln : lns) {
            if (!at2(ln.face, i, side)) continue;
            LevelsView& T = tables[ln.field];
            if (lane_emit(ln, T, K, s, i, side, S) && !ln.empty.data()[i]) lane_receive(ln, T, K, L, i, side, S);
        }
    }
}

}  // namespace

NB_MODULE(_pass_impl, m) {
    m.doc() = "One directional pass of the composition-transfer policy, on the backbone's tables.";
    m.def("transfer_pass", &transfer_pass, nb::arg("lam"), nb::arg("seq"), nb::arg("nbr"),
          nb::arg("terminal"), nb::arg("own"), nb::arg("own_mask"), nb::arg("f_kind"), nb::arg("f_row"),
          nb::arg("f_row2"), nb::arg("f_n_u"), nb::arg("f_n_s"), nb::arg("f_a_b"), nb::arg("f_a_x"),
          nb::arg("f_width"), nb::arg("f_var"), nb::arg("f_rows"), nb::arg("marginal_nodes"),
          nb::arg("lanes"), nb::arg("has_composition"), nb::arg("composition"), nb::arg("levels"),
          "Run one pass in chain order: for every destination in `seq` with a neighbour `nbr[i] >= 0` "
          "that is not a terminal, apply the face's composition rule and carry each lane's level, "
          "writing the Received tables in place.");
    m.def("trigamma", &trigamma, nb::arg("x"), "zeta(2, x), the counting variance's one home in C++.");
}
