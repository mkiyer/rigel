// prepare_kernel.cpp — the BUILDERS of the composition-transfer policy for one block, as one native call.
//
// `TransferPolicy.prepare` (`calibration/messages/transfer.py`) builds, per block, every node's own claim
// (`_claims`), the recipient's rule per directed face (`_splice_faces`, `_edge_level`, `_terminus_rules`,
// `_alternative_splice_site`, written into the `Faces` tables) and the level lanes (`lanes.gdna_lane`,
// `lanes.rna_lanes`: each lane's faces, own levels, junction flux levels and flux witnesses). This file is
// that arithmetic on the same tables, written in place, with the per-node Python removed: the Python
// builders are the executable specification and the reference the gate compares against
// (`tests/calibration/test_prepare_kernel.py`: every table of both on the toy's captured sweep; the wiring by
// a spy), with `profiling/sweep_replay.py replay --tolerance` on a captured sweep as the second verdict. The
// row constructors are `transfer_rows.h`, shared with the pass kernel.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "transfer_rows.h"

namespace nb = nanobind;
using namespace transfer_rows;

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

// ---- the boundary flags (calibration/splice_graph.py) ----------------------------------------------

constexpr int TSS_POS = 1 << 0, TSS_NEG = 1 << 1, TES_POS = 1 << 2, TES_NEG = 1 << 3;
constexpr int DON_POS = 1 << 4, DON_NEG = 1 << 5, ACC_POS = 1 << 6, ACC_NEG = 1 << 7;
constexpr int TERMINUS = TSS_POS | TSS_NEG | TES_POS | TES_NEG;
constexpr int SJ_FLAGS = DON_POS | DON_NEG | ACC_POS | ACC_NEG;
//: a transcript body that extends genomic-RIGHT from its terminus leaves the OUTSIDE flank on the left
constexpr int BODY_RIGHT = TSS_POS | TES_NEG, BODY_LEFT = TES_POS | TSS_NEG;
//: a DONOR bit marks the intron's LOW end (the intron lies right), an ACCEPTOR its HIGH end
constexpr int INTRON_RIGHT = DON_POS | DON_NEG, INTRON_LEFT = ACC_POS | ACC_NEG;
//: a strand's four boundary bits and its terminus bits (transfer_rows.strand_bits), by strand 0 = +, 1 = −
constexpr int ALL_BITS[2] = {TSS_POS | TES_POS | DON_POS | ACC_POS, TSS_NEG | TES_NEG | DON_NEG | ACC_NEG};
constexpr int TERM_BITS[2] = {TSS_POS | TES_POS, TSS_NEG | TES_NEG};

// face_is_licensed: no terminus on the face and both flanks admit the same strand set
inline bool face_is_licensed(int f, bool fp_e, bool fn_e, bool fp_i, bool fn_i) {
    return !(f & TERMINUS) && fp_e == fp_i && fn_e == fn_i;
}
// boundary_shares_strand: the same strand set on both, and that set a single strand
inline bool boundary_shares_strand(bool fp_b, bool fn_b, bool fp_i, bool fn_i) {
    return fp_b == fp_i && fn_b == fn_i && fp_b != fn_b;
}
// outside_flank: (outside, inside) of a terminus boundary, or (-1, -1)
inline void outside_flank(int f, int64_t left, int64_t right, int64_t& o, int64_t& i) {
    o = i = -1;
    if (!(f & TERMINUS)) return;
    const bool to_right = f & BODY_RIGHT, to_left = f & BODY_LEFT;
    if (to_right && !to_left) { o = left; i = right; }
    else if (to_left && !to_right) { o = right; i = left; }
}
// junction_exon_side: the flank on a junction's exon side, or -1
inline int64_t junction_exon_side(int f, int64_t left, int64_t right) {
    const bool don = f & (DON_POS | DON_NEG), acc = f & (ACC_POS | ACC_NEG);
    if (don && !acc) return left;
    if (acc && !don) return right;
    return -1;
}
// junction_flanks: (C, E) at an exon|exon junction with no terminus, or (-1, -1)
inline void junction_flanks(int f, int64_t left, int64_t right, int64_t& c, int64_t& e) {
    c = e = -1;
    if ((f & TERMINUS) || !(f & SJ_FLAGS)) return;
    const bool to_right = f & INTRON_RIGHT, to_left = f & INTRON_LEFT;
    if (to_right && !to_left) { c = right; e = left; }
    else if (to_left && !to_right) { c = left; e = right; }
}
// read_column: the genome-strand column strand `col`'s RNA reads on
inline int read_column(int col, bool has_strand, double kappa) {
    return (!has_strand || kappa >= 0.5) ? col : 1 - col;
}

// ---- the chain as the builders read it (transfer._Chain) --------------------------------------------

struct Chain {
    int n, K;
    const double* lam;
    const bool *is_bnd, *is_exon, *fp, *fn, *exon_pos, *exon_neg, *has_own;
    const int64_t *left, *right;
    const uint16_t* flags;
    const double *n_u, *n_s, *a_g, *a_r, *belief, *cnt, *src;
    const double *route_lo, *route_hi, *sj_lo, *sj_hi;  // (n, 2) by transcript strand
    const double* flux;                                 // the sj count over both strands
    std::vector<char> is_intron, is_intergenic;
    bool has_strand;
    double kappa, od_g, od_r;

    bool intron(int i) const { return is_intron[i]; }
    bool intergenic(int i) const { return is_intergenic[i]; }
    int64_t other_flank(int64_t b, int64_t e) const { return right[b] == e ? left[b] : right[b]; }
    // (exon, intron) across boundary b, or (-1, -1)
    void intron_exon_pair(int64_t b, int64_t& e, int64_t& i) const {
        e = i = -1;
        const int64_t lo = left[b], hi = right[b];
        if (lo < 0 || hi < 0) return;
        if (is_exon[lo] && intron(hi)) { e = lo; i = hi; }
        else if (is_exon[hi] && intron(lo)) { e = hi; i = lo; }
    }
    // a slot's own strand profile, the variance frozen at its incoming belief
    void strand_profile(int64_t x, double* out) const {
        const double f_ref = std::isfinite(belief[x]) ? belief[x] : 0.5;
        strand_row(lam, K, cnt[2 * x], cnt[2 * x + 1], fp[x], kappa, od_g, od_r, f_ref, out);
    }
    // a node's own strand MODE and its log-odds variance; false at a vertex mode
    bool strand_mode(int64_t y, double& f, double& v) const {
        const double c0 = cnt[2 * y], nn = c0 + cnt[2 * y + 1];
        const double p = c0 / nn;
        const double ks = fp[y] ? kappa : 1.0 - kappa;
        f = (p - ks) / (0.5 - ks);
        if (!(0.0 < f && f < 1.0)) return false;
        const double v_log = p * (1.0 - p) / nn / ((p - ks) * (p - ks));
        v = v_log / ((1.0 - f) * (1.0 - f));
        return true;
    }
};

// ---- the tables the builders write ----------------------------------------------------------------

struct RowsOut {  // faces.RowTable: rows (n, K) and mask (n,)
    double* rows; bool* mask; int K;
    void set(int64_t i, const double* row) { std::copy(row, row + K, rows + i * K); mask[i] = true; }
    void set_zero(int64_t i) { std::fill(rows + i * K, rows + (i + 1) * K, 0.0); mask[i] = true; }
    const double* get(int64_t i) const { return mask[i] ? rows + i * K : nullptr; }
};

struct FacesOut {  // faces.Faces
    const Chain* c;
    int8_t* kind; int32_t* row; int32_t* row2;
    double *n_u, *n_s, *a_b, *a_x, *width, *var;
    double* rows; int cap; int n_rows = 0;

    int keep(const double* r) {
        if (r == nullptr) return -1;
        if (n_rows == cap) throw std::runtime_error("transfer_prepare: the face row store is full");
        std::copy(r, r + c->K, rows + static_cast<size_t>(n_rows) * c->K);
        return n_rows++;
    }
    void set(int64_t s, int64_t i, int k, const double* r = nullptr, const double* r2 = nullptr, double nu = 0.0,
             double ns = 0.0, double ab = 0.0, double ax = 0.0, double w = 0.0, double v = 0.0) {
        const int side = s < i ? 0 : 1;
        const int64_t nbr = side == 0 ? c->left[i] : c->right[i];
        if (nbr != s)
            throw std::runtime_error("transfer_prepare: no face into " + std::to_string(i) + " from " +
                                     std::to_string(s));
        const size_t at = static_cast<size_t>(i) * 2 + side;
        if (kind[at] != NONE)
            throw std::runtime_error("transfer_prepare: the face into " + std::to_string(i) + " from " +
                                     std::to_string(s) + " already carries a rule");
        kind[at] = static_cast<int8_t>(k);
        row[at] = keep(r);
        row2[at] = keep(r2);
        n_u[at] = nu; n_s[at] = ns; a_b[at] = ab; a_x[at] = ax; width[at] = w; var[at] = v;
    }
};

struct LaneOut {  // lanes.LevelLane's tables: faces, two-sided faces, own levels, flux levels, witnesses
    bool* face; bool* two_sided; RowsOut own_level;
    double* flux_rows; bool* flux_mask;      // (n, 2, K), (n, 2)
    double* witness; bool* witness_mask;     // (n, 2), (n,)
};

// ---- the builders — one per shipped message (transfer.py, lanes.py) ---------------------------------

void claims(const Chain& c, RowsOut& own, Scratch& S) {
    for (int i = 0; i < c.n; ++i) {
        if (!c.intron(i)) continue;
        const double* r = c.src + static_cast<size_t>(i) * c.K;
        if (ptp(r, c.K) > EPS) {
            const double m = vmax(r, c.K);
            for (int j = 0; j < c.K; ++j) S.a[j] = r[j] - m;
            own.set(i, S.a.data());
        }
    }
    if (!c.has_strand) return;
    for (int x = 0; x < c.n; ++x) {
        const bool single = c.fp[x] != c.fn[x];
        if (!(single && c.has_own[x] && (c.is_exon[x] || c.is_bnd[x]))) continue;
        c.strand_profile(x, S.a.data());
        own.set(x, S.a.data());
    }
}

void splice_faces(const Chain& c, FacesOut& F, Scratch& S) {
    for (int64_t b = 0; b < c.n; ++b) {
        if (!c.is_bnd[b]) continue;
        int64_t e, i;
        c.intron_exon_pair(b, e, i);
        if (e < 0) continue;
        F.set(i, b, FORWARD);
        if (boundary_shares_strand(c.fp[b], c.fn[b], c.fp[i], c.fn[i])) F.set(b, i, FORWARD);
        const int hi = c.left[e] == b ? 1 : 0;
        if (!face_is_licensed(c.flags[b], c.fp[e], c.fn[e], c.fp[i], c.fn[i])) continue;
        if (!(c.n_u[b] > 0 && c.a_g[b] > 0 && c.a_r[b] > 0 && c.a_g[e] > 0 && c.a_r[e] > 0)) continue;
        const double* rr = hi ? c.route_hi : c.route_lo;
        const double* sc = hi ? c.sj_hi : c.sj_lo;
        const double rate = rr[2 * b] + rr[2 * b + 1], s = sc[2 * b] + sc[2 * b + 1];
        face_map_lambda(c.lam, c.K, c.n_u[b], c.a_g[b], c.a_r[b], c.a_g[e], c.a_r[e], rate, S.a.data());
        F.set(b, e, TRANSPORT, S.a.data(), nullptr, c.n_u[b], s);
        F.set(e, b, SPLICE_OUT, nullptr, nullptr, c.n_u[b], s, c.a_g[b], c.a_g[e]);
    }
}

void edge_level(const Chain& c, RowsOut& own, FacesOut& F, Scratch& S) {
    for (int64_t e = 0; e < c.n; ++e) {
        if (!c.is_exon[e]) continue;
        for (const int64_t b : {c.left[e], c.right[e]}) {
            if (b < 0 || !c.is_bnd[b]) continue;
            const int64_t o = c.other_flank(b, e);
            if (!(o >= 0 && c.intergenic(o) && c.a_g[b] > 0 && c.a_g[e] > 0 && c.n_u[e] > 0)) continue;
            own.set_zero(b);  // the level claim (the rule reads the count, not the row)
            edge_level_row(c.lam, c.K, c.n_u[b], c.n_u[e], c.a_g[b], c.a_g[e], S.a.data());
            F.set(b, e, EDGE, S.a.data());
        }
    }
}

void terminus_rules(const Chain& c, FacesOut& F, Scratch& S) {
    for (int64_t b = 0; b < c.n; ++b) {
        if (!(c.is_bnd[b] && c.left[b] >= 0 && c.right[b] >= 0)) continue;
        const int64_t lo = c.left[b], hi = c.right[b];
        int64_t o, i;
        outside_flank(c.flags[b], lo, hi, o, i);
        if (o < 0) continue;
        const int64_t ex_side = (c.flags[b] & SJ_FLAGS) ? junction_exon_side(c.flags[b], lo, hi) : -1;
        const double flux_b = ex_side >= 0 ? c.flux[b] : 0.0;
        if (c.is_exon[lo] && c.is_exon[hi] && boundary_shares_strand(c.fp[b], c.fn[b], c.fp[o], c.fn[o]) &&
            c.n_u[b] > 0 && c.a_g[b] > 0 && c.a_g[o] > 0) {
            const double s_out = c.n_s[b] + (ex_side == o ? flux_b : 0.0);
            face_map_lambda(c.lam, c.K, c.n_u[b], c.a_g[b], c.a_g[b], c.a_g[o], c.a_g[o], s_out / c.a_g[b],
                            S.a.data());
            F.set(o, b, SPLICE_OUT, nullptr, nullptr, c.n_u[b], s_out, c.a_g[b], c.a_g[o]);
            F.set(b, o, TRANSPORT, S.a.data(), nullptr, c.n_u[b], s_out);
        }
        // THE LEVEL RULE: the boundary's gDNA level into the inside region
        if (!(c.is_exon[i] && (c.is_exon[o] || c.intron(o)))) continue;
        if (!(c.n_u[b] > 0 && c.n_u[i] > 0 && c.a_g[b] > 0 && c.a_g[i] > 0)) continue;
        if (!boundary_shares_strand(c.fp[b], c.fn[b], c.fp[i], c.fn[i])) continue;
        const double density_b = c.n_u[b] / c.a_g[b];
        const double total_b = c.n_u[b] + c.n_s[b];
        const double total_i = total_b + (ex_side == i ? flux_b : 0.0);
        const double r = (c.n_u[i] / c.a_g[i]) / (total_i / c.a_g[b]);
        const double lr = std::log(r);
        double v_pair = std::max(0.0, lr * lr - (1.0 / c.n_u[i] + 1.0 / total_i));
        const double r_mode = (c.n_u[i] / c.a_g[i]) / (total_b / c.a_g[b]);
        if (c.has_strand && c.has_own[b] && c.has_own[i]) {
            double f_b, v_b, f_i, v_i;
            if (c.strand_mode(b, f_b, v_b) && c.strand_mode(i, f_i, v_i)) {
                const double f_pred = std::min(f_b / r_mode, 1.0 - 1e-9);
                const double d = std::log(f_i / (1.0 - f_i)) - std::log(f_pred / (1.0 - f_pred));
                v_pair += std::max(0.0, d * d - (v_b + v_i + 1.0 / c.n_u[i] + 1.0 / total_b));
            }
        }
        const double v_level = (count_logvar(c.n_u[b]) + count_logvar(c.n_u[i])) + v_pair;
        level_map_lambda(c.lam, c.K, density_b, c.a_g[i], c.n_u[i], S.a.data());
        level_bound_row(c.lam, c.K, density_b, c.a_g[i], c.n_u[i], v_level, S.b.data());
        F.set(b, i, LEVEL, S.a.data(), S.b.data(), 0.0, 0.0, 0.0, 0.0, 0.0, v_level);
    }
}

void alternative_splice_site(const Chain& c, FacesOut& F, Scratch& S) {
    if (!c.has_strand) return;
    for (int64_t b = 0; b < c.n; ++b) {
        if (!(c.is_bnd[b] && c.left[b] >= 0 && c.right[b] >= 0)) continue;
        const int64_t lo = c.left[b], hi = c.right[b];
        if (!(c.is_exon[lo] && c.is_exon[hi])) continue;
        int64_t c_side, e_side;
        junction_flanks(c.flags[b], lo, hi, c_side, e_side);
        if (c_side < 0 || !(c.n_u[b] > 0 && c.a_g[b] > 0)) continue;
        const int64_t flank[2] = {e_side, c_side};
        const double leaving[2] = {c.n_s[b] + c.flux[b], c.n_s[b]};
        for (int t = 0; t < 2; ++t) {
            const int64_t x = flank[t];
            const double s_out = leaving[t];
            if (!(boundary_shares_strand(c.fp[b], c.fn[b], c.fp[x], c.fn[x]) && c.a_g[x] > 0)) continue;
            double width = 0.0;
            if (c.has_own[b] && c.has_own[x]) {
                double f_b, v_b, f_x, v_x;
                if (c.strand_mode(b, f_b, v_b) && c.strand_mode(x, f_x, v_x)) {
                    const double lo_b = std::log(f_b / (1.0 - f_b)), lo_x = std::log(f_x / (1.0 - f_x));
                    const double v_ratio = s_out > 0 ? s_out / (c.n_u[b] * (c.n_u[b] + s_out)) : 0.0;
                    const double d = lo_b - lo_x - std::log((c.n_u[b] + s_out) / c.n_u[b]);
                    width = std::max(0.0, d * d - (v_b + v_x + v_ratio));
                }
            }
            face_map_lambda(c.lam, c.K, c.n_u[b], c.a_g[b], c.a_g[b], c.a_g[x], c.a_g[x], s_out / c.a_g[b],
                            S.a.data());
            F.set(x, b, SPLICE_OUT, nullptr, nullptr, c.n_u[b], s_out, c.a_g[b], c.a_g[x], width);
            F.set(b, x, TRANSPORT, S.a.data(), nullptr, c.n_u[b], s_out, 0.0, 0.0, width);
        }
    }
}

// lanes.gdna_lane: the lane's faces and every full node's own level
void gdna_lane(const Chain& c, const RowsOut& own, const FacesOut& F, double rho_ref, LaneOut& L, Scratch& S) {
    const double* u = c.lam;
    for (int64_t x = 0; x < c.n; ++x) {
        const bool empty = !(c.n_u[x] > 0.0) || !(c.a_g[x] > 0.0);
        if (empty || c.intergenic(x)) continue;
        bool gene_edge = false;
        if (c.is_bnd[x]) {
            const int64_t lo = c.left[x], hi = c.right[x];
            gene_edge = (lo >= 0 && c.intergenic(lo)) || (hi >= 0 && c.intergenic(hi));
        }
        if (gene_edge) {
            poisson_level(u, c.K, c.n_u[x], c.a_g[x], rho_ref, S.a.data());
            L.own_level.set(x, S.a.data());
        } else if (const double* o = own.get(x); o != nullptr && ptp(o, c.K) > EPS) {
            level_of_profile(o, c.lam, u, c.K, c.n_u[x], c.a_g[x], rho_ref, S, S.c.data());
            L.own_level.set(x, S.c.data());
        }
    }
    for (int64_t i = 0; i < c.n; ++i) {
        for (int side = 0; side < 2; ++side) {
            const int64_t s = side == 0 ? c.left[i] : c.right[i];
            L.face[i * 2 + side] = s >= 0 && !c.intergenic(i) && !c.intergenic(s) &&
                                   F.kind[static_cast<size_t>(i) * 2 + side] == NONE;
        }
    }
}

// lanes.rna_lanes: one strand's lane — its faces from the flag bits, its sources, its flux levels
void rna_lane(const Chain& c, const RowsOut& own, int strand, double rho_ref, LaneOut& L, Scratch& S) {
    const bool* free = strand == 0 ? c.fp : c.fn;
    const bool* exon_s = strand == 0 ? c.exon_pos : c.exon_neg;
    const int col = strand;
    const int col_read = read_column(col, c.has_strand, c.kappa);
    const double kappa_read = c.has_strand ? std::max(c.kappa, 1.0 - c.kappa) : 0.5;
    const int all_bits = ALL_BITS[strand], term_bits = TERM_BITS[strand];
    auto intron_s = [&](int64_t r) { return !c.is_bnd[r] && free[r] && !exon_s[r]; };
    for (int64_t i = 0; i < c.n; ++i) {
        for (int side = 0; side < 2; ++side) {
            const int64_t nbr = side == 0 ? c.left[i] : c.right[i];
            bool face = false, two = false;
            if (nbr >= 0 && !c.intergenic(i) && free[nbr] && free[i]) {
                const int64_t b = c.is_bnd[nbr] ? nbr : i, reg = c.is_bnd[nbr] ? i : nbr;
                const int f = c.flags[b];
                const bool crossing = (f & all_bits) == 0;
                const bool into_own_intron = !crossing && (f & term_bits) == 0 && intron_s(reg);
                face = crossing || into_own_intron;
                two = (crossing && intron_s(reg)) || into_own_intron;
            }
            L.face[i * 2 + side] = face;
            L.two_sided[i * 2 + side] = two;
        }
    }
    if (!(rho_ref > 0.0)) return;
    for (int64_t x = 0; x < c.n; ++x) {
        if (!free[x]) continue;
        const bool empty = !(c.n_u[x] > 0.0) || !(c.a_r[x] > 0.0);
        const bool single = !(c.fp[x] && c.fn[x]);
        int n_parts = 0;
        double* acc = S.d.data();  // the intersection of the parts, folded as they come
        if (const double* o = own.get(x); !empty && single && o != nullptr && ptp(o, c.K) > EPS) {
            rna_level_of_profile(o, c.lam, c.lam, c.K, c.n_u[x], c.a_r[x], rho_ref, S, acc);
            n_parts = 1;
        }
        double c_sum = 0.0, a_sum = 0.0;
        if (c.is_exon[x]) {
            for (const int64_t b : {c.left[x], c.right[x]}) {
                if (b < 0 || !c.is_bnd[b]) continue;
                const int hi = c.left[x] == b ? 1 : 0;
                const double c_j = (hi ? c.sj_hi : c.sj_lo)[2 * b + col];
                const double r_j = (hi ? c.route_hi : c.route_lo)[2 * b + col];
                if (!(c_j > 0.0 && r_j > 0.0)) continue;
                const double v = hop_price(c_j, c_j / r_j, c.cnt[2 * x + col_read], kappa_read * c.a_r[x]);
                double* fl = S.e.data();
                flux_level(c.lam, c.K, c_j, r_j, rho_ref, v, S, fl);
                if (n_parts == 0) std::copy(fl, fl + c.K, acc); else fold_min(fl, c.K, acc);
                ++n_parts;
                const int side = b < x ? 0 : 1;
                std::copy(fl, fl + c.K, L.flux_rows + (static_cast<size_t>(x) * 2 + side) * c.K);
                L.flux_mask[x * 2 + side] = true;
                c_sum += c_j;
                a_sum += c_j / r_j;
            }
        }
        if (n_parts) {
            norm_inplace(acc, c.K);
            L.own_level.set(x, acc);
            if (empty) {
                L.witness[2 * x] = c_sum;
                L.witness[2 * x + 1] = a_sum;
                L.witness_mask[x] = true;
            }
        }
    }
}

// ---- the call --------------------------------------------------------------------------------------------

LaneOut lane_view(nb::tuple t, int K) {
    LaneOut L;
    L.face = nb::cast<BoolMat>(t[0]).data();
    L.two_sided = nb::cast<BoolMat>(t[1]).data();
    L.own_level = RowsOut{nb::cast<Mat>(t[2]).data(), nb::cast<BoolVec>(t[3]).data(), K};
    L.flux_rows = nb::cast<Cube>(t[4]).data();
    L.flux_mask = nb::cast<BoolMat>(t[5]).data();
    L.witness = nb::cast<Mat>(t[6]).data();
    L.witness_mask = nb::cast<BoolVec>(t[7]).data();
    return L;
}

int transfer_prepare(Vec lam, BoolVec is_boundary, BoolVec is_exon, BoolVec free_pos, BoolVec free_neg,
                     BoolVec exon_pos, BoolVec exon_neg, IdxVec left, IdxVec right, FlagVec flags, Vec n_u,
                     Vec n_s, Vec a_g, Vec a_r, Mat cnt, Vec belief, BoolVec has_own_composition, Vec flux,
                     Mat route_rate_lo, Mat route_rate_hi, Mat sj_count_lo, Mat sj_count_hi, Mat src,
                     bool has_strand, double kappa, double od_g, double od_r, double rho_gdna, double rho_rna,
                     Mat own, BoolVec own_mask, nb::tuple faces, nb::object gdna, nb::tuple pos, nb::tuple neg) {
    Chain c;
    c.n = static_cast<int>(n_u.shape(0));
    c.K = static_cast<int>(lam.shape(0));
    if (static_cast<int>(src.shape(1)) != c.K || static_cast<int>(own.shape(1)) != c.K)
        throw std::invalid_argument("transfer_prepare: the rows and the grid disagree");
    c.lam = lam.data();
    c.is_bnd = is_boundary.data(); c.is_exon = is_exon.data();
    c.fp = free_pos.data(); c.fn = free_neg.data();
    c.exon_pos = exon_pos.data(); c.exon_neg = exon_neg.data();
    c.has_own = has_own_composition.data();
    c.left = left.data(); c.right = right.data();
    c.flags = flags.data();
    c.n_u = n_u.data(); c.n_s = n_s.data(); c.a_g = a_g.data(); c.a_r = a_r.data();
    c.belief = belief.data(); c.cnt = cnt.data(); c.src = src.data();
    c.route_lo = route_rate_lo.data(); c.route_hi = route_rate_hi.data();
    c.sj_lo = sj_count_lo.data(); c.sj_hi = sj_count_hi.data();
    c.flux = flux.data();
    c.has_strand = has_strand; c.kappa = kappa; c.od_g = od_g; c.od_r = od_r;
    c.is_intron.resize(c.n); c.is_intergenic.resize(c.n);
    for (int i = 0; i < c.n; ++i) {
        const bool region = !c.is_bnd[i] && !c.is_exon[i];
        c.is_intron[i] = region && (c.fp[i] || c.fn[i]);
        c.is_intergenic[i] = region && !c.fp[i] && !c.fn[i];
    }
    Scratch S(c.K);
    RowsOut claims_out{own.data(), own_mask.data(), c.K};
    Mat f_rows = nb::cast<Mat>(faces[9]);
    FacesOut F{&c, nb::cast<KindMat>(faces[0]).data(), nb::cast<IdxMat>(faces[1]).data(),
               nb::cast<IdxMat>(faces[2]).data(), nb::cast<Mat>(faces[3]).data(), nb::cast<Mat>(faces[4]).data(),
               nb::cast<Mat>(faces[5]).data(), nb::cast<Mat>(faces[6]).data(), nb::cast<Mat>(faces[7]).data(),
               nb::cast<Mat>(faces[8]).data(), f_rows.data(), static_cast<int>(f_rows.shape(0))};
    claims(c, claims_out, S);
    splice_faces(c, F, S);
    edge_level(c, claims_out, F, S);
    terminus_rules(c, F, S);
    alternative_splice_site(c, F, S);
    if (!gdna.is_none()) {
        LaneOut G = lane_view(nb::cast<nb::tuple>(gdna), c.K);
        gdna_lane(c, claims_out, F, rho_gdna, G, S);
    }
    LaneOut P = lane_view(pos, c.K), N = lane_view(neg, c.K);
    rna_lane(c, claims_out, 0, rho_rna, P, S);
    rna_lane(c, claims_out, 1, rho_rna, N, S);
    return F.n_rows;
}

}  // namespace

NB_MODULE(_prepare_impl, m) {
    m.doc() = "The composition-transfer policy's builders for one block, on the policy's tables.";
    m.def("transfer_prepare", &transfer_prepare, nb::arg("lam"), nb::arg("is_boundary"), nb::arg("is_exon"),
          nb::arg("free_pos"), nb::arg("free_neg"), nb::arg("exon_pos"), nb::arg("exon_neg"), nb::arg("left"),
          nb::arg("right"), nb::arg("flags"), nb::arg("n_u"), nb::arg("n_s"), nb::arg("a_g"), nb::arg("a_r"),
          nb::arg("cnt"), nb::arg("belief"), nb::arg("has_own_composition"), nb::arg("flux"),
          nb::arg("route_rate_lo"), nb::arg("route_rate_hi"), nb::arg("sj_count_lo"), nb::arg("sj_count_hi"),
          nb::arg("src"), nb::arg("has_strand"), nb::arg("kappa"), nb::arg("od_g"), nb::arg("od_r"),
          nb::arg("rho_gdna"), nb::arg("rho_rna"), nb::arg("own"), nb::arg("own_mask"), nb::arg("faces"),
          nb::arg("gdna").none(), nb::arg("pos"), nb::arg("neg"),
          "Build one block's own claims, face rules and level lanes into the tables given: `own` / "
          "`own_mask` (the claims), `faces` (kind, row, row2, n_u, n_s, a_b, a_x, width, var, rows), and per "
          "lane (face, two_sided, own rows, own mask, flux rows, flux mask, witness rows, witness mask) — the "
          "gDNA lane's `None` when the library has no gDNA coordinate. Returns the number of face rows written.");
}
