// transfer_kernel.cpp — the COMPOSITION TRANSFER policy's three kernels, one native module (`_transfer_impl`):
//
//   transfer_prepare  the BUILDERS for one block (phase 0): every node's own claim, the recipient's rule per
//                     directed face and the level lanes, written into the tables `TransferPolicy.prepare`
//                     allocates (`calibration/messages/transfer.py`, `faces.py`, `lanes.py`);
//   transfer_pass     ONE DIRECTIONAL PASS (phase 1): in chain order the recipient receives what its neighbour
//                     sends — the face's composition rule and each lane's level — writing the `Received`
//                     tables in place (`calibration/sweep.py::_pass`);
//   transfer_solve    THE SOLVE (phase 2, the policy's half): the two held tables into ψ's two channels — the
//                     fused λ rows and the cube delivery at the AMBIG nodes (`_PreparedTransfer.solve`).
//
// The row constructors are `transfer_rows.h`, shared with ψ's kernel. This file is the ONE implementation
// of the policy's arithmetic — there is no Python kernel beside it — and its gates are the transfer gates
// (`tests/calibration/test_transfer_faces.py`, `test_transfer_policy.py`, `test_transfer_rna_lanes.py`: the
// tables, single hops of the pass and the delivered channels against independent recomputes, the row
// constructors and flag predicates held to their analytic properties through the `rows` submodule bound at
// the end of this file) and `profiling/sweep_replay.py replay` on a captured sweep.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>

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

// ═══ THE BUILDERS (phase 0) ═══════════════════════════════════════════════════════════════════════════

// ---- the boundary flags (calibration/splice_graph.py) ----------------------------------------------

constexpr int TSS_POS = 1 << 0, TSS_NEG = 1 << 1, TES_POS = 1 << 2, TES_NEG = 1 << 3;
constexpr int DON_POS = 1 << 4, DON_NEG = 1 << 5, ACC_POS = 1 << 6, ACC_NEG = 1 << 7;
constexpr int TERMINUS = TSS_POS | TSS_NEG | TES_POS | TES_NEG;
constexpr int SJ_FLAGS = DON_POS | DON_NEG | ACC_POS | ACC_NEG;
//: a transcript body that extends genomic-RIGHT from its terminus leaves the OUTSIDE flank on the left
constexpr int BODY_RIGHT = TSS_POS | TES_NEG, BODY_LEFT = TES_POS | TSS_NEG;
//: a DONOR bit marks the intron's LOW end (the intron lies right), an ACCEPTOR its HIGH end
constexpr int INTRON_RIGHT = DON_POS | DON_NEG, INTRON_LEFT = ACC_POS | ACC_NEG;
//: a strand's four boundary bits and its terminus bits, by strand 0 = +, 1 = −
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

LaneOut lane_tables(nb::tuple t, int K) {
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
        LaneOut G = lane_tables(nb::cast<nb::tuple>(gdna), c.K);
        gdna_lane(c, claims_out, F, rho_gdna, G, S);
    }
    LaneOut P = lane_tables(pos, c.K), N = lane_tables(neg, c.K);
    rna_lane(c, claims_out, 0, rho_rna, P, S);
    rna_lane(c, claims_out, 1, rho_rna, N, S);
    return F.n_rows;
}


// ═══ THE PASS (phase 1) ═══════════════════════════════════════════════════════════════════════════════

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


// ═══ THE SOLVE (phase 2) ══════════════════════════════════════════════════════════════════════════════

// ---- THE SOLVE (phase 2, the policy's half): the two held tables into ψ's channels ----------------------
//
// `_PreparedTransfer.solve` (`calibration/messages/transfer.py`): at every node the two held compositions
// add (independent witnesses about one slot); a held gDNA level on either side is read as the node's
// composition row through its own total (two bounds on one density intersect) and joins them; the fused
// row is re-normalised. THE CEILINGS at single-strand nodes: an RNA level of the node's live strand, read
// ONLY from a face that sent no composition (a licensed face's map already carries the flux and a
// composition already carries its sender's witnesses), intersected with the node's own junction flux at
// that face, says "at most this much gDNA" and joins the row. THE CUBE at AMBIG nodes: per strand the
// held levels intersected with the node's own level's lower side, delivered with the node's total, RNA
// opportunity and the lanes' reference density — the ingredients ψ evaluates at its own θ nodes
// (`psi_kernel.cpp`). Gates: `tests/calibration/test_transfer_policy.py`, `test_transfer_rna_lanes.py`.

struct Held {  // one side's Received table, read only
    BoolVec has_neighbour, has_composition; Mat composition;
    std::array<BoolVec, 3> present; std::array<Mat, 3> profile;  // the three level lanes, in Received.LANES order
};

struct SolveLane {  // a level lane as the solve reads it
    bool built = false; int field = 0;
    BoolVec empty; Vec total; Vec a; double rho_ref = 0.0;
    Mat own_rows; BoolVec own_mask;      // the node's own level (an RNA lane's sources)
    Cube flux_rows; BoolMat flux_mask;   // the junction flux levels per (node, side)
};

Held held_view(nb::tuple t) {
    Held h{nb::cast<BoolVec>(t[0]), nb::cast<BoolVec>(t[1]), nb::cast<Mat>(t[2]), {}, {}};
    for (int l = 0; l < 3; ++l) {
        nb::tuple lv = nb::cast<nb::tuple>(t[3 + l]);
        h.present[l] = nb::cast<BoolVec>(lv[0]);
        h.profile[l] = nb::cast<Mat>(lv[1]);
    }
    return h;
}

SolveLane solve_lane(nb::object o) {
    SolveLane L;
    if (o.is_none()) return L;
    nb::tuple t = nb::cast<nb::tuple>(o);
    L.built = true;
    L.field = nb::cast<int>(t[0]);
    L.empty = nb::cast<BoolVec>(t[1]); L.total = nb::cast<Vec>(t[2]); L.a = nb::cast<Vec>(t[3]);
    L.rho_ref = nb::cast<double>(t[4]);
    L.own_rows = nb::cast<Mat>(t[5]); L.own_mask = nb::cast<BoolVec>(t[6]);
    L.flux_rows = nb::cast<Cube>(t[7]); L.flux_mask = nb::cast<BoolMat>(t[8]);
    return L;
}

// the pointwise-minimum accumulator of bounds on one density: the first part is copied, the rest folded
struct Intersection {
    double* acc; int K; int n = 0;
    void add(const double* p) {
        if (n++ == 0) std::copy(p, p + K, acc); else fold_min(p, K, acc);
    }
    bool any() const { return n > 0; }
    void finish() { norm_inplace(acc, K); }
};

std::pair<bool, int> transfer_solve(Vec lam, nb::tuple from_left, nb::tuple from_right, nb::object gdna,
                                    nb::object pos, nb::object neg, BoolVec ambig, BoolVec free_pos,
                                    BoolVec free_neg, Mat out_rows, IdxVec cube_slot, Mat cube_pos,
                                    BoolVec cube_has_pos, Mat cube_neg, BoolVec cube_has_neg, Vec cube_total,
                                    Vec cube_opportunity, Vec cube_rho) {
    const int K = static_cast<int>(lam.shape(0)), n = static_cast<int>(out_rows.shape(0));
    if (static_cast<int>(out_rows.shape(1)) != K) throw std::invalid_argument("transfer_solve: the rows and the grid disagree");
    const double* u = lam.data();  // a level's coordinate is the solve grid
    Held sides[2] = {held_view(from_left), held_view(from_right)};
    SolveLane G = solve_lane(gdna);
    SolveLane rna[2] = {solve_lane(pos), solve_lane(neg)};
    const bool* free_of[2] = {free_pos.data(), free_neg.data()};
    Scratch S(K);
    std::vector<double> bound(K), row(K);
    const int cube_capacity = static_cast<int>(cube_slot.shape(0));
    bool live = false;
    int d = 0;
    for (int i = 0; i < n; ++i) {
        double* r = out_rows.data() + static_cast<size_t>(i) * K;
        // the two held compositions add
        bool fused = false;
        for (const Held& h : sides) {
            if (!h.has_composition.data()[i]) continue;
            const double* c = h.composition.data() + static_cast<size_t>(i) * K;
            if (!fused) std::copy(c, c + K, r); else for (int j = 0; j < K; ++j) r[j] += c[j];
            fused = true;
        }
        // a held gDNA level on either side, read through the node's own total, joins as one more witness
        if (G.built && !G.empty.data()[i]) {
            Intersection inter{bound.data(), K};
            for (const Held& h : sides) {
                if (!h.present[G.field].data()[i]) continue;
                profile_of_level(h.profile[G.field].data() + static_cast<size_t>(i) * K, u, u, K, G.total.data()[i],
                                 G.a.data()[i], G.rho_ref, S, row.data());
                inter.add(row.data());
            }
            if (inter.any()) {
                inter.finish();
                if (fused) for (int j = 0; j < K; ++j) r[j] += bound[j]; else std::copy(bound.begin(), bound.end(), r);
                fused = true;
            }
        }
        if (fused) { norm_inplace(r, K); live = true; }
        // THE CEILINGS at a single-strand node, from the faces that sent no composition
        if (!ambig.data()[i]) {
            for (int s = 0; s < 2; ++s) {
                const SolveLane& L = rna[s];
                if (!L.built || !free_of[s][i] || L.empty.data()[i]) continue;
                Intersection inter{bound.data(), K};
                for (int side = 0; side < 2; ++side) {
                    const Held& h = sides[side];
                    if (!h.has_neighbour.data()[i] || h.has_composition.data()[i]) continue;
                    if (h.present[L.field].data()[i]) inter.add(h.profile[L.field].data() + static_cast<size_t>(i) * K);
                    if (L.flux_mask.data()[static_cast<size_t>(i) * 2 + side])
                        inter.add(L.flux_rows.data() + (static_cast<size_t>(i) * 2 + side) * K);
                }
                if (!inter.any()) continue;
                inter.finish();
                rna_row_of_level(bound.data(), u, u, K, L.total.data()[i], L.a.data()[i], L.rho_ref, S, row.data());
                if (ptp(row.data(), K) <= EPS) continue;
                if (ptp(r, K) > EPS) { for (int j = 0; j < K; ++j) r[j] += row[j]; norm_inplace(r, K); }
                else std::copy(row.begin(), row.end(), r);
                live = true;
            }
        }
        // THE CUBE at an AMBIG node: per strand, the held levels and the own level's lower side intersected
        if (ambig.data()[i] && rna[0].built && rna[1].built && !rna[0].empty.data()[i]) {
            if (d >= cube_capacity) throw std::runtime_error("transfer_solve: the cube table is full");
            bool any_profile = false;
            for (int s = 0; s < 2; ++s) {
                const SolveLane& L = rna[s];
                double* dest = (s == 0 ? cube_pos.data() : cube_neg.data()) + static_cast<size_t>(d) * K;
                Intersection inter{dest, K};
                for (const Held& h : sides)
                    if (h.present[L.field].data()[i]) inter.add(h.profile[L.field].data() + static_cast<size_t>(i) * K);
                if (L.own_mask.data()[i]) {
                    lower_side(L.own_rows.data() + static_cast<size_t>(i) * K, K, row.data());
                    inter.add(row.data());
                }
                if (inter.any()) inter.finish();
                (s == 0 ? cube_has_pos : cube_has_neg).data()[d] = inter.any();
                any_profile |= inter.any();
            }
            if (any_profile) {
                cube_slot.data()[d] = i;
                cube_total.data()[d] = rna[0].total.data()[i];
                cube_opportunity.data()[d] = rna[0].a.data()[i];
                cube_rho.data()[d] = rna[0].rho_ref;
                ++d;
            } else {
                cube_has_pos.data()[d] = false; cube_has_neg.data()[d] = false;
            }
        }
    }
    return {live, d};
}


// ═══ THE ROW CONSTRUCTORS AND FLAG PREDICATES, BOUND FOR THE GATES ═══════════════════════════════════════
// The one implementation of every row constructor (`transfer_rows.h`) and flag predicate (above) is read by
// the gates through these bindings — as ψ's is through `psi_cube` — each a fresh array from the constructor's
// own arguments, on a scratch of its own. Nothing in `src/` reads them.

using Row = nb::ndarray<nb::numpy, double, nb::ndim<1>>;

Row make_row(std::vector<double>&& v) {
    auto* p = new std::vector<double>(std::move(v));
    nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<double>*>(q); });
    const size_t shape[1] = {p->size()};
    return Row(p->data(), 1, shape, std::move(owner));
}

int K_of(const Vec& v) { return static_cast<int>(v.shape(0)); }

nb::object flank(int64_t x) { return x < 0 ? nb::none() : nb::cast(x); }

void bind_rows(nb::module_& m) {
    nb::module_ r = m.def_submodule(
        "rows", "The row constructors of transfer_rows.h and the builders' flag predicates, bound for the gates.");
    r.attr("EPS") = EPS;
    r.def("trigamma", &trigamma, nb::arg("x"), "zeta(2, x), the counting variance's one home.");
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
}

}  // namespace

NB_MODULE(_transfer_impl, m) {
    m.doc() = "The composition-transfer policy's builders, directional pass and solve, on the policy's tables.";
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
    m.def("transfer_pass", &transfer_pass, nb::arg("lam"), nb::arg("seq"), nb::arg("nbr"),
          nb::arg("terminal"), nb::arg("own"), nb::arg("own_mask"), nb::arg("f_kind"), nb::arg("f_row"),
          nb::arg("f_row2"), nb::arg("f_n_u"), nb::arg("f_n_s"), nb::arg("f_a_b"), nb::arg("f_a_x"),
          nb::arg("f_width"), nb::arg("f_var"), nb::arg("f_rows"), nb::arg("marginal_nodes"),
          nb::arg("lanes"), nb::arg("has_composition"), nb::arg("composition"), nb::arg("levels"),
          "Run one pass in chain order: for every destination in `seq` with a neighbour `nbr[i] >= 0` "
          "that is not a terminal, apply the face's composition rule and carry each lane's level, "
          "writing the Received tables in place.");
    m.def("transfer_solve", &transfer_solve, nb::arg("lam"), nb::arg("from_left"), nb::arg("from_right"),
          nb::arg("gdna").none(), nb::arg("pos").none(), nb::arg("neg").none(), nb::arg("ambig"),
          nb::arg("free_pos"), nb::arg("free_neg"), nb::arg("out_rows"), nb::arg("cube_slot"), nb::arg("cube_pos"),
          nb::arg("cube_has_pos"), nb::arg("cube_neg"), nb::arg("cube_has_neg"), nb::arg("cube_total"),
          nb::arg("cube_opportunity"), nb::arg("cube_rho"),
          "The policy's solve for one block: the two held tables into the fused λ rows (written in place) and the "
          "cube delivery at the AMBIG nodes (the cube arrays, written in place). Returns (live, delivered rows).");
    bind_rows(m);
}
