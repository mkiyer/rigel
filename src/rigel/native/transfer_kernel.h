// transfer_kernel.h — the COMPOSITION TRANSFER policy's three kernels, on plain views of their tables:
//
//   prepare_block   the BUILDERS for one block (phase 0): every node's own claim, the recipient's rule per directed
//                   face and the level lanes, written into the tables the caller provides;
//   pass_block      ONE DIRECTIONAL PASS (phase 1): in chain order the recipient receives what its neighbour sends —
//                   the face's composition rule and each lane's level — writing a received table in place;
//   solve_block     THE SOLVE (phase 2, the policy's half): the two received tables into ψ's two channels — the fused
//                   λ rows and the cube delivery at the AMBIG nodes.
//
// Every table is a VIEW of raw pointers with its sizes, so the same kernels run on the block pipeline's arena
// (`solve_kernel.cpp`, the production path: one native call per sweep, a pool of threads over the blocks) and on the
// numpy arrays the gates allocate through the bindings (the same file). The row constructors are `transfer_rows.h`,
// shared with ψ's kernel. This is the ONE implementation of the policy's arithmetic — there is no Python kernel
// beside it — and its gates are the transfer gates (`tests/calibration/test_transfer_faces.py`,
// `test_transfer_policy.py`, `test_transfer_rna_lanes.py`: the tables, single hops of the pass and the delivered
// channels against independent recomputes) and `profiling/sweep_replay.py replay` on a captured sweep.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "transfer_rows.h"

namespace transfer_kernel {

using namespace transfer_rows;

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

//: the splice-out row's marginal over log rho is taken on equal-probability nodes of the standard normal —
//: quadrature resolution, like ψ's tilt nodes, not a model constant: Φ⁻¹((t + ½)/9), t = 0..8
//: (`scipy.stats.norm.ppf`, every digit numpy prints — a hand-typed digit off at 1e-8 moved f_g by 1e-7)
constexpr int N_MARGINAL_NODES = 9;
constexpr double MARGINAL_NODES[N_MARGINAL_NODES] = {
    -1.5932188180230507, -0.9674215661017012, -0.5894557978497783, -0.28221614706250814, 0.0,
    0.28221614706250825, 0.5894557978497783, 0.9674215661017012, 1.59321881802305};

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
// read_column: the genome-strand column strand `col`'s RNA reads on — its own when the library reads sense
// (kappa >= 1/2, or no fitted strand model), the other under an antisense protocol
inline int read_column(int col, bool has_strand, double kappa) {
    return (!has_strand || kappa >= 0.5) ? col : 1 - col;
}

// ---- the chain as the builders read it ------------------------------------------------------------

struct Chain {
    int n, K;
    const double* lam;
    const bool *is_bnd, *is_exon, *fp, *fn, *exon_pos, *exon_neg, *has_own;
    const int64_t *left, *right;  // block-local, -1 off the block
    const uint16_t* flags;
    const double *n_u, *n_s, *a_g, *a_r, *belief, *cnt;
    const double *route_lo, *route_hi, *sj_lo, *sj_hi;  // (n, 2) by transcript strand
    const double* flux;                                 // the sj count over both strands
    //: the factory's row per slot, or nullptr where the slot has none (an intron with a row is the only reader)
    const double* const* src;
    const char *is_intron, *is_intergenic;
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
    // the two structural classes of a region, from its bits: an intron admits a strand, intergenic none
    static void classify(int n, const bool* is_bnd, const bool* is_exon, const bool* fp, const bool* fn, char* intron,
                         char* intergenic) {
        for (int i = 0; i < n; ++i) {
            const bool region = !is_bnd[i] && !is_exon[i];
            intron[i] = region && (fp[i] || fn[i]);
            intergenic[i] = region && !fp[i] && !fn[i];
        }
    }
};

// ---- the tables the builders write ----------------------------------------------------------------

struct RowsOut {  // OPTIONAL ROWS over the nodes: a matrix (n, K), unfilled, and a presence mask (n,)
    double* rows; bool* mask; int K;
    void set(int64_t i, const double* row) { std::copy(row, row + K, rows + i * K); mask[i] = true; }
    void set_zero(int64_t i) { std::fill(rows + i * K, rows + (i + 1) * K, 0.0); mask[i] = true; }
    const double* get(int64_t i) const { return mask[i] ? rows + i * K : nullptr; }
};

struct RowStore {  // a store of (K,) rows written in order — the face maps, the flux levels — read by index
    std::vector<double>* buf; int K; int n_rows = 0;
    int keep(const double* r) {
        if (r == nullptr) return -1;
        const size_t need = static_cast<size_t>(n_rows + 1) * K;
        if (buf->size() < need) buf->resize(std::max(need, buf->size() * 2));
        std::copy(r, r + K, buf->data() + static_cast<size_t>(n_rows) * K);
        return n_rows++;
    }
    const double* row(int r) const { return buf->data() + static_cast<size_t>(r) * K; }
};

struct FacesOut {  // the composition rules as typed tables over (destination, side)
    const Chain* c;
    int8_t* kind; int32_t* row; int32_t* row2;
    double *n_u, *n_s, *a_b, *a_x, *width, *var;
    RowStore store;

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
        row[at] = store.keep(r);
        row2[at] = store.keep(r2);
        n_u[at] = nu; n_s[at] = ns; a_b[at] = ab; a_x[at] = ax; width[at] = w; var[at] = v;
    }
};

struct LaneOut {  // a level lane's tables: faces, two-sided faces, own levels, junction flux levels, witnesses
    bool* face; bool* two_sided; RowsOut own_level;
    RowStore flux; int32_t* flux_index;      // the flux level per (exon, side): a row of the store, or -1
    double* witness; bool* witness_mask;     // (n, 2), (n,)
};

// ---- the builders — one per shipped message ---------------------------------------------------------

inline void claims(const Chain& c, RowsOut& own, Scratch& S) {
    for (int i = 0; i < c.n; ++i) {
        if (!c.intron(i)) continue;
        const double* r = c.src[i];
        if (r != nullptr && ptp(r, c.K) > EPS) {
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

inline void splice_faces(const Chain& c, FacesOut& F, Scratch& S) {
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

inline void edge_level(const Chain& c, RowsOut& own, FacesOut& F, Scratch& S) {
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

inline void terminus_rules(const Chain& c, FacesOut& F, Scratch& S) {
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

inline void alternative_splice_site(const Chain& c, FacesOut& F, Scratch& S) {
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

// the gDNA lane: the lane's faces and every full node's own level
inline void gdna_lane(const Chain& c, const RowsOut& own, const FacesOut& F, double rho_ref, LaneOut& L, Scratch& S) {
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

// one strand's RNA lane — its faces from the flag bits, its sources, its flux levels
inline void rna_lane(const Chain& c, const RowsOut& own, int strand, double rho_ref, LaneOut& L, Scratch& S) {
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
                L.flux_index[x * 2 + side] = L.flux.keep(fl);
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

// THE BUILDERS IN ONE CALL: every node's own claim, the face rules, the lanes — the gDNA lane only where the library
// has a gDNA coordinate (`gdna` nullptr otherwise; a gDNA-free library still has its RNA lanes)
inline void prepare_block(const Chain& c, RowsOut& own, FacesOut& F, LaneOut* gdna, LaneOut& pos, LaneOut& neg,
                          double rho_gdna, double rho_rna, Scratch& S) {
    claims(c, own, S);
    splice_faces(c, F, S);
    edge_level(c, own, F, S);
    terminus_rules(c, F, S);
    alternative_splice_site(c, F, S);
    if (gdna != nullptr) gdna_lane(c, own, F, rho_gdna, *gdna, S);
    rna_lane(c, own, 0, rho_rna, pos, S);
    rna_lane(c, own, 1, rho_rna, neg, S);
}


// ═══ THE PASS (phase 1) ═══════════════════════════════════════════════════════════════════════════════

// ---- the tables as the pass reads them ------------------------------------------------------------

struct FacesView {
    const int8_t* kind; const int32_t* row; const int32_t* row2;
    const double *n_u, *n_s, *a_b, *a_x, *width, *var;
    const double* rows; int K;
    const double* map(int r) const { return rows + static_cast<size_t>(r) * K; }
};

struct LevelsView {  // one lane's received table, written in place
    bool* present; double* profile; double* count; double* opportunity; bool* has_witness; double* rna_count;
    double* rna_count_var; int K;
    double* row(int i) { return profile + static_cast<size_t>(i) * K; }
    const double* row(int i) const { return profile + static_cast<size_t>(i) * K; }
};

struct ReceivedView {  // what every node holds from one side: the backbone's neighbour bit, the composition, three lanes
    bool* has_neighbour; bool* has_composition; double* composition; LevelsView lanes[3]; int K;
    double* comp(int i) { return composition + static_cast<size_t>(i) * K; }
    const double* comp(int i) const { return composition + static_cast<size_t>(i) * K; }
};

struct LaneView {  // a level lane as the pass reads it
    int field;  // 0 gdna, 1 rna_pos, 2 rna_neg — which received lane it writes
    const bool* face; const bool* two_sided; const bool* empty;
    const double* own_level; const bool* own_mask;
    const double* count; const double* a;
    const double* other;                   // the other column, or nullptr
    const double* witness; const bool* witness_mask;  // (n, 2): count, opportunity
};

inline bool at2(const bool* m, int i, int side) { return m[static_cast<size_t>(i) * 2 + side]; }

inline void levels_write(LevelsView& L, int i, const double* profile, double count, double opportunity, bool witness,
                         double rna_count, double rna_var) {
    L.present[i] = true;
    double* dest = L.row(i);
    if (dest != profile) std::copy(profile, profile + L.K, dest);
    L.count[i] = count;
    L.opportunity[i] = opportunity;
    L.has_witness[i] = witness;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    L.rna_count[i] = witness ? rna_count : nan;
    L.rna_count_var[i] = witness ? rna_var : nan;
}

inline void levels_forward(LevelsView& L, int s, int i) {
    L.present[i] = L.present[s];
    std::copy(L.row(s), L.row(s) + L.K, L.row(i));
    L.count[i] = L.count[s];
    L.opportunity[i] = L.opportunity[s];
    L.has_witness[i] = L.has_witness[s];
    L.rna_count[i] = L.rna_count[s];
    L.rna_count_var[i] = L.rna_count_var[s];
}

// intersect the present parts (pointwise minimum), max-normalised, into out; whether any part
inline bool intersect2(const double* p, const double* q, int K, double* out) {
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

// what s sends toward x on a lane, written into row x; whether anything was sent
inline bool lane_emit(const LaneView& ln, LevelsView& L, int K, int s, int x, int side, Scratch& S) {
    const double* own = ln.own_mask[s] ? ln.own_level + static_cast<size_t>(s) * K : nullptr;
    const bool far = L.present[s];
    if (ln.empty[s] && own == nullptr) {
        if (far) levels_forward(L, s, x);
        return far;
    }
    const double* own_used = own;
    if (own != nullptr && !at2(ln.two_sided, x, side)) {
        lower_side(own, K, S.a.data());
        own_used = S.a.data();
    }
    const double* held = far ? L.row(s) : nullptr;
    if (!intersect2(own_used, held, K, S.b.data())) return false;
    if (ln.empty[s]) {
        if (!ln.witness_mask[s]) throw std::runtime_error("an empty node with an own level has no flux witness");
        const double* w = ln.witness + static_cast<size_t>(s) * 2;
        levels_write(L, x, S.b.data(), w[0], w[1], false, 0.0, 0.0);
        return true;
    }
    if (ln.other != nullptr) {
        const double c_r = ln.count[s], c_o = ln.other[s];
        levels_write(L, x, S.b.data(), c_r, ln.a[s], true, c_r - c_o, c_r + c_o);
    } else {
        levels_write(L, x, S.b.data(), ln.count[s], ln.a[s], false, 0.0, 0.0);
    }
    return true;
}

// a FULL recipient re-prices row x in place
inline void lane_receive(const LaneView& ln, LevelsView& L, int K, const double* lam, int x, int side, Scratch& S) {
    const double n_sent = L.count[x], a_sent = L.opportunity[x];
    const double cnt = ln.count[x], a_x = ln.a[x];
    double v;
    if (ln.other == nullptr || !L.has_witness[x]) {
        v = hop_price(n_sent, a_sent, cnt, a_x);
    } else {
        v = count_logvar(n_sent) + count_logvar(cnt);
        const double n_s = L.rna_count[x], v_s = L.rna_count_var[x];
        const double n_x = cnt - ln.other[x], v_x = cnt + ln.other[x];
        if (n_s > 0.0 && n_x > 0.0) {
            const double l = std::log((n_x / a_x) / (n_s / a_sent));
            v += std::max(0.0, l * l - (v_s / (n_s * n_s) + v_x / (n_x * n_x)));
        }
    }
    const double* held = L.row(x);
    const double* p = held;
    if (!at2(ln.two_sided, x, side)) {
        lower_side(held, K, S.a.data());
        p = S.a.data();
    }
    if (v > 0.0) {
        blur_row(p, K, lam, v, S.b.data(), S.blur);
        p = S.b.data();
    }
    if (ln.other != nullptr) {
        const double c_o = ln.other[x];
        levels_write(L, x, p, cnt, a_x, true, cnt - c_o, cnt + c_o);
    } else {
        levels_write(L, x, p, cnt, a_x, false, 0.0, 0.0);
    }
}

// the composition rule at the face into i from side, on what the sender sends (its own claim, what it holds):
// the row for i in out, or nullptr for no claim
inline const double* faces_apply(const FacesView& F, int i, int side, const double* own, const double* held,
                                 const double* lam, int K, Scratch& S, double* out) {
    const size_t at = static_cast<size_t>(i) * 2 + side;
    const int k = F.kind[at];
    if (k == FORWARD) {
        if (own == nullptr && held == nullptr) return nullptr;
        for (int j = 0; j < K; ++j) out[j] = (own ? own[j] : 0.0) + (held ? held[j] : 0.0);
        norm_inplace(out, K);
        return out;
    }
    if (k == EDGE) {
        const int r = F.row[at];
        if (r < 0) return nullptr;
        const double* row = F.map(r);
        if (ptp(row, K) <= EPS) return nullptr;
        std::copy(row, row + K, out);
        return out;
    }
    if (k == LEVEL) {
        if (own == nullptr) {
            const int r2 = F.row2[at];
            if (r2 < 0) return nullptr;
            std::copy(F.map(r2), F.map(r2) + K, out);
            return out;
        }
        level_row(own, lam, K, F.map(F.row[at]), F.var[at], S, out);
        return out;
    }
    if (k != TRANSPORT && k != SPLICE_OUT) return nullptr;
    if (own == nullptr && held == nullptr) return nullptr;
    double* sending = S.e.data();
    for (int j = 0; j < K; ++j) sending[j] = (own ? own[j] : 0.0) + (held ? held[j] : 0.0);
    norm_inplace(sending, K);
    const double n_u = F.n_u[at], n_s = F.n_s[at];
    double* tmp = S.f.data();
    if (k == TRANSPORT) {
        transport_row(sending, lam, K, F.map(F.row[at]), n_u, n_s, S, tmp);
    } else {
        splice_out_row(sending, lam, K, n_u, n_s, F.a_b[at], F.a_x[at], MARGINAL_NODES, N_MARGINAL_NODES, S, tmp);
    }
    const double w = F.width[at];
    if (w > 0.0) {
        blur_row(tmp, K, lam, w, out, S.blur);
    } else {
        std::copy(tmp, tmp + K, out);
    }
    return out;
}

// ONE DIRECTIONAL PASS: for every destination in `seq` (chain order) with a neighbour `nbr[i] >= 0` that is not a
// terminal, the face's composition rule and each lane's level, written into the received table in place. The
// backbone's `has_neighbour` bit is the caller's to set.
inline void pass_block(const double* lam, int K, const int64_t* seq, int n, const int64_t* nbr, const bool* term,
                       const double* own, const bool* own_mask, const FacesView& F, const std::vector<LaneView>& lanes,
                       ReceivedView& R, Scratch& S, std::vector<double>& out) {
    out.resize(K);
    for (int idx = 0; idx < n; ++idx) {
        const int i = static_cast<int>(seq[idx]);
        const int s = static_cast<int>(nbr[i]);
        if (s < 0 || term[i]) continue;
        const int side = s < i ? 0 : 1;
        if (F.kind[static_cast<size_t>(i) * 2 + side] != NONE) {
            const double* own_s = own_mask[s] ? own + static_cast<size_t>(s) * K : nullptr;
            const double* held = R.has_composition[s] ? R.comp(s) : nullptr;
            const double* r = faces_apply(F, i, side, own_s, held, lam, K, S, out.data());
            if (r != nullptr && ptp(r, K) > EPS) {
                double* dest = R.comp(i);
                std::copy(r, r + K, dest);
                norm_inplace(dest, K);
                R.has_composition[i] = true;
            }
        }
        for (const LaneView& ln : lanes) {
            if (!at2(ln.face, i, side)) continue;
            LevelsView& T = R.lanes[ln.field];
            if (lane_emit(ln, T, K, s, i, side, S) && !ln.empty[i]) lane_receive(ln, T, K, lam, i, side, S);
        }
    }
}


// ═══ THE SOLVE (phase 2) ══════════════════════════════════════════════════════════════════════════════
//
// At every node the two held compositions add (independent witnesses about one slot); a held gDNA level on either
// side is read as the node's composition row through its own total (two bounds on one density intersect) and joins
// them; the fused row is re-normalised. THE CEILINGS at single-strand nodes: an RNA level of the node's live strand,
// read ONLY from a face that sent no composition (a licensed face's map already carries the flux and a composition
// already carries its sender's witnesses), intersected with the node's own junction flux at that face, says "at most
// this much gDNA" and joins the row. THE CUBE at AMBIG nodes: per strand the held levels intersected with the node's
// own level's lower side, delivered with the node's total, RNA opportunity and the lanes' reference density — the
// ingredients ψ evaluates at its own θ nodes (`psi_kernel.h`).

struct SolveLane {  // a level lane as the solve reads it
    bool built = false; int field = 0;
    const bool* empty = nullptr; const double* total = nullptr; const double* a = nullptr; double rho_ref = 0.0;
    const double* own_rows = nullptr; const bool* own_mask = nullptr;   // the node's own level (an RNA lane's sources)
    const double* flux_store = nullptr; const int32_t* flux_index = nullptr;  // the junction flux levels per (node, side)
};

struct CubeOut {  // the cube delivery, one row per delivered AMBIG node, written in order
    int64_t* slot; double* pos; bool* has_pos; double* neg; bool* has_neg; double* total; double* opportunity;
    double* rho; int capacity;
};

// the pointwise-minimum accumulator of bounds on one density: the first part is copied, the rest folded
struct Intersection {
    double* acc; int K; int n = 0;
    void add(const double* p) {
        if (n++ == 0) std::copy(p, p + K, acc); else fold_min(p, K, acc);
    }
    bool any() const { return n > 0; }
    void finish() { norm_inplace(acc, K); }
};

// the two received tables into ψ's channels: the fused λ rows (written in place, zero where nothing fused; `written`
// marks the rows that were) and the cube delivery. Returns (live, the number of cube rows delivered).
inline std::pair<bool, int> solve_block(const double* lam, int K, int n, const ReceivedView sides[2], const SolveLane& G,
                                        const SolveLane rna[2], const bool* free_pos, const bool* free_neg,
                                        double* out_rows, bool* written, CubeOut& cube, Scratch& S,
                                        std::vector<double>& bound, std::vector<double>& row) {
    const double* u = lam;  // a level's coordinate is the solve grid
    const bool* free_of[2] = {free_pos, free_neg};
    bound.resize(K); row.resize(K);
    bool live = false;
    int d = 0;
    for (int i = 0; i < n; ++i) {
        double* r = out_rows + static_cast<size_t>(i) * K;
        const bool ambig = free_pos[i] && free_neg[i];
        // the two held compositions add
        bool fused = false;
        for (int side = 0; side < 2; ++side) {
            const ReceivedView& h = sides[side];
            if (!h.has_composition[i]) continue;
            const double* c = h.comp(i);
            if (!fused) std::copy(c, c + K, r); else for (int j = 0; j < K; ++j) r[j] += c[j];
            fused = true;
        }
        // a held gDNA level on either side, read through the node's own total, joins as one more witness
        if (G.built && !G.empty[i]) {
            Intersection inter{bound.data(), K};
            for (int side = 0; side < 2; ++side) {
                const LevelsView& L = sides[side].lanes[G.field];
                if (!L.present[i]) continue;
                profile_of_level(L.row(i), u, u, K, G.total[i], G.a[i], G.rho_ref, S, row.data());
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
        if (!ambig) {
            for (int s = 0; s < 2; ++s) {
                const SolveLane& L = rna[s];
                if (!L.built || !free_of[s][i] || L.empty[i]) continue;
                Intersection inter{bound.data(), K};
                for (int side = 0; side < 2; ++side) {
                    const ReceivedView& h = sides[side];
                    if (!h.has_neighbour[i] || h.has_composition[i]) continue;
                    const LevelsView& lv = h.lanes[L.field];
                    if (lv.present[i]) inter.add(lv.row(i));
                    const int32_t fr = L.flux_index[static_cast<size_t>(i) * 2 + side];
                    if (fr >= 0) inter.add(L.flux_store + static_cast<size_t>(fr) * K);
                }
                if (!inter.any()) continue;
                inter.finish();
                rna_row_of_level(bound.data(), u, u, K, L.total[i], L.a[i], L.rho_ref, S, row.data());
                if (ptp(row.data(), K) <= EPS) continue;
                if (ptp(r, K) > EPS) { for (int j = 0; j < K; ++j) r[j] += row[j]; norm_inplace(r, K); }
                else std::copy(row.begin(), row.end(), r);
                live = true;
            }
        }
        written[i] = fused || ptp(r, K) > EPS;
        // THE CUBE at an AMBIG node: per strand, the held levels and the own level's lower side intersected
        if (ambig && rna[0].built && rna[1].built && !rna[0].empty[i]) {
            if (d >= cube.capacity) throw std::runtime_error("transfer_solve: the cube table is full");
            bool any_profile = false;
            for (int s = 0; s < 2; ++s) {
                const SolveLane& L = rna[s];
                double* dest = (s == 0 ? cube.pos : cube.neg) + static_cast<size_t>(d) * K;
                Intersection inter{dest, K};
                for (int side = 0; side < 2; ++side) {
                    const LevelsView& lv = sides[side].lanes[L.field];
                    if (lv.present[i]) inter.add(lv.row(i));
                }
                if (L.own_mask[i]) {
                    lower_side(L.own_rows + static_cast<size_t>(i) * K, K, row.data());
                    inter.add(row.data());
                }
                if (inter.any()) inter.finish();
                (s == 0 ? cube.has_pos : cube.has_neg)[d] = inter.any();
                any_profile |= inter.any();
            }
            if (any_profile) {
                cube.slot[d] = i;
                cube.total[d] = rna[0].total[i];
                cube.opportunity[d] = rna[0].a[i];
                cube.rho[d] = rna[0].rho_ref;
                ++d;
            } else {
                cube.has_pos[d] = false; cube.has_neg[d] = false;
            }
        }
    }
    return {live, d};
}

}  // namespace transfer_kernel
