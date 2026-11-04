// transfer_store.h — the transfer kernels' tables over storage they own.
//
// The kernels (`transfer_kernel.h`) read and write views of raw pointers. These structs own the arrays behind the
// views for one block: the chain's per-node sums, the builders' tables, a received table and the cube delivery. Each
// is sized per block and keeps its capacity, so the block pipeline's arena reuses it across blocks and the gates'
// bindings (`testing/solve.cpp`) build it fresh. Storage the kernels read before writing is cleared here.
#pragma once

#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "transfer_kernel.h"

namespace transfer_kernel {

static_assert(sizeof(bool) == 1, "the bit tables are byte arrays");

inline bool* as_bool(std::vector<uint8_t>& v) { return reinterpret_cast<bool*>(v.data()); }
inline const bool* as_bool(const std::vector<uint8_t>& v) { return reinterpret_cast<const bool*>(v.data()); }

template <class V>
void grow(V& v, size_t m) {
    if (v.size() < m) v.resize(m);
}

// The chain's observation arrays as a call receives them, per node: (n,), or (n, 2) by strand column for the counts and
// the rates.
struct ChainArrays {
    const bool *is_bnd, *is_exon, *fp, *fn, *exon_pos, *exon_neg, *terminal;
    const uint16_t* flags;
    const double *cnt, *spl, *sj, *sj_lo, *sj_hi, *route_lo, *route_hi;  // (n, 2)
    const double *a_g, *a_r;

    // The same arrays from node `start` on.
    ChainArrays from(int64_t start) const {
        return ChainArrays{is_bnd + start,      is_exon + start,      fp + start,           fn + start,
                           exon_pos + start,    exon_neg + start,     terminal + start,     flags + start,
                           cnt + 2 * start,     spl + 2 * start,      sj + 2 * start,       sj_lo + 2 * start,
                           sj_hi + 2 * start,   route_lo + 2 * start, route_hi + 2 * start, a_g + start,
                           a_r + start};
    }
};

// The per-node sums the builders and the solve read: the unspliced and spliced totals and the junction flux over both
// strand columns, and the two unspliced columns apart.
struct ChainSums {
    std::vector<double> n_u, n_s, flux, cnt_col[2];

    void fill(const ChainArrays& c, int n) {
        for (auto* v : {&n_u, &n_s, &flux, &cnt_col[0], &cnt_col[1]}) grow(*v, n);
        for (int i = 0; i < n; ++i) {
            n_u[i] = c.cnt[2 * i] + c.cnt[2 * i + 1];
            n_s[i] = c.spl[2 * i] + c.spl[2 * i + 1];
            flux[i] = c.sj[2 * i] + c.sj[2 * i + 1];
            cnt_col[0][i] = c.cnt[2 * i];
            cnt_col[1][i] = c.cnt[2 * i + 1];
        }
    }
};

// The chain as the builders read it: the observations, their sums, the block-local links (-1 off the block) and each
// node's own-evidence bit.
inline Chain chain_of(const ChainArrays& a, const ChainSums& s, int n, const int64_t* left, const int64_t* right,
                      const bool* has_own, const double* lam, int K, double kappa) {
    Chain c;
    c.n = n; c.K = K; c.lam = lam;
    c.is_bnd = a.is_bnd; c.is_exon = a.is_exon; c.fp = a.fp; c.fn = a.fn; c.exon_pos = a.exon_pos;
    c.exon_neg = a.exon_neg; c.has_own = has_own; c.left = left; c.right = right; c.flags = a.flags;
    c.n_u = s.n_u.data(); c.n_s = s.n_s.data(); c.a_g = a.a_g; c.a_r = a.a_r; c.cnt = a.cnt;
    c.route_lo = a.route_lo; c.route_hi = a.route_hi; c.sj_lo = a.sj_lo; c.sj_hi = a.sj_hi; c.flux = s.flux.data();
    c.terminal = a.terminal;
    c.kappa = kappa;
    return c;
}

// The builders' tables for one block: every node's own claim, the face rules over (destination, side), and the three
// level lanes (0 gDNA, 1 RNA+, 2 RNA−).
struct BuilderTables {
    struct Lane {
        std::vector<uint8_t> face, two_sided, own_mask, witness_mask;   // (n, 2), (n, 2), (n,), (n,)
        std::vector<double> own_rows, witness, flux_store;              // (n, K), (n, 2), the flux level rows
        std::vector<int32_t> flux_index;                                // (n, 2): a row of flux_store, or -1
        int n_flux_rows = 0;
    };
    int K = 0;
    std::vector<double> own_rows, f_par[9], face_store;   // (n, K), (n, 2) each, the face map rows
    std::vector<uint8_t> own_mask;
    std::vector<int8_t> f_kind;
    std::vector<int32_t> f_row, f_row2;
    int n_face_rows = 0;
    Lane lanes[3];

    // The builders on chain c into the tables; the gDNA lane is built only where the library has a gDNA coordinate.
    void build(const Chain& c, double rho_gdna, double rho_rna, Scratch& S) {
        size(c.n, c.K);
        RowsOut own{own_rows.data(), as_bool(own_mask), K};
        FacesOut F{&c, f_kind.data(), f_row.data(), f_row2.data(), f_par[0].data(), f_par[1].data(), f_par[2].data(),
                   f_par[3].data(), f_par[4].data(), f_par[5].data(), f_par[6].data(), f_par[7].data(),
                   f_par[8].data(), RowStore{&face_store, K}};
        LaneOut out[3];
        for (int l = 0; l < 3; ++l) {
            Lane& L = lanes[l];
            out[l] = LaneOut{as_bool(L.face), as_bool(L.two_sided), RowsOut{L.own_rows.data(), as_bool(L.own_mask), K},
                             RowStore{&L.flux_store, K}, L.flux_index.data(), L.witness.data(), as_bool(L.witness_mask)};
        }
        prepare_block(c, own, F, rho_gdna > 0.0 ? &out[0] : nullptr, out[1], out[2], rho_gdna, rho_rna, S);
        n_face_rows = F.store.n_rows;
        for (int l = 0; l < 3; ++l) lanes[l].n_flux_rows = out[l].flux.n_rows;
    }

    FacesView faces() const {
        return FacesView{f_kind.data(), f_row.data(), f_row2.data(), f_par[0].data(), f_par[1].data(),
                         f_par[2].data(), f_par[3].data(), f_par[4].data(), f_par[5].data(), f_par[6].data(),
                         f_par[7].data(), f_par[8].data(), face_store.data(), K};
    }

private:
    // Sized for n nodes, and cleared where the builders read before they write: the presence masks, the rule kinds,
    // the rows a rule names and the flux index.
    void size(int n, int K_) {
        K = K_;
        const size_t nk = static_cast<size_t>(n) * K;
        grow(own_rows, nk);
        grow(own_mask, n);
        for (auto& v : f_par) grow(v, 2 * n);
        grow(f_kind, 2 * n); grow(f_row, 2 * n); grow(f_row2, 2 * n);
        std::fill_n(own_mask.begin(), n, 0);
        std::fill_n(f_kind.begin(), 2 * n, static_cast<int8_t>(NONE));
        std::fill_n(f_row.begin(), 2 * n, -1);
        std::fill_n(f_row2.begin(), 2 * n, -1);
        for (Lane& L : lanes) {
            grow(L.face, 2 * n); grow(L.two_sided, 2 * n); grow(L.own_mask, n); grow(L.witness_mask, n);
            grow(L.own_rows, nk); grow(L.witness, 2 * n); grow(L.flux_index, 2 * n);
            std::fill_n(L.own_mask.begin(), n, 0);
            std::fill_n(L.witness_mask.begin(), n, 0);
            std::fill_n(L.flux_index.begin(), 2 * n, -1);
        }
    }
};

// One side's received table (`ReceivedView`'s arrays): the neighbour and composition bits, the composition, and the
// three level lanes.
struct ReceivedStore {
    struct Lane {
        std::vector<uint8_t> present, has_witness;
        std::vector<double> profile, count, opportunity, rna_count, rna_count_var;
    };
    std::vector<uint8_t> has_neighbour, has_composition;
    std::vector<double> composition;
    Lane lanes[3];

    // The table's view over n nodes, its presence bits cleared and its neighbour bits read from nbr (-1: none).
    ReceivedView view(int n, int K, const int64_t* nbr) {
        const size_t nk = static_cast<size_t>(n) * K;
        grow(has_neighbour, n); grow(has_composition, n); grow(composition, nk);
        std::fill_n(has_composition.begin(), n, 0);
        ReceivedView R;
        R.has_neighbour = as_bool(has_neighbour);
        R.has_composition = as_bool(has_composition);
        R.composition = composition.data();
        R.K = K;
        for (int l = 0; l < 3; ++l) {
            Lane& L = lanes[l];
            grow(L.present, n); grow(L.has_witness, n); grow(L.profile, nk);
            for (auto* v : {&L.count, &L.opportunity, &L.rna_count, &L.rna_count_var}) grow(*v, n);
            std::fill_n(L.present.begin(), n, 0);
            R.lanes[l] = LevelsView{as_bool(L.present), L.profile.data(), L.count.data(), L.opportunity.data(),
                                    as_bool(L.has_witness), L.rna_count.data(), L.rna_count_var.data(), K};
        }
        for (int i = 0; i < n; ++i) R.has_neighbour[i] = nbr[i] >= 0;
        return R;
    }

    // The table's first n nodes, as the diagnostics capture keeps them.
    ReceivedStore prefix(int n, int K) const {
        auto head = [](const auto& v, size_t m) { return std::decay_t<decltype(v)>(v.begin(), v.begin() + m); };
        const size_t nk = static_cast<size_t>(n) * K;
        ReceivedStore out;
        out.has_neighbour = head(has_neighbour, n);
        out.has_composition = head(has_composition, n);
        out.composition = head(composition, nk);
        for (int l = 0; l < 3; ++l) {
            const Lane& L = lanes[l];
            out.lanes[l] = Lane{head(L.present, n), head(L.has_witness, n), head(L.profile, nk), head(L.count, n),
                                head(L.opportunity, n), head(L.rna_count, n), head(L.rna_count_var, n)};
        }
        return out;
    }
};

// The cube delivery (`CubeOut`'s arrays): one row per delivered AMBIG node, at most one per node.
struct CubeStore {
    std::vector<int64_t> slot;
    std::vector<double> pos, neg, total, opportunity, rho;   // the two profiles (d, K), the rest (d,)
    std::vector<uint8_t> has_pos, has_neg;

    CubeOut out(int n, int K) {
        const size_t nk = static_cast<size_t>(n) * K;
        grow(slot, n); grow(pos, nk); grow(neg, nk); grow(has_pos, n); grow(has_neg, n);
        for (auto* v : {&total, &opportunity, &rho}) grow(*v, n);
        return CubeOut{slot.data(), pos.data(), as_bool(has_pos), neg.data(), as_bool(has_neg), total.data(),
                       opportunity.data(), rho.data(), n};
    }

    // The first d rows.
    CubeStore prefix(int d, int K) const {
        auto head = [](const auto& v, size_t m) { return std::decay_t<decltype(v)>(v.begin(), v.begin() + m); };
        const size_t dk = static_cast<size_t>(d) * K;
        return CubeStore{head(slot, d), head(pos, dk), head(neg, dk), head(total, d), head(opportunity, d),
                         head(rho, d), head(has_pos, d), head(has_neg, d)};
    }
};

}  // namespace transfer_kernel
