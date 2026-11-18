// testing/solve.cpp — `_solve_impl.testing`: the calibration solve's kernels bound for the gates, which production never
// calls. Each runs the production kernel on arrays the call receives or allocates and returns.
//
//   psi_cube, posterior_median, compose, gdna_arm
//                    ψ's pieces (`psi_kernel.h`).
//   transfer_prepare, transfer_pass, transfer_solve
//                    the transfer policy's three kernels (`transfer_kernel.h`) on tables the call allocates and returns,
//                    or the harness builds by hand (`tests/calibration/_transfer_harness.py`).
//   the row constructors of `transfer_rows.h` and the blur's radius, the builders' flag predicates and the strand
//   evidence, each returning a fresh array; and the constants the gates read (EPS, MARGINAL_NODES, OWN_EVIDENCE_EPS,
//   JEFFREYS_REF, and RULE_CODES, the composition rules' codes by name).

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "constants.h"
#include "ndarray_util.h"
#include "psi_kernel.h"
#include "solve_bindings.h"
#include "transfer_kernel.h"
#include "transfer_store.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;
using namespace transfer_kernel;
using namespace solve_bindings;
using rigel::bool_ndarray;
using rigel::vec_to_ndarray;

namespace {

size_t sz(int x) { return static_cast<size_t>(x); }
int K_of(const Vec& v) { return static_cast<int>(v.shape(0)); }

// ---- ψ's pieces ----

// ψ itself: the cube (m, K, C) with the two strand-fraction grids and the tilt beside it, for slots of one class
// (ambig: n_tilt + 2 columns; else one).
void psi_cube(Vec u_pos, Vec u_neg, BoolVec allow_pos, BoolVec allow_neg, Vec fg_ref, Vec fpos_ref, Vec fneg_ref,
              double kappa, Vec lam, nb::object gdna, nb::object row_logprior, IdxVec cube_slot, Mat cube_pos,
              BoolVec cube_has_pos, Mat cube_neg, BoolVec cube_has_neg, Vec cube_u, Vec cube_total,
              Vec cube_opportunity, Vec cube_rho, int n_tilt, bool ambig, Cube out_psi, Cube out_fpos, Cube out_fneg,
              Cube out_tau) {
    const PsiInputs in(u_pos, u_neg, allow_pos, allow_neg, fg_ref, fpos_ref, fneg_ref, kappa, lam, gdna, row_logprior,
                       cube_slot, cube_pos, cube_has_pos, cube_neg, cube_has_neg, cube_u, cube_total,
                       cube_opportunity, cube_rho, n_tilt);
    const int K = in.K(), C = in.grid.columns(ambig);
    if (static_cast<int>(out_psi.shape(1)) != K || static_cast<int>(out_psi.shape(2)) != C)
        throw std::invalid_argument("psi_cube: the output is not (m, K, columns)");
    std::vector<double> scratch;
    for (int i = 0; i < in.slots(); ++i) {
        const size_t off = sz(i) * K * C;
        P::slot_cube(in.grid, in.slot(i), ambig, out_psi.data() + off, out_fpos.data() + off, out_fneg.data() + off,
                     out_tau.data() + off, scratch);
    }
}

void posterior_median(Mat post, Vec lam, Vec out) {
    const int m = static_cast<int>(post.shape(0)), K = static_cast<int>(post.shape(1));
    std::vector<double> cdf;
    for (int i = 0; i < m; ++i) out.data()[i] = P::posterior_median(post.data() + sz(i) * K, lam.data(), K, cdf);
}

// The fitted gDNA arm as an (m, K) matrix: the kernel's own construction at every (slot, cell).
void gdna_arm(Vec log_rho, Vec logP, Vec lam, Vec mass, Vec eff, Mat out) {
    const int m = static_cast<int>(mass.shape(0)), K = K_of(lam);
    if (static_cast<int>(out.shape(0)) != m || static_cast<int>(out.shape(1)) != K)
        throw std::invalid_argument("gdna_arm: the output is not (m, K)");
    if (logP.shape(0) != log_rho.shape(0) || log_rho.shape(0) < 1 || static_cast<int>(eff.shape(0)) != m)
        throw std::invalid_argument("gdna_arm: the curve is two (G,) arrays and the support one (mass, eff) per slot");
    const P::Arm arm{log_rho.data(), logP.data(), static_cast<int>(log_rho.shape(0)), mass.data(), eff.data()};
    std::vector<double> log_frac(K);
    for (int k = 0; k < K; ++k) log_frac[k] = P::arm_abscissa(sigmoid(lam.data()[k]));
    for (int i = 0; i < m; ++i) {
        const double sh = arm.shift(i);
        double* o = out.data() + sz(i) * K;
        for (int k = 0; k < K; ++k) o[k] = arm.at(log_frac[k], sh);
    }
}

void compose(Vec f_g, Vec w_pos, BoolVec allow_pos, BoolVec allow_neg, Vec out_fpos, Vec out_fneg) {
    for (size_t i = 0; i < f_g.shape(0); ++i)
        P::compose(f_g.data()[i], w_pos.data()[i], allow_pos.data()[i], allow_neg.data()[i], out_fpos.data()[i],
                   out_fneg.data()[i]);
}

// ---- the transfer kernels, on tables the call allocates and returns ----

const char* const FACE_PARAMS[9] = {"n_u", "n_s", "a_b", "a_x", "width", "var", "a_r_b", "a_r_x", "splice_rate"};
const char* const LANE_NAMES[3] = {"gdna", "pos", "neg"};

// One lane of the builders' tables as the harness reads it.
nb::dict lane_tables(BuilderTables::Lane&& L, int l, int n, int K, const ChainSums& sums, const ChainArrays& obs,
                     nb::object n_u, double rho_ref, double kappa, bool split_live) {
    nb::dict d;
    d["field"] = l;
    d["face"] = bool_ndarray(std::move(L.face), {sz(n), 2});
    d["two_sided"] = bool_ndarray(std::move(L.two_sided), {sz(n), 2});
    const double* a = l == 0 ? obs.a_g : obs.a_r;
    std::vector<uint8_t> empty(n);
    for (int i = 0; i < n; ++i) empty[i] = !(sums.n_u[i] > 0.0) || !(a[i] > 0.0);
    d["empty"] = bool_ndarray(std::move(empty));
    d["own_rows"] = vec_to_ndarray(std::move(L.own_rows), {sz(n), sz(K)});
    d["own_mask"] = bool_ndarray(std::move(L.own_mask));
    d["rho_ref"] = rho_ref;
    d["total"] = n_u;
    if (l == 0) {
        d["count"] = n_u;
        d["a"] = vec_to_ndarray(std::vector<double>(a, a + n));
        d["other"] = nb::none();
    } else {
        const int col_read = read_column(l - 1, kappa);
        d["count"] = vec_to_ndarray(std::vector<double>(sums.cnt_col[col_read]));
        d["a"] = vec_to_ndarray(std::vector<double>(a, a + n));
        d["other"] = split_live ? vec_to_ndarray(std::vector<double>(sums.cnt_col[1 - col_read])) : nb::none();
    }
    L.flux_store.resize(sz(L.n_flux_rows) * K);
    d["flux_rows"] = vec_to_ndarray(std::move(L.flux_store), {sz(L.n_flux_rows), sz(K)});
    d["flux_index"] = vec_to_ndarray(std::move(L.flux_index), {sz(n), 2});
    d["witness"] = vec_to_ndarray(std::move(L.witness), {sz(n), 2});
    d["witness_mask"] = bool_ndarray(std::move(L.witness_mask));
    return d;
}

// One block's tables as the gates read them: a dict of fresh arrays.
nb::dict transfer_prepare(Vec lam, BoolVec is_boundary, BoolVec is_exon, BoolVec free_pos, BoolVec free_neg,
                          BoolVec exon_pos, BoolVec exon_neg, BoolVec terminal, IdxVec left, IdxVec right,
                          FlagVec flags, Mat cnt, Mat spliced, Mat sj_count, Mat sj_count_lo, Mat sj_count_hi,
                          Mat route_rate_lo, Mat route_rate_hi, Vec eff_gdna, Vec eff_rna,
                          BoolVec has_own_composition, double kappa, double rho_gdna, double rho_rna,
                          bool split_live) {
    const int n = static_cast<int>(free_pos.shape(0)), K = K_of(lam);
    const ChainArrays obs{is_boundary.data(), is_exon.data(), free_pos.data(), free_neg.data(), exon_pos.data(),
                          exon_neg.data(), terminal.data(), flags.data(), cnt.data(), spliced.data(), sj_count.data(),
                          sj_count_lo.data(), sj_count_hi.data(), route_rate_lo.data(), route_rate_hi.data(),
                          eff_gdna.data(), eff_rna.data()};
    ChainSums sums;
    sums.fill(obs, n);
    const Chain c = chain_of(obs, sums, n, left.data(), right.data(), has_own_composition.data(), lam.data(), K, kappa);
    BuilderTables T;
    Scratch S(K);
    T.build(c, rho_gdna, rho_rna, S);

    nb::dict out;
    out["lam"] = vec_to_ndarray(std::vector<double>(lam.data(), lam.data() + K));
    nb::object n_u = vec_to_ndarray(std::vector<double>(sums.n_u));
    out["n_u"] = n_u;
    nb::dict own;
    own["rows"] = vec_to_ndarray(std::move(T.own_rows), {sz(n), sz(K)});
    own["mask"] = bool_ndarray(std::move(T.own_mask));
    out["own"] = own;
    nb::dict faces;
    faces["kind"] = vec_to_ndarray(std::move(T.f_kind), {sz(n), 2});
    faces["row"] = vec_to_ndarray(std::move(T.f_row), {sz(n), 2});
    faces["row2"] = vec_to_ndarray(std::move(T.f_row2), {sz(n), 2});
    for (int q = 0; q < 9; ++q) faces[FACE_PARAMS[q]] = vec_to_ndarray(std::move(T.f_par[q]), {sz(n), 2});
    T.face_store.resize(sz(T.n_face_rows) * K);
    faces["rows"] = vec_to_ndarray(std::move(T.face_store), {sz(T.n_face_rows), sz(K)});
    out["faces"] = faces;
    nb::dict lanes;
    for (int l = rho_gdna > 0.0 ? 0 : 1; l < 3; ++l)
        lanes[LANE_NAMES[l]] = lane_tables(std::move(T.lanes[l]), l, n, K, sums, obs, n_u,
                                           l == 0 ? rho_gdna : rho_rna, kappa, split_live);
    out["lanes"] = lanes;
    return out;
}

// A received table's dict (the harness's `Received.arrays()`), as the pass writes it in place.
struct ReceivedArrays {
    BoolVec has_nbr, has_comp;
    Mat comp;
    struct L { BoolVec present, has_wit; Mat profile; Vec count, opp, rna, rnav; } lanes[3];
    ReceivedView view;
    explicit ReceivedArrays(nb::dict d) {
        has_nbr = exact<BoolVec>(d["has_neighbour"]);
        has_comp = exact<BoolVec>(d["has_composition"]);
        comp = exact<Mat>(d["composition"]);
        const char* names[3] = {"level_gdna", "level_rna_pos", "level_rna_neg"};
        view.has_neighbour = has_nbr.data(); view.has_composition = has_comp.data(); view.composition = comp.data();
        view.K = static_cast<int>(comp.shape(1));
        for (int l = 0; l < 3; ++l) {
            nb::dict ld = nb::cast<nb::dict>(d[names[l]]);
            L& x = lanes[l];
            x.present = exact<BoolVec>(ld["present"]); x.profile = exact<Mat>(ld["profile"]);
            x.count = exact<Vec>(ld["count"]); x.opp = exact<Vec>(ld["opportunity"]);
            x.has_wit = exact<BoolVec>(ld["has_witness"]); x.rna = exact<Vec>(ld["rna_count"]);
            x.rnav = exact<Vec>(ld["rna_count_var"]);
            view.lanes[l] = LevelsView{x.present.data(), x.profile.data(), x.count.data(), x.opp.data(),
                                       x.has_wit.data(), x.rna.data(), x.rnav.data(), view.K};
        }
    }
};

// The tables dict `transfer_prepare` returns, or the harness builds by hand, as the pass and the solve read it.
struct TablesArrays {
    Vec lam; Mat own_rows; BoolVec own_mask;
    KindMat f_kind; IdxMat f_row, f_row2; Mat f_par[9]; Mat f_rows;
    struct L {
        bool built = false; int field = 0; BoolMat face, two; BoolVec empty, own_mask, wit_mask;
        Mat own_rows, flux_rows, wit; IdxMat flux_index; Vec count, a, other, total; bool has_other = false;
        double rho_ref = 0.0;
    } lanes[3];
    int n = 0, K = 0;
    explicit TablesArrays(nb::dict t) {
        lam = nb::cast<Vec>(t["lam"]); K = K_of(lam);
        nb::dict o = nb::cast<nb::dict>(t["own"]);
        own_rows = nb::cast<Mat>(o["rows"]); own_mask = nb::cast<BoolVec>(o["mask"]);
        n = static_cast<int>(own_mask.shape(0));
        nb::dict f = nb::cast<nb::dict>(t["faces"]);
        f_kind = nb::cast<KindMat>(f["kind"]); f_row = nb::cast<IdxMat>(f["row"]); f_row2 = nb::cast<IdxMat>(f["row2"]);
        for (int q = 0; q < 9; ++q) f_par[q] = nb::cast<Mat>(f[FACE_PARAMS[q]]);
        f_rows = nb::cast<Mat>(f["rows"]);
        nb::dict ls = nb::cast<nb::dict>(t["lanes"]);
        for (int l = 0; l < 3; ++l) {
            if (!ls.contains(LANE_NAMES[l])) continue;
            nb::dict d = nb::cast<nb::dict>(ls[LANE_NAMES[l]]);
            L& x = lanes[l];
            x.built = true; x.field = nb::cast<int>(d["field"]);
            x.face = nb::cast<BoolMat>(d["face"]); x.two = nb::cast<BoolMat>(d["two_sided"]);
            x.empty = nb::cast<BoolVec>(d["empty"]);
            x.own_rows = nb::cast<Mat>(d["own_rows"]); x.own_mask = nb::cast<BoolVec>(d["own_mask"]);
            x.count = nb::cast<Vec>(d["count"]); x.a = nb::cast<Vec>(d["a"]); x.total = nb::cast<Vec>(d["total"]);
            x.has_other = !d["other"].is_none();
            if (x.has_other) x.other = nb::cast<Vec>(d["other"]);
            x.flux_rows = nb::cast<Mat>(d["flux_rows"]); x.flux_index = nb::cast<IdxMat>(d["flux_index"]);
            x.wit = nb::cast<Mat>(d["witness"]); x.wit_mask = nb::cast<BoolVec>(d["witness_mask"]);
            x.rho_ref = nb::cast<double>(d["rho_ref"]);
        }
    }
    FacesView faces() const {
        return FacesView{f_kind.data(), f_row.data(), f_row2.data(), f_par[0].data(), f_par[1].data(), f_par[2].data(),
                         f_par[3].data(), f_par[4].data(), f_par[5].data(), f_par[6].data(), f_par[7].data(),
                         f_par[8].data(), f_rows.data(), K};
    }
    std::vector<LaneView> lane_views() const {
        std::vector<LaneView> v;
        for (const L& x : lanes) {
            if (!x.built) continue;
            v.push_back(LaneView{x.field, x.face.data(), x.two.data(), x.empty.data(), x.own_rows.data(),
                                 x.own_mask.data(), x.count.data(), x.a.data(), x.has_other ? x.other.data() : nullptr,
                                 x.wit.data(), x.wit_mask.data()});
        }
        return v;
    }
    SolveLane solve_lane(int l) const {
        const L& x = lanes[l];
        SolveLane s;
        if (!x.built) return s;
        s.built = true; s.field = x.field; s.empty = x.empty.data(); s.total = x.total.data(); s.a = x.a.data();
        s.rho_ref = x.rho_ref; s.own_rows = x.own_rows.data(); s.own_mask = x.own_mask.data();
        s.flux_store = x.flux_rows.data(); s.flux_index = x.flux_index.data();
        return s;
    }
};

// One directional pass on the tables, writing the received table in place: the gates' whole passes and single hops.
void transfer_pass(nb::dict tables, nb::dict received, IdxVec seq, IdxVec nbr, BoolVec terminal) {
    TablesArrays T(tables);
    ReceivedArrays R(received);
    if (R.view.K != T.K)
        throw std::invalid_argument("transfer_pass: the tables and the received table disagree on the grid");
    Scratch S(T.K);
    std::vector<double> out;
    pass_block(T.lam.data(), T.K, seq.data(), static_cast<int>(seq.shape(0)), nbr.data(), terminal.data(),
               T.own_rows.data(), T.own_mask.data(), T.faces(), T.lane_views(), R.view, S, out);
}

// The two received tables into ψ's channels: (live, rows (n, K), cube dict or None).
nb::tuple transfer_solve(nb::dict tables, nb::dict from_left, nb::dict from_right, BoolVec free_pos, BoolVec free_neg) {
    TablesArrays T(tables);
    ReceivedArrays Lft(from_left), Rgt(from_right);
    const int n = T.n, K = T.K;
    const ReceivedView sides[2] = {Lft.view, Rgt.view};
    const SolveLane G = T.solve_lane(0), rna[2] = {T.solve_lane(1), T.solve_lane(2)};
    std::vector<double> rows(sz(n) * K, 0.0), bound, row;
    CubeStore store;
    CubeOut cube = store.out(n, K);
    Scratch S(K);
    const auto [live, d] = solve_block(T.lam.data(), K, n, sides, G, rna, free_pos.data(), free_neg.data(),
                                       rows.data(), cube, S, bound, row);
    nb::object cube_out = nb::none();
    if (d) {
        CubeStore c = store.prefix(d, K);
        nb::dict cd;
        cd["slot"] = vec_to_ndarray(std::move(c.slot));
        cd["profile_pos"] = vec_to_ndarray(std::move(c.pos), {sz(d), sz(K)});
        cd["has_pos"] = bool_ndarray(std::move(c.has_pos));
        cd["profile_neg"] = vec_to_ndarray(std::move(c.neg), {sz(d), sz(K)});
        cd["has_neg"] = bool_ndarray(std::move(c.has_neg));
        cd["u"] = vec_to_ndarray(std::vector<double>(T.lam.data(), T.lam.data() + K));
        cd["total"] = vec_to_ndarray(std::move(c.total));
        cd["opportunity"] = vec_to_ndarray(std::move(c.opportunity));
        cd["rho_ref"] = vec_to_ndarray(std::move(c.rho));
        cube_out = cd;
    }
    return nb::make_tuple(live, vec_to_ndarray(std::move(rows), {sz(n), sz(K)}), cube_out);
}

// ---- the row constructors, flag predicates and strand evidence, each on a fresh array ----

// A row constructor's output: a fresh (K,) row, written with its own scratch.
template <class Fill>
nb::object row_of(int K, Fill&& fill) {
    Scratch S(K);
    std::vector<double> out(K);
    fill(S, out.data());
    return vec_to_ndarray(std::move(out));
}

void bind_rows(nb::module_& t) {
    t.def("trigamma", &trigamma, nb::arg("x"), "zeta(2, x), the counting variance's one home.");
    t.def("count_logvar", &count_logvar, nb::arg("n"));
    t.def("hop_price", &hop_price, nb::arg("n_s"), nb::arg("a_s"), nb::arg("n_x"), nb::arg("a_x"));
    t.def("blur_radius", &blur_radius, nb::arg("v"), nb::arg("step"),
          "The blur kernel's half-width in grid steps at variance v: the footprint blur_row convolves over.");
    t.def("blur_row", [](Vec row, Vec lam, double v) {
        return row_of(K_of(row), [&](Scratch& S, double* out) { blur_row(row.data(), K_of(row), lam.data(), v, out, S.blur); });
    }, nb::arg("row"), nb::arg("lam"), nb::arg("v"));
    t.def("lower_side", [](Vec p) {
        return row_of(K_of(p), [&](Scratch&, double* out) { lower_side(p.data(), K_of(p), out); });
    }, nb::arg("profile"));
    t.def("face_map_lambda", [](Vec lam, double n_u, double a_g_b, double a_r_b, double e_g_e, double e_r_e, double s) {
        return row_of(K_of(lam), [&](Scratch&, double* out) {
            face_map_lambda(lam.data(), K_of(lam), n_u, a_g_b, a_r_b, e_g_e, e_r_e, s, out);
        });
    }, nb::arg("lam"), nb::arg("n_u"), nb::arg("a_g_b"), nb::arg("a_r_b"), nb::arg("e_g_e"), nb::arg("e_r_e"), nb::arg("s"));
    t.def("transport_row", [](Vec row, Vec lam, Vec map, double n_u, double n_s) {
        return row_of(K_of(lam), [&](Scratch& S, double* out) {
            transport_row(row.data(), lam.data(), K_of(lam), map.data(), n_u, n_s, S, out);
        });
    }, nb::arg("row"), nb::arg("lam"), nb::arg("lam_e_of_u"), nb::arg("n_u"), nb::arg("n_s"));
    t.def("splice_out_row", [](Vec row_e, Vec lam, double n_u, double n_s, double a_g_b, double a_g_e, double a_r_b,
                               double a_r_e, double splice_rate, Vec nodes) {
        return row_of(K_of(lam), [&](Scratch& S, double* out) {
            splice_out_row(row_e.data(), lam.data(), K_of(lam), n_u, n_s, a_g_b, a_g_e, a_r_b, a_r_e, splice_rate,
                           nodes.data(), K_of(nodes), S, out);
        });
    }, nb::arg("row_e"), nb::arg("lam"), nb::arg("n_u"), nb::arg("n_s"), nb::arg("a_g_b"), nb::arg("a_g_e"),
       nb::arg("a_r_b"), nb::arg("a_r_e"), nb::arg("splice_rate"), nb::arg("nodes"));
    t.def("level_map_lambda", [](Vec lam, double density_b, double opportunity_i, double total_i) {
        return row_of(K_of(lam), [&](Scratch&, double* out) {
            level_map_lambda(lam.data(), K_of(lam), density_b, opportunity_i, total_i, out);
        });
    }, nb::arg("lam"), nb::arg("density_b"), nb::arg("opportunity_i"), nb::arg("total_i"));
    t.def("level_row", [](Vec row_b, Vec lam, Vec map, double v) {
        return row_of(K_of(lam), [&](Scratch& S, double* out) {
            level_row(row_b.data(), lam.data(), K_of(lam), map.data(), v, S, out);
        });
    }, nb::arg("row_b"), nb::arg("lam"), nb::arg("lam_i_of_b"), nb::arg("v"));
    t.def("level_bound_row", [](Vec lam, double density_b, double opportunity_i, double total_i, double v) {
        return row_of(K_of(lam), [&](Scratch&, double* out) {
            level_bound_row(lam.data(), K_of(lam), density_b, opportunity_i, total_i, v, out);
        });
    }, nb::arg("lam"), nb::arg("density_b"), nb::arg("opportunity_i"), nb::arg("total_i"), nb::arg("v"));
    t.def("edge_level_row", [](Vec lam, double n_b, double n_e, double a_g_b, double a_g_e) {
        return row_of(K_of(lam), [&](Scratch&, double* out) {
            edge_level_row(lam.data(), K_of(lam), n_b, n_e, a_g_b, a_g_e, out);
        });
    }, nb::arg("lam"), nb::arg("n_b"), nb::arg("n_e"), nb::arg("a_g_b"), nb::arg("a_g_e"));
    t.def("poisson_level", [](Vec u, double n, double a, double rho_ref) {
        return row_of(K_of(u), [&](Scratch&, double* out) { poisson_level(u.data(), K_of(u), n, a, rho_ref, out); });
    }, nb::arg("u"), nb::arg("n"), nb::arg("a"), nb::arg("rho_ref"));
    t.def("level_of_profile", [](Vec row, Vec lam, Vec u, double n, double a, double rho_ref) {
        return row_of(K_of(lam), [&](Scratch& S, double* out) {
            level_of_profile(row.data(), lam.data(), u.data(), K_of(lam), n, a, rho_ref, S, out);
        });
    }, nb::arg("row"), nb::arg("lam"), nb::arg("u"), nb::arg("n"), nb::arg("a"), nb::arg("rho_ref"));
    t.def("rna_level_of_profile", [](Vec row, Vec lam, Vec u, double n, double a_r, double rho_ref) {
        return row_of(K_of(lam), [&](Scratch& S, double* out) {
            rna_level_of_profile(row.data(), lam.data(), u.data(), K_of(lam), n, a_r, rho_ref, S, out);
        });
    }, nb::arg("row"), nb::arg("lam"), nb::arg("u"), nb::arg("n"), nb::arg("a_r"), nb::arg("rho_ref"));
    t.def("profile_of_level", [](Vec p, Vec u, Vec lam, double n, double a, double rho_ref) {
        return row_of(K_of(lam), [&](Scratch& S, double* out) {
            profile_of_level(p.data(), u.data(), lam.data(), K_of(lam), n, a, rho_ref, S, out);
        });
    }, nb::arg("profile"), nb::arg("u"), nb::arg("lam"), nb::arg("n"), nb::arg("a"), nb::arg("rho_ref"));
    t.def("rna_row_of_level", [](Vec p, Vec u, Vec lam, double n, double a_r, double rho_ref) {
        return row_of(K_of(lam), [&](Scratch& S, double* out) {
            rna_row_of_level(p.data(), u.data(), lam.data(), K_of(lam), n, a_r, rho_ref, S, out);
        });
    }, nb::arg("profile"), nb::arg("u"), nb::arg("lam"), nb::arg("n"), nb::arg("a_r"), nb::arg("rho_ref"));
    t.def("flux_level", [](Vec u, double count, double rate, double rho_ref, double v) -> nb::object {
        bool claimed = false;
        nb::object row = row_of(K_of(u), [&](Scratch& S, double* out) {
            claimed = flux_level(u.data(), K_of(u), count, rate, rho_ref, v, S, out);
        });
        return claimed ? row : nb::none();
    }, nb::arg("u"), nb::arg("count"), nb::arg("rate"), nb::arg("rho_ref"), nb::arg("v") = 0.0,
       "The certified flux as a lower-sided RNA level, or None for a zero count.");
}

nb::object flank(int64_t x) { return x < 0 ? nb::none() : nb::cast(x); }

// The builders' flag predicates on the boundary flags, the strand evidence, and the constants the gates read.
void bind_predicates(nb::module_& t) {
    t.attr("EPS") = EPS;
    nb::dict rule_codes;
    rule_codes["none"] = transfer_rows::NONE;
    rule_codes["forward"] = transfer_rows::FORWARD;
    rule_codes["transport"] = transfer_rows::TRANSPORT;
    rule_codes["splice_out"] = transfer_rows::SPLICE_OUT;
    rule_codes["edge"] = transfer_rows::EDGE;
    rule_codes["level"] = transfer_rows::LEVEL;
    t.attr("RULE_CODES") = rule_codes;
    t.attr("MARGINAL_NODES") = vec_to_ndarray(std::vector<double>(MARGINAL_NODES, MARGINAL_NODES + N_MARGINAL_NODES));
    t.attr("OWN_EVIDENCE_EPS") = rigel::OWN_EVIDENCE_EPS;
    t.attr("JEFFREYS_REF") = P::JEFFREYS_REF;
    t.def("face_is_licensed", [](int64_t f, bool fp_e, bool fn_e, bool fp_i, bool fn_i) {
        return face_is_licensed(static_cast<int>(f), fp_e, fn_e, fp_i, fn_i);
    }, nb::arg("flags"), nb::arg("fp_e"), nb::arg("fn_e"), nb::arg("fp_i"), nb::arg("fn_i"));
    t.def("boundary_shares_strand", &boundary_shares_strand, nb::arg("fp_b"), nb::arg("fn_b"), nb::arg("fp_i"),
          nb::arg("fn_i"));
    t.def("outside_flank", [](int64_t f, int64_t left, int64_t right) {
        int64_t o, i;
        outside_flank(static_cast<int>(f), left, right, o, i);
        return nb::make_tuple(flank(o), flank(i));
    }, nb::arg("flags"), nb::arg("left"), nb::arg("right"), "(outside, inside) of a terminus boundary, or (None, None).");
    t.def("junction_exon_side", [](int64_t f, int64_t left, int64_t right) {
        return flank(junction_exon_side(static_cast<int>(f), left, right));
    }, nb::arg("flags"), nb::arg("left"), nb::arg("right"), "The flank on a junction's exon side, or None.");
    t.def("junction_flanks", [](int64_t f, int64_t left, int64_t right) {
        int64_t c, e;
        junction_flanks(static_cast<int>(f), left, right, c, e);
        return nb::make_tuple(flank(c), flank(e));
    }, nb::arg("flags"), nb::arg("left"), nb::arg("right"),
       "(C, E) at an exon|exon junction with no terminus, or (None, None).");
    // the strand evidence the block pipeline computes inline
    t.def("strand_evidence", [](Vec u_pos, Vec u_neg, Vec fg_loc, double kappa, double disc) {
        std::vector<double> out(u_pos.shape(0));
        for (size_t i = 0; i < out.size(); ++i)
            out[i] = strand_evidence(u_pos.data()[i], u_neg.data()[i], fg_loc.data()[i], kappa, disc);
        return vec_to_ndarray(std::move(out));
    }, nb::arg("u_pos"), nb::arg("u_neg"), nb::arg("fg_loc"), nb::arg("kappa"), nb::arg("disc"),
       "The reference-free strand composition evidence I_strand at fg_loc, `disc` the protocol's discriminability "
       "(0: the channel is dead).");
}

}  // namespace

void rigel::bind_testing(nb::module_ t) {
    t.def("psi_cube", &psi_cube, nb::arg("u_pos"), nb::arg("u_neg"), nb::arg("allow_pos"), nb::arg("allow_neg"),
          nb::arg("fg_ref"), nb::arg("fpos_ref"), nb::arg("fneg_ref"), nb::arg("kappa"), nb::arg("lam"),
          nb::arg("gdna").none(), nb::arg("row_logprior").none(), nb::arg("cube_slot"), nb::arg("cube_pos"),
          nb::arg("cube_has_pos"), nb::arg("cube_neg"), nb::arg("cube_has_neg"), nb::arg("cube_u"),
          nb::arg("cube_total"), nb::arg("cube_opportunity"), nb::arg("cube_rho"), nb::arg("n_tilt"), nb::arg("ambig"),
          nb::arg("out_psi").noconvert(), nb::arg("out_fpos").noconvert(), nb::arg("out_fneg").noconvert(),
          nb::arg("out_tau").noconvert(),
          "ψ over the (m, K, columns) cube for slots of one class, with f_pos, f_neg and the tilt beside it.");
    t.def("posterior_median", &posterior_median, nb::arg("post"), nb::arg("lam"), nb::arg("out").noconvert(),
          "The continuous ½-quantile of each row's λ posterior, read on λ's midpoint edges and mapped through σ.");
    t.def("compose", &compose, nb::arg("f_g"), nb::arg("w_pos"), nb::arg("allow_pos"), nb::arg("allow_neg"),
          nb::arg("out_fpos").noconvert(), nb::arg("out_fneg").noconvert(),
          "The composition as the image of (f_g, w_pos) on the admissible strands.");
    t.def("gdna_arm", &gdna_arm, nb::arg("log_rho"), nb::arg("logP"), nb::arg("lam"), nb::arg("mass"), nb::arg("eff"),
          nb::arg("out").noconvert(),
          "The fitted gDNA arm as (m, K): the landscape's curve at log σ(λ) + log M − log E per slot and cell, numpy's "
          "interpolation with the ends held — the kernel's own construction, written in place.");
    t.def("transfer_prepare", &transfer_prepare, nb::arg("lam"), nb::arg("is_boundary"), nb::arg("is_exon"),
          nb::arg("free_pos"), nb::arg("free_neg"), nb::arg("exon_pos"), nb::arg("exon_neg"), nb::arg("terminal"),
          nb::arg("left"), nb::arg("right"), nb::arg("flags"), nb::arg("cnt"), nb::arg("spliced"), nb::arg("sj_count"),
          nb::arg("sj_count_lo"), nb::arg("sj_count_hi"), nb::arg("route_rate_lo"), nb::arg("route_rate_hi"),
          nb::arg("eff_gdna"), nb::arg("eff_rna"), nb::arg("has_own_composition"), nb::arg("kappa"),
          nb::arg("rho_gdna"), nb::arg("rho_rna"), nb::arg("split_live"),
          "Build one block's own claims, face rules and level lanes into fresh tables, reading the caller's `terminal` "
          "as the intergenic class, and return them as a dict — `own` (rows, mask), `faces` (kind, row, row2, n_u, "
          "n_s, a_b, a_x, width, var, rows), `lanes` (gdna where the library has a gDNA coordinate, pos, neg: each its "
          "faces, two-sided faces, emptiness, own rows and mask, count, a, other, flux rows and index, witnesses), "
          "`lam`, `n_u`.");
    t.def("transfer_pass", &transfer_pass, nb::arg("tables"), nb::arg("received"), nb::arg("seq"), nb::arg("nbr"),
          nb::arg("terminal"),
          "Run one pass in chain order on the tables: for every destination in `seq` with a neighbour `nbr[i] >= 0` that "
          "is not a terminal, apply the face's composition rule and carry each lane's level, writing the received "
          "table's arrays in place.");
    t.def("transfer_solve", &transfer_solve, nb::arg("tables"), nb::arg("from_left"), nb::arg("from_right"),
          nb::arg("free_pos"), nb::arg("free_neg"),
          "The policy's solve for one block: the two received tables into the fused λ rows and the cube delivery at the "
          "AMBIG nodes. Returns (live, rows, cube dict or None).");
    bind_rows(t);
    bind_predicates(t);
}
