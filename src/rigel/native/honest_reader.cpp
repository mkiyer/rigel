/**
 * honest_reader.cpp — the standalone binding of the honest capture reader (honest_reader.h): the per-slot factors
 * arrive packed as arrays (the transfer harness produces them, the same inputs the NumPy prototype consumed) and
 * the kernel returns the modes. This is the gates' entry point; the block solve (solve_kernel.cpp) runs the same
 * evaluator in place after the last refit.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <thread>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "honest_reader.h"

namespace nb = nanobind;

namespace {

using Vec = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using Mat = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using U8Vec = nb::ndarray<uint8_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using namespace rigel::honest;

struct Packed {
    int n_slots, n_lam, n_dna, n_pos, n_neg;
    const double *u, *v, *eg, *er, *q;
    const uint8_t *both, *admits, *has_dna, *has_pos, *has_neg;
    const double *lam, *comp;
    const double *dna_grid, *dna_vals, *dna_origin, *pos_grid, *pos_vals, *pos_origin, *neg_grid, *neg_vals, *neg_origin;

    Slot slot(int i) const {
        Slot s;
        s.u = u[i]; s.v = v[i]; s.eg = eg[i]; s.er = er[i]; s.q = q[i];
        s.both = both[i] != 0; s.admits = admits[i] != 0;
        s.comp = comp + size_t(i) * n_lam;
        auto level = [&](const uint8_t* has, const double* grid, const double* vals, const double* origin, int width) {
            Level L;
            L.present = has[i] != 0;
            if (L.present) {
                L.grid = grid + size_t(i) * width; L.vals = vals + size_t(i) * width; L.origin = origin[i];
                int n = 0;
                while (n < width && !std::isnan(L.grid[n])) ++n;  // NaN-padded rows
                L.n = n;
            }
            return L;
        };
        s.dna = level(has_dna, dna_grid, dna_vals, dna_origin, n_dna);
        s.pos = level(has_pos, pos_grid, pos_vals, pos_origin, n_pos);
        s.neg = level(has_neg, neg_grid, neg_vals, neg_origin, n_neg);
        return s;
    }
};

Packed unpack(Vec& u, Vec& v, Vec& eg, Vec& er, Vec& q, U8Vec& both, U8Vec& admits, Vec& lam, Mat& comp,
              Mat& dna_grid, Mat& dna_vals, Vec& dna_origin, U8Vec& has_dna,
              Mat& pos_grid, Mat& pos_vals, Vec& pos_origin, U8Vec& has_pos,
              Mat& neg_grid, Mat& neg_vals, Vec& neg_origin, U8Vec& has_neg) {
    Packed P;
    P.n_slots = int(u.shape(0)); P.n_lam = int(lam.shape(0));
    if (int(comp.shape(0)) != P.n_slots || int(comp.shape(1)) != P.n_lam) throw std::invalid_argument("comp must be (n_slots, n_lam)");
    P.n_dna = int(dna_grid.shape(1)); P.n_pos = int(pos_grid.shape(1)); P.n_neg = int(neg_grid.shape(1));
    P.u = u.data(); P.v = v.data(); P.eg = eg.data(); P.er = er.data(); P.q = q.data();
    P.both = both.data(); P.admits = admits.data(); P.has_dna = has_dna.data(); P.has_pos = has_pos.data(); P.has_neg = has_neg.data();
    P.lam = lam.data(); P.comp = comp.data();
    P.dna_grid = dna_grid.data(); P.dna_vals = dna_vals.data(); P.dna_origin = dna_origin.data();
    P.pos_grid = pos_grid.data(); P.pos_vals = pos_vals.data(); P.pos_origin = pos_origin.data();
    P.neg_grid = neg_grid.data(); P.neg_vals = neg_vals.data(); P.neg_origin = neg_origin.data();
    return P;
}

}  // namespace

/// Every slot's posterior mode of log density under the landscape (log_rho, logP); writes `mode` and
/// `max_post` (n_slots each). Threaded over slots.
void honest_modes(Vec u, Vec v, Vec eg, Vec er, Vec q, U8Vec both, U8Vec admits, Vec lam, Mat comp,
                  Mat dna_grid, Mat dna_vals, Vec dna_origin, U8Vec has_dna,
                  Mat pos_grid, Mat pos_vals, Vec pos_origin, U8Vec has_pos,
                  Mat neg_grid, Mat neg_vals, Vec neg_origin, U8Vec has_neg,
                  Vec land_logrho, Vec land_logP, Vec mode, Vec max_post,
                  double coarse_step, double window_sd, int tilt_nodes, bool use_coarse, bool ystar_slot_q, double refine,
                  bool derived_tilt, bool skip_far_tilt, int search_stride, int n_threads) {
    Packed P = unpack(u, v, eg, er, q, both, admits, lam, comp, dna_grid, dna_vals, dna_origin, has_dna,
                      pos_grid, pos_vals, pos_origin, has_pos, neg_grid, neg_vals, neg_origin, has_neg);
    int nx = int(land_logrho.shape(0));
    if (int(land_logP.shape(0)) != nx) throw std::invalid_argument("landscape grid and logP differ in length");
    if (int(mode.shape(0)) != P.n_slots || int(max_post.shape(0)) != P.n_slots) throw std::invalid_argument("mode/max_post must be (n_slots,)");
    Rules R{coarse_step, window_sd, tilt_nodes, use_coarse, ystar_slot_q, refine, derived_tilt, skip_far_tilt, search_stride};
    std::vector<double> rho(nx);
    const double* x = land_logrho.data();
    for (int i = 0; i < nx; ++i) rho[i] = std::exp(x[i]);
    const double* logP = land_logP.data();
    double* out_mode = mode.data();
    double* out_max = max_post.data();
    std::vector<int> prior_peaks = prior_peaks_of(logP, nx);
    int t = n_threads > 0 ? n_threads : int(std::thread::hardware_concurrency());
    t = std::max(1, std::min(t, P.n_slots));
    {
        nb::gil_scoped_release release;
        std::vector<std::thread> pool;
        for (int w = 0; w < t; ++w) {
            pool.emplace_back([&, w]() {
                std::vector<double> logL(nx);
                for (int i = w; i < P.n_slots; i += t) {
                    Slot s = P.slot(i);
                    out_mode[i] = slot_mode(s, P.lam, P.n_lam, R, rho.data(), x, logP, nx, prior_peaks, logL, &out_max[i]);
                }
            });
        }
        for (auto& th : pool) th.join();
    }
}

/// One slot's log L curve on the landscape grid (diagnostics and the evaluator gate).
void honest_logL(Vec u, Vec v, Vec eg, Vec er, Vec q, U8Vec both, U8Vec admits, Vec lam, Mat comp,
                 Mat dna_grid, Mat dna_vals, Vec dna_origin, U8Vec has_dna,
                 Mat pos_grid, Mat pos_vals, Vec pos_origin, U8Vec has_pos,
                 Mat neg_grid, Mat neg_vals, Vec neg_origin, U8Vec has_neg,
                 Vec land_logrho, int slot, Vec out,
                 double coarse_step, double window_sd, int tilt_nodes, bool use_coarse, bool ystar_slot_q, double refine,
                 bool derived_tilt, bool skip_far_tilt) {
    Packed P = unpack(u, v, eg, er, q, both, admits, lam, comp, dna_grid, dna_vals, dna_origin, has_dna,
                      pos_grid, pos_vals, pos_origin, has_pos, neg_grid, neg_vals, neg_origin, has_neg);
    int nx = int(land_logrho.shape(0));
    if (int(out.shape(0)) != nx) throw std::invalid_argument("out must match the grid");
    if (slot < 0 || slot >= P.n_slots) throw std::invalid_argument("slot out of range");
    Rules R{coarse_step, window_sd, tilt_nodes, use_coarse, ystar_slot_q, refine, derived_tilt, skip_far_tilt, 0};
    std::vector<double> rho(nx);
    for (int i = 0; i < nx; ++i) rho[i] = std::exp(land_logrho.data()[i]);
    Slot s = P.slot(slot);
    Evaluator ev(s, P.lam, P.n_lam, R);
    ev.log_L(rho.data(), nx, out.data());
}

NB_MODULE(_honest_reader_impl, m) {
    m.doc() = "The honest capture reader: per-slot posterior modes of gDNA density under the landscape (standalone, the gates' entry point).";
    m.def("honest_modes", &honest_modes, nb::arg("u"), nb::arg("v"), nb::arg("eg"), nb::arg("er"), nb::arg("q"), nb::arg("both"), nb::arg("admits"),
          nb::arg("lam"), nb::arg("comp"), nb::arg("dna_grid"), nb::arg("dna_vals"), nb::arg("dna_origin"), nb::arg("has_dna"),
          nb::arg("pos_grid"), nb::arg("pos_vals"), nb::arg("pos_origin"), nb::arg("has_pos"),
          nb::arg("neg_grid"), nb::arg("neg_vals"), nb::arg("neg_origin"), nb::arg("has_neg"),
          nb::arg("land_logrho"), nb::arg("land_logP"), nb::arg("mode"), nb::arg("max_post"),
          nb::arg("coarse_step"), nb::arg("window_sd"), nb::arg("tilt_nodes"), nb::arg("use_coarse"), nb::arg("ystar_slot_q"), nb::arg("refine"),
          nb::arg("derived_tilt"), nb::arg("skip_far_tilt"), nb::arg("search_stride"), nb::arg("n_threads"));
    m.def("honest_logL", &honest_logL, nb::arg("u"), nb::arg("v"), nb::arg("eg"), nb::arg("er"), nb::arg("q"), nb::arg("both"), nb::arg("admits"),
          nb::arg("lam"), nb::arg("comp"), nb::arg("dna_grid"), nb::arg("dna_vals"), nb::arg("dna_origin"), nb::arg("has_dna"),
          nb::arg("pos_grid"), nb::arg("pos_vals"), nb::arg("pos_origin"), nb::arg("has_pos"),
          nb::arg("neg_grid"), nb::arg("neg_vals"), nb::arg("neg_origin"), nb::arg("has_neg"),
          nb::arg("land_logrho"), nb::arg("slot"), nb::arg("out"),
          nb::arg("coarse_step"), nb::arg("window_sd"), nb::arg("tilt_nodes"), nb::arg("use_coarse"), nb::arg("ystar_slot_q"), nb::arg("refine"),
          nb::arg("derived_tilt"), nb::arg("skip_far_tilt"));
}
