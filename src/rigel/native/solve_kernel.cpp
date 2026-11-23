// solve_kernel.cpp — the calibration solve's native module, `_solve_impl`. One call runs a whole sweep over the locus
// blocks on a thread pool, with ψ's kernel (`psi_kernel.h`) and the composition transfer's (`transfer_kernel.h`) as its
// pieces. The gates read the same kernels through `_solve_impl.testing` (testing/solve.cpp).
//
//   solve_blocks     the sweep. Per locus block: the introns' abundance rows, the self-solve ψ, the own-evidence
//                    precision, the message layer
//                    (builders, two passes, solve; skipped under the silent policy), the final ψ, the write-back on
//                    the owned slots, `has_composition`, and the counts of the backbone's assertions. Each block runs
//                    on one pool thread with that thread's arena; the belief is written in place and the diagnostics
//                    captured on request. Bit-identical at every thread count: a block's arithmetic does not depend on
//                    its thread, and only integer counts are reduced across blocks (`sweep.solve_chain`).
//   psi_solve        ψ over a slot list, one slot at a time on the same pool: the pre-sweep solve
//                    (`simplex_logodds.solve_regions_logodds`) and the gates (`tests/calibration/_psi_reference.py`).
//   rows             `read_column` and `ARM_EPS`, the kernel's values its Python readers share.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "constants.h"
#include "ndarray_util.h"
#include "psi_kernel.h"
#include "solve_bindings.h"
#include "thread_pool.h"
#include "transfer_kernel.h"
#include "transfer_store.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;
using namespace transfer_kernel;
using namespace solve_bindings;
using rigel::bool_ndarray;
using rigel::vec_to_ndarray;

namespace {

// ---- the module's thread pool ----

// One pool, shared by the block loop and ψ's slot loop. It persists across a run's calls, is rebuilt only when the
// thread count changes, and serialises the calls.
rigel::ThreadPool& pool_of(int n_threads) {
    static std::unique_ptr<rigel::ThreadPool> pool;
    if (!pool || pool->n_threads() != n_threads) pool = std::make_unique<rigel::ThreadPool>(n_threads);
    return *pool;
}
std::mutex& pool_mutex() {
    static std::mutex m;
    return m;
}

// The threads a call runs on: its budget (0: every core), at most one per task.
int pool_threads(int n_threads, int64_t n_tasks) {
    return static_cast<int>(std::min<int64_t>(rigel::thread_budget(n_threads), std::max<int64_t>(n_tasks, 1)));
}

// task(q, tid) for every q in [0, n_tasks), pulled one at a time by `threads` threads without the GIL. Pulling single
// tasks, with no chunk size, balances the load on asymmetric cores; `tid` names the thread's own state.
template <class Task>
void run_tasks(int threads, int64_t n_tasks, Task&& task) {
    std::atomic<int64_t> next{0};
    auto worker = [&](int tid) {
        for (int64_t q; (q = next.fetch_add(1, std::memory_order_relaxed)) < n_tasks;) task(q, tid);
    };
    nb::gil_scoped_release release;
    if (threads == 1) {
        worker(0);
        return;
    }
    std::lock_guard<std::mutex> lock(pool_mutex());
    pool_of(threads).run_parallel(worker);
}

// ---- ψ for the pre-sweep solve ----

// Each slot is solved on its own cube, sharing only the read-only grid and delivered rows, and writes its own four
// outputs, so the result is bit-identical to the serial loop at every thread count.
void psi_solve(IdxVec slots, Vec u_pos, Vec u_neg, BoolVec allow_pos, BoolVec allow_neg, Vec fg_ref, Vec fpos_ref,
               Vec fneg_ref, double kappa, Vec lam, nb::object gdna, nb::object row_logprior, IdxVec cube_slot,
               Mat cube_pos, BoolVec cube_has_pos, Mat cube_neg, BoolVec cube_has_neg, Vec cube_u, Vec cube_total,
               Vec cube_opportunity, Vec cube_rho, int n_tilt, Vec out_fg, Vec out_fpos, Vec out_fneg, Vec out_var,
               int n_threads) {
    if (n_tilt < 2) throw std::invalid_argument("psi_solve: the tilt needs at least two nodes");
    const PsiInputs in(u_pos, u_neg, allow_pos, allow_neg, fg_ref, fpos_ref, fneg_ref, kappa, lam, gdna, row_logprior,
                       cube_slot, cube_pos, cube_has_pos, cube_neg, cube_has_neg, cube_u, cube_total,
                       cube_opportunity, cube_rho, n_tilt);
    const int64_t n_sel = static_cast<int64_t>(slots.shape(0));
    const int T = pool_threads(n_threads, n_sel);
    std::vector<P::Scratch> scratch(T);
    run_tasks(T, n_sel, [&](int64_t q, int tid) {
        const int i = static_cast<int>(slots.data()[q]);
        P::solve_slot(in.grid, in.slot(i), scratch[tid], out_fg.data()[i], out_fpos.data()[i], out_fneg.data()[i],
                      out_var.data()[i]);
    });
}

// ---- the intron abundance model ----
// An intron's gDNA is one draw from the intergenic background's law at its own gDNA opportunity E, except with prior
// probability ε, when it is not at the background and any g in [0, C] is as likely. So its unspliced count C carries a
// row on the λ grid, log((1 − ε)·NegBinom(f_g·C; rate·E, size) + ε/(C + 1)), offset so its maximum is 0: the peak is
// the background's where the background is right, and the tail is capped, so the intron's own decisive strand evidence
// outweighs a background it contradicts. The row is read at the λ nodes, except that the node whose lattice cell holds
// the background's mode is read at the mode: a deep intron's peak, narrower than a cell, falls between two nodes, and
// read at the nodes alone it would sit under the departure floor at both and the row would be flat. The law is fitted
// once per calibration (`gdna_density.background_law`). The row joins ψ, its precision joins the intron's own
// composition evidence, and it is the intron's own claim to the message layer.

struct Background {
    double log_rate;   // the background rate's log, fragments per unit of gDNA opportunity
    double size;       // the negative binomial size an intron's gDNA is drawn with (∞: Poisson)
    double departure;  // ε, the prior probability that the intron's gDNA is not at the background
};

// log NegBinom(g; mean μ, size r) at a continuous g ≥ 0; r = ∞ is the Poisson limit
double log_negbinom(double g, double mu, double size) {
    mu = std::max(mu, P::ARM_EPS);
    const double lg_g1 = std::lgamma(g + 1.0);
    if (!std::isfinite(size)) {
        const double t = g * std::log(mu);
        return t - mu - lg_g1;
    }
    const double r = std::max(size, P::ARM_EPS);
    const double rpm = r + mu;
    const double t_r = r * (std::log(r) - std::log(rpm));
    const double t_g = g * (std::log(mu) - std::log(rpm));
    return std::lgamma(g + r) - std::lgamma(r) - lg_g1 + t_r + t_g;
}

// one intron's row over the uniform λ lattice: its gDNA count f_g·C against the background at its opportunity, at each
// node, and at the background's mode for the node whose cell holds it
void intron_row(const Background& bg, double count, double eff, const double* lam, int K, double* out) {
    const double mu = std::exp(bg.log_rate) * std::max(eff, P::ARM_EPS);
    // the continuous count's mode, where ψ(g + size) − ψ(g + 1) = log((μ + size)/μ), to O(1/g²); 0 below size 1
    const double mode = std::max(0.0, (std::isfinite(bg.size) ? mu * (1.0 - 1.0 / bg.size) : mu) - 0.5);
    const double half = K > 1 ? 0.5 * (lam[K - 1] - lam[0]) / (K - 1) : 0.0;
    for (int k = 0; k < K; ++k) {
        const double lo = P::sigmoid(lam[k] - half) * count, hi = P::sigmoid(lam[k] + half) * count;
        const double node = std::clamp(P::sigmoid(lam[k]), P::ARM_EPS, 1.0 - P::ARM_EPS) * count;
        out[k] = log_negbinom(mode >= lo && mode < hi ? mode : node, mu, bg.size);
    }
    // the departure: log((1 − ε)·NB + ε/(C + 1)), summed in the log domain
    const double a = std::log1p(-bg.departure), b = std::log(bg.departure) - std::log(count + 1.0);
    for (int k = 0; k < K; ++k) {
        const double x = a + out[k], hi = std::max(x, b);
        out[k] = hi + std::log(std::exp(x - hi) + std::exp(b - hi));
    }
    const double m = vmax(out, K);
    for (int k = 0; k < K; ++k) out[k] -= m;
}

// the composition evidence a row carries, 1/Var_λ under the row normalised; a flat row carries none
double row_precision(const double* row, const double* lam, int K) {
    if (!(ptp(row, K) > P::ARM_EPS)) return 0.0;
    const double m = vmax(row, K);
    double z = 0.0, s1 = 0.0, s2 = 0.0;
    for (int k = 0; k < K; ++k) {
        const double w = std::exp(row[k] - m);
        z += w;
        s1 += w * lam[k];
        s2 += w * lam[k] * lam[k];
    }
    const double mean = s1 / z, var = s2 / z - mean * mean;
    return var > P::ARM_EPS ? 1.0 / var : 0.0;
}

// ---- the block pipeline ----

// The backbone's assertions, counted per block on the owned slots as (violations, eligible); the backbone raises.
enum Assertion { LAM_ROWS_FINITE, CUBE_ROWS_FINITE, N_ASSERT };
const char* const ASSERTION_NAMES[N_ASSERT] = {"lam_rows_finite", "cube_rows_finite"};

struct Belief {  // the incoming belief, read-only
    const double *fpos, *fneg, *fg;
    Belief from(int64_t start) const { return Belief{fpos + start, fneg + start, fg + start}; }
};

struct Sweep {  // the whole chain and the call's settings, read-only
    ChainArrays obs;
    const int64_t *left, *right;  // global links
    Belief belief;
    double kappa, disc;  // ψ's strand model; the protocol's discriminability (0: the channel is dead)
    const double* lam; int K;
    // the policy, the library's gDNA and RNA lane coordinates, and whether the protocol permits strand witnesses
    int policy; double rho_gdna, rho_rna; bool split_live;
    const P::Grid* grid;
    bool want_capture;
    const Background* background;  // the intergenic background's law, or nullptr: no intron abundance rows
};

struct Outputs {  // written in place, the owned slots of each block
    double *fpos, *fneg, *fg, *var; bool* has_comp;
    // the diagnostics capture, or nullptr
    double *fg_loc = nullptr, *tau_lam = nullptr, *lam_rows = nullptr;
};

struct BlockResult {
    int64_t counts[N_ASSERT][2] = {};
    bool rows_delivered = false, cube_delivered = false, captured = false;
    // the capture: the cube rows the solve delivered, in the block's local slots, and the two received tables on the
    // owned slots
    CubeStore cube;
    ReceivedStore received[2];
};

// One thread's tables, sized per block; the capacity is kept across blocks.
struct Arena {
    int K = 0;
    ChainSums sums;
    std::vector<int64_t> left, right, seq_f, seq_b;   // the block's links (-1 off the block) and its two pass orders
    std::vector<double> fg_loc, tau_lam, rows, out, bound, row;
    std::vector<uint8_t> solvable, has_own, own_ev, held, empty_g, empty_r;
    std::vector<const double*> row_ptr;
    std::vector<double> intron;               // the block's intron abundance rows, (n, K)
    std::vector<const double*> intron_ptr;    // each slot's row, or nullptr
    BuilderTables builders;
    ReceivedStore received[2];
    CubeStore cube;
    std::unique_ptr<Scratch> S;  // transfer_rows' scratch, built for K
    P::Scratch psi;

    void size(int n, int K_) {
        if (K_ != K || !S) S = std::make_unique<Scratch>(K_);
        K = K_;
        for (auto* v : {&left, &right, &seq_f, &seq_b}) grow(*v, n);
        for (auto* v : {&fg_loc, &tau_lam}) grow(*v, n);
        grow(rows, static_cast<size_t>(n) * K);
        for (auto* v : {&solvable, &has_own, &own_ev, &held, &empty_g, &empty_r}) grow(*v, n);
        grow(row_ptr, n);
        grow(intron, static_cast<size_t>(n) * K);
        grow(intron_ptr, n);
    }
};

// One block's slice of the sweep: its observations, belief and size.
struct Block {
    ChainArrays obs;
    Belief belief;
    int64_t start;
    int n, n_owned;
};

// The block's links re-based (a neighbour outside the block is no neighbour), its two pass orders and its sums.
void link_block(const Sweep& W, const Block& b, Arena& A) {
    const int n = b.n;
    for (int i = 0; i < n; ++i) {
        const int64_t l = W.left[b.start + i] - b.start, r = W.right[b.start + i] - b.start;
        A.left[i] = (l >= 0 && l < n) ? l : -1;
        A.right[i] = (r >= 0 && r < n) ? r : -1;
        A.seq_f[i] = i;
        A.seq_b[i] = n - 1 - i;
    }
    A.sums.fill(b.obs, n);
}

// Slot i's inputs to ψ at the incoming belief, with the delivered composition and RNA level rows (or none).
P::SlotInputs slot_inputs(const Sweep& W, const Block& b, const Arena& A, int i, const double* row_prior,
                          const P::Delivered* delivered) {
    const bool built = W.grid->arm->built();  // the fitted gDNA arm
    return P::SlotInputs{b.obs.cnt[2 * i], b.obs.cnt[2 * i + 1], b.belief.fg[i], b.belief.fpos[i], b.belief.fneg[i],
                         b.obs.fp[i], b.obs.fn[i], built,
                         built ? P::Arm::shift_of(A.sums.n_u[i], b.obs.a_g[i]) : 0.0, row_prior, delivered,
                         A.intron_ptr[i]};
}

// Each intron's abundance row on the grid: a non-exon region that admits a strand, read against the background. No
// row where the background is off, nor at any other slot.
void intron_rows(const Sweep& W, const Block& b, Arena& A) {
    const int K = W.K;
    for (int i = 0; i < b.n; ++i) {
        const double* r = nullptr;
        if (W.background != nullptr && !b.obs.is_bnd[i] && !b.obs.is_exon[i] && (b.obs.fp[i] || b.obs.fn[i])) {
            double* dest = A.intron.data() + static_cast<size_t>(i) * K;
            intron_row(*W.background, b.obs.cnt[2 * i] + b.obs.cnt[2 * i + 1], b.obs.a_g[i], W.lam, K, dest);
            r = dest;
        }
        A.intron_ptr[i] = r;
    }
}

// The self-solve: every solvable slot (a live strand and an unspliced fragment) on its own cube, the incoming belief
// frozen as its reference, gives fg_loc; any other slot keeps the incoming belief.
void self_solve(const Sweep& W, const Block& b, Arena& A) {
    for (int i = 0; i < b.n; ++i) {
        as_bool(A.solvable)[i] = (b.obs.fp[i] || b.obs.fn[i]) && A.sums.n_u[i] > 0.0;
        A.fg_loc[i] = b.belief.fg[i];
        if (!as_bool(A.solvable)[i]) continue;
        double f_g, f_p, f_n, v_g;
        P::solve_slot(*W.grid, slot_inputs(W, b, A, i, nullptr, nullptr), A.psi, f_g, f_p, f_n, v_g);
        A.fg_loc[i] = f_g;
    }
}

// The own evidence: a single-strand slot's strand precision at fg_loc, and an intron's abundance-row precision. The
// builders' bit is the strand claim's; the abundance row reaches them as the intron's own row (`Chain::intron_row`).
void own_evidence(const Sweep& W, const Block& b, Arena& A) {
    for (int i = 0; i < b.n; ++i) {
        const bool single = b.obs.fp[i] != b.obs.fn[i];
        const double i_strand =
            strand_evidence(b.obs.cnt[2 * i], b.obs.cnt[2 * i + 1], A.fg_loc[i], W.kappa, W.disc);
        const double tau_strand = single ? i_strand : 0.0;
        const double tau_row = A.intron_ptr[i] != nullptr ? row_precision(A.intron_ptr[i], W.lam, W.K) : 0.0;
        const double tau = tau_strand + tau_row;
        A.tau_lam[i] = tau;
        as_bool(A.has_own)[i] = tau_strand > 0.0;
        as_bool(A.own_ev)[i] = tau > rigel::OWN_EVIDENCE_EPS;
    }
}

// The level lanes as the pass reads them: the gDNA lane where the library has a gDNA coordinate, and each strand's RNA
// lane on the column its reads land on, the other column its witness where the protocol permits.
std::vector<LaneView> pass_lanes(const Sweep& W, const Block& b, Arena& A) {
    const BuilderTables& T = A.builders;
    std::vector<LaneView> lanes;
    if (W.rho_gdna > 0.0) {
        const BuilderTables::Lane& L = T.lanes[0];
        lanes.push_back(LaneView{0, as_bool(L.face), as_bool(L.two_sided), as_bool(A.empty_g), L.own_rows.data(),
                                 as_bool(L.own_mask), A.sums.n_u.data(), b.obs.a_g, nullptr, L.witness.data(),
                                 as_bool(L.witness_mask)});
    }
    for (int s = 0; s < 2; ++s) {
        const BuilderTables::Lane& L = T.lanes[1 + s];
        const int col_read = read_column(s, W.kappa);
        lanes.push_back(LaneView{1 + s, as_bool(L.face), as_bool(L.two_sided), as_bool(A.empty_r), L.own_rows.data(),
                                 as_bool(L.own_mask), A.sums.cnt_col[col_read].data(), b.obs.a_r,
                                 W.split_live ? A.sums.cnt_col[1 - col_read].data() : nullptr, L.witness.data(),
                                 as_bool(L.witness_mask)});
    }
    return lanes;
}

// The lanes as the solve reads them: the gDNA lane (unbuilt without a gDNA coordinate) and the two RNA lanes with
// their own levels and junction flux levels.
void solve_lanes(const Sweep& W, const Block& b, Arena& A, SolveLane& G, SolveLane rna[2]) {
    if (W.rho_gdna > 0.0) {
        G.built = true; G.field = 0; G.empty = as_bool(A.empty_g); G.total = A.sums.n_u.data(); G.a = b.obs.a_g;
        G.rho_ref = W.rho_gdna;
    }
    for (int s = 0; s < 2; ++s) {
        const BuilderTables::Lane& L = A.builders.lanes[1 + s];
        rna[s].built = true; rna[s].field = 1 + s; rna[s].empty = as_bool(A.empty_r); rna[s].total = A.sums.n_u.data();
        rna[s].a = b.obs.a_r; rna[s].rho_ref = W.rho_rna; rna[s].own_rows = L.own_rows.data();
        rna[s].own_mask = as_bool(L.own_mask); rna[s].flux_store = L.flux_store.data();
        rna[s].flux_index = L.flux_index.data();
    }
}

// The message layer: the builders, the two passes into the two received tables (the forward pass reads each node's
// low neighbour, the backward its high one), and the solve into the fused λ rows and the cube delivery. Marks the
// slots that hold a composition from either side. Returns (rows live, cube rows delivered).
std::pair<bool, int> message_layer(const Sweep& W, const Block& b, Arena& A) {
    const int n = b.n, K = W.K;
    Scratch& S = *A.S;
    Chain c = chain_of(b.obs, A.sums, n, A.left.data(), A.right.data(), as_bool(A.has_own), W.lam, K, W.kappa);
    c.intron_row = A.intron_ptr.data();
    A.builders.build(c, W.rho_gdna, W.rho_rna, S);
    for (int i = 0; i < n; ++i) {
        as_bool(A.empty_g)[i] = !(A.sums.n_u[i] > 0.0) || !(b.obs.a_g[i] > 0.0);
        as_bool(A.empty_r)[i] = !(A.sums.n_u[i] > 0.0) || !(b.obs.a_r[i] > 0.0);
    }
    const std::vector<LaneView> lanes = pass_lanes(W, b, A);
    const FacesView faces = A.builders.faces();
    ReceivedView R[2];
    for (int s = 0; s < 2; ++s) {
        const int64_t* nbr = s == 0 ? A.left.data() : A.right.data();
        R[s] = A.received[s].view(n, K, nbr);
        pass_block(W.lam, K, (s == 0 ? A.seq_f : A.seq_b).data(), n, nbr, b.obs.terminal, A.builders.own_rows.data(),
                   as_bool(A.builders.own_mask), faces, lanes, R[s], S, A.out);
    }
    SolveLane G, rna[2];
    solve_lanes(W, b, A, G, rna);
    std::fill_n(A.rows.data(), static_cast<size_t>(n) * K, 0.0);
    CubeOut cube = A.cube.out(n, K);
    const auto delivered = solve_block(W.lam, K, n, R, G, rna, b.obs.fp, b.obs.fn, A.rows.data(), cube, S, A.bound,
                                       A.row);
    for (int i = 0; i < n; ++i) as_bool(A.held)[i] = R[0].has_composition[i] || R[1].has_composition[i];
    return delivered;
}

// The counts of the backbone's assertions on what was delivered, over the owned slots: every λ row finite, and every
// delivered cube row's profiles finite.
void count_assertions(const Block& b, const Arena& A, int K, bool has_rows, int d, const P::Rows& delivered,
                      BlockResult& res) {
    res.rows_delivered = has_rows;
    if (has_rows) {
        for (int i = 0; i < b.n_owned; ++i) {
            const double* r = A.rows.data() + static_cast<size_t>(i) * K;
            bool finite = true;
            for (int k = 0; k < K; ++k) finite = finite && std::isfinite(r[k]);
            res.counts[LAM_ROWS_FINITE][0] += !finite;
        }
        res.counts[LAM_ROWS_FINITE][1] += b.n_owned;
    }
    res.cube_delivered = d > 0;
    for (int q = 0; q < d; ++q) {
        if (A.cube.slot[q] >= b.n_owned) continue;
        const P::Delivered& row = delivered.rows[q];
        bool bad = false;
        for (const double* prof : {row.pos, row.neg}) {
            if (prof == nullptr) continue;
            for (int k = 0; k < K; ++k) bad = bad || !std::isfinite(prof[k]);
        }
        res.counts[CUBE_ROWS_FINITE][0] += bad;
        res.counts[CUBE_ROWS_FINITE][1] += 1;
    }
}

// The diagnostics on the owned slots: fg_loc, the own-evidence precision and, where asked, the fused λ rows (zero
// where none were delivered).
void capture_slots(const Block& b, const Arena& A, int K, bool has_rows, Outputs& out) {
    for (int i = 0; i < b.n_owned; ++i) {
        out.fg_loc[b.start + i] = A.fg_loc[i];
        out.tau_lam[b.start + i] = A.tau_lam[i];
        if (!out.lam_rows) continue;
        double* dest = out.lam_rows + static_cast<size_t>(b.start + i) * K;
        const double* src = A.rows.data() + static_cast<size_t>(i) * K;
        if (has_rows) std::copy(src, src + K, dest);
        else std::fill(dest, dest + K, 0.0);
    }
}

// The final solve of every solvable owned slot, with the fused λ row and its delivered cube row, written back
// clamped; a locked or empty slot keeps the incoming belief. Then `has_composition`: an own composition channel,
// structural certainty, or a composition held from either side.
void final_solve(const Sweep& W, const Block& b, Arena& A, const P::Rows& delivered, Outputs& out) {
    for (int i = 0; i < b.n_owned; ++i) {
        if (!as_bool(A.solvable)[i]) continue;
        const P::Delivered* row =
            delivered.index.empty() || delivered.index[i] < 0 ? nullptr : &delivered.rows[delivered.index[i]];
        double f_g, f_p, f_n, v_g;
        P::solve_slot(*W.grid, slot_inputs(W, b, A, i, A.row_ptr[i], row), A.psi, f_g, f_p, f_n, v_g);
        const int64_t at = b.start + i;
        out.fpos[at] = std::clamp(f_p, 0.0, 1.0);
        out.fneg[at] = std::clamp(f_n, 0.0, 1.0);
        out.fg[at] = std::clamp(f_g, 0.0, 1.0);
        out.var[at] = v_g;
    }
    for (int i = 0; i < b.n_owned; ++i)
        out.has_comp[b.start + i] = as_bool(A.own_ev)[i] || (!b.obs.fp[i] && !b.obs.fn[i]) || as_bool(A.held)[i];
}

// One locus block, end to end, on one thread: its slots [start, end), of which [start, stop) are its own.
void solve_one_block(const Sweep& W, int64_t start, int64_t stop, int64_t end, Outputs& out, Arena& A,
                     BlockResult& res) {
    const Block b{W.obs.from(start), W.belief.from(start), start, static_cast<int>(end - start),
                  static_cast<int>(stop - start)};
    const int K = W.K;
    A.size(b.n, K);
    link_block(W, b, A);
    intron_rows(W, b, A);
    self_solve(W, b, A);
    own_evidence(W, b, A);
    bool has_rows = false;
    int d = 0;
    P::Rows delivered;
    std::fill_n(as_bool(A.held), b.n, false);
    if (W.policy == rigel::POLICY_TRANSFER) {
        std::tie(has_rows, d) = message_layer(W, b, A);
        if (d)
            delivered = P::unpack_rows(b.n, K, d, A.cube.slot.data(), A.cube.pos.data(), as_bool(A.cube.has_pos),
                                       A.cube.neg.data(), as_bool(A.cube.has_neg), W.lam, A.cube.total.data(),
                                       A.cube.opportunity.data(), A.cube.rho.data());
        if (W.want_capture) {
            res.captured = true;
            res.cube = A.cube.prefix(d, K);
            for (int s = 0; s < 2; ++s) res.received[s] = A.received[s].prefix(b.n_owned, K);
        }
    }
    for (int i = 0; i < b.n; ++i) A.row_ptr[i] = has_rows ? A.rows.data() + static_cast<size_t>(i) * K : nullptr;
    count_assertions(b, A, K, has_rows, d, delivered, res);
    if (W.want_capture) capture_slots(b, A, K, has_rows, out);
    final_solve(W, b, A, delivered, out);
}

// ---- the sweep's result ----

nb::dict lane_capture(ReceivedStore::Lane&& ln, int n_owned, int K) {
    nb::dict d;
    d["present"] = bool_ndarray(std::move(ln.present));
    d["profile"] = vec_to_ndarray(std::move(ln.profile), {static_cast<size_t>(n_owned), static_cast<size_t>(K)});
    d["count"] = vec_to_ndarray(std::move(ln.count));
    d["opportunity"] = vec_to_ndarray(std::move(ln.opportunity));
    d["has_witness"] = bool_ndarray(std::move(ln.has_witness));
    d["rna_count"] = vec_to_ndarray(std::move(ln.rna_count));
    d["rna_count_var"] = vec_to_ndarray(std::move(ln.rna_count_var));
    return d;
}

// A block's two received tables as dicts, on its owned slots.
nb::list received_capture(ReceivedStore (&received)[2], int n_owned, int K) {
    nb::list sides;
    for (ReceivedStore& rc : received) {
        nb::dict side;
        side["has_neighbour"] = bool_ndarray(std::move(rc.has_neighbour));
        side["has_composition"] = bool_ndarray(std::move(rc.has_composition));
        side["composition"] = vec_to_ndarray(std::move(rc.composition), {static_cast<size_t>(n_owned), static_cast<size_t>(K)});
        side["level_gdna"] = lane_capture(std::move(rc.lanes[0]), n_owned, K);
        side["level_rna_pos"] = lane_capture(std::move(rc.lanes[1]), n_owned, K);
        side["level_rna_neg"] = lane_capture(std::move(rc.lanes[2]), n_owned, K);
        sides.append(side);
    }
    return sides;
}

// A block's delivered cube rows as (slot, pos, has_pos, neg, has_neg, total, opportunity, rho), or None.
nb::object cube_capture(CubeStore&& c, int K) {
    const size_t d = c.slot.size(), k = static_cast<size_t>(K);
    if (!d) return nb::none();
    return nb::make_tuple(vec_to_ndarray(std::move(c.slot)), vec_to_ndarray(std::move(c.pos), {d, k}),
                          bool_ndarray(std::move(c.has_pos)), vec_to_ndarray(std::move(c.neg), {d, k}),
                          bool_ndarray(std::move(c.has_neg)), vec_to_ndarray(std::move(c.total)),
                          vec_to_ndarray(std::move(c.opportunity)), vec_to_ndarray(std::move(c.rho)));
}

// The counts per block and assertion as (violations, eligible), the delivery bits, and under a capture the cubes and
// received tables.
nb::dict sweep_result(std::vector<BlockResult>& results, const int64_t* blocks, int K) {
    const size_t B = results.size();
    std::vector<int64_t> counts(B * N_ASSERT * 2);
    std::vector<uint8_t> rows_del(B), cube_del(B);
    nb::list cubes, received;
    for (size_t b = 0; b < B; ++b) {
        BlockResult& r = results[b];
        for (int a = 0; a < N_ASSERT; ++a) {
            counts[(b * N_ASSERT + a) * 2] = r.counts[a][0];
            counts[(b * N_ASSERT + a) * 2 + 1] = r.counts[a][1];
        }
        rows_del[b] = r.rows_delivered;
        cube_del[b] = r.cube_delivered;
        if (!r.captured) {
            cubes.append(nb::none());
            received.append(nb::none());
            continue;
        }
        cubes.append(cube_capture(std::move(r.cube), K));
        received.append(received_capture(r.received, static_cast<int>(blocks[3 * b + 1] - blocks[3 * b]), K));
    }
    nb::dict result;
    result["counts"] = vec_to_ndarray(std::move(counts), {B, static_cast<size_t>(N_ASSERT), 2});
    result["rows_delivered"] = bool_ndarray(std::move(rows_del));
    result["cube_delivered"] = bool_ndarray(std::move(cube_del));
    result["cubes"] = cubes;
    result["received"] = received;
    nb::list names;
    for (int a = 0; a < N_ASSERT; ++a) names.append(nb::str(ASSERTION_NAMES[a]));
    result["assertions"] = names;
    return result;
}

// The sweep. `blocks` is (B, 3) int64: start, stop, end per block. `gdna` is None or (log_rho, logP). The out arrays
// arrive as copies of the incoming belief and are written in place on the owned slots; `diagnostics` is None or a dict
// of capture arrays to fill. Returns a dict: `counts` (B, 2, 2) int64, per assertion (violations, eligible);
// `assertions`, their names; `rows_delivered` and `cube_delivered`, (B,) bool; and under a capture `cubes`, per block
// None or the delivered cube rows as (slot, pos, has_pos, neg, has_neg, total, opportunity, rho) in the block's local
// slots, and `received`, per block None or the two received tables as dicts.
nb::dict solve_blocks(BoolVec is_boundary, BoolVec is_exon, BoolVec free_pos, BoolVec free_neg, BoolVec exon_pos,
                      BoolVec exon_neg, BoolVec terminal, IdxVec left, IdxVec right, FlagVec flags, Mat cnt,
                      Mat spliced, Mat sj_count, Mat sj_count_lo, Mat sj_count_hi, Mat route_rate_lo,
                      Mat route_rate_hi, Vec eff_gdna, Vec eff_rna, Vec belief_fpos, Vec belief_fneg, Vec belief_fg,
                      nb::ndarray<int64_t, nb::ndim<2>, nb::c_contig> blocks, double kappa, double disc, Vec lam,
                      int n_tilt, int policy, double rho_gdna, double rho_rna, bool split_live, nb::object gdna,
                      nb::object background, Vec out_fpos, Vec out_fneg, Vec out_fg, Vec out_var, BoolVec out_has_composition,
                      nb::object diagnostics, int n_threads) {
    const int64_t n = static_cast<int64_t>(free_pos.shape(0));
    const int K = static_cast<int>(lam.shape(0)), B = static_cast<int>(blocks.shape(0));
    if (n_tilt < 2) throw std::invalid_argument("solve_blocks: the tilt needs at least two nodes");
    if (policy != rigel::POLICY_SILENT && policy != rigel::POLICY_TRANSFER)
        throw std::invalid_argument("solve_blocks: unknown policy code " + std::to_string(policy));
    if (static_cast<int>(cnt.shape(1)) != 2 || static_cast<int64_t>(cnt.shape(0)) != n)
        throw std::invalid_argument("solve_blocks: the counts are (n, 2)");
    P::Arm curve;  // the fitted gDNA arm's curve, or none
    Vec log_rho, logP;
    if (!gdna.is_none()) {
        nb::tuple t = nb::cast<nb::tuple>(gdna);
        log_rho = exact<Vec>(t[0]);
        logP = exact<Vec>(t[1]);
        curve.log_rho = log_rho.data(); curve.logP = logP.data();
        curve.G = static_cast<int>(log_rho.shape(0));
        if (curve.G < 1 || static_cast<int>(logP.shape(0)) != curve.G)
            throw std::invalid_argument("solve_blocks: the gDNA arm's curve is two (G,) arrays");
    }
    Background bg{};  // the intergenic background's law, or none
    const Background* bg_ptr = nullptr;
    if (!background.is_none()) {
        nb::tuple t = nb::cast<nb::tuple>(background);
        if (nb::len(t) != 3)
            throw std::invalid_argument("solve_blocks: the background is None or (log_rate, size, departure)");
        bg = Background{nb::cast<double>(t[0]), nb::cast<double>(t[1]), nb::cast<double>(t[2])};
        if (!(bg.departure > 0.0 && bg.departure < 1.0))
            throw std::invalid_argument("solve_blocks: the background's departure prior is in (0, 1)");
        bg_ptr = &bg;
    }
    Outputs out{out_fpos.data(), out_fneg.data(), out_fg.data(), out_var.data(), out_has_composition.data()};
    const bool want_capture = !diagnostics.is_none();
    if (want_capture) {
        nb::dict cd = nb::cast<nb::dict>(diagnostics);
        out.fg_loc = exact<Vec>(cd["fg_loc"]).data();
        out.tau_lam = exact<Vec>(cd["tau_lam"]).data();
        out.lam_rows = cd["lam_rows"].is_none() ? nullptr : exact<Mat>(cd["lam_rows"]).data();
    }
    const P::Grid grid(lam.data(), K, kappa, n_tilt, &curve);
    const Sweep W{ChainArrays{is_boundary.data(), is_exon.data(), free_pos.data(), free_neg.data(), exon_pos.data(),
                              exon_neg.data(), terminal.data(), flags.data(), cnt.data(), spliced.data(),
                              sj_count.data(), sj_count_lo.data(), sj_count_hi.data(), route_rate_lo.data(),
                              route_rate_hi.data(), eff_gdna.data(), eff_rna.data()},
                  left.data(), right.data(), Belief{belief_fpos.data(), belief_fneg.data(), belief_fg.data()},
                  kappa, disc, lam.data(), K, policy, rho_gdna, rho_rna, split_live, &grid, want_capture, bg_ptr};
    const int64_t* bl = blocks.data();
    std::vector<BlockResult> results(B);
    const int T = pool_threads(n_threads, B);
    std::vector<Arena> arenas(T);
    run_tasks(T, B, [&](int64_t q, int tid) {
        solve_one_block(W, bl[3 * q], bl[3 * q + 1], bl[3 * q + 2], out, arenas[tid], results[q]);
    });
    return sweep_result(results, bl, K);
}

}  // namespace

void bind_cut_shares(nb::module_& m);  // cut_shares.cpp

NB_MODULE(_solve_impl, m) {
    m.doc() = "The calibration solve: the block pipeline in one call per sweep, and ψ over a slot list.";
    // the policy codes `solve_blocks` switches on; the silent policy runs no message layer
    m.attr("POLICY_SILENT") = rigel::POLICY_SILENT;
    m.attr("POLICY_TRANSFER") = rigel::POLICY_TRANSFER;
    bind_cut_shares(m);
    m.def("solve_blocks", &solve_blocks, nb::arg("is_boundary"), nb::arg("is_exon"), nb::arg("free_pos"), nb::arg("free_neg"),
          nb::arg("exon_pos"), nb::arg("exon_neg"), nb::arg("terminal"), nb::arg("left"), nb::arg("right"), nb::arg("flags"),
          nb::arg("cnt"), nb::arg("spliced"), nb::arg("sj_count"), nb::arg("sj_count_lo"), nb::arg("sj_count_hi"),
          nb::arg("route_rate_lo"), nb::arg("route_rate_hi"), nb::arg("eff_gdna"), nb::arg("eff_rna"), nb::arg("belief_fpos"),
          nb::arg("belief_fneg"), nb::arg("belief_fg"), nb::arg("blocks"), nb::arg("kappa"), nb::arg("disc"), nb::arg("lam"),
          nb::arg("n_tilt"), nb::arg("policy"), nb::arg("rho_gdna"), nb::arg("rho_rna"), nb::arg("split_live"),
          nb::arg("gdna").none(), nb::arg("background").none(), nb::arg("out_fpos").noconvert(),
          nb::arg("out_fneg").noconvert(), nb::arg("out_fg").noconvert(), nb::arg("out_var").noconvert(),
          nb::arg("out_has_composition").noconvert(), nb::arg("diagnostics").none(), nb::arg("n_threads"),
          "Solve every locus block of the chain — the introns' abundance rows, the self-solve, the own evidence, the "
          "message layer, the "
          "final solve, the write-back, the assertions' counts — on a pool of threads, one block at a time; the belief and "
          "`has_composition` written in place on the owned slots; the counts and the capture returned. "
          "Bit-identical at every thread count.");
    m.def("psi_solve", &psi_solve, nb::arg("slots"), nb::arg("u_pos"), nb::arg("u_neg"), nb::arg("allow_pos"),
          nb::arg("allow_neg"), nb::arg("fg_ref"), nb::arg("fpos_ref"), nb::arg("fneg_ref"), nb::arg("kappa"),
          nb::arg("lam"), nb::arg("gdna").none(), nb::arg("row_logprior").none(), nb::arg("cube_slot"),
          nb::arg("cube_pos"), nb::arg("cube_has_pos"), nb::arg("cube_neg"), nb::arg("cube_has_neg"), nb::arg("cube_u"),
          nb::arg("cube_total"), nb::arg("cube_opportunity"), nb::arg("cube_rho"), nb::arg("n_tilt"),
          nb::arg("out_fg").noconvert(), nb::arg("out_fpos").noconvert(), nb::arg("out_fneg").noconvert(),
          nb::arg("out_var").noconvert(), nb::arg("n_threads"),
          "Solve every slot in `slots` on its own cube and write f_g, f_pos, f_neg and Var(log f_g) in place, on "
          "`n_threads` threads (0: every core) — bit-identical at every thread count. `gdna` is the fitted arm as "
          "(log_rho, logP, mass, eff) or None; `row_logprior` the delivered composition rows or None.");

    nb::module_ rows = m.def_submodule("rows", "read_column and ARM_EPS, the kernel's values its Python readers share.");
    rows.def("read_column", [](int64_t col, double kappa) { return read_column(static_cast<int>(col), kappa); },
             nb::arg("col"), nb::arg("kappa"),
             "The genome-strand column strand `col`'s RNA reads on: its own when the library reads sense, else the other.");
    rows.attr("ARM_EPS") = P::ARM_EPS;

    rigel::bind_testing(m.def_submodule("testing", "ψ's pieces, the transfer's kernels and its row constructors, bound for "
                                            "the gates."));
}
