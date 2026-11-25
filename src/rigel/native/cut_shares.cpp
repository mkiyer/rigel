// cut_shares.cpp — the per-cut arithmetic of a conserved cut share (`effective_length._side_share`): each window of
// lengths is a pair of rows of the pmf's cumulative tables, and each sum over it a difference of two entries.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "ndarray_util.h"
#include "thread_pool.h"

namespace nb = nanobind;

namespace {

using Vec = rigel::Array1<double>;

// Below this many cuts per thread a thread costs more to start than it saves.
constexpr size_t kMinCutsPerThread = 16384;

// The integer lengths lo < w <= hi as rows of a cumulative table of n + 1 entries, and whether the window holds one.
struct Window {
    int64_t i0, i1;
    bool live;
};
inline int64_t table_row(double x, int64_t n) {
    const int64_t i = static_cast<int64_t>(std::floor(x)) + 1;
    return i < 0 ? 0 : (i > n ? n : i);
}
inline Window window(double lo, double hi, int64_t n) {
    const int64_t i0 = table_row(lo, n), i1 = table_row(hi, n);
    return {i0, i1 > i0 ? i1 : i0, i1 > i0};
}
// the table's summand summed over the window's lengths
inline double over(const double* cum, const Window& w) { return w.live ? cum[w.i1] - cum[w.i0] : 0.0; }

// One side's share at one cut (`_side_share`'s docstring derives it): what the cut keeps of the placements on this side,
// less those the far side's end rules out.
inline double side_share(const double* cum0, const double* cum1, const double* cumr, int64_t n, double a, double near,
                         double far) {
    const double top = near + far;
    const double starts_inside = std::min(a + 1.0, top);
    const double runs_through = std::min(near + 1.0, top);
    const Window inside = window(-1.0, starts_inside, n);
    const double kept = 0.5 * (over(cum1, inside) - over(cum0, inside)) +
                        0.5 * a * over(cum0, window(starts_inside, runs_through, n)) +
                        0.5 * a * (near + 1.0) * over(cumr, window(runs_through, top, n));
    const double hi = std::min(far + 1.0 + a, top);
    const Window near_far = window(far + 1.0, hi, n), beyond = window(hi, top, n);
    const double removed = 0.5 * (over(cum1, near_far) - (2.0 * far + 1.0) * over(cum0, near_far) +
                                  far * (far + 1.0) * over(cumr, near_far)) +
                           0.5 * a * (over(cum0, beyond) - far * over(cumr, beyond));
    return std::max(kept - removed, 0.0);
}

nb::object side_shares(Vec cum0, Vec cum1, Vec cumr, Vec a, Vec near, Vec far, int n_threads) {
    const size_t rows = cum0.shape(0);
    if (rows < 1 || cum1.shape(0) != rows || cumr.shape(0) != rows)
        throw std::invalid_argument("side_shares: three cumulative tables of one length, at least one entry");
    const size_t m = a.shape(0);
    if (near.shape(0) != m || far.shape(0) != m)
        throw std::invalid_argument("side_shares: a, near and far are aligned per cut");
    const int64_t n = static_cast<int64_t>(rows) - 1;
    std::vector<double> out(m);
    {
        nb::gil_scoped_release release;
        const double *c0 = cum0.data(), *c1 = cum1.data(), *cr = cumr.data();
        const double *pa = a.data(), *pn = near.data(), *pf = far.data();
        // every cut is its own, so the cuts split into one contiguous range per thread
        const size_t threads = std::min<size_t>(rigel::thread_budget(n_threads), std::max<size_t>(m / kMinCutsPerThread, 1));
        const auto range = [&](size_t t) {
            for (size_t i = m * t / threads, end = m * (t + 1) / threads; i < end; ++i)
                out[i] = side_share(c0, c1, cr, n, pa[i], pn[i], pf[i]);
        };
        // on the pool, so a thread's exception reaches the caller rather than ending the process
        rigel::ThreadPool(static_cast<int>(threads)).run_parallel([&](int t) { range(static_cast<size_t>(t)); });
    }
    return rigel::vec_to_ndarray(std::move(out));
}

}  // namespace

void bind_cut_shares(nb::module_& m) {
    m.def("side_shares", &side_shares, nb::arg("cum0"), nb::arg("cum1"), nb::arg("cumr"), nb::arg("a"),
          nb::arg("near"), nb::arg("far"), nb::arg("n_threads") = 1,
          "One side's conserved share at every cut (`effective_length._side_share`), from the pmf's three cumulative "
          "tables (f, w·f and f/w, each with a leading 0), on n_threads threads (0: every core).");
}
