// cut_shares.cpp — the per-cut arithmetic of a conserved cut share (`effective_length._side_share`), element by
// element in numpy's own operation order. This file is built with -ffp-contract=off: no product is fused into a sum,
// so every value is the one numpy's element-wise operations give, to the bit. The pmf's cumulative tables are numpy's.
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "ndarray_util.h"

namespace nb = nanobind;

namespace {

using Vec = nb::ndarray<const double, nb::ndim<1>, nb::c_contig>;

// numpy's np.minimum and np.maximum on one pair: a NaN in either propagates
inline double np_minimum(double x, double y) { return (x <= y || std::isnan(x)) ? x : y; }
inline double np_maximum(double x, double y) { return (x >= y || std::isnan(x)) ? x : y; }

// The integer lengths lo < w <= hi as rows of a cumulative table of n + 1 entries: np.clip(np.floor(x) + 1, 0, n)
// at each end, the upper end held at the lower, and whether the window holds a length.
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

// One side's share at one cut: `_side_share`'s expression, each operation in numpy's order.
inline double side_share(const double* cum0, const double* cum1, const double* cumr, int64_t n, double a, double near,
                         double far) {
    const double top = near + far;
    const double starts_inside = np_minimum(a + 1.0, top);
    const double runs_through = np_minimum(near + 1.0, top);
    const Window inside = window(-1.0, starts_inside, n);
    const double half_a = 0.5 * a;
    const double kept_inside = 0.5 * (over(cum1, inside) - over(cum0, inside));
    const double kept_through = half_a * over(cum0, window(starts_inside, runs_through, n));
    const double kept_beyond = (half_a * (near + 1.0)) * over(cumr, window(runs_through, top, n));
    const double kept = (kept_inside + kept_through) + kept_beyond;
    const double lo = far + 1.0;
    const double hi = np_minimum((far + 1.0) + a, top);
    const Window near_far = window(lo, hi, n), beyond = window(hi, top, n);
    const double near_far_sum = (over(cum1, near_far) - ((2.0 * far) + 1.0) * over(cum0, near_far)) +
                                (far * (far + 1.0)) * over(cumr, near_far);
    const double removed = 0.5 * near_far_sum + half_a * (over(cum0, beyond) - far * over(cumr, beyond));
    return np_maximum(kept - removed, 0.0);
}

nb::object side_shares(Vec cum0, Vec cum1, Vec cumr, Vec a, Vec near, Vec far) {
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
        for (size_t i = 0; i < m; ++i) out[i] = side_share(c0, c1, cr, n, pa[i], pn[i], pf[i]);
    }
    return rigel::vec_to_ndarray(std::move(out));
}

}  // namespace

void bind_cut_shares(nb::module_& m) {
    m.def("side_shares", &side_shares, nb::arg("cum0"), nb::arg("cum1"), nb::arg("cumr"), nb::arg("a"),
          nb::arg("near"), nb::arg("far"),
          "One side's conserved share at every cut (`effective_length._side_share`), from the pmf's three cumulative "
          "tables (f, w·f and f/w, each with a leading 0): numpy's element-wise arithmetic, to the bit.");
}
