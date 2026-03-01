/**
 * cgranges_bind.cpp — nanobind wrapper for the vendored cgranges C library.
 *
 * Exposes a Python class ``cgranges`` with the same API surface used by
 * hulkrna's index module:
 *
 *   g = cgranges()
 *   g.add(ctg, start, end, label)
 *   g.index()
 *   for start, end, label in g.overlap(ctg, start, end):
 *       ...
 *
 * MIT License — see cgranges.h and khash.h for original copyright notices.
 */

#include <cstdlib>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>

extern "C" {
#include "cgranges.h"
}

namespace nb = nanobind;

class CGRanges {
    cgranges_t *cr_;
    int64_t *buf_;      // reusable output buffer for overlap/contain queries
    int64_t buf_cap_;

public:
    CGRanges() : cr_(cr_init()), buf_(nullptr), buf_cap_(0) {}

    ~CGRanges() {
        if (cr_) cr_destroy(cr_);
        free(buf_);
    }

    // Disable copy (the C struct owns heap memory)
    CGRanges(const CGRanges &) = delete;
    CGRanges &operator=(const CGRanges &) = delete;

    void add(const char *ctg, int32_t start, int32_t end, int32_t label) {
        cr_add(cr_, ctg, start, end, label);
    }

    void index() {
        cr_index(cr_);
    }

    /// Return list of (start, end, label) tuples for intervals overlapping
    /// [start, end).
    nb::list overlap(const char *ctg, int32_t start, int32_t end) {
        int64_t n = cr_overlap(cr_, ctg, start, end, &buf_, &buf_cap_);
        nb::list result;
        for (int64_t i = 0; i < n; i++) {
            int64_t idx = buf_[i];
            result.append(nb::make_tuple(
                cr_start(cr_, idx),
                cr_end(cr_, idx),
                cr_label(cr_, idx)
            ));
        }
        return result;
    }

    /// Return list of (start, end, label) tuples for intervals contained
    /// within [start, end).
    nb::list contain(const char *ctg, int32_t start, int32_t end) {
        int64_t n = cr_contain(cr_, ctg, start, end, &buf_, &buf_cap_);
        nb::list result;
        for (int64_t i = 0; i < n; i++) {
            int64_t idx = buf_[i];
            result.append(nb::make_tuple(
                cr_start(cr_, idx),
                cr_end(cr_, idx),
                cr_label(cr_, idx)
            ));
        }
        return result;
    }
};

NB_MODULE(_cgranges_impl, m) {
    m.doc() = "Vendored cgranges interval-overlap library (nanobind binding)";

    nb::class_<CGRanges>(m, "cgranges")
        .def(nb::init<>(), "Create a new empty interval collection.")
        .def("add", &CGRanges::add,
             nb::arg("ctg"), nb::arg("start"), nb::arg("end"), nb::arg("label"),
             "Add interval [start, end) on contig *ctg* with integer *label*.")
        .def("index", &CGRanges::index,
             "Sort and index all added intervals. Must be called before queries.")
        .def("overlap", &CGRanges::overlap,
             nb::arg("ctg"), nb::arg("start"), nb::arg("end"),
             "Return list of (start, end, label) for intervals overlapping [start, end).")
        .def("contain", &CGRanges::contain,
             nb::arg("ctg"), nb::arg("start"), nb::arg("end"),
             "Return list of (start, end, label) for intervals contained in [start, end).");
}
