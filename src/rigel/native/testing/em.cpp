// testing/em.cpp — `_em_impl.testing`: the EM's kernels bound for the tests, which production never calls.
//
//   apply_grouped_prior_update   one grouped prior update (`em_prior.h`), for its conservation identity
//   fast_exp                     the E-step's exponential (`fast_exp.h`) over an array, its NEON lanes and scalar tail
//   EM_LOG_EPSILON               the EM's log floor

#include <cstddef>
#include <stdexcept>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "constants.h"
#include "em_prior.h"
#include "fast_exp.h"
#include "ndarray_util.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;

namespace {

using f64_1d = rigel::Array1<double>;

nb::object fast_exp(f64_1d inputs) {
    std::vector<double> out(inputs.data(), inputs.data() + inputs.shape(0));
    rigel::fast_exp_inplace(out.data(), out.size());
    return rigel::vec_to_ndarray(std::move(out));
}

}  // namespace

void rigel::bind_testing(nb::module_ m) {
    // One update on raw_counts; an empty carried_state or rna_prior_weight is nullptr, as on the EM's lane. The
    // identity holds per call, not end to end: a different rna_prior changes theta, hence the next E-step's raw counts
    // and the converged gDNA total.
    m.def("apply_grouped_prior_update",
          [](f64_1d raw_counts, f64_1d carried_state, f64_1d rna_prior_weight, double gdna_prior_fragments,
             double rna_prior_fragments, int gdna_index, bool has_gdna_candidate) {
              const size_t n = raw_counts.shape(0);
              if (carried_state.size() != 0 && carried_state.shape(0) != n)
                  throw std::runtime_error(
                      "apply_grouped_prior_update: carried_state must be empty or the same length as raw_counts");
              if (rna_prior_weight.size() != 0 && rna_prior_weight.shape(0) != n)
                  throw std::runtime_error(
                      "apply_grouped_prior_update: rna_prior_weight must be empty or the same length as raw_counts");
              const rigel::AggregatePrior prior{gdna_prior_fragments, rna_prior_fragments, gdna_index,
                                                has_gdna_candidate,
                                                rna_prior_weight.size() == 0 ? nullptr : rna_prior_weight.data()};
              std::vector<double> out(n);
              rigel::apply_grouped_prior_update(raw_counts.data(),
                                                carried_state.size() == 0 ? nullptr : carried_state.data(), prior,
                                                out.data(), static_cast<int>(n));
              return rigel::vec_to_ndarray(std::move(out));
          },
          nb::arg("raw_counts"), nb::arg("carried_state"), nb::arg("rna_prior_weight"),
          nb::arg("gdna_prior_fragments"), nb::arg("rna_prior_fragments"), nb::arg("gdna_index"),
          nb::arg("has_gdna_candidate"),
          "One grouped prior update; returns out_counts. Empty carried_state / rna_prior_weight mean nullptr.");
    m.def("fast_exp", &fast_exp, nb::arg("inputs"), "fast_exp over an array: the NEON lanes, then the scalar tail.");
    m.attr("EM_LOG_EPSILON") = rigel::EM_LOG_EPSILON;
}
