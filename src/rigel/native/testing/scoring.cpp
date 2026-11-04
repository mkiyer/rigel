// testing/scoring.cpp — `_scoring_impl.testing`: the scorer's constants, bound for the parity tests that hold them to
// Python's own values. Production never reads them from here.

#include <nanobind/nanobind.h>

#include "constants.h"
#include "testing/bind_testing.h"

namespace nb = nanobind;

void rigel::bind_testing(nb::module_ m) {
    m.attr("LOG_HALF") = rigel::LOG_HALF;
    m.attr("TAIL_DECAY_LP") = rigel::TAIL_DECAY_LP;
}
