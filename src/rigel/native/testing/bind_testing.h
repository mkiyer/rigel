// bind_testing.h — an extension's test seams, bound as its `testing` submodule by the last line of its NB_MODULE.
// Each extension defines its own in testing/: bam.cpp, resolve.cpp, scoring.cpp, em.cpp and solve.cpp.
#pragma once

#include <nanobind/nanobind.h>

namespace rigel {

void bind_testing(nanobind::module_ m);

}  // namespace rigel
