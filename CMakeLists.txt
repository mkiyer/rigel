# Set the minimum CMake version and policies for highest tested version
cmake_minimum_required(VERSION 3.15...3.30)

# Set up the project and ensure there is a working C++ compiler
project(rigel_ext LANGUAGES C CXX)

# Warn if the user invokes CMake directly
if (NOT SKBUILD)
  message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build-core'. \
  Running it directly will almost certainly not produce the desired \
  result. If you are a user trying to install this package, use: \
  ===================================================================== \
   $ pip install . \
  ===================================================================== \
  For development with incremental rebuilds: \
  ===================================================================== \
   $ pip install nanobind scikit-build-core[pyproject] \
   $ pip install --no-build-isolation -ve . \
  =====================================================================")
endif()

# Don't use compiler extensions (e.g. -std=gnu++17 -> -std=c++17)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release build by default
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# --- Find Python and nanobind ------------------------------------------------

# Try to import all Python components potentially needed by nanobind
if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()

find_package(Python 3.12
  REQUIRED COMPONENTS Interpreter ${DEV_MODULE}
  OPTIONAL_COMPONENTS Development.SABIModule)

# Import nanobind through CMake's find_package mechanism
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# --- Build options ------------------------------------------------------------

# RIGEL_PORTABLE (default ON): omit -march=native, so a build runs on every machine of its architecture. A native
# build, tuned to the build host's CPU, dies with an illegal instruction on an older one; build one only for the
# machine it runs on, with -DRIGEL_PORTABLE=OFF (e.g. CMAKE_ARGS="-DRIGEL_PORTABLE=OFF" pip install ...).
option(RIGEL_PORTABLE "Build for every machine of the architecture (no -march=native)" ON)
# RIGEL_PROFILE_NATIVE (default OFF): symbols and frame pointers, no LTO, for a sampling profiler. pip strips the
# installed modules too unless told not to: -C cmake.define.RIGEL_PROFILE_NATIVE=ON -C install.strip=false.
option(RIGEL_PROFILE_NATIVE "Build with debug symbols and frame pointers, and without LTO, for a profiler" OFF)

# The flags are GCC's and Clang's; another compiler builds with its defaults. A Debug build keeps its own flags: no
# optimisation, assertions on, no link-time optimisation.
set(RIGEL_DEBUG_BUILD OFF)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(RIGEL_DEBUG_BUILD ON)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  set(RIGEL_OPT_FLAGS "")
  if(NOT RIGEL_DEBUG_BUILD)
    list(APPEND RIGEL_OPT_FLAGS -O3 -DNDEBUG)
  endif()
  if(NOT RIGEL_PORTABLE)
    list(APPEND RIGEL_OPT_FLAGS -march=native)
  endif()
  if(RIGEL_PROFILE_NATIVE)
    list(APPEND RIGEL_OPT_FLAGS -g -fno-omit-frame-pointer)
  endif()
  set(RIGEL_FAST_MATH -ffast-math -fno-finite-math-only)
  set(RIGEL_FMA_ONLY -ffp-contract=fast)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT RIGEL_IPO_SUPPORTED)

# One native module, built and installed like every other: optimised, link-time optimised unless profiling, and
# importable as rigel.<name>. The headers in native/ and the vendored cgranges are visible to every module.
function(rigel_extension name)
  # nanobind strips a Release module's symbols unless told not to, which a profiler needs
  set(rigel_symbols "")
  if(RIGEL_PROFILE_NATIVE)
    set(rigel_symbols NOSTRIP)
  endif()
  nanobind_add_module(${name} NOMINSIZE STABLE_ABI ${rigel_symbols} ${ARGN})
  target_include_directories(${name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/rigel/native
    ${CMAKE_CURRENT_SOURCE_DIR}/src/rigel/cgranges)
  target_compile_options(${name} PRIVATE ${RIGEL_OPT_FLAGS})
  if(RIGEL_IPO_SUPPORTED AND NOT RIGEL_PROFILE_NATIVE AND NOT RIGEL_DEBUG_BUILD)
    set_property(TARGET ${name} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
  install(TARGETS ${name} LIBRARY DESTINATION rigel)
endfunction()

# --- htslib, which the BAM scanner links ----------------------------------------

# pkg-config first, without IMPORTED_TARGET: that resolves htslib's own dependencies through their .pc files (bzip2.pc),
# which a conda cross-compilation build may lack. The shared libhts is linked directly, so its dependencies need no
# link flags.
find_package(PkgConfig)
if(PkgConfig_FOUND)
  pkg_check_modules(HTSLIB htslib)
endif()

# Otherwise an htslib installed under a prefix
macro(_find_htslib_in_prefix _prefix _label)
  if(NOT HTSLIB_FOUND AND EXISTS "${_prefix}/include/htslib/hts.h")
    foreach(_lib libhts.dylib libhts.so libhts.a)
      if(EXISTS "${_prefix}/lib/${_lib}")
        set(HTSLIB_FOUND TRUE)
        set(HTSLIB_INCLUDE_DIRS "${_prefix}/include")
        set(HTSLIB_LIBRARY_DIRS "${_prefix}/lib")
        set(HTSLIB_LIBRARIES "hts")
        message(STATUS "htslib found via ${_label}: ${_prefix}")
        break()
      endif()
    endforeach()
  endif()
endmacro()

# conda-build's host environment ($PREFIX), then the active environment ($CONDA_PREFIX)
if(DEFINED ENV{PREFIX})
  _find_htslib_in_prefix("$ENV{PREFIX}" "conda PREFIX")
endif()
if(DEFINED ENV{CONDA_PREFIX})
  _find_htslib_in_prefix("$ENV{CONDA_PREFIX}" "CONDA_PREFIX")
endif()

if(NOT HTSLIB_FOUND)
  message(FATAL_ERROR
    "htslib NOT found. The BAM scanner needs it: install htslib (conda install -c conda-forge htslib) and re-run.")
endif()

# --- The native modules ---------------------------------------------------------

# _bam_impl: the BAM scan, the annotated-BAM pass, and the calibration accumulator the scan fills
rigel_extension(_bam_impl
  src/rigel/native/scan/bam_module.cpp
  src/rigel/native/scan/bam_record.cpp
  src/rigel/native/scan/template_hits.cpp
  src/rigel/native/scan/bam_scanner.cpp
  src/rigel/native/scan/scan_result.cpp
  src/rigel/native/scan/annotation_writer.cpp
  src/rigel/native/calibration/accumulator.cpp
  src/rigel/native/calibration/accumulator_set.cpp
  src/rigel/native/calibration/deferred_fragments.cpp
  src/rigel/native/resolve/fragment_resolver.cpp
  src/rigel/native/testing/bam.cpp
  src/rigel/cgranges/cgranges.c)
target_include_directories(_bam_impl PRIVATE ${HTSLIB_INCLUDE_DIRS})
target_link_directories(_bam_impl PRIVATE ${HTSLIB_LIBRARY_DIRS})
target_link_libraries(_bam_impl PRIVATE ${HTSLIB_LIBRARIES})

# _resolve_impl: a fragment's overlap resolution against the index
rigel_extension(_resolve_impl
  src/rigel/native/resolve/resolve_module.cpp
  src/rigel/native/resolve/fragment_resolver.cpp
  src/rigel/native/testing/resolve.cpp
  src/rigel/cgranges/cgranges.c)

# _scoring_impl: per-fragment likelihoods. Full fast-math, since nothing there is compensated, but with infinities,
# since the scorer's sentinels are -inf.
rigel_extension(_scoring_impl src/rigel/native/scoring.cpp src/rigel/native/testing/scoring.cpp)
target_compile_options(_scoring_impl PRIVATE ${RIGEL_FAST_MATH})

# _em_impl: the locus EM. FMA contraction only: fast-math would let the compiler reassociate the E-step's Kahan
# compensated sums away.
rigel_extension(_em_impl src/rigel/native/em_solver.cpp src/rigel/native/testing/em.cpp)
target_compile_options(_em_impl PRIVATE ${RIGEL_FMA_ONLY})

# _solve_impl: the calibration solve, the sweep's block pipeline in one call per sweep
rigel_extension(_solve_impl src/rigel/native/solve_kernel.cpp src/rigel/native/testing/solve.cpp
                src/rigel/native/cut_shares.cpp)
