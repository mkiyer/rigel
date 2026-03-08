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

# --- Build the _cgranges_impl extension ---------------------------------------

# Compile the nanobind extension with optimizations for speed (NOMINSIZE)
# since cgranges overlap queries are a hot path.
nanobind_add_module(
  _cgranges_impl
  NOMINSIZE
  STABLE_ABI
  src/rigel/cgranges/cgranges_bind.cpp
  src/rigel/cgranges/cgranges.c
)

# The C wrapper needs to find cgranges.h and khash.h
target_include_directories(_cgranges_impl PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rigel/cgranges
)

# Install into the rigel package directory so it's importable as
# rigel._cgranges_impl
install(TARGETS _cgranges_impl LIBRARY DESTINATION rigel)

# --- Build the _scoring_impl extension ----------------------------------------

# C++ scoring hot-path kernels (compute_fragment_weight, etc.)
nanobind_add_module(
  _scoring_impl
  NOMINSIZE
  STABLE_ABI
  src/rigel/native/scoring.cpp
)

# Install into the rigel package directory so it's importable as
# rigel._scoring_impl
install(TARGETS _scoring_impl LIBRARY DESTINATION rigel)

# --- Build the _resolve_impl extension ----------------------------------------

# C++ fragment resolution kernel (resolve_fragment + cgranges queries)
nanobind_add_module(
  _resolve_impl
  NOMINSIZE
  STABLE_ABI
  src/rigel/native/resolve.cpp
  src/rigel/cgranges/cgranges.c
)

# Needs cgranges.h, khash.h, and shared native headers
target_include_directories(_resolve_impl PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rigel/cgranges
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rigel/native
)

# Install into the rigel package directory so it's importable as
# rigel._resolve_impl
install(TARGETS _resolve_impl LIBRARY DESTINATION rigel)

# --- Build the _bam_impl extension --------------------------------------------

# C++ BAM scanner using htslib — replaces the Python BAM parsing hot path.
# Reads BAM → resolves fragments → buffers results in one C++ call.

# Find htslib (pkg-config or cmake find)
find_package(PkgConfig)
if(PkgConfig_FOUND)
  pkg_check_modules(HTSLIB IMPORTED_TARGET htslib)
endif()

# Fallback: look for htslib in conda prefix or CONDA_PREFIX
if(NOT HTSLIB_FOUND AND DEFINED ENV{CONDA_PREFIX})
  set(_conda "$ENV{CONDA_PREFIX}")
  if(EXISTS "${_conda}/include/htslib/hts.h" AND EXISTS "${_conda}/lib/libhts.dylib")
    set(HTSLIB_FOUND TRUE)
    set(HTSLIB_INCLUDE_DIRS "${_conda}/include")
    set(HTSLIB_LIBRARY_DIRS "${_conda}/lib")
    set(HTSLIB_LIBRARIES "hts")
    message(STATUS "htslib found via CONDA_PREFIX: ${_conda}")
  elseif(EXISTS "${_conda}/include/htslib/hts.h" AND EXISTS "${_conda}/lib/libhts.so")
    set(HTSLIB_FOUND TRUE)
    set(HTSLIB_INCLUDE_DIRS "${_conda}/include")
    set(HTSLIB_LIBRARY_DIRS "${_conda}/lib")
    set(HTSLIB_LIBRARIES "hts")
    message(STATUS "htslib found via CONDA_PREFIX: ${_conda}")
  endif()
endif()

if(HTSLIB_FOUND)
  nanobind_add_module(
    _bam_impl
    NOMINSIZE
    STABLE_ABI
    src/rigel/native/bam_scanner.cpp
    src/rigel/cgranges/cgranges.c
  )

  target_include_directories(_bam_impl PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/rigel/cgranges
    ${CMAKE_CURRENT_SOURCE_DIR}/src/rigel/native
  )

  # Link htslib: prefer PkgConfig target, fall back to manual link
  if(TARGET PkgConfig::HTSLIB)
    target_link_libraries(_bam_impl PRIVATE PkgConfig::HTSLIB)
  else()
    target_include_directories(_bam_impl PRIVATE ${HTSLIB_INCLUDE_DIRS})
    target_link_directories(_bam_impl PRIVATE ${HTSLIB_LIBRARY_DIRS})
    target_link_libraries(_bam_impl PRIVATE ${HTSLIB_LIBRARIES})
  endif()

  install(TARGETS _bam_impl LIBRARY DESTINATION rigel)

  message(STATUS "htslib found — building _bam_impl C++ BAM scanner")
else()
  message(STATUS "htslib NOT found — _bam_impl will not be built. "
                 "Install htslib (conda install -c conda-forge htslib) "
                 "for native BAM scanning.")
endif()

# --- Build the _em_impl extension ---------------------------------------------

# C++ EM solver — replaces the Python EM hot path (_em_step, _vbem_step,
# _build_equiv_classes, SQUAREM loop) with a single native call.
# Pure C++17 + nanobind + numpy, no external dependencies.
nanobind_add_module(
  _em_impl
  NOMINSIZE
  STABLE_ABI
  src/rigel/native/em_solver.cpp
)
# Keep debug symbols for profiling (dsym/DWARF) even in Release builds
target_compile_options(_em_impl PRIVATE -g)

install(TARGETS _em_impl LIBRARY DESTINATION rigel)

# --- SIMD / Architecture Detection -------------------------------------------
#
# Detect CPU architecture and report available SIMD instruction sets.
# The C++ header simd_detect.h uses standard compiler-defined macros
# (__ARM_NEON, __AVX2__, etc.) for compile-time detection, plus CPUID
# for runtime dispatch on x86.  With -march=native the compiler sets
# these macros automatically; for portable builds the C++ code uses
# function-level __attribute__((target(...))) for multi-versioned SIMD.
#
# This block just reports what's available and sets RIGEL_ARCH for
# any downstream CMake logic.

if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
    set(RIGEL_ARCH "ARM64")
    message(STATUS "SIMD: ARM64 — NEON is baseline (always available)")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
    set(RIGEL_ARCH "X86_64")
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-mavx2 -mfma" _rigel_has_avx2_fma)
    check_cxx_compiler_flag("-mavx512f" _rigel_has_avx512f)
    if(_rigel_has_avx512f)
        message(STATUS "SIMD: x86_64 — SSE2 + AVX2/FMA + AVX-512F available")
    elseif(_rigel_has_avx2_fma)
        message(STATUS "SIMD: x86_64 — SSE2 + AVX2/FMA available")
    else()
        message(STATUS "SIMD: x86_64 — SSE2 only")
    endif()
else()
    set(RIGEL_ARCH "UNKNOWN")
    message(STATUS "SIMD: ${CMAKE_SYSTEM_PROCESSOR} — no explicit SIMD detection")
endif()

# --- Compiler optimizations ---------------------------------------------------
#
# Phase 2 performance flags: -O3, -march=native, FMA contraction, and LTO.
# Applied per-target to control fast-math scope.
#
# _em_impl gets -ffp-contract=fast (enables FMA) but NOT full -ffast-math,
# because -ffast-math allows the compiler to reassociate floating-point
# additions, which would break the Kahan compensated summation used in
# the E-step kernel for numerical stability.
#
# _scoring_impl gets full -ffast-math (no compensated summation there).
# All other targets get -O3 -march=native without fast-math.

# RIGEL_PORTABLE: when ON, omit -march=native so wheels run on any machine
# of the same architecture. CI wheel builds set this; local dev leaves it OFF.
option(RIGEL_PORTABLE "Build portable wheels (no -march=native)" OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    if(RIGEL_PORTABLE)
        set(RIGEL_OPT_FLAGS -O3 -DNDEBUG)
    else()
        set(RIGEL_OPT_FLAGS -O3 -march=native -DNDEBUG)
    endif()

    target_compile_options(_em_impl      PRIVATE ${RIGEL_OPT_FLAGS} -ffp-contract=fast)
    target_compile_options(_scoring_impl PRIVATE ${RIGEL_OPT_FLAGS} -ffast-math)
    target_compile_options(_resolve_impl PRIVATE ${RIGEL_OPT_FLAGS})
    target_compile_options(_cgranges_impl PRIVATE ${RIGEL_OPT_FLAGS})
    if(TARGET _bam_impl)
        target_compile_options(_bam_impl PRIVATE ${RIGEL_OPT_FLAGS})
    endif()
endif()

# Link-Time Optimization (LTO): enables cross-TU inlining and dead code
# elimination across compilation units.
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set_property(TARGET _em_impl      PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET _scoring_impl PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET _resolve_impl PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET _cgranges_impl PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    if(TARGET _bam_impl)
        set_property(TARGET _bam_impl PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif()
