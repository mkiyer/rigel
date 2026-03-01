# Set the minimum CMake version and policies for highest tested version
cmake_minimum_required(VERSION 3.15...3.30)

# Set up the project and ensure there is a working C++ compiler
project(hulkrna_ext LANGUAGES C CXX)

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
  src/hulkrna/cgranges/cgranges_bind.cpp
  src/hulkrna/cgranges/cgranges.c
)

# The C wrapper needs to find cgranges.h and khash.h
target_include_directories(_cgranges_impl PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/hulkrna/cgranges
)

# Install into the hulkrna package directory so it's importable as
# hulkrna._cgranges_impl
install(TARGETS _cgranges_impl LIBRARY DESTINATION hulkrna)

# --- Build the _scoring_impl extension ----------------------------------------

# C++ scoring hot-path kernels (compute_fragment_weight, etc.)
nanobind_add_module(
  _scoring_impl
  NOMINSIZE
  STABLE_ABI
  src/hulkrna/native/scoring.cpp
)

# Install into the hulkrna package directory so it's importable as
# hulkrna._scoring_impl
install(TARGETS _scoring_impl LIBRARY DESTINATION hulkrna)

# --- Build the _resolve_impl extension ----------------------------------------

# C++ fragment resolution kernel (resolve_fragment + cgranges queries)
nanobind_add_module(
  _resolve_impl
  NOMINSIZE
  STABLE_ABI
  src/hulkrna/native/resolve.cpp
  src/hulkrna/cgranges/cgranges.c
)

# Needs cgranges.h and khash.h for direct C cgranges calls
target_include_directories(_resolve_impl PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/hulkrna/cgranges
)

# Install into the hulkrna package directory so it's importable as
# hulkrna._resolve_impl
install(TARGETS _resolve_impl LIBRARY DESTINATION hulkrna)
