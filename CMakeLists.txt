cmake_minimum_required(VERSION 3.20)

# Version detection via gitversion.py
find_package(Python3 COMPONENTS Interpreter REQUIRED)
execute_process(
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/gitversion.py
  OUTPUT_VARIABLE EON_VERSION_RAW
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT EON_VERSION_RAW OR EON_VERSION_RAW STREQUAL "")
  set(EON_VERSION_RAW "0.0.0")
endif()
# Extract just the base version (first 3 dotted numbers)
string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+" EON_VERSION "${EON_VERSION_RAW}")
if(NOT EON_VERSION)
  set(EON_VERSION "0.0.0")
endif()

project(
  eon
  VERSION ${EON_VERSION}
  DESCRIPTION "eOn saddle point finding framework"
  LANGUAGES C CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------- Options (matching meson_options.txt)
option(WITH_WATER "Build Water potentials" OFF)
option(WITH_AMS "Build AMS potentials" OFF)
# LAMMPS, ARTn, IRA, XTB, and VASP have no build-time toggle: they
# are unconditionally compiled and dlopen / fork their runtime
# library / binary on first use. Mirrors meson_options.txt.
option(WITH_MPI "Build with MPI support (via MPItrampoline)" OFF)
option(WITH_GPRD "Build with GPR dimer" OFF)
option(WITH_FORTRAN "Build Fortran potentials" ON)
option(WITH_CUH2 "Build CuH2 potential" ON)
option(WITH_TESTS "Build tests" ON)
option(WITH_GP_SURROGATE "Build GP surrogate" OFF)
option(WITH_CATLEARN "Build CatLearn potential" OFF)
option(WITH_ASE_ORCA "Build ASE ORCA potential" OFF)
option(WITH_ASE_NWCHEM "Build ASE NWChem potential" OFF)
option(WITH_ASE "Build ASE potential" OFF)
option(WITH_QSC "Build QSC potential" OFF)
option(USE_MKL "Enable Intel MKL support" OFF)
option(WITH_METATOMIC "Build Metatomic potential" OFF)
option(WITH_PYTHON "Build Python-embedding potentials" OFF)
set(TORCH_PATH "" CACHE STRING "Path to Torch installation")
set(TORCH_VERSION "2.9" CACHE STRING "Torch version")
option(PIP_METATOMIC "Use pip versions of torch/metatensor/metatomic" OFF)

# ---------------------- Fortran (optional)
if(WITH_FORTRAN OR WITH_CUH2)
  enable_language(Fortran)
endif()

# ---------------------- Python embedding
set(EON_NEED_PYTHON_EMBED FALSE)
if(WITH_ASE OR WITH_ASE_ORCA OR WITH_ASE_NWCHEM OR WITH_CATLEARN)
  set(EON_NEED_PYTHON_EMBED TRUE)
endif()
if(WITH_QSC)
  set(EON_NEED_PYTHON_EMBED TRUE)
endif()
if(WITH_MPI)
  set(EON_NEED_PYTHON_EMBED TRUE)
endif()

if(EON_NEED_PYTHON_EMBED)
  find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
else()
  find_package(Python3 COMPONENTS Interpreter)
endif()

# ---------------------- MPI (via MPItrampoline)
#
# MPI integration always goes through MPItrampoline -- the
# FlexiBLAS-of-MPI. We link against its stable wrapper ABI and one
# eonclient binary forwards every MPI_* call to any spec-compliant
# libmpi at run time, picked via MPITRAMPOLINE_LIB. Mirrors the
# meson MPI block in client/meson.build.
#
# The fall-through here matches meson's wrap-driven behaviour:
#   1. find_package(MPItrampoline) -- system or
#      $CMAKE_PREFIX_PATH-resolvable install.
#   2. FetchContent from upstream v5.5.1 -- analogue of the
#      subprojects/mpitrampoline.wrap. There is never a
#      "no MPI library" failure mode; we deliberately do NOT fall
#      back to FindMPI / a direct system MPI link, since that
#      ABI-locks the binary to one flavour and defeats the point.
if(WITH_MPI)
  find_package(MPItrampoline QUIET)
  if(NOT MPItrampoline_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      MPItrampoline
      GIT_REPOSITORY https://github.com/eschnett/MPItrampoline.git
      GIT_TAG v5.5.1
      GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(MPItrampoline)
  endif()
endif()

# ---------------------- Subdirectories
add_subdirectory(client)
