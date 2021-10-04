cmake_minimum_required(VERSION 3.11...3.19)

project(
  eonclient
  VERSION 1.0
  DESCRIPTION "EON client"
  LANGUAGES C CXX Fortran)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

set(CMAKE_COLOR_MAKEFILE ON)
cmake_policy(SET CMP0079 NEW)

# # CXX Compiler Flags
# string(APPEND CMAKE_CXX_FLAGS
#        " -Wno-unused-local-typedefs -Wno-unused-function -Wfatal-errors ")
# string(APPEND CMAKE_CXX_FLAGS_RELEASE
#        " -O3 -march=native -ffast-math -funroll-loops ")
# # string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -Warray-bounds -Wall -Wno-unused -g3 -ggdb ")

# # Fortran Compiler Flags
# string(APPEND CMAKE_Fortran_FLAGS
#        " -fbounds-check -ffpe-trap=invalid,zero,overflow ")
# string(APPEND CMAKE_Fortran_FLAGS_RELEASE
#        " -O3 -march=native -ffast-math -funroll-loops -std=legacy ")
# string(APPEND CMAKE_Fortran_FLAGS_DEBUG
#        " -O0 -Warray-bounds -std=legacy -Wall -Wno-unused ")

# Other flags
# set(DOUBLE "-r8")

if(NO_WARN)
  string(APPEND CMAKE_CXX_FLAGS " -w ")
  string(APPEND CMAKE_Fortran_FLAGS " -w ")
endif()

# add_compile_options("$<$<CONFIG:Debug>:-O0;-g3;-ggdb;-Warray-bounds;-Wall;-Wno-unused>")

# add_compile_options(# "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
# "$<$<CONFIG:Debug>:-O0;-g3;-ggdb>")
# -g -Wall -Wno-unused-local-typedefs -Wno-unused-function -Wfatal-errors
# -lpython3 -O0  -w -DLINUX

# Messages
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message(STATUS "C compiler: " ${CMAKE_C_COMPILER})
message(STATUS "CXX compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS "Python: " ${Python_EXECUTABLE})

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
  find_package(BLAS)
  find_package(LAPACK)
  find_package(OpenMP)
  message(STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS}
                 ${CMAKE_Fortran_FLAGS_RELEASE} ${COMPILE_OPTIONS})
  message(STATUS "CXX compiler flags:" ${CMAKE_CXX_FLAGS}
                 ${CMAKE_CXX_FLAGS_RELEASE} ${COMPILE_OPTIONS})
elseif(${CMAKE_BUILD_TYPE} MATCHES "Debug")
  message(STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS}
                 ${CMAKE_Fortran_FLAGS_DEBUG} ${COMPILE_OPTIONS})
  message(STATUS "CXX compiler flags:" ${CMAKE_CXX_FLAGS}
                 ${CMAKE_CXX_FLAGS_DEBUG} ${COMPILE_OPTIONS})
endif()

# Build
add_executable(eonclient ClientEON.cpp CommandLine.cpp)

# Optional, for nix
if(FIND_EIGEN)
  find_package(Eigen3 3.3 REQUIRED NO_MODULE)
  include_directories(${EIGEN3_INCLUDE_DIR})
  message(STATUS "EIGEN: " ${EIGEN3_INCLUDE_DIR})
  target_link_libraries(eonclient Eigen3::Eigen)
endif()

# Set version
target_compile_features(eonclient PUBLIC cxx_std_17)
set_target_properties(
  eonclient
  PROPERTIES CXX_STANDARD 17
             CXX_STANDARD_REQUIRED YES
             CXX_EXTENSIONS NO)

# Versioning
add_custom_target(
  version
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/version.sh >
          ${CMAKE_CURRENT_SOURCE_DIR}/version.h
  COMMENT "run version update command in ${CMAKE_CURRENT_SOURCE_DIR}")

# Build order
add_dependencies(eonclient potentials eonlib version gprdimer)

# Libraries
target_link_libraries(eonclient potentials eonlib gprdimer)

# Sub-directories
add_subdirectory(subprojects/gprdimer)
add_subdirectory(potentials)

# Dependencies
add_library(
  eonlib SHARED
  MinModeSaddleSearch.cpp
  Dimer.cpp
  EpiCenters.cpp
  Hessian.cpp
  ConjugateGradients.cpp
  HelperFunctions.cpp
  Matter.cpp
  Parameters.cpp
  INIFile.cpp
  Potential.cpp
  Quickmin.cpp
  ProcessSearchJob.cpp
  PointJob.cpp
  MinimizationJob.cpp
  HessianJob.cpp
  ParallelReplicaJob.cpp
  SafeHyperJob.cpp
  TADJob.cpp
  ReplicaExchangeJob.cpp
  Dynamics.cpp
  BondBoost.cpp
  FiniteDifferenceJob.cpp
  NudgedElasticBandJob.cpp
  # TestJob.cpp
  GPRHelpers.cpp
  BasinHoppingJob.cpp
  SaddleSearchJob.cpp
  ImprovedDimer.cpp
  NudgedElasticBand.cpp
  Lanczos.cpp
  Bundling.cpp
  Job.cpp
  CommandLine.cpp
  DynamicsJob.cpp
  Log.cpp
  LBFGS.cpp
  LowestEigenmode.cpp
  Optimizer.cpp
  Prefactor.cpp
  DynamicsSaddleSearch.cpp
  PrefactorJob.cpp
  FIRE.cpp
  GlobalOptimizationJob.cpp
  GlobalOptimization.cpp
  StructureComparisonJob.cpp
  MonteCarloJob.cpp
  MonteCarlo.cpp
  SteepestDescent.cpp
  BasinHoppingSaddleSearch.cpp
  BiasedGradientSquaredDescent.cpp
  AtomicGPDimer.cpp)

# Options
if(MPI)
  target_compile_definitions(eonclient PUBLIC -DEONMPI)
  # Python is used with MPI
  find_package(Python COMPONENTS Interpreter Development)
  include_directories(${Python_INCLUDE_DIRS})
  # Need the MPI libraries
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  target_link_libraries(eonclient ${Python_LIBRARIES} MPI::MPI_CXX)
endif()

# https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/How-To-Write-Platform-Checks
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_compile_definitions(eonlib PUBLIC -DOSX)
  target_compile_definitions(eonclient PUBLIC -DOSX)
  elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  target_compile_definitions(eonlib PUBLIC -DLINUX)
    target_compile_definitions(eonclient PUBLIC -DLINUX)
  elseif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_compile_definitions(eonlib PUBLIC -DWIN32) # No VASP
    target_compile_definitions(eonclient PUBLIC -DWIN32)
endif()

if(WITH_WATER)
  target_compile_definitions(eonlib PUBLIC -DWITH_WATER)
endif()

if(WITH_AMS)
  target_compile_definitions(eonlib PUBLIC -DWITH_AMS)
endif()

if(WITH_GPRD)
  target_compile_definitions(eonlib PUBLIC -DWITH_GPRD)
endif()

if(WITH_LAMMPS)
  target_compile_definitions(eonlib PUBLIC -DLAMMPS_POT)
  target_compile_definitions(eonclient PUBLIC -DLAMMPS_POT)
endif()

if(WITH_FORTRAN)
  target_compile_definitions(eonlib PUBLIC -DWITH_FORTRAN)
endif()

if(WITH_PYTHON)
  # Python is used for potentials
  find_package(Python COMPONENTS Interpreter Development)
  include_directories(${Python_INCLUDE_DIRS})
  target_link_libraries(eonclient ${Python_LIBRARIES})
endif()

# Tests
option(PACKAGE_TESTS "Build the tests" OFF)
if(PACKAGE_TESTS)
  find_package(GTest REQUIRED)
  enable_testing()
  include(GoogleTest)
  add_subdirectory(gtests)
endif()

# Install
install(TARGETS eonclient)
install(TARGETS eonlib)
