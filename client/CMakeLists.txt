cmake_minimum_required(VERSION 3.10)
set(CMAKE_MODULE_PATH /usr/share/cmake-3.19/Modules/)

# set the project name
project(
  eonclient
  VERSION 1.0
  DESCRIPTION "EON client"
  LANGUAGES C CXX Fortran)

set(CMAKE_COLOR_MAKEFILE ON)

add_compile_options(# "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
                    "$<$<CONFIG:Debug>:-O0;-g3;-ggdb>")
# Messages
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message(STATUS "Fortran compiler: " ${CMAKE_C_COMPILER})

# Compiler gfortran Flags
set(CMAKE_Fortran_FLAGS "-fbounds-check ") # -qopenmp
set(CMAKE_Fortran_FLAGS_RELEASE
    "-O3 -march=native -ffast-math -funroll-loops -std=legacy -fopenmp")
set(CMAKE_Fortran_FLAGS_DEBUG
    "-O0 -Warray-bounds  -std=legacy -Wall -Wno-unused ")
set(DOUBLE "-r8")

# Add the new cpp files here
add_executable(eonclient ClientEON.cpp CommandLine.cpp INIFile.cpp)

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
  message(STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_RELEASE})
elseif(${CMAKE_BUILD_TYPE} MATCHES "Debug")
  message(STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_DEBUG})
endif()

add_custom_target(
  version
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/version.sh >
          ${CMAKE_CURRENT_SOURCE_DIR}/version.h
  COMMENT "run version update command in ${CMAKE_CURRENT_SOURCE_DIR}")
add_dependencies(eonclient version)

target_link_libraries(
  eonclient
  Aluminum
  EAM
  EDIP
  EMT
  ExtPot
  FeHe
  GPRPotential
  # IMD LAMMPS
  Lenosky
  LJ
  LJCluster
  Morse
  # MPIPot
  NewPot
  # PyAMFF
  QSC
  SW
  Tersoff
  VASP
  Water
  Water_H
  Water_Pt
  eonlib)

set(basePotentials
    Aluminum
    EAM
    EDIP
    EMT
    ExtPot
    FeHe
    GPRPotential
    Lenosky
    LJ
    LJCluster
    Morse
    NewPot
    QSC
    SW
    Tersoff
    VASP
    Water
    Water_H
    Water_Pt)

foreach(pot ${basePotentials})
  add_subdirectory(potentials/${pot})
endforeach()

# add_library(potentials SHARED ${basePotentials}) set_property(TARGET
# potentials PROPERTY LINKER_LANGUAGE Fortran)

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
  TestJob.cpp
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
  AtomicGPDimer.cpp
  gprdimer/gpr/AtomicDimer.cpp
  gprdimer/backend/DistributionFunctions.cpp
  gprdimer/gpr/auxiliary/ProblemSetUp.cpp
  gprdimer/gpr/auxiliary/Gradient.cpp
  gprdimer/gpr/auxiliary/Distance.cpp
  gprdimer/gpr/prior/PriorGaussian.cpp
  gprdimer/gpr/prior/PriorLogUnif.cpp
  gprdimer/gpr/prior/PriorSqrtt.cpp
  gprdimer/gpr/prior/PriorT.cpp
  gprdimer/gpr/dimer/LBFGS.cpp
  gprdimer/gpr/covariance_functions/ConstantCF.cpp
  gprdimer/gpr/covariance_functions/SexpatCF.cpp
  gprdimer/gpr/observation_models/LikGaussian.cpp
  gprdimer/gpr/ml/GaussianProcessRegression.cpp
  gprdimer/gpr/ml/SCG.cpp)
