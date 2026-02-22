# ---------------------- Compiler setup
set(EON_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty
)

# Fortran compiler flags
if(WITH_FORTRAN OR WITH_CUH2)
  include(CheckFortranCompilerFlag)
  set(EON_FORTRAN_FLAGS "")
  foreach(flag -fno-implicit-none -Wno-compare-reals -Wno-conversion
               -Wno-line-truncation -Wno-surprising -Wno-unused-dummy-argument
               -Wno-deprecated -Wno-unused-label -w -Wno-maybe-uninitialized
               -Wno-uninitialized -Wno-unused-but-set-parameter
               -Wno-unused-but-set-variable -Wno-unused-parameter
               -Wno-unused-variable)
    string(MAKE_C_IDENTIFIER "FORTRAN_SUPPORTS${flag}" flag_var)
    check_fortran_compiler_flag("${flag}" ${flag_var})
    if(${flag_var})
      list(APPEND EON_FORTRAN_FLAGS "${flag}")
    endif()
  endforeach()
endif()

# MKL support
if(USE_MKL)
  find_package(MKL CONFIG REQUIRED)
  add_compile_definitions(EIGEN_USE_MKL_ALL)
  link_libraries(mkl_intel_lp64 mkl_core mkl_sequential)
endif()

# ---------------------- Platform detection
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  add_compile_definitions(OSX=TRUE)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-faligned-allocation>)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  add_compile_definitions(IS_WINDOWS)
endif()

# Compiler warning suppression (matching Meson)
include(CheckCXXCompilerFlag)
foreach(flag -Wno-class-memaccess -Wno-unused-function -Wno-reorder
             -Wno-sign-compare -Wno-unused-but-set-parameter
             -Wno-unused-but-set-variable -Wno-unused-parameter
             -Wno-unused-result -Wno-unused-variable)
  string(MAKE_C_IDENTIFIER "CXX_SUPPORTS${flag}" flag_var)
  check_cxx_compiler_flag("${flag}" ${flag_var})
  if(${flag_var})
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${flag}>)
  endif()
endforeach()

# Safety checks in non-release builds
if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_definitions(EON_CHECKS=TRUE)
endif()

# ---------------------- Version detection
# Python3 already found by top-level CMakeLists.txt
execute_process(
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/get_version.py
  OUTPUT_VARIABLE VERSION_OUTPUT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  ERROR_QUIET)
if(NOT VERSION_OUTPUT)
  set(VERSION_OUTPUT "unknown,unknown,unknown")
endif()
string(REPLACE "," ";" VERSION_LIST "${VERSION_OUTPUT}")
list(GET VERSION_LIST 0 VERSION)
list(GET VERSION_LIST 1 BUILD_DATE)
list(GET VERSION_LIST 2 GIT_HASH)
set(OS_INFO ${CMAKE_SYSTEM_NAME})
set(ARCH ${CMAKE_SYSTEM_PROCESSOR})

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  @ONLY)

# Messages
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "eOn version: ${VERSION}")
message(STATUS "Platform: ${OS_INFO} ${ARCH}")

# ---------------------- Dependencies
find_package(Threads REQUIRED)
find_package(Eigen3 3.3 REQUIRED)
find_package(spdlog REQUIRED)

set(EON_DEPS Eigen3::Eigen spdlog::spdlog Threads::Threads)

if(NOT WIN32)
  find_library(MATH_LIB m)
  if(MATH_LIB)
    list(APPEND EON_DEPS ${MATH_LIB})
  endif()
endif()

if(EON_NEED_PYTHON_EMBED)
  list(APPEND EON_DEPS Python3::Python Python3::NumPy)
endif()

# ---------------------- eoncbase library
add_library(eoncbase SHARED INIFile.cpp fpe_handler.cpp)
target_include_directories(eoncbase PUBLIC
  ${EON_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(eoncbase PUBLIC ${EON_DEPS})
install(TARGETS eoncbase)

# ---------------------- Potentials
add_subdirectory(potentials)

# ---------------------- eonclib library
set(EONCLIB_SOURCES
  Parameters.cpp
  Optimizer.cpp
  IDPPObjectiveFunction.cpp
  PrefactorJob.cpp
  LBFGS.cpp
  ReplicaExchangeJob.cpp
  BondBoost.cpp
  Job.cpp
  GlobalOptimization.cpp
  LowestEigenmode.cpp
  MinModeSaddleSearch.cpp
  StructureComparisonJob.cpp
  SteepestDescent.cpp
  ImprovedDimer.cpp
  PointJob.cpp
  Prefactor.cpp
  ConjugateGradients.cpp
  Matter.cpp
  FiniteDifferenceJob.cpp
  Lanczos.cpp
  HessianJob.cpp
  TADJob.cpp
  ProcessSearchJob.cpp
  Bundling.cpp
  NEBInitialPaths.cpp
  NudgedElasticBand.cpp
  MonteCarloJob.cpp
  DynamicsJob.cpp
  MonteCarlo.cpp
  Hessian.cpp
  NudgedElasticBandJob.cpp
  DynamicsSaddleSearch.cpp
  HelperFunctions.cpp
  MatrixHelpers.hpp
  Dimer.cpp
  Dynamics.cpp
  GlobalOptimizationJob.cpp
  BiasedGradientSquaredDescent.cpp
  SafeHyperJob.cpp
  MinimizationJob.cpp
  Quickmin.cpp
  ParallelReplicaJob.cpp
  Potential.cpp
  SurrogatePotential.cpp
  BasinHoppingJob.cpp
  FIRE.cpp
  EpiCenters.cpp
  SaddleSearchJob.cpp
  BasinHoppingSaddleSearch.cpp)

# Conditional sources
if(WITH_GPRD)
  list(APPEND EONCLIB_SOURCES AtomicGPDimer.cpp GPRHelpers.cpp)
endif()
if(WITH_GP_SURROGATE)
  list(APPEND EONCLIB_SOURCES GPSurrogateJob.cpp)
endif()

add_library(eonclib SHARED ${EONCLIB_SOURCES})
target_include_directories(eonclib PUBLIC
  ${EON_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(eonclib PUBLIC eoncbase potentials ${EON_DEPS})
install(TARGETS eonclib)

# ---------------------- Compile definitions for options
if(WITH_GPRD)
  target_compile_definitions(eonclib PUBLIC WITH_GPRD)
endif()
if(WITH_GP_SURROGATE)
  target_compile_definitions(eonclib PUBLIC WITH_GP_SURROGATE)
  if(WITH_CATLEARN)
    target_compile_definitions(eonclib PUBLIC WITH_CATLEARN)
  endif()
endif()
if(WITH_VASP)
  target_compile_definitions(eonclib PUBLIC WITH_VASP)
endif()
if(WITH_WATER)
  target_compile_definitions(eonclib PUBLIC WITH_WATER)
endif()
if(WITH_AMS)
  target_compile_definitions(eonclib PUBLIC WITH_AMS)
endif()
if(WITH_LAMMPS)
  target_compile_definitions(eonclib PUBLIC LAMMPS_POT)
endif()
if(WITH_FORTRAN)
  target_compile_definitions(eonclib PUBLIC WITH_FORTRAN)
endif()
if(WITH_CUH2)
  target_compile_definitions(eonclib PUBLIC CUH2_POT WITH_FORTRAN)
endif()
if(WITH_XTB)
  target_compile_definitions(eonclib PUBLIC WITH_XTB)
endif()
if(WITH_METATOMIC)
  target_compile_definitions(eonclib PUBLIC WITH_METATOMIC)
endif()
if(EON_NEED_PYTHON_EMBED)
  target_compile_definitions(eonclib PUBLIC EMBED_PYTHON)
endif()
if(WITH_ASE)
  target_compile_definitions(eonclib PUBLIC ASE_POT)
endif()
if(WITH_ASE_ORCA)
  target_compile_definitions(eonclib PUBLIC WITH_ASE_ORCA)
endif()
if(WITH_ASE_NWCHEM)
  target_compile_definitions(eonclib PUBLIC WITH_ASE_NWCHEM)
endif()
if(WITH_MPI)
  target_compile_definitions(eonclib PUBLIC EONMPI)
endif()

# ---------------------- GPRD subproject
if(WITH_GPRD)
  find_package(gpr_optim QUIET)
  if(gpr_optim_FOUND)
    target_link_libraries(eonclib PUBLIC gpr_optim::gpr_optim)
  else()
    message(WARNING "gpr_optim not found. WITH_GPRD is ON but the library is missing.")
  endif()
endif()

# ---------------------- XTB dependency
if(WITH_XTB)
  find_package(xtb QUIET)
  if(xtb_FOUND)
    target_link_libraries(eonclib PUBLIC xtb::xtb)
  else()
    find_library(XTB_LIBRARY xtb)
    find_path(XTB_INCLUDE_DIR xtb.h)
    if(XTB_LIBRARY)
      target_link_libraries(eonclib PUBLIC ${XTB_LIBRARY})
      if(XTB_INCLUDE_DIR)
        target_include_directories(eonclib PUBLIC ${XTB_INCLUDE_DIR})
      endif()
    else()
      message(FATAL_ERROR "XTB not found. Set xtb_DIR or install xtb.")
    endif()
  endif()
endif()

# ---------------------- Metatomic / Torch dependencies
if(WITH_METATOMIC)
  if(TORCH_PATH)
    set(LIB_TORCH_PATH "${TORCH_PATH}")
  else()
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "import torch; print(torch.utils.cmake_prefix_path)"
      OUTPUT_VARIABLE TORCH_CMAKE_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(TORCH_CMAKE_PREFIX)
      list(APPEND CMAKE_PREFIX_PATH "${TORCH_CMAKE_PREFIX}")
    endif()
  endif()

  find_package(Torch QUIET)
  if(Torch_FOUND)
    target_link_libraries(eonclib PUBLIC ${TORCH_LIBRARIES})
    target_include_directories(eonclib SYSTEM PUBLIC ${TORCH_INCLUDE_DIRS})
  else()
    # Manual fallback: find libraries directly
    if(TORCH_PATH)
      set(LIB_TORCH_LIB_PATH "${TORCH_PATH}/lib")
      set(LIB_TORCH_INCL_PATH "${TORCH_PATH}/include")
      foreach(lib c10 torch torch_cpu torch_global_deps)
        find_library(${lib}_LIB ${lib} PATHS ${LIB_TORCH_LIB_PATH} NO_DEFAULT_PATH)
        if(${lib}_LIB)
          target_link_libraries(eonclib PUBLIC ${${lib}_LIB})
        endif()
      endforeach()
      target_include_directories(eonclib SYSTEM PUBLIC
        ${LIB_TORCH_INCL_PATH}
        ${LIB_TORCH_INCL_PATH}/torch/csrc/api/include)
    else()
      message(FATAL_ERROR "Metatomic requires Torch. Set TORCH_PATH or install PyTorch.")
    endif()
  endif()

  # metatensor and metatomic
  find_library(METATENSOR_LIB metatensor)
  find_library(METATENSOR_TORCH_LIB metatensor_torch)
  find_library(METATOMIC_TORCH_LIB metatomic_torch)
  if(METATENSOR_LIB AND METATENSOR_TORCH_LIB AND METATOMIC_TORCH_LIB)
    target_link_libraries(eonclib PUBLIC
      ${METATENSOR_LIB} ${METATENSOR_TORCH_LIB} ${METATOMIC_TORCH_LIB})
  else()
    message(FATAL_ERROR "metatensor/metatomic libraries not found.")
  endif()

  # vesin
  find_library(VESIN_LIB vesin QUIET)
  if(VESIN_LIB)
    target_link_libraries(eonclib PUBLIC ${VESIN_LIB})
  else()
    # Build vesin from thirdparty
    add_library(vesin_internal STATIC thirdparty/vesin/vesin-single-build.cpp)
    target_include_directories(vesin_internal PUBLIC thirdparty/vesin)
    target_link_libraries(eonclib PUBLIC vesin_internal)
  endif()
endif()

# ---------------------- AMS dependencies
if(WITH_AMS)
  find_package(absl CONFIG REQUIRED)
  find_package(fmt CONFIG REQUIRED)
  find_package(Boost CONFIG REQUIRED COMPONENTS headers iostreams)
  target_link_libraries(eonclib PUBLIC
    absl::base absl::strings absl::numeric absl::utility
    Boost::headers Boost::iostreams
    fmt::fmt)
endif()

# ---------------------- MPI
if(WITH_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(eonclib PUBLIC MPI::MPI_CXX)
endif()

# ---------------------- Executable
add_executable(eonclient ClientEON.cpp CommandLine.cpp)
target_include_directories(eonclient PRIVATE
  ${EON_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(eonclient PRIVATE eonclib)
install(TARGETS eonclient)

# ---------------------- Tests
if(WITH_TESTS)
  enable_testing()
  add_subdirectory(gtests)

  add_executable(testWriteR testWriteR.cpp)
  target_include_directories(testWriteR PRIVATE
    ${EON_INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR})
  target_link_libraries(testWriteR PRIVATE eonclib)
endif()
