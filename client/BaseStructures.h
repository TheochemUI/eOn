#pragma once
#include <memory>
#include <string>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <fmt/ostream.h>
#include <fmt/os.h> // To write the R style data frame

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#include <spdlog/spdlog.h>
#include <spdlog/cfg/env.h>  // support for loading levels from the environment variable
#include <spdlog/fmt/ostr.h> // support for user defined types
#include <spdlog/logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

using namespace std::string_literals; // For ""s

// This file contains forward declarations and enum classes

/* Don't guard with compiler directives anymore because that will break ABI for
 * any of these */
enum class PotType {
  // Only add to the end of this!!!
  UNKNOWN = 0,
  EMT,
  EXT,
  LJ,
  LJCLUSTER,
  MORSE_PT,
  NEW,
  CUH2,
  IMD,
  TIP4P,
  TIP4P_PT,
  TIP4P_H,
  SPCE,
  EAM_AL,
  EDIP,
  FEHE,
  LENOSKY_SI,
  SW_SI,
  TERSOFF_SI,
  VASP,
  LAMMPS,
  MPI,
  PYAMFF,
  QSC,
  BOPFOX, // unused?
  BOP,    // unused?
  // Add newer entries here
  AMS,
  AMS_IO,
  GPR,
  PYTHON,
  PYSURROGATE
};

enum class JobType {
  // Only add to the end of this!!!
  Unknown = 0,
  ProcessSearch,
  SaddleSearch,
  Minimization,
  Point,
  ParallelReplica,
  SafeHyperdynamics,
  TAD,
  ReplicaExchange,
  BasinHopping,
  Hessian,
  FiniteDifference,
  NEB,
  Dynamics,
  Prefactor,
  GlobalOptimization,
  StructureComparison,
  MonteCarlo,
  Test,
  GPSurrogate
};

enum class RunStatus { GOOD = 0, MAX_ITERATIONS, POTENTIAL_FAILED };

namespace helper_functions {
PotType getPotentialType(std::string pname);
std::string getPotentialName(PotType ptype);
JobType getJobType(std::string jname);
std::string getJobName(JobType jtype);
std::string getRunStatusName(RunStatus rstype);
} // namespace helper_functions
