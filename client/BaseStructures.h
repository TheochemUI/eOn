#pragma once
#include <string>

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
  PYTHON
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
  Test
};

namespace helper_functions {
PotType getPotentialType(std::string pname);
std::string getPotentialName(PotType ptype);
JobType getJobType(std::string jname);
std::string getJobName(JobType jtype);
} // namespace helper_functions
