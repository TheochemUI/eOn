/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once

#include <magic_enum/magic_enum_all.hpp>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/os.h> // To write the R style data frame
#include <fmt/ostream.h>
#include <fmt/printf.h>

#include "client/BaseTypes.hpp"

namespace eonc {
using namespace std::string_literals; // For ""s

// TODO(rg): Maybe put these in a namespace
// This file contains forward declarations and enum classes
// NOTE(rg):
// We use magic_enum for converting <-> strings so the names have to match what
// the user will input, though we are case-insensitive

/* Don't guard with compiler directives anymore because that will break ABI for
 * any of these */
enum class PotType {
  // Only add to the end of this!!!
  UNKNOWN = 0,
  EMT,
  EMT_RAS,
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
  EAM_STANDALONE,
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
  CatLearn,
  XTB,
  ASE_ORCA,
  ASE_POT,
  ASE_NWCHEM,
  METATOMIC
};

enum class JobType {
  // Only add to the end of this!!!
  Unknown = 0,
  Process_Search,
  Saddle_Search,
  Minimization,
  Point,
  Parallel_Replica,
  Safe_Hyperdynamics,
  TAD,
  Replica_Exchange,
  Basin_Hopping,
  Hessian,
  Finite_Difference,
  Nudged_Elastic_Band,
  Dynamics,
  Prefactor,
  Global_Optimization,
  Structure_Comparison,
  Monte_Carlo,
  Test,
  GP_Surrogate
};

enum class OptType {
  // Only add to the end of this!!!
  Unknown = -1, // an error case
  None = 0,
  QM,
  CG,
  LBFGS,
  FIRE,
  SD
};

enum class RunStatus {
  NOT_BEGUN = 0,
  GOOD,
  FAIL_MAX_ITERATIONS,
  FAIL_POTENTIAL_FAILED
};

enum class ConvergenceMeasure {
  UNKNOWN = -1, // an error case
  NORM,
  MAX_ATOM,
  MAX_COMPONENT
};

struct BaseOptParams {                                // [Optimizer]
  OptType optM{OptType::CG};                          // opt_method
  ConvergenceMeasure optCM{ConvergenceMeasure::NORM}; // convergence_metric
  ScalarType optConvergedForce{1e-3};                 // converged_force
  size_t optMaxIter{1000L};                           // max_iterations
  ScalarType optMaxMove{0.2};                         // max_move
};

namespace Prefactor {
enum class TYPE { REACTANT, SADDLE, PRODUCT };
enum class RATE { HTST, QQHTST };
enum class FILTER { CUTOFF, FRACTION };
} // namespace Prefactor

} // namespace eonc
