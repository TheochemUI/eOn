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

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#include <spdlog/cfg/env.h> // support for loading levels from the environment variable
#include <spdlog/fmt/ostr.h> // support for user defined types
#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

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

enum class RunStatus { GOOD = 0, FAIL_MAX_ITERATIONS, FAIL_POTENTIAL_FAILED };

namespace Prefactor {
enum class TYPE { REACTANT, SADDLE, PRODUCT };
enum class RATE { HTST, QQHTST };
enum class FILTER { CUTOFF, FRACTION };
} // namespace Prefactor
