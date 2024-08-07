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

#define HAVE_ALIGNED_ALLOC 1
#define HAVE_POSIX_MEMALIGN 1
#define CACHELOT_PLATFORM_BITS 64
#include <cachelot/cache.h>
#include <cachelot/common.h>

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#include <spdlog/cfg/env.h> // support for loading levels from the environment variable
#include <spdlog/fmt/ostr.h> // support for user defined types
#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

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

namespace Prefactor {
enum class TYPE { REACTANT, SADDLE, PRODUCT };
enum class RATE { HTST, QQHTST };
enum class FILTER { CUTOFF, FRACTION };
} // namespace Prefactor

template <typename T>
inline void hash_combine(std::size_t &seed, const T &value) {
  std::hash<T> hasher;
  seed ^= hasher(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

constexpr size_t cache_memory = 64 * cachelot::Megabyte;
constexpr size_t page_size = 4 * cachelot::Megabyte;
constexpr size_t hash_initial = 131072;
} // namespace eonc
