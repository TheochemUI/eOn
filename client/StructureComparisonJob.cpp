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
#include "eon/StructureComparisonJob.h"
#include "eon/BaseStructures.h"
#include "eon/EonLogger.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "magic_enum/magic_enum.hpp"

#include <cmath>
#include <format>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace eonc {

std::vector<std::string> StructureComparisonJob::run() {
  std::vector<std::string> returnFiles;
  const std::string file1 = eonc::helpers::getRelevantFile("matter1.con");
  const std::string file2 = eonc::helpers::getRelevantFile("matter2.con");

  auto matter1 = std::make_unique<Matter>(pot, params);
  if (!eonc::io::io_ok(matter1->con2matter(file1))) {
    EONC_LOG_CRITICAL("Failed to load {}", file1);
    throw std::runtime_error("failed to load " + file1);
  }
  auto matter2 = std::make_unique<Matter>(pot, params);
  if (!eonc::io::io_ok(matter2->con2matter(file2))) {
    EONC_LOG_CRITICAL("Failed to load {}", file2);
    throw std::runtime_error("failed to load " + file2);
  }

  // Matter::compare can translate the left operand. Probe on a copy.
  Matter probe(*matter1);
  const bool match = probe.compare(
      *matter2, params.structure_comparison_options.indistinguishable_atoms);

  double distance = std::numeric_limits<double>::quiet_NaN();
  double perAtom = std::numeric_limits<double>::quiet_NaN();
  if (matter1->numberOfAtoms() == matter2->numberOfAtoms()) {
    distance = matter1->distanceTo(*matter2);
    perAtom = matter1->perAtomNorm(*matter2);
  }

  const double e1 = matter1->getPotentialEnergy();
  const double e2 = matter2->getPotentialEnergy();

  const std::string resultsFilename("results.dat");
  std::ofstream out(resultsFilename, std::ios::binary);
  if (!out) {
    EONC_LOG_CRITICAL("Failed to open {}", resultsFilename);
    throw std::runtime_error("failed to open " + resultsFilename);
  }
  out << std::format("{} termination_reason\n",
                     static_cast<int>(RunStatus::GOOD));
  out << std::format("{} termination_reason_text\n",
                     magic_enum::enum_name<RunStatus>(RunStatus::GOOD));
  out << "structure_comparison job_type\n";
  out << std::format("{} match\n", match ? 1 : 0);
  out << std::format("{:.12f} distance\n", distance);
  out << std::format("{:.12f} per_atom_norm\n", perAtom);
  out << std::format("{:.12f} energy_1\n", e1);
  out << std::format("{:.12f} energy_2\n", e2);
  out << std::format("{:.12f} energy_abs_diff\n", std::abs(e1 - e2));
  out.close();
  if (!out) {
    throw std::runtime_error("failed to write " + resultsFilename);
  }
  returnFiles.push_back(resultsFilename);
  return returnFiles;
}

} // namespace eonc
