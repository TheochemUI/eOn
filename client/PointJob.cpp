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
#include "eon/PointJob.h"
#include "eon/BaseStructures.h"
#include "eon/Matter.h"
#include "eon/PotRegistry.h"
#include "magic_enum/magic_enum.hpp"

#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<std::string> PointJob::run() {
  std::vector<std::string> returnFiles;
  std::string posInFilename("pos.con");
  std::string resultsFilename("results.dat");

  auto pos = std::make_unique<Matter>(pot, params);
  if (!eonc::io::io_ok(pos->con2matter(posInFilename))) {
    QUILL_LOG_CRITICAL(log, "Failed to load {}", posInFilename);
    throw std::runtime_error("failed to load " + posInFilename);
  }

  QUILL_LOG_DEBUG(log, "Energy:         {:.12f}", pos->getPotentialEnergy());
  std::stringstream freeForcesStream;
  freeForcesStream << pos->getForcesFree();
  QUILL_LOG_DEBUG(log, "(free) Forces:\n{}", freeForcesStream.str());
  QUILL_LOG_DEBUG(log, "Max atom force: {:.12f}", pos->maxForce());

  std::ofstream outFile(resultsFilename, std::ios::binary);
  if (!outFile) {
    QUILL_LOG_CRITICAL(log, "Failed to open {}", resultsFilename);
    throw std::runtime_error("failed to open " + resultsFilename);
  }
  // Energy and Max_Force are the SVN reference format (see
  // data/reference/point_*.dat); the rest is the key set every other job
  // writes and eon.explorer reads.
  outFile << std::format("{} termination_reason\n",
                         static_cast<int>(RunStatus::GOOD));
  outFile << std::format("{} termination_reason_text\n",
                         magic_enum::enum_name<RunStatus>(RunStatus::GOOD));
  outFile << "point job_type\n";
  outFile << std::format(
      "{} potential_type\n",
      magic_enum::enum_name<PotType>(params.potential_options.potential));
  outFile << std::format("{} total_force_calls\n",
                         PotRegistry::get().total_force_calls());
  outFile << std::format("{:.12f} potential_energy\n",
                         pos->getPotentialEnergy());
  outFile << std::format("{:.12f} Energy\n", pos->getPotentialEnergy());
  outFile << std::format("{:.12f} Max_Force\n", pos->maxForce());
  outFile.close();
  if (!outFile) {
    QUILL_LOG_CRITICAL(log, "Failed to write {}", resultsFilename);
    throw std::runtime_error("failed to write " + resultsFilename);
  }
  returnFiles.push_back(resultsFilename);

  return returnFiles;
}
