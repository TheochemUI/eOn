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
#include "eon/HelperFunctions.h"
#include "eon/JobResult.h"
#include "eon/Matter.h"
#include "eon/PotRegistry.h"

#include <sstream>
#include <stdexcept>

namespace eonc {

std::vector<std::string> PointJob::run() {
  std::vector<std::string> returnFiles;
  std::string posInFilename = eonc::helpers::getRelevantFile("pos.con");
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

  // Energy and Max_Force are the SVN reference format (see
  // data/reference/point_*.dat); the rest is the key set every other job
  // writes and eon.explorer reads.
  auto env = JobResultEnvelope::fromMinimization(
      RunStatus::GOOD, params.potential_options().potential,
      PotRegistry::get().total_force_calls(), true, pos->getPotentialEnergy());
  env.job_type = "point";
  env.extras.emplace_back("Energy", pos->getPotentialEnergy());
  env.extras.emplace_back("Max_Force", pos->maxForce());
  env.writeResultsDat(resultsFilename);
  returnFiles.push_back(resultsFilename);

  return returnFiles;
}

} // namespace eonc
