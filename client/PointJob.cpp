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

#include <sstream>
#include <stdexcept>

namespace eonc {

std::vector<std::string> PointJob::run() {
  std::vector<std::string> returnFiles;
  std::string posInFilename = eonc::helpers::getRelevantFile("pos.con");
  std::string resultsFilename("results.dat");

  std::unique_ptr<Matter> pos;
  if (ownedPot) {
    pos = std::make_unique<Matter>(std::move(ownedPot), params);
  } else if (pot) {
    pos = std::make_unique<Matter>(pot, params);
  } else {
    throw std::runtime_error("PointJob: no potential");
  }
  if (!eonc::io::io_ok(pos->con2matter(posInFilename))) {
    QUILL_LOG_CRITICAL(log, "Failed to load {}", posInFilename);
    throw std::runtime_error("failed to load " + posInFilename);
  }

  // A point job evaluates the selected potential even when the input
  // contains a saved energy and force evaluation.
  pos->setPositions(pos->getPositions());
  const double energy = pos->getPotentialEnergy();
  QUILL_LOG_DEBUG(log, "Energy:         {:.12f}", energy);
  std::stringstream freeForcesStream;
  freeForcesStream << pos->getForcesFree();
  QUILL_LOG_DEBUG(log, "(free) Forces:\n{}", freeForcesStream.str());
  const double maxForce = pos->maxForce();
  QUILL_LOG_DEBUG(log, "Max atom force: {:.12f}", maxForce);
  const size_t calls = pos->getPotential()->registryForceCalls();

  // Energy and Max_Force are the SVN reference format (see
  // data/reference/point_*.dat); the rest is the key set every other job
  // writes and eon.explorer reads.
  auto env = JobResultEnvelope::fromMinimization(
      RunStatus::GOOD, params.potential_options().potential, calls, true,
      energy);
  env.provenance = provenanceForJob(params);
  env.job_type = "point";
  env.extras.emplace_back("Energy", energy);
  env.extras.emplace_back("Max_Force", maxForce);
  env.writeResultsDat(resultsFilename);
  returnFiles.push_back(resultsFilename);

  return returnFiles;
}

} // namespace eonc
