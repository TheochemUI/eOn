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
#include <stdexcept>
#include <string>

#include "eon/BaseStructures.h"
#include "eon/Dynamics.h"
#include "eon/DynamicsJob.h"
#include "eon/EonLogger.h"
#include "eon/HelperFunctions.h"
#include "eon/JobResult.h"
#include "eon/Parameters.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"

namespace eonc {

std::vector<std::string> DynamicsJob::run(void) {
  auto R = std::make_shared<Matter>(pot, params);
  const std::string posFile = eonc::helpers::getRelevantFile("pos.con");
  if (!eonc::io::io_ok(R->con2matter(posFile))) {
    EONC_LOG_CRITICAL("Failed to load {}", posFile);
    throw std::runtime_error("failed to load " + posFile);
  }

  auto d = std::make_unique<Dynamics>(R.get(), params);
  d->run();

  std::string productFilename("final.con");
  if (!eonc::io::io_ok(R->matter2con(productFilename))) {
    EONC_LOG_ERROR("Failed to write {}", productFilename);
  }

  const std::string resultsFilename("results.dat");
  auto env = JobResultEnvelope::fromMinimization(
      RunStatus::GOOD, params.potential_options().potential,
      PotRegistry::get().total_force_calls(), true, R->getPotentialEnergy());
  env.job_type = "dynamics";
  env.writeResultsDat(resultsFilename);

  std::vector<std::string> returnFiles;
  returnFiles.push_back(productFilename);
  returnFiles.push_back(resultsFilename);
  return returnFiles;
}

} // namespace eonc
