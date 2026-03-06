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
#include "PointJob.h"
#include "Matter.h"
#include <format>
#include <fstream>

std::vector<std::string> PointJob::run(void) {
  std::vector<std::string> returnFiles;
  string posInFilename("pos.con");
  string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  auto pos = std::make_unique<Matter>(pot, params);
  pos->con2matter(posInFilename);

  LOG_DEBUG(log, "Energy:         {:.12f}", pos->getPotentialEnergy());
  std::stringstream freeForcesStream;
  freeForcesStream << pos->getForcesFree();
  LOG_DEBUG(log, "(free) Forces:\n{}", freeForcesStream.str());
  LOG_DEBUG(log, "Max atom force: {:.12f}", pos->maxForce());

  std::ofstream outFile("results.dat");
  if (outFile.is_open()) {
    outFile << std::format("{:.12f} Energy\n", pos->getPotentialEnergy());
    outFile << std::format("{:.12f} Max_Force\n", pos->maxForce());
  } else {
    LOG_ERROR(log, "Failed to open results.dat for writing");
  }

  return returnFiles;
}
