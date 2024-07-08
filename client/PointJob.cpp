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
namespace eonc {
std::vector<std::string> PointJob::run(void) {
  std::vector<std::string> returnFiles;
  std::string posInFilename("pos.con");
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  auto pos = std::make_unique<Matter>(pot, params);
  pos->con2matter(posInFilename);

  SPDLOG_LOGGER_DEBUG(log, "Energy:         {:.12f}",
                      pos->getPotentialEnergy());
  SPDLOG_LOGGER_DEBUG(log, "Max atom force: {:.12f}", pos->maxForce());

  std::shared_ptr<spdlog::logger> fileLogger;
  fileLogger = spdlog::basic_logger_mt("point", "results.dat", true);

  fileLogger->set_pattern("%v");
  fileLogger->info("{:.12f} Energy", pos->getPotentialEnergy());
  fileLogger->info("{:.12f} Max_Force", pos->maxForce());

  spdlog::drop("point");
  fileLogger.reset();
  return returnFiles;
}

} // namespace eonc
