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

  auto fileLogger = quill::Frontend::create_or_get_logger(
      "point",
      quill::Frontend::create_or_get_sink<quill::FileSink>(
          "results.dat",
          []() {
            quill::FileSinkConfig cfg;
            cfg.set_open_mode('w');
            return cfg;
          }(),
          quill::FileEventNotifier{}),
      quill::PatternFormatterOptions{
          quill::PatternFormatterOptions{quill::PatternFormatterOptions{
              quill::PatternFormatterOptions{"%(message)"}}}},
      quill::ClockSourceType::System);
  LOG_INFO(fileLogger, "{:.12f} Energy", pos->getPotentialEnergy());
  LOG_INFO(fileLogger, "{:.12f} Max_Force", pos->maxForce());
  fileLogger->flush_log();

  return returnFiles;
}
