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
namespace eonc {
bool PointJob::runImpl(void) {
  double energy = _mat.getPotentialEnergy();
  auto maxForce = _mat.maxForce();

  // SPDLOG_LOGGER_DEBUG(log, "(free) Forces:         {:.12}",
  //                     fmt::streamed(_mat.getForcesFree()));
  SPDLOG_LOGGER_DEBUG(log, "Energy:         {:.12f}", energy);
  SPDLOG_LOGGER_DEBUG(log, "Max atom force: {:.12e}", maxForce);

  std::shared_ptr<spdlog::logger> fileLogger;
  fileLogger = spdlog::basic_logger_mt("point", "results.dat", true);

  fileLogger->set_pattern("%v");
  fileLogger->info("{:.12f} Energy", energy);
  fileLogger->info("{:.12e} Max_Force", maxForce);

  spdlog::drop("point");
  fileLogger.reset();
  return true;
}

} // namespace eonc
