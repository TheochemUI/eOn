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
#include "MonteCarloJob.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "MonteCarlo.h"

#include <filesystem>

std::vector<std::string> MonteCarloJob::run(void) {
  string posInFilename("pos.con");
  string posOutFilename("out.con");

  if (params.main_options.checkpoint) {
    if (std::filesystem::exists("pos_cp.con")) {
      posInFilename = "pos_cp.con";
      SPDLOG_LOGGER_DEBUG(log, "Resuming from checkpoint\n");
    } else {
      SPDLOG_LOGGER_DEBUG(log, "No checkpoint files found\n");
    }
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back(posOutFilename);

  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(posInFilename);

  // code will go
  MonteCarlo mc = MonteCarlo(matter, params);
  mc.run(params.monte_carlo_options.steps, params.main_options.temperature,
         params.monte_carlo_options.step_size);

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  return returnFiles;
}
