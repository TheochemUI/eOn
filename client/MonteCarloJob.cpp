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
namespace eonc {
std::vector<std::string> MonteCarloJob::run(void) {
  std::string posInFilename("pos.con");
  std::string posOutFilename("out.con");

  if (params->main.checkpoint) {
    FILE *pos;
    pos = fopen("pos_cp.con", "r");
    if (pos != NULL) {
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
  mc.run(params->monte_carlo.steps, params->main.temperature,
         params->monte_carlo.stepSize);

  // FILE *fileResults;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  // fileResults = fopen(resultsFilename.c_str(), "wb");

  // fprintf(fileResults, "%d termination_reason\n", status);
  // fprintf(fileResults, "minimization job_type\n");
  // fprintf(fileResults, "%s potential_type\n",
  // helper_functions::getPotentialName(params->potential).c_str());
  // fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
  // if (status != STATUS_POTENTIAL_FAILED) {
  //     fprintf(fileResults, "%f potential_energy\n",
  //     pos->getPotentialEnergy());
  // }
  // fclose(fileResults);

  return returnFiles;
}

} // namespace eonc
