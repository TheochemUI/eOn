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
using namespace std;

std::vector<std::string> MonteCarloJob::run(void) {
  string posInFilename("pos.con");
  string posOutFilename("out.con");

  if (params.main_options.checkpoint) {
    FILE *pos;
    pos = fopen("pos_cp.con", "r");
    if (pos != NULL) {
      posInFilename = "pos_cp.con";
      QUILL_LOG_DEBUG(log, "Resuming from checkpoint\n");
    } else {
      QUILL_LOG_DEBUG(log, "No checkpoint files found\n");
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
  FILE *fileResults = fopen(resultsFilename.c_str(), "wb");
  fprintf(fileResults, "%s potential_type\n",
          std::string{magic_enum::enum_name<PotType>(
                          params.potential_options.potential)}
              .c_str());
  fprintf(fileResults, "%zu total_force_calls\n",
          PotRegistry::get().total_force_calls());
  fprintf(fileResults, "%f potential_energy\n", matter->getPotentialEnergy());
  fclose(fileResults);

  return returnFiles;
}
