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
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

std::vector<std::string> MonteCarloJob::run(void) {
  std::string posInFilename("pos.con");
  std::string posOutFilename("out.con");

  if (params.main_options.checkpoint) {
    if (std::filesystem::exists("pos_cp.con")) {
      posInFilename = "pos_cp.con";
      QUILL_LOG_DEBUG(log, "Resuming from checkpoint\n");
    } else {
      QUILL_LOG_DEBUG(log, "No checkpoint files found\n");
    }
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back(posOutFilename);

  auto matter = std::make_shared<Matter>(pot, params);
  if (!eonc::io::io_ok(matter->con2matter(posInFilename))) {
    QUILL_LOG_CRITICAL(log, "Failed to load {}", posInFilename);
    throw std::runtime_error("failed to load " + posInFilename);
  }

  MonteCarlo mc = MonteCarlo(matter, params);
  mc.run(params.monte_carlo_options.steps, params.main_options.temperature,
         params.monte_carlo_options.step_size);

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  std::ofstream out(resultsFilename, std::ios::binary);
  if (out) {
    out << std::format(
        "{} potential_type\n",
        magic_enum::enum_name<PotType>(params.potential_options.potential));
    out << std::format("{} total_force_calls\n",
                       PotRegistry::get().total_force_calls());
    out << std::format("{:f} potential_energy\n", matter->getPotentialEnergy());
  }

  return returnFiles;
}
