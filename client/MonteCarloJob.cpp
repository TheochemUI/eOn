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
#include "eon/MonteCarloJob.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/MonteCarlo.h"

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

  auto matter = std::make_shared<Matter>(pot, params);
  if (!eonc::io::io_ok(matter->con2matter(posInFilename))) {
    QUILL_LOG_CRITICAL(log, "Failed to load {}", posInFilename);
    throw std::runtime_error("failed to load " + posInFilename);
  }

  MonteCarlo mc = MonteCarlo(matter, params);
  mc.run(params.monte_carlo_options.steps, params.main_options.temperature,
         params.monte_carlo_options.step_size);

  QUILL_LOG_DEBUG(log, "Saving result to {}", posOutFilename);
  if (eonc::io::io_ok(matter->matter2con(posOutFilename))) {
    returnFiles.push_back(posOutFilename);
  } else {
    QUILL_LOG_ERROR(log, "Failed to write {}", posOutFilename);
  }

  std::string resultsFilename("results.dat");

  std::ofstream out(resultsFilename, std::ios::binary);
  if (!out) {
    QUILL_LOG_CRITICAL(log, "Failed to open {}", resultsFilename);
    throw std::runtime_error("failed to open " + resultsFilename);
  }
  out << std::format(
      "{} potential_type\n",
      magic_enum::enum_name<PotType>(params.potential_options.potential));
  out << std::format("{} total_force_calls\n",
                     PotRegistry::get().total_force_calls());
  out << std::format("{:f} potential_energy\n", matter->getPotentialEnergy());
  out.close();
  if (!out) {
    QUILL_LOG_CRITICAL(log, "Failed to write {}", resultsFilename);
    throw std::runtime_error("failed to write " + resultsFilename);
  }
  returnFiles.push_back(resultsFilename);

  return returnFiles;
}
