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
#include "eon/MinimizationJob.h"
#include "eon/BaseStructures.h"
#include "eon/HelperFunctions.h"
#include "eon/JobResult.h"
#include "eon/Matter.h"
#include "eon/Optimizer.h"

#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace eonc {

std::vector<std::string> MinimizationJob::run() {
  std::string posInFilename("pos.con");
  std::string posOutFilename("min.con");

  if (params.main_options.checkpoint) {
    if (std::filesystem::exists("pos_cp.con")) {
      posInFilename = "pos_cp.con";
      QUILL_LOG_DEBUG(log, "[Minimization] Resuming from checkpoint");
    } else {
      QUILL_LOG_DEBUG(log, "[Minimization] No checkpoint files found");
    }
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back(posOutFilename);

  auto pos = std::make_shared<Matter>(pot, params);
  if (!eonc::io::io_ok(pos->con2matter(posInFilename))) {
    QUILL_LOG_CRITICAL(log, "Failed to load {}", posInFilename);
    throw std::runtime_error("failed to load " + posInFilename);
  }

  QUILL_LOG_DEBUG(log, "\nBeginning minimization of {}", posInFilename);

  bool converged;
  try {
    converged =
        pos->relax(false, params.debug_options.write_movies,
                   params.main_options.checkpoint, "minimization", "pos");
    if (converged) {
      status = RunStatus::GOOD;
      QUILL_LOG_DEBUG(log, "Minimization converged within tolerence");
    } else {
      status = RunStatus::FAIL_MAX_ITERATIONS;
      QUILL_LOG_DEBUG(log, "Minimization did not converge to tolerence!"
                           "Maybe try to increase max_iterations?");
    }
  } catch (int e) {
    if (e == 100) {
      status = RunStatus::FAIL_POTENTIAL_FAILED;
    } else {
      throw e;
    }
  } catch (const std::exception &e) {
    QUILL_LOG_ERROR(log, "Minimization potential failed: {}", e.what());
    status = RunStatus::FAIL_POTENTIAL_FAILED;
  }

  QUILL_LOG_DEBUG(log, "Saving result to {}", posOutFilename);
  if (!eonc::io::io_ok(pos->matter2con(posOutFilename))) {
    QUILL_LOG_ERROR(log, "Failed to write {}", posOutFilename);
  }
  if (status != RunStatus::FAIL_POTENTIAL_FAILED) {
    QUILL_LOG_DEBUG(log, "Final Energy: {}", pos->getPotentialEnergy());
  }

  std::filesystem::path resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename.string());

  const bool hasE = status != RunStatus::FAIL_POTENTIAL_FAILED;
  const double energy = hasE ? pos->getPotentialEnergy() : 0.0;
  JobResultEnvelope::fromMinimization(
      status, params.potential_options.potential,
      this->pot->forceCallCounter.load(), hasE, energy)
      .writeResultsDat(resultsFilename.string());

  return returnFiles;
}

} // namespace eonc
