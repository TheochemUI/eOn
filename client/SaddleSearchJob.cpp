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
#include "SaddleSearchJob.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "Potential.h"

#include <filesystem>
#include <string>

using namespace std;

std::vector<std::string> SaddleSearchJob::run(void) {
  string reactantFilename("pos.con");
  string displacementFilename("displacement.con");
  string modeFilename("direction.dat");

  if (params.main_options.checkpoint) {
    namespace fs = std::filesystem;
    if (fs::exists("displacement_cp.con") && fs::exists("mode_cp.dat")) {
      displacementFilename = "displacement_cp.con";
      modeFilename = "mode_cp.dat";
      SPDLOG_LOGGER_DEBUG(log, "Resuming from checkpoint");
    } else {
      SPDLOG_LOGGER_DEBUG(log, "No checkpoint files found");
    }
  }

  initial = std::make_shared<Matter>(pot, params);
  displacement = std::make_shared<Matter>(pot, params);
  saddle = std::make_shared<Matter>(pot, params);

  initial->con2matter(reactantFilename);

  if (params.saddle_search_options.displace_type == EpiCenters::DISP_LOAD) {
    // displacement was passed from the server
    saddle->con2matter(displacementFilename);
  } else {
    // displacement and mode will be made on the client
    // in saddleSearch->initialize(...)
    *saddle = *initial;
  }
  AtomMatrix mode;
  if (params.saddle_search_options.displace_type == EpiCenters::DISP_LOAD) {
    // mode was passed from the server
    mode = helper_functions::loadMode(modeFilename, initial->numberOfAtoms());
  }

  saddleSearch = std::make_unique<MinModeSaddleSearch>(
      saddle, mode, initial->getPotentialEnergy(), params, pot);

  int status;
  status = doSaddleSearch();
  printEndState(status);
  saveData(status);

  return returnFiles;
}

int SaddleSearchJob::doSaddleSearch() {
  Matter matterTemp(pot, params);
  long status;
  int f1{0};
  f1 = this->pot->forceCallCounter;
  try {
    status = saddleSearch->run();
  } catch (int e) {
    if (e == 100) {
      status = MinModeSaddleSearch::STATUS_POTENTIAL_FAILED;
    } else {
      printf("unknown exception: %i\n", e);
      throw e;
    }
  }

  if (params.saddle_search_options.minmode_method ==
      LowestEigenmode::MINMODE_GPRDIMER) {
    fCallsSaddle = saddleSearch->forcecalls - f1; // TODO: Check if this
    // works
  } else {
    fCallsSaddle += this->pot->forceCallCounter - f1;
  }

  return status;
}

void SaddleSearchJob::saveData(int status) {
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  {
    auto out = fmt::output_file(resultsFilename);
    out.print("{} termination_reason\n", status);
    out.print("saddle_search job_type\n");
    out.print("{} random_seed\n", params.main_options.randomSeed);
    out.print("{} potential_type\n", std::string{magic_enum::enum_name<PotType>(
                                         params.potential_options.potential)});
    out.print("{} total_force_calls\n", this->pot->forceCallCounter);
    out.print("{} force_calls_saddle\n", fCallsSaddle);
    out.print("{} iterations\n", saddleSearch->iteration);
    if (status != MinModeSaddleSearch::STATUS_POTENTIAL_FAILED) {
      out.print("{:f} potential_energy_saddle\n", saddle->getPotentialEnergy());
      out.print("{:f} final_eigenvalue\n", saddleSearch->getEigenvalue());
    }
    out.print("{:f} potential_energy_reactant\n",
              initial->getPotentialEnergy());
  }

  std::string modeFilename("mode.dat");
  returnFiles.push_back(modeFilename);
  helper_functions::saveMode(modeFilename, saddle,
                             saddleSearch->getEigenvector());

  std::string saddleFilename("saddle.con");
  returnFiles.push_back(saddleFilename);
  saddle->matter2con(saddleFilename);
}

void SaddleSearchJob::printEndState(int status) {
  SPDLOG_LOGGER_DEBUG(log, "[Saddle Search] Final status: ");

  if (status == MinModeSaddleSearch::STATUS_GOOD)
    SPDLOG_LOGGER_DEBUG(log, "Success");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_CONVEX)
    SPDLOG_LOGGER_DEBUG(log,
                        "Initial displacement unable to reach convex region");

  else if (status == MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY)
    SPDLOG_LOGGER_DEBUG(log, "Barrier too high");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MAX_CONCAVE_ITERATIONS)
    SPDLOG_LOGGER_DEBUG(log, "Too many iterations in concave region");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS)
    SPDLOG_LOGGER_DEBUG(log, "Too many iterations");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED)
    SPDLOG_LOGGER_DEBUG(log, "Saddle is not connected to initial state");

  else if (status == MinModeSaddleSearch::STATUS_BAD_PREFACTOR)
    SPDLOG_LOGGER_DEBUG(log, "Prefactors not within window");

  else if (status == MinModeSaddleSearch::STATUS_FAILED_PREFACTOR)
    SPDLOG_LOGGER_DEBUG(log, "Hessian calculation failed");

  else if (status == MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER)
    SPDLOG_LOGGER_DEBUG(log, "Energy barrier not within window");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MINIMA)
    SPDLOG_LOGGER_DEBUG(log, "Minimizations from saddle did not converge");

  else if (status == MinModeSaddleSearch::STATUS_NONNEGATIVE_ABORT)
    SPDLOG_LOGGER_DEBUG(log, "Nonnegative initial mode, aborting");

  else if (status == MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER)
    SPDLOG_LOGGER_DEBUG(log, "Negative barrier detected");

  else if (status == MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT)
    SPDLOG_LOGGER_DEBUG(log, "No reaction found during MD trajectory");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE)
    SPDLOG_LOGGER_DEBUG(
        log, "Converged to stationary point with zero negative modes");

  else if (status == MinModeSaddleSearch::STATUS_BAD_NO_BARRIER)
    SPDLOG_LOGGER_DEBUG(log,
                        "No forward barrier was found along minimized band");

  else if (status == MinModeSaddleSearch::STATUS_ZEROMODE_ABORT)
    SPDLOG_LOGGER_DEBUG(log, "Zero mode abort.");

  else if (status == MinModeSaddleSearch::STATUS_OPTIMIZER_ERROR)
    SPDLOG_LOGGER_DEBUG(log, "Optimizer error.");

  else
    SPDLOG_LOGGER_DEBUG(log, "Unknown status: {}!", status);

  return;
}
