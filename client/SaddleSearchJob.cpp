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
#ifdef WITH_ARTN
#include "ARTnSaddleSearch.h"
#endif
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "Potential.h"

#include <filesystem>
#include <format>
#include <fstream>
#include <string>

std::vector<std::string> SaddleSearchJob::run() {
  std::string reactantFilename("pos.con");
  std::string displacementFilename("displacement.con");
  std::string modeFilename("direction.dat");

  if (params.main_options.checkpoint) {
    if (std::filesystem::exists("displacement_cp.con") &&
        std::filesystem::exists("mode_cp.dat")) {
      displacementFilename = "displacement_cp.con";
      modeFilename = "mode_cp.dat";
      QUILL_LOG_DEBUG(log, "Resuming from checkpoint");
    } else {
      QUILL_LOG_DEBUG(log, "No checkpoint files found");
    }
  }

  initial = std::make_shared<Matter>(pot, params);
  displacement = std::make_shared<Matter>(pot, params);
  saddle = std::make_shared<Matter>(pot, params);

  initial->con2matter(reactantFilename);

  const bool standaloneARTn = params.saddle_search_options.method == "artn";

  if (!standaloneARTn &&
      params.saddle_search_options.displace_type == eonc::EpiCenters::DISP_LOAD) {
    saddle->con2matter(displacementFilename);
  } else {
    *saddle = *initial;
  }

  AtomMatrix mode = AtomMatrix::Zero(initial->numberOfAtoms(), 3);
  const bool canLoadMode =
      params.saddle_search_options.displace_type == eonc::EpiCenters::DISP_LOAD;
  if (canLoadMode && std::filesystem::exists(modeFilename)) {
    mode = eonc::helpers::loadMode(modeFilename, initial->numberOfAtoms());
  }

#ifdef WITH_ARTN
  if (params.saddle_search_options.method == "artn" ||
      params.saddle_search_options.minmode_method == "artn") {
    saddleSearch =
        std::make_unique<ARTnSaddleSearch>(saddle, pot, mode, params);
  } else
#endif
  {
    saddleSearch = std::make_unique<MinModeSaddleSearch>(
        saddle, mode, initial->getPotentialEnergy(), params, pot);
  }

  int status = doSaddleSearch();
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
          LowestEigenmode::MINMODE_GPRDIMER ||
      params.saddle_search_options.method == "artn" ||
      params.saddle_search_options.minmode_method == "artn") {
    fCallsSaddle = saddleSearch->forcecalls;
  } else {
    fCallsSaddle += this->pot->forceCallCounter - f1;
  }

  return status;
}

void SaddleSearchJob::saveData(int status) {
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  std::ofstream out(resultsFilename, std::ios::binary);
  if (out) {
    out << std::format("{} termination_reason\n", status);
    out << std::format("{} termination_reason_text\n",
                       saddleSearch->describeStatus(status));
    out << "saddle_search job_type\n";
    out << std::format("{} random_seed\n", params.main_options.randomSeed);
    out << std::format(
        "{} potential_type\n",
        magic_enum::enum_name<PotType>(params.potential_options.potential));
    out << std::format("{} total_force_calls\n",
                       this->pot->forceCallCounter.load());
    out << std::format("{} force_calls_saddle\n", fCallsSaddle);
    out << std::format("{} iterations\n", saddleSearch->iteration);
    if (status != MinModeSaddleSearch::STATUS_POTENTIAL_FAILED) {
      out << std::format("{:f} potential_energy_saddle\n",
                         saddle->getPotentialEnergy());
      out << std::format("{:f} final_eigenvalue\n",
                         saddleSearch->getEigenvalue());
    }
    out << std::format("{:f} potential_energy_reactant\n",
                       initial->getPotentialEnergy());
  }

  std::string modeFilename("mode.dat");
  returnFiles.push_back(modeFilename);
  {
    std::ofstream modeOut(modeFilename, std::ios::binary);
    if (modeOut) {
      auto eigenvec = saddleSearch->getEigenvector();
      for (long row = 0; row < eigenvec.rows(); ++row) {
        modeOut << std::format("{:12.6f} {:12.6f} {:12.6f}\n", eigenvec(row, 0),
                               eigenvec(row, 1), eigenvec(row, 2));
      }
    }
  }

  std::string saddleFilename("saddle.con");
  returnFiles.push_back(saddleFilename);
  saddle->matter2con(saddleFilename);
}

void SaddleSearchJob::printEndState(int status) {
  auto msg = saddleSearch->describeStatus(status);
  if (status == MinModeSaddleSearch::STATUS_GOOD) {
    QUILL_LOG_DEBUG(log, "[Saddle Search] {}", msg);
  } else {
    QUILL_LOG_WARNING(log, "[Saddle Search] {}", msg);
  }
}
