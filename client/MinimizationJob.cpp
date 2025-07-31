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
#include "MinimizationJob.h"
#include "BaseStructures.h"
#include "magic_enum/magic_enum.hpp"
namespace eonc {

double MinObjF::getEnergy() const { return _mat.getPotentialEnergy(); }

VectorType MinObjF::getGradient(bool fdstep) const {
  return -_mat.getForcesFreeV();
}

VectorType MinObjF::getPositions() const { return _mat.getPositionsFreeV(); }

int MinObjF::degreesOfFreedom() const { return 3 * _mat.numberOfFreeAtoms(); }

bool MinObjF::isConverged() const { return getConvergence() < _convForce; }

double MinObjF::getConvergence() const {
  switch (_metric) {
  case ObjectiveFunction::Conv::NORM:
    return _mat.getForcesFreeV().norm();
  case ObjectiveFunction::Conv::MAX_ATOM:
    return _mat.maxForce();
  case ObjectiveFunction::Conv::MAX_COMPONENT:
    return _mat.getForces().maxCoeff();
  default:
    SPDLOG_CRITICAL("{} Unknown opt_convergence_metric: {}", "[Minimize]"s,
                    magic_enum::enum_name(_metric));
    std::exit(1);
  }
}

void MinObjF::setPositions(VectorType x) { _mat.setPositionsFreeV(x); }

VectorType MinObjF::difference(VectorType a, VectorType b) {
  return _mat.pbcV(a - b);
}

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

std::vector<std::string> MinimizationJob::run(void) {
  std::string posInFilename("pos.con");
  std::string posOutFilename("min.con");

  if (params->main.checkpoint) {
    FILE *pos_file;
    pos_file = fopen("pos_cp.con", "r");
    if (pos_file != NULL) {
      posInFilename = "pos_cp.con";
      SPDLOG_LOGGER_DEBUG(log, "[Minimization] Resuming from checkpoint");
    } else {
      SPDLOG_LOGGER_DEBUG(log, "[Minimization] No checkpoint files found");
    }
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back(posOutFilename);

  auto pos = Matter(pot, params);
  pos.con2matter(posInFilename);

  SPDLOG_LOGGER_DEBUG(log, "\nBeginning minimization of {}", posInFilename);

  bool converged;
  try {
    converged = pos.relax(false, params->debug.writeMovies,
                          params->main.checkpoint, "minimization", "pos");
    if (converged) {
      status = RunStatus::GOOD;
      SPDLOG_LOGGER_DEBUG(log, "Minimization converged within tolerence");
    } else {
      status = RunStatus::FAIL_MAX_ITERATIONS;
      SPDLOG_LOGGER_DEBUG(log, "Minimization did not converge to tolerence!"
                               "Maybe try to increase max_iterations?");
    }
  } catch (int e) {
    if (e == 100) {
      status = RunStatus::FAIL_POTENTIAL_FAILED;
    } else {
      throw e;
    }
  }

  SPDLOG_LOGGER_DEBUG(log, "Saving result to {}", posOutFilename);
  pos->matter2con(posOutFilename);
  if (status != RunStatus::FAIL_POTENTIAL_FAILED) {
    SPDLOG_LOGGER_DEBUG(log, "Final Energy: {}", pos->getPotentialEnergy());
  }

  std::filesystem::path resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename.string());

  std::ofstream fileResults(resultsFilename, std::ios::binary);

  if (!fileResults.is_open()) {
    std::cerr << "Error opening file " << resultsFilename << ": "
              << std::strerror(errno) << std::endl;
    throw std::runtime_error("Failed to open results file: " +
                             std::string(std::strerror(errno)));
    return returnFiles;
  }

  fileResults << magic_enum::enum_name<RunStatus>(status)
              << " termination_reason\n";
  fileResults << "minimization job_type\n";
  fileResults << magic_enum::enum_name<PotType>(params->potential)
              << " potential_type\n";
  fileResults << this->pot->forceCallCounter << " total_force_calls\n";

  if (status != RunStatus::FAIL_POTENTIAL_FAILED) {
    fileResults << pos->getPotentialEnergy() << " potential_energy\n";
  }

  // No explicit fclose needed; RAII handles it.
  if (!fileResults.good()) {
    std::cerr << "Error writing to file " << resultsFilename
              << ": May be incomplete." << std::endl;
    // Consider throwing, depending on the severity.
  }

  return returnFiles;
}

} // namespace eonc
