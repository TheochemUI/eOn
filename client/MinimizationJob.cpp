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

std::vector<std::string> MinimizationJob::run(void) {
  std::string posInFilename("pos.con");
  std::string posOutFilename("min.con");

  if (checkpoint) {
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

  _mat.con2matter(posInFilename);

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

  FILE *fileResults;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");

  fprintf(fileResults, "%s termination_reason\n",
          (std::string{magic_enum::enum_name<RunStatus>(status)}).c_str());
  fprintf(fileResults, "minimization job_type\n");
  fprintf(fileResults, "%s potential_type\n",
          std::string{magic_enum::enum_name<PotType>(params->pot.potential)}
              .c_str());
  // fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
  if (status != RunStatus::FAIL_POTENTIAL_FAILED) {
    fprintf(fileResults, "%f potential_energy\n", pos->getPotentialEnergy());
  }
  fclose(fileResults);

  return returnFiles;
}

bool MinimizationJob::relax(bool quiet, bool writeMovie, bool checkpoint,
                            string prefixMovie, string prefixCheckpoint) {
  auto objf = std::make_shared<MinObjF>(_mat);
  auto optim = helpers::create::mkOptim(objf, OptType::CG, parameters);

  ostringstream min;
  min << prefixMovie;
  if (writeMovie) {
    matter2con(min.str(), false);
  }

  int iteration = 0;
  if (!quiet) {
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10s}  {:14s}  {:18s}  {:13s}\n",
                        "[Matter]", "Iter", "Step size",
                        parameters->optim.convergenceMetricLabel, "Energy");
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}\n",
                        "[Matter]", iteration, 0.0, objf->getConvergence(),
                        getPotentialEnergy());
  }

  while (!objf->isConverged() && iteration < parameters->optim.maxIterations) {

    AtomMatrix pos = _mat.getPositions();

    optim->step(parameters->optim.maxMove);
    iteration++;
    setPositionsFreeV(objf->getPositions());

    double stepSize =
        helper_functions::maxAtomMotion(pbc(getPositions() - pos));

    if (!quiet) {
      SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                          "[Matter]", iteration, stepSize,
                          objf->getConvergence(), getPotentialEnergy());
    }

    if (writeMovie) {
      matter2con(min.str(), true);
    }

    if (checkpoint) {
      ostringstream chk;
      chk << prefixCheckpoint << "_cp";
      matter2con(chk.str(), false);
    }
  }

  if (iteration == 0) {
    if (!quiet) {
      SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                          "[Matter]", iteration, 0.0, objf->getConvergence(),
                          getPotentialEnergy());
    }
  }
  return objf->isConverged();
}
} // namespace eonc
