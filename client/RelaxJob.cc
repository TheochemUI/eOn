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
#include "client/RelaxJob.hpp"
#include "client/HelperFunctions.h"
#include "client/Optimizer.h"

namespace eonc {
double MatterObjectiveFunction::getEnergy() const {
  return matter.getPotentialEnergy();
}

VectorType MatterObjectiveFunction::getGradient(bool fdstep) const {
  return -matter.getForcesFreeV();
}

void MatterObjectiveFunction::setPositions(const VectorType &x) {
  // Only mutable point
  const_cast<Matter &>(matter).setPositionsFreeV(x);
}

VectorType MatterObjectiveFunction::getPositions() const {
  return matter.getPositionsFreeV();
}

int MatterObjectiveFunction::degreesOfFreedom() const {
  return 3 * matter.numberOfFreeAtoms();
}

bool MatterObjectiveFunction::isConverged() const {
  return getConvergence() < convForce;
}

double MatterObjectiveFunction::getConvergence() const {
  if (convMetric == "norm") {
    return matter.getForcesFreeV().norm();
  } else if (convMetric == "max_atom") {
    return matter.getMaxForce();
  } else if (convMetric == "max_component") {
    return matter.getForcesFree().maxCoeff();
  } else {
    SPDLOG_CRITICAL("{} Unknown opt_convergence_metric: {}", "[Matter]",
                    convMetric);
    std::exit(1);
  }
}

VectorType MatterObjectiveFunction::difference(const VectorType &a,
                                               const VectorType &b) const {
  return matter.pbcV(a - b);
}

bool RelaxJob::runImpl(Matter &mat) {
  // TODO(rg): params are from toml
  eonc::MatterObjectiveFunction objf({"norm", 1e-3}, mat);

  const auto config = toml::table{{"Optimizer", toml::table{{"method", "cg"}}}};
  auto optim = helpers::create::mkOptim(objf, config);

  // std::ostringstream min;
  // min << prefixMovie;
  // if (writeMovie) {
  //   matter2con(min.str(), false);
  // }
  bool quiet = false;
  const size_t maxIter = 1000;
  const double maxMove = 0.2;
  int iteration = 0;
  if (!quiet) {
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10s}  {:14s}  {:18s}  {:13s}", "[Matter]",
                        "Iter", "Step size", "norm", "Energy");
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                        "[Matter]", iteration, 0.0, objf.getConvergence(),
                        objf.getEnergy());
  }

  while (!objf.isConverged() && iteration < maxIter) {

    VectorType pos = objf.getPositions();

    optim->step(maxMove);
    iteration++;
    mat.setPositionsFreeV(objf.getPositions());

    double stepSize = helper_functions::maxAtomMotionV(
        mat.pbcV(mat.getPositionsFreeV() - pos));

    if (!quiet) {
      SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                          "[Matter]", iteration, stepSize,
                          objf.getConvergence(), objf.getEnergy());
    }

    // if (writeMovie) {
    //   matter2con(min.str(), true);
    // }

    // if (checkpoint) {
    //   ostringstream chk;
    //   chk << prefixCheckpoint << "_cp";
    //   matter2con(chk.str(), false);
    // }
  }

  if (iteration == 0) {
    if (!quiet) {
      SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                          "[Matter]", iteration, 0.0, objf.getConvergence(),
                          objf.getEnergy());
    }
  }
  //    bool converged = optimizer->run(parameters->optMaxIterations,
  //    parameters->optMaxMove);
  return objf.isConverged();
  return true;
}

} // namespace eonc
