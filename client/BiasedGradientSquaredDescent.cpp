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
#include "eon/BiasedGradientSquaredDescent.h"
#include "eon/EigenmodeStrategy.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/ObjectiveFunction.h"
#include "eon/Optimizer.h"
#include "eon/SaddleSearchMethod.h"
#include "eon/SafeMath.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>

namespace eonc {

class BGSDObjectiveFunction : public ObjectiveFunction {
  Matter &matter;

public:
  BGSDObjectiveFunction(Matter &matterRef, double reactantEnergyPassed,
                        double bgsdAlphaPassed,
                        const Parameters &parametersPassed)
      : ObjectiveFunction(parametersPassed),
        matter{matterRef} {
    bgsdAlpha = bgsdAlphaPassed;
    reactantEnergy = reactantEnergyPassed;
  }

  ~BGSDObjectiveFunction() = default;

  double getEnergy() {
    VectorXd Vforce = matter.getForcesFreeV();
    double Henergy = 0.5 * Vforce.dot(Vforce) +
                     0.5 * bgsdAlpha *
                         (matter.getPotentialEnergy() -
                          (reactantEnergy + params.bgsd_options.beta)) *
                         (matter.getPotentialEnergy() -
                          (reactantEnergy + params.bgsd_options.beta));
    return Henergy;
  }

  VectorXd getGradient(bool fdstep = false) {
    (void)fdstep;
    VectorXd Vforce = matter.getForcesFreeV();
    const double magVforce = Vforce.norm();
    const double fd = params.bgsd_options.gradient_finite_difference;
    if (!(magVforce > 0.0) || !std::isfinite(magVforce) || !(fd > 0.0)) {
      return VectorXd::Zero(Vforce.size());
    }
    const VectorXd normVforce = Vforce / magVforce;
    const VectorXd Vpositions = matter.getPositionsFreeV();
    matter.setPositionsFreeV(Vpositions - normVforce * fd);
    VectorXd Vforcenew = matter.getForcesFreeV();
    matter.setPositionsFreeV(Vpositions);
    VectorXd Hforce = magVforce * (Vforcenew - Vforce) / fd +
                      bgsdAlpha *
                          (matter.getPotentialEnergy() -
                           (reactantEnergy + params.bgsd_options.beta)) *
                          Vforce;
    return -Hforce;
  }

  double getGradientnorm() {
    VectorXd Hforce = getGradient();
    double Hnorm = Hforce.norm();
    return Hnorm;
  }

  void setPositions(const VectorXd &x) { matter.setPositionsFreeV(x); }
  VectorXd getPositions() { return matter.getPositionsFreeV(); }
  int degreesOfFreedom() { return 3 * matter.numberOfFreeAtoms(); }
  bool isConverged() { return isConvergedH() && isConvergedV(); }
  bool isConvergedH() {
    return getConvergenceH() < params.bgsd_options.h_force_convergence;
  }
  bool isConvergedV() {
    return getConvergenceV() < params.bgsd_options.grad2energy_convergence;
  }
  bool isConvergedIP() {
    return getConvergenceH() < params.bgsd_options.grad2force_convergence;
  }

  double getConvergence() { return getGradient().norm(); }
  double getConvergenceH() { return getGradient().norm(); }
  double getConvergenceV() { return getEnergy(); }
  VectorXd difference(const VectorXd &a, const VectorXd &b) {
    return matter.pbcV(a - b);
  }

private:
  double reactantEnergy;
  double bgsdAlpha;
};

int BiasedGradientSquaredDescent::run() {
  auto objf = std::make_shared<BGSDObjectiveFunction>(
      *saddle, reactantEnergy, params.bgsd_options.alpha, params);
  auto optim = eonc::helpers::create::mkOptim(
      objf, params.optimizer_options.method, params);
  int iteration = 0;
  const int max_iter = params.optimizer_options.max_iterations;
  QUILL_LOG_DEBUG(
      log,
      "starting optimization of H with params alpha and beta: {:.2f} {:.2f}",
      params.bgsd_options.alpha, params.bgsd_options.beta);
  while (iteration < max_iter && (!objf->isConvergedH() || iteration == 0)) {
    if (!std::isfinite(objf->getEnergy())) {
      break;
    }
    optim->step(params.optimizer_options.max_move);
    QUILL_LOG_DEBUG(log,
                    "iteration {} Henergy, gradientHnorm, and Venergy: "
                    "{:.8f} {:.8f} {:.8f}",
                    iteration, objf->getEnergy(), objf->getGradientnorm(),
                    saddle->getPotentialEnergy());
    iteration++;
  }
  auto objf2 = std::make_shared<BGSDObjectiveFunction>(*saddle, reactantEnergy,
                                                       0.0, params);
  auto optim2 = eonc::helpers::create::mkOptim(
      objf2, params.optimizer_options.method, params);
  int iter2 = 0;
  while (iter2 < max_iter && (!objf2->isConvergedV() || iter2 == 0)) {
    if (objf2->isConvergedIP() || !std::isfinite(objf2->getEnergy())) {
      break;
    }
    optim2->step(params.optimizer_options.max_move);
    QUILL_LOG_DEBUG(log,
                    "gradient squared iteration {} Henergy, gradientHnorm, "
                    "and Venergy: {:.8f} {:.8f} {:.8f}",
                    iteration, objf2->getEnergy(), objf2->getGradientnorm(),
                    saddle->getPotentialEnergy());
    ++iteration;
    ++iter2;
  }

  auto minModeMethod = eonc::buildEigenmodeStrategy(saddle, params, pot);

  eigenvector.setRandom();
  for (int i = 0; i < saddle->numberOfAtoms(); i++) {
    for (int j = 0; j < 3; j++) {
      if (saddle->getFixed(i)) {
        eigenvector(i, j) = 0.0;
      };
    }
  }
  eonc::safemath::safe_normalize_inplace(eigenvector);
  eonc::eigenmodeCompute(*minModeMethod, saddle, eigenvector);
  eigenvector = eonc::eigenmodeGetEigenvector(*minModeMethod);
  eigenvalue = eonc::eigenmodeGetEigenvalue(*minModeMethod);
  QUILL_LOG_DEBUG(log, "lowest eigenvalue {:.8f}", eigenvalue);
  if (objf2->isConvergedV()) {
    return 0;
  } else if (objf2->isConvergedIP()) {
    return 1;
  } else {
    return 1;
  };
}

double BiasedGradientSquaredDescent::getEigenvalue() { return eigenvalue; }

AtomMatrix BiasedGradientSquaredDescent::getEigenvector() {
  return eigenvector;
}

} // namespace eonc
