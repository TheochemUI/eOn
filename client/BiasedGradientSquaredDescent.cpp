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

#include "BiasedGradientSquaredDescent.h"

#include "EigenmodeStrategy.h"

#include "HelperFunctions.h"

#include "Matter.h"

#include "ObjectiveFunction.h"

#include "Optimizer.h"

#include "SaddleSearchMethod.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>

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

  double getEnergy() override {
    VectorXd Vforce = matter.getForcesFreeV();
    double Henergy = 0.5 * Vforce.dot(Vforce) +
                     0.5 * bgsdAlpha *
                         (matter.getPotentialEnergy() -
                          (reactantEnergy + params.bgsd_options.beta)) *
                         (matter.getPotentialEnergy() -
                          (reactantEnergy + params.bgsd_options.beta));
    return Henergy;
  }

  VectorXd getGradient(bool fdstep = false) override {
    VectorXd Vforce = matter.getForcesFreeV();
    double magVforce = Vforce.norm();
    VectorXd normVforce = Vforce / magVforce;
    VectorXd Vpositions = matter.getPositionsFreeV();
    matter.setPositionsFreeV(
        matter.getPositionsFreeV() -
        normVforce * params.bgsd_options.gradient_finite_difference);
    VectorXd Vforcenew = matter.getForcesFreeV();
    matter.setPositionsFreeV(Vpositions);
    VectorXd Hforce = magVforce * (Vforcenew - Vforce) /
                          params.bgsd_options.gradient_finite_difference +
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

  void setPositions(const VectorXd &x) override { matter.setPositionsFreeV(x); }
  VectorXd getPositions() override { return matter.getPositionsFreeV(); }
  int degreesOfFreedom() override { return 3 * matter.numberOfFreeAtoms(); }
  bool isConverged() override { return isConvergedH() && isConvergedV(); }
  bool isConvergedH() {
    return getConvergenceH() < params.bgsd_options.h_force_convergence;
  }
  bool isConvergedV() {
    return getConvergenceV() < params.bgsd_options.grad2energy_convergence;
  }
  bool isConvergedIP() {
    return getConvergenceH() < params.bgsd_options.grad2force_convergence;
  }

  double getConvergence() override {
    return getEnergy() && getGradient().norm();
  }
  double getConvergenceH() { return getGradient().norm(); }
  double getConvergenceV() { return getEnergy(); }
  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return matter.pbcV(a - b);
  }

private:
  double reactantEnergy;
  double bgsdAlpha;
};

SaddleStatus BiasedGradientSquaredDescent::run() {
  auto objf = std::make_shared<BGSDObjectiveFunction>(
      *saddle, reactantEnergy, params.bgsd_options.alpha, params);
  auto optim = eonc::helpers::create::mkOptim(
      objf, params.optimizer_options.method, params);
  int iteration = 0;
  QUILL_LOG_DEBUG(
      log,
      "starting optimization of H with params alpha and beta: {:.2f} {:.2f}",
      params.bgsd_options.alpha, params.bgsd_options.beta);
  while (!objf->isConvergedH() || iteration == 0) {
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
  while (!objf2->isConvergedV() || iteration == 0) {
    if (objf2->isConvergedIP()) {
      break;
    };
    optim2->step(params.optimizer_options.max_move);
    QUILL_LOG_DEBUG(log,
                    "gradient squared iteration {} Henergy, gradientHnorm, "
                    "and Venergy: {:.8f} {:.8f} {:.8f}",
                    iteration, objf2->getEnergy(), objf2->getGradientnorm(),
                    saddle->getPotentialEnergy());
    iteration++;
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
  eigenvector.normalize();
  eonc::eigenmodeCompute(*minModeMethod, saddle, eigenvector);
  eigenvector = eonc::eigenmodeGetEigenvector(*minModeMethod);
  eigenvalue = eonc::eigenmodeGetEigenvalue(*minModeMethod);
  QUILL_LOG_DEBUG(log, "lowest eigenvalue {:.8f}", eigenvalue);
  // Two convergence flavours: V (true convergence on a saddle) and
  // IP (inflection-point fallback). isConvergedIP() and the trailing
  // else both signal "not-V converged" -- preserve the pre-typed-
  // status mapping (0 -> Good, 1 -> Init).
  return objf2->isConvergedV() ? SaddleStatus::Good : SaddleStatus::Init;
}

double BiasedGradientSquaredDescent::getEigenvalue() { return eigenvalue; }

AtomMatrix BiasedGradientSquaredDescent::getEigenvector() {
  return eigenvector;
}
