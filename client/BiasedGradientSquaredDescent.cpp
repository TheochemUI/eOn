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
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
namespace eonc {
class BGSDObjectiveFunction : public ObjectiveFunction {
public:
  BGSDObjectiveFunction(std::shared_ptr<Matter> matterPassed,
                        double reactantEnergyPassed, double bgsdAlphaPassed,
                        std::shared_ptr<Parameters> parametersPassed)
      : ObjectiveFunction(matterPassed, parametersPassed) {
    bgsdAlpha = bgsdAlphaPassed;
    reactantEnergy = reactantEnergyPassed;
  }

  ~BGSDObjectiveFunction(void) = default;

  double getEnergy() {
    VectorType Vforce = matter->getForcesFreeV();
    double Henergy = 0.5 * Vforce.dot(Vforce) +
                     0.5 * bgsdAlpha *
                         (matter->getPotentialEnergy() -
                          (reactantEnergy + parameters->bgsd.beta)) *
                         (matter->getPotentialEnergy() -
                          (reactantEnergy + parameters->bgsd.beta));
    return Henergy;
  }

  VectorType getGradient(bool fdstep = false) {
    VectorType Vforce = matter->getForcesFreeV();
    double magVforce = Vforce.norm();
    VectorType normVforce = Vforce / magVforce;
    VectorType Vpositions = matter->getPositionsFreeV();
    matter->setPositionsFreeV(matter->getPositionsFreeV() -
                              normVforce *
                                  parameters->bgsd.gradientfinitedifference);
    VectorType Vforcenew = matter->getForcesFreeV();
    matter->setPositionsFreeV(Vpositions);
    VectorType Hforce = magVforce * (Vforcenew - Vforce) /
                            parameters->bgsd.gradientfinitedifference +
                        bgsdAlpha *
                            (matter->getPotentialEnergy() -
                             (reactantEnergy + parameters->bgsd.beta)) *
                            Vforce;
    return -Hforce;
  }

  double getGradientnorm() {
    VectorType Hforce = getGradient();
    double Hnorm = Hforce.norm();
    return Hnorm;
  }

  void setPositions(VectorType x) { matter->setPositionsFreeV(x); }
  VectorType getPositions() { return matter->getPositionsFreeV(); }
  int degreesOfFreedom() { return 3 * matter->numberOfFreeAtoms(); }
  bool isConverged() { return isConvergedH() && isConvergedV(); }
  bool isConvergedH() {
    return getConvergenceH() < parameters->bgsd.Hforceconvergence;
  }
  bool isConvergedV() {
    return getConvergenceV() < parameters->bgsd.grad2energyconvergence;
  }
  bool isConvergedIP() {
    return getConvergenceH() < parameters->bgsd.grad2forceconvergence;
  }

  double getConvergence() { return getEnergy() && getGradient().norm(); }
  double getConvergenceH() { return getGradient().norm(); }
  double getConvergenceV() { return getEnergy(); }
  VectorType difference(VectorType a, VectorType b) {
    return matter->pbcV(a - b);
  }

private:
  Matter *matter;
  Parameters *parameters;
  double reactantEnergy;
  double bgsdAlpha;
};

int BiasedGradientSquaredDescent::run() {
  auto objf = std::make_shared<BGSDObjectiveFunction>(
      saddle, reactantEnergy, params->bgsd.alpha, params);
  auto optim = helpers::create::mkOptim(objf, params->optim.method, params);
  int iteration = 0;
  SPDLOG_LOGGER_DEBUG(
      log,
      "starting optimization of H with params alpha and beta: {:.2f} {:.2f}",
      params->bgsd.alpha, params->bgsd.beta);
  while (!objf->isConvergedH() || iteration == 0) {
    optim->step(params->optim.maxMove);
    SPDLOG_LOGGER_DEBUG(log,
                        "iteration {} Henergy, gradientHnorm, and Venergy: "
                        "{:.8f} {:.8f} {:.8f}",
                        iteration, objf->getEnergy(), objf->getGradientnorm(),
                        saddle->getPotentialEnergy());
    iteration++;
  }
  auto objf2 = std::make_shared<BGSDObjectiveFunction>(saddle, reactantEnergy,
                                                       0.0, params);
  auto optim2 = helpers::create::mkOptim(objf2, params->optim.method, params);
  while (!objf2->isConvergedV() || iteration == 0) {
    if (objf2->isConvergedIP()) {
      break;
    };
    optim2->step(params->optim.maxMove);
    SPDLOG_LOGGER_DEBUG(log,
                        "gradient squared iteration {} Henergy, gradientHnorm, "
                        "and Venergy: {:.8f} {:.8f} {:.8f}",
                        iteration, objf2->getEnergy(), objf2->getGradientnorm(),
                        saddle->getPotentialEnergy());
    iteration++;
  }

  LowestEigenmode *minModeMethod;
  if (params->saddle.minmodeMethod == LowestEigenmode::MINMODE_DIMER) {
    if (params->dimer.improved) {
      minModeMethod = new ImprovedDimer(saddle, params, pot);
    } else {
      minModeMethod = new Dimer(saddle, params, pot);
    }
  } else if (params->saddle.minmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
    minModeMethod = new Lanczos(saddle, params, pot);
  }

  //   eigenvector.setZero();
  //   eigenvector(384,0) = 1.;
  eigenvector.setRandom();
  for (int i = 0; i < saddle->numberOfAtoms(); i++) {
    for (int j = 0; j < 3; j++) {
      if (saddle->getFixed(i)) {
        eigenvector(i, j) = 0.0;
      };
    }
  }
  eigenvector.normalize();
  minModeMethod->compute(saddle, eigenvector);
  eigenvector = minModeMethod->getEigenvector();
  eigenvalue = minModeMethod->getEigenvalue();
  SPDLOG_LOGGER_DEBUG(log, "lowest eigenvalue {:.8f}", eigenvalue);
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
