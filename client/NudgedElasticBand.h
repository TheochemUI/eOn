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
#pragma once

#include <filesystem>

#include "Eigen.h"
#include "EonLogger.h"

#include "EigenmodeStrategy.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "NEBProjection.h"
#include "NEBTangent.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"

namespace eonc::neb {
class OCINEBController;
}

namespace eonc {

// NEB method for determining a minimum energy path between two matter objects
class NudgedElasticBand {
  friend class eonc::neb::OCINEBController;

public:
  enum class NEBStatus {
    GOOD = 0,
    INIT = 1,
    BAD_MAX_ITERATIONS = 2,
    RUNNING,
    MAX_UNCERTAINTY
  };
  NudgedElasticBand(std::shared_ptr<Matter> initialPassed,
                    std::shared_ptr<Matter> finalPassed,
                    const Parameters &parametersPassed,
                    std::shared_ptr<Potential> potPassed);
  NudgedElasticBand(std::vector<Matter> initPath,
                    const Parameters &parametersPassed,
                    std::shared_ptr<Potential> potPassed);
  ~NudgedElasticBand() = default;

  NudgedElasticBand::NEBStatus compute(void);
  NudgedElasticBand::NEBStatus getStatus() { return this->status; };
  void updateForces(bool ci_active);
  void updateForces(void) { updateForces(ci_enabled_); }
  void setCIEnabled(bool enabled) { ci_enabled_ = enabled; }
  double convergenceForce(void);
  void findExtrema(void);
  void printImageData(bool writeToFile = false, size_t idx = 0);
  std::vector<std::shared_ptr<EigenmodeStrategy>> eigenmode_solvers;

  int atoms{0};
  long numImages{0}, climbingImage{0}, numExtrema{0};
  std::vector<std::shared_ptr<Matter>> path; // NEB images
  std::vector<std::shared_ptr<AtomMatrix>> tangent;
  std::vector<std::shared_ptr<AtomMatrix>> projectedForce;
  std::vector<double> extremumEnergy;
  std::vector<double> extremumPosition;
  std::vector<double> extremumCurvature;

  std::size_t maxEnergyImage{0};
  bool movedAfterForceCall{false};
  bool perImagePotentials_{false}; ///< Whether per-image potential instances exist
  double ksp{0.0};
  double k_u{0.0}; // Upper-bound value for the spring constant
  double k_l{0.0}; // Lower-bound value for the spring constant
  double E_ref; // Reference energy chosen to be equal to the max energy of the
                // reactant or product energy minimum

private:
  bool ci_enabled_{false}; // runtime CI state, set by compute()
  double baseline_force{-1.0};
  Parameters params;
  std::shared_ptr<Potential> pot;
  NEBStatus status;
  eonc::log::Scoped log;

  // Cached strategies (constant across iterations)
  neb::TangentStrategy tangentStrat_;
  neb::ProjectionStrategy projectionStrat_;
};

class NEBObjectiveFunction : public ObjectiveFunction {
public:
  NEBObjectiveFunction(NudgedElasticBand *nebPassed,
                       const Parameters &parametersPassed)
      : ObjectiveFunction(parametersPassed),
        neb{nebPassed} {}
  // This is the odd one out, doesn't take a Matter so we null it

  ~NEBObjectiveFunction(void) {};

  VectorXd getGradient(bool fdstep = false);
  double getEnergy();
  void setPositions(const VectorXd &x);
  VectorXd getPositions();
  int degreesOfFreedom();
  bool isConverged();
  bool isUncertain();
  double getConvergence();
  VectorXd difference(const VectorXd &a, const VectorXd &b);
  NudgedElasticBand::NEBStatus status;

private:
  NudgedElasticBand *neb;
};

} // namespace eonc

using eonc::NEBObjectiveFunction;
using eonc::NudgedElasticBand;
