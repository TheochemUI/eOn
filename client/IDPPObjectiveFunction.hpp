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

#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"
#include <Eigen/Dense>

class IDPPObjectiveFunction : public ObjectiveFunction {
public:
  IDPPObjectiveFunction(std::shared_ptr<Matter> matterPassed,
                        std::shared_ptr<Parameters> paramsPassed,
                        const Eigen::MatrixXd &targetDistances)
      : ObjectiveFunction(matterPassed, paramsPassed),
        d_target(targetDistances) {

    // Initialize working variables to avoid re-allocation
    int natoms = matter->numberOfAtoms();
    d_current.resize(natoms, natoms);
  }

  // IDPP Energy: E = 0.5 * sum( w * (r_ij - d_target_ij)^2 )
  // w = 1 / r_ij^4
  double getEnergy() override;

  // IDPP Gradient
  VectorXd getGradient(bool fdstep = false) override;

  // Standard Interface Plumbing
  void setPositions(VectorXd x) override {
    // Map 3N vector back to Matter
    matter->setPositions(MatrixXd::Map(x.data(), matter->numberOfAtoms(), 3));
  }

  VectorXd getPositions() override {
    // Map Matter positions to 3N vector
    return VectorXd::Map(matter->getPositions().data(),
                         3 * matter->numberOfAtoms());
  }

  int degreesOfFreedom() override { return 3 * matter->numberOfAtoms(); }

  // Use the force convergence criteria from params, or a specific IDPP one if
  // added
  bool isConverged() override {
    return getConvergence() < params->nebConvergedForce;
  }

  double getConvergence() override {
    // Return max force component or norm depending on preference
    // Using norm here for simplicity in path generation
    return getGradient().norm();
  }

  // Handles PBC difference correctly using the Matter object
  VectorXd difference(VectorXd a, VectorXd b) override {
    return matter->pbcV(a - b);
  }

private:
  Eigen::MatrixXd d_target;  // The interpolated "ideal" distances
  Eigen::MatrixXd d_current; // Workspace for current distances
};
