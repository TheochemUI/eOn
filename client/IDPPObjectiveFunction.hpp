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
#include <cmath>
#include <vector>

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

class CollectiveIDPPObjectiveFunction : public ObjectiveFunction {
public:
  CollectiveIDPPObjectiveFunction(std::vector<Matter> &pathRef,
                                  std::shared_ptr<Parameters> paramsPassed)
      : ObjectiveFunction(
            nullptr, paramsPassed), // We manage a vector, not single Matter
        path(pathRef) {

    // Initialize distances for endpoints
    dInit = getDistanceMatrix(path.front());
    dFinal = getDistanceMatrix(path.back());
  }

  // Return total energy (IDPP + Springs) - Optional for optimization but good
  // for debugging
  double getEnergy() override { return 0.0; }

  // THE CORE LOGIC: mimicking ORCA's NEB force projection
  VectorXd getGradient(bool fdstep = false) override;

  // Plumbing to map the entire path (all images) to one vector
  void setPositions(VectorXd x) override {
    int atoms = path[0].numberOfAtoms();
    // Skip endpoints (0 and N+1)
    for (size_t i = 1; i < path.size() - 1; ++i) {
      path[i].setPositions(MatrixXd::Map(
          x.segment(3 * atoms * (i - 1), 3 * atoms).data(), atoms, 3));
    }
  }

  VectorXd getPositions() override {
    int atoms = path[0].numberOfAtoms();
    int n_free_images = path.size() - 2;
    VectorXd pos(3 * atoms * n_free_images);

    for (size_t i = 1; i < path.size() - 1; ++i) {
      pos.segment(3 * atoms * (i - 1), 3 * atoms) =
          VectorXd::Map(path[i].getPositions().data(), 3 * atoms);
    }
    return pos;
  }

  int degreesOfFreedom() override {
    return 3 * path[0].numberOfAtoms() * (path.size() - 2);
  }

  // Check convergence of the IDPP-NEB
  bool isConverged() override {
      // TODO(rg): parameterize, but tighter tolerance here works better
    return getConvergence() < 0.001;
  }

  double getConvergence() override { return lastMaxForce; }

  VectorXd difference(VectorXd a, VectorXd b) override {
    // Simple difference for this purpose, assuming pre-aligned or handling PBC
    // inside
    return a - b;
  }

private:
  std::vector<Matter> &path;
  Eigen::MatrixXd dInit, dFinal;
  double lastMaxForce = 100.0;

  Eigen::MatrixXd getDistanceMatrix(const Matter &m);
  Eigen::MatrixXd getIDPPForces(const Matter &m,
                                const Eigen::MatrixXd &dTarget);
};
