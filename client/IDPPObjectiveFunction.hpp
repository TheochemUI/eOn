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

  bool isConverged() override {
    return getConvergence() <
           params->neb_options.initialization.force_tolerance;
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
  Eigen::MatrixXd d_target; // The interpolated "ideal" distances
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
    return getConvergence() <
           params->neb_options.initialization.force_tolerance;
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

class ZBLRepulsiveIDPPObjective : public ObjectiveFunction {
public:
  std::shared_ptr<ObjectiveFunction> idpp_obj;
  std::shared_ptr<Potential> zbl_pot;
  std::vector<Matter> &path; // Reference to the actual path vector
  double zbl_weight;

  ZBLRepulsiveIDPPObjective(std::shared_ptr<ObjectiveFunction> idpp,
                            std::shared_ptr<Potential> zbl,
                            std::vector<Matter> &p,
                            std::shared_ptr<Parameters> params,
                            double weight = 1.0)
      : ObjectiveFunction(nullptr, params),
        idpp_obj(idpp),
        zbl_pot(zbl),
        path(p),
        zbl_weight(weight) {}

  double getEnergy() override {
    // IDPP "Energy" (Residual) + ZBL Energy
    return idpp_obj->getEnergy();
  }

  VectorXd getGradient(bool fdstep = false) override {
    // 1. Get IDPP Gradient (forces atoms towards interpolated distances)
    VectorXd grad = idpp_obj->getGradient(fdstep);

    // 2. Calculate ZBL Forces for every image
    int n_images = path.size();
    int atoms_per_image = path[0].numberOfAtoms();

    // ZBL calculation loop
    for (int i = 1; i < n_images - 1; ++i) { // Skip endpoints
      AtomMatrix forces = MatrixXd::Zero(atoms_per_image, 3);
      double energy = 0;

      // Calculate ZBL forces for this image
      zbl_pot->force(atoms_per_image, path[i].getPositions().data(),
                     path[i].getAtomicNrs().data(), forces.data(), &energy,
                     nullptr, path[i].getCell().data());

      int segment_start = (i - 1) * 3 * atoms_per_image;
      VectorXd zbl_grad_vec = VectorXd::Map(forces.data(), 3 * atoms_per_image);

      // Add repulsive push (negate force to get gradient)
      grad.segment(segment_start, 3 * atoms_per_image) -=
          (zbl_grad_vec * zbl_weight);
    }

    return grad;
  }

  // Delegate other methods
  void setPositions(VectorXd x) override { idpp_obj->setPositions(x); }
  VectorXd getPositions() override { return idpp_obj->getPositions(); }
  int degreesOfFreedom() override { return idpp_obj->degreesOfFreedom(); }
  bool isConverged() override { return idpp_obj->isConverged(); }
  double getConvergence() override { return idpp_obj->getConvergence(); }

  VectorXd difference(VectorXd a, VectorXd b) override {
    return idpp_obj->difference(a, b);
  }
};
