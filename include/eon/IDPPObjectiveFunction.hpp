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

namespace eonc {

class IDPPObjectiveFunction : public ObjectiveFunction {
  std::shared_ptr<Matter> matter;

public:
  IDPPObjectiveFunction(std::shared_ptr<Matter> matterPassed,
                        const Parameters &paramsPassed,
                        const MatrixXd &targetDistances)
      : ObjectiveFunction(paramsPassed),
        matter{std::move(matterPassed)},
        d_target(targetDistances) {

    // Initialize working variables to avoid re-allocation
    int natoms = matter->numberOfAtoms();
  }

  // IDPP Energy: E = 0.5 * sum( w * (r_ij - d_target_ij)^2 )
  // w = 1 / r_ij^4
  double getEnergy() override;

  // IDPP Gradient
  VectorXd getGradient(bool fdstep = false) override;

  // Free-atom DOF only. Full-3N writes let nearby movers drag frozen atoms
  // (TheochemUI/eOn#410). Same contract as MatterObjectiveFunction.
  void setPositions(const VectorXd &x) override {
    matter->setPositionsFreeV(x);
  }

  VectorXd getPositions() override { return matter->getPositionsFreeV(); }

  int degreesOfFreedom() override {
    return 3 * static_cast<int>(matter->numberOfFreeAtoms());
  }

  bool isConverged() override {
    return getConvergence() < params.neb_options.initialization.force_tolerance;
  }

  double getConvergence() override {
    // Return max force component or norm depending on preference
    // Using norm here for simplicity in path generation
    return getGradient().norm();
  }

  // Handles PBC difference correctly using the Matter object
  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return matter->pbcV(a - b);
  }

private:
  MatrixXd d_target; // The interpolated "ideal" distances
};

class CollectiveIDPPObjectiveFunction : public ObjectiveFunction {
public:
  CollectiveIDPPObjectiveFunction(std::vector<Matter> &pathRef,
                                  const Parameters &paramsPassed)
      : ObjectiveFunction(paramsPassed),
        path(pathRef) {

    // Initialize distances for endpoints
    dInit = getDistanceMatrix(path.front());
    dFinal = getDistanceMatrix(path.back());
  }

  // Return total energy (IDPP + Springs) - Optional for optimization but good
  // for debugging
  double getEnergy() override { return 0.0; }

  VectorXd getGradient(bool fdstep = false) override;

  void setPositions(const VectorXd &x) override {
    const int nfree = static_cast<int>(path[0].numberOfFreeAtoms());
    const int seg = 3 * nfree;
    for (size_t i = 1; i < path.size() - 1; ++i) {
      path[i].setPositionsFreeV(x.segment(seg * static_cast<int>(i - 1), seg));
    }
  }

  VectorXd getPositions() override {
    const int nfree = static_cast<int>(path[0].numberOfFreeAtoms());
    const int seg = 3 * nfree;
    const int n_free_images = static_cast<int>(path.size()) - 2;
    VectorXd pos(seg * n_free_images);
    for (size_t i = 1; i < path.size() - 1; ++i) {
      pos.segment(seg * static_cast<int>(i - 1), seg) =
          path[i].getPositionsFreeV();
    }
    return pos;
  }

  int degreesOfFreedom() override {
    return 3 * static_cast<int>(path[0].numberOfFreeAtoms()) *
           (static_cast<int>(path.size()) - 2);
  }

  // Check convergence of the IDPP-NEB
  bool isConverged() override {
    return getConvergence() < params.neb_options.initialization.force_tolerance;
  }

  double getConvergence() override { return lastMaxForce; }

  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    // Simple difference for this purpose, assuming pre-aligned or handling PBC
    // inside
    return a - b;
  }

private:
  std::vector<Matter> &path;
  MatrixXd dInit, dFinal;
  double lastMaxForce = 100.0;

  MatrixXd getDistanceMatrix(const Matter &m);
  MatrixXd getIDPPForces(const Matter &m, const MatrixXd &dTarget);
};

class ZBLRepulsiveIDPPObjective : public ObjectiveFunction {
public:
  std::shared_ptr<ObjectiveFunction> idpp_obj;
  std::shared_ptr<Potential> zbl_pot;
  std::vector<Matter> &path; // Reference to the actual path vector
  double zbl_weight;

  ZBLRepulsiveIDPPObjective(std::shared_ptr<ObjectiveFunction> idpp,
                            std::shared_ptr<Potential> zbl,
                            std::vector<Matter> &p, const Parameters &params,
                            double weight = 1.0)
      : ObjectiveFunction(params),
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

    const int nfree = static_cast<int>(path[0].numberOfFreeAtoms());
    const int seg = 3 * nfree;
    for (int i = 1; i < n_images - 1; ++i) {
      AtomMatrix forces = MatrixXd::Zero(atoms_per_image, 3);
      double energy = 0;

      zbl_pot->force(atoms_per_image, path[i].getPositions().data(),
                     path[i].getAtomicNrs().data(), forces.data(), &energy,
                     nullptr, path[i].getCell().data());

      VectorXd zbl_free(seg);
      long k = 0;
      for (int a = 0; a < atoms_per_image; ++a) {
        if (!path[i].getFixed(a)) {
          zbl_free[k++] = forces(a, 0);
          zbl_free[k++] = forces(a, 1);
          zbl_free[k++] = forces(a, 2);
        }
      }
      grad.segment((i - 1) * seg, seg) -= (zbl_free * zbl_weight);
    }

    return grad;
  }

  // Delegate other methods
  void setPositions(const VectorXd &x) override { idpp_obj->setPositions(x); }
  VectorXd getPositions() override { return idpp_obj->getPositions(); }
  int degreesOfFreedom() override { return idpp_obj->degreesOfFreedom(); }
  bool isConverged() override { return idpp_obj->isConverged(); }
  double getConvergence() override { return idpp_obj->getConvergence(); }

  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return idpp_obj->difference(a, b);
  }
};

} // namespace eonc

using eonc::CollectiveIDPPObjectiveFunction;
using eonc::IDPPObjectiveFunction;
using eonc::ZBLRepulsiveIDPPObjective;
