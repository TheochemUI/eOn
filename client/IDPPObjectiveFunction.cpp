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

#include "IDPPObjectiveFunction.hpp"

double IDPPObjectiveFunction::getEnergy() {
  double energy = 0.0;
  int natoms = matter->numberOfAtoms();
  AtomMatrix pos = matter->getPositions();

  // Loop over unique pairs
  for (int i = 0; i < natoms; ++i) {
    for (int j = i + 1; j < natoms; ++j) {
      // Respect PBC
      double r = matter->pbc(pos.row(i) - pos.row(j)).norm();

      // Weight function w = 1 / r^4
      // Avoid division by zero if atoms overlap perfectly (unlikely in IDPP but
      // possible)
      if (r < 1e-4)
        r = 1e-4;

      double diff = r - d_target(i, j);
      double weight = 1.0 / std::pow(r, 4);

      energy += 0.5 * weight * diff * diff;
    }
  }
  return energy;
}

VectorXd IDPPObjectiveFunction::getGradient(bool fdstep) {
  int natoms = matter->numberOfAtoms();
  AtomMatrix pos = matter->getPositions();
  AtomMatrix forces = AtomMatrix::Zero(natoms, 3);

  for (int i = 0; i < natoms; ++i) {
    for (int j = i + 1; j < natoms; ++j) {
      // Vector pointing from j to i
      Eigen::RowVector3d dr_vec = matter->pbc(pos.row(i) - pos.row(j));
      double r = dr_vec.norm();

      if (r < 1e-4)
        r = 1e-4;

      double diff = r - d_target(i, j);
      double r2 = r * r;
      double r5 = r2 * r2 * r;

      // Derivative of E_pair = 0.5 * (1/r^4) * (r - d_target)^2
      // dE/dr = (r - d_target)/r^4 - 2(r - d_target)^2 / r^5
      // Simplified: (r - d_target) * (1 - 2(r - d_target)/r) / r^4

      double dEdr = (diff * (1.0 - 2.0 * diff / r)) / std::pow(r, 4);

      // Force contribution: F = -dE/dr * (dr_vec / r)
      Eigen::RowVector3d f_contribution = -dEdr * (dr_vec / r);

      forces.row(i) += f_contribution;
      forces.row(j) -= f_contribution; // Newton's 3rd law
    }
  }

  // Convert N x 3 matrix to 3N vector and return negative gradient (force)
  // BUT getGradient expects the Gradient (positive derivative), so we return
  // -Forces Actually, typical eOn getGradient returns dV/dx.
  return VectorXd::Map(forces.data(), 3 * natoms) * -1.0;
}

Eigen::MatrixXd
CollectiveIDPPObjectiveFunction::getDistanceMatrix(const Matter &m) {
  int natoms = m.numberOfAtoms();
  Eigen::MatrixXd d(natoms, natoms);
  auto pos = m.getPositions();
  for (int i = 0; i < natoms; ++i) {
    for (int j = 0; j < natoms; ++j) {
      d(i, j) = m.pbc(pos.row(i) - pos.row(j)).norm();
    }
  }
  return d;
}

Eigen::MatrixXd
CollectiveIDPPObjectiveFunction::getIDPPForces(const Matter &m,
                                               const Eigen::MatrixXd &dTarget) {
  int natoms = m.numberOfAtoms();
  Eigen::MatrixXd forces = Eigen::MatrixXd::Zero(natoms, 3);
  auto pos = m.getPositions();

  for (int i = 0; i < natoms; ++i) {
    for (int j = i + 1; j < natoms; ++j) {
      Eigen::RowVector3d dr_vec = m.pbc(pos.row(i) - pos.row(j));
      double r = dr_vec.norm();
      if (r < 1e-4)
        r = 1e-4;

      double diff = r - dTarget(i, j);
      // SOTA Weighting: w = 1/r^4. Gradient logic matches Smidstrup/ASE
      double dEdr = (diff / std::pow(r, 5)) * (2.0 * dTarget(i, j) - r);

      Eigen::RowVector3d f = -dEdr * (dr_vec / r); // -dE/dr * r_hat
      forces.row(i) += f;
      forces.row(j) -= f;
    }
  }
  return forces;
}

VectorXd CollectiveIDPPObjectiveFunction::getGradient(bool fdstep) {
  int nImgs = path.size() - 2; // Exclude fixed endpoints
  int natoms = path[0].numberOfAtoms();
  VectorXd totalGradient(3 * natoms * nImgs);
  double maxForce = 0.0;

  // 1. Compute Raw IDPP Forces and Tangents
  std::vector<Eigen::MatrixXd> rawForces(path.size());
  std::vector<Eigen::MatrixXd> tangents(path.size());

  // We compute for 1..N (moving images)
  for (size_t i = 1; i <= nImgs; ++i) {
    // Interpolate Target
    double xi = static_cast<double>(i) / (nImgs + 1);
    Eigen::MatrixXd dTarget = (1.0 - xi) * dInit + xi * dFinal;

    rawForces[i] = getIDPPForces(path[i], dTarget);

    // Simple Tangent: Next - Prev
    Eigen::MatrixXd nextPos = path[i + 1].getPositions();
    Eigen::MatrixXd prevPos = path[i - 1].getPositions();
    tangents[i] = path[i].pbc(nextPos - prevPos);
    tangents[i].normalize(); // Unit tangent
  }

  // 2. Project Forces and Add Springs (The "NEB" part of IDPP-NEB)
  double k = params->neb_options.spring.constant;

  for (size_t i = 1; i <= nImgs; ++i) {
    Eigen::MatrixXd f = rawForces[i];
    Eigen::MatrixXd t = tangents[i];

    // Perpendicular Force (IDPP optimization)
    double f_dot_t = (f.array() * t.array()).sum();
    Eigen::MatrixXd f_perp = f - f_dot_t * t;

    // Spring Force (Spacing optimization)
    double distNext = path[i].distanceTo(path[i + 1]);
    double distPrev = path[i].distanceTo(path[i - 1]);
    Eigen::MatrixXd f_spring = k * (distNext - distPrev) * t;

    // Total NEB Force
    Eigen::MatrixXd f_neb = f_perp + f_spring;

    // Store as Gradient (-Force)
    totalGradient.segment(3 * natoms * (i - 1), 3 * natoms) =
        VectorXd::Map(f_neb.data(), 3 * natoms) * -1.0;

    // Tracking convergence
    maxForce = std::max(maxForce, f_neb.template lpNorm<Eigen::Infinity>());
  }

  lastMaxForce = maxForce;
  return totalGradient;
}
