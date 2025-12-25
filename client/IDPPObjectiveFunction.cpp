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
