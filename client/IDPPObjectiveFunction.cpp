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

#include "eon/IDPPObjectiveFunction.hpp"

#include <algorithm>
#include <cmath>

namespace eonc {

namespace {
VectorXd packFree(const Matter &m, const AtomMatrix &forces) {
  const long nfree = m.numberOfFreeAtoms();
  AtomMatrix freeF(nfree, 3);
  long k = 0;
  for (long i = 0; i < m.numberOfAtoms(); ++i) {
    if (!m.getFixed(i)) {
      freeF.row(k++) = forces.row(i);
    }
  }
  return VectorXd(VectorXd::Map(freeF.data(), 3 * nfree));
}
} // namespace

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
      double r2 = r * r;
      double weight = 1.0 / (r2 * r2); // 1/r^4

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

      double r4 = r2 * r2;
      double dEdr = (diff * (1.0 - 2.0 * diff / r)) / r4;

      // Force contribution: F = -dE/dr * (dr_vec / r)
      Eigen::RowVector3d f_contribution = -dEdr * (dr_vec / r);

      forces.row(i) += f_contribution;
      forces.row(j) -= f_contribution; // Newton's 3rd law
    }
  }

  // dV/dx on free atoms only. Frozen rows stay in `forces` for Newton's
  // third law during the pair loop, then are dropped.
  return packFree(*matter, forces) * -1.0;
}

MatrixXd CollectiveIDPPObjectiveFunction::getDistanceMatrix(const Matter &m) {
  int natoms = m.numberOfAtoms();
  MatrixXd d(natoms, natoms);
  auto pos = m.getPositions();
  for (int i = 0; i < natoms; ++i) {
    for (int j = 0; j < natoms; ++j) {
      d(i, j) = m.pbc(pos.row(i) - pos.row(j)).norm();
    }
  }
  return d;
}

MatrixXd
CollectiveIDPPObjectiveFunction::getIDPPForces(const Matter &m,
                                               const MatrixXd &dTarget) {
  int natoms = m.numberOfAtoms();
  AtomMatrix forces = AtomMatrix::Zero(natoms, 3);
  auto pos = m.getPositions();

  for (int i = 0; i < natoms; ++i) {
    for (int j = i + 1; j < natoms; ++j) {
      Eigen::RowVector3d dr_vec = m.pbc(pos.row(i) - pos.row(j));
      double r = dr_vec.norm();
      if (r < 1e-4)
        r = 1e-4;

      double diff = r - dTarget(i, j);
      // SOTA Weighting: w = 1/r^4. Gradient logic matches Smidstrup/ASE
      double r2 = r * r;
      double r5 = r2 * r2 * r;
      double dEdr = (diff / r5) * (2.0 * dTarget(i, j) - r);

      Eigen::RowVector3d f = -dEdr * (dr_vec / r); // -dE/dr * r_hat
      forces.row(i) += f;
      forces.row(j) -= f;
    }
  }
  return forces;
}

VectorXd CollectiveIDPPObjectiveFunction::getGradient(bool fdstep) {
  int nImgs = path.size() - 2; // Exclude fixed endpoints
  int nfree = static_cast<int>(path[0].numberOfFreeAtoms());
  VectorXd totalGradient(3 * nfree * nImgs);
  double maxForce = 0.0;

  // 1. Compute Raw IDPP Forces and Tangents
  std::vector<AtomMatrix> rawForces(path.size());
  std::vector<AtomMatrix> tangents(path.size());

  // We compute for 1..N (moving images)
  for (size_t i = 1; i <= nImgs; ++i) {
    // Interpolate Target
    double xi = static_cast<double>(i) / (nImgs + 1);
    MatrixXd dTarget = (1.0 - xi) * dInit + xi * dFinal;

    rawForces[i] = getIDPPForces(path[i], dTarget);

    // Simple Tangent: Next - Prev
    AtomMatrix nextPos = path[i + 1].getPositions();
    AtomMatrix prevPos = path[i - 1].getPositions();
    tangents[i] = path[i].pbc(nextPos - prevPos);
    const double tnorm = tangents[i].norm();
    if (tnorm > 1e-10) {
      tangents[i] /= tnorm;
    }
  }

  // 2. Project Forces and Add Springs (The "NEB" part of IDPP-NEB)
  double k = params.neb_options.spring.constant;

  for (size_t i = 1; i <= nImgs; ++i) {
    AtomMatrix f = rawForces[i];
    AtomMatrix t = tangents[i];

    // Perpendicular Force (IDPP optimization)
    double f_dot_t = matDot(f, t);
    AtomMatrix f_perp = f - f_dot_t * t;

    // Spring Force (Spacing optimization)
    double distNext = path[i].distanceTo(path[i + 1]);
    double distPrev = path[i].distanceTo(path[i - 1]);
    AtomMatrix f_spring = k * (distNext - distPrev) * t;

    AtomMatrix f_neb = f_perp + f_spring;

    VectorXd freeForce = packFree(path[i], f_neb);
    totalGradient.segment(3 * nfree * static_cast<int>(i - 1), 3 * nfree) =
        freeForce * -1.0;

    // Free-atom residuals only. Frozen pair rows stay in f_neb for
    // Newton's third law and must not pin lastMaxForce.
    maxForce = std::max(maxForce, freeForce.lpNorm<Eigen::Infinity>());
  }

  lastMaxForce = maxForce;
  return totalGradient;
}

} // namespace eonc
