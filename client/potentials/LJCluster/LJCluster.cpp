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

#include "eon/potentials/LJCluster/LJCluster.h"
#include "eon/VesinNeighbors.h"
#include <cmath>

void LJCluster::cleanMemory(void) { return; }

void LJCluster::setParameters(double u0Recieved, double cuttOffRRecieved,
                              double psiRecieved) {
  u0 = u0Recieved;
  psi = psiRecieved;

  cuttOffR = cuttOffRRecieved;
  cuttOffU = 4 * u0 * (pow(psi / cuttOffR, 12) - pow(psi / cuttOffR, 6));
  return;
}

void LJCluster::force(long N, const double *R, const int * /*atomicNrs*/,
                      double *F, double *U, double *variance,
                      const double *box) {
  variance = nullptr;
  *U = 0;
  for (int i = 0; i < N; i++) {
    F[3 * i] = 0;
    F[3 * i + 1] = 0;
    F[3 * i + 2] = 0;
  }
  if (N < 2) {
    return;
  }

  // Cluster: free boundary (no PBC). Large cutoff if cuttOffR unused
  // historically.
  thread_local eonc::VesinNeighbors nl;
  eonc::VesinNeighbors::Options opt;
  opt.cutoff = (cuttOffR > 0.0) ? cuttOffR : 1.0e6;
  opt.full = false;
  opt.return_distances = true;
  opt.return_vectors = true;
  opt.periodic = {{false, false, false}};
  // Dummy orthorhombic box (ignored when non-periodic)
  double free_box[9] = {1e6, 0, 0, 0, 1e6, 0, 0, 0, 1e6};
  const double *box_use = (box != nullptr) ? box : free_box;
  nl.compute(R, static_cast<std::size_t>(N), box_use, opt);

  for (std::size_t p = 0; p < nl.size(); ++p) {
    const int i = static_cast<int>(nl.i(p));
    const int j = static_cast<int>(nl.j(p));
    if (i == j) {
      continue;
    }
    const double diffR = nl.distance(p);
    if (diffR <= 0.0) {
      continue;
    }
    const double *v = nl.vector(p); // r_j - r_i
    const double diffRX = -v[0];
    const double diffRY = -v[1];
    const double diffRZ = -v[2];

    const double a = pow(psi / diffR, 6);
    const double b = 4 * u0 * a;
    *U = *U + b * (a - 1) - cuttOffU;

    const double dU = -6 * b / diffR * (2 * a - 1);
    F[3 * i] = F[3 * i] - dU * diffRX / diffR;
    F[3 * i + 1] = F[3 * i + 1] - dU * diffRY / diffR;
    F[3 * i + 2] = F[3 * i + 2] - dU * diffRZ / diffR;
    F[3 * j] = F[3 * j] + dU * diffRX / diffR;
    F[3 * j + 1] = F[3 * j + 1] + dU * diffRY / diffR;
    F[3 * j + 2] = F[3 * j + 2] + dU * diffRZ / diffR;
  }
}

LJCluster::~LJCluster() { cleanMemory(); }
