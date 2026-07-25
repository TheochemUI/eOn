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
  eonc::CachedPairList::Options opt;
  opt.cutoff = (cuttOffR > 0.0) ? cuttOffR : 1.0e6;
  opt.periodic = {{false, false, false}};
  // Dummy orthorhombic box (ignored when non-periodic)
  double free_box[9] = {1e6, 0, 0, 0, 1e6, 0, 0, 0, 1e6};
  const double *box_use = (box != nullptr) ? box : free_box;
  const double psi2 = psi * psi;
  double energyAcc = 0.0;
  eonc::PairListCache::global().ensureVisit(
      R, static_cast<std::size_t>(N), box_use, opt,
      [&](int32_t i, int32_t j, double dx, double dy, double dz, double r2) {
    const double invR2 = 1.0 / r2;
    const double sr2 = psi2 * invR2;
    const double a = sr2 * sr2 * sr2; // (psi/r)^6 without pow()
    const double b = 4 * u0 * a;
    energyAcc += b * (a - 1) - cuttOffU;

    // -dU/dr / r along d = r_i - r_j, same sign convention as the historical
    // double loop.
    const double fscale = 6 * b * invR2 * (2 * a - 1);
    const double fx = fscale * dx;
    const double fy = fscale * dy;
    const double fz = fscale * dz;

    F[3 * i] += fx;
    F[3 * i + 1] += fy;
    F[3 * i + 2] += fz;
    F[3 * j] -= fx;
    F[3 * j + 1] -= fy;
    F[3 * j + 2] -= fz;
  });
  *U = energyAcc;
}

LJCluster::~LJCluster() { cleanMemory(); }
