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

#include "eon/potentials/LJ/LJ.h"
#include "eon/VesinNeighbors.h"
#include <cmath>

void LJ::setParameters(double u0In, double cutoffIn, double psiIn) {
  u0 = u0In;
  psi = psiIn;
  cuttOffR = cutoffIn;
  double sr = psi / cuttOffR;
  double sr2 = sr * sr;
  double r6 = sr2 * sr2 * sr2;
  cuttOffU = 4.0 * u0 * r6 * (r6 - 1.0);
}

void LJ::force(long N, const double *R, const int * /*atomicNrs*/, double *F,
               double *U, double *variance, const double *box) {
  variance = nullptr;
  *U = 0.0;
  for (long i = 0; i < 3 * N; i++) {
    F[i] = 0.0;
  }
  if (N < 2 || cuttOffR <= 0.0) {
    return;
  }

  eonc::CachedPairList::Options opt;
  opt.cutoff = cuttOffR;
  const auto &nl = eonc::PairListCache::local().ensure(
      R, static_cast<std::size_t>(N), box, opt);

  const double psi2 = psi * psi;
  double energyAcc = 0.0;
  nl.forEach(R, [&](int32_t i, int32_t j, double dx, double dy, double dz,
                    double r2) {
    const double invR2 = 1.0 / r2;
    const double sr2 = psi2 * invR2;
    const double sr6 = sr2 * sr2 * sr2;
    const double e = 4.0 * u0 * sr6;
    energyAcc += e * (sr6 - 1.0) - cuttOffU;

    const double fscale = 6.0 * e * invR2 * (2.0 * sr6 - 1.0);
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
