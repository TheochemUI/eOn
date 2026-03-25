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

#include "LJ.h"
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

  // Pre-compute reciprocals for PBC and cutoff
  const double invBox0 = 1.0 / box[0];
  const double invBox4 = 1.0 / box[4];
  const double invBox8 = 1.0 / box[8];
  const double cutoff2 = cuttOffR * cuttOffR;
  const double psi2 = psi * psi;

  for (long i = 0; i < N - 1; i++) {
    const double xi = R[3 * i];
    const double yi = R[3 * i + 1];
    const double zi = R[3 * i + 2];

    for (long j = i + 1; j < N; j++) {
      double dx = xi - R[3 * j];
      double dy = yi - R[3 * j + 1];
      double dz = zi - R[3 * j + 2];

      // Minimum image convention (hoisted reciprocals)
      dx -= box[0] * std::floor(dx * invBox0 + 0.5);
      dy -= box[4] * std::floor(dy * invBox4 + 0.5);
      dz -= box[8] * std::floor(dz * invBox8 + 0.5);

      double r2 = dx * dx + dy * dy + dz * dz;

      // r^2 cutoff avoids sqrt for pairs outside cutoff
      if (r2 < cutoff2) {
        // Compute (sigma/r)^6 without pow(): sr2^3
        double invR2 = 1.0 / r2;
        double sr2 = psi2 * invR2;
        double sr6 = sr2 * sr2 * sr2;
        double e = 4.0 * u0 * sr6;

        *U += e * (sr6 - 1.0) - cuttOffU;

        // Force: -dU/dr * (1/r) * dr_vec = fscale * dr_vec
        double fscale = 6.0 * e * invR2 * (2.0 * sr6 - 1.0);
        double fx = fscale * dx;
        double fy = fscale * dy;
        double fz = fscale * dz;

        F[3 * i] += fx;
        F[3 * i + 1] += fy;
        F[3 * i + 2] += fz;
        F[3 * j] -= fx;
        F[3 * j + 1] -= fy;
        F[3 * j + 2] -= fz;
      }
    }
  }
}
