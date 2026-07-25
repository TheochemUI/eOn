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

  eonc::VesinNeighbors::Options opt;
  opt.cutoff = cuttOffR;
  opt.full = false; // half list: one pair entry per unordered bond
  opt.return_distances = true;
  opt.return_vectors = true;
  nl_.compute(R, static_cast<std::size_t>(N), box, opt);

  const double psi2 = psi * psi;
  for (std::size_t p = 0; p < nl_.size(); ++p) {
    const long i = static_cast<long>(nl_.i(p));
    const long j = static_cast<long>(nl_.j(p));
    if (i == j) {
      continue;
    }
    const double r = nl_.distance(p);
    if (r <= 0.0) {
      continue;
    }
    const double r2 = r * r;
    const double invR2 = 1.0 / r2;
    // vesin vector is r_j - r_i; LJ force uses r_i - r_j
    const double *v = nl_.vector(p);
    const double dx = -v[0];
    const double dy = -v[1];
    const double dz = -v[2];

    const double sr2 = psi2 * invR2;
    const double sr6 = sr2 * sr2 * sr2;
    const double e = 4.0 * u0 * sr6;
    *U += e * (sr6 - 1.0) - cuttOffU;

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
  }
}
