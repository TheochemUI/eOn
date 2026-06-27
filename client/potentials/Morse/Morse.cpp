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

#include "Morse.h"
#include <cassert>
#include <cmath>

/** @file
      @brief Morse potential for platinum
      @author Anonymous (possibly A. Pedersen or G. Henkelman), revision: Jean
   Claude C. Berthet
      @date Unknown, revision: 2010, University of Iceland
      */

void Morse::setParameters(double De, double a, double re, double cutoff) {
  re_ = re;
  De_ = De;
  a_ = a;
  cutoff_ = cutoff;
  double f;
  morse(cutoff, energyCutoff_, f);
}

void Morse::force(long N, const double *R, const int * /*atomicNrs*/, double *F,
                  double *U, double *variance, const double *box) {
  variance = nullptr;
  assert(box[0] > 0.0 && box[4] > 0.0 && box[8] > 0.0);

  *U = 0.0;
  for (long k = 0; k < 3 * N; k++) {
    F[k] = 0.0;
  }

  // Orthorhombic MIC (box diagonal). Precompute inverses and hot constants.
  const double invBox0 = 1.0 / box[0];
  const double invBox4 = 1.0 / box[4];
  const double invBox8 = 1.0 / box[8];
  const double cutoff2 = cutoff_ * cutoff_;
  const double a = a_;
  const double re = re_;
  const double De = De_;
  const double twoDeA = 2.0 * De * a;
  const double eCut = energyCutoff_;
  double energyAcc = 0.0;

  for (long i = 0; i < N; i++) {
    const double xi = R[3 * i];
    const double yi = R[3 * i + 1];
    const double zi = R[3 * i + 2];
    double *Fi = F + 3 * i;

    for (long j = i + 1; j < N; j++) {
      double dx = xi - R[3 * j];
      double dy = yi - R[3 * j + 1];
      double dz = zi - R[3 * j + 2];

      // Minimum image: floor(x+0.5) avoids fmod; branchless on modern CPUs.
      dx -= box[0] * std::floor(dx * invBox0 + 0.5);
      dy -= box[4] * std::floor(dy * invBox4 + 0.5);
      dz -= box[8] * std::floor(dz * invBox8 + 0.5);

      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 >= cutoff2) {
        continue;
      }

      const double r = std::sqrt(r2);
      // Inline morse() to keep the pair loop in one TU for inlining / LTO.
      const double d = 1.0 - std::exp(-a * (r - re));
      const double energy = De * d * d - De;
      const double fmag = twoDeA * d * (d - 1.0);

      energyAcc += energy - eCut;

      // Pair force magnitude along MIC vector (keep /r for bit-stability vs
      // main)
      const double fscale = fmag / r;
      const double fx = fscale * dx;
      const double fy = fscale * dy;
      const double fz = fscale * dz;

      Fi[0] += fx;
      Fi[1] += fy;
      Fi[2] += fz;
      F[3 * j] -= fx;
      F[3 * j + 1] -= fy;
      F[3 * j + 2] -= fz;
    }
  }
  *U = energyAcc;
}

void Morse::morse(double r, double &energy, double &force) const {
  double const d = 1.0 - std::exp(-a_ * (r - re_));
  energy = De_ * d * d - De_;
  force = 2.0 * De_ * d * (d - 1.0) * a_;
}
