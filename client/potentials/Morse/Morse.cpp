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

#include "eon/potentials/Morse/Morse.h"
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
  if (N < 2 || cutoff_ <= 0.0) {
    return;
  }

  eonc::VesinNeighbors::Options opt;
  opt.cutoff = cutoff_;
  opt.full = false;
  opt.return_distances = true;
  opt.return_vectors = true;
  nl_.compute(R, static_cast<std::size_t>(N), box, opt);

  const double a = a_;
  const double re = re_;
  const double De = De_;
  const double twoDeA = 2.0 * De * a;
  const double eCut = energyCutoff_;
  double energyAcc = 0.0;

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
    // vesin: r_j - r_i; Morse force used r_i - r_j
    const double *v = nl_.vector(p);
    const double dx = -v[0];
    const double dy = -v[1];
    const double dz = -v[2];

    const double d = 1.0 - std::exp(-a * (r - re));
    const double energy = De * d * d - De;
    const double fmag = twoDeA * d * (d - 1.0);
    energyAcc += energy - eCut;

    const double fscale = fmag / r;
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
  *U = energyAcc;
}

void Morse::morse(double r, double &energy, double &force) const {
  double const d = 1.0 - std::exp(-a_ * (r - re_));
  energy = De_ * d * d - De_;
  force = 2.0 * De_ * d * (d - 1.0) * a_;
}
