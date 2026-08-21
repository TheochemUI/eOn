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
#pragma once

#include "BaseStructures.h"
#include "Eigen.h"
#include "ObjectiveFunction.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace eonc {
namespace pairhess {

inline void addPairBlock(Eigen::MatrixXd &P, int i, int j,
                         const Eigen::Vector3d &rij, double k_par,
                         double k_perp) {
  const double r = rij.norm();
  if (r < 1.0e-14 || (!std::isfinite(k_par) && !std::isfinite(k_perp))) {
    return;
  }
  if (k_par == 0.0 && k_perp == 0.0) {
    return;
  }
  const Eigen::Vector3d u = rij / r;
  const Eigen::Matrix3d H = k_perp * Eigen::Matrix3d::Identity() +
                            (k_par - k_perp) * (u * u.transpose());
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      const int ia = 3 * i + a;
      const int ib = 3 * i + b;
      const int ja = 3 * j + a;
      const int jb = 3 * j + b;
      P(ia, ib) += H(a, b);
      P(ja, jb) += H(a, b);
      P(ia, jb) -= H(a, b);
      P(ja, ib) -= H(a, b);
    }
  }
}

inline void ljDerivs(double r, double eps, double sigma, double &vp,
                     double &vpp) {
  const double inv = sigma / r;
  const double x6 = inv * inv * inv * inv * inv * inv;
  const double x12 = x6 * x6;
  vp = 4.0 * eps * (-12.0 * x12 + 6.0 * x6) / r;
  vpp = 4.0 * eps * (156.0 * x12 - 42.0 * x6) / (r * r);
}

inline void morseDerivs(double r, double D, double a, double re, double &vp,
                        double &vpp) {
  const double e1 = std::exp(-a * (r - re));
  const double e2 = e1 * e1;
  vp = 2.0 * D * a * (e1 - e2);
  vpp = 2.0 * D * a * a * (2.0 * e2 - e1);
}

/// Analytic pair Hessian. `kind` is pair / pair_abs / pair_full / lindh / exp / c1.
inline Eigen::MatrixXd
build(const Eigen::VectorXd &pos, const std::string &kind, PotType pot,
      double A, double mu, double rcut_in, ObjectiveFunction &objf) {
  const int n = static_cast<int>(pos.size());
  const int nat = n / 3;
  std::vector<double> nn(static_cast<size_t>(nat),
                         std::numeric_limits<double>::infinity());
  auto mic = [&](int i, int j) {
    Eigen::Vector3d dr = pos.segment<3>(3 * i) - pos.segment<3>(3 * j);
    objf.minimumImage(dr);
    return dr;
  };
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < nat; ++j) {
      if (i == j) {
        continue;
      }
      const double r = mic(i, j).norm();
      if (r > 1.0e-12 && r < nn[static_cast<size_t>(i)]) {
        nn[static_cast<size_t>(i)] = r;
      }
    }
  }
  double r_nn = 0.0;
  for (double v : nn) {
    if (std::isfinite(v)) {
      r_nn = std::max(r_nn, v);
    }
  }
  if (r_nn <= 0.0) {
    r_nn = 1.0;
  }
  const double rcut = rcut_in > 0.0 ? rcut_in : 2.5 * r_nn;
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < nat; ++i) {
    for (int j = i + 1; j < nat; ++j) {
      const Eigen::Vector3d rij = mic(i, j);
      const double r = rij.norm();
      if (r >= rcut || r < 1.0e-14) {
        continue;
      }
      if (kind == "pair" || kind == "pair_abs" || kind == "pair_full") {
        double vp = 0.0;
        double vpp = 0.0;
        if (pot == PotType::MORSE_PT) {
          morseDerivs(r, 0.7102, 1.6047, 2.897, vp, vpp);
        } else {
          ljDerivs(r, 1.0, 1.0, vp, vpp);
        }
        double k_par = std::max(vpp, 0.0);
        double k_perp = std::max(vp / r, 0.0);
        if (kind == "pair_abs") {
          k_par = std::abs(vpp);
          k_perp = std::abs(vp / r);
        } else if (kind == "pair_full") {
          k_par = vpp;
          k_perp = vp / r;
        }
        addPairBlock(P, i, j, rij, k_par, k_perp);
      } else if (kind == "lindh") {
        const double alpha = A > 0.0 ? A : 1.0;
        const double k = mu * std::exp(alpha * (r_nn * r_nn - r * r));
        addPairBlock(P, i, j, rij, k, 0.0);
      } else if (kind == "c1") {
        const double x = r / rcut;
        const double w = mu * (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
        addPairBlock(P, i, j, rij, w, w);
      } else {
        const double w = mu * std::exp(-A * (r / r_nn - 1.0));
        addPairBlock(P, i, j, rij, w, w);
      }
    }
  }
  const double shift = 1.0e-6 * std::max(1.0, P.diagonal().mean());
  P.diagonal().array() += shift;
  return P;
}

} // namespace pairhess
} // namespace eonc
