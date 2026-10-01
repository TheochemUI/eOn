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

#include "eon/Eigen.h"
#include "eon/Tunneling.h"

#include <cmath>
#include <numbers>
#include <vector>

namespace eonc::testing {

using tunneling::BatchPotential;
using tunneling::kHbar;

// Symmetric Eckart barrier V = V0 / cosh^2(x / a) at unit mass, the H + H2
// like barrier every instanton paper tests on: V0 = 0.425 eV, a = 0.734
// amu^0.5 Angstrom, T_c = 150 K.
struct Eckart {
  double v0 = 0.425, a = 0.734;
  double value(double x) const { return v0 / std::pow(std::cosh(x / a), 2); }
  double slope(double x) const {
    return -2.0 * v0 * std::tanh(x / a) / (a * std::pow(std::cosh(x / a), 2));
  }
  double curvature(double x) const {
    const double s = std::sinh(x / a), c = std::cosh(x / a);
    return 2.0 * v0 * (2.0 * s * s - 1.0) / (a * a * std::pow(c, 4));
  }
  MatrixXd hessian_at(const VectorXd &q) const {
    return MatrixXd::Constant(1, 1, curvature(q(0)));
  }
  MatrixXd hessian_at_top() const {
    return MatrixXd::Constant(1, 1, -2.0 * v0 / (a * a));
  }
  BatchPotential batch() const {
    return [this](const std::vector<VectorXd> &q, std::vector<double> &v,
                  std::vector<VectorXd> &grad) {
      v.resize(q.size());
      grad.resize(q.size());
      for (size_t j = 0; j < q.size(); ++j) {
        v[j] = value(q[j](0));
        grad[j] = VectorXd::Constant(1, slope(q[j](0)));
      }
    };
  }
  // Exact transmission probability (Eckart 1930).
  double transmission(double e) const {
    const double alpha = a * std::sqrt(2.0 * e) / kHbar;
    const double d2 = 8.0 * v0 * a * a / (kHbar * kHbar) - 1.0;
    const double ch = d2 > 0.0 ? std::cosh(std::numbers::pi * std::sqrt(d2))
                               : std::cos(std::numbers::pi * std::sqrt(-d2));
    const double ca = std::cosh(2.0 * std::numbers::pi * alpha);
    return (ca - 1.0) / (ca + ch);
  }
  // ln of the exact thermal flux (1 / 2 pi hbar) int P(E) e^{-beta E} dE,
  // k Z_r for a free reactant per unit mass-weighted length.
  double logExactFlux(double beta) const {
    const double eMax = v0 + 60.0 / beta;
    const int n = 200000;
    double sum = 0.0;
    for (int k = 0; k <= n; ++k) {
      const double e = eMax * k / n + 1e-9;
      const double w = (k == 0 || k == n) ? 0.5 : 1.0;
      sum += w * transmission(e) * std::exp(-beta * e);
    }
    return std::log(sum * eMax / n / (2.0 * std::numbers::pi * kHbar));
  }
};

} // namespace eonc::testing
