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

#include <algorithm>
#include <cmath>

namespace eonc::safemath {

inline constexpr double eps = 1e-300;

inline constexpr double safe_div(double num, double denom,
                                 double fallback = 0.0) {
  return (std::abs(denom) < eps) ? fallback : num / denom;
}

inline constexpr double safe_recip(double x, double fallback = 0.0) {
  return safe_div(1.0, x, fallback);
}

inline double safe_acos(double x) {
  return std::acos(std::clamp(x, -1.0, 1.0));
}

inline double safe_sqrt(double x) { return std::sqrt(std::max(0.0, x)); }

inline double safe_atan_ratio(double num, double denom,
                              double fallback = 0.0) {
  return (std::abs(denom) < eps) ? fallback : std::atan(num / denom);
}

} // namespace eonc::safemath

// Eigen-dependent utilities, available only when Eigen is already included
#ifdef EIGEN_CORE_H
namespace eonc::safemath {

template <typename Derived>
auto safe_normalized(const Eigen::MatrixBase<Derived> &v,
                     double min_norm = eps) {
  using ResultType = typename Derived::PlainObject;
  double n = v.norm();
  if (n < min_norm) {
    return ResultType::Zero(v.rows(), v.cols()).eval();
  }
  return (v / n).eval();
}

} // namespace eonc::safemath
#endif
