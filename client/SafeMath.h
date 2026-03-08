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

[[nodiscard]] inline constexpr double safe_div(double num, double denom,
                                               double fallback = 0.0) {
  if (std::abs(denom) < eps) [[unlikely]] {
    return fallback;
  }
  return num / denom;
}

[[nodiscard]] inline constexpr double safe_recip(double x,
                                                 double fallback = 0.0) {
  return safe_div(1.0, x, fallback);
}

[[nodiscard]] inline double safe_acos(double x) {
  return std::acos(std::clamp(x, -1.0, 1.0));
}

[[nodiscard]] inline double safe_sqrt(double x) {
  return std::sqrt(std::max(0.0, x));
}

[[nodiscard]] inline double safe_atan_ratio(double num, double denom,
                                            double fallback = 0.0) {
  if (std::abs(denom) < eps) [[unlikely]] {
    return fallback;
  }
  return std::atan(num / denom);
}

} // namespace eonc::safemath

// Eigen-dependent utilities, available only when Eigen is already included.
// Include order: Eigen headers first, then SafeMath.h.
#ifdef EIGEN_CORE_H
namespace eonc::safemath {

template <typename Derived>
[[nodiscard]] auto safe_normalized(const Eigen::MatrixBase<Derived> &v,
                                   double min_norm = eps) {
  using ResultType = typename Derived::PlainObject;
  double n = v.norm();
  if (n < min_norm) [[unlikely]] {
    return ResultType::Zero(v.rows(), v.cols()).eval();
  }
  return (v / n).eval();
}

template <typename Derived>
void safe_normalize_inplace(Eigen::MatrixBase<Derived> &v,
                            double min_norm = eps) {
  double n = v.norm();
  if (n < min_norm) [[unlikely]] {
    v.derived().setZero();
  } else {
    v.derived() /= n;
  }
}

} // namespace eonc::safemath
#endif
