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

// Correctly rounded dot products and norms.
//
// Every finite double is an integer times a power of two, so x_i * y_i is
// the 106-bit integer product of the two significands times 2^(e_x + e_y),
// with no rounding. The products are added into a fixed-point integer
// accumulator that spans the whole exponent range (Kulisch), where addition
// is associative, and the exact total is rounded once to nearest, ties to
// even. The result is therefore a function of the multiset of products
// alone: the same bits for any summation order, vector width, contraction
// setting, compiler or Eigen release. This is the integer formulation of
// the error-free splitting in Ozaki, Ogita, Oishi and Rump, Numer.
// Algorithms 59 (2012), doi:10.1007/s11075-011-9478-1, and of Ozaki,
// Uchino and Imamura's scheme II (arXiv:2504.08009), taken to a single
// correctly rounded scalar.
//
// Use these where a reduction decides something (step lengths, curvature
// estimates, convergence tests, acceptance); Eigen's own reductions change
// their association between releases (3.4 sums two packet accumulators,
// 5.0 four) and with the SIMD width.

#include <cmath>
#include <cstddef>

#include <Eigen/Core>

namespace eonc::exact {

/// Round-to-nearest-even value of the exact sum of x[i] * y[i]. An exact
/// zero returns +0. A non-finite input returns the plain IEEE sum, so NaN
/// and infinity propagate as they would without this function.
[[nodiscard]] double dot(const double *x, const double *y, std::size_t n);

/// Round-to-nearest-even value of the exact sum of x[i].
[[nodiscard]] double sum(const double *x, std::size_t n);

/// A fixed summation tree: four lanes (term i goes to lane i mod 4, each
/// lane a left-to-right sum of rounded products), then (l0 + l1) + (l2 +
/// l3). Not exact, but every operation is a single IEEE multiply or add in
/// an order this translation unit fixes (it is compiled without
/// contraction), so the bits match on every compiler, Eigen release and
/// CPU. For bulk products that set a vector, not a branch.
[[nodiscard]] double orderedDot(const double *x, const double *y,
                                std::size_t n);

// Eigen front ends. A plain object or a contiguous block is read in place;
// any other expression is evaluated into a plain temporary first. The
// result does not depend on the storage order or the stride.
template <typename A, typename B>
[[nodiscard]] double dot(const Eigen::MatrixBase<A> &x,
                         const Eigen::MatrixBase<B> &y) {
  eigen_assert(x.size() == y.size());
  const auto &xe = x.derived().eval();
  const auto &ye = y.derived().eval();
  return dot(xe.data(), ye.data(), static_cast<std::size_t>(xe.size()));
}

template <typename A, typename B>
[[nodiscard]] double orderedDot(const Eigen::MatrixBase<A> &x,
                                const Eigen::MatrixBase<B> &y) {
  eigen_assert(x.size() == y.size());
  const auto &xe = x.derived().eval();
  const auto &ye = y.derived().eval();
  return orderedDot(xe.data(), ye.data(), static_cast<std::size_t>(xe.size()));
}

template <typename A>
[[nodiscard]] double squaredNorm(const Eigen::MatrixBase<A> &x) {
  const auto &xe = x.derived().eval();
  return dot(xe.data(), xe.data(), static_cast<std::size_t>(xe.size()));
}

/// sqrt of the correctly rounded squared norm; IEEE sqrt is itself
/// correctly rounded, so this is reproducible as well.
template <typename A>
[[nodiscard]] double norm(const Eigen::MatrixBase<A> &x) {
  return std::sqrt(squaredNorm(x));
}

/// x / norm(x), or zero below min_norm (the safemath convention).
template <typename A>
[[nodiscard]] typename A::PlainObject
normalized(const Eigen::MatrixBase<A> &x, double min_norm = 1e-300) {
  const double n = norm(x);
  if (n < min_norm) [[unlikely]] {
    return A::PlainObject::Zero(x.rows(), x.cols());
  }
  return (x / n).eval();
}

} // namespace eonc::exact
