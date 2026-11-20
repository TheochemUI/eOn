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

#include "Eigen.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace eonc::tunneling {

/// out[j] = H(q_j) v[j]: the mass-weighted Hessian of every ring bead
/// applied to that bead's part of a ring vector, from exact products.
using RingHessianProducts = std::function<void(const std::vector<VectorXd> &v,
                                               std::vector<VectorXd> &out)>;

/// H v for the mass-weighted Hessian at one structure.
using HessianProduct = std::function<VectorXd(const VectorXd &v)>;

/// A log-determinant from stochastic Lanczos quadrature. The error terms
/// are kept apart: the standard error of the probe mean, the Gauss
/// quadrature tail averaged over the probes, and the error of the removed
/// near-zero eigenvalue.
struct StochasticLogDet {
  double value{0.0};
  double standardError{0.0};
  double quadratureError{0.0};
  double modeError{0.0};
  long probes{0};
  /// Krylov directions whose block enters exactly.
  long deflated{0};
  /// Operator applications, each one product per bead.
  long products{0};
  /// Negative eigenvalues, the removed near-zero one left out.
  long negativeModes{0};
  /// Lowest eigenvalue apart from the removed one.
  double lowestEigenvalue{0.0};
  /// The removed eigenvalue nearest the cycle direction.
  double zeroEigenvalue{0.0};
  /// tr(H) / dimension of the well operator.
  double meanCurvature{0.0};
  /// The Krylov space filled the whole space, so the value is exact.
  bool exhausted{false};
  [[nodiscard]] double error() const {
    return standardError + quadratureError + modeError;
  }
};

/// ln |det' J| for the ring Hessian J = blockdiag(H_j) + c L, L the cyclic
/// Laplacian over beads, with its rigid null space and the eigenvalue
/// along the cycle direction removed. Only products with J are used.
///
/// J is congruent to M^{-1/2} J M^{-1/2} with M = c L + shift I, which the
/// real discrete Fourier basis over beads diagonalises, so ln det M is a
/// closed sum. A two-column block Lanczos run seeded with the cycle
/// direction and one Gaussian direction resolves the near-zero and
/// negative directions, its block enters exactly, and the Schur complement
/// on the rest is estimated with Gaussian probes. The Krylov block doubles
/// while the probes it saves cost more than the block itself, while the
/// removed eigenvalue is not resolved to a quarter of `tolerance`, and
/// while a probe encounters a direction the block has not captured. Each
/// quadrature uses a tail-estimate target of a quarter of `tolerance`; probes
/// stop at an estimated standard error of half `tolerance`. The reported
/// sum of standard, quadrature-tail and mode estimates is diagnostic, not a
/// certified bound.
///
/// `cycle` has one vector per bead; `rigid` holds Euclidean orthonormal null
/// vectors of J, each one vector per bead. `shift` must be positive.
StochasticLogDet
ringLogDetPrime(const RingHessianProducts &products, double c, double shift,
                const std::vector<VectorXd> &cycle,
                const std::vector<std::vector<VectorXd>> &rigid,
                double tolerance, std::uint64_t seed);

/// sum_k ln det(H + shifts[k]) over the complement of the orthonormal
/// columns of `rigid`, from products with H; one Krylov space serves every
/// shift. A probe whose Krylov space fills the complement gives the
/// spectrum exactly, and the next probe is run to that end once the probes
/// still needed would cost more steps than the complement has directions.
/// Otherwise the probes stop at an estimated standard error of half
/// `tolerance`, each quadrature at an estimated tail of a quarter. These
/// terms are diagnostics, not certified bounds.
StochasticLogDet wellLogDet(const HessianProduct &product, long freedom,
                            const MatrixXd &rigid,
                            const std::vector<double> &shifts, double tolerance,
                            std::uint64_t seed);

} // namespace eonc::tunneling
