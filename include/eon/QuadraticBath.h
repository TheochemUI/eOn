/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace eonc::tunneling {

/// Exact Gaussian elimination of a stable physical mass-weighted harmonic bath.
/// All nonpositive complementary curvatures are promoted to reaction
/// coordinates. The measured Hessian defines a local quadratic model; this
/// class does not assert that an anharmonic potential has a globally quadratic
/// bath.
class QuadraticBath {
public:
  QuadraticBath(const Eigen::MatrixXd &hessian,
                const Eigen::MatrixXd &reactionBasis) {
    const Eigen::Index n = hessian.rows();
    if (n < 1 || hessian.cols() != n || !hessian.allFinite() ||
        !hessian.isApprox(hessian.transpose(),
                          32 * std::numeric_limits<double>::epsilon()) ||
        reactionBasis.rows() != n || reactionBasis.cols() < 1 ||
        reactionBasis.cols() > n || !reactionBasis.allFinite() ||
        !(reactionBasis.transpose() * reactionBasis)
             .isApprox(Eigen::MatrixXd::Identity(reactionBasis.cols(),
                                                 reactionBasis.cols()),
                       1e-12)) {
      throw std::invalid_argument("quadratic bath requires symmetric finite H "
                                  "and orthonormal reaction basis");
    }
    const Eigen::MatrixXd h = 0.5 * hessian + 0.5 * hessian.transpose();
    Eigen::MatrixXd retained = reactionBasis;
    for (;;) {
      Eigen::HouseholderQR<Eigen::MatrixXd> qr(retained);
      const Eigen::MatrixXd orthogonal =
          qr.householderQ() * Eigen::MatrixXd::Identity(n, n);
      reaction_ = orthogonal.leftCols(retained.cols());
      bath_ = orthogonal.rightCols(n - retained.cols());
      if (bath_.cols() == 0) {
        frequenciesSquared_.resize(0);
        break;
      }
      const Eigen::MatrixXd e = bath_.transpose() * h * bath_;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(e);
      if (eig.info() != Eigen::Success || !eig.eigenvalues().allFinite()) {
        throw std::runtime_error("quadratic bath eigensolver failed");
      }
      Eigen::Index promote = 0;
      while (promote < eig.eigenvalues().size() &&
             eig.eigenvalues()[promote] <= 0.0)
        ++promote;
      if (promote == 0) {
        bath_ = bath_ * eig.eigenvectors();
        frequenciesSquared_ = eig.eigenvalues();
        break;
      }
      retained.resize(n, reaction_.cols() + promote);
      retained.leftCols(reaction_.cols()) = reaction_;
      retained.rightCols(promote) =
          bath_ * eig.eigenvectors().leftCols(promote);
    }
    a_ = reaction_.transpose() * h * reaction_;
    c_ = reaction_.transpose() * h * bath_;
  }

  const Eigen::MatrixXd &reactionBasis() const { return reaction_; }
  const Eigen::MatrixXd &bathBasis() const { return bath_; }
  const Eigen::MatrixXd &reactionCurvature() const { return a_; }
  const Eigen::MatrixXd &coupling() const { return c_; }
  const Eigen::VectorXd &bathCurvatures() const { return frequenciesSquared_; }

  /// Schur complement of H+nu^2 I, including reaction kinetic curvature.
  Eigen::MatrixXd kernel(double nu) const {
    const double nu2 = squaredFrequency(nu);
    Eigen::MatrixXd result = a_;
    result.diagonal().array() += nu2;
    for (Eigen::Index k = 0; k < frequenciesSquared_.size(); ++k)
      result.noalias() -=
          c_.col(k) * c_.col(k).transpose() / (frequenciesSquared_[k] + nu2);
    return result;
  }

  /// Laplace transform of the physical harmonic memory kernel for s >= 0.
  Eigen::MatrixXd memory(double s) const {
    const double s2 = squaredFrequency(s);
    if (s < 0)
      throw std::invalid_argument("Laplace frequency must be nonnegative");
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(a_.rows(), a_.cols());
    for (Eigen::Index k = 0; k < frequenciesSquared_.size(); ++k) {
      const double e = frequenciesSquared_[k];
      result.noalias() +=
          (s / e / (e + s2)) * c_.col(k) * c_.col(k).transpose();
    }
    return result;
  }

  double logBathDeterminant(double nu) const {
    const double nu2 = squaredFrequency(nu);
    return (frequenciesSquared_.array() + nu2).log().sum();
  }

  /// Normalized quantum partition of the uncoupled conditional harmonic bath.
  /// Ratios cancel only for identical conditional baths and normalization.
  double logBathPartition(double betaHbar) const {
    if (!(betaHbar > 0) || !std::isfinite(betaHbar))
      throw std::invalid_argument("beta*hbar must be positive and finite");
    double value = 0;
    for (Eigen::Index k = 0; k < frequenciesSquared_.size(); ++k) {
      const double z = betaHbar * std::sqrt(frequenciesSquared_[k]);
      if (!(z > 0) || !std::isfinite(z))
        throw std::invalid_argument(
            "thermal bath frequency overflows or underflows");
      value -= z / 2 + std::log(-std::expm1(-z));
    }
    return value;
  }

private:
  static double squaredFrequency(double nu) {
    if (!std::isfinite(nu) || !std::isfinite(nu * nu))
      throw std::invalid_argument("frequency must have a finite square");
    return nu * nu;
  }
  Eigen::MatrixXd reaction_, bath_, a_, c_;
  Eigen::VectorXd frequenciesSquared_;
};
} // namespace eonc::tunneling
