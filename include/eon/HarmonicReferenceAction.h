/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace eonc::tunneling {

/// Exact nonnegative harmonic-reference links in mass-weighted coordinates.
/// tau = beta*hbar/N. Add correction(q) to the sum of physical bead energies;
/// the physical energy zero is retained. The residual-potential splitting is
/// second order, and is exact only when the residual commutes with the reference.
/// This action does not define a real-time transmission factor or rate prefactor.
class HarmonicReferenceAction {
public:
  HarmonicReferenceAction(const Eigen::MatrixXd &h0,
                          const Eigen::VectorXd &center, double tau)
      : h0_(h0), center_(center) {
    if (!(tau > 0.0) || !std::isfinite(tau) || center.size() == 0 ||
        !center.allFinite() || h0.rows() != center.size() ||
        h0.cols() != center.size() || !h0.allFinite() ||
        !h0.isApprox(h0.transpose(), 32 * std::numeric_limits<double>::epsilon())) {
      throw std::invalid_argument("harmonic reference requires finite symmetric H0 and positive tau");
    }
    h0_ = 0.5*h0+0.5*h0.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(h0_);
    if (eig.info() != Eigen::Success || eig.eigenvalues().minCoeff() < 0.0) {
      throw std::invalid_argument("harmonic reference requires nonnegative H0");
    }
    Eigen::VectorXd b(center.size()), d(center.size());
    logMeasurePerBead_ = 0.0;
    for (Eigen::Index i = 0; i < center.size(); ++i) {
      const double omega = std::sqrt(eig.eigenvalues()[i]);
      const double z = tau * omega;
      if (!std::isfinite(z)) {
        throw std::invalid_argument("harmonic reference frequency-time product overflows");
      }
      double logRatio;
      if (z < 1e-4) {
        const double z2 = z*z;
        const double ratio = 1.0-z2/6.0+7.0*z2*z2/360.0;
        b[i] = ratio / tau / tau;
        logRatio = -z2/6.0+z2*z2/180.0;
      } else {
        logRatio = std::log(2.0)+std::log(z)-z-std::log(-std::expm1(-2*z));
        b[i] = std::exp(logRatio) / tau / tau;
      }
      d[i] = 2*omega*std::tanh(z/2)/tau;
      logMeasurePerBead_ += 0.5*logRatio;
    }
    b_ = eig.eigenvectors()*b.asDiagonal()*eig.eigenvectors().transpose();
    d_ = eig.eigenvectors()*d.asDiagonal()*eig.eigenvectors().transpose();
    if (!b_.allFinite() || !d_.allFinite() || !std::isfinite(logMeasurePerBead_)) {
      throw std::invalid_argument("harmonic reference coefficients overflow");
    }
  }

  const Eigen::MatrixXd &link() const { return b_; }
  const Eigen::MatrixXd &onsite() const { return d_; }
  double logMeasureCorrection(long beads) const {
    if (beads < 1) throw std::invalid_argument("ring requires at least one bead");
    return beads*logMeasurePerBead_;
  }

  /// Rows are beads; columns are ordered mass-weighted coordinates.
  double correction(const Eigen::MatrixXd &q) const {
    validate(q);
    double value = 0.0;
    for (Eigen::Index j = 0; j < q.rows(); ++j) {
      const Eigen::VectorXd x = q.row(j).transpose()-center_;
      const Eigen::VectorXd dx = q.row((j+1)%q.rows()).transpose()-q.row(j).transpose();
      value += 0.5*x.dot((d_-h0_)*x)+0.5*dx.dot(b_*dx);
    }
    return value;
  }

  Eigen::MatrixXd gradientCorrection(const Eigen::MatrixXd &q) const {
    validate(q);
    Eigen::MatrixXd result(q.rows(), q.cols());
    for (Eigen::Index j = 0; j < q.rows(); ++j) {
      const Eigen::VectorXd x = q.row(j).transpose()-center_;
      const Eigen::VectorXd lap = 2*q.row(j).transpose()
          - q.row((j+1)%q.rows()).transpose()
          - q.row((j+q.rows()-1)%q.rows()).transpose();
      result.row(j) = ((d_-h0_)*x+b_*lap).transpose();
    }
    return result;
  }

  /// Flattening order is bead-major, including the N=1 and N=2 cycles.
  Eigen::MatrixXd hessianCorrection(long beads) const {
    if (beads < 1) throw std::invalid_argument("ring requires at least one bead");
    const Eigen::Index f = center_.size();
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(beads*f, beads*f);
    for (long j = 0; j < beads; ++j) {
      const long k = (j+1)%beads;
      result.block(j*f,j*f,f,f) += d_-h0_;
      result.block(j*f,j*f,f,f) += b_;
      result.block(k*f,k*f,f,f) += b_;
      result.block(j*f,k*f,f,f) -= b_;
      result.block(k*f,j*f,f,f) -= b_;
    }
    return result;
  }

private:
  void validate(const Eigen::MatrixXd &q) const {
    if (q.rows() < 1 || q.cols() != center_.size() || !q.allFinite()) {
      throw std::invalid_argument("ring coordinates have invalid shape or values");
    }
  }
  Eigen::MatrixXd h0_, b_, d_;
  Eigen::VectorXd center_;
  double logMeasurePerBead_;
};
} // namespace eonc::tunneling
