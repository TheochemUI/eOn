/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#pragma once

#include "CyclicBlockBath.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace eonc::tunneling {

/// V(q,b)=V_fixed(q)+source(q)^T b+b^T E(q)b/2.
/// The bath represents physical coordinates, not an Ohmic ring thermostat.
struct PhysicalQuadraticBath {
  Eigen::MatrixXd coupling;   ///< retained x bath, source=C^T q by default
  Eigen::VectorXd curvatures; ///< strictly positive physical bath eigenvalues
  /// With V_relaxed=V_fixed-source^T E^-1 source/2, restore its counterterm.
  bool relaxedPotential = false;
  /// Source, bath x retained Jacobian, one retained Hessian per bath mode.
  using Source =
      std::function<void(const Eigen::VectorXd &, Eigen::VectorXd &,
                         Eigen::MatrixXd &, std::vector<Eigen::MatrixXd> &)>;
  Source source;

  struct CurvatureDerivatives {
    Eigen::MatrixXd value;
    std::vector<Eigen::MatrixXd> first;
    std::vector<std::vector<Eigen::MatrixXd>> second;
  };
  /// Full bath curvature in a fixed basis, with derivatives through order 0--2.
  using Curvature =
      std::function<CurvatureDerivatives(const Eigen::VectorXd &, int)>;
  Curvature curvature;

  struct Evaluation {
    double energy = 0.0;
    std::vector<Eigen::VectorXd> gradient;
    Eigen::MatrixXd hessian;
    std::vector<Eigen::VectorXd> bathPositions;
    std::vector<Eigen::VectorXd> inverseBath;
    /// Mixed bath/retained Hessians, including E_i*y for variable curvature.
    std::vector<Eigen::MatrixXd> sourceJacobians;
    std::shared_ptr<const CyclicBlockBath> blockBath;
    double bathLogDeterminant = 0.0;

    Eigen::MatrixXd solveBath(const Eigen::MatrixXd &rhs) const {
      if (blockBath)
        return blockBath->solve(rhs);
      if (inverseBath.empty())
        throw std::logic_error("bath inverse requires a Hessian evaluation");
      const long n = inverseBath.size(), b = inverseBath.front().size();
      if (rhs.rows() != n * b || !rhs.allFinite())
        throw std::invalid_argument("physical bath right hand side mismatch");
      Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n * b, rhs.cols());
      for (long j = 0; j < n; ++j)
        for (long l = 0; l < n; ++l)
          result.middleRows(j * b, b).noalias() +=
              inverseBath[(j - l + n) % n].asDiagonal() *
              rhs.middleRows(l * b, b);
      return result;
    }
  };

  void validate(long dof) const {
    if (dof < 1 || coupling.rows() != dof ||
        coupling.cols() != curvatures.size() || !coupling.allFinite() ||
        !curvatures.allFinite() || (curvatures.array() <= 0).any())
      throw std::invalid_argument(
          "physical bath requires finite coupling and positive bath curvature");
    if (curvature && relaxedPotential)
      throw std::invalid_argument(
          "variable physical bath requires a frozen potential callback");
  }

  CurvatureDerivatives curvatureAt(const Eigen::VectorXd &q, int order) const {
    const long r = q.size(), b = curvatures.size();
    if (!curvature || order < 0 || order > 2 || b < 1)
      throw std::invalid_argument("invalid physical bath curvature request");
    auto out = curvature(q, order);
    auto valid = [&](const Eigen::MatrixXd &matrix) {
      return matrix.rows() == b && matrix.cols() == b && matrix.allFinite() &&
             matrix.isApprox(matrix.transpose(), 1e-12);
    };
    if (!valid(out.value) ||
        (order >= 1 && out.first.size() != static_cast<size_t>(r)) ||
        (order >= 2 && out.second.size() != static_cast<size_t>(r)))
      throw std::invalid_argument(
          "physical bath curvature derivative mismatch");
    if (order >= 1)
      for (const auto &matrix : out.first)
        if (!valid(matrix))
          throw std::invalid_argument(
              "invalid physical bath curvature gradient");
    if (order >= 2) {
      for (long i = 0; i < r; ++i) {
        if (out.second[i].size() != static_cast<size_t>(r))
          throw std::invalid_argument(
              "physical bath curvature Hessian mismatch");
        for (long j = 0; j < r; ++j)
          if (!valid(out.second[i][j]))
            throw std::invalid_argument(
                "invalid physical bath curvature Hessian");
      }
      for (long i = 0; i < r; ++i)
        for (long j = 0; j < i; ++j)
          if (!out.second[i][j].isApprox(out.second[j][i], 1e-12))
            throw std::invalid_argument(
                "asymmetric physical bath curvature Hessian");
    }
    return out;
  }

  /// Physical curvature at the conditional bath equilibrium. The negative
  /// eigenvalue of this full operator sets the dynamical crossover frequency.
  Eigen::MatrixXd fullHessian(const Eigen::VectorXd &q,
                              const Eigen::MatrixXd &potentialHessian) const {
    const long r = q.size(), b = curvatures.size();
    validate(r);
    if (potentialHessian.rows() != r || potentialHessian.cols() != r ||
        !potentialHessian.allFinite())
      throw std::invalid_argument("physical bath potential Hessian mismatch");
    Eigen::VectorXd s;
    Eigen::MatrixXd jac;
    std::vector<Eigen::MatrixXd> second;
    if (source)
      source(q, s, jac, second);
    else {
      s = coupling.transpose() * q;
      jac = coupling.transpose();
      second.assign(b, Eigen::MatrixXd::Zero(r, r));
    }
    if (s.size() != b || jac.rows() != b || jac.cols() != r || !s.allFinite() ||
        !jac.allFinite() || second.size() != static_cast<size_t>(b))
      throw std::invalid_argument("physical bath source derivative mismatch");
    for (const auto &h : second)
      if (h.rows() != r || h.cols() != r || !h.allFinite() ||
          !h.isApprox(h.transpose(), 1e-12))
        throw std::invalid_argument("physical bath source Hessian mismatch");
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(r + b, r + b);
    out.topLeftCorner(r, r) = potentialHessian;
    if (curvature) {
      const auto state = curvatureAt(q, 2);
      const CyclicBlockBath factor({state.value}, 1.0);
      const Eigen::VectorXd y = -factor.solve(s).col(0);
      Eigen::MatrixXd mixed = jac;
      for (long a = 0; a < b; ++a)
        out.topLeftCorner(r, r).noalias() += y[a] * second[a];
      for (long i = 0; i < r; ++i) {
        mixed.col(i).noalias() += state.first[i] * y;
        for (long j = 0; j < r; ++j)
          out(i, j) += 0.5 * y.dot(state.second[i][j] * y);
      }
      out.topRightCorner(r, b) = mixed.transpose();
      out.bottomLeftCorner(b, r) = mixed;
      out.bottomRightCorner(b, b) = state.value;
      if (!out.allFinite() || !out.isApprox(out.transpose(), 1e-12))
        throw std::invalid_argument(
            "physical bath full Hessian is not finite and symmetric");
      return out;
    }
    if (relaxedPotential)
      out.topLeftCorner(r, r).noalias() +=
          jac.transpose() * curvatures.cwiseInverse().asDiagonal() * jac;
    else
      for (long k = 0; k < b; ++k)
        out.topLeftCorner(r, r).noalias() -= (s[k] / curvatures[k]) * second[k];
    out.topRightCorner(r, b) = jac.transpose();
    out.bottomLeftCorner(b, r) = jac;
    out.bottomRightCorner(b, b) = curvatures.asDiagonal();
    if (!out.allFinite() || !out.isApprox(out.transpose(), 1e-12))
      throw std::invalid_argument(
          "physical bath full Hessian is not finite and symmetric");
    return out;
  }

  /// Exact finite primitive-ring Gaussian elimination, with c=(beta_N*hbar)^-2.
  Evaluation evaluate(const std::vector<Eigen::VectorXd> &q, double c,
                      bool needHessian = false) const {
    if (q.empty() || !(c > 0) || !std::isfinite(c))
      throw std::invalid_argument(
          "physical bath requires a nonempty ring and positive spring");
    const long n = q.size(), r = q.front().size(), b = curvatures.size();
    validate(r);
    Evaluation out;
    out.gradient.assign(n, Eigen::VectorXd::Zero(r));
    out.bathPositions.assign(n, Eigen::VectorXd::Zero(b));
    if (needHessian)
      out.hessian = Eigen::MatrixXd::Zero(n * r, n * r);
    std::vector<Eigen::VectorXd> s(n);
    std::vector<Eigen::MatrixXd> jac(n);
    std::vector<std::vector<Eigen::MatrixXd>> second(n);
    for (long j = 0; j < n; ++j) {
      if (q[j].size() != r || !q[j].allFinite())
        throw std::invalid_argument("physical bath ring coordinate mismatch");
      if (source)
        source(q[j], s[j], jac[j], second[j]);
      else {
        s[j] = coupling.transpose() * q[j];
        jac[j] = coupling.transpose();
        second[j].assign(b, Eigen::MatrixXd::Zero(r, r));
      }
      if (s[j].size() != b || jac[j].rows() != b || jac[j].cols() != r ||
          !s[j].allFinite() || !jac[j].allFinite() ||
          second[j].size() != static_cast<size_t>(b))
        throw std::invalid_argument("physical bath source derivative mismatch");
      for (const auto &h : second[j])
        if (h.rows() != r || h.cols() != r || !h.allFinite() ||
            !h.isApprox(h.transpose(), 1e-12))
          throw std::invalid_argument("physical bath source Hessian mismatch");
    }
    if (curvature) {
      std::vector<Eigen::MatrixXd> blocks;
      blocks.reserve(n);
      Eigen::MatrixXd flatSource(n * b, 1);
      for (long j = 0; j < n; ++j) {
        blocks.push_back(curvatureAt(q[j], 0).value);
        flatSource.middleRows(j * b, b) = s[j];
      }
      auto factor = std::make_shared<CyclicBlockBath>(blocks, c);
      const Eigen::MatrixXd flatY = -factor->solve(flatSource);
      out.bathLogDeterminant = factor->logDeterminant();
      if (needHessian) {
        out.blockBath = factor;
        out.sourceJacobians.resize(n);
      }
      for (long j = 0; j < n; ++j) {
        const Eigen::VectorXd y = flatY.middleRows(j * b, b).col(0);
        out.bathPositions[j] = y;
        out.energy += 0.5 * s[j].dot(y);
        out.gradient[j].noalias() += jac[j].transpose() * y;
        const auto state = curvatureAt(q[j], needHessian ? 2 : 1);
        Eigen::MatrixXd mixed = jac[j];
        for (long i = 0; i < r; ++i) {
          out.gradient[j][i] += 0.5 * y.dot(state.first[i] * y);
          if (needHessian) {
            mixed.col(i).noalias() += state.first[i] * y;
            for (long k = 0; k < r; ++k)
              out.hessian(j * r + i, j * r + k) +=
                  0.5 * y.dot(state.second[i][k] * y);
          }
        }
        if (needHessian) {
          for (long a = 0; a < b; ++a)
            out.hessian.block(j * r, j * r, r, r).noalias() +=
                y[a] * second[j][a];
          out.sourceJacobians[j] = std::move(mixed);
        }
      }
      if (needHessian) {
        for (long l = 0; l < n; ++l) {
          Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(n * b, r);
          rhs.middleRows(l * b, b) = out.sourceJacobians[l];
          const Eigen::MatrixXd response = factor->solve(rhs);
          for (long j = 0; j < n; ++j)
            out.hessian.block(j * r, l * r, r, r).noalias() -=
                out.sourceJacobians[j].transpose() *
                response.middleRows(j * b, b);
        }
      }
      return out;
    }
    // inverse[lag] is diagonal in the fixed physical bath eigenbasis.
    std::vector<Eigen::VectorXd> inverse(n, Eigen::VectorXd::Zero(b));
    // The positive cyclic Green function solves (E+c*L)G=I. Its two
    // decaying branches sum periodic images without Fourier cancellation.
    const double rootSpring = std::sqrt(c);
    for (long a = 0; a < b; ++a) {
      const double frequency = std::sqrt(curvatures[a]);
      double alpha = 2 * std::asinh(frequency / (2 * rootSpring));
      if (!std::isfinite(alpha))
        alpha = std::log(curvatures[a]) - std::log(c);
      const double periodic = -std::expm1(-n * alpha);
      const double scale =
          1 / frequency / std::hypot(frequency, 2 * rootSpring) / periodic;
      out.bathLogDeterminant +=
          n * (std::log(c) + alpha) + 2 * std::log(periodic);
      for (long lag = 0; lag < n; ++lag)
        inverse[lag][a] =
            scale * (std::exp(-lag * alpha) + std::exp(-(n - lag) * alpha));
    }
    for (long j = 0; j < n; ++j) {
      Eigen::VectorXd y = Eigen::VectorXd::Zero(b);
      for (long l = 0; l < n; ++l)
        y.array() += inverse[(j - l + n) % n].array() * s[l].array();
      out.bathPositions[j] = -y;
      if (relaxedPotential)
        y.array() -= s[j].array() / curvatures.array();
      out.energy -= 0.5 * s[j].dot(y);
      out.gradient[j].noalias() -= jac[j].transpose() * y;
      if (needHessian) {
        for (long l = 0; l < n; ++l) {
          Eigen::VectorXd weights = inverse[(j - l + n) % n];
          if (relaxedPotential && j == l)
            weights.array() -= curvatures.cwiseInverse().array();
          out.hessian.block(j * r, l * r, r, r).noalias() -=
              jac[j].transpose() * weights.asDiagonal() * jac[l];
        }
        for (long a = 0; a < b; ++a)
          out.hessian.block(j * r, j * r, r, r).noalias() -=
              y[a] * second[j][a];
      }
    }
    if (needHessian) {
      out.inverseBath = inverse;
      out.sourceJacobians = jac;
    }
    return out;
  }
};
} // namespace eonc::tunneling
