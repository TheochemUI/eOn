/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#include "eon/HarmonicReferenceAction.h"
#include "catch2/catch_amalgamated.hpp"

using eonc::tunneling::HarmonicReferenceAction;

TEST_CASE("Mehler reference has the exact harmonic trace at every bead count",
          "[instanton][harmonic-reference]") {
  const double beta = 3.2, hbar = 0.7, offset = -2.3;
  Eigen::Matrix2d h;
  h << 2.0, 0.4, 0.4, 0.8;
  const Eigen::Vector2d center(1.1, -0.3);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(h);
  const double exact =
      -beta * offset -
      (2 * (beta * hbar * eig.eigenvalues().array().sqrt() / 2).sinh())
          .log()
          .sum();
  for (long n : {1L, 2L, 3L, 8L, 17L}) {
    INFO(n);
    const double betaN = beta / n;
    HarmonicReferenceAction action(h, center, betaN * hbar);
    Eigen::MatrixXd k = action.hessianCorrection(n);
    for (long j = 0; j < n; ++j)
      k.block<2, 2>(2 * j, 2 * j) += h;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> modes(k);
    REQUIRE(modes.eigenvalues().minCoeff() > 0);
    const double logZ = -beta * offset + action.logMeasureCorrection(n) -
                        n * std::log(betaN * hbar * hbar) -
                        0.5 * modes.eigenvalues().array().log().sum();
    REQUIRE(std::abs(logZ - exact) < 2e-12);
  }
}

TEST_CASE("Mehler correction derivatives include all cyclic links",
          "[instanton][harmonic-reference]") {
  Eigen::Matrix2d h;
  h << 2.0, 0.4, 0.4, 0.8;
  HarmonicReferenceAction action(h, Eigen::Vector2d(0.3, -0.2), 0.4);
  for (long n : {1L, 2L, 5L}) {
    INFO(n);
    Eigen::MatrixXd q(n, 2);
    for (long j = 0; j < n; ++j)
      q.row(j) << 0.1 * j, std::cos(0.7 * j);
    const Eigen::MatrixXd g = action.gradientCorrection(q);
    const Eigen::MatrixXd k = action.hessianCorrection(n);
    REQUIRE(k.isApprox(k.transpose(), 1e-14));
    const double step = 1e-5;
    for (long j = 0; j < n; ++j)
      for (long a = 0; a < 2; ++a) {
        Eigen::MatrixXd plus = q, minus = q;
        plus(j, a) += step;
        minus(j, a) -= step;
        REQUIRE(std::abs((action.correction(plus) - action.correction(minus)) /
                             (2 * step) -
                         g(j, a)) < 2e-9);
        const Eigen::MatrixXd dg = (action.gradientCorrection(plus) -
                                    action.gradientCorrection(minus)) /
                                   (2 * step);
        for (long l = 0; l < n; ++l)
          for (long b = 0; b < 2; ++b)
            REQUIRE(std::abs(dg(l, b) - k(2 * l + b, 2 * j + a)) < 2e-9);
      }
  }
}

TEST_CASE("Mehler reference is covariant and retains the primitive limit",
          "[instanton][harmonic-reference]") {
  Eigen::Matrix2d h;
  h << 2.0, 0.4, 0.4, 0.8;
  Eigen::Matrix2d rotation;
  rotation << std::cos(0.4), -std::sin(0.4), std::sin(0.4), std::cos(0.4);
  Eigen::Vector2d center(0.1, -0.3);
  HarmonicReferenceAction original(h, center, 0.4);
  HarmonicReferenceAction rotated(rotation * h * rotation.transpose(),
                                  rotation * center, 0.4);
  Eigen::MatrixXd q(3, 2);
  q << 0.2, 0.7, 0.8, -0.1, -0.5, 0.4;
  REQUIRE(std::abs(original.correction(q) -
                   rotated.correction(q * rotation.transpose())) < 1e-12);
  REQUIRE((original.gradientCorrection(q) * rotation.transpose())
              .isApprox(rotated.gradientCorrection(q * rotation.transpose()),
                        1e-13));
  HarmonicReferenceAction small(h, center, 1e-7);
  REQUIRE((small.link() * 1e-14 - Eigen::Matrix2d::Identity()).norm() < 1e-13);
  REQUIRE((small.onsite() - h).norm() < 1e-13);
  REQUIRE(std::abs(small.logMeasureCorrection(1) + 1e-14 * h.trace() / 12) <
          1e-28);
  HarmonicReferenceAction large(h, center, 1000);
  REQUIRE(large.link().allFinite());
  REQUIRE(std::isfinite(large.logMeasureCorrection(8)));
}

TEST_CASE("Mehler reference retains exact free modes without spectral flooring",
          "[instanton][harmonic-reference]") {
  Eigen::Matrix2d h = Eigen::Matrix2d::Zero();
  h(1, 1) = 4.0;
  const double tau = 0.4;
  HarmonicReferenceAction action(h, Eigen::Vector2d::Zero(), tau);
  REQUIRE(action.link()(0, 0) == 1.0 / tau / tau);
  REQUIRE(action.onsite()(0, 0) == 0.0);
  REQUIRE(action.link()(0, 1) == 0.0);
  const double z = 2 * tau;
  REQUIRE(std::abs(action.logMeasureCorrection(3) -
                   1.5 * std::log(z / std::sinh(z))) < 1e-14);
  HarmonicReferenceAction free(Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(),
                               tau);
  REQUIRE(free.logMeasureCorrection(7) == 0.0);
  REQUIRE(free.onsite().isZero());
  REQUIRE(free.link().isApprox(Eigen::Matrix2d::Identity() / tau / tau, 1e-15));
  h(0, 0) = -1e-15;
  REQUIRE_THROWS_AS(HarmonicReferenceAction(h, Eigen::Vector2d::Zero(), tau),
                    std::invalid_argument);
}

TEST_CASE("Mehler reference rejects negative or malformed input",
          "[instanton][harmonic-reference]") {
  Eigen::Matrix2d h = Eigen::Matrix2d::Identity();
  const Eigen::Vector2d center = Eigen::Vector2d::Zero();
  REQUIRE_THROWS_AS(HarmonicReferenceAction(h, center, 0),
                    std::invalid_argument);
  h(0, 0) = -1;
  REQUIRE_THROWS_AS(HarmonicReferenceAction(h, center, 1),
                    std::invalid_argument);
  h.setIdentity();
  h(0, 1) = 0.1;
  REQUIRE_THROWS_AS(HarmonicReferenceAction(h, center, 1),
                    std::invalid_argument);
  h.setIdentity();
  h(0, 0) = std::numeric_limits<double>::infinity();
  REQUIRE_THROWS_AS(HarmonicReferenceAction(h, center, 1),
                    std::invalid_argument);
  HarmonicReferenceAction action(Eigen::Matrix2d::Identity(), center, 1);
  REQUIRE_THROWS_AS(action.correction(Eigen::MatrixXd::Zero(2, 3)),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(action.hessianCorrection(0), std::invalid_argument);
}
