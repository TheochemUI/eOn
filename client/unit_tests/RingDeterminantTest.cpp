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
#include "eon/RingDeterminant.h"

#include "catch2/catch_amalgamated.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <utility>
#include <vector>

using namespace eonc;
using namespace eonc::tunneling;

namespace {

/// A Gaussian estimate is checked at four standard errors (two-sided
/// p ~ 6e-5) on top of its deterministic error terms.
constexpr double kSigmas = 4.0;

MatrixXd randomSymmetric(long n, std::mt19937_64 &engine, double low,
                         double high) {
  std::normal_distribution<double> normal;
  MatrixXd a(n, n);
  for (long i = 0; i < n; ++i) {
    for (long j = 0; j < n; ++j) {
      a(i, j) = normal(engine);
    }
  }
  const Eigen::HouseholderQR<MatrixXd> qr(a);
  const MatrixXd q = qr.householderQ();
  VectorXd values(n);
  for (long i = 0; i < n; ++i) {
    values(i) = low + (high - low) * static_cast<double>(i) /
                          static_cast<double>(n - 1);
  }
  return q * values.asDiagonal() * q.transpose();
}

/// The dense ring Hessian, as the test oracle only.
MatrixXd denseRing(const std::vector<MatrixXd> &hessians, double c) {
  const long beads = static_cast<long>(hessians.size());
  const long f = hessians.front().rows();
  MatrixXd j = MatrixXd::Zero(beads * f, beads * f);
  const MatrixXd eye = MatrixXd::Identity(f, f);
  for (long b = 0; b < beads; ++b) {
    j.block(b * f, b * f, f, f) = hessians[static_cast<size_t>(b)] + 2.0 * c * eye;
    const long next = (b + 1) % beads;
    j.block(b * f, next * f, f, f) -= c * eye;
    j.block(next * f, b * f, f, f) -= c * eye;
  }
  return j;
}

struct PlantedRing {
  std::vector<MatrixXd> hessians;
  std::vector<VectorXd> cycle;
  std::vector<std::vector<VectorXd>> rigid;
  double c{0.0};
  double exact{0.0};
  /// Rounding in the dense oracle and in an exhausted Krylov space.
  double roundoff{0.0};
  double zero{0.0};
  double negative{0.0};
  long negatives{0};
};

/// Bead Hessians with one rigid direction e shared by every bead, shifted
/// so the ring has one negative eigenvalue and one of size `zero`; the
/// cycle is that eigenvector, perturbed the way a central difference is.
PlantedRing plantedRing(long beads, long f, double c, double zero,
                        std::uint64_t seed) {
  std::mt19937_64 engine(seed);
  VectorXd e = VectorXd::Zero(f);
  e(0) = 1.0;
  const MatrixXd p = MatrixXd::Identity(f, f) - e * e.transpose();
  PlantedRing ring;
  ring.c = c;
  // An orbit's Hessians vary smoothly with imaginary time.
  const MatrixXd base = randomSymmetric(f, engine, 0.5, 4.0);
  const MatrixXd swing = randomSymmetric(f, engine, -0.3, 0.3);
  for (long b = 0; b < beads; ++b) {
    const double phase = 2.0 * std::numbers::pi * static_cast<double>(b) /
                         static_cast<double>(beads);
    MatrixXd h = base + std::sin(phase) * swing;
    h(1, 1) -= 1.5 * (1.0 + std::cos(phase));
    ring.hessians.push_back(p * h * p);
  }
  std::vector<VectorXd> rigid(static_cast<size_t>(beads),
                              e / std::sqrt(static_cast<double>(beads)));
  ring.rigid.push_back(rigid);

  const long n = beads * f;
  MatrixXd basis(n, 1);
  for (long b = 0; b < beads; ++b) {
    basis.block(b * f, 0, f, 1) = rigid[static_cast<size_t>(b)];
  }
  const MatrixXd complement =
      Eigen::HouseholderQR<MatrixXd>(basis).householderQ() *
      MatrixXd::Identity(n, n).rightCols(n - 1);
  auto spectrum = [&]() {
    const MatrixXd reduced =
        complement.transpose() * denseRing(ring.hessians, c) * complement;
    return Eigen::SelfAdjointEigenSolver<MatrixXd>(reduced);
  };
  const double shift = spectrum().eigenvalues()(1) - zero;
  for (auto &h : ring.hessians) {
    h -= shift * p;
  }
  const auto es = spectrum();
  const VectorXd lambda = es.eigenvalues();
  ring.zero = lambda(1);
  ring.negative = lambda(0);
  ring.negatives = 1;
  ring.exact = 0.0;
  const double eigenvalueError = std::numeric_limits<double>::epsilon() *
                                 static_cast<double>(n) *
                                 lambda.cwiseAbs().maxCoeff();
  for (long i = 0; i < lambda.size(); ++i) {
    if (i != 1) {
      ring.exact += std::log(std::abs(lambda(i)));
    }
    ring.roundoff += eigenvalueError / std::abs(lambda(i));
  }
  // A central difference errs by mixing in neighbouring smooth modes.
  VectorXd tau = complement * (es.eigenvectors().col(1) +
                               0.05 * es.eigenvectors().col(2));
  for (long b = 0; b < beads; ++b) {
    ring.cycle.push_back(tau.segment(b * f, f));
  }
  return ring;
}

RingHessianProducts productsOf(const std::vector<MatrixXd> &hessians,
                               long &calls) {
  return [&hessians, &calls](const std::vector<VectorXd> &v,
                             std::vector<VectorXd> &out) {
    out.resize(v.size());
    for (size_t j = 0; j < v.size(); ++j) {
      out[j] = hessians[j] * v[j];
      ++calls;
    }
  };
}

} // namespace

TEST_CASE("ring log-determinant matches the dense spectrum",
          "[ring_determinant]") {
  const double tolerance = 0.05;
  // The small ring fills its Krylov space; the large one has an orbit's
  // spring, nu_1 = 4 c sin^2(pi / N) about the largest bead curvature, so
  // the kinetic preconditioner leaves a remainder for the probes.
  for (const auto &[beads, c] : {std::pair{8L, 3.0}, std::pair{64L, 400.0}}) {
    const PlantedRing ring = plantedRing(beads, 6, c, 0.02, 11);
    long calls = 0;
    const StochasticLogDet result =
        ringLogDetPrime(productsOf(ring.hessians, calls), ring.c, 2.25,
                        ring.cycle, ring.rigid, tolerance, 5);
    INFO("beads " << beads << " value " << result.value << " exact "
                  << ring.exact << " error " << result.error() << " probes "
                  << result.probes << " deflated " << result.deflated);
    CHECK(std::abs(result.value - ring.exact) <=
          kSigmas * result.standardError + result.quadratureError +
              result.modeError + ring.roundoff);
    CHECK(result.error() <= tolerance);
    CHECK(result.negativeModes == ring.negatives);
    // modeError bounds |ln(theta / lambda_0)|.
    CHECK(std::abs(std::log(result.zeroEigenvalue / ring.zero)) <=
          result.modeError + ring.roundoff);
    CHECK(result.lowestEigenvalue < 0.0);
    CHECK(calls == result.products * beads);
    if (beads == 64) {
      CHECK_FALSE(result.exhausted);
      CHECK(result.probes > 0);
    }
    if (result.exhausted) {
      CHECK(result.error() == 0.0);
    }
  }
}

TEST_CASE("ring log-determinant estimates are reproducible from a seed",
          "[ring_determinant]") {
  const PlantedRing ring = plantedRing(32, 6, 3.0, 0.02, 3);
  long calls = 0;
  const auto products = productsOf(ring.hessians, calls);
  const StochasticLogDet a =
      ringLogDetPrime(products, ring.c, 2.0, ring.cycle, ring.rigid, 0.05, 9);
  const StochasticLogDet b =
      ringLogDetPrime(products, ring.c, 2.0, ring.cycle, ring.rigid, 0.05, 9);
  CHECK(a.value == b.value);
  CHECK(a.products == b.products);
}

TEST_CASE("ring log-determinant rejects bad inputs", "[ring_determinant]") {
  const PlantedRing ring = plantedRing(8, 4, 3.0, 0.02, 2);
  long calls = 0;
  const auto products = productsOf(ring.hessians, calls);
  CHECK_THROWS_AS(
      ringLogDetPrime(products, ring.c, 2.0, ring.cycle, ring.rigid, 0.0, 1),
      std::invalid_argument);
  CHECK_THROWS_AS(
      ringLogDetPrime(products, ring.c, 0.0, ring.cycle, ring.rigid, 0.1, 1),
      std::invalid_argument);
  const RingHessianProducts broken = [](const std::vector<VectorXd> &v,
                                        std::vector<VectorXd> &out) {
    out.assign(v.size(), VectorXd());
  };
  CHECK_THROWS_AS(
      ringLogDetPrime(broken, ring.c, 2.0, ring.cycle, ring.rigid, 0.1, 1),
      std::runtime_error);
}

TEST_CASE("well log-determinant sums the shifted spectrum",
          "[ring_determinant]") {
  std::mt19937_64 engine(4);
  const long f = 6;
  VectorXd e = VectorXd::Zero(f);
  e(2) = 1.0;
  const MatrixXd p = MatrixXd::Identity(f, f) - e * e.transpose();
  const MatrixXd h = p * randomSymmetric(f, engine, 0.3, 5.0) * p;
  const std::vector<double> shifts{0.0, 0.7, 0.7, 2.8};
  long calls = 0;
  const StochasticLogDet result = wellLogDet(
      [&](const VectorXd &v) {
        ++calls;
        return VectorXd(h * v);
      },
      f, e, shifts, 0.05, 1);
  const MatrixXd complement =
      Eigen::HouseholderQR<MatrixXd>(MatrixXd(e)).householderQ() *
      MatrixXd::Identity(f, f).rightCols(f - 1);
  const VectorXd lambda = Eigen::SelfAdjointEigenSolver<MatrixXd>(
                              complement.transpose() * h * complement)
                              .eigenvalues();
  // Backward-stable eigenvalues move by about f eps |H|, and ln(l + s)
  // by that over l + s.
  const double eigenvalueError = std::numeric_limits<double>::epsilon() *
                                 static_cast<double>(f) *
                                 lambda.cwiseAbs().maxCoeff();
  double exact = 0.0;
  double bound = 0.0;
  for (long i = 0; i < lambda.size(); ++i) {
    for (const double s : shifts) {
      exact += std::log(lambda(i) + s);
      bound += eigenvalueError / (lambda(i) + s);
    }
  }
  CHECK(result.exhausted);
  CHECK(result.error() == 0.0);
  CHECK(std::abs(result.value - exact) <= bound);
  CHECK(std::abs(result.meanCurvature - lambda.mean()) <= eigenvalueError);
  CHECK(calls == result.products);

  const MatrixXd indefinite = h - 1.0 * p;
  CHECK_THROWS_AS(wellLogDet([&](const VectorXd &v) {
                    return VectorXd(indefinite * v);
                  },
                             f, e, shifts, 0.05, 1),
                  std::runtime_error);
}
