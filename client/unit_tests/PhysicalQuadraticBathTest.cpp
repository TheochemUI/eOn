/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#include "catch2/catch_amalgamated.hpp"
#include "eon/DerivativeRbfModel.h"
#include "eon/Tunneling.h"

using namespace eonc::tunneling;

TEST_CASE("derivative RBF site signs and physical derivatives agree with "
          "Gaussian calculus",
          "[instanton][physical-bath][interpolation]") {
  const double x = 0.37, site = -0.2, length = 0.8;
  for (int order : {0, 1, 2}) {
    const nlohmann::json data = {
        {"kind", "derivative_rbf"}, {"length_scales", {length}},
        {"constant", {0.0}},        {"sites", {{site}}},
        {"weights", {{1.0}}},       {"derivative_orders", {{order}}}};
    const auto model = DerivativeRbfModel::fromJson(data);
    const auto value = model.evaluate(VectorXd::Constant(1, x));
    const double dx = x - site, inv = 1 / (length * length),
                 gaussian = std::exp(-0.5 * dx * dx * inv);
    const double expected =
        order == 0 ? gaussian
                   : (order == 1 ? dx * inv * gaussian
                                 : (dx * dx * inv * inv - inv) * gaussian);
    REQUIRE(std::abs(value.value[0] - expected) < 1e-14);
    const double step = 1e-5;
    const auto plus = model.evaluate(VectorXd::Constant(1, x + step));
    const auto minus = model.evaluate(VectorXd::Constant(1, x - step));
    REQUIRE(std::abs((plus.value[0] - minus.value[0]) / (2 * step) -
                     value.jacobian(0, 0)) < 1e-8);
    REQUIRE(std::abs((plus.jacobian(0, 0) - minus.jacobian(0, 0)) / (2 * step) -
                     value.hessians[0](0, 0)) < 1e-8);
  }
}

TEST_CASE("finite physical bath action equals the full primitive Schur model",
          "[instanton][physical-bath]") {
  PhysicalQuadraticBath bath;
  bath.coupling = MatrixXd::Constant(1, 2, 0.3);
  bath.curvatures = VectorXd(2);
  bath.curvatures << 2.0, 5.0;
  const long n = 7;
  const double c = 1.7;
  std::vector<VectorXd> q(n, VectorXd::Zero(1));
  for (long j = 0; j < n; ++j)
    q[j][0] = 0.3 + std::cos(2 * std::numbers::pi * j / n);
  const auto ev = bath.evaluate(q, c, true);
  MatrixXd db = MatrixXd::Zero(2 * n, 2 * n), cross = MatrixXd::Zero(n, 2 * n);
  VectorXd source(2 * n), bp(2 * n);
  for (long j = 0; j < n; ++j) {
    const long k = (j + 1) % n;
    db.block<2, 2>(2 * j, 2 * j) += bath.curvatures.asDiagonal();
    for (long a = 0; a < 2; ++a) {
      db(2 * j + a, 2 * j + a) += c;
      db(2 * k + a, 2 * k + a) += c;
      db(2 * j + a, 2 * k + a) -= c;
      db(2 * k + a, 2 * j + a) -= c;
    }
    source.segment<2>(2 * j) = bath.coupling.transpose() * q[j];
    bp.segment<2>(2 * j) = ev.bathPositions[j];
    cross.block(j, 2 * j, 1, 2) = bath.coupling;
  }
  REQUIRE((db * bp + source).norm() < 1e-13);
  REQUIRE(std::abs(ev.energy - (0.5 * bp.dot(db * bp) + source.dot(bp))) <
          1e-13);
  REQUIRE(
      ev.hessian.isApprox(-cross * db.inverse() * cross.transpose(), 1e-13));
  REQUIRE(std::abs(ev.bathLogDeterminant - std::log(db.determinant())) < 1e-12);
}

TEST_CASE(
    "nonlinear physical bath source supplies the complete retained Hessian",
    "[instanton][physical-bath]") {
  PhysicalQuadraticBath bath;
  bath.coupling = MatrixXd::Zero(1, 1);
  bath.curvatures = VectorXd::Constant(1, 2.3);
  bath.source = [](const VectorXd &q, VectorXd &s, MatrixXd &j,
                   std::vector<MatrixXd> &h) {
    s = VectorXd::Constant(1, q[0] * q[0] - 0.4);
    j = MatrixXd::Constant(1, 1, 2 * q[0]);
    h.assign(1, MatrixXd::Constant(1, 1, 2));
  };
  for (bool relaxed : {false, true}) {
    bath.relaxedPotential = relaxed;
    std::vector<VectorXd> q{VectorXd::Constant(1, -0.7),
                            VectorXd::Constant(1, 0.2),
                            VectorXd::Constant(1, 0.8)};
    const auto ev = bath.evaluate(q, 1.2, true);
    for (long k = 0; k < 3; ++k) {
      const double step = 1e-5;
      auto plus = q, minus = q;
      plus[k][0] += step;
      minus[k][0] -= step;
      const auto ep = bath.evaluate(plus, 1.2), em = bath.evaluate(minus, 1.2);
      REQUIRE(std::abs((ep.energy - em.energy) / (2 * step) -
                       ev.gradient[k][0]) < 1e-9);
      for (long j = 0; j < 3; ++j)
        REQUIRE(std::abs((ep.gradient[j][0] - em.gradient[j][0]) / (2 * step) -
                         ev.hessian(j, k)) < 1e-9);
    }
  }
}

TEST_CASE("physical bath phase coarea preserves the full Gaussian rate weight",
          "[instanton][physical-bath]") {
  MatrixXd d(2, 2);
  d << 3.0, 0.2, 0.2, 2.0;
  MatrixXd c(3, 2);
  c << 0.3, 0.2, -0.1, 0.4, 0.2, -0.2;
  MatrixXd s = MatrixXd::Zero(3, 3);
  s.diagonal() << -2.0, 0.0, 4.0;
  MatrixXd full(5, 5);
  full << s + c * d.inverse() * c.transpose(), c, c.transpose(), d;
  VectorXd ur = VectorXd::Unit(3, 1), uf(5);
  uf << ur, -d.inverse() * c.transpose() * ur;
  REQUIRE((full * uf).norm() < 1e-14);
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(full);
  double fullPrime = 1;
  long zeros = 0;
  for (double value : es.eigenvalues()) {
    if (std::abs(value) < 1e-13)
      ++zeros;
    else
      fullPrime *= value;
  }
  REQUIRE(zeros == 1);
  const double reducedPrime = -8.0;
  REQUIRE(std::abs(fullPrime - d.determinant() * reducedPrime *
                                   uf.squaredNorm() / ur.squaredNorm()) <
          1e-12);
  const double fullWeight = uf.norm() / std::sqrt(std::abs(fullPrime));
  const double reducedWeight =
      ur.norm() / std::sqrt(std::abs(d.determinant() * reducedPrime));
  REQUIRE(std::abs(fullWeight - reducedWeight) < 1e-14);
}

TEST_CASE("physical bath optimizes and normalizes a coupled quartic rate ring",
          "[instanton][physical-bath]") {
  const double coupling = 0.35, e = 4.0,
               xr = -std::sqrt(1 + coupling * coupling / e),
               yr = -coupling * xr / e;
  const double vR =
      0.25 * std::pow(xr * xr - 1, 2) + 0.5 * e * yr * yr + coupling * xr * yr;
  BatchPotential frozen = [=](const std::vector<VectorXd> &q,
                              std::vector<double> &v,
                              std::vector<VectorXd> &g) {
    v.resize(q.size());
    g.resize(q.size());
    for (size_t j = 0; j < q.size(); ++j) {
      const double x = xr + q[j][0];
      v[j] = 0.25 * std::pow(x * x - 1, 2) + 0.5 * e * yr * yr +
             coupling * x * yr - vR;
      g[j] = VectorXd::Constant(1, x * (x * x - 1) + coupling * yr);
    }
  };
  auto bath = std::make_shared<PhysicalQuadraticBath>();
  bath->coupling = MatrixXd::Constant(1, 1, coupling);
  bath->curvatures = VectorXd::Constant(1, e);
  RateInstantonOptions options;
  options.beads = 64;
  options.maxIterations = 500;
  options.forceTolerance = 1e-9;
  options.physicalBath = bath;
  options.halfRing = false;
  options.checkOddSector = false;
  options.maxStep = 0.2;
  options.collapseFraction = 0;
  const double beta = 10 / kHbar;
  std::vector<VectorXd> guess(options.beads, VectorXd::Zero(1));
  for (long j = 0; j < options.beads; ++j)
    guess[j][0] =
        -xr + 0.9 * std::cos(2 * std::numbers::pi * j / options.beads);
  const auto ring = optimizeRateInstanton(VectorXd::Constant(1, -xr),
                                          MatrixXd::Constant(1, 1, -1), beta,
                                          guess, frozen, options);
  REQUIRE(ring.converged);
  MatrixXd physicalSaddle(2,2);
  physicalSaddle<<-1.0,coupling,coupling,e;
  REQUIRE(std::abs(ring.crossover-crossoverTemperature(physicalSaddle))<1e-12);
  REQUIRE(ring.crossover<crossoverTemperature(MatrixXd::Constant(1,1,-1-coupling*coupling/e)));
  const double spring = 1 / std::pow(beta / options.beads * kHbar, 2);
  const auto eliminated = bath->evaluate(ring.beads, spring);
  double fullEnergy = 0;
  for (long j = 0; j < options.beads; ++j) {
    const long next = (j + 1) % options.beads;
    const double r = ring.beads[j][0], b = eliminated.bathPositions[j][0],
                 x = xr + r, y = yr + b;
    fullEnergy +=
        0.25 * std::pow(x * x - 1, 2) + 0.5 * e * y * y + coupling * x * y - vR;
    fullEnergy += 0.5 * spring *
                  (std::pow(ring.beads[next][0] - r, 2) +
                   std::pow(eliminated.bathPositions[next][0] - b, 2));
  }
  REQUIRE(std::abs(fullEnergy - ring.ringPotential) < 1e-11);
  auto rate = ring;
  RingBeadHessian hessian = [=](long, const VectorXd &q) {
    return MatrixXd::Constant(1, 1, 3 * std::pow(xr + q[0], 2) - 1);
  };
  instantonRate(rate, hessian, MatrixXd::Constant(1, 1, 3 * xr * xr - 1), 0);
  REQUIRE(rate.negativeModes == 1);
  REQUIRE(rate.phaseZeroVerified);
  REQUIRE(rate.rateValid);
  REQUIRE(std::isfinite(rate.logRate));
}
