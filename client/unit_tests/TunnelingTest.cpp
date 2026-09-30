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
#include "eon/Tunneling.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

using namespace eonc;
using namespace eonc::tunneling;
using Catch::Matchers::WithinRel;

namespace {

// V(x) = V0 (x^2 / a^2 - 1)^2 on a mass-weighted coordinate, sampled as a
// band of `images` points from one minimum to the other.
Profile quarticBand(double v0, double a, int images) {
  std::vector<double> s, v;
  for (int i = 0; i < images; ++i) {
    const double si = 2.0 * a * i / (images - 1);
    const double x = si - a;
    const double q = x * x / (a * a) - 1.0;
    s.push_back(si);
    v.push_back(v0 * q * q);
  }
  return Profile(s, v);
}

// The exact splitting: the two lowest levels of H = -hbar^2/2 d2/ds2 + V on
// a finite-difference grid wide enough that both states decay to zero.
double exactSplitting(double v0, double a, int n = 900) {
  const double half = 2.2 * a;
  const double h = 2.0 * half / (n - 1);
  const double t = kHbar * kHbar / (2.0 * h * h);
  Eigen::MatrixXd hmat = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    const double x = -half + i * h;
    const double q = x * x / (a * a) - 1.0;
    hmat(i, i) = 2.0 * t + v0 * q * q;
    if (i + 1 < n) {
      hmat(i, i + 1) = -t;
      hmat(i + 1, i) = -t;
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(hmat,
                                                    Eigen::EigenvaluesOnly);
  return es.eigenvalues()(1) - es.eigenvalues()(0);
}

} // namespace

TEST_CASE("hbar in eV amu Angstrom units", "[Tunneling]") {
  const double expected =
      1.054571817e-34 /
      (std::sqrt(1.602176634e-19 * 1.66053906660e-27) * 1e-10);
  REQUIRE_THAT(kHbar, WithinRel(expected, 1e-12));
}

TEST_CASE("The profile passes through the band and never overshoots it",
          "[Tunneling]") {
  const Profile p = quarticBand(0.1, 1.0, 11);
  for (size_t i = 0; i < p.s().size(); ++i) {
    REQUIRE_THAT(p(p.s()[i]), WithinRel(p.v()[i], 1e-12) ||
                                  Catch::Matchers::WithinAbs(p.v()[i], 1e-14));
  }
  // Monotone between the band points: no value outside the bracketing pair.
  for (size_t k = 0; k + 1 < p.s().size(); ++k) {
    const double lo = std::min(p.v()[k], p.v()[k + 1]);
    const double hi = std::max(p.v()[k], p.v()[k + 1]);
    for (int j = 1; j < 10; ++j) {
      const double x = p.s()[k] + (p.s()[k + 1] - p.s()[k]) * j / 10.0;
      REQUIRE(p(x) >= lo - 1e-12);
      REQUIRE(p(x) <= hi + 1e-12);
    }
  }
  REQUIRE_THROWS(Profile({0.0, 0.0}, {0.0, 1.0}));
  REQUIRE_THROWS(Profile({0.0}, {0.0}));
}

TEST_CASE("The well curvature of a 21-image band recovers hbar omega",
          "[Tunneling]") {
  for (double v0 : {0.02, 0.08, 0.3}) {
    const Profile p = quarticBand(v0, 1.0, 21);
    const double exact = kHbar * std::sqrt(8.0 * v0);
    REQUIRE_THAT(hbarOmega(wellCurvature(p, true)), WithinRel(exact, 0.02));
    REQUIRE_THAT(hbarOmega(wellCurvature(p, false)), WithinRel(exact, 0.02));
  }
  REQUIRE_THROWS(hbarOmega(0.0));
}

TEST_CASE("WKB splittings match the exact double well as it deepens",
          "[Tunneling]") {
  // V0 / hbar omega of about 1.5, 2.1 and 3.0.
  const struct {
    double v0;
    double tolerance;
  } cases[] = {{0.08, 0.15}, {0.15, 0.10}, {0.3, 0.05}};
  for (const auto &c : cases) {
    const Profile p = quarticBand(c.v0, 1.0, 21);
    const double hw = hbarOmega(wellCurvature(p, true));
    const Splitting sp = wkbSplitting(p, hw, hw);
    const double exact = exactSplitting(c.v0, 1.0);
    INFO("V0 = " << c.v0 << ", WKB " << sp.delta0 << ", exact " << exact);
    REQUIRE_THAT(sp.delta0, WithinRel(exact, c.tolerance));
    REQUIRE(sp.deepWells);
    REQUIRE_THAT(sp.delta, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(sp.tlsEnergy(), WithinRel(sp.delta0, 1e-12));
  }
}

TEST_CASE("An asymmetric pair tunnels from the higher ground level",
          "[Tunneling]") {
  // A tilt of 5 meV on a 0.15 eV double well.
  std::vector<double> s, v;
  for (int i = 0; i < 21; ++i) {
    const double si = 2.0 * i / 20.0;
    const double x = si - 1.0;
    s.push_back(si);
    v.push_back(0.15 * std::pow(x * x - 1.0, 2) + 0.0025 * (x + 1.0));
  }
  const Profile p(s, v);
  const double hwL = hbarOmega(wellCurvature(p, true));
  const double hwR = hbarOmega(wellCurvature(p, false));
  const Splitting sp = wkbSplitting(p, hwL, hwR);
  REQUIRE_THAT(sp.delta, WithinRel(0.005, 1e-9));
  REQUIRE_THAT(
      sp.referenceEnergy,
      WithinRel(std::max(v.front() + 0.5 * hwL, v.back() + 0.5 * hwR), 1e-12));
  REQUIRE(sp.tlsEnergy() > sp.delta);
  REQUIRE_THAT(sp.tlsEnergy(),
               WithinRel(std::hypot(sp.delta, sp.delta0), 1e-12));
}

TEST_CASE("A shallow double well is flagged", "[Tunneling]") {
  const Profile p = quarticBand(0.005, 1.0, 21);
  const double hw = hbarOmega(wellCurvature(p, true));
  REQUIRE(hw > 0.005);
  REQUIRE_FALSE(wkbSplitting(p, hw, hw).deepWells);
}

TEST_CASE("Mass-weighted distance takes the masses and the minimum image",
          "[Tunneling]") {
  auto params = std::make_shared<Parameters>();
  ParametersLoadAccess::potential_options(*params).potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, *params);
  Matter a(pot, *params);
  a.resize(2);
  Matrix3d cell = Matrix3d::Identity() * 10.0;
  a.setCell(cell);
  AtomMatrix r(2, 3);
  r << 0.5, 5.0, 5.0, //
      5.0, 5.0, 5.0;
  a.setPositions(r);
  VectorXd m(2);
  m << 1.0, 16.0;
  a.setMasses(m);
  Matter b(a);
  AtomMatrix rb = r;
  rb(0, 0) = 9.5; // one Angstrom across the boundary, not nine
  rb(1, 1) = 5.5; // half an Angstrom
  b.setPositions(rb);
  const double expected = std::sqrt(1.0 * 1.0 + 16.0 * 0.25);
  REQUIRE_THAT(massWeightedDistance(a, b), WithinRel(expected, 1e-12));
  VectorXd zero(2);
  zero << 1.0, 0.0;
  a.setMasses(zero);
  REQUIRE_THROWS(massWeightedDistance(a, b));
}
