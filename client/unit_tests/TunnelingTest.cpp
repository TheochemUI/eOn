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
  // h / (2 pi) with the exact SI h, over sqrt(eV * dalton) * Angstrom with the
  // CODATA 2022 dalton; kHbar is sollya's correctly rounded double of this.
  const double expected =
      6.62607015e-34 / (2.0 * std::numbers::pi) /
      (std::sqrt(1.602176634e-19 * 1.66053906892e-27) * 1e-10);
  REQUIRE_THAT(kHbar, WithinRel(expected, 1e-14));
  REQUIRE_THAT(kBoltzmann, WithinRel(1.380649e-23 / 1.602176634e-19, 1e-14));
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

TEST_CASE("The WKB action lands in sollya's certified enclosures",
          "[Tunneling]") {
  // Enclosures from data/tunneling/wkb_quartic.sollya (diam = 1e-5).
  const struct {
    double v0, lo, hi;
  } refs[] = {{0.08, 5.063598115650122, 5.063808157855117},
              {0.15, 7.961002615858806, 7.961336193876939},
              {0.3, 12.473664501321133, 12.474192881817101}};
  for (const auto &r : refs) {
    const double level = 0.5 * kHbar * std::sqrt(8.0 * r.v0);
    // A dense band: the quadrature and the interpolant, nothing else.
    const double dense = wkbAction(quarticBand(r.v0, 1.0, 801), level);
    INFO("V0 = " << r.v0 << ", dense action " << dense);
    REQUIRE(dense > r.lo - 1e-5);
    REQUIRE(dense < r.hi + 1e-5);
    // The 21 images a band has: within 3e-4 of the certified midpoint,
    // which moves the splitting by under 0.4 percent.
    const double sparse = wkbAction(quarticBand(r.v0, 1.0, 21), level);
    REQUIRE_THAT(sparse, WithinRel(0.5 * (r.lo + r.hi), 3e-4));
  }
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

namespace {

// V0 (x^2 - 1)^2 + (K / 2) (y - C (1 - x^2))^2 at unit mass: minima at
// (-1, 0) and (1, 0), a valley that bows out to y = C at the barrier.
struct CurvedValley {
  double v0, k, c;
  double value(const VectorXd &q) const {
    const double x = q(0), t = q(1) - c * (1.0 - x * x);
    return v0 * (x * x - 1.0) * (x * x - 1.0) + 0.5 * k * t * t;
  }
  VectorXd gradient(const VectorXd &q) const {
    const double x = q(0), t = q(1) - c * (1.0 - x * x);
    VectorXd g(2);
    g << 4.0 * v0 * x * (x * x - 1.0) + 2.0 * k * c * x * t, k * t;
    return g;
  }
  MatrixXd hessian(const VectorXd &q) const {
    const double x = q(0), t = q(1) - c * (1.0 - x * x);
    MatrixXd h(2, 2);
    h(0, 0) = 4.0 * v0 * (3.0 * x * x - 1.0) + 2.0 * k * c * t +
              4.0 * k * c * c * x * x;
    h(0, 1) = h(1, 0) = 2.0 * k * c * x;
    h(1, 1) = k;
    return h;
  }
  BatchPotential batch() const {
    return [this](const std::vector<VectorXd> &q, std::vector<double> &v,
                  std::vector<VectorXd> &g) {
      v.resize(q.size());
      g.resize(q.size());
      for (size_t i = 0; i < q.size(); ++i) {
        v[i] = value(q[i]);
        g[i] = gradient(q[i]);
      }
    };
  }
};

Instanton valleyInstanton(const CurvedValley &pes, long beads,
                          double betaHbarOmega) {
  VectorXd a(2), b(2);
  a << -1.0, 0.0;
  b << 1.0, 0.0;
  const double omega = pathOmega(pes.hessian(a), pes.hessian(b), a, b);
  InstantonOptions opt;
  opt.beads = beads;
  opt.betaHbarOmega = betaHbarOmega;
  opt.forceTolerance = 1e-9;
  opt.maxIterations = 20000;
  Instanton inst =
      optimizeInstanton(a, b, betaHbarOmega / omega, {}, pes.batch(), opt);
  instantonSplitting(
      inst, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(a), pes.hessian(b));
  return inst;
}

} // namespace

// Exact splittings of the two-dimensional valley from the lowest two
// eigenvalues of a fourth-order finite-difference Hamiltonian on a 361 x 181
// grid over [-2.2, 2.2] x [-1.2, 1.6] (converged to 1e-5 against 241 x 121),
// K = 4 eV / Angstrom^2, unit mass. The instanton values beside them are an
// independent implementation of the same discretisation (scipy L-BFGS to a
// gradient of 1e-12, dense eigenvalues of the chain Hessian) at P = 512 and
// beta hbar omega = 40.
TEST_CASE("Instanton splitting in a curved valley matches the exact gap",
          "[Tunneling][Instanton]") {
  struct Case {
    double v0, c, exact, reference;
  };
  const Case cases[] = {
      {0.12, 0.35, 2.5195793978e-05, 2.949267463340166e-05},
      {0.30, 0.35, 9.8178409097e-08, 1.0623251042029217e-07},
      {0.12, 0.00, 2.0396238276e-05, 2.2784988974931514e-05},
      {0.30, 0.00, 1.1953758103e-07, 1.2776094740454857e-07},
  };
  for (const auto &cs : cases) {
    CAPTURE(cs.v0, cs.c);
    const CurvedValley pes{cs.v0, 4.0, cs.c};
    const Instanton inst = valleyInstanton(pes, 512, 40.0);
    REQUIRE(inst.converged);
    REQUIRE(inst.symmetricEnough);
    REQUIRE(inst.modeSeparation > 1e3);
    REQUIRE_THAT(inst.delta0, WithinRel(cs.reference, 2e-3));
    // Semiclassical error: 17 and 12 percent at V0 = 0.12 eV (V0 about
    // 2.2 hbar omega), 8 and 7 percent at 0.3 eV.
    REQUIRE_THAT(inst.delta0, WithinRel(cs.exact, cs.v0 < 0.2 ? 0.2 : 0.09));
  }
}

TEST_CASE("The instanton cuts the corner the minimum energy path takes",
          "[Tunneling][Instanton]") {
  const CurvedValley pes{0.12, 4.0, 0.35};
  const Instanton inst = valleyInstanton(pes, 256, 30.0);
  REQUIRE(inst.converged);
  // Halfway in imaginary time the instanton sits inside the valley's bow.
  const VectorXd &mid = inst.path[inst.path.size() / 2];
  REQUIRE(std::abs(mid(0)) < 1e-4);
  REQUIRE(mid(1) > 0.0);
  REQUIRE(mid(1) < 0.35);
  // One-dimensional WKB along the valley floor misses the gap by a factor
  // of 2.8; the instanton lands within 20 percent.
  std::vector<double> s, v;
  VectorXd prev(2);
  for (int i = 0; i <= 400; ++i) {
    const double x = -1.0 + 2.0 * i / 400.0;
    VectorXd q(2);
    q << x, 0.35 * (1.0 - x * x);
    s.push_back(i == 0 ? 0.0 : s.back() + (q - prev).norm());
    v.push_back(pes.value(q));
    prev = q;
  }
  const Profile floor(s, v);
  const double hw = hbarOmega(wellCurvature(floor, true));
  const Splitting wkb = wkbSplitting(floor, hw, hw);
  const double exact = 2.5195793978e-05;
  REQUIRE(wkb.delta0 < 0.5 * exact);
  REQUIRE_THAT(inst.delta0, WithinRel(exact, 0.2));
}

TEST_CASE("Instanton inputs are checked", "[Tunneling][Instanton]") {
  const CurvedValley pes{0.12, 4.0, 0.0};
  VectorXd a(2);
  a << -1.0, 0.0;
  InstantonOptions opt;
  opt.beads = 2;
  REQUIRE_THROWS_AS(optimizeInstanton(a, a, 1.0, {}, pes.batch(), opt),
                    std::invalid_argument);
  opt.beads = 16;
  REQUIRE_THROWS_AS(optimizeInstanton(a, a, 1.0, {}, pes.batch(), opt),
                    std::invalid_argument);
  MatrixXd flat = MatrixXd::Zero(2, 2);
  REQUIRE_THROWS_AS(pathOmega(flat, flat, a, -a), std::invalid_argument);
}
