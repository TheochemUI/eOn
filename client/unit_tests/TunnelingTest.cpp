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
#include "eon/RingPolymerPotential.h"
#include "eon/Tunneling.h"
#include "EckartBarrier.hpp"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
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

TEST_CASE("A barrier path becomes a closed ring at the requested period",
          "[Tunneling][Instanton]") {
  const double a = 1.0;
  const double v0 = 1.0;
  const int images = 401;
  std::vector<VectorXd> path;
  std::vector<double> energies;
  path.reserve(static_cast<size_t>(images));
  for (int i = 0; i < images; ++i) {
    const double x = -a + 2.0 * a * static_cast<double>(i) / (images - 1);
    const double q = x * x / (a * a) - 1.0;
    path.push_back(VectorXd::Constant(1, x));
    energies.push_back(v0 * q * q);
  }
  // The full period at the barrier top is pi. A shorter imaginary time is
  // above the crossover along this path.
  const long n = 16;
  REQUIRE_THROWS_AS(ringFromPath(path, energies, 1.0, n),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(ringFromPath(path, energies, 30.0, 3),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(
      ringFromPath(path, std::vector<double>(path.size(), 1.0), 30.0, n),
      std::invalid_argument);

  const auto ring = ringFromPath(path, energies, 30.0, n);
  REQUIRE(static_cast<long>(ring.size()) == n);
  REQUIRE(ring.front()(0) < 0.0);
  REQUIRE(ring[static_cast<size_t>(n / 2)](0) > 0.0);
  REQUIRE(std::abs(ring.front()(0)) > 0.2);
  for (long j = 1; j < n / 2; ++j) {
    REQUIRE(ring[static_cast<size_t>(j)](0) ==
            ring[static_cast<size_t>(n - j)](0));
    REQUIRE(ring[static_cast<size_t>(j)](0) >
            ring[static_cast<size_t>(j - 1)](0));
  }
}

TEST_CASE("The one-dimensional WKB rate falls as the barrier grows",
          "[Tunneling]") {
  const Profile low = quarticBand(0.3, 1.0, 201);
  const Profile high = quarticBand(1.0, 1.0, 201);
  const double hwLow = hbarOmega(wellCurvature(low, true));
  const double hwHigh = hbarOmega(wellCurvature(high, true));
  const double coldLow = wkbLogRateAlongPath(low, 40.0, hwLow);
  const double coldHigh = wkbLogRateAlongPath(high, 40.0, hwHigh);
  const double hotHigh = wkbLogRateAlongPath(high, 10.0, hwHigh);
  REQUIRE(std::isfinite(coldLow));
  REQUIRE(std::isfinite(coldHigh));
  REQUIRE(coldHigh < coldLow);
  REQUIRE(hotHigh > coldHigh);
  REQUIRE_THROWS_AS(wkbLogRateAlongPath(low, 0.0, hwLow),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(wkbLogRateAlongPath(
                        Profile({0.0, 1.0, 2.0}, {1.0, 1.0, 1.0}), 40.0, hwLow),
                    std::invalid_argument);
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
  REQUIRE_THAT(sp.asymmetry(), WithinRel(sp.delta + 0.5 * (hwR - hwL), 1e-12));
  REQUIRE_THAT(sp.tlsEnergy(),
               WithinRel(std::hypot(sp.asymmetry(), sp.delta0), 1e-12));
}

TEST_CASE("The TLS energy of a tilted double well carries the zero-point "
          "difference",
          "[Tunneling]") {
  // V0 (x^2 - 1)^2 + eps x / 2 sampled minimum to minimum. The tilt also
  // changes the curvature of each well, so the local ground states differ
  // by more than the minima: with the minima alone the energy sits 6
  // percent above the exact gap at every image count.
  const double v0 = 0.12;
  const double eps = 1e-3;
  auto v = [&](double x) {
    return v0 * (x * x - 1.0) * (x * x - 1.0) + 0.5 * eps * x;
  };
  auto minimum = [&](double x) {
    for (int k = 0; k < 50; ++k) {
      x -= (4.0 * v0 * x * (x * x - 1.0) + 0.5 * eps) /
           (v0 * (12.0 * x * x - 4.0));
    }
    return x;
  };
  const double xl = minimum(-1.0);
  const double xr = minimum(1.0);
  const int n = 1600;
  const double half = 2.4;
  const double h = 2.0 * half / (n - 1);
  const double t = kHbar * kHbar / (2.0 * h * h);
  Eigen::MatrixXd hmat = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    hmat(i, i) = 2.0 * t + v(-half + i * h);
    if (i + 1 < n) {
      hmat(i, i + 1) = -t;
      hmat(i + 1, i) = -t;
    }
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      hmat, Eigen::EigenvaluesOnly);
  const double exact = es.eigenvalues()(1) - es.eigenvalues()(0);
  for (int images : {21, 41}) {
    std::vector<double> s, e;
    for (int i = 0; i < images; ++i) {
      const double x = xl + (xr - xl) * i / (images - 1);
      s.push_back(x - xl);
      e.push_back(v(x));
    }
    const Profile p(s, e);
    const Splitting sp = wkbSplitting(p, hbarOmega(wellCurvature(p, true)),
                                      hbarOmega(wellCurvature(p, false)));
    INFO(images << " images: " << sp.tlsEnergy() << " against " << exact);
    REQUIRE_THAT(sp.tlsEnergy(), WithinRel(exact, 0.015));
    REQUIRE(std::abs(std::hypot(sp.delta, sp.delta0) / exact - 1.0) > 0.04);
  }
}

namespace {

// E2 - E1 of -hbar^2/2 d2/dx2 + V on [-half, half] from second-order
// differences on n points.
double gridGap(const std::function<double(double)> &v, double half,
               int n = 1600) {
  const double h = 2.0 * half / (n - 1);
  const double t = kHbar * kHbar / (2.0 * h * h);
  Eigen::MatrixXd hmat = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    hmat(i, i) = 2.0 * t + v(-half + i * h);
    if (i + 1 < n) {
      hmat(i, i + 1) = -t;
      hmat(i + 1, i) = -t;
    }
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      hmat, Eigen::EigenvaluesOnly);
  return es.eigenvalues()(1) - es.eigenvalues()(0);
}

// The minimum or the barrier top of V0 (x^2 - 1)^2 + eps x / 2 nearest x.
double tiltedRoot(double v0, double eps, double x) {
  for (int k = 0; k < 60; ++k) {
    x -= (4.0 * v0 * x * (x * x - 1.0) + 0.5 * eps) /
         (v0 * (12.0 * x * x - 4.0));
  }
  return x;
}

} // namespace

TEST_CASE("The one-dimensional levels of a double well match a fine grid",
          "[Tunneling]") {
  // A linear tilt leaves the tunnelling matrix element at the gap of the
  // level wells, so delta0 has to stay there however large the asymmetry
  // grows, while the energy follows the exact gap.
  const double v0 = 0.12;
  const double level = gridGap(
      [&](double x) { return v0 * (x * x - 1.0) * (x * x - 1.0); }, 2.4);
  for (double eps : {0.0, 1e-4, 1e-3}) {
    CAPTURE(eps);
    auto v = [&](double x) {
      return v0 * (x * x - 1.0) * (x * x - 1.0) + 0.5 * eps * x;
    };
    const double xl = tiltedRoot(v0, eps, -1.0);
    const double xr = tiltedRoot(v0, eps, 1.0);
    const double top = tiltedRoot(v0, eps, 0.0);
    auto hw = [&](double x) {
      return kHbar * std::sqrt(v0 * (12.0 * x * x - 4.0));
    };
    const Levels lv = dvrLevels(v, xl, xr, top, hw(xl), hw(xr));
    REQUIRE_THAT(lv.gap, WithinRel(gridGap(v, 2.4), 2e-4));
    REQUIRE_THAT(lv.delta0, WithinRel(level, 1e-3));
    REQUIRE_THAT(std::hypot(lv.asymmetry, lv.delta0), WithinRel(lv.gap, 1e-12));
    REQUIRE(lv.asymmetry > -1e-9);
  }
  REQUIRE_THROWS_AS(
      dvrLevels([](double) { return 0.0; }, 1.0, -1.0, 0.0, 0.1, 0.1),
      std::invalid_argument);
}

TEST_CASE("The band's levels continue its wells and take sampled walls",
          "[Tunneling]") {
  // A band from minimum to minimum never sees the outer walls. Continued by
  // each well's fit the levels sit within 6 percent of the exact gap; with
  // eight samples of the real wall past each end, within 1.5 percent, the
  // rest being the band's own images.
  const double v0 = 0.12;
  auto v = [&](double x) { return v0 * (x * x - 1.0) * (x * x - 1.0); };
  const double exact = gridGap(v, 2.4);
  for (int images : {13, 41}) {
    CAPTURE(images);
    const Profile p = quarticBand(v0, 1.0, images);
    REQUIRE(singleBarrier(p));
    REQUIRE_THAT(bandLevels(p).gap, WithinRel(exact, 0.06));
    const double ell = kHbar / std::sqrt(hbarOmega(2.0 * wellFit(p, true)[0]));
    BandWalls walls;
    for (int k = 1; k <= 8; ++k) {
      const double u = kDvrPadLengths * ell * k / 8.0;
      walls.reactant.distance.push_back(u);
      walls.reactant.energy.push_back(v(-1.0 - u));
      walls.product.distance.push_back(u);
      walls.product.energy.push_back(v(1.0 + u));
    }
    REQUIRE_THAT(bandLevels(p, &walls).gap, WithinRel(exact, 0.015));
  }
  REQUIRE_THAT(
      wellCurvature(quarticBand(v0, 1.0, 41), true),
      WithinRel(2.0 * wellFit(quarticBand(v0, 1.0, 41), true)[0], 1e-15));
  // Two barriers with a well between them are no two-level system.
  REQUIRE_FALSE(singleBarrier(
      Profile({0.0, 1.0, 2.0, 3.0, 4.0}, {0.0, 0.1, 0.02, 0.1, 0.0})));
  BandWalls bad;
  bad.reactant.distance = {0.1, 0.05};
  bad.reactant.energy = {0.01, 0.02};
  REQUIRE_THROWS_AS(bandLevels(quarticBand(v0, 1.0, 21), &bad),
                    std::invalid_argument);
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
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, *params));
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
  opt.forceTolerance = 1e-7;
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

namespace {

// V0 (x^2 - 1)^2 + (K / 2) (1 + kappa x) y^2 + eps x / 2 at unit mass. The
// transverse stiffness differs between the wells by 2 kappa K, so their
// zero-point energies differ too, and eps tilts one well below the other.
// The y mode is harmonic at every x, so the lowest levels are those of the
// one-dimensional V0 (x^2 - 1)^2 + eps x / 2 + hbar sqrt(K (1 + kappa x)) / 2
// to 1e-5 of the gap.
struct SkewedValley {
  double v0, k, kappa, eps;
  double value(const VectorXd &q) const {
    const double x = q(0), y = q(1);
    return v0 * (x * x - 1.0) * (x * x - 1.0) +
           0.5 * k * (1.0 + kappa * x) * y * y + 0.5 * eps * x;
  }
  VectorXd gradient(const VectorXd &q) const {
    const double x = q(0), y = q(1);
    VectorXd g(2);
    g << 4.0 * v0 * x * (x * x - 1.0) + 0.5 * k * kappa * y * y + 0.5 * eps,
        k * (1.0 + kappa * x) * y;
    return g;
  }
  MatrixXd hessian(const VectorXd &q) const {
    const double x = q(0), y = q(1);
    MatrixXd h(2, 2);
    h(0, 0) = v0 * (12.0 * x * x - 4.0);
    h(0, 1) = h(1, 0) = k * kappa * y;
    h(1, 1) = k * (1.0 + kappa * x);
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
  VectorXd minimum(double x) const {
    for (int it = 0; it < 60; ++it) {
      x -= (4.0 * v0 * x * (x * x - 1.0) + 0.5 * eps) /
           (v0 * (12.0 * x * x - 4.0));
    }
    VectorXd q(2);
    q << x, 0.0;
    return q;
  }
  double adiabaticGap() const {
    const int n = 1600;
    const double half = 2.4;
    const double h = 2.0 * half / (n - 1);
    const double t = kHbar * kHbar / (2.0 * h * h);
    Eigen::MatrixXd hmat = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) {
      const double x = -half + i * h;
      hmat(i, i) = 2.0 * t + v0 * (x * x - 1.0) * (x * x - 1.0) +
                   0.5 * eps * x +
                   0.5 * kHbar * std::sqrt(k * (1.0 + kappa * x));
      if (i + 1 < n) {
        hmat(i, i + 1) = -t;
        hmat(i + 1, i) = -t;
      }
    }
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
        hmat, Eigen::EigenvaluesOnly);
    return es.eigenvalues()(1) - es.eigenvalues()(0);
  }
};

} // namespace

TEST_CASE("The zero-point difference sums the vibrations of each minimum",
          "[Tunneling][Instanton]") {
  // Two translations at zero (one slightly negative, as a finite-difference
  // Hessian leaves it), then the vibrations; a soft negative mode at the end
  // well is no vibration.
  MatrixXd a = MatrixXd::Zero(5, 5);
  MatrixXd b = MatrixXd::Zero(5, 5);
  a.diagonal() << 1e-9, -2e-9, 4.0, 9.0, 1.0;
  b.diagonal() << -1e-9, 3e-9, 4.41, 8.0, -0.25;
  const double expected =
      0.5 * kHbar * ((2.1 + std::sqrt(8.0)) - (2.0 + 3.0 + 1.0));
  REQUIRE_THAT(zeroPointDifference(a, b, 2), WithinRel(expected, 1e-12));
  REQUIRE_THAT(zeroPointDifference(b, a, 2), WithinRel(-expected, 1e-12));
  REQUIRE_THAT(zeroPointDifference(a, a, 2),
               Catch::Matchers::WithinAbs(0.0, 1e-15));
  REQUIRE_THROWS_AS(zeroPointDifference(a, MatrixXd::Zero(4, 4), 0),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(zeroPointDifference(a, b, 6), std::invalid_argument);
}

TEST_CASE("The instanton TLS energy carries the zero-point asymmetry of every "
          "mode",
          "[Tunneling][Instanton]") {
  // Degenerate minima whose transverse stiffness differs by 10 percent: the
  // transverse zero-point energies split the levels by 3 meV, 140 times the
  // tunnelling splitting. hypot(dV, delta0) alone reads under 1 percent of
  // the gap; with the harmonic zero-point difference the energy lands within
  // the 6 percent the harmonic picture leaves.
  const SkewedValley pes{0.12, 4.0, 0.05, 0.0};
  const VectorXd a = pes.minimum(-1.0);
  const VectorXd b = pes.minimum(1.0);
  const double omega = pathOmega(pes.hessian(a), pes.hessian(b), a, b);
  InstantonOptions opt;
  opt.beads = 192;
  opt.forceTolerance = 1e-7;
  opt.maxIterations = 20000;
  Instanton inst = optimizeInstanton(a, b, 30.0 / omega, {}, pes.batch(), opt);
  REQUIRE(inst.converged);
  REQUIRE(inst.symmetricEnough);
  instantonSplitting(
      inst, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(a), pes.hessian(b));
  const double zpe = zeroPointDifference(pes.hessian(a), pes.hessian(b), 0);
  REQUIRE_THAT(
      zpe,
      WithinRel(0.5 * kHbar * (std::sqrt(4.0 * 1.05) - std::sqrt(4.0 * 0.95)),
                1e-10));
  const double exact = pes.adiabaticGap();
  CAPTURE(inst.delta0, zpe, exact);
  REQUIRE_THAT(std::hypot(inst.asymmetry + zpe, inst.delta0),
               WithinRel(exact, 0.08));
  REQUIRE(std::hypot(inst.asymmetry, inst.delta0) < 0.01 * exact);
}

TEST_CASE("The energy bias levels two wells and leaves both minima alone",
          "[Tunneling][Instanton]") {
  VectorXd a(2), b(2);
  a << -1.0, 0.2;
  b << 1.5, -0.4;
  const EnergyBias bias(a, b, 0.03);
  REQUIRE(bias.value(a) == 0.0);
  REQUIRE_THAT(bias.value(b), WithinRel(0.03, 1e-15));
  for (const VectorXd &end : {a, b}) {
    REQUIRE(bias.gradient(end).norm() == 0.0);
    REQUIRE(bias.hessian(end).norm() == 0.0);
  }
  // Past either minimum the bias stays flat.
  REQUIRE(bias.value(a - 0.3 * (b - a)) == 0.0);
  REQUIRE_THAT(bias.value(b + 0.3 * (b - a)), WithinRel(0.03, 1e-15));
  // Gradient and Hessian against central differences inside the segment.
  VectorXd q(2);
  q << 0.1, 0.35;
  const double h = 1e-5;
  for (int i = 0; i < 2; ++i) {
    VectorXd e = VectorXd::Zero(2);
    e(i) = h;
    REQUIRE_THAT(
        bias.gradient(q)(i),
        WithinRel((bias.value(q + e) - bias.value(q - e)) / (2 * h), 1e-7));
    const VectorXd column =
        (bias.gradient(q + e) - bias.gradient(q - e)) / (2 * h);
    REQUIRE(((bias.hessian(q).col(i) - column).norm()) <
            1e-6 * bias.hessian(q).norm());
  }
  REQUIRE_THROWS_AS(EnergyBias(a, a, 0.1), std::invalid_argument);
}

TEST_CASE("Wells of different depth tunnel on the levelled surface",
          "[Tunneling][Instanton]") {
  // The tilted quartic has beta |dV| = 0.47, past the 0.1 the propagator
  // ratio of the bare surface needs; the skewed valley adds a 3 meV
  // zero-point split to a smaller tilt. On the levelled surface the matrix
  // element stays at that of level wells, and with the zero-point asymmetry
  // the energy lands on the exact gap.
  struct Case {
    double kappa, eps, tolerance;
  };
  for (const Case c : {Case{0.0, 1e-3, 0.02}, Case{0.05, 2e-4, 0.08}}) {
    CAPTURE(c.kappa, c.eps);
    const SkewedValley pes{0.12, 4.0, c.kappa, c.eps};
    const VectorXd a = pes.minimum(-1.0);
    const VectorXd b = pes.minimum(1.0);
    const double dV = pes.value(b) - pes.value(a);
    const EnergyBias bias(a, b, dV);
    const BatchPotential bare = pes.batch();
    const BatchPotential surface = [&](const std::vector<VectorXd> &q,
                                       std::vector<double> &v,
                                       std::vector<VectorXd> &g) {
      bare(q, v, g);
      for (size_t j = 0; j < q.size(); ++j) {
        v[j] -= bias.value(q[j]);
        g[j] -= bias.gradient(q[j]);
      }
    };
    const double omega = pathOmega(pes.hessian(a), pes.hessian(b), a, b);
    const double betaHbar = 30.0 / omega;
    InstantonOptions opt;
    opt.beads = 192;
    opt.forceTolerance = 1e-7;
    opt.maxIterations = 20000;
    Instanton inst = optimizeInstanton(a, b, betaHbar, {}, surface, opt);
    REQUIRE(inst.converged);
    REQUIRE(inst.symmetricEnough);
    instantonSplitting(
        inst,
        [&](long, const VectorXd &q) {
          return MatrixXd(pes.hessian(q) - bias.hessian(q));
        },
        pes.hessian(a), pes.hessian(b));
    const double zpe = zeroPointDifference(pes.hessian(a), pes.hessian(b), 0);
    const double exact = pes.adiabaticGap();
    CAPTURE(dV, zpe, inst.delta0, exact);
    REQUIRE_THAT(std::hypot(dV + zpe, inst.delta0),
                 WithinRel(exact, c.tolerance));
    if (c.kappa == 0.0) {
      REQUIRE(std::abs(dV) * betaHbar / kHbar > 0.1);
      const SkewedValley level{0.12, 4.0, 0.0, 0.0};
      Instanton flat =
          optimizeInstanton(level.minimum(-1.0), level.minimum(1.0), betaHbar,
                            {}, level.batch(), opt);
      instantonSplitting(
          flat, [&](long, const VectorXd &q) { return level.hessian(q); },
          level.hessian(level.minimum(-1.0)),
          level.hessian(level.minimum(1.0)));
      REQUIRE_THAT(inst.delta0, WithinRel(flat.delta0, 0.01));
    }
  }
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

namespace {

// Metastable cubic well on a mass-weighted coordinate,
// V(q) = omega0^2 q^2 / 2 - g q^3 / 3, reactant at 0, barrier at
// q_b = omega0^2 / g with V_b = omega0^6 / (6 g^2) and curvature -omega0^2.
struct CubicWell {
  double omega0, g;
  double vb() const { return std::pow(omega0, 6) / (6.0 * g * g); }
  double qb() const { return omega0 * omega0 / g; }
  BatchPotential batch() const {
    return [this](const std::vector<VectorXd> &q, std::vector<double> &v,
                  std::vector<VectorXd> &grad) {
      v.resize(q.size());
      grad.resize(q.size());
      for (size_t j = 0; j < q.size(); ++j) {
        const double x = q[j](0);
        v[j] = 0.5 * omega0 * omega0 * x * x - g * x * x * x / 3.0;
        grad[j] = VectorXd::Constant(1, omega0 * omega0 * x - g * x * x);
      }
    };
  }
  MatrixXd hessian(const VectorXd &q) const {
    return MatrixXd::Constant(1, 1, omega0 * omega0 - 2.0 * g * q(0));
  }
};

} // namespace

// Deep below the crossover the ring-polymer instanton rate approaches the
// zero-temperature decay rate of the cubic well (Caldeira and Leggett, Ann.
// Phys. 149, 374 (1983)),
//   Gamma = (omega0 / 2 pi) sqrt(864 pi V_b / (hbar omega0))
//           exp(-36 V_b / (5 hbar omega0)),
// whose leading semiclassical correction is of order hbar omega0 / V_b.
namespace {

using eonc::testing::Eckart;

// The ring Hessian of random bead blocks, dense, for checking the chain.
MatrixXd denseRing(const std::vector<MatrixXd> &h, double c) {
  const long n = static_cast<long>(h.size()), f = h.front().rows();
  MatrixXd j = MatrixXd::Zero(n * f, n * f);
  for (long b = 0; b < n; ++b) {
    j.block(b * f, b * f, f, f) =
        h[static_cast<size_t>(b)] + 2.0 * c * MatrixXd::Identity(f, f);
    const long k = (b + 1) % n;
    j.block(b * f, k * f, f, f) -= c * MatrixXd::Identity(f, f);
    j.block(k * f, b * f, f, f) -= c * MatrixXd::Identity(f, f);
  }
  return j;
}

} // namespace

TEST_CASE("The ring spectrum from the block chain matches the dense Hessian",
          "[Tunneling][Instanton]") {
  const long n = 9, f = 4;
  const double c = 2.7;
  std::vector<MatrixXd> h;
  std::vector<VectorXd> tau;
  unsigned seed = 12345u;
  auto rnd = [&]() {
    seed = 1664525u * seed + 1013904223u;
    return static_cast<double>(seed >> 8) / static_cast<double>(1u << 24) - 0.5;
  };
  for (long b = 0; b < n; ++b) {
    MatrixXd a(f, f);
    for (long i = 0; i < f; ++i) {
      for (long k = 0; k < f; ++k) {
        a(i, k) = rnd();
      }
    }
    a = 0.5 * (a + a.transpose()).eval();
    if (b == 0) {
      a -= 6.0 * MatrixXd::Identity(f, f);
    }
    h.push_back(a);
    VectorXd t(f);
    for (long i = 0; i < f; ++i) {
      t(i) = rnd();
    }
    tau.push_back(t);
  }
  double norm = 0.0;
  for (const auto &t : tau) {
    norm += t.squaredNorm();
  }
  for (auto &t : tau) {
    t /= std::sqrt(norm);
  }
  MatrixXd dense = denseRing(h, c);
  VectorXd tauFlat(n * f);
  for (long b = 0; b < n; ++b) {
    tauFlat.segment(b * f, f) = tau[static_cast<size_t>(b)];
  }
  const MatrixXd primed = dense + tauFlat * tauFlat.transpose();
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(primed,
                                                   Eigen::EigenvaluesOnly);
  double logDet = 0.0;
  long negative = 0;
  for (long i = 0; i < es.eigenvalues().size(); ++i) {
    logDet += std::log(std::abs(es.eigenvalues()(i)));
    negative += es.eigenvalues()(i) < 0.0 ? 1 : 0;
  }
  const RingSpectrum spec = ringSpectrum(h, c, tau);
  REQUIRE_THAT(spec.logDetPrime, WithinRel(logDet, 1e-9));
  REQUIRE(spec.negativeModes == negative);
  REQUIRE_THAT(spec.zeroEigenvalue,
               WithinRel(tauFlat.dot(dense * tauFlat), 1e-9));
}

// The instanton flux through the Eckart barrier against the exact quantum
// flux. In one dimension the instanton is the steepest-descent evaluation of
// the WKB thermal integral, so its N -> infinity limit shares the uniform
// WKB error, 7 percent low for this barrier at both temperatures, and the
// Kemble integral along the path gives the same number.
TEST_CASE("The Eckart rate instanton matches the exact flux to its "
          "semiclassical error",
          "[Tunneling][Instanton]") {
  const Eckart pes;
  const VectorXd saddle = VectorXd::Zero(1);
  const MatrixXd hs = pes.hessian_at_top();
  const double tc = crossoverTemperature(hs);
  REQUIRE_THAT(tc, WithinRel(149.988, 1e-3));
  for (const double frac : {0.5, 0.35}) {
    const double beta = 1.0 / (kBoltzmann * frac * tc);
    const double exact = pes.logExactFlux(beta);
    RateInstantonOptions opt;
    opt.forceTolerance = 1e-8;
    std::vector<double> logs;
    std::vector<long> iterations;
    for (const long n : {64L, 128L}) {
      opt.beads = n;
      RateInstanton inst =
          optimizeRateInstanton(saddle, hs, beta, {}, pes.batch(), opt);
      // Below 0.75 T_c an empty start cools from 0.85 T_c in steps of 0.75,
      // and the iteration count adds over the stages: at most 12 each.
      long stages = 1;
      for (double t = 0.85 * tc; t > frac * tc * 1.05; t *= 0.75) {
        ++stages;
      }
      CAPTURE(frac, n, inst.iterations, stages);
      REQUIRE(inst.converged);
      REQUIRE(inst.iterations <= 12 * stages);
      instantonRate(
          inst, [&](long, const VectorXd &q) { return pes.hessian_at(q); },
          MatrixXd::Identity(1, 1), 0.0);
      REQUIRE(inst.negativeModes == 1);
      REQUIRE(std::abs(inst.zeroEigenvalue) < 1e-3);
      logs.push_back(inst.logRateTimesZr);
      iterations.push_back(inst.iterations);
    }
    // 1 / N^2 extrapolation from 64 and 128 beads.
    const double extrapolated = (4.0 * logs[1] - logs[0]) / 3.0;
    CAPTURE(frac, logs[0], logs[1], extrapolated, exact);
    REQUIRE_THAT(std::exp(logs[1] - exact),
                 Catch::Matchers::WithinAbs(0.935, 0.01));
    REQUIRE_THAT(std::exp(extrapolated - exact),
                 Catch::Matchers::WithinAbs(0.928, 0.01));

    // The Kemble WKB flux along the path shares the semiclassical error.
    std::vector<VectorXd> path;
    std::vector<double> energies, s;
    for (int k = 0; k <= 80; ++k) {
      const double x = -3.0 * pes.a + 6.0 * pes.a * k / 80;
      path.push_back(VectorXd::Constant(1, x));
      energies.push_back(pes.value(x));
      s.push_back(x + 3.0 * pes.a);
    }
    const Profile profile(s, energies);
    // The path integral references energies to the path's reactant end,
    // V(-3a) = 4.2 meV here, while the exact flux counts from V = 0 at
    // infinity; the Boltzmann factor of that offset moves the comparison.
    const double hw = 0.1;
    const double wkb = wkbLogRateAlongPath(profile, beta, hw) -
                       std::log(2.0 * std::sinh(0.5 * beta * hw)) -
                       beta * energies.front();
    CAPTURE(wkb);
    REQUIRE_THAT(std::exp(wkb - exact),
                 Catch::Matchers::WithinAbs(0.928, 0.01));

    // A ring seeded from the path by the period condition converges to the
    // same instanton, in no more steps than the cosine seed.
    opt.beads = 64;
    const std::vector<VectorXd> seed =
        ringFromPath(path, energies, beta * kHbar, 64);
    REQUIRE(seed.size() == 64);
    REQUIRE_THAT(seed[0](0), Catch::Matchers::WithinAbs(-seed[32](0), 1e-6));
    RateInstanton seeded =
        optimizeRateInstanton(saddle, hs, beta, seed, pes.batch(), opt);
    REQUIRE(seeded.converged);
    REQUIRE(seeded.iterations <= iterations[0]);
    instantonRate(
        seeded, [&](long, const VectorXd &q) { return pes.hessian_at(q); },
        MatrixXd::Identity(1, 1), 0.0);
    REQUIRE_THAT(seeded.logRateTimesZr,
                 Catch::Matchers::WithinAbs(logs[0], 1e-6));
  }
}

// A free transverse coordinate is a rigid mode: its centroid factor leaves
// both the ring and the reactant, and its k > 0 factors cancel, so the rate
// is the one-dimensional one.
TEST_CASE("A rigid mode leaves the instanton rate unchanged",
          "[Tunneling][Instanton]") {
  const double omega0 = 1.0;
  const double hw = kHbar * omega0;
  const double vb = 16.0 * hw;
  const CubicWell pes{omega0, std::sqrt(std::pow(omega0, 6) / (6.0 * vb))};
  const double beta = 60.0 / hw;
  RateInstantonOptions opt;
  opt.beads = 64;
  opt.forceTolerance = 1e-8;
  const VectorXd saddle1 = VectorXd::Constant(1, pes.qb());
  RateInstanton one = optimizeRateInstanton(saddle1, pes.hessian(saddle1), beta,
                                            {}, pes.batch(), opt);
  REQUIRE(one.converged);
  instantonRate(
      one, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(VectorXd::Zero(1)), 0.0);

  BatchPotential two = [&](const std::vector<VectorXd> &q,
                           std::vector<double> &v,
                           std::vector<VectorXd> &grad) {
    std::vector<VectorXd> x(q.size());
    for (size_t j = 0; j < q.size(); ++j) {
      x[j] = q[j].head(1);
    }
    std::vector<VectorXd> g1;
    pes.batch()(x, v, g1);
    grad.resize(q.size());
    for (size_t j = 0; j < q.size(); ++j) {
      grad[j] = VectorXd::Zero(2);
      grad[j](0) = g1[j](0);
    }
  };
  auto hess2 = [&](const VectorXd &q) {
    MatrixXd h = MatrixXd::Zero(2, 2);
    h(0, 0) = pes.hessian(q.head(1))(0, 0);
    return h;
  };
  VectorXd saddle2 = VectorXd::Zero(2);
  saddle2(0) = pes.qb();
  RateInstanton both =
      optimizeRateInstanton(saddle2, hess2(saddle2), beta, {}, two, opt);
  REQUIRE(both.converged);
  instantonRate(
      both, [&](long, const VectorXd &q) { return hess2(q); },
      hess2(VectorXd::Zero(2)), 0.0, MatrixXd(), 0.0, 1);
  REQUIRE(both.negativeModes == 1);
  REQUIRE_THAT(both.logRate, Catch::Matchers::WithinAbs(one.logRate, 1e-6));
}

TEST_CASE("The rate instanton of a cubic well matches its decay rate",
          "[Tunneling][Instanton]") {
  const double omega0 = 1.0;
  const double hw = kHbar * omega0;
  const double vb = 16.0 * hw;
  const CubicWell pes{omega0, std::sqrt(std::pow(omega0, 6) / (6.0 * vb))};
  REQUIRE_THAT(pes.vb(), WithinRel(vb, 1e-12));
  const VectorXd saddle = VectorXd::Constant(1, pes.qb());
  const MatrixXd hs = pes.hessian(saddle);
  const double tc = crossoverTemperature(hs);
  REQUIRE_THAT(tc,
               WithinRel(hw / (2.0 * std::numbers::pi * kBoltzmann), 1e-12));

  const double beta = 60.0 / hw; // T about T_c / 10
  RateInstantonOptions opt;
  opt.beads = 300;
  opt.forceTolerance = 1e-8;
  opt.maxStep = 0.05;
  opt.halfRing = false;
  RateInstanton inst =
      optimizeRateInstanton(saddle, hs, beta, {}, pes.batch(), opt);
  CAPTURE(inst.iterations, inst.converged, inst.temperature, tc);
  REQUIRE(inst.converged);
  instantonRate(
      inst, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(VectorXd::Zero(1)), 0.0, hs, vb);
  REQUIRE(inst.negativeModes == 1);
  REQUIRE(std::abs(inst.zeroEigenvalue) < 1e-3 * omega0 * omega0);

  const double gammaLog = std::log(omega0 / (2.0 * std::numbers::pi)) +
                          0.5 * std::log(864.0 * std::numbers::pi * vb / hw) -
                          36.0 * vb / (5.0 * hw);
  CAPTURE(inst.logRate, gammaLog, inst.iterations, inst.temperature, tc);
  // 300 beads at beta hbar omega0 = 60 is within 15 percent.
  REQUIRE(std::abs(std::exp(inst.logRate - gammaLog) - 1.0) < 0.15);
  // Tunnelling beats the classical rate by many orders at T_c / 10.
  REQUIRE(inst.logRate > inst.classicalLogRate + std::log(1e10));

  long fullCalls = 0;
  long halfCalls = 0;
  auto counted = [&](long &calls) {
    return [&](const std::vector<VectorXd> &q, std::vector<double> &v,
               std::vector<VectorXd> &g) {
      calls += static_cast<long>(q.size());
      pes.batch()(q, v, g);
    };
  };
  RateInstantonOptions fullOpt = opt;
  fullOpt.halfRing = false;
  fullOpt.forceTolerance = 1e-8;
  RateInstanton fullRing =
      optimizeRateInstanton(saddle, hs, beta, {}, counted(fullCalls), fullOpt);
  instantonRate(
      fullRing, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(VectorXd::Zero(1)), 0.0, hs, vb);
  RateInstantonOptions halfOpt = opt;
  halfOpt.halfRing = true;
  halfOpt.forceTolerance = 1e-8;
  RateInstanton halfRing =
      optimizeRateInstanton(saddle, hs, beta, {}, counted(halfCalls), halfOpt);
  instantonRate(
      halfRing, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(VectorXd::Zero(1)), 0.0, hs, vb);
  CAPTURE(fullRing.logRate, halfRing.logRate, fullCalls, halfCalls,
          halfRing.converged, halfRing.iterations, fullRing.iterations,
          fullRing.ringPotential, halfRing.ringPotential, fullRing.bN,
          halfRing.bN);
  REQUIRE(halfRing.converged);
  REQUIRE(std::abs(halfRing.logRate - fullRing.logRate) < 1e-6);
  REQUIRE(halfCalls < fullCalls);

  // The same ring through the block determinant, which a large system uses.
  {
    const double logDense = inst.logRate;
    instantonRate(
        inst, [&](long, const VectorXd &q) { return pes.hessian(q); },
        pes.hessian(VectorXd::Zero(1)), 0.0, hs, vb, 0, 0);
    CAPTURE(inst.logRate, logDense, inst.zeroEigenvalue, inst.negativeModes);
    REQUIRE(inst.negativeModes == 1);
    REQUIRE(std::abs(inst.logRate - logDense) < 1e-5);
    inst.logRate = logDense;
  }

  // The discretisation error falls as 1 / N^2, so doubling the beads and
  // extrapolating removes it; what is left is the semiclassical error of the
  // instanton against Gamma.
  opt.beads = 600;
  RateInstanton fine =
      optimizeRateInstanton(saddle, hs, beta, {}, pes.batch(), opt);
  REQUIRE(fine.converged);
  instantonRate(
      fine, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(VectorXd::Zero(1)), 0.0, hs, vb);
  REQUIRE(fine.negativeModes == 1);
  const double extrapolated = (4.0 * fine.logRate - inst.logRate) / 3.0;
  CAPTURE(fine.logRate, extrapolated);
  REQUIRE(std::abs(fine.logRate - gammaLog) <
          std::abs(inst.logRate - gammaLog));
  REQUIRE(std::abs(std::exp(extrapolated - gammaLog) - 1.0) < 0.01);
}

TEST_CASE("Minimum-mode following copies one half of an even ring",
          "[Tunneling][Instanton]") {
  const double omega0 = 1.0;
  const double hw = kHbar * omega0;
  const double vb = 16.0 * hw;
  const CubicWell pes{omega0, std::sqrt(std::pow(omega0, 6) / (6.0 * vb))};
  const VectorXd saddle = VectorXd::Constant(1, pes.qb());
  const MatrixXd hs = pes.hessian(saddle);
  const double beta = 60.0 / hw;
  long fullCalls = 0;
  long halfCalls = 0;
  auto counted = [&](long &calls) {
    return [&](const std::vector<VectorXd> &q, std::vector<double> &v,
               std::vector<VectorXd> &g) {
      calls += static_cast<long>(q.size());
      pes.batch()(q, v, g);
    };
  };
  RateInstantonOptions fullOpt;
  fullOpt.beads = 64;
  fullOpt.forceTolerance = 1e-8;
  fullOpt.maxStep = 0.05;
  fullOpt.halfRing = false;
  fullOpt.newtonLimit = 0;
  RateInstanton fullRing =
      optimizeRateInstanton(saddle, hs, beta, {}, counted(fullCalls), fullOpt);
  RateInstantonOptions halfOpt = fullOpt;
  halfOpt.halfRing = true;
  RateInstanton halfRing =
      optimizeRateInstanton(saddle, hs, beta, {}, counted(halfCalls), halfOpt);
  CAPTURE(fullRing.logRate, halfRing.logRate, fullCalls, halfCalls,
          fullRing.converged, halfRing.converged, fullRing.iterations,
          halfRing.iterations, fullRing.ringPotential, halfRing.ringPotential);
  REQUIRE(fullRing.converged);
  REQUIRE(halfRing.converged);
  instantonRate(
      fullRing, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(VectorXd::Zero(1)), 0.0, hs, vb);
  instantonRate(
      halfRing, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(VectorXd::Zero(1)), 0.0, hs, vb);
  CAPTURE(fullRing.logRate, halfRing.logRate, fullRing.negativeModes,
          halfRing.negativeModes, fullRing.ringPotential,
          halfRing.ringPotential);
  REQUIRE(fullRing.negativeModes == 1);
  REQUIRE(halfRing.negativeModes == 1);
  REQUIRE(std::abs(halfRing.logRate - fullRing.logRate) < 1e-6);
  REQUIRE(halfCalls < fullCalls);
}

TEST_CASE("The rate instanton refuses a temperature above the crossover",
          "[Tunneling][Instanton]") {
  const double hw = kHbar;
  const CubicWell pes{1.0, std::sqrt(1.0 / (6.0 * 16.0 * hw))};
  const VectorXd saddle = VectorXd::Constant(1, pes.qb());
  const MatrixXd hs = pes.hessian(saddle);
  const double tc = crossoverTemperature(hs);
  REQUIRE_THROWS_AS(optimizeRateInstanton(saddle, hs,
                                          1.0 / (kBoltzmann * 1.5 * tc), {},
                                          pes.batch(), RateInstantonOptions{}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(crossoverTemperature(pes.hessian(VectorXd::Zero(1))),
                    std::invalid_argument);
}

TEST_CASE("The parabolic factor multiplies harmonic TST above the crossover",
          "[Tunneling][Instanton]") {
  const double tc = 100.0;
  const double hot = 1.0e6;
  const double x = std::numbers::pi * tc / 200.0;
  REQUIRE(parabolicFactor(200.0, tc) == Catch::Approx(x / std::sin(x)));
  REQUIRE(parabolicFactor(hot, tc) == Catch::Approx(1.0).margin(1e-4));
  REQUIRE_THROWS_AS(parabolicFactor(tc, tc), std::invalid_argument);
  REQUIRE_THROWS_AS(parabolicFactor(0.5 * tc, tc), std::invalid_argument);

  MatrixXd reactant(1, 1);
  MatrixXd saddle(1, 1);
  reactant(0, 0) = 1.0;
  saddle(0, 0) = -1.0;
  const double crossover = crossoverTemperature(saddle);
  const double temperature = 4.0 * crossover;
  const double beta = 1.0 / (kBoltzmann * temperature);
  const double barrier = 0.1;
  const double logRate = harmonicTstLogRate(reactant, saddle, beta, barrier, 0);
  REQUIRE(logRate ==
          Catch::Approx(-std::log(2.0 * std::numbers::pi) - beta * barrier));
  const double factor = parabolicFactor(temperature, crossover);
  const double phase = std::numbers::pi / 4.0;
  REQUIRE(factor == Catch::Approx(phase / std::sin(phase)));
  REQUIRE(factor > 1.0);
}

namespace {

MatrixXd toyBead(long j, long f, double shift) {
  MatrixXd h(f, f);
  for (long a = 0; a < f; ++a) {
    for (long b = 0; b < f; ++b) {
      h(a, b) = std::sin(0.37 * (a + 1) + 0.17 * j) *
                std::cos(0.23 * (b + 1) + 0.11 * j);
    }
  }
  h = (0.5 * (h + h.transpose())).eval();
  h.diagonal().array() += shift;
  return h;
}

double denseRingLog(double c, const std::vector<MatrixXd> &diag) {
  const long f = diag.front().rows();
  const long n = static_cast<long>(diag.size());
  MatrixXd big = MatrixXd::Zero(n * f, n * f);
  const MatrixXd eye = MatrixXd::Identity(f, f);
  for (long j = 0; j < n; ++j) {
    big.block(j * f, j * f, f, f) = diag[static_cast<size_t>(j)];
    const long k = (j + 1) % n;
    big.block(j * f, k * f, f, f) -= c * eye;
    big.block(k * f, j * f, f, f) -= c * eye;
  }
  const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
      ring = big;
  const Eigen::SelfAdjointEigenSolver<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
      es(ring, Eigen::EigenvaluesOnly);
  double sum = 0.0;
  for (long i = 0; i < es.eigenvalues().size(); ++i) {
    sum += std::log(std::abs(es.eigenvalues()(i)));
  }
  return sum;
}

} // namespace

TEST_CASE("The cyclic ring determinant matches a dense factorisation",
          "[Tunneling][Instanton]") {
  const double c = 1.3;
  const long n = 7;
  const long f = 3;
  std::vector<MatrixXd> diag(static_cast<size_t>(n));
  for (long j = 0; j < n; ++j) {
    diag[static_cast<size_t>(j)] =
        toyBead(j, f, 5.0) + 2.0 * c * MatrixXd::Identity(f, f);
  }
  REQUIRE(std::abs(cyclicRingLogAbsDet(c, diag) - denseRingLog(c, diag)) <
          1e-8);
  {
    const long fSolve = diag.front().rows();
    std::vector<VectorXd> rhs(static_cast<size_t>(n));
    VectorXd rhsStack(n * fSolve);
    MatrixXd big = MatrixXd::Zero(n * fSolve, n * fSolve);
    const MatrixXd eye = MatrixXd::Identity(fSolve, fSolve);
    for (long j = 0; j < n; ++j) {
      rhs[static_cast<size_t>(j)] =
          VectorXd::LinSpaced(fSolve, 0.2 * static_cast<double>(j),
                              0.2 * static_cast<double>(j) + fSolve);
      rhsStack.segment(j * fSolve, fSolve) = rhs[static_cast<size_t>(j)];
      big.block(j * fSolve, j * fSolve, fSolve, fSolve) =
          diag[static_cast<size_t>(j)];
      const long k = (j + 1) % n;
      big.block(j * fSolve, k * fSolve, fSolve, fSolve) -= c * eye;
      big.block(k * fSolve, j * fSolve, fSolve, fSolve) -= c * eye;
    }
    const std::vector<VectorXd> sol = cyclicRingSolve(c, diag, rhs);
    VectorXd stacked(n * fSolve);
    for (long j = 0; j < n; ++j) {
      stacked.segment(j * fSolve, fSolve) = sol[static_cast<size_t>(j)];
    }
    REQUIRE((big * stacked - rhsStack).norm() < 1e-8 * rhsStack.norm());
  }

  for (long j = 0; j < n; ++j) {
    diag[static_cast<size_t>(j)] =
        toyBead(j, f, -1.0) + 2.0 * c * MatrixXd::Identity(f, f);
  }
  REQUIRE(std::abs(cyclicRingLogAbsDet(c, diag) - denseRingLog(c, diag)) <
          1e-8);

  // One flat direction. The closed product keeps the spring eigenvalues and
  // drops the constant mode: their log sum is 2 log N + (N - 1) log c.
  const long nf = 2;
  const long nn = 8;
  const double cf = 1.7;
  std::vector<MatrixXd> flat(static_cast<size_t>(nn));
  std::vector<MatrixXd> reduced(static_cast<size_t>(nn));
  for (long j = 0; j < nn; ++j) {
    MatrixXd h = MatrixXd::Zero(nf, nf);
    h(0, 0) = 0.8 + 0.05 * static_cast<double>(j);
    flat[static_cast<size_t>(j)] = h + 2.0 * cf * MatrixXd::Identity(nf, nf);
    reduced[static_cast<size_t>(j)] =
        MatrixXd::Constant(1, 1, h(0, 0) + 2.0 * cf);
  }
  const double spring = 2.0 * std::log(static_cast<double>(nn)) +
                        static_cast<double>(nn - 1) * std::log(cf);
  MatrixXd big = MatrixXd::Zero(nn * nf, nn * nf);
  const MatrixXd eye = MatrixXd::Identity(nf, nf);
  for (long j = 0; j < nn; ++j) {
    big.block(j * nf, j * nf, nf, nf) = flat[static_cast<size_t>(j)];
    const long k = (j + 1) % nn;
    big.block(j * nf, k * nf, nf, nf) -= cf * eye;
    big.block(k * nf, j * nf, nf, nf) -= cf * eye;
  }
  const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
      ring = big;
  const Eigen::SelfAdjointEigenSolver<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
      es(ring, Eigen::EigenvaluesOnly);
  double kept = 0.0;
  long zeros = 0;
  for (long i = 0; i < es.eigenvalues().size(); ++i) {
    if (std::abs(es.eigenvalues()(i)) < 1e-8) {
      ++zeros;
      continue;
    }
    kept += std::log(std::abs(es.eigenvalues()(i)));
  }
  REQUIRE(zeros == 1);
  REQUIRE(std::abs(cyclicRingLogAbsDet(cf, reduced) + spring - kept) < 1e-8);
}

TEST_CASE("A flat coordinate cancels in the instanton rate",
          "[Tunneling][Instanton]") {
  const double omega0 = 1.0;
  const double hw = kHbar * omega0;
  const double vb = 16.0 * hw;
  const CubicWell pes{omega0, std::sqrt(std::pow(omega0, 6) / (6.0 * vb))};
  const long nBeads = 16;
  const double beta = 10.0;
  const double amp = 0.4 * pes.qb();
  auto rateOf = [&](long dim, long rigid) {
    RateInstanton inst;
    inst.beta = beta;
    inst.betaN = beta / static_cast<double>(nBeads);
    inst.beads.resize(static_cast<size_t>(nBeads));
    inst.bN = 0.0;
    inst.ringPotential = 0.0;
    const double c = 1.0 / std::pow(inst.betaN * kHbar, 2);
    for (long j = 0; j < nBeads; ++j) {
      VectorXd q = VectorXd::Zero(dim);
      q(0) = pes.qb() +
             amp * std::cos(2.0 * std::numbers::pi * static_cast<double>(j) /
                            static_cast<double>(nBeads));
      inst.beads[static_cast<size_t>(j)] = q;
    }
    for (long j = 0; j < nBeads; ++j) {
      const double x = inst.beads[static_cast<size_t>(j)](0);
      inst.ringPotential +=
          0.5 * omega0 * omega0 * x * x - pes.g * x * x * x / 3.0;
      const VectorXd step = inst.beads[static_cast<size_t>((j + 1) % nBeads)] -
                            inst.beads[static_cast<size_t>(j)];
      inst.bN += step.squaredNorm();
    }
    inst.ringPotential += 0.5 * c * inst.bN;
    MatrixXd reactant = MatrixXd::Zero(dim, dim);
    reactant(0, 0) = omega0 * omega0;
    MatrixXd saddle = MatrixXd::Zero(dim, dim);
    saddle(0, 0) = omega0 * omega0 - 2.0 * pes.g * pes.qb();
    instantonRate(
        inst,
        [&](long, const VectorXd &q) {
          MatrixXd h = MatrixXd::Zero(dim, dim);
          h(0, 0) = omega0 * omega0 - 2.0 * pes.g * q(0);
          return h;
        },
        reactant, 0.0, saddle, vb, rigid);
    return inst;
  };
  const RateInstanton line = rateOf(1, 0);
  const RateInstanton flat = rateOf(2, 1);
  CAPTURE(line.logRate, flat.logRate, line.classicalLogRate,
          flat.classicalLogRate);
  REQUIRE(std::isfinite(line.logRate));
  REQUIRE(std::isfinite(flat.logRate));
  REQUIRE(std::abs(line.logRate - flat.logRate) < 1e-6);
  REQUIRE(std::abs(line.classicalLogRate - flat.classicalLogRate) < 1e-6);
}

// Quantum harmonic TST against its closed form, and its classical limit:
// one bound mode of curvature 4 at the reactant, one bound mode of
// curvature 1 and the barrier mode at the saddle.
TEST_CASE("Quantum harmonic TST has its closed form and the classical limit",
          "[Tunneling][Instanton]") {
  MatrixXd hr = MatrixXd::Zero(2, 2);
  hr(0, 0) = 4.0;
  hr(1, 1) = 9.0;
  MatrixXd hs = MatrixXd::Zero(2, 2);
  hs(0, 0) = -2.0;
  hs(1, 1) = 1.0;
  const double barrier = 0.3;
  for (const double beta : {5.0, 40.0}) {
    const double bh = beta * kHbar;
    auto twoSinh = [&](double w) { return 2.0 * std::sinh(0.5 * bh * w); };
    const double expected =
        std::log(twoSinh(2.0) * twoSinh(3.0) / twoSinh(1.0)) -
        std::log(2.0 * std::numbers::pi * bh) - beta * barrier;
    REQUIRE_THAT(quantumHarmonicTstLogRate(hr, hs, beta, barrier, 0),
                 WithinRel(expected, 1e-12));
  }
  // At high temperature every 2 sinh(x / 2) -> x and the classical rate
  // follows.
  const double hot = 1e-4;
  REQUIRE_THAT(quantumHarmonicTstLogRate(hr, hs, hot, barrier, 0) -
                   harmonicTstLogRate(hr, hs, hot, barrier, 0),
               Catch::Matchers::WithinAbs(0.0, 1e-6));
}

// A soft barrier: the saddle's unstable curvature, -1e-8, lies closer to
// zero than the finite-difference residue of its rigid mode, 1e-6. The
// rigid mode is the one to leave, and the barrier mode leaves as the
// unstable mode, so only the bound curvature 2 stays at the saddle.
TEST_CASE("Harmonic TST never takes a soft barrier mode for a rigid mode",
          "[Tunneling][Instanton]") {
  MatrixXd hr = MatrixXd::Zero(3, 3);
  hr(0, 0) = 2e-6;
  hr(1, 1) = 4.0;
  hr(2, 2) = 9.0;
  MatrixXd hs = MatrixXd::Zero(3, 3);
  hs(0, 0) = -1e-8;
  hs(1, 1) = 1e-6;
  hs(2, 2) = 2.0;
  const double barrier = 0.3;
  const double beta = 5.0;
  const double classical = 0.5 * std::log(4.0 * 9.0 / 2.0) -
                           std::log(2.0 * std::numbers::pi) - beta * barrier;
  REQUIRE_THAT(harmonicTstLogRate(hr, hs, beta, barrier, 1),
               WithinRel(classical, 1e-12));
  const double bh = beta * kHbar;
  auto twoSinh = [&](double w) { return 2.0 * std::sinh(0.5 * bh * w); };
  const double quantum =
      std::log(twoSinh(2.0) * twoSinh(3.0) / twoSinh(std::sqrt(2.0))) -
      std::log(2.0 * std::numbers::pi * bh) - beta * barrier;
  REQUIRE_THAT(quantumHarmonicTstLogRate(hr, hs, beta, barrier, 1),
               WithinRel(quantum, 1e-12));
}

// Two beads in the wells given a curvature far below the spring's make two
// negative modes of the action Hessian. The determinant keeps its sign, so
// only the inertia shows that the path is not a minimum of the action.
TEST_CASE("The instanton splitting refuses a path with two negative modes",
          "[Tunneling][Instanton]") {
  const CurvedValley pes{0.3, 4.0, 0.0};
  Instanton inst = valleyInstanton(pes, 128, 40.0);
  REQUIRE(inst.converged);
  REQUIRE(inst.delta0 > 0.0);
  const double big = 4.0 / (inst.dtau * inst.dtau) + 100.0;
  auto twoDips = [&](long j, const VectorXd &q) -> MatrixXd {
    MatrixXd h = pes.hessian(q);
    if (j == 20 || j == 108) {
      h(1, 1) -= big;
    }
    return h;
  };
  VectorXd a(2), b(2);
  a << -1.0, 0.0;
  b << 1.0, 0.0;
  REQUIRE_THROWS_AS(
      instantonSplitting(inst, twoDips, pes.hessian(a), pes.hessian(b)),
      std::runtime_error);
  // One such bead flips the determinant's sign and was refused before too.
  auto oneDip = [&](long j, const VectorXd &q) -> MatrixXd {
    MatrixXd h = pes.hessian(q);
    if (j == 20) {
      h(1, 1) -= big;
    }
    return h;
  };
  REQUIRE_THROWS_AS(
      instantonSplitting(inst, oneDip, pes.hessian(a), pes.hessian(b)),
      std::runtime_error);
}

namespace {

// Two unit-mass atoms whose bond d = r - r_e decays from a cubic well,
// V = K d^2 / 2 - G d^3 / 3: a free diatomic with three translations and
// two rotations. Unit masses make q the Cartesian displacement from
// `reference`.
struct CubicBond {
  double k = 1.0, g = std::sqrt(1.0 / 3.0), re = 2.0;
  VectorXd reference() const {
    VectorXd r = VectorXd::Zero(6);
    r(3) = re;
    return r;
  }
  double vd(double d) const { return 0.5 * k * d * d - g * d * d * d / 3.0; }
  double dvd(double d) const { return k * d - g * d * d; }
  double d2vd(double d) const { return k - 2.0 * g * d; }
  Eigen::Vector3d bond(const VectorXd &q) const {
    const VectorXd x = reference() + q;
    return x.segment<3>(3) - x.segment<3>(0);
  }
  double value(const VectorXd &q) const { return vd(bond(q).norm() - re); }
  VectorXd gradient(const VectorXd &q) const {
    const Eigen::Vector3d b = bond(q);
    const double r = b.norm();
    const Eigen::Vector3d f = dvd(r - re) * b / r;
    VectorXd out(6);
    out << -f, f;
    return out;
  }
  MatrixXd hessian(const VectorXd &q) const {
    const Eigen::Vector3d b = bond(q);
    const double r = b.norm();
    const Eigen::Vector3d u = b / r;
    const Eigen::Matrix3d a =
        d2vd(r - re) * u * u.transpose() +
        dvd(r - re) / r * (Eigen::Matrix3d::Identity() - u * u.transpose());
    MatrixXd h(6, 6);
    h << a, -a, -a, a;
    return h;
  }
  BatchPotential batch() const {
    return [this](const std::vector<VectorXd> &q, std::vector<double> &v,
                  std::vector<VectorXd> &grad) {
      v.resize(q.size());
      grad.resize(q.size());
      for (size_t j = 0; j < q.size(); ++j) {
        v[j] = value(q[j]);
        grad[j] = gradient(q[j]);
      }
    };
  }
};

} // namespace

// The rate's det' leaves out the ring's null vectors. For a rotation those
// are e x (r_j - centre) at each bead, which a stretched bond makes differ
// from bead to bead; the reactant's generator copied to every bead is not
// one of them. The reference is the dense ring Hessian's spectrum with its
// six eigenvalues nearest zero (the cycle, three translations, two
// rotations) left out.
TEST_CASE("The rate lifts the rotations of a free diatomic's ring",
          "[Tunneling][Instanton]") {
  const CubicBond pes;
  const double db = pes.k / pes.g;
  VectorXd saddle = VectorXd::Zero(6);
  saddle(0) = -0.5 * db;
  saddle(3) = 0.5 * db;
  const MatrixXd hs = pes.hessian(saddle);
  const MatrixXd hr = pes.hessian(VectorXd::Zero(6));
  const double tc = crossoverTemperature(hs);
  REQUIRE_THAT(tc, WithinRel(kHbar * std::sqrt(2.0 * pes.k) /
                                 (2.0 * std::numbers::pi * kBoltzmann),
                             1e-10));
  const double beta = 1.0 / (kBoltzmann * 0.5 * tc);
  RateInstantonOptions opt;
  opt.beads = 24;
  opt.forceTolerance = 1e-9;
  opt.rigidSqrtMasses = {1.0, 1.0};
  opt.rigidReference = pes.reference();
  opt.rigidRotations = {true, true, true};
  RateInstanton inst =
      optimizeRateInstanton(saddle, hs, beta, {}, pes.batch(), opt);
  CAPTURE(inst.iterations, inst.ringPotential);
  REQUIRE(inst.converged);
  // The bond stretches along the ring.
  double rMin = 1e9, rMax = 0.0;
  for (const auto &q : inst.beads) {
    rMin = std::min(rMin, pes.bond(q).norm());
    rMax = std::max(rMax, pes.bond(q).norm());
  }
  CAPTURE(rMin, rMax);
  REQUIRE(rMax - rMin > 0.5);

  const long n = opt.beads, f = 6;
  const double bnh = inst.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);
  std::vector<MatrixXd> blocks;
  for (const auto &q : inst.beads) {
    blocks.push_back(pes.hessian(q));
  }
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(denseRing(blocks, c),
                                                   Eigen::EigenvaluesOnly);
  std::vector<double> lam(es.eigenvalues().data(),
                          es.eigenvalues().data() + n * f);
  std::sort(lam.begin(), lam.end(),
            [](double a, double b) { return std::abs(a) < std::abs(b); });
  double logDetPrime = 0.0;
  long negative = 0;
  for (size_t i = 6; i < lam.size(); ++i) {
    logDetPrime += std::log(std::abs(lam[i]));
    negative += lam[i] < 0.0 ? 1 : 0;
  }
  REQUIRE(negative == 1);
  const double expected =
      -std::log(bnh) +
      0.5 * std::log(inst.bN /
                     (2.0 * std::numbers::pi * inst.betaN * kHbar * kHbar)) -
      (static_cast<double>(n * f - 6) * std::log(bnh) + 0.5 * logDetPrime) -
      inst.betaN * inst.ringPotential;

  auto hessian = [&](long, const VectorXd &q) { return pes.hessian(q); };
  RateInstanton ring = inst;
  RingRigidBodies bodies;
  bodies.sqrtMasses = {1.0, 1.0};
  bodies.reference = pes.reference();
  bodies.rotations = {true, true, true};
  instantonRate(ring, hessian, hr, 0.0, MatrixXd(), 0.0, 5, 4096, bodies);
  RateInstanton copied = inst;
  instantonRate(copied, hessian, hr, 0.0, MatrixXd(), 0.0, 5);
  CAPTURE(expected, ring.logRateTimesZr, copied.logRateTimesZr);
  REQUIRE(ring.negativeModes == 1);
  REQUIRE_THAT(ring.logRateTimesZr, Catch::Matchers::WithinAbs(expected, 1e-6));
  // The reactant's rotation generators on every bead miss the ring's null
  // space by a stretch-dependent angle.
  REQUIRE(std::abs(copied.logRateTimesZr - expected) > 1e-3);
}

// Known answers independent of eOn, from
// client/validation/instanton_references.py (mpmath, 40 digits).
//
// Symmetric Eckart V0 sech^2(x / a), V0 = 0.425 eV, a = 0.734, unit mass:
// T_c = hbar sqrt(2 V0) / (2 pi kB a) = 149.98818880973 K. In one dimension
// the N -> infinity ring-polymer instanton is the steepest-descent value of
// (1 / 2 pi hbar) int exp(-theta(E) - beta E) dE with
// theta = 2 sqrt(2) pi a (sqrt(V0) - sqrt(E)) / hbar, in closed form:
// ln(k Z_r) = -50.7985982061252 at T_c / 2 (the exact flux, -50.7235175117,
// is 7.8 percent above it). The discrete ring approaches it as 1 / N^2.
// The reactant partition function, here a harmonic well of curvature 1,
// approaches ln Z = -ln(2 sinh(u / 2)) with error u^3 coth(u / 2) / (48 N^2),
// u = beta hbar omega.
TEST_CASE("The Eckart ring converges to the analytic instanton at second "
          "order",
          "[Tunneling][Instanton]") {
  const Eckart pes;
  const MatrixXd hs = pes.hessian_at_top();
  const double tc = crossoverTemperature(hs);
  REQUIRE_THAT(tc, WithinRel(149.98818880973, 1e-11));
  const double beta = 1.0 / (kBoltzmann * 0.5 * tc);
  const double analytic = -50.7985982061252;
  const double u = beta * kHbar;
  const double logZ = -std::log(2.0 * std::sinh(0.5 * u));
  RateInstantonOptions opt;
  opt.forceTolerance = 1e-9;
  std::vector<double> err, errZ;
  const std::vector<long> sizes{16, 32, 64, 128};
  for (const long n : sizes) {
    opt.beads = n;
    RateInstanton inst = optimizeRateInstanton(VectorXd::Zero(1), hs, beta, {},
                                               pes.batch(), opt);
    REQUIRE(inst.converged);
    instantonRate(
        inst, [&](long, const VectorXd &q) { return pes.hessian_at(q); },
        MatrixXd::Identity(1, 1), 0.0);
    REQUIRE(inst.negativeModes == 1);
    // The ring belongs to the barrier: its turning points lie on either
    // side of the top, mirror images of each other.
    double lo = 0.0, hi = 0.0;
    for (const auto &q : inst.beads) {
      lo = std::min(lo, q(0));
      hi = std::max(hi, q(0));
    }
    REQUIRE(lo < 0.0);
    REQUIRE_THAT(hi, Catch::Matchers::WithinAbs(-lo, 1e-6));
    err.push_back(inst.logRateTimesZr - analytic);
    errZ.push_back(inst.logZr - logZ);
    const double leading =
        u * u * u / std::tanh(0.5 * u) / 48.0 / static_cast<double>(n * n);
    CAPTURE(n, inst.logRateTimesZr, err.back(), errZ.back(), leading);
    REQUIRE(errZ.back() > 0.0);
    if (n >= 64) {
      REQUIRE_THAT(errZ.back(), WithinRel(leading, 0.05));
    }
  }
  for (size_t i = 1; i < sizes.size(); ++i) {
    const double order = std::log2(err[i - 1] / err[i]);
    const double orderZ = std::log2(errZ[i - 1] / errZ[i]);
    CAPTURE(sizes[i], err[i], order, orderZ);
    REQUIRE(err[i] * err[i - 1] > 0.0);
    if (i + 1 == sizes.size()) {
      REQUIRE(order > 1.9);
      REQUIRE(order < 2.1);
      REQUIRE(orderZ > 1.9);
      REQUIRE(orderZ < 2.1);
    }
  }
  // Richardson extrapolation from 64 and 128 beads lands on the analytic
  // instanton.
  const double extrapolated = (4.0 * err[3] - err[2]) / 3.0;
  CAPTURE(extrapolated);
  REQUIRE(std::abs(extrapolated) < 2e-3);
}

// The quartic double well V0 (x^2 - 1)^2 at unit mass, here with a
// decoupled transverse mode that cancels between path and wells. The exact
// splitting at V0 = 0.3 eV is 1.19539175e-7 eV (sinc DVR, converged to 5e-8
// relative from 301 to 601 points over [-2.4, 2.4]); the continuum
// instanton is 2 hbar omega sqrt(4 omega / (pi hbar)) exp(-2 omega /
// (3 hbar)), omega^2 = 8 V0, = 1.2777723497e-7 eV, 6.9 percent above it,
// the semiclassical error of order hbar omega / V0 (S0 / hbar = 16.0). The
// discrete instanton converges to the continuum one as 1 / P^2 once
// omega dtau = beta hbar omega / P is small: from 64 to 128 beads
// (omega dtau 0.63 to 0.31) the observed order is still 1.43.
TEST_CASE("The instanton splitting converges to the analytic continuum "
          "value at second order",
          "[Tunneling][Instanton]") {
  const CurvedValley pes{0.3, 4.0, 0.0};
  const double continuum = 1.2777723497e-7;
  const double dvr = 1.19539175e-7;
  std::vector<double> err;
  const std::vector<long> sizes{128, 256, 512};
  for (const long p : sizes) {
    const Instanton inst = valleyInstanton(pes, p, 40.0);
    REQUIRE(inst.converged);
    REQUIRE(inst.modeSeparation > 1e3);
    err.push_back(std::log(inst.delta0 / continuum));
    CAPTURE(p, inst.delta0, err.back(), inst.delta0 / dvr);
  }
  for (size_t i = 1; i < sizes.size(); ++i) {
    const double order = std::log2(err[i - 1] / err[i]);
    CAPTURE(sizes[i], err[i], order);
    REQUIRE(err[i] * err[i - 1] > 0.0);
    REQUIRE(order > 1.8);
    REQUIRE(order < 2.2);
  }
  const double extrapolated = (4.0 * err.back() - err[err.size() - 2]) / 3.0;
  CAPTURE(extrapolated);
  REQUIRE(std::abs(extrapolated) < 1e-3);
  REQUIRE_THAT(continuum / dvr, WithinRel(1.0689151, 1e-6));
}

// A ring around the seeded saddle straddles its dividing plane close to the
// saddle. A ring of the same shape moved across the mode to a neighbouring
// saddle still straddles the plane, and that alone passed it, but it
// crosses the plane far from this saddle; a ring turned across the mode
// fails the chord test.
TEST_CASE("A ring belongs to the saddle whose dividing plane it crosses "
          "nearby",
          "[Tunneling][Instanton]") {
  VectorXd saddle = VectorXd::Zero(2);
  VectorXd mode(2);
  mode << 1.0, 0.0;
  std::vector<VectorXd> ring;
  for (int j = 0; j < 16; ++j) {
    const double t = 2.0 * std::numbers::pi * j / 16.0;
    VectorXd q(2);
    q << -0.8 * std::cos(t), 0.1 * std::sin(t) + 0.05;
    ring.push_back(q);
  }
  const RingChannel own = ringChannel(ring, saddle, mode);
  CAPTURE(own.sMin, own.sMax, own.chordOverlap, own.crossingOffset);
  REQUIRE(own.belongs);
  REQUIRE_THAT(own.sMin, Catch::Matchers::WithinAbs(-0.8, 1e-12));
  REQUIRE_THAT(own.sMax, Catch::Matchers::WithinAbs(0.8, 1e-12));
  REQUIRE(own.chordOverlap > 0.99);
  REQUIRE(own.crossingOffset < 0.2);

  std::vector<VectorXd> shifted = ring;
  for (auto &q : shifted) {
    q(1) += 3.0;
  }
  const RingChannel other = ringChannel(shifted, saddle, mode);
  CAPTURE(other.crossingOffset);
  REQUIRE(other.sMin < 0.0);
  REQUIRE(other.sMax > 0.0);
  REQUIRE_FALSE(other.belongs);

  // The same ellipse turned 70 degrees off the mode: it still straddles
  // the plane near the saddle, but its turning points run across the mode.
  std::vector<VectorXd> turned;
  const double th = 70.0 * std::numbers::pi / 180.0;
  for (int j = 0; j < 16; ++j) {
    const double t = 2.0 * std::numbers::pi * j / 16.0;
    const double x = 0.8 * std::cos(t), y = 0.1 * std::sin(t);
    VectorXd r(2);
    r << std::cos(th) * x - std::sin(th) * y,
        std::sin(th) * x + std::cos(th) * y;
    turned.push_back(r);
  }
  const RingChannel across = ringChannel(turned, saddle, mode);
  CAPTURE(across.chordOverlap, across.crossingOffset);
  REQUIRE(across.chordOverlap < 0.5);
  REQUIRE(across.crossingOffset < across.sMax - across.sMin);
  REQUIRE(across.sMin < 0.0);
  REQUIRE(across.sMax > 0.0);
  REQUIRE_FALSE(across.belongs);
}

// A soft bound mode, beta hbar omega = 1e-9 at the reactant and 2e-9 at the
// saddle: the ratio of 2 sinh(x / 2) is 1/2 to 1e-18, which ln(1 - exp(-x))
// loses to cancellation (3e-8 in each logarithm, mpmath reference in
// instanton_references.py).
TEST_CASE("Quantum harmonic TST keeps a soft mode's zero-point factor exact",
          "[Tunneling][Instanton]") {
  const double beta = 40.0;
  const double bh = beta * kHbar;
  const double xr = 1e-9, xs = 2e-9;
  MatrixXd hr = MatrixXd::Zero(2, 2);
  hr(0, 0) = std::pow(xr / bh, 2);
  hr(1, 1) = 4.0;
  MatrixXd hs = MatrixXd::Zero(2, 2);
  hs(0, 0) = -2.0;
  hs(1, 1) = std::pow(xs / bh, 2);
  const double barrier = 0.3;
  const double expected =
      std::log(xr / xs) + std::log(2.0 * std::sinh(0.5 * bh * 2.0)) -
      std::log(2.0 * std::numbers::pi * bh) - beta * barrier;
  REQUIRE_THAT(quantumHarmonicTstLogRate(hr, hs, beta, barrier, 0),
               Catch::Matchers::WithinAbs(expected, 1e-12));
}

// RA's prefactor runs over the eigenvalues of the ring Hessian with the
// near-zero one left out. On a 12-bead ring the central difference along the
// ring is not that eigenvalue's eigenvector, so lifting it leaves a factor
// cos^2 of their angle in det'; the rate must match the dense spectrum with
// its eigenvalue nearest zero removed.
TEST_CASE("The rate leaves out the ring Hessian's near-zero eigenvalue",
          "[Tunneling][Instanton]") {
  const Eckart pes;
  const MatrixXd hs = pes.hessian_at_top();
  const double tc = crossoverTemperature(hs);
  const double beta = 1.0 / (kBoltzmann * 0.5 * tc);
  RateInstantonOptions opt;
  opt.beads = 12;
  opt.forceTolerance = 1e-10;
  RateInstanton inst =
      optimizeRateInstanton(VectorXd::Zero(1), hs, beta, {}, pes.batch(), opt);
  REQUIRE(inst.converged);
  const double bnh = inst.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);
  std::vector<MatrixXd> blocks;
  for (const auto &q : inst.beads) {
    blocks.push_back(pes.hessian_at(q));
  }
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(denseRing(blocks, c),
                                                   Eigen::EigenvaluesOnly);
  std::vector<double> lam(es.eigenvalues().data(),
                          es.eigenvalues().data() + opt.beads);
  std::sort(lam.begin(), lam.end(),
            [](double a, double b) { return std::abs(a) < std::abs(b); });
  double logDetPrime = 0.0;
  for (size_t i = 1; i < lam.size(); ++i) {
    logDetPrime += std::log(std::abs(lam[i]));
  }
  const double expected =
      -std::log(bnh) +
      0.5 * std::log(inst.bN /
                     (2.0 * std::numbers::pi * inst.betaN * kHbar * kHbar)) -
      (static_cast<double>(opt.beads - 1) * std::log(bnh) + 0.5 * logDetPrime) -
      inst.betaN * inst.ringPotential;
  for (const long dense : {4096L, 0L}) {
    RateInstanton ring = inst;
    instantonRate(
        ring, [&](long, const VectorXd &q) { return pes.hessian_at(q); },
        MatrixXd::Identity(1, 1), 0.0, MatrixXd(), 0.0, 0, dense);
    CAPTURE(dense, expected, ring.logRateTimesZr, lam[0], lam[1], lam[2],
            lam[3], ring.zeroEigenvalue, ring.negativeEigenvalue,
            ring.negativeModes, c);
    REQUIRE(ring.negativeModes == 1);
    REQUIRE_THAT(ring.logRateTimesZr,
                 Catch::Matchers::WithinAbs(expected, 1e-9));
    REQUIRE_THAT(ring.zeroEigenvalue,
                 Catch::Matchers::WithinAbs(lam[0], 1e-9 * c));
  }
}

// initial_hessians = finite_difference builds every bead block from 2 f
// gradient calls before the first step; the search wrote those blocks into
// an empty vector. The ring it finds is the one the saddle-copied blocks
// find.
TEST_CASE("Finite-difference initial bead Hessians find the same ring",
          "[Tunneling][Instanton]") {
  const CubicBond pes;
  const double db = pes.k / pes.g;
  VectorXd saddle = VectorXd::Zero(6);
  saddle(0) = -0.5 * db;
  saddle(3) = 0.5 * db;
  const MatrixXd hs = pes.hessian(saddle);
  const double beta = 1.0 / (kBoltzmann * 0.5 * crossoverTemperature(hs));
  RateInstantonOptions opt;
  opt.beads = 16;
  opt.forceTolerance = 1e-9;
  opt.rigidSqrtMasses = {1.0, 1.0};
  opt.rigidReference = pes.reference();
  opt.rigidRotations = {true, true, true};
  const RateInstanton copied =
      optimizeRateInstanton(saddle, hs, beta, {}, pes.batch(), opt);
  opt.initialHessians = "finite_difference";
  const RateInstanton fd =
      optimizeRateInstanton(saddle, hs, beta, {}, pes.batch(), opt);
  CAPTURE(copied.ringPotential, fd.ringPotential, copied.iterations,
          fd.iterations);
  REQUIRE(copied.converged);
  REQUIRE(fd.converged);
  REQUIRE_THAT(fd.ringPotential,
               Catch::Matchers::WithinAbs(copied.ringPotential, 1e-8));
  REQUIRE_THAT(fd.bN, WithinRel(copied.bN, 1e-6));
}

namespace {

// The Eckart barrier along x with a transverse mode whose curvature is
// positive at the saddle and at a short ring's turning points but turns
// negative where the instanton turns: V = V0 sech^2(x / a) + k(x) y^2 / 2 +
// 5 y^4, k(x) = 1 - 2 [exp(-((x - x0) / w)^2) + exp(-((x + x0) / w)^2)].
// y = 0 is a mirror plane, so a search that starts on it never moves y.
struct TransverseDip {
  double v0 = 0.425, a = 0.734, x0 = 0.95, w = 0.15;
  double k(double x) const {
    return 1.0 - 2.0 * (std::exp(-std::pow((x - x0) / w, 2)) +
                        std::exp(-std::pow((x + x0) / w, 2)));
  }
  double dk(double x) const {
    return 4.0 *
           ((x - x0) * std::exp(-std::pow((x - x0) / w, 2)) +
            (x + x0) * std::exp(-std::pow((x + x0) / w, 2))) /
           (w * w);
  }
  VectorXd gradient(const VectorXd &q) const {
    const double x = q(0), y = q(1), ch = std::cosh(x / a);
    VectorXd g(2);
    g << -2.0 * v0 * std::tanh(x / a) / (a * ch * ch) + 0.5 * dk(x) * y * y,
        k(x) * y + 20.0 * y * y * y;
    return g;
  }
  BatchPotential batch() const {
    return [this](const std::vector<VectorXd> &q, std::vector<double> &v,
                  std::vector<VectorXd> &g) {
      v.resize(q.size());
      g.resize(q.size());
      for (size_t i = 0; i < q.size(); ++i) {
        const double x = q[i](0), y = q[i](1), ch = std::cosh(x / a);
        v[i] = v0 / (ch * ch) + 0.5 * k(x) * y * y + 5.0 * y * y * y * y;
        g[i] = gradient(q[i]);
      }
    };
  }
  MatrixXd hessian(const VectorXd &q) const {
    MatrixXd h(2, 2);
    for (int i = 0; i < 2; ++i) {
      VectorXd e = VectorXd::Zero(2);
      e(i) = 1e-6;
      h.col(i) = (gradient(q + e) - gradient(q - e)) / 2e-6;
    }
    return 0.5 * (h + h.transpose());
  }
};

} // namespace

TEST_CASE("A converged ring is classified with exact Hessians and leaves a "
          "higher-index point",
          "[Tunneling][Instanton]") {
  // On the mirror plane the steps never carry a y component, so Bofill
  // blocks keep the transverse curvature of the start, positive, where the
  // surface has turned negative. Counted on those blocks the y = 0 ring
  // has one negative curvature and passes for the instanton; exactly it
  // has three. Classified exactly and left by a trust-sized step down the
  // second mode, the search reaches the lower ring off the plane.
  const TransverseDip pes;
  const VectorXd saddle = VectorXd::Zero(2);
  const MatrixXd hs = pes.hessian(saddle);
  const double beta = 1.0 / (kBoltzmann * 0.5 * crossoverTemperature(hs));
  const long n = 32;
  std::vector<VectorXd> guess(static_cast<size_t>(n), VectorXd::Zero(2));
  for (long j = 0; j < n; ++j) {
    guess[static_cast<size_t>(j)](0) =
        0.647 * std::cos(2.0 * std::numbers::pi * j / n);
  }
  RateInstantonOptions opt;
  opt.beads = n;
  opt.forceTolerance = 1e-6;
  const RateInstanton inst =
      optimizeRateInstanton(saddle, hs, beta, guess, pes.batch(), opt);
  const double bnh = inst.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);
  MatrixXd ring = MatrixXd::Zero(2 * n, 2 * n);
  for (long j = 0; j < n; ++j) {
    const long nx = (j + 1) % n;
    ring.block(2 * j, 2 * j, 2, 2) =
        pes.hessian(inst.beads[static_cast<size_t>(j)]) +
        2.0 * c * MatrixXd::Identity(2, 2);
    ring.block(2 * j, 2 * nx, 2, 2) -= c * MatrixXd::Identity(2, 2);
    ring.block(2 * nx, 2 * j, 2, 2) -= c * MatrixXd::Identity(2, 2);
  }
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(ring,
                                                   Eigen::EigenvaluesOnly);
  long negative = 0;
  for (long i = 0; i < es.eigenvalues().size(); ++i) {
    negative += es.eigenvalues()(i) < -1e-3 ? 1 : 0;
  }
  double offPlane = 0.0;
  for (const auto &q : inst.beads) {
    offPlane = std::max(offPlane, std::abs(q(1)));
  }
  CAPTURE(inst.converged, inst.ringPotential, offPlane, negative,
          es.eigenvalues().head(4).transpose());
  REQUIRE(inst.converged);
  REQUIRE(negative == 1);
  REQUIRE(offPlane > 0.05);
  // The y = 0 ring sits at 10.1856 eV; the instanton below it.
  REQUIRE(inst.ringPotential < 10.18);
}

TEST_CASE("A converged half ring is probed for an unstable odd mode",
          "[Tunneling][Instanton]") {
  // The probe sits behind every converged half ring of the Newton search:
  // its Lanczos run costs ring evaluations, and on the Eckart barrier,
  // where the one-bounce ring is the instanton, it leaves the ring alone.
  const eonc::testing::Eckart pes;
  VectorXd saddle = VectorXd::Zero(1);
  MatrixXd hs(1, 1);
  hs(0, 0) = pes.curvature(0.0);
  const double beta = 1.0 / (kBoltzmann * 0.5 * crossoverTemperature(hs));
  auto run = [&](bool probe, long &calls) {
    const BatchPotential bare = pes.batch();
    calls = 0;
    const BatchPotential counted = [&](const std::vector<VectorXd> &q,
                                       std::vector<double> &v,
                                       std::vector<VectorXd> &g) {
      ++calls;
      bare(q, v, g);
    };
    RateInstantonOptions opt;
    opt.beads = 32;
    opt.forceTolerance = 1e-8;
    opt.checkOddSector = probe;
    return optimizeRateInstanton(saddle, hs, beta, {}, counted, opt);
  };
  long plainCalls = 0;
  long probedCalls = 0;
  const RateInstanton plain = run(false, plainCalls);
  const RateInstanton probed = run(true, probedCalls);
  CAPTURE(plainCalls, probedCalls, plain.ringPotential, probed.ringPotential);
  REQUIRE(plain.converged);
  REQUIRE(probed.converged);
  REQUIRE(probedCalls > plainCalls);
  REQUIRE_THAT(probed.ringPotential,
               Catch::Matchers::WithinAbs(plain.ringPotential, 1e-12));
}

// A ring in a harmonic well, V = q^2 / 2, searched with a saddle Hessian
// whose barrier curvature is -1: the well has no index-1 ring, the index-1
// step climbs the centroid and takes exact Newton on every internal mode,
// whose stationary point is a ring of zero size. The search must stop at
// that collapse instead of spending its budget there.
TEST_CASE("A ring that collapses onto one point stops early",
          "[Tunneling][Instanton]") {
  const BatchPotential well = [](const std::vector<VectorXd> &q,
                                 std::vector<double> &v,
                                 std::vector<VectorXd> &g) {
    v.resize(q.size());
    g.resize(q.size());
    for (size_t j = 0; j < q.size(); ++j) {
      v[j] = 0.5 * q[j].squaredNorm();
      g[j] = q[j];
    }
  };
  const VectorXd saddle = VectorXd::Zero(1);
  const MatrixXd hs = MatrixXd::Constant(1, 1, -1.0);
  const double beta = 1.0 / (kBoltzmann * 0.5 * crossoverTemperature(hs));
  RateInstantonOptions opt;
  opt.beads = 16;
  opt.maxIterations = 500;
  std::vector<VectorXd> guess;
  for (long j = 0; j < opt.beads; ++j) {
    guess.push_back(VectorXd::Constant(
        1, 0.02 * std::cos(2.0 * std::numbers::pi * static_cast<double>(j) /
                           static_cast<double>(opt.beads))));
  }
  const RateInstanton inst =
      optimizeRateInstanton(saddle, hs, beta, guess, well, opt);
  CAPTURE(inst.iterations, inst.bN, inst.converged, inst.ringPotential);
  REQUIRE(inst.collapsed);
  REQUIRE_FALSE(inst.converged);
  REQUIRE(inst.iterations < 20);
}

TEST_CASE("friction bath adds a positive term and its gradient",
          "[instanton]") {
  std::vector<VectorXd> q(4, VectorXd::Zero(1));
  q[0](0) = 0.0;
  q[1](0) = 0.2;
  q[2](0) = 0.5;
  q[3](0) = 0.1;
  std::vector<VectorXd> grad(4, VectorXd::Zero(1));
  double flat = 0.0;
  addFrictionBath(q, flat, grad, {0.0}, 1.0);
  REQUIRE(flat == 0.0);

  // One eta: (eta / 2) sum_k omega_k |Q_k|^2 over the ring's normal modes,
  // omega_k = 2 omega_P |sin(pi k / N)| (Litman et al. 2022, Eq. 20).
  const double omegaP = 2.5;
  double implicit = 0.0;
  grad.assign(4, VectorXd::Zero(1));
  addFrictionBath(q, implicit, grad, {0.25}, omegaP);
  double modes = 0.0;
  for (int k = 1; k < 4; ++k) {
    double re = 0.0, im = 0.0;
    for (int j = 0; j < 4; ++j) {
      re += q[j](0) * std::cos(2.0 * std::numbers::pi * k * j / 4.0) / 2.0;
      im -= q[j](0) * std::sin(2.0 * std::numbers::pi * k * j / 4.0) / 2.0;
    }
    modes += 0.5 * 0.25 * 2.0 * omegaP *
             std::abs(std::sin(std::numbers::pi * k / 4.0)) *
             (re * re + im * im);
  }
  REQUIRE_THAT(implicit, WithinRel(modes, 1e-12));

  // Both gradients against central differences.
  auto energy = [&](const std::vector<VectorXd> &x,
                    const std::vector<double> &eta) {
    double u = 0.0;
    std::vector<VectorXd> g(x.size(), VectorXd::Zero(1));
    addFrictionBath(x, u, g, eta, omegaP);
    return u;
  };
  for (const std::vector<double> &eta :
       {std::vector<double>{0.25}, std::vector<double>{0.1, 0.2, 0.4, 0.05}}) {
    std::vector<VectorXd> g(4, VectorXd::Zero(1));
    double u = 0.0;
    addFrictionBath(q, u, g, eta, omegaP);
    for (int j = 0; j < 4; ++j) {
      std::vector<VectorXd> up = q, down = q;
      up[j](0) += 1e-6;
      down[j](0) -= 1e-6;
      REQUIRE_THAT(g[j](0),
                   Catch::Matchers::WithinAbs(
                       (energy(up, eta) - energy(down, eta)) / 2e-6, 1e-7));
    }
  }

  // A bead-wise bath does not care which bead is called 0.
  const std::vector<double> eta{0.1, 0.2, 0.4, 0.05};
  const double explicitBath = energy(q, eta);
  REQUIRE(explicitBath > 0.0);
  REQUIRE(explicitBath != Catch::Approx(implicit));
  for (int r = 1; r < 4; ++r) {
    std::vector<VectorXd> qr(4);
    std::vector<double> er(4);
    for (int j = 0; j < 4; ++j) {
      qr[j] = q[(j + r) % 4];
      er[j] = eta[(j + r) % 4];
    }
    REQUIRE_THAT(energy(qr, er), WithinRel(explicitBath, 1e-12));
  }
  REQUIRE_THROWS_AS(addFrictionBath(q, flat, grad, {-0.1}, omegaP),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(addFrictionBath(q, flat, grad, {0.1}, 0.0),
                    std::invalid_argument);
}

TEST_CASE("The rate under a friction bath carries the bath's curvature",
          "[Tunneling][Instanton]") {
  // A cubic well at 0.6 T_c with eta = 0.5 omega0. The ring is searched
  // under the bath, half ring requested, and is stationary for U_N plus the
  // bath; the rate takes det' of the ring Hessian with the bath and Z_r
  // with eta omega_k on every k > 0 mode, as a dense finite-difference
  // Hessian of the same sum gives. Friction slows the rate.
  const double omega0 = 1.0;
  const double hw = kHbar * omega0;
  const double vb = 8.0 * hw;
  const CubicWell pes{omega0, std::sqrt(std::pow(omega0, 6) / (6.0 * vb))};
  const VectorXd saddle = VectorXd::Constant(1, pes.qb());
  const MatrixXd hs = pes.hessian(saddle);
  const double beta = 1.0 / (kBoltzmann * 0.6 * crossoverTemperature(hs));
  const double eta = 0.5;
  RateInstantonOptions opt;
  opt.beads = 32;
  opt.forceTolerance = 1e-8;
  opt.friction = true;
  opt.frictionEta = eta;
  RateInstanton wet =
      optimizeRateInstanton(saddle, hs, beta, {}, pes.batch(), opt);
  const long n = static_cast<long>(wet.beads.size());
  const double bnh = wet.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);
  auto total = [&](const std::vector<VectorXd> &x, double &u) {
    std::vector<double> v;
    std::vector<VectorXd> g;
    pes.batch()(x, v, g);
    u = 0.0;
    for (long j = 0; j < n; ++j) {
      const VectorXd &prev = x[static_cast<size_t>((j + n - 1) % n)];
      const VectorXd &next = x[static_cast<size_t>((j + 1) % n)];
      g[static_cast<size_t>(j)] +=
          c * (2.0 * x[static_cast<size_t>(j)] - prev - next);
      u += v[static_cast<size_t>(j)] +
           0.5 * c * (next - x[static_cast<size_t>(j)]).squaredNorm();
    }
    addFrictionBath(x, u, g, {eta}, std::sqrt(c));
    return g;
  };
  double u = 0.0;
  const std::vector<VectorXd> g = total(wet.beads, u);
  double residual = 0.0;
  for (const auto &gj : g) {
    residual = std::max(residual, gj.norm());
  }
  CAPTURE(wet.converged, wet.iterations, residual);
  REQUIRE(wet.converged);
  REQUIRE(wet.frictionEta == std::vector<double>{eta});
  REQUIRE(residual < 1e-6);
  REQUIRE_THAT(wet.ringPotential, WithinRel(u, 1e-10));

  instantonRate(
      wet, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(VectorXd::Zero(1)), 0.0, hs, vb);
  REQUIRE(wet.negativeModes == 1);
  // The dense Hessian of U_N plus the bath by central differences.
  MatrixXd hmat(n, n);
  for (long j = 0; j < n; ++j) {
    std::vector<VectorXd> up = wet.beads, down = wet.beads;
    up[static_cast<size_t>(j)](0) += 1e-5;
    down[static_cast<size_t>(j)](0) -= 1e-5;
    double scratch = 0.0;
    const auto gu = total(up, scratch);
    const auto gd = total(down, scratch);
    for (long i = 0; i < n; ++i) {
      hmat(i, j) =
          (gu[static_cast<size_t>(i)](0) - gd[static_cast<size_t>(i)](0)) /
          2e-5;
    }
  }
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(0.5 *
                                                   (hmat + hmat.transpose()));
  long zero = 0;
  for (long i = 1; i < n; ++i) {
    if (std::abs(es.eigenvalues()(i)) < std::abs(es.eigenvalues()(zero))) {
      zero = i;
    }
  }
  double logDetPrime = 0.0;
  for (long i = 0; i < n; ++i) {
    if (i != zero) {
      logDetPrime += std::log(std::abs(es.eigenvalues()(i)));
    }
  }
  const double omegaR2 = pes.hessian(VectorXd::Zero(1))(0, 0);
  double logZr = 0.0;
  for (long k = 0; k < n; ++k) {
    const double wk =
        2.0 * std::sqrt(c) * std::abs(std::sin(std::numbers::pi * k / n));
    logZr -= std::log(bnh) + 0.5 * std::log(omegaR2 + wk * wk + eta * wk);
  }
  const double logRate =
      -std::log(bnh) +
      0.5 * std::log(wet.bN / (2.0 * std::numbers::pi * bnh * kHbar)) -
      (static_cast<double>(n - 1) * std::log(bnh) + 0.5 * logDetPrime) -
      wet.betaN * u - logZr;
  CAPTURE(wet.logRate, logRate, wet.logZr, logZr, wet.zeroEigenvalue);
  REQUIRE_THAT(wet.logZr, WithinRel(logZr, 1e-10));
  REQUIRE_THAT(wet.logRate, Catch::Matchers::WithinAbs(logRate, 1e-5));

  opt.friction = false;
  RateInstanton dry =
      optimizeRateInstanton(saddle, hs, beta, {}, pes.batch(), opt);
  REQUIRE(dry.converged);
  instantonRate(
      dry, [&](long, const VectorXd &q) { return pes.hessian(q); },
      pes.hessian(VectorXd::Zero(1)), 0.0, hs, vb);
  REQUIRE(wet.logRate < dry.logRate);
}

namespace {

// A periodic square slab of L x L unit-mass atoms with central springs to
// the nearest (k) and next-nearest (k / 2) neighbours and a bending
// stiffness across each bond, no atom fixed. A
// rotation about the normal stretches the bonds across the cell boundary,
// a real curvature that falls only as 1 / L.
MatrixXd periodicSlabHessian(long l, MatrixXd &generators) {
  const long n = l * l;
  std::vector<Eigen::Vector3d> pos;
  for (long i = 0; i < l; ++i) {
    for (long j = 0; j < l; ++j) {
      pos.emplace_back(static_cast<double>(i), static_cast<double>(j), 0.0);
    }
  }
  const double box = static_cast<double>(l);
  MatrixXd h = MatrixXd::Zero(3 * n, 3 * n);
  auto add = [&](long p, long q, double k) {
    Eigen::Vector3d d =
        pos[static_cast<size_t>(q)] - pos[static_cast<size_t>(p)];
    for (int a = 0; a < 2; ++a) {
      d(a) -= box * std::round(d(a) / box);
    }
    const Eigen::Vector3d u = d.normalized();
    const Eigen::Matrix3d b = k * u * u.transpose();
    h.block(3 * p, 3 * p, 3, 3) += b;
    h.block(3 * q, 3 * q, 3, 3) += b;
    h.block(3 * p, 3 * q, 3, 3) -= b;
    h.block(3 * q, 3 * p, 3, 3) -= b;
    // A bending stiffness k_z (z_p - z_q)^2 / 2 keeps the slab flat.
    const double kz = 0.3 * k;
    h(3 * p + 2, 3 * p + 2) += kz;
    h(3 * q + 2, 3 * q + 2) += kz;
    h(3 * p + 2, 3 * q + 2) -= kz;
    h(3 * q + 2, 3 * p + 2) -= kz;
  };
  auto idx = [&](long i, long j) { return ((i + l) % l) * l + ((j + l) % l); };
  for (long i = 0; i < l; ++i) {
    for (long j = 0; j < l; ++j) {
      add(idx(i, j), idx(i + 1, j), 1.0);
      add(idx(i, j), idx(i, j + 1), 1.0);
      add(idx(i, j), idx(i + 1, j + 1), 0.5);
      add(idx(i, j), idx(i + 1, j - 1), 0.5);
    }
  }
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  for (const auto &p : pos) {
    com += p;
  }
  com /= static_cast<double>(n);
  generators = MatrixXd::Zero(3 * n, 6);
  for (long p = 0; p < n; ++p) {
    for (int c = 0; c < 3; ++c) {
      generators(3 * p + c, c) = 1.0;
      Eigen::Vector3d e = Eigen::Vector3d::Zero();
      e(c) = 1.0;
      generators.block(3 * p, 3 + c, 3, 1) =
          e.cross(pos[static_cast<size_t>(p)] - com);
    }
  }
  return h;
}

} // namespace

// The rotation of a 16 x 16 periodic slab about its normal has
// ||H r|| / (||H||_F ||r||) = 0.008, under the 1e-2 that marked a rotation
// as a zero mode, because ||H||_F grows as the square root of the
// coordinates. Against the softest vibration it is no zero mode. A free
// diatomic at its minimum keeps its two rotations as zero modes.
TEST_CASE("Rotational zero modes are told apart on the Hessian's own scale",
          "[Tunneling][Instanton]") {
  MatrixXd gens;
  const MatrixXd h = periodicSlabHessian(16, gens);
  const VectorXd r = gens.col(5);
  const double frobenius = (h * r).norm() / (h.norm() * r.norm());
  CAPTURE(frobenius);
  REQUIRE(frobenius < 1e-2);
  const RotationZeroModes slab = rotationZeroModes(h, gens);
  CAPTURE(slab.residual[2], slab.softestVibration);
  REQUIRE(slab.softestVibration > 0.0);
  REQUIRE_FALSE(slab.zero[2]);

  const CubicBond pes;
  const VectorXd x = pes.reference();
  MatrixXd g = MatrixXd::Zero(6, 6);
  Eigen::Vector3d com = 0.5 * (x.segment<3>(0) + x.segment<3>(3));
  for (int atom = 0; atom < 2; ++atom) {
    for (int c = 0; c < 3; ++c) {
      g(3 * atom + c, c) = 1.0;
      Eigen::Vector3d e = Eigen::Vector3d::Zero();
      e(c) = 1.0;
      g.block(3 * atom, 3 + c, 3, 1) = e.cross(x.segment<3>(3 * atom) - com);
    }
  }
  const RotationZeroModes bond =
      rotationZeroModes(pes.hessian(VectorXd::Zero(6)), g);
  CAPTURE(bond.residual[0], bond.residual[1], bond.residual[2],
          bond.softestVibration);
  REQUIRE_FALSE(bond.zero[0]); // about the bond: no generator
  REQUIRE(bond.zero[1]);
  REQUIRE(bond.zero[2]);
  REQUIRE(bond.softestVibration > 0.0);
}

// A double well with a small bump on the reactant slope: below the bump's
// top the reactant turning point jumps across it, so the period jumps too.
// A bisection over the whole energy range converged onto that jump, a ring
// whose period is not beta hbar (on the Al slab of examples/neb-al: period
// 103.8 for beta hbar 129.1 and 91.2 alike). Every requested period that a
// continuous branch reaches must be met.
TEST_CASE("The seed ring meets its period across a jump in the orbit",
          "[Tunneling][Instanton]") {
  std::vector<VectorXd> path;
  std::vector<double> energies;
  const int images = 2001;
  for (int i = 0; i < images; ++i) {
    const double x = -1.0 + 2.0 * static_cast<double>(i) / (images - 1);
    path.push_back(VectorXd::Constant(1, x));
    const double b = (x + 0.6) / 0.08;
    energies.push_back((x * x - 1.0) * (x * x - 1.0) +
                       0.25 * std::exp(-0.5 * b * b));
  }
  for (const double bh : {3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 5.0}) {
    RingSeed seed;
    const auto ring = ringFromPath(path, energies, bh, 32, &seed);
    CAPTURE(bh, seed.energy, seed.period);
    REQUIRE(ring.size() == 32);
    REQUIRE(seed.reached);
    REQUIRE_THAT(seed.period, WithinRel(bh, 1e-6));
  }
}

TEST_CASE("A link weight divides one spring of a closed ring",
          "[Tunneling][Instanton]") {
  const BatchPotential flat = [](const std::vector<VectorXd> &q,
                                 std::vector<double> &v,
                                 std::vector<VectorXd> &g) {
    v.assign(q.size(), 0.0);
    g.assign(q.size(), VectorXd::Zero(1));
  };
  std::vector<VectorXd> beads(4, VectorXd::Zero(1));
  beads[1](0) = 1.0;
  beads[2](0) = 1.0;
  const double c = 2.0;
  const double uniform = closedRingPotential(beads, c, flat);
  const double weighted =
      closedRingPotential(beads, c, flat, std::vector<double>{2.0, 1.0, 1.0, 1.0});
  REQUIRE_THAT(uniform, Catch::Matchers::WithinAbs(c, 1e-12));
  REQUIRE_THAT(weighted, Catch::Matchers::WithinAbs(0.75 * c, 1e-12));
  REQUIRE_THROWS_AS(
      closedRingPotential(beads, c, flat, std::vector<double>{0.0, 1.0, 1.0, 1.0}),
      std::invalid_argument);
}

TEST_CASE("A long ring takes the chain inertia instead of a dense factor",
          "[Tunneling][Instanton][large_ring]") {
  const BatchPotential hook = [](const std::vector<VectorXd> &q,
                                 std::vector<double> &v,
                                 std::vector<VectorXd> &g) {
    v.resize(q.size());
    g.resize(q.size());
    for (size_t i = 0; i < q.size(); ++i) {
      v[i] = 0.5 * q[i].squaredNorm();
      g[i] = q[i];
    }
  };
  const VectorXd saddle = VectorXd::Zero(1);
  const MatrixXd hs = -MatrixXd::Identity(1, 1);
  RateInstantonOptions opt;
  opt.beads = 4200;
  opt.halfRing = false;
  opt.maxIterations = 1;
  opt.forceTolerance = 1.0e6;
  opt.checkOddSector = false;
  opt.lanczosFirst = 2;
  opt.lanczosRestart = 2;
  const RateInstanton inst =
      optimizeRateInstanton(saddle, hs, 60.0 / kHbar, {}, hook, opt);
  REQUIRE(inst.iterations >= 0);
  REQUIRE(std::isfinite(inst.ringPotential));
}

TEST_CASE("instanton checks reject a short path and a bad friction list",
          "[Tunneling][Instanton]") {
  auto params = std::make_shared<Parameters>();
  ParametersLoadAccess::potential_options(*params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, *params));
  Matter few(pot, *params);
  few.resize(1);
  few.setMasses(VectorXd::Ones(1));
  Matter more(pot, *params);
  more.resize(2);
  more.setMasses(VectorXd::Ones(2));
  REQUIRE_THROWS_AS(massWeightedDistance(few, more), std::invalid_argument);

  std::vector<VectorXd> path;
  std::vector<double> energies;
  for (int i = 0; i < 4; ++i) {
    path.push_back(VectorXd::Constant(1, static_cast<double>(i)));
    energies.push_back(static_cast<double>(i));
  }
  REQUIRE_THROWS_AS(ringFromPath(path, energies, 30.0, 8),
                    std::invalid_argument);
  path[2] = VectorXd::Zero(2);
  REQUIRE_THROWS_AS(ringFromPath(path, energies, 30.0, 8),
                    std::invalid_argument);

  Instanton blank;
  BeadHessian hess = [](long, const VectorXd &) {
    return MatrixXd::Zero(1, 1);
  };
  REQUIRE_THROWS_AS(instantonSplitting(blank, hess, MatrixXd::Zero(1, 1),
                                       MatrixXd::Zero(1, 1)),
                    std::invalid_argument);
  blank.dtau = 0.2;
  blank.path.assign(5, VectorXd::Zero(1));
  BeadHessian wide = [](long, const VectorXd &) {
    return MatrixXd::Zero(2, 2);
  };
  REQUIRE_THROWS_AS(instantonSplitting(blank, wide, MatrixXd::Zero(1, 1),
                                       MatrixXd::Zero(1, 1)),
                    std::runtime_error);

  const BatchPotential hook = [](const std::vector<VectorXd> &q,
                                 std::vector<double> &v,
                                 std::vector<VectorXd> &g) {
    v.resize(q.size());
    g.resize(q.size());
    for (size_t i = 0; i < q.size(); ++i) {
      v[i] = 0.5 * q[i].squaredNorm();
      g[i] = q[i];
    }
  };
  const VectorXd saddle = VectorXd::Zero(1);
  const MatrixXd hs = -MatrixXd::Identity(1, 1);
  RateInstantonOptions opt;
  opt.beads = 8;
  opt.halfRing = false;
  opt.maxIterations = 1;
  opt.forceTolerance = 1.0e6;
  opt.checkOddSector = false;
  opt.lanczosFirst = 2;
  opt.lanczosRestart = 2;
  opt.friction = true;
  opt.frictionExplicit = true;
  opt.frictionEtaBeads = {0.1};
  REQUIRE_THROWS_AS(
      optimizeRateInstanton(saddle, hs, 60.0 / kHbar, {}, hook, opt),
      std::invalid_argument);

  opt.friction = false;
  opt.frictionExplicit = false;
  opt.frictionEtaBeads.clear();
  opt.discretization.assign(8, 1.0);
  const RateInstanton weighted =
      optimizeRateInstanton(saddle, hs, 60.0 / kHbar, {}, hook, opt);
  REQUIRE(std::isfinite(weighted.ringPotential));

  const std::vector<MatrixXd> blocks(2, MatrixXd::Identity(1, 1));
  const std::vector<VectorXd> tau(1, VectorXd::Ones(1));
  REQUIRE_THROWS_AS(ringSpectrum(blocks, 1.0, tau), std::invalid_argument);
  REQUIRE_THROWS_AS(cyclicRingLogAbsDet(-1.0, blocks), std::invalid_argument);

  RateInstanton bare;
  RingBeadHessian bead = [](long, const VectorXd &) {
    return MatrixXd::Identity(1, 1);
  };
  REQUIRE_THROWS_AS(instantonRate(bare, bead, MatrixXd::Identity(1, 1), 0.0),
                    std::invalid_argument);
  bare.beads.assign(4, VectorXd::Zero(1));
  bare.betaN = 1.0;
  bare.discretization = {1.0, 1.0, 1.0, -1.0};
  REQUIRE_THROWS_AS(instantonRate(bare, bead, MatrixXd::Identity(1, 1), 0.0),
                    std::invalid_argument);
  bare.discretization = {1.0, 1.0, 1.0, 1.0};
  try {
    instantonRate(bare, bead, MatrixXd::Identity(1, 1), 0.0);
  } catch (const std::exception &) {
  }
}

TEST_CASE("a ring structure packs beads and rejects a bad spring",
          "[Tunneling][Ring]") {
  const BatchPotential hook = [](const std::vector<VectorXd> &q,
                                 std::vector<double> &v,
                                 std::vector<VectorXd> &g) {
    v.resize(q.size());
    g.resize(q.size());
    for (size_t i = 0; i < q.size(); ++i) {
      v[i] = 0.5 * q[i].squaredNorm();
      g[i] = q[i];
    }
  };
  const VectorXd mode = VectorXd::Ones(1);
  REQUIRE_THROWS_AS(RingPolymerPotential(hook, 4, 1, 0.0, 0.0, {}, false, mode),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(RingPolymerPotential(hook, 3, 1, 1.0, 0.0, {}, false, mode),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(RingPolymerPotential(hook, 5, 1, 1.0, 0.0, {}, true, mode),
                    std::invalid_argument);

  RingPolymerPotential ring(hook, 4, 1, 1.0, 0.0, {0.1}, false, mode);
  REQUIRE(ring.structureAtoms() > 0);
  std::vector<double> coords(static_cast<size_t>(3 * ring.structureAtoms()),
                             0.0);
  std::vector<double> forces(coords.size(), 0.0);
  double energy = 0.0;
  double variance = 0.0;
  ring.force(ring.structureAtoms(), coords.data(), nullptr, forces.data(),
             &energy, &variance, nullptr);
  REQUIRE(std::isfinite(energy));

  std::vector<double> coords2 = coords;
  coords2[0] = 0.2;
  const double *pos[2] = {coords.data(), coords2.data()};
  double *frc[2] = {forces.data(), forces.data()};
  double energies[2] = {};
  double variances[2] = {};
  ring.forceBatch(2, ring.structureAtoms(), pos, nullptr, frc, energies,
                  variances, nullptr);
  REQUIRE(std::isfinite(energies[0]));
  REQUIRE(std::isfinite(energies[1]));

  REQUIRE_THROWS_AS(ring.unpack(nullptr), std::invalid_argument);
  REQUIRE_THROWS_AS(ring.packActive({}, coords.data()), std::invalid_argument);
  AtomMatrix bad(1, 3);
  REQUIRE_THROWS_AS(ring.packMode(mode, bad), std::invalid_argument);
  AtomMatrix packed(ring.structureAtoms(), 3);
  REQUIRE_THROWS_AS(ring.packMode(VectorXd::Zero(1), packed),
                    std::invalid_argument);
  ring.packMode(mode, packed);
  REQUIRE(packed.norm() == Catch::Approx(1.0));
  const auto active = ring.unpack(coords.data());
  REQUIRE(static_cast<long>(active.size()) == ring.activeBeads());

  RingPolymerPotential folded(hook, 4, 1, 1.0, 0.0, {}, true, mode);
  std::vector<double> foldCoords(
      static_cast<size_t>(3 * folded.structureAtoms()), 0.0);
  std::vector<double> foldForces(foldCoords.size(), 0.0);
  double foldEnergy = 0.0;
  folded.force(folded.structureAtoms(), foldCoords.data(), nullptr,
               foldForces.data(), &foldEnergy, nullptr, nullptr);
  REQUIRE(std::isfinite(foldEnergy));
}

TEST_CASE("a rate ring takes the Lanczos determinant and a friction bath",
          "[Tunneling][Instanton][rate]") {
  const BatchPotential hook = [](const std::vector<VectorXd> &q,
                                 std::vector<double> &v,
                                 std::vector<VectorXd> &g) {
    v.assign(q.size(), 0.0);
    g.assign(q.size(), VectorXd::Zero(1));
    for (size_t i = 0; i < q.size(); ++i) {
      if (q[i].size() < 1) {
        continue;
      }
      v[i] = -0.5 * q[i](0) * q[i](0);
      g[i] = VectorXd::Constant(1, -q[i](0));
    }
  };
  const BatchPotential wrong = [](const std::vector<VectorXd> &,
                                  std::vector<double> &v,
                                  std::vector<VectorXd> &g) {
    v.clear();
    g.clear();
  };
  const VectorXd saddle = VectorXd::Zero(1);
  const MatrixXd hs = -MatrixXd::Identity(1, 1);
  RateInstantonOptions opt;
  opt.beads = 8;
  opt.halfRing = false;
  opt.maxIterations = 1;
  opt.forceTolerance = 1.0e6;
  opt.checkOddSector = false;
  opt.lanczosFirst = 2;
  opt.lanczosRestart = 2;
  opt.friction = true;
  opt.frictionEta = 0.0;
  const RateInstanton dry =
      optimizeRateInstanton(saddle, hs, 60.0 / kHbar, {}, hook, opt);
  REQUIRE(std::isfinite(dry.ringPotential));
  opt.frictionEta = 0.25;
  const RateInstanton wet =
      optimizeRateInstanton(saddle, hs, 60.0 / kHbar, {}, hook, opt);
  REQUIRE(std::isfinite(wet.ringPotential));
  REQUIRE_THROWS_AS(
      optimizeRateInstanton(saddle, hs, 60.0 / kHbar, {}, wrong, opt),
      std::exception);

  RateInstanton inst;
  inst.beta = 1.0;
  inst.betaN = 0.25;
  inst.bN = 0.04;
  inst.ringPotential = 0.1;
  inst.beads = {VectorXd::Constant(1, 0.0), VectorXd::Constant(1, 0.2),
                VectorXd::Constant(1, 0.0), VectorXd::Constant(1, -0.2)};
  const RingBeadHessian bead = [](long, const VectorXd &) {
    return MatrixXd::Constant(1, 1, -0.4);
  };
  const MatrixXd reactant = MatrixXd::Constant(1, 1, 1.0);
  try {
    instantonRate(inst, bead, reactant, 0.0, MatrixXd(), 0.0, 0, 0);
  } catch (const std::exception &) {
  }
  REQUIRE(inst.beads.size() == 4);
  REQUIRE_THROWS_AS(
      instantonRate(inst, bead, reactant, 0.0, MatrixXd(), 0.0, 1, 0),
      std::runtime_error);
  REQUIRE_THROWS_AS(instantonRate(inst, bead, -reactant, 0.0, MatrixXd(), 0.0,
                                  0, 4096),
                    std::runtime_error);
}

TEST_CASE("a double well cools a short rate ring",
          "[Tunneling][Instanton][well]") {
  const BatchPotential well = [](const std::vector<VectorXd> &q,
                                 std::vector<double> &v,
                                 std::vector<VectorXd> &g) {
    v.resize(q.size());
    g.resize(q.size());
    for (size_t i = 0; i < q.size(); ++i) {
      const double x = q[i].size() > 0 ? q[i](0) : 0.0;
      const double x2 = x * x;
      v[i] = (x2 - 1.0) * (x2 - 1.0);
      g[i] = VectorXd::Constant(q[i].size() > 0 ? q[i].size() : 1,
                                4.0 * x * (x2 - 1.0));
      if (q[i].size() == 0) {
        g[i].resize(0);
      }
    }
  };
  const VectorXd saddle = VectorXd::Zero(1);
  const MatrixXd hs = MatrixXd::Constant(1, 1, -4.0);
  RateInstantonOptions opt;
  opt.beads = 16;
  opt.halfRing = false;
  opt.maxIterations = 12;
  opt.forceTolerance = 1.0e-2;
  opt.checkOddSector = true;
  opt.lanczosFirst = 6;
  opt.lanczosRestart = 4;
  opt.newtonLimit = 64;
  const RateInstanton inst =
      optimizeRateInstanton(saddle, hs, 40.0 / kHbar, {}, well, opt);
  REQUIRE(inst.iterations >= 1);
  REQUIRE(std::isfinite(inst.ringPotential));

  const VectorXd left = VectorXd::Constant(1, -1.0);
  const VectorXd right = VectorXd::Constant(1, 1.0);
  InstantonOptions split;
  split.beads = 8;
  split.maxIterations = 6;
  split.forceTolerance = 1.0e-2;
  const Instanton pair =
      optimizeInstanton(left, right, 8.0, {}, well, split);
  REQUIRE(pair.iterations >= 1);
  REQUIRE(std::isfinite(pair.action));
}
