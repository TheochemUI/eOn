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
#include "EckartBarrier.hpp"
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
    const Eigen::Matrix3d a = d2vd(r - re) * u * u.transpose() +
                              dvd(r - re) / r *
                                  (Eigen::Matrix3d::Identity() -
                                   u * u.transpose());
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
  std::sort(lam.begin(), lam.end(), [](double a, double b) {
    return std::abs(a) < std::abs(b);
  });
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
  REQUIRE_THAT(ring.logRateTimesZr,
               Catch::Matchers::WithinAbs(expected, 1e-6));
  // The reactant's rotation generators on every bead miss the ring's null
  // space by a stretch-dependent angle.
  REQUIRE(std::abs(copied.logRateTimesZr - expected) > 1e-3);
}
