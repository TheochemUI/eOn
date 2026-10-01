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
#include "eon/PIQTST.h"

#include "EckartBarrier.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Tunneling.h"


#include <cmath>
#include <numbers>
#include <vector>

using namespace eonc;
using eonc::testing::Eckart;

namespace {

// V(x, y) = V0 sech^2(x / a) + k0 (1 + c sech^2(x / a)) y^2 / 2 at unit
// mass. With c = 0 and y fixed it is the one-dimensional Eckart barrier.
struct EckartPot final : Potential {
  Eckart pes;
  double k0{0.0};
  double c{0.0};

  EckartPot()
      : Potential(PotType::LJ) {}

  void evaluate(const double *x, double *f, double *e) const {
    const double sech2 = 1.0 / std::pow(std::cosh(x[0] / pes.a), 2);
    const double k = k0 * (1.0 + c * sech2);
    const double dsech2 = -2.0 * sech2 * std::tanh(x[0] / pes.a) / pes.a;
    *e = pes.value(x[0]) + 0.5 * k * x[1] * x[1];
    f[0] = -pes.slope(x[0]) - 0.5 * k0 * c * dsech2 * x[1] * x[1];
    f[1] = -k * x[1];
    f[2] = 0.0;
  }

  void force(long, const double *positions, const int *, double *forces,
             double *energy, double *variance, const double *) override {
    evaluate(positions, forces, energy);
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }

  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return true;
  }

  void forceBatch(long nSystems, long, const double *const *positions,
                  const int *const *, double *const *forces, double *energies,
                  double *variances, const double *const *) override {
    for (long s = 0; s < nSystems; ++s) {
      evaluate(positions[s], forces[s], &energies[s]);
      if (variances != nullptr) {
        variances[s] = 0.0;
      }
    }
  }
};

piqtst::Coordinate line(bool transverse) {
  piqtst::Coordinate c;
  c.atoms = 1;
  c.masses = {1.0};
  c.numbers = {1};
  c.free = {1, static_cast<char>(transverse ? 1 : 0), 0};
  c.reference = VectorXd::Zero(3);
  c.direction = VectorXd::Zero(3);
  c.direction(0) = 1.0;
  return c;
}

// Planes from -5a to the top: spacing 0.1 Angstrom out to -1.5a, where the
// mean force is small and nearly noiseless, then nTop intervals up to the
// top, where the ring's fluctuations carry the variance.
std::vector<double> planes(double a, long nTop) {
  std::vector<double> s;
  const double knee = -1.5 * a;
  const long nFar = std::lround((knee + 5.0 * a) / 0.1);
  for (long j = 0; j < nFar; ++j) {
    s.push_back(-5.0 * a + (knee + 5.0 * a) * static_cast<double>(j) /
                               static_cast<double>(nFar));
  }
  for (long j = 0; j <= nTop; ++j) {
    s.push_back(knee -
                knee * static_cast<double>(j) / static_cast<double>(nTop));
  }
  return s;
}

pathintegral::Options ring(double temperature, long beads, double dt) {
  pathintegral::Options o;
  o.beads = beads;
  o.temperature = temperature;
  o.kB = tunneling::kBoltzmann;
  o.hbar = tunneling::kHbar;
  o.dt = dt;
  o.pileTau = 1.0;
  o.seed = 20261001;
  return o;
}

// With F = 0 far on the reactant side, the PI-QTST rate per unit reactant
// density is k Z_r = exp(-beta (F(s*) - F(-inf))) / (2 pi beta hbar), the
// one-dimensional limit of the reactant-integral formula in piqtst::rate.
// The planes start at -5a, where V is 8e-5 eV.
double logFlux(double barrier, double beta) {
  return -beta * barrier -
         std::log(2.0 * std::numbers::pi * beta * tunneling::kHbar);
}

piqtst::ScanOptions eckartScan(double temperature, long beads) {
  piqtst::ScanOptions o;
  o.planes = planes(Eckart{}.a, 40);
  o.equilibration = 100;
  o.production = 2000;
  o.blocks = 10;
  o.ring = ring(temperature, beads, 0.1);
  o.ring.pileScale = 0.5;
  return o;
}

} // namespace

// The reference barriers are the same centroid free energy at the same
// bead count by coupling-constant integration over free-ring normal modes
// with Metropolis sampling (an independent sampler), to 5e-5 eV. The
// allowance beside the sampling error covers the trapezoid rule on these
// planes (2e-4 eV from the end curvatures) and the 0.1 time step.
//
// Against the exact flux: for a symmetric Eckart barrier PI-QTST lies below
// the exact rate, by about 10 percent at T_c and 30 percent at 0.7 T_c (the
// 64-bead reference here: 0.89 and 0.71). Ring-polymer MD, PI-QTST times a
// transmission factor near 1 for this barrier, is 15 percent low just above
// T_c, 23 percent at 0.79 T_c and 45 percent at 0.53 T_c for the H + H2
// Eckart barrier (Suleimanov, Aoiz and Guo, J. Phys. Chem. A 120, 8488
// (2016), table 1). The bands below bracket those, widened by three
// standard errors of the sampled rate.
TEST_CASE("PI-QTST matches the exact Eckart flux near and below the "
          "crossover",
          "[PIQTST]") {
  const Eckart pes;
  const double tc = tunneling::crossoverTemperature(pes.hessian_at_top());
  EckartPot pot;
  struct Case {
    double fraction;
    long beads;
    double reference;
    double low;
    double high;
  };
  for (const Case cs : {Case{1.0, 16, 0.39653, 0.80, 1.0},
                        Case{0.7, 24, 0.36276, 0.60, 0.95}}) {
    const double t = cs.fraction * tc;
    const double beta = 1.0 / (tunneling::kBoltzmann * t);
    const auto result = piqtst::scan(pot, line(false), eckartScan(t, cs.beads));
    const double barrier = result.back().freeEnergy;
    const double error = result.back().freeEnergyError;
    const double ratio =
        std::exp(logFlux(barrier, beta) - pes.logExactFlux(beta));
    CAPTURE(cs.fraction, cs.beads, barrier, error, ratio);
    REQUIRE(error > 0.0);
    REQUIRE(beta * error < 0.15);
    REQUIRE(std::abs(barrier - cs.reference) < 3.0 * error + 5e-4);
    REQUIRE(ratio > cs.low * std::exp(-3.0 * beta * error));
    REQUIRE(ratio < cs.high * std::exp(3.0 * beta * error));
    // The top plane is the dividing surface; F' vanishes there.
    REQUIRE(std::abs(result.back().meanForce) <
            3.0 * result.back().meanForceError + 1e-3);
    // The ring spreads along the free coordinate only, more at the top.
    REQUIRE(result.back().spread[1] == 0.0);
    REQUIRE(result.back().spread[0] > result.front().spread[0]);
  }
}

// One bead is classical TST along s. A transverse mode whose stiffness
// rises by 1 + c at the top adds (kT / 2) ln((1 + c sech^2) / ...) to the
// free energy, sampled by the centroid thermostat within each plane; c = 0
// leaves the bare barrier and the rate exp(-beta V0) / (2 pi beta hbar).
TEST_CASE("One bead gives classical TST on the free-energy profile",
          "[PIQTST]") {
  const Eckart pes;
  const double t = 300.0;
  const double beta = 1.0 / (tunneling::kBoltzmann * t);
  const double s0 = -5.0 * pes.a;
  const double sech0 = 1.0 / std::pow(std::cosh(s0 / pes.a), 2);
  for (const double c : {0.0, 1.0}) {
    EckartPot pot;
    pot.k0 = 1.0;
    pot.c = c;
    piqtst::ScanOptions o;
    o.planes = planes(pes.a, 40);
    o.equilibration = 100;
    o.production = 1000;
    o.blocks = 10;
    o.ring = ring(t, 1, 0.1);
    const auto result = piqtst::scan(pot, line(true), o);
    const double barrier = result.back().freeEnergy;
    const double error = result.back().freeEnergyError;
    const double exact = pes.v0 - pes.value(s0) +
                         0.5 / beta * std::log((1.0 + c) / (1.0 + c * sech0));
    CAPTURE(c, barrier, error, exact);
    REQUIRE(std::abs(barrier - exact) < 3.0 * error + 5e-4);
    const double logRate = logFlux(barrier, beta);
    const double classical = -beta * pes.v0 - std::log(2.0 * std::numbers::pi *
                                                       beta * tunneling::kHbar);
    if (c == 0.0) {
      // No transverse sampling: the error is the trapezoid rule alone.
      REQUIRE(error < 1e-12);
      REQUIRE(std::abs(logRate - classical) < beta * 5e-4);
    } else {
      REQUIRE(error > 0.0);
    }
  }
}

// A harmonic well F = omega^2 s^2 / 2 cut by a plane at s* gives classical
// TST, omega / (2 pi) exp(-beta F(s*)), when the planes span the well.
TEST_CASE("The PI-QTST rate of a harmonic profile is harmonic TST",
          "[PIQTST]") {
  const double omega = 0.8;
  const double beta = 50.0;
  const double sStar = 0.6;
  std::vector<piqtst::Plane> planes;
  const long n = 2001;
  for (long j = 0; j < n; ++j) {
    piqtst::Plane p;
    p.s = -sStar +
          2.0 * sStar * static_cast<double>(j) / static_cast<double>(n - 1);
    p.meanForce = omega * omega * p.s;
    p.meanForceError = 0.0;
    planes.push_back(p);
  }
  piqtst::integrate(planes);
  const piqtst::Rate r = piqtst::rate(planes, beta);
  const double fStar = 0.5 * omega * omega * sStar * sStar;
  REQUIRE(r.reactant == (n - 1) / 2);
  REQUIRE_THAT(r.barrier, Catch::Matchers::WithinAbs(fStar, 1e-6));
  REQUIRE(r.firstPlaneHeight > 5.0);
  // The reactant integral is cut at -s*, a 1e-3 relative tail here.
  const double expected =
      std::log(omega / (2.0 * std::numbers::pi)) - beta * fStar +
      std::log(1.0 / std::erf(sStar * omega * std::sqrt(0.5 * beta)));
  REQUIRE_THAT(r.logRate, Catch::Matchers::WithinAbs(expected, 1e-5));
  REQUIRE(r.logRateError == 0.0);
}

TEST_CASE("PI-QTST errors propagate through the trapezoid rule", "[PIQTST]") {
  std::vector<piqtst::Plane> planes(3);
  for (long j = 0; j < 3; ++j) {
    planes[static_cast<size_t>(j)].s = 0.5 * static_cast<double>(j);
    planes[static_cast<size_t>(j)].meanForce = 1.0 - static_cast<double>(j);
    planes[static_cast<size_t>(j)].meanForceError = 0.1;
  }
  piqtst::integrate(planes);
  // F_2 = 0.25 (F'_0 + 2 F'_1 + F'_2), error 0.25 sqrt(1 + 4 + 1) * 0.1.
  REQUIRE_THAT(planes[2].freeEnergy, Catch::Matchers::WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(planes[1].freeEnergy, Catch::Matchers::WithinAbs(0.25, 1e-14));
  REQUIRE_THAT(planes[2].freeEnergyError,
               Catch::Matchers::WithinRel(0.025 * std::sqrt(6.0), 1e-12));
  const piqtst::Rate r = piqtst::rate(planes, 10.0);
  REQUIRE(r.reactant == 0);
  REQUIRE_THAT(r.barrierError,
               Catch::Matchers::WithinRel(0.025 * std::sqrt(6.0), 1e-12));
  REQUIRE(r.logRateError > 0.0);
}
