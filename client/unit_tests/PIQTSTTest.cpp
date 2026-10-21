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
#include "eon/Parameters.h"
#include "eon/ParametersJSON.h"
#include "eon/Tunneling.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <stdexcept>
#include <vector>

using namespace eonc;
using eonc::testing::Eckart;

namespace {

// In mass-weighted coordinates q = sqrt(mass) (x, y), rotated by theta to
// xi along u = (cos theta, sin theta) and eta across it,
// V = V0 sech^2(xi / a) + k0 (1 + c sech^2(xi / a)) eta^2 / 2.
// With unit mass, theta = 0, c = 0 and y fixed it is the one-dimensional
// Eckart barrier.
struct EckartPot final : Potential {
  Eckart pes;
  double k0{0.0};
  double c{0.0};
  double mass{1.0};
  double theta{0.0};

  EckartPot()
      : Potential(PotType::LJ) {}

  void evaluate(const double *x, double *f, double *e) const {
    const double sm = std::sqrt(mass);
    const double ct = std::cos(theta), st = std::sin(theta);
    const double xi = sm * (ct * x[0] + st * x[1]);
    const double eta = sm * (-st * x[0] + ct * x[1]);
    const double sech2 = 1.0 / std::pow(std::cosh(xi / pes.a), 2);
    const double k = k0 * (1.0 + c * sech2);
    const double dsech2 = -2.0 * sech2 * std::tanh(xi / pes.a) / pes.a;
    *e = pes.value(xi) + 0.5 * k * eta * eta;
    const double dxi = pes.slope(xi) + 0.5 * k0 * c * dsech2 * eta * eta;
    const double deta = k * eta;
    f[0] = -sm * (ct * dxi - st * deta);
    f[1] = -sm * (st * dxi + ct * deta);
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

piqtst::Coordinate line(bool transverse, double mass = 1.0,
                        double theta = 0.0) {
  piqtst::Coordinate c;
  c.atoms = 1;
  c.masses = {mass};
  c.numbers = {1};
  c.free = {1, static_cast<char>(transverse ? 1 : 0), 0};
  c.reference = VectorXd::Zero(3);
  c.direction = VectorXd::Zero(3);
  c.direction(0) = std::cos(theta);
  c.direction(1) = std::sin(theta);
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
// rises by 1 + c at the top adds (kT / 2) ln((1 + c sech^2(xi / a)) / ...)
// to the free energy, sampled by the centroid thermostat within each
// plane; c = 0 leaves the bare barrier and the rate
// exp(-beta V0) / (2 pi beta hbar). A mass of 4 amu with the barrier along
// a rotated mass-weighted direction is the same profile in s.
TEST_CASE("One bead gives classical TST on the free-energy profile",
          "[PIQTST]") {
  const Eckart pes;
  const double t = 300.0;
  const double beta = 1.0 / (tunneling::kBoltzmann * t);
  const double s0 = -5.0 * pes.a;
  const double sech0 = 1.0 / std::pow(std::cosh(s0 / pes.a), 2);
  struct Case {
    double c, mass, theta;
  };
  for (const Case cs :
       {Case{0.0, 1.0, 0.0}, Case{1.0, 1.0, 0.0}, Case{1.0, 4.0, 0.6}}) {
    EckartPot pot;
    pot.k0 = 1.0;
    pot.c = cs.c;
    pot.mass = cs.mass;
    pot.theta = cs.theta;
    piqtst::ScanOptions o;
    o.planes = planes(pes.a, 40);
    o.equilibration = 100;
    o.production = 1000;
    o.blocks = 10;
    o.ring = ring(t, 1, 0.1);
    const auto result = piqtst::scan(pot, line(true, cs.mass, cs.theta), o);
    const double barrier = result.back().freeEnergy;
    const double error = result.back().freeEnergyError;
    const double exact =
        pes.v0 - pes.value(s0) +
        0.5 / beta * std::log((1.0 + cs.c) / (1.0 + cs.c * sech0));
    CAPTURE(cs.c, cs.mass, cs.theta, barrier, error, exact);
    REQUIRE(std::abs(barrier - exact) < 3.0 * error + 5e-4);
    const double logRate = logFlux(barrier, beta);
    const double classical = -beta * pes.v0 - std::log(2.0 * std::numbers::pi *
                                                       beta * tunneling::kHbar);
    if (cs.c == 0.0) {
      // No transverse coupling: the error is the trapezoid rule alone.
      REQUIRE(error < 1e-12);
      REQUIRE(std::abs(logRate - classical) < beta * 5e-4);
    } else {
      REQUIRE(error > 0.0);
    }
    // The centroid stays on the line it was seeded on, within three
    // transverse thermal widths sqrt(kT / k0) in mass-weighted units.
    const VectorXd mid = result[result.size() / 2].centroid;
    const double across = std::sqrt(cs.mass) * (-std::sin(cs.theta) * mid(0) +
                                                std::cos(cs.theta) * mid(1));
    REQUIRE(std::abs(across) < 3.0 * std::sqrt(1.0 / (beta * pot.k0)));
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

// F falling all the way to s*: the reactant is the lowest plane before s*,
// so the barrier comes out negative, where the minimum over every plane
// took s* itself and reported a barrier of zero. The rate does not depend
// on the reference.
TEST_CASE("PI-QTST takes the reactant before the dividing plane", "[PIQTST]") {
  std::vector<piqtst::Plane> planes(5);
  for (long j = 0; j < 5; ++j) {
    planes[static_cast<size_t>(j)].s = 0.25 * static_cast<double>(j);
    planes[static_cast<size_t>(j)].meanForce = -0.4;
  }
  piqtst::integrate(planes);
  const piqtst::Rate r = piqtst::rate(planes, 10.0);
  REQUIRE(r.reactant == 3);
  REQUIRE_THAT(r.barrier, Catch::Matchers::WithinAbs(-0.1, 1e-14));
  double z = 0.0;
  for (long j = 0; j < 5; ++j) {
    const double width = j == 0 || j == 4 ? 0.125 : 0.25;
    z += width * std::exp(-10.0 * planes[static_cast<size_t>(j)].freeEnergy);
  }
  const double expected =
      std::log(0.5 * std::sqrt(2.0 / (std::numbers::pi * 10.0))) -
      10.0 * planes.back().freeEnergy - std::log(z);
  REQUIRE_THAT(r.logRate, Catch::Matchers::WithinAbs(expected, 1e-12));
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

TEST_CASE("PI-QTST options round-trip through JSON and are checked",
          "[PIQTST][params]") {
  nlohmann::json j = {{"Instanton",
                       {{"mode", "rate"},
                        {"pi_planes", 9},
                        {"pi_beads", 12},
                        {"pi_equilibration_steps", 30},
                        {"pi_sampling_steps", 400},
                        {"pi_time_step", 0.25},
                        {"pi_thermostat", "PILE"},
                        {"pi_pile_tau", 50.0},
                        {"pi_pile_scale", 0.5},
                        {"pi_seed", 7},
                        {"pi_direction", "Line"},
                        {"pi_reactant_extent", 0.75}}}};
  Parameters p;
  eonc::config::from_json(j, p);
  const auto &o = p.instanton_options();
  REQUIRE(o.pi_planes == 9);
  REQUIRE(o.pi_beads == 12);
  REQUIRE(o.pi_equilibration_steps == 30);
  REQUIRE(o.pi_sampling_steps == 400);
  REQUIRE(o.pi_time_step == Catch::Approx(0.25));
  REQUIRE(o.pi_thermostat == "pile");
  REQUIRE(o.pi_pile_tau == Catch::Approx(50.0));
  REQUIRE(o.pi_pile_scale == Catch::Approx(0.5));
  REQUIRE(o.pi_seed == 7);
  REQUIRE(o.pi_direction == "line");
  REQUIRE(o.pi_reactant_extent == Catch::Approx(0.75));
  const nlohmann::json back = eonc::config::to_json(p);
  Parameters again;
  eonc::config::from_json(back, again);
  REQUIRE(again.instanton_options().pi_planes == 9);
  REQUIRE(again.instanton_options().pi_direction == "line");
  REQUIRE(again.instanton_options().pi_pile_scale == Catch::Approx(0.5));

  for (const nlohmann::json &bad :
       {nlohmann::json{{"Instanton", {{"pi_planes", 4}}}},
        nlohmann::json{{"Instanton", {{"mode", "rate"}, {"pi_planes", 1}}}},
        nlohmann::json{{"Instanton",
                        {{"mode", "rate"},
                         {"pi_planes", 4},
                         {"pi_thermostat", "piglet"}}}},
        nlohmann::json{{"Instanton",
                        {{"mode", "rate"},
                         {"pi_planes", 4},
                         {"pi_direction", "sideways"}}}}}) {
    Parameters q;
    CAPTURE(bad.dump());
    REQUIRE_THROWS_AS(eonc::config::from_json(bad, q), std::invalid_argument);
  }
}

namespace {

piqtst::RecrossingOptions recrossingAtTop(double temperature, long beads,
                                          long parents, long children) {
  piqtst::RecrossingOptions o;
  o.s = 0.0;
  o.equilibration = 200;
  o.parents = parents;
  o.spacing = 20;
  o.children = children;
  o.steps = 100;
  o.ring = ring(temperature, beads, 0.1);
  return o;
}

struct Estimate {
  double kappa{0.0};
  double error{0.0};
};

// Classical flux-side transmission through x = 0 on the rotated barrier by
// a separate sampler: Metropolis on y along the line x = 0, Gaussian
// velocities at kT, plain velocity Verlet, and the same estimator as the
// RPMD children (momentum-reversed pairs, mean of kappa(t) over the last
// quarter, jackknife over blocks of pairs).
Estimate classicalTransmission(const EckartPot &pot, double beta, long pairs,
                               long blocks, long steps, double dt) {
  std::mt19937_64 rng(4242);
  std::normal_distribution<double> gauss(0.0, 1.0);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  auto energy = [&](double x, double y) {
    const double q[3] = {x, y, 0.0};
    double f[3];
    double e = 0.0;
    pot.evaluate(q, f, &e);
    return e;
  };
  double y = 0.0;
  double ey = energy(0.0, y);
  const double width = 0.5 / std::sqrt(beta * pot.k0);
  auto metropolis = [&]() {
    for (int sweep = 0; sweep < 10; ++sweep) {
      const double trial = y + width * (2.0 * uniform(rng) - 1.0);
      const double et = energy(0.0, trial);
      if (uniform(rng) < std::exp(-beta * (et - ey))) {
        y = trial;
        ey = et;
      }
    }
  };
  for (int burn = 0; burn < 200; ++burn) {
    metropolis();
  }
  const long first = steps - steps / 4;
  const double span = static_cast<double>(steps + 1 - first);
  const double sigma = 1.0 / std::sqrt(beta * pot.mass);
  std::vector<double> num(static_cast<size_t>(blocks), 0.0);
  std::vector<double> den(static_cast<size_t>(blocks), 0.0);
  const long perBlock = pairs / blocks;
  for (long block = 0; block < blocks; ++block) {
    for (long pair = 0; pair < perBlock; ++pair) {
      metropolis();
      const double vx0 = sigma * gauss(rng);
      const double vy0 = sigma * gauss(rng);
      for (const double sign : {1.0, -1.0}) {
        double q[3] = {0.0, y, 0.0};
        double v[2] = {sign * vx0, sign * vy0};
        double f[3];
        double e = 0.0;
        pot.evaluate(q, f, &e);
        const double sdot = std::sqrt(pot.mass) * v[0];
        den[static_cast<size_t>(block)] += std::max(sdot, 0.0);
        for (long i = 1; i <= steps; ++i) {
          for (int a = 0; a < 2; ++a) {
            v[a] += 0.5 * dt * f[a] / pot.mass;
            q[a] += dt * v[a];
          }
          pot.evaluate(q, f, &e);
          for (int a = 0; a < 2; ++a) {
            v[a] += 0.5 * dt * f[a] / pot.mass;
          }
          if (i >= first && q[0] > 0.0) {
            num[static_cast<size_t>(block)] += sdot / span;
          }
        }
      }
    }
  }
  double sumNum = 0.0, sumDen = 0.0;
  for (long b = 0; b < blocks; ++b) {
    sumNum += num[static_cast<size_t>(b)];
    sumDen += den[static_cast<size_t>(b)];
  }
  Estimate out;
  out.kappa = sumNum / sumDen;
  double mean = 0.0;
  std::vector<double> leave(static_cast<size_t>(blocks));
  for (long b = 0; b < blocks; ++b) {
    leave[static_cast<size_t>(b)] = (sumNum - num[static_cast<size_t>(b)]) /
                                    (sumDen - den[static_cast<size_t>(b)]);
    mean += leave[static_cast<size_t>(b)];
  }
  mean /= static_cast<double>(blocks);
  double var = 0.0;
  for (const double l : leave) {
    var += (l - mean) * (l - mean);
  }
  const double nb = static_cast<double>(blocks);
  out.error = std::sqrt((nb - 1.0) / nb * var);
  return out;
}

} // namespace

// One classical particle leaving the top of a one-dimensional barrier never
// returns: energy conservation keeps its kinetic energy above zero on the
// product side. Every forward child ends on the product side and every
// reversed child on the reactant side, so kappa(t) = 1 at every t.
TEST_CASE("The classical transmission at the top of a 1D barrier is one",
          "[PIQTST][recrossing]") {
  EckartPot pot;
  const auto o = recrossingAtTop(300.0, 1, 20, 10);
  const auto k = piqtst::recrossing(pot, line(false), o);
  CAPTURE(k.plateau, k.plateauError, k.trajectories);
  REQUIRE(k.trajectories == 2 * o.parents * o.children);
  REQUIRE(k.kappa.size() == static_cast<size_t>(o.steps + 1));
  REQUIRE(std::abs(k.plateau - 1.0) <= 3.0 * k.plateauError + 1e-12);
  for (const double v : k.kappa) {
    REQUIRE_THAT(v, Catch::Matchers::WithinAbs(1.0, 1e-12));
  }
}

// The barrier rotated by theta from the dividing line x = 0, with transverse
// stiffness k0 = 4 omega_b^2: the harmonic saddle gives
// kappa = sqrt(cos^2 theta - sin^2 theta omega_b^2 / k0) = 0.7755 for
// theta = 0.6. One-bead RPMD is classical MD and must agree with plain
// velocity Verlet from an independent sampler of the same line, and both
// with the harmonic value, which the Eckart anharmonicity moves by 1e-3 at
// 300 K.
TEST_CASE("One-bead RPMD recrossing matches classical trajectories on a "
          "rotated barrier",
          "[PIQTST][recrossing]") {
  const Eckart pes;
  const double t = 300.0;
  const double beta = 1.0 / (tunneling::kBoltzmann * t);
  const double omegaB2 = 2.0 * pes.v0 / (pes.a * pes.a);
  const double theta = 0.6;
  EckartPot pot;
  pot.k0 = 4.0 * omegaB2;
  pot.theta = theta;
  // Parents 40 steps apart, past the transverse period of 2.5 time units,
  // so the jackknife over parents sees nearly independent samples.
  auto o = recrossingAtTop(t, 1, 800, 20);
  o.spacing = 40;
  const auto k = piqtst::recrossing(pot, line(true), o);
  const Estimate md = classicalTransmission(pot, beta, 20000, 50, o.steps, 0.1);
  const double harmonic =
      std::sqrt(std::pow(std::cos(theta), 2) -
                std::pow(std::sin(theta), 2) * omegaB2 / pot.k0);
  const double combined = std::hypot(k.plateauError, md.error);
  CAPTURE(k.plateau, k.plateauError, md.kappa, md.error, harmonic);
  REQUIRE(k.plateauError > 0.0);
  REQUIRE(k.plateauError < 0.01);
  REQUIRE(k.plateau + 5.0 * k.plateauError < 0.9);
  REQUIRE(std::abs(k.plateau - md.kappa) < 3.0 * combined);
  REQUIRE(std::abs(k.plateau - harmonic) < 3.0 * k.plateauError + 0.005);
  REQUIRE(std::abs(md.kappa - harmonic) < 3.0 * md.error + 0.005);
  REQUIRE_THAT(k.kappa.front(), Catch::Matchers::WithinAbs(1.0, 1e-12));
  // Recrossing only lowers the flux-side correlation from its t = 0 value.
  REQUIRE(k.kappa.back() < 1.0);
}

// kappa(0) takes the side of each child from the sign of sdot(0), so the
// numerator and denominator are the same sum. A quantum ring below the
// crossover, where the internal modes and the centroid exchange energy,
// still starts at 1, and the time axis is steps of the ring time step.
TEST_CASE("The transmission curve starts at one", "[PIQTST][recrossing]") {
  const Eckart pes;
  const double tc = tunneling::crossoverTemperature(pes.hessian_at_top());
  EckartPot pot;
  pot.k0 = 1.0;
  auto o = recrossingAtTop(0.8 * tc, 8, 4, 4);
  o.steps = 20;
  const auto k = piqtst::recrossing(pot, line(true), o);
  REQUIRE(k.time.front() == 0.0);
  REQUIRE_THAT(k.time.back(), Catch::Matchers::WithinRel(20 * 0.1, 1e-12));
  REQUIRE_THAT(k.kappa.front(), Catch::Matchers::WithinAbs(1.0, 1e-14));
  for (const double v : k.kappa) {
    REQUIRE(std::isfinite(v));
  }
}

// The parents sample the ring polymer's own distribution under PILE even
// when the scan asked for PIGLET, so a PIGLET request without a GLE matrix
// still yields a transmission curve.
TEST_CASE("Recrossing parents take PILE whatever the scan's thermostat",
          "[PIQTST][recrossing]") {
  const Eckart pes;
  const double tc = tunneling::crossoverTemperature(pes.hessian_at_top());
  EckartPot pot;
  pot.k0 = 1.0;
  auto o = recrossingAtTop(0.8 * tc, 8, 4, 4);
  o.steps = 20;
  const auto pile = piqtst::recrossing(pot, line(true), o);
  o.ring.thermostat = pathintegral::Thermostat::Piglet;
  o.ring.gleFile.clear();
  const auto asked = piqtst::recrossing(pot, line(true), o);
  REQUIRE(asked.kappa.size() == pile.kappa.size());
  for (size_t i = 0; i < pile.kappa.size(); ++i) {
    REQUIRE(asked.kappa[i] == pile.kappa[i]);
  }
}

TEST_CASE("PI-QTST recrossing options round-trip and are checked",
          "[PIQTST][params]") {
  nlohmann::json j = {{"Instanton",
                       {{"mode", "rate"},
                        {"pi_planes", 9},
                        {"pi_recrossing_parents", 12},
                        {"pi_recrossing_children", 6},
                        {"pi_recrossing_time", 80.0},
                        {"pi_recrossing_spacing", 25}}}};
  Parameters p;
  eonc::config::from_json(j, p);
  const auto &o = p.instanton_options();
  REQUIRE(o.pi_recrossing_parents == 12);
  REQUIRE(o.pi_recrossing_children == 6);
  REQUIRE(o.pi_recrossing_time == Catch::Approx(80.0));
  REQUIRE(o.pi_recrossing_spacing == 25);
  Parameters again;
  eonc::config::from_json(eonc::config::to_json(p), again);
  REQUIRE(again.instanton_options().pi_recrossing_parents == 12);
  REQUIRE(again.instanton_options().pi_recrossing_time == Catch::Approx(80.0));
  REQUIRE(Parameters{}.instanton_options().pi_recrossing_parents == 0);

  auto rate = [](nlohmann::json extra) {
    nlohmann::json s = {{"mode", "rate"}, {"pi_planes", 4}};
    s.update(extra);
    return nlohmann::json{{"Instanton", s}};
  };
  for (const nlohmann::json &bad :
       {rate({{"pi_recrossing_parents", 1}}),
        rate({{"pi_recrossing_parents", -2}}),
        rate({{"pi_recrossing_parents", 4}, {"pi_recrossing_children", 0}}),
        rate({{"pi_recrossing_parents", 4}, {"pi_recrossing_spacing", 0}}),
        rate({{"pi_recrossing_parents", 4}, {"pi_recrossing_time", 1.0}})}) {
    Parameters q;
    CAPTURE(bad.dump());
    REQUIRE_THROWS_AS(eonc::config::from_json(bad, q), std::invalid_argument);
  }
}

// Ring-polymer MD on the H + H2 model barrier of Craig and Manolopoulos
// (J. Chem. Phys. 122, 084106 (2005)): V0 = 0.425 eV, a = 0.734 bohr,
// m = 1061 electron masses, so a = 0.296329 amu^0.5 Angstrom in mass-weighted
// coordinates and T_c = 371.5 K (kB beta_c = 2.69e-3 / K). The RPMD rate is
// the PI-QTST rate on the centroid plane through the top times the
// recrossing plateau; the reference is the exact quantum flux. RPMD lies
// below the exact rate for a symmetric barrier (Richardson and Althorpe,
// J. Chem. Phys. 131, 214106 (2009)); for the Eckart barrier with T_c =
// 239 K the deviation runs from -10 percent at 1.6 T_c to -45 percent at
// 0.53 T_c (Suleimanov, Aoiz and Guo, J. Phys. Chem. A 120, 8488 (2016),
// table 1). A validation run, minutes long, not part of the suite:
//   test_piqtst "[validation]"
TEST_CASE("RPMD rates of the Craig-Manolopoulos Eckart barrier",
          "[.][validation][PIQTST]") {
  const double bohr = 0.529177210544;
  const double me = 5.485799090441e-4;
  const Eckart pes{0.425, 0.734 * bohr * std::sqrt(1061.0 * me)};
  const double tc = tunneling::crossoverTemperature(pes.hessian_at_top());
  REQUIRE_THAT(1.0 / tc, Catch::Matchers::WithinRel(2.69e-3, 2e-3));
  EckartPot pot;
  pot.pes = pes;
  struct Case {
    double kBeta; // 1e-3 / K
    long beads;
  };
  // Every temperature, then 0.54 T_c (kB beta 5e-3 / K) at 16, 32, 64 and
  // 128 beads for the convergence in N.
  for (const Case cs :
       {Case{2.0, 16}, Case{3.0, 32}, Case{5.0, 48}, Case{7.0, 64},
        Case{5.0, 16}, Case{5.0, 32}, Case{5.0, 64}, Case{5.0, 128}}) {
    const double t = 1e3 / cs.kBeta;
    const double beta = 1.0 / (tunneling::kBoltzmann * t);
    piqtst::ScanOptions so;
    so.planes = planes(pes.a, 40);
    so.equilibration = 200;
    so.production = 4000;
    so.blocks = 10;
    so.ring = ring(t, cs.beads, 0.05);
    so.ring.pileScale = 0.5;
    const auto scan = piqtst::scan(pot, line(false), so);
    const double barrier = scan.back().freeEnergy;
    const double error = scan.back().freeEnergyError;
    auto ro = recrossingAtTop(t, cs.beads, 200, 10);
    ro.ring.dt = 0.05;
    ro.steps = 200;
    const auto k = piqtst::recrossing(pot, line(false), ro);
    const double exact = pes.logExactFlux(beta);
    const double qtst = std::exp(logFlux(barrier, beta) - exact);
    const double rpmd = qtst * k.plateau;
    const double rpmdError =
        rpmd * std::hypot(beta * error, k.plateauError / k.plateau);
    WARN("T = " << t << " K (T / T_c = " << t / tc << "), N = " << cs.beads
                << ": F = " << barrier << " +- " << error
                << " eV, PI-QTST / exact = " << qtst
                << ", kappa = " << k.plateau << " +- " << k.plateauError
                << ", RPMD / exact = " << rpmd << " +- " << rpmdError);
    REQUIRE(k.plateau <= 1.0 + 3.0 * k.plateauError);
    REQUIRE(rpmd < 1.0);
    REQUIRE(rpmd > 0.3);
  }
}
