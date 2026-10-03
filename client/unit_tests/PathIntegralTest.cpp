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

#include "eon/PathIntegral.h"

#include "catch2/catch_amalgamated.hpp"

#include <Eigen/Core>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <string>
#include <vector>

using eonc::pathintegral::Options;
using eonc::pathintegral::RingPolymer;
using eonc::pathintegral::Springs;
using eonc::pathintegral::Thermostat;

namespace {

constexpr double kX = 20.0;

double exactKinetic(double x) { return 0.25 * x / std::tanh(0.5 * x); }

double trotterKinetic(long beads, double x) {
  double sum = 0.0;
  for (long k = 0; k < beads; ++k) {
    const double ratio =
        k == 0 ? 0.0
               : 2.0 * (static_cast<double>(beads) / x) *
                     std::sin(std::numbers::pi * static_cast<double>(k) /
                              static_cast<double>(beads));
    sum += 1.0 / (1.0 + ratio * ratio);
  }
  return 0.5 * sum;
}

struct Harmonic final : eonc::Potential {
  double mass{1.0};
  double omega{kX};
  bool sawForce{false};
  long calls{0};
  long minSystems{-1};
  long maxSystems{-1};

  Harmonic()
      : Potential(eonc::PotType::LJ) {}

  void force(long nAtoms, const double * /*positions*/,
             const int * /*atomicNrs*/, double *forces, double *energy,
             double *variance, const double * /*box*/) override {
    sawForce = true;
    *energy = 0.0;
    if (variance != nullptr) {
      *variance = 0.0;
    }
    for (long i = 0; i < nAtoms * 3; ++i) {
      forces[i] = 0.0;
    }
  }

  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return true;
  }

  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const * /*atomicNrs*/, double *const *forces,
                  double *energies, double *variances,
                  const double *const * /*boxes*/) override {
    REQUIRE(nAtoms == 1);
    ++calls;
    if (minSystems < 0 || nSystems < minSystems) {
      minSystems = nSystems;
    }
    if (nSystems > maxSystems) {
      maxSystems = nSystems;
    }
    const double spring = mass * omega * omega;
    for (long s = 0; s < nSystems; ++s) {
      const double x = positions[s][0];
      energies[s] = 0.5 * spring * x * x;
      forces[s][0] = -spring * x;
      forces[s][1] = 0.0;
      forces[s][2] = 0.0;
      if (variances != nullptr) {
        variances[s] = 0.0;
      }
    }
  }
};

Options baseOptions(long beads) {
  Options opt;
  opt.beads = beads;
  opt.temperature = 1.0;
  opt.kB = 1.0;
  opt.hbar = 1.0;
  opt.dt = 0.005;
  opt.pileTau = 0.2;
  opt.pileScale = 1.0;
  opt.seed = 11;
  return opt;
}

RingPolymer makePolymer(const Options &opt) {
  return RingPolymer(1, std::vector<double>{1.0}, std::vector<int>{1},
                     std::vector<char>{1, 0, 0}, opt);
}

void writeGle(const std::filesystem::path &path, long nModes, double tk) {
  std::ofstream out(path);
  REQUIRE(out);
  out.setf(std::ios::scientific);
  out.precision(16);
  out << nModes << " 2\n";
  for (long mode = 0; mode < nModes; ++mode) {
    out << "30 8\n-8 30\n";
    out << tk << " 0\n0 " << tk << "\n";
  }
}

} // namespace

TEST_CASE("Normal mode matrix is orthogonal", "[path-integral]") {
  for (long n : {1, 2, 7, 8}) {
    const Eigen::MatrixXd c = eonc::pathintegral::normalModeMatrix(n);
    const Eigen::MatrixXd eye = c * c.transpose();
    for (long i = 0; i < n; ++i) {
      for (long j = 0; j < n; ++j) {
        const double expect = i == j ? 1.0 : 0.0;
        REQUIRE(eye(i, j) == Catch::Approx(expect).margin(1e-12));
      }
    }
  }
}

TEST_CASE("Economised eigenvalues match the fit", "[path-integral]") {
  const Eigen::VectorXd small = eonc::pathintegral::ecoEigenvalues(8, 20.0);
  REQUIRE(small.size() == 8);
  REQUIRE(small[0] == Catch::Approx(0.0).margin(1e-15));
  for (long k = 1; k < 8; ++k) {
    REQUIRE(small[k] == Catch::Approx(1.126555635687).epsilon(1e-6));
  }
  const Eigen::VectorXd wide = eonc::pathintegral::ecoEigenvalues(48, 20.0);
  REQUIRE(wide[1] == Catch::Approx(0.130893053639).epsilon(1e-5));
  // Stationary Nyquist sample of the fit at 48 beads and xmax 20.
  REQUIRE(wide[24] == Catch::Approx(1.268436246082).epsilon(1e-5));
  REQUIRE(wide[47] == Catch::Approx(wide[1]).margin(1e-12));
}

TEST_CASE("Economised springs are refused with PIGLET and the instanton",
          "[path-integral]") {
  Options opt = baseOptions(4);
  opt.springs = Springs::Eco;
  opt.thermostat = Thermostat::Piglet;
  opt.ecoOmegaMax = 20.0;
  opt.gleFile = "unused";
  REQUIRE_THROWS_AS(makePolymer(opt), std::invalid_argument);
  REQUIRE_THROWS_AS(
      eonc::pathintegral::requireTrotterSprings("eco", "instanton"),
      std::invalid_argument);
  REQUIRE_NOTHROW(
      eonc::pathintegral::requireTrotterSprings("trotter", "instanton"));
}

TEST_CASE("Centroid hyperplane mean force of a harmonic oscillator",
          "[path-integral]") {
  Options opt = baseOptions(4);
  opt.seed = 2;
  Harmonic pot;
  pot.omega = 2.0;
  RingPolymer ring = makePolymer(opt);
  const double s = 0.25;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd origin = q;
  Eigen::VectorXd normal = q;
  q[0] = s;
  origin[0] = s;
  normal[0] = 1.0;
  ring.setAllBeads(q.data());
  ring.setHyperplane(normal, origin);
  const auto sample = ring.sample(pot, nullptr, 0, 8);
  REQUIRE(sample.meanForce == Catch::Approx(-1.0).margin(1e-8));
  // One bead batch per step, plus the first step's.
  REQUIRE(sample.batches == 9);
  REQUIRE(pot.calls == 9);
  REQUIRE(pot.minSystems == 4);
  REQUIRE(pot.maxSystems == 4);
  REQUIRE_FALSE(pot.sawForce);
  REQUIRE(ring.centroid()[0] == Catch::Approx(s).margin(1e-10));
}

TEST_CASE(
    "PIGLET kinetic energy of a harmonic oscillator at beta hbar omega 20",
    "[path-integral]") {
  const double exact = exactKinetic(kX);
  const long beads = 6;
  double internal = 0.0;
  for (long k = 1; k < beads; ++k) {
    const double ratio = 2.0 * (static_cast<double>(beads) / kX) *
                         std::sin(std::numbers::pi * static_cast<double>(k) /
                                  static_cast<double>(beads));
    internal += 1.0 / (1.0 + ratio * ratio);
  }
  const double alpha = (exact - 0.5) / (0.5 * internal);
  const auto path =
      std::filesystem::temp_directory_path() / "eon-path-integral-gle.txt";
  writeGle(path, beads - 1, alpha * static_cast<double>(beads));

  Options opt = baseOptions(beads);
  opt.thermostat = Thermostat::Piglet;
  opt.gleFile = path.string();
  opt.dt = 0.001;
  opt.seed = 5;
  Harmonic pot;
  RingPolymer ring = makePolymer(opt);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(3);
  ring.setAllBeads(q.data());
  const long production = 6000000;
  const auto sample = ring.sample(pot, nullptr, 200000, production);
  REQUIRE(sample.kineticCv == Catch::Approx(exact).epsilon(0.01));
  REQUIRE(sample.batches == production);
  REQUIRE(pot.minSystems == beads);
  REQUIRE(pot.maxSystems == beads);
  REQUIRE_FALSE(pot.sawForce);
  std::filesystem::remove(path);
}

TEST_CASE("Economised springs reach the Trotter error at half the beads",
          "[path-integral]") {
  const double exact = exactKinetic(kX);
  long trotterBeads = 0;
  for (long beads : {16L, 32L, 48L, 64L, 96L, 128L}) {
    const double rel = std::fabs(trotterKinetic(beads, kX) - exact) / exact;
    if (rel < 0.01) {
      trotterBeads = beads;
      break;
    }
  }
  REQUIRE(trotterBeads == 96);
  const long ecoBeads = trotterBeads / 2;
  REQUIRE(std::fabs(trotterKinetic(ecoBeads, kX) - exact) / exact > 0.01);

  Options opt = baseOptions(ecoBeads);
  opt.springs = Springs::Eco;
  opt.ecoOmegaMax = kX;
  opt.dt = 0.002;
  opt.seed = 9;
  Harmonic pot;
  RingPolymer ring = makePolymer(opt);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(3);
  ring.setAllBeads(q.data());
  const auto sample = ring.sample(pot, nullptr, 50000, 800000);
  REQUIRE(sample.kineticCv == Catch::Approx(exact).epsilon(0.01));
  REQUIRE(pot.minSystems == ecoBeads);
  REQUIRE(pot.maxSystems == ecoBeads);
  REQUIRE_FALSE(pot.sawForce);
}
