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

#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/HelperFunctions.h"
#include "eon/ImprovedDimer.h"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/Parameters.h"

#include <atomic>
#include <chrono>
#include <cmath>
#include <memory>
#include <thread>
#include <vector>

namespace tests {

using eonc::ImprovedDimer;
using eonc::LowestEigenmode;
using eonc::Matter;
using eonc::OptType;
using eonc::Parameters;
using eonc::ParametersLoadAccess;
using eonc::PotType;

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("ImprovedDimer computes eigenvalue on displaced cluster",
          "[ImprovedDimer]") {
  Parameters parameters;
  ParametersLoadAccess::potential_options(parameters).potential = PotType::LJ;
  ParametersLoadAccess::optimizer_options(parameters).method = OptType::CG;
  ParametersLoadAccess::optimizer_options(parameters).converged_force = 0.001;
  ParametersLoadAccess::dimer_options(parameters).converged_angle = 0.001;
  ParametersLoadAccess::saddle_search_options(parameters).minmode_method =
      LowestEigenmode::MINMODE_DIMER;

  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(parameters));
  auto initial = std::make_shared<Matter>(pot, parameters);
  auto displacement = std::make_shared<Matter>(pot, parameters);
  auto saddle = std::make_shared<Matter>(pot, parameters);

  initial->con2matter(std::string("pos.con"));
  saddle->con2matter(std::string("displacement.con"));

  AtomMatrix mode =
      eonc::helpers::loadMode("direction.dat", initial->numberOfAtoms());

  auto minModeMethod = std::make_unique<ImprovedDimer>(saddle, parameters,
                                                       saddle->getPotential());
  minModeMethod->compute(saddle, mode);

  double eigenvalue = minModeMethod->getEigenvalue();
  REQUIRE(std::isfinite(eigenvalue));
  // On a displaced cluster the dimer should find a negative eigenvalue
  REQUIRE(eigenvalue < 0.0);
}

// A per-image backend (like ext_pot) that cannot be shared across threads.
// Each force() call records whether another thread was already inside the
// same instance.
class ExclusivePerImagePot final : public eonc::Potential {
public:
  struct Log {
    std::atomic<int> overlaps{0};
    std::atomic<int> instancesUsed{0};
  };
  explicit ExclusivePerImagePot(std::shared_ptr<Log> log)
      : Potential(PotType::UNKNOWN),
        log_(std::move(log)) {}

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }
  [[nodiscard]] std::shared_ptr<eonc::Potential>
  clonePotential() const override {
    return std::make_shared<ExclusivePerImagePot>(log_);
  }

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    static_cast<void>(atomicNrs);
    static_cast<void>(box);
    if (inside_.fetch_add(1) != 0) {
      log_->overlaps.fetch_add(1);
    }
    if (!used_.exchange(true)) {
      log_->instancesUsed.fetch_add(1);
    }
    // Long enough for the two dimer images to meet if they share this
    // instance.
    std::this_thread::sleep_for(std::chrono::milliseconds(2));
    const double curvature[3] = {-10.0, 0.1, 5.0};
    double e = 0.0;
    for (long a = 0; a < nAtoms; ++a) {
      for (int c = 0; c < 3; ++c) {
        const double x = positions[3 * a + c];
        forces[3 * a + c] = -curvature[c] * x;
        e += 0.5 * curvature[c] * x * x;
      }
    }
    *energy = e;
    if (variance != nullptr) {
      *variance = 0.0;
    }
    inside_.fetch_sub(1);
  }

private:
  std::shared_ptr<Log> log_;
  std::atomic<int> inside_{0};
  std::atomic<bool> used_{false};
};

TEST_CASE("ImprovedDimer evaluates per-image potentials on separate instances",
          "[ImprovedDimer][parallel]") {
  Parameters parameters;
  REQUIRE(parameters.main_options().parallel);
  auto log = std::make_shared<ExclusivePerImagePot::Log>();
  auto pot = std::make_shared<ExclusivePerImagePot>(log);
  auto matter = std::make_shared<Matter>(pot, parameters);
  matter->resize(1);
  matter->setAtomicNr(0, 1);
  matter->setMass(0, 1.0);
  matter->setPeriodic(false);
  matter->setCell(Matrix3d::Identity() * 20.0);
  AtomMatrix start(1, 3);
  start << 0.1, 0.1, 0.1;
  matter->setPositions(start);

  AtomMatrix mode(1, 3);
  mode << 1.0, 1.0, 1.0;
  mode /= std::sqrt(3.0);

  ImprovedDimer dimer(matter, parameters, pot);
  dimer.compute(matter, mode);

  // The two images ran at the same time on two instances, never on one.
  REQUIRE(log->instancesUsed.load() >= 2);
  REQUIRE(log->overlaps.load() == 0);
  REQUIRE(dimer.getEigenvalue() == Catch::Approx(-10.0).margin(0.5));
}

// Diagonal Hessian, E = 0.5 (c_x x^2 + c_y y^2 + c_z z^2), one or more atoms.
class DiagHessianPot final : public eonc::Potential {
public:
  DiagHessianPot(double cx, double cy, double cz)
      : Potential(PotType::UNKNOWN),
        c_{cx, cy, cz} {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    static_cast<void>(atomicNrs);
    static_cast<void>(box);
    double e = 0.0;
    for (long a = 0; a < nAtoms; ++a) {
      for (int k = 0; k < 3; ++k) {
        const double x = positions[3 * a + k];
        forces[3 * a + k] = -c_[k] * x;
        e += 0.5 * c_[k] * x * x;
      }
    }
    *energy = e;
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }

private:
  double c_[3];
};

static std::shared_ptr<Matter>
oneAtom(const std::shared_ptr<eonc::Potential> &pot,
        const Parameters &parameters, const Eigen::Vector3d &r) {
  auto matter = std::make_shared<Matter>(pot, parameters);
  matter->resize(1);
  matter->setAtomicNr(0, 1);
  matter->setMass(0, 1.0);
  matter->setPeriodic(false);
  matter->setCell(Matrix3d::Identity() * 20.0);
  AtomMatrix start(1, 3);
  start.row(0) = r.transpose();
  matter->setPositions(start);
  return matter;
}

TEST_CASE("ImprovedDimer CG direction persists across translations",
          "[ImprovedDimer][cg]") {
  // EON r2549 keeps init_cg, F_R_Old and thetaOld on the dimer object.
  // A per-compute reset forces gamma = 0 on the single rotation a loose
  // angle tolerance allows, which is steepest descent.
  Parameters parameters;
  ParametersLoadAccess::main_options(parameters).parallel = false;
  ParametersLoadAccess::main_options(parameters).finiteDifference = 0.01;
  ParametersLoadAccess::dimer_options(parameters).opt_method = OptType::CG;
  ParametersLoadAccess::dimer_options(parameters).improved = true;
  ParametersLoadAccess::dimer_options(parameters).converged_angle = 180.0;
  ParametersLoadAccess::dimer_options(parameters).rotations_max = 1;

  auto pot = std::make_shared<DiagHessianPot>(-8.0, 3.0, 6.0);
  auto matter = oneAtom(pot, parameters, Eigen::Vector3d::Zero());
  ImprovedDimer dimer(matter, parameters, pot);

  AtomMatrix mode(1, 3);
  mode << 1.0, -0.4, 0.7;
  dimer.compute(matter, mode);
  REQUIRE_FALSE(dimer.init_cg);
  REQUIRE(dimer.gamma == Catch::Approx(0.0));
  const VectorXd thetaBefore = dimer.thetaOld;
  const VectorXd frBefore = dimer.F_R_Old;
  REQUIRE(frBefore.norm() > 1.0e-6);

  mode << 0.5, 0.2, -1.0;
  dimer.compute(matter, mode);

  const double b = frBefore.squaredNorm();
  const double overlap = std::abs(dimer.F_R.dot(frBefore));
  REQUIRE(overlap < 0.5 * b);
  const double expectGamma = dimer.F_R.dot(dimer.F_R - frBefore) / b;
  REQUIRE(std::abs(expectGamma) > 1.0e-3);
  REQUIRE(dimer.gamma == Catch::Approx(expectGamma).epsilon(1e-12));
  REQUIRE_FALSE(dimer.init_cg);

  // theta = F_R + thetaOld * gamma, then drop the component along tau.
  // thetaOld is that vector before it is normalised (r2549). The r2006
  // form multiplies gamma by |F_R_Old| and does not satisfy this.
  VectorXd expectTheta = dimer.F_R + thetaBefore * dimer.gamma;
  expectTheta -= expectTheta.dot(dimer.tau) * dimer.tau;
  REQUIRE((dimer.thetaOld - expectTheta).norm() ==
          Catch::Approx(0.0).margin(1e-10));
}

TEST_CASE("ImprovedDimer CG direction persists after a multi-rotation estimate",
          "[ImprovedDimer][cg]") {
  // An estimate that rotated more than once still leaves F_R_Old and
  // thetaOld for the next centre. The following rotation uses them, so
  // gamma is the Polak-Ribiere value rather than the steepest-descent 0
  // a per-compute reset would force.
  Parameters parameters;
  ParametersLoadAccess::main_options(parameters).parallel = false;
  ParametersLoadAccess::main_options(parameters).finiteDifference = 0.01;
  ParametersLoadAccess::dimer_options(parameters).opt_method = OptType::CG;
  ParametersLoadAccess::dimer_options(parameters).improved = true;
  ParametersLoadAccess::dimer_options(parameters).converged_angle = 1.0e-6;
  ParametersLoadAccess::dimer_options(parameters).rotations_max = 2;

  auto pot = std::make_shared<DiagHessianPot>(-8.0, 3.0, 6.0);
  auto matter = oneAtom(pot, parameters, Eigen::Vector3d::Zero());
  ImprovedDimer dimer(matter, parameters, pot);

  AtomMatrix mode(1, 3);
  mode << 1.0, 1.0, 0.3;
  dimer.compute(matter, mode);
  REQUIRE(dimer.statsRotations > 1);
  REQUIRE_FALSE(dimer.init_cg);

  // The converged rotational force is tiny. Plant a direction the next
  // F_R does not already contain, so a reset (gamma = 0) cannot pass.
  dimer.F_R_Old = VectorXd::Zero(3);
  dimer.F_R_Old << 0.0, 4.0, -1.0;
  dimer.thetaOld = VectorXd::Zero(3);
  dimer.thetaOld << 0.2, -0.5, 0.1;
  const VectorXd thetaBefore = dimer.thetaOld;
  const VectorXd frBefore = dimer.F_R_Old;

  ParametersLoadAccess::dimer_options(parameters).converged_angle = 180.0;
  ParametersLoadAccess::dimer_options(parameters).rotations_max = 1;
  mode << 0.5, 0.2, -1.0;
  dimer.compute(matter, mode);

  const double b = frBefore.squaredNorm();
  const double overlap = std::abs(dimer.F_R.dot(frBefore));
  REQUIRE(overlap < 0.5 * b);
  const double expectGamma = dimer.F_R.dot(dimer.F_R - frBefore) / b;
  REQUIRE(std::abs(expectGamma) > 1.0e-3);
  REQUIRE(dimer.gamma == Catch::Approx(expectGamma).epsilon(1e-12));
  REQUIRE_FALSE(dimer.init_cg);
  VectorXd expectTheta = dimer.F_R + thetaBefore * dimer.gamma;
  expectTheta -= expectTheta.dot(dimer.tau) * dimer.tau;
  REQUIRE((dimer.thetaOld - expectTheta).norm() ==
          Catch::Approx(0.0).margin(1e-10));
}

static void prepareMinMode(Parameters &parameters) {
  ParametersLoadAccess::main_options(parameters).parallel = false;
  ParametersLoadAccess::main_options(parameters).finiteDifference = 0.01;
  ParametersLoadAccess::main_options(parameters).removeNetForce = false;
  ParametersLoadAccess::optimizer_options(parameters).method = OptType::LBFGS;
  ParametersLoadAccess::optimizer_options(parameters).max_move = 0.2;
  ParametersLoadAccess::optimizer_options(parameters).lbfgs.auto_scale = true;
  ParametersLoadAccess::dimer_options(parameters).opt_method = OptType::CG;
  ParametersLoadAccess::dimer_options(parameters).improved = true;
  ParametersLoadAccess::dimer_options(parameters).converged_angle = 180.0;
  ParametersLoadAccess::dimer_options(parameters).rotations_max = 1;
  ParametersLoadAccess::saddle_search_options(parameters).minmode_method =
      LowestEigenmode::MINMODE_DIMER;
  ParametersLoadAccess::saddle_search_options(parameters).converged_force =
      1.0e-12;
  ParametersLoadAccess::saddle_search_options(parameters).max_iterations = 1;
}

TEST_CASE("Steep mode curvature yields to the effective-force probe",
          "[ImprovedDimer][lbfgs]") {
  // Eigenvector of H = diag(-40, 100, 100), so |C_τ| = 40. H0 = 1/40
  // times |F_eff| is past max_move. The finite-difference curvature of
  // F_eff is stiffer, and the accepted step is that shorter step.
  // A capped |C_τ| step would have length 0.2 and would skip the probe.
  Parameters parameters;
  prepareMinMode(parameters);
  const Eigen::Vector3d r0(0.05, 0.1, 0.1);
  const Eigen::Vector3d tau(1.0, 0.0, 0.0);
  const Eigen::Vector3d feff(-40.0 * r0.x(), -100.0 * r0.y(), -100.0 * r0.z());
  const Eigen::Vector3d fhat = feff.normalized();
  const double cFd = 40.0 * fhat.x() * fhat.x() + 100.0 * fhat.y() * fhat.y() +
                     100.0 * fhat.z() * fhat.z();
  const Eigen::Vector3d expect = feff / cFd;
  REQUIRE(expect.norm() < 0.2);
  REQUIRE((feff / 40.0).norm() > 0.2);

  auto pot = std::make_shared<DiagHessianPot>(-40.0, 100.0, 100.0);
  auto matter = oneAtom(pot, parameters, r0);
  AtomMatrix mode(1, 3);
  mode.row(0) = tau.transpose();
  const long before = static_cast<long>(pot->forceCallCounter.load());
  eonc::MinModeSaddleSearch search(matter, mode, 0.0, parameters, pot);
  const int status = search.run(1);
  const long calls = static_cast<long>(pot->forceCallCounter.load()) - before;
  const Eigen::Vector3d step = matter->getPositions().row(0).transpose() - r0;
  INFO("status " << status << " force calls " << calls << " step "
                 << step.norm());
  REQUIRE(step.norm() == Catch::Approx(expect.norm()).epsilon(1e-6));
  REQUIRE(step.normalized().dot(fhat) == Catch::Approx(1.0).margin(1e-6));
  // Centre, forward image, the probe, and the accepted step.
  REQUIRE(calls == 4);
}

TEST_CASE("Mild mode curvature yields to the effective-force probe",
          "[ImprovedDimer][lbfgs]") {
  // |C_τ| = 2 is mild, and H0 = 1/|C_τ| times |F_eff| is still past
  // max_move. The finite-difference curvature of F_eff is used anyway:
  // the probe is not reserved for steep modes.
  Parameters parameters;
  prepareMinMode(parameters);
  const Eigen::Vector3d r0(0.01, 0.05, 0.05);
  const Eigen::Vector3d tau(1.0, 0.0, 0.0);
  const Eigen::Vector3d feff(-2.0 * r0.x(), -50.0 * r0.y(), -50.0 * r0.z());
  const Eigen::Vector3d fhat = feff.normalized();
  const double cFd = 2.0 * fhat.x() * fhat.x() + 50.0 * fhat.y() * fhat.y() +
                     50.0 * fhat.z() * fhat.z();
  const Eigen::Vector3d expect = feff / cFd;
  REQUIRE((feff / 2.0).norm() > 0.2);
  REQUIRE(expect.norm() < 0.2);

  auto pot = std::make_shared<DiagHessianPot>(-2.0, 50.0, 50.0);
  auto matter = oneAtom(pot, parameters, r0);
  AtomMatrix mode(1, 3);
  mode.row(0) = tau.transpose();
  const long before = static_cast<long>(pot->forceCallCounter.load());
  eonc::MinModeSaddleSearch search(matter, mode, 0.0, parameters, pot);
  const int status = search.run(1);
  const long calls = static_cast<long>(pot->forceCallCounter.load()) - before;
  const Eigen::Vector3d step = matter->getPositions().row(0).transpose() - r0;
  INFO("status " << status << " force calls " << calls << " step "
                 << step.norm());
  REQUIRE(step.norm() == Catch::Approx(expect.norm()).epsilon(1e-6));
  REQUIRE(step.normalized().dot(fhat) == Catch::Approx(1.0).margin(1e-6));
  // Centre, forward image, the probe, and the accepted step.
  REQUIRE(calls == 4);
}

// Saddle along x, with a wall on the far side of the saddle. A long
// steepest-descent step crosses the wall. max_energy stays 20.
class FarWallPot final : public eonc::Potential {
public:
  FarWallPot()
      : Potential(PotType::UNKNOWN) {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    static_cast<void>(nAtoms);
    static_cast<void>(atomicNrs);
    static_cast<void>(box);
    const double x = positions[0];
    const double y = positions[1];
    const double z = positions[2];
    double e = -2.0 * x * x + 20.0 * y * y + 20.0 * z * z;
    double fx = 4.0 * x;
    double fy = -40.0 * y;
    double fz = -40.0 * z;
    if (x < minX) {
      minX = x;
    }
    if (x < -0.05) {
      ++belowWall;
      const double w = x + 0.05;
      e += 400.0 + 2000.0 * w * w;
      fx += -4000.0 * w;
    }
    forces[0] = fx;
    forces[1] = fy;
    forces[2] = fz;
    *energy = e;
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }

  mutable double minX{1.0e9};
  mutable int belowWall{0};
};

class AlwaysHotPot final : public eonc::Potential {
public:
  AlwaysHotPot()
      : Potential(PotType::UNKNOWN) {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    static_cast<void>(nAtoms);
    static_cast<void>(atomicNrs);
    static_cast<void>(box);
    const double x = positions[0];
    forces[0] = -x;
    forces[1] = 0.0;
    forces[2] = 0.0;
    *energy = 100.0 + 0.5 * x * x;
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }
};

static void prepareWallSearch(Parameters &parameters) {
  ParametersLoadAccess::main_options(parameters).parallel = false;
  ParametersLoadAccess::main_options(parameters).finiteDifference = 0.001;
  ParametersLoadAccess::main_options(parameters).removeNetForce = false;
  ParametersLoadAccess::optimizer_options(parameters).method = OptType::SD;
  ParametersLoadAccess::optimizer_options(parameters).max_move = 0.5;
  ParametersLoadAccess::optimizer_options(parameters).sd.alpha = 0.5;
  ParametersLoadAccess::optimizer_options(parameters).sd.two_point = false;
  ParametersLoadAccess::dimer_options(parameters).opt_method = OptType::CG;
  ParametersLoadAccess::dimer_options(parameters).improved = true;
  ParametersLoadAccess::dimer_options(parameters).converged_angle = 5.0;
  ParametersLoadAccess::dimer_options(parameters).rotations_max = 8;
  ParametersLoadAccess::dimer_options(parameters).rotations_min = 1;
  ParametersLoadAccess::saddle_search_options(parameters).minmode_method =
      LowestEigenmode::MINMODE_DIMER;
  ParametersLoadAccess::saddle_search_options(parameters).converged_force =
      1.0e-3;
  ParametersLoadAccess::saddle_search_options(parameters).max_energy = 20.0;
}

TEST_CASE("Min-mode shortens a translation that crosses max_energy",
          "[MinMode][max_energy]") {
  Parameters parameters;
  prepareWallSearch(parameters);
  auto pot = std::make_shared<FarWallPot>();
  auto matter = oneAtom(pot, parameters, Eigen::Vector3d(0.2, 0.0, 0.0));
  AtomMatrix mode(1, 3);
  mode << 1.0, 0.0, 0.0;
  eonc::MinModeSaddleSearch search(matter, mode, 0.0, parameters, pot);
  const int status = search.run(40);
  const double x = matter->getPositions()(0, 0);
  const double energy = matter->getPotentialEnergy();
  INFO("status " << status << " x " << x << " energy " << energy);
  REQUIRE(status == eonc::MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(std::abs(x) < 0.02);
  REQUIRE(energy < 20.0);
  // The first halving accepts. That probe is not in the eigenmode
  // counter, and getForceCalls() has to include it.
  REQUIRE(pot->belowWall >= 1);
  REQUIRE(search.getForceCalls() >= 1);
}

TEST_CASE("Min-mode still aborts when every point is past max_energy",
          "[MinMode][max_energy]") {
  Parameters parameters;
  prepareWallSearch(parameters);
  auto pot = std::make_shared<AlwaysHotPot>();
  auto matter = oneAtom(pot, parameters, Eigen::Vector3d(0.1, 0.0, 0.0));
  AtomMatrix mode(1, 3);
  mode << 1.0, 0.0, 0.0;
  eonc::MinModeSaddleSearch search(matter, mode, 0.0, parameters, pot);
  const int status = search.run(8);
  REQUIRE(status == eonc::MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY);
}

// Energy is under the limit only at the starting point. Every probe
// misses, so the stall counter climbs without a reset.
class CentreColdPot final : public eonc::Potential {
public:
  explicit CentreColdPot(double x0)
      : Potential(PotType::UNKNOWN),
        x0_{x0} {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    static_cast<void>(nAtoms);
    static_cast<void>(atomicNrs);
    static_cast<void>(box);
    const double x = positions[0];
    const double y = positions[1];
    const double z = positions[2];
    forces[0] = 4.0 * x + 2.0;
    forces[1] = -80.0 * y;
    forces[2] = -80.0 * z;
    const bool atCentre = std::abs(x - x0_) < 1.0e-12 &&
                          std::abs(y) < 1.0e-12 && std::abs(z) < 1.0e-12;
    *energy = atCentre ? 0.0 : 100.0;
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }

private:
  double x0_;
};

// Legal for x >= 0. Above x = 0.12 the climb points at the wall, and
// at or below 0.12 it points back into the legal region. A crossing
// therefore accepts half the step, the next step stays under the
// limit, and the cap doubles. The stall counter resets on both, so
// six consecutive stalls never happen. The curvature along x stays
// negative on each piece, and y and z are stiff.
class HalfLinePot final : public eonc::Potential {
public:
  HalfLinePot()
      : Potential(PotType::UNKNOWN) {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    static_cast<void>(nAtoms);
    static_cast<void>(atomicNrs);
    static_cast<void>(box);
    const double x = positions[0];
    const double y = positions[1];
    const double z = positions[2];
    const double push = x > 0.12 ? 2.0 : -2.0;
    forces[0] = push + 0.1 * x;
    forces[1] = -80.0 * y;
    forces[2] = -80.0 * z;
    *energy = x >= 0.0 ? 0.0 : 100.0;
    if (x >= 0.0 && x < 0.14) {
      ++shortAccepted;
    }
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }

  mutable int shortAccepted{0};
};

// Legal for x >= 0. Above x = 0.1 the climb points at the wall.
// At or below 0.1 it points back into the legal region, so the next
// step stays under the limit and the cap doubles.
class RegrowthPot final : public eonc::Potential {
public:
  RegrowthPot()
      : Potential(PotType::UNKNOWN) {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    static_cast<void>(nAtoms);
    static_cast<void>(atomicNrs);
    static_cast<void>(box);
    const double x = positions[0];
    const double y = positions[1];
    const double z = positions[2];
    const double push = x > 0.1 ? 2.0 : -2.0;
    forces[0] = push + 0.1 * x;
    forces[1] = -200.0 * y;
    forces[2] = -200.0 * z;
    *energy = x >= 0.0 ? 0.0 : 100.0;
    if (x < minX) {
      minX = x;
    }
    if (x > maxX) {
      maxX = x;
    }
    if (x < -0.05) {
      ++hot;
    }
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }

  mutable double minX{1.0e9};
  mutable double maxX{-1.0e9};
  mutable int hot{0};
};

// Saddle at x = 0, with a softening x^4 term. A curvature step from
// x = 0.2 is longer than the distance to the saddle, so it crosses
// the wall at x < -0.05. Steepest descent with a fixed step does too.
class SoftWallPot final : public eonc::Potential {
public:
  SoftWallPot()
      : Potential(PotType::UNKNOWN) {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    static_cast<void>(nAtoms);
    static_cast<void>(atomicNrs);
    static_cast<void>(box);
    const double x = positions[0];
    const double y = positions[1];
    const double z = positions[2];
    double e = -2.0 * x * x + 5.0 * x * x * x * x + 20.0 * y * y + 20.0 * z * z;
    double fx = 4.0 * x - 20.0 * x * x * x;
    double fy = -40.0 * y;
    double fz = -40.0 * z;
    if (x < minX) {
      minX = x;
    }
    if (x < -0.05) {
      const double w = x + 0.05;
      e += 400.0 + 2000.0 * w * w;
      fx += -4000.0 * w;
    }
    forces[0] = fx;
    forces[1] = fy;
    forces[2] = fz;
    *energy = e;
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }

  mutable double minX{1.0e9};
};

// Reported energy is hot for x below wallX. An x^2 y term tilts the
// soft mode, so one dimer rotation does not finish the estimate.
class XWallPot final : public eonc::Potential {
public:
  explicit XWallPot(double wallX)
      : Potential(PotType::UNKNOWN),
        wallX_{wallX} {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    static_cast<void>(nAtoms);
    static_cast<void>(atomicNrs);
    static_cast<void>(box);
    const double x = positions[0];
    const double y = positions[1];
    const double z = positions[2];
    // E = -2 x^2 + 8 y^2 + 8 z^2 + 20 x^2 y, plus the wall.
    forces[0] = 4.0 * x - 40.0 * x * y;
    forces[1] = -16.0 * y - 20.0 * x * x;
    forces[2] = -16.0 * z;
    const double mechanical =
        -2.0 * x * x + 8.0 * y * y + 8.0 * z * z + 20.0 * x * x * y;
    *energy = x < wallX_ ? 100.0 : mechanical;
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }

private:
  double wallX_;
};

TEST_CASE("Min-mode aborts after six consecutive energy stalls",
          "[MinMode][max_energy]") {
  Parameters parameters;
  prepareWallSearch(parameters);
  ParametersLoadAccess::optimizer_options(parameters).sd.alpha = 1.0;
  ParametersLoadAccess::optimizer_options(parameters).max_move = 0.2;
  auto pot = std::make_shared<CentreColdPot>(0.2);
  auto matter = oneAtom(pot, parameters, Eigen::Vector3d(0.2, 0.0, 0.0));
  AtomMatrix mode(1, 3);
  mode << 1.0, 0.0, 0.0;
  eonc::MinModeSaddleSearch search(matter, mode, 0.0, parameters, pot);
  const int status = search.run(40);
  INFO("status " << status << " iterations " << search.getIterationCount());
  REQUIRE(status == eonc::MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY);
  // Six restores, not the iteration cap.
  REQUIRE(search.getIterationCount() == 6);
  REQUIRE(search.getForceCalls() >= 6);
}

TEST_CASE("Min-mode stops a pinned search after the wall-event bound",
          "[MinMode][max_energy]") {
  Parameters parameters;
  prepareWallSearch(parameters);
  ParametersLoadAccess::optimizer_options(parameters).sd.alpha = 1.0;
  ParametersLoadAccess::optimizer_options(parameters).max_move = 0.2;
  auto pot = std::make_shared<HalfLinePot>();
  auto matter = oneAtom(pot, parameters, Eigen::Vector3d(0.35, 0.0, 0.0));
  AtomMatrix mode(1, 3);
  mode << 1.0, 0.0, 0.0;
  eonc::MinModeSaddleSearch search(matter, mode, 0.0, parameters, pot);
  const int status = search.run(200);
  INFO("status " << status << " iterations " << search.getIterationCount()
                 << " short " << pot->shortAccepted);
  REQUIRE(status == eonc::MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY);
  // One legal step, then pairs of (found probe, under-limit regrowth).
  // The 16th crossing is the 32nd step. Six stalls would stop at 6,
  // and the old path would run to the 200-iteration cap.
  REQUIRE(search.getIterationCount() == 32);
  REQUIRE(search.getIterationCount() < 200);
  REQUIRE(pot->shortAccepted >= 16);
  // Improved-dimer force accounting does not include these probes.
  REQUIRE(search.getForceCalls() >= 16);
}

TEST_CASE("Min-mode regrows the translation cap after an under-limit step",
          "[MinMode][max_energy]") {
  Parameters parameters;
  prepareWallSearch(parameters);
  ParametersLoadAccess::optimizer_options(parameters).sd.alpha = 1.0;
  ParametersLoadAccess::optimizer_options(parameters).max_move = 0.4;
  auto pot = std::make_shared<RegrowthPot>();
  auto matter = oneAtom(pot, parameters, Eigen::Vector3d(0.25, 0.0, 0.0));
  AtomMatrix mode(1, 3);
  mode << 1.0, 0.0, 0.0;
  eonc::MinModeSaddleSearch search(matter, mode, 0.0, parameters, pot);
  const int status = search.run(3);
  INFO("status " << status << " min " << pot->minX << " max " << pot->maxX
                 << " hot " << pot->hot);
  // The away step is still on the shortened cap, so x stays below 0.3.
  // The doubled cap then crosses the wall a second time.
  REQUIRE(pot->maxX < 0.35);
  REQUIRE(pot->hot >= 2);
  REQUIRE(search.getIterationCount() == 3);
  static_cast<void>(status);
}

TEST_CASE("A very short accepted step re-initialises the mode",
          "[MinMode][max_energy]") {
  const Eigen::Vector3d tau(1.0, 0.3, 0.0);
  const Eigen::Vector3d tauHat = tau.normalized();
  auto dotMode = [&](double wallX) {
    Parameters parameters;
    prepareWallSearch(parameters);
    ParametersLoadAccess::optimizer_options(parameters).sd.alpha = 1.0;
    ParametersLoadAccess::optimizer_options(parameters).max_move = 0.4;
    // CG keeps rotational memory across translations. A reseed drops it
    // and puts the original tangent back. SD would rebuild the same
    // step from the force, so the two paths would match.
    ParametersLoadAccess::dimer_options(parameters).opt_method = OptType::CG;
    ParametersLoadAccess::dimer_options(parameters).rotations_max = 1;
    ParametersLoadAccess::dimer_options(parameters).converged_angle = 1.0e-8;
    ParametersLoadAccess::saddle_search_options(parameters).converged_force =
        1.0e-12;
    auto pot = std::make_shared<XWallPot>(wallX);
    auto matter = oneAtom(pot, parameters, Eigen::Vector3d(0.5, 0.0, 0.0));
    AtomMatrix mode(1, 3);
    mode.row(0) = tauHat.transpose();
    eonc::MinModeSaddleSearch search(matter, mode, 0.0, parameters, pot);
    const int status = search.run(2);
    const AtomMatrix ev = search.getEigenvector();
    const double dot =
        ev(0, 0) * tauHat.x() + ev(0, 1) * tauHat.y() + ev(0, 2) * tauHat.z();
    INFO("wall " << wallX << " status " << status << " dot " << dot << " ev "
                 << ev(0, 0) << " " << ev(0, 1));
    REQUIRE(status != eonc::MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY);
    return std::abs(dot);
  };
  // x < 0.499 needs more than six halvings, so the mode is seeded
  // again. x < 0.2 accepts at 1/2 and keeps the rotational CG state.
  // The fresh rotation moves farther from the original tangent, toward
  // the soft mode. The stale CG step does not.
  const double reseeded = dotMode(0.499);
  const double kept = dotMode(0.2);
  REQUIRE(kept > reseeded + 0.05);
}

static void requireRecovered(OptType method, bool lanczos) {
  Parameters parameters;
  prepareWallSearch(parameters);
  ParametersLoadAccess::optimizer_options(parameters).method = method;
  if (method == OptType::FIRE) {
    // |F| is about 0.6. dt^2 times that force, capped by max_move,
    // has to pass the saddle and the wall.
    ParametersLoadAccess::optimizer_options(parameters).time_step = 1.0;
    ParametersLoadAccess::optimizer_options(parameters).max_time_step = 2.0;
  }
  if (lanczos) {
    ParametersLoadAccess::saddle_search_options(parameters).minmode_method =
        LowestEigenmode::MINMODE_LANCZOS;
    ParametersLoadAccess::lanczos_options(parameters).max_iterations = 6;
    ParametersLoadAccess::lanczos_options(parameters).tolerance = 0.05;
  }
  // The softening well makes a curvature step overshoot the saddle.
  // Lanczos reports linear dependence on that well, so it climbs the
  // quadratic wall instead. Steepest descent still overshoots there.
  std::shared_ptr<eonc::Potential> pot;
  double *minX = nullptr;
  if (lanczos) {
    auto far = std::make_shared<FarWallPot>();
    minX = &far->minX;
    pot = far;
  } else {
    auto soft = std::make_shared<SoftWallPot>();
    minX = &soft->minX;
    pot = soft;
  }
  auto matter = oneAtom(pot, parameters, Eigen::Vector3d(0.2, 0.0, 0.0));
  AtomMatrix mode(1, 3);
  mode << 1.0, 0.0, 0.0;
  eonc::MinModeSaddleSearch search(matter, mode, 0.0, parameters, pot);
  const int status = search.run(40);
  const double x = matter->getPositions()(0, 0);
  const double energy = matter->getPotentialEnergy();
  INFO("method " << static_cast<int>(method) << " lanczos " << lanczos
                 << " status " << status << " x " << x << " energy " << energy
                 << " min " << *minX << " calls " << search.getForceCalls());
  REQUIRE(status == eonc::MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(energy < 20.0);
  REQUIRE(*minX < -0.05);
  REQUIRE(search.getForceCalls() >= 1);
}

TEST_CASE("LBFGS translation recovers a step that crosses max_energy",
          "[MinMode][max_energy]") {
  requireRecovered(OptType::LBFGS, false);
}

TEST_CASE("CG translation recovers a step that crosses max_energy",
          "[MinMode][max_energy]") {
  requireRecovered(OptType::CG, false);
}

TEST_CASE("FIRE translation recovers a step that crosses max_energy",
          "[MinMode][max_energy]") {
  requireRecovered(OptType::FIRE, false);
}

TEST_CASE("Lanczos translation recovers a step that crosses max_energy",
          "[MinMode][max_energy]") {
  requireRecovered(OptType::SD, true);
}

} /* namespace tests */
