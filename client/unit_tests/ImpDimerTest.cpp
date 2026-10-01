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
    (void)atomicNrs;
    (void)box;
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

} /* namespace tests */
