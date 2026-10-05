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
#include "eon/Matter.h"

#include <thread>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("EAM_AL potential returns finite energy on Al FCC cluster",
          "[pot][eam][al]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::EAM_AL;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  // SVN reference (data/reference/point_eam_al.dat):
  double energy = matter->getPotentialEnergy();
  REQUIRE(energy == Catch::Approx(-5.217864).epsilon(1e-4));

  double maxForce = matter->getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(0.968647).epsilon(1e-3));
}

TEST_CASE("EAM_AL needs per-image instances, which run on threads",
          "[pot][eam][al][thread_safety]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::EAM_AL;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));

  // One instance is not shared across threads; each image gets its own,
  // and each instance owns its neighbour workspace.
  REQUIRE_FALSE(pot->isSharedInstanceThreadSafe());
  REQUIRE(pot->needsPerImageInstance());

  auto instance = [&]() -> std::shared_ptr<Potential> {
    auto cloned = pot->clonePotential();
    return cloned ? cloned
                  : eonc::helpers::sharePotential(
                        eonc::helpers::makePotential(params));
  };
  Matter a(instance(), params);
  Matter b(instance(), params);
  a.con2matter(std::string("pos.con"));
  b.con2matter(std::string("pos.con"));
  REQUIRE(a.getPotential().get() != b.getPotential().get());
  // Two geometries, so the instances hold different neighbour tables.
  AtomMatrix moved = b.getPositions();
  moved(0, 0) += 0.05;
  b.setPositions(moved);

  // Serial references, each on its own fresh instance.
  Matter refA(instance(), params);
  Matter refB(instance(), params);
  refA.con2matter(std::string("pos.con"));
  refB.con2matter(std::string("pos.con"));
  refB.setPositions(moved);
  const AtomMatrix fA = refA.getForces();
  const AtomMatrix fB = refB.getForces();
  const double eA = refA.getPotentialEnergy();
  const double eB = refB.getPotentialEnergy();

  // a and b are never evaluated, so each round's copies are, and a copy
  // keeps its source's instance: both instances run at once, every round,
  // and must give the serial values bit for bit.
  REQUIRE(a.needsForceUpdate());
  REQUIRE(b.needsForceUpdate());
  for (int round = 0; round < 20; ++round) {
    Matter ta(a);
    Matter tb(b);
    AtomMatrix ra;
    double ea = 0.0;
    std::thread worker([&] {
      ra = ta.getForces();
      ea = ta.getPotentialEnergy();
    });
    const AtomMatrix rb = tb.getForces();
    const double eb = tb.getPotentialEnergy();
    worker.join();
    REQUIRE(ea == eA);
    REQUIRE(eb == eB);
    REQUIRE((ra.array() == fA.array()).all());
    REQUIRE((rb.array() == fB.array()).all());
  }
}

TEST_CASE("EAM_AL minimization converges on Al FCC",
          "[pot][eam][al][minimization]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::EAM_AL;
  ParametersLoadAccess::optimizer_options(params).method = OptType::LBFGS;
  ParametersLoadAccess::optimizer_options(params).converged_force = 0.01;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 50;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  double e_before = matter->getPotentialEnergy();
  bool converged = matter->relax(false, false, false, "eam_test", "eam_test");
  double e_after = matter->getPotentialEnergy();

  // SVN reference: minimized energy = -5.562432, 7 force calls
  REQUIRE(e_after == Catch::Approx(-5.562432).epsilon(1e-4));
  REQUIRE(e_after <= e_before + 1e-10);
}

TEST_CASE("EAM_AL FIRE minimization matches SVN",
          "[pot][eam][al][minimization][fire]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::EAM_AL;
  ParametersLoadAccess::optimizer_options(params).method = OptType::FIRE;
  ParametersLoadAccess::optimizer_options(params).converged_force = 0.01;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 200;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  matter->relax(false, false, false, "eam_fire_test", "eam_fire_test");

  // SVN reference (data/reference/minimization_eam_fire.dat):
  // energy = -5.562431, 26 force calls
  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-5.562431).epsilon(1e-4));
}

} /* namespace tests */
