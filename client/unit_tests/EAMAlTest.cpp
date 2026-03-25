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

#include "Matter.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("EAM_AL potential returns finite energy on Al FCC cluster",
          "[pot][eam][al]") {
  Parameters params;
  params.potential_options.potential = PotType::EAM_AL;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  // SVN reference (data/reference/point_eam_al.dat):
  double energy = matter->getPotentialEnergy();
  REQUIRE(energy == Catch::Approx(-5.217864).epsilon(1e-4));

  double maxForce = matter->getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(0.968647).epsilon(1e-3));
}

TEST_CASE("EAM_AL minimization converges on Al FCC",
          "[pot][eam][al][minimization]") {
  Parameters params;
  params.potential_options.potential = PotType::EAM_AL;
  params.optimizer_options.method = OptType::LBFGS;
  params.optimizer_options.converged_force = 0.01;
  params.optimizer_options.max_iterations = 50;
  auto pot = eonc::helpers::makePotential(params);
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
  params.potential_options.potential = PotType::EAM_AL;
  params.optimizer_options.method = OptType::FIRE;
  params.optimizer_options.converged_force = 0.01;
  params.optimizer_options.max_iterations = 200;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  matter->relax(false, false, false, "eam_fire_test", "eam_fire_test");

  // SVN reference (data/reference/minimization_eam_fire.dat):
  // energy = -5.562431, 26 force calls
  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-5.562431).epsilon(1e-4));
}

} /* namespace tests */
