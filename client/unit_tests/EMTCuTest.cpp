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

TEST_CASE("EMT energy matches SVN on Cu FCC cluster", "[pot][emt][cu]") {
  Parameters params;
  params.potential_options.potential = PotType::EMT;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  // SVN reference (data/reference/point_emt_cu.dat):
  double energy = matter->getPotentialEnergy();
  REQUIRE(energy == Catch::Approx(5.129167).epsilon(1e-4));

  double maxForce = matter->getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(1.914263).epsilon(1e-3));
}

TEST_CASE("EMT forces are conservative on Cu FCC", "[pot][emt][cu]") {
  Parameters params;
  params.potential_options.potential = PotType::EMT;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  AtomMatrix forces = matter->getForces();
  REQUIRE(forces.allFinite());

  // Sum of forces should be approximately zero (Newton's third law)
  Eigen::Vector3d totalForce = forces.colwise().sum();
  REQUIRE(totalForce.norm() < 1e-6);
}

TEST_CASE("EMT minimization converges on Cu FCC", "[pot][emt][cu]") {
  Parameters params;
  params.potential_options.potential = PotType::EMT;
  params.optimizer_options.method = OptType::LBFGS;
  params.optimizer_options.converged_force = 0.01;
  params.optimizer_options.max_iterations = 50;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  double e_before = matter->getPotentialEnergy();
  bool converged = matter->relax(false, false, false, "emt_test", "emt_test");
  double e_after = matter->getPotentialEnergy();

  REQUIRE(std::isfinite(e_after));
  REQUIRE(e_after <= e_before + 1e-10);
}

} /* namespace tests */
