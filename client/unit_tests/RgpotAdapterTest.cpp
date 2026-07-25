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

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// Guards the RgpotAdapter + parameter mapping for every migrated arm: the
// same pinned numbers as the pre-migration in-tree kernels, evaluated
// through helpers::makePotential exactly as jobs do. Kernel numerics
// themselves are pinned upstream in rgpot's own test suite.

TEST_CASE("Adapter LJ matches the pinned reference", "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-39.965351).epsilon(1e-4));
  REQUIRE(matter->getForces().rowwise().norm().maxCoeff() ==
          Catch::Approx(0.004704).epsilon(1e-2));
}

TEST_CASE("Adapter LJCluster matches the pinned reference",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::LJCLUSTER;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-39.965379).epsilon(1e-4));
}

TEST_CASE("Adapter Morse evaluates and reports shared-safe",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(params);

  REQUIRE(pot->isThreadSafe());
  REQUIRE(pot->isSharedInstanceThreadSafe());
  REQUIRE_FALSE(pot->needsPerImageInstance());

  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
  REQUIRE(matter->getForces().allFinite());
}

TEST_CASE("Adapter ZBL maps config.ini cutoffs through",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::ZBL;
  params.zbl_options.cut_inner = 2.0;
  params.zbl_options.cut_global = 2.5;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
}

} /* namespace tests */
