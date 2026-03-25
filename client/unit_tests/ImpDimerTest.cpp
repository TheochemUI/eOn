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

#include "HelperFunctions.h"
#include "ImprovedDimer.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("ImprovedDimer computes eigenvalue on displaced cluster",
          "[ImprovedDimer]") {
  Parameters parameters;
  parameters.potential_options.potential = PotType::LJ;
  parameters.optimizer_options.method = OptType::CG;
  parameters.optimizer_options.converged_force = 0.001;
  parameters.dimer_options.converged_angle = 0.001;
  parameters.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_DIMER;

  auto pot = eonc::helpers::makePotential(parameters);
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

} /* namespace tests */
