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

TEST_CASE("FeHe potential returns finite energy on Fe BCC cluster",
          "[pot][fehe][fe]") {
  Parameters params;
  params.potential_options.potential = PotType::FEHE;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  // SVN reference (data/reference/point_fehe.dat):
  double energy = matter->getPotentialEnergy();
  REQUIRE(energy == Catch::Approx(-43.959774).epsilon(1e-4));

  double maxForce = matter->getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(1.245094).epsilon(1e-3));
}

} /* namespace tests */
