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
#include "eon/NudgedElasticBand.h"
#include "eon/Potential.h"
#include <memory>

#include <cmath>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("rgsaddle band steps a short LJ path", "[neb][rgsaddle]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::neb_options(params).opt_method = OptType::XTSCI;
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).max_iterations = 3;
  ParametersLoadAccess::neb_options(params).force_tolerance = 0.01;
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
  ParametersLoadAccess::neb_options(params).climbing_image.ocineb.use_mmf =
      false;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::LINEAR;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  ParametersLoadAccess::optimizer_options(params).xtsci.method = "fire";

  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  reactant->con2matter(std::string("reactant.con"));
  product->con2matter(std::string("reactant.con"));
  auto pos = product->getPositions();
  pos(0, 0) += 0.5;
  product->setPositions(pos);

  NudgedElasticBand neb(reactant, product, params, pot);
  const auto status = neb.compute();
  REQUIRE(status != NudgedElasticBand::NEBStatus::INIT);
  for (const auto &image : neb.path) {
    REQUIRE(std::isfinite(image->getPotentialEnergy()));
  }
}

} // namespace tests
