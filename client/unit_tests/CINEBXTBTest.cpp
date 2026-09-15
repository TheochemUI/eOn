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

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// Regression test: CI-NEB with XTB on a small (9-atom) molecule.
// Reproduces a bug where removing EIGEN_DEFAULT_TO_ROW_MAJOR silently changed
// the storage order of bare MatrixXd types, corrupting force projections and
// causing the NEB to diverge from the first step (issue introduced in
// 6e8461c3).
TEST_CASE("CI-NEB XTB regression", "[neb][xtb]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::XTB;
  ParametersLoadAccess::xtb_options(params).paramset = "GFN2xTB";
  ParametersLoadAccess::xtb_options(params).acc = 1.0;
  ParametersLoadAccess::xtb_options(params).elec_temperature = 300.0;
  ParametersLoadAccess::xtb_options(params).maxiter = 250;

  ParametersLoadAccess::neb_options(params).image_count = 10;
  ParametersLoadAccess::neb_options(params).spring.weighting.enabled = true;
  ParametersLoadAccess::neb_options(params).spring.weighting.k_min = 0.972;
  ParametersLoadAccess::neb_options(params).spring.weighting.k_max = 9.72;
  ParametersLoadAccess::neb_options(params).spring.weighting.trigger = 0.5;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::LINEAR;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = true;
  ParametersLoadAccess::neb_options(params).climbing_image.converged_only =
      true;
  ParametersLoadAccess::neb_options(params).climbing_image.trigger_force = 0.5;
  ParametersLoadAccess::neb_options(params).climbing_image.trigger_factor = 0.8;
  ParametersLoadAccess::neb_options(params).force_tolerance = 0.0514221;
  ParametersLoadAccess::optimizer_options(params).method = OptType::LBFGS;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 100;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.1;

  auto pot = eonc::helpers::makePotential(params.potential_options().potential,
                                          params);
  auto initial = std::make_shared<Matter>(pot, params);
  auto final_state = std::make_shared<Matter>(pot, params);

  std::string reactFile("reactant.con");
  std::string prodFile("product.con");
  initial->con2matter(reactFile);
  final_state->con2matter(prodFile);

  auto neb =
      std::make_unique<NudgedElasticBand>(initial, final_state, params, pot);
  auto status = neb->compute();

  REQUIRE(static_cast<int>(status) ==
          static_cast<int>(NudgedElasticBand::NEBStatus::GOOD));

  neb->findExtrema();
  REQUIRE(neb->numExtrema >= 1);
}

} /* namespace tests */
