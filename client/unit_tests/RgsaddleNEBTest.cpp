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
#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("rgsaddle CI-NEB reaches the LJ13 saddle and counts every call",
          "[neb][rgsaddle][force_calls][climbing_image]") {
  // lbfgs and fire band steps run in the rgsaddle session in this build;
  // cg stays on the host.
  const OptType method = GENERATE(OptType::LBFGS, OptType::FIRE, OptType::CG);
  CAPTURE(static_cast<int>(method));
  Parameters params;
  ParametersLoadAccess::main_options(params).parallel = false;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::neb_options(params).opt_method = method;
  ParametersLoadAccess::neb_options(params).image_count = 7;
  ParametersLoadAccess::neb_options(params).max_iterations = 1000;
  ParametersLoadAccess::neb_options(params).force_tolerance = 0.01;
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = true;
  ParametersLoadAccess::neb_options(params).climbing_image.ocineb.use_mmf =
      false;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::LINEAR;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.1;
  auto pot = std::make_shared<eonc::helpers::test::CountingPotential>(
      eonc::helpers::sharePotential(
          eonc::helpers::makePotential(PotType::LJ, params)));
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  Matter saddle(pot, params);
  reactant->con2matter(std::string("reactant.con"));
  product->con2matter(std::string("product.con"));
  saddle.con2matter(std::string("saddle.con"));
  const double saddleEnergy = saddle.getPotentialEnergy();

  const size_t counterBefore = pot->forceCallCounter;
  const size_t kernelBefore = pot->kernelCalls;
  NudgedElasticBand neb(reactant, product, params, pot);
  const auto status = neb.compute();
  const size_t counted = pot->forceCallCounter - counterBefore;
  const size_t kernel = pot->kernelCalls - kernelBefore;
  CAPTURE(counted, kernel, neb.lastIteration(), saddleEnergy);
  REQUIRE(status == NudgedElasticBand::NEBStatus::GOOD);
  REQUIRE(kernel > 0);
  REQUIRE(counted == kernel);
  REQUIRE(neb.climbingImage > 0);
  REQUIRE(
      neb.path[static_cast<size_t>(neb.climbingImage)]->getPotentialEnergy() ==
      Catch::Approx(saddleEnergy).margin(1e-2));
}

} // namespace tests
